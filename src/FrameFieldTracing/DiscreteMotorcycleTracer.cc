#include "DiscreteMotorcycleTracer.hh"
#include <math.h>
#include <algorithm>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <iostream>

// The following functions are some Utilifty functions.
/// Conversion
///  from Eigen Vector to OpenMesh Vector
Eigen::Vector3d toEigen(OpenMesh::Vec3d vec)
{
    return Eigen::Vector3d(vec[0], vec[1], vec[2]);
}
/// Conversion from OpenMesh Vector to Eigen Vector
OpenMesh::Vec3d toOM(Eigen::Vector3d &vec)
{
    return OpenMesh::Vec3d(vec.x(), vec.y(), vec.z());
}
/// Returns the sign of a number.
template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
/// Modulo operation (differs to %, which is remainder.)
int modulo(int x, int N)
{
    return (x % N + N) % N;
}

///	Returns the vector between the two points of the halfedge.
OpenMesh::Vec3d halfedge_vector(OpenMesh::TriMesh& mesh, OpenMesh::SmartHalfedgeHandle h)
{
    return mesh.point(h.to()) - mesh.point(h.from());
}

/// Returns the angle between two vectors.
double vector_angle(OpenMesh::Vec3d h, OpenMesh::Vec3d u, OpenMesh::Vec3d n)
{
    h.normalize();
    u.normalize();
    n.normalize();
    double sign = sgn(OpenMesh::dot(OpenMesh::cross(h, u), n)) < 0. ? -1. : 1.;
    auto res = sign * acos(OpenMesh::dot(h, u) * (.9999999));
    return res;
}


double angle_average(double alpha, double beta){
    auto y = sin(alpha) + sin(beta);
    auto x = cos(alpha) + cos(beta);
    return atan2(y,x);
}
double difference_angle(double alpha,double beta){
    OpenMesh::Vec2d vec_alpha(cos(alpha), sin(alpha));
    OpenMesh::Vec2d vec_beta(cos(beta), sin(beta));
    return acos(vec_alpha.dot(vec_beta));
}

/// Returns the frame oriented to the direction of Travel.
/// Given the frame and a direction, the function returns a pair of two vectors where the first vector is facing the direction,
/// while the second vector is the next direction in ccw order.
template <class T, class U>
std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> get_frame_by_direction(OpenMesh::TriMesh& mesh, T u, T v, U face, int matching)
{
    switch (matching)
    {
    case 0:
        return std::make_pair(mesh.property(u, face), mesh.property(v, face));
        break;
    case 1:
        return std::make_pair(mesh.property(v, face), -mesh.property(u, face));
        break;
    case 2:
        return std::make_pair(-mesh.property(u, face), -mesh.property(v, face));
        break;
    case 3:
        return std::make_pair(-mesh.property(v, face), mesh.property(u, face));
        break;
    default:
        throw std::logic_error("Matching is not 0,1,2 or 3.");
        break;
    }
}

/// Returns the frame of a face corrected to a neighbouring face's orientation.
/// Similar to 'get__frame_by_direction()', but rotates the frame in cw direction.
std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> get_frame_for_matching(OpenMesh::TriMesh& mesh, OpenMesh::FPropHandleT<OpenMesh::Vec3d> u, OpenMesh::FPropHandleT<OpenMesh::Vec3d> v, OpenMesh::SmartFaceHandle face, int matching)
{
    switch (matching)
    {
    case 0:
        return std::make_pair(mesh.property(u, face), mesh.property(v, face));
        break;
    case 1:
        return std::make_pair(-mesh.property(v, face), mesh.property(u, face));
        break;
    case 2:
        return std::make_pair(-mesh.property(u, face), -mesh.property(v, face));
        break;
    case 3:
        return std::make_pair(mesh.property(v, face), -mesh.property(u, face));
        break;
    default:
        throw std::logic_error("Matching is not 0,1,2 or 3.");
        break;
    }
}

/// Constructor.
DiscreteMotorcycleTracer::DiscreteMotorcycleTracer(OpenMesh::TriMesh& mesh, double angle_threshold): mesh(mesh), total_angle_threshold(angle_threshold)
{
    if(!mesh.get_property_handle(u, "dir_u"))
        std::cerr << "dir_u not found." << std::endl;
    if(!mesh.get_property_handle(v, "dir_v"))
        std::cerr <<  "dir_v not found." << std::endl;
    if(!mesh.get_property_handle(ep_pjump, "DF_EP_pjump"))
        std::cerr << "pjump not found." << std::endl;
    if(!mesh.get_property_handle(vp_index, "DF_VP_idx"))
        std::cerr << "Field index not found." << std::endl;

    //Add some properties.
    mesh.add_property(matching, "FF_EP_matching");
    mesh.add_property(HP_motorcycle_trail_index, "edge");
    mesh.add_property(VP_motorcycle_trail_index, "visited");
    mesh.add_property(edgenormal, "edgenormal");
    mesh.add_property(travel_direction, "travel_direction");
    mesh.add_property(ep_total_travel_distance, "total_travel_distance");
    mesh.add_property(ep_travel_distance_edge, "edge_travel_distance");
    mesh.add_property(vp_crash_vertex, "vp_crash_vertex");
    mesh.add_property(vp_original,"vp_original");
    mesh.add_property(ep_original,"ep_original");

}

int DiscreteMotorcycleTracer::get_motorcycle_number(){
    return motorcycle_index;
}

/// Returns a vector of outgoing halfeges, sorted to start at a given halfedge.
std::vector<OpenMesh::SmartHalfedgeHandle> halfedge_vector_from_start(const OpenMesh::SmartVertexHandle v, const OpenMesh::SmartHalfedgeHandle &h)
{
    auto vec = v.outgoing_halfedges().to_vector();
    std::rotate(vec.begin(), std::find(vec.begin(), vec.end(), h), vec.end());
    return vec;
}

/// Returns whether or not there is a seperatix between two faces.
bool check_for_seperatrix(OpenMesh::Vec3d f0n, OpenMesh::Vec3d f1n, OpenMesh::Vec3d v0, OpenMesh::Vec3d v1, OpenMesh::Vec3d d0, OpenMesh::Vec3d d1)
{
    return 0 > v0.cross(d0).dot(f0n) * v1.cross(d1).dot(f1n);
}

/// Finds possible starting directions for a given halfedge by checking for seperatrices in different directions.
std::vector<int> DiscreteMotorcycleTracer::get_starting_directions(OpenMesh::SmartHalfedgeHandle eh, OpenMesh::SmartVertexHandle vertex)
{

    auto facehandle = eh.face();
    auto ofacehandle = eh.opp().face();
    const auto &vpoint = mesh.point(vertex);
    //find which edges should have motorcycles send down.
    OpenMesh::Vec3d b0 = mesh.calc_centroid(facehandle);
    OpenMesh::Vec3d b1 = mesh.calc_centroid(ofacehandle);
    OpenMesh::Vec3d v0 = b0 - vpoint;
    OpenMesh::Vec3d v1 = b1 - vpoint;
    int match = mesh.property(matching, eh.edge());

    v0.normalize();
    v1.normalize();

    auto d0 = get_frame_for_matching(mesh, u, v, facehandle, 0);
    auto d1 = get_frame_for_matching(mesh, u, v, ofacehandle, eh == eh.edge().h0() ? match : modulo(-match,4));

    auto f0n = d0.first.cross(d0.second).normalized();
    auto f1n = d1.first.cross(d1.second).normalized();

    auto left_bisector = (halfedge_vector(mesh,eh).normalized() + halfedge_vector(mesh,eh.prev().opp()).normalized()).normalized();
    auto right_bisector = (halfedge_vector(mesh,eh).normalized() + halfedge_vector(mesh,eh.opp().next()).normalized()).normalized();

    std::vector<int> ret;
    if (check_for_seperatrix(f0n, f1n, left_bisector, right_bisector, d0.first, d1.first))
    {
        ret.push_back(left_bisector.dot(d0.first.normalized()) + right_bisector.dot(d1.first.normalized())> 0. ? 0 : 2);
    }
    if (check_for_seperatrix(f0n, f1n, left_bisector, right_bisector, d0.second, d1.second))
    {
        ret.push_back(left_bisector.dot(d0.second.normalized()) + right_bisector.dot(d1.second.normalized())> 0. ? 1 : 3);
    }
    return ret;
}

/// Evaluates a vector of halfedges for their cost. Creates startedge candidates and pushes them to the vector.
void DiscreteMotorcycleTracer::get_startedge_candidates(int starting_direction, std::vector<OpenMesh::SmartHalfedgeHandle> halfedges, std::vector<TravelCandidate> &startedge_candidates, bool clockwise)
{
    //we have a halfedge and a direction.
    auto direction = starting_direction;
    for (size_t i = 1; i < halfedges.size(); i++)
    {

        direction = clockwise ? get_direction_in_neighbour(direction,halfedges[i-1]) : get_direction_in_neighbour(direction,halfedges[i-1].prev());
        if (/*abs(resulting_total_angle(0,0, direction, halfedges[i])) < total_angle_threshold  &&*/ !halfedges[i].is_boundary() && abs(get_angle_to_direction(halfedges[i],direction)) < M_PI_4)
        { //Passes test: Angle <= Threshhold.
            bool head_on = false;
            auto cost = get_cost_halfedge(halfedges[i], direction, 0.);

            if (check_for_head_on_collision(halfedges[i], direction))
            {
                cost = -INFINITY;
                head_on = true;
            }
            startedge_candidates.push_back(TravelCandidate(halfedges[i], direction, cost, head_on));
        }
        else
        {
            break;
        }
    }
}

DiscreteMotorcycleTracer::TravelCandidate::TravelCandidate(OpenMesh::SmartHalfedgeHandle he, int dir, double cost, int head_on) : halfedge(he), direction(dir), cost(cost), head_on(head_on){};

/// Refine the mesh around singularities to increase resolution.
/// Currently Simply 2-4 splits the outgoing halfedges in clockwise direction.
void DiscreteMotorcycleTracer::refine_around_singularites()
{
    std::vector<OpenMesh::SmartVertexHandle> singularities;
    for (auto vertex : mesh.vertices())
    {
        if (mesh.property(vp_index, vertex) != 0 || vertex.feature())
        {
            singularities.push_back(vertex);
        }
    }
    for (auto singularity : singularities)
    {

        auto hes = singularity.outgoing_halfedges().to_vector();

        for (auto halfedge : hes)
        {
            if(halfedge.is_boundary())
                split_edge(halfedge.opp(), 0.5);
            else
                split_edge(halfedge,0.5);
        }
        hes = singularity.outgoing_halfedges().to_vector();
        for(size_t i = 0; i < hes.size(); i++){
            if(!(hes[modulo(i-1,hes.size())].is_boundary() && hes[i].is_boundary())){
                if( acos(halfedge_vector(mesh,hes[modulo(i-1,hes.size())]).normalized().dot(halfedge_vector(mesh,hes[i]).normalized())) > M_PI_4)
                    split_edge(hes[i].next());
            }
        }
    }
    number_of_original_vertices = mesh.n_vertices();
}

//Returns the direction on the other side of the halfedge. (The face adjacent to he.opp())
int DiscreteMotorcycleTracer::get_direction_in_neighbour(int direction, OpenMesh::SmartHalfedgeHandle he){
    return he == he.edge().h0() ? ( direction + modulo(-mesh.property(matching, he.edge()),4) ) % 4 : ( direction + mesh.property(matching,he.edge())) % 4;
}


int DiscreteMotorcycleTracer::get_best_direction(OpenMesh::SmartHalfedgeHandle he){
    std::vector<std::pair<double,double>> direction_angles;
    for(int i= 0 ; i<=3;i++){
        direction_angles.push_back(std::make_pair(abs(get_angle_to_direction(he,i)),i));
    }
    std::sort(direction_angles.begin(),direction_angles.end());
    return direction_angles[0].second;
}

/// Starts motorcycles.
void DiscreteMotorcycleTracer::start_motorcycles()
{
    // First refine around singularites, specifically to increase the amount of discrete steps between close lying singularities.
    refine_around_singularites();

    motorcycle_index = 1;
    std::vector<DiscreteMotorcycle> boundary_cycles;
    std::vector<DiscreteMotorcycle> feature_cycles;
    std::vector<OpenMesh::SmartVertexHandle> singularities;
    for (auto vertex : mesh.vertices())
    {
        if (mesh.property(vp_index, vertex) != 0 || vertex.feature())
        {
            singularities.push_back(vertex);
        }
    }

    for (auto singularity : singularities)
    { //If vertex is a singularity
        int index = mesh.property(vp_index, singularity);
        std::vector<MotorcycleStart> possible_starts;
        if(!(singularity.is_boundary() && index == 1)){
            for (auto heh : singularity.outgoing_halfedges())
            {
                std::vector<int> starting_directions;
                auto halfedges = halfedge_vector_from_start(singularity, heh);
                auto reverse_halfedges = halfedges;
                std::reverse(reverse_halfedges.begin() + 1, reverse_halfedges.end());

                //Seperatrix might pass between barycenter of face next to boundary and boundary. The ordinary test doesnt catch that.
                if(heh.edge().is_boundary()){
                    starting_directions.empty();
                    auto best_direction = get_best_direction(heh);
                    auto facehandle = heh.is_boundary() ? heh.opp().face() : heh.face();
                    auto frame = std::make_pair(mesh.property(u,facehandle  ),mesh.property(v,facehandle));
                    OpenMesh::Vec3d b0 = mesh.calc_centroid(facehandle);
                    OpenMesh::Vec3d v0 = b0 - mesh.point(heh.from());
                    auto f0n = mesh.calc_face_normal(facehandle);

                    int boundary_orthogonal_start_direction;
                    bool found_seperatrix = false;

                    if(best_direction == 0 || best_direction == 2){
                        if(check_for_seperatrix(f0n,f0n,v0,halfedge_vector(mesh,heh),frame.second,frame.second)){
                            boundary_orthogonal_start_direction = (v0.dot(frame.second) > 0. ? 1 : 3);
                            found_seperatrix = true;
                        }
                    } else {
                        if(check_for_seperatrix(f0n,f0n,v0,halfedge_vector(mesh,heh),frame.first,frame.first)){
                            boundary_orthogonal_start_direction = (v0.dot(frame.first) > 0. ? 0 : 2);
                            found_seperatrix = true;
                        }
                    }
                    if(!found_seperatrix)
                        continue;
                    auto he = heh.is_boundary() ? heh.opp().next() : heh;
                    MotorcycleStart startedge_candidates(he, boundary_orthogonal_start_direction);
                    if(heh.is_boundary()){
                        get_startedge_candidates(boundary_orthogonal_start_direction, halfedges, startedge_candidates.edge_candidates, true);
                    }
                    else{
                        get_startedge_candidates(boundary_orthogonal_start_direction, reverse_halfedges, startedge_candidates.edge_candidates, false);
                    }
                    possible_starts.push_back(startedge_candidates);
                    continue;

                } else{
                    starting_directions = get_starting_directions(heh, singularity);
                }

                for (int starting_direction : starting_directions)
                {

                    MotorcycleStart startedge_candidates(heh, starting_direction);
                    //TODO: Put this somewhere else.
                    if (abs(resulting_total_angle(0,0,starting_direction, heh)) < total_angle_threshold && abs(get_angle_to_direction(heh,starting_direction)) < M_PI_2)
                    { //Passes test: Angle <= Threshhold.
                        auto cost = get_cost_halfedge(heh, starting_direction, 0.);
                        int head_on = check_for_head_on_collision(heh, starting_direction);
                        if (head_on)
                        {
                            cost = -INFINITY;
                        }
                        startedge_candidates.edge_candidates.push_back(TravelCandidate(heh, starting_direction, cost, head_on));
                    }
                    get_startedge_candidates(starting_direction, halfedges, startedge_candidates.edge_candidates, true);
                    get_startedge_candidates(starting_direction, reverse_halfedges, startedge_candidates.edge_candidates, false);
                    possible_starts.push_back(startedge_candidates);

                }


                //IF there is any without edge -> treat those first. Split.
                //now we have 3 or 5 maps with keys and costs. Solve via best matching.
            }
        }

        //Treat boundaries:
        for(auto he : singularity.outgoing_halfedges()){

            if(he.edge().is_boundary()){

                std::vector<std::pair<double,double>> direction_angles;
                for(int i= 0 ; i<=3;i++){
                    direction_angles.push_back(std::make_pair(abs(get_angle_to_direction(he,i)),i));
                }
                std::sort(direction_angles.begin(),direction_angles.end());
                int direction = direction_angles[0].second;
                auto cycle = DiscreteMotorcycle(motorcycle_index,direction,he,singularity);

                mesh.property(HP_motorcycle_trail_index, he) = motorcycle_index;
                mesh.property(travel_direction, he) = direction;
                auto travel = decompose_travel_direction_halfedge(he,direction);
                cycle.length = travel.first;
                cycle.current_vertex = he.to();

                boundary_cycles.push_back(cycle);

                mesh.property(ep_total_travel_distance,he.edge()) = travel.first;
                mesh.property(ep_travel_distance_edge,he.edge()) = travel.first;

                motorcycle_index++;
            }
        }

        if (possible_starts.size() != (size_t)(4 - mesh.property(vp_index, singularity)))
            //std::cerr << "Warning: Vertex " << singularity.idx() << "with Index: " << mesh.property(vp_index, singularity) << " has " << possible_startedges.size() << " starting edges.");
            std::sort(possible_starts.begin(), possible_starts.end(), [](MotorcycleStart &c0, MotorcycleStart &c1) { return c0.edge_candidates.size() < c1.edge_candidates.size(); });


        std::vector<MotorcycleStart> duplicate_free_starts;
        std::vector<std::pair<OpenMesh::SmartHalfedgeHandle, int>> duplicate_check_vec;

        for(auto start : possible_starts){
            //IF NOT DUPLIACTE OF ANY ON DUPLICATE FREE && NOT DUPLICATE OF A BOUNDARY
            bool is_duplicate = false;
            if(std::find(duplicate_check_vec.begin(),duplicate_check_vec.end(), std::make_pair(start.original_halfedge,start.original_direction)) != duplicate_check_vec.end()) {
                is_duplicate = true;
            }
            //ADD TO DUPLICATE FREE
            for(auto ec : start.edge_candidates){
                if(std::find(duplicate_check_vec.begin(),duplicate_check_vec.end(), std::make_pair(ec.halfedge,ec.direction)) != duplicate_check_vec.end()) {
                    is_duplicate = true;
                }
            }
            if(!is_duplicate){
                duplicate_free_starts.push_back(start);
                duplicate_check_vec.push_back( std::make_pair(start.original_halfedge,start.original_direction));
                std::transform(start.edge_candidates.begin(),start.edge_candidates.end(),std::back_inserter(duplicate_check_vec),[](TravelCandidate& c){return std::make_pair(c.halfedge,c.direction);});
            }
        }

        possible_starts = duplicate_free_starts;

        for (auto it = possible_starts.begin(); it != possible_starts.end(); it++)
        {
            auto candidate = *it;
            std::sort(candidate.edge_candidates.begin(), candidate.edge_candidates.end(), [](const TravelCandidate &c0, const TravelCandidate &c1) { return c0.cost < c1.cost; });
            if (candidate.edge_candidates.size() != 0)
            {
                auto best_pick = candidate.edge_candidates[0];
                assert(best_pick.halfedge.is_valid());
                auto cycle = DiscreteMotorcycle(motorcycle_index, best_pick.direction, best_pick.halfedge, singularity);

                mesh.property(HP_motorcycle_trail_index, best_pick.halfedge) = motorcycle_index;
                mesh.property(travel_direction, best_pick.halfedge) = best_pick.direction;
                auto travel = get_travel_distance_halfedge(best_pick.halfedge,best_pick.direction);
                mesh.property(ep_total_travel_distance,best_pick.halfedge.edge()) = travel;
                mesh.property(ep_travel_distance_edge,best_pick.halfedge.edge()) = travel;
                cycle.length = travel;
                if (best_pick.head_on)
                {
                    mesh.property(VP_motorcycle_trail_index, best_pick.halfedge.to()) = motorcycle_index;
                    head_on_crashes.push_back(std::make_pair(motorcycle_index, best_pick.head_on));
                }

                else
                {
                    motorcycles.push_back(cycle);
                }
                motorcycle_index++;
                for (auto erase_it = it + 1; erase_it != possible_starts.end(); erase_it++)
                {
                    erase_it->edge_candidates.erase(remove_if(erase_it->edge_candidates.begin(), erase_it->edge_candidates.end(), [best_pick](TravelCandidate &c) { return c.halfedge == best_pick.halfedge; }), erase_it->edge_candidates.end());
                }
            }
            else
            {
                //Try to find the best split within the two edges to the left and two to the right.

                double min_cost = INFINITY;
                double min_split;
                OpenMesh::SmartHalfedgeHandle min_edge;
                int min_direction = 0;

                int steps; //= candidate.original_halfedge.prev().edge().is_boundary() ? 3 : 4:
                OpenMesh::SmartHalfedgeHandle current_edge;
                int current_direction;

                if(candidate.original_halfedge.prev().edge().is_boundary()){
                    steps = 3;
                    current_edge = candidate.original_halfedge;
                    current_direction = candidate.original_direction;
                } else {
                    steps = 4;
                    current_edge = candidate.original_halfedge.prev().opp();
                    current_direction = get_direction_in_neighbour(candidate.original_direction,current_edge.opp());
                }

                double current_cost;
                double current_split;
                if(!candidate.original_halfedge.edge().is_boundary()){
                    for(int i=0;i<steps;i++){
                        if(evaluate_split(current_cost,current_edge,current_direction,0.,current_split)){
                            if(current_cost < min_cost)
                            {
                                min_cost = current_cost;
                                min_direction = current_direction;
                                min_edge = current_edge;
                                min_split = current_split;
                            }
                        }
                        if(current_edge.edge().is_boundary())
                            break;
                        current_direction = get_direction_in_neighbour(current_direction,current_edge);
                        current_edge = current_edge.opp().next();
                    }
                } else {
                    if(candidate.original_halfedge.is_boundary()){
                        min_edge = candidate.original_halfedge.opp().next();
                        evaluate_split(min_cost,min_edge,candidate.original_direction,0,min_split);
                        min_direction = candidate.original_direction;
                    } else {
                        min_edge = candidate.original_halfedge;
                        evaluate_split(min_cost,min_edge,candidate.original_direction,0,min_split);
                        min_direction = candidate.original_direction;
                    }
                }

                OpenMesh::SmartHalfedgeHandle next_edge;

                int direction = min_direction;
                split_edge(min_edge.next(), min_split);
                next_edge = min_edge.prev().opp();



                // Occupy next edge like a normal cycle.
                // next_edge, direction.

                auto cycle = DiscreteMotorcycle(motorcycle_index, direction, next_edge, singularity);
                mesh.property(HP_motorcycle_trail_index, next_edge) = motorcycle_index;
                mesh.property(travel_direction, next_edge) = direction;
                auto travel = get_travel_distance_halfedge(next_edge,direction);
                mesh.property(ep_total_travel_distance,next_edge.edge()) = travel;
                mesh.property(ep_travel_distance_edge,next_edge.edge()) = travel;
                cycle.length = travel;
                motorcycles.push_back(cycle);
                motorcycle_index++;
            }
        }
    }

    while(!boundary_cycles.empty()){
        auto next_cycle = std::min_element(boundary_cycles.begin(), boundary_cycles.end(),
                                           [](DiscreteMotorcycle &m0, DiscreteMotorcycle &m1) { return m0.GetLength() < m1.GetLength(); });
        auto vertex = next_cycle->current_vertex;
        OpenMesh::SmartHalfedgeHandle entry_halfedge;
        for(auto he : vertex.incoming_halfedges()){
            if(mesh.property(HP_motorcycle_trail_index, he) == next_cycle->index){
                entry_halfedge = he;
            }
        }

        //Crash, handle it.
        if(mesh.property(VP_motorcycle_trail_index,vertex) != 0){
            mesh.property(vp_crash_vertex,vertex) = true;
            boundary_cycles.erase(next_cycle);
            continue;
        }


        // If ingoing trail, then crash, add head on,
        int collisioncycle = 0;
        for(auto incoming : vertex.incoming_halfedges()){
            int edgetrail = mesh.property(HP_motorcycle_trail_index, incoming);
            if(edgetrail != 0 && edgetrail != next_cycle->index)
                collisioncycle = edgetrail;
        }

        if(collisioncycle){
            head_on_crashes.push_back(std::make_pair(next_cycle->index, collisioncycle));
            mesh.property(VP_motorcycle_trail_index, vertex) = next_cycle->index;
            boundary_cycles.erase(next_cycle);
            continue;
        }
        OpenMesh::SmartHalfedgeHandle next_edge;
        for(OpenMesh::TriMesh::VertexOHalfedgeIter he_iter = mesh.voh_iter(vertex); he_iter.is_valid(); he_iter++){

            if(he_iter->edge() != entry_halfedge.edge()){
                if(he_iter->edge().is_boundary()){
                    next_edge = *he_iter;
                }
            }
        }

        //Get direction: Lowest cost.
        std::vector<std::pair<double,double>> direction_angles;
        for(int i= 0 ; i<=3;i++){
            direction_angles.push_back(std::make_pair(abs(get_angle_to_direction(next_edge,i)),i));
        }
        std::sort(direction_angles.begin(),direction_angles.end());
        int direction = direction_angles[0].second;

        auto distance = get_travel_distance_halfedge(next_edge,direction);

        next_cycle->length += distance;
        next_cycle->direction = direction;
        next_cycle->current_vertex= next_edge.to();
        mesh.property(ep_travel_distance_edge, next_edge.edge()) = distance;
        mesh.property(ep_total_travel_distance, next_edge.edge()) = next_cycle->length;

        // Leave a trail for other cycles to crash in.
        mesh.property(HP_motorcycle_trail_index, next_edge) = next_cycle->index;
        mesh.property(VP_motorcycle_trail_index, vertex) = next_cycle->index;
        mesh.property(travel_direction, next_edge) = direction;




    }

    while(!feature_cycles.empty()){
        auto next_cycle = std::min_element(feature_cycles.begin(), feature_cycles.end(),
                                           [](DiscreteMotorcycle &m0, DiscreteMotorcycle &m1) { return m0.GetLength() < m1.GetLength(); });
        auto vertex = next_cycle->current_vertex;
        OpenMesh::SmartHalfedgeHandle entry_halfedge;
        for(auto he : vertex.incoming_halfedges()){
            if(mesh.property(HP_motorcycle_trail_index, he) == next_cycle->index){
                entry_halfedge = he;
            }
        }

        //Crash, handle it.
        if(mesh.property(VP_motorcycle_trail_index,vertex) != 0){
            mesh.property(vp_crash_vertex,vertex) = true;
            feature_cycles.erase(next_cycle);
            continue;
        }

        // If ingoing trail, then crash, add head on,
        int collisioncycle = 0;
        for(auto incoming : vertex.incoming_halfedges()){
            int edgetrail = mesh.property(HP_motorcycle_trail_index, incoming);
            if(edgetrail != 0 && edgetrail != next_cycle->index)
                collisioncycle = edgetrail;
        }

        if(collisioncycle){
            head_on_crashes.push_back(std::make_pair(next_cycle->index, collisioncycle));
            mesh.property(VP_motorcycle_trail_index, vertex) = next_cycle->index;
            feature_cycles.erase(next_cycle);
            continue;
        }

        OpenMesh::SmartHalfedgeHandle next_edge;
        for(auto outgoing : vertex.outgoing_halfedges()){

            if(outgoing.edge() != entry_halfedge.edge()){
                if(outgoing.feature()){
                    next_edge = outgoing;
                }
            }
        }
        //Get direction: Lowest cost.
        std::vector<std::pair<double,double>> direction_angles;
        for(int i= 0 ; i<=3;i++){
            direction_angles.push_back(std::make_pair(abs(get_angle_to_direction(next_edge,i)),i));
        }
        std::sort(direction_angles.begin(),direction_angles.end());
        int direction = direction_angles[0].second;

        auto distance = get_travel_distance_halfedge(next_edge,direction);

        next_cycle->length += distance;
        next_cycle->direction = direction;
        next_cycle->current_vertex = next_edge.to();
        mesh.property(ep_travel_distance_edge, next_edge.edge()) = distance;
        mesh.property(ep_total_travel_distance, next_edge.edge()) = next_cycle->length;

        // Leave a trail for other cycles to crash in.
        mesh.property(HP_motorcycle_trail_index, next_edge) = next_cycle->index;
        mesh.property(VP_motorcycle_trail_index, vertex) = next_cycle->index;
        mesh.property(travel_direction, next_edge) = direction;

    }


    //These get treated first.
};

/// For a Halfedge, get the angle to the given direction in reference to the edges first halfedges face.
double DiscreteMotorcycleTracer::get_angle_to_direction(OpenMesh::SmartHalfedgeHandle h, int dir)
{
    std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame0 = h.is_boundary() ? std::make_pair(mesh.property(u, h.opp().face()),mesh.property(v, h.opp().face())) : std::make_pair(mesh.property(u, h.face()),mesh.property(v, h.face()));
    std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame1 = h.opp().is_boundary() ?
                                                   std::make_pair(mesh.property(u, h.face()),mesh.property(v, h.face())) :
                                                   get_frame_for_matching(mesh,u,v,h.opp().face(), h == h.edge().h0() ? mesh.property(matching,h.edge()) : modulo(-mesh.property(matching,h.edge()),4));


    double res;
    double v0, v1;
    switch (dir)
    {
    case 0:

        v0 = vector_angle(halfedge_vector(mesh, h), frame0.first, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), frame1.first, frame1.first.cross(frame1.second));
        break;
    case 1:
        v0 = vector_angle(halfedge_vector(mesh, h), frame0.second, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), frame1.second, frame1.first.cross(frame1.second));        ;
        break;
    case 2:
        v0 = vector_angle(halfedge_vector(mesh, h), -frame0.first, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), -frame1.first, frame1.first.cross(frame1.second));
        break;
    case 3:
        v0 = vector_angle(halfedge_vector(mesh, h), -frame0.second, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), -frame1.second, frame1.first.cross(frame1.second));
        break;

    default:
        throw std::logic_error("Direction is not 0,1,2 or 3.");
        break;
    }

    res = angle_average(v0,v1);
    return res;
}


bool DiscreteMotorcycleTracer::is_imploding_or_exploding(OpenMesh::SmartHalfedgeHandle h, int dir)
{
    std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame0 = h.is_boundary() ? std::make_pair(mesh.property(u, h.opp().face()),mesh.property(v, h.opp().face())) : std::make_pair(mesh.property(u, h.face()),mesh.property(v, h.face()));
    std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame1 = h.opp().is_boundary() ?
                                                   std::make_pair(mesh.property(u, h.face()),mesh.property(v, h.face())) :
                                                   get_frame_for_matching(mesh,u,v,h.opp().face(), h == h.edge().h0() ? mesh.property(matching,h.edge()) : modulo(-mesh.property(matching,h.edge()),4));


    double res;
    double v0, v1;
    switch (dir)
    {
    case 0:

        v0 = vector_angle(halfedge_vector(mesh, h), frame0.first, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), frame1.first, frame1.first.cross(frame1.second));
        break;
    case 1:
        v0 = vector_angle(halfedge_vector(mesh, h), frame0.second, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), frame1.second, frame1.first.cross(frame1.second));        ;
        break;
    case 2:
        v0 = vector_angle(halfedge_vector(mesh, h), -frame0.first, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), -frame1.first, frame1.first.cross(frame1.second));
        break;
    case 3:
        v0 = vector_angle(halfedge_vector(mesh, h), -frame0.second, frame0.first.cross(frame0.second));
        v1 = vector_angle(halfedge_vector(mesh, h), -frame1.second, frame1.first.cross(frame1.second));
        break;

    default:
        throw std::logic_error("Direction is not 0,1,2 or 3.");
        break;
    }
    return v0 * v1 < 0.;
}

DiscreteMotorcycle::DiscreteMotorcycle(int idx, int dir, OpenMesh::SmartHalfedgeHandle& starting_edge, OpenMesh::SmartVertexHandle& start_vertex) : index(idx), direction(dir), start_vertex(start_vertex){
    current_vertex = starting_edge.to();
    drift = 0.;
    crash = false;
}

double DiscreteMotorcycleTracer::calculate_drift_halfedge(OpenMesh::SmartHalfedgeHandle h, int dir)
{
    return decompose_travel_direction_halfedge(h, dir).second;
}
/// Calculate distance traveled for a given frame and vector.
double DiscreteMotorcycleTracer::get_travel_distance(std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame, OpenMesh::Vec3d travel_vector, OpenMesh::Vec3d normal_vector)

{
    return decompose_travel_direction(frame, travel_vector, normal_vector).first;
}
double DiscreteMotorcycleTracer::get_travel_distance_halfedge(OpenMesh::SmartHalfedgeHandle h, int dir)

{
    return decompose_travel_direction_halfedge(h, dir).first;
}
std::pair<double, double> DiscreteMotorcycleTracer::decompose_travel_direction(std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame, OpenMesh::Vec3d travel_vector, OpenMesh::Vec3d normal_vector)
{
    auto b0 = frame.first.normalized();
    auto b1 = b0.cross(normal_vector);
    Eigen::Matrix2d f;
    f << frame.first.dot(b0), b0.dot(frame.second), 0., b1.dot(frame.second);
    auto he_vec = travel_vector;
    Eigen::Vector2d he_plane(b0.dot(he_vec), b1.dot(he_vec));
    auto finv = f.inverse();
    auto res = finv * he_plane;
    return std::make_pair(res[0], res[1]);
}
std::pair<double, double> DiscreteMotorcycleTracer::decompose_travel_direction_halfedge(OpenMesh::SmartHalfedgeHandle h, int dir)
{
    auto frame0 =  h.is_boundary() ? get_frame_by_direction(mesh,u,v,h.opp().face(),dir) : get_frame_by_direction(mesh,u,v,h.face(),dir);
    auto frame1 =  h.opp().is_boundary() ?  get_frame_by_direction(mesh,u,v,h.face(),dir)  : get_frame_by_direction(mesh,u,v,h.opp().face(),get_direction_in_neighbour(dir,h));
    auto travel_vector = halfedge_vector(mesh, h);
    auto normal0 = frame0.first.cross(frame0.second);
    auto normal1 = frame1.first.cross(frame1.second);
    auto decomp0 = decompose_travel_direction(frame0, travel_vector, normal0);
    auto decomp1 = decompose_travel_direction(frame1, travel_vector, normal1);

    return std::make_pair((decomp0.first + decomp1.first) * 0.5 , (decomp0.second + decomp1.second) * 0.5);
}

double DiscreteMotorcycleTracer::get_cost_halfedge(OpenMesh::SmartHalfedgeHandle h, int direction, double drift)
{
    auto res = decompose_travel_direction_halfedge(h,direction);
    return abs(res.second + drift) - res.first;
}


//Checks which direction around the vertex can be taken.
DiscreteMotorcycleTracer::CollisionAvoidance DiscreteMotorcycleTracer::check_for_parallel_collision(OpenMesh::SmartHalfedgeHandle h, int direction)
{
    auto vertex = h.to();
    auto halfedges = halfedge_vector_from_start(vertex, h.opp());
    auto originaldirection = direction;
    auto prev_angle = get_angle_to_direction(halfedges[0].opp(),direction);

    //Find two Halfedges that should not be crossed while iterating around the vertex to avoid flipping over due to traingles with a very large opening angle, where both edges have an angle < M_PI_2
    double highest = -M_PI;
    OpenMesh::SmartHalfedgeHandle highestHE;
    bool positive = true;

    double first_direction = get_direction_in_neighbour(direction,h);
    for (size_t i = 0; i < halfedges.size(); i++)
    {        //go through the edges, accounting for direction.
        auto incoming_halfedge = halfedges[i].opp();

        first_direction = get_direction_in_neighbour(first_direction,halfedges[i]);
        auto angle = get_angle_to_direction(incoming_halfedge, first_direction);
        if(abs(angle) > highest)
        {
            highestHE = incoming_halfedge;
            highest = abs(angle);
            if(angle > 0)
                positive = true;
            else
                positive = false;
        }
    }

    bool singularityOrBoundary = h.to().is_boundary() || mesh.property(vp_index,h.to()) != 0;

    for (size_t i = 1; i < halfedges.size(); i++)
    {        //go through the edges, accounting for direction.
        auto incoming_halfedge = halfedges[i].opp();

        direction = get_direction_in_neighbour(direction,halfedges[i]);
        auto angle = get_angle_to_direction(incoming_halfedge, direction);

        if(abs(angle) > M_PI_2 || halfedges[i-1].edge().is_boundary() || (!singularityOrBoundary && halfedges[i-((int)positive)].opp() == highestHE)){
            break;
        }
        prev_angle = angle;
        if(mesh.property(HP_motorcycle_trail_index, halfedges[i]) != 0 && (get_direction_in_neighbour(mesh.property(travel_direction, halfedges[i]),halfedges[i]) + direction) % 2 == 0){
            return CollisionAvoidance::RIGHT;
        }
        if(mesh.property(HP_motorcycle_trail_index, incoming_halfedge) != 0 && (direction + mesh.property(travel_direction, incoming_halfedge))%2 == 0){
            return CollisionAvoidance::RIGHT;
        }
    }
    auto reverse_halfedges = halfedges;
    std::reverse(reverse_halfedges.begin() + 1, reverse_halfedges.end());
    direction = originaldirection;
    prev_angle = get_angle_to_direction(halfedges[0].opp(), direction);
    for (size_t i = 1; i < reverse_halfedges.size(); i++)
    {
        //go through the edges, accounting for direction.
        direction = get_direction_in_neighbour(direction,reverse_halfedges[i-1].opp());
        auto incoming_halfedge = reverse_halfedges[i].opp();
        auto angle = get_angle_to_direction(incoming_halfedge, direction);
        if(abs(angle) > M_PI_2 || reverse_halfedges[i-1].edge().is_boundary() ||
            (!singularityOrBoundary && reverse_halfedges[i-(int)(!positive)].opp() == highestHE))
        {
            break;
        }
        prev_angle  =angle;
        if(mesh.property(HP_motorcycle_trail_index, reverse_halfedges[i]) != 0 && (get_direction_in_neighbour(mesh.property(travel_direction, reverse_halfedges[i]),reverse_halfedges[i]) + direction) % 2 == 0){
            return CollisionAvoidance::LEFT;
        }
        if(mesh.property(HP_motorcycle_trail_index, incoming_halfedge) != 0 && (direction + mesh.property(travel_direction, incoming_halfedge))%2 == 0){
            return CollisionAvoidance::LEFT;
        }
    }
    return CollisionAvoidance::BOTH;
}

int DiscreteMotorcycleTracer::check_for_head_on_collision(OpenMesh::SmartHalfedgeHandle h, int direction)
{
    int ret = 0;
    auto vertex = h.to();
    auto halfedges = halfedge_vector_from_start(vertex, h.opp());
    if (mesh.property(VP_motorcycle_trail_index, vertex) == 0)
    {
        for (size_t i = 1; i < halfedges.size(); i++)
        {
            direction = get_direction_in_neighbour(direction,halfedges[i]);
            if (mesh.property(HP_motorcycle_trail_index, halfedges[i].opp()) != 0)
            {
                int othertraveldirection = mesh.property(travel_direction, halfedges[i].opp());
                // Parallel traveling motorcycles
                if (((othertraveldirection + direction) % 2) == 0)
                {
                    ret = mesh.property(HP_motorcycle_trail_index, halfedges[i].opp());
                }
            }
        }
    }
    return ret;
}


std::vector<DiscreteMotorcycleTracer::TravelCandidate>
    DiscreteMotorcycleTracer::
    get_tracing_candidates(OpenMesh::SmartVertexHandle vertex, OpenMesh::SmartHalfedgeHandle entry_halfedge, DiscreteMotorcycle& cycle, int direction, double local_angle_treshold)
{
    std::vector<DiscreteMotorcycleTracer::TravelCandidate> candidates;

    //Get outgoing halfedges.
    auto halfedges = halfedge_vector_from_start(vertex, entry_halfedge.opp());

    //This iterates around the vertex in clockwise diretion starting at the halfedge that was entered.
    double existing_total_angle = atan((cycle.drift)/(cycle.length));
    double allowed_total_angle = std::max(existing_total_angle * 1.001,total_angle_threshold);

    int prev_direction = get_direction_in_neighbour(direction,halfedges[0].opp());
    for (size_t j = 1; j <= halfedges.size(); j++)
    {
        //Correct direction for new halfedge.
        int i =  j%halfedges.size();

        if(j>1)
        {
            prev_direction = direction;
            direction = get_direction_in_neighbour(direction,halfedges[j-1]);
        }

        auto ang0 = get_angle_to_direction(halfedges[i], direction);
        auto ang1 = get_angle_to_direction(halfedges[j-1], prev_direction);


        CollisionAvoidance evade;
        if (i==0){
            evade = CollisionAvoidance::LEFT;
        }
        else if(mesh.property(HP_motorcycle_trail_index, halfedges[i].opp()) != 0){
            int otherdirection = get_direction_in_neighbour(mesh.property(travel_direction, halfedges[i].opp()),halfedges[i].opp());
            evade = (direction + 1) % 4 == otherdirection ? CollisionAvoidance::LEFT : CollisionAvoidance::RIGHT;
        }else {
            evade = check_for_parallel_collision(halfedges[i],direction);
        }

        if((abs(ang0) > local_angle_treshold && abs(ang1) >  local_angle_treshold && !(sgn(ang0) > 0. &&  sgn(ang1) < 0)) || (abs(ang0) > M_PI * .7 && abs(ang1) > M_PI * .7)){
            if(abs(ang0) < M_PI_2){
                if(evade == CollisionAvoidance::LEFT){
                    break;
                }
                if(evade == CollisionAvoidance::RIGHT){
                    candidates.clear();
                }
            }
            continue;
        }

        if(evade != CollisionAvoidance::RIGHT){
            double splitcost;
            double splitpoint;
            if (evaluate_split(splitcost, halfedges[i], direction, cycle.drift, splitpoint))
            {
                auto h  =halfedges[i];
                auto candidate = TravelCandidate(halfedges[i], direction, splitcost, false);
                candidate.cost = splitcost;
                candidate.type = TravelCandidate::TravelType::SPLIT;
                candidate.splitparam = splitpoint;
                candidates.push_back(candidate);
            } else {
            }

        }

        if(evade == CollisionAvoidance::BOTH && !mesh.property(vp_index, halfedges[i].to()))
        {
            //Check wether halfedge is okay for traveling
            auto tmp = decompose_travel_direction_halfedge(halfedges[i], direction);

            if ((abs(resulting_total_angle(cycle.length, cycle.drift, direction,halfedges[i])) < allowed_total_angle && abs(get_angle_to_direction(halfedges[i],direction)) < local_angle_treshold) || is_imploding_or_exploding(halfedges[i],direction))
            { //Passes test: Angle <= Threshhold.c

                int headon = check_for_head_on_collision(halfedges[i], direction);
                auto cost = get_cost_halfedge(halfedges[i], direction, cycle.drift);
                //If there is a possible head on collision, provoke it.
                if (headon)
                {
                    cost = -INFINITY;
                }
                auto candidate = TravelCandidate(halfedges[i], direction, cost, headon);
                candidate.type = TravelCandidate::TravelType::EGDE;
                candidates.push_back(candidate);

            }

        }
        if(evade == CollisionAvoidance::LEFT){
            break;
        }
        if(evade == CollisionAvoidance::RIGHT){
            candidates.clear();

        }
    }

    return candidates;
}

//trace the motorcycles until no more cycles exist
void DiscreteMotorcycleTracer::trace_motorcycles()
{
    int current_vertex;
    auto active_motorcycles = motorcycles;
    std::vector<DiscreteMotorcycle> inactive_motorcycles;
    std::vector<DiscreteMotorcycle> non_crash_cycles;
    int hoc = 0;
    int hoc_ac = 0;
    int count = 0;
    int singularity_crash = 0;
    size_t crash = 0;
    int splits = 0;
    bool go_again = false;
    std::vector<DiscreteMotorcycle>::iterator prev_cycle;

    try
    {
        while (!active_motorcycles.empty())
        {
            //Find out which cycle to start with.

            auto next_cycle = go_again ? prev_cycle : std::min_element(active_motorcycles.begin(), active_motorcycles.end(), [](DiscreteMotorcycle &m0, DiscreteMotorcycle &m1) { return m0.GetLength() < m1.GetLength(); });
            go_again = false;

            auto vertex = next_cycle->current_vertex;
            OpenMesh::SmartHalfedgeHandle entry_halfedge;
            for(auto he : vertex.incoming_halfedges()){
                if(mesh.property(HP_motorcycle_trail_index, he) == next_cycle->index){
                    entry_halfedge = he;
                }
            }

            current_vertex = vertex.idx();
            //preliminary collision detection.
            if (mesh.property(vp_index, vertex) != 0)
            {
                // Crash into singularity, should not happen.
                inactive_motorcycles.push_back(*next_cycle);
                non_crash_cycles.push_back(*next_cycle);
                active_motorcycles.erase(next_cycle);
                singularity_crash++;
                continue;
            }
            if (mesh.property(VP_motorcycle_trail_index, vertex))
            {
                //Crashed into a trail.
                mesh.property(vp_crash_vertex,vertex) = true;
                next_cycle->crash =true;
                inactive_motorcycles.push_back(*next_cycle);
                active_motorcycles.erase(next_cycle);
                crash++;
                continue;
            }
            //Get the direction to travel down.
            int direction = next_cycle->direction;
            std::vector<TravelCandidate> candidates = get_tracing_candidates(vertex, entry_halfedge, *next_cycle,direction, M_PI_4);
            if(candidates.empty()){
                candidates = get_tracing_candidates(vertex, entry_halfedge, *next_cycle,direction, M_PI_2);
            }
            //choose edge
            if (candidates.empty())
            {
                non_crash_cycles.push_back(*next_cycle);
                active_motorcycles.erase(next_cycle);
                OpenMesh::EPropHandleT<int> edgemarker;
                mesh.add_property(edgemarker, "trailedge");

                for(auto edge : mesh.edges())
                    mesh.property(edgemarker, edge) = std::max(mesh.property(HP_motorcycle_trail_index,edge.h0()),mesh.property(HP_motorcycle_trail_index,edge.h1()));

            }
            //Sort by Travel type, and then cost. Such that an edge is always prefered.
            std::sort(candidates.begin(), candidates.end(), [](const TravelCandidate &c0, const TravelCandidate &c1) { return (c0.type == c1.type)? c0.cost < c1.cost : c0.type == TravelCandidate::TravelType::EGDE; });

            //If there is no edge now, after the splits, then somethings horribly wrong and breaking is okay.
            auto best_candidate = candidates[0];
            if(best_candidate.type == TravelCandidate::TravelType::SPLIT){
                //Split the edge and set the new edge as travel edge.
                split_edge(best_candidate.halfedge.next(),best_candidate.splitparam);
                best_candidate.halfedge = best_candidate.halfedge.prev().opp();
            }




            // update cycle.
            next_cycle->direction = best_candidate.direction;
            next_cycle->current_vertex = best_candidate.halfedge.to();
            auto traveldistance = get_travel_distance_halfedge(best_candidate.halfedge, best_candidate.direction);
            mesh.property(ep_travel_distance_edge, best_candidate.halfedge.edge()) = traveldistance;
            next_cycle->length += traveldistance;
            mesh.property(ep_total_travel_distance, best_candidate.halfedge.edge()) = next_cycle->length;
            next_cycle->drift += calculate_drift_halfedge(best_candidate.halfedge, best_candidate.direction);

            //Not an original vertex, go again until an original vertex is reached.
            if(!mesh.property(vp_original,vertex)){
                prev_cycle = next_cycle;
                go_again = false;
            }

            // Leave a trail for other cycles to crash in.
            mesh.property(HP_motorcycle_trail_index, best_candidate.halfedge) = next_cycle->index;
            mesh.property(VP_motorcycle_trail_index, vertex) = next_cycle->index;
            mesh.property(travel_direction, best_candidate.halfedge) = best_candidate.direction;
            if (best_candidate.head_on)
            {
                //Check if this was a head on collision. If yes, note it down, the information will be used to set corners properly when creating the tmesh.
                mesh.property(VP_motorcycle_trail_index, best_candidate.halfedge.to()) = next_cycle->index;
                head_on_crashes.push_back(std::make_pair(next_cycle->index, best_candidate.head_on));
                hoc_ac++;
            }

            count++;
        }
    }
    catch (std::exception &e)
    {
        //std::cerr << "Exception in vertex " << current_vertex << ": " << e.what());
    }

}

void DiscreteMotorcycleTracer::update_cycle(DiscreteMotorcycle &cycle, int direction, OpenMesh::SmartHalfedgeHandle &halfedge)
{
    // update cycle.
    cycle.direction = direction;
    cycle.current_vertex = halfedge.to();
    cycle.length += get_travel_distance_halfedge(halfedge, direction);
    cycle.drift += calculate_drift_halfedge(halfedge, direction);

    mesh.property(HP_motorcycle_trail_index, halfedge) = 1; //next_cycle->index;
    mesh.property(VP_motorcycle_trail_index, halfedge.to()) = cycle.index;
}

std::vector<std::pair<int,int>> DiscreteMotorcycleTracer::get_head_on_collisions()
{
    return head_on_crashes;
}

void DiscreteMotorcycleTracer::generate()
{
    for (auto edge : mesh.edges())
    {
        mesh.property(matching, edge) = modulo(mesh.property(ep_pjump, edge),4);
    }
    start_motorcycles();
    for (auto edge : mesh.edges())
    {
        mesh.property(ep_original,edge) = true;
    }
    for (auto vert : mesh.vertices()){
        mesh.property(vp_original, vert) = true;
    }
    trace_motorcycles();

}

double DiscreteMotorcycleTracer::resulting_total_angle(double length, double drift, int direction, OpenMesh::SmartHalfedgeHandle h){
    auto travel =  decompose_travel_direction_halfedge(h, direction);
    return atan((drift + travel.second)/(length   + travel.first));
}

bool DiscreteMotorcycleTracer::evaluate_split(double &cost, OpenMesh::SmartHalfedgeHandle h, int direction, double drift, double &splitpoint)
{
    int direction_on_this_face = direction;// h == h.edge().h0() ? direction : modulo(direction - mesh.property(matching, h.edge()), 4);
    auto frame = get_frame_by_direction(mesh, u, v, h.face(), direction_on_this_face);
    auto normal = frame.first.cross(frame.second);


    //If both halfedge have to much angle to the desired direction, then there can be no valid split on the edge between them.
    auto he0vec = halfedge_vector(mesh, h).normalized();
    auto he1vec = halfedge_vector(mesh, h.prev().opp()).normalized();
    auto ang0 = vector_angle(he0vec, frame.first, normal);
    auto ang1 = vector_angle(he1vec, frame.first, normal);

    if (abs(ang0) > M_PI_2 && abs(ang1) > M_PI_2 && !(sgn(ang0) > 0. &&  sgn(ang1) < 0))
    {
        return false;
    }



    //Valid edge, now get the best cost.
    //There is two cases: Integral Line lies ON the triangle: SPLIT THAT WAY
    //Integral Line DOES NOT lie on the triangle: Best split is as close as possible to the closer edge.

    auto pt0 = mesh.point(h.to());
    auto pt1 = mesh.point(h.next().to());


    auto v = mesh.point(h.from());

    auto pt0_decomposed = decompose_travel_direction(frame,(pt0 - v),normal);
    auto pt1_decomposed = decompose_travel_direction(frame,(pt1 - v),normal);
    auto pt0dabs = abs(pt0_decomposed.first);
    auto pt1dabs = abs(pt1_decomposed.first);

    if(pt0_decomposed.first < 0){
        pt0 = (pt0*(pt1dabs/(pt0dabs+pt1dabs))) +  (pt1 * (pt0dabs/(pt0dabs+pt1dabs)));//set new pt0, set pt0_decomposed.
        pt0_decomposed = decompose_travel_direction(frame,(pt0 - v),normal);

    }
    if(pt1_decomposed.first < 0){
        //set new pt1, set pt1_decomposed.
        pt1 = (pt0*(pt1dabs/(pt0dabs+pt1dabs))) +  (pt1 * (pt0dabs/(pt0dabs+pt1dabs)));//set new pt0, set pt0_decomposed.
        pt1_decomposed = decompose_travel_direction(frame,(pt1 - v),normal);
    }

    auto d0 = pt0_decomposed.second;//((pt0 - v) % frame.first.normalized()).length();
    auto d1 = pt1_decomposed.second;//((pt1 - v) % frame.first.normalized()).length();
    auto d0_abs = abs(d0);
    auto d1_abs = abs(d1);
    OpenMesh::Vec3d pt;
    pt = (d1_abs / (d0_abs + d1_abs)) * pt0 + (d0_abs / (d0_abs + d1_abs)) * pt1;
    if (d0 * d1 < 0.)
    { // case 1: Integral line lies in triangle.
        pt = (d1_abs / (d0_abs + d1_abs)) * pt0 + (d0_abs / (d0_abs + d1_abs)) * pt1;
        auto vec = pt - mesh.point(h.from());
        auto traveldirections = decompose_travel_direction(frame,vec,normal);
        cost = -INFINITY;
    }
    else
    {
        //Swap so that pt0 is always closer.
        if ( (2. * d0_abs) - pt0_decomposed.first > (2. * d1_abs) - pt1_decomposed.first)
        {
            std::swap(pt0, pt1);
        }
        pt = .99 * pt0 + .01 * pt1;
        auto vec = pt - mesh.point(h.from());
        auto traveldirections = decompose_travel_direction(frame,vec,normal);
        cost = abs(vector_angle(vec,frame.first,normal));//(abs(drift + traveldirections.second) * 2.) -traveldirections.first;
    }
    auto dist_pt0 = (pt-mesh.point(h.to())).length();
    auto dist_pt1 = (pt-mesh.point(h.next().to())).length();
    splitpoint = dist_pt0 / (dist_pt0 + dist_pt1);

    return true;
}
void DiscreteMotorcycleTracer::split_edge(OpenMesh::SmartHalfedgeHandle h)
{
    auto pt = .5;
    split_edge(h, pt);
}

void DiscreteMotorcycleTracer::split_edge(OpenMesh::SmartHalfedgeHandle h, double splitpoint)
{
    if (h.is_boundary())
    {
        split_edge(h.opp(), 1.0 - splitpoint);
        return;
    }

    bool boundary = h.edge().is_boundary();
    bool original = mesh.property(ep_original, h.edge());
    auto to_vh = h.to();
    auto from_vh = h.from();

    auto leftu = mesh.property(u, h.face());
    auto leftv = mesh.property(v, h.face());
    auto rightu = boundary ? OpenMesh::Vec3d(0,0,0) : mesh.property(u, h.opp().face());
    auto rightv = boundary ? OpenMesh::Vec3d(0,0,0) : mesh.property(v, h.opp().face());
    auto lefttop = h.next();
    auto leftbottom = lefttop.next();
    auto rightbottom = h.opp().next();
    auto righttop = rightbottom.next();


    auto edge_distance = mesh.property(ep_travel_distance_edge,h.edge());
    auto total_distance = mesh.property(ep_total_travel_distance, h.edge());

    int left_trail = mesh.property(HP_motorcycle_trail_index, h);
    int right_trail = mesh.property(HP_motorcycle_trail_index, h.opp());
    int left_direction = mesh.property(travel_direction, h);
    int right_direction = mesh.property(travel_direction, h.opp());


    bool going_left = left_trail > 0;

    int prev_matching = mesh.property(matching, h.edge());
    bool h_first = h == h.edge().h0();
    bool feature = h.edge().feature();

    auto vertexProp = std::max(mesh.property(HP_motorcycle_trail_index, h), mesh.property(HP_motorcycle_trail_index, h.opp()));

    const OpenMesh::Vec3d p_new = splitpoint * mesh.point(h.to()) + (1. - splitpoint) * mesh.point(h.from());
    const auto pjump_before = mesh.property(ep_pjump, h.edge());
    auto res = mesh.split(h.edge(), p_new);
    const auto heh_from = mesh.find_halfedge(from_vh, res);
    const auto heh_to = mesh.find_halfedge(res, to_vh);
    mesh.property(ep_pjump, heh_from.edge()) = pjump_before;
    mesh.property(ep_pjump, heh_to.edge()) = pjump_before;

    mesh.property(HP_motorcycle_trail_index, leftbottom.next()) = left_trail;
    mesh.property(HP_motorcycle_trail_index, lefttop.prev()) = left_trail;
    mesh.property(HP_motorcycle_trail_index, leftbottom.next().opp()) = right_trail;
    mesh.property(HP_motorcycle_trail_index, lefttop.prev().opp()) = right_trail;

    mesh.property(travel_direction, leftbottom.next()) = left_direction;
    mesh.property(travel_direction, lefttop.prev()) = left_direction;
    mesh.property(travel_direction, leftbottom.next().opp()) = right_direction;
    mesh.property(travel_direction, lefttop.prev().opp()) = right_direction;

    mesh.property(VP_motorcycle_trail_index, res) = vertexProp;
    mesh.property(u, lefttop.face()) = leftu;
    mesh.property(v, lefttop.face()) = leftv;
    mesh.property(u, leftbottom.face()) = leftu;
    mesh.property(v, leftbottom.face()) = leftv;

    if(!boundary){
        mesh.property(u, righttop.face()) = rightu;
        mesh.property(v, righttop.face()) = rightv;
        mesh.property(u, rightbottom.face()) = rightu;
        mesh.property(v, rightbottom.face()) = rightv;
    }

    auto newoutgoing = leftbottom.next();
    auto otheroutgoing = lefttop.prev();
    bool otheroutgoing_first = otheroutgoing.edge().h0() == otheroutgoing;
    bool newoutgoing_first = newoutgoing.edge().h0() == newoutgoing;

    double left_length = halfedge_vector(mesh,lefttop.prev()).length();
    double right_length = halfedge_vector(mesh,leftbottom.next().opp()).length();
    auto left_distance = edge_distance * (left_length / (left_length +right_length));
    auto right_distance = edge_distance * (right_length / (left_length +right_length));

    if(going_left){
        mesh.property(ep_total_travel_distance, leftbottom.next().edge()) = total_distance - left_distance;
        mesh.property(ep_travel_distance_edge, leftbottom.next().edge()) = right_distance;
        mesh.property(ep_total_travel_distance, lefttop.prev().edge()) = total_distance;
        mesh.property(ep_travel_distance_edge, lefttop.prev().edge()) = left_distance;
    } else {
        mesh.property(ep_total_travel_distance, leftbottom.next().edge()) = total_distance;
        mesh.property(ep_travel_distance_edge, leftbottom.next().edge()) = right_distance;
        mesh.property(ep_total_travel_distance, lefttop.prev().edge()) = total_distance - right_distance;
        mesh.property(ep_travel_distance_edge, lefttop.prev().edge()) = left_distance;
    }

    mesh.property(matching, newoutgoing.edge()) = newoutgoing_first == h_first ? prev_matching : modulo(-prev_matching, 4);
    mesh.property(matching, otheroutgoing.edge()) = otheroutgoing_first == h_first ? prev_matching : modulo(-prev_matching, 4);

    auto splithalfedge0 = mesh.find_halfedge(res, to_vh);
    auto splithalfedge1 = mesh.find_halfedge(from_vh, res);
    for(auto e : res.edges()){
        mesh.status(e).set_feature(false);
        mesh.property(ep_original,e) = false;
    }
    mesh.status(splithalfedge0.edge()).set_feature(feature);
    mesh.status(splithalfedge1.edge()).set_feature(feature);
    mesh.property(ep_original,splithalfedge0.edge()) = original;
    mesh.property(ep_original,splithalfedge1.edge()) = original;
    mesh.property(vp_original, res) = original;
    mesh.property(matching, lefttop.next().edge()) = 0;


    if(!boundary){
        mesh.property(matching, rightbottom.next().edge()) = 0;
    }
}
