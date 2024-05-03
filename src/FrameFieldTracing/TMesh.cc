#include "TMesh.hh"
#include <FrameFieldTracing/DiscreteMotorcycleTracer.hh>

enum BoundaryState{
    NONE,
    H0,
    H1
};

struct PotentialArc
{
    PotentialArc(bool it, bool nc, int m, double td, BoundaryState boundary, bool feature, OpenMesh::EdgeHandle eh)
        : incoming_trace(it), is_next_corner(nc), motorcycle(m), travel_distance(td), boundary_state(boundary), is_feature(feature), eh(eh) { }
    PotentialArc(bool it, bool nc, int m, double td, BoundaryState boundary, bool feature, OpenMesh::EdgeHandle eh, pm::edge_handle ah)
        : incoming_trace(it), is_next_corner(nc), motorcycle(m), travel_distance(td), boundary_state(boundary), is_feature(feature), eh(eh), arc_handle(ah) { }

    BoundaryState boundary_state;
    bool incoming_trace;
    bool is_next_corner;
    int motorcycle;
    double travel_distance;
    bool is_feature;

    OpenMesh::EdgeHandle eh;
    pm::edge_handle arc_handle;
};

struct PotentialNode
{
    std::vector<PotentialArc> potential_arcs;
    pm::vertex_handle node_handle;
    bool init = false;

    void add_arc(const PotentialArc& pa)
    {
        potential_arcs.push_back(pa);
    }

    std::vector<PotentialArc>::iterator get_arc(int motorcycle, bool incoming)
    {
        auto relevant_arc = [motorcycle, incoming](const PotentialArc &a) {
            return a.motorcycle == motorcycle && a.incoming_trace == incoming;
        };
        return find_if(potential_arcs.begin(), potential_arcs.end(), relevant_arc);
    }
};

std::vector<OpenMesh::VertexHandle> list_vh_between_nodes(int cycle_index, OpenMesh::VertexHandle start_vh, OpenMesh::VertexHandle end_vh, OpenMesh::TriMesh &tri_mesh_, OpenMesh::HPropHandleT<int>& prop)
{
    std::vector<OpenMesh::VertexHandle> results;
    OpenMesh::VertexHandle current_vh = start_vh;
    bool found = false;

    while(current_vh != end_vh){
        for(auto outgoing_heh = tri_mesh_.voh_begin(current_vh); outgoing_heh != tri_mesh_.voh_end(current_vh); outgoing_heh++){
            if(tri_mesh_.property(prop, *outgoing_heh) == cycle_index)
            {
                current_vh = tri_mesh_.to_vertex_handle(*outgoing_heh);
                break;
            }
        }
        if(current_vh != end_vh)
            results.push_back(current_vh);
    }
    return results;
};

void compute_t_mesh(OpenMesh::TriMesh _mesh, pm::Mesh& _tmesh, pm::vertex_attribute<tg::pos3>& _tmesh_pos)
{
    OpenMesh::VPropHandleT<int> vp_index;

    _mesh.get_property_handle(vp_index, "DF_VP_idx");

    _tmesh_pos = _tmesh.vertices().make_attribute<tg::pos3>();
    DiscreteMotorcycleTracer discreteTracer(_mesh, 0.02);
    discreteTracer.generate();

    auto head_on_collisions = discreteTracer.get_head_on_collisions();
    size_t number_motorcycles = discreteTracer.get_motorcycle_number();

    //Collect all intersections of traces:
    //Singularities, Boundary Starts, Featureedge Starts, and of course Collisions. Create tmeshnode for all of them.
    OpenMesh::VPropHandleT<bool> vp_collision;
    OpenMesh::HPropHandleT<int> hp_trail;
    OpenMesh::EPropHandleT<double> ep_travel_distance;
    OpenMesh::EPropHandleT<double> ep_edge_travel_distance;

    _mesh.get_property_handle(hp_trail, "edge");
    _mesh.get_property_handle(vp_collision, "vp_crash_vertex");
    _mesh.get_property_handle(ep_travel_distance, "total_travel_distance");
    _mesh.get_property_handle(ep_edge_travel_distance, "edge_travel_distance");

    OpenMesh::EPropHandleT<int> edgemarker;
    _mesh.add_property(edgemarker, "trailedge");

    for(auto edge : _mesh.edges())
        _mesh.property(edgemarker, edge) = std::max(_mesh.property(hp_trail,edge.h0()), _mesh.property(hp_trail,edge.h1()));

    _tmesh_pos = _tmesh.vertices().make_attribute<tg::pos3>();
    pm::vertex_attribute<OpenMesh::VertexHandle> original_vertex(_tmesh);

    std::vector<PotentialNode *> potential_nodes;
    // iterate over vertices
    for (auto vertex : _mesh.vertices())
    {
        // find singularities
        if (_mesh.property(vp_index, vertex) != 0 || vertex.feature())
        {
            //This is a singularity
            PotentialNode* pn = new PotentialNode();
            for (auto outgoing_he : vertex.outgoing_halfedges())
            {
                int cycle_val = _mesh.property(hp_trail, outgoing_he);
                bool is_corner = outgoing_he.is_boundary() || !outgoing_he.edge().is_boundary();
                BoundaryState boundary = outgoing_he.edge().is_boundary() ? (outgoing_he.is_boundary() ?  BoundaryState::H0 : BoundaryState::H1): BoundaryState::NONE;
                if (cycle_val != 0)
                {
                    pn->add_arc(PotentialArc(false, is_corner, cycle_val, 0, boundary, outgoing_he.feature(), outgoing_he.edge()));
                }
            }

            const auto nh = _tmesh.vertices().add();
            pn->node_handle = nh;
            _tmesh_pos[nh] = tg::pos3(_mesh.point(vertex));
            original_vertex[nh] = vertex;

            potential_nodes.push_back(pn);
        }
        else if (_mesh.property(vp_collision, vertex))
        {
            //This is a collision vertex. Create a node for it.
            PotentialNode* pn = new PotentialNode();
            for (OpenMesh::TriMesh::VertexIHalfedgeIter he_iter = _mesh.vih_iter(vertex); he_iter.is_valid(); he_iter++)
            {
                auto incoming_he = *he_iter;
                //For each halfedge that has a motorcycle trail, create an appropriate potential_arc.
                if (_mesh.property(hp_trail, incoming_he))
                {
                    pn->add_arc(PotentialArc(true, false, _mesh.property(hp_trail, incoming_he), _mesh.property(ep_travel_distance, incoming_he.edge()), BoundaryState::NONE, incoming_he.feature(), incoming_he.edge()));
                }
                else if (_mesh.property(hp_trail, incoming_he.opp()))
                {
                    BoundaryState boundary = incoming_he.opp().edge().is_boundary() ? (incoming_he.opp().is_boundary() ?  BoundaryState::H0 : BoundaryState::H1): BoundaryState::NONE;
                    pn->add_arc(PotentialArc(false, false, _mesh.property(hp_trail, incoming_he.opp()), _mesh.property(ep_travel_distance, incoming_he.edge()) - _mesh.property(ep_edge_travel_distance, incoming_he.edge()),boundary, incoming_he.feature(), incoming_he.edge()));
                }

                //For each neighbouring pair of edges, check wether they are the same, or have had a head on collision. If not, then add a corner.
                for (size_t i = 0; i < pn->potential_arcs.size(); i++)
                {
                    auto collisionpair = std::make_pair(pn->potential_arcs[i].motorcycle, pn->potential_arcs[(i + 1) % pn->potential_arcs.size()].motorcycle);
                    bool sameedge = collisionpair.first == collisionpair.second;
                    auto& a0 = pn->potential_arcs[i];
                    auto& a1 = pn->potential_arcs[(i + 1) % pn->potential_arcs.size()];
                    if(sameedge){
                        //Self intersection
                        if(a0.incoming_trace && a1.incoming_trace)
                            sameedge = false;
                        if(a0.incoming_trace && !a1.incoming_trace && a0.travel_distance * .999 > a1.travel_distance)
                            sameedge = false;
                        if(!a0.incoming_trace && a1.incoming_trace && a0.travel_distance < .999 * a1.travel_distance)
                            sameedge = false;
                    }
                    if (!sameedge)
                    {
                        if (std::find(head_on_collisions.begin(), head_on_collisions.end(), collisionpair) != head_on_collisions.end())
                            sameedge = true;

                        std::swap(collisionpair.first, collisionpair.second);
                        if (std::find(head_on_collisions.begin(), head_on_collisions.end(), collisionpair) != head_on_collisions.end())
                            sameedge = true;
                    }
                    pn->potential_arcs[i].is_next_corner = !sameedge;
                }
            }

            auto nh = _tmesh.vertices().add();
            pn->node_handle = nh;

            _tmesh_pos[nh] = tg::pos3(_mesh.point(vertex));
            original_vertex[nh] = vertex;

            potential_nodes.push_back(pn);
        }
    }

    // Start building Graph.
    for (size_t motorcycle_index = 1; motorcycle_index <= number_motorcycles; motorcycle_index++)
    {
        // Get all touched nodes.
        auto contains_relevant_arc = [motorcycle_index](const PotentialArc &a) {
            return a.motorcycle == motorcycle_index;
        };

        std::vector<PotentialNode *> touched_nodes;
        for (auto node : potential_nodes)
        {
            if (node->potential_arcs.end() != std::find_if(node->potential_arcs.begin(), node->potential_arcs.end(), contains_relevant_arc))
                touched_nodes.push_back(node);
        }
        if(touched_nodes.empty())
            continue;

        auto touched_nodes_copy = touched_nodes;
        for(const auto& pn : touched_nodes_copy){
            for(auto pnit = pn->potential_arcs.begin(); pnit != pn->potential_arcs.end();){
                auto cmpit = pnit;
                pnit++;
                auto other_it = std::find_if(pnit, pn->potential_arcs.end(), [cmpit](const PotentialArc& a0){return (a0.incoming_trace && cmpit->incoming_trace && a0.motorcycle == cmpit->motorcycle);});
                if(other_it != pn->potential_arcs.end()){
                    PotentialNode* dummy = new PotentialNode();
                    dummy->node_handle = pn->node_handle;

                    auto move_to_dummy = cmpit->travel_distance < other_it->travel_distance ? other_it : cmpit;
                    for(auto it = pn->potential_arcs.begin(); it != pn->potential_arcs.end(); it++){
                        dummy->add_arc(PotentialArc(it->incoming_trace,it->is_next_corner,move_to_dummy == it ? it->motorcycle : -1,it->travel_distance,it->boundary_state,it->is_feature, it->eh, it->arc_handle));
                    }
                    move_to_dummy->motorcycle = -1;
                    touched_nodes.push_back(dummy);
                    break;
                }
            }
        }

        //Sort by distance traveled of relevant motorcycle.
        std::sort(touched_nodes.begin(), touched_nodes.end(), [&](PotentialNode *n0, PotentialNode *n1) {
            auto ait0 = std::find_if(n0->potential_arcs.begin(), n0->potential_arcs.end(), contains_relevant_arc);
            auto ait1 = std::find_if(n1->potential_arcs.begin(), n1->potential_arcs.end(), contains_relevant_arc);
            return ait0->travel_distance < ait1->travel_distance;
        });

        for (size_t n = 0; n < touched_nodes.size() - 1; n++)
        {
            auto start_node = touched_nodes[n];
            auto end_node = touched_nodes[n + 1];
            auto arc_start = start_node->get_arc(motorcycle_index, false);
            auto arc_end = end_node->get_arc(motorcycle_index, true);

            start_node->init = true;
            end_node->init = true;

            auto vh_from = original_vertex[start_node->node_handle];
            auto vh_to = original_vertex[end_node->node_handle];

            std::vector<OpenMesh::VertexHandle> vh_vec = list_vh_between_nodes(motorcycle_index,vh_from,vh_to,_mesh,hp_trail);

            if (vh_vec.empty())
            {
                _tmesh.edges().add_or_get(start_node->node_handle, end_node->node_handle);
                continue;
            }

            auto vh_tmesh_prev = _tmesh.vertices().add();
            _tmesh_pos[vh_tmesh_prev] = tg::pos3(_mesh.point(vh_vec.at(0)));

            auto ah_new = _tmesh.edges().add_or_get(start_node->node_handle, vh_tmesh_prev);

            for (int i = 1; i < vh_vec.size(); ++i)
            {
                auto vh_tmesh = _tmesh.vertices().add();
                _tmesh_pos[vh_tmesh] = tg::pos3(_mesh.point(vh_vec.at(i)));

                const auto ah_new = _tmesh.edges().add_or_get(vh_tmesh_prev, vh_tmesh);
                vh_tmesh_prev = vh_tmesh;
            }

            ah_new = _tmesh.edges().add_or_get(vh_tmesh_prev, end_node->node_handle);
        }
    }
}
