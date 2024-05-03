#pragma once

#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>

class DiscreteMotorcycle
{
	public:
    DiscreteMotorcycle(int idx, int dir, OpenMesh::SmartHalfedgeHandle& starting_edge, OpenMesh::SmartVertexHandle& start_vertex);
	double GetLength() {return length;};
	void SetLength(double l){length = l;};

	// To identify motorcycle tracks.
	int index;
	// Either 0: u, 1: v, 2: -u, 3: -v;
	int direction;
	double length;
	double drift;
	//Maybe only the last halfedge is needed.
    OpenMesh::SmartVertexHandle current_vertex;
	OpenMesh::SmartVertexHandle start_vertex;
    bool crash;
};

class DiscreteMotorcycleTracer
{
public:
    DiscreteMotorcycleTracer(OpenMesh::TriMesh& mesh,double angle_threshold);
    void generate();
	std::vector<std::pair<int,int>> get_head_on_collisions();
		int get_motorcycle_number();
private:
struct TravelCandidate
{
    enum TravelType{
        EGDE,
        SPLIT
    };
	TravelCandidate(OpenMesh::SmartHalfedgeHandle he, int dir, double cost, int head_on);
	OpenMesh::SmartHalfedgeHandle halfedge;
	int direction;
	double cost;
	int head_on;
    TravelType type;
    double splitparam;

    bool operator==(const TravelCandidate& cmp)
    {
        return this->halfedge == cmp.halfedge && this->direction == cmp.direction;
    }
};

struct MotorcycleStart{
        MotorcycleStart(OpenMesh::SmartHalfedgeHandle oh, int od) : original_halfedge(oh), original_direction(od) {};
        OpenMesh::SmartHalfedgeHandle original_halfedge;
        int original_direction;
        std::vector<TravelCandidate> edge_candidates;
};
enum CollisionAvoidance{
    BOTH,
    LEFT,
    RIGHT
};
	double get_cost(std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> frame, OpenMesh::Vec3d travel_vector, OpenMesh::Vec3d normal_vector, double drift);
	double get_cost_halfedge(OpenMesh::SmartHalfedgeHandle h, int direction, double cost);
	void start_motorcycles();
	void trace_motorcycles();
	double calculate_drift(std::pair<OpenMesh::Vec3d,OpenMesh::Vec3d> frame, OpenMesh::Vec3d travel_vector, OpenMesh::Vec3d normal_vector);
	double calculate_drift_halfedge(OpenMesh::SmartHalfedgeHandle h, int dir);
	double get_travel_distance(std::pair<OpenMesh::Vec3d,OpenMesh::Vec3d> frame, OpenMesh::Vec3d travel_vector, OpenMesh::Vec3d normal_vector);
	double get_travel_distance_halfedge(OpenMesh::SmartHalfedgeHandle h,int dir);
	double get_angle_to_direction(OpenMesh::SmartHalfedgeHandle h, int dir);
	std::vector<int> get_starting_directions(OpenMesh::SmartHalfedgeHandle eh, OpenMesh::SmartVertexHandle vertex);
	void update_cycle(DiscreteMotorcycle& cycle, int direction,OpenMesh::SmartHalfedgeHandle& halfedge);
	void get_startedge_candidates(int starting_direction, std::vector<OpenMesh::SmartHalfedgeHandle> halfedges,std::vector<TravelCandidate>& startedge_candidates,bool clockwise);
    bool evaluate_split(double& cost, OpenMesh::SmartHalfedgeHandle h, int direction, double drift, double& splitpoint);
	void split_edge(OpenMesh::SmartHalfedgeHandle h);
    void split_edge(OpenMesh::SmartHalfedgeHandle h, double splitpoint);
	void refine_around_singularites();
	std::pair<double,double> decompose_travel_direction(std::pair<OpenMesh::Vec3d,OpenMesh::Vec3d> frame, OpenMesh::Vec3d travel_vector, OpenMesh::Vec3d normal_vector);
	std::pair<double, double> decompose_travel_direction_halfedge(OpenMesh::SmartHalfedgeHandle h, int dir);
    CollisionAvoidance check_for_parallel_collision(OpenMesh::SmartHalfedgeHandle h, int direction);
    int check_for_head_on_collision(OpenMesh::SmartHalfedgeHandle h, int direction);
    bool is_imploding_or_exploding(OpenMesh::SmartHalfedgeHandle h, int dir);

    OpenMesh::TriMesh& mesh;

	std::vector<std::pair<int,int>> head_on_crashes;
	double total_angle_threshold;

	std::vector<DiscreteMotorcycle> motorcycles;
    int motorcycle_index;
    int number_of_original_vertices;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> u;
	OpenMesh::FPropHandleT<OpenMesh::Vec3d> v;
	OpenMesh::EPropHandleT<int> matching;
	OpenMesh::VPropHandleT<int> vp_index;
	OpenMesh::VPropHandleT<int> VP_motorcycle_trail_index;
	OpenMesh::VPropHandleT<int> VP_motorcycle_active_index;
	OpenMesh::HPropHandleT<int> HP_motorcycle_trail_index;
    OpenMesh::EPropHandleT<double> ep_pjump;
	OpenMesh::EPropHandleT<double> ep_travel_distance_edge;
	OpenMesh::EPropHandleT<double> ep_total_travel_distance;
    OpenMesh::VPropHandleT<bool> vp_crash_vertex;
    OpenMesh::VPropHandleT<bool> vp_original;
    OpenMesh::EPropHandleT<bool> ep_original;

	OpenMesh::HPropHandleT<int> travel_direction;
	OpenMesh::EPropHandleT<OpenMesh::Vec3d> edgenormal;

    std::vector<std::pair<OpenMesh::SmartHalfedgeHandle, int> > get_outgoing_directions(OpenMesh::SmartHalfedgeHandle h, int dir);
    double resulting_total_angle(double length, double drift, int direction, OpenMesh::SmartHalfedgeHandle h);
    int get_direction_in_neighbour(int direction, OpenMesh::SmartHalfedgeHandle he);
    int get_best_direction(OpenMesh::SmartHalfedgeHandle he);
    std::vector<TravelCandidate> get_tracing_candidates(OpenMesh::SmartVertexHandle vertex, OpenMesh::SmartHalfedgeHandle entry_halfedge, DiscreteMotorcycle &cycle, int direction, double angle_treshold);
};
