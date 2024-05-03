#include "FrameFieldTracing/TMesh.hh"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/DefaultTriMesh.hh"
#include "OpenMesh/Core/Utils/Property.hh"
#include "glow-extras/glfw/GlfwContext.hh"
#include <glow-extras/viewer/view.hh>
#include <filesystem>
#include <string>
#include <fstream>

#include <FrameFieldTracing/DiscreteMotorcycleTracer.hh>


void read_vertex_property(OpenMesh::VPropHandleT<int>& _ph, OpenMesh::TriMesh& _mesh, const std::filesystem::path& _path)
{
    std::string line;
    std::ifstream file(_path);

    int curr_idx = 0;
    if (file.is_open())
    {
        while (getline(file, line))
        {
            _mesh.property(_ph, OpenMesh::VertexHandle(curr_idx++)) = std::stoi(line);
        }

        file.close();
    }
    else
    {
        std::cerr << "Could not read file " << _path << std::endl;
    }

    if (curr_idx < _mesh.n_vertices())
        std::cerr << "Problem with input file" << std::endl;
}

void read_edge_property(OpenMesh::EPropHandleT<double>& _ph, OpenMesh::TriMesh& _mesh, const std::filesystem::path& _path)
{
    std::string line;
    std::ifstream file(_path);

    int curr_idx = 0;
    if (file.is_open())
    {
        while (getline(file, line))
        {
            _mesh.property(_ph, OpenMesh::EdgeHandle(curr_idx++)) = std::stod(line);
        }
        file.close();
    }
    else
    {
        std::cerr << "Could not read file " << _path << std::endl;
    }

    if (curr_idx < _mesh.n_edges())
        std::cerr << "Problem with input file" << std::endl;
}

void read_face_property(OpenMesh::FPropHandleT<OpenMesh::Vec3d>& _ph, OpenMesh::TriMesh& _mesh, const std::filesystem::path& _path)
{
    std::string line;
    std::ifstream file(_path);

    int curr_idx = 0;
    if (file.is_open())
    {
        while (getline(file, line))
        {
            OpenMesh::Vec3d vec;
            std::istringstream iss(line);
            iss >> vec[0] >> vec[1] >> vec[2];
            _mesh.property(_ph, OpenMesh::FaceHandle(curr_idx++)) = vec;
        }
        file.close();
    }
    else
    {
        std::cerr << "Could not read file " << _path << std::endl;
    }

    if (curr_idx < _mesh.n_faces())
        std::cerr << "Problem with input file" << std::endl;
}

void read_mesh(OpenMesh::TriMesh& _mesh, const std::string& _path)
{
    _mesh.request_edge_status();
    _mesh.request_vertex_status();
    _mesh.request_halfedge_status();
    _mesh.request_face_status();

    //Read Cross-Field Information
    OpenMesh::FPropHandleT<OpenMesh::Vec3d> u;
    OpenMesh::FPropHandleT<OpenMesh::Vec3d> v;
    OpenMesh::EPropHandleT<double> ep_pjump;
    OpenMesh::VPropHandleT<int> vp_index;

    _mesh.add_property(u, "dir_u");
    _mesh.add_property(v, "dir_v");
    _mesh.add_property(vp_index, "DF_VP_idx");
    _mesh.add_property(ep_pjump, "DF_EP_pjump");
    _mesh.property(u).set_persistent(true);
    _mesh.property(v).set_persistent(true);
    _mesh.property(vp_index).set_persistent(true);
    _mesh.property(ep_pjump).set_persistent(true);


    OpenMesh::IO::Options opts = OpenMesh::IO::Options::Custom;
    OpenMesh::IO::read_mesh(_mesh, _path, opts);
}

pm::vertex_attribute<tg::pos3> to_poly_mesh(const OpenMesh::TriMesh& _mesh_om, pm::Mesh& _mesh_pm)
{
    pm::vertex_attribute<tg::pos3> pos(_mesh_pm);
    for (auto vh : _mesh_om.vertices())
    {
        pos[_mesh_pm.vertices().add()] = tg::pos3(_mesh_om.point(vh));
    }
    for (auto fh : _mesh_om.faces())
    {
        auto vertcies = fh.vertices_ccw();
        std::vector<pm::vertex_handle> vercies_pm = {_mesh_pm.vertices()[vertcies.to_vector()[0].idx()],
                                                     _mesh_pm.vertices()[vertcies.to_vector()[1].idx()],
                                                     _mesh_pm.vertices()[vertcies.to_vector()[2].idx()]};
        _mesh_pm.faces().add(vercies_pm);
    }

    return pos;
}

struct TMesh
{
    pm::Mesh mesh;
    pm::vertex_attribute<tg::pos3> pos;
};

int main()
{
    glow::glfw::GlfwContext ctx;
    auto style = gv::config(gv::dark_ui, gv::no_grid, gv::no_outline, gv::background_color(tg::color3::white), gv::ssao_power(0.5f),
                            gv::shadow_screen_fadeout_distance(50.f), gv::no_shadow);

    std::vector<OpenMesh::TriMesh> meshes(3);
    read_mesh(meshes[0], DATA_PATH "/spot.om");
    read_mesh(meshes[1], DATA_PATH "/fertility.om");
    read_mesh(meshes[2], DATA_PATH "/rocker_arm.om");

    std::vector<TMesh> tmeshes(meshes.size());

    #pragma omp parallel for
    for (int i = 0; i < meshes.size(); ++i)
        compute_t_mesh(meshes[i], tmeshes[i].mesh, tmeshes[i].pos);

    auto viewer = gv::grid();
    for (int i = 0; i < meshes.size(); ++i)
    {
        pm::Mesh m;
        auto pos = to_poly_mesh(meshes[i], m);

        auto c = gv::CameraController::create();
        auto v = gv::view(gv::lines(tmeshes[i].pos).line_width_px(2.), tg::color3::red, c);
        gv::view(pos);
        gv::view(gv::lines(pos).line_width_px(.5));
    }
}
