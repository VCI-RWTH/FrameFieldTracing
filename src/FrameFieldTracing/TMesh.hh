#pragma once

#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>
#include <polymesh/Mesh.hh>
#include <typed-geometry/types/pos.hh>

void compute_t_mesh(OpenMesh::TriMesh _mesh, polymesh::Mesh &_tmesh, pm::vertex_attribute<tg::pos3> &_tmesh_pos, double _angle_treshold = 0.02);
