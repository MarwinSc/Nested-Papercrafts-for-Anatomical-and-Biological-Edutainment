// boolean.h : Include file for standard system include files,
// or project specific include files.

#ifndef BOOLEAN_H_
#define BOOLEAN_H_

#pragma once

#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/squared_distance_3.h>
#include <fstream>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;

typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::Vertex_handle Vertex_handle;

namespace ae {
	class boolean_interface {
	public:
		boolean_interface();
		~boolean_interface();
		void boolean(std::string first, std::string second);
		void boolUnion(std::string first, std::string second);

		Mesh merge_vertices_by_distance(Mesh first_mesh, double threshold, K::Plane_3 cut_plane);
		void triangulateCut(std::string first, float n_x, float n_y, float n_z, float o_x, float o_y, float o_z);

	private:
		Polyhedron boolean_interface::load(std::string file);
		void boolean_interface::remesh(Polyhedron mesh);

		Mesh connectBoundaries(Mesh mesh, K::Vector_3 compareNormal);
		std::vector<float> getEdgeLengths(halfedge_descriptor hd, Mesh mesh);
		float triangleAspectRatio(face_descriptor fd, Mesh mesh);
		Mesh removeDegenFaces(Mesh mesh, float threshold, K::Plane_3 cut_plane);
	};
}

extern "C" 	__declspec(dllexport) ae::boolean_interface * __stdcall _boolean_interface();

extern "C" 	__declspec(dllexport) void __stdcall _boolean(ae::boolean_interface * g, char* first, char* second);

extern "C" 	__declspec(dllexport) void __stdcall _boolUnion(ae::boolean_interface * g, char* first, char* second);

extern "C" 	__declspec(dllexport) void __stdcall _triangulateCut(ae::boolean_interface * g, char* first, float x, float y, float z, float o_x, float o_y, float o_z);

#endif