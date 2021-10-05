// boolean.cpp : Defines the entry point for the application.
//

#include "boolean.h"
#include <vector>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <CGAL/Vector_3.h>

#include <fstream>

namespace ae{

	namespace PMP = CGAL::Polygon_mesh_processing;
	namespace params = PMP::parameters;
	namespace SMS = CGAL::Surface_mesh_simplification;

	struct Vector_pmap_wrapper {
		std::vector<bool>& vect;
		Vector_pmap_wrapper(std::vector<bool>& v) : vect(v) {}
		friend bool get(const Vector_pmap_wrapper& m, face_descriptor f)
		{
			return m.vect[f];
		}
		friend void put(const Vector_pmap_wrapper& m, face_descriptor f, bool b)
		{
			m.vect[f] = b;
		}
	};

	boolean_interface::boolean_interface(){}

	boolean_interface::~boolean_interface()
	{
	}
	

	void boolean_interface::merge(std::string first, float threshold, float n_x, float n_y, float n_z, float o_x, float o_y, float o_z) {
	
		std::cout << first << std::endl;
		const char* filename1 = first.c_str();
		std::ifstream input(filename1);
		Mesh mesh;
		if (!input || !(input >> mesh))
		{
			std::cerr << "First mesh is not a valid off file." << std::endl;
		}
		input.close();

		K::Plane_3 cut_plane = K::Plane_3(K::Point_3(o_x, o_y, o_z), K::Vector_3(n_x, n_y, n_z));
		mesh = boolean_interface::merge_vertices_by_distance(mesh, threshold, cut_plane);
		mesh = boolean_interface::merge_vertices_by_distance(mesh, threshold, cut_plane);

		mesh.collect_garbage();
		std::ofstream output("../out/3D/merged.off");
		output.precision(17);
		output << mesh;
		output.close();
		output.clear();
	}

	void boolean_interface::triangulateCut(std::string first, float threshold, float n_x, float n_y, float n_z, float o_x, float o_y, float o_z) {
		std::cout << first << std::endl;
		const char* filename1 = first.c_str();
		std::ifstream input(filename1);
		Mesh mesh;
		if (!input || !(input >> mesh))
		{
			std::cerr << "First mesh is not a valid off file." << std::endl;
		}
		input.close();

		K::Plane_3 cut_plane = K::Plane_3(K::Point_3(o_x, o_y, o_z), K::Vector_3(n_x, n_y, n_z));
		mesh = boolean_interface::merge_vertices_by_distance(mesh, threshold, cut_plane);
		mesh = boolean_interface::merge_vertices_by_distance(mesh, threshold, cut_plane);

		mesh = boolean_interface::removeDegenFaces(mesh,0.2f,cut_plane);
		mesh = boolean_interface::connectBoundaries(mesh, K::Vector_3(n_x, n_y, n_z));

	}

	Mesh boolean_interface::connectBoundaries(Mesh mesh, K::Vector_3 compareNormal) {
		std::vector<halfedge_descriptor> boundary_halfedges;
		std::vector<vertex_descriptor> boundary_vertices;

		//collect all boundary halfedges and vertices
		for (halfedge_descriptor hd : mesh.halfedges()) {
			if (mesh.face(hd) == Mesh::null_face()) {
				boundary_halfedges.push_back(hd);
				boundary_vertices.push_back(mesh.target(hd));
			}
		}

		//split them into the two boundary loops
		std::vector<halfedge_descriptor> boundary_halfedges_second;
		halfedge_descriptor startPoint = boundary_halfedges.at(0);
		halfedge_descriptor current = boundary_halfedges.at(0);
		boundary_halfedges_second.push_back(startPoint);
		do {
			//iterate over one boundary loop 
			for (halfedge_descriptor hd : halfedges_around_source(mesh.target(current), mesh)) {
				if (mesh.face(hd) == Mesh::null_face() && hd != current) {
					current = hd;
					boundary_halfedges_second.push_back(hd);
				}
			}
		} while (startPoint != current);
		//remove from other vector
		for (int i = 0; i < boundary_halfedges.size(); i++) {
			if (std::count(boundary_halfedges_second.begin(), boundary_halfedges_second.end(), boundary_halfedges.at(i))) {
				boundary_halfedges.erase(boundary_halfedges.begin() + i);
				i--;
			}
		}

		//pick one point at random and find the closest in the other loop
		bool found = false;
		int i = 0;
		halfedge_descriptor secondPoint;
		while (!found) {
			startPoint = boundary_halfedges_second.at(i++);
			float smallestDis = 1000000.0f;
			for (halfedge_descriptor hd : boundary_halfedges) {
				float distance = CGAL::squared_distance(mesh.point(mesh.target(startPoint)), mesh.point(mesh.target(hd)));

				//check that there isn't a closer vertex in loop "boundary_halfedges_second" so the first added triangle won't self intersect 
				float sameloop_smallestDis = 1000000.0f;
				for (halfedge_descriptor sameloop_hd : boundary_halfedges_second) {
					float sameloop_distance = CGAL::squared_distance(mesh.point(mesh.target(sameloop_hd)), mesh.point(mesh.target(hd)));
					if (sameloop_distance < sameloop_smallestDis) {
						sameloop_smallestDis = sameloop_distance;
					}
				}

				if (distance < smallestDis && distance <= sameloop_smallestDis) {
					smallestDis = distance;
					secondPoint = hd;
					found = true;
				}
			}
		}

		//order other vector
		std::vector<halfedge_descriptor> boundary_halfedges_first;
		startPoint = secondPoint;
		current = secondPoint;
		boundary_halfedges_first.push_back(startPoint);
		do {
			//iterate over one boundary loop 
			for (halfedge_descriptor hd : halfedges_around_target(mesh.source(current), mesh)) {
				if (mesh.face(hd) == Mesh::null_face() && hd != current) {
					current = hd;
					boundary_halfedges_first.push_back(hd);
				}
			}
		} while (startPoint != current);

		//iterate over length and add triangles

		Mesh m;

		int size_first = boundary_halfedges_first.size();
		int size_second = boundary_halfedges_second.size();

		int indexFirstLoop = 0;
		int indexSecondLoop = 0;

		while (indexFirstLoop < size_first - 1 || indexSecondLoop < size_second - 1) {

			Mesh tempMesh = Mesh(mesh);

			bool firstTriValid = false;
			float aspectRatio1 = 0.0;
			bool secondTriValid = false;
			float aspectRatio2 = 0.0;

			if (indexFirstLoop < size_first - 1) {

				vertex_descriptor u = tempMesh.target(boundary_halfedges_second.at(indexSecondLoop));
				vertex_descriptor v = tempMesh.target(boundary_halfedges_first.at(indexFirstLoop));
				vertex_descriptor w = tempMesh.target(boundary_halfedges_first.at(indexFirstLoop + 1));

				face_descriptor newFace = tempMesh.add_face(w, v, u);
				K::Vector_3 normal = PMP::compute_face_normal(newFace, tempMesh);
				float s = CGAL::scalar_product(normal, compareNormal);
				aspectRatio1 = std::abs(1.33333 - boolean_interface::triangleAspectRatio(newFace, tempMesh));
				
				//offset edgeloops
				
				if (indexSecondLoop < size_second - 1){
					vertex_descriptor vd = tempMesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));
					K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) / 100.0;
					K::Point_3 pt = tempMesh.point(vd);
					K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
					tempMesh.point(vd) = new_pt;
				}
				/*
				for (vertex_descriptor vd : boundary_vertices) {
					if (vd != u && vd != v && vd != w) {
						K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) / 3.0;
						K::Point_3 pt = tempMesh.point(vd);
						K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
						tempMesh.point(vd) = new_pt;
					}
				}
				*/
				tempMesh.collect_garbage();
				std::ofstream temp("temp.off");
				temp.precision(17);
				temp << tempMesh;
				temp.close();
				temp.clear();

				bool intersection = PMP::does_self_intersect(tempMesh);
				std::cout << "First Loop Scalarproduct: " << s << std::endl;
				std::cout << "Second Loop Intersection: " << intersection << std::endl;

				if (s > 0.0 && !intersection) {
					firstTriValid = true;
				}
			}

			Mesh tempMesh2 = Mesh(mesh);

			if (indexSecondLoop < size_second - 1) {

				vertex_descriptor u = tempMesh.target(boundary_halfedges_second.at(indexSecondLoop));
				vertex_descriptor v = tempMesh.target(boundary_halfedges_first.at(indexFirstLoop));
				vertex_descriptor w = tempMesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));

				face_descriptor newFace = tempMesh2.add_face(w, v, u);
				K::Vector_3 normal = PMP::compute_face_normal(newFace, tempMesh2);
				float s = CGAL::scalar_product(normal, compareNormal);
				aspectRatio2 = std::abs(1.33333 - boolean_interface::triangleAspectRatio(newFace, tempMesh2));
				
				//offset edgeloops
				
				if (indexFirstLoop < size_first - 1) {
					vertex_descriptor vd = tempMesh2.target(boundary_halfedges_first.at(indexFirstLoop + 1));
					K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) / 100.0;
					K::Point_3 pt = tempMesh2.point(vd);
					K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
					tempMesh2.point(vd) = new_pt;
				}
				/*
				for (vertex_descriptor vd : boundary_vertices) {
					if (vd != u && vd != v && vd != w) {
						K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) / 3.0;
						K::Point_3 pt = tempMesh2.point(vd);
						K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
						tempMesh2.point(vd) = new_pt;
					}
				}
				*/
				tempMesh2.collect_garbage();
				std::ofstream temp2("temp2.off");
				temp2.precision(17);
				temp2 << tempMesh2;
				temp2.close();
				temp2.clear();

				bool intersection = PMP::does_self_intersect(tempMesh2);
				std::cout << "Second Loop Scalarproduct: " << s << std::endl;
				std::cout << "Second Loop Intersection: " << intersection << std::endl;
				
				if (s > 0.0 && !intersection) {
					secondTriValid = true;
				}
			}

			if (firstTriValid && secondTriValid) {
				if (aspectRatio1 < aspectRatio2) {
					vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
					vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
					vertex_descriptor w = mesh.target(boundary_halfedges_first.at(indexFirstLoop + 1));
					face_descriptor newFace = mesh.add_face(w, v, u);
					indexFirstLoop++;
				}
				else {
					vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
					vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
					vertex_descriptor w = mesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));
					face_descriptor newFace = mesh.add_face(w, v, u);
					indexSecondLoop++;
				}
			}
			else if (firstTriValid) {
				vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
				vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
				vertex_descriptor w = mesh.target(boundary_halfedges_first.at(indexFirstLoop + 1));
				face_descriptor newFace = mesh.add_face(w, v, u);
				indexFirstLoop++;
			}
			else if (secondTriValid) {
				vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
				vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
				vertex_descriptor w = mesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));
				face_descriptor newFace = mesh.add_face(w, v, u);
				indexSecondLoop++;
			}

			mesh.collect_garbage();
			std::cout << "connected boundaries\n";
			std::ofstream output("../out/3D/connected.off");
			output.precision(17);
			output << mesh;
			output.close();
			output.clear();
		}
		return mesh;
	}

	Mesh boolean_interface::merge_vertices_by_distance(Mesh first_mesh, double threshold, K::Plane_3 cut_plane) {

		int count = 0;

		Mesh::Property_map<vertex_descriptor, K::Point_3> location = first_mesh.points();
		Mesh::Vertex_range range = first_mesh.vertices();
		Mesh::Vertex_range::iterator  vb, ve;
		vb = boost::begin(range);
		ve = boost::end(range);
		for (boost::tie(vb, ve) = first_mesh.vertices(); vb != ve; ++vb) {
			vertex_descriptor vd = *vb;

			CGAL::Vertex_around_target_circulator<Mesh> vbegin(first_mesh.halfedge(vd), first_mesh), done(vbegin);

			do {
				vertex_descriptor vd2 = *vbegin++;
				//std::cout << vd << " to " << vd2 << std::endl;
				float distance = CGAL::squared_distance(location[vd], location[vd2]);
				//check distance, all other conditions most likely got obsolete
				if (distance < threshold && (vd != vd2) && !(first_mesh.is_removed(vd)) && !(first_mesh.is_removed(vd2))) {

					K::Point_3 loc = location[vd];
					K::Point_3 loc2 = location[vd2];
					//std::cout << vd << " at " << loc[2] << std::endl;
					//std::cout << " to " << std::endl;
					//std::cout << vd2 << " at " << loc2[2] << std::endl;
					//std::cout << "has length " << distance << std::endl;

					//vertex merge
					vertex_descriptor new_vd = CGAL::Euler::collapse_edge(first_mesh.edge(first_mesh.halfedge(vd, vd2)), first_mesh);

					count++;
					//CGAL::Euler::join_vertex(first_mesh.halfedge(vd2, vd), first_mesh);


					K::Point_3 point = first_mesh.point(new_vd);
					if (!cut_plane.has_on(point)) {
						first_mesh.point(new_vd) = cut_plane.projection(point);
					}
					//if the source is lower than the target assign its z coordinate to target.
					//if (loc[2] > loc2[2]) {
					//	first_mesh.point(new_vd) = K::Point_3(loc[0], loc[1], loc2[2]);
					//}
					//repeat once more for the current vertex because of merge
					vb--;
					break;
				}
			} while (vbegin != done);
		}

		first_mesh.collect_garbage();
		std::cout << "Merge was successfully computed\n";
		std::cout << "Removed Vertices: " << count << std::endl;
		std::ofstream output("../out/3D/merged_vertices.off");
		output.precision(17);
		output << first_mesh;
		output.close();
		output.clear();

		return first_mesh;
	}

	//removes needle triangles based on aspect ratio and area threshold
	Mesh boolean_interface::removeDegenFaces(Mesh mesh, float threshold, K::Plane_3 cut_plane) {
		for (face_descriptor fd : mesh.faces()) {
			std::vector<float> lengths = boolean_interface::getEdgeLengths(mesh.halfedge(fd), mesh);
			if (lengths.size() == 3)
			{
				float aspectRatio = boolean_interface::triangleAspectRatio(fd, mesh);

				float area = CGAL::Polygon_mesh_processing::face_area(fd, mesh);

				halfedge_descriptor current = mesh.halfedge(fd);
				K::Point_3 p1 = mesh.point(mesh.target(current));
				current = mesh.next(current);
				K::Point_3 p2 = mesh.point(mesh.target(current));
				current = mesh.next(current);
				K::Point_3 p3 = mesh.point(mesh.target(current));

				bool sameZ = (cut_plane.has_on(p1) && cut_plane.has_on(p2) && cut_plane.has_on(p3));

				std::cout << "aspectRatio: " << aspectRatio << "area: " << area << "\n";
				if (std::abs(aspectRatio) > 50000 && area < 80.0f && sameZ) {

					CGAL::Euler::remove_face(mesh.halfedge(fd), mesh);
					std::cout << "remove face\n";

				}
			}
			else {
				std::cerr << "removeDegenFaces(): Not a valid triangle" << std::endl;
			}
		}
		mesh.collect_garbage();
		std::cout << "Degenerate triangles removed\n";
		std::ofstream output("../out/3D/removed_degeneratedTris.off");
		output.precision(17);
		output << mesh;
		output.close();
		output.clear();

		return mesh;
	}

	float boolean_interface::triangleAspectRatio(face_descriptor fd, Mesh mesh) {
		float aspectRatio = 0.0;
		std::vector<float> lengths = boolean_interface::getEdgeLengths(mesh.halfedge(fd), mesh);
		if (lengths.size() == 3)
		{
			float a = lengths.at(0);
			float b = lengths.at(1);
			float c = lengths.at(2);

			float longestEdge = 0.0;

			if (a >= b) {
				longestEdge = a;
			}
			else {
				longestEdge = b;
			}
			if (c > longestEdge) {
				longestEdge = c;
			}
			//float s = (a + b + c) / 2.0f;
			//aspectRatio = (a * b * c) / (8 * (s - a) * (s - b) * (s - c));

			float area = CGAL::Polygon_mesh_processing::face_area(fd, mesh);

			aspectRatio = (std::pow(longestEdge, 2.0) / area);
		}
		else {
			std::cerr << "triangleAspectRatio(): Not a valid triangle" << std::endl;
		}
		return aspectRatio;
	}

	std::vector<float> boolean_interface::getEdgeLengths(halfedge_descriptor hd, Mesh mesh) {
		halfedge_descriptor current = hd;
		std::vector<float> lengths;
		do {
			float temp = (CGAL::squared_distance(mesh.point(mesh.target(current)), mesh.point(mesh.target(mesh.next(current)))));
			lengths.push_back(temp);
			current = mesh.next(current);
			std::cout << "edgelength: " << temp << "\n";
		} while (current != hd);
		return lengths;
	}
}

ae::boolean_interface* __stdcall _boolean_interface()
{
	return new ae::boolean_interface();
}

void __stdcall _triangulateCut(ae::boolean_interface* g, char* first, float t, float x, float y, float z, float o_x, float o_y, float o_z)
{
	g->triangulateCut(std::string(first), t, x, y, z, o_x, o_y, o_z);
}

void __stdcall _merge(ae::boolean_interface * g, char* first, float t, float x, float y, float z, float o_x, float o_y, float o_z)
{
	g->merge(std::string(first), t, x, y, z, o_x, o_y, o_z);
}