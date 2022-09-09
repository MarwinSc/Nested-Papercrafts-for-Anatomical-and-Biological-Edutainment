// boolean.cpp : Defines the entry point for the application.
//

#include "boolean.h"
#include <vector>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

#include <CGAL/Vector_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <fstream>

namespace ae{

	namespace PMP = CGAL::Polygon_mesh_processing;
	namespace params = PMP::parameters;
	namespace SMS = CGAL::Surface_mesh_simplification;
	namespace VSA = CGAL::Surface_mesh_approximation;

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

	void boolean_interface::simplify(std::string first, double stop_ratio) {
		std::cout << first << std::endl;
		const char* filename1 = first.c_str();
		std::ifstream input(filename1);
		Mesh mesh;
		if (!input || !(input >> mesh))
		{
			std::cerr << "First mesh is not a valid off file." << std::endl;
		}
		input.close();

		
		if (stop_ratio <= 1.0) {
			SMS::Count_ratio_stop_predicate<Mesh> stop(stop_ratio);
			int r = SMS::edge_collapse(mesh, stop);
		}
		else {
			//stop_ratio equals the number of edges in the resulting mesh
			SMS::Count_stop_predicate<Mesh> stop(stop_ratio);
			int r = SMS::edge_collapse(mesh, stop);
		}

		mesh.collect_garbage();
		std::ofstream output("../out/3D/simplified.off");
		output.precision(17);
		output << mesh;
		output.close();
		output.clear();
	}


	void boolean_interface::approximate(std::string first) {
		std::cout << first << std::endl;
		const char* filename1 = first.c_str();
		std::ifstream input(filename1);
		Mesh mesh;
		if (!input || !(input >> mesh))
		{
			std::cerr << "First mesh is not a valid off file." << std::endl;
		}
		input.close();

		boolean_interface::approximate_mesh(mesh);

	}

	typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;
	typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> Mesh_approximation;
	typedef Mesh_approximation::Error_metric L21_metric;

	typedef std::vector<std::size_t> Polygon;

	void boolean_interface::approximate_mesh(Mesh mesh)
	{
		Vertex_point_map vpmap = get(boost::vertex_point, const_cast<Mesh&>(mesh));
		// error metric and fitting function
		L21_metric error_metric(mesh, vpmap);
		// creates VSA algorithm instance
		Mesh_approximation approx(mesh, vpmap, error_metric);
		// seeds 100 random proxies
		approx.initialize_seeds(CGAL::parameters::seeding_method(VSA::HIERARCHICAL)
			//.max_number_of_proxies(15)
			.min_error_drop(0.05));

		// runs 30 iterations
		approx.run(5);
		// adds 3 proxies to the one with the maximum fitting error,
		// running 5 iterations between each addition
		//approx.add_to_furthest_proxies(3, 5);
		// runs 10 iterations
		//approx.run(5);
		// teleports 2 proxies to tunnel out of local minima,
		// running 5 iterations between each teleport
		//approx.teleport_proxies(2, 5);
		// runs 10 iterations
		//approx.run(10);
		// extract approximated mesh with default parameters
		approx.extract_mesh(CGAL::parameters::all_default());
		// get approximated triangle soup
		std::vector<K::Point_3> anchors;
		std::vector<std::array<std::size_t, 3> > triangles;

		approx.output(CGAL::parameters::anchors(std::back_inserter(anchors)).
			triangles(std::back_inserter(triangles)));

		std::vector<Polygon> polygons;

		Polygon p;
		for (auto &tri : triangles) {
			p.push_back(tri.at(0));
			p.push_back(tri.at(1));
			p.push_back(tri.at(2));
			polygons.push_back(p);
			p.clear();
		}

		PMP::repair_polygon_soup(anchors, polygons);

		// convert from soup to surface mesh
		PMP::orient_polygon_soup(anchors, polygons);

		if (!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
			std::cerr << "Not a valid polygon mesh!" << std::endl;
		}

		Mesh output;
		PMP::polygon_soup_to_polygon_mesh(anchors, polygons, output);
		if (CGAL::is_closed(output) && (!PMP::is_outward_oriented(output))) {
			PMP::reverse_face_orientations(output);
		}

		PMP::orient_to_bound_a_volume(output);

		// Fix manifoldness by splitting non-manifold vertices
		PMP::duplicate_non_manifold_vertices(output);

		for (vertex_descriptor v : vertices(output))
		{
			if (PMP::is_non_manifold_vertex(v, output))
			{
				std::cout << "vertex " << v << " is non-manifold" << std::endl;
			}
		}

		PMP::stitch_borders(output);

		std::ofstream out("../out/3D/dump.off");
		out << output;
		out.close();
	}

	void boolean_interface::triangulateCut(std::string first, float threshold, float n_x, float n_y, float n_z, float o_x, float o_y, float o_z) 
	{
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
		//mesh = boolean_interface::merge_vertices_by_distance(mesh, threshold, cut_plane);
		mesh = boolean_interface::merge_vertices_by_distance(mesh, 2.0, cut_plane);
		
		mesh = boolean_interface::connectBoundaries(mesh, K::Vector_3(n_x, n_y, n_z), threshold);

		//mesh = boolean_interface::removeNeedleTriangles(mesh, 100000.0f, cut_plane);

	}
	/*
	Mesh boolean_interface::connectBoundaries(Mesh mesh, K::Vector_3 compareNormal) 
	{
		std::vector<halfedge_descriptor> border_cycles;
		PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

		for (halfedge_descriptor hd : border_cycles)
		{
			std::vector<face_descriptor>  patch_facets;
			std::vector<vertex_descriptor> patch_vertices;
			PMP::triangulate_and_refine_hole(mesh, hd, std::back_inserter(patch_facets), std::back_inserter(patch_vertices));
		}

		mesh.collect_garbage();
		std::cout << "connected boundaries\n";
		std::ofstream output("../out/3D/connected.off");
		output.precision(17);
		output << mesh;
		output.close();
		output.clear();

		return mesh;
	}
	
	*/
	Mesh boolean_interface::connectBoundaries(Mesh mesh, K::Vector_3 compareNormal, float aspectRatioThreshold) {
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

		std::cout << "boundary1 size: " << boundary_halfedges.size() << std::endl;
		std::cout << "boundary2 size: " << boundary_halfedges_second.size() << std::endl;


		if (boundary_halfedges.size() > 0) {


			//pick the vertices with the highest coordinates in both loops 
			halfedge_descriptor highest_firstLoop;
			K::Point_3 highest_pt = K::Point_3(-1000.0, -1000.0, -1000.0);
			double sumHighest_first = highest_pt.x() + highest_pt.y() + highest_pt.z();
			for (halfedge_descriptor hd : boundary_halfedges) {
				K::Point_3 pt = mesh.point(mesh.target(hd));

				sumHighest_first = highest_pt.x() + highest_pt.y() + highest_pt.z();
				double sumCurrent = pt.x() + pt.y() + pt.z();

				std::cout << highest_firstLoop << " " << sumHighest_first << std::endl;
				std::cout << hd << " " << sumCurrent << std::endl;

				if (sumCurrent > sumHighest_first)
				{
					highest_firstLoop = hd;
					highest_pt = pt;
				}
			}

			halfedge_descriptor highest_secondLoop;
			highest_pt = K::Point_3(-1000.0, -1000.0, -1000.0);
			double sumHighest_second = highest_pt.x() + highest_pt.y() + highest_pt.z();
			for (halfedge_descriptor hd : boundary_halfedges_second) {
				K::Point_3 pt = mesh.point(mesh.target(hd));

				sumHighest_second = highest_pt.x() + highest_pt.y() + highest_pt.z();
				double sumCurrent = pt.x() + pt.y() + pt.z();
				if (sumCurrent > sumHighest_second)
				{
					highest_secondLoop = hd;
					highest_pt = pt;
				}
			}

			halfedge_descriptor startPoint_first;
			halfedge_descriptor startPoint_second;
			//to ensures firstLoop is the outer loop
			if (sumHighest_first > sumHighest_second) {
				startPoint_first = highest_firstLoop;
				startPoint_second = highest_secondLoop;
			}
			else {
				startPoint_first = highest_secondLoop;
				startPoint_second = highest_firstLoop;
			}

			std::cout << "start1: " << startPoint_first << std::endl;
			std::cout << "start2: " << startPoint_second << std::endl;

			//order both vectors
			std::vector<halfedge_descriptor> boundary_halfedges_first;
			current = startPoint_first;
			boundary_halfedges_first.push_back(startPoint_first);
			do {
				//iterate over one boundary loop 
				for (halfedge_descriptor hd : halfedges_around_source(mesh.target(current), mesh)) {
					if (mesh.face(hd) == Mesh::null_face() && hd != current) {
						current = hd;
						boundary_halfedges_first.push_back(hd);
					}
				}
			} while (startPoint_first != current);

			boundary_halfedges_second.clear();
			current = startPoint_second;
			boundary_halfedges_second.push_back(startPoint_second);
			do {
				//iterate over one boundary loop 
				for (halfedge_descriptor hd : halfedges_around_target(mesh.source(current), mesh)) {
					if (mesh.face(hd) == Mesh::null_face() && hd != current) {
						current = hd;
						boundary_halfedges_second.push_back(hd);
					}
				}
			} while (startPoint_second != current);

			std::cout << "triangulate cut" << std::endl;

			//iterate over length and add triangles

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

					face_descriptor newFace = tempMesh.add_face(u, v, w);
					if (!(newFace == tempMesh.null_face())) {
						K::Vector_3 normal = PMP::compute_face_normal(newFace, tempMesh);
						float s = CGAL::scalar_product(normal, compareNormal);
						//aspectRatio1 = std::abs(1.33333 - boolean_interface::triangleAspectRatio(newFace, tempMesh));

						//offset edgeloops to provoke intersection
						
						for (halfedge_descriptor hd : boundary_halfedges_second) {
							vertex_descriptor vd = tempMesh.target(hd);
							if (vd != u) {
								K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) ;
								K::Point_3 pt = tempMesh.point(vd);
								K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
								tempMesh.point(vd) = new_pt;
							}
						}
						for (halfedge_descriptor hd : boundary_halfedges_first) {
							vertex_descriptor vd = tempMesh.target(hd);
							if (vd != v && vd != w) {
								K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length()));
								K::Point_3 pt = tempMesh.point(vd);
								K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
								tempMesh.point(vd) = new_pt;
							}
						}
							/*
							* if (indexSecondLoop < size_second - 1) {
							vertex_descriptor vd = tempMesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));
							K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) / 100.0;
							K::Point_3 pt = tempMesh.point(vd);
							K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
							tempMesh.point(vd) = new_pt;
							}
							*/
						

						tempMesh.collect_garbage();
						std::ofstream temp("../out/3D/temp.off");
						temp.precision(17);
						temp << tempMesh;
						temp.close();
						temp.clear();

						bool intersection = PMP::does_self_intersect(tempMesh);
						std::cout << "First Loop Scalarproduct: " << s << std::endl;
						std::cout << "First Loop Intersection: " << intersection << std::endl;

						float aspectRatio1 = PMP::face_aspect_ratio(newFace, tempMesh);
						std::cout << "First Loop aspectratio: " << aspectRatio1 << std::endl;

						if (s > 0.0 && !intersection && aspectRatio1 < aspectRatioThreshold) {
							firstTriValid = true;
						}
					}
				}

				Mesh tempMesh2 = Mesh(mesh);

				if (indexSecondLoop < size_second - 1) {

					vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
					vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
					vertex_descriptor w = mesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));

					face_descriptor newFace = tempMesh2.add_face(u, v, w);
					if (!(newFace == tempMesh2.null_face())) {
						K::Vector_3 normal = PMP::compute_face_normal(newFace, tempMesh2);
						float s = CGAL::scalar_product(normal, compareNormal);
						//aspectRatio2 = std::abs(1.33333 - boolean_interface::triangleAspectRatio(newFace, tempMesh2));

						//offset edgeloops to provoke intersection
						for (halfedge_descriptor hd : boundary_halfedges_second) {
							vertex_descriptor vd = tempMesh2.target(hd);
							if (vd != u && vd != w) {
								K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) ;
								K::Point_3 pt = tempMesh2.point(vd);
								K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
								tempMesh2.point(vd) = new_pt;
							}
						}
						for (halfedge_descriptor hd : boundary_halfedges_first) {
							vertex_descriptor vd = tempMesh2.target(hd);
							if (vd != v) {
								K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length()));
								K::Point_3 pt = tempMesh2.point(vd);
								K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
								tempMesh2.point(vd) = new_pt;
							}
						}
						//if (indexFirstLoop < size_first - 1) {

							//vertex_descriptor vd = tempMesh2.target(boundary_halfedges_first.at(indexFirstLoop + 1));
							//K::Vector_3 tempNormal = (compareNormal / CGAL::sqrt(compareNormal.squared_length())) / 100.0;
							//K::Point_3 pt = tempMesh2.point(vd);
							//K::Point_3 new_pt = K::Point_3(pt.x() + tempNormal.x(), pt.y() + tempNormal.y(), pt.z() + tempNormal.z());
							//tempMesh2.point(vd) = new_pt;
							//}


						tempMesh2.collect_garbage();
						std::ofstream temp2("../out/3D/temp2.off");
						temp2.precision(17);
						temp2 << tempMesh2;
						temp2.close();
						temp2.clear();

						bool intersection = PMP::does_self_intersect(tempMesh2);
						std::cout << "Second Loop Scalarproduct: " << s << std::endl;
						std::cout << "Second Loop Intersection: " << intersection << std::endl;

						float aspectRatio2 = PMP::face_aspect_ratio(newFace, tempMesh2);
						std::cout << "Second Loop aspectratio: " << aspectRatio2 << std::endl;

						if (s > 0.0 && !intersection && aspectRatio2 < aspectRatioThreshold) {
							secondTriValid = true;
						}
					}

				}

				if (firstTriValid && secondTriValid) {
					if (aspectRatio1 < aspectRatio2) {
						secondTriValid = false;
					}
					else {
						firstTriValid = false;
					}
				}

				if (firstTriValid) {
					vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
					vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
					vertex_descriptor w = mesh.target(boundary_halfedges_first.at(indexFirstLoop + 1));
					face_descriptor newFace = mesh.add_face(u, v, w);
					indexFirstLoop++;
				}
				else if (secondTriValid) {
					vertex_descriptor u = mesh.target(boundary_halfedges_second.at(indexSecondLoop));
					vertex_descriptor v = mesh.target(boundary_halfedges_first.at(indexFirstLoop));
					vertex_descriptor w = mesh.target(boundary_halfedges_second.at(indexSecondLoop + 1));
					face_descriptor newFace = mesh.add_face(u, v, w);
					indexSecondLoop++;
				}
				else {

					float ratioOuter = (indexFirstLoop + 1.0f) / static_cast<float>(boundary_halfedges_first.size());
					float ratioInner = (indexSecondLoop + 1.0f) / static_cast<float>(boundary_halfedges_second.size());
					std::cout << "Ratio Inner: " << ratioInner << std::endl;
					std::cout << "Ratio Outer: " << ratioOuter << std::endl;
					if (ratioInner >= ratioOuter) {
						indexFirstLoop++;
					}
					else {
						indexSecondLoop++;
					}
				}

				mesh.collect_garbage();
				std::ofstream output("../out/3D/connected.off");
				output.precision(17);
				output << mesh;
				output.close();
				output.clear();
			}

		}
		std::cout << "fill boundaries\n";


		//fill all remaining holes
		std::vector<halfedge_descriptor> border_cycles;
		PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

		std::vector<face_descriptor> allFacets;

		for (halfedge_descriptor hd : border_cycles)
		{
			std::vector<face_descriptor>  patch_facets;
			std::vector<vertex_descriptor> patch_vertices;
			PMP::triangulate_and_refine_hole(mesh, hd, std::back_inserter(patch_facets), std::back_inserter(patch_vertices));
			allFacets.insert(allFacets.end(), patch_facets.begin(), patch_facets.end());
		}

		if (true){

			std::cout << "1st fill\n";

			int j = 0;

			for (int i = 0; i < 3; i++) {

				for (face_descriptor fd : mesh.faces()) {


					if (std::count(mesh.faces().begin(), mesh.faces().end(), fd)) {

						float aspectRatio = PMP::face_aspect_ratio(fd, mesh);
						if (aspectRatio > 10) {

							for (halfedge_descriptor hd : CGAL::halfedges_around_face(mesh.halfedge(fd), mesh)) {

								std::cout << "checking halfedge\n";

								K::Vector_3 normal = PMP::compute_face_normal(mesh.face(mesh.opposite(hd)), mesh);
							
								//float s = CGAL::scalar_product(normal, compareNormal);
								float s = normal.x() * compareNormal.x() + normal.y() * compareNormal.y() + normal.z() * compareNormal.z();
								std::cout << "scalar Product: " << s << "\n";
								float len_normal = std::sqrt(normal.x() * normal.x() + normal.y() * normal.y() + normal.z() * normal.z());
								float len_compare_normal = std::sqrt(compareNormal.x() * compareNormal.x() + compareNormal.y() * compareNormal.y() + compareNormal.z() * compareNormal.z());
								aspectRatio = PMP::face_aspect_ratio(mesh.face(mesh.opposite(hd)), mesh);

								std::cout << "difference: " << s - (len_normal * len_compare_normal) << "\n";
								if (std::abs(s - (len_normal*len_compare_normal)) < 0.001 && mesh.face(mesh.opposite(hd)) != mesh.null_face() && aspectRatio > 4) {
									std::cout << "removed an additional face\n";

									//remove_faces.insert(mesh.face(mesh.opposite(hd)));
									CGAL::Euler::remove_face(mesh.opposite(hd), mesh);
									j++;
								}
							}
							if (fd != mesh.null_face()) {
								std::cout << "removing original bad \n";
								//remove_faces.insert(fd);
								CGAL::Euler::remove_face(mesh.halfedge(fd), mesh);
								j++;
							}

							mesh.collect_garbage();
							std::ofstream output("../out/3D/temp3.off");
							output.precision(17);
							output << mesh;
							output.close();
							output.clear();

							std::cout << "triangulate and refine \n";
							border_cycles.clear();
							PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

							for (halfedge_descriptor hd : border_cycles)
							{
								std::vector<face_descriptor>  patch_facets;
								std::vector<vertex_descriptor> patch_vertices;
								PMP::triangulate_and_refine_hole(mesh, hd, std::back_inserter(patch_facets), std::back_inserter(patch_vertices));
							}
						}
					}
				}
			}
			std::cout << "removed faces: " << std::to_string(j) << "\n";
		}

		mesh.collect_garbage();
		std::cout << "connected boundaries\n";
		std::ofstream output("../out/3D/connected.off");
		output.precision(17);
		output << mesh;
		output.close();
		output.clear();

		return mesh;
	}
	

	Mesh boolean_interface::merge_vertices_by_distance(Mesh first_mesh, double threshold, K::Plane_3 cut_plane) {

		int count = 0;

		Mesh::Property_map<vertex_descriptor, K::Point_3> location = first_mesh.points();
		bool breakOut = false;
		bool repeat = true;
		while (repeat)
		{
			for (vertex_descriptor vd : first_mesh.vertices()) {

				CGAL::Vertex_around_target_circulator<Mesh> vbegin(first_mesh.halfedge(vd), first_mesh), done(vbegin);

				do {
					vertex_descriptor vd2 = *vbegin++;
					//std::cout << vd << " to " << vd2 << std::endl;
					float distance = CGAL::squared_distance(location[vd], location[vd2]);
		
					if (distance < threshold) {

						K::Point_3 loc = location[vd];
						K::Point_3 loc2 = location[vd2];

						bool was_on_plane = cut_plane.has_on(loc2) || cut_plane.has_on(loc);

						//vertex merge
						vertex_descriptor new_vd = CGAL::Euler::collapse_edge(first_mesh.edge(first_mesh.halfedge(vd2, vd)), first_mesh);

						count++;

						K::Point_3 point = first_mesh.point(new_vd);
						if (!cut_plane.has_on(point)) {
							first_mesh.point(new_vd) = cut_plane.projection(point);
						}
						breakOut = true;
						break;
					}
				} while (vbegin != done);

			}
			if (breakOut) {
				breakOut = false;
			}
			else {
				break;
			}
		}
		
		first_mesh.collect_garbage();

		//remove triangles where two edges are on the boundary 
		for (face_descriptor fd : first_mesh.faces()) {
			halfedge_descriptor hd = first_mesh.halfedge(fd);
			halfedge_descriptor current = hd;
			bool isOnBoundary = false;
			do {
				if (first_mesh.null_face() == first_mesh.face(first_mesh.opposite(current)))
				{
					if (isOnBoundary) {
						//first_mesh.remove_face(fd);
						CGAL::Euler::remove_face(current, first_mesh);
						break;
					}
					else {
						isOnBoundary = true;
					}
				}
				current = first_mesh.next(current);
			} while (current != hd);
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
	Mesh boolean_interface::removeNeedleTriangles(Mesh mesh, float threshold, K::Plane_3 cut_plane) {
		for (face_descriptor fd : mesh.faces()) {
			
			float aspectRatio = boolean_interface::triangleAspectRatio(fd, mesh);

			halfedge_descriptor current = mesh.halfedge(fd);
			K::Point_3 p1 = mesh.point(mesh.target(current));
			//current = mesh.next(current);
			//K::Point_3 p2 = mesh.point(mesh.target(current));
			//current = mesh.next(current);
			//K::Point_3 p3 = mesh.point(mesh.target(current));
			K::Point_3 p2 = mesh.point(mesh.source(current));

			bool onPlane = (cut_plane.has_on(p1) && cut_plane.has_on(p2));

			std::cout << "aspectRatio: " << aspectRatio << std::endl;

			if (std::abs(aspectRatio) > threshold && !onPlane) {

				CGAL::Euler::join_face(mesh.halfedge(fd), mesh);
				std::cout << "removed face\n";

			}
			
		}
		mesh.collect_garbage();
		std::cout << "Needle triangles removed\n";
		std::ofstream output("../out/3D/removed_needleTris.off");
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

void __stdcall _approximate(ae::boolean_interface* g, char* first)
{
	g->approximate(std::string(first));
}

void __stdcall _simplify(ae::boolean_interface* g, char* first, double stop_ratio)
{
	g->simplify(std::string(first), stop_ratio);
}