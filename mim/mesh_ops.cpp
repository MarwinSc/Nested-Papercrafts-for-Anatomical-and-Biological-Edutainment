#include "mesh_ops.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>
#include <limits>
#include "mesh_ops.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polygon_2<K> Polygon_2;

Polyhedron get_polyhedron(char* mesh)
{
	std::ifstream input(mesh);
	Polyhedron poly;
	if (!input || !(input >> poly) || poly.empty()
		|| !CGAL::is_triangle_mesh(poly))
	{
		throw std::exception("not a valid input file.");
	}

	input.close();
	return poly;
}

bool __stdcall meshB_inside_of_meshA(char* meshA, char* meshB)
{
	Polyhedron polyA = get_polyhedron(meshA);
	CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(polyA);
	Polyhedron polyB = get_polyhedron(meshB);

	int nb_inside = 0;
	int nb_boundary = 0;
	int point_count = 0;
	auto iterator = polyB.points_begin();
	while (iterator != polyB.points_end())
	{
		CGAL::Bounded_side res = inside(*iterator);
		if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
		if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }

		iterator++;
		point_count++;
	}

	std::vector<std::pair<size_t, size_t>> intersections(1);
	std::vector<Polyhedron> triangleMeshRange = { polyA, polyB };
	auto iti = CGAL::Polygon_mesh_processing::intersecting_meshes<std::vector<Polyhedron>, std::vector<std::pair<size_t, size_t>>::iterator>(triangleMeshRange, intersections.begin());
	bool intersects = false;
	if (iti->first == 1 || iti->second == 1)
	{
		intersects = true;
	}

	std::cout << "Checked " << point_count << " points" << std::endl;
	std::cout << "  " << nb_inside << " points inside " << std::endl;
	std::cout << "  " << nb_boundary << " points on boundary " << std::endl;
	std::cout << "  " << point_count - nb_inside - nb_boundary << " points outside " << std::endl;
	std::cout << "  " << (intersects ? "meshes intersect" : "meshes don't intersect") << std::endl;

	return !intersects && point_count - nb_inside == 0;
}

bool __stdcall convex_hull_of_mesh(char* mesh, char* hull_file) {
	Polyhedron poly;
	try
	{
		poly = get_polyhedron(mesh);
	}
	catch (...)
	{
		std::cout << "failed to get polyhedron" << std::endl;
		return false;
	}
	Polyhedron hull;

	try
	{
		CGAL::convex_hull_3(poly.points_begin(), poly.points_end(), hull);
	}
	catch (...)
	{
		std::cout << "failed to compute convex hull" << std::endl;
		return false;
	}

	try
	{
		std::ofstream off(hull_file);
		CGAL::set_ascii_mode(off);
		off << hull;
	}
	catch (...)
	{
		std::cout << "failed to save" << std::endl;
		return false;
	}

	return true;
}