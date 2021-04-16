#include "mesh_inside.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>

#include <vector>
#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

double max_coordinate(const Polyhedron& poly)
{
	double max_coord = -std::numeric_limits<double>::infinity();
	for (Polyhedron::Vertex_handle v : vertices(poly))
	{
		Point p = v->point();
		max_coord = (std::max)(max_coord, p.x());
		max_coord = (std::max)(max_coord, p.y());
		max_coord = (std::max)(max_coord, p.z());
	}
	return max_coord;
}

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