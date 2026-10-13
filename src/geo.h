#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Random.h>

namespace CGAL
    {
    typedef CGAL::Simple_cartesian<double> K;
    //typedef CGAL::Exact_predicates_exact_constructions_kernel K;

    typedef K::FT FT;
    typedef K::Ray_3 Ray;
    typedef K::Line_3 Line;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;
    typedef K::Triangle_3 Triangle;

    typedef std::list<Triangle>::iterator Iterator;
    typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
    typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

    typedef Tree::Object_and_primitive_id Object_and_primitive_id;
    typedef Tree::Primitive_id Primitive_id;

    typedef CGAL::Search_traits_3<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Search_Tree;

    using Facet2Pt   = std::map<int, Point>;
    using Tet2LastPt = std::map<int, Point>;
    using Pt2Facet   = std::map<Point, int>;
    using Pt2Tet     = std::map<Point, int>;
    using Point_List = std::list<Point>;
    using DistanceFromSource_to_Point = std::map<double, Point>;
}

