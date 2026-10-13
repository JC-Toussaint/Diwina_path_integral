#ifndef pathIntegral_h
#define pathIntegral_h

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

#include "fem.h"
#include "fem2d.h"
#include "tetra.h"

namespace PathInt
    {
    using Triplet         = std::array<int,   3>;
    using Quadruplet      = std::array<int,   4>;

    /** hash function for Quadruplet, to be used in an unordered_map */
    struct QuadrupletHash
        {
        std::size_t operator()(const Quadruplet &q) const noexcept
            {
            std::size_t h = 0;
            for (int i : q)
                { h ^= std::hash<int>{}(i) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); }
            return h;
            }
        };

    /** intersection point of a ray with a facet, t is the abscissa of the point along the ray */
    struct Hit
        {
        CGAL::Point point;
        Triplet facet;
        double t;
        };

    class pathIntegral
        {
        public:

    /** constructor
    populate CGAL triangles and the facet (node indices) associated to each of them
    populate quadruplet to tetrahedron dictionary
    build CGAL tree
    */
            pathIntegral(Mesh::mesh & _msh,bool _filled):refMsh(_msh),filled(_filled)
                {
                const int NOD = refMsh.getNbNodes();
                pt.resize(NOD);
                for (int i=0; i<NOD; i++)
                    {
                    Eigen::Vector3d coords = refMsh.getNode_p(i);
                    pt[i] = CGAL::Point(coords[0], coords[1], coords[2]);
                    }

                populate_triangles(refMsh.tet);
                populate_quads(refMsh.tet);

                // constructs AABB tree
                tree = std::make_unique<CGAL::Tree>(triangles.cbegin(), triangles.cend());
                tree->build();
                }

        /** return tetrahedron number if the two facets are adjacent else -1 */
        int identify_tet(const Triplet &facet_curr, const Triplet &facet_prev) const;

    void processNode(int nod, Fem2d& fem2d, const Settings& mySettings) const;

    void integrate_u(const std::vector<Hit> &hits, double &path_length, CGAL::Vector &u_integral) const;

    /** fill hits with the ray-facets intersection points, sorted along the ray, return the number of intersected facets */
    int search_intersections(const CGAL::Ray &ray_query, std::vector<Hit> &hits) const;

    void path_integrals(Settings const &mySets, const CGAL::Point &ori, const CGAL::Point &ext, const std::vector<Hit> &hits, double &path_length, CGAL::Vector &M_integral, double &contrast) const;

        private:
            /** ref to the mesh */
            Mesh::mesh & refMsh;

            /** if true the mesh is hollow, filled with non magnetic material */
            const bool filled;

            /** CGAL AABB-Tree */
            std::unique_ptr<CGAL::Tree> tree;

            /** Node coordinates stores into a vector of CGAL::Point */
            std::vector<CGAL::Point> pt;

            /** tetrahedron mesh splitted into triangles */
            std::vector<CGAL::Triangle> triangles;

            /** facets[i] are the sorted node indices of triangles[i] */
            std::vector<Triplet> facets;

            /** quadruplet to its tetrahedron number */
            std::unordered_map<Quadruplet, int, QuadrupletHash> quad_to_tetId;

    /**
    Populate a vector of CGAL triangles from the tetrahedral mesh, and the vector of the
    corresponding triplets of node indices.
    */
        void populate_triangles(const std::vector<Tetra::Tet> &tet);

        void populate_quads(const std::vector<Tetra::Tet> &tet);

        std::array<double, Facette::N> facet_lagrange_polynomials(const Triplet &facet_indices, const CGAL::Point &point_curr) const;

        CGAL::Vector interpolate_u(const Triplet &facet_indices, const CGAL::Point &point) const;
        }; // end class
    } // end namespace PathInt
#endif
