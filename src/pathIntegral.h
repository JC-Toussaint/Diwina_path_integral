#ifndef pathIntegral_h
#define pathIntegral_h

#include <array>
#include <memory>
#include <vector>

#include "fem.h"
#include "fem2d.h"
#include "tetra.h"

namespace PathInt
    {
    using Triplet         = std::array<int,   3>;
    using Quadruplet      = std::array<int,   4>;

    /** intersection point of a ray with a facet, t is the abscissa of the point along the ray */
    struct Hit
        {
        CGAL::Point point;
        Triplet facet;
        int facetIdx;  /**< index of the facet in facets */
        double t;
        int cluster;   /**< hits of the same cluster are the same intersection point */
        };

    /** relative tolerance, with respect to the ray length, under which intersection points are
    considered as the same point (the ray goes through an edge or a vertex of the mesh) */
    const double CLUSTER_TOL = 1e-9;

    /** tolerance on barycentric coordinates to decide that a point is inside a tetrahedron */
    const double BARY_TOL = 1e-6;

    class pathIntegral
        {
        public:

    /** constructor
    populate CGAL triangles and the facet (node indices) associated to each of them, and the
    tetrahedrons sharing each facet
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

                // constructs AABB tree
                tree = std::make_unique<CGAL::Tree>(triangles.cbegin(), triangles.cend());
                tree->build();
                }

    void processNode(int nod, Fem2d& fem2d, const Settings& mySettings) const;

    void integrate_u(const std::vector<Hit> &hits, double &path_length, CGAL::Vector &u_integral) const;

    /** fill hits with the ray-facets intersection points, sorted along the ray, return the number of intersected facets.
    When the ray goes through an edge or a vertex, several hits share (nearly) the same point:
    they are grouped into the same cluster. */
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

            /** facet_tets[i] are the (at most two) tetrahedrons sharing facets[i], -1 if none */
            std::vector<std::array<int, 2>> facet_tets;

    /**
    Populate a vector of CGAL triangles from the tetrahedral mesh, the vector of the
    corresponding triplets of node indices, and the tetrahedrons sharing each of them.
    */
        void populate_triangles(const std::vector<Tetra::Tet> &tet);

        /** return the index of the end of the cluster of hits starting at begin */
        size_t cluster_end(const std::vector<Hit> &hits, size_t begin) const;

        /** smallest barycentric coordinate of point p in tetrahedron t */
        double min_barycentric(int t, const CGAL::Point &p) const;

        /** return the tetrahedron containing the segment between the clusters of hits [a_begin, a_end)
        and [b_begin, b_end), else -1. The candidates are the tetrahedrons sharing the facets of
        the two clusters, the one containing the middle of the segment is chosen. */
        int locate_tet(const std::vector<Hit> &hits, size_t a_begin, size_t a_end,
                       size_t b_begin, size_t b_end) const;

        std::array<double, Facette::N> facet_lagrange_polynomials(const Triplet &facet_indices, const CGAL::Point &point_curr) const;

        CGAL::Vector interpolate_u(const Triplet &facet_indices, const CGAL::Point &point) const;
        }; // end class
    } // end namespace PathInt
#endif
