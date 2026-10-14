#ifndef tetScan_h
#define tetScan_h

#include <vector>

#include "fem.h"
#include "fem2d.h"
#include "tetra.h"

namespace PathInt
    {
    /**
    \class tetScan
    Computes the integrals along the beam by scanning the tetrahedrons rather than throwing rays.
    The beam is along z, so the integration lines are vertical lines through the nodes of the
    regular square sensor grid. Each tetrahedron is projected on the grid: the nodes inside the
    bounding box of its projection are its candidate nodes. For each node, the chord of the
    vertical line in each candidate tetrahedron is computed from the barycentric coordinates,
    which are linear along the line.
    */
    class tetScan
        {
        public:
            /** constructor: computes the barycentric gradients of all the tetrahedrons */
            tetScan(Mesh::mesh & _msh, bool _filled);

            /** compute path length, magnetization integrals and contrast on all the nodes of fem2d */
            void run(Fem2d &fem2d, const Settings &mySets) const;

        private:
            /** ref to the mesh */
            Mesh::mesh & refMsh;

            /** if true the mesh is hollow, filled with non magnetic material */
            const bool filled;

            /** chord of a vertical line in a tetrahedron, between z_in and z_out */
            struct Chord
                {
                double z_in, z_out;
                int tet;
                };

            /** data of a tetrahedron needed to compute the chords */
            struct TetData
                {
                Eigen::Vector3d p0;  /**< first vertex */
                /** gradients of the barycentric coordinates: lambda_i(X) = (i == 0) + grad.row(i).dot(X - p0) */
                Eigen::Matrix<double, Tetra::N, Nodes::DIM> grad;
                double xmin, xmax, ymin, ymax;  /**< bounding box of the projection in the xy plane */
                };

            std::vector<TetData> tetData;

            /** relative tolerance, with respect to the mesh height, under which chords are
            considered of zero length, or as overlapping */
            static constexpr double Z_TOL = 1e-9;

            /** relative tolerance on a barycentric gradient to consider it orthogonal to z */
            static constexpr double GRAD_TOL = 1e-12;

            /** tolerance on barycentric coordinates to decide that a point is inside a tetrahedron */
            static constexpr double BARY_TOL = 1e-10;

            /** compute the chord of the vertical line through (x, y) in tetrahedron t, return false if
            the line misses the tetrahedron */
            bool chord(int t, double x, double y, Chord &c) const;

            /** magnetization interpolated at the point (x, y, z) of tetrahedron t */
            Eigen::Vector3d interpolate_u(int t, double x, double y, double z) const;
        };
    } // end namespace PathInt
#endif
