#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

#include "tetScan.h"

using namespace PathInt;

tetScan::tetScan(Mesh::mesh & _msh, bool _filled) : refMsh(_msh), filled(_filled)
    {
    const int TET = refMsh.getNbTets();
    tetData.resize(TET);
    for (int t = 0; t < TET; t++)
        {
        const std::vector<int> &ind = refMsh.tet[t].ind;
        TetData &d = tetData[t];
        d.p0 = refMsh.getNode_p(ind[0]);
        Eigen::Matrix3d J;
        for (int k = 1; k < Tetra::N; k++)
            { J.col(k-1) = refMsh.getNode_p(ind[k]) - d.p0; }
        // rows of J^-1 are the gradients of lambda_1, lambda_2, lambda_3
        const Eigen::Matrix3d Jinv = J.inverse();
        d.grad.bottomRows<3>() = Jinv;
        d.grad.row(0) = -Jinv.colwise().sum();

        d.xmin = d.xmax = d.p0.x();
        d.ymin = d.ymax = d.p0.y();
        for (int k = 1; k < Tetra::N; k++)
            {
            const Eigen::Vector3d p = refMsh.getNode_p(ind[k]);
            d.xmin = std::min(d.xmin, p.x());
            d.xmax = std::max(d.xmax, p.x());
            d.ymin = std::min(d.ymin, p.y());
            d.ymax = std::max(d.ymax, p.y());
            }
        } // endfor t
    }

bool tetScan::chord(int t, double x, double y, Chord &c) const
    {
    const TetData &d = tetData[t];
    const double dx = x - d.p0.x();
    const double dy = y - d.p0.y();
    // along the line, lambda_i = a + b*w with w = z - p0.z
    double w_in = -std::numeric_limits<double>::infinity();
    double w_out = +std::numeric_limits<double>::infinity();
    for (int i = 0; i < Tetra::N; i++)
        {
        const double a = (i == 0 ? 1.0 : 0.0) + d.grad(i, 0)*dx + d.grad(i, 1)*dy;
        const double b = d.grad(i, 2);
        if (std::abs(b) <= GRAD_TOL*d.grad.row(i).norm())
            { // the line is parallel to the face opposite to vertex i
            if (a < -BARY_TOL) return false;
            continue;
            }
        const double w = -a/b;
        if (b > 0)
            { w_in = std::max(w_in, w); }
        else
            { w_out = std::min(w_out, w); }
        }
    if (w_in > w_out) return false;
    c = Chord{d.p0.z() + w_in, d.p0.z() + w_out, t};
    return true;
    }

Eigen::Vector3d tetScan::interpolate_u(int t, double x, double y, double z) const
    {
    const TetData &d = tetData[t];
    const Eigen::Vector3d dX(x - d.p0.x(), y - d.p0.y(), z - d.p0.z());
    const std::vector<int> &ind = refMsh.tet[t].ind;
    Eigen::Vector3d u = Eigen::Vector3d::Zero();
    for (int i = 0; i < Tetra::N; i++)
        {
        const double lambda = (i == 0 ? 1.0 : 0.0) + d.grad.row(i).dot(dX);
        u += lambda*refMsh.getNode_u(ind[i]);
        }
    return u;
    }

void tetScan::run(Fem2d &fem2d, const Settings &mySets) const
    {
    std::vector<Node2d> &node = fem2d.getNodes();
    const int NOD = node.size();
    const int N = static_cast<int>(std::lround(std::sqrt(static_cast<double>(NOD))));
    const double pixel = fem2d.pixel_size;
    const double x0 = node[0].p[0];
    const double y0 = node[0].p[1];

    // candidate tetrahedrons of each node (compressed sparse row): the nodes of the grid inside
    // the bounding box of the projection of the tetrahedron. Node (ix, iy) has index ix*N + iy.
    auto index_range = [N, pixel](double lo, double hi, double origin, int &i0, int &i1)
        {
        i0 = std::max(0, static_cast<int>(std::floor((lo - origin)/pixel)));
        i1 = std::min(N - 1, static_cast<int>(std::ceil((hi - origin)/pixel)));
        };
    const int TET = tetData.size();
    std::vector<int> first(NOD + 1, 0);
    for (int t = 0; t < TET; t++)
        {
        int ix0, ix1, iy0, iy1;
        index_range(tetData[t].xmin, tetData[t].xmax, x0, ix0, ix1);
        index_range(tetData[t].ymin, tetData[t].ymax, y0, iy0, iy1);
        for (int ix = ix0; ix <= ix1; ix++)
            for (int iy = iy0; iy <= iy1; iy++)
                { first[ix*N + iy + 1]++; }
        }
    for (int n = 0; n < NOD; n++)
        { first[n + 1] += first[n]; }
    std::vector<int> candidates(first[NOD]);
    std::vector<int> next(first.begin(), first.end() - 1);
    for (int t = 0; t < TET; t++)
        {
        int ix0, ix1, iy0, iy1;
        index_range(tetData[t].xmin, tetData[t].xmax, x0, ix0, ix1);
        index_range(tetData[t].ymin, tetData[t].ymax, y0, iy0, iy1);
        for (int ix = ix0; ix <= ix1; ix++)
            for (int iy = iy0; iy <= iy1; iy++)
                { candidates[next[ix*N + iy]++] = t; }
        }

    const double tol_z = Z_TOL*refMsh.l[2];

    #pragma omp parallel
        {
        std::vector<Chord> chords;

        // the cost of a node depends on the number of tetrahedrons above it: dynamic scheduling
        #pragma omp for schedule(dynamic, 256)
        for (int n = 0; n < NOD; n++)
            {
            Node2d &node2d = node[n];
            const double x = node2d.p[0];
            const double y = node2d.p[1];

            bool touched = false;
            chords.clear();
            for (int k = first[n]; k < first[n + 1]; k++)
                {
                Chord c;
                if (chord(candidates[k], x, y, c))
                    {
                    touched = true;
                    if (c.z_out - c.z_in > tol_z) chords.push_back(c);
                    }
                }
            std::sort(chords.begin(), chords.end(), [](const Chord &a, const Chord &b)
                      { return std::tie(a.z_in, a.z_out, a.tet) < std::tie(b.z_in, b.z_out, b.tet); });

            double path_length = 0.;
            Eigen::Vector3d M_integral = Eigen::Vector3d::Zero();
            double arg_plus = 0, arg_minus = 0;
            double z_first = 0., z_last = -std::numeric_limits<double>::infinity();
            for (const Chord &c : chords)
                {
                // a chord overlapping the previous one is a duplicate: the line goes along a
                // vertical face shared by two tetrahedrons
                if (c.z_in < z_last - tol_z) continue;
                if (z_last == -std::numeric_limits<double>::infinity()) z_first = c.z_in;
                z_last = c.z_out;

                const Tetra::Tet &tet_ = refMsh.tet[c.tet];
                const double J = mySets.paramTetra[tet_.idxPrm].J;
                const double inm = 1e+9;
                const double lbp = mySets.paramTetra[tet_.idxPrm].lbp*inm;
                const double lbm = mySets.paramTetra[tet_.idxPrm].lbm*inm;

                const Eigen::Vector3d u_in = interpolate_u(c.tet, x, y, c.z_in);
                const Eigen::Vector3d u_out = interpolate_u(c.tet, x, y, c.z_out);
                const double delta = c.z_out - c.z_in;
                if (!filled) path_length += delta;
                M_integral += 0.5*(u_in + u_out)*delta*J/mu0;
                const double scm = 0.5*(u_in.z() + u_out.z());  // in [-1, 1], the beam is along z
                arg_plus  += 0.5*(lbp + lbm + (lbp - lbm)*scm)*delta;  // sigma=+1
                arg_minus += 0.5*(lbp + lbm - (lbp - lbm)*scm)*delta;  // sigma=-1
                } // endfor c
            if (filled && !chords.empty()) path_length = z_last - z_first;

            const double Iph_plus  = exp(-arg_plus);  // sigma=+1
            const double Iph_minus = exp(-arg_minus); // sigma=-1
            node2d.flag_inside = touched ? 1 : 0;
            node2d.path_length = path_length;
            node2d.Mx_integral = M_integral.x();
            node2d.My_integral = M_integral.y();
            node2d.Mz_integral = M_integral.z();
            node2d.contrast = (Iph_plus - Iph_minus)/(Iph_plus + Iph_minus);
            } // endfor n
        } // end omp parallel
    }
