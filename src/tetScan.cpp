#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>

#include <omp.h>

#include "tetScan.h"

using namespace PathInt;

tetScan::tetScan(Mesh::mesh & _msh, bool _filled) : refMsh(_msh), filled(_filled)
    {
    nbTets = refMsh.getNbTets();
    tetData.reset(new TetData[nbTets]);
    #pragma omp parallel for
    for (int t = 0; t < nbTets; t++)
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

bool tetScan::chord(int t, double x, double y, double tol_z, Chord &c) const
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
    // a line touching the tetrahedron (edge or vertex) gives w_in = w_out up to rounding errors
    if (w_in > w_out + tol_z) return false;
    w_out = std::max(w_in, w_out);
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
    // Each grid row ix is handled by a single thread, so that no atomic operation is needed. For
    // each node, the candidates are sorted by increasing tetrahedron index.
    const int TET = nbTets;
    std::vector<std::array<int, 4>> ranges(TET);  // ix0, ix1, iy0, iy1

    // first step: list of the tetrahedrons overlapping each row, built with per thread counters
    std::vector<int> rowFirst(N + 1);
    std::vector<int> rowTets;
    std::vector<int> rowCount;
    #pragma omp parallel
        {
        #pragma omp single
        { rowCount.assign(static_cast<size_t>(omp_get_num_threads())*N, 0); }
        int *count = rowCount.data() + static_cast<size_t>(omp_get_thread_num())*N;

        #pragma omp for schedule(static)
        for (int t = 0; t < TET; t++)
            {
            auto &[ix0, ix1, iy0, iy1] = ranges[t];
            index_range(tetData[t].xmin, tetData[t].xmax, x0, ix0, ix1);
            index_range(tetData[t].ymin, tetData[t].ymax, y0, iy0, iy1);
            for (int ix = ix0; ix <= ix1; ix++)
                { count[ix]++; }
            }

        #pragma omp single
            { // counters are turned into the first position of each thread in each row
            const int nbThreads = rowCount.size()/N;
            int pos = 0;
            for (int ix = 0; ix < N; ix++)
                {
                rowFirst[ix] = pos;
                for (int th = 0; th < nbThreads; th++)
                    {
                    const int c = rowCount[static_cast<size_t>(th)*N + ix];
                    rowCount[static_cast<size_t>(th)*N + ix] = pos;
                    pos += c;
                    }
                }
            rowFirst[N] = pos;
            rowTets.resize(pos);
            }

        // same static distribution as the loop above: each thread fills the positions it counted
        #pragma omp for schedule(static)
        for (int t = 0; t < TET; t++)
            {
            for (int ix = ranges[t][0]; ix <= ranges[t][1]; ix++)
                { rowTets[count[ix]++] = t; }
            }
        } // end omp parallel

    // second step: number of candidates of each node, then first position of each node
    std::vector<int> first(NOD + 1);
    std::vector<int> rowSum(N + 1, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int ix = 0; ix < N; ix++)
        {
        int *count = first.data() + static_cast<size_t>(ix)*N;
        std::fill(count, count + N, 0);
        for (int k = rowFirst[ix]; k < rowFirst[ix + 1]; k++)
            {
            const auto &r = ranges[rowTets[k]];
            for (int iy = r[2]; iy <= r[3]; iy++)
                { count[iy]++; }
            }
        rowSum[ix + 1] = std::accumulate(count, count + N, 0);
        }
    for (int ix = 0; ix < N; ix++)
        { rowSum[ix + 1] += rowSum[ix]; }
    #pragma omp parallel for
    for (int ix = 0; ix < N; ix++)
        {
        int pos = rowSum[ix];
        for (int iy = 0; iy < N; iy++)
            {
            const int c = first[ix*N + iy];
            first[ix*N + iy] = pos;
            pos += c;
            }
        }
    first[NOD] = rowSum[N];

    // third step: candidates, not value initialized, filled row by row
    std::unique_ptr<int[]> candidates(new int[first[NOD]]);
    #pragma omp parallel
        {
        std::vector<int> next(N);
        #pragma omp for schedule(dynamic, 1)
        for (int ix = 0; ix < N; ix++)
            {
            std::copy(first.begin() + static_cast<size_t>(ix)*N, first.begin() + static_cast<size_t>(ix + 1)*N, next.begin());
            for (int k = rowFirst[ix]; k < rowFirst[ix + 1]; k++)
                {
                const int t = rowTets[k];
                const auto &r = ranges[t];
                for (int iy = r[2]; iy <= r[3]; iy++)
                    { candidates[next[iy]++] = t; }
                }
            }
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
                if (chord(candidates[k], x, y, tol_z, c))
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
