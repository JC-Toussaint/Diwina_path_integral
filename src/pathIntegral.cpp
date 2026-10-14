#include "pathIntegral.h"

using namespace PathInt;

void pathIntegral::processNode(int nod, Fem2d& fem2d, const Settings& mySettings) const
    {
	Node2d& node2d = fem2d.getNode(nod);

	const double x = node2d.p[0];
	const double y = node2d.p[1];
	const CGAL::Point ori = CGAL::Point(x, y, +refMsh.c[2] - 2 * refMsh.l[2]);
	const CGAL::Point ext = CGAL::Point(x, y, +refMsh.c[2] + 2 * refMsh.l[2]);

	double path_length = 0.;
	double contrast = 0.;
	CGAL::Vector M_integral(0., 0., 0.);

	const CGAL::Ray ray_query(ori, ext);
	/** intersection points and their associated facets, sorted along the ray */
	thread_local std::vector<Hit> hits;

	int nb_intersect = search_intersections(ray_query, hits);

	if (nb_intersect)
	    {
		node2d.flag_inside = 1; // intersection
		if (!hits.empty())
			path_integrals(mySettings, ori, ext, hits, path_length, M_integral, contrast);
	    }
	node2d.path_length = path_length;
	node2d.Mx_integral = M_integral.cartesian(0);
	node2d.My_integral = M_integral.cartesian(1);
	node2d.Mz_integral = M_integral.cartesian(2);
	node2d.contrast = contrast;
    }

void pathIntegral::integrate_u(const std::vector<Hit> &hits, double &path_length, CGAL::Vector &u_integral) const
    {
    path_length = 0.;
    u_integral = CGAL::Vector (0, 0, 0);
    if (hits.empty()) return;
    int inside(0);
    size_t begin_prev = 0;
    size_t end_prev = cluster_end(hits, begin_prev);
    CGAL::Point point_prev = hits[begin_prev].point;
    CGAL::Vector u_prev = interpolate_u(hits[begin_prev].facet, point_prev);

    while (end_prev < hits.size()) {
        const size_t begin_curr = end_prev;
        const size_t end_curr = cluster_end(hits, begin_curr);
        const CGAL::Point &point_curr = hits[begin_curr].point;
        CGAL::Vector u = interpolate_u(hits[begin_curr].facet, point_curr);
        if (locate_tet(hits, begin_prev, end_prev, begin_curr, end_curr) == -1) {inside = 0;} else {inside = 1;}
        double delta = sqrt((point_curr-point_prev).squared_length());
        if (filled)
           path_length += delta;
        else
           path_length += delta*inside;

        u_integral += 0.5*(u_prev+u)*delta*inside; 
        begin_prev = begin_curr;
        end_prev = end_curr;
        point_prev = point_curr;
        u_prev = u;
        } // endwhile
    }

/**
The intersection points are sorted along the ray. If the ray goes through an edge or a vertex of
the mesh, all the facets sharing it are kept: they give (nearly) the same point, and are grouped
later on into a cluster by cluster_end().
*/
int pathIntegral::search_intersections(const CGAL::Ray &ray_query, std::vector<Hit> &hits) const
    {
    using Intersection = CGAL::Tree::Intersection_and_primitive_id<CGAL::Ray>::Type;
    thread_local std::vector<Intersection> intersections;

    hits.clear();
    intersections.clear();
    tree->all_intersections(ray_query, std::back_inserter(intersections));
    if (intersections.empty()) return 0;

    const CGAL::Point &ori = ray_query.source();
    const CGAL::Vector dir = ray_query.second_point() - ori;
// the ray passes through a triangular facet at a point of intersection called iPoint
    for (const Intersection &inter : intersections)
        {
        if (const CGAL::Point *point = boost::get<CGAL::Point>(&inter.first))
            { // intersection object is a point
            const int facetIdx = std::distance(triangles.cbegin(), inter.second);
            hits.push_back(Hit{*point, facets[facetIdx], facetIdx, dir * (*point - ori), 0});
            }
        }//endfor

    // sort along the ray
    std::sort(hits.begin(), hits.end(),
              [](const Hit &a, const Hit &b) { return (a.t < b.t) || (a.t == b.t && a.facetIdx < b.facetIdx); });

    // hits closer than tol_t along the ray are the same point (t is scaled by the ray length squared)
    const double tol_t = CLUSTER_TOL * dir.squared_length();
    size_t first = 0;
    for (size_t k = 1; k < hits.size(); ++k)
        {
        if (hits[k].t - hits[first].t <= tol_t)
            { hits[k].cluster = hits[first].cluster; }
        else
            {
            hits[k].cluster = hits[first].cluster + 1;
            first = k;
            }
        }

    return intersections.size();
    }

size_t pathIntegral::cluster_end(const std::vector<Hit> &hits, size_t begin) const
    {
    size_t end = begin + 1;
    while (end < hits.size() && hits[end].cluster == hits[begin].cluster)
        { ++end; }
    return end;
    }

double pathIntegral::min_barycentric(int t, const CGAL::Point &p) const
    {
    const std::vector<int> &ind = refMsh.tet[t].ind;
    const CGAL::Point *q[Tetra::N] = {&pt[ind[0]], &pt[ind[1]], &pt[ind[2]], &pt[ind[3]]};
    // six times the signed volume of the tetrahedron (a, b, c, d)
    auto vol = [](const CGAL::Point &a, const CGAL::Point &b, const CGAL::Point &c, const CGAL::Point &d)
        { return CGAL::cross_product(b - a, c - a) * (d - a); };

    const double V = vol(*q[0], *q[1], *q[2], *q[3]);
    double lambda_min = vol(p, *q[1], *q[2], *q[3]) / V;
    lambda_min = std::min(lambda_min, vol(*q[0], p, *q[2], *q[3]) / V);
    lambda_min = std::min(lambda_min, vol(*q[0], *q[1], p, *q[3]) / V);
    lambda_min = std::min(lambda_min, vol(*q[0], *q[1], *q[2], p) / V);
    return lambda_min;
    }

int pathIntegral::locate_tet(const std::vector<Hit> &hits, size_t a_begin, size_t a_end,
                             size_t b_begin, size_t b_end) const
    {
    const CGAL::Point mid = CGAL::midpoint(hits[a_begin].point, hits[b_begin].point);
    int best_tet = -1;
    double best_lambda = -BARY_TOL;
    for (size_t k : {a_begin, b_begin})
        {
        const size_t end = (k == a_begin) ? a_end : b_end;
        for (size_t h = k; h < end; ++h)
            {
            for (int t : facet_tets[hits[h].facetIdx])
                {
                if (t < 0 || t == best_tet) continue;
                const double lambda = min_barycentric(t, mid);
                if (lambda > best_lambda)
                    {
                    best_lambda = lambda;
                    best_tet = t;
                    }
                }
            }
        }
    return best_tet;
    }

void pathIntegral::path_integrals(Settings const &mySets, const CGAL::Point &ori, const CGAL::Point &ext, const std::vector<Hit> &hits, double &path_length, CGAL::Vector &M_integral, double &contrast) const
    {
	CGAL::Vector uk = ext-ori;
	uk /= sqrt(uk.squared_length());
	path_length = 0.;
	M_integral = CGAL::Vector (0, 0, 0);
	int inside = 0;
	double J = 0;
	double lbp = 0, lbm = 0;
	double arg_plus = 0, arg_minus = 0;
    size_t begin_prev = 0;
    size_t end_prev = cluster_end(hits, begin_prev);
    CGAL::Point   point_prev = hits[begin_prev].point;
    CGAL::Vector  u_prev = interpolate_u(hits[begin_prev].facet, point_prev);

    while (end_prev < hits.size()) {
            const size_t begin_curr = end_prev;
            const size_t end_curr = cluster_end(hits, begin_curr);
            const CGAL::Point &point_curr = hits[begin_curr].point;
            CGAL::Vector u = interpolate_u(hits[begin_curr].facet, point_curr);
            int tet_id = locate_tet(hits, begin_prev, end_prev, begin_curr, end_curr);

			if (tet_id == -1) {
				inside = 0;
				J   = 0.;
				lbp = 0.;
				lbm = 0.;
			}
			else {
				inside = 1;
                const Tetra::Tet &tet_ = refMsh.tet[tet_id];
                J = mySets.paramTetra[tet_.idxPrm].J;
                const double inm=1e+9;
                lbp = mySets.paramTetra[tet_.idxPrm].lbp *inm; // lbp et lbm lus correctement
                lbm = mySets.paramTetra[tet_.idxPrm].lbm *inm;
                }

		double delta = sqrt((point_curr-point_prev).squared_length());
        if (filled)
           path_length += delta;
        else
           path_length += delta*inside;
		M_integral += 0.5*(u_prev+u)*delta*inside*J/mu0;
		double scm = 0.5* uk*(u_prev+u); // in [-1, 1]
		arg_plus  +=  0.5*(lbp+lbm + (lbp-lbm)*scm)*delta*inside;  // sigma=+1
		arg_minus +=  0.5*(lbp+lbm - (lbp-lbm)*scm)*delta*inside; // sigma=-1
		begin_prev = begin_curr;
		end_prev = end_curr;
		point_prev = point_curr;
		u_prev = u;
		} // endwhile

	double Iph_plus  = exp(-arg_plus ); // sigma=+1
	double Iph_minus = exp(-arg_minus); // sigma=-1
	assert(Iph_plus <1.);
	assert(Iph_minus<1.);
	contrast = (Iph_plus-Iph_minus)/(Iph_plus+Iph_minus);
    }

    /**
    Populate a vector of CGAL triangles from the tetrahedral mesh, sorted by their triplet of node
    indices, and the vector of the corresponding triplets.
    */
    void pathIntegral::populate_triangles(const std::vector<Tetra::Tet> &tet)
        {
        facets.clear();
        facets.reserve(4*tet.size());
        for (const auto &tet_ : tet)
            {
            const int ia = tet_.ind[0];
            const int ib = tet_.ind[1];
            const int ic = tet_.ind[2];
            const int id = tet_.ind[3];
            for (Triplet facet : {Triplet{ia, ib, ic}, Triplet{ia, ib, id}, Triplet{ia, ic, id}, Triplet{ib, ic, id}})
                {
                std::sort(facet.begin(), facet.end());
                facets.push_back(facet);
                }
            } // endfor tet_
        std::sort(facets.begin(), facets.end());
        facets.erase(std::unique(facets.begin(), facets.end()), facets.end());

        triangles.clear();
        triangles.reserve(facets.size());
        for (const Triplet &facet : facets)
            { triangles.emplace_back(pt[facet[0]], pt[facet[1]], pt[facet[2]]); }

        facet_tets.assign(facets.size(), {-1, -1});
        for (unsigned int t = 0; t < tet.size(); t++)
            {
            const std::vector<int> &ind = tet[t].ind;
            for (Triplet facet : {Triplet{ind[0], ind[1], ind[2]}, Triplet{ind[0], ind[1], ind[3]},
                                  Triplet{ind[0], ind[2], ind[3]}, Triplet{ind[1], ind[2], ind[3]}})
                {
                std::sort(facet.begin(), facet.end());
                const int f = std::distance(facets.begin(), std::lower_bound(facets.begin(), facets.end(), facet));
                std::array<int, 2> &owners = facet_tets[f];
                if (owners[0] == -1) owners[0] = t;
                else owners[1] = t;
                }
            } // endfor t
        }

/** computes Lagrange polynoms on facette */
std::array<double, Facette::N> pathIntegral::facet_lagrange_polynomials(const Triplet &facet_indices, const CGAL::Point &point_curr) const
    {
    std::array<double, Facette::N> alpha;
    CGAL::Vector vec[Facette::N];

    for (int ie=0; ie<Facette::N; ie++)
        { vec[ie] = CGAL::Vector(point_curr, pt[facet_indices[ie]]); }
             
    double subArea[Facette::N];
    double area(0); //pour calculer la surface des petits triangles
    for (int ie=0; ie<Facette::N; ie++)
        {
        CGAL::Vector crossProduct = CGAL::cross_product(vec[(ie+1)%Facette::N], vec[(ie+2)%Facette::N]);
        subArea[ie] = 0.5*sqrt(crossProduct.squared_length());
        area += std::abs(subArea[ie]); //aire totale du triangle
        }
    for (int ie=0; ie<Facette::N; ie++)
        { alpha[ie] = subArea[ie]/area; }   
    return alpha;              
    } 

/** computes magnetization interpolation */
CGAL::Vector pathIntegral::interpolate_u(const Triplet &facet_indices, const CGAL::Point &point) const
    {
    std::array<double, Facette::N> alpha = facet_lagrange_polynomials(facet_indices, point);
    CGAL::Vector u(0., 0., 0.);
    for (int ie=0; ie<Facette::N; ie++)
        {
        int i = facet_indices[ie];
        Eigen::Vector3d u_nod= refMsh.getNode_u(i);
        u += alpha[ie]*CGAL::Vector(u_nod[0], u_nod[1], u_nod[2]); 
        }
    return u;
    }
