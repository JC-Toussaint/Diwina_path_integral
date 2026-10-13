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
    CGAL::Point point_prev = hits[0].point;
    Triplet facet_prev = hits[0].facet;
    CGAL::Vector u_prev = interpolate_u(facet_prev, point_prev);

    for (size_t k = 1; k < hits.size(); ++k) {
        const CGAL::Point &point_curr = hits[k].point;
        const Triplet &facet_curr = hits[k].facet;
        CGAL::Vector u = interpolate_u(facet_curr, point_curr);
        if (identify_tet(facet_curr, facet_prev) == -1) {inside = 0;} else {inside = 1;}
        double delta = sqrt((point_curr-point_prev).squared_length());
        if (filled)
           path_length += delta;
        else
           path_length += delta*inside;

        u_integral += 0.5*(u_prev+u)*delta*inside; 
        point_prev = point_curr;
        facet_prev = facet_curr;
        u_prev = u;
        } // endfor k
    }

/**
The intersection points are sorted along the ray. If the ray hits several facets at the same point
(edge or vertex), only one of them is kept: for identical points the last facet found, and for
distinct points with the same abscissa along the ray the lexicographically smallest point.
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
            const Triplet &facet = facets[std::distance(triangles.cbegin(), inter.second)];
            hits.push_back(Hit{*point, facet, dir * (*point - ori)});
            }
        }//endfor

    // identical points: keep the last one found
    std::stable_sort(hits.begin(), hits.end(),
                     [](const Hit &a, const Hit &b) { return a.point < b.point; });
    auto last = std::unique(hits.rbegin(), hits.rend(),
                            [](const Hit &a, const Hit &b) { return a.point == b.point; });
    hits.erase(hits.begin(), last.base());

    // sort along the ray: keep the first one for a given abscissa
    std::stable_sort(hits.begin(), hits.end(),
                     [](const Hit &a, const Hit &b) { return a.t < b.t; });
    hits.erase(std::unique(hits.begin(), hits.end(),
                           [](const Hit &a, const Hit &b) { return !(a.t < b.t) && !(b.t < a.t); }),
               hits.end());

    return intersections.size();
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
    CGAL::Point   point_prev = hits[0].point;
    Triplet facet_prev = hits[0].facet;
    CGAL::Vector  u_prev = interpolate_u(facet_prev, point_prev);

    for (size_t k = 1; k < hits.size(); ++k) {
            const CGAL::Point &point_curr = hits[k].point;
            const Triplet &facet_curr = hits[k].facet;
            CGAL::Vector u = interpolate_u(facet_curr, point_curr);
            int tet_id = identify_tet(facet_curr, facet_prev);

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
		point_prev = point_curr;
		facet_prev = facet_curr;
		u_prev = u;
		} // endfor k

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
        }

    void pathIntegral::populate_quads(const std::vector<Tetra::Tet> &tet)
        {  
        quad_to_tetId.reserve(tet.size());
        for (unsigned int t=0; t<tet.size(); t++)
            {
            const Tetra::Tet &tet_ = tet[t];
            Quadruplet quad ={tet_.ind[0], tet_.ind[1], tet_.ind[2], tet_.ind[3]};
            std::sort(quad.begin(), quad.end());
            quad_to_tetId[quad] = t;       
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


    /** return tetrahedron number if the two facets are adjacent else -1 
    * In the ideal case, knowing the triplets of 2 triangles sharing an edge, we deduce the tetrahedron containing these 2 triangles
         3
        /|\
       / | \
      /  2  \
     / .   . \
    0---------1
    
    One case, of which there are certainly many, that poses difficulties is when the ray crosses 2 triangles belonging to different tetrahedrons. 
    We store only the first 4 indices encountered in the quadruplet quad. 
         3         2
        / \       / \
       /   \     /   \
      /     \   /     \
     /       \ /       \
    0---------1--------4 
    */
int pathIntegral::identify_tet(const Triplet &facet_curr, const Triplet &facet_prev) const
    {
    // sorted union of the indices of the two facets, only the 4 smallest ones are kept
    std::array<int, 2*Facette::N> indices{facet_curr[0], facet_curr[1], facet_curr[2],
                                          facet_prev[0], facet_prev[1], facet_prev[2]};
    std::sort(indices.begin(), indices.end());
    auto last = std::unique(indices.begin(), indices.end());

    Quadruplet quad{};
    std::copy(indices.begin(), std::min(last, indices.begin() + quad.size()), quad.begin());

    auto it = quad_to_tetId.find(quad);
    return (it != quad_to_tetId.end()) ? it->second : -1;
    }
