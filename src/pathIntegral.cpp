#include "pathIntegral.h"

using namespace PathInt;

void pathIntegral::processNode(int nod, Fem2d& fem2d, const Settings& mySettings)
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
	/** intersection point to its associated facet */
	std::map<CGAL::Point, Triplet>  iPoint_to_facet;

	int nb_intersect = search_intersections(ray_query, iPoint_to_facet);

	if (nb_intersect)
	    {
		node2d.flag_inside = 1; // intersection
		path_integrals(mySettings, ori, ext, iPoint_to_facet, path_length, M_integral, contrast);
	    }
	node2d.path_length = path_length;
	node2d.Mx_integral = M_integral.cartesian(0);
	node2d.My_integral = M_integral.cartesian(1);
	node2d.Mz_integral = M_integral.cartesian(2);
	node2d.contrast = contrast;
    }

void pathIntegral::integrate_u(CGAL::Point &ori, CGAL::Point &ext, std::map<CGAL::Point, Triplet> &iPoint_to_facet, double &path_length, CGAL::Vector &u_integral)
    {
	auto compare = [ori, ext](const auto& lhs, const auto& rhs) -> bool
	                         { return ((ext-ori) * (lhs-ori)) < ((ext-ori) * (rhs-ori)); };

    std::map<CGAL::Point, Triplet, decltype(compare)> sorted_iPoint_to_facet(iPoint_to_facet.begin(), iPoint_to_facet.end(), compare);
    path_length = 0.;
    u_integral = CGAL::Vector (0, 0, 0);
    int inside(0);
    auto it = sorted_iPoint_to_facet.begin();
    const auto& entry = *it;
    CGAL::Point point_prev = entry.first;
    Triplet facet_prev = entry.second;
    CGAL::Vector u_prev = interpolate_u(facet_prev, point_prev);

    if (it != sorted_iPoint_to_facet.end()) {
        ++it; // Move the iterator to the second element
        for (; it != sorted_iPoint_to_facet.end(); ++it) {
            const auto& entry = *it;
            CGAL::Point point_curr = entry.first;
            Triplet facet_curr = entry.second;
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
            } // endfor it
       }//endif      
    }

int pathIntegral::search_intersections(const CGAL::Ray &ray_query, std::map<CGAL::Point, Triplet> &iPoint_to_facet)
    {
    if (!tree->do_intersect(ray_query)) return 0;

    // computes all intersections with segment query (as pairs object - primitive_id)
    std::list<CGAL::Object_and_primitive_id> intersections;
    tree->all_intersections(ray_query, std::back_inserter(intersections));

// the ray passes through a triangular facet at a point of intersection called iPoint
    for (std::list<CGAL::Object_and_primitive_id>::const_iterator it = intersections.begin(); it != intersections.end(); ++it)
        {
        CGAL::Object_and_primitive_id op = *it;
        CGAL::Object object = op.first;
        std::list<CGAL::Triangle>::iterator itv=op.second;

        CGAL::Point point;
        if (CGAL::assign(point, object)){ // intersection object is a point
           Triplet facet = {point_to_id[itv->vertex(0)], point_to_id[itv->vertex(1)], point_to_id[itv->vertex(2)]};
           std::sort(facet.begin(), facet.end());
           iPoint_to_facet[point] = facet;
           }
	    }//endfor
    return tree->number_of_intersected_primitives(ray_query);
    }

void pathIntegral::path_integrals(Settings const &mySets, const CGAL::Point &ori, const CGAL::Point &ext, std::map<CGAL::Point, Triplet> &iPoint_to_facet, double &path_length, CGAL::Vector &M_integral, double &contrast)
    {
	auto compare = [ori, ext](const auto& lhs, const auto& rhs) -> bool
	                         { return ((ext-ori) * (lhs-ori)) < ((ext-ori) * (rhs-ori)); };

    std::map<CGAL::Point, Triplet, decltype(compare)> sorted_iPoint_to_facet(iPoint_to_facet.begin(), iPoint_to_facet.end(), compare);
	CGAL::Vector uk = ext-ori;
	uk /= sqrt(uk.squared_length());
	path_length = 0.;
	M_integral = CGAL::Vector (0, 0, 0);
	int inside = 0;
	double J = 0;
	double lbp = 0, lbm = 0;
	double arg_plus = 0, arg_minus = 0;
    auto it = sorted_iPoint_to_facet.begin();
    const auto& entry = *it;
    CGAL::Point   point_prev = entry.first;
    Triplet facet_prev = entry.second;
    CGAL::Vector  u_prev = interpolate_u(facet_prev, point_prev);

    if (it != sorted_iPoint_to_facet.end()) {
        ++it; // Move the iterator to the second element
        for (; it != sorted_iPoint_to_facet.end(); ++it) {
            const auto& entry = *it;
            CGAL::Point   point_curr = entry.first;
            Triplet facet_curr = entry.second;
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
                Tetra::Tet &tet_ = refMsh.tet[tet_id];
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
		} // endfor it
	}//endif

	double Iph_plus  = exp(-arg_plus ); // sigma=+1
	double Iph_minus = exp(-arg_minus); // sigma=-1
	assert(Iph_plus <1.);
	assert(Iph_minus<1.);
	contrast = (Iph_plus-Iph_minus)/(Iph_plus+Iph_minus);
    }

    /**
    Populate a list of CGAL triangles from the tetrahedral mesh. Associate each triplet of points 
    that forms a triangle with the number of the tetrahedron from which it originates.
    */
    void pathIntegral::populate_triangles(std::vector<Tetra::Tet> &tet)
        {
        std::map<Triplet, CGAL::Triangle> facet_to_triangle;
        for (auto tet_ : tet)
            {
            int ia = tet_.ind[0];
            int ib = tet_.ind[1];
            int ic = tet_.ind[2];
            int id = tet_.ind[3];
                {  
                Triplet facet = {ia, ib, ic};
                std::sort(facet.begin(), facet.end());
                facet_to_triangle[facet] = CGAL::Triangle(pt[facet[0]], pt[facet[1]], pt[facet[2]]);
                }
                {  
                Triplet facet = {ia, ib, id};    
                std::sort(facet.begin(), facet.end());
                facet_to_triangle[facet] = CGAL::Triangle(pt[facet[0]], pt[facet[1]], pt[facet[2]]);
                }
                {  
                Triplet facet = {ia, ic, id};    
                std::sort(facet.begin(), facet.end());
                facet_to_triangle[facet] = CGAL::Triangle(pt[facet[0]], pt[facet[1]], pt[facet[2]]);
                }
                {  
                Triplet facet = {ib, ic, id};    
                std::sort(facet.begin(), facet.end());
                facet_to_triangle[facet] = CGAL::Triangle(pt[facet[0]], pt[facet[1]], pt[facet[2]]);
                }
            } // endfor tet_
        for (const auto& entry : facet_to_triangle) { triangles.push_back(entry.second); }
        }

    void pathIntegral::populate_quads(std::vector<Tetra::Tet> &tet)
        {  
        for (unsigned int t=0; t<tet.size(); t++)
            {
            const Tetra::Tet &tet_ = tet[t];
            Quadruplet quad ={tet_.ind[0], tet_.ind[1], tet_.ind[2], tet_.ind[3]};
            std::sort(quad.begin(), quad.end());
            quad_to_tetId[quad] = t;       
            } // endfor t                
        }

/** computes Lagrange polynoms on facette */
std::vector<double> pathIntegral::facet_lagrange_polynomials(const Triplet &facet_indices, const CGAL::Point &point_curr)
    {
    std::vector<double> alpha(Facette::N);
    CGAL::Vector vec[Facette::N];

    for (int ie=0; ie<Facette::N; ie++)
        {
        int i = facet_indices[ie];
        Eigen::Vector3d xyz = refMsh.getNode_p(i);
        CGAL::Point tmp(xyz[0], xyz[1], xyz[2]);
        vec[ie] = CGAL::Vector(point_curr, tmp);   
        }
             
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
CGAL::Vector pathIntegral::interpolate_u(const Triplet &facet_indices, const CGAL::Point &point)
    {
    std::vector<double> alpha = facet_lagrange_polynomials(facet_indices, point);
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
int pathIntegral::identify_tet(Triplet &facet_curr, Triplet &facet_prev)
    {
    auto compare = [](const int a,const int b) -> bool { return a < b; };
    std::set<int, decltype(compare)> setOfQuads(compare);
    setOfQuads.insert(facet_curr[0]); setOfQuads.insert(facet_curr[1]); setOfQuads.insert(facet_curr[2]);
    setOfQuads.insert(facet_prev[0]); setOfQuads.insert(facet_prev[1]); setOfQuads.insert(facet_prev[2]);
        
    Quadruplet quad{};
    // Ensures that we only copy as many elements as the size of quad
    std::copy_n(setOfQuads.begin(), quad.size(), quad.begin());
 
    if (quad_to_tetId.find(quad) != quad_to_tetId.end())
       return quad_to_tetId[quad];
    else
       return -1;
    }               
