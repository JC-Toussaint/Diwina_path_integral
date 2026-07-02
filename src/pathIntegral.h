#ifndef pathIntegral_h
#define pathIntegral_h

#include <map>
#include <set>

#include "fem.h"
#include "fem2d.h"
#include "tetra.h"

namespace PathInt
    {    
    using Triplet         = std::array<int,   3>;
    using Quadruplet      = std::array<int,   4>;

    class pathIntegral
        {
        public:
            
    /** constructor
    populate a vector of CGAL::Point coordinates and a dictionary numbering each mesh node from 
    its coordinates stored in a CGAL::Point class
    populate CGAL triangles and quadrangles
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
                    point_to_id[pt[i]]= i;
                    }
                
                populate_triangles(refMsh.tet);
                populate_quads(refMsh.tet);
                        
                // constructs AABB tree
                tree=new CGAL::Tree(triangles.begin(), triangles.end());
                }
            
            ~pathIntegral() { delete tree; }

        /** return tetrahedron number if the two facets are adjacent else -1 */
        int identify_tet(Triplet &facet_curr, Triplet &facet_prev);

    void processNode(int nod, Fem2d& fem2d, const Settings& mySettings);
    
    void integrate_u(CGAL::Point &ori, CGAL::Point &ext, std::map<CGAL::Point, Triplet> &iPoint_to_facet, double &path_length, CGAL::Vector &u_integral);

    /** number of ray-facets intersections */
    int search_intersections(const CGAL::Ray &ray_query, std::map<CGAL::Point, Triplet> &iPoint_to_facet);

    void path_integrals(Settings const &mySets, const CGAL::Point &ori, const CGAL::Point &ext, std::map<CGAL::Point, Triplet> &iPoint_to_facet, double &path_length, CGAL::Vector &M_integral, double &contrast);

        private:
            /** ref to the mesh */
            Mesh::mesh & refMsh;

            /** if true the mesh is hollow, filled with non magnetic material */
            const bool filled;

            /** CGAL AABB-Tree */
            CGAL::Tree* tree;  

            /** Node coordinates stores into a vector of CGAL::Point */
            std::vector<CGAL::Point> pt;

            /** tetrahedron mesh splitted into triangles */
            std::list<CGAL::Triangle> triangles;

            /** point to its associted number */
            std::map<CGAL::Point,     int>  point_to_id;                  

            /** quadruplet to its tetrahedron number */
            std::map<Quadruplet, int> quad_to_tetId;

    /**
    Populate a list of CGAL triangles from the tetrahedral mesh. Associate each triplet of points 
    that forms a triangle with the number of the tetrahedron from which it originates.
    */
        void populate_triangles(std::vector<Tetra::Tet> &tet);

        void populate_quads(std::vector<Tetra::Tet> &tet);
        
        std::vector<double> facet_lagrange_polynomials(const Triplet &facet_indices, const CGAL::Point &point_curr);
 
        CGAL::Vector interpolate_u(const Triplet &facet_indices, const CGAL::Point &point);
        }; // end class
    } // end namespace PathInt
#endif
