#ifndef fem2d_h
#define fem2d_h

#include <vector>
#include <initializer_list>
#include <stdexcept>
#include <gmsh.h>

#include <eigen3/Eigen/Dense>

#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "settings.h"

const double VACUUM_PERMEABILITY = 1.25663706144e-6; /* kg m / A^2 s^2   */
const double PLANCKS_H          = 6.62606896e-34;   /* kg m^2 / s */
const double PLANCKS_HBAR       = 1.05457162825e-34; /* kg m^2 / s */
const double CHARGE_ELECTRON    = 1.602176487e-19;   /* C */

enum class ExportType { CONTRAST=0, MZ_INTEGRAL, PATH_LENGTH, HOLO_PHASE };

struct Node2d
    {
	Eigen::Vector2d p;  /**< Physical position p=(x,y)  of the node in the z plane*/
	/** the flag_inside attribute is set to 1 if the node is inside the projected geometry
            of the 3D micromagnetic system in the Oxy plane; otherwise 0 */
	int flag_inside;
        /** path length through the magnetic materials */ 
	double path_length;
	/** magnetization components integrated along beam */
	double Mx_integral, My_integral, Mz_integral;
	/** STXM contrast */
	double contrast;
	double Mx, My; // for holography phase calculations
	/** Magnetic part of the Holography phase */
    double sol;
    };

namespace Triangle {
	static const int NBN = 3, NPI = 3;
	static constexpr double u[NPI] {1/6., 2/3., 1/6.};
	static constexpr double v[NPI] {1/6., 1/6., 2/3.};
	static constexpr double pds[NPI]   {1/6., 1/6., 1/6.};

    /** Interpolation polynomials */
    static constexpr double a[NBN][NPI] = {
        		{1. - u[0] - v[0], 1. - u[1] - v[1], 1. - u[2] - v[2]},
        		{u[0], u[1], u[2]},
        		{v[0], v[1], v[2]}};
	class Tri
		{
		public:
        
		/** Constructor */
		Tri(std::initializer_list<int> _i): surf(0), detJ(0), overlap(false)
			{
			if (_i.size() != NBN)
				{ throw std::invalid_argument("Initializer list must contain exactly 3 elements."); }
			ind.assign(_i.begin(), _i.end());
			//zeroBasing(); // Passage de la convention Matlab/msh à C++
			}

		/** Getter for indices */
		const std::vector<int>& getIndices() const { return ind; }

		int getIndice(int i) const { 
			if (i < 0 || i >= NBN) {
			   throw std::out_of_range("Index out of bounds in getIndice");
		       }
		    return ind[i]; 
		    }

        /** surface of the element */
        double surf;

        /** determinant(J) */
        double detJ;
		    
		/** points de gauss */
		double x[NPI], y[NPI];

		/** poids de gauss x detJ */
		double weight[NPI];

		void chapeaux(std::vector<Node2d> &node);

        /** initialize surf to the surface of the triangle, re-orientate if needed */
        void calc_surf(std::vector<Node2d> &node);

        /** The overlap attribute is set to true if the triangle overlaps the projected geometry 
        of the 3D micromagnetic system in the Oxy plane; otherwise false */
        bool overlap;

		private:
			std::vector<int> ind;

		/** zeroBasing: index convention Matlab/msh (one-based) -> C++ (zero-based) */
			inline void zeroBasing() { for (int& index : ind) { --index; } }
		};
    }

class Fem2d
{
public:
	/** constructor */
	Fem2d(Eigen::Vector3d &_c, Eigen::Vector3d &_l, double _CE, double _V, double _zoomFactor, double _meshSize):
		zoomFactor(_zoomFactor), meshSize(_meshSize), CE(_CE), V(_V)
	    {
	    double xmin = _c[0] - _l[0]/(2.0*zoomFactor);
	    double xmax = _c[0] + _l[0]/(2.0*zoomFactor);
	    double ymin = _c[1] - _l[1]/(2.0*zoomFactor);
	    double ymax = _c[1] + _l[1]/(2.0*zoomFactor);

	    double xymin = std::min(xmin, ymin);
	    double xymax = std::max(xmax, ymax);

		l = Eigen::Vector2d(xymax-xymin, xymax-xymin);
		diam = l.maxCoeff();
		c = Eigen::Vector2d(0.5 * (xymax + xymin), 0.5 * (xymax + xymin));
		grid_generator(xymin, xymax, meshSize);
		chapeaux();
	    }

    inline void infos(void) const
        {
        std::cout << "Sensor : " << std::endl;
        std::cout << "-- sensor mesh -- number of nodes : " << node.size() << std::endl;
        std::cout << "-- meshSize : " << meshSize << std::endl;
	    std::cout << "-- zoom     : " << zoomFactor << std::endl;
        }

    inline int getNbNodes(void) const { return node.size(); }
    inline int getNbTriangles(void) const { return tri.size(); }
    
	int exportRatioGrayScaleImage(const Settings &settings, ExportType eType);
	int exportRatioRGBscaleImage(const Settings &settings, ExportType eType);
	
    int exportMagIntegrals(const std::string &simName);
    int exportHoloPhase(const std::string &simName);
    int exportMagIntegrals_nd_HoloPhase(const std::string &simName);

	/** getter : return node.p */
	inline const Eigen::Vector2d getNode_p(const int i) const { return node[i].p; }

	// Getter for the entire vector of nodes
	inline std::vector<Node2d>& getNodes() { return node; }

	// Getter to get a non constant reference 
	inline Node2d& getNode(const int i) {
		if (i < 0 || i >= static_cast<int>(node.size())) {
			throw std::out_of_range("Index out of bounds in getNode");
		}
		return node[i];
	}

	// Getter for the entire vector of triangles
	inline std::vector<Triangle::Tri>& getTriangles() { return tri; }
	
	// Getter to get a non constant reference 
	inline Triangle::Tri& getTri(const int t) {
		if (t < 0 || t >= static_cast<int>(tri.size())) {
			throw std::out_of_range("Index out of bounds in getTri");
		}
		return tri[t];
	}
	
    /** fix member sol to zero for all nodes in vector node */
    void zero_node_sol(void);

    /** Calculation of the number of active sources */
    int calc_nb_active_sources(void);

    /** initialize many data structures in both Fem2d struct but also triangles */
    void util(void);
    
    /** call chapeaux member function for all triangles */
    void chapeaux(void);

	/** isobarycenter */
	Eigen::Vector2d c;

	/** lengths along x,y,z axis */
	Eigen::Vector2d l;
	
	/** diameter */
	double diam;

    /** a zoom factor relative to the plane where images are computed */
	double zoomFactor;
	
	/** triangular mesh size */
	double meshSize;
	
	double surf;
	
    /** electrostatic constant */
    double CE;

    /** inner potential */
    double V;

private:
	/** 2d nodes */
	std::vector<Node2d> node;

	/** 2d triangles */
	std::vector<Triangle::Tri> tri;

	/**
	 * Generate a structured square mesh in 3D space.
	 * The mesh consists of nodes placed on a regular grid in the XY-plane (Z=0),
	 * and triangles connecting the nodes to form a plane.
	 * Nodes and triangle indices are stored in the corresponding containers.
	 */
	void grid_generator(double xymin, double xymax, double meshSize)
	    {
		int Nx = static_cast<int>((xymax - xymin) / meshSize + 0.5) + 1;
		int Ny = static_cast<int>((xymax - xymin) / meshSize + 0.5) + 1;		
		auto index = [&Nx](int i, int j) -> int { return j * Nx + i; };  // zero based numbering
		
		// Creation de nodes
		for (int ix = 0; ix < Nx; ++ix) {
			for (int iy = 0; iy < Ny; ++iy) {
				Node2d node_{};
				node_.p[0] = xymin + meshSize * ix;
				node_.p[1] = xymin + meshSize * iy;
				node_.p[2] = 0.0;
				node.push_back(node_);
			}
		}

		// Creation of triangles
		for (int ix = 0; ix < Nx - 1; ++ix) {       
			for (int iy = 0; iy < Ny - 1; ++iy) {   
				tri.push_back(Triangle::Tri{ { 
				    index(ix    , iy    ), 
				    index(ix + 1, iy    ), 
				    index(ix + 1, iy + 1) 
				} });

				tri.push_back(Triangle::Tri{ { 
				    index(ix    , iy    ), 
				    index(ix + 1, iy + 1), 
				    index(ix    , iy + 1) 
				} });
			}
		}
		assert (node.size() == Nx*Ny);
		assert (tri.size() == 2*(Nx-1)*(Ny-1));
		}


	/**
	The function handleExport is designed to handle the export of a grayscale image based on a specific export type.
	It takes two parameters:
    - ExportType eType: An enumeration value representing the type of data to export (e.g., MZ_INTEGRAL or CONTRAST).
    - const std::function<double(const Node2d&)>& valueExtractor: A callable function (or lambda) that
      specifies how to extract the desired data from a Node2d object.
	*/
	int handleExport(const Settings &settings, ExportType eType, const std::function<double(const Node2d&)>& valueExtractor);
	int handleRGBexport(const Settings &settings, ExportType eType, const std::function<double(const Node2d&)>& valueExtractor);
};

/** return square */
inline constexpr double sq(const double x) {return x*x;}

namespace pot2D {
    const int scalFMM_height= 8;
    const int scalFMM_order = 7;
    const int NB_PTS_INTEGRATION = 4;// ATTENTION 4 points d'integration != Tri::NPI
    constexpr double ksi = 1.0/sqrt(3.0);
    constexpr double u[] {0.50*(1.0-ksi), 0.50*(1.0-ksi), 0.50*(1.0+ksi), 0.50*(1.0+ksi)};
    constexpr double v[] {0.25*sq(1.0-ksi), 0.25*(1.0-sq(ksi)), 0.25*(1.0-sq(ksi)), 0.25*sq(1.0+ksi) };
    constexpr double w[] {0.125*(1.0-ksi), 0.125*(1.0-ksi), 0.125*(1.0+ksi), 0.125*(1.0+ksi)};
    const int perm[][Triangle::NBN]={{0, 1, 2}, {1, 2, 0}, {2, 0, 1}}; // toutes les permutations circulaires

    /** probably some tree depth that defines the precision */
    const int IPREC = 4;

    /** computation using fmm in 2D */ 
    int fmm2d_sum(Fem2d &fem);     // Greengard's approach
    int scalfmm2d_sum(Fem2d &fem); // Coulaud's   approach
    int direct2d_sum(Fem2d &fem);  // direct   approach
    
    /** correction function for triangle number t */
    void correction(Fem2d &fem,const Triangle::Tri &t,double xk,double yk,double Mxk,double Myk,double wk_detJk);

    /** compute the whole correction */
    void integre_correction(Fem2d &fem,const Triangle::Tri &t);
    }
#endif
