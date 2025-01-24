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
			zeroBasing(); // Passage de la convention Matlab/msh à C++
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
	const int DIM_OBJ_1D = 1;
	const int DIM_OBJ_2D = 2;
	const int DIM_OBJ_3D = 3;
	const int SIZE_TRIANGLE = 3;

	/** constructor */
	//Fem2d(double xymin, double xymax, double meshSize_):meshSize(meshSize_)
	Fem2d(Eigen::Vector3d &_c, Eigen::Vector3d &_l, double _zoomFactor, double meshSize_):zoomFactor(_zoomFactor), meshSize(meshSize_)
	    {
	    double xmin = _c[0] - _l[0]/(2.0*zoomFactor);
	    double xmax = _c[0] + _l[0]/(2.0*zoomFactor);
	    double ymin = _c[1] - _l[1]/(2.0*zoomFactor);
	    double ymax = _c[1] + _l[1]/(2.0*zoomFactor);
	    double xymin = std::min(xmin, ymin);
	    double xymax = std::min(xmax, ymax);
		l = Eigen::Vector2d(xymax-xymin, xymax-xymin);
		c = Eigen::Vector2d(0.5 * (xymax + xymin), 0.5 * (xymax + xymin));
		grid_generator(xymin, xymax, meshSize);
		chapeaux();
		
		CE=0; /* A REVOIR */
		V=0;
	    }

    inline void infos(void) const
        {
        std::cout << "Detector : " << std::endl;
        std::cout << "-- detector mesh -- number of nodes : " << node.size() << std::endl;
        std::cout << "-- meshSize : " << meshSize << std::endl;
	    std::cout << "-- zoom     : " << zoomFactor << std::endl;
        }

    inline int getNbNodes(void) const { return node.size(); }
    inline int getNbTriangles(void) const { return tri.size(); }
    inline int getNbTrianglesWithSources(void) const { return triWithSources.size(); }
    
	int exportRatioGrayScaleImage(const Settings &settings, ExportType eType);
	int exportRatioRGBscaleImage(const Settings &settings, ExportType eType);
	
    int exportMagIntegrals(const std::string &simName);
    int exportHoloPhase(const std::string &simName);

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

	// Getter to get a non constant reference 
	inline Triangle::Tri& getTri(const int t) {
		if (t < 0 || t >= static_cast<int>(tri.size())) {
			throw std::out_of_range("Index out of bounds in getTri");
		}
		return tri[t];
	}

	// Getter for the entire vector of triangles
	inline std::vector<Triangle::Tri>& getTriangles() { return tri; }

	// Getter to get a non constant reference to the triangle t with sources
	inline Triangle::Tri& getTriWithSources(const int t) {
		if (t < 0 || t >= static_cast<int>(triWithSources.size())) {
			throw std::out_of_range("Index out of bounds in getTri");
		}
		return triWithSources[t];
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

    /** a zoom factor relative to the plane where images are computed */
	double zoomFactor;
	
	/** triangular mesh size */
	double meshSize;
	
	double surf;
	
    /** electric constant */
    double CE;

    /** applied potential */
    double V;

private:
	/** 2d nodes */
	std::vector<Node2d> node;

	/** 2d triangles */
	std::vector<Triangle::Tri> tri;

	/** 2d triangles in which are located sources */
	std::vector<Triangle::Tri> triWithSources;

	/**
	using gmsh geo, rectangle method buids a triangular mesh of a rectangle ((xmin,ymin),(xmax,ymax)) in 3D space,
	in a xOy plane at altitude z, with meshSize carac length
	if outputFile is true a file named rectangle.msh is written (text format 2.2)
	if structured is true output mesh file is structured, meaning that nodes are on a regular grid
	 */
	int rectangle(bool outputFile,bool structured,double xmin,double xmax,double ymin,double ymax,double z,double meshSize)
	    {
		using namespace gmsh::model;
		const int idx_rectangle = 1;

		geo::addPoint(xmin,ymin,z,meshSize,1);
		geo::addPoint(xmax,ymin,z,meshSize,2);
		geo::addPoint(xmax,ymax,z,meshSize,3);
		geo::addPoint(xmin,ymax,z,meshSize,4);
        geo::addLine(1,2,1);
		geo::addLine(2,3,2);
		geo::addLine(3,4,3);
		geo::addLine(4,1,4);
        geo::addCurveLoop({1,2,3,4},1);
		geo::addPlaneSurface({1},idx_rectangle);

		if(structured)
		    { geo::mesh::setTransfiniteSurface(idx_rectangle,"Left",{1,2,3,4}); }

		geo::synchronize();
        addPhysicalGroup(DIM_OBJ_1D,{1,2,3,4},5);
		setPhysicalName(DIM_OBJ_1D,5,"frontier");
		addPhysicalGroup(DIM_OBJ_2D,{1},10);
		setPhysicalName(DIM_OBJ_2D,10,"Rectangle");
        mesh::generate(DIM_OBJ_3D);

		if(outputFile)
		    {
			gmsh::option::setNumber("Mesh.MshFileVersion",2.2);
			gmsh::write("detector.msh");
		    }
		return idx_rectangle;
	    }

	/**
	build a mesh of a rectangle in 3D space(calling above rectangle function),
	 and print back to terminal nodes and triangle indices of the rectangle mesh
	 */
	void grid_generator(double xymin, double xymax, double meshSize)
	    {
		const bool OUTPUT_FILE = true;
    	const bool STRUCTURED = true;

		gmsh::initialize();
		gmsh::option::setNumber("General.ExpertMode", 1);
		gmsh::model::add("Rectangle");

		int idx_mesh = rectangle(OUTPUT_FILE,STRUCTURED,xymin,xymax,xymin,xymax,0.0,meshSize);
		std::vector<std::size_t> nodeT;
		std::vector<double> coord;
		gmsh::model::mesh::getNodesForPhysicalGroup(DIM_OBJ_2D,10,nodeT,coord);

		std::for_each(nodeT.begin(),nodeT.end(),[this, &coord](std::size_t &idx)
			{
			long j=idx-1;
			Node2d node_;
			node_.p[0] = coord[3*j+0];
			node_.p[1] = coord[3*j+1];
			node_.p[2] = coord[3*j+2];
			node.push_back(node_);
			});

		std::vector<int> elemTypes;
		std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
		gmsh::model::mesh::getElements(elemTypes,elemTags,elemNodeTags,DIM_OBJ_2D,idx_mesh);

		for(unsigned int i=0;i<elemNodeTags[0].size();i+=SIZE_TRIANGLE)
		    {
			int i0 = elemNodeTags[0][i+0];
			int i1 = elemNodeTags[0][i+1];
			int i2 = elemNodeTags[0][i+2];
			Triangle::Tri tri_({i0, i1, i2});
			tri.push_back(tri_);
		    }
		gmsh::finalize();
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
    const int NB_PTS_INTEGRATION = 4;// ATTENTION 4 points d'integration != Tri::NPI
    constexpr double ksi = 1.0/sqrt(3.0);
    constexpr double u[] {0.50*(1.0-ksi), 0.50*(1.0-ksi), 0.50*(1.0+ksi), 0.50*(1.0+ksi)};
    constexpr double v[] {0.25*sq(1.0-ksi), 0.25*(1.0-sq(ksi)), 0.25*(1.0-sq(ksi)), 0.25*sq(1.0+ksi) };
    constexpr double w[] {0.125*(1.0-ksi), 0.125*(1.0-ksi), 0.125*(1.0+ksi), 0.125*(1.0+ksi)};
    const int perm[][Triangle::NBN]={{0, 1, 2}, {1, 2, 0}, {2, 0, 1}}; // toutes les permutations circulaires

    /** probably some tree depth that defines the precision */
    const int IPREC = 4;

    /** computation using fmm in 2D */ 
    int fmm2d_sum(Fem2d &fem);

    /** computation using direct in 2D */
    int direct2d_sum(Fem2d &fem);

    /** correction function for triangle number t */
    void correction(Fem2d &fem,const Triangle::Tri &t,double xk,double yk,double Mxk,double Myk,double wk_detJk);

    /** compute the whole correction */
    void integre_correction(Fem2d &fem,const Triangle::Tri &t);
    }
#endif
