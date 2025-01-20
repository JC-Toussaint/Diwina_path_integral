#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, readMesh is expecting a mesh file in gmsh text format 2.2, with first order tetraedrons and triangular facettes.
 */

#include <algorithm>
#include <execution>
#include <map>
#include <numeric>
#include <set>

#include "facette.h"
#include "node.h"
#include "tetra.h"
#include "geo.h" // JCT
#include "settings.h"

namespace Mesh
{


/** \class mesh
class for storing the mesh, including mesh geometry values, containers for the nodes, triangular
faces and tetrahedrons. nodes data are private. They are accessible only through getter and setter.
 */
class mesh
{
public:
	/** constructor : read mesh file, reorder indices and computes some values related to the mesh :
     center and length along coordinates and diameter = max(l(x|y|z)), volume and surface */
	inline mesh(Settings const &mySets /**< [in] */)
	{
		readMesh(mySets);
		indexReorder(mySets.paramTetra);  // reordering of index nodes for facette orientation
		if (mySets.verbose)
		{ std::cout << "  reindexed\n"; }

		double xmin = minNodes(Nodes::IDX_X);
		double xmax = maxNodes(Nodes::IDX_X);

		double ymin = minNodes(Nodes::IDX_Y);
		double ymax = maxNodes(Nodes::IDX_Y);

		double zmin = minNodes(Nodes::IDX_Z);
		double zmax = maxNodes(Nodes::IDX_Z);

		l = Eigen::Vector3d(xmax - xmin, ymax - ymin, zmax - zmin);
		diam = l.maxCoeff();
		c = Eigen::Vector3d(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));

		vol = std::transform_reduce(EXEC_POL, tet.begin(), tet.end(), 0.0, std::plus{},
				[](Tetra::Tet const &te) { return te.calc_vol(); });

		surf = std::transform_reduce(EXEC_POL, fac.begin(), fac.end(), 0.0, std::plus{},
				[](Facette::Fac const &fa) { return fa.surf; });
	}

	//    ~mesh() { delete tree; }

	/** return lengths, used after rotation */
	void update_lengths(){
		double xmin = minNodes(Nodes::IDX_X);
		double xmax = maxNodes(Nodes::IDX_X);

		double ymin = minNodes(Nodes::IDX_Y);
		double ymax = maxNodes(Nodes::IDX_Y);

		double zmin = minNodes(Nodes::IDX_Z);
		double zmax = maxNodes(Nodes::IDX_Z);

		l = Eigen::Vector3d(xmax - xmin, ymax - ymin, zmax - zmin);
		diam = l.maxCoeff();
		c = Eigen::Vector3d(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));
	}

	/** return number of nodes  */
	inline int getNbNodes(void) const { return node.size(); }

	/** return number of triangular fac */
	inline int getNbFacs(void) const { return fac.size(); }

	/** return number of tetrahedrons */
	inline int getNbTets(void) const { return tet.size(); }

	/** getter : return node.p */
	inline const Eigen::Vector3d getNode_p(const int i) const { return node[i].p; }

	/** getter : return node.u */
	inline const Eigen::Vector3d getNode_u(const int i) const { return node[i].get_u(Nodes::NEXT); }

	/** setter for u0 */
	inline void set_node_u0(const int i, Eigen::Vector3d const &val)
	{
		node[i].d[Nodes::CURRENT].u = val;
	}

	/** setter for p */ // JCT
	inline void setNode_p(const int i, Eigen::Vector3d const &val)
	{
		node[i].p = val;
	}

	/** setter for u */
	inline void setNode_u(const int i, Eigen::Vector3d const &val)
	{
		node[i].d[Nodes::NEXT].u = val;
	}

	/** basic informations on the mesh */
	void infos(void) const
	{
		std::cout << "mesh:\n";
		std::cout << "  bounding box diam:  " << diam << '\n';
		std::cout << "  nodes:              " << getNbNodes() << '\n';
		std::cout << "  faces:              " << fac.size() << '\n';
		std::cout << "  tetraedrons:        " << tet.size() << '\n';
		std::cout << "  total surface:      " << surf << '\n';
		std::cout << "  total volume:       " << vol << '\n';
	}

	/** isobarycenter */
	Eigen::Vector3d c;

	/** lengths along x,y,z axis */
	Eigen::Vector3d l;

	/** max of l coordinates, to define a bounding box */
	double diam;

	/** total surface */
	double surf;

	/** total volume of the mesh */
	double vol;

	/** face container */
	std::vector<Facette::Fac> fac;

	/** tetrahedron container */
	std::vector<Tetra::Tet> tet;

	/** read a solution from a file (tsv formated) and initialize fem struct to restart computation
	 * from that distribution, return time
	 */
	double readSol(bool VERBOSE /**< [in] */,
			const std::string fileName /**< [in] input .sol text file */);

	/** computes an analytical initial magnetization distribution as a starting point for the
	 * simulation */
	inline void init_distrib(Settings const &mySets /**< [in] */)
	{
		std::for_each(node.begin(),node.end(),[&mySets](Nodes::Node &n)
				{
			n.d[Nodes::CURRENT].u = mySets.getMagnetization(n.p);
			n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
				} );
	}

	/** setter for node[i]; what_to_set will fix what is the part of the node struct to set (usefull
	 * for fmm_demag.h) */
	inline void set(const int i, std::function<void(Nodes::Node &, const double)> what_to_set,
			const double val)
	{
		what_to_set(node[i], val);
	}

private:
	/** node container: not initialized by constructor, but later while reading the mesh by member
	 * function init_node */
	std::vector<Nodes::Node> node;

	/** map of the surface region physical names from mesh file */
	std::map<int, std::string> surfRegNames;

	/** map of the volume region physical names from mesh file */
	std::map<int, std::string> volRegNames;

	/** test if mesh file contains surfaces and regions mentionned in yaml settings and their dimensions */
	void checkMeshFile(Settings const &mySets /**< [in] */);

	/** read Nodes from mesh file */
	void readNodes(Settings const &mySets /**< [in] */);

	/** read tetraedrons of the settings volume regions */
	void readTetraedrons(Settings const &mySets /**< [in] */);

	/** read facettes of the settings surface regions */
	void readTriangles(Settings const &mySets /**< [in] */);

	/** reading mesh format 2.2 text file function */
	void readMesh(Settings const &mySets /**< [in] */);

	/** memory allocation for the nodes */
	inline void init_node(const int Nb) { node.resize(Nb); }

	/** loop on nodes to apply predicate 'whatTodo'  */
	double doOnNodes(const double init_val, const Nodes::index coord,
			std::function<bool(double, double)> whatToDo) const
	{
		double result(init_val);
		std::for_each(node.begin(), node.end(),
				[&result, coord, whatToDo](Nodes::Node const &n)
				{
			double val = n.p(coord);
			if (whatToDo(val, result)) result = val;
				});
		return result;
	}

	/** return the minimum of all nodes coordinate along coord axis */
	inline double minNodes(const Nodes::index coord) const
	{
		return doOnNodes(__DBL_MAX__, coord, [](double a, double b) { return a < b; });
	}

	/** return the maximum of all nodes coordinate along coord axis */
	inline double maxNodes(const Nodes::index coord) const
	{
		return doOnNodes(__DBL_MIN__, coord, [](double a, double b) { return a > b; });
	}

	/** redefine orientation of triangular faces in accordance with the tetrahedron
	 * reorientation of the tetrahedrons if needed; definition of Ms on facette elements
    Indices and orientation convention :

                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w

	 */
	void indexReorder(std::vector<Tetra::prm> const &prmTetra)
	{
		std::set<Facette::Fac> sf;  // implicit use of operator< redefined in class Fac

		std::for_each(tet.begin(), tet.end(), [this,&sf](Tetra::Tet const &te)
				{
			const int ia = te.ind[0];
			const int ib = te.ind[1];
			const int ic = te.ind[2];
			const int id = te.ind[3];

			sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ic, ib} ));
			sf.insert(Facette::Fac(node, 0, te.idxPrm, {ib, ic, id} ));
			sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, id, ic} ));
			sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ib, id} ));
				});  // end for_each

				std::for_each(fac.begin(), fac.end(), [this, &prmTetra, &sf](Facette::Fac &fa)
						{
					int i0 = fa.ind[0], i1 = fa.ind[1], i2 = fa.ind[2];
					std::set<Facette::Fac>::iterator it = sf.end();
					for (int perm = 0; perm < 2; perm++)
					{
						for (int nrot = 0; nrot < 3; nrot++)
						{
							Facette::Fac fc(node, 0, 0, {0, 0, 0} );
							fc.ind[(0 + nrot) % 3] = i0;
							fc.ind[(1 + nrot) % 3] = i1;
							fc.ind[(2 + nrot) % 3] = i2;
							it = sf.find(fc);
							if (it != sf.end()) break;
						}
						std::swap(i1, i2);  // it seems from ref archive we do not want to swap
						// inner fac indices but local i1 and i2
						fa.n = fa.calc_norm(); // update normal vector
					}                   // end perm
						});  // end for_each
	}
};

}  // end namespace Mesh

#endif
