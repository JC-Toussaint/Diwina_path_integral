#include<cmath>
#include<limits>

#include "fem2d.h"

using namespace std;

void Fem2d::zero_node_sol(void)
{
	std::for_each(node.begin(),node.end(),
			[](Node2d &n){ n.sol = 0.0; } );
}

int Fem2d::calc_nb_active_sources(void)
{
	int nsource(0);
	std::for_each(tri.begin(),tri.end(),[this,&nsource](Triangle::Tri &tri)
			{
		if (tri.overlap){
			nsource += Triangle::NPI;
		}
			});
	//    std::cout << "DEBUG femutil nsource : " << nsource << std::endl;
	return nsource;
}

void Triangle::Tri::calc_surf(std::vector<Node2d> &node)
{
	int i0 = ind[0];
	int i1 = ind[1];
	int i2 = ind[2];
	double x0,y0, x1,y1, x2,y2;
	x0 = node[i0].p[0];   y0 = node[i0].p[1];
	x1 = node[i1].p[0];   y1 = node[i1].p[1];
	x2 = node[i2].p[0];   y2 = node[i2].p[1];

	surf = 0.5*((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0));
	if (surf < 0.)
	{
		ind[2]=i1;
		ind[1]=i2;
		surf *= -1.;
	}
}

void Fem2d::util(void)
{
	surf = 0.0;
	// calcul des surfaces et reorientation des segments si necessaire
	std::for_each(tri.begin(),tri.end(),
			[this](Triangle::Tri &t) { t.calc_surf(node); surf += t.surf; } );

	/*

 Throwing a ray vertically from a star node intercepts the 3D micromagnetic mesh, while the others do not.
 The “overlap” variable of a triangle with one or more star nodes is set to true otherwise false.

	 * inside
   o outside         triangles whose overlap variable is set to true

	 *---*---o       *---*   o
     |\  |\  |       |\  |\
     | \ | \ |       | \ | \
     |  \|  \|       |  \|  \
     o---*---o       *---*---*
     |\  |\  |        \  |\  |
     | \ | \ |         \ | \ |
     |  \|  \|          \|  \|
     o---o---o       o   *---*

	 */
	std::for_each(tri.begin(),tri.end(),[this](Triangle::Tri &tri)
			{
		bool inside=false;
		tri.overlap = false;
		for (int ie=0; ie < Triangle::NBN; ie++){
			int i=tri.getIndice(ie);
			Node2d &node_ = node[i];
			inside |= node_.flag_inside;
		}
		if (inside){
			tri.overlap = true;
		}
			});
}
