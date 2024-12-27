#include "fem2d.h"

using namespace std;

void Triangle::Tri::chapeaux(std::vector<Node2d> &node)
{
	vector<double> xp(NBN), yp(NBN);    
	for (size_t ie=0; ie<NBN; ie++)
	{
		size_t i = ind[ie];
		xp[ie] = node[i].p[0];
		yp[ie] = node[i].p[1];
	}

    double detJ = (xp[1]-xp[0])*(yp[2]-yp[0])-(xp[2]-xp[0])*(yp[1]-yp[0]);
    assert(detJ > 0);
    
	for (size_t npi=0; npi<NPI; npi++)
	{
	    vector<double> _a{1. - u[npi] - v[npi], u[npi], v[npi] };                
        /* Gauss point locations */
        x[npi] = _a[0]*xp[0]+_a[1]*xp[1]+_a[2]*xp[2];
        y[npi] = _a[0]*yp[0]+_a[1]*yp[1]+_a[2]*yp[2];
        /* Gauss weights multiplied by detJ */
		weight[npi] = detJ * pds[npi];
	}
}

void Fem2d::chapeaux(void)
{
	/****************** TRIANGLES *****************/
	std::for_each(tri.begin(),tri.end(),
			[this](Triangle::Tri &t) { t.chapeaux(node); } );
}
