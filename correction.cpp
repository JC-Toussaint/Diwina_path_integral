#include "fem2d.h"

void pot2D::correction(Fem2d &fem,const Triangle::Tri &t,double xk,double yk,double Mxk,double Myk,double wk_detJk)
    {
    for (int ie=0; ie<Triangle::NBN; ie++)
	    {
	    int i = t.getIndice(ie);
        Node2d &node = fem.getNode(i);
        double xi = node.p[0];
        double yi = node.p[1];
        double d2 = sq(xk-xi)+sq(yk-yi);
        node.sol -= VACUUM_PERMEABILITY/(2.0*M_PI)*(Mxk*(xi-xk) + Myk*(yi-yk))/d2*wk_detJk;
        }//endfor ie
    }

/*
     v
     ^
   1 |       .
     |     ./|
     |   ./  |
     | ./    |
     |/      |
     |---------> u   ATTENTION  triangle toujours pointu en 0 => permuter circulairement la table de connectivite
     0       1
*/

void pot2D::integre_correction(Fem2d &fem,const Triangle::Tri &tri)
    {
    for (int k=0; k<NB_PTS_INTEGRATION; k++)
        {
        for (int iperm=0; iperm<3; iperm++)
            { // iperm permutation
	        int ie0 = perm[iperm][0], ie1 =  perm[iperm][1], ie2 =  perm[iperm][2];
            int i0 = tri.getIndice(ie0), i1 = tri.getIndice(ie1), i2 = tri.getIndice(ie2);
	        double x0 = fem.getNode(i0).p[0], x1 = fem.getNode(i1).p[0], x2 = fem.getNode(i2).p[0];
            double y0 = fem.getNode(i0).p[1], y1 = fem.getNode(i1).p[1], y2 = fem.getNode(i2).p[1];
            double alpha0 = 1.0-u[k], alpha1 = u[k]-v[k], alpha2 = v[k];
            double xk = alpha0*x0 + alpha1*x1 + alpha2*x2;
	        double yk = alpha0*y0 + alpha1*y1 + alpha2*y2;
	        double wk_detJk = w[k]*tri.detJ;
	        
            double Mxk = 0;
            double Myk = 0;
            for (int ie=0; ie<Triangle::NBN; ie++){
                int i=tri.getIndice(ie);
                Node2d &node=fem.getNode(i);
                Mxk += Triangle::a[ie][k]*node.Mx;
                Myk += Triangle::a[ie][k]*node.My;          
                }

	        Node2d &node = fem.getNode(i0);
	        double d2 = sq(xk-x0)+sq(yk-y0);
	        node.sol += VACUUM_PERMEABILITY/(2.0*M_PI)*(Mxk*(x0-xk) + Myk*(y0-yk))/d2*wk_detJk;
	        }
        } // endfor k
    }
