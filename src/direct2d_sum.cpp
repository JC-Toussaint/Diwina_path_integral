#include <cmath>
#include <chrono>
#include <execution>
#include <string>

#include "fem2d.h"

// Surcharge de l'opérateur << pour afficher les composantes d'un vecteur
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
	os << "[ ";
	for (const auto& element : vec) {
		os << element << " ";
	}
	os << "]";
	return os;
}

int pot2D::direct2d_sum(Fem2d &fem) {
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double, std::milli> micros;
	std::cout << " direct2d_sum start\n "; 

	start = std::chrono::high_resolution_clock::now();

	fem.zero_node_sol();
	
    std::for_each(fem.getTriangles().begin(),fem.getTriangles().end(),[&fem](const Triangle::Tri &tri)
        {
        if (!tri.overlap) return;//continue;

        for (int k=0; k<Triangle::NPI; k++)
            {
            double xk = tri.x[k];
            double yk = tri.y[k];
            double wk_detJk = tri.weight[k];
			double Mxk = 0;
			double Myk = 0;
			for (int ie=0; ie<Triangle::NBN; ie++){
				int i=tri.getIndice(ie);
				Node2d &node = fem.getNode(i);
				Mxk += Triangle::a[ie][k]*node.Mx;
				Myk += Triangle::a[ie][k]*node.My;
			}

            std::for_each(fem.getNodes().begin(),fem.getNodes().end(),
                          [&xk,&yk,&Mxk,&Myk,&wk_detJk](Node2d &node)
                                {
                                double xi = node.p[0];
                                double yi = node.p[1];
                                double d2 = sq(xk-xi)+sq(yk-yi);
                                node.sol += VACUUM_PERMEABILITY/(2.0*M_PI)*(Mxk*(xi-xk) + Myk*(yi-yk))/d2*wk_detJk;
                                });

            correction(fem, tri, xk, yk, Mxk, Myk, wk_detJk);
            } // endfor k

        integre_correction(fem, tri);
        } );// end for_each on triangles

	end = std::chrono::high_resolution_clock::now();
	micros = end - start;

	std::cout << std::endl << boost::format("%5t fast multipole time %50T. ");
	std::cout << micros.count() << " [ms]\n";
			
	return 0;
}

