#include <cmath>
#include <chrono>
#include <execution>
#include <string>

#include "fem2d.h"
#include "scalfmm3.hpp"
#include "scalfmm_2d.hpp"
#include "scalfmm/utils/sort.hpp"
#include "scalfmm/tools/fma_loader.hpp"

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

int pot2D::scalfmm2d_sum(Fem2d &fem) {
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double, std::milli> micros;
	std::cout << " scalfmm2d_sum start\n "; const int NDIM = 2;
	const int NOD = fem.getNbNodes();
	const int TRI = fem.getNbTriangles();
	const int NBN = Triangle::NBN;
	const int NPI = Triangle::NPI;

	start = std::chrono::high_resolution_clock::now();

	int nsource = fem.calc_nb_active_sources();
	std::vector<particle_source_type> source(nsource);
	// poids
	std::vector<double> dipstr(nsource);

	int ntarget = fem.getNbNodes();
	std::vector<particle_target_type> target(ntarget);
	std::vector<double> pottarg( ntarget); // pottarg[i] is the potential at the ith target

	std::cout << "nb of sources & targets : " << nsource << " " << ntarget << std::endl;

	fem.zero_node_sol();
	point_type S_max(-100000.0), S_min(100000.0);

	std::cout << "Calcul de la boite englobante \n";
	// On remplit la structure cible = noeuds du maillage
	//  positon en espace + indice initiale (avant renumerotation)
	for (int i = 0; i < ntarget; ++i) {
		auto &pos = target[i].position();
		pos[0] = fem.getNode(i).p[0];
		pos[1] = fem.getNode(i).p[1];
		target[i].variables(i);
		S_max[0] = std::max(S_max[0], pos[0]);
		S_max[1] = std::max(S_max[1], pos[1]);
		S_min[0] = std::min(S_min[0], pos[0]);
		S_min[1] = std::min(S_min[1], pos[1]);
	}

	// Construct the Gauss points (sources)
	int ns = 0;
	// Boucle sur les triangles
	for (int t = 0; t < TRI; t++) {
		Triangle::Tri &tri = fem.getTri(t);
		if (!tri.overlap) continue;

		// Boucle sur les points de Gauss du triangle
		for (int k = 0; k < NPI; k++) {
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
			dipstr[ns] = VACUUM_PERMEABILITY/(2.0*M_PI) * wk_detJk;

			// On met les informations dans le format scalfmm
			auto &pos = source[ns].position();
			pos[0] = xk;
			pos[1] = yk;
			S_max[0] = std::max(S_max[0], xk);
			S_max[1] = std::max(S_max[1], yk);
			S_min[0] = std::min(S_min[0], xk);
			S_min[1] = std::min(S_min[1], yk);
			// garde la position initale de la source
			source[ns].variables(ns);
			//
			auto &inputs = source[ns].inputs();
			inputs[0] = Mxk * dipstr[ns];
			inputs[1] = Myk * dipstr[ns];

			++ns;
		    correction(fem, tri, xk, yk, Mxk, Myk, wk_detJk);
		} // endfor k

		integre_correction(fem, tri);
	} // endfor t
	
	S_min -= 0.1 ;
	S_max += 0.1;

	assert(ns == nsource);
	// Nous trions les particules en fonction de leur indice de Morton  !.
	int level_max{ pot2D::scalFMM_height + 1} ;
	box_type box(S_min, S_max);
	std::cout << "  box_type box(S_min, S_max)  " << box<< std::endl;

	//   box
	scalfmm::utils::sort_container(box, level_max, source);
	scalfmm::utils::sort_container(box, level_max, target);

	const int group_size = 10;
	bool sorted = true;
	// the container is now sorted
	int height = pot2D::scalFMM_height;
	int order  = pot2D::scalFMM_order;
	tree_source_type tree_source(
			height, order, box, group_size, group_size, source, sorted);

	tree_target_type tree_target(height, order, box, group_size, group_size,
			target, sorted);
	end = std::chrono::high_resolution_clock::now();
	micros = end - start;
	std::cout << std::endl << boost::format("%5t initialization time %50T. ");
	std::cout << micros.count() << " [ms]\n";
	start = std::chrono::high_resolution_clock::now();
	tree_source.statistics("tree source" ,std::cout) ;
	tree_target.statistics("tree trarget" ,std::cout) ;
	int ierr{0};
	int iprec(IPREC);
	scalfmm_execute(ierr, iprec, tree_source, tree_target, pottarg);

	end = std::chrono::high_resolution_clock::now();
	micros = end - start;

	std::cout << std::endl << boost::format("%5t fast multipole time %50T. ");
	std::cout << micros.count() << " [ms]\n";

	for (int i = 0; i < NOD; i++) {
	    fem.getNode(i).sol -= pottarg[i];
	} // minus sign see potential definition given by Gimbutas-Greengard

	return ierr;
}

