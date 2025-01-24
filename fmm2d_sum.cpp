#include<cmath>
#include <chrono>
#include "fem2d.h"

using namespace std;

extern "C"{
void rfmm2d(int *ierr, int *iprec, int *nsource, double *source, double* dipstr, double *dipvec, int *ntarget, double *target, double *pottarg);
}

int pot2D::fmm2d_sum(Fem2d &fem)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double,std::milli> micros;

	const int NDIM = 2;
	const int NOD = fem.getNbNodes();
	const int TRI = fem.getNbTriangles();
	const int NPI = Triangle::NPI;

	std::cout << "DEBUG " << TRI << " " << NPI << std::endl;

	start = std::chrono::high_resolution_clock::now();

	int nsource = fem.calc_nb_active_sources();
	std::vector<double> source(NDIM*nsource);
	std::vector<double> dipstr(nsource);
	std::vector<double> dipvec(NDIM*nsource);

	int ntarget = NOD;
	std::vector<double> target(NDIM*ntarget);
	std::vector<double> pottarg(ntarget); // pottarg[i] is the potential at the ith target

	std::cout << "nb of sources & targets : " << nsource << " " << ntarget << std::endl;

	fem.zero_node_sol();

	for (int i=0; i<NOD; i++)
	{
		target[NDIM*i+0] = fem.getNode(i).p[0];
		target[NDIM*i+1] = fem.getNode(i).p[1];
	}

	int ns=0;
	for (int t=0; t<TRI; t++)
	{
		const Triangle::Tri &tri = fem.getTri(t);
		if (!tri.overlap) continue;

		for (int k=0; k<NPI; k++)
		{
			double xk = tri.x[k];
			double yk = tri.y[k];
			double wk_detJk = tri.weight[k];

			double Mxk = 0;
			double Myk = 0;
			for (int ie=0; ie < Triangle::NBN; ie++){
				int i=tri.getIndice(ie);
				Node2d &node=fem.getNode(i);
				Mxk += Triangle::a[ie][k]*node.Mx;
				Myk += Triangle::a[ie][k]*node.My;
			}

			source[NDIM*ns+0] = xk;
			source[NDIM*ns+1] = yk;

			dipstr[ns] = VACUUM_PERMEABILITY/(2.0*M_PI)*wk_detJk;
			dipvec[NDIM*ns+0] = Mxk;
			dipvec[NDIM*ns+1] = Myk;
			ns++;

			correction(fem, tri, xk, yk, Mxk, Myk, wk_detJk);
		} // endfor k

		integre_correction(fem, tri);
	} // endfor t

	std::cout << "DEBUG ns & nsource : " << ns << " " << nsource << std::endl;
	assert(ns==nsource);

	end = std::chrono::high_resolution_clock::now();
	micros = end-start;
	std::cout << std::endl << boost::format("%5t initialization time %50T. ");
	std::cout << micros.count() << " [ms]\n";
	double cputime = micros.count();

	start = std::chrono::high_resolution_clock::now();

	int ierr;
	int iprec{pot2D::IPREC};
	rfmm2d(&ierr, &iprec, &nsource, source.data(), dipstr.data(), dipvec.data(), &ntarget, target.data(), pottarg.data());

	end = std::chrono::high_resolution_clock::now();
	micros = end-start;

	std::cout << std::endl << boost::format("%5t fast multipole time %50T. ");
	std::cout << micros.count() << " [ms]\n";
	cputime += micros.count();

	start = std::chrono::high_resolution_clock::now();
	for (int i=0; i<NOD; i++){
	    fem.getNode(i).sol -= pottarg[i];
	} // minus sign see potential definition given by Gimbutas-Greengard
	end = std::chrono::high_resolution_clock::now();
	micros = end-start;
	std::cout << std::endl << boost::format("%5t get pot %50T. ");
	std::cout << micros.count() << " [ms]\n";
	cputime += micros.count();

	std::cout << std::endl << boost::format("%5t Greengard's FMM time %50T. ");
	std::cout << cputime << " [ms]\n";

	return ierr;
}

