#include <cstdint>
#include <limits>
#include <vector>
#include <filesystem>
#include <iostream>

#include <cmath>
#include <string>
#include <fstream>

#include <stdexcept>
#include <sstream>

#include "fem2d.h"
#include <boost/format.hpp>

int Fem2d::exportMagIntegrals_nd_HoloPhase(const std::string &simName)
    {
	std::string filename = simName+"_STXM_HOLO.out";
	std::ofstream fout(filename);
	if (fout.fail()) exit(1);

	fout << "## columns" << std::endl;
	fout << boost::format("## %3s %20s %20s %2s %20s %20s %20s %20s %20s %20s")
     	% "1" % "2" % "3" % "4" % "5" % "6" % "7" % "8" % "9" % "10" << std::endl;
	fout << boost::format("## %3s %20s %20s %2s %20s %20s %20s %20s %20s %20s")
     	% "idx" % "x" % "y" % "in" % "path_length" % "Mx_integral" % "My_integral" % "Mz_integral" % "contrast" % "phase"<< std::endl;
     
	for (unsigned int nod=0; nod<node.size(); nod++)
	    {
		Node2d& node2d = getNode(nod);
		double x = node2d.p[0];
		double y = node2d.p[1];
		int flag_inside = node2d.flag_inside;
		double path_length = node2d.path_length;
		double Mx_integral = node2d.Mx_integral;
		double My_integral = node2d.My_integral;
		double Mz_integral = node2d.Mz_integral;
		double contrast    = node2d.contrast;
		double phase = CE*V*path_length-CHARGE_ELECTRON/PLANCKS_HBAR*node2d.sol;

		fout << boost::format("%6d %+20.10e %+20.10e %2d %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e")
                    				% nod % x % y % flag_inside % path_length % Mx_integral % My_integral % Mz_integral % contrast % phase << std::endl;
	    }

	fout.close();
	std::cout << "Integrals & Holo Phase saved in " << filename << std::endl;
	return 0;
    }

