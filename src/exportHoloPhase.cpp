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

int Fem2d::exportHoloPhase(const std::string &simName)
    {
	std::string filename = simName+"_Holo.out";
	std::ofstream fout(filename);
	if (fout.fail()) exit(1);

	fout << "## columns" << std::endl;
	fout << boost::format("## %3s %20s %20s %2s %20s %20s")
     	% "1" % "2" % "3" % "4" % "5" % "6" << std::endl;
	fout << boost::format("## %3s %20s %20s %2s %20s %20s")
     	% "idx" % "x" % "y" % "in" % "path_length" % "phase" << std::endl;
     	
	for (unsigned int nod=0; nod<node.size(); nod++)
	    {
		Node2d& node = getNode(nod);
		double x = node.p[0];
		double y = node.p[1];
		int flag_inside = node.flag_inside;
		double path_length = node.path_length;
                double phase = CE*V*path_length-CHARGE_ELECTRON/PLANCKS_HBAR*node.sol;

		fout << boost::format("%6d %+20.10e %+20.10e %2d %+20.10e %+20.10e")
                    				% nod % x % y % flag_inside % path_length % phase << std::endl;
	    }

	fout.close();
	std::cout << "Holography phase saved in " << filename << std::endl;
	return 0;
    }

