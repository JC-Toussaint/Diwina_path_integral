#include <cstdint>
#include <limits>
#include <vector>
#include <filesystem>
#include <iostream>

#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>

#include <stdexcept>
#include <sstream>

#include "fem2d.h"
#include <boost/format.hpp>

int Fem2d::exportMagIntegrals(const std::string &simName)
    {
	std::string filename = simName+"_M_integrals.out";
	std::ofstream fout(filename);
	if (fout.fail()) exit(1);

	fout << "## columns" << std::endl;
	fout << boost::format("## %3s %20s %20s %2s %20s %20s %20s %20s %20s")
     	% "1" % "2" % "3" % "4" % "5" % "6" % "7" % "8" % "9" << std::endl;
	fout << boost::format("## %3s %20s %20s %2s %20s %20s %20s %20s %20s")
     	% "idx" % "x" % "y" % "in" % "path_length" % "Mx_integral" % "My_integral" % "Mz_integral" % "contrast" << std::endl;
     
	// lines are formatted in parallel by chunks, then written in order
	const int NOD = node.size();
	const int CHUNK = 4096;
	const int nbChunks = (NOD + CHUNK - 1) / CHUNK;
	std::vector<std::ostringstream> buffers(nbChunks);
	#pragma omp parallel for schedule(dynamic)
	for (int ch = 0; ch < nbChunks; ch++)
	    {
		std::ostringstream &buf = buffers[ch];
		boost::format line("%6d %+20.10e %+20.10e %2d %+20.10e %+20.10e %+20.10e %+20.10e %+20.10e");
		for (int nod = ch*CHUNK; nod < std::min(NOD, (ch+1)*CHUNK); nod++)
		    {
			const Node2d& node2d = node[nod];
			buf << line % nod % node2d.p[0] % node2d.p[1] % node2d.flag_inside % node2d.path_length
			            % node2d.Mx_integral % node2d.My_integral % node2d.Mz_integral % node2d.contrast << '\n';
		    }
	    }
	for (const std::ostringstream &buf : buffers)
		{ fout << buf.str(); }

	fout.close();
	std::cout << "Integrals saved in " << filename << std::endl;
	return 0;
    }

