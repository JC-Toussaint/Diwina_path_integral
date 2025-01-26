#ifndef fem_h
#define fem_h

/** \file fem.h
\brief principal header, contains the struct fem
This file is called by many source files in feellgood.
It does also contains the definition of many constants for the solver, and for scalfmm
*/
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include <time.h>
#include <unistd.h>

#include "mesh.h"

#include "ANN.h"
#include "settings.h"

#include "geo.h"

#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/math/quaternion.hpp>

using Triplet         = std::array<int,   3>;
using Quadruplet      = std::array<int,   4>;

/** \class Fem
class container to grab altogether all parameters of a simulation, including mesh geometry,
containers for the mesh
*/
class Fem
    {
public:
    /** constructor: call mesh constructor, initialize pts,kdtree and many inner variables */
    inline Fem(Settings const &mySets) : msh(mySets)
        {
        if (mySets.restoreFileName == "")
            {
            msh.init_distrib(mySets);
            }
        else
            {
            msh.readSol(mySets.verbose, mySets.restoreFileName);
            }
        }

    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh */
    Mesh::mesh msh;

    };  // end class

#endif
