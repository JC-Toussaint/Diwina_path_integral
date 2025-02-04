# pathIntegral (from Diwina project)
# based on FeeLLGood – A micromagnetic solver ![Build Status](https://github.com/feellgood/FeeLLGood/actions/workflows/tests.yml/badge.svg)

pathIntegral is a STXM image simulator using feeLLGood simulations inputs. feeLLGood is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation, developped by JC Toussaint & al. The code is being modified without any warranty it works. A dedicated website can be found [here][]  

### Dependencies

* C++17 and the STL
* Fortran 2018 (gfortran)
* [GMSH][]
* [TBB][]
* [libpng][]
* [yaml-cpp][]
* [Duktape][] 2.7.0
* [Eigen][] ≥ 3.3
* BLAS
* LAPACK
* [fmmlib2d][] 1.2.4

### License

Copyright (C) 2024  Jean-Christophe Toussaint, C. Thirion and E. Bonet.

pathIntegral is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Additional permission under GNU GPL version 3 section 7: If you modify this Program, or any covered work, by linking or combining it with the Intel® MKL library (or a modified version of that library), containing parts covered by the terms of Intel Simplified Software License, the licensors of this Program grant you additional permission to convey the resulting work.

The libraries used by pathIntegral are distributed under different licenses, and this is documented in their respective Web sites.

[here]: https://feellgood.neel.cnrs.fr/
[GMSH]: https://gmsh.info/
[TBB]: https://www.threadingbuildingblocks.org/
[libpng]: http://libpng.org/pub/png/libpng.html
[yaml-cpp]: https://github.com/jbeder/yaml-cpp
[Duktape]: https://duktape.org/
[Eigen]: https://eigen.tuxfamily.org/
[fmmlib2d]: https://github.com/zgimbutas/fmmlib2d

A configuration template file can be generated, giving you all the fields you need to fill in to run a simulation.  To do this, enter :
./pathIntegral --print-defaults > settings.yml

The resulting file is in YAML format. It has lots of comments (everything after a # sign) documenting the file structure 
and the meaning of the individual settings.
With the comments stripped, the file looks like this:

outputs:
  directory: .
  file_basename: sim
mesh:
  filename: sphere.msh
  length_unit: 1e-9
  volume_regions:
    sphere_volume:
      Js: 1
      lbp: 0.01
      lbm: 0.018
  surface_regions:
    sphere_surface:
initial_magnetization:
  - 1
  - 0
  - 0
rotations:
  angle1: 0
  axe1: [1, 0, 0]
  angle2: 0
  axe2: [0, 0, 1]
filled: false
electrostatics:
  CE: 0
  V : 0
detector:
  zoom: 0.5
  meshSize: 1.0e-9

The previous settings file corresponds to the case of a sphere uniformly magnetized along the Ox axis, with a radius of r=100nm. 
The Gmsh mesh generator has already been used to generate the mesh of the sphere. We define two regions: sphere_volume and sphere_surface. 
For the sphere_volume region, you need to specify the material’s magnetization in Tesla, as well as the photon absorption lengths 
for left- and right-circular polarizations in nanometers.
The micromagnetic system can be rotated along two axes (i.e., axis1 and axis2) by specifying the rotation angles (i.e., angle1 and angle2) in degrees.
The filled flag is not relevant for a sphere entirely made of magnetic material. However, it is relevant for tubes filled with a non-magnetic material, such as copper (A VERIFIER).
The CE and V parameters account for electrostatic effects on the holographic phase. CE depends on the type of microscope and the electron energy, while the inner potential V depends on the type of material.
Finally, the holographic phase image is calculated at any point on a grid called the detector. Since its dimensions must be larger than those of the projected system, the zoom parameter must be less than 1 in the current version. The resolution is determined by the meshSize parameter.

 The program generates the following ASCII files:

    sim_M_integrals.out – Columns represent the node number, its x and y coordinates, a boolean value indicating whether the beam propagates through materials, the length of material traversed, the integrals of Mx, My, and Mz over the beam path, and the STXM intensity contrast.
    sim_Holo.out – Columns represent the node number, its x and y coordinates, a boolean value indicating whether the beam propagates through materials, the length of material traversed, and the holographic phase in radians, taking into account both electrostatic and magnetic contributions.

The program also generates the following image files:

    sim_MZ.png – A 16-bit grayscale map of the integral of Mz.
    sim_PATH_LENGTH.png – A 16-bit grayscale map of the material path length.
    sim_STXM_XMCD.png – A 16-bit grayscale map of the STXM contrast.
    sim_HOLO_PHASE.png – An 8-bit indexed RGB map of the holographic phase.
 
 
