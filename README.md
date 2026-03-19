# pathIntegral from Diwina project
## based on FeeLLGood – A micromagnetic solver

pathIntegral is a STXM-XMCD and electron holography image simulator using feeLLGood simulations inputs. feeLLGood is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation, developped by JC Toussaint & al. The code is being modified without any warranty it works. A dedicated website can be found [here][]  

### Dependencies

* C++17 and the STL
* [GMSH][]
* [TBB][]
* [libpng][]
* [yaml-cpp][]
* [Duktape][] 2.7.0
* [Eigen][] 3.4
* BLAS
* LAPACK

### License

Copyright (C) 2024  Jean-Christophe Toussaint, C. Thirion and E. Bonet.

pathIntegral is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Additional permission under GNU GPL version 3 section 7: If you modify this program, or any covered work, by linking or combining it with the Intel® MKL library (or a modified version of that library), containing parts covered by the terms of Intel Simplified Software License, the licensors of this Program grant you additional permission to convey the resulting work.

The libraries used by pathIntegral are distributed under different licenses, and this is documented in their respective Web sites.

[here]: https://feellgood.neel.cnrs.fr/
[GMSH]: https://gmsh.info/
[TBB]: https://www.threadingbuildingblocks.org/
[libpng]: http://libpng.org/pub/png/libpng.html
[yaml-cpp]: https://github.com/jbeder/yaml-cpp
[Duktape]: https://duktape.org/
[Eigen]: https://eigen.tuxfamily.org/
[fmmlib2d]: https://github.com/zgimbutas/fmmlib2d
 
# Detailed Installation Protocol

## For Windows only, installation of WSL2 (Windows Subsystem for Linux)

** Install Windows Subsystem for Linux (WSL) - This enables running a Linux environment directly on Windows **
```bash
wsl --install
```

## For Linux or WSL2, initial Setup 

** Update the list of available packages - This refreshes the package database to ensure you have access to the latest versions of software available in the repositories **
```bash
sudo apt update
```

** Upgrade all installed packages to their latest versions - This ensures system security and compatibility by installing the most recent updates and patches **
```bash
sudo apt upgrade
```

## System Dependencies Installation


```bash
sudo apt install wget git cmake pkgconf g++ gfortran python3-venv
sudo apt install libgmsh-dev libeigen3-dev libtbb-dev libyaml-cpp-dev duktape-dev libfftw3-dev libopenblas-dev libpng-dev libcgal-dev
```
** ParaView (scientific visualization software) - A powerful, open-source application developed by Kitware, Sandia National Laboratories, and Los Alamos National Laboratory for visualizing large datasets, commonly used in scientific and engineering applications for 3D data analysis **

** FFTW library for fast Fourier transforms - A high-performance C library developed by Matteo Frigo and Steven G. Johnson at MIT for computing discrete Fourier transforms, essential for signal processing and numerical analysis **

** OpenBLAS (optimized linear algebra library) - A highly optimized implementation of Basic Linear Algebra Subprograms (BLAS) developed by Zhang Xianyi at the Chinese Academy of Sciences, providing fast matrix operations crucial for scientific computing **

** GNU Fortran compiler - Part of the GNU Compiler Collection (GCC) developed by the Free Software Foundation, required for compiling Fortran source code commonly used in scientific computing applications **

** PNG development library - The reference implementation of the PNG specification, developed by the PNG Development Group, providing support for reading and writing PNG image files **

** CGAL (computational geometry algorithms library) - A comprehensive C++ library developed by INRIA, ETH Zurich, and other European research institutions, providing efficient and reliable geometric algorithms for computational geometry applications **

## ScalFMM Installation

** Clone ScalFMM repository with all its submodules - ScalFMM is a Fast Multipole Method library developed by INRIA (French National Institute for Research in Digital Science and Technology) for N-body simulations, widely used in computational physics and engineering. The --recursive flag ensures all dependent submodules are also downloaded **
```bash
git clone --recursive https://gitlab.inria.fr/solverstack/ScalFMM.git
cd ScalFMM
```

** Create a build directory - This follows CMake best practices (CMake developed by Kitware) by keeping build files separate from source code, making it easier to clean builds and manage multiple build configurations **
```bash
mkdir build
```

** Navigate to the build directory - All compilation will happen in this isolated directory **
```bash
cd build
```

** Configure compilation in Release mode (optimized) - CMake generates the necessary build files with optimization flags enabled for maximum performance in production use **
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
```

** Compile the ScalFMM project - This builds all the libraries and executables using GNU Make (part of the GNU Project by the Free Software Foundation), which may take several minutes depending on your system's performance **
```bash
make -j $(nproc)
```

** Install ScalFMM headers **
```bash
sudo make install
cd ../..
```

## Install ANN from sources

If ANN is not already installed, download it and install it as follows:

```shell
wget https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
tar xzf ann_1.1.2.tar.gz
cd ann_1.1.2/
sed -i 's/CFLAGS =.* -O3/& -std=c++98/' Make-config
make linux-g++
sudo cp lib/libANN.a /usr/local/lib/
sudo cp include/ANN/ANN.h /usr/local/include/
cd ..
```

## Diwina Path Integral Project Installation

** Clone the main branch of the Diwina Path Integral project specifically adapted to work with the ScalFMM library from INRIA **
```bash
git clone https://github.com/JC-Toussaint/Diwina_path_integral.git
```

** Navigate to the Diwina_path directory - Enter the project directory to begin configuration and compilation **
```bash
cd Diwina_path_integral
```

** Configure compilation of the Diwina project in Release mode - This sets up the build system using CMake for optimal performance, linking against the previously compiled ScalFMM libraries **
```bash
cmake . -DCMAKE_BUILD_TYPE=Release
```
## Compile Diwina_path_integral Project

```bash
make -j $(nproc)
```
### installation
```bash
sudo make install
```

## Python Environment Setup

** Create a Python virtual environment - This uses Python's built-in venv module (part of the Python standard library developed by the Python Software Foundation) to isolate Python packages for this project, preventing conflicts with system Python packages **
```bash
python3 -m venv .venv
```

** Activate the virtual environment - This modifies your shell environment to use the isolated Python installation and package directory **
```bash
source .venv/bin/activate
```

** Install meshio pyyaml scipy PyQt5 pyvista pyvistaqt python modules **
```bash
pip3 install gmsh meshio pyyaml scipy PyQt5 pyvista pyvistaqt
```

** meshio (library for reading/writing mesh files) - Developed by Nico Schlömer, this library is essential for handling various mesh file formats commonly used in finite element analysis **

** PyYAML (for processing YAML files) - Developed by Kirill Simonov and Ingy döt Net, this library enables parsing and generation of YAML configuration files, commonly used for application settings and data serialization **

** SciPy (scientific computing library) - Developed by Travis Oliphant, Pearu Peterson, and Eric Jones, built on NumPy, this library provides algorithms for optimization, linear algebra, integration, and other mathematical tasks **

** PyQt5 (graphical user interface) - Developed by Riverbank Computing, Python bindings for the Qt5 framework (originally developed by Trolltech, now The Qt Company), enabling creation of desktop applications **

** PyVista and PyVistaQt (3D visualization) - Developed by Bane Sullivan and Alex Kaszynski, Python
library for 3D plotting and mesh analysis, built on top of VTK (Visualization Toolkit developed by
Kitware) and its Qt GUI **


# Getting started

A configuration template file can be generated, giving you all the fields you need to fill in to run a simulation.  To do this, enter :
./pathIntegral --print-defaults > settings.yml

The resulting file is in YAML format. It has lots of comments (everything after a # sign) documenting the file structure 
and the meaning of the individual settings.
With the comments stripped, the file looks like this:

```yaml
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
```

The previous settings file corresponds to the case of a sphere uniformly magnetized along the Ox axis, with a radius of r=100nm. 
The Gmsh mesh generator has already been used to generate the mesh of the sphere. We define two regions: sphere_volume and sphere_surface. 
For the sphere_volume region, you need to specify the material’s magnetization in Tesla, as well as the photon absorption lengths 
for left- and right-circular polarizations in nanometers.
The micromagnetic system can be rotated along two axes (i.e., axis1 and axis2) by specifying the rotation angles (i.e., angle1 and angle2) in degrees.
The filled flag is not relevant for a sphere entirely made of magnetic material. However, it is
relevant for tubes filled with a non-magnetic material, such as copper (TO CHECK).
The CE and V parameters account for electrostatic effects on the holographic phase. CE depends on the type of microscope and the electron energy, while the inner potential V depends on the type of material.
Finally, the holographic phase image is calculated at any point on a grid called the detector. Since its dimensions must be larger than those of the projected system, the zoom parameter must be less than 1 in the current version. The resolution is determined by the meshSize parameter.

## Launch the Path Integral Software

In a feeLLGood directory, launch path_integral GUI:

```bash
cd examples/uni_sphere
yml-editor-extended.py settings.yml
```

yml-editor-extended.py reads the yaml settings in the input file settings.yml and generates seven files:
```bash
sim_M_integrals.out
sim_Holo.out
sim_MZ.png
sim_PATH_LENGTH.png
sim_STXM_XMCD.png
sim_HOLO_PHASE.png
sim_HOLO_PHASE_RGB.png
```

The .out files are text:

    sim_M_integrals.out – Columns represent the node number, its x and y coordinates, a boolean value indicating whether the beam propagates through materials, the length of material traversed, the integrals of Mx, My, and Mz over the beam path, and the STXM intensity contrast.
    sim_Holo.out – Columns represent the node number, its x and y coordinates, a boolean value indicating whether the beam propagates through materials, the length of material traversed, and the holographic phase in radians, taking into account both electrostatic and magnetic contributions.

The .png files are image files:

    sim_MZ.png – A 16-bit grayscale map of the integral of Mz.
    sim_PATH_LENGTH.png – A 16-bit grayscale map of the material path length.
    sim_STXM_XMCD.png – A 16-bit grayscale map of the STXM contrast.
    sim_HOLO_PHASE.png – A 16-bit grayscale map of the holographic phase.
    sim_HOLO_PHASE_RGB.png – An 8-bit indexed RGB map of the holographic phase.

The user might visualize them typing:
```bash
png_viewer.py
```

> **Warning**: run_path_integral.py uses the default values from the `default-settings.yml` file if it does not find them in `settings.yml`. In particular, the absorption coefficients are set to 0.01 and 0.018 (nm⁻¹). 
