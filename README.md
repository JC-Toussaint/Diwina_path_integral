# pathIntegral from Diwina project
## based on FeeLLGood – A micromagnetic solver

pathIntegral is a STXM-XMCD and electron holography magnetic phase image simulator using feeLLGood simulations inputs. feeLLGood is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation, developped by JC Toussaint & al. The code is being modified without any warranty it works. A dedicated website can be found [here][]  

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

Install Windows Subsystem for Linux (WSL) - This enables running a Linux environment directly on Windows
```bash
wsl --install
```

## For Linux or WSL2, initial Setup 

Update the list of available packages and upggrade to their latest versions
```bash
sudo apt update
sudo apt upgrade
```

## System Dependencies Installation


```bash
sudo apt install wget git cmake pkgconf g++ gfortran python3-venv
sudo apt install libgmsh-dev libeigen3-dev libtbb-dev libyaml-cpp-dev duktape-dev
sudo apt install libfftw3-dev libopenblas-dev libpng-dev libcgal-dev
```

## ScalFMM Installation

Clone ScalFMM repository with all its submodules - ScalFMM is a Fast Multipole Method library developed by INRIA (French National Institute for Research in Digital Science and Technology) for N-body simulations, widely used in computational physics and engineering. The --recursive flag ensures all dependent submodules are also downloaded
```bash
git clone --recursive https://gitlab.inria.fr/solverstack/ScalFMM.git
cd ScalFMM
```

Create a build directory and navigate to it - All compilation will happen in this isolated directory
```bash
mkdir build
cd build
```

Configure compilation in Release mode (optimized) - CMake generates the necessary build files with optimization flags enabled for maximum performance in production use
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
```

Compile the ScalFMM project - This builds all the libraries and executables using GNU Make (part of the GNU Project by the Free Software Foundation)
```bash
make -j $(nproc)
```

Install ScalFMM headers
```bash
sudo make install
cd ../..
```

## Diwina Path Integral Project Installation

Clone the main branch of the Diwina Path Integral project specifically adapted to work with the ScalFMM library from INRIA
```bash
git clone https://github.com/JC-Toussaint/Diwina_path_integral.git
```

Navigate to the Diwina_path directory - Enter the project directory to begin configuration and compilation
```bash
cd Diwina_path_integral
```

Configure compilation of the Diwina project in Release mode - This sets up the build system using CMake for optimal performance, linking against the previously compiled ScalFMM libraries
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

Create a Python virtual environment. This uses Python's built-in venv module
```bash
python3 -m venv .venv
```

Activate the virtual environment
```bash
source .venv/bin/activate
```

Install meshio pyyaml scipy PyQt5 pyvista pyvistaqt python modules
```bash
pip3 install meshio pyyaml scipy PyQt5 pyvista pyvistaqt pypng
```

# Getting started

A configuration template file can be generated, giving you all the fields you need to fill in to run a simulation.
To do this, from command line, enter :
```bash
./pathIntegral --print-defaults > settings.yml
```
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
initial_magnetization: [1, 0, 0]
rotations:
  angle1: 0
  axe1: [1, 0, 0]
  angle2: 0
  axe2: [0, 0, 1]
filled: false
sensor:
  relative_size: 2
  pixel_size: 1.0e-9
```

The previous settings file corresponds to the case of a sphere uniformly magnetized along the Ox axis, with a radius of r=100nm. 
The mesh generator has already been used to generate the mesh of the sphere. We define two regions: sphere_volume and sphere_surface. 
For the sphere_volume region, you need to specify the material’s magnetization in Tesla, as well as the photon absorption lengths 
for left- and right-circular polarizations in nanometers.
The micromagnetic system can be rotated along two axes (i.e., axis1 and axis2) by specifying the rotation angles (i.e., angle1 and angle2) in degrees.
The filled flag is not relevant for a sphere entirely made of magnetic material. However, it is
relevant for tubes filled with a non-magnetic material, such as copper (TO CHECK).
Finally, the holographic phase image is calculated at any point on a grid called the sensor. Since its dimensions must be larger than those of the projected system, the sensor's relative size must be larger than 1. The resolution is determined by the pixel\_size parameter.

## Launch pathIntegral graphic interface

In a feeLLGood directory, launch fmit-calculator, the graphic interface to pathIntegral:

```bash
cd examples/uni_sphere
fmit-calculator unisphere.yml
```
In that example, fmit-calculator reads the yaml settings in the input file unisphere.yml and generates six files:
```bash
sim_M_integrals.out
sim_STXM_HOLO.out
sim_MZ.png
sim_PATH_LENGTH.png
sim_STXM_XMCD.png
sim_HOLO_PHASE.png
```

The .out files are text (tsv) files:

    sim_M_integrals.out – Columns represent the node number, its x and y coordinates, a boolean value indicating whether the beam propagates through materials, the length of material traversed, the integrals of Mx, My, and Mz over the beam path, and the STXM intensity contrast.
    sim_STXM_HOLO.out – Columns represent the node number, its x and y coordinates, a boolean value,
    crossed thickness, integrated Mx, integrated My, integrated Mz, XMCD conrast, holographic phase
    indicating whether the beam propagates through materials, the length of material traversed, and
    the holographic phase in radians, taking into account magnetic contributions only.

The .png files are image files:

    sim_MZ.png – A 16-bit grayscale map of the integral of Mz.
    sim_PATH_LENGTH.png – A 16-bit grayscale map of the material path length.
    sim_STXM_XMCD.png – A 16-bit grayscale map of the STXM contrast.
    sim_HOLO_PHASE.png – A 16-bit grayscale map of the holographic phase.

The user might visualize them typing:
```bash
fmit-png-viewer
```

Basic information about the PNG files can be obtained with the standard `file` utility:
```bash
file sim_HOLO_PHASE.png
```
This should print a line similar to the following:
```
sim_HOLO_PHASE.png: PNG image data, 487 x 487, 16-bit grayscale, non-interlaced
```
The metadata within these files is show by `fmit-calculator`. It can also be retrieved from the command line using the `pngmeta` utility, which has to be installed first:
```bash
sudo apt install pngmeta
```
Then run:
```
pngmeta sim_HOLO_PHASE.png
```
This prints the information needed for getting the physical quantities from the image (pixel size and mapping from pixel values to the imaged quantity), as well as the YAML file that was fed to `pathIntegral` in order to build the images.

**Warning**: fmit-calculator uses the default values from the `default-settings.yml` file if it does not
 find them in the yaml input file. The absorption coefficient values are set to 0.01 and 0.018 (nm⁻¹). 
