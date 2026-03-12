# Detailed Installation Protocol

## For Windows only, installation of WSL2 (Windows Subsystem for Linux)

**Install Windows Subsystem for Linux (WSL) - This enables running a Linux environment directly on Windows without requiring a separate virtual machine or dual-boot setup**
```bash
wsl --install
```

## For Linux or WSL2, initial Setup 

**Update the list of available packages - This refreshes the package database to ensure you have access to the latest versions of software available in the repositories**
```bash
sudo apt update
```

**Upgrade all installed packages to their latest versions - This ensures system security and compatibility by installing the most recent updates and patches**
```bash
sudo apt upgrade
```

## System Dependencies Installation

**Install ParaView (scientific visualization software) - A powerful, open-source application developed by Kitware, Sandia National Laboratories, and Los Alamos National Laboratory for visualizing large datasets, commonly used in scientific and engineering applications for 3D data analysis**
```bash
sudo apt install paraview
```

**Install FFTW library for fast Fourier transforms - A high-performance C library developed by Matteo Frigo and Steven G. Johnson at MIT for computing discrete Fourier transforms, essential for signal processing and numerical analysis**
```bash
sudo apt install libfftw3-dev
```

**Install OpenBLAS (optimized linear algebra library) - A highly optimized implementation of Basic Linear Algebra Subprograms (BLAS) developed by Zhang Xianyi at the Chinese Academy of Sciences, providing fast matrix operations crucial for scientific computing**
```bash
sudo apt install libopenblas-dev
```

**Install GNU Fortran compiler - Part of the GNU Compiler Collection (GCC) developed by the Free Software Foundation, required for compiling Fortran source code commonly used in scientific computing applications**
```bash
sudo apt install gfortran
```

**Install PNG development library - The reference implementation of the PNG specification, developed by the PNG Development Group, providing support for reading and writing PNG image files**
```bash
sudo apt install libpng-dev
```

**Install CGAL (computational geometry algorithms library) - A comprehensive C++ library developed by INRIA, ETH Zurich, and other European research institutions, providing efficient and reliable geometric algorithms for computational geometry applications**
```bash
sudo apt install libcgal-dev
```

## ScalFMM Installation

**Clone ScalFMM repository with all its submodules - ScalFMM is a Fast Multipole Method library developed by INRIA (French National Institute for Research in Digital Science and Technology) for N-body simulations, widely used in computational physics and engineering. The --recursive flag ensures all dependent submodules are also downloaded**
```bash
git clone --recursive https://gitlab.inria.fr/solverstack/ScalFMM.git
```

**Create a build directory - This follows CMake best practices (CMake developed by Kitware) by keeping build files separate from source code, making it easier to clean builds and manage multiple build configurations**
```bash
mkdir build
```

**Navigate to the build directory - All compilation will happen in this isolated directory**
```bash
cd build
```

**Configure compilation in Release mode (optimized) - CMake generates the necessary build files with optimization flags enabled for maximum performance in production use**
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
```

**Compile the ScalFMM project - This builds all the libraries and executables using GNU Make (part of the GNU Project by the Free Software Foundation), which may take several minutes depending on your system's performance**
```bash
make
```

**Install ScalFMM headers**
```bash
sudo make install
```

## Diwina Path Integral Project Installation

**Clone the main branch of the Diwina Path Integral project specifically adapted to work with the ScalFMM library from INRIA**
```bash
git clone git@github.com:JC-Toussaint/Diwina_path_integral.git
```

**Navigate to the Diwina_path directory - Enter the project directory to begin configuration and compilation**
```bash
cd Diwina_path
```

**Configure compilation of the Diwina project in Release mode - This sets up the build system using CMake for optimal performance, linking against the previously compiled ScalFMM libraries**
```bash
cmake . -DCMAKE_BUILD_TYPE=Release
```
## Compile Diwina_path_integral Project

```bash
cmake .
make -j 8
```

## Python Environment Setup

**Create a Python virtual environment - This uses Python's built-in venv module (part of the Python standard library developed by the Python Software Foundation) to isolate Python packages for this project, preventing conflicts with system Python packages**
```bash
python3 -m venv .venv
```

**Activate the virtual environment - This modifies your shell environment to use the isolated Python installation and package directory**
```bash
source .venv/bin/activate
```

**Install meshio (library for reading/writing mesh files) - Developed by Nico Schlömer, this library is essential for handling various mesh file formats commonly used in finite element analysis**
```bash
pip3 install meshio
```

**Install PyYAML (for processing YAML files) - Developed by Kirill Simonov and Ingy döt Net, this library enables parsing and generation of YAML configuration files, commonly used for application settings and data serialization**
```bash
pip3 install pyyaml
```

**Install SciPy (scientific computing library) - Developed by Travis Oliphant, Pearu Peterson, and Eric Jones, built on NumPy, this library provides algorithms for optimization, linear algebra, integration, and other mathematical tasks**
```bash
pip3 install scipy
```

**Install PyQt5 (graphical user interface) - Developed by Riverbank Computing, this provides comprehensive Python bindings for the Qt5 framework (originally developed by Trolltech, now The Qt Company), enabling creation of desktop applications**
```bash
pip3 install PyQt5
```

**Install PyVista (3D visualization) - Developed by Bane Sullivan and Alex Kaszynski, this is a powerful Python library for 3D plotting and mesh analysis, built on top of VTK (Visualization Toolkit developed by Kitware)**
```bash
pip3 install pyvista
```
**Install PyVistaQt - GUI for PyVista**
```bash
pip3 install pyvistaqt
```

## Launch the Path Integral Software

In a feeLLGood directory, launch the path_integral software:

```bash
../../tools/run_path_integral.py settings.yml
```

> **Warning**: The program uses the default values from the `default-settings.yml` file if it does not find them in `settings.yml`. In particular, the absorption coefficients are set to 0.01 and 0.018 (nm⁻¹).
