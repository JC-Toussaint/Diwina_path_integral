# FeeLLGood Micromagnetic Imaging Toolbox (fmit)

## fmit-calculator
    python executable that launch the main application with a graphic interface. It builds some
settings, save them in yaml with some 3D visualization for simulating phase images and XMCD images.
It is a front-end for a command line executable pathIntegral. The GUI uses PyQt5 and pyvista(VTK) for 3D rendering.

### How To Force Nvidia GPU Usage when using VTK
    When working with Ubuntu systems that have both integrated graphics and discrete Nvidia GPUs
(commonly found in laptops with Nvidia Optimus technology), the VTK OpenGL error often occurs
because the application is trying to use the integrated GPU instead of the dedicated Nvidia GPU.
This is particularly common in hybrid graphics setups.

### Understanding the Problem
    Some drivers, such as NVidia's, are generously fault tolerant and run cleanly on code that's not
following the GL spec to the letter VTK/OpenGL Driver Information - KitwarePublic
but the issue here is usually related to GPU selection rather than OpenGL compliance.

### Solution: Environment Variables Method
    The most effective approach is to use Nvidia's PRIME render offload system through environment variables.
Set these environment variables before running your Python script:
```bash
export __NV_PRIME_RENDER_OFFLOAD=1
export __GLX_VENDOR_LIBRARY_NAME=nvidia
export __VK_LAYER_NV_optimus=NVIDIA_only
```
You might optionaly enable GPU debugging information
```bash
export VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/nvidia_icd.json
```

And finally run your pyvista(VTK) Python script

## phase_contours.py
    iso values of the holographic phase (using matplotlib)

## fresnel_from_png_phase.py
    post-processing tool to compute Fresnel intensities (Work In Progress)

## fg-yml.py
    'tree like' text display of a yaml file (command line)

## msh-explorer.py
    informations on a mesh file: print the content of a mesh file, and smallest tetrahedron and its
    volume (command line tool)

## fmit-png-extractor
    tool to get some metadatas from an png file from pathIntegral (command line tool)

## fmit-png-viewer
    tool to display pngs output from pathIntegral (using matplotlib)

