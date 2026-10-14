# Examples
In this drectory there are four different examples to illustrate basic use of pathIntegral and
fmit-calculator. Each example comes with a python script to generate a mesh, and a yaml input file to feed either directly pathIntegral from command line, or fmit-calculator if graphic interface is needed to find adequate rotations.

## uni_hollow_sphere
this example is a hollow sphere with uniform magnetization.

## uni_sphere
this example is sphere with uniform magnetization.

## tube
this example is a tube with an analytic ansatz mimicking a rolled Néel domain wall onto the tube.

## curling
this example is a cylinder (diameter 200 nm, length 1000 nm) with easy-plane anisotropy, in the curling mode (see `curling/README.md`). The angle between the magnetization and the cylinder axis is the solution of the 1D Euler-Lagrange equation, solved by `fg-gensol_ai`; `init_cfg.py` maps it onto the mesh nodes.
```
pip install gmsh pyyaml
python3 fg_cyl_sample.py --diameter 200 --length 1000 --mesh_size 5 --settings settings.yml
./fg-gensol_ai settings.yml
pathIntegral settings.yml
```
`settings.yml` is the single input file of feeLLGood, `fg-gensol_ai`, pathIntegral and fmit-calculator: it holds the material constants, the mesh file name (written by `fg_cyl_sample.py --settings`) and the region names. `fg-gensol_ai` deduces the radius of the cylinder from the mesh, and writes the initial magnetization `config.sol`.
