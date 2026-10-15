# python_gmsh

Mesh generation scripts using the gmsh Python API.

## fg-mesh-cylinder

`fg-mesh-cylinder` meshes a cylinder. The script is in this directory;
`examples/curling/fg-mesh-cylinder` is a symbolic link to it, so it can be run from `python_gmsh`
or from `examples/curling`.

### Requirements

- the gmsh Python module: `pip install gmsh`,
- pyyaml, with the `--settings` option only: `pip install pyyaml`.

### Usage

Mesh of a cylinder of diameter 200 nm and length 2000 nm, with a mesh size of 5 nm:

```bash
./fg-mesh-cylinder --diameter 200 --length 2000 --mesh_size 5 --version 4
```

This generates the mesh `cyl_D200_L2000_h5_fg.msh` in the gmsh ASCII format, version 4.1. This is
the default format, so `--version 4` can be omitted.

To generate a mesh in the gmsh ASCII format, version 2.2:

```bash
./fg-mesh-cylinder --diameter 200 --length 2000 --mesh_size 5 --version 2
```

To also write the name of the mesh file in a settings file (`mesh: filename:`), add
`--settings settings.yml`. Only this line of the settings file is modified, its comments are kept.
The script also warns about the regions of the settings file that are not physical groups of the
mesh (`cyl_volume`, `inlet`, `outlet`, `side_surface`).

Without any argument, the script prints its help.
