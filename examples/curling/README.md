# Curling mode in a cylinder

This example simulates a magnetic cylinder (diameter 200 nm, length 1000 nm) with easy-plane
anisotropy, starting from the curling mode, then computes its holographic phase and XMCD images
with pathIntegral.

The workflow has four steps:

1. generate the mesh with `fg_cyl_sample.py`,
2. complete `settings.yml`,
3. compute the initial magnetization `config.sol` with `fg-gensol_ai`,
4. relax the magnetization with feeLLGood, then image it with pathIntegral.

Each parameter is set in a single place, so that the files cannot become inconsistent:

- the geometry (diameter, length, mesh size) on the command line of `fg_cyl_sample.py`: it is
  then stored in the mesh, and `fg-gensol_ai` deduces the radius of the cylinder from the mesh;
- the material constants, the name of the mesh file and the region names in `settings.yml`,
  which is read by feeLLGood, `fg-gensol_ai`, pathIntegral and fmit-calculator.

| File | Role |
|------|------|
| `fg_cyl_sample.py` | generates the cylinder mesh (gmsh Python API), and writes its name in `settings.yml` |
| `settings.yml` | input file of feeLLGood, `fg-gensol_ai`, pathIntegral and fmit-calculator |
| `fg-gensol_ai` | solves the curling profile and writes `config.sol` |
| `init_cfg.py` | magnetization field of each volume region, called by `fg-gensol_ai` |
| `config.sol` | initial magnetization (feeLLGood `.sol` format), one line per mesh node |

## Prerequisites

- `feellgood` and `pathIntegral` installed (see the main README),
- the Python modules `numpy`, `scipy`, `matplotlib`, `meshio`, `gmsh` and `pyyaml`:

```bash
pip install numpy scipy matplotlib meshio gmsh pyyaml
```

## 1. Generate the mesh

```bash
cd examples/curling
python3 fg_cyl_sample.py --diameter 200 --length 1000 --mesh_size 5 --settings settings.yml
```

This writes `cyl_D200_L1000_h5_fg.msh` (about 190 000 nodes and 1 090 000 tetrahedra, in about
one minute). The mesh file is not stored in git, so this step is required.

- With `--settings settings.yml`, the name of the mesh file is written in `settings.yml`
  (`mesh: filename:`); only this line is modified, the comments are kept. The script also warns
  if a region of `settings.yml` is not a physical group of the mesh. Without this option,
  `settings.yml` is not modified.

- The mesh is written in the gmsh ASCII format 4.1 by default, or 2.2 with `--version 2`.
  feeLLGood, pathIntegral and `fg-gensol_ai` read both formats and give the same results.
- Lengths are in nanometers; this is why `length_unit: 1e-9` in the settings files.
- The disk lies in the plane z = -L/2 and is extruded along z: the cylinder axis is z and
  z ranges from -L/2 to +L/2.
- The mesh size should be smaller than the exchange length
  l_ex = sqrt(2A / (µ0 Ms²)) ≈ 5.7 nm for this material.

The mesh contains four physical groups, whose names are used in the settings files:

| Physical group | Dimension | Description |
|----------------|-----------|-------------|
| `cyl_volume` | 3 | volume of the cylinder |
| `side_surface` | 2 | lateral surface |
| `inlet` | 2 | end face at z = -L/2 |
| `outlet` | 2 | end face at z = +L/2 |

Running the script without argument prints its help.

## 2. Complete `settings.yml`

`settings.yml` is a commented feeLLGood template. The settings to adapt to the mesh are:

```yaml
mesh:
  filename: cyl_D200_L1000_h5_fg.msh   # written by fg_cyl_sample.py --settings
  length_unit: 1e-9                    # mesh lengths are in nm
  volume_regions:
    cyl_volume:                        # one entry per physical volume
      Ae: 1.3e-11                      # exchange constant (J/m)
      Js: 1.00531                      # µ0 Ms (T)
      K: -1.e5                         # uniaxial anisotropy (J/m³), K < 0: easy plane
      uk: [0, 0, 1]                    # anisotropy axis = cylinder axis
      alpha_LLG: 0.5
  surface_regions:                     # one entry per physical surface
    side_surface:
      suppress_charges: false
    inlet:
      suppress_charges: false
    outlet:
      suppress_charges: false

initial_magnetization: config.sol      # written by fg-gensol_ai
```

- The region names must be exactly those of the physical groups of the mesh. `fg-gensol_ai` stops
  if the volume regions of the mesh and of `settings.yml` differ.
- `suppress_charges: true` removes the magnetic charges of a surface, to mimic an infinite
  medium. Set it on `inlet` and `outlet` to simulate an infinite cylinder.
- `fg-gensol_ai` reads the material constants `Ae` and `K` from this file (see step 3).
- Check also the `outputs` section: `final_time` (the time at which feeLLGood stops),
  `evol_time_step`, and `mag_config_every`, which sets how often a magnetization file is saved.
  On this mesh, a time step took about 30 s on our test machine, so choose `final_time` accordingly:
  the value 1e-3 s of the template is far too long.
- `demagnetizing_field_solver.nb_threads` and `finite_element_solver.nb_threads` should match
  the number of cores of your computer (0 means all of them for the finite element solver).

Check the file before running a simulation:

```bash
feellgood --verify settings.yml
```

## 3. Compute the initial magnetization

```bash
./fg-gensol_ai settings.yml
```

The script:

1. reads in `settings.yml` the mesh file name, its length unit, and the material constants of the
   volume region: A = `Ae` and Ku = `K`. It deduces the radius R of the cylinder from the mesh, as
   the largest distance of a node to the z axis, and prints these values;
2. solves the 1D Euler-Lagrange equation of the curling mode for the angle ω(ρ) between the
   magnetization and the cylinder axis, with ρ = r/R:

   ω'' = -ω'/ρ + (1/(2ρ²) + κ) sin(2ω),  κ = Ku R² / (2A)

   with ω = 0 on the axis and ω' = 0 on the lateral surface;
3. saves the profile ω(ρ) in `curling_mode_adimensionne.png`;
4. evaluates, at every node of the region `cyl_volume`, the field defined in `init_cfg.py`:

   m = (-sin ω · y/r, sin ω · x/r, cos ω)

   while printing a progress percentage. `define_region_functions(sol_sol, R)` receives the radius
   R in mesh units (nm);
5. writes `config.sol`.

No physical parameter is set in `fg-gensol_ai` or `init_cfg.py`: to change the diameter, generate
a new mesh; to change the material, edit `settings.yml`. The curling profile is computed for a
single material, so all the volume regions must have the same `Ae` and `K`. The anisotropy axis
`uk` should be the cylinder axis z, otherwise a warning is printed.

`config.sol` must be regenerated each time the mesh changes: it has one line per mesh node, and
feeLLGood stops if the node indices of the mesh and of the `.sol` file do not match.

## 4. Run the simulation and compute the images

### Micromagnetic relaxation with feeLLGood

```bash
feellgood settings.yml
```

feeLLGood writes, with the basename `sim` set in `outputs.file_basename`:

- `sim.evol`: evolution of the quantities listed in `evol_columns` (average magnetization,
  energies),
- `sim_iter<N>.sol`: magnetization configurations, saved every `mag_config_every` steps.

### Phase and XMCD images with pathIntegral

pathIntegral reads `settings.yml` too: it ignores the feeLLGood settings, and uses the mesh, the
volume regions and the last sections of the file:

```yaml
initial_magnetization: config.sol     # magnetization to image
rotations:                            # rotations of the sample, in degrees
  angle1: 90
  axe1: [1, 0, 0]                     # cylinder axis perpendicular to the beam
  angle2: 0
  axe2: [0, 0, 1]
filled: false
sensor:
  relative_size: 2                    # sensor size / sample size (≥ 1)
  pixel_size: 0.5e-9                  # in m
```

The number of nodes of the sensor grid grows as 1/pixel_size²: with the cylinder perpendicular to
the beam, `pixel_size: 0.5e-9` gives a 4001 × 4001 grid, and `2e-9` a 1001 × 1001 grid.

- To image the relaxed state instead of the initial curling mode, set `initial_magnetization` to
  one of the `sim_iter<N>.sol` files written by feeLLGood.
- The XMCD absorption lengths `lbp` and `lbm` are not set in this file: pathIntegral uses its
  default values (0.01 and 0.018 nm⁻¹). Add them under `cyl_volume` to change them.

Then run pathIntegral from the command line:

```bash
pathIntegral settings.yml
```

or with the graphic interface, which helps choosing the rotations:

```bash
fmit-calculator settings.yml
```

fmit-calculator saves the parameters chosen in its interface in `settings_ray.yml` before running
pathIntegral. This file is generated, it is not stored in git.

pathIntegral writes `sim_M_integrals.out`, `sim_STXM_HOLO.out`, `sim_MZ.png`,
`sim_PATH_LENGTH.png`, `sim_STXM_XMCD.png`, `sim_HOLO_PHASE.png` and `sim_HOLO_PHASE_RGB.png`.
The `.out` and grayscale `.png` files are described in the main README; `sim_HOLO_PHASE_RGB.png`
is a color version of the holographic phase. The images can be displayed with `fmit-png-viewer`.
