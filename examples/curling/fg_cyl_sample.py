"""Mesh of a cylinder for st-feeLLGood (gmsh Python API).

Disk of diameter D in the plane z = -L/2, extruded along z over the length L;
lengths in nm (length_unit: 1e-9 in settings.yml). Physical groups: inlet (z = -L/2),
outlet (z = +L/2), side_surface, cyl_volume. Output: cyl_D<D>_L<L>_h<h>_fg.msh.

    python3 fg_cyl_sample.py --diameter 200 --length 2000 --mesh_size 5 [--version 2|4]
"""
import gmsh
import sys
import argparse

# -----------------------------
# command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Generation of a meshed cylinder")

parser.add_argument("--diameter", type=float, default=200.0,
                    help="Diameter of the cylinder (nm), default = 200")
parser.add_argument("--length", type=float, default=2000.0,
                    help="Length of the cylinder (nm), default = 2000")
parser.add_argument("--mesh_size", type=float, default=5.0,
                    help="Characteristic mesh size (nm), default = 5")
parser.add_argument("--version", type=int, choices=[2, 4], default=4,
                    help="Version of the Gmsh ASCII mesh format (2 or 4), default = 4")

# Print the help if no argument is given
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

# Initialisation of Gmsh
gmsh.initialize()
gmsh.model.add("cylinder")

# -----------------------------
# geometrical parameters
# -----------------------------
D = args.diameter
R = D / 2.0  # radius
L = args.length

# wished mesh size
lc = args.mesh_size

# base disc in the plane z = -L/2
z_base = -L / 2.0

# centre of the disc
p1 = gmsh.model.geo.addPoint(0, 0, z_base, lc)

# points on the circle (4 points to define the circle)
p2 = gmsh.model.geo.addPoint(R, 0, z_base, lc)
p3 = gmsh.model.geo.addPoint(0, R, z_base, lc)
p4 = gmsh.model.geo.addPoint(-R, 0, z_base, lc)
p5 = gmsh.model.geo.addPoint(0, -R, z_base, lc)

# circular arcs
c1 = gmsh.model.geo.addCircleArc(p2, p1, p3)
c2 = gmsh.model.geo.addCircleArc(p3, p1, p4)
c3 = gmsh.model.geo.addCircleArc(p4, p1, p5)
c4 = gmsh.model.geo.addCircleArc(p5, p1, p2)

# curve loop forming the circle
loop = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])

# surface of the disc (inlet)
inlet_tag = gmsh.model.geo.addPlaneSurface([loop])

# synchronisation required before the extrusion
gmsh.model.geo.synchronize()

# extrusion of the disc along the z axis
# extrude returns a list of (dimension, tag) pairs
extrusion = gmsh.model.geo.extrude([(2, inlet_tag)], 0, 0, L)

# final synchronisation
gmsh.model.geo.synchronize()

# printout for debugging
print("Result of the extrusion:")
for i, entity in enumerate(extrusion):
    print(f"  extrusion[{i}] = {entity} (dim={entity[0]}, tag={entity[1]})")

# tags of the entities
# extrusion[0] = (dim=2, tag) -> outlet surface (top)
# extrusion[1] = (dim=3, tag) -> volume
# extrusion[2:] = (dim=2, tag) -> side surfaces

outlet_tag = extrusion[0][1]
volume_tag = extrusion[1][1]
side_surfaces = [extrusion[i][1] for i in range(2, len(extrusion))]

print(f"\ninlet_tag = {inlet_tag}")
print(f"outlet_tag = {outlet_tag}")
print(f"volume_tag = {volume_tag}")
print(f"side_surfaces = {side_surfaces}")

# named physical groups (regions of settings.yml)
gmsh.model.addPhysicalGroup(2, [inlet_tag], name="inlet")
gmsh.model.addPhysicalGroup(2, [outlet_tag], name="outlet")
gmsh.model.addPhysicalGroup(2, side_surfaces, name="side_surface")
gmsh.model.addPhysicalGroup(3, [volume_tag], name="cyl_volume")

# 3D mesh generation
gmsh.model.mesh.generate(3)

# output format: ASCII, chosen version
if args.version == 2:
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
else:  # version 4
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.Binary", 0)  # 0 = ASCII, 1 = binary

# file name with D, L and h as integers
filename = f"cyl_D{int(D)}_L{int(L)}_h{int(lc)}_fg.msh"

# save the mesh
gmsh.write(filename)

# open the graphical interface (optional)
#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

# finalisation
gmsh.finalize()
