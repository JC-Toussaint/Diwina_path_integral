#!/usr/bin/python3

import gmsh
import os
import sys
import subprocess

# Initialiser Gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# Créer un fichier .geo temporaire
geo_filename = "hollow_sphere.geo"
geo_content = """\
	SetFactory("OpenCASCADE");
	Mesh.CharacteristicLengthMax = 10.0;
	Mesh.CharacteristicLengthMin = 10.0;
	R2 = 200;
	R1 = 100;

	Sphere(1) = {0,0,0, R2};
	Sphere(2) = {0,0,0, R1};
	BooleanDifference{ Volume{1}; Delete;}{ Volume{2}; Delete;}

	Physical Volume("hollow_volume", 10) = {1};
	Physical Surface("hollow_surface", 11) = {2, 3};
"""

# Écrire le fichier .geo
with open(geo_filename, "w") as geo_file:
    geo_file.write(geo_content)

# Lire et générer le maillage avec Gmsh
gmsh.open(geo_filename)
gmsh.model.mesh.generate(3)  # Maillage 3D

# Sauvegarder le maillage au format MSH
msh_filename = "hollow_sphere.msh"
gmsh.write(msh_filename)

# Afficher dans l'interface graphique de Gmsh (optionnel)
#if "-nopopup" not in os.sys.argv:
#    gmsh.fltk.run()

# Nettoyage et terminaison
gmsh.finalize()

print(f"Maillage sauvegardé dans : {msh_filename}")

val = subprocess.run(["../../pathIntegral", "uni_hollow_sphere.yml"], text=True)


# Use ANSI colors if printing to a terminal.
if sys.stdout.isatty():
    green  = "\x1b[1;32m"
    red    = "\x1b[1;31m"
    normal = "\x1b[m"
else:
    green  = ""
    red    = ""
    normal = ""

print(f"{green}FINISHED{normal}")
sys.exit(0)

