#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess

# Définir la variable d'environnement DIWINA_PATH
diwina_path = os.path.expanduser("~/Devel/Diwina_path_integral")
os.environ["DIWINA_PATH"] = diwina_path

# Copier la colormap pour l'holographie
holo_scale_src = os.path.join(diwina_path, "tools", "HoloScale.tsv")
holo_scale_dst = os.path.join(os.getcwd(), "HoloScale.tsv")
shutil.copy(holo_scale_src, holo_scale_dst)

# Vérifie que l'argument est fourni
if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} settings.yml")
    sys.exit(1)

settings_file = sys.argv[1]

# Vérifie que le fichier indiqué existe
if not os.path.isfile(settings_file):
    print(f"Error: file '{settings_file}' does not exist.")
    sys.exit(1)

# Lancer l'éditeur YAML
yml_editor_path = os.path.join(diwina_path, "tools", "yml-editor.py")
subprocess.run([yml_editor_path, settings_file], check=True)

# Extraire le nom sans l'extension
basename = os.path.splitext(settings_file)[0]

# Construire le nouveau nom
filename = f"{basename}_ray.yml"

# Exécuter path_integral
path_integral_exe = os.path.join(diwina_path, "pathIntegral")
subprocess.run([path_integral_exe, filename], check=True)

# Exécuter png_viewer
png_viewer = os.path.join(diwina_path, "tools", "png_viewer.py")
subprocess.run([png_viewer], check=True)

