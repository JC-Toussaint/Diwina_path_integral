#!/usr/bin/env python3

import os
import sys
import subprocess

# Définir la variable d'environnement DIWINA_PATH
#diwina_path = os.path.expanduser("~/Devel/Diwina_path_integral")
#os.environ["DIWINA_PATH"] = diwina_path

# Vérifie que l'argument est fourni
if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} settings.yml")
    sys.exit(1)

settings_file = sys.argv[1]

# Lancer l'éditeur YAML
yml_editor_path = "yml-editor-extended.py"   #os.path.join(diwina_path, "tools", "yml-editor-extended.py")
subprocess.run([yml_editor_path, settings_file], check=True)

# Extraire le nom sans l'extension
basename = os.path.splitext(settings_file)[0]

# Construire le nouveau nom
filename = f"{basename}_ray.yml"

# Exécuter path_integral
path_integral_exe = "pathIntegral" #os.path.join(diwina_path, "pathIntegral")
subprocess.run([path_integral_exe, filename], check=True)

# Exécuter png_viewer
png_viewer = "png_viewer.py" #os.path.join(diwina_path, "tools", "png_viewer.py")
subprocess.run([png_viewer], check=True)

