#!/usr/bin/env python3

import os
import sys
import subprocess


# check that argument is valid
if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} settings.yml")
    sys.exit(1)

settings_file = sys.argv[1]

# Launch pathIntegral GUI
subprocess.run(["yml-editor-extended.py", settings_file], check=True)

# Extract file base name
basename = os.path.splitext(settings_file)[0]

# Build new name
filename = f"{basename}_ray.yml"

# Launch pathIntegral
subprocess.run(["pathIntegral", filename], check=True)
