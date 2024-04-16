#!/usr/bin/python3

import os
import sys
import subprocess
from feellgood import meshMaker # you have to install this module with sudo make install from feellgood public directory

r = 100 
meshSize = 4

sphere = meshMaker.Ellipsoid(r,r,meshSize,"sphere_surface","sphere_volume")
sphere.make("smoke-test.msh")

val = subprocess.run(["./pathIntegral","smoke-test.yml"], text=True)
success = not(val.returncode)

# Use ANSI colors if printing to a terminal.
if sys.stdout.isatty():
    green  = "\x1b[1;32m"
    red    = "\x1b[1;31m"
    normal = "\x1b[m"
else:
    green  = ""
    red    = ""
    normal = ""

# Report success status.
if success:
    print(f"test {green}PASSED{normal}")
    sys.exit(0)
else:
    print(f"test {red}FAILED{normal}")
    sys.exit(1)
