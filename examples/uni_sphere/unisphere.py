#!/usr/bin/env python3

import os
import sys
import subprocess

'''
you have to install meshMaker feellgood module with pip install in your venv
typing pip install ../FeeLLGood/python-modules
with appropriate path
'''

import meshMaker
r = 100 
meshSize = 4

sphere = meshMaker.Ellipsoid(r,r,meshSize,"sphere_surface","sphere_volume")
sphere.make("sphere.msh")

val = subprocess.run(["../../pathIntegral","unisphere.yml"], text=True)

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

