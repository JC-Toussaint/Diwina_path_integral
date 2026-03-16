#!/usr/bin/python3

import os
import sys
import subprocess

# you have to install feellgood module with sudo make install from feellgood public directory
from feellgood import meshMaker
r = 100
meshSize = 4

sphere = meshMaker.Ellipsoid(r,r,meshSize,"sphere_surface","sphere_volume")
sphere.make("smoke-test.msh")

val = subprocess.run(["./pathIntegral","smoke-test.yml"], text=True)

# Use ANSI colors if printing to a terminal.
if sys.stdout.isatty():
    green  = "\x1b[1;32m"
    red    = "\x1b[1;31m"
    normal = "\x1b[m"
else:
    green  = ""
    red    = ""
    normal = ""

num_ref_line = 70000
values = [num_ref_line, -2.80e-08, -2.80e-08, 1, 1.836e-07, 0, 0, 0.14613901, 0.6258561806]
okVals = False
tol = 1e-8
if(val.returncode == 0):
    with open("smoke-test_M_integrals.out",'r') as f:
        nb_lig = 0;
        for line in f:
            if (nb_lig <num_ref_line):
                nb_lig += 1
                pass
            else:
                break
        data = line.split()
    okVals = True
    for i in range(len(data)):
        okVals = ( okVals and ( abs(values[i] - float(data[i]) ) < tol ) )

    # Report success status.
if ((val.returncode == 0) and okVals):
    print(f"test {green}PASSED{normal}")
    sys.exit(0)
else:
    print(f"test {red}FAILED{normal}")
    sys.exit(1)
