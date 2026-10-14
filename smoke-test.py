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

val = subprocess.run(["src/pathIntegral","smoke-test.yml"], text=True)

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
values = [num_ref_line, -2.6e-8, 2.6e-8, 1, 1.86e-7, 0, 0, 0.1479698160, 0.6314219671]
okVals = False
tol = 1e-8
if(val.returncode == 0):
    with open("smoke-test_M_integrals.out",'r') as f:
        nb_lig = 0;
        for line in f:
            if line.startswith("#"):
                continue
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
    if (val.returncode != 0):
        print(f"pathIntegral returned {val.returncode}")
    else:
        names = ["idx", "x", "y", "in", "path_length", "Mx_integral", "My_integral", "Mz_integral", "contrast"]
        print(f"line {num_ref_line} of smoke-test_M_integrals.out (tolerance {tol}):")
        for i in range(len(values)):
            ok = (i < len(data)) and (abs(values[i] - float(data[i])) < tol)
            got = data[i] if i < len(data) else "missing"
            print(f"  {names[i]:12s} expected {values[i]:<14} got {got:<20} {'ok' if ok else 'MISMATCH'}")
    print(f"test {red}FAILED{normal}")
    sys.exit(1)
