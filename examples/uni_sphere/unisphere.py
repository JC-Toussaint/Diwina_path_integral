#!/usr/bin/env python3

import os
import sys
import subprocess

'''
you have to install meshMaker feellgood module with pip install in your venv
typing pip install ../FeeLLGood/python-modules
with appropriate path
'''

from feellgood import meshMaker
r = 100 
meshSize = 4

sphere = meshMaker.Ellipsoid(r,r,meshSize,"sphere_surface","sphere_volume")
sphere.make("sphere.msh")

