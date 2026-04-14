#!/usr/bin/python3

from feellgood.meshMaker import Tube

Tube(
    radius1 = 50,  # inner radius [nm]
    radius2 = 70,  # outer radius [nm]
    length = 200,  # [nm]
    mesh_size = 4,
    surfName = "surface",
    volName = "volume").make("tube.msh")
