#!/usr/bin/env python3
# Auto-generated configuration module
# Edit the lambda functions below to define magnetization vector field for each region

import numpy as np

def define_region_functions(sol_sol, R):
    """Define magnetization vector field for each region
    sol_sol: continuous solution omega(rho) of the curling profile, rho = r/R
    R: radius of the cylinder, deduced from the mesh, in mesh units (nm)
    """
    import numpy as np

    def curling(x, y, z):
        r = np.sqrt(x**2 + y**2)
        # Sur l'axe, omega(0) = 0 : m est le long de z (evite la division par r = 0)
        if r == 0:
            return [0.0, 0.0, 1.0]

        # 1. On calcule le rayon adimensionnel du point courant
        rho_point = r / R

        # 2. On interpole l'angle omega grâce à la fonction continue sol_sol
        # sol_sol(rho_point)[0] donne la valeur de omega
        omega = sol_sol(rho_point)[0]

        # 3. On génère le vecteur [mx, my, mz] basé sur cet angle
        mx = -np.sin(omega) * (y / r)
        my = +np.sin(omega) * (x / r)
        mz = +np.cos(omega)

        return [mx, my, mz]

    return {
        'cyl_volume': lambda x, y, z: curling(x, y, z),
    }
