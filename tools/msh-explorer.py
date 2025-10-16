#!/usr/bin/env python3

import sys
import numpy as np
import meshio

def load_mesh(filename):
    try:
        # Lecture du maillage avec meshio
        mesh = meshio.read(filename)

        # Coordonnées des nœuds (x, y, z)
        points = mesh.points[:, :3]

        # Initialisation des listes pour triangles et tétraèdres
        triangles = []
        tetrahedra = []

        # Parcours des blocs de cellules
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                triangles.append(cell_block.data)
            elif cell_block.type == "tetra":
                tetrahedra.append(cell_block.data)

        # Fusionner les blocs s'ils existent
        triangles  = np.vstack(triangles)  if triangles  else np.empty((0, 3), dtype=int)
        tetrahedra = np.vstack(tetrahedra) if tetrahedra else np.empty((0, 4), dtype=int)

        return points, triangles, tetrahedra

    except Exception as e:
        print(f"[ERREUR] Échec de la lecture du fichier '{filename}': {e}", file=sys.stderr)
        sys.exit(1)

# Exemple d’utilisation
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python msh-explorer.py fichier.msh", file=sys.stderr)
        sys.exit(1)

    filename = sys.argv[1]
    points, triangles, tetrahedra = load_mesh(filename)

    print(f"[INFO] {len(points)} nœuds chargés.")
    print(f"[INFO] {len(triangles)} triangles trouvés.")
    print(f"[INFO] {len(tetrahedra)} tétraèdres trouvés.")
    
    for i, p in enumerate(points):
    	print(f"noeud {i} : {p}")
    	
    # 1-based gmsh 3895 3963 4018 4047 
	# Définition des points
    A = points[3894, :]
    B = points[3962, :]
    C = points[4017, :]
    D = points[4046, :]

	# Calcul des vecteurs
    AB = B - A
    AC = C - A
    AD = D - A

	# Calcul du volume via le produit mixte
    volume = abs(np.dot(AB, np.cross(AC, AD))) / 6

    print(f"Le volume du tétraèdre est : {volume:.6f} unités cubiques")


