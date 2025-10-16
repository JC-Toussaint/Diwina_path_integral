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

def vol_tetra(points, tet):
	
	# Définition des points
    A = points[tet[0], :]
    B = points[tet[1], :]
    C = points[tet[2], :]
    D = points[tet[3], :]

	# Calcul des vecteurs
    AB = B - A
    AC = C - A
    AD = D - A

	# Calcul du volume via le produit mixte
    return abs(np.dot(AB, np.cross(AC, AD))) / 6

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
    	
    tet_id = 0
    tet_min = tetrahedra[0]
    vol_min	= vol_tetra(points, tet_min)

    for i, t in enumerate(tetrahedra):
        vol = vol_tetra(points, t)
        print(f"tet {i} : {t} -> vol {vol}")

        # Mise à jour du plus petit volume
        if vol < vol_min:
           vol_min = vol
           tet_min = t  
           tet_id  = i

	# Affichage final du plus petit volume
    print(f"\ntet {tet_id} avec le + petit volume {tet_min} -> vol {vol_min}")
     	   	

