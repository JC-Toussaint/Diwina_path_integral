#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Lecture du fichier
def read_holo_file(filename):
    """
    Lit le fichier sim_Holo.out et retourne les données
    """
    data = []
    
    with open(filename, 'r') as f:
        for line in f:
            # Ignorer les lignes de commentaires
            if line.startswith('#'):
                continue
            
            # Lire les données
            values = line.split()
            if len(values) == 6:
                idx = int(values[0])
                x = float(values[1])
                y = float(values[2])
                in_val = int(values[3])
                path_length = float(values[4])
                phase = float(values[5])
                data.append([idx, x, y, in_val, path_length, phase])
    
    return np.array(data)

# Lire le fichier
data = read_holo_file('sim_Holo.out')

# Extraire les colonnes
x = data[:, 1]
y = data[:, 2]
phase = data[:, 5]

# Déterminer la grille
nx = int(np.sqrt(len(x) + 0.5))
ny = nx

print(f"Grille détectée: {nx} x {ny} points")
print(f"Phase : [ {np.min(phase)}  {np.max(phase)} ]")

# Restructurer les données sur une grille 2D
X = x.reshape(ny, nx)
Y = y.reshape(ny, nx)
Phase = phase.reshape(ny, nx)

# Exemple : tableau fourni par ton collègue
# colonnes : x, R, G, B
cmap_data = np.array([
    [0.000000, 0.000000, 0.000000, 0.000000], 
    [0.200000, 0.000000, 0.000000, 1.000000],
    [0.400000, 0.000000, 0.666666, 0.000000],  
    [0.600000, 1.000000, 0.000000, 0.000000],  
    [0.800000, 1.000000, 1.000000, 0.000000], 
    [1.000000, 1.000000, 1.000000, 1.000000]
])

# Construire la colormap
cdict = {'red': [], 'green': [], 'blue': []}
for x, r, g, b in cmap_data:
    cdict['red'].append((x, r, r))
    cdict['green'].append((x, g, g))
    cdict['blue'].append((x, b, b))

cmap = LinearSegmentedColormap('my_cmap', segmentdata=cdict, N=256)

plt.figure(figsize=(10, 10))
cs = plt.contour(X, Y, Phase, levels=20, linewidths=1, cmap='jet')
plt.clabel(cs)
plt.show()

# Créer la visualisation des isovaleurs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Graphique 1: Contours remplis (filled contours)
levels = 20  # Nombre de niveaux
contourf = ax1.contourf(X, Y, Phase, levels=levels, cmap=cmap)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Phase - Contours remplis')
ax1.set_aspect('equal')
plt.colorbar(contourf, ax=ax1, label='Phase (rad)')

# Graphique 2: Lignes de contour (isovaleurs)
levels=np.linspace(np.min(phase), np.max(phase), 10);
contour = ax2.contour(X, Y, Phase, levels=levels, colors='black', linewidths=0.5)
ax2.clabel(contour, inline=True, fontsize=8)
#contourf2 = ax2.contourf(X, Y, Phase, levels=levels, cmap='RdBu_r', alpha=0.6)
ax2.set_xlabel('x (m)')
ax2.set_ylabel('y (m)')
ax2.set_title('Phase - Lignes isovaleurs')
ax2.set_aspect('equal')
#plt.colorbar(contourf2, ax=ax2, label='Phase (rad)')

plt.tight_layout()
plt.savefig('phase_isovaleurs.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"\nStatistiques de la phase:")
print(f"  Min: {phase.min():.6f}")
print(f"  Max: {phase.max():.6f}")
print(f"  Moyenne: {phase.mean():.6f}")
print(f"  Écart-type: {phase.std():.6f}")
