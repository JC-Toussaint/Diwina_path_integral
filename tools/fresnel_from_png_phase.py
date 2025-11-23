#!/usr/bin/env python3
"""
charge un PNG 16-bit contenant une phase,
extrait uniquement les métadonnées utiles, simule Fresnel,
et sauvegarde l'image de sortie dans fresnel.png.
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import argparse
import re

# Constantes
h = 6.62607015e-34
e = 1.602176634e-19
m0 = 9.10938356e-31
c = 299792458.0

def relativistic_wavelength(eV):
    E = eV * e
    p = np.sqrt(E**2 + 2*E*m0*c**2) / c
    return h / p

def find_meta_value(meta, patterns):
    """Cherche une clé correspondant à une regex dans meta."""
    for key in meta:
        key_norm = key.lower().strip()
        for pat in patterns:
            if re.search(pat, key_norm):
                return key, meta[key]
    return None, None

def safe_float(x):
    try:
        if isinstance(x, bytes):
            x = x.decode('utf-8', errors='ignore')
        return float(str(x).strip())
    except Exception:
        return None

def load_phase_png_with_metadata(path_png):
    im = Image.open(path_png)
    meta = dict(im.info)  # métadonnées PNG

    # Recherche des valeurs utiles
    pxw_key, pxw_val = find_meta_value(meta, [r'pixel *width', r'pixel.*x', r'physical.*pixel.*width'])
    pxh_key, pxh_val = find_meta_value(meta, [r'pixel *height', r'pixel.*y', r'physical.*pixel.*height'])
    pxunit_key, pxunit_val = find_meta_value(meta, [r'pixel *unit', r'pixel.*unit', r'unit *pixel'])

    vmin_key, vmin_val = find_meta_value(meta, [r'value *min', r'min *value', r'phase *min'])
    vmax_key, vmax_val = find_meta_value(meta, [r'value *max', r'max *value', r'phase *max'])
    vunit_key, vunit_val = find_meta_value(meta, [r'value *unit', r'unit *value', r'phase *unit'])

    dx = safe_float(pxw_val)
    dy = safe_float(pxh_val)
    px_unit = str(pxunit_val).strip() if pxunit_val else None

    phase_min = safe_float(vmin_val)
    phase_max = safe_float(vmax_val)
    value_unit = str(vunit_val).strip() if vunit_val else None

    # Chargement image
    arr = np.array(im)
    if arr.dtype != np.uint16:
        arr = arr.astype(np.uint16)
    arr = arr.astype(np.float64)

    # Mapping en phase
    if (phase_min is not None) and (phase_max is not None):
        phase = phase_min + (arr / 65535.0) * (phase_max - phase_min)
        note = f"Mapping linéaire 0..65535 -> [{phase_min}, {phase_max}] {value_unit or ''}"
    else:
        phase = (arr / 65535.0) * (2*np.pi) - np.pi
        note = "Aucune Value Min/Max : mapping 0..65535 -> [-π, π]"

    # Affichage résumé
    print("→ Métadonnées utilisées :")
    if dx is not None and dy is not None:
        print(f"   Pixel Width  (dx) = {dx} {px_unit or ''}")
        print(f"   Pixel Height (dy) = {dy} {px_unit or ''}")
    else:
        print("   Pixel Width/Height non trouvés")

    if phase_min is not None and phase_max is not None:
        print(f"   Value Min = {phase_min} {value_unit or ''}")
        print(f"   Value Max = {phase_max} {value_unit or ''}")
    else:
        print("   Value Min/Max non trouvés")

    print(f"   Note : {note}\n")

    return phase, dx, dy

def apply_window(phase, window_type='hann', width_fraction=0.1):
    """Applique une fenêtre sur les bords de l'image"""
    Ny, Nx = phase.shape
    
    # Largeur de la zone de transition
    wx = int(Nx * width_fraction)
    wy = int(Ny * width_fraction)
    
    # Création des fenêtres 1D
    win_x = np.ones(Nx)
    win_y = np.ones(Ny)
    
    if window_type == 'hann':
        win_x[:wx] = np.sin(np.linspace(0, np.pi/2, wx))**2
        win_x[-wx:] = np.sin(np.linspace(np.pi/2, 0, wx))**2
        win_y[:wy] = np.sin(np.linspace(0, np.pi/2, wy))**2
        win_y[-wy:] = np.sin(np.linspace(np.pi/2, 0, wy))**2
    
    # Fenêtre 2D
    window = np.outer(win_y, win_x)
    
    return phase * window

# Zero-padding
# Augmenter la taille du tableau en ajoutant des zéros autour 
def simulate_fresnel_padded(phase, dx, dy, E0_eV=200e3, defocus=2e-6, Cs=1.2e-3, pad_factor=2):
    Ny, Nx = phase.shape
    lam = relativistic_wavelength(E0_eV)
    
    # Dimensions avec padding
    Ny_pad = Ny * pad_factor
    Nx_pad = Nx * pad_factor
    
    # Créer psi0 avec padding
    psi0 = np.exp(1j * phase)
    psi0_padded = np.zeros((Ny_pad, Nx_pad), dtype=complex)
    
    # Centrer l'image originale
    y_start = (Ny_pad - Ny) // 2
    x_start = (Nx_pad - Nx) // 2
    psi0_padded[y_start:y_start+Ny, x_start:x_start+Nx] = psi0
    
    # FFT sur le tableau paddé
    kx = np.fft.fftfreq(Nx_pad, d=dx)
    ky = np.fft.fftfreq(Ny_pad, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2
    
    chi = np.pi*lam*defocus*k2 + 0.5*np.pi*Cs*(lam**3)*(k2**2)
    H = np.exp(-1j * chi)
    
    Psi0 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psi0_padded)))
    Psi_def = Psi0 * H
    psi_def = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Psi_def)))
    
    # Extraire la région centrale
    psi_def_crop = psi_def[y_start:y_start+Ny, x_start:x_start+Nx]
    
    return np.abs(psi_def_crop)**2

# Soustraction de tendance (Detrending)
# Retirer la pente moyenne de la phase avant FFT 
def detrend_phase(phase):
    """Retire les tendances linéaires de la phase"""
    Ny, Nx = phase.shape
    y, x = np.mgrid[0:Ny, 0:Nx]
    
    # Ajustement d'un plan
    A = np.c_[x.ravel(), y.ravel(), np.ones(Nx*Ny)]
    coef, _, _, _ = np.linalg.lstsq(A, phase.ravel(), rcond=None)
    
    plane = coef[0]*x + coef[1]*y + coef[2]
    return phase - plane

# Approche combinée (recommandée)
def simulate_fresnel_improved(phase, dx, dy, E0_eV=200e3, defocus=2e-6, Cs=1.2e-3, 
                              use_window=True, pad_factor=1.5):
    Ny, Nx = phase.shape
    lam = relativistic_wavelength(E0_eV)
    
    # 1. Fenêtrage optionnel
    if use_window:
        phase_proc = apply_window(phase, width_fraction=0.05)
    else:
        phase_proc = phase
    
    psi0 = np.exp(1j * phase_proc)
    
    # 2. Zero-padding
    Ny_pad = int(Ny * pad_factor)
    Nx_pad = int(Nx * pad_factor)
    psi0_padded = np.zeros((Ny_pad, Nx_pad), dtype=complex)
    
    y_start = (Ny_pad - Ny) // 2
    x_start = (Nx_pad - Nx) // 2
    psi0_padded[y_start:y_start+Ny, x_start:x_start+Nx] = psi0
    
    # 3. FFT et propagation
    kx = np.fft.fftfreq(Nx_pad, d=dx)
    ky = np.fft.fftfreq(Ny_pad, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2
    
    chi = np.pi*lam*defocus*k2 + 0.5*np.pi*Cs*(lam**3)*(k2**2)
    H = np.exp(-1j * chi)
    
    Psi0 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psi0_padded)))
    Psi_def = Psi0 * H
    psi_def = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Psi_def)))
    
    # 4. Extraction de la région d'intérêt
    psi_def_crop = psi_def[y_start:y_start+Ny, x_start:x_start+Nx]
    
    return np.abs(psi_def_crop)**2

def simulate_fresnel(phase, dx, dy, E0_eV=200e3, defocus=2e-6, Cs=1.2e-3):
    Ny, Nx = phase.shape
    lam = relativistic_wavelength(E0_eV)

    psi0 = np.exp(1j * phase)
    kx = np.fft.fftfreq(Nx, d=dx)
    ky = np.fft.fftfreq(Ny, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2

    chi = np.pi*lam*defocus*k2 + 0.5*np.pi*Cs*(lam**3)*(k2**2)
    H = np.exp(-1j * chi)

    Psi0 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psi0)))
    Psi_def = Psi0 * H
    psi_def = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Psi_def)))

    return np.abs(psi_def)**2

def make_periodic(phase):
    """Force la phase à être périodique en soustrayant une surface bilinéaire"""
    Ny, Nx = phase.shape
    
    # Valeurs moyennes sur chaque bord
    top = np.mean(phase[0, :])
    bottom = np.mean(phase[-1, :])
    left = np.mean(phase[:, 0])
    right = np.mean(phase[:, -1])
    
    # Créer une grille de correction
    y, x = np.mgrid[0:Ny, 0:Nx]
    
    # Interpolation bilinéaire des bords
    correction = (top * (Ny-1-y) + bottom * y) / (Ny-1)
    correction += (left * (Nx-1-x) + right * x) / (Nx-1)
    correction -= (top + bottom + left + right) / 4  # normalisation
    
    return phase - correction

def cosine_edge_apodization(phase, edge_width=0.15):
    """Applique une apodisation cosinus uniquement sur les bords"""
    Ny, Nx = phase.shape
    
    # Largeur de transition en pixels
    wx = int(Nx * edge_width)
    wy = int(Ny * edge_width)
    
    # Fonction de transition cosinus
    def cosine_window(n, w):
        win = np.ones(n)
        if w > 0:
            # Bord gauche/haut
            win[:w] = 0.5 * (1 - np.cos(np.pi * np.arange(w) / w))
            # Bord droit/bas
            win[-w:] = 0.5 * (1 - np.cos(np.pi * np.arange(w-1, -1, -1) / w))
        return win
    
    win_x = cosine_window(Nx, wx)
    win_y = cosine_window(Ny, wy)
    window = np.outer(win_y, win_x)
    
    # Ramener progressivement vers la moyenne sur les bords
    phase_mean = np.mean(phase)
    return phase_mean + (phase - phase_mean) * window

def simulate_fresnel_mirrored(phase, dx, dy, E0_eV=200e3, defocus=2e-6, Cs=1.2e-3):
    Ny, Nx = phase.shape
    lam = relativistic_wavelength(E0_eV)
    
    # Padding par symétrie miroir
    # Haut-bas
    phase_pad = np.vstack([
        np.flipud(phase),
        phase,
        np.flipud(phase)
    ])
    
    # Gauche-droite
    phase_pad = np.hstack([
        np.fliplr(phase_pad),
        phase_pad,
        np.fliplr(phase_pad)
    ])
    
    Ny_pad, Nx_pad = phase_pad.shape
    psi0 = np.exp(1j * phase_pad)
    
    # FFT sur l'image étendue
    kx = np.fft.fftfreq(Nx_pad, d=dx)
    ky = np.fft.fftfreq(Ny_pad, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2
    
    chi = np.pi*lam*defocus*k2 + 0.5*np.pi*Cs*(lam**3)*(k2**2)
    H = np.exp(-1j * chi)
    
    Psi0 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psi0)))
    Psi_def = Psi0 * H
    psi_def = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Psi_def)))
    
    # Extraire la région centrale originale
    psi_def_crop = psi_def[Ny:2*Ny, Nx:2*Nx]
    
    return np.abs(psi_def_crop)**2

def remove_edge_discontinuity(phase):
    """Retire les discontinuités de bord par résolution de Poisson"""
    from scipy.ndimage import laplace
    from scipy.sparse import diags, csr_matrix
    from scipy.sparse.linalg import spsolve
    
    Ny, Nx = phase.shape
    
    # Calculer le Laplacien de la phase
    lap = laplace(phase)
    
    # Résoudre l'équation de Poisson avec conditions aux limites périodiques
    # Pour simplification, on utilise une approximation
    phase_corrected = np.copy(phase)
    
    # Alternative simple : soustraire un plan qui minimise les sauts aux bords
    edges = np.concatenate([
        phase[0, :],   # haut
        phase[-1, :],  # bas
        phase[:, 0],   # gauche
        phase[:, -1]   # droite
    ])
    edge_mean = np.mean(edges)
    
    phase_corrected -= edge_mean
    
    return phase_corrected

def simulate_fresnel_robust(phase, dx, dy, E0_eV=200e3, defocus=2e-6, Cs=1.2e-3):
    """Version robuste avec gestion des discontinuités de bord + padding par symétrie miroir"""
    Ny, Nx = phase.shape
    lam = relativistic_wavelength(E0_eV)
    
    # Étape 1 : Forcer la périodicité
    phase_periodic = make_periodic(phase)
    
    # Étape 2 : Apodisation douce des bords restants
    phase_apo = cosine_edge_apodization(phase_periodic, edge_width=0.1)
    
    # Étape 3 : Padding par symétrie miroir pour assurer la continuité
    # Haut-bas
    phase_pad = np.vstack([
        np.flipud(phase_apo),
        phase_apo,
        np.flipud(phase_apo)
    ])
    
    # Gauche-droite
    phase_pad = np.hstack([
        np.fliplr(phase_pad),
        phase_pad,
        np.fliplr(phase_pad)
    ])
    
    # Étape 4 : Créer psi0 sur l'image étendue
    Ny_pad, Nx_pad = phase_pad.shape
    psi0 = np.exp(1j * phase_pad)
    
    # Étape 5 : Propagation de Fresnel
    kx = np.fft.fftfreq(Nx_pad, d=dx)
    ky = np.fft.fftfreq(Ny_pad, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2
    
    chi = np.pi*lam*defocus*k2 + 0.5*np.pi*Cs*(lam**3)*(k2**2)
    H = np.exp(-1j * chi)
    
    Psi0 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psi0)))
    Psi_def = Psi0 * H
    psi_def = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Psi_def)))
    
    # Étape 6 : Extraire la région centrale originale
    psi_def_crop = psi_def[Ny:2*Ny, Nx:2*Nx]
    
    return np.abs(psi_def_crop)**2
    
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--png", required=True, help="PNG 16-bit contenant la phase")
    p.add_argument("--defocus", type=float, default=2e-6)
    p.add_argument("--E0", type=float, default=200e3)
    p.add_argument("--Cs", type=float, default=1.2e-3)
    return p.parse_args()

def main():
    args = parse_args()

    # Charger phase + métadonnées utiles
    phase, dx, dy = load_phase_png_with_metadata(args.png)

    # Fallback si dx/dy absents
    if dx is None or dy is None:
        dx = dy = 1e-9
        print(f"⚠ dx/dy manquants → utilisation de dx = dy = {dx} m\n")

    # Affichage des paramètres de simulation
    print("→ Paramètres de simulation :")
    print(f"   Tension E0 = {args.E0} eV")
    print(f"   Défocus    = {args.defocus} m")
    print(f"   Constante Cs = {args.Cs} m\n")

    # Simulation Fresnel
    I = simulate_fresnel_robust(phase, dx, dy, E0_eV=args.E0, defocus=args.defocus, Cs=args.Cs)

    Ny, Nx = I.shape
    extent_nm = (-(Nx*dx)/2*1e9, (Nx*dx)/2*1e9, -(Ny*dy)/2*1e9, (Ny*dy)/2*1e9)

    # Sauvegarde de l’image
    plt.figure(figsize=(7,6))
    plt.imshow(I, cmap="gray", extent=extent_nm)
    plt.xlabel("x (nm)")
    plt.ylabel("y (nm)")
    plt.title("Contraste de Fresnel")
    plt.colorbar(label="Intensité")
    plt.savefig("fresnel.png", dpi=150)
    plt.close()

    print("✔  Image sauvegardée dans : fresnel.png")

if __name__ == "__main__":
    main()

