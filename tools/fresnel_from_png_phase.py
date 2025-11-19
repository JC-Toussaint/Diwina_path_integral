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
    I = simulate_fresnel(phase, dx, dy, E0_eV=args.E0, defocus=args.defocus, Cs=args.Cs)

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

