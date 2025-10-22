#!/usr/bin/env python3

def define_region_functions():
    """Define vector field for each region."""
    return {
        'V1': lambda x, y, z: np.array([0, 0, 1]),
        'V2': lambda x, y, z: np.array([x, y, 0]),
        'V3': lambda x, y, z: np.array([np.sin(x), np.cos(y), z]),
        'V4': lambda x, y, z: np.array([np.sin(x), np.cos(y), z]),
        'V5': lambda x, y, z: np.array([np.sin(x), np.cos(y), z])
    }
    
import os
import sys
import numpy as np
import meshio

from datetime import datetime

import warnings

warnings.simplefilter("always")  # pour forcer l'affichage
warnings.formatwarning = lambda message, category, filename, lineno, line=None: f"{message}\n"

def get_region_elements(mesh):
    """
    Return a dictionary mapping region name -> list of element descriptions.
    Each element description is a dict:
      {
        "cell_type": "triangle" / "tetra" / ...,
        "block_index": i,            # index of the cell block in mesh.cells
        "local_element_index": j,    # element index inside its cell block
        "global_element_index": g,   # element index counting across all blocks
        "nodes": array([...], dtype=int)  # node indices for that element
      }

    Requires mesh.field_data to map names -> [physical_id, dim].
    Uses mesh.cell_data to find the physical tag per element (common key is "gmsh:physical").
    """
    if not hasattr(mesh, "cells"):
        raise ValueError("mesh has no attribute 'cells'")

    # 1) build mapping physical_id -> region_name using mesh.field_data
    id_to_name = {}
    for name, arr in mesh.field_data.items():            # <-- mesh.field_data.items() as requested
        try:
            phys_id = int(np.asarray(arr).flatten()[0])
        except Exception:
            continue
        id_to_name[phys_id] = name

    # 2) find the appropriate key in mesh.cell_data that holds physical tags
    phys_key = None
    if hasattr(mesh, "cell_data") and mesh.cell_data:
        # mesh.cell_data is a dict: key -> list_of_arrays (one per cell block)
        # prefer the exact "gmsh:physical" key, else any key containing "physical"
        for k in mesh.cell_data.keys():
            if k == "gmsh:physical":
                phys_key = k
                break
        if phys_key is None:
            for k in mesh.cell_data.keys():
                if "physical" in k.lower():
                    phys_key = k
                    break

    if phys_key is None:
        raise ValueError("No physical tags found in mesh.cell_data. Look for key 'gmsh:physical' or similar.")

    # 3) iterate cell blocks and attach elements to regions
    regions = {}     # region_name -> list of element dicts
    global_idx = 0
    # mesh.cell_data[phys_key] should be a list/sequence with same length as mesh.cells
    phys_tags_list = mesh.cell_data[phys_key]
    if len(phys_tags_list) != len(mesh.cells):
        # sometimes meshio stores cell_data as dict of lists; if shapes mismatch, try to handle gracefully
        raise ValueError("Mismatch between number of cell blocks and number of cell_data arrays for key "
                         f"'{phys_key}' (found {len(phys_tags_list)} vs {len(mesh.cells)}).")

    for block_index, (cell_block, phys_tags) in enumerate(zip(mesh.cells, phys_tags_list)):
        cell_type = cell_block.type
        connectivity = np.asarray(cell_block.data, dtype=int)  # shape (n_elems_block, n_nodes_per_elem)
        phys_tags = np.asarray(phys_tags, dtype=int).flatten() # length n_elems_block

        if connectivity.shape[0] != phys_tags.shape[0]:
            raise ValueError(f"Number of elements in cell block {block_index} ({connectivity.shape[0]}) "
                             f"does not match number of physical tags ({phys_tags.shape[0]}).")

        for local_idx in range(connectivity.shape[0]):
            phys_id = int(phys_tags[local_idx])
            region_name = id_to_name.get(phys_id, f"phys_{phys_id}")  # fallback name if not found in field_data
            elem = {
                "cell_type": cell_type,
                "block_index": block_index,
                "local_element_index": local_idx,
                "global_element_index": global_idx,
                "nodes": connectivity[local_idx].tolist()
            }
            regions.setdefault(region_name, []).append(elem)
            global_idx += 1

    return regions

def check_gmsh_version(filename):
    """
    Check the Gmsh .msh file version before reading it.
    Raises a ValueError if the version is below 4.0.
    """
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() == "$MeshFormat":
                version_line = next(f).strip()
                try:
                    version = float(version_line.split()[0])
                    if version < 4.0:
                        raise ValueError(f"Gmsh format version {version} is not supported. Use version 4.0 or higher.")
                    print(f"Gmsh file format version: {version}")
                except Exception:
                    raise ValueError("Unable to read Gmsh file version.")
                break
    return True


def load_mesh(filename):
    """
    Load a mesh from a file using meshio.
    Returns the point coordinates and the triangular connectivity.
    Raises an error if the file does not exist or contains no triangles.
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"The file '{filename}' does not exist.")

    try:
        # Check Gmsh version > 4.0
        # check_gmsh_version(filename)

        # Load mesh using meshio
        mesh = meshio.read(filename)
        # For each cell type (triangle, tetra, etc.)
        for cell_block, cell_data in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
            cell_type = cell_block.type
            print(f"\nCell type: {cell_type}")
            for region_name, (region_id, dim) in mesh.field_data.items():
                element_indices = np.where(cell_data == region_id)[0]
                if len(element_indices) > 0:
                    print(f"  Region '{region_name}' (ID={region_id}, dim={dim}): {len(element_indices)} elements")

        regions = get_region_elements(mesh)

        # afficher le nombre d'éléments par région
        # for region_name, elems in regions.items():
        #     print(f"Region '{region_name}': {len(elems)} elements")
        #     # afficher le premier élément (node indices)
        #     print("  first element nodes:", elems[0]["nodes"])

        # Display all element node indices for each region
        # for region_name, elems in regions.items():
        #     print(f"\nRegion '{region_name}' — {len(elems)} elements:")
        #     for elem in elems:
        #         print(f"  Element {elem['global_element_index']:4d} ({elem['cell_type']}): nodes {elem['nodes']}")
 
        points = mesh.points[:, :3]

        triangles = []
        tetrahedrons = []
        unsupported_types = []

        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                triangles.append(cell_block.data)
            elif cell_block.type == "tetra":
                tetrahedrons.append(cell_block.data)
            else:
                unsupported_types.append(cell_block.type)

        # Raise an error if any unsupported element types are found
        if unsupported_types:
            unique_types = ", ".join(sorted(set(unsupported_types)))
            raise ValueError(
                f"Unsupported element types found in mesh: {unique_types}. "
                f"Only 'triangle' and 'tetra' elements are supported."
            )

        if not triangles:
            raise ValueError(f"No triangular cells found in mesh '{filename}'.")

        if not tetrahedrons:
            raise ValueError(f"No tetrahedron cells found in mesh '{filename}'.")

        # Concatenate all triangle blocks into a single array
        triangles = np.vstack(triangles)

        # Concatenate all tetrahedron blocks into a single array
        tetrahedrons = np.vstack(tetrahedrons)

        return points, triangles, tetrahedrons, regions

    except Exception as e:
        raise RuntimeError(f"Error while loading mesh '{filename}': {e}")


def apply_region_functions(points, regions, region_functions):
    """Apply region defined functions."""
    m = np.full((len(points), 3), np.nan, dtype=np.float64)
    for region_name, elems in regions.items():
        if region_name not in region_functions:
            warnings.warn(f"⚠️ No function defined for region '{region_name}'. Region ignored.")
            continue

        func = region_functions[region_name]
        for elem in elems:
            for node_id in elem['nodes']:
                idx = node_id - 1
                xyz = points[idx, :]
                m[idx] = func(*xyz)

    if np.isnan(m).any():
        warnings.warn("m contains non-initialized values")
        sys.exit("❌ abort")

    return m


def write_solution_file(m, filename="sol.in"):
    """Écrit le champ m et les métadonnées dans le fichier sol.in."""
    now = datetime.now().isoformat()
    header_lines = [
        "## feeLLGood version: 0.10.3-24-g3103a8d",
        "## hostname: jc-Precision-5570",
        f"## real-world time: {now}",
        "## settings file: settings.yml",
        "## time: 0.000000e+00",
        "## columns: idx\tmx\tmy\tmz\tphi"
    ]

    idx = np.arange(m.shape[0])[:, np.newaxis]
    phi = np.zeros((m.shape[0], 1))
    tab = np.hstack((idx, m, phi))

    with open(filename, "w") as f:
        for line in header_lines:
            f.write(line + "\n")
        np.savetxt(f, tab, fmt=["%d", "%.7e", "%.7e", "%.7e", "%.7e"], delimiter="\t")

    print(f"✅ File '{filename}' successfully written ({m.shape[0]} lines).")

def main(mesh_file):
    """Routine principale."""
    points, triangles, tetrahedrons, regions = load_mesh(mesh_file)
    region_functions = define_region_functions()
    m = apply_region_functions(points, regions, region_functions)
    write_solution_file(m)
    
if __name__ == "__main__":
    # Check for required argument
    if len(sys.argv) < 2:
        print("Usage: python fg-gensol.py mesh_file.msh")
        sys.exit(1)

    mesh_file = sys.argv[1]
    main(mesh_file)
