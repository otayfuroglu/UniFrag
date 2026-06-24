#!/usr/bin/env python3
"""
prepare_linker4qm.py

Post-processing script that converts a folder of linker structure files
(.xyz or .cif, produced by UniFrag/moffragmentor helper fragmentation) into
a single QM-ready ExtXYZ collection file.

Each linker is processed through the SAME QM-preparation pipeline used inside
fragmentation_oop.py:

  Step 1  Read structure (.xyz or .cif)
  Step 2  Strip metal atoms  (Zn, Cu, Mg, Fe, Co, Ni, Mn, Zr, Ti, V, Cr, Al, …)
  Step 3  Cap open terminal atoms with H
            - O with one C/P/S neighbor and no existing H  → O-H cap
            - C with one heavy neighbor (e.g. broken linker stub)  → C-H cap
            - N with one heavy neighbor                             → N-H cap
          Direction: away from the bonding neighbor (same as _cap_path_j_open_oxygens)
  Step 4  Optimise capped-H geometry only (sp2 C-H planarity, O-H tetrahedral angle)
  Step 5  Fix odd electron count (remove lowest-priority H if needed)
  Step 6  Deduplicate by heavy-atom formula (only unique chemistries written)
  Step 7  Write ExtXYZ frame with label=<REFCODE>LinkerMof

Usage:
    python runUniFrag/prepare_linker4qm.py <linker_folder> [output_extxyz]

Arguments:
    linker_folder   Path to folder containing .xyz and/or .cif linker files
    output_extxyz   (Optional) Output ExtXYZ filepath.
                    Defaults to <linker_folder>/linkers_collection.extxyz

Notes:
    - Both .xyz and .cif files in the same folder are processed together.
    - CIF files are read with ASE (fractional → Cartesian conversion automatic).
    - Processing order: all .xyz files first (alphabetical), then .cif files.
    - The capping/QM logic is borrowed directly from fragmentation_oop.py and
      is self-contained here so no import from that module is needed.
"""

import os
import sys
import re
import glob
import numpy as np
from collections import Counter, deque

# ---------------------------------------------------------------------------
# Constants (identical to BaseFragmenter in fragmentation_oop.py)
# ---------------------------------------------------------------------------

METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi",
}

LARGE_NON_METALS = {"Br", "I", "S", "P", "Cl"}

_ATOMIC_NUMBERS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33,
    "Se": 34, "Br": 35, "Kr": 36, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52,
    "I": 53, "Xe": 54, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62,
    "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}

# Bond-length defaults for H caps (identical to cap_bond_length in fragmentation_oop.py)
def _cap_bond_length(element):
    if element == "C":  return 1.09
    if element == "N":  return 1.01
    if element == "Si": return 1.48
    if element == "B":  return 1.19
    return 0.96   # O, S, P, and everything else


# ---------------------------------------------------------------------------
# Bond validity (identical to is_valid_bond in fragmentation_oop.py)
# ---------------------------------------------------------------------------

def _is_valid_bond(s1, s2, dist):
    if s1 in METALS or s2 in METALS:
        return dist < 2.6
    if "H" in (s1, s2):
        return dist < 1.2
    if s1 in LARGE_NON_METALS or s2 in LARGE_NON_METALS:
        return dist < 2.2
    return dist < 1.8


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _read_xyz(filepath):
    """
    Read a plain .xyz file.
    Returns (species: list[str], coords: list[[x,y,z]]).
    """
    with open(filepath, "r") as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError("XYZ file has fewer than 3 lines")
    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(f"Cannot parse atom count: {lines[0]!r}")
    species, coords = [], []
    for i in range(2, 2 + n_atoms):
        if i >= len(lines):
            break
        parts = lines[i].strip().split()
        if len(parts) >= 4:
            species.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if len(species) != n_atoms:
        raise ValueError(f"Expected {n_atoms} atoms but parsed {len(species)}")
    return species, coords


def _read_cif(filepath):
    """
    Read a CIF file using ASE (handles fractional → Cartesian automatically).
    Returns (species: list[str], coords: list[[x,y,z]]).
    """
    try:
        import ase.io
    except ImportError:
        raise ValueError("ASE is required to read CIF files. Install with: pip install ase")
    try:
        atoms = ase.io.read(filepath, index=0)
    except Exception as exc:
        raise ValueError(f"ASE could not read CIF: {exc}")
    species = list(atoms.get_chemical_symbols())
    coords  = [list(pos) for pos in atoms.get_positions()]
    if not species:
        raise ValueError("CIF file contains no atoms")
    return species, coords


def _read_structure(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".cif":
        return _read_cif(filepath)
    elif ext == ".xyz":
        return _read_xyz(filepath)
    else:
        raise ValueError(f"Unsupported file extension: {ext!r}")


def _write_extxyz_frame(fh, species, coords, label):
    """Append one ExtXYZ frame — format matches the UniFrag pipeline."""
    fh.write(f"{len(species)}\n")
    fh.write(f'Properties=species:S:1:pos:R:3 label={label} pbc="F F F"\n')
    for sym, pos in zip(species, coords):
        fh.write(f"{sym:<8}{pos[0]:16.8f}{pos[1]:16.8f}{pos[2]:16.8f}\n")


# ---------------------------------------------------------------------------
# Step 2: strip metal atoms
# ---------------------------------------------------------------------------

def _strip_metals(species, coords):
    """Remove all metal-element atoms and return (species, coords, n_stripped)."""
    keep = [i for i, s in enumerate(species) if s not in METALS]
    n_stripped = len(species) - len(keep)
    return [species[i] for i in keep], [coords[i] for i in keep], n_stripped


# ---------------------------------------------------------------------------
# Step 3: cap open terminal atoms
# (mirrors _cap_path_j_open_oxygens + stub-C/N capping from fragmentation_oop.py)
# ---------------------------------------------------------------------------

def _orthonormal_basis(u):
    trial = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(u, trial)) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, trial)
    n1 = np.linalg.norm(e1)
    if n1 < 1e-12:
        e1 = np.array([0.0, 0.0, 1.0])
        n1 = np.linalg.norm(e1)
    e1 /= n1
    e2 = np.cross(u, e1)
    e2 /= np.linalg.norm(e2)
    return e1, e2


def _score_h_candidate(cand, species, coords, parent_idx):
    min_h, min_heavy = float("inf"), float("inf")
    for i, s in enumerate(species):
        d = np.linalg.norm(cand - np.array(coords[i], dtype=float))
        if s == "H":
            if d < min_h: min_h = d
        else:
            if i != parent_idx and d < min_heavy: min_heavy = d
    return min_h, min_heavy


def _place_capping_h(parent_idx, base_vec, bl, species, coords,
                     min_hh=1.5, min_heavy=0.9, min_o_contact=1.5,
                     capped_h_flags=None):
    """Place one capping H on parent_idx along base_vec direction at distance bl."""
    parent_pos = np.array(coords[parent_idx], dtype=float)
    vnorm = np.linalg.norm(base_vec)
    if vnorm < 1e-12:
        return False
    u = base_vec / vnorm
    e1, e2 = _orthonormal_basis(u)

    theta_list = [0.0, 20.0, 35.0, 50.0, 65.0, 80.0, 110.0, 140.0, 170.0]
    phi_list   = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
    directions = []
    for th in theta_list:
        th_r = np.deg2rad(th)
        ct, st = np.cos(th_r), np.sin(th_r)
        if th == 0.0:
            directions.append(u)
        else:
            for ph in phi_list:
                ph_r = np.deg2rad(ph)
                dv = ct * u + st * (np.cos(ph_r) * e1 + np.sin(ph_r) * e2)
                directions.append(dv / np.linalg.norm(dv))

    parent_sp = species[parent_idx]
    for dvec in directions:
        cand = parent_pos + dvec * bl
        too_close_o = False
        for oi, os_ in enumerate(species):
            if os_ != "O":
                continue
            if parent_sp == "O" and oi == parent_idx:
                continue
            if np.linalg.norm(cand - np.array(coords[oi], dtype=float)) < min_o_contact:
                too_close_o = True
                break
        if too_close_o:
            continue
        mh, mheavy = _score_h_candidate(cand, species, coords, parent_idx)
        if mh >= min_hh and mheavy >= min_heavy:
            species.append("H")
            coords.append(list(cand))
            if capped_h_flags is not None:
                capped_h_flags.append(True)
            return True
    return False


def _oxygen_already_protonated(i, species, coords, oh_cutoff=1.5):
    if species[i] != "O":
        return False
    opos = np.array(coords[i], dtype=float)
    for j, s in enumerate(species):
        if s != "H":
            continue
        if np.linalg.norm(opos - np.array(coords[j], dtype=float)) <= oh_cutoff:
            return True
    return False


def _cap_open_atoms(species, coords, capped_h_flags):
    """
    Cap open terminal atoms that lost their metal neighbour during metal stripping.

    Rules (same as fragmentation_oop.py Path J capping):
      O with exactly one C/P/S neighbour and no H  → O-H (away from C)
      C with exactly one heavy neighbour            → C-H (away from neighbour)
      N with exactly one or two heavy neighbours    → N-H (away from neighbours avg)
    """
    n_capped = 0

    # Track C/P/S centres that already have one O-H cap, to avoid double-capping
    # the same carboxylate carbon (mirrors capped_central_atoms logic)
    capped_central = set()
    for i, sp in enumerate(species):
        if sp != "O":
            continue
        if not _oxygen_already_protonated(i, species, coords):
            continue
        opos = np.array(coords[i], dtype=float)
        for j, spj in enumerate(species):
            if j == i or spj == "H":
                continue
            d = np.linalg.norm(opos - np.array(coords[j], dtype=float))
            if _is_valid_bond("O", spj, d) and spj in ("C", "P", "S"):
                capped_central.add(j)

    heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]

    for i in list(heavy_idx):
        if i >= len(species):
            continue
        sp = species[i]

        # ── Oxygen capping ──────────────────────────────────────────────
        if sp == "O":
            if _oxygen_already_protonated(i, species, coords):
                continue
            opos = np.array(coords[i], dtype=float)
            heavy_nbs = []
            has_metal = False
            for j, spj in enumerate(species):
                if j == i or spj == "H":
                    continue
                d = np.linalg.norm(opos - np.array(coords[j], dtype=float))
                if _is_valid_bond("O", spj, d):
                    heavy_nbs.append(j)
                    if spj in METALS:
                        has_metal = True
            # Only cap terminal oxygens still bonded to exactly one C/P/S
            if has_metal or len(heavy_nbs) != 1 or species[heavy_nbs[0]] not in ("C", "P", "S"):
                continue
            central = heavy_nbs[0]
            if central in capped_central:
                continue
            base = opos - np.array(coords[central], dtype=float)
            bl = _cap_bond_length("O")
            before = len(species)
            _place_capping_h(i, base, bl, species, coords, min_hh=1.5,
                             capped_h_flags=capped_h_flags)
            if len(species) == before and np.linalg.norm(base) > 1e-12:
                _place_capping_h(i, -base, bl, species, coords, min_hh=1.5,
                                 capped_h_flags=capped_h_flags)
            if len(species) > before:
                n_capped += 1
                capped_central.add(central)

        # ── Carbon capping (stub from broken bond) ───────────────────────
        elif sp == "C":
            cpos = np.array(coords[i], dtype=float)
            heavy_nbs = []
            for j, spj in enumerate(species):
                if j == i or spj == "H":
                    continue
                d = np.linalg.norm(cpos - np.array(coords[j], dtype=float))
                if _is_valid_bond("C", spj, d):
                    heavy_nbs.append(j)
            # Only cap if exactly one heavy neighbour remains (broken stub end)
            if len(heavy_nbs) != 1:
                continue
            base = cpos - np.array(coords[heavy_nbs[0]], dtype=float)
            bl = _cap_bond_length("C")
            before = len(species)
            _place_capping_h(i, base, bl, species, coords, min_hh=1.5,
                             capped_h_flags=capped_h_flags)
            if len(species) > before:
                n_capped += 1

        # ── Nitrogen capping ─────────────────────────────────────────────
        elif sp == "N":
            npos = np.array(coords[i], dtype=float)
            heavy_nbs = []
            h_nbs = 0
            for j, spj in enumerate(species):
                if j == i:
                    continue
                d = np.linalg.norm(npos - np.array(coords[j], dtype=float))
                if spj == "H" and d < 1.2:
                    h_nbs += 1
                elif spj != "H" and _is_valid_bond("N", spj, d):
                    heavy_nbs.append(j)
            # Cap N with 1–2 heavy neighbours and no H yet
            if h_nbs > 0 or len(heavy_nbs) == 0 or len(heavy_nbs) > 2:
                continue
            if len(heavy_nbs) == 1:
                base = npos - np.array(coords[heavy_nbs[0]], dtype=float)
            else:
                vecs = [npos - np.array(coords[j], dtype=float) for j in heavy_nbs]
                base = np.mean(vecs, axis=0)
            bl = _cap_bond_length("N")
            before = len(species)
            _place_capping_h(i, base, bl, species, coords, min_hh=1.5,
                             capped_h_flags=capped_h_flags)
            if len(species) > before:
                n_capped += 1

    return n_capped


# ---------------------------------------------------------------------------
# Step 4: optimise capped H geometry
# (mirrors optimize_capped_h_geometry_only in fragmentation_oop.py)
# ---------------------------------------------------------------------------

def _heavy_neighbors_of(idx, species, coords):
    out = []
    ci = np.array(coords[idx], dtype=float)
    for j, sp in enumerate(species):
        if j == idx or sp == "H":
            continue
        d = np.linalg.norm(ci - np.array(coords[j], dtype=float))
        if _is_valid_bond(species[idx], sp, d):
            out.append((j, d))
    out.sort(key=lambda x: x[1])
    return [j for j, _ in out]


def _is_in_ring(idx, species, coords):
    nbs = lambda i: _heavy_neighbors_of(i, species, coords)
    q = deque([(idx, [idx])])
    while q:
        curr, path = q.popleft()
        if len(path) > 8:
            continue
        for nb in nbs(curr):
            if nb == idx and len(path) > 2:
                return True
            if nb not in path:
                q.append((nb, path + [nb]))
    return False


def _enforce_sp2_capped_h(species, coords, capped_h_indices):
    """Place aromatic/sp2 C-H and N-H caps in correct planar direction."""
    for hidx in capped_h_indices:
        if hidx < 0 or hidx >= len(species) or species[hidx] != "H":
            continue
        hpos = np.array(coords[hidx], dtype=float)
        parent, best = None, 1e9
        for i, sp in enumerate(species):
            if sp == "H":
                continue
            d = np.linalg.norm(hpos - np.array(coords[i], dtype=float))
            if _is_valid_bond("H", sp, d) and d < best:
                best = d
                parent = i
        if parent is None or species[parent] not in ("C", "N"):
            continue
        nbs = _heavy_neighbors_of(parent, species, coords)
        if len(nbs) < 2:
            continue
        h_nbs = sum(
            1 for j, spj in enumerate(species)
            if j != parent and spj == "H"
            and np.linalg.norm(np.array(coords[parent]) - np.array(coords[j])) < 1.2
        )
        if len(nbs) == 2 and h_nbs >= 2:
            continue

        p = np.array(coords[parent], dtype=float)
        u = []
        for nb in nbs[:2]:
            v = np.array(coords[nb], dtype=float) - p
            nv = np.linalg.norm(v)
            if nv > 1e-10:
                u.append(v / nv)
        if len(u) < 2:
            continue
        dvec = -(u[0] + u[1])
        nd = np.linalg.norm(dvec)
        if nd < 1e-8:
            continue
        dvec /= nd

        # Enforce global ring planarity
        if _is_in_ring(parent, species, coords):
            bfs = deque([(parent, 0)])
            seen = {parent}
            ring_atoms = []
            while bfs:
                curr, depth = bfs.popleft()
                ring_atoms.append(curr)
                if depth < 3:
                    for nb in _heavy_neighbors_of(curr, species, coords):
                        if nb not in seen:
                            seen.add(nb)
                            bfs.append((nb, depth + 1))
            if len(ring_atoms) >= 5:
                pts = np.array([coords[a] for a in ring_atoms], dtype=float)
                centroid = np.mean(pts, axis=0)
                _, _, vh = np.linalg.svd(pts - centroid, full_matrices=False)
                normal = vh[-1, :]
                dvec = dvec - np.dot(dvec, normal) * normal
                nd = np.linalg.norm(dvec)
                if nd > 1e-8:
                    dvec /= nd

        coords[hidx] = list(p + _cap_bond_length(species[parent]) * dvec)


def _enforce_capped_oh(species, coords, capped_h_indices):
    """Place O-H cap at correct tetrahedral angle (109.5°)."""
    target_angle = np.deg2rad(109.5)
    oh_len = _cap_bond_length("O")

    for hidx in capped_h_indices:
        if hidx < 0 or hidx >= len(species) or species[hidx] != "H":
            continue
        hpos = np.array(coords[hidx], dtype=float)
        parent_o, best = None, float("inf")
        for i, sp in enumerate(species):
            if sp != "O":
                continue
            d = np.linalg.norm(hpos - np.array(coords[i], dtype=float))
            if d < best and _is_valid_bond("H", "O", d):
                best = d
                parent_o = i
        if parent_o is None:
            continue

        opos = np.array(coords[parent_o], dtype=float)
        heavy_nbs = []
        for j, sp in enumerate(species):
            if j == parent_o or sp == "H":
                continue
            d = np.linalg.norm(opos - np.array(coords[j], dtype=float))
            if _is_valid_bond("O", sp, d):
                heavy_nbs.append((j, d))
        if not heavy_nbs:
            continue
        heavy_nbs.sort(key=lambda x: (species[x[0]] != "C", x[1]))

        anchor = heavy_nbs[0][0]
        anchor_pos = np.array(coords[anchor], dtype=float)
        oa = anchor_pos - opos
        oa_norm = np.linalg.norm(oa)
        if oa_norm < 1e-12:
            continue
        e_a = oa / oa_norm

        plane_ref = None
        if species[anchor] == "C":
            best_ref = float("inf")
            for j, sp in enumerate(species):
                if j in {parent_o, anchor} or sp == "H":
                    continue
                d = np.linalg.norm(anchor_pos - np.array(coords[j], dtype=float))
                if _is_valid_bond("C", sp, d) and d < best_ref:
                    plane_ref = j
                    best_ref = d

        if plane_ref is not None:
            ref_vec = np.array(coords[plane_ref], dtype=float) - anchor_pos
            perp = ref_vec - np.dot(ref_vec, e_a) * e_a
            if np.linalg.norm(perp) > 1e-12:
                e_p = -perp / np.linalg.norm(perp)
            else:
                e_p, _ = _orthonormal_basis(e_a)
        else:
            e_p, _ = _orthonormal_basis(e_a)

        def _score(pos):
            mh, mheavy = float("inf"), float("inf")
            for k, sp in enumerate(species):
                if k in {hidx, parent_o}:
                    continue
                d = np.linalg.norm(pos - np.array(coords[k], dtype=float))
                if sp == "H":
                    mh = min(mh, d)
                else:
                    mheavy = min(mheavy, d)
            return mh, mheavy

        dir1 = np.cos(target_angle) * e_a + np.sin(target_angle) * e_p
        dir1 /= np.linalg.norm(dir1)
        c1 = opos + dir1 * oh_len

        dir2 = np.cos(target_angle) * e_a - np.sin(target_angle) * e_p
        n2 = np.linalg.norm(dir2)
        if n2 > 1e-12:
            c2 = opos + (dir2 / n2) * oh_len
            coords[hidx] = list(c1 if _score(c1) >= _score(c2) else c2)
        else:
            coords[hidx] = list(c1)


def _optimise_capped_h_geometry(species, coords, capped_h_indices):
    """Run the full geometry optimisation on capped H atoms only."""
    cap_idx = [i for i in capped_h_indices if 0 <= i < len(species) and species[i] == "H"]
    if not cap_idx:
        return
    _enforce_sp2_capped_h(species, coords, cap_idx)
    _enforce_capped_oh(species, coords, cap_idx)
    # Second pass: re-enforce sp2 after O-H may have slightly shifted ring
    _enforce_sp2_capped_h(species, coords, cap_idx)


# ---------------------------------------------------------------------------
# Step 5: QM even-electron fix
# (mirrors fix_odd_electron_multiplicity in fragmentation_oop.py)
# ---------------------------------------------------------------------------

def _classify_cap_removal_priority(h_idx, species, coords_arr):
    """Return (priority_score, group_label) for removing h_idx."""
    n = len(species)
    h_pos = coords_arr[h_idx]
    best_d, parent_idx = float("inf"), -1
    for j in range(n):
        if j == h_idx or species[j] == "H":
            continue
        d = float(np.linalg.norm(coords_arr[j] - h_pos))
        if d < best_d:
            best_d, parent_idx = d, j
    if parent_idx < 0 or best_d >= 1.5:
        return 0, "unknown"
    parent_sym = species[parent_idx]
    if parent_sym == "C":
        return 0, "carbon"

    h_neighbors = sum(
        1 for j in range(n)
        if species[j] == "H"
        and float(np.linalg.norm(coords_arr[parent_idx] - coords_arr[j])) < 1.3
    )
    if parent_sym == "O" and h_neighbors >= 3:
        return 100, "overcoordinated-O"

    grandparents = [
        (j, species[j], float(np.linalg.norm(coords_arr[j] - coords_arr[parent_idx])))
        for j in range(n)
        if j not in (h_idx, parent_idx) and species[j] != "H"
        and float(np.linalg.norm(coords_arr[j] - coords_arr[parent_idx])) <= 2.2
    ]

    def _gp_o_count(gp_idx):
        return sum(
            1 for j in range(n)
            if j not in (gp_idx, parent_idx) and species[j] == "O"
            and float(np.linalg.norm(coords_arr[j] - coords_arr[gp_idx])) <= 1.5
        )
    def _gp_heavy_count(gp_idx):
        return sum(
            1 for j in range(n)
            if j not in (gp_idx, parent_idx) and species[j] != "H"
            and float(np.linalg.norm(coords_arr[j] - coords_arr[gp_idx])) <= 1.8
        )

    if parent_sym == "O":
        if not grandparents:
            return 30, "bare-O"
        gp_idx, gp_sym, _ = grandparents[0]
        if gp_sym == "C":
            if _gp_o_count(gp_idx) >= 1: return 80, "carboxylate"
            if _gp_heavy_count(gp_idx) >= 3: return 60, "phenol"
            return 50, "alcohol"
        if gp_sym == "S":
            return (35, "sulfonate") if _gp_o_count(gp_idx) >= 2 else (70, "thiol")
        if gp_sym == "P":
            return 75, "phosphonate"
        return 40, f"O-on-{gp_sym}"

    if parent_sym == "N":
        if not grandparents:
            return 25, "bare-N"
        gp_idx, gp_sym, _ = grandparents[0]
        if gp_sym == "S": return 45, "sulfonamide"
        if gp_sym == "C" and _gp_o_count(gp_idx) >= 1: return 10, "amide"
        return (20, "primary-amine") if len(grandparents) == 1 else (15, "secondary-amine")

    if parent_sym == "B":
        return 5, "boron"
    return 8, f"{parent_sym}-group"


def _fix_odd_electron(species, coords, capped_h_indices, label):
    """
    Ensure even total electron count.
    Prefers removing a capped H; if none available tries any H.
    Returns (species, coords, capped_h_indices, was_fixed).
    """
    z_sum = sum(_ATOMIC_NUMBERS.get(s, 6) for s in species)
    if z_sum % 2 == 0:
        return species, coords, capped_h_indices, False

    n = len(species)
    coords_arr = np.array(coords, dtype=float)

    # Prefer capped H atoms for removal (same priority logic as fragmentation_oop.py)
    candidates = []
    for h_idx in capped_h_indices:
        if h_idx < 0 or h_idx >= n or species[h_idx] != "H":
            continue
        score, group = _classify_cap_removal_priority(h_idx, species, coords_arr)
        if score > 0:
            candidates.append((score, h_idx, group))

    if not candidates:
        # Fall back: any H atom
        for h_idx, sp in enumerate(species):
            if sp != "H":
                continue
            score, group = _classify_cap_removal_priority(h_idx, species, coords_arr)
            candidates.append((score, h_idx, group))

    if not candidates:
        print(f"  QM-Fix WARNING: Odd electron count for '{label}' but no H to remove.")
        return species, coords, capped_h_indices, False

    candidates.sort(key=lambda x: x[0], reverse=True)
    score, target_idx, group = candidates[0]

    parent_sym = "?"
    for j in range(n):
        if j != target_idx and species[j] != "H":
            if float(np.linalg.norm(coords_arr[j] - coords_arr[target_idx])) < 1.5:
                parent_sym = species[j]
                break

    print(f"  QM-Fix [{group}]: Removed capping H[{target_idx}] from {parent_sym} in '{label}'.")
    keep = [i for i in range(n) if i != target_idx]
    new_sp    = [species[i] for i in keep]
    new_co    = [coords[i]  for i in keep]
    new_caps  = [i for i in capped_h_indices if i != target_idx]
    new_caps  = [i - 1 if i > target_idx else i for i in new_caps]
    return new_sp, new_co, new_caps, True


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def _heavy_formula_key(species):
    cnt = Counter(s for s in species if s not in METALS and s != "H")
    return tuple(sorted(cnt.items()))


# ---------------------------------------------------------------------------
# Filename → label
# ---------------------------------------------------------------------------

def _parse_filename(filename):
    base = os.path.splitext(filename)[0]
    m = re.match(r'^(.+?)_(\d+)$', base)
    if m:
        return m.group(1), int(m.group(2))
    return base, 0


def _stem_to_label_base(stem):
    stem = re.sub(r'[\[\](){}]', '', stem)
    parts = re.split(r'[_\-]+', stem)
    result = []
    for p in parts:
        if not p:
            continue
        if p.isdigit():
            result.append(p)
        elif p.isalpha() and p.isupper():
            result.append(p)
        else:
            result.append(p[0].upper() + p[1:])
    return ''.join(result)


def _make_label(stem, extra_index):
    base = _stem_to_label_base(stem)
    if extra_index == 0:
        return f"{base}LinkerMof"
    return f"{base}{extra_index:02d}LinkerMof"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def prepare_linker4qm(linker_folder, output_extxyz=None):
    """Full QM-ready linker preparation pipeline."""

    if output_extxyz is None:
        output_extxyz = os.path.join(linker_folder, "linkers_collection.extxyz")

    xyz_paths = sorted(glob.glob(os.path.join(linker_folder, "*.xyz")))
    cif_paths = sorted(glob.glob(os.path.join(linker_folder, "*.cif")))
    xyz_paths = [p for p in xyz_paths
                 if os.path.abspath(p) != os.path.abspath(output_extxyz)]
    all_paths = xyz_paths + cif_paths

    if not all_paths:
        print(f"No .xyz or .cif files found in: {linker_folder}")
        return

    print("=" * 65)
    print("UniFrag – Linker ExtXYZ QM-Ready Preparation")
    print("=" * 65)
    print(f"Input folder : {linker_folder}")
    print(f"  .xyz files : {len(xyz_paths)}")
    print(f"  .cif files : {len(cif_paths)}")
    print(f"  Total      : {len(all_paths)}")
    print(f"Output file  : {output_extxyz}")
    print()
    print("Pipeline per structure:")
    print("  Read → Strip metals → Cap open atoms → Optimise H → QM-fix → Write")
    print()

    if os.path.exists(output_extxyz):
        os.remove(output_extxyz)

    seen_formula_keys = set()
    stem_extra_index  = {}

    n_written  = 0
    n_dup      = 0
    n_err      = 0
    n_qm_fixed = 0

    csv_rows = []

    with open(output_extxyz, "w") as out_fh:
        for path in all_paths:
            filename = os.path.basename(path)
            ext      = os.path.splitext(filename)[1].lower()
            stem, _  = _parse_filename(filename)

            # 1. Read
            try:
                species, coords = _read_structure(path)
            except Exception as exc:
                print(f"  [ERROR ] {filename}: {exc}")
                n_err += 1
                csv_rows.append({
                    "filename": filename,
                    "label": "",
                    "num_atoms": 0,
                    "num_heavy_atoms": 0,
                    "num_hydrogens": 0,
                    "formula": "",
                    "metals_stripped": 0,
                    "atoms_capped": 0,
                    "qm_fixed": "False",
                    "status": "ERROR",
                    "reason": str(exc)
                })
                continue

            # 2. Strip metals
            species, coords, n_metals = _strip_metals(species, coords)
            if len(species) == 0:
                print(f"  [SKIP  ] {filename}: no non-metal atoms left after stripping")
                n_err += 1
                csv_rows.append({
                    "filename": filename,
                    "label": "",
                    "num_atoms": 0,
                    "num_heavy_atoms": 0,
                    "num_hydrogens": 0,
                    "formula": "",
                    "metals_stripped": n_metals,
                    "atoms_capped": 0,
                    "qm_fixed": "False",
                    "status": "ERROR",
                    "reason": "no non-metal atoms left after stripping"
                })
                continue

            # 3. Cap open terminal atoms — track which H are caps
            capped_h_flags  = [False] * len(species)
            n_capped = _cap_open_atoms(species, coords, capped_h_flags)
            capped_h_indices = [i for i, f in enumerate(capped_h_flags) if f]

            # 4. Optimise capped H geometry
            _optimise_capped_h_geometry(species, coords, capped_h_indices)

            # 5. Deduplication (after capping, before QM-fix)
            fkey = _heavy_formula_key(species)
            if fkey in seen_formula_keys:
                print(f"  [SKIP  ] {filename}  (duplicate heavy-atom formula)")
                n_dup += 1
                formula_str = "".join(f"{el}{cnt}" for el, cnt in sorted(Counter(species).items()))
                num_heavy = sum(1 for s in species if s != "H")
                num_h = sum(1 for s in species if s == "H")
                csv_rows.append({
                    "filename": filename,
                    "label": "",
                    "num_atoms": len(species),
                    "num_heavy_atoms": num_heavy,
                    "num_hydrogens": num_h,
                    "formula": formula_str,
                    "metals_stripped": n_metals,
                    "atoms_capped": n_capped,
                    "qm_fixed": "False",
                    "status": "SKIP_DUP",
                    "reason": "duplicate heavy-atom formula"
                })
                continue
            seen_formula_keys.add(fkey)

            # 6. Build label
            extra_idx = stem_extra_index.get(stem, 0)
            label = _make_label(stem, extra_idx)
            stem_extra_index[stem] = extra_idx + 1

            # 7. QM-fix: even electron count
            species, coords, capped_h_indices, was_fixed = _fix_odd_electron(
                species, coords, capped_h_indices, label
            )
            if was_fixed:
                n_qm_fixed += 1

            # 8. Write
            _write_extxyz_frame(out_fh, species, coords, label)
            n_written += 1

            formula_str = "".join(
                f"{el}{cnt}" for el, cnt in sorted(Counter(species).items())
            )
            metal_note = f" [stripped {n_metals} metals]" if n_metals else ""
            cap_note   = f" [capped {n_capped} atoms]" if n_capped else ""
            tag = ext.lstrip(".")
            print(
                f"  [OK {tag:>3}] {filename:<50}  "
                f"label={label:<40}  {formula_str}{metal_note}{cap_note}"
            )

            num_heavy = sum(1 for s in species if s != "H")
            num_h = sum(1 for s in species if s == "H")
            csv_rows.append({
                "filename": filename,
                "label": label,
                "num_atoms": len(species),
                "num_heavy_atoms": num_heavy,
                "num_hydrogens": num_h,
                "formula": formula_str,
                "metals_stripped": n_metals,
                "atoms_capped": n_capped,
                "qm_fixed": str(was_fixed),
                "status": "OK",
                "reason": ""
            })

    # Write CSV summary
    import csv
    csv_path = os.path.splitext(output_extxyz)[0] + "_summary.csv"
    fieldnames = [
        "filename", "label", "num_atoms", "num_heavy_atoms", "num_hydrogens",
        "formula", "metals_stripped", "atoms_capped", "qm_fixed", "status", "reason"
    ]
    try:
        with open(csv_path, "w", newline="") as csv_fh:
            writer = csv.DictWriter(csv_fh, fieldnames=fieldnames)
            writer.writeheader()
            for r in csv_rows:
                writer.writerow(r)
    except Exception as exc:
        print(f"  [ERROR ] Failed to write CSV summary: {exc}")

    print()
    print("=" * 65)
    print("Summary")
    print("-" * 65)
    print(f"  Input files (.xyz)      : {len(xyz_paths)}")
    print(f"  Input files (.cif)      : {len(cif_paths)}")
    print(f"  Written (unique) frames : {n_written}")
    print(f"  Skipped (duplicates)    : {n_dup}")
    print(f"  Skipped (errors)        : {n_err}")
    print(f"  QM-fixed (H removed)    : {n_qm_fixed}")
    print(f"  Output file             : {output_extxyz}")
    print(f"  CSV summary file        : {csv_path}")
    print("=" * 65)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    folder = sys.argv[1]
    out    = sys.argv[2] if len(sys.argv) > 2 else None
    prepare_linker4qm(folder, out)
