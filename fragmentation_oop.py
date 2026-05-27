import argparse
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pymatgen.core import Molecule, Structure
from pymatgen.io.cif import CifParser


@dataclass
class FragmentResult:
    species: list
    coords: list


def _write_extxyz(filepath, species, coords, label):
    """
    Writes a molecular fragment to an ExtXYZ file.
    Does not require any external library (like ASE) to ensure zero-dependency robustness.
    """
    with open(filepath, "a") as f:
        f.write(f"{len(species)}\n")
        f.write(f"Properties=species:S:1:pos:R:3 label={label} pbc=\"F F F\"\n")
        for sym, pos in zip(species, coords):
            f.write(f"{sym:<8}{pos[0]:16.8f}{pos[1]:16.8f}{pos[2]:16.8f}\n")


def _parse_extxyz(filepath):
    """
    Parses a multi-frame ExtXYZ file and returns a list of dictionaries:
    [{'label': label, 'species': [...], 'coords': [...]}, ...]
    """
    import os
    if not os.path.exists(filepath):
        return []
    frames = []
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n_atoms = int(line)
        except ValueError:
            i += 1
            continue
        
        if i + 1 >= n_lines:
            break
        comment = lines[i + 1].strip()
        
        label = ""
        if "label=" in comment:
            parts = comment.split("label=")
            if len(parts) > 1:
                val = parts[1].split()[0]
                if val.startswith('"') and val.endswith('"'):
                    val = val[1:-1]
                elif val.startswith("'") and val.endswith("'"):
                    val = val[1:-1]
                label = val
        
        species = []
        coords = []
        for j in range(n_atoms):
            idx = i + 2 + j
            if idx >= n_lines:
                break
            atom_line = lines[idx].strip()
            if not atom_line:
                continue
            parts = atom_line.split()
            if len(parts) >= 4:
                species.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        
        if len(species) == n_atoms:
            frames.append({
                "label": label,
                "species": species,
                "coords": coords
            })
        i += 2 + n_atoms
    return frames


def _load_seen_keys_from_extxyz(filepath):
    seen_keys = set()
    frames = _parse_extxyz(filepath)
    for frame in frames:
        key = MOFFragmenter._chemical_identity_key(frame["species"], frame["coords"])
        seen_keys.add(key)
    return seen_keys


def _update_csv_rows(csv_path, file_basename, new_rows_data, header):
    import os
    existing_lines = []
    if os.path.exists(csv_path):
        with open(csv_path, "r") as f:
            existing_lines = f.readlines()
            
    filtered_lines = []
    header_found = False
    for idx, line in enumerate(existing_lines):
        if idx == 0 and ("cif_file" in line or "pdb_file" in line):
            filtered_lines.append(line)
            header_found = True
            continue
        parts = line.split(",")
        if parts and parts[0] == file_basename:
            continue
        filtered_lines.append(line)
        
    if not header_found and header:
        filtered_lines.insert(0, header if header.endswith("\n") else header + "\n")
        
    temp_path = csv_path + ".tmp"
    try:
        with open(temp_path, "w") as f:
            f.writelines(filtered_lines)
            f.writelines(new_rows_data)
        if os.path.exists(csv_path):
            os.remove(csv_path)
        os.rename(temp_path, csv_path)
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise e


def _update_extxyz_collection(extxyz_path, clean_base, new_fragments):
    import os
    existing_frames = []
    if os.path.exists(extxyz_path):
        existing_frames = _parse_extxyz(extxyz_path)
    
    remaining_frames = [f for f in existing_frames if not f["label"].startswith(clean_base)]
    
    temp_path = extxyz_path + ".tmp"
    try:
        if os.path.exists(temp_path):
            os.remove(temp_path)
            
        for f in remaining_frames:
            _write_extxyz(temp_path, f["species"], f["coords"], f["label"])
            
        for species, coords, label in new_fragments:
            _write_extxyz(temp_path, species, coords, label)
            
        if os.path.exists(temp_path):
            if os.path.exists(extxyz_path):
                os.remove(extxyz_path)
            os.rename(temp_path, extxyz_path)
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise e


class BaseFragmenter:
    def __init__(self, radius=6.0):
        self.radius = float(radius)

    @staticmethod
    def cap_bond_length(site_species):
        if site_species == "C":
            return 1.09
        if site_species == "N":
            return 1.01
        if site_species == "Si":
            return 1.48
        if site_species == "B":
            return 1.19
        return 0.96

    @staticmethod
    def _orthonormal_basis(u):
        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(u, trial)) > 0.9:
            trial = np.array([0.0, 1.0, 0.0])
        e1 = np.cross(u, trial)
        n1 = np.linalg.norm(e1)
        if n1 < 1e-12:
            e1 = np.array([0.0, 0.0, 1.0])
            n1 = np.linalg.norm(e1)
        e1 = e1 / n1
        e2 = np.cross(u, e1)
        e2 = e2 / np.linalg.norm(e2)
        return e1, e2

    @staticmethod
    def _score_candidate_h(candidate_pos, species, coords, parent_idx):
        min_h = float("inf")
        min_heavy = float("inf")
        for i, s in enumerate(species):
            d = np.linalg.norm(candidate_pos - np.array(coords[i]))
            if s == "H":
                if d < min_h:
                    min_h = d
            else:
                if i != parent_idx and d < min_heavy:
                    min_heavy = d
        return min_h, min_heavy

    def place_capping_h(self, parent_idx, base_vec, bl, species, coords, min_hh=1.5, min_heavy=0.9, min_o_contact=1.5, capped_h_flags=None):
        parent_pos = np.array(coords[parent_idx])
        vnorm = np.linalg.norm(base_vec)
        if vnorm < 1e-12:
            return
        u = base_vec / vnorm
        e1, e2 = self._orthonormal_basis(u)

        theta_list = [0.0, 20.0, 35.0, 50.0, 65.0, 80.0, 110.0, 140.0, 170.0]
        phi_list = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        directions = []
        for th in theta_list:
            th_r = np.deg2rad(th)
            ct, st = np.cos(th_r), np.sin(th_r)
            if th == 0.0:
                directions.append(u)
            else:
                for ph in phi_list:
                    ph_r = np.deg2rad(ph)
                    dir_vec = ct * u + st * (np.cos(ph_r) * e1 + np.sin(ph_r) * e2)
                    directions.append(dir_vec / np.linalg.norm(dir_vec))

        parent_species = species[parent_idx]
        for dvec in directions:
            cand = parent_pos + dvec * bl
            too_close_to_o = False
            for oi, os in enumerate(species):
                if os != "O":
                    continue
                if parent_species == "O" and oi == parent_idx:
                    continue
                if np.linalg.norm(cand - np.array(coords[oi])) < min_o_contact:
                    too_close_to_o = True
                    break
            if too_close_to_o:
                continue
            mh, mheavy = self._score_candidate_h(cand, species, coords, parent_idx)
            if mh >= min_hh and mheavy >= min_heavy:
                species.append("H")
                coords.append(cand)
                if capped_h_flags is not None:
                    capped_h_flags.append(True)
                return

    @staticmethod
    def oxygen_already_protonated(parent_idx, species, coords, oh_cutoff=1.5):
        if species[parent_idx] != "O":
            return False
        opos = np.array(coords[parent_idx])
        for i, s in enumerate(species):
            if s != "H":
                continue
            if np.linalg.norm(opos - np.array(coords[i])) <= oh_cutoff:
                return True
        return False

    def refine_h_geometry_with_rdkit(self, species, coords, capped_h_indices=None, max_iters=300):
        if not species or "H" not in species:
            return
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit.Geometry import Point3D
        except Exception:
            return

        try:
            rw = Chem.RWMol()
            for sp in species:
                rw.AddAtom(Chem.Atom(sp))

            n = len(species)
            for i in range(n):
                ci = np.array(coords[i], dtype=float)
                for j in range(i + 1, n):
                    d = np.linalg.norm(ci - np.array(coords[j], dtype=float))
                    if self.is_valid_bond(species[i], species[j], d):
                        rw.AddBond(i, j, Chem.BondType.SINGLE)

            mol = rw.GetMol()
            mol.UpdatePropertyCache(strict=False)
            try:
                Chem.SanitizeMol(
                    mol,
                    sanitizeOps=(
                        Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                        | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                        | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                    ),
                )
            except Exception:
                pass

            conf = Chem.Conformer(n)
            for i, c in enumerate(coords):
                x, y, z = [float(v) for v in c]
                conf.SetAtomPosition(i, Point3D(x, y, z))
            mol.AddConformer(conf, assignId=True)

            if not AllChem.UFFHasAllMoleculeParams(mol):
                return
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
            if ff is None:
                return

            if capped_h_indices is None:
                movable_h = {i for i, sp in enumerate(species) if sp == "H"}
            else:
                movable_h = {i for i in capped_h_indices if 0 <= i < len(species) and species[i] == "H"}

            if not movable_h:
                return

            for i in range(n):
                if i not in movable_h:
                    ff.AddFixedPoint(i)

            ff.Initialize()
            ff.Minimize(maxIts=max_iters)

            conf2 = mol.GetConformer(0)
            for i in range(n):
                p = conf2.GetAtomPosition(i)
                coords[i] = np.array([p.x, p.y, p.z], dtype=float)
        except Exception:
            return

    def enforce_sp2_capped_h_geometry(self, species, coords, capped_h_indices=None):
        if not species or not capped_h_indices:
            return

        def heavy_neighbors(idx):
            out = []
            ci = np.array(coords[idx], dtype=float)
            for j, sp in enumerate(species):
                if j == idx or sp == "H":
                    continue
                d = np.linalg.norm(ci - np.array(coords[j], dtype=float))
                if self.is_valid_bond(species[idx], sp, d):
                    out.append((j, d))
            out.sort(key=lambda x: x[1])
            return [j for j, _ in out]

        for hidx in capped_h_indices:
            if hidx < 0 or hidx >= len(species) or species[hidx] != "H":
                continue

            hpos = np.array(coords[hidx], dtype=float)
            parent = None
            best = 1e9
            for i, sp in enumerate(species):
                if sp == "H":
                    continue
                d = np.linalg.norm(hpos - np.array(coords[i], dtype=float))
                if self.is_valid_bond("H", sp, d) and d < best:
                    best = d
                    parent = i
            if parent is None or species[parent] not in ("C", "N"):
                continue

            nbs = heavy_neighbors(parent)
            if len(nbs) < 2:
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

            # sp2 direction at aromatic/phenyl carbon: opposite to sum of two sigma bonds
            dvec = -(u[0] + u[1])
            nd = np.linalg.norm(dvec)
            if nd < 1e-8:
                continue
            dvec = dvec / nd

            # --- GLOBAL AROMATIC PLANARITY ENFORCEMENT ---
            # To ensure the H is perfectly planar with the entire ring (e.g. phenyl),
            # we find all heavy atoms within 3 bonds to fit a best plane.
            def is_in_ring(idx):
                from collections import deque
                q = deque([(idx, [idx])])
                while q:
                    curr, path = q.popleft()
                    if len(path) > 8:
                        continue
                    for nb in heavy_neighbors(curr):
                        if len(path) > 2 and nb == idx:
                            return True
                        if nb not in path:
                            q.append((nb, path + [nb]))
                return False

            if is_in_ring(parent):
                from collections import deque
                q = deque([(parent, 0)])
                seen_heavy = {parent}
                ring_atoms = []
                while q:
                    curr, depth = q.popleft()
                    ring_atoms.append(curr)
                    if depth < 3:
                        for nb in heavy_neighbors(curr):
                            if nb not in seen_heavy:
                                seen_heavy.add(nb)
                                q.append((nb, depth + 1))
                
                if len(ring_atoms) >= 5:
                    pts = np.array([coords[idx] for idx in ring_atoms], dtype=float)
                    centroid = np.mean(pts, axis=0)
                    pts_centered = pts - centroid
                    _, _, vh = np.linalg.svd(pts_centered, full_matrices=False)
                    normal = vh[-1, :]
                    dvec = dvec - np.dot(dvec, normal) * normal
                    nd = np.linalg.norm(dvec)
                    if nd > 1e-8:
                        dvec = dvec / nd
            # ---------------------------------------------

            bl = self.cap_bond_length(species[parent])
            coords[hidx] = p + bl * dvec


    def enforce_capped_oh_geometry(self, species, coords, capped_h_indices=None):
        if not species or not capped_h_indices:
            return
        target_angle = np.deg2rad(109.5)
        oh_len = self.cap_bond_length("O")

        for hidx in capped_h_indices:
            if hidx < 0 or hidx >= len(species) or species[hidx] != "H":
                continue

            hpos = np.array(coords[hidx], dtype=float)
            parent_o = None
            best = float("inf")
            for i, sp in enumerate(species):
                if sp != "O":
                    continue
                d = np.linalg.norm(hpos - np.array(coords[i], dtype=float))
                if d < best and self.is_valid_bond("H", "O", d):
                    best = d
                    parent_o = i
            if parent_o is None:
                continue

            opos = np.array(coords[parent_o], dtype=float)
            heavy_neighbors = []
            for j, sp in enumerate(species):
                if j == parent_o or sp == "H":
                    continue
                d = np.linalg.norm(opos - np.array(coords[j], dtype=float))
                if self.is_valid_bond("O", sp, d):
                    heavy_neighbors.append((j, d))
            if not heavy_neighbors:
                continue
            heavy_neighbors.sort(key=lambda x: (species[x[0]] != "C", x[1]))

            anchor = heavy_neighbors[0][0]
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
                    if self.is_valid_bond("C", sp, d) and d < best_ref:
                        plane_ref = j
                        best_ref = d

            if plane_ref is not None:
                ref_vec = np.array(coords[plane_ref], dtype=float) - anchor_pos
                perp = ref_vec - np.dot(ref_vec, e_a) * e_a
                if np.linalg.norm(perp) > 1e-12:
                    e_p = -perp / np.linalg.norm(perp)
                else:
                    e_p, _ = self._orthonormal_basis(e_a)
            else:
                e_p, _ = self._orthonormal_basis(e_a)

            direction = np.cos(target_angle) * e_a + np.sin(target_angle) * e_p
            norm = np.linalg.norm(direction)
            if norm < 1e-12:
                continue
            candidate = opos + (direction / norm) * oh_len

            alt_direction = np.cos(target_angle) * e_a - np.sin(target_angle) * e_p
            alt_norm = np.linalg.norm(alt_direction)
            if alt_norm > 1e-12:
                alt = opos + (alt_direction / alt_norm) * oh_len

                def score(pos):
                    min_h = float("inf")
                    min_heavy = float("inf")
                    for k, sp in enumerate(species):
                        if k in {hidx, parent_o}:
                            continue
                        d = np.linalg.norm(pos - np.array(coords[k], dtype=float))
                        if sp == "H":
                            min_h = min(min_h, d)
                        else:
                            min_heavy = min(min_heavy, d)
                    return min_h, min_heavy

                if score(alt) > score(candidate):
                    candidate = alt

            coords[hidx] = candidate

    def enforce_single_molecule(self, species, coords, flags=None):
        """Keeps only the largest connected component in the molecule graph,
        removing any disjoint/floating solvent or sub-fragments.
        """
        import numpy as np
        from collections import deque
        n = len(species)
        if n <= 1:
            return (species, coords, flags) if flags is not None else (species, coords)
            
        adj = {i: [] for i in range(n)}
        for i in range(n):
            ci = np.array(coords[i], dtype=float)
            for j in range(i + 1, n):
                d = float(np.linalg.norm(ci - np.array(coords[j], dtype=float)))
                if self.is_valid_bond(species[i], species[j], d):
                    adj[i].append(j)
                    adj[j].append(i)
                    
        visited = set()
        components = []
        for i in range(n):
            if i not in visited:
                comp = []
                q = deque([i])
                visited.add(i)
                while q:
                    curr = q.popleft()
                    comp.append(curr)
                    for nb in adj[curr]:
                        if nb not in visited:
                            visited.add(nb)
                            q.append(nb)
                components.append(comp)
                
        if len(components) <= 1:
            return (species, coords, flags) if flags is not None else (species, coords)
            
        largest = set(max(components, key=len))
        
        new_sp = [s for i, s in enumerate(species) if i in largest]
        new_co = [c for i, c in enumerate(coords) if i in largest]
        
        if flags is not None:
            new_fl = [f for i, f in enumerate(flags) if i in largest]
            return new_sp, new_co, new_fl
        return new_sp, new_co

    def optimize_capped_h_geometry_only(self, species, coords, capped_h_indices=None):
        if not capped_h_indices:
            return
        capped_h_indices = [
            i for i in capped_h_indices
            if 0 <= i < len(species) and species[i] == "H"
        ]
        if not capped_h_indices:
            return
        self.enforce_sp2_capped_h_geometry(species, coords, capped_h_indices)
        self.enforce_capped_oh_geometry(species, coords, capped_h_indices)
        # Global rule: RDKit/UFF relaxation for capped H only, with all
        # non-capped atoms fixed.
        # RDKit can slightly pull aromatic C-H caps out of the phenyl plane;
        # enforce the final sp2 direction again before writing coordinates.
        self.enforce_sp2_capped_h_geometry(species, coords, capped_h_indices)


    @staticmethod
    def _classify_cap_removal_priority(h_idx, species, coords_arr):
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
            
        h_neighbors = 0
        p_pos = coords_arr[parent_idx]
        for j in range(n):
            if species[j] == "H":
                if float(np.linalg.norm(coords_arr[j] - p_pos)) < 1.3:
                    h_neighbors += 1
        if parent_sym == "O" and h_neighbors >= 3:
            return 100, "overcoordinated-O"
            
        grandparents = []
        for j in range(n):
            if j == h_idx or j == parent_idx or species[j] == "H":
                continue
            d = float(np.linalg.norm(coords_arr[j] - p_pos))
            if d <= 2.2:
                grandparents.append((j, species[j], d))
        
        def _gp_o_count(gp_idx):
            gp_pos = coords_arr[gp_idx]
            count = 0
            for j in range(n):
                if j == gp_idx or j == parent_idx or species[j] != "O":
                    continue
                if float(np.linalg.norm(coords_arr[j] - gp_pos)) <= 1.5:
                    count += 1
            return count

        def _gp_heavy_count(gp_idx):
            gp_pos = coords_arr[gp_idx]
            count = 0
            for j in range(n):
                if j == gp_idx or j == parent_idx or species[j] == "H":
                    continue
                if float(np.linalg.norm(coords_arr[j] - gp_pos)) <= 1.8:
                    count += 1
            return count

        if parent_sym == "O":
            if not grandparents:
                return 30, "bare-O"
            gp_idx, gp_sym, _ = grandparents[0]
            if gp_sym == "C":
                n_o = _gp_o_count(gp_idx)
                if n_o >= 1:
                    return 80, "carboxylate"
                if _gp_heavy_count(gp_idx) >= 3:
                    return 60, "phenol"
                return 50, "alcohol"
            if gp_sym == "S":
                n_o = _gp_o_count(gp_idx)
                if n_o >= 2:
                    return 35, "sulfonate"
                return 70, "thiol"
            if gp_sym == "P":
                return 75, "phosphonate"
            return 40, f"O-on-{gp_sym}"
            
        if parent_sym == "N":
            if not grandparents:
                return 25, "bare-N"
            gp_idx, gp_sym, _ = grandparents[0]
            if gp_sym == "S":
                return 45, "sulfonamide"
            if gp_sym == "C":
                n_o = _gp_o_count(gp_idx)
                if n_o >= 1:
                    return 10, "amide"
            n_heavy = len(grandparents)
            if n_heavy == 1:
                return 20, "primary-amine"
            return 15, "secondary-amine"
            
        if parent_sym == "B":
            return 5, "boron"
            
        return 8, f"{parent_sym}-group"

    def fix_odd_electron_multiplicity(self, species, coords, capped_h_indices, label):
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
            "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83
        }
        z_sum = sum(_ATOMIC_NUMBERS.get(s, 6) for s in species)
        if z_sum % 2 == 0:
            return species, coords, capped_h_indices
            
        n = len(species)
        coords_arr = np.array(coords, dtype=float)
        
        removal_candidates = []
        for h_idx in capped_h_indices:
            score, group = self._classify_cap_removal_priority(h_idx, species, coords_arr)
            if score > 0:
                removal_candidates.append((score, h_idx, group))
                
        if removal_candidates:
            removal_candidates.sort(key=lambda x: x[0], reverse=True)
            score, target_idx, group = removal_candidates[0]
            
            p_sym = "unknown"
            for j in range(n):
                if species[j] != "H" and float(np.linalg.norm(coords_arr[j] - coords_arr[target_idx])) < 1.5:
                    p_sym = species[j]
                    break
                    
            print(f"QM-Fix [{group}]: Removed capping H from {p_sym}[{target_idx}] (priority={score}) to achieve even electron count for '{label}'.")
            
            keep = [i for i in range(n) if i != target_idx]
            species = [species[i] for i in keep]
            coords = [coords[i] for i in keep]
            new_caps = [i for i in capped_h_indices if i != target_idx]
            capped_h_indices = [i - 1 if i > target_idx else i for i in new_caps]
            return species, coords, capped_h_indices
            
        def _heavy_and_h_nbs(atom_idx):
            pos = coords_arr[atom_idx]
            heavy, h_nb = [], []
            for j in range(n):
                if j == atom_idx:
                    continue
                d = float(np.linalg.norm(coords_arr[j] - pos))
                if species[j] == "H":
                    if d < 1.3:
                        h_nb.append(j)
                elif d < 2.2:
                    heavy.append(j)
            return heavy, h_nb

        _ADD_ELIGIBLE = {"O", "N"}
        add_candidates = []
        for i, sym in enumerate(species):
            if sym not in _ADD_ELIGIBLE:
                continue
            heavy_nbs, h_nbs = _heavy_and_h_nbs(i)
            if len(h_nbs) > 0:
                continue
            if sym == "O" and len(heavy_nbs) != 1:
                continue
            if sym == "N" and len(heavy_nbs) != 2:
                continue
                
            simulated_coords = np.vstack([coords_arr, coords_arr[i] + [0, 0, 1.0]])
            score, group = self._classify_cap_removal_priority(n, species + ["H"], simulated_coords)
            add_candidates.append((score, i, group, heavy_nbs))
            
        if add_candidates:
            add_candidates.sort(key=lambda x: x[0], reverse=True)
            score, target_idx, group, heavy_nbs = add_candidates[0]
            target_sym = species[target_idx]
            
            parent_pos = coords_arr[target_idx]
            if len(heavy_nbs) == 1:
                vec = parent_pos - coords_arr[heavy_nbs[0]]
                nn = np.linalg.norm(vec)
                base_vec = vec / nn if nn > 1e-9 else np.array([0.0, 0.0, 1.0])
            else:
                vecs = [parent_pos - coords_arr[nb] for nb in heavy_nbs]
                avg = np.mean(vecs, axis=0)
                nn = np.linalg.norm(avg)
                base_vec = avg / nn if nn > 1e-9 else np.array([0.0, 0.0, 1.0])
                
            bl = self.cap_bond_length(target_sym)
            before = len(species)
            self.place_capping_h(target_idx, base_vec, bl, species, coords, min_hh=1.5, capped_h_flags=None)
            if len(species) > before:
                capped_h_indices.append(len(species) - 1)
                print(f"QM-Fix [{group}]: Added H to {target_sym}[{target_idx}] (no removable non-C capping H found) to achieve even electron count for '{label}'.")
                return species, coords, capped_h_indices
                
        print(f"[QM WARNING] Fragment '{label}': could not automatically fix odd electron count (Z_sum={z_sum}). Manual inspection needed.")
        return species, coords, capped_h_indices

class MOFFragmenter(BaseFragmenter):
    METALS = {"Mg", "Zn", "Cu", "Fe", "Co", "Ni", "Mn", "Zr", "Ti", "V", "Cr", "Al"}
    LARGE_NON_METALS = {"Br", "I", "S", "P", "Cl"}

    def is_valid_bond(self, s1_str, s2_str, dist):
        if s1_str in self.METALS or s2_str in self.METALS:
            return dist < 2.6
        if "H" in (s1_str, s2_str):
            return dist < 1.2
        if s1_str in self.LARGE_NON_METALS or s2_str in self.LARGE_NON_METALS:
            return dist < 2.2
        return dist < 1.8

    @staticmethod
    def _safe_name(text):
        safe = []
        for ch in str(text):
            safe.append(ch if ch.isalnum() or ch in {"-", "_"} else "_")
        return "".join(safe).strip("_") or "fragment"

    @staticmethod
    def _species_coords_unique_key(species, coords, decimals=3):
        species = [str(sp) for sp in species]
        coords = [np.array(c, dtype=float) for c in coords]
        counts = tuple(sorted((sp, species.count(sp)) for sp in set(species)))
        if len(coords) <= 1:
            return counts, ()

        distances = []
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                pair = tuple(sorted((species[i], species[j])))
                d = round(float(np.linalg.norm(coords[i] - coords[j])), decimals)
                distances.append((pair[0], pair[1], d))
        return counts, tuple(sorted(distances))

    @staticmethod
    def _chemical_identity_key(species, coords, decimals=1):
        """Chemical duplicate key based on heavy-atom formula only.
        Hydrogen atoms are excluded so that minor H-capping differences
        do not make otherwise identical skeletons appear unique.
        Pairwise distances are intentionally not used because small
        lattice parameter shifts across CIF files cause systematic
        distance deviations that defeat histogram binning."""
        species = [str(sp) for sp in species]
        # Filter to heavy atoms only
        heavy_species = [sp for sp in species if sp != "H"]
        if not heavy_species:
            return (("H", len(species)),)
        counts = tuple(sorted((sp, heavy_species.count(sp)) for sp in set(heavy_species)))
        return counts

    @classmethod
    def _molecule_unique_key(cls, mol, decimals=3):
        return cls._species_coords_unique_key(mol.species, mol.cart_coords, decimals=decimals)

    @classmethod
    def _xyz_unique_key(cls, path, decimals=3):
        species = []
        coords = []
        for line in Path(path).read_text().splitlines()[2:]:
            parts = line.split()
            if len(parts) < 4:
                continue
            species.append(parts[0])
            coords.append(np.array([float(x) for x in parts[1:4]], dtype=float))
        if not species:
            return None
        return cls._species_coords_unique_key(species, coords, decimals=decimals)



    def _append_merged_atoms(self, species, coords, add_species, add_coords, tol=0.08):
        for sp, coord in zip(add_species, add_coords):
            c = np.array(coord, dtype=float)
            duplicate = False
            for i, old_sp in enumerate(species):
                if old_sp != sp:
                    continue
                if np.linalg.norm(np.array(coords[i], dtype=float) - c) <= tol:
                    duplicate = True
                    break
            if not duplicate:
                species.append(sp)
                coords.append(c)

    def _cap_path_j_open_oxygens(self, species, coords, capped_h_flags):
        heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
        for i in list(heavy_idx):
            if i >= len(species) or species[i] != "O":
                continue
            if self.oxygen_already_protonated(i, species, coords):
                continue
            opos = np.array(coords[i], dtype=float)
            heavy_neighbors = []
            has_metal = False
            for j, spj in enumerate(species):
                if i == j or spj == "H":
                    continue
                d = np.linalg.norm(opos - np.array(coords[j], dtype=float))
                if self.is_valid_bond("O", spj, d):
                    heavy_neighbors.append(j)
                    if spj in self.METALS:
                        has_metal = True
            # Open carboxylate/terminal oxygens from helper linkers have only
            # the carbonyl/carboxyl carbon left after the missing node was cut.
            # Cap them with H using the same placement/refinement machinery as
            # the legacy MOF path.
            if has_metal or len(heavy_neighbors) != 1 or species[heavy_neighbors[0]] not in ("C", "P"):
                continue
            base = opos - np.array(coords[heavy_neighbors[0]], dtype=float)
            cap_len = self.cap_bond_length("O")
            before = len(species)
            self.place_capping_h(i, base, cap_len, species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)
            if len(species) == before and np.linalg.norm(base) > 1e-12:
                self.place_capping_h(i, -base, cap_len, species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

    def _first_connected_ring_fragment(self, linker_species, linker_coords, node_species, node_coords):
        heavy = [i for i, sp in enumerate(linker_species) if sp != "H"]
        if not heavy:
            return [], [], []

        hadj = {i: [] for i in heavy}
        for a in range(len(heavy)):
            i = heavy[a]
            ci = np.array(linker_coords[i], dtype=float)
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                d = np.linalg.norm(ci - np.array(linker_coords[j], dtype=float))
                if self.is_valid_bond(linker_species[i], linker_species[j], d):
                    hadj[i].append(j)
                    hadj[j].append(i)

        bridge_atoms = set()
        for i in heavy:
            ci = np.array(linker_coords[i], dtype=float)
            for nsp, nco in zip(node_species, node_coords):
                d = np.linalg.norm(ci - np.array(nco, dtype=float))
                if self.is_valid_bond(linker_species[i], nsp, d):
                    bridge_atoms.add(i)
                    break
        if not bridge_atoms:
            return [], [], []

        # 1. Find the 2-core (all rings and their interconnects)
        core_atoms = set(heavy)
        changed = True
        while changed:
            changed = False
            to_remove = set()
            for u in core_atoms:
                degree = sum(1 for v in hadj.get(u, []) if v in core_atoms)
                if degree <= 1:
                    to_remove.add(u)
            if to_remove:
                core_atoms -= to_remove
                changed = True

        keep_heavy = set(bridge_atoms)

        # 2. Find shortest path from bridge_atoms to the nearest 2-core component
        if core_atoms:
            q = deque([(ba, [ba]) for ba in bridge_atoms])
            seen = set(bridge_atoms)
            path_to_core = []
            while q:
                curr, path = q.popleft()
                if curr in core_atoms:
                    path_to_core = path
                    break
                for nb in hadj.get(curr, []):
                    if nb not in seen:
                        seen.add(nb)
                        q.append((nb, path + [nb]))
            
            if path_to_core:
                keep_heavy.update(path_to_core)
                core_target = path_to_core[-1]
                
                def is_bridge(u, v, adj_dict, valid_nodes):
                    q_b = deque([u])
                    seen_b = {u}
                    while q_b:
                        curr_b = q_b.popleft()
                        if curr_b == v: return False
                        for nb_b in adj_dict.get(curr_b, []):
                            if nb_b in valid_nodes and nb_b not in seen_b:
                                if (curr_b == u and nb_b == v) or (curr_b == v and nb_b == u):
                                    continue
                                seen_b.add(nb_b)
                                q_b.append(nb_b)
                    return True

                core_comp = {core_target}
                cq = deque([core_target])
                found_cycle = False
                
                while cq:
                    c = cq.popleft()
                    
                    curr_in_cycle = any(not is_bridge(c, nb, hadj, core_atoms) for nb in hadj.get(c, []) if nb in core_atoms)
                    if curr_in_cycle:
                        found_cycle = True
                        
                    for nb in hadj.get(c, []):
                        if nb in core_atoms and nb not in core_comp:
                            edge_is_bridge = is_bridge(c, nb, hadj, core_atoms)
                            
                            if found_cycle:
                                if not edge_is_bridge:
                                    core_comp.add(nb)
                                    cq.append(nb)
                            else:
                                core_comp.add(nb)
                                cq.append(nb)
                                
                keep_heavy.update(core_comp)

        removed_nbr_count = {i: sum(1 for nb in hadj.get(i, []) if nb not in keep_heavy) for i in keep_heavy}

        keep = set(keep_heavy)
        for i, sp in enumerate(linker_species):
            if sp != "H":
                continue
            hpos = np.array(linker_coords[i], dtype=float)
            for j in keep_heavy:
                d = np.linalg.norm(hpos - np.array(linker_coords[j], dtype=float))
                if self.is_valid_bond("H", linker_species[j], d):
                    keep.add(i)
                    break

        ordered = [i for i in range(len(linker_species)) if i in keep]
        out_species = [linker_species[i] for i in ordered]
        out_coords = [np.array(linker_coords[i], dtype=float) for i in ordered]
        out_flags = [False] * len(out_species)
        old_to_new = {old_i: new_i for new_i, old_i in enumerate(ordered)}

        # Cap aromatic/connector carbons that became unsaturated when the
        # partial linker branch was trimmed to the first connected ring.
        for old_i, n_removed in removed_nbr_count.items():
            if n_removed <= 0 or linker_species[old_i] != "C":
                continue
            new_i = old_to_new.get(old_i)
            if new_i is None:
                continue

            cpos = np.array(out_coords[new_i], dtype=float)
            heavy_nbs = []
            h_nbs = 0
            for j, sp in enumerate(out_species):
                if j == new_i:
                    continue
                d = np.linalg.norm(cpos - np.array(out_coords[j], dtype=float))
                if sp == "H" and self.is_valid_bond("C", "H", d):
                    h_nbs += 1
                elif sp != "H" and self.is_valid_bond("C", sp, d):
                    heavy_nbs.append(j)

            # Aromatic/linker carbons should not be left below valence 3 in
            # these first-ring partial branches. Add at most the number of
            # heavy neighbors lost during trimming.
            n_cap = min(n_removed, max(0, 3 - len(heavy_nbs) - h_nbs))
            for _ in range(n_cap):
                base = np.zeros(3)
                if heavy_nbs:
                    for nb in heavy_nbs:
                        base -= np.array(out_coords[nb], dtype=float) - cpos
                else:
                    base = np.array([1.0, 0.0, 0.0])
                before = len(out_species)
                self.place_capping_h(new_i, base, self.cap_bond_length("C"), out_species, out_coords, min_hh=1.5, capped_h_flags=out_flags)
                if len(out_species) == before:
                    break
                h_nbs += 1

        return out_species, out_coords, out_flags

    def _recover_open_node_linkers_from_structure(self, struct, node_species, node_coords, species, coords, minimize=False):
        recovered = []
        if struct is None or not species:
            return recovered

        def current_heavy_neighbors(idx):
            out = []
            ci = np.array(coords[idx], dtype=float)
            for j, sp in enumerate(species):
                if j == idx or sp == "H":
                    continue
                d = np.linalg.norm(ci - np.array(coords[j], dtype=float))
                if self.is_valid_bond(species[idx], sp, d):
                    out.append(j)
            return out

        open_carbons = []
        for i, sp in enumerate(species):
            if sp != "C":
                continue
            nbs = current_heavy_neighbors(i)
            o_ct = sum(1 for nb in nbs if species[nb] == "O")
            c_ct = sum(1 for nb in nbs if species[nb] == "C")
            if o_ct >= 2 and c_ct == 0:
                open_carbons.append(i)
        if not open_carbons:
            return recovered

        lattice = np.array(struct.lattice.matrix, dtype=float)
        image_atoms = []
        for si, site in enumerate(struct):
            sp = site.species_string
            for ia in (-1, 0, 1):
                for ib in (-1, 0, 1):
                    for ic in (-1, 0, 1):
                        shift = ia * lattice[0] + ib * lattice[1] + ic * lattice[2]
                        image_atoms.append((sp, np.array(site.coords, dtype=float) + shift, si, (ia, ib, ic)))

        existing_heavy = [np.array(c, dtype=float) for sp, c in zip(species, coords) if sp != "H"]

        for cidx in open_carbons:
            target = np.array(coords[cidx], dtype=float)
            matches = [
                (np.linalg.norm(pos - target), ai)
                for ai, (sp, pos, _, _) in enumerate(image_atoms)
                if sp == "C" and np.linalg.norm(pos - target) <= 0.35
            ]
            if not matches:
                continue
            _, match_ai = min(matches, key=lambda x: x[0])
            match_pos = image_atoms[match_ai][1]

            local = []
            for ai, (sp, pos, si, img) in enumerate(image_atoms):
                if np.linalg.norm(pos - target) <= 35.0:
                    local.append((ai, sp, pos, si, img))
            local_ids = {ai: li for li, (ai, _, _, _, _) in enumerate(local)}
            match_li = local_ids.get(match_ai)
            if match_li is None:
                continue

            start_candidates = []
            for li, (ai, sp, pos, _, _) in enumerate(local):
                if li == match_li or sp != "C":
                    continue
                d = np.linalg.norm(pos - match_pos)
                if self.is_valid_bond("C", "C", d):
                    # Prefer the carbon that is not already present in the assembled fragment.
                    nearest_existing = min((np.linalg.norm(pos - old) for old in existing_heavy), default=float("inf"))
                    start_candidates.append((nearest_existing, d, li))
            start_candidates = [x for x in start_candidates if x[0] > 0.35]
            if not start_candidates:
                continue
            _, _, start_li = min(start_candidates, key=lambda x: (x[1], -x[0]))

            adj = {li: [] for li in range(len(local))}
            for a in range(len(local)):
                _, spa, posa, _, _ = local[a]
                if spa in self.METALS:
                    continue
                for b in range(a + 1, len(local)):
                    _, spb, posb, _, _ = local[b]
                    if spb in self.METALS:
                        continue
                    d = np.linalg.norm(posa - posb)
                    if self.is_valid_bond(spa, spb, d):
                        adj[a].append(b)
                        adj[b].append(a)

            q = deque([start_li, match_li])
            seen = {start_li, match_li}
            comp = []
            while q and len(comp) <= 120:
                u = q.popleft()
                comp.append(u)
                for v in adj.get(u, []):
                    if v not in seen:
                        seen.add(v)
                        q.append(v)
            if len(comp) > 120:
                continue

            comp_species = [local[li][1] for li in comp]
            comp_coords = [local[li][2] for li in comp]
            if minimize:
                psp, pco, pflags = self._first_connected_ring_fragment(comp_species, comp_coords, node_species, node_coords)
                if psp:
                    recovered.append((psp, pco, pflags))
            else:
                recovered.append((comp_species, comp_coords, [False] * len(comp_species)))
        return recovered

    def _path_j_second_node_family(self, stem):
        name = self._safe_name(stem).upper()
        return name.startswith("PCN-") or name.startswith("PCN_") or name.startswith("NU-") or name.startswith("NU_")

    def _one_node_paddlewheel_family(self, stem):
        name = self._safe_name(stem).upper()
        return name.startswith("CU-BTC") or name.startswith("DUT-49")

    def _export_moffragmentor_library(self, result, stem):
        node_dir = Path("mof_nodes_lib")
        linker_dir = Path("mof_linkers_lib")
        node_dir.mkdir(exist_ok=True)
        linker_dir.mkdir(exist_ok=True)
        file_stem = self._safe_name(stem)

        for collection, out_dir in ((getattr(result, "nodes", []), node_dir), (getattr(result, "linkers", []), linker_dir)):
            if any(out_dir.glob(f"{file_stem}_*.xyz")):
                continue

            existing_keys = set()
            for old in out_dir.glob("*.xyz"):
                try:
                    key = self._xyz_unique_key(old)
                except Exception:
                    key = None
                if key is not None:
                    existing_keys.add(key)

            sbus = getattr(collection, "sbus", list(collection) if collection is not None else [])
            seen = set(existing_keys)
            out_idx = 0
            for sbu in sbus:
                mol = getattr(sbu, "molecule", None)
                if mol is None:
                    continue
                key = self._molecule_unique_key(mol)
                if key in seen:
                    continue
                seen.add(key)
                out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"
                while out.exists():
                    out_idx += 1
                    out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"
                mol.to(filename=str(out), fmt="xyz")
                out_idx += 1

    def _export_mof_component_library(self, species, coords, stem, kind):
        out_dir = Path("mof_nodes_lib") if kind == "node" else Path("mof_linkers_lib")
        out_dir.mkdir(exist_ok=True)
        file_stem = self._safe_name(stem)

        if any(out_dir.glob(f"{file_stem}_*.xyz")):
            return False

        existing_keys = set()
        for old in out_dir.glob("*.xyz"):
            try:
                key = self._xyz_unique_key(old)
            except Exception:
                key = None
            if key is not None:
                existing_keys.add(key)

        key = self._species_coords_unique_key(species, coords)
        if key in existing_keys:
            return False

        out_idx = 0
        out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"
        while out.exists():
            out_idx += 1
            out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"

        Molecule(species, coords).to(filename=str(out), fmt="xyz")
        return True

    def _fallback_export_mof_node_linker(self, mof_path):
        try:
            struct = Structure.from_file(mof_path)
        except Exception:
            return False

        symbols = [str(site.species_string) for site in struct]
        neighs = struct.get_all_neighbors(r=3.6)
        graph = [[] for _ in range(len(struct))]
        for i, ns in enumerate(neighs):
            si = symbols[i]
            ci = struct[i].coords
            for n in ns:
                j = n.index
                if j <= i:
                    continue
                d = float(np.linalg.norm(ci - n.coords))
                if self.is_valid_bond(si, symbols[j], d):
                    graph[i].append(j)
                    graph[j].append(i)

        metals = [i for i, s in enumerate(symbols) if s in self.METALS]
        if not metals:
            return False

        center = struct.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
        center_m = min(metals, key=lambda i: float(np.linalg.norm(struct[i].coords - center)))

        node = {center_m}
        q = deque([center_m])
        while q:
            u = q.popleft()
            for v in graph[u]:
                if symbols[v] in self.METALS and v not in node:
                    node.add(v)
                    q.append(v)

        expanded_node = set(node)
        for u in list(node):
            for v in graph[u]:
                if symbols[v] not in self.METALS:
                    expanded_node.add(v)

        linker_comp = set()
        vis = set(expanded_node)
        seeds = []
        for u in expanded_node:
            for v in graph[u]:
                if v not in expanded_node and symbols[v] != "H":
                    seeds.append(v)

        for st in seeds:
            if st in vis:
                continue
            comp = {st}
            qq = deque([st])
            vis.add(st)
            while qq:
                u = qq.popleft()
                for v in graph[u]:
                    if v in vis or v in expanded_node or symbols[v] in self.METALS:
                        continue
                    vis.add(v)
                    comp.add(v)
                    qq.append(v)
            if len(comp) > len(linker_comp):
                linker_comp = comp

        stem = Path(mof_path).stem
        node_sp = [symbols[i] for i in sorted(expanded_node)]
        node_co = [np.array(struct[i].coords, dtype=float) for i in sorted(expanded_node)]
        node_written = self._export_mof_component_library(node_sp, node_co, stem, "node")

        linker_written = False
        if linker_comp:
            lk = sorted(linker_comp)
            lsp = [symbols[i] for i in lk]
            lco = [np.array(struct[i].coords, dtype=float) for i in lk]
            linker_written = self._export_mof_component_library(lsp, lco, stem, "linker")

        if node_written or linker_written:
            print("  -> Helper fragments available in mof_nodes_lib/ and mof_linkers_lib/ (fallback export).")
            return True
        return False

    def _try_moffragmentor_node_linker_fragment(self, mof_path, output_path, minimize=False):
        try:
            from moffragmentor import MOF
        except Exception as exc:
            print(f"  moffragmentor unavailable: {exc}")
            return None
        try:
            struct = Structure.from_file(mof_path)
            result = MOF.from_cif(mof_path).fragment()
        except Exception as exc:
            print(f"  moffragmentor node+linker failed: {exc}")
            return None

        if getattr(result, "has_1d_sbu", False):
            print("  moffragmentor detected 1D SBU. Skipping moffragmentor node extraction to allow --nmetals slicing via fallback path.")
            return None

        nodes = list(getattr(getattr(result, "nodes", []), "sbus", []))
        linkers = list(getattr(getattr(result, "linkers", []), "sbus", []))
        if not nodes or not linkers:
            return None

        stem = Path(mof_path).stem
        self._export_moffragmentor_library(result, stem)

        center = struct.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
        node = min(nodes, key=lambda s: np.linalg.norm(np.array(getattr(s, "center", np.mean(s.molecule.cart_coords, axis=0))) - center))
        node_sp = [str(x) for x in node.molecule.species]
        node_co = [np.array(c, dtype=float) for c in node.molecule.cart_coords]

        # --- Merge nearby nodes that belong to the same polynuclear SBU ---
        # moffragmentor may split a dinuclear/polynuclear SBU into individual
        # 1-atom metal nodes.  Detect this by checking metal-metal distances
        # (with periodic images) and merge when distance < 3.5 Å.
        _merge_cutoff = 3.5
        _lat = struct.lattice.matrix
        merged_ids = {id(node)}
        changed = True
        while changed:
            changed = False
            node_metals = [node_co[k] for k in range(len(node_sp)) if node_sp[k] in self.METALS]
            for other in nodes:
                if id(other) in merged_ids:
                    continue
                osp = [str(x) for x in other.molecule.species]
                oco = [np.array(c, dtype=float) for c in other.molecule.cart_coords]
                should_merge = False
                best_shift = np.zeros(3)
                best_d = float("inf")
                for ia in (-1, 0, 1):
                    for ib in (-1, 0, 1):
                        for ic in (-1, 0, 1):
                            shift = ia * _lat[0] + ib * _lat[1] + ic * _lat[2]
                            for om_pos in [oco[k] + shift for k in range(len(osp)) if osp[k] in self.METALS]:
                                for nm_pos in node_metals:
                                    d = float(np.linalg.norm(om_pos - nm_pos))
                                    if 0.1 < d < _merge_cutoff and d < best_d:
                                        should_merge = True
                                        best_shift = shift
                                        best_d = d
                if should_merge:
                    for k in range(len(osp)):
                        node_sp.append(osp[k])
                        node_co.append(oco[k] + best_shift)
                    merged_ids.add(id(other))
                    changed = True  # re-scan in case of chains of bridging metals

        if len(merged_ids) > 1:
            print(f"  -> Merged {len(merged_ids)} moffragmentor nodes into one polynuclear SBU.")
        node_ctr = np.mean(node_co, axis=0)

        image_vectors = []
        for ia in (-1, 0, 1):
            for ib in (-1, 0, 1):
                for ic in (-1, 0, 1):
                    image_vectors.append(
                        ia * np.array(struct.lattice.matrix[0], dtype=float)
                        + ib * np.array(struct.lattice.matrix[1], dtype=float)
                        + ic * np.array(struct.lattice.matrix[2], dtype=float)
                    )

        def linker_image_score(linker, image_shift):
            lsp = [str(x) for x in linker.molecule.species]
            lco = [np.array(c, dtype=float) + image_shift for c in linker.molecule.cart_coords]
            metal_bonds = 0
            min_d = float("inf")
            for i, si in enumerate(node_sp):
                for j, sj in enumerate(lsp):
                    d = float(np.linalg.norm(node_co[i] - lco[j]))
                    min_d = min(min_d, d)
                    if (si in self.METALS) != (sj in self.METALS) and self.is_valid_bond(si, sj, d):
                        metal_bonds += 1
            lctr = np.mean(lco, axis=0)
            return (metal_bonds, -min_d, -float(np.linalg.norm(lctr - node_ctr)))

        scored_images = []
        for linker in linkers:
            for image_shift in image_vectors:
                score = linker_image_score(linker, image_shift)
                if score[0] > 0:
                    scored_images.append((score, linker, image_shift))
        if not scored_images:
            return None
        scored_images.sort(key=lambda x: x[0], reverse=True)
        normal_size = len(node_sp) + sum(len(l.molecule) for _, l, _ in scored_images)
        should_minimize = minimize and normal_size >= 50

        selected_images = scored_images[:1] if should_minimize else scored_images
        partial_images = scored_images[1:] if should_minimize else []

        species = []
        coords = []
        partial_capped_h = []
        self._append_merged_atoms(species, coords, node_sp, node_co)
        for _, linker, image_shift in selected_images:
            lsp = [str(x) for x in linker.molecule.species]
            lco = [np.array(c, dtype=float) + image_shift for c in linker.molecule.cart_coords]
            self._append_merged_atoms(species, coords, lsp, lco)
        for _, linker, image_shift in partial_images:
            lsp = [str(x) for x in linker.molecule.species]
            lco = [np.array(c, dtype=float) + image_shift for c in linker.molecule.cart_coords]
            psp, pco, pflags = self._first_connected_ring_fragment(lsp, lco, node_sp, node_co)
            for p_sp, p_coord, p_flag in zip(psp, pco, pflags):
                before = len(species)
                self._append_merged_atoms(species, coords, [p_sp], [p_coord])
                if p_flag and len(species) > before and species[-1] == "H":
                    partial_capped_h.append(len(species) - 1)
        # --- Coordination completion: Ensure full coordination sphere by adding any stray coordinating atoms ---
        all_species_coords = list(zip(species, coords))
        stray_found = False
        
        for i, sp_m in enumerate(node_sp):
            if sp_m not in self.METALS:
                continue
            m_pos = node_co[i]
            # Find all neighbors in struct (using images)
            for ia in (-1, 0, 1):
                for ib in (-1, 0, 1):
                    for ic in (-1, 0, 1):
                        shift = ia * struct.lattice.matrix[0] + ib * struct.lattice.matrix[1] + ic * struct.lattice.matrix[2]
                        for s_idx, site in enumerate(struct):
                            d = np.linalg.norm(site.coords + shift - m_pos)
                            if 0.1 < d < 2.6: # Bonded
                                # Check if this atom is already in species/coords
                                is_missing = True
                                for ex_sp, ex_co in all_species_coords:
                                    if ex_sp == site.species_string and np.linalg.norm(ex_co - (site.coords + shift)) < 0.5:
                                        is_missing = False
                                        break
                                if is_missing:
                                    # BFS to find the full connected organic component for this stray atom
                                    q_stray = deque([(site.coords + shift, site.species_string, s_idx, (ia, ib, ic))])
                                    visited_stray = {(s_idx, (ia, ib, ic))}
                                    
                                    while q_stray:
                                        curr_pos, curr_sp, curr_si, curr_img = q_stray.popleft()
                                        
                                        # Add to result if not a duplicate
                                        dup_in_result = False
                                        curr_idx_in_res = -1
                                        for ex_idx, (ex_sp, ex_co) in enumerate(all_species_coords):
                                            if ex_sp == curr_sp and np.linalg.norm(ex_co - curr_pos) < 0.5:
                                                dup_in_result = True
                                                curr_idx_in_res = ex_idx
                                                break
                                        if not dup_in_result:
                                            before_len = len(species)
                                            self._append_merged_atoms(species, coords, [curr_sp], [curr_pos])
                                            if len(species) > before_len:
                                                curr_idx_in_res = len(species) - 1
                                                all_species_coords.append((curr_sp, curr_pos))
                                                stray_found = True
                                            else:
                                                # somehow it was a duplicate after all or append failed
                                                continue
                                        capping_hs_added = 0
                                        # Find neighbors of this stray atom to continue BFS or cap
                                        for nia in (-1, 0, 1):
                                            for nib in (-1, 0, 1):
                                                for nic in (-1, 0, 1):
                                                    n_shift = (curr_img[0]+nia)*struct.lattice.matrix[0] + (curr_img[1]+nib)*struct.lattice.matrix[1] + (curr_img[2]+nic)*struct.lattice.matrix[2]
                                                    for n_s_idx, n_site in enumerate(struct):
                                                        n_pos = n_site.coords + n_shift
                                                        nd = np.linalg.norm(n_pos - curr_pos)
                                                        if not self.is_valid_bond(curr_sp, n_site.species_string, nd) or nd < 0.1:
                                                            continue

                                                        if n_site.species_string in self.METALS:
                                                            # If it's a metal NOT in our core node, we need to cap this bond
                                                            is_core_metal = False
                                                            for m_c in node_co:
                                                                if np.linalg.norm(m_c - n_pos) < 0.5:
                                                                    is_core_metal = True
                                                                    break
                                                            if not is_core_metal:
                                                                if curr_sp == "O":
                                                                    h_count_struct = sum(1 for nn in struct.get_neighbors(struct[curr_si], 1.3) if nn.species_string == "H")
                                                                    # Check if this O has a heavy non-metal neighbor (C, P, S, etc.)
                                                                    # If so, it's a coordinating O (e.g. P-O-Metal, C-O-Metal),
                                                                    # and should get at most 1 H total to replace the metal bond.
                                                                    # Only pure water-like O (bonded only to metals/H) can have 2 H.
                                                                    has_heavy_nonmetal = any(
                                                                        nn.species_string not in self.METALS and nn.species_string != "H"
                                                                        for nn in struct.get_neighbors(struct[curr_si], 1.8)
                                                                    )
                                                                    max_h = 1 if has_heavy_nonmetal else 2
                                                                    if (h_count_struct + capping_hs_added) >= max_h:
                                                                        continue
                                                                # CAP IT
                                                                base = n_pos - curr_pos
                                                                before = len(species)
                                                                capped_h_flags_temp = [False] * (len(species) + 1) # dummy
                                                                self.place_capping_h(curr_idx_in_res, base, self.cap_bond_length(curr_sp), 
                                                                                      species, coords, min_hh=1.2, min_o_contact=1.2, capped_h_flags=capped_h_flags_temp)
                                                                if len(species) > before:
                                                                    # The new H index is before
                                                                    partial_capped_h.append(before)
                                                                    all_species_coords.append((species[before], coords[before]))
                                                                    capping_hs_added += 1
                                                            continue

                                                        if (n_s_idx, (curr_img[0]+nia, curr_img[1]+nib, curr_img[2]+nic)) not in visited_stray:
                                                            visited_stray.add((n_s_idx, (curr_img[0]+nia, curr_img[1]+nib, curr_img[2]+nic)))
                                                            q_stray.append((n_pos, n_site.species_string, n_s_idx, (curr_img[0]+nia, curr_img[1]+nib, curr_img[2]+nic)))

        if stray_found:
            print("  -> Coordination completion: Added stray coordinating groups.")

        recovered_branches = self._recover_open_node_linkers_from_structure(struct, node_sp, node_co, species, coords, minimize=minimize)
        for rsp, rco, rflags in recovered_branches:
            for r_sp, r_coord, r_flag in zip(rsp, rco, rflags):
                before = len(species)
                self._append_merged_atoms(species, coords, [r_sp], [r_coord])
                if r_flag and len(species) > before and species[-1] == "H":
                    partial_capped_h.append(len(species) - 1)

        capped_h_flags = [False] * len(species)
        for hidx in partial_capped_h:
            if 0 <= hidx < len(capped_h_flags):
                capped_h_flags[hidx] = True
        self._cap_path_j_open_oxygens(species, coords, capped_h_flags)
        capped_h_indices = [i for i, is_cap in enumerate(capped_h_flags) if is_cap and species[i] == "H"]
        # Path J keeps helper node/linker heavy atoms fixed. Only capped H
        # coordinates are locally adjusted.
        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)
        _label = str(output_path or 'mof_fragment')
        species, coords, capped_h_indices = self.fix_odd_electron_multiplicity(species, coords, capped_h_indices, label=_label)

        mol = Molecule(species, coords)
        if output_path:
            mol.to(filename=output_path, fmt="xyz")
            print("  -> MOF Path J (moffragmentor node+linker combine).")
            print(f"  -> Helper fragments available in mof_nodes_lib/ and mof_linkers_lib/.")
            print(f"Final size: {len(species)} atoms.")
            print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)

    def extract(self, mof_path, center_idx=-1, nmetals=3, output_path="fragment.xyz", minimize=False):
        print(f"Loading '{mof_path}'...")
        combined = self._try_moffragmentor_node_linker_fragment(mof_path, output_path, minimize=minimize)
        if combined is not None:
            return combined
        self._fallback_export_mof_node_linker(mof_path)
        struct = Structure.from_file(mof_path)
        structure_stem = Path(mof_path).stem
        if nmetals < 1:
            raise ValueError(f"--nmetals must be >= 1 (got {nmetals}).")

        metals = self.METALS
        user_center = center_idx != -1
        if user_center:
            if center_idx < 0 or center_idx >= len(struct):
                raise IndexError(f"--center index {center_idx} is out of range [0, {len(struct)-1}]")
            if struct[center_idx].species_string not in metals:
                raise ValueError(
                    f"--center index {center_idx} is {struct[center_idx].species_string}, not a metal in {sorted(metals)}"
                )

        print("Creating supercell...")
        dims = [max(1, int(np.ceil(30.0 / a))) for a in struct.lattice.abc]
        dims = [max(3, d) if a < 15.0 else d for d, a in zip(dims, struct.lattice.abc)]
        supercell = struct * dims

        sc_center_cart = supercell.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
        best_dist = float("inf")
        sc_center_idx = -1

        if not user_center:
            for i, site in enumerate(supercell):
                if site.species_string not in metals:
                    continue
                d = np.linalg.norm(site.coords - sc_center_cart)
                if d < best_dist:
                    best_dist = d
                    sc_center_idx = i
        else:
            target_site = struct[center_idx]
            target_frac = np.array(target_site.frac_coords)
            dims_arr = np.array(dims, dtype=float)
            for i, site in enumerate(supercell):
                if site.species_string != target_site.species_string:
                    continue
                parent_frac = np.mod(np.array(site.frac_coords) * dims_arr, 1.0)
                frac_delta = np.abs(np.mod(parent_frac - target_frac + 0.5, 1.0) - 0.5)
                if np.max(frac_delta) > 1e-3:
                    continue
                d = np.linalg.norm(site.coords - sc_center_cart)
                if d < best_dist:
                    best_dist = d
                    sc_center_idx = i

        if sc_center_idx == -1:
            if user_center:
                raise ValueError(f"Could not map --center index {center_idx} into the generated supercell.")
            raise ValueError("No metal found in the input structure.")
        sc_center_site = supercell[sc_center_idx]

        print("Detecting topology...")
        sbu_all_neighs = supercell.get_all_neighbors(r=3.6)
        sbu_metals = {sc_center_idx}
        sbu_queue = deque([sc_center_idx])
        is_infinite_sbu = False
        while sbu_queue:
            curr = sbu_queue.popleft()
            for n in sbu_all_neighs[curr]:
                if n.species_string in metals and n.index not in sbu_metals:
                    sbu_metals.add(n.index)
                    sbu_queue.append(n.index)
            if len(sbu_metals) > 20:
                is_infinite_sbu = True
                break

        if not is_infinite_sbu and len(sbu_metals) > 1:
            coords = np.array([supercell[i].coords for i in sbu_metals])
            dists = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=-1))
            if dists.max() > min(struct.lattice.abc) * 0.5:
                is_infinite_sbu = True

        zif_mode = False
        if not is_infinite_sbu and len(sbu_metals) == 1:
            first_shell = supercell.get_neighbors(sc_center_site, 2.3)
            n_neighbors = sum(1 for n in first_shell if n.species_string == "N")
            if n_neighbors >= 3:
                zif_mode = True

        sc_neighbors = supercell.get_neighbors(sc_center_site, self.radius)
        initial_indices = {sc_center_idx}
        for n in sc_neighbors:
            if zif_mode:
                if n.species_string not in metals:
                    initial_indices.add(n.index)
            elif not is_infinite_sbu:
                initial_indices.add(n.index)
            else:
                if n.species_string not in metals:
                    initial_indices.add(n.index)

        is_one_node_paddlewheel = self._one_node_paddlewheel_family(structure_stem)
        is_cu_paddlewheel_pair = (
            is_one_node_paddlewheel
            and len(sbu_metals) == 2
            and all(supercell[m].species_string == "Cu" for m in sbu_metals)
        )

        if is_infinite_sbu:
            print(f"  -> Path B (Infinite). Metals: {nmetals}")
            all_sc_m = []
            for i, site in enumerate(supercell):
                if site.species_string in metals:
                    all_sc_m.append((np.linalg.norm(site.coords - sc_center_site.coords), i))
            all_sc_m.sort()
            core_metals = {idx for _, idx in all_sc_m[:nmetals]}
            if len(core_metals) < nmetals:
                raise ValueError(
                    f"Requested --nmetals={nmetals}, but only found {len(core_metals)} metals in generated supercell."
                )
            initial_indices.update(core_metals)
            is_path_c = False
        elif len(sbu_metals) == 2 and not minimize and not is_cu_paddlewheel_pair:
            print(f"  -> Path C (Discrete, 2 SBUs). Auto-detected small SBU size: {len(sbu_metals)}")
            c0 = np.mean([supercell[m].coords for m in sbu_metals], axis=0)
            q = deque([(list(sbu_metals)[0], 0)])
            visited = set(sbu_metals)
            struct_all_neighs = supercell.get_all_neighbors(r=3.6)

            found_sbus = []
            while q:
                curr_idx, dist = q.popleft()
                curr_site = supercell[curr_idx]
                for n in struct_all_neighs[curr_idx]:
                    n_idx = n.index
                    n_dist = np.linalg.norm(n.coords - curr_site.coords)
                    if not self.is_valid_bond(n.species_string, curr_site.species_string, n_dist):
                        continue
                    if n.species_string in metals and n_idx not in sbu_metals:
                        m_comp = {n_idx}
                        mq = deque([n_idx])
                        mv = {n_idx}
                        while mq:
                            mc = mq.popleft()
                            for mn in struct_all_neighs[mc]:
                                if mn.species_string in metals and mn.index not in mv:
                                    md = np.linalg.norm(mn.coords - supercell[mc].coords)
                                    if md < 3.6:
                                        mv.add(mn.index)
                                        m_comp.add(mn.index)
                                        mq.append(mn.index)
                        if len(m_comp) == len(sbu_metals):
                            c1 = np.mean([supercell[m].coords for m in m_comp], axis=0)
                            if not any(np.linalg.norm(x[2] - c1) < 1.0 for x in found_sbus):
                                found_sbus.append((m_comp, dist + 1, c1))
                    elif n_idx not in visited:
                        visited.add(n_idx)
                        q.append((n_idx, dist + 1))

            if not found_sbus:
                all_metal_indices = [i for i, s in enumerate(supercell) if s.species_string in metals]
                seen_m = set()
                metal_components = []
                for m in all_metal_indices:
                    if m in seen_m:
                        continue
                    comp = {m}
                    mq = deque([m])
                    seen_m.add(m)
                    while mq:
                        mc = mq.popleft()
                        for mn in struct_all_neighs[mc]:
                            if mn.species_string in metals and mn.index not in seen_m:
                                md = np.linalg.norm(mn.coords - supercell[mc].coords)
                                if md < 3.6:
                                    seen_m.add(mn.index)
                                    comp.add(mn.index)
                                    mq.append(mn.index)
                    metal_components.append(comp)

                c_ref = np.mean([supercell[m].coords for m in sbu_metals], axis=0)
                fallback_candidates = []
                for comp in metal_components:
                    if len(comp) != len(sbu_metals):
                        continue
                    c_comp = np.mean([supercell[m].coords for m in comp], axis=0)
                    if np.linalg.norm(c_comp - c_ref) < 1.0:
                        continue
                    fallback_candidates.append((np.linalg.norm(c_comp - c_ref), comp, c_comp))

                if fallback_candidates:
                    fallback_candidates.sort(key=lambda x: x[0])
                    best_dist, best_comp, best_center = fallback_candidates[0]
                    print(f"     Fallback Path C: selected nearest matching SBU at {best_dist:.2f} A")
                    found_sbus.append((best_comp, 1, best_center))

            if found_sbus:
                c1 = np.mean([supercell[m].coords for m in sbu_metals], axis=0)
                found_sbus.sort(key=lambda x: np.linalg.norm(x[2] - c1))
                max_candidates = min(20, len(found_sbus))
                candidates_to_test = found_sbus[:max_candidates]
                print(f"     Testing {len(candidates_to_test)} SBU candidates to minimize atom count...")

                best_size = float("inf")
                best_result = None
                best_dist = float("inf")
                best_valid_size = float("inf")
                best_valid_result = None
                best_valid_dist = float("inf")
                expected_metals = len(sbu_metals) * 2
                for idx, candidate in enumerate(candidates_to_test):
                    test_core = set(sbu_metals) | candidate[0]
                    test_init = set(initial_indices) | test_core
                    sp, co = self._get_fragment(
                        test_core,
                        test_init,
                        supercell,
                        sc_center_idx,
                        is_infinite_sbu,
                        nmetals,
                        minimize,
                        zif_mode,
                    )
                    sz = len(sp)
                    metal_ct = sum(1 for x in sp if x in metals)
                    dist_val = np.linalg.norm(candidate[2] - c1)
                    print(f"       Candidate {idx + 1} (dist {dist_val:.2f} A) -> {sz} atoms, {metal_ct} metals")
                    if sz < best_size or (sz == best_size and dist_val < best_dist):
                        best_size = sz
                        best_result = (sp, co)
                        best_dist = dist_val
                    if metal_ct >= expected_metals:
                        if sz < best_valid_size or (sz == best_valid_size and dist_val < best_valid_dist):
                            best_valid_size = sz
                            best_valid_result = (sp, co)
                            best_valid_dist = dist_val

                if best_valid_result is not None:
                    print(f"     Selected metal-complete SBU candidate yielding {best_valid_size} atoms.")
                    final_species, final_coords = best_valid_result
                else:
                    print(f"     No metal-complete candidate found. Selected smallest candidate yielding {best_size} atoms.")
                    final_species, final_coords = best_result
            else:
                print("     Could not find adjacent SBU. Reverting to 1 SBU.")
                core_metals = set(sbu_metals)
                initial_indices.update(core_metals)
                final_species, final_coords = self._get_fragment(
                    core_metals,
                    initial_indices,
                    supercell,
                    sc_center_idx,
                    is_infinite_sbu,
                    nmetals,
                    minimize,
                    zif_mode,
                )

            is_path_c = True
        else:
            if zif_mode:
                print(f"  -> Path D (ZIF-like). SBU size: {len(sbu_metals)}")
            else:
                print(f"  -> Path A (Discrete). SBU size: {len(sbu_metals)}")
            core_metals = sbu_metals
            initial_indices.update(core_metals)
            is_path_c = False

        if not is_path_c:
            final_species, final_coords = self._get_fragment(
                core_metals,
                initial_indices,
                supercell,
                sc_center_idx,
                is_infinite_sbu,
                nmetals,
                minimize,
                zif_mode,
            )

        bool_flags = [i in getattr(self, "_last_capped_h_indices", []) for i in range(len(final_species))]
        final_species, final_coords, bool_flags = self.enforce_single_molecule(final_species, final_coords, bool_flags)
        self._last_capped_h_indices = [i for i, f in enumerate(bool_flags) if f]

        self.optimize_capped_h_geometry_only(final_species, final_coords, self._last_capped_h_indices)
        _label = str(output_path or 'mof_fragment')
        final_species, final_coords, self._last_capped_h_indices = self.fix_odd_electron_multiplicity(final_species, final_coords, self._last_capped_h_indices, label=_label)
        if output_path:
            mol = Molecule(final_species, final_coords)
            mol.to(filename=output_path, fmt="xyz")
            print(f"Saved: {output_path}")
        return FragmentResult(species=final_species, coords=final_coords)

    def _get_fragment(self, core_metals, initial_indices, supercell, sc_center_idx, is_infinite_sbu, nmetals, minimize, zif_mode=False):
        metals = self.METALS
        sc_all_neighbors = supercell.get_all_neighbors(r=3.0)

        seeded = set()
        for idx in initial_indices:
            if supercell[idx].species_string in metals and idx not in core_metals:
                continue
            seeded.add(idx)
        final_indices = set(seeded)
        queue = deque(seeded)
        visited = set(seeded)
        broken_bonds = []

        unwrapped_coords = {idx: supercell[idx].coords for idx in initial_indices}

        while queue:
            curr_idx = queue.popleft()
            curr_site = supercell[curr_idx]

            if is_infinite_sbu and curr_site.species_string in metals and curr_idx != sc_center_idx:
                continue
            if (not is_infinite_sbu) and curr_site.species_string in metals and curr_idx not in core_metals:
                continue

            for n in sc_all_neighbors[curr_idx]:
                n_idx = n.index
                n_dist = np.linalg.norm(n.coords - curr_site.coords)
                if not self.is_valid_bond(n.species_string, curr_site.species_string, n_dist):
                    continue

                if n.species_string in metals:
                    if n_idx not in core_metals:
                        unwrapped_nb_pos = unwrapped_coords[curr_idx] + (n.coords - curr_site.coords)
                        broken_bonds.append((curr_idx, unwrapped_nb_pos))
                    elif n_idx not in visited:
                        final_indices.add(n_idx)
                        visited.add(n_idx)
                        queue.append(n_idx)
                        unwrapped_coords[n_idx] = unwrapped_coords[curr_idx] + (n.coords - curr_site.coords)
                elif n_idx not in visited:
                    final_indices.add(n_idx)
                    visited.add(n_idx)
                    queue.append(n_idx)
                    unwrapped_coords[n_idx] = unwrapped_coords[curr_idx] + (n.coords - curr_site.coords)

        if is_infinite_sbu and nmetals > 1:
            print("Completing coordination for edge metals...")
            edge_metals = core_metals - {sc_center_idx}
            comp_queue = deque(edge_metals)

            while comp_queue:
                idx = comp_queue.popleft()
                site = supercell[idx]

                if site.species_string in metals and idx not in core_metals:
                    continue

                for n in sc_all_neighbors[idx]:
                    n_idx = n.index
                    n_dist = np.linalg.norm(n.coords - site.coords)
                    if not self.is_valid_bond(n.species_string, site.species_string, n_dist):
                        continue

                    if n.species_string in metals:
                        if n_idx not in core_metals:
                            unwrapped_nb_pos = unwrapped_coords[idx] + (n.coords - site.coords)
                            broken_bonds.append((idx, unwrapped_nb_pos))
                    elif n_idx not in visited:
                        final_indices.add(n_idx)
                        visited.add(n_idx)
                        comp_queue.append(n_idx)
                        unwrapped_coords[n_idx] = unwrapped_coords[idx] + (n.coords - site.coords)

        print("Pruning partial linkers...")
        organic_indices = {idx for idx in final_indices if supercell[idx].species_string not in metals}
        local_adj = {idx: [] for idx in final_indices}
        for idx in final_indices:
            s1 = supercell[idx]
            for n in sc_all_neighbors[idx]:
                if n.index in final_indices:
                    nd = np.linalg.norm(n.coords - s1.coords)
                    if self.is_valid_bond(s1.species_string, n.species_string, nd):
                        local_adj[idx].append(n.index)

        org_visited = set()
        components = []
        for seed in organic_indices:
            if seed in org_visited:
                continue
            component = {seed}
            org_visited.add(seed)
            q = deque([seed])
            while q:
                cur = q.popleft()
                for nb in local_adj[cur]:
                    if nb not in org_visited and supercell[nb].species_string not in metals:
                        org_visited.add(nb)
                        component.add(nb)
                        q.append(nb)
            if component:
                components.append(component)

        normal_size_estimate = len(final_indices) + len(broken_bonds)
        should_minimize = minimize and normal_size_estimate >= 50

        if should_minimize:
            partial_to_remove = set()
            bridge_atoms_to_cap = []
            linker_evaluations = []

            for comp in components:
                touching_metals = set()
                bridge_atoms = set()
                for atom_idx in comp:
                    for nb in local_adj[atom_idx]:
                        if nb in core_metals:
                            touching_metals.add(nb)
                            bridge_atoms.add(atom_idx)

                remove_cond = len(touching_metals) < 1
                if remove_cond:
                    atoms_to_cut = comp - bridge_atoms
                    partial_to_remove.update(atoms_to_cut)
                    for ba in bridge_atoms:
                        removed_neighbor_coords = []
                        for nb in local_adj[ba]:
                            if nb in atoms_to_cut:
                                removed_neighbor_coords.append(unwrapped_coords[nb])
                        if removed_neighbor_coords:
                            bridge_atoms_to_cap.append((ba, removed_neighbor_coords))
                    continue

                keep_linker = False
                touched_coords = [supercell[m].coords for m in touching_metals]
                max_dist = 0
                for i in range(len(touched_coords)):
                    for j in range(i + 1, len(touched_coords)):
                        d = np.linalg.norm(touched_coords[i] - touched_coords[j])
                        if d > max_dist:
                            max_dist = d
                if max_dist > 4.5:
                    keep_linker = True

                linker_evaluations.append((comp, bridge_atoms, keep_linker))

            kept_any_full_linker = any(keep for _, _, keep in linker_evaluations)
            if not kept_any_full_linker and linker_evaluations:
                linker_evaluations.sort(key=lambda x: len(x[0]), reverse=True)
                best_comp, best_ba, _ = linker_evaluations[0]
                linker_evaluations[0] = (best_comp, best_ba, True)

            for comp, bridge_atoms, keep_linker in linker_evaluations:
                if not keep_linker:
                    keep_atoms = set(comp)

                    changed = True
                    while changed:
                        changed = False
                        to_remove = set()
                        for ka in keep_atoms - bridge_atoms:
                            internal_bonds = sum(1 for nb in local_adj[ka] if nb in keep_atoms)
                            if internal_bonds <= 1:
                                to_remove.add(ka)
                        if to_remove:
                            keep_atoms -= to_remove
                            changed = True

                    atoms_to_cut = comp - keep_atoms
                    partial_to_remove.update(atoms_to_cut)

                    for ka in keep_atoms:
                        removed_neighbor_coords = []
                        for nb in local_adj[ka]:
                            if nb in atoms_to_cut:
                                removed_neighbor_coords.append(unwrapped_coords[nb])
                        if removed_neighbor_coords:
                            bridge_atoms_to_cap.append((ka, removed_neighbor_coords))
        else:
            partial_to_remove = set()
            bridge_atoms_to_cap = []
            for comp in components:
                touching_metals = set()
                bridge_atoms = set()
                for atom_idx in comp:
                    for nb in local_adj[atom_idx]:
                        if nb in core_metals:
                            touching_metals.add(nb)
                            bridge_atoms.add(atom_idx)
                remove_cond = len(touching_metals) < 1
                if remove_cond:
                    atoms_to_cut = comp - bridge_atoms
                    partial_to_remove.update(atoms_to_cut)
                    for ba in bridge_atoms:
                        removed_neighbor_coords = []
                        for nb in local_adj[ba]:
                            if nb in atoms_to_cut:
                                removed_neighbor_coords.append(unwrapped_coords[nb])
                        if removed_neighbor_coords:
                            bridge_atoms_to_cap.append((ba, removed_neighbor_coords))

        if partial_to_remove:
            final_indices -= partial_to_remove
            broken_bonds = [(l, n) for l, n in broken_bonds if l not in partial_to_remove]

        heavy_indices = [idx for idx in final_indices if supercell[idx].species_string != "H"]
        if heavy_indices:
            heavy_adj = {idx: [] for idx in heavy_indices}
            heavy_set = set(heavy_indices)
            for idx in heavy_indices:
                s1 = supercell[idx]
                for n in sc_all_neighbors[idx]:
                    if n.index in heavy_set:
                        nd = np.linalg.norm(n.coords - s1.coords)
                        if self.is_valid_bond(s1.species_string, n.species_string, nd):
                            heavy_adj[idx].append(n.index)

            comps = []
            vis = set()
            for idx in heavy_indices:
                if idx in vis:
                    continue
                q = deque([idx])
                vis.add(idx)
                comp = {idx}
                while q:
                    cur = q.popleft()
                    for nb in heavy_adj[cur]:
                        if nb not in vis:
                            vis.add(nb)
                            comp.add(nb)
                            q.append(nb)
                comps.append(comp)

            if len(comps) > 1:
                comp_core = [set(i for i in comp if i in core_metals) for comp in comps]
                comp_metal_ct = [sum(1 for i in comp if supercell[i].species_string in metals) for comp in comps]

                full_cover = [k for k, cset in enumerate(comp_core) if len(cset) == len(core_metals)]
                if full_cover:
                    best_k = max(full_cover, key=lambda k: (comp_metal_ct[k], -len(comps[k])))
                    keep_heavy = set(comps[best_k])
                else:
                    uncovered = set(core_metals)
                    chosen = []
                    remaining = set(range(len(comps)))
                    while uncovered and remaining:
                        best_k = max(
                            remaining,
                            key=lambda k: (len(comp_core[k] & uncovered), comp_metal_ct[k], -len(comps[k])),
                        )
                        gain = comp_core[best_k] & uncovered
                        if not gain:
                            break
                        chosen.append(best_k)
                        uncovered -= gain
                        remaining.remove(best_k)
                    if chosen:
                        keep_heavy = set().union(*[comps[k] for k in chosen])
                    else:
                        def comp_key(comp):
                            core_count = sum(1 for i in comp if i in core_metals)
                            metal_count = sum(1 for i in comp if supercell[i].species_string in metals)
                            return (core_count, metal_count, len(comp))
                        keep_heavy = max(comps, key=comp_key)
                final_indices = {idx for idx in final_indices if (idx in keep_heavy) or (supercell[idx].species_string == "H")}
                broken_bonds = [(l, n) for l, n in broken_bonds if l in final_indices]

        anchor_idx = sc_center_idx if sc_center_idx in final_indices else next(iter(final_indices))
        old_unwrapped = {idx: np.array(unwrapped_coords[idx]) for idx in final_indices if idx in unwrapped_coords}
        rebuilt_unwrapped = {anchor_idx: np.array(supercell[anchor_idx].coords)}
        q = deque([anchor_idx])
        while q:
            cur = q.popleft()
            cur_site = supercell[cur]
            cur_frac = np.array(cur_site.frac_coords)
            cur_pos = rebuilt_unwrapped[cur]
            for n in sc_all_neighbors[cur]:
                nb = n.index
                if nb not in final_indices or nb in rebuilt_unwrapped:
                    continue
                d = np.linalg.norm(n.coords - cur_site.coords)
                if not self.is_valid_bond(cur_site.species_string, n.species_string, d):
                    continue
                nb_frac = np.array(supercell[nb].frac_coords)
                dfrac = nb_frac - cur_frac
                dfrac -= np.round(dfrac)
                rebuilt_unwrapped[nb] = cur_pos + supercell.lattice.get_cartesian_coords(dfrac)
                q.append(nb)
        for idx in final_indices:
            if idx not in rebuilt_unwrapped:
                rebuilt_unwrapped[idx] = np.array(unwrapped_coords[idx])
        unwrapped_coords.update(rebuilt_unwrapped)

        delta_by_idx = {}
        for idx in final_indices:
            if idx in old_unwrapped:
                delta_by_idx[idx] = unwrapped_coords[idx] - old_unwrapped[idx]
        if delta_by_idx:
            new_broken = []
            for l, npos in broken_bonds:
                shift = delta_by_idx.get(l, np.zeros(3))
                new_broken.append((l, np.array(npos) + shift))
            broken_bonds = new_broken

            shifted_bridge = []
            for ba, removed_coords in bridge_atoms_to_cap:
                shift = delta_by_idx.get(ba, np.zeros(3))
                shifted_bridge.append((ba, [np.array(rc) + shift for rc in removed_coords]))
            bridge_atoms_to_cap = shifted_bridge

        ordered_indices = sorted(final_indices)
        sites = [supercell[idx] for idx in ordered_indices]
        species = [s.species_string for s in sites]
        coords = [unwrapped_coords[idx] for idx in ordered_indices]
        capped_h_flags = [False] * len(species)
        local_index = {sc_idx: i for i, sc_idx in enumerate(ordered_indices)}

        bonds_by_ligand = {}
        for lid, unwrapped_m_pos in broken_bonds:
            if lid not in bonds_by_ligand:
                bonds_by_ligand[lid] = []
            bonds_by_ligand[lid].append(unwrapped_m_pos)

        bridge_cap_map = {}
        for ba_idx, removed_coords in bridge_atoms_to_cap:
            bridge_cap_map.setdefault(ba_idx, []).extend(removed_coords)

        for ba_idx, removed_coords in bridge_cap_map.items():
            ba_site = supercell[ba_idx]
            ba_pos = np.array(unwrapped_coords[ba_idx])
            bl = self.cap_bond_length(ba_site.species_string)
            parent_local_idx = local_index.get(ba_idx)
            if parent_local_idx is None:
                continue

            if ba_site.species_string == "O":
                if self.oxygen_already_protonated(parent_local_idx, species, coords):
                    continue
                avg_vec = np.zeros(3)
                for rc in removed_coords:
                    avg_vec = avg_vec + (rc - ba_pos)
                self.place_capping_h(parent_local_idx, avg_vec, bl, species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)
            else:
                for rc in removed_coords:
                    vec = rc - ba_pos
                    self.place_capping_h(parent_local_idx, vec, bl, species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        already_capped = {ba_idx for ba_idx, _ in bridge_atoms_to_cap}

        for lid, missing_coords in bonds_by_ligand.items():
            if lid in already_capped:
                continue
            lsite = supercell[lid]
            lpos = unwrapped_coords[lid]
            kept = 0
            for n in sc_all_neighbors[lid]:
                if n.index in final_indices:
                    d = np.linalg.norm(n.coords - lsite.coords)
                    if self.is_valid_bond(n.species_string, lsite.species_string, d):
                        kept = kept + 1
            if lsite.species_string == "O" and kept >= 2:
                continue
            avg_vec = np.zeros(3)
            for c in missing_coords:
                avg_vec = avg_vec + (c - lpos)
            if np.linalg.norm(avg_vec) > 0:
                bl = self.cap_bond_length(lsite.species_string)
                parent_local_idx = local_index.get(lid)
                if parent_local_idx is not None:
                    if species[parent_local_idx] == "O" and self.oxygen_already_protonated(parent_local_idx, species, coords):
                        continue
                    self.place_capping_h(parent_local_idx, avg_vec, bl, species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        if not minimize:
            heavy_idx = [i for i, s in enumerate(species) if s != "H"]
            if heavy_idx:
                heavy_adj = {i: [] for i in heavy_idx}
                for a in range(len(heavy_idx)):
                    i = heavy_idx[a]
                    ci = np.array(coords[i])
                    for b in range(a + 1, len(heavy_idx)):
                        j = heavy_idx[b]
                        d = np.linalg.norm(ci - np.array(coords[j]))
                        if self.is_valid_bond(species[i], species[j], d):
                            heavy_adj[i].append(j)
                            heavy_adj[j].append(i)

                comps = []
                vis = set()
                for i in heavy_idx:
                    if i in vis:
                        continue
                    q = deque([i])
                    vis.add(i)
                    comp = {i}
                    while q:
                        cur = q.popleft()
                        for nb in heavy_adj[cur]:
                            if nb not in vis:
                                vis.add(nb)
                                comp.add(nb)
                                q.append(nb)
                    comps.append(comp)

                if len(comps) > 1:
                    metal_counts = [sum(1 for i in comp if species[i] in metals) for comp in comps]
                    max_m = max(metal_counts)
                    if max_m > 0:
                        keep_heavy = set().union(*[comp for comp, m in zip(comps, metal_counts) if m == max_m])
                    else:
                        keep_heavy = max(comps, key=len)
                    keep = set(keep_heavy)
                    for i, s in enumerate(species):
                        if s != "H":
                            continue
                        ci = np.array(coords[i])
                        for j in keep_heavy:
                            d = np.linalg.norm(ci - np.array(coords[j]))
                            if self.is_valid_bond("H", species[j], d):
                                keep.add(i)
                                break
                    species = [s for i, s in enumerate(species) if i in keep]
                    coords = [c for i, c in enumerate(coords) if i in keep]
                    capped_h_flags = [f for i, f in enumerate(capped_h_flags) if i in keep]

        self._last_capped_h_indices = [i for i, f in enumerate(capped_h_flags) if f and species[i] == "H"]
        return species, coords


class COFFragmenter(BaseFragmenter):
    COV_RAD = {
        "H": 0.31,
        "B": 0.84,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "F": 0.57,
        "Si": 1.11,
        "P": 1.07,
        "S": 1.05,
        "Cl": 1.02,
        "Br": 1.20,
        "I": 1.39,
    }

    def __init__(self, radius=6.0, layer_mode="auto"):
        super().__init__(radius=radius)
        allowed = {"auto", "monomer", "dimer"}
        if layer_mode not in allowed:
            raise ValueError(f"layer_mode must be one of {sorted(allowed)}")
        self.layer_mode = layer_mode

    def _rad(self, sym):
        return self.COV_RAD.get(sym, 0.77)

    def is_valid_bond(self, s1, s2, dist):
        if s1 == "H" and s2 == "H":
            return dist < 0.9
        cutoff = 1.25 * (self._rad(s1) + self._rad(s2))
        cutoff = min(2.2, max(1.1, cutoff))
        return dist <= cutoff

    def _cap_open_oxygens(self, species, coords, capped_h_flags):
        heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
        for i in list(heavy_idx):
            if i >= len(species) or species[i] not in {"O", "N", "C"}:
                continue
            sp = species[i]
            pos = np.array(coords[i], dtype=float)
            h_cut = 1.25 if sp in {"N", "O"} else 1.20
            has_h = any(
                spj == "H" and np.linalg.norm(pos - np.array(coords[j], dtype=float)) <= h_cut
                for j, spj in enumerate(species)
            )
            if has_h:
                continue
            if sp == "O" and self.oxygen_already_protonated(i, species, coords):
                continue

            heavy_neighbors = []
            for j, spj in enumerate(species):
                if i == j or spj == "H":
                    continue
                d = np.linalg.norm(pos - np.array(coords[j], dtype=float))
                if self.is_valid_bond(sp, spj, d):
                    heavy_neighbors.append(j)

            if sp == "O":
                if len(heavy_neighbors) != 1 or species[heavy_neighbors[0]] not in {"C", "B", "Si", "P", "S", "N"}:
                    continue
                base = pos - np.array(coords[heavy_neighbors[0]], dtype=float)
            elif sp == "N":
                if len(heavy_neighbors) != 1 or species[heavy_neighbors[0]] not in {"C", "N"}:
                    continue
                base = pos - np.array(coords[heavy_neighbors[0]], dtype=float)
            else:
                # Phenyl/aromatic edge C: two retained heavy neighbors but no H.
                # Avoid carbonyl-like C by requiring C/N neighbors only.
                if len(heavy_neighbors) != 2:
                    continue
                if not all(species[j] in {"C", "N"} for j in heavy_neighbors):
                    continue
                base = np.zeros(3)
                for nb in heavy_neighbors:
                    base -= np.array(coords[nb], dtype=float) - pos
                if np.linalg.norm(base) < 1e-12:
                    base = pos - np.mean([np.array(coords[nb], dtype=float) for nb in heavy_neighbors], axis=0)

            before = len(species)
            self.place_capping_h(i, base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)
            if len(species) == before and np.linalg.norm(base) > 1e-12:
                self.place_capping_h(i, -base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

    def _prune_duplicate_cof_helper_files(self, out_dir, decimals=1):
        seen = {}
        removed = 0
        for old in sorted(Path(out_dir).glob("*.xyz")):
            try:
                key = MOFFragmenter._xyz_unique_key(old, decimals=decimals)
            except Exception:
                key = None
            if key is None:
                continue
            if key in seen:
                old.unlink()
                removed += 1
            else:
                seen[key] = old
        return removed

    def _export_coffragmentor_library(self, result, stem):
        node_dir = Path("cof_nodes_lib")
        linker_dir = Path("cof_linkers_lib")
        node_dir.mkdir(exist_ok=True)
        linker_dir.mkdir(exist_ok=True)
        file_stem = MOFFragmenter._safe_name(stem)

        for collection, out_dir in ((getattr(result, "nodes", []), node_dir), (getattr(result, "linkers", []), linker_dir)):
            self._prune_duplicate_cof_helper_files(out_dir)
            if any(out_dir.glob(f"{file_stem}_*.xyz")):
                continue

            existing_keys = set()
            for old in out_dir.glob("*.xyz"):
                try:
                    key = MOFFragmenter._xyz_unique_key(old, decimals=1)
                except Exception:
                    key = None
                if key is not None:
                    existing_keys.add(key)

            sbus = list(collection) if collection is not None else []
            seen = set(existing_keys)
            out_idx = 0
            for sbu in sbus:
                mol = getattr(sbu, "molecule", None)
                if mol is None:
                    continue
                key = MOFFragmenter._molecule_unique_key(mol, decimals=1)
                if key in seen:
                    continue
                seen.add(key)
                out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"
                while out.exists():
                    out_idx += 1
                    out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"
                mol.to(filename=str(out), fmt="xyz")
                out_idx += 1

    def _export_cof_component_library(self, species, coords, stem, kind):
        out_dir = Path("cof_nodes_lib") if kind == "node" else Path("cof_linkers_lib")
        out_dir.mkdir(exist_ok=True)
        file_stem = MOFFragmenter._safe_name(stem)

        self._prune_duplicate_cof_helper_files(out_dir)
        if any(out_dir.glob(f"{file_stem}_*.xyz")):
            return False

        existing_keys = set()
        for old in out_dir.glob("*.xyz"):
            try:
                key = MOFFragmenter._xyz_unique_key(old, decimals=1)
            except Exception:
                key = None
            if key is not None:
                existing_keys.add(key)

        key = MOFFragmenter._species_coords_unique_key(species, coords, decimals=1)
        if key in existing_keys:
            return False

        out_idx = 0
        out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"
        while out.exists():
            out_idx += 1
            out = out_dir / f"{file_stem}_{out_idx:02d}.xyz"

        Molecule(species, coords).to(filename=str(out), fmt="xyz")
        return True

    def _export_cof6_fallback_library(self, stem, comps, comp_id, center_node, comp_graph, unwrapped, sc_sym):
        try:
            node_comp = sorted(comps[center_node])
            node_sp = [sc_sym[i] for i in node_comp]
            node_co = [np.array(unwrapped[i], dtype=float) for i in node_comp]
            node_written = self._export_cof_component_library(node_sp, node_co, stem, "node")

            linker_written = False
            neigh = sorted(comp_graph.get(center_node, set()))
            for lk in neigh:
                if lk == center_node:
                    continue
                comp = sorted(comps[lk])
                lsp = [sc_sym[i] for i in comp]
                lco = [np.array(unwrapped[i], dtype=float) for i in comp]
                if self._export_cof_component_library(lsp, lco, stem, "linker"):
                    linker_written = True
                    break

            if node_written or linker_written:
                print("  -> COF helper fragments available in cof_nodes_lib/ and cof_linkers_lib/ (Path J6 fallback).")
        except Exception as exc:
            print(f"  COF helper fallback export failed: {exc}")

    def _try_export_coffragmentor_library(self, cif_path):
        try:
            from coffragmentor import COF
        except Exception as exc:
            print(f"  coffragmentor unavailable: {exc}")
            return False
        try:
            result = COF.from_cif(cif_path).fragment()
        except Exception as exc:
            print(f"  coffragmentor library export failed: {exc}")
            return False

        nodes = getattr(result, "nodes", [])
        linkers = getattr(result, "linkers", [])
        if not nodes or not linkers:
            print("  coffragmentor found no node/linker set to export.")
            return False

        stem = Path(cif_path).stem
        self._export_coffragmentor_library(result, stem)
        print("  -> COF helper fragments available in cof_nodes_lib/ and cof_linkers_lib/.")
        return True

    def _try_cof_graph_node_linker_fragment(self, cif_path, output_path, minimize=False):
        try:
            struct = Structure.from_file(cif_path)
        except Exception:
            return None

        dims = [max(1, int(np.ceil(28.0 / a))) for a in struct.lattice.abc]
        dims = [max(3, d) if a < 15.0 else d for d, a in zip(dims, struct.lattice.abc)]
        supercell = struct * dims

        def _site_symbol(site):
            try:
                return site.specie.symbol
            except Exception:
                return max(site.species.items(), key=lambda kv: kv[1])[0].symbol

        sc_sym = [_site_symbol(site) for site in supercell]
        all_neigh = supercell.get_all_neighbors(r=2.4)
        graph = [[] for _ in range(len(supercell))]
        for i, neighs in enumerate(all_neigh):
            si = sc_sym[i]
            ci = supercell[i].coords
            for n in neighs:
                j = n.index
                if j <= i:
                    continue
                d = np.linalg.norm(ci - n.coords)
                if self.is_valid_bond(si, sc_sym[j], d):
                    graph[i].append(j)
                    graph[j].append(i)

        cut_edges = []
        for i in range(len(supercell)):
            si = sc_sym[i]
            for j in graph[i]:
                if j <= i:
                    continue
                if {si, sc_sym[j]} == {"B", "O"}:
                    cut_edges.append((i, j))
        if not cut_edges:
            return None

        g2 = [set(nbs) for nbs in graph]
        for i, j in cut_edges:
            g2[i].discard(j)
            g2[j].discard(i)

        comp_id = [-1] * len(supercell)
        comps = []
        cid = 0
        for i in range(len(supercell)):
            if comp_id[i] != -1:
                continue
            q = deque([i])
            comp_id[i] = cid
            comp = {i}
            while q:
                u = q.popleft()
                for v in g2[u]:
                    if comp_id[v] == -1:
                        comp_id[v] = cid
                        comp.add(v)
                        q.append(v)
            comps.append(comp)
            cid += 1

        node_candidates = []
        for k, comp in enumerate(comps):
            ccount = sum(1 for x in comp if sc_sym[x] == "C")
            bcount = sum(1 for x in comp if sc_sym[x] == "B")
            if bcount >= 1 and ccount >= 6:
                node_candidates.append(k)
        if not node_candidates:
            return None

        comp_graph = {k: set() for k in range(len(comps))}
        for i, j in cut_edges:
            ci, cj = comp_id[i], comp_id[j]
            if ci != cj:
                comp_graph[ci].add(cj)
                comp_graph[cj].add(ci)

        ctr = supercell.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
        def comp_center_cart(k):
            return np.mean([supercell[i].coords for i in comps[k]], axis=0)
        def comp_center_frac_mod1(k):
            fr = np.mean([np.array(supercell[i].frac_coords, dtype=float) for i in comps[k]], axis=0)
            return np.mod(fr, 1.0)

        center_node = min(node_candidates, key=lambda k: float(np.linalg.norm(comp_center_cart(k) - ctr)))

        # Topology split: keep boroxine-like nodes for the stricter path below,
        # and use this graph fallback for other COF topologies (e.g., COF-66).
        cc0 = sum(1 for x in comps[center_node] if sc_sym[x] == "C")
        bb0 = sum(1 for x in comps[center_node] if sc_sym[x] == "B")
        oo0 = sum(1 for x in comps[center_node] if sc_sym[x] == "O")
        is_boroxine_like = (bb0 >= 2 and 4 <= cc0 <= 12 and oo0 <= 2)
        if is_boroxine_like:
            return None

        def linker_neighbors(node_comp_idx):
            return [k for k in comp_graph.get(node_comp_idx, set()) if k not in node_candidates]

        def lk_key(k):
            ccount = sum(1 for x in comps[k] if sc_sym[x] == "C")
            ocount = sum(1 for x in comps[k] if sc_sym[x] == "O")
            bcount = sum(1 for x in comps[k] if sc_sym[x] == "B")
            return (ocount, -bcount, ccount, len(comps[k]))

        center_linkers = linker_neighbors(center_node)
        if not center_linkers:
            return None

        layer_mode = getattr(self, "layer_mode", "auto")
        selected_nodes = [center_node]
        partner_vec = None
        if layer_mode != "monomer" and len(node_candidates) > 1:
            # choose second layer along shortest lattice axis using minimum-image
            # fractional center difference, to avoid in-plane partner mistakes.
            axis = int(np.argmin(np.array(struct.lattice.abc, dtype=float)))
            f0 = comp_center_frac_mod1(center_node)
            best = None
            best_key = None
            center_sig = (
                sum(1 for x in comps[center_node] if sc_sym[x] == "B"),
                sum(1 for x in comps[center_node] if sc_sym[x] == "C"),
            )
            for k in node_candidates:
                if k == center_node:
                    continue
                sig_k = (
                    sum(1 for x in comps[k] if sc_sym[x] == "B"),
                    sum(1 for x in comps[k] if sc_sym[x] == "C"),
                )
                if sig_k != center_sig:
                    continue
                fk = comp_center_frac_mod1(k)
                dfrac = fk - f0
                dfrac -= np.round(dfrac)
                dz = abs(float(dfrac[axis]))
                if dz < 1e-3:
                    continue
                dperp = float(np.linalg.norm(np.delete(dfrac, axis)))
                # Prefer same in-plane phase first, then nearest adjacent
                # layer along stacking axis (smallest non-zero dz).
                nlk = len([x for x in comp_graph.get(k, set()) if x not in node_candidates])
                cdist = float(np.linalg.norm(comp_center_cart(k) - ctr))
                key = (dperp, dz, cdist, -nlk)
                if best_key is None or key < best_key:
                    best_key = key
                    best = k
            if best is not None:
                selected_nodes.append(best)
                f_best = comp_center_frac_mod1(best)
                f_cent = comp_center_frac_mod1(center_node)
                dfrac = f_best - f_cent
                dfrac -= np.round(dfrac)
                partner_vec = np.array(supercell.lattice.get_cartesian_coords(dfrac), dtype=float)

        # For normal dimer in graph fallback, build one complete monomer
        # (center node + all attached linkers), then place second monomer by
        # crystal layer vector. This avoids asymmetric second-layer growth.
        growth_nodes = selected_nodes
        if layer_mode != "monomer" and partner_vec is not None:
            growth_nodes = [center_node]

        keep_comps = set(growth_nodes)
        selected_linkers = []
        for nidx in growth_nodes:
            n_linkers = linker_neighbors(nidx)
            if not n_linkers:
                continue
            if minimize:
                picked = max(n_linkers, key=lk_key)
                selected_linkers.append(picked)
                keep_comps.add(picked)
            else:
                selected_linkers.extend(sorted(n_linkers))
                keep_comps |= set(n_linkers)

        if not selected_linkers:
            return None

        keep_atoms = set()
        for k in keep_comps:
            keep_atoms |= comps[k]

        # unwrap coordinates on full graph, then align kept components via cut edges
        sc_center_idx = min(comps[center_node], key=lambda i: np.linalg.norm(supercell[i].coords - ctr))
        unwrapped = [None] * len(supercell)
        unwrapped[sc_center_idx] = np.array(supercell[sc_center_idx].coords, dtype=float)
        q = deque([sc_center_idx])
        while q:
            i = q.popleft()
            ui = unwrapped[i]
            fi = np.array(supercell[i].frac_coords, dtype=float)
            for j in graph[i]:
                if unwrapped[j] is not None:
                    continue
                fj = np.array(supercell[j].frac_coords, dtype=float)
                dfrac = fj - fi
                dfrac -= np.round(dfrac)
                unwrapped[j] = ui + supercell.lattice.get_cartesian_coords(dfrac)
                q.append(j)
        for i in range(len(unwrapped)):
            if unwrapped[i] is None:
                unwrapped[i] = np.array(supercell[i].coords, dtype=float)

        lattice = np.array(supercell.lattice.matrix, dtype=float)
        comp_shift = {k: np.zeros(3) for k in keep_comps}
        comp_adj = {k: [] for k in keep_comps}
        for i, j in cut_edges:
            ci, cj = comp_id[i], comp_id[j]
            if ci in keep_comps and cj in keep_comps and ci != cj:
                comp_adj[ci].append((cj, i, j))
                comp_adj[cj].append((ci, j, i))

        seen_comp = set()
        # Align each selected node-layer subgraph from its own root so linker
        # images are correct for both dimer layers.
        root_order = [r for r in selected_nodes if r in keep_comps]
        if center_node in keep_comps and center_node not in root_order:
            root_order.insert(0, center_node)
        for root in root_order:
            if root in seen_comp:
                continue
            qcomp = deque([root])
            seen_comp.add(root)
            while qcomp:
                ca = qcomp.popleft()
                shift_a = comp_shift[ca]
                for cb, ia, ib in comp_adj.get(ca, []):
                    if cb in seen_comp:
                        continue
                    pa = np.array(unwrapped[ia], dtype=float) + shift_a
                    pb0 = np.array(unwrapped[ib], dtype=float)
                    best = None
                    bestd = float('inf')
                    for aa in (-1, 0, 1):
                        for bb in (-1, 0, 1):
                            for cc in (-1, 0, 1):
                                sh = aa * lattice[0] + bb * lattice[1] + cc * lattice[2]
                                d = float(np.linalg.norm(pa - (pb0 + sh)))
                                if d < bestd:
                                    bestd = d
                                    best = sh
                    comp_shift[cb] = best
                    seen_comp.add(cb)
                    qcomp.append(cb)

        for k in keep_comps:
            sh = comp_shift.get(k)
            if sh is None:
                continue
            for gi in comps[k]:
                if gi in keep_atoms:
                    unwrapped[gi] = np.array(unwrapped[gi], dtype=float) + sh

        # export helper node/linker from fallback decomposition
        stem = Path(cif_path).stem
        node_atoms = sorted(comps[center_node])
        node_sp = [sc_sym[i] for i in node_atoms]
        node_co = [np.array(unwrapped[i], dtype=float) for i in node_atoms]
        self._export_cof_component_library(node_sp, node_co, stem, "node")
        lk = sorted(selected_linkers)[0]
        lk_atoms = sorted(comps[lk])
        lk_sp = [sc_sym[i] for i in lk_atoms]
        lk_co = [np.array(unwrapped[i], dtype=float) for i in lk_atoms]
        self._export_cof_component_library(lk_sp, lk_co, stem, "linker")

        ordered = sorted(keep_atoms)
        local = {gi: li for li, gi in enumerate(ordered)}
        species = [sc_sym[i] for i in ordered]
        coords = [np.array(unwrapped[i], dtype=float) for i in ordered]
        capped_h_flags = [False] * len(species)

        # cap boundary cuts to excluded components
        for u in ordered:
            li = local[u]
            su = species[li]
            if su == "H":
                continue
            for v in graph[u]:
                if v in keep_atoms:
                    continue
                vec = np.array(unwrapped[v], dtype=float) - np.array(unwrapped[u], dtype=float)
                self.place_capping_h(li, vec, self.cap_bond_length(su), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        # Minimized mode: ensure edge O/B stay and are capped if open.
        if minimize:
            heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
            hadj = {i: [] for i in heavy_idx}
            for a in range(len(heavy_idx)):
                i = heavy_idx[a]
                ci = np.array(coords[i], dtype=float)
                for b in range(a + 1, len(heavy_idx)):
                    j = heavy_idx[b]
                    d = float(np.linalg.norm(ci - np.array(coords[j], dtype=float)))
                    if self.is_valid_bond(species[i], species[j], d):
                        hadj[i].append(j)
                        hadj[j].append(i)
            target_valence = {"O": 2, "B": 3}
            for i in heavy_idx:
                sp = species[i]
                if sp not in target_valence:
                    continue
                if sp == "O" and self.oxygen_already_protonated(i, species, coords):
                    continue
                deficit = max(0, target_valence[sp] - len(hadj.get(i, [])))
                if sp == "B":
                    deficit = min(deficit, 2)
                if deficit <= 0:
                    continue
                base = np.zeros(3)
                if hadj.get(i):
                    c0 = np.array(coords[i], dtype=float)
                    for nb in hadj[i]:
                        base -= (np.array(coords[nb], dtype=float) - c0)
                else:
                    base = np.array([1.0, 0.0, 0.0])
                for _ in range(deficit):
                    self.place_capping_h(i, base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        # Enforce boron cap limit in this fallback path: at most 2 bonded H per B.
        remove_h = set()
        for i, sp in enumerate(species):
            if sp != "B":
                continue
            bpos = np.array(coords[i], dtype=float)
            h_nbs = []
            for j, sj in enumerate(species):
                if sj != "H" or j in remove_h:
                    continue
                d = float(np.linalg.norm(bpos - np.array(coords[j], dtype=float)))
                if self.is_valid_bond("B", "H", d):
                    h_nbs.append((d, j))
            if len(h_nbs) > 2:
                h_nbs.sort(key=lambda x: x[0])
                for _, hj in h_nbs[2:]:
                    remove_h.add(hj)

        if remove_h:
            keep_idx = [i for i in range(len(species)) if i not in remove_h]
            species = [species[i] for i in keep_idx]
            coords = [coords[i] for i in keep_idx]
            capped_h_flags = [capped_h_flags[i] for i in keep_idx]

        capped_h_indices = [i for i, f in enumerate(capped_h_flags) if f and species[i] == "H"]
        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)
        _label = str(output_path or 'cof_fragment')
        species, coords, capped_h_indices = self.fix_odd_electron_multiplicity(species, coords, capped_h_indices, label=_label)

        # Normal dimer: preserve a complete monomer then place second copy at
        # the actual crystal layer vector to avoid one-sided linker loss.
        if layer_mode != "monomer" and partner_vec is not None:
            base_species = list(species)
            base_coords = [np.array(c, dtype=float) for c in coords]
            species = base_species + base_species
            coords = base_coords + [c + partner_vec for c in base_coords]

        # Keep principal connected layers only; remove stray far islands.
        if species:
            adj = [[] for _ in range(len(species))]
            for i in range(len(species)):
                ci = np.array(coords[i], dtype=float)
                for j in range(i + 1, len(species)):
                    d = float(np.linalg.norm(ci - np.array(coords[j], dtype=float)))
                    if self.is_valid_bond(species[i], species[j], d):
                        adj[i].append(j)
                        adj[j].append(i)
            vis = set()
            comps_local = []
            for i in range(len(species)):
                if i in vis:
                    continue
                q = deque([i])
                vis.add(i)
                comp = {i}
                while q:
                    u = q.popleft()
                    for v in adj[u]:
                        if v not in vis:
                            vis.add(v)
                            comp.add(v)
                            q.append(v)
                comps_local.append(comp)
            if len(comps_local) > 1:
                def ckey(comp):
                    node_ct = sum(1 for idx in comp if species[idx] in {"B", "O"})
                    return (node_ct, len(comp))
                ranked = sorted(comps_local, key=ckey, reverse=True)
                keep_n = 1 if layer_mode == "monomer" else 2
                keep = set().union(*ranked[:keep_n])
                species = [x for i, x in enumerate(species) if i in keep]
                coords = [x for i, x in enumerate(coords) if i in keep]

        if output_path:
            Molecule(species, coords).to(filename=output_path, fmt="xyz")
            print(f"  -> COF Path J (graph node+linker fallback combine).")
            print("  -> Helper fragments available in cof_nodes_lib/ and cof_linkers_lib/.")
            print(f"Final size: {len(species)} atoms")
            print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)

    def _try_coffragmentor_node_linker_fragment(self, cif_path, output_path, minimize=False):
        try:
            from coffragmentor import COF
        except Exception as exc:
            print(f"  coffragmentor unavailable: {exc}")
            return None
        try:
            struct = Structure.from_file(cif_path)
            result = COF.from_cif(cif_path).fragment()
        except Exception as exc:
            print(f"  coffragmentor node+linker failed: {exc}")
            return None

        nodes = list(getattr(result, "nodes", []))
        linkers = list(getattr(result, "linkers", []))
        if not nodes or not linkers:
            return None

        self._export_coffragmentor_library(result, Path(cif_path).stem)

        center = struct.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
        node = min(nodes, key=lambda sbu: float(np.linalg.norm(np.mean(sbu.molecule.cart_coords, axis=0) - center)))
        node_sp = [str(x) for x in node.molecule.species]
        node_co = [np.array(c, dtype=float) for c in node.molecule.cart_coords]
        node_ctr = np.mean(node_co, axis=0)

        image_vectors = []
        for ia in (-1, 0, 1):
            for ib in (-1, 0, 1):
                for ic in (-1, 0, 1):
                    image_vectors.append(
                        ia * np.array(struct.lattice.matrix[0], dtype=float)
                        + ib * np.array(struct.lattice.matrix[1], dtype=float)
                        + ic * np.array(struct.lattice.matrix[2], dtype=float)
                    )

        def linker_image_score(linker, image_shift):
            lsp = [str(x) for x in linker.molecule.species]
            lco = [np.array(c, dtype=float) + image_shift for c in linker.molecule.cart_coords]
            attach_bonds = 0
            min_d = float("inf")
            for i, si in enumerate(node_sp):
                if si == "H":
                    continue
                for j, sj in enumerate(lsp):
                    if sj == "H":
                        continue
                    d = float(np.linalg.norm(node_co[i] - lco[j]))
                    min_d = min(min_d, d)
                    if self.is_valid_bond(si, sj, d):
                        attach_bonds += 1
            lctr = np.mean(lco, axis=0)
            return (attach_bonds, -min_d, -float(np.linalg.norm(lctr - node_ctr)))

        scored_images = []
        for linker in linkers:
            for image_shift in image_vectors:
                score = linker_image_score(linker, image_shift)
                if score[0] > 0:
                    scored_images.append((score, linker, image_shift))
        if not scored_images:
            return None
        scored_images.sort(key=lambda x: x[0], reverse=True)
        selected_images = scored_images[:1] if minimize else scored_images

        def append_merged(species, coords, add_species, add_coords, tol=0.08):
            for sp, coord in zip(add_species, add_coords):
                c = np.array(coord, dtype=float)
                duplicate = False
                for i, old_sp in enumerate(species):
                    if old_sp != sp:
                        continue
                    if np.linalg.norm(np.array(coords[i], dtype=float) - c) <= tol:
                        duplicate = True
                        break
                if not duplicate:
                    species.append(sp)
                    coords.append(c)

        species = []
        coords = []
        append_merged(species, coords, node_sp, node_co)
        for _, linker, image_shift in selected_images:
            lsp = [str(x) for x in linker.molecule.species]
            lco = [np.array(c, dtype=float) + image_shift for c in linker.molecule.cart_coords]
            append_merged(species, coords, lsp, lco)

        capped_h_flags = [False] * len(species)

        # COF Path J edge completion for minimized mode:
        # keep edge O/B atoms and cap open valences by H.
        if minimize and species:
            heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
            hadj = {i: [] for i in heavy_idx}
            for a in range(len(heavy_idx)):
                i = heavy_idx[a]
                ci = np.array(coords[i], dtype=float)
                for b in range(a + 1, len(heavy_idx)):
                    j = heavy_idx[b]
                    d = float(np.linalg.norm(ci - np.array(coords[j], dtype=float)))
                    if self.is_valid_bond(species[i], species[j], d):
                        hadj[i].append(j)
                        hadj[j].append(i)

            target_valence = {"O": 2, "B": 3}
            for i in heavy_idx:
                sp = species[i]
                if sp not in target_valence:
                    continue
                if sp == "O" and self.oxygen_already_protonated(i, species, coords):
                    continue
                cur_deg = len(hadj.get(i, []))
                deficit = max(0, target_valence[sp] - cur_deg)
                if deficit <= 0:
                    continue

                base = np.zeros(3)
                if hadj.get(i):
                    c0 = np.array(coords[i], dtype=float)
                    for nb in hadj[i]:
                        base -= (np.array(coords[nb], dtype=float) - c0)
                else:
                    base = np.array([1.0, 0.0, 0.0])

                for _ in range(deficit):
                    self.place_capping_h(i, base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        # Keep helper heavy atoms fixed; only adjust capped H atoms.
        self._cap_open_oxygens(species, coords, capped_h_flags)
        capped_h_indices = [i for i, is_cap in enumerate(capped_h_flags) if is_cap and species[i] == "H"]
        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)
        _label = str(output_path or 'cof_fragment')
        species, coords, capped_h_indices = self.fix_odd_electron_multiplicity(species, coords, capped_h_indices, label=_label)

        # Global layered-COF rule for Path J: if a face-to-face layer spacing is
        # present and dimer output is requested, duplicate the completed fragment
        # by the nearest stacking vector from the crystal.
        layer_mode = getattr(self, "layer_mode", "auto")
        abc = np.array(struct.lattice.abc, dtype=float)
        stack_axis = int(np.argmin(abc))
        stack_len = float(abc[stack_axis])
        if layer_mode != "monomer" and 2.5 <= stack_len <= 5.0 and len(species) > 0:
            layer_vec = np.array(struct.lattice.matrix[stack_axis], dtype=float)
            base_species = list(species)
            base_coords = [np.array(c, dtype=float) for c in coords]
            base_flags = list(capped_h_flags)
            species = base_species + base_species
            coords = base_coords + [c + layer_vec for c in base_coords]
            capped_h_flags = base_flags + base_flags

        if output_path:
            mol = Molecule(species, coords)
            mol.to(filename=output_path, fmt="xyz")
            print("  -> COF Path J (coffragmentor node+linker combine).")
            print("  -> Helper fragments available in cof_nodes_lib/ and cof_linkers_lib/.")
            print(f"Final size: {len(species)} atoms")
            print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)

    def extract(self, cif_path, output_path="cof_fragment.xyz", center_idx=-1, minimize=False):
        print(f"Loading '{cif_path}'...")
        try:
            struct = Structure.from_file(cif_path)
        except Exception as exc:
            print(f"  Standard CIF load failed: {exc}")
            print("  Retrying with tolerant CIF parser...")
            parser = CifParser(cif_path, occupancy_tolerance=2.1)
            structs = parser.get_structures(primitive=False)
            if not structs:
                raise
            struct = structs[0]
        combined = self._try_coffragmentor_node_linker_fragment(cif_path, output_path, minimize=minimize)
        if combined is not None:
            return combined
        combined = self._try_cof_graph_node_linker_fragment(cif_path, output_path, minimize=minimize)
        if combined is not None:
            return combined
        print("Creating supercell...")
        dims = [max(1, int(np.ceil(28.0 / a))) for a in struct.lattice.abc]
        dims = [max(3, d) if a < 15.0 else d for d, a in zip(dims, struct.lattice.abc)]
        supercell = struct * dims

        def _site_symbol(site):
            try:
                return site.specie.symbol
            except Exception:
                # Disordered site: choose highest-occupancy species.
                return max(site.species.items(), key=lambda kv: kv[1])[0].symbol

        sc_sym = [_site_symbol(site) for site in supercell]

        print("Building bond graph...")
        all_neigh = supercell.get_all_neighbors(r=2.4)
        graph = [[] for _ in range(len(supercell))]
        for i, neighs in enumerate(all_neigh):
            si = sc_sym[i]
            ci = supercell[i].coords
            for n in neighs:
                j = n.index
                if j <= i:
                    continue
                d = np.linalg.norm(ci - n.coords)
                if self.is_valid_bond(si, sc_sym[j], d):
                    graph[i].append(j)
                    graph[j].append(i)

        structure_stem = Path(cif_path).stem

        # Boroxine-like node+linker strict path (previously COF-6-specific):
        # detect building blocks by cutting B-O bridges, then assemble node +
        # attached linkers from topology similarity (not filename).
        if True:
            cut_edges = []
            for i in range(len(supercell)):
                si = sc_sym[i]
                for j in graph[i]:
                    if j <= i:
                        continue
                    sj = sc_sym[j]
                    if {si, sj} == {"B", "O"}:
                        cut_edges.append((i, j))

            if cut_edges:
                g2 = [set(nbs) for nbs in graph]
                for i, j in cut_edges:
                    if j in g2[i]:
                        g2[i].remove(j)
                    if i in g2[j]:
                        g2[j].remove(i)

                comp_id = [-1] * len(supercell)
                comps = []
                cid = 0
                for i in range(len(supercell)):
                    if comp_id[i] != -1:
                        continue
                    q = deque([i])
                    comp_id[i] = cid
                    comp = {i}
                    while q:
                        u = q.popleft()
                        for v in g2[u]:
                            if comp_id[v] == -1:
                                comp_id[v] = cid
                                comp.add(v)
                                q.append(v)
                    comps.append(comp)
                    cid += 1

                # classify components; swap node/linker assignment for COF-6
                # per project convention: node is B-containing unit.
                node_candidates = []
                for k, comp in enumerate(comps):
                    ccount = sum(1 for x in comp if sc_sym[x] == "C")
                    bcount = sum(1 for x in comp if sc_sym[x] == "B")
                    ocount = sum(1 for x in comp if sc_sym[x] == "O")
                    if bcount >= 1 and ccount >= 6:
                        node_candidates.append(k)

                if node_candidates:
                    ctr = supercell.lattice.get_cartesian_coords([0.5, 0.5, 0.5])

                    # Topology-similarity gate: require a boroxine-like node
                    # signature (B-rich aromatic node) before entering this
                    # strict node+linker path.
                    def _node_sig(comp):
                        ccount = sum(1 for x in comp if sc_sym[x] == "C")
                        bcount = sum(1 for x in comp if sc_sym[x] == "B")
                        ocount = sum(1 for x in comp if sc_sym[x] == "O")
                        return ccount, bcount, ocount
                    center_probe = min(node_candidates, key=lambda k: float(np.linalg.norm(np.mean([supercell[i].coords for i in comps[k]], axis=0) - ctr)))
                    cc, bb, oo = _node_sig(comps[center_probe])
                    is_boroxine_like = (bb >= 2 and 4 <= cc <= 12 and oo <= 2)
                    if not is_boroxine_like:
                        pass
                    else:
                        def comp_center(k):
                            return np.mean([supercell[i].coords for i in comps[k]], axis=0)

                        # Build component adjacency from cut edges.
                        comp_graph = {k: set() for k in range(len(comps))}
                        for i, j in cut_edges:
                            ci, cj = comp_id[i], comp_id[j]
                            if ci != cj:
                                comp_graph[ci].add(cj)
                                comp_graph[cj].add(ci)

                        def flood_keep(node_set):
                            forbidden_nodes = set(node_candidates) - set(node_set)
                            keep = set(node_set)
                            qcg = deque(node_set)
                            while qcg:
                                ca = qcg.popleft()
                                for cb in comp_graph.get(ca, set()):
                                    if cb in keep or cb in forbidden_nodes:
                                        continue
                                    keep.add(cb)
                                    qcg.append(cb)
                            return keep

                        center_node = min(node_candidates, key=lambda k: float(np.linalg.norm(comp_center(k) - ctr)))
                        # Build one chemically clean monomer first; for dimer mode
                        # we duplicate this validated monomer with a layer shift.
                        do_stack_dimer = (getattr(self, "layer_mode", "auto") != "monomer")
                        chosen_nodes = [center_node]

                        node_set = set(chosen_nodes)
                        if minimize:
                            # COF-6 minimum: keep center node + one O-rich linker
                            # component only (no extra terminal B-rich components).
                            neighbor_candidates = [k for k in comp_graph.get(center_node, set()) if k not in node_candidates]
                            if not neighbor_candidates:
                                neighbor_candidates = [k for k in comp_graph.get(center_node, set()) if k != center_node]

                            branch_start = None
                            if neighbor_candidates:
                                def lk_key(k):
                                    ccount = sum(1 for x in comps[k] if sc_sym[x] == "C")
                                    bcount = sum(1 for x in comps[k] if sc_sym[x] == "B")
                                    ocount = sum(1 for x in comps[k] if sc_sym[x] == "O")
                                    return (ocount, -bcount, ccount, len(comps[k]))
                                branch_start = max(neighbor_candidates, key=lk_key)

                            keep_comps = set(node_set)
                            if branch_start is not None:
                                keep_comps.add(branch_start)
                        else:
                            keep_comps = flood_keep(node_set)
                        keep_atoms = set()
                        for k in keep_comps:
                            keep_atoms |= comps[k]

                        # Build coordinates with unwrapping (existing routine below)
                        sc_center_idx = min(comps[center_node], key=lambda i: np.linalg.norm(supercell[i].coords - ctr))
                        unwrapped = [None] * len(supercell)
                        unwrapped[sc_center_idx] = np.array(supercell[sc_center_idx].coords)
                        q = deque([sc_center_idx])
                        while q:
                            i = q.popleft()
                            ui = unwrapped[i]
                            fi = np.array(supercell[i].frac_coords)
                            for j in graph[i]:
                                if unwrapped[j] is not None:
                                    continue
                                fj = np.array(supercell[j].frac_coords)
                                dfrac = fj - fi
                                dfrac -= np.round(dfrac)
                                unwrapped[j] = ui + supercell.lattice.get_cartesian_coords(dfrac)
                                q.append(j)
                        for i in range(len(unwrapped)):
                            if unwrapped[i] is None:
                                unwrapped[i] = np.array(supercell[i].coords)

                        # Align kept cut-components to nearest periodic images to
                        # avoid wrong wrapping / apparent dangling atoms in dimer.
                        lattice = np.array(supercell.lattice.matrix, dtype=float)
                        comp_shift = {k: np.zeros(3) for k in keep_comps}
                        node_root = center_node
                        comp_adj = {k: [] for k in keep_comps}
                        for i, j in cut_edges:
                            ci, cj = comp_id[i], comp_id[j]
                            if ci in keep_comps and cj in keep_comps and ci != cj:
                                comp_adj[ci].append((cj, i, j))
                                comp_adj[cj].append((ci, j, i))

                        qcomp = deque([node_root])
                        seen_comp = {node_root}
                        while qcomp:
                            ca = qcomp.popleft()
                            shift_a = comp_shift[ca]
                            for cb, ia, ib in comp_adj.get(ca, []):
                                if cb in seen_comp:
                                    continue
                                pa = np.array(unwrapped[ia], dtype=float) + shift_a
                                pb0 = np.array(unwrapped[ib], dtype=float)
                                best = None
                                bestd = float('inf')
                                for aa in (-1, 0, 1):
                                    for bb in (-1, 0, 1):
                                        for cc in (-1, 0, 1):
                                            sh = aa * lattice[0] + bb * lattice[1] + cc * lattice[2]
                                            d = float(np.linalg.norm(pa - (pb0 + sh)))
                                            if d < bestd:
                                                bestd = d
                                                best = sh
                                comp_shift[cb] = best
                                seen_comp.add(cb)
                                qcomp.append(cb)

                        # Apply component shifts to kept atoms.
                        for k in keep_comps:
                            sh = comp_shift.get(k)
                            if sh is None:
                                continue
                            for gi in comps[k]:
                                if gi in keep_atoms:
                                    unwrapped[gi] = np.array(unwrapped[gi], dtype=float) + sh

                        self._export_cof6_fallback_library(structure_stem, comps, comp_id, center_node, comp_graph, unwrapped, sc_sym)

                        ordered = sorted(keep_atoms)
                        local = {gi: li for li, gi in enumerate(ordered)}
                        species = [sc_sym[i] for i in ordered]
                        coords = [np.array(unwrapped[i]) for i in ordered]
                        capped_h_flags = [False] * len(species)

                        # cap edges cut to excluded components
                        for u in ordered:
                            li = local[u]
                            su = species[li]
                            if su == "H":
                                continue
                            for v in graph[u]:
                                if v in keep_atoms:
                                    continue
                                vec = np.array(unwrapped[v]) - np.array(unwrapped[u])
                                self.place_capping_h(li, vec, self.cap_bond_length(su), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

                        capped_h_indices = [i for i, f in enumerate(capped_h_flags) if f and species[i] == "H"]
                        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)
                        _label = str(output_path or 'cof_fragment')
                        species, coords, capped_h_indices = self.fix_odd_electron_multiplicity(species, coords, capped_h_indices, label=_label)

                        if do_stack_dimer:
                            # Use geometric layer shift from nearest node component
                            # (minimum-image), not shortest lattice vector.
                            stack_vec = None
                            others = [k for k in node_candidates if k != center_node]
                            if others:
                                c0 = np.mean([supercell[i].coords for i in comps[center_node]], axis=0)
                                fi0 = np.mean([supercell[i].frac_coords for i in comps[center_node]], axis=0)
                                best_d = float('inf')
                                best_v = None
                                for k in others:
                                    ck = np.mean([supercell[i].coords for i in comps[k]], axis=0)
                                    fik = np.mean([supercell[i].frac_coords for i in comps[k]], axis=0)
                                    dfrac = np.array(fik) - np.array(fi0)
                                    dfrac -= np.round(dfrac)
                                    v = np.array(supercell.lattice.get_cartesian_coords(dfrac), dtype=float)
                                    d = float(np.linalg.norm(v))
                                    if d < best_d and d > 1e-6:
                                        best_d = d
                                        best_v = v
                                # Expected interlayer spacing for COF-6 is around 3.6 A.
                                if best_v is not None and 2.0 <= best_d <= 6.0:
                                    stack_vec = best_v

                            if stack_vec is None:
                                lat = np.array(supercell.lattice.matrix, dtype=float)
                                stack_vec = min(lat, key=lambda v: float(np.linalg.norm(v)))

                            base_species = list(species)
                            base_coords = [np.array(c, dtype=float) for c in coords]
                            species = base_species + base_species
                            coords = base_coords + [c + stack_vec for c in base_coords]

                        print("  -> COF Path J6 (boroxine-like node+all attached linkers).")
                        print(f"Final size: {len(species)} atoms")
                        mol = Molecule(species, coords)
                        if output_path:
                            mol.to(filename=output_path, fmt="xyz")
                            print(f"Saved: {output_path}")
                        return FragmentResult(species=species, coords=coords)

        bo_node_species = {"B", "O"}
        bo_node_atoms = {i for i in range(len(supercell)) if sc_sym[i] in bo_node_species}

        ctr = supercell.lattice.get_cartesian_coords([0.5, 0.5, 0.5])

        path_mode = "A"
        node_atoms = set()
        core_nodes = set()
        node_species = set()

        # Path H (COF-105-like only): large-cell Si-node SBU
        # SBU = Si + 4 phenyl; linker branch is B-ended 3-edge side.
        si_nodes = [i for i, el in enumerate(sc_sym) if el == "Si"]
        if si_nodes and max(struct.lattice.abc) > 40.0:
            h_cores = []
            for si in si_nodes:
                core = {si}
                qh = deque([(si, 0)])
                seen_h = {si}
                while qh:
                    u, dep = qh.popleft()
                    if dep >= 5:
                        continue
                    for v in graph[u]:
                        if v in seen_h or sc_sym[v] == "H":
                            continue
                        sv = sc_sym[v]
                        if sv in {"B", "O"}:
                            continue
                        if sv not in {"Si", "C"}:
                            continue
                        seen_h.add(v)
                        core.add(v)
                        qh.append((v, dep + 1))
                ccount = sum(1 for i in core if sc_sym[i] == "C")
                if ccount >= 12:
                    h_cores.append(core)

            if h_cores:
                def _ctr(core):
                    return np.mean([supercell[i].coords for i in core], axis=0)
                center_core = min(h_cores, key=lambda c: np.linalg.norm(_ctr(c) - ctr))
                core_nodes = set(center_core)
                node_atoms = set().union(*h_cores)
                node_species = {"Si"}
                if center_idx >= 0:
                    sc_center_idx = center_idx
                else:
                    si_in_core = [i for i in center_core if sc_sym[i] == "Si"]
                    pool = si_in_core if si_in_core else list(center_core)
                    sc_center_idx = min(pool, key=lambda i: np.linalg.norm(supercell[i].coords - ctr))
                path_mode = "H"
                print(f"  -> COF Path H (COF-105 Si-node SBU). cores: {len(h_cores)}, core size: {len(center_core)}")

        if path_mode == "A" and bo_node_atoms:
            node_species = bo_node_species
            node_atoms = set(bo_node_atoms)

            if center_idx >= 0:
                sc_center_idx = center_idx
            else:
                sc_center_idx = min(node_atoms, key=lambda i: np.linalg.norm(supercell[i].coords - ctr))

            node_comps = []
            seen = set()
            for seed in node_atoms:
                if seed in seen:
                    continue
                comp = {seed}
                q = deque([seed])
                seen.add(seed)
                while q:
                    u = q.popleft()
                    for v in graph[u]:
                        if v in node_atoms and v not in seen:
                            seen.add(v)
                            comp.add(v)
                            q.append(v)
                node_comps.append(comp)

            center_comp = None
            for comp in node_comps:
                if sc_center_idx in comp:
                    center_comp = comp
                    break
            if center_comp is None:
                center_comp = {sc_center_idx}

            core_nodes = set(center_comp)
            center_sig = {}
            for i in center_comp:
                el = sc_sym[i]
                center_sig[el] = center_sig.get(el, 0) + 1

            # Path I (COF-108-like only): tetra-C-centered aryl node in B/O frameworks.
            if "Si" not in sc_sym and center_sig.get("B", 0) == 1 and center_sig.get("O", 0) == 2 and len(center_comp) == 3:
                c_nodes = []
                for ci in range(len(supercell)):
                    if sc_sym[ci] != "C":
                        continue
                    c_nbs = [j for j in graph[ci] if sc_sym[j] == "C"]
                    h_nbs = [j for j in graph[ci] if sc_sym[j] == "H"]
                    if len(c_nbs) == 4 and len(h_nbs) == 0:
                        c_nodes.append(ci)

                i_cores = []
                for c0 in c_nodes:
                    core = {c0}
                    qh = deque([(c0, 0)])
                    seen_h = {c0}
                    while qh:
                        u, dep = qh.popleft()
                        if dep >= 5:
                            continue
                        for v in graph[u]:
                            if v in seen_h or sc_sym[v] == "H":
                                continue
                            sv = sc_sym[v]
                            if sv in {"B", "O"}:
                                continue
                            if sv != "C":
                                continue
                            seen_h.add(v)
                            core.add(v)
                            qh.append((v, dep + 1))
                    ccount = sum(1 for ii in core if sc_sym[ii] == "C")
                    if ccount >= 12:
                        i_cores.append(core)

                if i_cores:
                    def _ctr(core):
                        return np.mean([supercell[ii].coords for ii in core], axis=0)
                    center_core = min(i_cores, key=lambda c: np.linalg.norm(_ctr(c) - ctr))
                    core_nodes = set(center_core)
                    node_atoms = set().union(*i_cores)
                    node_species = {"C"}
                    if center_idx >= 0:
                        sc_center_idx = center_idx
                    else:
                        c_in_core = [ii for ii in center_core if sc_sym[ii] == "C"]
                        pool = c_in_core if c_in_core else list(center_core)
                        sc_center_idx = min(pool, key=lambda ii: np.linalg.norm(supercell[ii].coords - ctr))
                    path_mode = "I"
                    print(f"  -> COF Path I (COF-108 tetra-C node SBU). cores: {len(i_cores)}, core size: {len(center_core)}")

            layered_candidates = []
            single_block_layer_vec = None
            center_list = list(center_comp)
            for comp in node_comps:
                if comp is center_comp:
                    continue
                sig = {}
                for i in comp:
                    el = sc_sym[i]
                    sig[el] = sig.get(el, 0) + 1
                if sig != center_sig:
                    continue

                min_d = float("inf")
                min_pair = None
                for i in center_list:
                    ci = supercell[i].coords
                    for j in comp:
                        d = np.linalg.norm(ci - supercell[j].coords)
                        if d < min_d:
                            min_d = d
                            min_pair = (i, j)
                layered_candidates.append((min_d, comp, min_pair))

            layered_component_count = 1
            if path_mode == "A" and layered_candidates:
                layered_candidates.sort(key=lambda x: x[0])
                best_d, best_comp, best_pair = layered_candidates[0]
                if 1.0 <= best_d <= 4.5:
                    if best_pair is not None:
                        single_block_layer_vec = np.array(supercell[best_pair[1]].coords) - np.array(supercell[best_pair[0]].coords)
                    layer_mode = getattr(self, "layer_mode", "auto")
                    is_single_block_topology = len(center_comp) <= 2
                    if layer_mode == "monomer" or is_single_block_topology:
                        # Keep single-layer node set. Single-building-block COFs
                        # build dimers later by duplicating the finished fragment
                        # with the closest-layer vector.
                        pass
                    else:
                        # Keep exactly one adjacent layer (dimer-like), not all
                        # symmetry-tied neighbors, to avoid unintended 3-layer
                        # normal fragments in structures like COF-6.
                        core_nodes |= set(best_comp)
                        layered_component_count = 2
                        path_mode = "B"

            if path_mode == "B":
                print(f"  -> COF Path B (Layered set). Node components: {layered_component_count}, spacing: {best_d:.2f} A")
            elif path_mode == "A":
                print(f"  -> COF Path A (Single node). Node component size: {len(center_comp)}")

        elif path_mode == "A":
            # Path C for COF-3xx-like families: use tetra-connected central
            # carbon nodes (carbon bonded to 4 heavy neighbors, mostly carbons).
            c4_nodes = []
            for i in range(len(supercell)):
                if sc_sym[i] != "C":
                    continue
                heavy_nbs = [v for v in graph[i] if sc_sym[v] != "H"]
                if len(heavy_nbs) < 4:
                    continue
                c_nbs = sum(1 for v in heavy_nbs if sc_sym[v] == "C")
                if len(heavy_nbs) >= 4 and c_nbs >= 3:
                    c4_nodes.append(i)

            if c4_nodes:
                node_atoms = set(c4_nodes)
                if center_idx >= 0:
                    sc_center_idx = center_idx
                else:
                    sc_center_idx = min(c4_nodes, key=lambda i: np.linalg.norm(supercell[i].coords - ctr))
                core_nodes = {sc_center_idx}
                path_mode = "C"
                print(f"  -> COF Path C (Tetra-C node). Node candidates: {len(c4_nodes)}")
            else:
                # Path D for porphyrinic COFs (e.g., COF-366): detect N4 core.
                n_atoms = [i for i in range(len(supercell)) if sc_sym[i] == "N"]
                porph_cores = []
                if n_atoms:
                    # N-N graph by heavy-atom shortest-path proximity.
                    n_adj = {i: set() for i in n_atoms}
                    heavy_ok = lambda x: sc_sym[x] != "H"
                    for ni in n_atoms:
                        q = deque([(ni, 0)])
                        seen = {ni}
                        while q:
                            u, dep = q.popleft()
                            if dep >= 6:
                                continue
                            for v in graph[u]:
                                if v in seen or (not heavy_ok(v)):
                                    continue
                                seen.add(v)
                                if sc_sym[v] == "N" and v != ni:
                                    n_adj[ni].add(v)
                                q.append((v, dep + 1))

                    # connected components on N-proximity graph
                    n_vis = set()
                    n_comps = []
                    for ni in n_atoms:
                        if ni in n_vis:
                            continue
                        q = deque([ni])
                        n_vis.add(ni)
                        comp = {ni}
                        while q:
                            u = q.popleft()
                            for v in n_adj[u]:
                                if v not in n_vis:
                                    n_vis.add(v)
                                    comp.add(v)
                                    q.append(v)
                        n_comps.append(comp)

                    for ncomp in n_comps:
                        if len(ncomp) < 4:
                            continue
                        # porphyrin core atoms: within 2 bonds of N set (heavy only)
                        core = set(ncomp)
                        q = deque([(u, 0) for u in ncomp])
                        seen = set(ncomp)
                        while q:
                            u, dep = q.popleft()
                            if dep >= 2:
                                continue
                            for v in graph[u]:
                                if v in seen or (not heavy_ok(v)):
                                    continue
                                seen.add(v)
                                core.add(v)
                                q.append((v, dep + 1))
                        porph_cores.append(core)

                if porph_cores:
                    # pick core nearest supercell center
                    def core_ctr(core):
                        return np.mean([supercell[i].coords for i in core], axis=0)
                    core_nodes = min(porph_cores, key=lambda c: np.linalg.norm(core_ctr(c) - ctr))
                    node_atoms = set().union(*porph_cores)
                    if center_idx >= 0:
                        sc_center_idx = center_idx
                    else:
                        sc_center_idx = min(core_nodes, key=lambda i: np.linalg.norm(supercell[i].coords - ctr))
                    path_mode = "D"
                    print(f"  -> COF Path D (Porphyrin core). N-rich cores: {len(porph_cores)}")
                else:
                    node_atoms = set()
                    if center_idx >= 0:
                        sc_center_idx = center_idx
                    else:
                        sc_center_idx = min(
                            (i for i in range(len(supercell)) if sc_sym[i] != "H"),
                            key=lambda i: np.linalg.norm(supercell[i].coords - ctr),
                        )
                    core_nodes = {sc_center_idx}
                    print("  -> COF Path A (Fallback single-center mode).")
                    print("Warning: Rare fallback COF topology detected (e.g., COF-505-like helix).")
                    print("  -> This topology is not implemented yet in UniFrag.")
        unwrapped = [None] * len(supercell)
        unwrapped[sc_center_idx] = np.array(supercell[sc_center_idx].coords)
        q = deque([sc_center_idx])
        while q:
            i = q.popleft()
            ui = unwrapped[i]
            fi = np.array(supercell[i].frac_coords)
            for j in graph[i]:
                if unwrapped[j] is not None:
                    continue
                fj = np.array(supercell[j].frac_coords)
                dfrac = fj - fi
                dfrac -= np.round(dfrac)
                unwrapped[j] = ui + supercell.lattice.get_cartesian_coords(dfrac)
                q.append(j)
        for i in range(len(unwrapped)):
            if unwrapped[i] is None:
                unwrapped[i] = np.array(supercell[i].coords)

        cpos = np.array(unwrapped[sc_center_idx])
        sphere = {i for i in range(len(supercell)) if np.linalg.norm(np.array(unwrapped[i]) - cpos) <= self.radius}

        # Normal B/O-COF rule: build as node + all fully attached linker
        # components around the kept node set.
        if (not minimize) and (path_mode in {"A", "B"}) and node_atoms:
            final = set(core_nodes)
            broken = []

            non_node = [i for i in range(len(supercell)) if i not in node_atoms]
            non_node_set = set(non_node)
            seen_non = set()
            comps = []
            for seed in non_node:
                if seed in seen_non:
                    continue
                qn = deque([seed])
                seen_non.add(seed)
                comp = {seed}
                while qn:
                    u = qn.popleft()
                    for v in graph[u]:
                        if v in non_node_set and v not in seen_non:
                            seen_non.add(v)
                            comp.add(v)
                            qn.append(v)
                comps.append(comp)

            for comp in comps:
                touches_core = False
                for u in comp:
                    for v in graph[u]:
                        if v in core_nodes:
                            touches_core = True
                            break
                    if touches_core:
                        break
                if touches_core:
                    final |= comp

            # Boundary cuts for capping: any kept atom bonded to excluded atom.
            for u in list(final):
                for v in graph[u]:
                    if v in final:
                        continue
                    broken.append((u, np.array(unwrapped[v]) - np.array(unwrapped[u])))
        else:
            # Grow from core B/O nodes only; this avoids detached islands from
            # direct sphere-seeding while still keeping local linker chemistry.
            final = set(core_nodes)
            visited = set(core_nodes)
            queue = deque(core_nodes)
            broken = []
            while queue:
                u = queue.popleft()
                su = sc_sym[u]
                if u in node_atoms and u not in core_nodes:
                    continue
                for v in graph[u]:
                    if v in node_atoms and v not in core_nodes:
                        broken.append((u, np.array(unwrapped[v]) - np.array(unwrapped[u])))
                        continue
                    if v in visited:
                        continue
                    visited.add(v)
                    final.add(v)
                    queue.append(v)

        single_block_keep_heavy = None

        # Helper export for single-building-block Path-A COFs (e.g., COF-JLU2):
        # keep one ring-centered motif with C/N arms and two O groups.
        try:
            if path_mode == "A" and len(core_nodes) <= 2:
                center_cart = np.array(unwrapped[sc_center_idx], dtype=float)
                node_set = set(node_atoms)
                non_node = [i for i in final if i not in node_set and sc_sym[i] != "H"]
                # Single-block COFs may classify nearly all atoms as node; in that
                # case, fall back to all heavy atoms in final for motif extraction.
                if not non_node:
                    non_node = [i for i in final if sc_sym[i] != "H"]
                non_node_set = set(non_node)

                # pick nearest non-node component touching core
                seen_nn = set()
                best_comp = set()
                best_key = None
                for seed in non_node:
                    if seed in seen_nn:
                        continue
                    qn = deque([seed])
                    seen_nn.add(seed)
                    comp = {seed}
                    touches = False
                    while qn:
                        u = qn.popleft()
                        for v in graph[u]:
                            if v in core_nodes:
                                touches = True
                            if v in non_node_set and v not in seen_nn:
                                seen_nn.add(v)
                                comp.add(v)
                                qn.append(v)
                    if not touches:
                        continue
                    c_count = sum(1 for i in comp if sc_sym[i] == "C")
                    if c_count < 6:
                        continue
                    cctr = np.mean([np.array(unwrapped[i], dtype=float) for i in comp], axis=0)
                    key = (-c_count, float(np.linalg.norm(cctr - center_cart)), -len(comp))
                    if best_key is None or key < best_key:
                        best_key = key
                        best_comp = comp

                if best_comp:
                    carb = [i for i in best_comp if sc_sym[i] == "C"]
                    carb.sort(key=lambda i: float(np.linalg.norm(np.array(unwrapped[i], dtype=float) - center_cart)))
                    ring = set(carb[:6])

                    # extend ring to connected carbons if needed
                    qext = deque(list(ring))
                    while qext and len(ring) < 6:
                        u = qext.popleft()
                        for v in graph[u]:
                            if v in best_comp and sc_sym[v] == "C" and v not in ring:
                                ring.add(v)
                                qext.append(v)
                                if len(ring) >= 6:
                                    break

                    # Keep a local motif around the aromatic ring seeds by graph
                    # distance. This preserves nearby functional groups (N/O arms)
                    # without extending into the full periodic chain.
                    keep_heavy = set(ring)
                    qloc = deque([(i, 0) for i in ring])
                    seen_loc = set(ring)
                    max_hops = 2
                    allowed = {"C", "N", "O", "B", "Si", "P", "S", "F", "Cl", "Br", "I"}
                    while qloc:
                        u, depth = qloc.popleft()
                        if depth >= max_hops:
                            continue
                        for v in graph[u]:
                            if v in seen_loc or v not in best_comp:
                                continue
                            sv = sc_sym[v]
                            if sv not in allowed:
                                continue
                            seen_loc.add(v)
                            keep_heavy.add(v)
                            qloc.append((v, depth + 1))

                    # Ensure terminal hetero atoms on this local block are kept
                    # even when component partitioning classifies them as node-side.
                    extra_hetero = set()
                    for u in list(keep_heavy):
                        if sc_sym[u] != "C":
                            continue
                        for v in graph[u]:
                            if sc_sym[v] in {"N", "O", "B"}:
                                extra_hetero.add(v)
                    keep_heavy |= extra_hetero

                    # Prune duplicated N-N terminal pairs on side arms: keep the N
                    # connected to the carbon framework, drop the distal extra N.
                    drop_n = set()
                    for nidx in list(keep_heavy):
                        if sc_sym[nidx] != "N":
                            continue
                        nn_in_keep = [nb for nb in graph[nidx] if nb in keep_heavy and sc_sym[nb] == "N"]
                        if not nn_in_keep:
                            continue
                        has_c_anchor = any(nb in keep_heavy and sc_sym[nb] == "C" for nb in graph[nidx])
                        if not has_c_anchor:
                            drop_n.add(nidx)
                    keep_heavy -= drop_n
                    if keep_heavy:
                        single_block_keep_heavy = set(keep_heavy)

                    keep = set(keep_heavy)
                    # hydrogens bonded to retained heavy atoms
                    for i in range(len(supercell)):
                        if sc_sym[i] != "H":
                            continue
                        for nb in graph[i]:
                            if nb in keep_heavy:
                                keep.add(i)
                                break

                    blk = sorted(keep)
                    if blk:
                        blk_sp = [sc_sym[i] for i in blk]
                        blk_co = [np.array(unwrapped[i], dtype=float) for i in blk]
                        self._export_cof_component_library(blk_sp, blk_co, structure_stem, "node")
        except Exception:
            pass

        # For single-building-block COFs, construct the actual fragment directly
        # from extracted block(s), then proceed with standard edge capping.
        if single_block_keep_heavy:
            final = set(single_block_keep_heavy)

            def _extract_neighbor_single_block(v0, base_final):
                # Extract a full neighboring block from an attachment seed,
                # mirroring the center ring + functional-group rule.
                outside_heavy = {i for i in range(len(supercell)) if sc_sym[i] != "H" and i not in base_final}
                comp2 = set()
                qn2 = deque([v0])
                seen2 = {v0}
                while qn2:
                    u = qn2.popleft()
                    if u not in outside_heavy:
                        continue
                    comp2.add(u)
                    for w in graph[u]:
                        if w in seen2:
                            continue
                        if w in outside_heavy:
                            seen2.add(w)
                            qn2.append(w)

                if not comp2:
                    return set()

                v0_pos = np.array(unwrapped[v0], dtype=float)
                carb2 = [i for i in comp2 if sc_sym[i] == "C"]
                carb2.sort(key=lambda i: float(np.linalg.norm(np.array(unwrapped[i], dtype=float) - v0_pos)))
                ring2 = set(carb2[:6])

                qext2 = deque(list(ring2))
                while qext2 and len(ring2) < 6:
                    u = qext2.popleft()
                    for w in graph[u]:
                        if w in comp2 and sc_sym[w] == "C" and w not in ring2:
                            ring2.add(w)
                            qext2.append(w)
                            if len(ring2) >= 6:
                                break

                keep2 = set(ring2)
                qloc2 = deque([(i, 0) for i in ring2])
                seen_loc2 = set(ring2)
                allowed2 = {"C", "N", "O", "B", "Si", "P", "S", "F", "Cl", "Br", "I"}
                max_hops2 = 2
                while qloc2:
                    u, depth = qloc2.popleft()
                    if depth >= max_hops2:
                        continue
                    for w in graph[u]:
                        if w in seen_loc2 or w not in comp2:
                            continue
                        sw = sc_sym[w]
                        if sw not in allowed2:
                            continue
                        seen_loc2.add(w)
                        keep2.add(w)
                        qloc2.append((w, depth + 1))

                extra2 = set()
                for u in list(keep2):
                    if sc_sym[u] != "C":
                        continue
                    for w in graph[u]:
                        if sc_sym[w] in {"N", "O", "B"} and w in comp2:
                            extra2.add(w)
                keep2 |= extra2

                drop_n2 = set()
                for nidx in list(keep2):
                    if sc_sym[nidx] != "N":
                        continue
                    nn_in = [nb for nb in graph[nidx] if nb in keep2 and sc_sym[nb] == "N"]
                    if not nn_in:
                        continue
                    has_c_anchor = any(nb in keep2 and sc_sym[nb] == "C" for nb in graph[nidx])
                    if not has_c_anchor:
                        drop_n2.add(nidx)
                keep2 -= drop_n2
                return keep2

            # Min: center + 1 attached block. Normal: center + 3 attached blocks.
            boundary_pairs = []
            center_final = set(final)
            for u in center_final:
                for v in graph[u]:
                    if sc_sym[v] == "H" or v in center_final:
                        continue
                    boundary_pairs.append((u, v))

            if boundary_pairs:
                center_cart = np.array(unwrapped[sc_center_idx], dtype=float)
                nn_pairs = [(u, v) for (u, v) in boundary_pairs if sc_sym[u] == "N" and sc_sym[v] == "N"]
                cand_pairs = nn_pairs if nn_pairs else boundary_pairs
                cand_pairs = sorted(
                    cand_pairs,
                    key=lambda uv: float(np.linalg.norm(np.array(unwrapped[uv[1]], dtype=float) - center_cart)),
                )
                target_count = 1 if minimize else 3
                selected_pairs = []
                used_center_atoms = set()
                used_neighbor_atoms = set()
                for u, v in cand_pairs:
                    if u in used_center_atoms or v in used_neighbor_atoms:
                        continue
                    selected_pairs.append((u, v))
                    used_center_atoms.add(u)
                    used_neighbor_atoms.add(v)
                    if len(selected_pairs) >= target_count:
                        break

                for _, v0 in selected_pairs:
                    final |= _extract_neighbor_single_block(v0, center_final)

            broken = []
            for u in list(final):
                for v in graph[u]:
                    if sc_sym[v] == "H":
                        continue
                    if v in final:
                        continue
                    broken.append((u, np.array(unwrapped[v]) - np.array(unwrapped[u])))

        ordered = sorted(final)
        local = {gi: li for li, gi in enumerate(ordered)}
        species = [sc_sym[i] for i in ordered]
        coords = [np.array(unwrapped[i]) for i in ordered]
        capped_h_flags = [False] * len(species)

        broken_by_parent = {}
        for gi, vec in broken:
            broken_by_parent.setdefault(gi, []).append(vec)

        # Build kept-neighbor map in final fragment to improve capping direction
        # at linker termini (place H away from preserved bonded atoms).
        final_set = set(final)
        kept_neighbor_vectors = {}
        for gi in broken_by_parent:
            neigh_vecs = []
            gpos = np.array(unwrapped[gi])
            for nb in graph[gi]:
                if nb in final_set:
                    neigh_vecs.append(np.array(unwrapped[nb]) - gpos)
            kept_neighbor_vectors[gi] = neigh_vecs

        for gi, vecs in broken_by_parent.items():
            li = local.get(gi)
            if li is None:
                continue
            if species[li] == "H":
                continue
            if species[li] == "O" and self.oxygen_already_protonated(li, species, coords):
                continue

            # Primary direction: opposite of vectors to kept neighbors.
            dir_vec = np.zeros(3)
            kv = kept_neighbor_vectors.get(gi, [])
            if kv:
                for v in kv:
                    dir_vec -= v
            else:
                # Fallback: old behavior from removed-side average.
                for v in vecs:
                    dir_vec += v

            # Cap only true cut sites. Number of H is computed from local
            # valence deficit after cut, and limited by number of broken bonds.
            target_valence = {"C": 4, "B": 3, "Si": 4, "N": 3, "O": 2}
            sp = species[li]
            n_broken = len(vecs)
            if sp in target_valence and n_broken > 0:
                # local coordination in current fragment (after cut)
                cur_deg = 0
                li_pos = np.array(coords[li])
                for jj in range(len(species)):
                    if jj == li:
                        continue
                    d = np.linalg.norm(li_pos - np.array(coords[jj]))
                    if self.is_valid_bond(sp, species[jj], d):
                        cur_deg += 1
                deficit = max(0, target_valence[sp] - cur_deg)
                n_cap = min(n_broken, deficit)

                # oxygen special-case: avoid protonating already protonated O
                if sp == "O" and self.oxygen_already_protonated(li, species, coords):
                    n_cap = 0

                for _ in range(n_cap):
                    self.place_capping_h(li, dir_vec, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        # COF fragments can leave terminal O atoms without a broken heavy-atom
        # edge marker; cap those O sites with H before geometry refinement.
        self._cap_open_oxygens(species, coords, capped_h_flags)

        if minimize and species and path_mode == "H":
            heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
            if heavy_idx:
                hadj = {i: [] for i in heavy_idx}
                for a in range(len(heavy_idx)):
                    i = heavy_idx[a]
                    ci = np.array(coords[i])
                    for b in range(a + 1, len(heavy_idx)):
                        j = heavy_idx[b]
                        d = np.linalg.norm(ci - np.array(coords[j]))
                        if self.is_valid_bond(species[i], species[j], d):
                            hadj[i].append(j)
                            hadj[j].append(i)

                local_nodes = {local[g] for g in core_nodes if g in local}
                organic = [i for i in heavy_idx if i not in local_nodes]
                organic_set = set(organic)

                comps = []
                vis = set()
                for i in organic:
                    if i in vis:
                        continue
                    q = deque([i])
                    vis.add(i)
                    comp = {i}
                    while q:
                        u = q.popleft()
                        for v in hadj[u]:
                            if v in organic_set and v not in vis:
                                vis.add(v)
                                comp.add(v)
                                q.append(v)
                    comps.append(comp)

                best = None
                best_key = None
                for comp in comps:
                    touch = 0
                    bterm = 0
                    for u in comp:
                        if any(v in local_nodes for v in hadj[u]):
                            touch += 1
                        if species[u] == "B":
                            bterm += 1
                    key = (touch, bterm, len(comp))
                    if best_key is None or key > best_key:
                        best_key = key
                        best = comp

                keep_heavy = set(local_nodes)
                if best is not None:
                    keep_heavy |= set(best)

                # Preserve Si/B attached to retained phenyl carbons.
                for i in heavy_idx:
                    if species[i] not in {"Si", "B"}:
                        continue
                    if i in keep_heavy or any((nb in keep_heavy and species[nb] == "C") for nb in hadj.get(i, [])):
                        keep_heavy.add(i)

                keep = set(keep_heavy)
                for i, sp in enumerate(species):
                    if sp != "H":
                        continue
                    ci = np.array(coords[i])
                    for j in keep_heavy:
                        if j >= len(coords):
                            continue
                        d = np.linalg.norm(ci - np.array(coords[j]))
                        if self.is_valid_bond("H", species[j], d):
                            keep.add(i)
                            break

                species = [x for i, x in enumerate(species) if i in keep]
                coords = [x for i, x in enumerate(coords) if i in keep]
                capped_h_flags = [x for i, x in enumerate(capped_h_flags) if i in keep]

        if minimize and species and path_mode != "H":
            # COF minimize target:
            # 1) keep one whole linker branch
            # 2) for other branches keep only first benzene ring near node
            heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
            if heavy_idx:
                hset = set(heavy_idx)
                hadj = {i: [] for i in heavy_idx}
                for a in range(len(heavy_idx)):
                    i = heavy_idx[a]
                    ci = np.array(coords[i])
                    for b in range(a + 1, len(heavy_idx)):
                        j = heavy_idx[b]
                        d = np.linalg.norm(ci - np.array(coords[j]))
                        if self.is_valid_bond(species[i], species[j], d):
                            hadj[i].append(j)
                            hadj[j].append(i)

                local_nodes = {local[g] for g in core_nodes if g in local}
                organic = [i for i in heavy_idx if i not in local_nodes]
                organic_set = set(organic)

                # Organic connected components (excluding node atoms).
                comps = []
                vis = set()
                for i in organic:
                    if i in vis:
                        continue
                    q = deque([i])
                    vis.add(i)
                    comp = {i}
                    while q:
                        u = q.popleft()
                        for v in hadj[u]:
                            if v in organic_set and v not in vis:
                                vis.add(v)
                                comp.add(v)
                                q.append(v)
                    comps.append(comp)

                # Score each component by node attachments.
                comp_info = []
                for comp in comps:
                    touching_nodes = set()
                    bridge_atoms = set()
                    for u in comp:
                        for v in hadj[u]:
                            if v in local_nodes:
                                touching_nodes.add(v)
                                bridge_atoms.add(u)
                    comp_info.append((comp, touching_nodes, bridge_atoms))

                # Choose whole linker branch(es) to keep fully.
                # For layered Path B systems, keep one full branch per heavy
                # layer component so one dimer side is not over-trimmed.
                primary_indices = set()
                if comp_info:
                    if path_mode == "B":
                        heavy_comp_id = {}
                        vis_h = set()
                        comp_id = 0
                        for seed in heavy_idx:
                            if seed in vis_h:
                                continue
                            qh = deque([seed])
                            vis_h.add(seed)
                            while qh:
                                u = qh.popleft()
                                heavy_comp_id[u] = comp_id
                                for v in hadj.get(u, []):
                                    if v not in vis_h:
                                        vis_h.add(v)
                                        qh.append(v)
                            comp_id += 1

                        best_by_layer = {}
                        for idx, (comp, touching_nodes, bridge_atoms) in enumerate(comp_info):
                            layer_ids = [heavy_comp_id.get(u) for u in comp]
                            layer_ids = [x for x in layer_ids if x is not None]
                            layer_id = layer_ids[0] if layer_ids else idx
                            key = (len(touching_nodes), len(comp))
                            if layer_id not in best_by_layer or key > best_by_layer[layer_id][0]:
                                best_by_layer[layer_id] = (key, idx)
                        primary_indices = {idx for _, idx in best_by_layer.values()}
                    else:
                        primary_idx = max(range(len(comp_info)), key=lambda idx: (len(comp_info[idx][1]), len(comp_info[idx][0])))
                        primary_indices = {primary_idx}

                    keep_heavy = set(local_nodes)
                    for idx in primary_indices:
                        keep_heavy |= set(comp_info[idx][0])
                else:
                    keep_heavy = set(local_nodes) if local_nodes else {heavy_idx[0]}

                def find_six_cycle(comp, bridge_atoms):
                    carbon = {u for u in comp if species[u] == "C"}
                    if len(carbon) < 6:
                        return set()
                    cadj = {u: [v for v in hadj[u] if v in carbon] for u in carbon}

                    seeds = [u for u in bridge_atoms if u in carbon]
                    if not seeds:
                        seeds = list(carbon)

                    def dfs(start, cur, path, used):
                        if len(path) == 6:
                            if start in cadj[cur]:
                                return list(path)
                            return None
                        for nb in cadj[cur]:
                            if nb in used:
                                continue
                            used.add(nb)
                            path.append(nb)
                            got = dfs(start, nb, path, used)
                            if got is not None:
                                return got
                            path.pop()
                            used.remove(nb)
                        return None

                    for st in seeds:
                        got = dfs(st, st, [st], {st})
                        if got is not None:
                            return set(got)
                    return set()

                # For non-primary branches keep first benzene ring only.
                for comp_idx, (comp, touching_nodes, bridge_atoms) in enumerate(comp_info):
                    if comp_idx in primary_indices:
                        continue
                    ring = find_six_cycle(comp, bridge_atoms)
                    if ring:
                        keep_heavy |= ring
                    else:
                        # Fallback if no clear benzene ring found: keep up to
                        # six nearest carbons from bridge atoms.
                        seeds = list(bridge_atoms) if bridge_atoms else list(comp)[:1]
                        q = deque(seeds)
                        dist = {u: 0 for u in seeds}
                        seen = set(seeds)
                        carbons = []
                        while q and len(carbons) < 6:
                            u = q.popleft()
                            if u in comp and species[u] == "C":
                                carbons.append(u)
                            for v in hadj[u]:
                                if v in comp and v not in seen:
                                    seen.add(v)
                                    dist[v] = dist[u] + 1
                                    q.append(v)
                        keep_heavy |= set(carbons)

                # Preserve adjacent C atoms for Path A COFs (e.g., COF-102) so
                # linker carbons are not over-pruned before capping.
                if path_mode == "A":
                    for cidx in heavy_idx:
                        if species[cidx] != "C" or cidx in keep_heavy:
                            continue
                        if any((nb in keep_heavy and species[nb] in {"B", "O", "C"}) for nb in hadj.get(cidx, [])):
                            keep_heavy.add(cidx)

                # Preserve N atoms adjacent to kept carbon fragments (important
                # for porphyrinic / imine motifs in COF-3xx), then add capping H
                # if those N become terminal after trimming.
                for nidx in heavy_idx:
                    if species[nidx] != "N" or nidx in keep_heavy:
                        continue
                    if any((nb in keep_heavy and species[nb] in {"C", "N"}) for nb in hadj[nidx]):
                        keep_heavy.add(nidx)

                # Preserve Si attached to retained phenyl carbons.
                for i in heavy_idx:
                    if i >= len(species) or species[i] != "Si":
                        continue
                    if i in keep_heavy or any((nb in keep_heavy and species[nb] == "C") for nb in hadj.get(i, [])):
                        keep_heavy.add(i)

                # Boundary chemistry rule for minimized COFs:
                # 1) Keep O on linker edge.
                # 2) Keep B on node-open side.
                node_region = set(local_nodes)
                # linker-edge O: O attached to kept linker carbon-like atom.
                for i in heavy_idx:
                    if i in keep_heavy or species[i] != "O":
                        continue
                    if any((nb in keep_heavy and nb not in node_region and species[nb] in {"C", "N", "O"}) for nb in hadj.get(i, [])):
                        keep_heavy.add(i)

                # node-open B: B attached to kept node-side atoms.
                for i in heavy_idx:
                    if i in keep_heavy or species[i] != "B":
                        continue
                    if any((nb in keep_heavy and nb in node_region and species[nb] in {"O", "B", "C"}) for nb in hadj.get(i, [])):
                        keep_heavy.add(i)

                # General hetero retention near kept chemistry.
                for i in heavy_idx:
                    if i in keep_heavy:
                        continue
                    if species[i] not in {"O", "B"}:
                        continue
                    if any((nb in keep_heavy and species[nb] in {"C", "B", "O", "N", "Si"}) for nb in hadj.get(i, [])):
                        keep_heavy.add(i)

                # Track how many heavy neighbors each kept atom loses during minimize trim.
                removed_nbr_count = {}
                for ii in keep_heavy:
                    removed_nbr_count[ii] = sum(1 for nb in hadj.get(ii, []) if nb not in keep_heavy)

                # Keep hydrogens bonded to retained heavy atoms.
                keep = set(keep_heavy)
                for i, sp in enumerate(species):
                    if sp != "H":
                        continue
                    ci = np.array(coords[i])
                    for j in keep_heavy:
                        d = np.linalg.norm(ci - np.array(coords[j]))
                        if self.is_valid_bond("H", species[j], d):
                            keep.add(i)
                            break

                kept_idx = [i for i in range(len(species)) if i in keep]
                old_to_new = {old_i: new_i for new_i, old_i in enumerate(kept_idx)}
                species = [x for i, x in enumerate(species) if i in keep]
                coords = [x for i, x in enumerate(coords) if i in keep]
                capped_h_flags = [x for i, x in enumerate(capped_h_flags) if i in keep]

                # Cap only carbons that became dangling due to minimize trim.
                for old_i, n_removed in removed_nbr_count.items():
                    if n_removed <= 0:
                        continue
                    new_i = old_to_new.get(old_i)
                    if new_i is None:
                        continue
                    if species[new_i] != "C":
                        continue

                    # Current coordination in trimmed fragment
                    cpos = np.array(coords[new_i])
                    nbs = []
                    for j in range(len(species)):
                        if j == new_i:
                            continue
                        d = np.linalg.norm(cpos - np.array(coords[j]))
                        if self.is_valid_bond("C", species[j], d):
                            nbs.append(j)

                    cur_deg = len(nbs)
                    deficit = max(0, 4 - cur_deg)
                    n_cap = min(n_removed, deficit)
                    if n_cap <= 0:
                        continue

                    base = np.zeros(3)
                    if nbs:
                        for nb in nbs:
                            base -= (np.array(coords[nb]) - cpos)
                    else:
                        base = np.array([1.0, 0.0, 0.0])

                    for _ in range(n_cap):
                        self.place_capping_h(new_i, base, self.cap_bond_length("C"), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

                # Post-trim N capping: if retained N has low retained valence,
                # place one H opposite to kept neighbors.
                hset2 = [i for i, sp in enumerate(species) if sp != "H"]
                hadj2 = {i: [] for i in hset2}
                for a in range(len(hset2)):
                    i = hset2[a]
                    ci = np.array(coords[i])
                    for b in range(a + 1, len(hset2)):
                        j = hset2[b]
                        d = np.linalg.norm(ci - np.array(coords[j]))
                        if self.is_valid_bond(species[i], species[j], d):
                            hadj2[i].append(j)
                            hadj2[j].append(i)

                for i, sp in list(enumerate(species)):
                    if sp != "N":
                        continue
                    kept_nbs = hadj2.get(i, [])
                    if len(kept_nbs) >= 2:
                        continue
                    base = np.zeros(3)
                    if kept_nbs:
                        for nb in kept_nbs:
                            base -= (np.array(coords[nb]) - np.array(coords[i]))
                    else:
                        base = np.array([1.0, 0.0, 0.0])
                    self.place_capping_h(i, base, self.cap_bond_length("N"), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        # Final minimize Si/B capping: cap undercoordinated Si/B with one H.
        if minimize and species:
            hidx = [i for i, sp in enumerate(species) if sp != "H"]
            hadj_si = {i: [] for i in hidx}
            for a in range(len(hidx)):
                i = hidx[a]
                ci = np.array(coords[i])
                for b in range(a + 1, len(hidx)):
                    j = hidx[b]
                    d = np.linalg.norm(ci - np.array(coords[j]))
                    if self.is_valid_bond(species[i], species[j], d):
                        hadj_si[i].append(j)
                        hadj_si[j].append(i)

            for i, sp in list(enumerate(species)):
                if sp not in {"Si", "B"}:
                    continue
                hcut = 1.7 if sp == "Si" else 1.5
                has_h = any(species[j] == "H" and np.linalg.norm(np.array(coords[i]) - np.array(coords[j])) <= hcut for j in range(len(species)))
                if has_h:
                    continue
                kept_nbs = hadj_si.get(i, [])
                target_val = 4 if sp == "Si" else 3
                if len(kept_nbs) >= target_val:
                    continue
                base = np.zeros(3)
                if kept_nbs:
                    for nb in kept_nbs:
                        base -= (np.array(coords[nb]) - np.array(coords[i]))
                else:
                    base = np.array([1.0, 0.0, 0.0])
                self.place_capping_h(i, base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        # Keep a single connected component for Path A.
        # For Path B (stacked dimer-like layers), keep both principal layers.
        if species:
            adj = [[] for _ in range(len(species))]
            for i in range(len(species)):
                ci = np.array(coords[i])
                for j in range(i + 1, len(species)):
                    d = np.linalg.norm(ci - np.array(coords[j]))
                    if self.is_valid_bond(species[i], species[j], d):
                        adj[i].append(j)
                        adj[j].append(i)
            comps = []
            vis = set()
            for i in range(len(species)):
                if i in vis:
                    continue
                q = deque([i])
                vis.add(i)
                comp = {i}
                while q:
                    u = q.popleft()
                    for v in adj[u]:
                        if v not in vis:
                            vis.add(v)
                            comp.add(v)
                            q.append(v)
                comps.append(comp)
            if len(comps) > 1:
                def comp_key(comp):
                    ns = node_species if node_species else {"B", "O"}
                    node_ct = sum(1 for k in comp if species[k] in ns)
                    return (node_ct, len(comp))
                if path_mode == "B":
                    ranked = sorted(comps, key=comp_key, reverse=True)
                    layer_mode = getattr(self, "layer_mode", "auto")
                    keep_n = 1 if layer_mode == "monomer" else 2
                    keep = set().union(*ranked[:keep_n])
                else:
                    keep = max(comps, key=comp_key)
                species = [x for i, x in enumerate(species) if i in keep]
                coords = [x for i, x in enumerate(coords) if i in keep]
                capped_h_flags = [x for i, x in enumerate(capped_h_flags) if i in keep]

        if single_block_keep_heavy and getattr(self, "layer_mode", "auto") == "dimer":
            layer_vec = locals().get("single_block_layer_vec")
            if layer_vec is not None:
                base_species = list(species)
                base_coords = [np.array(c, dtype=float) for c in coords]
                base_flags = list(capped_h_flags)
                species = base_species + base_species
                coords = base_coords + [c + np.array(layer_vec, dtype=float) for c in base_coords]
                capped_h_flags = base_flags + base_flags

        self._cap_open_oxygens(species, coords, capped_h_flags)
        capped_h_indices = [i for i, is_cap in enumerate(capped_h_flags) if is_cap and species[i] == "H"]
        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)
        _label = str(output_path or 'cof_fragment')
        species, coords, capped_h_indices = self.fix_odd_electron_multiplicity(species, coords, capped_h_indices, label=_label)
        print(f"Final size: {len(species)} atoms")
        mol = Molecule(species, coords)
        if output_path:
            mol.to(filename=output_path, fmt="xyz")
            print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)


# ---------------------------------------------------------------------------
# BioMolFragmenter — sliding-window fragmenter for single-chain PDB structures
# ---------------------------------------------------------------------------

class BioMolFragmenter(BaseFragmenter):
    """Sliding-window fragmenter for biological macromolecules (PDB format).

    Strategy
    --------
    - Parse residues from a PDB ATOM record chain.
    - Slide a window of ``window_size`` residues along the sequence with
      ``stride`` residue steps.
    - For each window:
        * Collect all heavy atoms from the selected residues.
        * Add geometry-estimated H atoms for every heavy atom using standard
          bond lengths and idealized angles (no external force-field needed).
        * Cap the N-terminal cut with an ACE group (CH3-C(=O)-) attached to
          the window's first backbone N.
        * Cap the C-terminal cut with an NME group (-NH-CH3) attached to the
          window's last backbone C.
        * Write the window as an XYZ file.
    - Write a summary CSV with residue range, sequence, and atom counts.

    Parameters
    ----------
    window_size : int
        Number of residues per fragment window (default 5).
    stride : int
        Step size in residues between consecutive windows (default 1).

    H atoms are added **only** at the two severed peptide bonds (a single H cap
    on the N-terminal cut, and a single H cap on the C-terminal cut). Interior
    heavy atoms are taken directly from the PDB coordinates and are not modified.
    """

    # Standard covalent radii (Å) for bio elements
    COV_RAD = {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
        "S": 1.05, "P": 1.07, "Se": 1.20,
    }

    # Ideal bond lengths (Å) used when adding H atoms
    _BL = {"C": 1.09, "N": 1.01, "O": 0.96, "S": 1.34}

    # Expected heavy-atom valence for H-addition counting
    _VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2}

    # Backbone atom names (PDB naming)
    _BACKBONE = {"N", "CA", "C", "O", "OXT"}

    # Three-letter → one-letter amino acid code
    _AA1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "HSD": "H", "HSE": "H", "HSP": "H", "HIE": "H", "HID": "H",
        "HIP": "H", "MSE": "M", "SEC": "U",
    }

    def __init__(self, window_size=5, stride=1, use_pdbfixer=True, ph=7.0):
        """Initialise the fragmenter.

        Parameters
        ----------
        window_size : int
            Number of residues per sliding window (default 5).
        stride : int
            Step size in residues between consecutive windows (default 1).
        use_pdbfixer : bool
            If True (default), run PDBFixer on the input PDB before parsing
            to add missing heavy atoms and hydrogens.  Requires the
            ``pdbfixer`` and ``openmm`` packages.
        ph : float
            pH used by PDBFixer when adding missing hydrogens (default 7.0).
        """
        super().__init__(radius=0.0)
        if window_size < 1:
            raise ValueError("window_size must be >= 1")
        if stride < 1:
            raise ValueError("stride must be >= 1")
        self.window_size = int(window_size)
        self.stride = int(stride)
        self.use_pdbfixer = bool(use_pdbfixer)
        self.ph = float(ph)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def extract(self, pdb_path, output_dir=None, chain=None, write_files=True):
        """Fragment a PDB file into sliding-window XYZ files.

        Parameters
        ----------
        pdb_path : str | Path
            Path to the input PDB file.
        output_dir : str | Path | None
            Directory for output XYZ files and summary CSV.
            Defaults to a ``bio_fragments/`` sub-folder next to the PDB.
        chain : str | None
            Chain ID to extract (e.g. ``"B"``).  If None, the first chain
            with ATOM records is used.
        write_files : bool
            If True, write out individual XYZ files and a summary CSV.
            If False, only return the FragmentResult objects without writing files.

        Returns
        -------
        list[FragmentResult]
            One entry per window.
        """
        pdb_path = Path(pdb_path)
        if write_files:
            if output_dir is None:
                output_dir = pdb_path.parent / "bio_fragments"
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # --- optional PDBFixer pre-processing ----------------------------
        parse_path = pdb_path
        pdb_has_h = False          # True once PDBFixer has added H atoms
        temp_fixed_path = None
        if self.use_pdbfixer:
            fixed_pdb_str = self._fix_with_pdbfixer(pdb_path, self.ph)
            if fixed_pdb_str is not None:
                if write_files:
                    fixed_path = output_dir / (pdb_path.stem + "_fixed.pdb")
                    fixed_path.write_text(fixed_pdb_str)
                    parse_path = fixed_path
                    print(f"PDBFixer: fixed PDB written to {fixed_path.name}")
                else:
                    temp_fixed_path = pdb_path.parent / f".{pdb_path.stem}_fixed.pdb"
                    temp_fixed_path.write_text(fixed_pdb_str)
                    parse_path = temp_fixed_path
                pdb_has_h = True
            else:
                print("PDBFixer unavailable or failed — using original PDB.")
        # -----------------------------------------------------------------

        try:
            residues = self._parse_pdb(parse_path, chain=chain, keep_h=pdb_has_h)
            if not residues:
                raise ValueError(f"No ATOM records found in {pdb_path} (chain={chain!r})")

            stem = pdb_path.stem
            results = []
            csv_rows = [["window", "res_start", "res_end", "sequence",
                         "n_residues", "n_heavy", "n_total", "xyz_file"]]

            n_res = len(residues)
            win_idx = 0
            start = 0
            seen_keys = set()
            while start < n_res:
                end = min(start + self.window_size, n_res)
                window_res = residues[start:end]

                species, coords, capped_h_flags = self._build_window(
                    window_res,
                    is_nterm=(start == 0),
                    is_cterm=(end == n_res),
                )

                res_ids   = [r["res_seq"] for r in window_res]
                res_names = [r["res_name"] for r in window_res]
                seq1 = "".join(self._AA1.get(rn, "X") for rn in res_names)
                n_heavy = sum(1 for s in species if s != "H")

                key = MOFFragmenter._species_coords_unique_key(species, coords, decimals=1)
                if key in seen_keys:
                    print(
                        f"Window {win_idx:3d}  res {res_ids[0]:3d}-{res_ids[-1]:<3d}"
                        f"  seq={seq1}  (Skipped, duplicate)"
                    )
                    csv_rows.append([
                        win_idx, res_ids[0], res_ids[-1], seq1,
                        len(window_res), n_heavy, len(species), "duplicate",
                    ])
                else:
                    seen_keys.add(key)
                    fname = f"{stem}_w{win_idx:03d}_r{res_ids[0]}-{res_ids[-1]}.xyz"

                    if write_files:
                        out_path = output_dir / fname
                        # write XYZ
                        mol = Molecule(species, coords)
                        mol.to(filename=str(out_path), fmt="xyz")
                        print(
                            f"Window {win_idx:3d}  res {res_ids[0]:3d}-{res_ids[-1]:3d}"
                            f"  seq={seq1}  heavy={n_heavy}  total={len(species)}"
                            f"  → {out_path.name}"
                        )
                    else:
                        print(
                            f"Window {win_idx:3d}  res {res_ids[0]:3d}-{res_ids[-1]:3d}"
                            f"  seq={seq1}  heavy={n_heavy}  total={len(species)}"
                        )

                    csv_rows.append([
                        win_idx, res_ids[0], res_ids[-1], seq1,
                        len(window_res), n_heavy, len(species), fname,
                    ])
                    results.append(FragmentResult(species=species, coords=coords))

                win_idx += 1
                start += self.stride
                if end == n_res:
                    break  # last window already processed

            if write_files:
                # write summary CSV
                csv_path = output_dir / f"{stem}_windows.csv"
                with open(csv_path, "w") as fh:
                    for row in csv_rows:
                        fh.write(",".join(str(v) for v in row) + "\n")
                print(f"\nSummary CSV: {csv_path}")
                print(f"Total windows: {len(results)}")
            
            return results

        finally:
            if temp_fixed_path and temp_fixed_path.exists():
                try:
                    temp_fixed_path.unlink()
                except Exception:
                    pass

    # ------------------------------------------------------------------
    # PDBFixer pre-processing
    # ------------------------------------------------------------------

    @staticmethod
    def _fix_with_pdbfixer(pdb_path, ph=7.0):
        """Run PDBFixer on *pdb_path* and return the fixed PDB as a string.

        Performs, in order:
        1. Find and fill missing residues in loops/termini.
        2. Find and add missing heavy atoms (side-chain and backbone).
        3. Add missing hydrogen atoms at the requested pH.

        Returns the fixed PDB text as a ``str``, or ``None`` if PDBFixer or
        OpenMM is not installed or if fixing fails for any reason.
        """
        try:
            import io
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
        except ImportError:
            return None

        try:
            fixer = PDBFixer(filename=str(pdb_path))

            # Fill in missing residues (inserts template coords for gaps)
            fixer.findMissingResidues()
            # Find and add missing heavy atoms (e.g. incomplete side chains)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            # Add H at the requested pH
            fixer.addMissingHydrogens(pH=ph)

            buf = io.StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, buf)
            return buf.getvalue()

        except Exception as exc:
            print(f"PDBFixer warning: {exc}  — continuing without fixing.")
            return None

    # ------------------------------------------------------------------
    # PDB parser
    # ------------------------------------------------------------------

    def _parse_pdb(self, pdb_path, chain=None, keep_h=False):
        """Return ordered list of residue dicts from ATOM records.

        Each dict has keys: ``res_seq``, ``res_name``, ``chain_id``,
        ``atoms``  (list of ``{"name", "element", "coord"}``).

        Parameters
        ----------
        keep_h : bool
            If True, include H atoms found in the PDB (e.g. from a PDBFixer
            fixed file).  Default False (skip all H).
        """
        raw = {}        # (chain_id, res_seq, ins_code) → residue dict
        order = []      # insertion-ordered keys

        with open(pdb_path) as fh:
            for line in fh:
                rec = line[:6].strip()
                if rec not in ("ATOM", "HETATM"):
                    continue
                # PDB column layout (1-indexed in spec, 0-indexed here)
                chain_id  = line[21]
                res_name  = line[17:20].strip()
                try:
                    res_seq = int(line[22:26].strip())
                except ValueError:
                    continue
                ins_code  = line[26]
                atom_name = line[12:16].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                # element column (cols 77-78, 0-indexed 76-77)
                element = line[76:78].strip() if len(line) >= 78 else ""
                if not element:
                    # fall back: first character of atom name that is alphabetic
                    element = "".join(c for c in atom_name if c.isalpha())[:1]
                element = element.capitalize()

                # Skip H atoms from original PDB (no H) unless the file was
                # produced by PDBFixer, in which case keep_h=True is passed.
                if not keep_h and (element == "H" or atom_name.startswith("H")):
                    continue
                # Skip water, metals, HETATM non-amino-acids
                if rec == "HETATM" and res_name not in self._AA1:
                    continue

                key = (chain_id, res_seq, ins_code)
                if key not in raw:
                    raw[key] = {
                        "chain_id": chain_id,
                        "res_seq": res_seq,
                        "res_name": res_name,
                        "atoms": [],
                    }
                    order.append(key)
                raw[key]["atoms"].append({
                    "name": atom_name,
                    "element": element,
                    "coord": np.array([x, y, z], dtype=float),
                })

        # select chain
        if chain is None:
            # pick chain with the most residues
            from collections import Counter
            cnt = Counter(k[0] for k in order)
            if not cnt:
                return []
            chain = cnt.most_common(1)[0][0]

        residues = [raw[k] for k in order if k[0] == chain]
        return residues

    # ------------------------------------------------------------------
    # Window assembly
    # ------------------------------------------------------------------

    def _build_window(self, window_res, is_nterm=False, is_cterm=False):
        """Collect atoms, add H, add ACE/NME caps.

        Returns (species, coords, capped_h_flags).
        """
        species: list = []
        coords:  list = []
        capped_h_flags: list = []

        # 1) collect heavy atoms from each residue (coordinates taken directly
        #    from the PDB — no H added to interior atoms)
        for res in window_res:
            for atom in res["atoms"]:
                species.append(atom["element"])
                coords.append(atom["coord"].copy())
                capped_h_flags.append(False)

        # 2) Simple H cap at N-terminus of the window (if not the real N-term).
        #    Only the bond that was severed between the previous residue's C
        #    and this residue's N is capped.
        if not is_nterm:
            self._add_n_term_h_cap(window_res[0], species, coords, capped_h_flags)

        # 3) Simple H cap at C-terminus of the window (if not the real C-term).
        #    Only the bond that was severed between this residue's C and the
        #    next residue's N is capped.
        if not is_cterm:
            self._add_c_term_h_cap(window_res[-1], species, coords, capped_h_flags)

        # 4) Neutralize the fragment (adjust side-chain protonation states)
        species, coords, capped_h_flags = self._neutralize_window(species, coords, capped_h_flags)

        # 5) Ensure only a single connected component remains (remove floating molecules)
        species, coords, capped_h_flags = self.enforce_single_molecule(species, coords, capped_h_flags)

        # 6) geometry-clean only the cap H atoms (per project decision)
        cap_h_indices = [i for i, f in enumerate(capped_h_flags) if f and species[i] == "H"]
        self.optimize_capped_h_geometry_only(species, coords, cap_h_indices)
        species, coords, cap_h_indices = self.fix_odd_electron_multiplicity(species, coords, cap_h_indices, label='bio_fragment')

        return species, coords, capped_h_flags

    # ------------------------------------------------------------------
    # Neutralization
    # ------------------------------------------------------------------

    def _neutralize_window(self, species, coords, capped_h_flags):
        """Ensure the fragment is formally neutral by adding/removing H atoms.
        
        Uses a distance-based bond graph to detect common charged functional groups
        at pH 7 (carboxylates, ammonium, guanidinium, imidazolium) and adjusts
        their protonation state to neutral.
        """
        import numpy as np

        n = len(species)
        adj = {i: [] for i in range(n)}
        for i in range(n):
            ci = coords[i]
            for j in range(i + 1, n):
                d = float(np.linalg.norm(ci - coords[j]))
                if self.is_valid_bond(species[i], species[j], d):
                    adj[i].append(j)
                    adj[j].append(i)

        to_remove = set()
        to_add_sp = []
        to_add_co = []
        to_add_flags = []

        # 1. Carboxylates (-COO-) -> add H
        for i in range(n):
            if species[i] == "C":
                o_nbs = [j for j in adj[i] if species[j] == "O"]
                if len(o_nbs) == 2 and len(adj[o_nbs[0]]) == 1 and len(adj[o_nbs[1]]) == 1:
                    o_idx = o_nbs[0]
                    c_pos = coords[i]
                    o_pos = coords[o_idx]
                    vec = o_pos - c_pos
                    if np.linalg.norm(vec) < 1e-6:
                        continue
                    vec /= np.linalg.norm(vec)
                    e1, e2 = self._orthonormal_basis(vec)
                    angle = np.deg2rad(109.5)
                    dv = np.cos(angle) * vec + np.sin(angle) * e1
                    dv /= np.linalg.norm(dv)
                    h_pos = o_pos + 0.96 * dv
                    
                    to_add_sp.append("H")
                    to_add_co.append(h_pos)
                    to_add_flags.append(True) # Flag as cap to optimize its angle
                    adj[o_idx].append(-1) # Prevent double processing

        # 2. Ammonium (-NH3+) -> remove one H
        for i in range(n):
            if species[i] == "N" and len(adj[i]) == 4:
                h_nbs = [j for j in adj[i] if species[j] == "H"]
                if len(h_nbs) >= 1:
                    to_remove.add(h_nbs[0])

        # 3. Guanidinium (Arg) -> remove one H
        for i in range(n):
            if species[i] == "C":
                n_nbs = [j for j in adj[i] if species[j] == "N"]
                if len(n_nbs) == 3:
                    n_h_nbs = []
                    for n_idx in n_nbs:
                        n_h_nbs.extend([j for j in adj[n_idx] if species[j] == "H"])
                    if len(n_h_nbs) == 5:
                        for n_idx in n_nbs:
                            hs = [j for j in adj[n_idx] if species[j] == "H"]
                            if len(hs) == 2:
                                to_remove.add(hs[0])
                                break

        # 4. Imidazolium (Protonated His) -> remove one H
        his_n = []
        for i in range(n):
            if species[i] == "N" and len(adj[i]) == 3:
                c_nbs = [j for j in adj[i] if species[j] == "C"]
                h_nbs = [j for j in adj[i] if species[j] == "H"]
                if len(c_nbs) == 2 and len(h_nbs) == 1:
                    his_n.append((i, h_nbs[0]))
        
        for idx1, (n1, h1) in enumerate(his_n):
            for idx2 in range(idx1 + 1, len(his_n)):
                n2, h2 = his_n[idx2]
                common_c = set(adj[n1]).intersection(adj[n2])
                if common_c and any(species[c] == "C" for c in common_c):
                    to_remove.add(h2)

        if to_remove:
            keep = [i for i in range(n) if i not in to_remove]
            species = [species[i] for i in keep]
            coords = [coords[i] for i in keep]
            capped_h_flags = [capped_h_flags[i] for i in keep]

        if to_add_sp:
            species.extend(to_add_sp)
            coords.extend(to_add_co)
            capped_h_flags.extend(to_add_flags)

        return species, coords, capped_h_flags

    # ------------------------------------------------------------------
    # Hydrogen addition
    # ------------------------------------------------------------------

    def _build_bond_graph(self, species, coords):
        """Return adjacency dict of heavy-atom indices."""
        n = len(species)
        adj = {i: [] for i in range(n)}
        for i in range(n):
            if species[i] == "H":
                continue
            ci = coords[i]
            for j in range(i + 1, n):
                if species[j] == "H":
                    continue
                d = float(np.linalg.norm(ci - coords[j]))
                if self.is_valid_bond(species[i], species[j], d):
                    adj[i].append(j)
                    adj[j].append(i)
        return adj

    def _add_heavy_hydrogens(self, species, coords, capped_h_flags):
        """Add geometry-estimated H atoms to all heavy atoms.

        Uses valence deficit: H_count = expected_valence − heavy_bond_count.
        H directions are placed using idealized tetrahedral/trigonal geometry.
        """
        # snapshot of heavy-atom-only list for adjacency
        heavy_idx = [i for i, s in enumerate(species) if s != "H"]
        # build adjacency among heavy atoms
        adj: dict[int, list] = {i: [] for i in heavy_idx}
        for a, i in enumerate(heavy_idx):
            ci = coords[i]
            for b in range(a + 1, len(heavy_idx)):
                j = heavy_idx[b]
                d = float(np.linalg.norm(ci - coords[j]))
                if self.is_valid_bond(species[i], species[j], d):
                    adj[i].append(j)
                    adj[j].append(i)

        new_sp:    list = []
        new_co:    list = []
        new_flags: list = []

        for i in heavy_idx:
            el = species[i]
            expected = self._VALENCE.get(el, 0)
            if expected == 0:
                continue
            nbs = adj[i]
            n_h_needed = expected - len(nbs)
            if n_h_needed <= 0:
                continue

            bl = self._BL.get(el, 1.09)
            pi = coords[i]

            if len(nbs) == 0:
                # isolated atom — place H along x
                dirs = [np.array([1.0, 0.0, 0.0])]
            elif len(nbs) == 1:
                # one heavy neighbor: use tetrahedral/trigonal directions
                pj = coords[nbs[0]]
                u = pi - pj
                u /= np.linalg.norm(u)
                e1, e2 = self._orthonormal_basis(u)
                if el == "C" and expected == 3:
                    # sp2 — one H in plane, bisecting the 120° gap
                    dirs = [u]
                elif n_h_needed == 1:
                    # sp3-like: tetrahedral off-axis
                    angle = np.deg2rad(109.5)
                    dirs = [np.cos(angle) * u + np.sin(angle) * e1]
                else:
                    # sp3 with multiple H: distribute around axis
                    angle = np.deg2rad(109.5)
                    phi_step = 2 * np.pi / n_h_needed
                    dirs = []
                    for k in range(n_h_needed):
                        phi = k * phi_step
                        d = (np.cos(angle) * u
                             + np.sin(angle) * (np.cos(phi) * e1 + np.sin(phi) * e2))
                        dirs.append(d / np.linalg.norm(d))
            elif len(nbs) == 2:
                # two heavy neighbors: bisector direction (sp3) or anti-bisector (sp2)
                u1 = coords[nbs[0]] - pi
                u2 = coords[nbs[1]] - pi
                u1 /= np.linalg.norm(u1); u2 /= np.linalg.norm(u2)
                bisect = -(u1 + u2)
                nb = np.linalg.norm(bisect)
                if nb < 1e-8:
                    bisect, _ = self._orthonormal_basis(u1)
                else:
                    bisect /= nb
                if n_h_needed == 1:
                    dirs = [bisect]
                else:
                    e1, _ = self._orthonormal_basis(bisect)
                    ang = np.deg2rad(109.5)
                    phi_step = 2 * np.pi / n_h_needed
                    dirs = []
                    for k in range(n_h_needed):
                        phi = k * phi_step
                        d = np.cos(ang) * bisect + np.sin(ang) * np.cos(phi) * e1
                        d /= np.linalg.norm(d)
                        dirs.append(d)
            else:
                # 3+ neighbors: use sum-of-bonds direction
                nb_sum = np.zeros(3)
                for j in nbs:
                    v = coords[j] - pi
                    nb_sum += v / np.linalg.norm(v)
                anti = -nb_sum
                na = np.linalg.norm(anti)
                if na < 1e-8:
                    anti = np.array([0.0, 0.0, 1.0])
                else:
                    anti /= na
                dirs = [anti]

            for k, dv in enumerate(dirs[:n_h_needed]):
                dv = np.array(dv, dtype=float)
                dv /= np.linalg.norm(dv)
                hp = pi + bl * dv
                new_sp.append("H")
                new_co.append(hp)
                new_flags.append(False)   # native H, not a cap

        species.extend(new_sp)
        coords.extend(new_co)
        capped_h_flags.extend(new_flags)

    # ------------------------------------------------------------------
    # Simple H capping
    # ------------------------------------------------------------------

    def _get_backbone(self, res, *names):
        """Return coordinates of backbone atoms by PDB name."""
        result = {}
        for atom in res["atoms"]:
            if atom["name"] in names:
                result[atom["name"]] = atom["coord"]
        return result

    def _add_n_term_h_cap(self, first_res, species, coords, capped_h_flags):
        """Add a simple H cap to the N-terminal cut of the window.
        
        The H is placed along the CA -> N bond extended backwards.
        """
        bb = self._get_backbone(first_res, "N", "CA")
        if "N" not in bb or "CA" not in bb:
            return

        n_pos  = bb["N"]
        ca_pos = bb["CA"]

        # Direction from CA → N
        bond_vec = n_pos - ca_pos
        bond_len = np.linalg.norm(bond_vec)
        if bond_len < 1e-6:
            return
        u = bond_vec / bond_len

        # H placed 1.01 Å from N along the extension
        h_pos = n_pos + u * 1.01

        species.append("H")
        coords.append(h_pos)
        capped_h_flags.append(True)

    def _add_c_term_h_cap(self, last_res, species, coords, capped_h_flags):
        """Add a simple H cap to the C-terminal cut of the window.
        
        The H is placed along the CA -> C bond extended outward.
        """
        bb = self._get_backbone(last_res, "C", "CA")
        if "C" not in bb or "CA" not in bb:
            return

        c_pos  = bb["C"]
        ca_pos = bb["CA"]

        # Direction from CA → C
        bond_vec = c_pos - ca_pos
        bond_len = np.linalg.norm(bond_vec)
        if bond_len < 1e-6:
            return
        u = bond_vec / bond_len

        # H placed 1.09 Å from C along extension (typical C-H bond length)
        h_pos = c_pos + u * 1.09

        species.append("H")
        coords.append(h_pos)
        capped_h_flags.append(True)

    # ------------------------------------------------------------------
    # Bond validity (organic molecules only, no metals)
    # ------------------------------------------------------------------

    def is_valid_bond(self, s1, s2, dist):
        r1 = self.COV_RAD.get(s1, 0.77)
        r2 = self.COV_RAD.get(s2, 0.77)
        if s1 == "H" and s2 == "H":
            return dist < 0.9
        cutoff = 1.3 * (r1 + r2)
        cutoff = min(2.0, max(0.9, cutoff))
        return dist <= cutoff


def _get_formula(species):
    from collections import Counter
    c = Counter(species)
    return " ".join(f"{el}{c[el]}" for el in sorted(c))


def _process_cof_file(args_tuple):
    cif_path, radius, center, layer_mode = args_tuple
    import os
    base = os.path.basename(cif_path)
    if base.endswith(".cif"):
        base = base[:-4]
    try:
        frag = COFFragmenter(radius=radius, layer_mode=layer_mode)
        res = frag.extract(cif_path, center_idx=center, output_path=None, minimize=False)
        norm_atoms = len(res.species) if res else 0
        norm_formula = _get_formula(res.species) if res else "N/A"
        min_atoms = "N/A"
        min_formula = "N/A"
        res_min = None
        if res and len(res.species) > 50:
            print(f"[{base}] Normal size > 50. Auto-generating minimize version...")
            res_min = frag.extract(cif_path, center_idx=center, output_path=None, minimize=True)
            if res_min:
                min_atoms = len(res_min.species)
                min_formula = _get_formula(res_min.species)
        return {
            "cif": os.path.basename(cif_path),
            "norm_atoms": norm_atoms,
            "norm_formula": norm_formula,
            "min_atoms": min_atoms,
            "min_formula": min_formula,
            "norm_res": res,
            "min_res": res_min,
        }
    except Exception as e:
        print(f"[{base}] Error: {e}")
        return {
            "cif": os.path.basename(cif_path),
            "norm_atoms": "ERROR", "norm_formula": "ERROR",
            "min_atoms": "ERROR", "min_formula": "ERROR",
            "norm_res": None, "min_res": None,
        }


def _process_bio_file(args_tuple):
    pdb_path, window_size, stride, use_pdbfixer, ph, chain = args_tuple
    import os
    base = os.path.basename(pdb_path)
    if base.lower().endswith(".pdb"):
        base = base[:-4]
    try:
        frag = BioMolFragmenter(
            window_size=window_size, stride=stride,
            use_pdbfixer=use_pdbfixer, ph=ph,
        )
        # Pass write_files=False to avoid saving individual window XYZs on disk
        results = frag.extract(pdb_path, output_dir=None, chain=chain, write_files=False)
        n_windows = len(results) if results else 0
        total_atoms = sum(len(r.species) for r in results) if results else 0
        return {
            "pdb": os.path.basename(pdb_path),
            "n_windows": n_windows,
            "total_atoms": total_atoms,
            "results": results or [],
        }
    except Exception as e:
        print(f"[{base}] Error: {e}")
        return {
            "pdb": os.path.basename(pdb_path),
            "n_windows": "ERROR",
            "total_atoms": "ERROR",
            "results": [],
        }

def _process_mof_file(args_tuple):
    cif_path, radius, center, nmetals = args_tuple
    import os
    base = os.path.basename(cif_path)
    if base.endswith(".cif"):
        base = base[:-4]
    
    # Passing None to skip individual file saving in folder mode
    out_norm = None
    out_min = None
    
    try:
        frag = MOFFragmenter(radius=radius)
        res = frag.extract(cif_path, center_idx=center, nmetals=nmetals,
                           output_path=out_norm, minimize=False)
        
        norm_atoms = len(res.species) if res else 0
        norm_formula = _get_formula(res.species) if res else "N/A"
        
        min_atoms = "N/A"
        min_formula = "N/A"
        res_min = None
        
        if res and len(res.species) > 50:
            print(f"[{base}] Normal size > 50. Auto-generating minimize version...")
            res_min = frag.extract(cif_path, center_idx=center, nmetals=nmetals,
                                   output_path=out_min, minimize=True)
            if res_min:
                min_atoms = len(res_min.species)
                min_formula = _get_formula(res_min.species)
                
        return {
            "cif": os.path.basename(cif_path),
            "norm_atoms": norm_atoms,
            "norm_formula": norm_formula,
            "min_atoms": min_atoms,
            "min_formula": min_formula,
            "norm_res": res,
            "min_res": res_min
        }
    except Exception as e:
        print(f"[{base}] Error: {e}")
        return {
            "cif": os.path.basename(cif_path),
            "norm_atoms": "ERROR",
            "norm_formula": "ERROR",
            "min_atoms": "ERROR",
            "min_formula": "ERROR",
            "norm_res": None,
            "min_res": None
        }


def main():
    parser = argparse.ArgumentParser(
        description="OOP fragmenter core for MOF / COF / Bio-macromolecule"
    )
    parser.add_argument("input_path", help="CIF/PDB file or folder containing CIF files")
    parser.add_argument("--kind", choices=["mof", "cof", "bio"], default="mof")
    parser.add_argument("--nproc", type=int, default=1, help="Number of parallel processes for folder mode.")
    # MOF / COF shared
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--nmetals", type=int, default=3)
    parser.add_argument("--output", default="fragment.xyz")
    parser.add_argument("--minimize", action="store_true")
    parser.add_argument(
        "--cof-layer",
        choices=["auto", "monomer", "dimer"],
        default="auto",
        help="COF layer output mode for layered COFs.",
    )
    # Bio-specific
    parser.add_argument(
        "--window-size", type=int, default=5,
        help="Number of residues per sliding window (bio mode, default 5).",
    )
    parser.add_argument(
        "--stride", type=int, default=1,
        help="Residue step between consecutive windows (bio mode, default 1).",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory for bio-fragment XYZ files and summary CSV.",
    )
    parser.add_argument(
        "--chain", default=None,
        help="PDB chain ID to extract (bio mode). Default: largest chain.",
    )
    parser.add_argument(
        "--ph", type=float, default=7.0,
        help="pH for PDBFixer hydrogen placement (bio mode, default 7.0).",
    )
    parser.add_argument(
        "--no-pdbfixer", dest="no_pdbfixer", action="store_true",
        help="Skip PDBFixer pre-processing and use the PDB as-is (bio mode).",
    )
    args = parser.parse_args()

    if args.kind == "mof":
        import os
        import glob
        import multiprocessing
        
        is_dir = os.path.isdir(args.input_path)
        if is_dir:
            cif_files = sorted(glob.glob(os.path.join(args.input_path, "*.cif")))
            out_dir = args.input_path
        else:
            if not os.path.exists(args.input_path):
                print(f"File not found: {args.input_path}")
                return
            cif_files = [args.input_path]
            out_dir = os.path.dirname(os.path.abspath(args.input_path))
            
        if not cif_files:
            print(f"No CIF files found to process.")
            return
            
        csv_path = os.path.join(out_dir, "fragmentation_summary.csv")
        extxyz_path = os.path.join(out_dir, "fragments_collection.extxyz")
        
        seen_keys = set()
        if is_dir:
            # Initialize files
            with open(csv_path, "w") as f:
                f.write("cif_file,normal_atoms,normal_formula,norm_duplicate,min_atoms,min_formula,min_duplicate\n")
            if os.path.exists(extxyz_path):
                os.remove(extxyz_path)
        else:
            if os.path.exists(extxyz_path):
                seen_keys = _load_seen_keys_from_extxyz(extxyz_path)
                
        pool_args = [(cif, args.radius, args.center, args.nmetals) for cif in cif_files]
        
        def _flush_mof_result(r, seen_keys, csv_path, extxyz_path, is_dir):
            norm_dup = "no"
            min_dup = "no"
            if r['norm_res']:
                ckey = MOFFragmenter._chemical_identity_key(r['norm_res'].species, r['norm_res'].coords)
                if ckey in seen_keys: norm_dup = "yes"
            if r['min_res']:
                ckey = MOFFragmenter._chemical_identity_key(r['min_res'].species, r['min_res'].coords)
                if ckey in seen_keys: min_dup = "yes"

            new_extxyz_frags = []
            for key in ["norm_res", "min_res"]:
                res_obj = r[key]
                if res_obj:
                    ckey = MOFFragmenter._chemical_identity_key(res_obj.species, res_obj.coords)
                    if ckey in seen_keys:
                        continue
                    seen_keys.add(ckey)
                    suffix = "Min" if key == "min_res" else ""
                    base_name = r["cif"]
                    if base_name.endswith(".cif"): base_name = base_name[:-4]
                    clean_base = base_name.replace("[", "").replace("]", "").replace("_", "")
                    frag_name = f"{clean_base}FragMof{suffix}"
                    new_extxyz_frags.append((res_obj.species, res_obj.coords, frag_name))

            csv_row = f"{r['cif']},{r['norm_atoms']},{r['norm_formula']},{norm_dup},{r['min_atoms']},{r['min_formula']},{min_dup}\n"
            
            if is_dir:
                with open(csv_path, "a") as f:
                    f.write(csv_row)
                for species, coords, frag_name in new_extxyz_frags:
                    _write_extxyz(extxyz_path, species, coords, frag_name)
            else:
                base_name = r["cif"]
                if base_name.endswith(".cif"): base_name = base_name[:-4]
                clean_base = base_name.replace("[", "").replace("]", "").replace("_", "")
                
                _update_csv_rows(
                    csv_path, 
                    r["cif"], 
                    [csv_row], 
                    "cif_file,normal_atoms,normal_formula,norm_duplicate,min_atoms,min_formula,min_duplicate\n"
                )
                _update_extxyz_collection(extxyz_path, clean_base, new_extxyz_frags)

        if is_dir and args.nproc > 1:
            print(f"Processing {len(cif_files)} CIF files using {args.nproc} processes...")
            with multiprocessing.Pool(args.nproc) as pool:
                for r in pool.imap_unordered(_process_mof_file, pool_args):
                    _flush_mof_result(r, seen_keys, csv_path, extxyz_path, is_dir=True)
        else:
            for pa in pool_args:
                _flush_mof_result(_process_mof_file(pa), seen_keys, csv_path, extxyz_path, is_dir=is_dir)
                
        print(f"Processing complete.")
        print(f"CSV summary updated at: {csv_path}")
        if os.path.exists(extxyz_path):
            print(f"ExtXYZ collection updated at: {extxyz_path}")

    elif args.kind == "cof":
        import os
        import glob
        import multiprocessing

        is_dir = os.path.isdir(args.input_path)
        if is_dir:
            cif_files = sorted(glob.glob(os.path.join(args.input_path, "*.cif")))
            out_dir = args.input_path
        else:
            if not os.path.exists(args.input_path):
                print(f"File not found: {args.input_path}")
                return
            cif_files = [args.input_path]
            out_dir = os.path.dirname(os.path.abspath(args.input_path))

        if not cif_files:
            print(f"No CIF files found to process.")
            return

        csv_path = os.path.join(out_dir, "fragmentation_summary.csv")
        extxyz_path = os.path.join(out_dir, "fragments_collection.extxyz")
        
        seen_keys = set()
        if is_dir:
            with open(csv_path, "w") as f:
                f.write("cif_file,normal_atoms,normal_formula,norm_duplicate,min_atoms,min_formula,min_duplicate\n")
            if os.path.exists(extxyz_path):
                os.remove(extxyz_path)
        else:
            if os.path.exists(extxyz_path):
                seen_keys = _load_seen_keys_from_extxyz(extxyz_path)

        pool_args = [(cif, args.radius, args.center, args.cof_layer) for cif in cif_files]

        def _flush_cof_result(r, seen_keys, csv_path, extxyz_path, is_dir):
            norm_dup = "no"
            min_dup = "no"
            if r['norm_res']:
                ckey = MOFFragmenter._chemical_identity_key(r['norm_res'].species, r['norm_res'].coords)
                if ckey in seen_keys: norm_dup = "yes"
            if r['min_res']:
                ckey = MOFFragmenter._chemical_identity_key(r['min_res'].species, r['min_res'].coords)
                if ckey in seen_keys: min_dup = "yes"

            new_extxyz_frags = []
            for key in ["norm_res", "min_res"]:
                res_obj = r[key]
                if res_obj:
                    ckey = MOFFragmenter._chemical_identity_key(res_obj.species, res_obj.coords)
                    if ckey in seen_keys:
                        continue
                    seen_keys.add(ckey)
                    suffix = "Min" if key == "min_res" else ""
                    base_name = r["cif"]
                    if base_name.endswith(".cif"): base_name = base_name[:-4]
                    clean_base = base_name.replace("[", "").replace("]", "").replace("_", "")
                    frag_name = f"{clean_base}FragCof{suffix}"
                    new_extxyz_frags.append((res_obj.species, res_obj.coords, frag_name))

            csv_row = f"{r['cif']},{r['norm_atoms']},{r['norm_formula']},{norm_dup},{r['min_atoms']},{r['min_formula']},{min_dup}\n"

            if is_dir:
                with open(csv_path, "a") as f:
                    f.write(csv_row)
                for species, coords, frag_name in new_extxyz_frags:
                    _write_extxyz(extxyz_path, species, coords, frag_name)
            else:
                base_name = r["cif"]
                if base_name.endswith(".cif"): base_name = base_name[:-4]
                clean_base = base_name.replace("[", "").replace("]", "").replace("_", "")
                
                _update_csv_rows(
                    csv_path, 
                    r["cif"], 
                    [csv_row], 
                    "cif_file,normal_atoms,normal_formula,norm_duplicate,min_atoms,min_formula,min_duplicate\n"
                )
                _update_extxyz_collection(extxyz_path, clean_base, new_extxyz_frags)

        if is_dir and args.nproc > 1:
            print(f"Processing {len(cif_files)} CIF files using {args.nproc} processes...")
            with multiprocessing.Pool(args.nproc) as pool:
                for r in pool.imap_unordered(_process_cof_file, pool_args):
                    _flush_cof_result(r, seen_keys, csv_path, extxyz_path, is_dir=True)
        else:
            for pa in pool_args:
                _flush_cof_result(_process_cof_file(pa), seen_keys, csv_path, extxyz_path, is_dir=is_dir)

        print(f"Processing complete.")
        print(f"CSV summary updated at: {csv_path}")
        if os.path.exists(extxyz_path):
            print(f"ExtXYZ collection updated at: {extxyz_path}")

    else:  # bio
        import os
        import glob
        import multiprocessing

        is_dir = os.path.isdir(args.input_path)
        if is_dir:
            pdb_files = sorted(glob.glob(os.path.join(args.input_path, "*.pdb")))
            out_dir = args.input_path
        else:
            if not os.path.exists(args.input_path):
                print(f"File not found: {args.input_path}")
                return
            pdb_files = [args.input_path]
            out_dir = os.path.dirname(os.path.abspath(args.input_path))

        if not pdb_files:
            print(f"No PDB files found to process.")
            return

        csv_path = os.path.join(out_dir, "bio_fragmentation_summary.csv")
        extxyz_path = os.path.join(out_dir, "bio_fragments_collection.extxyz")
        
        seen_keys = set()
        if is_dir:
            with open(csv_path, "w") as f:
                f.write("pdb_file,window_name,n_atoms,is_duplicate\n")
            if os.path.exists(extxyz_path):
                os.remove(extxyz_path)
        else:
            if os.path.exists(extxyz_path):
                seen_keys = _load_seen_keys_from_extxyz(extxyz_path)

        pool_args = [
            (pdb, args.window_size, args.stride, not args.no_pdbfixer, args.ph, args.chain)
            for pdb in pdb_files
        ]

        def _flush_bio_result(r, seen_keys, csv_path, extxyz_path, is_dir):
            base_name = r["pdb"]
            if base_name.lower().endswith(".pdb"): base_name = base_name[:-4]
            clean_base = base_name.replace("[", "").replace("]", "").replace("_", "")
            
            rows_to_write = []
            new_extxyz_frags = []
            
            for win_idx, res_obj in enumerate(r["results"]):
                win_name = f"{clean_base}W{win_idx:03d}"
                ckey = MOFFragmenter._chemical_identity_key(res_obj.species, res_obj.coords)
                if ckey in seen_keys:
                    rows_to_write.append(f"{r['pdb']},{win_name},{len(res_obj.species)},yes\n")
                else:
                    seen_keys.add(ckey)
                    rows_to_write.append(f"{r['pdb']},{win_name},{len(res_obj.species)},no\n")
                    new_extxyz_frags.append((res_obj.species, res_obj.coords, win_name))

            if is_dir:
                with open(csv_path, "a") as f:
                    for row in rows_to_write:
                        f.write(row)
                for species, coords, win_name in new_extxyz_frags:
                    _write_extxyz(extxyz_path, species, coords, win_name)
            else:
                _update_csv_rows(
                    csv_path, 
                    r["pdb"], 
                    rows_to_write, 
                    "pdb_file,window_name,n_atoms,is_duplicate\n"
                )
                _update_extxyz_collection(extxyz_path, clean_base, new_extxyz_frags)

        if is_dir and args.nproc > 1:
            print(f"Processing {len(pdb_files)} PDB files using {args.nproc} processes...")
            with multiprocessing.Pool(args.nproc) as pool:
                for r in pool.imap_unordered(_process_bio_file, pool_args):
                    _flush_bio_result(r, seen_keys, csv_path, extxyz_path, is_dir=True)
        else:
            for pa in pool_args:
                _flush_bio_result(_process_bio_file(pa), seen_keys, csv_path, extxyz_path, is_dir=is_dir)

        print(f"Processing complete.")
        print(f"CSV summary updated at: {csv_path}")
        if os.path.exists(extxyz_path):
            print(f"ExtXYZ collection updated at: {extxyz_path}")


if __name__ == "__main__":
    main()
