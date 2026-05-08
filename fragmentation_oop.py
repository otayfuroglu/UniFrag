import argparse
import os
from collections import deque
from dataclasses import dataclass
import networkx as nx
from pathlib import Path

import numpy as np
from pymatgen.core import Molecule, Structure
from pymatgen.io.cif import CifParser


@dataclass
class FragmentResult:
    species: list
    coords: list


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

    def _warn_possible_impurity_components(self, struct, framework_kind="framework"):
        """Warn when CIF appears to contain detached small ionic/solvent components."""
        try:
            symbols = []
            for site in struct:
                try:
                    symbols.append(site.specie.symbol)
                except Exception:
                    symbols.append(max(site.species.items(), key=lambda kv: kv[1])[0].symbol)

            if len(symbols) < 2:
                return

            all_neigh = struct.get_all_neighbors(r=2.6)
            graph = [[] for _ in range(len(struct))]
            for i, neighs in enumerate(all_neigh):
                si = symbols[i]
                ci = struct[i].coords
                for n in neighs:
                    j = n.index
                    if j <= i:
                        continue
                    d = float(np.linalg.norm(ci - n.coords))
                    if self.is_valid_bond(si, symbols[j], d):
                        graph[i].append(j)
                        graph[j].append(i)

            comps = []
            seen = set()
            for i in range(len(symbols)):
                if i in seen:
                    continue
                q = deque([i])
                seen.add(i)
                comp = {i}
                while q:
                    u = q.popleft()
                    for v in graph[u]:
                        if v not in seen:
                            seen.add(v)
                            comp.add(v)
                            q.append(v)
                comps.append(comp)

            if len(comps) <= 1:
                return

            metals = getattr(self, "METALS", set())

            def comp_meta(comp):
                idx = list(comp)
                sp = [symbols[k] for k in idx]
                heavy = sum(1 for x in sp if x != "H")
                metal = sum(1 for x in sp if x in metals)
                return {
                    "idx": idx,
                    "atoms": len(idx),
                    "heavy": heavy,
                    "metal": metal,
                    "sp": sp,
                }

            meta = [comp_meta(c) for c in comps]
            largest_heavy = max((m["heavy"] for m in meta), default=0)
            if largest_heavy <= 0:
                return

            suspicious = []
            for m in meta:
                if m["heavy"] == largest_heavy:
                    continue
                is_small = m["heavy"] <= 14 or m["heavy"] <= int(0.25 * largest_heavy)
                is_guest_like = m["metal"] == 0 and m["atoms"] <= 60
                if is_small and is_guest_like:
                    suspicious.append(m)

            if not suspicious:
                return

            preview = []
            for m in sorted(suspicious, key=lambda x: (x["heavy"], x["atoms"]))[:3]:
                from collections import Counter
                c = Counter(m["sp"])
                formula = " ".join(f"{el}{c[el]}" for el in sorted(c))
                preview.append(f"{m['atoms']} atoms ({formula})")

            print(f"WARNING: {framework_kind} CIF has {len(comps)} disconnected components; {len(suspicious)} look like solvent/ionic guests.")
            print("  -> Guest-like components detected: " + "; ".join(preview))
            print("  -> Fragment quality can be affected. Consider cleaning/removing guests before fragmentation.")
        except Exception:
            # Warning path must never break fragmentation.
            return

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
            if parent is None or species[parent] != "C":
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

            bl = self.cap_bond_length("C")
            new_h = p + bl * dvec

            # quick collision guard with other atoms
            clash = False
            for j, sp in enumerate(species):
                if j in (hidx, parent):
                    continue
                dj = np.linalg.norm(new_h - np.array(coords[j], dtype=float))
                if sp == "H" and dj < 1.4:
                    clash = True
                    break
                if sp != "H" and dj < 0.9:
                    clash = True
                    break
            if not clash:
                coords[hidx] = new_h


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
            if has_metal or len(heavy_neighbors) != 1 or species[heavy_neighbors[0]] != "C":
                continue
            base = opos - np.array(coords[heavy_neighbors[0]], dtype=float)
            before = len(species)
            self.place_capping_h(i, base, self.cap_bond_length("O"), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)
            if len(species) == before and np.linalg.norm(base) > 1e-12:
                self.place_capping_h(i, -base, self.cap_bond_length("O"), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

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

        def shortest_path_without_edge(start, target, banned_edge):
            q = deque([start])
            parent = {start: None}
            while q:
                u = q.popleft()
                if u == target:
                    break
                for v in hadj.get(u, []):
                    if {u, v} == banned_edge or v in parent:
                        continue
                    parent[v] = u
                    q.append(v)
            if target not in parent:
                return None
            path = []
            u = target
            while u is not None:
                path.append(u)
                u = parent[u]
            return list(reversed(path))

        def distance_from_bridge(atoms):
            q = deque((seed, 0) for seed in bridge_atoms)
            seen = set(bridge_atoms)
            while q:
                u, dist = q.popleft()
                if u in atoms:
                    return dist
                for v in hadj.get(u, []):
                    if v not in seen:
                        seen.add(v)
                        q.append((v, dist + 1))
            return 999

        def find_first_cycle(max_cycle=8):
            cycles = []
            seen_cycles = set()
            for u in heavy:
                for v in hadj.get(u, []):
                    if u >= v:
                        continue
                    path = shortest_path_without_edge(u, v, {u, v})
                    if not path:
                        continue
                    cycle = set(path)
                    if len(cycle) < 3 or len(cycle) > max_cycle:
                        continue
                    key = tuple(sorted(cycle))
                    if key in seen_cycles:
                        continue
                    seen_cycles.add(key)
                    cycles.append(cycle)
            if not cycles:
                return set()
            cycles.sort(key=lambda cyc: (distance_from_bridge(cyc), len(cyc), sum(sorted(cyc))))
            return cycles[0]

        ring = find_first_cycle()
        if not ring:
            q = deque(bridge_atoms)
            seen = set(bridge_atoms)
            fallback = []
            while q and len(fallback) < 6:
                u = q.popleft()
                fallback.append(u)
                for v in hadj.get(u, []):
                    if v not in seen:
                        seen.add(v)
                        q.append(v)
            ring = set(fallback)
        if not ring:
            return [], [], []

        keep_heavy = set(bridge_atoms) | set(ring)
        # Keep the shortest connector from node-bound atoms into the first ring.
        q = deque(bridge_atoms)
        parent = {i: None for i in bridge_atoms}
        target = None
        while q:
            u = q.popleft()
            if u in ring:
                target = u
                break
            for v in hadj.get(u, []):
                if v not in parent:
                    parent[v] = u
                    q.append(v)
        while target is not None:
            keep_heavy.add(target)
            target = parent[target]

        # Preserve hetero atoms attached to kept connector/ring atoms on the node side.
        for i in heavy:
            if i in keep_heavy:
                continue
            if linker_species[i] in {"O", "N", "S", "P", "B"}:
                if any(nb in keep_heavy for nb in hadj.get(i, [])) and any(nb in bridge_atoms for nb in hadj.get(i, []) + [i]):
                    keep_heavy.add(i)

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
        selected_images = scored_images[:1] if minimize else scored_images
        partial_images = scored_images[1:] if minimize else []

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

        mol = Molecule(species, coords)
        mol.to(filename=output_path, fmt="xyz")
        print("  -> MOF Path J (moffragmentor node+linker combine).")
        print(f"  -> Helper fragments available in mof_nodes_lib/ and mof_linkers_lib/.")
        print(f"Final size: {len(species)} atoms.")
        print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)

    def extract(self, mof_path, center_idx=-1, nmetals=3, output_path="fragment.xyz", minimize=False):
        print(f"Loading '{mof_path}'...")
        structure_stem = Path(mof_path).stem
        if minimize:
            probe_output = str(Path("/private/tmp") / f"{self._safe_name(structure_stem)}_normal_probe.xyz")
            normal_probe = self.extract(mof_path, center_idx=center_idx, nmetals=nmetals, output_path=probe_output, minimize=False)
            if normal_probe is not None and len(normal_probe.species) < 80:
                Molecule(normal_probe.species, normal_probe.coords).to(filename=output_path, fmt="xyz")
                print(f"  -> Normal fragment has {len(normal_probe.species)} atoms (<80); using normal as minimized fragment.")
                print(f"Final size: {len(normal_probe.species)} atoms.")
                print(f"Saved: {output_path}")
                return normal_probe

        use_candidate_count_path = (not minimize) and self._path_j_second_node_family(structure_stem)
        if use_candidate_count_path:
            print("  -> Skipping MOF Path J so normal PCN/NU can use two-node atom-count selection.")
        else:
            combined = self._try_moffragmentor_node_linker_fragment(mof_path, output_path, minimize=minimize)
            if combined is not None:
                return combined
        struct = Structure.from_file(mof_path)
        self._warn_possible_impurity_components(struct, framework_kind="MOF")
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

        self.optimize_capped_h_geometry_only(final_species, final_coords, getattr(self, "_last_capped_h_indices", None))
        print(f"Final size: {len(final_species)} atoms.")
        mol = Molecule(final_species, final_coords)
        mol.to(filename=output_path, fmt="xyz")
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

        if minimize:
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

                remove_cond = len(touching_metals) < 1 if zif_mode else len(touching_metals) < 2
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
                    keep_atoms = set(bridge_atoms)
                    q = deque([(ba, 0) for ba in bridge_atoms])
                    while q:
                        curr, depth = q.popleft()
                        if depth < 5:
                            for nb in local_adj[curr]:
                                if nb in comp and nb not in keep_atoms:
                                    keep_atoms.add(nb)
                                    q.append((nb, depth + 1))

                    changed = True
                    while changed:
                        changed = False
                        to_remove = set()
                        for ka in keep_atoms - bridge_atoms:
                            if supercell[ka].species_string == "C":
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
                remove_cond = len(touching_metals) < 1 if zif_mode else len(touching_metals) < 2
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

    def _limit_to_single_helical_period(self, scored_images, node_ctr, lattice):
        """
        For helical COFs, keep linker images within ~one lattice period along
        the dominant lattice axis to avoid over-extended Path J assemblies.
        """
        if len(scored_images) < 3:
            return scored_images

        centers = []
        for _, _, linker, shift in scored_images:
            lco = [np.array(c, dtype=float) + shift for c in linker.molecule.cart_coords]
            centers.append(np.mean(lco, axis=0))
        centers = np.array(centers, dtype=float)
        if centers.size == 0:
            return scored_images

        axis_idx = None
        best_span = 0.0
        for i in range(3):
            vec = np.array(lattice[i], dtype=float)
            nrm = float(np.linalg.norm(vec))
            if nrm < 1e-6:
                continue
            u = vec / nrm
            proj = centers @ u
            span = float(np.max(proj) - np.min(proj))
            if span > best_span:
                best_span = span
                axis_idx = i

        if axis_idx is None:
            return scored_images

        axis = np.array(lattice[axis_idx], dtype=float)
        axis_len = float(np.linalg.norm(axis))
        if axis_len < 1e-6:
            return scored_images
        if best_span <= 1.05 * axis_len:
            return scored_images

        u = axis / axis_len
        node_proj = float(np.dot(node_ctr, u))
        kept = []
        for entry, ctr in zip(scored_images, centers):
            dproj = float(abs(np.dot(ctr, u) - node_proj))
            if dproj <= 0.55 * axis_len:
                kept.append(entry)

        if len(kept) >= 1:
            print(f"  -> COF Path J helical limiter kept {len(kept)}/{len(scored_images)} linker images (~1 period).")
            return kept
        return scored_images

    def _coffragmentor_retry_clean_cif(self, cif_path):
        """Try to sanitize problematic CIF occupancy/disorder for coffragmentor."""
        try:
            parser = CifParser(cif_path, occupancy_tolerance=2.1)
            structs = parser.get_structures(primitive=False)
            if not structs:
                return None
            cleaned = structs[0]
            tmp_path = Path('/private/tmp') / f"{Path(cif_path).stem}_cof_clean_for_pathj.cif"
            cleaned.to(filename=str(tmp_path), fmt='cif')
            return str(tmp_path)
        except Exception:
            return None

    def _try_strict_helix_node_linker_path(self, struct, output_path, minimize=False):
        """Strict helix COF rule: detect exactly two repeating building blocks
        (node/linker) and build fragment only by node + attached linker combine."""
        n = len(struct)
        dbg = os.environ.get("UNIFRAG_DEBUG_HELIX", "0") == "1"
        if dbg:
            print(f"  [helix] atoms={n}")
        if n < 8:
            if dbg:
                print("  [helix] reject: too few atoms")
            return None

        sym = []
        for site in struct:
            try:
                sym.append(site.specie.symbol)
            except Exception:
                sym.append(max(site.species.items(), key=lambda kv: kv[1])[0].symbol)

        all_neigh = struct.get_all_neighbors(r=2.4)
        graph = [[] for _ in range(n)]
        for i, neighs in enumerate(all_neigh):
            ci = struct[i].coords
            si = sym[i]
            for nn in neighs:
                j = nn.index
                if j <= i:
                    continue
                d = float(np.linalg.norm(ci - nn.coords))
                if self.is_valid_bond(si, sym[j], d):
                    graph[i].append(j)
                    graph[j].append(i)

        heavy_deg = {}
        for i in range(n):
            heavy_deg[i] = sum(1 for j in graph[i] if sym[j] != 'H')

        def in_small_ring(u, v, max_len=8):
            q = deque([(u, 0)])
            seen = {u}
            while q:
                x, dep = q.popleft()
                if dep >= max_len:
                    continue
                for y in graph[x]:
                    if (x == u and y == v) or (x == v and y == u):
                        continue
                    if y == v:
                        return True
                    if y not in seen:
                        seen.add(y)
                        q.append((y, dep + 1))
            return False

        cut_edges = set()
        for i in range(n):
            if sym[i] != 'C':
                continue
            for j in graph[i]:
                if j <= i or sym[j] != 'C':
                    continue
                if heavy_deg[i] >= 3 and heavy_deg[j] >= 3 and (not in_small_ring(i, j)):
                    cut_edges.add((i, j))

        if dbg:
            print(f"  [helix] cut_edges={len(cut_edges)}")
        if not cut_edges:
            if dbg:
                print("  [helix] reject: no cut edges")
            return None

        kept_graph = [[] for _ in range(n)]
        for i in range(n):
            for j in graph[i]:
                a, b = (i, j) if i < j else (j, i)
                if (a, b) in cut_edges:
                    continue
                kept_graph[i].append(j)

        comp_id = [-1] * n
        comps = []
        cid = 0
        for i in range(n):
            if comp_id[i] != -1:
                continue
            q = deque([i])
            comp_id[i] = cid
            comp = {i}
            while q:
                u = q.popleft()
                for v in kept_graph[u]:
                    if comp_id[v] == -1:
                        comp_id[v] = cid
                        comp.add(v)
                        q.append(v)
            comps.append(comp)
            cid += 1

        if dbg:
            print(f"  [helix] comps={len(comps)}")
        if len(comps) < 3:
            if dbg:
                print("  [helix] reject: too few components")
            return None

        def comp_hash(comp):
            nodes = sorted(comp)
            idx = {g: k for k, g in enumerate(nodes)}
            g2 = nx.Graph()
            for g in nodes:
                g2.add_node(idx[g], sp=sym[g])
            for g in nodes:
                for h in kept_graph[g]:
                    if h in comp and g < h:
                        g2.add_edge(idx[g], idx[h])
            return nx.weisfeiler_lehman_graph_hash(g2, node_attr='sp')

        hash_to_comps = {}
        for k, comp in enumerate(comps):
            h = comp_hash(comp)
            hash_to_comps.setdefault(h, []).append(k)

        if dbg:
            print(f"  [helix] hash classes={len(hash_to_comps)}")
        if len(hash_to_comps) < 2:
            if dbg:
                print("  [helix] reject: <2 hash classes")
            return None

        # Helix strict rule: use the two dominant repeating building blocks
        # and ignore minor residual classes from disorder/symmetry artifacts.
        ranked_classes = sorted(
            hash_to_comps.items(),
            key=lambda kv: (len(kv[1]), np.mean([len(comps[c]) for c in kv[1]]) if kv[1] else 0.0),
            reverse=True,
        )
        major = ranked_classes[:2]
        hash_to_comps = {k: v for k, v in major}

        comp_graph = {k: set() for k in range(len(comps))}
        for i, j in cut_edges:
            ci, cj = comp_id[i], comp_id[j]
            if ci != cj:
                comp_graph[ci].add(cj)
                comp_graph[cj].add(ci)

        groups = list(hash_to_comps.values())
        g0_deg = np.mean([len(comp_graph[c]) for c in groups[0]]) if groups[0] else 0.0
        g1_deg = np.mean([len(comp_graph[c]) for c in groups[1]]) if groups[1] else 0.0
        node_group = groups[0] if g0_deg >= g1_deg else groups[1]
        linker_group = groups[1] if node_group is groups[0] else groups[0]

        ctr = np.array(struct.lattice.get_cartesian_coords([0.5, 0.5, 0.5]), dtype=float)
        def comp_center(comp_idx):
            pts = [np.array(struct[i].coords, dtype=float) for i in comps[comp_idx]]
            return np.mean(pts, axis=0)

        node_comp = min(node_group, key=lambda c: float(np.linalg.norm(comp_center(c) - ctr)))
        linker_set = set(linker_group)
        attached_linkers = [c for c in comp_graph[node_comp] if c in linker_set]
        if dbg:
            print(f"  [helix] node_group={len(node_group)} linker_group={len(linker_group)}")
            print(f"  [helix] attached before fallback={len(attached_linkers)}")
        if not attached_linkers:
            # If class filtering breaks direct adjacency, pick nearest linker block.
            candidates = list(linker_set)
            if not candidates:
                if dbg:
                    print("  [helix] reject: no linker candidates")
                return None
            nc = comp_center(node_comp)
            nearest = min(candidates, key=lambda c: float(np.linalg.norm(comp_center(c) - nc)))
            attached_linkers = [nearest]

        if minimize:
            attached_linkers = attached_linkers[:1]

        keep_comps = {node_comp} | set(attached_linkers)
        keep_atoms = set()
        for c in keep_comps:
            keep_atoms |= set(comps[c])

        species = [sym[i] for i in sorted(keep_atoms)]
        coords = [np.array(struct[i].coords, dtype=float) for i in sorted(keep_atoms)]

        # Cap open cut sites after node+linker combine.
        local = {g: k for k, g in enumerate(sorted(keep_atoms))}
        capped_h_flags = [False] * len(species)
        for i, j in cut_edges:
            gi = i if (i in keep_atoms and j not in keep_atoms) else (j if (j in keep_atoms and i not in keep_atoms) else None)
            gj = j if gi == i else (i if gi == j else None)
            if gi is None or gj is None:
                continue
            li = local.get(gi)
            if li is None:
                continue
            sp = species[li]
            if sp not in {'C', 'B', 'O', 'N'}:
                continue
            base = np.array(coords[li]) - np.array(struct[gj].coords, dtype=float)
            self.place_capping_h(li, base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        capped_h_indices = [i for i, f in enumerate(capped_h_flags) if f and species[i] == 'H']
        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)

        mol = Molecule(species, coords)
        mol.to(filename=output_path, fmt='xyz')
        print('  -> COF Path J (strict helix node+linker).')
        print(f'Final size: {len(species)} atoms.')
        print(f'Saved: {output_path}')
        return FragmentResult(species=species, coords=coords)

    def _try_coffragmentor_path_j(self, cif_path, output_path, minimize=False):
        try:
            from coffragmentor import COF
        except Exception as exc:
            print(f"  coffragmentor unavailable: {exc}")
            return None
        try:
            struct = Structure.from_file(cif_path)
            self._warn_possible_impurity_components(struct, framework_kind="COF")
            # First-class strict helix rule: if we can detect two repeating
            # building blocks (node/linker), always build by direct node+linker
            # combine and return immediately.
            strict = self._try_strict_helix_node_linker_path(struct, output_path, minimize=minimize)
            if strict is not None:
                return strict
            result = COF.from_cif(cif_path).fragment()
        except Exception as exc:
            clean_cif = self._coffragmentor_retry_clean_cif(cif_path)
            if clean_cif is None:
                print(f"  coffragmentor COF Path J failed: {exc}")
                return None
            try:
                struct = Structure.from_file(clean_cif)
                self._warn_possible_impurity_components(struct, framework_kind="COF")
                strict = self._try_strict_helix_node_linker_path(struct, output_path, minimize=minimize)
                if strict is not None:
                    return strict
                result = COF.from_cif(clean_cif).fragment()
                print("  -> COF Path J recovered via cleaned CIF retry.")
            except Exception as exc2:
                print(f"  coffragmentor COF Path J failed: {exc2}")
                return None

        nodes = list(getattr(result, "nodes", []))
        linkers = list(getattr(result, "linkers", []))
        if not nodes or not linkers:
            strict = self._try_strict_helix_node_linker_path(struct, output_path, minimize=minimize)
            if strict is not None:
                return strict
            return None

        metals = {"Zn", "Cu", "Fe", "Co", "Ni", "Mn"}
        cell_ctr = np.array(struct.lattice.get_cartesian_coords([0.5, 0.5, 0.5]), dtype=float)

        def node_score(sbu):
            sp = [str(x) for x in sbu.molecule.species]
            return (sum(1 for x in sp if x in metals), sp.count("N"), len(sp))

        pc_nodes = [sbu for sbu in nodes if node_score(sbu)[0] >= 1 and node_score(sbu)[1] >= 4]
        is_metallo_pc = bool(pc_nodes)
        if is_metallo_pc:
            node = max(pc_nodes, key=node_score)
            node_kind = "metallo-PC"
        else:
            node = min(
                nodes,
                key=lambda sbu: float(np.linalg.norm(np.mean(sbu.molecule.cart_coords, axis=0) - cell_ctr)),
            )
            node_kind = "general"

        node_sp = [str(x) for x in node.molecule.species]
        node_co = [np.array(c, dtype=float) for c in node.molecule.cart_coords]
        node_ctr = np.mean(node_co, axis=0)

        lattice = np.array(struct.lattice.matrix, dtype=float)
        image_vectors = []
        for ia in (-1, 0, 1):
            for ib in (-1, 0, 1):
                for ic in (-1, 0, 1):
                    image_vectors.append(ia * lattice[0] + ib * lattice[1] + ic * lattice[2])

        scored_images = []
        for linker_idx, linker in enumerate(linkers):
            lsp = [str(x) for x in linker.molecule.species]
            lco0 = [np.array(c, dtype=float) for c in linker.molecule.cart_coords]
            best = None
            for shift in image_vectors:
                lco = [c + shift for c in lco0]
                contact_pairs = 0
                bo_contacts = 0
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
                            contact_pairs += 1
                            if {si, sj} == {"B", "O"}:
                                bo_contacts += 1
                if contact_pairs <= 0 or min_d > 2.2:
                    continue
                lctr = np.mean(lco, axis=0)
                score = (contact_pairs, -min_d, -float(np.linalg.norm(lctr - node_ctr)))
                entry = (score, linker_idx, linker, shift)
                if is_metallo_pc:
                    if bo_contacts > 0:
                        scored_images.append(entry)
                else:
                    if best is None or score > best[0]:
                        best = entry
            if (not is_metallo_pc) and best is not None:
                scored_images.append(best)

        if not scored_images:
            return None

        scored_images.sort(key=lambda item: item[0], reverse=True)
        if (not is_metallo_pc) and (not minimize):
            scored_images = self._limit_to_single_helical_period(scored_images, node_ctr, lattice)
        selected_images = scored_images[:1] if minimize else scored_images

        species = list(node_sp)
        coords = [np.array(c, dtype=float) for c in node_co]
        for _, _, linker, shift in selected_images:
            lsp = [str(x) for x in linker.molecule.species]
            lco = [np.array(c, dtype=float) + shift for c in linker.molecule.cart_coords]
            species.extend(lsp)
            coords.extend(lco)

        base_species = list(species)
        base_coords = [np.array(c, dtype=float) for c in coords]
        layer_mode = self.layer_mode
        if layer_mode == "auto":
            layer_mode = "dimer" if is_metallo_pc else "monomer"

        if layer_mode == "dimer":
            stack_vec = min(lattice, key=lambda v: float(np.linalg.norm(v)))
            species = base_species + base_species
            coords = base_coords + [c + stack_vec for c in base_coords]
        else:
            species = base_species
            coords = base_coords

        mol = Molecule(species, coords)
        mol.to(filename=output_path, fmt="xyz")
        print(f"  -> COF Path J ({node_kind} node+linker {layer_mode}).")
        print(f"Final size: {len(species)} atoms.")
        print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)

    def extract(self, cif_path, output_path="cof_fragment.xyz", center_idx=-1, minimize=False):
        print(f"Loading '{cif_path}'...")
        direct = self._try_coffragmentor_path_j(cif_path, output_path, minimize=minimize)
        if direct is not None:
            return direct
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
        self._warn_possible_impurity_components(struct, framework_kind="COF")
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
                for i in center_list:
                    ci = supercell[i].coords
                    for j in comp:
                        d = np.linalg.norm(ci - supercell[j].coords)
                        if d < min_d:
                            min_d = d
                layered_candidates.append((min_d, comp))

            node_component_by_atom = {i: comp for comp in node_comps for i in comp}

            layered_component_count = 1
            if path_mode == "A" and layered_candidates:
                layered_candidates.sort(key=lambda x: x[0])
                best_d, best_comp = layered_candidates[0]
                if 1.0 <= best_d <= 4.5:
                    # Some stacked COFs have two symmetry-equivalent neighboring
                    # node components at the same interlayer spacing. Keep all
                    # tied neighbors so one equivalent linker is not dropped.
                    same_layer = [comp for d, comp in layered_candidates if abs(d - best_d) <= 0.05]
                    for comp in same_layer:
                        core_nodes |= set(comp)
                    layered_component_count += len(same_layer)
                    path_mode = "B"

            if path_mode == "B":
                layer_mode = "dimer" if self.layer_mode == "auto" else self.layer_mode
                print(f"  -> COF Path B (Layered set {layer_mode}). Node components: {layered_component_count}, spacing: {best_d:.2f} A")
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
                    if path_mode == "B" and not minimize:
                        terminal_node = set(node_component_by_atom.get(v, {v}))
                        final |= terminal_node
                        visited |= terminal_node
                    else:
                        broken.append((u, np.array(unwrapped[v]) - np.array(unwrapped[u])))
                    continue
                if v not in sphere:
                    broken.append((u, np.array(unwrapped[v]) - np.array(unwrapped[u])))
                    continue
                if v in visited:
                    continue
                visited.add(v)
                final.add(v)
                queue.append(v)

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

                # Choose one whole linker branch to keep fully. For Path B
                # dimers, do this once per disconnected layer so one layer does
                # not get the full linker while the other is only ring-trimmed.
                primary_indices = set()
                if comp_info:
                    layer_mode = "dimer" if self.layer_mode == "auto" else self.layer_mode
                    if path_mode == "B" and layer_mode != "monomer":
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

                # General COF rule: do not cut rings in minimized fragments.
                # Keep all non-primary linker components intact.
                for comp_idx, (comp, touching_nodes, bridge_atoms) in enumerate(comp_info):
                    if comp_idx in primary_indices:
                        continue
                    keep_heavy |= set(comp)

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
                    layer_mode = "dimer" if self.layer_mode == "auto" else self.layer_mode
                    keep_n = 1 if layer_mode == "monomer" else 2
                    keep = set().union(*ranked[:keep_n])
                else:
                    ranked = sorted(comps, key=comp_key, reverse=True)
                    # General COF fallback rule: keep two principal building
                    # blocks in normal mode when disconnected components exist.
                    keep_n = 2 if (not minimize and len(ranked) >= 2) else 1
                    keep = set().union(*ranked[:keep_n])
                species = [x for i, x in enumerate(species) if i in keep]
                coords = [x for i, x in enumerate(coords) if i in keep]
                capped_h_flags = [x for i, x in enumerate(capped_h_flags) if i in keep]

        # Final edge-site capping after all pruning: ensure terminal B/O sites
        # at fragment boundaries are protonated/capped when undercoordinated.
        if species:
            heavy_idx = [i for i, sp in enumerate(species) if sp != "H"]
            hadj_bo = {i: [] for i in heavy_idx}
            for a in range(len(heavy_idx)):
                i = heavy_idx[a]
                ci = np.array(coords[i])
                for b in range(a + 1, len(heavy_idx)):
                    j = heavy_idx[b]
                    d = np.linalg.norm(ci - np.array(coords[j]))
                    if self.is_valid_bond(species[i], species[j], d):
                        hadj_bo[i].append(j)
                        hadj_bo[j].append(i)

            for i, sp in list(enumerate(species)):
                if sp not in {"B", "O"}:
                    continue

                cpos = np.array(coords[i])
                heavy_nbs = hadj_bo.get(i, [])
                h_nbs = 0
                for j, sj in enumerate(species):
                    if sj != "H":
                        continue
                    d = np.linalg.norm(cpos - np.array(coords[j]))
                    hcut = 1.5 if sp == "O" else 1.7
                    if d <= hcut:
                        h_nbs += 1

                if sp == "O":
                    if self.oxygen_already_protonated(i, species, coords):
                        continue
                    if len(heavy_nbs) >= 2:
                        continue
                    n_cap = min(1, max(0, 2 - len(heavy_nbs) - h_nbs))
                else:
                    if len(heavy_nbs) >= 3:
                        continue
                    n_cap = min(1, max(0, 3 - len(heavy_nbs) - h_nbs))

                if n_cap <= 0:
                    continue

                base = np.zeros(3)
                if heavy_nbs:
                    for nb in heavy_nbs:
                        base -= (np.array(coords[nb]) - cpos)
                else:
                    base = np.array([1.0, 0.0, 0.0])

                for _ in range(n_cap):
                    self.place_capping_h(i, base, self.cap_bond_length(sp), species, coords, min_hh=1.5, capped_h_flags=capped_h_flags)

        capped_h_indices = [i for i, is_cap in enumerate(capped_h_flags) if is_cap and species[i] == "H"]
        self.optimize_capped_h_geometry_only(species, coords, capped_h_indices)
        print(f"Final size: {len(species)} atoms")
        mol = Molecule(species, coords)
        mol.to(filename=output_path, fmt="xyz")
        print(f"Saved: {output_path}")
        return FragmentResult(species=species, coords=coords)
def main():
    parser = argparse.ArgumentParser(description="OOP fragmenter core for MOF/COF")
    parser.add_argument("cif_path")
    parser.add_argument("--kind", choices=["mof", "cof"], default="mof")
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--nmetals", type=int, default=3)
    parser.add_argument("--output", default="fragment.xyz")
    parser.add_argument("--minimize", action="store_true")
    parser.add_argument(
        "--cof-layer",
        choices=["auto", "monomer", "dimer"],
        default="auto",
        help="COF layer output for layered COFs. auto preserves existing dimer behavior.",
    )
    args = parser.parse_args()

    if args.kind == "mof":
        frag = MOFFragmenter(radius=args.radius)
        frag.extract(args.cif_path, center_idx=args.center, nmetals=args.nmetals, output_path=args.output, minimize=args.minimize)
    else:
        frag = COFFragmenter(radius=args.radius, layer_mode=args.cof_layer)
        frag.extract(args.cif_path, center_idx=args.center, output_path=args.output, minimize=args.minimize)


if __name__ == "__main__":
    main()
