#!/usr/bin/env python3
"""Automated MOF fragment QA: charge and parity.

Usage:  python3 check_mof_fragments.py <run_dir>

<run_dir> holds fragments_collection.extxyz and the parent CIFs. Exits 0 when
every MUST passes, 1 otherwise.

MUST
  [M1] closed-shell formal charge == the recorded charge= tag (#9). The
       formal charge is metal oxidation states plus each non-metal molecule's
       charge, computed by UniFrag itself (MOFFragmenter.formal_charge_of), so
       this report and the fragmenter cannot disagree. A redox metal's state
       (Cu, Co, Mn, Fe, Cr, V) is read from the parent CIF in <run_dir>; a
       fragment that cannot be judged (None) fails too.
  [M2] even electron count at the recorded charge.
  [M3] no carboxyl C carrying two O-H, R-C(OH)2 (#10). Where the linker is
       conjugated, the formal-charge check alone reads two such ends as a
       closed-shell quinoid and passes them; 279 frames of the #9 output did.
  [M4] no H on a complete sp atom (#12): a C with exactly two non-metal heavy
       partners whose bonds are within 15 deg of linear (alkyne, nitrile C),
       or an N whose only heavy partner is such a C, at the shorter of its two
       bonds (nitrile N; the central N of a cut dicyanamide arm, N-C#N, is
       not one). No absolute length is used: disordered vinylene/imine CIFs
       draw C=C at 1.11-1.25 A with their own H at 115-134 deg, and those
       stay legal.
"""
import sys, os, re, glob, json, warnings
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentation_oop import MOFFragmenter, _parse_extxyz, _clean_base_of, _label_belongs_to


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    run = sys.argv[1]
    xyz = os.path.join(run, "fragments_collection.extxyz")
    if not os.path.exists(xyz):
        print(f"no fragments_collection.extxyz in {run}")
        sys.exit(1)
    cifs = {_clean_base_of(p): p for p in glob.glob(os.path.join(run, "*.cif"))}
    frames = _parse_extxyz(xyz)
    frag = MOFFragmenter()
    parent_of = {}
    charge_bad, odd, gem, sp_h = [], [], [], []
    cos_sp = np.cos(np.deg2rad(165.0))
    for fr in frames:
        lab, sp = fr["label"], fr["species"]
        rec = int(fr.get("charge", 0) or 0)
        base = next((b for b in cifs if _label_belongs_to(lab, b)), None)
        needs_parent = any(s in frag.METALS and s not in frag._OX_FIXED for s in sp)
        if needs_parent and base != parent_of.get("base"):
            frag.structure = None
            if base is not None:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        frag.structure = frag._load_clean_structure(cifs[base])
                except Exception:
                    frag.structure = None
            parent_of["base"] = base
        fq = frag.formal_charge_of(sp, np.asarray(fr["coords"], dtype=float))
        if fq != rec:
            charge_bad.append((lab, fq, rec))
        if (MOFFragmenter._parity_z_sum(sp) - rec) % 2:
            odd.append(lab)
        adj, _ = frag._nonmetal_graph(sp, np.asarray(fr["coords"], dtype=float))
        for c in adj:
            oxy = [j for j in adj[c] if sp[j] == "O"]
            if sp[c] == "C" and len(oxy) == 2 and all(
                    any(sp[k] == "H" for k in adj[o]) and sum(sp[k] != "H" for k in adj[o]) == 1
                    for o in oxy):
                gem.append(lab)
                break
        xyz_a = np.asarray(fr["coords"], dtype=float)
        heavy_nb = {a: [j for j in adj[a] if sp[j] != "H"] for a in adj}

        def linear_c(a):
            nb = heavy_nb.get(a, [])
            if sp[a] != "C" or len(nb) != 2:
                return False
            u, v = xyz_a[nb[0]] - xyz_a[a], xyz_a[nb[1]] - xyz_a[a]
            return round(float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v))), 9) <= cos_sp
        def nitrile_n(a):
            if sp[a] != "N" or len(heavy_nb[a]) != 1 or not linear_c(heavy_nb[a][0]):
                return False
            c = heavy_nb[a][0]
            other = [j for j in heavy_nb[c] if j != a][0]
            dist = lambda j: round(float(np.linalg.norm(xyz_a[j] - xyz_a[c])), 3)
            return dist(a) < dist(other)
        for a in adj:
            if sp[a] == "H" or not any(sp[k] == "H" for k in adj[a]):
                continue
            if linear_c(a) or nitrile_n(a):
                sp_h.append(lab)
                break

    report = {"frames": len(frames), "formal_charge_mismatch": len(charge_bad),
              "odd_electron": len(odd), "gem_diol_carboxyl": len(gem),
              "h_on_sp_atom": len(sp_h)}
    fails = []
    if charge_bad:
        fails.append(f"[M1] {len(charge_bad)} fragments whose formal charge is not the "
                     f"recorded charge (label, formal, recorded): {charge_bad[:8]}")
    if odd:
        fails.append(f"[M2] {len(odd)} odd-electron fragments: {odd[:8]}")
    if gem:
        fails.append(f"[M3] {len(gem)} fragments with a carboxyl C carrying two O-H: {gem[:8]}")
    if sp_h:
        fails.append(f"[M4] {len(sp_h)} fragments with H on an sp C or a nitrile N: {sp_h[:8]}")
    print("=" * 68)
    print("MOF FRAGMENT CHECK REPORT")
    print("=" * 68)
    print(json.dumps(report, indent=2))
    print("-" * 68)
    if fails:
        print(f"MUST FAILURES ({len(fails)}):")
        for f in fails:
            print("  FAIL " + f)
    else:
        print("MUST: all passed")
    print("=" * 68)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
