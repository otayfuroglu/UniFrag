#!/usr/bin/env python3
"""Automated QA for guest-loaded fragments (guest_insertion.py output).

Usage:  python3 check_guest_fragments.py <run_dir>

<run_dir> holds fragments_collection.extxyz, guests_config.json and the
fragments_guests_*.extxyz files. Exits 0 when every MUST passes, 1 otherwise.

MUST
  [G1] host atoms (species, coordinates) equal the parent frame's, to 1e-6 A.
  [G2] no guest atom closer to an atom outside its own molecule than
       r(Z_i) + r(Z_j); radii rebuilt from guests_config.json.
  [G3] every guest is rigid: its pair distances match its template to 1e-4 A.
  [G4] guest atoms are the set's unit repeated a whole number of times, in the
       counts the guests= tag gives (mixture counts = units x ratio; a single
       set holds only its own gas).
  [G5] charge= present and equal to the parent's; guest Z sum even, so the
       frame stays closed-shell (the host passed parity on its way into the
       collection).
  [G6] no "_" in any label; labels unique across all guest files; parent=
       names a frame of fragments_collection.extxyz.
"""
import os
import sys
import json

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

import guest_insertion as gi  # noqa: E402
from fragmentation_oop import _parse_extxyz  # noqa: E402

CHECKS = ("G1", "G2", "G3", "G4", "G5", "G6")


def _pd(x):
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
    return d[np.triu_indices(len(x), 1)]


def _frame_failures(f, gset, parents, table):
    """Check ids this frame fails."""
    from ase.data import atomic_numbers
    bad = set()
    t = gi.parse_header_tags(f["comment"])
    label = t.get("label", "")
    parent = parents.get(t.get("parent"))
    if "_" in label or parent is None:
        bad.add("G6")
    if "charge" not in t or parent is None or int(t["charge"]) != int(parent["charge"]):
        bad.add("G5")
    sp = list(f["species"])
    co = np.asarray(f["coords"], float)
    try:
        n = int(t["n_host"])
    except (KeyError, ValueError):
        return bad | {"G1", "G4"}
    if parent is None or sp[:n] != list(parent["species"]) or \
            not np.allclose(co[:n], np.asarray(parent["coords"], float), atol=1e-6, rtol=0):
        bad.add("G1")
    gsp = sp[n:]
    if sum(atomic_numbers[s] for s in gsp) % 2:
        bad.add("G5")

    # G4: the guests= counts are units x ratio, and the atoms follow that order
    try:
        counts = {k: int(v) for k, v in (x.split(":") for x in t.get("guests", "").split())}
    except ValueError:
        return bad | {"G4"}
    names = [g.name for g in gset.guests]
    if set(counts) - set(names):
        return bad | {"G4"}
    per = {counts.get(nm, 0) / r for nm, r in zip(names, gset.ratios)}
    units = per.pop() if len(per) == 1 else None
    if units is None or units != int(units) or units < 1:
        return bad | {"G4"}
    unit = gset.unit()
    order = [g for _ in range(int(units)) for g in unit]
    if [s for g in order for s in g.symbols] != gsp:
        return bad | {"G4"}

    # G2 / G3 per guest molecule
    r = gi.radii_for(sp, table)
    start = n
    for g in order:
        idx = np.arange(start, start + len(g.symbols))
        start += len(g.symbols)
        if len(idx) > 1 and not np.allclose(_pd(co[idx]), _pd(g.positions), atol=1e-4, rtol=0):
            bad.add("G3")
        others = np.setdiff1d(np.arange(len(sp)), idx)
        d = np.linalg.norm(co[idx][:, None] - co[others][None], axis=-1)
        if (d < r[idx][:, None] + r[others][None] - 1e-6).any():
            bad.add("G2")
    return bad


def check(run_dir):
    """{"frames": n, "G1": [labels], ...}"""
    with open(os.path.join(run_dir, gi.CONFIG_JSON)) as fh:
        conf = json.load(fh)
    sh = conf["shared"]
    table = gi.radius_table(sh["vdw_radii"], sh["vdw_shift"], sh["vdw_scale"])
    parents = {f["label"]: f for f in
               _parse_extxyz(os.path.join(run_dir, "fragments_collection.extxyz"))}
    report = {"frames": 0, **{c: [] for c in CHECKS}}
    seen = set()
    for gset in gi.sets_from_specs(conf["set_specs"]):
        path = os.path.join(run_dir, gset.output)
        for f in _parse_extxyz(path):
            report["frames"] += 1
            bad = _frame_failures(f, gset, parents, table)
            if f["label"] in seen:          # labels are OpenMLP keys: unique across files
                bad.add("G6")
            seen.add(f["label"])
            for c in bad:
                report[c].append(f["label"])
    return report


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    rep = check(sys.argv[1])
    fails = [c for c in CHECKS if rep[c]]
    print("=" * 68)
    print("GUEST FRAGMENT CHECK REPORT")
    print("=" * 68)
    print(json.dumps({"frames": rep["frames"], **{c: len(rep[c]) for c in CHECKS}}, indent=2))
    print("-" * 68)
    if fails:
        print(f"MUST FAILURES ({len(fails)}):")
        for c in fails:
            print(f"  FAIL [{c}] {len(rep[c])} frames: {rep[c][:8]}")
    else:
        print("MUST: all passed")
    print("=" * 68)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
