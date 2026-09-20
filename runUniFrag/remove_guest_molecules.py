#!/usr/bin/env python3
"""Strip guest / solvent molecules from COF CIFs, producing a clean collection.

A COF framework is periodic: following its bonds you can walk out of the unit
cell and come back to the same atom in a different lattice image. A guest -
acetone, DMF, water - is a finite molecule that closes on itself with zero net
translation. That difference, not size or composition, is what separates them,
and it is what this script tests.

Originals are never modified. Cleaned copies are written to --out.

Usage:
    python3 remove_guest_molecules.py <in_dir> --out <out_dir> [--report r.csv]
                                      [--csd-check N] [--jobs N]
"""
import argparse, os, glob, csv, collections, warnings, sys
import numpy as np

warnings.filterwarnings("ignore")


def _prune_overbonded_h(sg, struct):
    """Hydrogen is monovalent; JmolNN can score a short H-bond as a bond.

    526's ketoenamine H sits 1.09 A from its carbon and 1.28 A from a keto
    oxygen. Left in, such an edge can tie a guest to the framework and hide it.
    """
    removed = 0
    for i, site in enumerate(struct):
        if site.specie.symbol != "H":
            continue
        nbrs = list(sg.graph.to_undirected().neighbors(i))
        if len(nbrs) <= 1:
            continue
        # Nearest HEAVY neighbour, not nearest of any kind: a disordered H
        # cluster can put another hydrogen closer than the real parent atom.
        heavy = [j for j in nbrs if struct[j].specie.symbol != "H"]
        keep = min(heavy or nbrs, key=lambda j: struct.get_distance(i, j))
        for j in nbrs:
            if j == keep:
                continue
            try:
                sg.break_edge(i, j, allow_reverse=True)
                removed += 1
            except Exception:
                pass
    return removed


def periodic_components(sg, n):
    """Split the bond graph into components, flagging which extend periodically.

    BFS carries each atom's lattice image. Reaching an already-visited atom at a
    DIFFERENT image means the walk left the cell and returned to the same atom
    one cell over, so the component is infinite - framework. Otherwise finite.
    """
    g = sg.graph
    adj = collections.defaultdict(list)
    for u, v, d in g.edges(data=True):
        img = np.array(d.get("to_jimage", (0, 0, 0)), dtype=int)
        adj[u].append((v, img))
        adj[v].append((u, -img))

    seen = set()
    comps = []
    for start in range(n):
        if start in seen:
            continue
        image = {start: np.zeros(3, dtype=int)}
        seen.add(start)
        stack = [start]
        members = [start]
        infinite = False
        while stack:
            u = stack.pop()
            for v, img in adj[u]:
                nimg = image[u] + img
                if v not in image:
                    image[v] = nimg
                    seen.add(v)
                    members.append(v)
                    stack.append(v)
                elif not np.array_equal(image[v], nimg):
                    infinite = True
        comps.append((sorted(members), infinite))
    return comps


def clean_one(path, out_dir):
    from pymatgen.core import Structure
    from pymatgen.analysis.graphs import StructureGraph
    from pymatgen.analysis.local_env import JmolNN

    stem = os.path.splitext(os.path.basename(path))[0]
    row = {"cif": stem, "status": "", "atoms_before": 0, "atoms_after": 0,
           "guests_removed": 0, "guest_formulas": "", "n_framework_comps": 0}
    try:
        s = Structure.from_file(path)
    except Exception as e:
        row["status"] = f"read_error:{type(e).__name__}"
        return row, None
    row["atoms_before"] = len(s)
    try:
        sg = StructureGraph.from_local_env_strategy(s, JmolNN())
    except Exception as e:
        row["status"] = f"graph_error:{type(e).__name__}"
        return row, None
    _prune_overbonded_h(sg, s)

    comps = periodic_components(sg, len(s))
    framework = [c for c, inf in comps if inf]
    finite = [c for c, inf in comps if not inf]
    row["n_framework_comps"] = len(framework)

    if not framework:
        # Nothing periodic was found. Could be a 0D cage, or bond perception
        # failed. Either way, removing "all finite components" would delete the
        # whole structure, so leave it alone and say so.
        row["status"] = "no_framework_detected__left_unchanged"
        return row, s

    if not finite:
        row["status"] = "clean_already"
        row["atoms_after"] = len(s)
        return row, s

    guest_idx = sorted(i for c in finite for i in c)
    forms = collections.Counter()
    for c in finite:
        f = collections.Counter(s[i].specie.symbol for i in c)
        forms["".join(f"{e}{n}" for e, n in sorted(f.items()))] += 1
    keep = [i for i in range(len(s)) if i not in set(guest_idx)]
    cleaned = Structure.from_sites([s[i] for i in keep])
    row["status"] = "guests_removed"
    row["atoms_after"] = len(cleaned)
    row["guests_removed"] = len(guest_idx)
    row["guest_formulas"] = ";".join(f"{n}x{f}" for f, n in forms.most_common())
    return row, cleaned


def csd_components(path):
    """Independent second opinion from the CSD toolkit."""
    import ccdc.io
    c = ccdc.io.CrystalReader(path)[0]
    m = c.molecule
    poly = [x for x in m.components if x.is_polymeric]
    finite = [x for x in m.components if not x.is_polymeric]
    forms = collections.Counter(x.formula for x in finite)
    return len(poly), sum(len(x.atoms) for x in finite), forms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("in_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", default="guest_removal_report.csv")
    ap.add_argument("--csd-check", type=int, default=0,
                    help="cross-check the first N structures against the CSD toolkit")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.in_dir, "*.cif")))
    rows = []
    tally = collections.Counter()
    for k, fp in enumerate(files):
        row, cleaned = clean_one(fp, args.out)
        if cleaned is not None:
            try:
                cleaned.to(filename=os.path.join(args.out, os.path.basename(fp)), fmt="cif")
            except Exception as e:
                row["status"] += f"|write_error:{type(e).__name__}"
        rows.append(row)
        tally[row["status"].split("|")[0]] += 1
        if (k + 1) % 100 == 0:
            print(f"  ... {k+1}/{len(files)}", flush=True)

    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nprocessed {len(files)} structures -> {args.out}")
    for k, v in tally.most_common():
        print(f"  {k:38s} {v}")
    removed = [r for r in rows if r["guests_removed"]]
    print(f"  total guest atoms removed: {sum(r['guests_removed'] for r in removed)}")
    print(f"  report: {args.report}")

    if args.csd_check:
        print(f"\nCSD cross-check on first {args.csd_check}:")
        agree = disagree = 0
        for r in rows[:args.csd_check]:
            fp = os.path.join(args.in_dir, r["cif"] + ".cif")
            try:
                npoly, nguest_atoms, forms = csd_components(fp)
            except Exception as e:
                print(f"   {r['cif']}: CSD error {type(e).__name__}")
                continue
            ours = r["guests_removed"]
            ok = (ours > 0) == (nguest_atoms > 0)
            agree += ok
            disagree += (not ok)
            flag = "" if ok else "   <-- DISAGREE"
            print(f"   {r['cif']:>6}: ours={ours:4d} guest atoms | CSD={nguest_atoms:4d} "
                  f"({len(forms)} distinct){flag}")
        print(f"  agreement: {agree}/{agree+disagree}")


if __name__ == "__main__":
    main()
