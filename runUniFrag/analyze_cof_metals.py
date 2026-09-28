#!/usr/bin/env python
"""Analyze metallo-COF distribution across the CoRE-COF CIF database.

Scans every .cif file in a folder, determines which structures contain a
metal element (using pymatgen's is_metal on each species), and reports:
  - metallo-COF count vs organic-only count
  - per-metal structure counts (which metal, how many)
  - co-occurrence: structures with more than one metal
"""
import argparse
import warnings
from collections import Counter
from pathlib import Path

warnings.filterwarnings("ignore")

# Elements to treat as "metal" for metallo-COF classification (transition
# metals, alkali/alkaline-earth, post-transition metals). Metalloids like Si,
# B, Sb are excluded since they are common COF backbone/linkage elements, not
# framework metal centers.
METALLIC_ELEMENTS = {
    "Li", "Na", "K", "Rb", "Cs",
    "Mg", "Ca", "Sr", "Ba",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src_dir", default="CoRECOF/CoRE-COFs_DT1242-v7.0")
    ap.add_argument("--out_md", default="CoRECOF/cof_metal_analysis.md")
    args = ap.parse_args()

    from pymatgen.core import Structure

    src_dir = Path(args.src_dir)
    cif_files = sorted(src_dir.glob("*.cif"))

    metal_struct_count = Counter()  # metal -> number of structures containing it
    struct_metals = {}  # filename -> tuple of metals present
    n_metals_hist = Counter()  # number of distinct metals per structure -> count
    failed = []

    for f in cif_files:
        try:
            s = Structure.from_file(str(f))
            els = {str(sp.symbol) for sp in s.composition.elements}
        except Exception as e:
            failed.append((f.name, str(e)))
            continue
        metals = tuple(sorted(els & METALLIC_ELEMENTS))
        struct_metals[f.name] = metals
        n_metals_hist[len(metals)] += 1
        for m in metals:
            metal_struct_count[m] += 1

    total = len(struct_metals)
    n_metallo = sum(1 for m in struct_metals.values() if m)
    n_organic = total - n_metallo

    print(f"Parsed {total}/{len(cif_files)} structures ({len(failed)} failed)")
    print(f"Metallo-COFs: {n_metallo} ({100*n_metallo/total:.2f}%)")
    print(f"Organic-only COFs: {n_organic} ({100*n_organic/total:.2f}%)")

    lines = []
    lines.append("# CoRE-COF Metallo-COF Distribution")
    lines.append("")
    lines.append(f"Source: `{args.src_dir}` ({len(cif_files)} CIF files, {total} parsed, {len(failed)} failed)")
    lines.append("")
    lines.append("Metal set: alkali/alkaline-earth, transition metals, and post-transition "
                  "metals (Al, Ga, In, Sn, Tl, Pb, Bi). Metalloids/nonmetals commonly used as "
                  "COF linkages (B, Si, P, Sb) are NOT counted as metals here.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- **Metallo-COFs** (contain >=1 metal): **{n_metallo}** ({100*n_metallo/total:.2f}%)")
    lines.append(f"- **Organic-only COFs** (no metal): **{n_organic}** ({100*n_organic/total:.2f}%)")
    lines.append("")
    lines.append("## Metal Distribution (structures containing each metal / total parsed)")
    lines.append("")
    lines.append("| Metal | # Structures | % of DB | % of Metallo-COFs |")
    lines.append("|---|---|---|---|")
    for el, cnt in metal_struct_count.most_common():
        pct_db = 100 * cnt / total
        pct_metallo = 100 * cnt / n_metallo if n_metallo else 0.0
        lines.append(f"| {el} | {cnt} | {pct_db:.2f}% | {pct_metallo:.2f}% |")

    lines.append("")
    lines.append("## Distribution of Distinct Metal Count per Structure")
    lines.append("")
    lines.append("| # Distinct Metals | # Structures |")
    lines.append("|---|---|")
    for n, cnt in sorted(n_metals_hist.items()):
        lines.append(f"| {n} | {cnt} |")

    multi_metal = [f for f, m in struct_metals.items() if len(m) > 1]
    if multi_metal:
        lines.append("")
        lines.append("## Multi-Metal Structures (bimetallic or more)")
        lines.append("")
        lines.append("| File | Metals |")
        lines.append("|---|---|")
        for f in sorted(multi_metal):
            lines.append(f"| {f} | {', '.join(struct_metals[f])} |")

    if failed:
        lines.append("")
        lines.append("## Failed to Parse")
        lines.append("")
        lines.append("| File | Error |")
        lines.append("|---|---|")
        for f, e in failed:
            lines.append(f"| {f} | {e} |")

    out_path = Path(args.out_md)
    out_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
