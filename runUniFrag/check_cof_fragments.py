#!/usr/bin/env python3
"""Automated COF fragment QA against cof-fragment-checklist.md.

Usage:  python3 check_cof_fragments.py <run_dir> [--n-input N]

<run_dir> holds fragments_collection.extxyz, fragmentation_summary.csv and
(optionally) log.out.  Exits 0 when every MUST passes, 1 otherwise.
"""
import sys, os, csv, re, json, math, collections
import numpy as np

# Bond perception and bond-order tables are taken from UniFrag itself so this
# report and the fragmenter can never disagree about what is bonded.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentation_oop import COFFragmenter
_F = COFFragmenter()

VALENCE = {"H":1,"B":3,"C":4,"N":3,"O":2,"F":1,"Si":4,"P":3,"S":2,"Cl":1,"Br":1,"I":1}
Z = {"H":1,"B":5,"C":6,"N":7,"O":8,"F":9,"Si":14,"P":15,"S":16,"Cl":17,"Br":35,"I":53}

def read_extxyz(path):
    frames=[]
    with open(path) as f:
        while True:
            line=f.readline()
            if not line: break
            line=line.strip()
            if not line: continue
            try: n=int(line)
            except ValueError: continue
            comment=f.readline()
            sp=[]; co=[]
            for _ in range(n):
                parts=f.readline().split()
                sp.append(parts[0]); co.append([float(x) for x in parts[1:4]])
            m=re.search(r"label=(\S+)", comment)
            frames.append({"label": m.group(1) if m else "?",
                           "species": sp, "coords": np.array(co)})
    return frames

def bonds_of(sp, co):
    """Neighbour list using UniFrag's own is_valid_bond."""
    n=len(sp)
    adj={i:[] for i in range(n)}
    if n < 2: return adj
    d=np.linalg.norm(co[:,None,:]-co[None,:,:], axis=-1)
    for i in range(n):
        for j in range(i+1,n):
            if _F.is_valid_bond(sp[i], sp[j], float(d[i,j])):
                adj[i].append(j); adj[j].append(i)
    return adj

def analyse(fr):
    """Only checks that do not need aromatic perception.

    A naive valence count is meaningless on aromatics (every aromatic carbon
    looks under-valent, and a triazine C-N at 1.34 A looks like a double bond),
    so instead of a global valence audit this reports two unambiguous defects:

      over_coord  - an atom with MORE sigma bonds than its valence permits
                    (a carbon with 5 neighbours, a hydrogen with 2). No bond
                    order needed; always a real bug.
      bad_terminal- a heavy atom with exactly ONE heavy neighbour whose H count
                    does not match the valence deficit implied by that bond's
                    order. This is precisely the cut-site case UniFrag's own
                    _terminal_bond_order is defined for, so the test is sound.
    """
    sp, co = fr["species"], fr["coords"]
    n=len(sp); out={"n":n}
    adj = bonds_of(sp, co)

    # Same length scale the fragmenter uses, measured here on the fragment's
    # own conjugated C-C bonds. Some CoRE-COF entries are idealised models
    # whose bonds are uniformly short; judging their caps against a fixed
    # table would disagree with the code that placed them.
    _cc=[]
    for i in range(n):
        if sp[i]!="C": continue
        for j in adj[i]:
            if j>i and sp[j]=="C":
                d=float(np.linalg.norm(co[i]-co[j]))
                if d<=1.46: _cc.append(d)
    scale=1.0
    if len(_cc)>=4:   # a single ring is enough here; the fragmenter sees the whole parent
        r=float(np.median(_cc))/1.39
        if abs(r-1.0)>=0.03: scale=float(min(1.08, max(0.92, r)))
    out["length_scale"]=scale

    over=[]; bad_term=[]
    for i in range(n):
        tgt=VALENCE.get(sp[i])
        if tgt is None: continue
        nbrs=adj[i]
        if len(nbrs) > tgt:
            over.append((i, sp[i], len(nbrs)))
        heavy=[j for j in nbrs if sp[j]!="H"]
        nh=len(nbrs)-len(heavy)
        if sp[i]!="H" and len(heavy)==1:
            j=heavy[0]
            d=float(np.linalg.norm(co[i]-co[j]))
            order=_F._terminal_bond_order(sp[i], sp[j], d, scale)
            deficit=tgt-order
            if nh != deficit:
                bad_term.append((i, sp[i], f"H={nh} need={deficit}"))
    out["over_coord"]=over
    out["bad_terminal"]=bad_term

    zsum=sum(Z.get(s,0) for s in sp)
    out["zsum"]=zsum; out["odd_electron"]=(zsum % 2 == 1)

    clash=None
    if n >= 2:
        d=np.linalg.norm(co[:,None,:]-co[None,:,:], axis=-1)
        excl=set()
        for i in range(n):
            for j in adj[i]:
                excl.add((i,j)); excl.add((j,i))
                for k in adj[j]:
                    if k!=i: excl.add((i,k)); excl.add((k,i))
        best=1e9
        for i in range(n):
            for j in range(i+1,n):
                if (i,j) in excl: continue
                if d[i,j] < best: best=float(d[i,j])
        clash = None if best>=1e9 else best
    out["min_nonbonded"]=clash

    seen=set(); pieces=0
    for i in range(n):
        if i in seen: continue
        pieces+=1; stack=[i]; seen.add(i)
        while stack:
            u=stack.pop()
            for v in adj[u]:
                if v not in seen: seen.add(v); stack.append(v)
    out["pieces"]=pieces
    return out

def main():
    run=sys.argv[1]
    n_input=None
    if "--n-input" in sys.argv:
        n_input=int(sys.argv[sys.argv.index("--n-input")+1])
    fails=[]; warns=[]; report={}

    # ---------- section 0: run level ----------
    csv_path=os.path.join(run,"fragmentation_summary.csv")
    if os.path.exists(csv_path):
        rows=list(csv.DictReader(open(csv_path)))
        tally=collections.Counter()
        for r in rows:
            if r["normal_atoms"].strip()=="ERROR": tally["ERROR"]+=1
            elif r["norm_duplicate"].strip()=="yes": tally["duplicate"]+=1
            else: tally["produced"]+=1
        report["csv"]=dict(tally); report["csv_rows"]=len(rows)
        if tally["ERROR"]: fails.append(f"[S0] {tally['ERROR']} structures ERRORed")
        if n_input and len(rows)!=n_input:
            fails.append(f"[S0] CSV has {len(rows)} rows, expected {n_input}")
    else:
        warns.append("[S0] no fragmentation_summary.csv")

    log=os.path.join(run,"log.out")
    if os.path.exists(log):
        t=open(log, errors="replace").read()
        fnf=len(re.findall(r"No such file or directory: '(?:cof_nodes_lib|cof_linkers_lib)", t))
        to=len(re.findall(r"TimeoutError", t))
        nolink=len(re.findall(r"no known linkage chemistry", t))
        qm=len(re.findall(r"QM WARNING", t))
        report["log"]={"lib_file_not_found":fnf,"timeouts":to,
                       "no_linkage":nolink,"qm_warning":qm}
        if fnf: fails.append(f"[S0] {fnf} helper-library FileNotFoundError (parallel race)")
        if to: warns.append(f"[S0] {to} structures timed out")
        if nolink: warns.append(f"[S1] {nolink} structures with no recognised linkage")
        if qm: warns.append(f"[S5] {qm} QM WARNING (unfixable odd electron)")

    # ---------- per-fragment ----------
    xyz=os.path.join(run,"fragments_collection.extxyz")
    if not os.path.exists(xyz):
        fails.append("[S0] no fragments_collection.extxyz"); _emit(report,fails,warns); return
    frames=read_extxyz(xyz)
    report["frames"]=len(frames)

    degenerate=[]; odd=[]; overv=[]; underv=[]; clashing=[]; multi=[]
    stems=set()
    for fr in frames:
        m=re.match(r"(\d+)", fr["label"])
        if m: stems.add(m.group(1))
        a=analyse(fr)
        lab=fr["label"]
        # S2 degenerate helper fragments
        if ("OnlyLinker" in lab or "OnlyNode" in lab) and a["n"] < 6:
            degenerate.append((lab, a["n"], "".join(sorted(fr["species"]))))
        # S5 parity
        if a["odd_electron"]: odd.append((lab, a["zsum"]))
        # S3 valence
        if a["over_coord"]: overv.append((lab, a["over_coord"][:3]))
        if a["bad_terminal"]: underv.append((lab, a["bad_terminal"][:3]))
        # S4/S6 clash
        if a["min_nonbonded"] is not None and a["min_nonbonded"] < 2.0:
            clashing.append((lab, round(a["min_nonbonded"],2)))
        if a["pieces"] > 2: multi.append((lab, a["pieces"]))

    report["unique_stems"]=len(stems)
    report["degenerate"]=degenerate
    report["odd_electron"]=len(odd)
    report["over_coordinated"]=len(overv)
    report["bad_terminal_caps"]=len(underv)
    report["clashing"]=len(clashing)
    report["multi_piece"]=len(multi)

    if degenerate: fails.append(f"[S2] {len(degenerate)} degenerate helper fragments: {degenerate[:6]}")
    if odd: fails.append(f"[S5] {len(odd)} odd-electron fragments (multiplicity != 1): {odd[:6]}")
    if overv: fails.append(f"[S3] {len(overv)} fragments with over-coordinated atoms: {overv[:4]}")
    if clashing: warns.append(f"[S4/S6] {len(clashing)} fragments with non-bonded contact < 2.0 A: {clashing[:6]}")
    if underv: warns.append(f"[S3] {len(underv)} fragments with mis-capped terminal sites: {underv[:4]}")
    if multi: warns.append(f"[S6] {len(multi)} fragments in >2 disconnected pieces: {multi[:6]}")

    _emit(report, fails, warns)

def _emit(report, fails, warns):
    print("="*68); print("COF FRAGMENT CHECKLIST REPORT"); print("="*68)
    print(json.dumps(report, indent=2, default=str))
    print("-"*68)
    if fails:
        print(f"MUST FAILURES ({len(fails)}):")
        for f in fails: print("  FAIL " + f)
    else:
        print("MUST: all passed")
    if warns:
        print(f"WARNINGS ({len(warns)}):")
        for w in warns: print("  WARN " + w)
    print("="*68)
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
