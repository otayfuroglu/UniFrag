"""Quality checks for a metallo-COF fragment folder.

Usage: python runUniFrag/check_metal_cof_fragments.py runUniFrag/CoRECOF/cu_cofs Cu

Per structure: metal retained, metal coordination number vs parent, odd
ligand (non-metal) electron count, disconnected pieces, CH3/CH2 absent from the parent, bent
2-coordinate (under-coordinated) C, out-of-plane sp2 caps, and contacts of
added H only (0.9 A within a layer, 2.0 A between layers).
"""
import re, sys, csv, warnings, numpy as np, networkx as nx
from pathlib import Path
from collections import Counter, defaultdict
from ase.io import read
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
folder, METAL = Path(sys.argv[1]), sys.argv[2]
R={"H":0.31,"B":0.84,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"Si":1.11,"P":1.07,"S":1.05,"Cl":1.02,"Br":1.2,"I":1.39}
MET={"Li","Na","Mg","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Mo","Ru","Rh","Cd","Au","Hg"}
Z={"H":1,"B":5,"C":6,"N":7,"O":8,"F":9,"Si":14,"P":15,"S":16,"Cl":17,"Br":35,"I":53,"Co":27,"Ni":28,"Cu":29,"Zn":30,"Mo":42}
def bonded(a,b,d):
    if a=="H" and b=="H": return d<0.9
    if a in MET or b in MET:
        if a in MET and b in MET: return False
        return a not in "CH" and b not in "CH" and d<2.6 if not ({a,b}&{"C","H"}) else False
    return d<=min(2.2,max(1.1,1.25*(R.get(a,.77)+R.get(b,.77))))
def parent_info(cif):
    s=Structure.from_file(str(cif)); sym=[x.specie.symbol for x in s]
    cn=[]; me=ch2=0
    for i,site in enumerate(s):
        nb=[n for n in s.get_neighbors(site,2.7) if bonded(sym[i],n.specie.symbol,n.nn_distance)]
        if sym[i]==METAL: cn.append(len(nb))
        if sym[i]=="C":
            hs=sum(1 for n in nb if n.specie.symbol=="H")
            me+=hs==3; ch2+=hs==2
    return dict(cn=Counter(cn), methyl=me, ch2=ch2)
def frame_info(a):
    s=a.get_chemical_symbols(); p=a.get_positions(); n=len(a)
    D=np.linalg.norm(p[:,None]-p[None],axis=2)
    g=nx.Graph(); g.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i+1,n):
            if bonded(s[i],s[j],D[i,j]): g.add_edge(i,j)
    comp={}
    for k,c in enumerate(nx.connected_components(g)):
        for i in c: comp[i]=k
    caps=[int(x) for x in re.findall(r"-?\d+", str(np.asarray(a.info["capped_h"]).tolist()))] if "capped_h" in a.info else None
    # Ligand parity is what must be even (M(II) in porphyrin2-/Pc2-); an odd
    # total left over is the metal's own spin (Cu(II) d9, Co(II) d7).
    r=dict(n=n, ncomp=nx.number_connected_components(g),
           odd=sum(Z.get(x,6) for x in s if x not in MET)%2==1,
           metal_open_shell=sum(Z.get(x,6) for x in s)%2==1)
    r["metals"]=[(len([j for j in g[i]]), round(min([D[i,j] for j in g[i]] or [9]),2)) for i in range(n) if s[i]==METAL]
    r["other_metals"]=sum(1 for x in s if x in MET and x!=METAL)
    ch3=ch2=under=npl=0
    for i in range(n):
        hv=[j for j in g[i] if s[j]!="H" and s[j] not in MET]; hs=[j for j in g[i] if s[j]=="H"]
        if s[i]=="C":
            ch3+=len(hs)==3; ch2+=len(hs)==2
            nb=list(g[i])
            if len(nb)<2: under+=1
            elif len(nb)==2:
                u=p[nb[0]]-p[i]; v=p[nb[1]]-p[i]
                ang=np.degrees(np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))
                under+=ang<165   # linear = sp (alkyne/nitrile), complete
        if s[i]=="N" and not list(g[i]): under+=1
        if caps is not None and s[i] in "CN" and len(hv)==2 and len(hs)==1 and hs[0] in caps:
            nr=np.cross(p[hv[0]]-p[i],p[hv[1]]-p[i])
            if np.linalg.norm(nr)>1e-6:
                nr/=np.linalg.norm(nr); v=p[hs[0]]-p[i]
                npl+=np.degrees(np.arcsin(abs(nr@v)/np.linalg.norm(v)))>10
    r.update(ch3=ch3, ch2=ch2, under=under, nonplanar_cap=npl)
    bad=[]
    if caps is not None:
        for h in caps:
            if h>=n: continue
            for j in range(n):
                if j==h or j in g[h]: continue
                lim=0.9 if comp[j]==comp[h] else 2.0
                if D[h,j]<lim: bad.append((h,j,round(D[h,j],2)))
    r["contacts"]=len(bad); r["caps_known"]=caps is not None
    return r
rows={r["cif_file"][:-4]:r for r in csv.DictReader(open(folder/"fragmentation_summary.csv"))}
frames=defaultdict(dict)
for a in read(folder/"fragments_collection.extxyz",index=":"):
    l=a.info["label"]; stem=l.split("Frag")[0]; frames[stem][l[len(stem):]]=a
flags=defaultdict(list); tally=Counter()
cifs=sorted(folder.glob("*.cif"), key=lambda c:(not c.stem.isdigit(), int(c.stem) if c.stem.isdigit() else 0, c.stem))
timed=sorted((folder/"timed_out_structures").glob("*.cif")) if (folder/"timed_out_structures").exists() else []
print(f"# {METAL}: {len(cifs)} CIFs processed, {len(timed)} timed out {[t.stem for t in timed]}")
print(f"{'cif':>5} {'normal':>7} {'min':>5}  parentCN   frag metal CN(min dist)          defects")
for c in cifs:
    st=c.stem; P=parent_info(c); fr=frames.get(st,{}); row=rows.get(st)
    pcn=max(P["cn"],key=P["cn"].get) if P["cn"] else 0
    line=[]
    for kind in ("FragCof","FragCofMin"):
        if kind not in fr: continue
        f=frame_info(fr[kind]); tag="N" if kind=="FragCof" else "M"
        tally["frames"]+=1
        if not f["metals"]: flags[st].append(f"{tag}: METAL LOST"); tally["metal_lost"]+=1
        low=[m for m in f["metals"] if m[0]<pcn]
        if low: flags[st].append(f"{tag}: {len(low)}/{len(f['metals'])} {METAL} under-coordinated vs parent CN{pcn} {low[:3]}"); tally["metal_lowCN"]+=1
        if f["odd"]: flags[st].append(f"{tag}: odd ligand electrons"); tally["odd_ligand"]+=1
        if f["metal_open_shell"]: tally["open_shell_metal (doublet etc., info)"]+=1
        if f["ch3"]>0 and P["methyl"]==0: flags[st].append(f"{tag}: {f['ch3']} CH3 (parent has none)"); tally["ch3_artifact"]+=1
        if f["ch2"]>0 and P["ch2"]==0: flags[st].append(f"{tag}: {f['ch2']} CH2 (parent has none)"); tally["ch2_artifact"]+=1
        if f["under"]: flags[st].append(f"{tag}: {f['under']} under-coordinated C"); tally["under"]+=1
        if f["nonplanar_cap"]: flags[st].append(f"{tag}: {f['nonplanar_cap']} out-of-plane sp2 cap H"); tally["nonplanar"]+=1
        if f["contacts"]: flags[st].append(f"{tag}: {f['contacts']} added-H contacts"); tally["contacts"]+=1
        if f["ncomp"]>2: flags[st].append(f"{tag}: {f['ncomp']} disconnected pieces"); tally["pieces"]+=1
        if not f["caps_known"]: tally["no_caps_info"]+=1
        line.append(f"{tag}:{Counter(m[0] for m in f['metals'])}")
    only=[k for k in fr if "Only" in k]
    only_metal=sum(1 for k in only if METAL in fr[k].get_chemical_symbols())
    nn=row["normal_atoms"] if row else "-"; mm=row["min_atoms"] if row else "-"
    if not row or nn in ("N/A","TIMEOUT","0"): flags[st].append(f"no normal fragment ({nn})"); tally["no_normal"]+=1
    print(f"{st:>5} {nn:>7} {mm:>5}  CN{pcn} x{sum(P['cn'].values()):<3} {'  '.join(line):40} {'; '.join(flags[st]) or 'ok'}  [Only*: {only_metal}/{len(only)} with {METAL}]")
print("TALLY", dict(tally))
