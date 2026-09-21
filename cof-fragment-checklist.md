# COF Fragment QA Checklist

Standard acceptance checklist for COF cluster models produced by `COFFragmenter`.
Every threshold below is calibrated against real failures found on the CoRE-COF
884-structure set; the case IDs in parentheses are the CIF stems that motivated
the rule.

**How to use it.** Work top to bottom. Sections 0–2 are cheap and catch the
failures that invalidate everything downstream, so never skip them to get to the
chemistry. A fragment is QM-ready only when every **MUST** passes.

---

## 0. Run-level sanity — do this before looking at any fragment

The single most important lesson from this dataset: **a zero exit code does not
mean the run succeeded.**

- [ ] **MUST — count `ERROR` rows in the summary CSV, not the exit status.**
      A 112-core run exited `0:0` in 5m59s while silently losing 329 of 884
      structures (37%). Per-structure exceptions are caught and logged; the
      driver still exits clean.
      ```bash
      python3 -c "
      import csv,collections
      c=collections.Counter()
      for r in csv.DictReader(open('cifs/fragmentation_summary.csv')):
          if r['normal_atoms'].strip()=='ERROR': c['ERROR']+=1
          elif r['norm_duplicate'].strip()=='yes': c['duplicate']+=1
          else: c['produced']+=1
      print(dict(c), 'total', sum(c.values()))"
      ```
      Expected: `ERROR == 0`, and `produced + duplicate + ERROR == n_input_cifs`.

- [ ] **MUST — the CSV row count equals the number of input CIFs.** A short CSV
      means structures were dropped before they were even attempted.

- [ ] **MUST — no `FileNotFoundError` on `cof_nodes_lib/` or `cof_linkers_lib/`.**
      With `--nproc` all workers share those directories. Any unguarded
      filesystem op in the prune/export path is a race.
      ```bash
      grep -oE "Error: .*" log.out | sed -E 's/[0-9]+_[0-9]+\.xyz/<F>.xyz/' | sort | uniq -c | sort -rn
      ```

- [ ] **MUST — twin jobs write to separate directories.** UniFrag writes
      `fragments_collection.extxyz` and `fragmentation_summary.csv` *into the
      input directory*. Two jobs pointed at the same `cifs/` (including via a
      symlink) will interleave frames into one file and produce a
      plausible-looking but corrupt collection. Copy the input directory; do not
      symlink it.

- [ ] **SHOULD — check timeouts.** `grep -c TimeoutError log.out`. The default
      `--timeout 300` is too short for the largest frameworks (913 severs 192
      bonds); use `--timeout 1800` for full-set runs. Timed-out CIFs are moved to
      `cifs/timed_out_structures/` and must be restored before a rerun.

- [ ] **SHOULD — record counts of `QM WARNING` and `no known linkage chemistry`**
      as the run's quality baseline, so regressions are visible next time.

---

## 0b. Bond perception

- [ ] **MUST — no hydrogen has more than one heavy neighbour.** `JmolNN` scores
      by distance alone, so a short non-covalent contact becomes a bond. In 526
      each beta-ketoenamine H sits 1.09 A from its own carbon and **1.28 A from
      a neighbouring keto oxygen** - the resonance-assisted H-bond drawn as a
      dashed line in every picture of a TpPa COF. Those six phantom edges
      bridged node and linker, so severing all six genuine linkages still left
      ONE 114-atom component: no blocks at all, silent fallback to Path A, and
      a "node" that was neither building block. Dropping them recovers the
      textbook decomposition, 3 x C12H8N4O4 linkers and 2 x C9H3O3 nodes.
      `coffragmentor` now keeps only each hydrogen's nearest neighbour.
      **35 of the 1242 CoRE-COF structures are affected** (worst: 284 with 48
      such hydrogens, then 133, 134, 832, 286, 55, 843, 525).

- [ ] **A structure that logs cut bonds can still produce no blocks.** 526 logged
      `6 bonds severed (imine)` and then fell to Path A, because the cuts did
      not disconnect anything. Treat `linkages recognised` followed by a
      `Path A` line as a failure signal, not a success.

## 0c. Guest / solvent removal (pre-processing)

Run `runUniFrag/remove_guest_molecules.py` over the raw CIFs BEFORE fragmenting
and point UniFrag at the cleaned collection. Originals are never modified.

- [ ] **MUST — strip finite (non-periodic) components.** A framework is
      periodic: following its bonds you leave the unit cell and return to the
      same atom in a different lattice image. A guest closes on itself with zero
      net translation. That, not size or composition, is the discriminator, and
      it is what the script tests (BFS carrying each atom's lattice image).
- [ ] **MUST — leave a structure alone when no periodic component is found.**
      Otherwise "remove every finite component" deletes the whole thing. The
      script flags these as `no_framework_detected__left_unchanged`.
- [ ] Measured on the 884 HCNO set: **11 structures carry guests, 436 atoms
      removed**. 1180 holds **28 acetone molecules (280 atoms)**, 476 holds
      C12H18N3O and C8H8O2, 267 a C20H12, 929/930 water and hydroxyl, and five
      more carry stray single hydrogens.
- [ ] Cross-checked against the CSD toolkit (`component.is_polymeric`):
      **10/11 exact atom-count agreement, 0 false negatives in 40 sampled
      "clean" structures**. The one difference is 70, whose disordered H cluster
      puts a hydrogen 1.07 A from a carbon and 0.86 A from another hydrogen;
      CSD calls it isolated, we bond it. Ambiguous input, not a systematic
      disagreement.
- [ ] **CORRECTION — the acetone in 1180 was never the problem.** An earlier
      version of this file claimed that removing it turned two confused nodes
      into a clean node/linker pair. That was a false comparison across two code
      versions: the two-node result came from the RAW cif run BEFORE the
      monovalent-H fix, the clean result from the CLEANED cif run after it.
      Fragmenting 1180 raw and cleaned with the same current code gives
      byte-identical output - FragCof(179), Min(128), OnlyLinker(68:C32H30N2O4),
      OnlyNode(53:C25H24N4) - so the improvement belongs entirely to the
      monovalent-H fix. **When attributing a change, hold the code constant.**

- [ ] **But measured over the full 884, cleaning changed nothing.** A complete
      run on the cleaned collection reproduced the uncleaned numbers, and NONE
      of the 11 cleaned structures moved on any metric - 476 and 929 emit
      byte-identical fragments either way. The reason is that UniFrag already
      discards guests downstream: the zero-cut-component skip in
      `coffragmentor` drops any block with no severed bond, which is exactly
      what a guest is. Treat this pass as a guarantee and a hygiene step for
      per-structure work, not as something that will move the aggregate
      numbers.

## 1. Linkage recognition

- [ ] **MUST — at least one bond severed**, unless the framework is genuinely
      fully fused. `0 bonds severed` emits:
      `!! UniFrag COF WARNING [stem]: no known linkage chemistry recognised`.
      Baseline on the full set: **31/884 (3.5%)** unrecognised, nearly all fully
      fused frameworks where node/linker decomposition does not apply.

- [ ] **MUST — the recognised linkage matches the chemistry you expect.**
      Full-884 distribution:
      | Linkage | Count |
      |---|---|
      | imine | 727 |
      | biaryl | 53 |
      | imine + vinylene | 27 |
      | vinylene | 25 |
      | oxazole | 8 |
      | dioxin + imine | 7 |
      | dioxin | 4 |
      | oxazole + vinylene | 2 |
      | none | 31 |

- [ ] **MUST — linkage tags are per-bond, never per-type.** Tags are
      `frozenset((u, v))`. A shared per-type tag makes `_count_linkages` merge
      independent bonds and undercount.

- [ ] **Vinylene (C=C) rule** — both carbons `heavy_degree == 2`, `h_count == 1`,
      bond ≤ 1.42 Å, **and not in a small ring**. The ring test is the only thing
      separating a vinylene linker from an ordinary aromatic CH=CH edge; dropping
      it shreds every aromatic ring.

- [ ] **Imine (C=N) rule** — requires `heavy_degree(C) == 2`,
      `heavy_degree(N) >= 2`, not in a small ring. The terminal-N guard is what
      stops a terminal `=NH` being treated as a linkage (760).

- [ ] **DO NOT add a bond-length guard to the imine rule.** This was proposed and
      **rejected on measurement**: across 810 cut C–N bonds, p75 = 1.348 Å and
      p95 = 1.447 Å. A 1.36 Å cutoff makes **15 of 100 structures lose every cut**
      (543 and 585 among them). The rule's length-blindness is load-bearing.

- [ ] **Secondary amine (C-N) rule - TRIED AND REVERTED, do not re-add
      without a per-structure review.** N with two heavy neighbours, both
      carbon, and one H, outside a small ring. It looked good on aggregate
      metrics over the full 884 (`no_linkage` 43->40, over-coordinated 71->67,
      odd-electron 40->39, mis-capped 292->287, clashing 835->828) and it did
      fix 1015, moving it off Path A onto Path J with a clean dimer.
      It was still wrong. The real motif is `Ar-NH-CH2-R`, not `Ar-NH-Ar`, and
      when the aliphatic side is a dead-end tail the cut detaches a pendant
      substituent that the linkage-count classifier then calls a linker. On
      1180 it emitted an 11-atom `HN-CH2-C(=O)-CH3` scrap AND destroyed a
      previously correct decomposition - two proper nodes (61 and 81 atoms)
      vanished. On 525 it split fragments into five disconnected pieces.
      Two guards were tried and both failed: requiring both flanking carbons
      to be aromatic matched nothing (every bridging N has exactly one ring and
      one non-ring carbon), and rejecting components with a single attachment
      point did not catch the scrap either.
      **Lesson: aggregate valence/clash metrics cannot see "this is not
      chemically a strut."** Every linkage-chemistry change needs a
      per-structure eyeball on the blocks it creates, not just a metric diff.

- [ ] **MUST — a fragment is one connected molecule, or two for a dimer.** A
      dimer's halves are comparable in size, so anything below ~40% of the
      largest piece is detached junk, not a layer. Baseline: 5 of 1737 frames
      (1103 carries six lone H, 1179 three detached 14-atom pieces).

- [ ] **Biaryl (C–C) fallback** — runs *only* when no other edge was found. Fuse
      rings ≤7 into ring systems; a biaryl bond joins two different systems and
      lies outside any ring; sever only bonds touching a system of degree ≥3.
      Vinylene + biaryl together cover **107/884 (12%)** of the set.

---

## 2. Node / linker decomposition

- [ ] **MUST — both a node and a linker are exported** when the framework is not
      fully fused. "Linker but no node", or a node with no linker, means the
      classification collapsed (411, 612).

- [ ] **MUST — no degenerate helper fragments.** These are always wrong and were
      observed in production output:
      ```
      931FragCofOnlyLinker_0  -> [N]            (1 atom)
      929FragCofOnlyLinker_1  -> [H, O]         (2 atoms)
      645FragCofOnlyLinker    -> [N, N, N, H]
      1226FragCofOnlyLinker_1 -> [C, N, N, H, H]
      ```
      Flag any helper fragment with fewer than ~6 atoms for inspection.

- [ ] **MUST — no redundant carbon carried onto the node side of a severed
      linkage** (411). Carry-over is allowed only for `_CARRYOVER_ELEMENTS` or an
      explicit `vinylene_partner` pair.

- [ ] **MUST — no orphaned bridging unit.** A linkage rule matches a bond
      pattern, not a whole linkage, so it can sever EVERY bond around a small
      bridge and leave it floating: `Ar-NH-Ar` loses both C-N bonds (931), an
      `-N=N-N-` triazene is cut on all sides (645), a C bridging two N likewise
      (1226). The orphan guard in `coffragmentor` restores one cut per component
      with fewer than 4 heavy atoms, so the unit stays attached as a proper
      terminus. Do not raise that threshold past 4 without re-measuring: a
      boroxine `B3O3` node is carbon-free with 6 heavy atoms and must survive.

- [ ] **MUST — guest and solvent molecules are not exported as linkers.** A
      component with ZERO severed bonds was already disconnected in the parent
      (929 carries six water/hydroxyl molecules). It reaches the classifier with
      no linkages, falls below the `>= 3` node test and would otherwise be
      emitted as a linker.

- [ ] **k+k symmetric frameworks** — when classification yields no linkers but
      more than one node, the symmetric tie-break demotes all but the smallest
      `(smiles, n_atoms)` group to linkers. Without it, symmetric COFs return
      all-nodes and no linker.

- [ ] **MUST — linker images are chosen per attachment site, not per index.**
      Per-index selection drops symmetry-equivalent placements and yields
      linkers with missing arms (411, 704).

---

## 3. Capping chemistry — valence and bond order

- [ ] **MUST — no over-valent atom.** Every severed site is capped according to
      its *perceived bond order*, not blindly with one H.
      Double-bond ceilings (Å): C–N 1.36, C–O 1.30, C–C 1.35, N–N 1.30, N–O 1.30.
      Triple-bond ceilings: C–N 1.20, C–C 1.24.

- [ ] **MUST — nitriles and alkynes are never protonated.** The saturation guard
      (`_skip_saturated_h_cap`) must run **before** the per-element branches in
      `_cap_open_oxygens`, and must apply to C, N *and* O. Scoping it to carbon
      only drove over-valent sites from 46 to 89.

- [ ] **MUST — an aldimine terminates as `Ar–CH=NH`, not `Ar–CH₃`.** Note: an
      `Ar–CH₃` in the output is *usually not* a capping bug. Measured on this set,
      153 aldimines produced **zero** imine-derived methyls; the methyls found
      were native alkyl side chains present in the parent CIF (e.g. 502's `C[27]`
      carries 3 H in the raw block, an `–O–CH₂–CH₃` ethoxy). Check the parent
      before blaming the capper.

- [ ] **MUST — no superimposed atoms.** Fragment assembly can emit the same
      atom twice - a carried-over heteroatom the 0.1 A guard in `coffragmentor`
      misses, a linker image on an already-filled site, a dimer layer
      overlapping its source. Duplicates land 0.15-0.75 A apart and surface as
      hydrogens with two bonds and carbons with six; because each copy carries
      its own electrons they ALSO flip the parity that the multiplicity repair
      is about to fix. `_dedupe_superimposed_atoms` (COF-only hook, runs inside
      `fix_odd_electron_multiplicity` before the repair) drops anything within
      0.85 A - below the shortest real bond here, O-H at ~0.96 A. Baseline
      before the fix: 42 of 1695 frames affected.

- [ ] **Separate inherited damage from capping bugs before chasing one.** Check
      the PARENT for the same defect. Measured on the failing set: 284's parent
      has a 0.78 A minimum distance and **80 over-coordinated atoms of its own**;
      1091's is 0.75 A; 1180/526/671 carry genuine quaternary `N(CCCC)`. Those
      are inherited and unfixable downstream. Only structures whose parent is
      clean (1015, 1061, 1112, 1128, 652, 900, 901) indicate a real capping bug.

- [ ] **SHOULD — under-valent sites are reviewed.** They are less damaging than
      over-valent ones but still distort the electronic structure. Baseline on
      300 structures: 46 under-valent, 11 over-valent, 98.8% valence correct.

---

## 4. Capping geometry — planarity

- [ ] **MUST — capping H on a conjugated site lies in the ring plane.** UFF's
      pyramidal default puts an NH₂ hydrogen ~0.85 Å out of plane, which is both
      chemically wrong (aryl amines are near-planar from conjugation) and the
      direct cause of interlayer clashes in stacked dimers.

- [ ] Planarization applies when `sp == "N"` or `anchor_order >= 2`, fits a plane
      by SVD over the conjugated system (BFS depth ≤3 from the anchor, N
      included), and places H at ±120° sp² slots.
      Guards: **≥5 atoms in the ring system** (`len(ring_atoms) < 5` bails), **flatness ≤ 0.35 Å**.

- [ ] **Do not "fix" this by capping N with a single H.** It was considered and
      rejected: measured NH1 sites were already 0.85 Å out of plane — no better
      than NH2 — and a single H leaves a radical, breaking the multiplicity = 1
      requirement.

---

## 5. Electron parity and multiplicity

- [ ] **MUST — every fragment is closed-shell, multiplicity = 1.** This is a hard
      project requirement, not a preference.

- [ ] **Parity rule** — for a neutral closed-shell CHNO fragment, `h + n` must be
      even (h = hydrogen count, n = nitrogen count).

- [ ] **MUST — repair parity by ADDING an H, not removing one**, wherever a site
      can accept it (`_can_accept_extra_h`: `_local_valence_used(idx) < target`).
      Removal is gated off whenever an add-candidate exists. Removing H opens a
      valence that was chemically correct.

- [ ] **MUST — parity repair runs AFTER dimerization**, in both the Path J branch
      and the A/B/C/D branch, with `capped_h_indices` extended across the
      duplicated layer. Repairing the monomer first leaves the dimer odd (167).

- [ ] **MUST — the main collection contains no open-shell fragment.** Anything
      that still has an odd electron count after the repair is written to
      `fragments_quarantine.extxyz` instead of `fragments_collection.extxyz`,
      tagged in its own comment line (`quarantine=odd_electron zsum=N`) and
      announced per fragment in the log plus a closing `QUARANTINE:` summary.
      So `odd_electron` in a checklist report on the main collection should read
      0; the count to watch is the size of the quarantine file. Current residue
      is 5 of 1737: `1005FragCofOnlyLinker`, `1104FragCof`, `1104FragCofMin`,
      `615FragCofOnlyLinker`, `744FragCofMin`.

- [ ] Known residual: `142FragCofMin` is a true monomer with no benign parity
      site (1 of 78). A `[QM WARNING] could not automatically fix odd electron
      count` needs manual inspection — do not ship it silently.

---

## 6. Stacked dimers

- [ ] **MUST — the layer is genuinely stackable before a dimer is built.**
      Unified `_is_stackable_layer` criteria:
      - flatness ≤ `max(0.6, 0.30 × |T|)` Å (relative to the translation, not absolute)
      - alignment ≥ 0.8
      - **measured** closest contact ≥ 2.4 Å

- [ ] **MUST — dimer geometry within bounds:** interlayer separation 2.5–5.0 Å,
      lateral offset ≤ 2.0 Å, closest contact ≥ 2.4 Å.

- [ ] **MUST — the planarity guard covers every path, not just A/B/C/D.** Scoping
      it narrowly passed on 20 structures purely by luck; at 100 structures
      **48 of 63 clashing frames were non-flat Path J dimers**.

- [ ] **Partner ops** — lattice translation, crystal symmetry op, or local screw
      rotation, as `(R, t, local)`. Sort key `(|delta|, rotation_angle, -contact)`;
      the rotation-angle term prefers the primitive screw (585 picked 180° over
      the correct 60° without it).

- [ ] **MUST — `lattice_reduce` slides laterally only.** Minimizing *total*
      displacement collapses pure lattice translations to zero and silently
      removes **every** dimer from the set.

- [ ] **MUST — a local screw op is re-centred per block.** Applying it verbatim to
      node and linker blocks (which have different centroids) puts the halves
      33–38 Å apart.

- [ ] **MUST — Path B cannot dimerize twice.** Guarded by `dimer_already_built`.

- [ ] Expect roughly **1 dimer per 3 fragments** (331 dimers / 906 frames on the
      300-structure validation).

---

## 7. Duplicate detection

- [ ] **MUST — COF duplicates use the topology-aware key**: heavy-atom formula +
      element-labelled Weisfeiler–Lehman hash of the heavy-atom bond graph.
      Formula alone merges constitutional isomers.

- [ ] **MUST NOT change the MOF key.** MOF duplicate detection remains
      formula-only by deliberate decision. `_identity_key_fn` is selected by
      `args.kind`; `_flush_cof_result` and `_flush_mof_result` keep separate keys.

---

## 8. MOF isolation — regression guard

Standing project constraint: **COF work must never alter MOF results.**

- [ ] **MUST — every COF-specific behaviour enters through a hook**: a no-op on
      `BaseFragmenter`, overridden on `COFFragmenter`. Never inline a COF-only
      helper into a `BaseFragmenter` method — `_local_valence_used` inlined into
      `_cap_open_oxygens` would have raised `AttributeError` on every MOF run.

- [ ] **MUST — verify identity, not intent**, after touching shared code:
      ```python
      MOFFragmenter._skip_saturated_h_cap is BaseFragmenter._skip_saturated_h_cap  # True
      ```
      For a COF-only method, assert it is not reachable from the MOF class at all:
      ```python
      hasattr(MOFFragmenter, '_prune_duplicate_cof_helper_files')  # False
      ```

- [ ] **MUST NOT apply COF-only heuristics to MOFs** — e.g. the reduction < 20
      rule is COF-only by explicit decision.

---

## 9. Parent structure caveats

- [ ] **MUST — parent CIFs stay untouched.** Every structure remains original;
      all repair happens on the extracted fragment.

- [ ] **Know the parent's own quality before blaming the fragmenter.**
      **147 of 567 layered parents (26%)** have a genuine sub-2.0 Å *non-bonded*
      interlayer contact in the source CIF. Heavy-atom contacts are healthy
      (p50 = 3.27 Å) — it is hydrogen placement in the deposited structure.

- [ ] **When measuring clashes, exclude pairs inside the covalent cutoff.** A
      first attempt that counted covalent bonds crossing the cell boundary
      reported 333/567 "clashing" parents; the corrected figure is 147/567.

---

## Appendix — validation methodology

Anti-patterns that cost real time on this project:

0. **Clear the output directory before a re-run.** UniFrag deduplicates
   against the existing `fragments_collection.extxyz`, so re-running in place
   after a code change silently reproduces the OLD result: the log goes quiet
   and the frames come back byte-identical. Delete the extxyz, the CSV and both
   helper libraries first, or you will validate the fix you did not apply.
1. **Never justify a path-specific guard on 20 structures.** 20 is too small for
   path coverage to be meaningful; the Path J planarity gap was invisible there
   and obvious at 100.
2. **Measure before adding a threshold.** The imine bond-length guard looked
   chemically reasonable and would have destroyed 15% of the set. Percentiles
   over the actual cut-bond population settled it in minutes.
3. **Trust the failure count, not the exit code.** See §0.
4. **RDKit `DetermineBonds` is combinatorial** and hangs indefinitely on ~800-atom
   fragments. Any RDKit-level validation at set scale needs a per-fragment
   timeout *and* an atom-count cap. The 300-structure validation used the
   bond-order heuristic instead, not RDKit.
5. **Distinguish "deduplicated" from "failed."** Both shrink the output; only one
   is a problem. Always split the CSV three ways: produced / duplicate / ERROR.

---

## Reference baselines

Full 884-structure CoRE-COF runs, scored by `runUniFrag/check_cof_fragments.py`.
Each column is one complete run; `d2f2054` is current.

| metric | pre-fix | +dedup/orphan `e28ee72` | +amine `5fd407d` | +parity `371c362` | +clearance `d2f2054` |
|---|---|---|---|---|---|
| frames | 1695 | 1737 | 1737 | 1737 | 1737 |
| ERROR structures | 329 | 0 | 0 | 0 | 0 |
| degenerate helpers | 4 | 0 | 0 | 0 | 0 |
| **odd-electron** | 34 | 40 | 39 | **5** | **5** |
| QM warnings | 47 | 68 | 66 | 10 | 11 |
| over-coordinated | 96 | 71 | 67 | 69 | 70 |
| mis-capped terminals | - | 292 | 287 | 214 | 212 |
| clashing (<2.0 A) | 808 | 835 | 828 | 849 | 842 |
| detached-junk frames | - | 5 | 7 | 5 | 5 |
| no recognised linkage | 41 | 43 | 40 | 37 | 40 |

Read the columns as cumulative, not independent: each adds to the one before.
`no recognised linkage` moves by a few between otherwise identical runs because
parallel processing order decides which structures are deduplicated before the
linkage report prints - treat differences under about 5 there as noise.

The parity repair is the single biggest win (odd-electron 39 -> 5) and it paid
for itself in mis-capped sites too (287 -> 212). It cost clashes (828 -> 849);
the clearance pass recovered a third of that (-> 842) and no more, because the
remaining contacts are 1.8-2.0 A H...H pairs where no roomier site exists.

### Known-unfixable residue

These do not respond to fragmenter changes and should be excluded before judging
a run, not chased:

- **Broken parent CIFs.** 284's own structure has a 0.78 A minimum interatomic
  distance and 80 over-coordinated atoms; 1091's minimum is 0.75 A. Whatever the
  fragmenter does, their fragments inherit it.
- **Genuine quaternary N.** 1180, 526, 671 carry real `N(CCCC)` centres.
- **5 odd-electron fragments** (1005, 1104 x2, 615, 744): no site can accept a
  hydrogen at any clearance. All are flagged by `QM WARNING`, never silent.
