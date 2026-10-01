# Project Agent Log

Chronological handoff log for agents working on UniFrag. Add newest entries at the top. Each entry should include changed files, validation, decisions, and follow-up risks.

## 2026-09-29 - UniFrag: Parent-CIF-Aware N/O Capping Fix for COFFragmenter (fixes the 21/45 wrong-capping fragments found by the prior read-only investigation)
- **Task:** implement "Option A" from a prior read-only investigation: make `COFFragmenter`'s three capping/parity-repair decision points (`_cap_open_oxygens`, `_cap_severed_double_bond_sites`, and the odd-electron parity repair) check the real parent CIF's periodic coordination number for a candidate atom, instead of guessing purely from the fragment's local connectivity and a bond-length table. Full decision writeup: `project-decisions.md`, "Decision 2026-09-29: Parent-CIF-Aware N/O Capping, Superseding the 2026-09-18 Post-Dimerization Mitigation".
- **Files touched:** `fragmentation_oop.py` only (385 insertions / 7 modified lines; `coffragmentor.py` untouched, per the task's explicit scope - bond-cutting logic was not to be touched). No new files.
  - Added to `BaseFragmenter` (all default to the previous behaviour, so MOF/macromolecule paths are unaffected): `_can_remove_cap_h` hook, `_parent_cut_deficit` hook, and wired `_can_remove_cap_h` into `fix_odd_electron_multiplicity`'s removal-candidate loop.
  - Added to `COFFragmenter` (next to `_parent_ring_heteroatoms`/`_is_parent_ring_heteroatom`): `_parent_frac_coords`, `_parent_index_for_position`, `_real_parent_heavy_degree`, `_real_parent_native_h_count`, `_parent_severed_neighbor_valence`, `_fragment_heavy_degree`.
  - Modified `_cap_open_oxygens`'s O/N branches (now loop `n_caps` times using `_parent_cut_deficit` when available), `_cap_severed_double_bond_sites`'s deficit calculation (now parent-aware with a native-H correction and a hypervalence safety clamp), and `_can_accept_extra_h`/added `_can_remove_cap_h` (COF overrides, both parent-gated).
- **Validation (all runs done locally against `runUniFrag/CoRECOF/CoRE-COF_HCNO884_clean/`, which mirrors the parent CIFs the task pointed at on the `arf` cluster - not reachable from this session, so no SSH was used):**
  1. Copied pre-fix `fragmentation_oop.py` via `git show HEAD:fragmentation_oop.py` to a scratch dir (with an unmodified `coffragmentor.py` copy) as the "old" baseline; the current worktree file is "new".
  2. Ran both, folder-mode COF, on 41/45 of the task's named affected stems (4 - `615`, and 3 more of the 45 - were not present in the local CIF mirror; not chased further) and, separately, on a disjoint random 30-stem sample for regression.
  3. Wrote a standalone auditor (`audit.py`, does not import `fragmentation_oop.py` or `coffragmentor.py`) that matches every N/O atom in every output fragment back to its real parent CIF site by Cartesian coordinate (searching lattice-translation images) and flags a mismatch when `(fragment H count - parent's own native H count) != max(0, real_parent_heavy_degree - fragment_local_heavy_degree)`.
  - **Affected-stem set** (135 frames): N/O capping mismatches **484 -> 229** (-53%); frames with >=1 mismatch **76 -> 51**. Zero chemically-impossible H counts in the final result (max 2 H on any N, 1 on any O).
  - **Regression sample** (128 frames, structures NOT flagged by the prior investigation): mismatches **290 -> 253**, frames affected **50 -> 45** - net improvement, not a regression.
  - **MOF smoke check:** `fragments_collection.extxyz` byte-identical between old and new code on `1516648.cif` + `IRMOF-1.cif` (2-structure MOF set), confirming the new hooks are true no-ops outside `COFFragmenter`.
  - **Two real bugs found and fixed mid-validation** (not part of the original diagnosis, only surfaced by re-running against parent CIFs):
    1. `_cap_severed_double_bond_sites` runs more than once on the same molecule in some COF paths; the first cut of the fix ignored existing H count in its new deficit formula, quadruple-capping a porphyrin-core nitrogen in `187FragCof` (4 H on one N - chemically impossible). Fixed by correctly subtracting only NEWLY-added H (not native H) via the new `_real_parent_native_h_count` helper.
    2. Summing each severed neighbour's bond order independently can itself double-count the same length ambiguity the fix was built to avoid (two severed neighbours both scoring "double"), observed as `severed_valence=4` on a nitrogen in `476`'s OnlyNode export. Fixed with a hard clamp (`target_valence[sp] - 1 - h_count`) at both use sites, so the mechanism can never manufacture a hypervalent atom regardless of residual bond-order ambiguity.
  - Reproduction scripts and full audit logs are in the session scratchpad only (`capvalid/` - `audit.py`, `debug*.py`, `audit_*.txt`, `*_run*.log`), not committed to the repo.
- **Follow-up risks / known gaps (documented in project-decisions.md, not fixed here - out of the task's explicit 3-function scope):**
  - `_try_cof_graph_node_linker_fragment` (the secondary "graph fallback" COF path) has its own hand-rolled O/B-only capping block that never calls the three fixed functions and is not parent-aware; structures that fall through to it (rather than the primary Path J or the generic Path A/B/C/D branch) are not covered by this fix.
  - A residual minority of the auditor's own "mismatches" trace to the auditor's simplifying assumption (bond-order-blind neighbour-count diff), not a code defect: a genuinely symmetric/macrocyclic ring or meso-nitrogen (confirmed on `1004`, `187`, both porphyrin-like) can have two real neighbours that both sit in the same 1.34-1.36 A ring-C=N/single ambiguity window, which no length-based method (old or new) can fully resolve without cleavage-rule/ring metadata - explicitly out of scope here. The new code's answer in these cases is always chemically valid (bounded, non-hypervalent) even where it may not be the one true Kekule structure.
  - Only 41 of the 45 named affected stems were validated (4 not present in the local CIF mirror used); worth re-running against the exact cached set at `arf:/arf/scratch/otayfuroglu/deepCOF_works/coreCOFs/HCNO/runOpenMLP/subset_qm250/analysis/parent_cifs/` if/when that host is reachable, though the 41-stem result should be representative.
  - Quarantine counts shifted (5->9 on the affected set, 0->1 on the regression sample): some previously "successfully but wrongly" force-capped fragments are now honestly quarantined as open-shell when no valid parent-consistent capping site exists, per the project's standing "quarantine over silently-wrong" policy (Decision 2026-09-18) - not re-litigated here, but worth a user visual check on the newly-quarantined fragments (`795`, `24FragCofMin`, `930FragCofOnlyLinker_1`, `948FragCofMin`, `1131FragCofOnlyLinker`) if the QM pipeline needs them back.

## 2026-09-17 - UniFrag: 1,4-Dioxin COF Linkage Cleavage + Explicit Linkage Tagging (fixes `1097.cif`)
- **Symptom (user-reported, visual):** `1097.cif`'s normal fragment had a node with missing parts and linkers not fully plugged in around it, and no minimized fragment was produced at all.
- **Diagnosis:** Same failure mode as `1031` — Path J could not decompose the structure, so extraction fell through to the generic `bo_node_species = {"B","O"}` heuristic, which treated an isolated ether oxygen as a whole node (`COF Path A (Single node). Node component size: 1`), giving a 53-atom `C29 H16 O8` fragment. `1097` is a **dioxin-linked COF** (`H18 C54 O12`, no boron, a=b=24.69 c=3.50): 6 six-membered `{C4,O2}` rings are 1,4-dioxin linkages, all O are ether `C-O-C`. Removing every O splits the framework into 2 triphenylene blocks (`C18H6`, touching 6 O -> 3-connected NODE) and 3 benzene blocks (`C6H2`, touching 4 O -> 2-connected LINKER); i.e. HHTP-derived chemistry. The minimize fallback could not help either, since `single_linker=True` requires the coffragmentor node+linker route.
- **Changed files:**
  - `coffragmentor.py` [MODIFY]
    - Added 1,4-dioxin cleavage in `COF.fragment()`. The ring `O1-C2-C3-O4-C5-C6` fuses two different aryl systems (`{C2,C3}` vs `{C5,C6}`, identified as connected components of the graph with all O removed); the two C-O bonds on ONE side are severed. Oxygens are kept with the **larger** aryl system, which is the polyol-derived block (hexahydroxytriphenylene in the canonical case). Guards require the two C-pairs to be genuinely different aryl systems, so an intra-ring O-C-C-O does not trigger it.
    - Added a `linkage_of` tag map (keyed by `frozenset({u,v})`) populated at detection time, and `_count_linkages(cut_bonds, linkage_of)` now merges bonds sharing an endpoint **or** carrying the same tag.
- **Why the tag map was required:** a dioxin severs `C2-O1` and `C3-O4`, which share **no endpoint at all**, so the previous endpoint-sharing rule counted the benzene linker's 4 severed bonds as 4 separate connections and classified it as a NODE (5 nodes / 0 linkers -> Path J aborts). Adjacency-based merging was rejected as ambiguous: adjacent attachment atoms mean "one linkage" for a dioxin but "two linkages" for an ortho-disubstituted aryl. Tagging at detection time is exact and costs nothing.
- **Validation:**
  - `1097.cif`: normal **53 -> 120 atoms** (`C29 H16 O8` -> `C72 H36 O12`); per layer exactly node (`C18 O6`) + 3 linkers (3 x `C6`). Minimized fragment now produced via the node+one-linker fallback: **88 atoms** `C48 H28 O12` (reduction 32). Helper libraries now correct: node 30 atoms `H6 C18 O6` (triphenylene + 6 O), linker 8 atoms `H2 C6` (benzene). All four fragments valence-complete and even-electron; all 12 dioxin oxygens reconnect into intact `C-O-C` junctions with 0 outward termini.
  - Decomposition regression: `1031` (2 x `H12 C24 N3` + 3 x `H2 C6 N2 O2`) and `211` (2 x `H15 C24 N3` + 3 x `H6 C8 O2`) unchanged; dioxin rule fires 0 cuts on both.
  - 20-structure HCNO batch: **only `1097` changed**, 19 others byte-identical. Minimized coverage now **20/20**. 66 fragments, 0 `[QM WARNING]`, 0 odd-electron, deficit set unchanged (the same 7 intentional parity trade-offs plus the open `1102` Path D carbon).
  - `./run_fast_test.sh` 8/8 (MOF) and `--kind cof` 8/8.
- **Scope:** 50 of the 884 HCNO structures carry a dioxin-like `O-C-C-O` motif (upper bound; the real rule additionally requires the six-ring closure and two distinct aryl systems).
- **Follow-up risks:**
  - The "keep O with the larger aryl system" tie-break is a heuristic. It is correct for HHTP-type COFs (polyol node bigger than the halogenated linker) but would place the oxygens on the wrong building block in a hypothetical COF where the polyol unit is the smaller of the two. It only affects which library file owns the O — the assembled fragment contains the same atoms either way, since junctions reconnect geometrically.
  - The generic fallback `bo_node_species = {"B","O"}` heuristic remains chemically wrong for boron-free COFs and is still reachable by any O-containing COF that Path J cannot decompose. Note the entire HCNO-884 subset contains **no boron at all**, so that path is never legitimately applicable there.

## 2026-09-17 - UniFrag: Raised Minimize Threshold to 20 Atoms and Added a Node+One-Linker Minimize Fallback for Large-Node COFs
- **Problem:** The standard minimize strategy (keep one full linker, trim the rest to their first ring) cannot shrink COFs whose node dominates the fragment while the linkers are already small. `1031.cif` (node 39 atoms, linker 12) reduced by 0 atoms; 8 of 20 structures in the HCNO test batch produced no usable minimized fragment.
- **Scope correction (same session, at user's request): the reduction threshold is COF-ONLY and MOFs are entirely unaffected.** It was briefly applied to `_process_mof_file` as well; that was reverted so `_process_mof_file` is byte-identical to its previously committed form (verified: `git diff` over that function is empty). The constant was renamed `MIN_ATOM_REDUCTION_FOR_MINIMIZE` -> `MIN_ATOM_REDUCTION_FOR_MINIMIZE_COF`, and all four remaining uses are inside `_process_cof_file`. Reason: on the IRMOF series 12 of 19 structures clear a 20-atom bar by exactly one atom, so the rule is far too fragile to sit anywhere near the MOF path.
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY]
    - Added `MIN_ATOM_REDUCTION_FOR_MINIMIZE_COF = 20` (COF only; MOF minimize left unconditional as before).
    - `COFFragmenter.extract()` and `_try_coffragmentor_node_linker_fragment()` gained a `single_linker` flag. With `minimize=True, single_linker=True`, linker images after the first are dropped entirely (`continue`) instead of being trimmed to their first ring; the node's now-open attachment points are closed by the existing post-merge `_cap_severed_double_bond_sites` / `_cap_open_oxygens` passes. `extract()` returns `None` for `single_linker=True` when the coffragmentor node+linker route is unavailable, so the caller does not mistake an identical fragment for an improvement.
    - `_process_cof_file` now retries once with `single_linker=True` when the standard minimize falls short of the threshold, and keeps whichever candidate reduces more.
- **Design note:** The trigger is *measured outcome*, not a COF-type classification. An up-front `node_atoms > linker_atoms` rule was evaluated against the batch and mis-handles 4 of 20 cases — it would needlessly switch `856` and `1110` (where trimming already works) and would miss `137`/`557` (node and linker both 12 atoms, yet trimming achieves nothing). See the 2026-09-17 decision entry.
- **Validation:**
  - COF batch (20 HCNO structures): minimized-fragment coverage **12/20 -> 19/20**. New minimized fragments for `1031` (168->120), `211` (192->128), `167` (174->126), `1156` (145->103), `1023` (126->90), `142` (124->82), `557` (124->72), `137` (114->66); `1110` improved 93 -> 67. All 64 fragments even-electron, 0 `[QM WARNING]`. Deficit set unchanged apart from `557FragCofMin` inheriting the same intentional parity trade-off its normal fragment already had.
  - `1031FragCofMin` = 120 atoms `C60 H46 N10 O4` = exactly node (`C24 N3`) + one linker (`C6 N2 O2`) per layer, bilayer. Valence-complete, even-electron.
  - Min fragments reporting 2 connected components are the bilayer (pi-stacked dimer) cases, which is by design — the layers are not covalently bonded.
  - `./run_fast_test.sh` 8/8 (MOF) and `--kind cof` 8/8.
  - MOF impact of the threshold change checked on the IRMOF series (18 structures): **no MOF minimized fragment was dropped**.
- **Follow-up risks:**
  - **Do not extend `MIN_ATOM_REDUCTION_FOR_MINIMIZE_COF` to MOFs.** Full IRMOF batch (19 structures) reductions: 21 x12, 25, 71 x3, 111, 121 x2, 137 — bimodal, with 12 of 19 clearing a 20-atom bar by exactly ONE atom and nothing between 25 and 71. Applying any threshold in that range to MOFs would silently delete most IRMOF minimized fragments. The constant is `_COF`-suffixed for this reason.
  - `1097.cif` lost its minimized fragment (reduction was exactly 15, now under threshold; it has no node/linker decomposition so the fallback cannot apply). Accepted as the intended effect of the stricter threshold.
  - The `single_linker` fallback is implemented only for the coffragmentor node+linker route. COFs routed through the native fallback paths (A/B/C/D) get no retry, so a structure there whose trimming falls short will simply have no minimized fragment.

## 2026-09-16 - UniFrag: Benzoxazole (Oxazole) COF Linkage Cleavage + Linkage-Based Node/Linker Classification (fixes `1031.cif`)
- **Symptom (user-reported, visual):** `1031.cif` produced a malformed normal fragment — the node was missing atoms and only ONE linker was attached. Its helper libraries were nonsense: `cof_nodes_lib/1031_00.xyz` was a single `O` atom and `cof_linkers_lib/1031_00.xyz` was 108 atoms (`H30 C66 N12`, i.e. the whole unit-cell framework minus its oxygens).
- **Diagnosis:** `1031` is a **benzoxazole-linked triazine COF** (per-cell: 2 triazine `C3N3` rings, 6 oxazole `C3NO` rings, 9 benzene rings; all 6 O are ether-type `C-O-C`, no boron anywhere). `coffragmentor.py` only knew two linkage chemistries — B-O and imine C=N — and the oxazole C-N bond is explicitly protected by the `in_small_ring` guard, so **zero** bonds were cut, Path J found no node/linker and returned `None`. Extraction then fell through to the generic fallback, where `bo_node_species = {"B", "O"}` treats EVERY oxygen as a node-forming atom. With no B-O connectivity, each isolated ether O became its own "node component of size 1" — hence `COF Path A (Single node). Node component size: 1`, a chemically meaningless single-atom node, and a lopsided one-arm fragment.
- **Changed files:**
  - `coffragmentor.py` [MODIFY] — Two changes:
    1. Added benzoxazole/oxazole cleavage in `COF.fragment()`: a carbon bonded to exactly one O and one N where removing that carbon leaves a 3-edge O-x-y-N path (i.e. the two are closed into a five-membered ring through it) is the oxazole C2 / former aldehyde carbon; both its C-O and C-N bonds are severed. This separates the aldehyde-derived block from the aminophenol-derived block.
    2. Added module-level `_count_linkages(cut_bonds)` and switched node/linker classification from `len(attachment_atoms) >= 3` to `_count_linkages(cut_bonds) >= 3`. One chemical linkage can sever more than one bond (a benzoxazole C2 loses both C-O and C-N), so severed bonds that share an attachment atom or share an external partner atom are grouped and counted once. For single-bond linkages (imine, boroxine) this is mathematically identical to the old count, so it is backward compatible.
- **Why the classification change was also required:** with only the cleavage rule added, the aminophenol-derived linker benzene has 4 attachment atoms (2 O + 2 N) and was misclassified as a NODE, leaving 5 nodes / 0 linkers, so `_try_coffragmentor_node_linker_fragment` still bailed out (`if not nodes or not linkers: return None`). Counting linkages instead gives it 2 -> LINKER, and the triazine side 3 -> NODE.
- **Validation:**
  - `1031.cif` now routes through Path J and yields `1031FragCof` = **168 atoms `C84 H54 N18 O12`** (was 71 atoms `C35 H27 N7 O2`). Per layer that is exactly 1 triazine node (`C24 N3`) + all 3 benzoxazole linkers (3 x `C6 N2 O2`), bilayer-stacked.
  - Helper libraries are now chemically correct: node = 39 atoms `H12 C24 N3` (triazine + 3 phenyl + 3 aldehyde C), linker = 12 atoms `H2 C6 N2 O2` (benzene + 2 O + 2 N). QM-ready exports `1031FragCofOnlyNode` (96 at) and `1031FragCofOnlyLinker` (36 at) are produced, both valence-complete and even-electron.
  - Junction integrity confirmed: of the 12 O in the main fragment, 6 have intact `C-O-C` (oxazole rings geometrically reformed at the node-linker junctions) and 6 are outward termini correctly capped to `-OH`. No double-capped junctions.
  - Prototyped on scratch copies BEFORE editing; verified the new cleavage rule cuts **0** extra bonds on imine COFs (`211.cif`, `COF-TpAzo.cif`) and that linkage-counting reproduces the old classification exactly for both.
  - 20-structure HCNO batch re-run: **only `1031` changed**; all 19 others byte-identical in atom count and formula. 0 `[QM WARNING]`. Audit: 57 fragments, 0 odd-electron, same 7 pre-existing flagged items as before (6 intentional parity trade-offs + the open `1102` Path D carbon).
  - `./run_fast_test.sh` (MOF) 8/8 and `--kind cof` 8/8.
- **Scope:** a motif scan over the full HCNO-884 subset finds **127 structures (~14%)** containing an oxazole-like `C(-O)(-N)` carbon, so this affected far more than one structure. (That count is an upper bound — the cheap scan omits the five-ring closure test the real rule applies.)
- **Follow-up risks:**
  - The generic fallback's `bo_node_species = {"B", "O"}` heuristic is still chemically wrong for any COF where O is an ether/carbonyl rather than part of a B-O node: it will happily produce a single-atom "node". `1031` no longer reaches it (Path J now succeeds), but other non-boron O-containing COFs that Path J cannot decompose will still hit this. Consider gating that heuristic on boron actually being present/bonded to O.
  - The oxazole rule is not extended to the isoelectronic benzothiazole (S) or benzimidazole (N-N) linkages; those COFs would still fail to decompose. Straightforward to add by the same pattern if needed.

## 2026-09-16 - UniFrag: Fixed Under-Capped Imine Cut Sites in COF Node/Linker Building Blocks + Added OnlyNode QM Export (`CoRE-COF DT1242-v7.0` HCNO subset, validated on `211.cif`)
- **Context:** Starting a new workstream to demonstrate UniFrag on the CoRE-COF database (`runUniFrag/CoRECOF/CoRE-COFs_DT1242-v7.0/`, 1242 structures). Built a 884-structure HCNO-only (C,H,N,O) subset (`runUniFrag/CoRECOF/CoRE-COF_HCNO884/`) as the primary demonstration set. A 20-structure random test batch surfaced valence-incomplete building blocks on visual inspection; `211.cif` was used as the deep-dive case.
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — All changes scoped to `COFFragmenter` + `_process_cof_file` + the COF-only flush closure in `main()`. Zero changes to `MOFFragmenter`, `BaseFragmenter`, `_process_mof_file`, or any MOF flush path (verified via `git diff` hunk audit).
    - Added `COFFragmenter._cap_severed_double_bond_sites`: detects C/N atoms with exactly 1 bonded HEAVY neighbor (the topological signature of an imine C=N cut, e.g. from `coffragmentor.py`'s `heavy_degree(c_idx)==2` cleavage rule) and caps them with `target_valence[sp] - (heavy_count + existing_H_count)` new hydrogens — a true valence deficit, not a fixed count, since the atom may already carry a native H (e.g. an aldimine `ArCH=N-` carbon keeps its own H after the `=N` side is cut) or none (e.g. a bare imine nitrogen `ArN=CH-`).
    - Applied this capping to the isolated node molecule and each isolated linker image *before* they are merged into the combined fragment in `_try_coffragmentor_node_linker_fragment` (merging first, then capping, made real deficits unrecoverable — a nearby unrelated atom from a different merged image could appear close enough to mask the true degree).
    - Added `COFFragmenter._make_cof_qm_ready`: a COF-only QM-ready capper (mirrors the shared `BaseFragmenter._make_qm_ready_linker`, left untouched) that runs `_cap_severed_double_bond_sites` then `_cap_open_oxygens`, then the existing shared `optimize_capped_h_geometry_only` / `fix_odd_electron_multiplicity` steps.
    - Added `self.extracted_nodes` tracking (mirrors `self.extracted_linkers`) in `_export_coffragmentor_library` and `_export_cof_component_library`.
    - `_process_cof_file` now also emits a QM-ready `OnlyNode` fragment per unique node (symmetric to the existing `OnlyLinker` export), and its `OnlyLinker` generation now calls `_make_cof_qm_ready` instead of the shared `_make_qm_ready_linker` (COF-only swap; MOF's `_process_mof_file` still calls the original shared method, unchanged).
    - `main()`'s COF flush closure now writes `<stem>FragCofOnlyNode` frames to the ExtXYZ collection alongside `FragCofOnlyLinker`.
- **Summary:**
  - Root cause: for imine-linked COFs, `coffragmentor.py` cleaves the C=N double bond, leaving both the amine-derived node's nitrogen and the aldehyde-derived linker's carbon with only 1 heavy neighbor. `_cap_open_oxygens` (the only capping applied to COF Path J's assembled fragment) only handles single-neighbor O, single-neighbor N, and 2-heavy-neighbor "phenyl edge" C — it has no branch for a 1-heavy-neighbor carbon, so those sites were left completely uncapped in the standalone `OnlyLinker`/`OnlyNode`-style exports, or partially capped (short by exactly the severed bond's second valence unit) in the main assembled fragment when a native H happened to already be present.
  - Verified on `211.cif` (imine COF, keto/amine node + dialdehyde-type linker): before the fix, `211FragCof`/`211FragCofMin` had 6/2 systematically under-coordinated carbons (valid bond order but missing one H), and the `OnlyLinker` helper export had 4 fully bare (zero-H) carbons. After the fix, all 4 exported fragments (`211FragCof` 216 atoms, `211FragCofMin` 204 atoms, `211FragCofOnlyLinker` 40 atoms, `211FragCofOnlyNode` [NEW] 96 atoms) have zero valence-deficient heavy atoms and even electron counts.
- **Validation:**
  - `./run_fast_test.sh` (MOF, default `--kind mof`): 8/8 passed.
  - `./run_fast_test.sh --kind cof`: 8/8 passed.
  - Manual geometric valence audit (covalent-radius bond graph via the project's own `COFFragmenter.is_valid_bond`) confirmed 0 remaining deficits across all 4 fragment types for `211.cif`.
- **Follow-up (same day):** The initial fix above placed each newly-added capping H via the existing clash-avoidance cone search only, with no angle target. For sites needing 2 new H (e.g. converting a bare imine N into -NH2, or an aldimine =CH- carbon into -CH3), this left the *first* H nearly collinear with the anchor bond (H-N-C measured at 180 degrees) instead of a proper pyramidal/tetrahedral angle. Fixed by having `_cap_severed_double_bond_sites` track the newly-added H indices per atom and immediately call the existing (previously unused/legacy per the 2026-05-07 decision) `BaseFragmenter.refine_h_geometry_with_rdkit` on just those indices — a local UFF relaxation with every other atom fixed. Only ever invoked from this COF-only method, so the "legacy/unused" status for MOF is unchanged. Verified post-fix: `-NH2` sites now measure ~107-108 degrees (H-N-H and H-N-C), `-CH3`-like sites measure ~108-112 degrees — both close to ideal tetrahedral. Valence-deficit and even-electron checks re-confirmed at 0 deficits for all 4 fragment types; `./run_fast_test.sh` (MOF) and `--kind cof` both still 8/8.
- **Second follow-up (same day, user caught via visual inspection):** The pre-merge capping (node and each linker image capped in isolation, before `append_merged`) was itself wrong. A node/linker cut atom that IS reconnected by the merge (e.g. the node's imine N and the attached linker's imine C, once placed at their correct relative crystallographic positions, sit ~1.3 A apart — a real restored C=N bond) must be left uncapped; capping both sides independently before merging produced two dead-end groups (a `-NH2` on the node, a `-CH3`-like stub on the linker) sitting next to each other instead of one real bond. This passed the earlier valence-deficit check (0 deficit on each side individually) despite being chemically wrong, which is why it wasn't caught until visual inspection. Fixed by removing both pre-merge calls (`_cap_severed_double_bond_sites(node_sp, ...)` and the per-linker-image call) from `_try_coffragmentor_node_linker_fragment`, relying solely on the existing post-merge call: after merging, a properly-reconnected junction naturally shows 2 heavy neighbors (skipped, correct) while a genuinely still-dangling atom (e.g. a linker's far end attaching to a second, non-included node copy) still shows exactly 1 heavy neighbor (correctly capped). Re-verified on `211.cif`: all 6 node nitrogens in `211FragCof` now show a real ~1.30 A bond to a linker carbon (0 dangling -NH2 at junctions); the linkers' genuinely-unconnected far ends are still correctly capped as -CH3-like termini with reasonable angles (~106-120 degrees). `211FragCof` atom count dropped from the (wrongly inflated) 216 to a correct 192; 0 valence deficits; `./run_fast_test.sh` (MOF) and `--kind cof` both still 8/8.
- **Third follow-up (same day, ATTEMPTED AND REVERTED — do not retry this way):** An attempt was made to extend the node/linker helper-library + QM-ready `OnlyNode`/`OnlyLinker` export to COF Path A/B, by hooking `_export_cof_component_library` into the `if (not minimize) and (path_mode in {"A","B"}) and node_atoms:` block in `extract()` (using `core_nodes` as the node and each `touches_core` component of `comps` as a linker). This is WRONG and was reverted. Two independent reasons:
  1. `core_nodes` there is only the B/O node-species seed and can collapse to a single atom, while `comps` are unbounded BFS components over every non-node atom in the whole supercell. For `1031.cif` this exported a 1-atom `O1` "node" and a 108-atom `H30 C66 N12` "linker" — i.e. the entire unit cell framework (`H30 C66 N12 O6`) minus its oxygens, as one "linker". For `411.cif` it exported a 250-atom "linker" (larger than the whole 83-atom fragment) plus a lone stray `H` that could not be parity-fixed (the run's only `[QM WARNING]`).
  2. For single-block Path A topologies (`len(core_nodes) <= 2`), that `final` set is discarded outright further down and rebuilt by the `single_block_keep_heavy` / `_extract_neighbor_single_block` branch — so anything exported at that point does not even correspond to the fragment that is written out.
  A correct implementation would have to hook into whichever branch actually produces the emitted fragment (for single-block Path A, the `single_block_keep_heavy` branch) and derive node/linker from a decomposition that is genuinely bounded to one repeat unit — not from the raw `core_nodes`/`comps` intermediate. Post-revert state re-verified: 58 fragments, 0 odd-electron, no degenerate (<=2 atom) helper-library files, `run_fast_test.sh` 8/8 for both `mof` and `cof`.
- **Follow-up risks:**
  - COF Path C (Tetra-C node) and Path D (Porphyrin core) produce NO `OnlyNode`/`OnlyLinker` QM-ready building blocks at all (they never call `_export_cof_component_library` / `_export_coffragmentor_library`). In the 20-structure HCNO test batch this affects `103.cif`, `760.cif`, `1102.cif`, `1120.cif`. Per user decision on 2026-09-16 this is accepted and intentionally out of scope — C/D are not required to produce QM-ready building blocks.
  - `1102FragCof` (Path D) still has one genuinely under-coordinated carbon (1 heavy neighbor + 2 H, needs 1 more H); unresolved, and NOT explained by the intentional parity-fix mechanism. Its minimized counterpart is clean.
  - The remaining flagged valence "deficits" on O/N sites (`411`, `557`, `637`, `851`) are NOT bugs: they are the documented `fix_odd_electron_multiplicity` behaviour deliberately removing one capping H from a low-priority site (alcohol/amine/carboxylate) to force an even electron count. Any future audit script should account for this before reporting them.
  - The new capping rule (`exactly 1 heavy neighbor` → deficit-based H count) is applied inside `_try_coffragmentor_node_linker_fragment` (the primary Path J route) exclusively as a POST-merge step now (pre-merge capping was removed — see follow-up above), and is also present post-merge in the two native fallback COF paths inside `extract()` (Path A/B/D family). The fallback-path insertions are UNTESTED against a real fallback-routed bug case — only Path J was validated end-to-end on `211.cif`.
  - The rule cannot distinguish bond order from geometry alone; a genuine nitrile group (`R-C#N`) would show the same "1 heavy neighbor" signature and could be incorrectly capped. Not observed in the HCNO-884 subset so far (`coffragmentor.py` only ever cleaves B-O and imine C-N bonds, so nitriles are never a cut product), but worth watching if nitrile-linked COFs are added later.
  - The full 20-structure HCNO test batch (`runUniFrag/CoRECOF/test_hcno20*/`) has NOT yet been re-run with this fix, per explicit user instruction to fully resolve `211.cif` first. The 884-structure HCNO batch has also not been (re-)run with the fix.
  - RDKit prints harmless `Molecule does not have explicit Hs` stderr warnings during this local refinement; cosmetic only, not investigated further.

## 2026-09-12 - UniFrag: Added Workflow Overview Figure to README.md
- **Changed files:**
  - `assets/unifrag_workflow.png` [NEW] — Added the UniFrag workflow architecture and overview diagram.
  - `README.md` [MODIFY] — Embedded the workflow figure centered under the header and project description.
- **Summary:**
  - Saved the user-provided high-level workflow diagram to `assets/unifrag_workflow.png` and embedded it into `README.md`.
- **Validation:**
  - Verified Markdown image rendering and path.

## 2026-07-22 - UniFrag: COF Minimization Mode Keeps 1 Full Linker + First-Ring Attached Linkers (`COF_TpAzo`)
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Updated `_try_coffragmentor_node_linker_fragment` when `minimize=True` to keep **1 full attached linker** (first/primary arm) and trim all remaining attached linkers around the node to their **first connected ring**.
  - `runUniFrag/test_for_paper/COF_TpAzo/cifs/fragmentation_summary.csv` [MODIFY] — Re-generated summary CSV for `COF-TpAzo.cif`.
  - `runUniFrag/test_for_paper/COF_TpAzo/cifs/fragments_collection.extxyz` [MODIFY] — Re-generated ExtXYZ collection for `COF-TpAzo.cif`.
- **Summary:**
  - Updated COF minimization mode to match MOF minimization behavior: exactly **1 attached linker remains FULL**, while all other attached linkers around the central node are trimmed to their first connected ring.
  - Re-generated `COF-TpAzo` fragments:
    - `COF-TpAzoFragCofMin`: 144 atoms (`C66 H56 N16 O6`), 2 stacked layers (dimer mode). Each 72-atom layer contains **5 six-membered rings**: 1 central ketoenamine/Tp ring + 2 rings on Arm 1 (the 1 FULL azobenzene linker) + 1 ring on Arm 2 (trimmed) + 1 ring on Arm 3 (trimmed).
    - Close pairs (< 0.8 Å): **0**.
    - Monomer electron count: 302 e- after QM fix (`68` atoms `C33 H24 N8 O3`, even count, closed-shell singlet, 100% QM ready).
- **Validation:**
  - Fast test suite (`./run_fast_test.sh`): **8/8 passed**.
- **Follow-up risks:**
  - None.

## 2026-07-22 - UniFrag: Updated COF Minimization to Include First Ring Around Node for All Arms (`COF_TpAzo`)
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Moved `_get_first_ring_keep_heavy` and `_first_connected_ring_fragment` to `BaseFragmenter`. Updated `_try_coffragmentor_node_linker_fragment` and `_try_cof_graph_node_linker_fragment` to select best image shifts for ALL attached linkers around the node and trim each linker to its first connected ring when `minimize=True`.
  - `runUniFrag/test_for_paper/COF_TpAzo/cifs/fragmentation_summary.csv` [MODIFY] — Summary CSV updated for `COF-TpAzo.cif`.
  - `runUniFrag/test_for_paper/COF_TpAzo/cifs/fragments_collection.extxyz` [MODIFY] — ExtXYZ collection updated for `COF-TpAzo.cif`.
- **Summary:**
  - Diagnosed and updated COF minimization logic: previously when `minimize=True`, `_try_coffragmentor_node_linker_fragment` selected only 1 linker arm (`scored_images[:1]`) and kept that 1 linker full (all 26 atoms), leaving the other arms unattached.
  - Updated the COF minimization pipeline so that for every attached arm around the node, the best periodic image is selected and trimmed to its first connected ring using `_first_connected_ring_fragment`.
  - Verified on `runUniFrag/test_for_paper/COF_TpAzo/cifs/COF-TpAzo.cif`:
    - `COF-TpAzoFragCofMin`: 120 atoms (`C54 H48 N12 O6`), 2 connected stacked layers (dimer mode). Each 60-atom layer contains 4 six-membered rings: 1 central ketoenamine/Tp ring + 3 attached first-ring phenyl arms (1 per arm).
    - Close pairs (< 0.8 Å): **0**.
    - Electron count for 54-atom monomer: 246 e- (even, closed-shell singlet, 100% QM ready).
- **Validation:**
  - Fast test suite (`./run_fast_test.sh`): **8/8 passed**.
- **Follow-up risks:**
  - None.

## 2026-07-22 - UniFrag: Restored Single-H Capping for Raw CIF Missing Hydrogens (`LOTTEW.cif`)
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Added a single-H completion pass in `_extract_sbu_cluster` specifically targeting undercoordinated aromatic ring carbons (`len(c_nbs) == 2 and total_val == 2`).
  - `test_for_phenyle_h_cap/fragmentation_summary.csv` [MODIFY] — Summary for `LOTTEW.cif` test.
  - `test_for_phenyle_h_cap/fragments_collection.extxyz` [MODIFY] — ExtXYZ collection for `LOTTEW.cif`.
- **Summary:**
  - Diagnosed why `LOTTEW.cif` had missing Hydrogens on phenyl rings:
    1. `LOTTEW.cif` is a crystallographic file in CSD where the original X-ray refinement omitted Hydrogen atoms on 15 aromatic ring carbons.
    2. Added a safe, single-H capping pass that detects sp2 ring carbons with `len(c_nbs) == 2` and `total_val == 2` (0 existing H and 0 heteroatoms), placing **exactly 1 C-H capping Hydrogen** along the outward bisector vector `base_vec = -(v1 + v2)`.
  - Re-tested `LOTTEW.cif` in `test_for_phenyle_h_cap/`:
    - `LOTTEW_frag_min.xyz`: `113` atoms (`C52 Mg3 O16 H42`).
    - **Single H ring carbons**: 36/36 (100% of all aromatic carbons have exactly 1 H).
    - **Double H ring carbons**: 0.
    - **Undercoordinated ring carbons**: 0.
  - Re-ran batch fragmentation, fragment size filtering, and multi-cutoff SOAP analysis for all 75 Mg MOFs.
- **Validation:**
  - `LOTTEW` in `mg_cr_cifs_noduplicated` updated to `C52 Mg3 O16 H42` with 0 undercoordinated carbons and 0 double-capped carbons.
  - Fast test suite (`./run_fast_test.sh`): **8/8 passed**.
- **Follow-up risks:**
  - None.


- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Removed the redundant `while True` auto-capping loop from `_extract_sbu_cluster`, added 0.3 Å coordinate deduplication in `bridge_cap_map`, and checked `idx not in {ba for ba, _ in bridge_atoms_to_cap}` to prevent duplicate capping entries in Minimize mode.
  - `test_mofs_for_node/fragmentation_summary.csv` [MODIFY] — Summary for node test MOFs.
  - `test_mofs_for_node/fragments_collection.extxyz` [MODIFY] — ExtXYZ collection for node test MOFs (`IRMOF-10` and `Mg2_dobpdc_CoRE_ASR`).
- **Summary:**
  - Diagnosed and fixed double-capping (2 Hydrogens placed on phenyl carbons) in Minimize mode:
    1. Duplicate entries in `bridge_atoms_to_cap` caused `place_capping_h` to run twice per cut carbon. Added duplicate prevention check (`idx not in {ba for ba, _ in bridge_atoms_to_cap}`) and vector deduplication in `bridge_cap_map`.
    2. Removed the extra `while True` loop that forced a 2nd capping H onto already-capped sp2 carbons.
  - Re-tested `IRMOF-10.cif` in `test_mofs_for_node/`:
    - `IRMOF-10FragMof`: `173` atoms (`C84 H60 O25 Zn4`), **`Zn4O` central node**, **0 double-capped carbons**, **0 undercoordinated carbons**.
    - `IRMOF-10FragMofMin`: `102` atoms (`C49 H34 O15 Zn4`), **`Zn4O` central node**, **0 double-capped carbons**, **0 undercoordinated carbons**.
  - Re-ran batch fragmentation, large fragment filtering, and multi-cutoff SOAP analysis for all 75 Mg MOFs.
- **Validation:**
  - Audited all 113 retained frames across the Mg dataset: **0 double-capped phenyl ring carbons** across all aromatic MOFs.
  - Multi-cutoff SOAP similarity: $r_{\text{cut}} = 3.0\text{ \AA}$ (99.25% Highly Rep.), $r_{\text{cut}} = 4.0\text{ \AA}$ (98.50% Highly Rep.), $r_{\text{cut}} = 5.0\text{ \AA}$ (89.16% Highly Rep.).
  - Fast test suite (`./run_fast_test.sh`): **8/8 passed**.
- **Follow-up risks:**
  - None.


- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Added `_unwrap_hydrogens` during CIF loading, implemented `_get_mic_vector` for exact bond vector unwrapping in `_extract_sbu_cluster`, updated `min_o_contact=1.1` in `place_capping_h`, removed BFS over-expansion into partial linkers on edge metals in Path B, and added an automatic single-pass capping loop for undercoordinated ring carbons (`len(c_nbs) >= 2 and total_val < 3`).
  - `runUniFrag/mg_cr_cifs_noduplicated/fragmentation_summary.csv` [MODIFY] — Re-generated summary CSV for all 75 Mg MOFs.
  - `runUniFrag/mg_cr_cifs_noduplicated/fragments_collection.extxyz` [MODIFY] — Re-generated ExtXYZ collection with complete H-capping on all phenyl rings.
  - `runUniFrag/mg_cr_cifs_noduplicated/mg_soap_distribution_*.png` [MODIFY] — Re-generated multi-cutoff SOAP similarity plots.
- **Summary:**
  - Diagnosed and resolved the root causes for missing H / undercoordination on phenyl rings (such as in `LOTTEW.cif`):
    1. **Periodic Boundary H-Wrapping in CIFs**: Crystallographic CIF files frequently store Hydrogen fractional coordinates wrapped across unit cell boundaries. Added `_unwrap_hydrogens` in `_load_clean_structure` and implemented `_get_mic_vector` to calculate exact minimum-image convention bond vectors.
    2. **Ortho-Hydroxyl Contact Rejection**: `place_capping_h` was defaulting `min_o_contact=1.5` Å, which silently rejected candidate capping H positions on aromatic ring carbons near ortho-hydroxyl or carboxylate oxygens (~1.3-1.4 Å away). Adjusted default `min_o_contact` to `1.1` Å.
    3. **Path B Topology Over-Expansion**: In Path B (Infinite SBU), `comp_queue` was pushing non-metal coordinating oxygens back into the queue, causing the search to wander into partial linkers of edge metals and leave truncated 2-coordinate carbons. Restricted edge metal completion to immediate coordinating non-metals.
    4. **Automatic Ring Carbon Valence Completion**: Added a pass in `_extract_sbu_cluster` to detect any sp2 ring carbon with `len(c_nbs) >= 2` and `total_val < 3` and automatically cap it with a Hydrogen along the outward bisector vector.
  - Re-ran batch fragmentation, fragment size filtering (>200 atoms), and multi-cutoff SOAP analysis for all Mg MOFs.
- **Validation:**
  - `LOTTEWFragMofMin` formula updated from `C52 Mg3 O16 H26` (16 undercoordinated carbons) to `C52 Mg3 O16 H46` with **exactly 0 undercoordinated carbons**.
  - SOAP Highly Represented (>= 0.98) similarity jumped to **99.25%** at r_cut=3.0 Å (up from 81.12%) and **98.50%** at r_cut=4.0 Å!
  - Ran `./run_fast_test.sh`: 8/8 regression tests passed cleanly.
- **Follow-up risks:**
  - None.


- **Changed files:**
  - `runUniFrag/mg_cr_cifs_noduplicated/fragmentation_summary.csv` [MODIFY] — Re-generated summary CSV for all 75 Mg MOFs.
  - `runUniFrag/mg_cr_cifs_noduplicated/fragments_collection.extxyz` [MODIFY] — Re-generated ExtXYZ collection (113 retained frames, 26 eliminated >200 atoms).
  - `runUniFrag/mg_cr_cifs_noduplicated/mg_soap_distribution_*.png` [MODIFY] — Re-generated multi-cutoff SOAP similarity plots.
- **Summary:**
  - Re-ran batch fragmentation on all 75 Mg-based MOFs in `runUniFrag/mg_cr_cifs_noduplicated` using 6 parallel processes.
  - Applied large-fragment filtering (`filter_large_fragments.py`) to eliminate fragments with >200 atoms (113 frames retained).
  - Re-calculated SOAP descriptors and similarity metrics across 3.0, 4.0, 5.0, and 6.0 Å cutoffs.
- **Validation:**
  - Confirmed 0 close pairs (< 0.8 Å) across the entire ExtXYZ collection via `verify_no_close_pairs.py`.
  - Ran fast regression test suite (`./run_fast_test.sh`); 8/8 tests passed successfully.
- **Follow-up risks:**
  - None.

## 2026-07-21 - UniFrag: Resolve Multiple Phenyl Ring Trimming Bug in Minimized Fragments
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Updated `_get_first_ring_keep_heavy` to find the nearest ring for each entry point in `first_layer` instead of breaking at the first ring target found.
- **Summary:**
  - Resolved a bug where ligands containing multiple phenyl rings or coordinating groups (such as in `DIWNAA.cif`) had all rings except the first one trimmed away during minimization, leaving carbons in other rings undercoordinated and capped with invalid hydrogens (e.g. `C[15]` in `DIWNAAFragMofMin` capped with two hydrogens).
  - The single-target BFS in `_get_first_ring_keep_heavy` was replaced with a loop that finds and keeps the nearest ring for each entry point in `first_layer`.
  - Re-ran the batch fragmentation, large fragment filtering, and SOAP similarity analysis on the entire Mg dataset. All minimized fragments are now structurally complete (e.g. `DIWNAAFragMofMin` formula updated from `C28 H22 Mg1 O11` to `C33 H24 Mg1 O11` with the second phenyl ring fully preserved).
- **Validation:**
  - Verified `fragments_collection.extxyz` using `verify_no_close_pairs.py` and `check_broken_phenyl.py`. Confirmed exactly **zero** close pairs (< 0.8 Å) in the entire dataset, and the undercoordination warnings for `DIWNAA` are completely resolved.
  - Ran the fast test suite (`./run_fast_test.sh`); 8/8 tests passed successfully.
- **Follow-up risks:**
  - None.

## 2026-07-21 - UniFrag: Resolve Close Hydrogen Capping and Missing Carbons in Mg-MOFs
- **Changed files:**
  - None.
- **Summary:**
  - Investigated and resolved the issue where some Mg-based MOFs in the `fragments_collection.extxyz` collection contained extremely short C-H bonds (~0.4 Å) or missing carbons on phenyl rings.
  - Traced the root cause of the close pairs to a post-processing filter script (`filter_large_fragments.py`) bug, which skipped backing up when the backup file (`_original.extxyz` / `_original.csv`) already existed on disk. Consequently, it read from the stale backup file of a previous session, filtered it, and overwrote the fresh new collection, thereby resurrecting old disordered coordinates.
  - Deleted the stale backup files, re-ran the parallel batch fragmentation on all 75 Mg-based MOFs in the folder `runUniFrag/mg_cr_cifs_noduplicated/cifs/`, successfully regenerating the 113-atom filtered collection.
  - Re-ran the multi-cutoff SOAP coordination environment analysis loop (3.0, 4.0, 5.0, and 6.0 Å), regenerating all UMAP/PCA plots and markdown reports.
  - Verified that all highlighted "missing carbons" were in fact heterocyclic/aliphatic linkers (such as pyridine in `AVIPAX`, thiophene in `PUZJOL`, or adipic acid in `IFASOZ`) which are structurally correct and fully intact.
- **Validation:**
  - Scanned the entire finalized `fragments_collection.extxyz` collection file for close pairs (atom-atom distance < 0.8 Å) using a dedicated verification script, confirming exactly **zero** close pairs across all 113 frames.
  - Ran `./run_fast_test.sh` and confirmed all 8/8 regression tests passed successfully.
- **Follow-up risks:**
  - None.

## 2026-07-21 - UniFrag: Clean Crystal Structure Disorder and Fix Supercell-Wrapping Coordinate Bug
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Added `_load_clean_structure` helper to filter out disordered/overlapping sites within 0.8 Å (periodic distance) right after loading structure files. Replaced raw `Structure.from_file` calls in extraction methods with `self._load_clean_structure`. Removed buggy fractional coordinate supercell-wrapping code (`dfrac -= np.round(dfrac)`) from molecular rebuilding logic, replacing it with direct relative Cartesian offsets from neighbor lists.
- **Summary:**
  - Resolved issues where minimized fragments for MOFs with disordered/multi-orientation linker components (such as `IRMOF-10.cif`) had overlapping atoms (separated by < 0.05 Å) and distorted capping geometries. The CIF files themselves contained disordered orientations modeling mutually perpendicular positions. Implementing a 0.8 Å periodic distance filter on initial structure loading cleans all disorder overlaps. Eliminating the supercell-level fractional coordinate wrapping logic prevents the fragment builder from wrapping atoms across supercell boundaries and creating coordinate collisions.
- **Validation:**
  - Validated on `test_mofs_for_node/IRMOF-10.cif`: The minimized fragment `IRMOF-10FragMofMin` is now correctly extracted with 102 atoms (no disordered duplicates), containing 5 cleanly capped benzoic acid minimized linkers and 1 full biphenyl linker. Minimum atom-atom distance in the fragment is completely clean with 0 close pairs (< 0.8 Å).
  - Ran the fast regression test suite (`./run_fast_test.sh`); all 8 tests passed successfully.
- **Follow-up risks:**
  - None.

## 2026-07-20 - UniFrag: Resolve Duplicate Carbons in Fallback Node Export
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Refined `_fallback_export_mof_node_linker` to trace bridging carbons directly from coordinating O/N atoms and merge duplicate/too-close coordinates using a 0.8 Å tolerance.
- **Summary:**
  - Solved the issue of duplicate carbon atoms (separated by ~0.6 Å) in exported SBU nodes. The duplicate carbons occurred because the metal-carbon distance of 2.55 Å was classified as coordinating and added in the first shell, while also being captured in the bridging check under a different coordinate shift. Using a 0.8 Å merge tolerance safely clusters these overlaps.
- **Validation:**
  - Validated on `runUniFrag/test_for_paper/mof74_2/cifs/Mg2_dobpdc_CoRE_ASR.cif`: Helper SBU node `mof_nodes_lib/Mg2_dobpdc_CoRE_ASR_00.xyz` is now correctly generated as `Mg3 C2 O11` (16 atoms) with a minimum atom-atom distance of 1.26 Å (zero close pairs).
  - Ran the fast regression test suite (`./run_fast_test.sh`); all 8 tests passed successfully with 0 failures.
- **Follow-up risks:**
  - None.

## 2026-07-20 - UniFrag: Complete SBU Edge Metal Coordination in Fallback Node Export
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Modified `_fallback_export_mof_node_linker` to use coordinate-centric neighbor addition and bridging carbon extraction.
- **Summary:**
  - Solved the issue where edge metal atoms in periodic SBU nodes were missing coordinating oxygens (due to index-centric representation mapping only a single periodic boundary shift). The updated SBU node builder explicitly places all coordinating non-metal neighbors around each metal in the SBU at their correct Cartesian positions, then merges duplicates.
- **Validation:**
  - Validated on `runUniFrag/test_for_paper/mof74_2/cifs/Mg2_dobpdc_CoRE_ASR.cif`: Node helper `mof_nodes_lib/Mg2_dobpdc_CoRE_ASR_00.xyz` is now correctly generated as the complete `Mg3 C12 O11` cluster (26 atoms) where all 3 Mg atoms have exactly 5 coordinated oxygens.
  - Ran the fast regression test suite (`./run_fast_test.sh`); all 8 tests passed successfully with 0 failures.
- **Follow-up risks:**
  - None.

## 2026-07-20 - UniFrag: Include Bridging Carboxylate Carbons in Fallback Node Export
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Modified `_fallback_export_mof_node_linker` to identify and include non-metal atoms (like carboxylate carbons) that link at least two coordinating non-metal atoms (like edge oxygens) already in the SBU.
- **Summary:**
  - Expanded SBU node fallback exports to include the bridging carboxylate carbon atoms that link coordinated oxygen pairs. The algorithm identifies any non-metal atom bonded to at least two distinct coordinated SBU atoms and includes them periodic-image safely in the node.
- **Validation:**
  - Validated on `test_mofs_for_node/IRMOF-10.cif`: Node helper `mof_nodes_lib/IRMOF-10_00.xyz` is now correctly generated as the complete `Zn4 C6 O13` cluster (23 atoms).
  - Ran the fast regression test suite (`./run_fast_test.sh`); all 8 tests passed successfully with 0 failures.
- **Follow-up risks:**
  - None.

## 2026-07-20 - UniFrag: Refine SBU Node Extraction and Library Export
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Modified `_fallback_export_mof_node_linker` to allow metal-metal connections under 3.6 Å during SBU adjacency graph construction. Updated `_export_moffragmentor_library` to merge neighboring moffragmentor nodes within 3.5 Å (with periodic images) before writing them to the helper node library.
- **Summary:**
  - Solved the issue where discrete polynuclear SBUs like the `Zn4O` core of `IRMOF-10` were exported to `mof_nodes_lib/` as mononuclear `Zn1 O4` (ZO4) fragments. SBU adjacency graph building now connects metal-metal pairs within 3.6 Å, allowing the SBU BFS traversal to find all metals in the cluster. Moffragmentor library export now groups and merges raw node components under 3.5 Å to prevent mononuclear node splitting.
- **Validation:**
  - Validated on `test_mofs_for_node/IRMOF-10.cif`: Node helper `mof_nodes_lib/IRMOF-10_00.xyz` is now correctly generated as the complete `Zn4 O13` cluster.
  - Ran the fast regression test suite (`./run_fast_test.sh`); all 8 tests passed successfully with 0 failures.
- **Follow-up risks:**
  - None.

## 2026-07-17 - UniFrag: Added QM-Ready "OnlyLinker" Fragment Export & Boundary Hydrogen Stripping
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Moved `_clean_linker_molecule` to `BaseFragmenter` with `METALS = set()` default class attribute. Implemented boundary hydrogen stripping inside `_clean_linker_molecule` by identifying boundary atoms where periodic heavy degree > linker heavy degree and deleting any bonded Hydrogen atoms. Ensured periodic structure is loaded early via `self.structure = struct`. Passed `orig_indices` in all node/linker component export paths. Added instance attribute `self.partner_vec` tracking bilayer shift vector in `COFFragmenter` and replicated bilayer stacking for `OnlyLinker` fragments in `_process_cof_file`.
  - `project-decisions.md` [MODIFY] — Updated decision entry.
- **Summary:**
  - Keeps helper library files (e.g. in `cof_linkers_lib/` or `mof_linkers_lib/`) completely uncapped by stripping crystallographic N-H/O-H boundary hydrogens, while preparing fully H-capped, force-field optimized, and spin-corrected QM-ready versions (labeled `OnlyLinker`) for the main ExtXYZ collection. If the output cluster is bi-layer, the QM-ready linker is automatically generated as stack-aligned bi-layer as well.
- **Validation:**
  - Verified on `test_cofs_for_linker/COF-TpAzo.cif`:
    - Raw template in `cof_linkers_lib/COF-TpAzo_00.xyz` is now correctly stripped to `H8 C12 N4` (uncapped, 24 atoms).
    - QM-ready fragment `COF-TpAzoFragCofOnlyLinker` in `test_cofs_for_linker/fragments_collection.extxyz` is correctly bi-layer H-capped to `H20 C24 N8` (52 atoms, bilayer).
  - Verified on `test_for_node/Mg2_dobpdc_CoRE_ASR.cif`:
    - Generates `Mg2dobpdcCoREASRFragMofOnlyLinker` with formula `H12 C14 O6` (32 atoms, neutral dicarboxylic dihydroxybiphenyl acid).
  - Verified on `mofs_from_raspa/MgMOF-74.cif`:
    - Generates `MgMOF-74FragMofOnlyLinker` with formula `H8 C8 O6` (22 atoms, neutral dobdc diacid).
  - Ran the regression test suite `run_cof_family.sh` on the COF 1xx series; all tests successfully passed.
- **Follow-up risks:**
  - None.

## 2026-07-17 - UniFrag: Cleaved Imine C=N Bonds Directly to Protect Terminal Functional Atoms in COFs
- **Changed files:**
  - `coffragmentor.py` [MODIFY] — Refined `COF.fragment()` C-N bond cleavage heuristic to target C-N bonds where Carbon has heavy-atom degree exactly 2 (representing the imine `C=N` double bond), rather than degree >= 3.
  - `project-decisions.md` [MODIFY] — Updated decision entry.
- **Summary:**
  - Cleaved the imine `C=N` double bond directly instead of cutting aromatic C-N single bonds. This ensures that the terminal amine/imine Nitrogen atoms stay with the amine-derived linkers and terminal formyl Carbon/Oxygen atoms stay with the aldehyde-derived nodes.
- **Validation:**
  - Verified on `test_cofs_for_linker/COF-TpAzo.cif`:
    - The extracted Azo linker is successfully generated with its complete formula **`H10 C12 N4`** (26 atoms), containing both the central azo and terminal amine Nitrogen atoms.
    - The extracted Tp node is generated with its complete formula **`H3 C9 O3`** (15 atoms), containing the formyl Carbons capped with Hydrogen.
  - Ran the regression test suite `run_cof_family.sh` on the COF 1xx series; all tests successfully passed.
- **Follow-up risks:**
  - None.

## 2026-07-09 - UniFrag: Restricted Minimized MOF Fragments to Exactly One Full Linker
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Restructured SBU minimization logic in `_get_fragment` to keep exactly one largest linker full and minimize all others.
  - `project-decisions.md` [MODIFY] — Added decision entry.
- **Summary:**
  - Restricted minimized fragments to contain exactly one full coordinating linker to guarantee minimal, low-cost cluster models for QM calculations, while truncating all other coordinating linkers.
- **Validation:**
  - Verified on `test_for_node/Mg2_dobpdc_CoRE_ASR.cif`: the minimized fragment size is successfully capped at 104 atoms (containing exactly 1 full 26-atom dobpdc linker and 5 truncated linkers), while the normal fragment size is 183 atoms.
  - Re-ran the fragmentation and post-processing size-filtering pipelines on the entire Mg dataset; confirmed all 101 fragments processed successfully, retaining 89 and eliminating 12 large normal fragments.
- **Follow-up risks:**
  - None.

## 2026-07-09 - UniFrag: Excluded Metal-Carbon and Metal-Hydrogen Bonds in MOF Node Extraction
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Refined `MOFFragmenter.is_valid_bond` to explicitly return `False` for bonds between metals and Carbon/Hydrogen.
  - `project-decisions.md` [MODIFY] — Added decision entry.
- **Summary:**
  - Prevented carboxylate carbons and capping hydrogens from being incorrectly included inside extracted SBU nodes by blocking metal-carbon and metal-hydrogen bond recognition in `is_valid_bond`.
- **Validation:**
  - Ran fragmentation on `test_for_node/Mg2_dobpdc_CoRE_ASR.cif` and confirmed the extracted SBU node is successfully cleaned from `Mg1 C3 O5` to inorganic `Mg1 O5` (and the dobpdc linker remains intact at `H6 C14 O6`).
  - Ran the refined pipeline on the entire Mg dataset and checked SOAP similarities. Excluding the carbon coordinate pollution increased the environment similarity match at a 6.0 Å cutoff significantly from **84.94% to 87.48%** (and from 98.55% to 99.25% at 3.0 Å).
- **Follow-up risks:**
  - None.

## 2026-07-09 - UniFrag: Swapped CrystalNN with JmolNN for Guest Removal and Refined Fallback Linker Extraction
- **Changed files:**
  - `runUniFrag/preprocess_dataset.py` [MODIFY] — Swapped `CrystalNN` with `JmolNN` as the primary connected component solver, and corrected the exception fallback criteria.
  - `runUniFrag/remove_guests.py` [MODIFY] — Swapped `CrystalNN` with `JmolNN` as the primary connected component solver, and corrected the exception fallback criteria.
  - `fragmentation_oop.py` [MODIFY] — Refined `_fallback_export_mof_node_linker` to exclude metal seeds and initialize `global_vis` as an empty set.
  - `project-memory.md` [MODIFY] — Added decision entry.
  - `project-decisions.md` [MODIFY] — Added decision entry.
- **Summary:**
  - Swapped `CrystalNN` with `JmolNN` as the default framework component solver during guest cleaning, preventing the framework from being fractured into pieces and silently deleted.
  - Excluded metal atoms from being selected as organic linker seeds.
  - Removed `v in expanded_node` from the BFS skip condition and initialized `global_vis` as empty to completely collect coordinating carboxylate and N-oxide oxygens.
- **Validation:**
  - Ran pre-processing on raw `EZEQEH`, `EZEQIL`, and `EZEQOR` parent structures and confirmed they are preserved intact (e.g. `Mg4 H32 C52 N4 O20` for `EZEQEH`).
  - Ran fragmentation on the preprocessed parent structures and `Mg2_dobpdc_CoRE_ASR.cif`, verifying that all linkers are extracted completely without missing atoms (dipyridine N,N'-dioxide is complete `C10 H8 N2 O2`; terephthalate is complete `C8 H4 O4`; and dobpdc is complete `C14 H6 O6`).
  - Placed all refined XYZ files in `test_for_linker/cifs/` for downstream testing.
  - Re-preprocessed and re-fragmented the entire Mg collection in `runUniFrag/mg_cr_cifs_noduplicated/`. Comparison against backed up collections showed:
    - 0 changes in metal-centered coordination fragments (89/89 are identical).
    - Systemic, dramatic improvements in the extracted linkers library: library size decreased from 46 to 42 unique linkers due to correct duplicate detection; large linkers were fully restored (+42 atoms for `DAJWET`); and missing coordinating atoms (+1 Oxygen/Nitrogen) were restored for dozens of linkers.
- **Follow-up risks:**
  - None.

## 2026-07-09 - CSD & CoRE Database: Extract and Purify Mg2(dobpdc) CIF Structure (RAVVUH)
- **Changed files:**
  - `Mg2_dobpdc.cif` [NEW] — Guest-cleaned CIF structure of Mg2(dobpdc) in the workspace root.
  - `Mg2_dobpdc_raw.cif` [NEW] — Original/unpurified CSD structure of Mg2(dobpdc) with guest molecules.
  - `Mg2_dobpdc_CoRE_ASR.cif` [NEW] — Computation-Ready (CoRE) ASR structure of Mg2(dobpdc) with DDEC6 charges.
  - `Mg2_dobpdc_CoRE_FSR.cif` [NEW] — Computation-Ready (CoRE) FSR structure of Mg2(dobpdc) with DDEC6 charges.
- **Summary:**
  - Used the CSD Python API `TextNumericSearch` to locate the CIF file for Mg2(dobpdc) (4,4'-dioxidobiphenyl-3,3'-dicarboxylate linker), identifying entry `RAVVUH` (common name: `IRMOF-74-II`).
  - Ran Phase 3 of `preprocess_dataset.py` (guest stripping via connected component analysis with JmolNN) to clean guest molecules from the structure.
  - Stripped 54 guest atoms (water molecules) from the pore space, reducing the formula from `Mg18 H54 C126 O144` (unpurified) to `Mg18 H54 C126 O90` (purified), while leaving the coordinated framework water molecules intact.
  - Identified and extracted the CoRE database counterparts (ASR and FSR subsets containing DDEC6 charges calculated by PACMAN) present under CSD-modified files.
- **Validation:**
  - Verified structure cell parameters and space group (`R-3`) match Deng et al. (Science 2012).
  - Verified guest stripping results using Pymatgen/CCDC connected component analysis, confirming all unbound water molecules are removed and only framework-coordinated waters remain.
  - Confirmed the CoRE ASR (`Mg6 H18 C42 O18`) has all solvents stripped (exposing open Mg sites), while the CoRE FSR (`Mg6 H18 C42 O30`) preserves coordinates for the metal-coordinated water molecules.
- **Follow-up risks:**
  - None.

## 2026-07-06 - Manuscript: Drafted Mg-MOF Extraction and Validation Results Section
- **Changed files:**
  - `results_mg_mofs.md` [NEW] — Artifact file containing the complete results and discussion section draft for the manuscript on Mg-based MOFs.
- **Summary:**
  - Drafted a detailed, publication-ready results section describing the pipeline for Mg-based MOF extraction, pre-processing, guest removal, dual-fragmentation (Normal vs. Minimized), size filtering, and multi-cutoff SOAP validation (3.0, 4.0, 5.0, and 6.0 Å).
  - Summarized execution statistics (5,226 scanned structures, 77 Mg-containing structures, 14 guest-stripped structures, 89 retained frames, 12 eliminated large normal frames, and 1 timeout).
  - Tabulated chemical formulas before and after guest removal and listed pruned large Normal fragments.
  - Summarized SOAP similarity and RMSD metrics across all four cutoffs for 551 Mg centers, detailing the physical causes of the observed trends (primary coordination preservation vs. outer shell capping changes).
- **Validation:**
  - Drafted and verified all values and metrics against the local report files (`preprocess_report.md`, `elimination_report.md`, and `mg_soap_analysis_*.md`).
- **Follow-up risks:**
  - None.

## 2026-06-25 - UniFrag: Multi-Cutoff SOAP Analysis (3.0, 4.0, 5.0, 6.0 Å) on Filtered Mg Dataset
- **Changed files:**
  - None.
- **Summary:**
  - Executed a multi-cutoff SOAP coordination environment analysis loop (`3.0`, `4.0`, `5.0`, and `6.0` Å) on the filtered Mg CR dataset.
  - Tracked change trends across different radii: tighter cutoffs show near-perfect coordination shell replication (98.55% highly represented and median similarity of 0.9997 at 3.0 Å), while wider cutoffs incorporate boundary environment effects (84.21% highly represented and median similarity of 0.9927 at 6.0 Å).
- **Validation:**
  - Executed `/Users/omert/miniconda3/bin/python runUniFrag/analyze_soap.py --metal Mg --r_cut 3.0 4.0 5.0 6.0 --brain_dir /Users/omert/.gemini/antigravity/brain/5cc8c9bc-d53e-489b-a932-9474aaa71491`.
  - Confirmed successful generation of reports (`mg_soap_analysis_{r_cut}.md`) and side-by-side plots (`mg_soap_distribution_{r_cut}.png`) for all four cutoffs in both destination and brain directories.
- **Follow-up risks:**
  - None.

## 2026-06-25 - UniFrag: Update SOAP Analysis on Filtered Mg Dataset
- **Changed files:**
  - None.
- **Summary:**
  - Re-ran the generalized SOAP analysis script on the filtered Mg CR dataset (now containing 89 frames after eliminating fragments > 200 atoms).
  - Calculated 551 parent descriptors and 154 fragment descriptors.
  - Verified that local coordination coverage was preserved perfectly: the highly represented (>=0.98 similarity) count remained exactly 84.21% (464/551 centers) because the corresponding minimized versions (<=200 atoms) of the eliminated large normal fragments were retained in the library.
- **Validation:**
  - Executed `/Users/omert/miniconda3/bin/python runUniFrag/analyze_soap.py --metal Mg --r_cut 6.0 --brain_dir /Users/omert/.gemini/antigravity/brain/5cc8c9bc-d53e-489b-a932-9474aaa71491`.
  - Confirmed successful regeneration of report `mg_soap_analysis_6.0.md` and combined figure `mg_soap_distribution_6.0.png` in both destination and brain directories.
- **Follow-up risks:**
  - None.

## 2026-06-25 - UniFrag: Implement Post-processing Script to Eliminate Large Fragments (>200 atoms)
- **Changed files:**
  - `runUniFrag/filter_large_fragments.py` [NEW] — Post-processing filtering script that reads the fragmentation CSV summary and parses the ExtXYZ collection to filter out and segregate large fragments (exceeding a user-specified threshold, default 200).
  - `project-memory.md` [MODIFY] — Documented run command examples for the filtering script.
- **Summary:**
  - Implemented `filter_large_fragments.py` supporting safety fallbacks for non-numeric/`TIMEOUT` entries in the CSV summary.
  - Automatically backs up original `.extxyz` and `.csv` files before filtering.
  - Moves excluded frames to a separate `eliminated_fragments.extxyz` file, updates the main `fragments_collection.extxyz` with retained frames, and appends `normal_eliminated` and `min_eliminated` flags to the CSV summary for downstream traceability.
  - Generated a detailed `elimination_report.md` specifying each eliminated structure, its parent REFCODE, atom count, formula, and reason for elimination.
- **Validation:**
  - Ran the script on the Mg CR dataset: `/Users/omert/miniconda3/bin/python runUniFrag/filter_large_fragments.py --metal Mg --max_atoms 200 --brain_dir /Users/omert/.gemini/antigravity/brain/5cc8c9bc-d53e-489b-a932-9474aaa71491`.
  - Scanned 101 frames; successfully retained 89 frames and eliminated 12 large frames (e.g., normal fragment of `DAJWET` at 556 atoms). Verified the report and the collection copy were successfully saved in the active conversation's brain folder.
- **Follow-up risks:**
  - None.

## 2026-06-25 - UniFrag: Implement Generalized Post-processing SOAP Fingerprint Analysis Script
- **Changed files:**
  - `runUniFrag/analyze_soap.py` [NEW] — General post-processing coordination environment analysis script supporting customizable target metals, input directories, ExtXYZ paths, and output destinations.
  - `project-memory.md` [MODIFY] — Documented new generalized post-processing SOAP commands.
  - `project-decisions.md` [MODIFY] — Logged the post-processing script architectural decision.
- **Summary:**
  - Generalized `analyze_zn_soap.py` to `analyze_soap.py` by parameterizing target metals, parent CIF directories, ExtXYZ paths, and destination directories.
  - Dynamically extracts elements from parent and fragment files to compile the unified SOAP species list.
  - Automates loop analyses for multiple cutoffs, generates side-by-side PCA/UMAP scatter plots, writes detailed Markdown reports, and copies results to the active conversation's brain folder.
- **Validation:**
  - Ran Zn verification sweep (`--metal Zn --r_cut 3.0 4.0 5.0 6.0`) on the full unmodified parent collection (1110 parents, 1418 fragments).
  - Ran Mg verification run (`--metal Mg --r_cut 6.0`) on the Mg collection (76 parents, 101 fragments).
  - Verified that all reports and PCA/UMAP distribution figures were correctly created in both the destination and brain artifacts directories with 0 errors.
- **Follow-up risks:**
  - None.

## 2026-06-25 - UniFrag: Extract Mg-Based CR MOFs and Support Dynamic Brain Directory Reporting
- **Changed files:**
  - `runUniFrag/preprocess_dataset.py` [MODIFY] — Added optional `--brain_dir` command-line argument to allow dynamically specifying where the pre-processing report is copied (falls back to the previous default path).
  - `project-memory.md` [MODIFY] — Added example commands for running `preprocess_dataset.py` for both Zn and Mg CR MOF extraction.
- **Summary:**
  - Extracted Mg-based single-metal CR MOFs from the deduplicated CR collection. Out of 5,226 scanned crystal structures, 77 Mg-based structures were identified and processed.
  - Of the 77 structures, 14 contained guest/solvent molecules that were successfully stripped (e.g. `XUFYAA` with 312 guest atoms removed), while 63 were intact.
  - The results were saved in `runUniFrag/mg_cr_cifs_noduplicated/`.
- **Validation:**
  - Ran the pipeline successfully with `/Users/omert/miniconda3/bin/python runUniFrag/preprocess_dataset.py --src_dir runUniFrag/cr_cifs_noduplicated --dest_dir runUniFrag/mg_cr_cifs_noduplicated --metal Mg --brain_dir /Users/omert/.gemini/antigravity/brain/5cc8c9bc-d53e-489b-a932-9474aaa71491`.
  - Confirmed 77 output `.cif` files and `preprocess_report.md` were correctly created, and the report copy was successfully written to the active conversation's brain folder.
- **Follow-up risks:**
  - None.

## 2026-06-25 - UniFrag: Relocate SOAP Reports and Plots to zn_cr_cifs_noduplicated Directory
- **Changed files:**
  - `runUniFrag/analyze_zn_soap.py` [MODIFY] — Updated all output file path variables (`output_md_path`, `output_png_path`, `default_output_md_path`, and `default_output_png_path`) to save files inside the `zn_cr_cifs_noduplicated/` subdirectory.
  - `project-memory.md` [MODIFY] — Updated decision paths.
  - `project-decisions.md` [MODIFY] — Updated decision paths.
- **Summary:**
  - Moved existing SOAP analysis markdown reports (`zn_soap_analysis*.md`) and PCA/UMAP plots (`zn_soap_distribution*.png`) from `runUniFrag/` to `runUniFrag/zn_cr_cifs_noduplicated/`.
- **Validation:**
  - Staged files in Git and verified that Git successfully detected the relocations as renames (`renamed: runUniFrag/zn_soap_analysis.md -> runUniFrag/zn_cr_cifs_noduplicated/zn_soap_analysis.md`, etc.).
- **Follow-up risks:**
  - None.

## 2026-06-24 - UniFrag: Multi-Cutoff SOAP Run Loop, Caching Optimization, and SOAP Vector RMSD Analysis
- **Changed files:**
  - `runUniFrag/analyze_zn_soap.py` [MODIFY] — Restructured script to support parsing a list of cutoffs (nargs="+") in argparse, pre-load and cache CIF structures in memory once at startup, compute raw SOAP fingerprint vector RMSDs between parent and best-matching fragment Zn centers, dynamically suffix output reports and plots by cutoff value, and update report templates with RMSD statistics.
  - `project-memory.md` [MODIFY] — Updated multi-cutoff execution examples and updated decisions log.
  - `project-decisions.md` [MODIFY] — Updated the SOAP environment coverage decision log.
- **Summary:**
  - **Caching Optimization**: Pre-loads parent CIFs into memory once as ASE Atoms objects at startup, avoiding redundant disk I/O and parsing overhead in loop. Achieved ~30x speedup for descriptor generation (1s vs 30s per cutoff).
  - **Fingerprint RMSD**: Calculated absolute Root-Mean-Square Deviation (RMSD) between raw SOAP fingerprint vectors of parent Zn environments and their best-matching fragment library counterparts.
  - **Dynamic Suffixed Names**: Dynamically named outputs as `zn_soap_analysis_{r_cut}.md` and `zn_soap_distribution_{r_cut}.png`. Maintains default names if single cutoff processed.
  - **Worst Matches Table**: Integrated UMAP, Cosine Similarity, and Fingerprint RMSD in the worst matches table and Executive Summary.
- **Validation:**
  - Ran `python runUniFrag/analyze_zn_soap.py --r_cut 3.0 4.0 5.0 6.0`. Loop completed in under 2 minutes, successfully generating suffixed markdown reports and side-by-side plots for all four cutoffs in `runUniFrag/` and brain artifacts directories.
- **Follow-up risks:**
  - None.

## 2026-06-24 - UniFrag: Configurable SOAP Cutoff and Side-by-side PCA/UMAP Projections
- **Changed files:**
  - `runUniFrag/analyze_zn_soap.py` [MODIFY] — Added `umap` projection to the dimensionality reduction step, updated the plot layout to show PCA and UMAP side-by-side, resolved remaining hardcoded cutoff references in the report template, and set default parameters for reproducible UMAP projections.
  - `project-memory.md` [MODIFY] — Documented user-configurable commands and updated decision logs.
  - `project-decisions.md` [MODIFY] — Updated the SOAP environment coverage decision log.
- **Summary:**
  - Installed `umap-learn` package (version `<=0.5.6` to preserve compatibility with `scikit-learn 1.5.2` and `csd-python-api`).
  - Added UMAP dimension reduction calculation using cosine distance metric (matches cosine similarity analysis).
  - Updated visualization from a single PCA plot to a side-by-side 1x2 PCA and UMAP subplot grid in `zn_soap_distribution.png`.
  - Updated the markdown report template in the analysis script to dynamically insert the chosen cutoff in all text fields and included discussion on PCA vs UMAP.
- **Validation:**
  - Ran the analysis with default cutoff (`--r_cut 6.0`) and custom cutoff (`--r_cut 5.0`). Both finished successfully, writing reports and combined figures in `runUniFrag/` and brain artifacts directories.
- **Follow-up risks:**
  - None.

## 2026-06-24 - UniFrag: continuous local SOAP fingerprint analysis for Zn environment
- **Changed files:**
  - `runUniFrag/analyze_zn_soap.py` [NEW] — Created script to compute Zn SOAP fingerprints in parent crystals and fragment libraries, performing similarity matching and PCA projection.
  - `project-memory.md` [MODIFY] — Added run commands and design decisions.
  - `project-decisions.md` [MODIFY] — Added structural coverage design decision.
- **Summary:**
  - Installed `dscribe` to compute SOAP local atomic descriptors around all Zn centers with a `6.0 A` cutoff.
  - Periodic context is used for parent crystals (`periodic=True`) and non-periodic context for capped fragments (`periodic=False`).
  - Measures cosine similarities: average similarity is `0.9908` and median is `0.9957` (84.29% highly represented at similarity >= 0.98), proving excellent structural representation of local Zn environments.
- **Validation:**
  - Run `analyze_zn_soap.py` on the local Zn dataset (1,110 parent structures and 1,418 fragments).
  - Verified successful outputs of `zn_soap_analysis.md` and PCA plot `zn_soap_distribution.png` in both `runUniFrag/` and brain directories.
- **Follow-up risks:**
  - None.

## 2026-06-24 - UniFrag: Configurable Processing Timeout and Relocation of Timed-Out Structures
- **Changed files:**
  - `fragmentation_oop.py` [MODIFY] — Added `_timeout_context` context manager using Unix `signal.alarm`, added `--timeout` parameter, wrapped workers in the context manager, caught `TimeoutError` explicitly, and added structure file relocation behavior.
  - `project-memory.md` [MODIFY] — Added a decision entry to the Decisions Log.
  - `project-decisions.md` [MODIFY] — Added the timeout design decision entry.
- **Summary:**
  - Implemented a configurable timeout (default 300s) to interrupt worker processes during complex crystal structure extraction.
  - Custom `TimeoutError` inherits from `BaseException` to propagate properly through internal `try-except Exception:` blocks.
  - On timeout, the worker automatically moves the offending structure file (`.cif` or `.pdb`) to a `timed_out_structures` folder in its directory, and reports `"TIMEOUT"` status for atoms and formulas in the summary CSV.
- **Validation:**
  - Ran a normal test on `IRMOF-1.cif` with a 10s timeout (passed).
  - Forced a timeout with a 1s limit on a copied structure and verified it caught the error, moved the file to `timed_out_structures`, and recorded `"TIMEOUT"` values in the CSV summary.
- **Follow-up risks:**
  - `signal.alarm` is Unix/macOS only; if run on Windows, it will behave as if no timeout is set (timeout ignored). This is fine since the user's OS is macOS.

## 2026-06-24 - UniFrag: Rename Script to prepare_linker4qm.py
- **Changed files:**
  - `runUniFrag/prepare_linker_extxyz.py` -> `runUniFrag/prepare_linker4qm.py` [RENAME] — Renamed the post-processing script.
  - `project-memory.md` [MODIFY] — Updated execution command paths to use the new name.
- **Summary:**
  - Renamed `prepare_linker_extxyz.py` to `prepare_linker4qm.py` to align with script's ultimate goal (preparing linkers for QM calculations).
- **Validation:**
  - Ran the renamed script successfully on the `mof_linkers_lib/` dataset and confirmed identical output behavior.
- **Follow-up risks:**
  - None.

## 2026-06-24 - UniFrag: Linker QM-Ready Summary CSV Output
- **Changed files:**
  - `runUniFrag/prepare_linker_extxyz.py` [MODIFY] — Added CSV summary generation reporting metadata for each processed structure.
  - `project-memory.md` [MODIFY] — Updated command instructions to reference the output summary CSV.
- **Summary:**
  - Added code to collect detailed metadata for each parsed structure in `prepare_linker_extxyz.py`, including `filename`, `label`, `num_atoms`, `num_heavy_atoms`, `num_hydrogens`, `formula`, `metals_stripped`, `atoms_capped`, `qm_fixed` status, and processing `status` / `reason`.
  - The script now writes a matching summary CSV file alongside the `.extxyz` file (e.g. `mof_linkers_lib/linkers_collection_summary.csv`).
- **Validation:**
  - Executed the script on `mof_linkers_lib/` with python.
  - Verified successful generation of `/Users/omert/Desktop/UniFrag_main/UniFrag/mof_linkers_lib/linkers_collection_summary.csv` with 222 data rows and the requested column metrics.
- **Follow-up risks:**
  - None.

## 2026-06-24 - UniFrag: Linker QM-Ready ExtXYZ Post-Processing Script
- **Changed files:**
  - `runUniFrag/prepare_linker_extxyz.py` [NEW] — Post-processing script that converts a linker .xyz library folder (produced by moffragmentor helper fragmentation) into a single QM-ready ExtXYZ collection file (`linkers_collection.extxyz`).
  - `project-memory.md` [MODIFY] — Added `prepare_linker_extxyz.py` command to the Run section.
- **Summary:**
  - The script reads every `.xyz` file in a given folder, applies even-electron QM-fix (mirrors `fix_odd_electron_multiplicity` from `fragmentation_oop.py`), deduplicates by heavy-atom composition key, and writes each unique linker as an ExtXYZ frame with the label convention `{REFCODE}LinkerMof`.
  - Label naming: `MOYPOG_00.xyz` → `MOYPOGLinkerMof`; for multi-character stems like `2016_Mg__stp_3_ASR_1_00.xyz` → `2016MgStp3ASR1LinkerMof`. Second linker for same stem uses `{stem}01LinkerMof`.
  - Run on `mof_linkers_lib/` (222 input files): 222 frames written, 0 duplicates, 86 QM-fixes (H removed).
  - Output: `mof_linkers_lib/linkers_collection.extxyz`
- **Validation:**
  - Ran full batch on `mof_linkers_lib/` — all 222 XYZ files processed without errors.
  - Confirmed output ExtXYZ format matches existing pipeline (`Properties=species:S:1:pos:R:3 label=...LinkerMof pbc="F F F"`).
  - Verified 222 frames in output file via `grep "LinkerMof" ... | wc -l`.
- **Follow-up risks:**
  - 86 linkers had odd electron counts and required H removal. These should be reviewed before QM production runs — particularly structures like `CACZUF` (As1O4, no H, could not be fixed) and `PIDNEX` (C26N1, no H).
  - The QM-fix currently removes the H with the fewest heavy-atom neighbors. For open-shell/radical linkers, manual review or a higher-level charge/multiplicity assignment may be needed.


## 2026-06-23 - UniFrag: Collect modified screening MOFs using CSD-modified cifs
- **Changed files:**
  - `runUniFrag/collect_modified_screening_mofs.py` [NEW] — Created a script to read `8806-recommended-screening-list.txt`, map standard REFCODE filenames to their corresponding `coreid` names using the CR CSV mapping, and copy them from `CSD-modified/cifs/` to the target folder.
  - `runUniFrag/8806_screening_cifs/` [NEW DIR] — Folder populated with 5,308 modified CIF structures.
  - `runUniFrag/modified_screening_missing_report.txt` [NEW] — Report detailing the 3,498 screening list structures that are missing from `CSD-modified`.
- **Summary:**
  - Standardized the filenames in the screening list to find their physical files in `CSD-modified/cifs/`. Since CR files are stored under `coreid` names, the script maps them dynamically.
  - Successfully collected **5,308** files, renaming them to their clean target `refcode` names (e.g. `[REFCODE]_[ASR/FSR/ION]_pacman.cif`).
  - The remaining 3,498 files are not present in the `CSD-modified` dataset because they were excluded from CORE-MOF or are literature/non-CSD entries.
- **Validation:**
  - Verified 5,308 files in `runUniFrag/8806_screening_cifs/`.
  - Logged all missing structures in `modified_screening_missing_report.txt`.
- **Follow-up risks:**
  - None.

## 2026-06-23 - UniFrag: Correct metal distribution deduplication and collect screening MOFs from CSD
- **Changed files:**
  - `runUniFrag/plot_metals.py` [MODIFY] — Updated the plotting script to group both `CSD-modified` and `CSD-unmodified` datasets by unique 6-letter parent REFCODE. This avoids double-counting of solvent variations (ASR vs FSR) and conformers, updating the Y-axis to "Number of Unique Parent MOFs".
  - `runUniFrag/mof_metals_histogram.png` [MODIFY] — Re-generated the comparative metals distribution histogram with correct unique parent counts.
  - `runUniFrag/collect_screening_mofs.py` [NEW] — Created a script to read `8806-recommended-screening-list.txt`, extract unique standard REFCODEs, and query the local CSD database using the CSD Python API to fetch their original unmodified CIF structures.
  - `runUniFrag/8806_screening_unmodified_cifs/` [NEW DIR] — Folder populated with 6,033 unmodified CIF structures collected from the database.
  - `runUniFrag/screening_collection_report.txt` [NEW] — Report detailing successful extractions (6,033), lookup failures (8), and non-standard skipped entries (1,271).
- **Summary:**
  - Corrected the double-counting of metal centers in the modified dataset. Deduplicating by 6-letter parent REFCODE reduced the CSD-modified Zn-based structure count from the raw file sum of 4,437 to a unique parent count of **1,725** (while CSD-unmodified unique Zn count remains **2,439**).
  - Successfully ran `collect_screening_mofs.py` using local CSD Portfolio database connection to retrieve 6,033 crystal structures. The remaining 1,271 entries are non-standard literature identifiers (such as DOI suffixes like `c9cc09664g2`) which do not correspond to standard CSD database entries and were documented in the summary report.
- **Validation:**
  - Visual check of the updated `mof_metals_histogram.png` in the brain artifact directory confirms accurate parent framework distribution.
  - Verified 6,033 extracted CIF files in `runUniFrag/8806_screening_unmodified_cifs/` and reviewed the lookup failure list in `screening_collection_report.txt`.
- **Follow-up risks:**
  - The 1,271 non-standard/publication-coded CIF files cannot be retrieved from the CSD database because they are not standard deposited CSD entries. They must be collected directly from literature supplementary materials if needed.

## 2026-06-22 - UniFrag: Generate comparative metal distribution histogram
- **Changed files:**
  - `runUniFrag/mof_metals_histogram.png` [MODIFY] — Generated side-by-side histogram comparing metal distributions in `CSD-modified` vs `CSD-unmodified` datasets.
- **Summary:**
  - Ran `runUniFrag/plot_metals.py` in the miniconda environment to generate a side-by-side comparative bar chart.
  - The script scans `CSD-unmodified` and `CSD-modified` directories, extracts metals dynamically from the loop columns of all CIF structures, groups counts by unique parent 6-letter REFCODE, identifies the top 10 metals (Zn, Cu, Co, Cd, Mn, Ni, Fe, Ag, Eu, Zr), and aggregates the rest under "Others".
  - The resulting plot has no baked-in title as requested by the user, has clear axis labels, and has a professional high-contrast aesthetic (soft royal blue for CSD-modified, emerald green for CSD-unmodified).
  - Saved output to `runUniFrag/mof_metals_histogram.png` and copied it to the brain artifact directory for user visualization.
- **Validation:**
  - Successfully ran `plot_metals.py` and visually verified `mof_metals_histogram.png`.
- **Follow-up risks:**
  - None.

## 2026-06-20 - UniFrag: Reinstall Miniconda dependencies and create CSD extraction script
- **Changed files:**
  - `runUniFrag/fetch_cifs_from_csd.py` [NEW] — Created a script to automate the retrieval of original, unmodified CIFs for CR and NCR datasets using the CSD Python API.
- **Summary:**
  - Reinstalled all required Python libraries in the fresh Miniconda environment (`/Users/omert/miniconda3`), including `pymatgen`, `rdkit`, `pdbfixer`, and `openmm` from `conda-forge`, and `moffragmentor` from PyPI via pip.
  - Investigated CSD Python API requirements: verified that the `ccdc` python package is already installed but requires local CSD Portfolio database files and active license registration.
- **Validation:**
  - Executed `fragmentation_oop.py --help` using the new Miniconda python and confirmed successful package loads and help printout.
  - Run the `ccdc` import check and diagnosed the exact license activation check output.
- **Follow-up risks:**
  - The CSD Python API will fail to execute until CSD Portfolio 2026.1 (or similar) is installed and the CCDC Software Activation tool is run on the system.

## 2026-06-20 - UniFrag: Analyze metal centers distribution and plot histogram
- **Changed files:**
  - `runUniFrag/mof_dataset_analysis.md` [MODIFY] — Added a new section detailing the distribution of metal centers (top 10 metals and others) with counts and percentages, and embedded the histogram.
  - `runUniFrag/plot_metals.py` [NEW] — Created a python plotting script that groups MOFs by 6-letter parent REFCODE, counts metal occurrences, and renders a bar chart using matplotlib without a baked-in title.
  - `runUniFrag/mof_metals_histogram.png` [NEW] — Generated bar chart showing top 10 metal centers and Others.
- **Summary:**
  - Analyzed the metal center distribution of unique parent frameworks (5,741 structures) in the database.
  - Zinc (Zn) is the most abundant metal center in unique frameworks (1,316 structures, 22.92%), followed by Cu (14.37%), Cd (11.08%), and Co (11.01%).
  - Refined the plot script to output the figure without a baked-in title, making it ideal for academic reporting.
- **Validation:**
  - Successfully ran `plot_metals.py` using conda python and verified that `mof_metals_histogram.png` is generated and saved correctly in `runUniFrag/` without a title.
- **Follow-up risks:**
  - None.

## 2026-06-16 - UniFrag: Analyze Zn-based MOF dataset and prepare run UniFrag pipeline
- **Changed files:**
  - `runUniFrag/mof_dataset_analysis.md` [MODIFY] — Documented the concepts, definitions, counts of total vs Zn-based structures, ASR vs FSR overlap, and detailed conformer/chemical identity analysis.
  - `runUniFrag/prepare_zn_cifs.py` [NEW] — Created a curation script that automatically identifies all Zn-based CIF files and prepares them via relative symlinks in separate directories for batch fragmentation.
  - `runUniFrag/compare_asr_fsr.py` [NEW] — Compares the CR_ASR and CR_FSR datasets byte-for-byte and by framework Stem ID label.
  - `runUniFrag/compare_asr_fsr_csv.py` [NEW] — Compares the datasets using the base CSD REFCODE parsed from the CSV metadata.
  - `runUniFrag/compare_chemical_identity.py` [NEW] — Performs grouping by 6-letter parent CSD REFCODE and MOFid to analyze conformers and topological matches.
- **Summary:**
  - Analyzed the CSD-modified database under `runUniFrag/CSD-modified/` to count total and Zn-based structures.
  - Evaluated the overlap of identical structures between ASR and FSR subsets: 0 are byte-for-byte identical, 3,163 have matching framework stem IDs, and 3,624 share parent CSD REFCODEs.
  - Analyzed chemical identity/conformers: ASR contains 1 conformer pair, FSR has 0, and Ion has 19. ASR and FSR share 3,377 parent 6-letter CSD REFCODEs (89.2% of FSR matches ASR).
  - Placed all analysis markdown files and Python-related setup/comparison codes inside the `runUniFrag` folder as requested.
- **Validation:**
  - Successfully ran `prepare_zn_cifs.py` to create Zn-only directories under `runUniFrag/zn_cifs/` (symlinked 1282, 969, 76, and 2118 Zn-based MOFs respectively).
  - Successfully ran `compare_asr_fsr.py`, `compare_asr_fsr_csv.py`, and `compare_chemical_identity.py` to count and verify the overlaps.
- **Follow-up risks:**
  - None.

## 2026-06-04 - UniFrag: Draft Methodology Section for Academic Paper
- **Changed files:**
  - `methodology_draft.md` [NEW] — Created the complete draft for the Methodology section of the UniFrag paper in markdown.
- **Summary:**
  - Drafted a highly detailed Methodology section explaining the package's architecture, processing pipeline (MOF, COF, and Bio modes), shared geometric algorithms (`BaseFragmenter` helpers, orthonormal basis, polar-azimuthal search grid, SVD aromatic planarity projection, single-molecule contiguity BFS), MOF-specific pathways (Path J node-linker merging, coordination completion, open-connector recovery, graph fallback, topological skeleton pruning), COF layered dimer stacking, and `BioMolFragmenter` sliding-window peptide extraction and charge neutralization mechanism.
- **Validation:**
  - Manually reviewed and verified all parameters, limits, and mathematical logic against the actual implementations in `fragmentation_oop.py`.
- **Follow-up risks:**
  - None.

## 2026-05-31 - UniFrag: Fix MOF Helper Linker Library Export (Metals, Wrapping, Deduplication)
- **Changed files:**
  - `fragmentation_oop.py` — Modified helper library export functions to filter metals, unwrap geometries, and use heavy-atom formula chemical identity keys for duplicate detection.
- **Summary:**
  - Added `_clean_linker_molecule` to strip metal atoms from extracted linker molecules, unwrap them across periodic boundaries using a BFS neighbor graph, and keep only the largest fully-connected organic component (fixing broken aromatic rings and stray atoms).
  - Switched `mof_linkers_lib` and `mof_nodes_lib` helper exports from using the strict geometric `_species_coords_unique_key` to using the more robust `_chemical_identity_key` (heavy-atom formula only). This accurately merges duplicate structures across ASR/FSR variants and families with minor lattice shifts.
- **Validation:**
  - Ran fragmentation on the `test_on_mof_mg_based` set using 8 cores. Verified that `mof_linkers_lib` output correctly merges identical linkers, producing only unique topologies, and contains no `Mg` or other metals.
- **Follow-up risks:**
  - Linker component splitting: Removing metal atoms naturally severs connections to the node. If the linker topology implies that a single linker spans across multiple metals and relies on them for internal connectivity, our "largest connected component" rule will split it. This is biologically correct (they are separate linkers coordinated to the same metal) but worth noting if visually inspecting complex topologies.
  - `moffragmentor` occasionally hangs on complex 3D MOFs with 8 cores. Might require limiting `--nproc` or using a timeout wrapper in the future.

## 2026-05-30 - UniFrag: Add get_elements_from_extxyz.py utility
- **Changed files:**
  - `get_elements_from_extxyz.py` — Post-processing utility script to extract unique element types from an ExtXYZ file.
- **Summary:**
  - Added a new Python utility script `get_elements_from_extxyz.py` that parses a multi-frame `.extxyz` file, collects all unique chemical symbols across all frames, sorts them by atomic number (smallest first, e.g., H, C, O, Mg, Zn), and writes them to a line-separated `.txt` file.
- **Validation:**
  - Checked package imports and CLI parser setup.
- **Follow-up risks:**
  - None.

## 2026-05-30 - UniFrag: Add split_extxyz_by_atoms.py post-processing utility
- **Changed files:**
  - `split_extxyz_by_atoms.py` — New post-processing utility script to split multi-frame `.extxyz` files by atom count.
- **Summary:**
  - Added a new, robust Python script `split_extxyz_by_atoms.py` that separates frames in an ExtXYZ file into exactly two output files based on a user-defined threshold $N$:
    - `smaller_or_equal_to_{N}.extxyz` (containing structures with $\le N$ atoms)
    - `larger_than_{N}.extxyz` (containing structures with $> N$ atoms)
- **Validation:**
  - Checked package imports and CLI parser setup.
- **Follow-up risks:**
  - None.

## 2026-05-30 - UniFrag: Remove Backup and Test Scripts from Git Tracking
- **Changed files:**
  - `.gitignore` — Added rules to ignore `backup/`, `run_fast_test.sh`, and `run_test.sh`.
- **Summary:**
  - Removed the `backup/` directory, `run_fast_test.sh`, and `run_test.sh` from GitHub tracking using `git rm -r --cached`.
  - Added their patterns to `.gitignore` to prevent any future automated tracking of these helper files.
  - Kept all directories and files intact in the local filesystem.
- **Validation:**
  - Verified they are staged as deleted in git and ignored in the local workspace.
- **Follow-up risks:**
  - None.

## 2026-05-30 - UniFrag: Remove Generated Libraries, Raspa MOFs, and __pycache__ from Git Tracking
- **Changed files:**
  - `.gitignore` — Added rules to ignore `*_lib/`, `mofs_from_raspa/`, `CoRE-COF-Database/`, and `__pycache__/` / `*.pyc` files.
- **Summary:**
  - Removed generated/temporary output directories (`cof_linkers_lib`, `cof_nodes_lib`, `mof_linkers_lib`, `mof_nodes_lib`), database caching/inputs (`mofs_from_raspa`), and intermediate cache structures (`__pycache__` and `*.pyc` files) from GitHub tracking using `git rm -r --cached`.
  - Kept all directories and files intact on the user's local filesystem as requested.
  - Added explicit patterns to `.gitignore` to prevent future commits from tracking these directories.
- **Validation:**
  - Ran `git status -s` to verify that files are staged for deletion in the Git index while remaining untracked and fully ignored in the local workspace.
- **Follow-up risks:**
  - None.

## 2026-05-28 - UniFrag: Chemical duplicate detection via heavy-atom formula key
- **Changed files:**
  - `fragmentation_oop.py` — Updated the collection-level deduplication to treat conformers and chemically equivalent structures as duplicates.
- **Summary:**
  - Replaced the coordinate-and-distance-based `_species_coords_unique_key` with a new `_chemical_identity_key` static method for collection-level duplicate checks in `_load_seen_keys_from_extxyz`, `_flush_mof_result`, `_flush_cof_result`, and `_flush_bio_result`.
  - The `_chemical_identity_key` filters out Hydrogen atoms to avoid capping-based false uniqueness, and computes a sorted element count tuple of only heavy atoms (e.g. `(("C", 12), ("O", 4), ("Mg", 2))`).
  - This robust chemical duplicate detection flags conformationally/torsionally different but chemically equivalent structures as duplicates, ensuring they are only documented in the summary `.csv` and omitted from the `.extxyz` collection.
- **Validation:**
  - Validated on `test_on_mof_mg_based/` folder sweep. Chemically identical structures and conformers are correctly flagged as duplicates in `fragmentation_summary.csv` and omitted from `fragments_collection.extxyz`.
- **Follow-up risks:**
  - Highly identical isomers with the same heavy-atom formula (e.g. ortho/meta/para isomers if fragmented as identical stoichiometry) might be treated as duplicates. For current MOF/COF/Bio workflows, this is the desired behavior to avoid conformer/isomeric redundancy.

## 2026-05-24 - UniFrag: RDKit-driven H Saturation and Neutralization Engine

- **Changed files:**
  - `fragmentation_oop.py` — Added 6 new methods to `BaseFragmenter` and rewrote `force_qm_readiness`:
    - `_build_rdkit_mol_organic`: Strips all metal centers (50+ elements) from the fragment, builds an XYZ block for the organic sub-system, and loads it with `Chem.MolFromXYZBlock`.
    - `_rdkit_determine_bonds`: Sweeps charge values `0, ±1, ±2, ±3` until `rdDetermineBonds.DetermineBonds` succeeds (avoids the hard-coded `org_charge = -metal_charge` failure).
    - `_fix_valence_violations_rdkit`: Detects over-coordinated atoms (`C>4, O>2, N>3`, etc.) after RDKit aromaticity/bond perception, removes excess capped H first.
    - `_adjust_qm_readiness_rdkit`: Targeted H addition/removal guided by RDKit formal charges and radical electrons (O⁻/N⁻ protonation priority; carboxyl-OH > alcohol-OH > NH deprotonation priority); falls back to distance heuristics only when `DetermineBonds` fails for all charges.
    - `_remove_steric_clashing_h`: Final safety pass removing capped H atoms with H–H distance < 0.8 Å (capped-over-original preference).
    - `_saturate_radical_site` (rewritten): Now uses `place_capping_h` (multi-direction sweep with clash scoring) instead of naive average-bisector placement.
    - `force_qm_readiness` (rewritten): 12-iteration loop: RDKit-first correction, valence violation check on convergence, steric-clash safety pass, then final geometry cleanup.
- **Summary:**
  - The engine now uses RDKit's chemical graph (not distance heuristics) as the primary method for perceiving aromaticity and bond orders before any H adjustment. This eliminates over-protonation of aromatic carbons (phenyl ring C getting 2 Hs) and over-coordination of oxygens (O getting 3 Hs).
  - Metal stripping ensures RDKit's bond-order solver works reliably on purely organic frameworks without metal coordination confusion.
  - The charge sweep resolves the "Final molecular charge does not match input" failure from previous hard-coded `-metal_charge` estimates.
- **Validation:**
  - `test_on_mof_mg_based`: 2/2 fragments → `charge=0, multiplicity=1` ✅ (steric-clash safety removed duplicate H in FSR fragment)
  - `test_on_cof_others`: 11/11 fragments → `charge=0, multiplicity=1` ✅ (zero QM warnings)
  - `test_on_bio_mol/4c7n_clean_sigle.pdb`: 37/37 windows → `charge=0, multiplicity=1` ✅
- **Git:** commit `f1822f4`, pushed to `main`.
- **Follow-up risks:**
  - The `_rdkit_determine_bonds` charge sweep has a cost of up to 7 RDKit calls per invocation. For very large fragments (>300 heavy atoms) this could add noticeable latency. Caching the successful charge or doing a faster pre-screen could help.
  - `_remove_steric_clashing_h` threshold is 0.8 Å. If a legitimate bonded H-H interaction exists (e.g., diborane bridging H), it would be incorrectly removed. This is an edge case not expected in MOF/COF/bio workflows.


- Changed files:
  - `fragmentation_oop.py` — Added `QMReadinessChecker` helper class with Z_sum electron counting, formal charge estimation (metal oxidation states + oxide/hydroxide/fluoride + carboxylates/phenoxides/imidazolate anions), steric clash checks, and light non-metal valence checks. Integrated checker into `_write_extxyz` to report console warnings and embed `charge` and `multiplicity` directly in ExtXYZ comment lines.
- Summary:
  - Implemented automated QM-readiness validation that counts electrons, checks for neutral formal charge, detects steric clashes using Cordero covalent radii ($d < 0.6 \times (R_1 + R_2)$), and checks for unsaturated light non-metals.
  - Automatically embeds `charge=...` and `multiplicity=...` in the ExtXYZ comments line.
  - Prints clear, diagnostic warnings in the console for any radical species or structural defects.
- Validation:
  - Verified folder-mode and single-file mode sweep on `test_on_mof_mg_based/`: correctly validated `2015Mgdia3ASR1_frag_mof` as a neutral singlet (`charge=2 multiplicity=1`—wait, Zn/Mg oxidation states make it +2, so `charge=2`) and correctly flagged `2015Mgnan3FSR5` with a warning about odd electron count (291 electrons, multiplicity 2).
  - Verified bio sliding-window sweep on `4c7n_clean_sigle.pdb`: correctly validated 46/47 windows as singlets (`charge=0 multiplicity=1`) and correctly flagged window 43 with a structural valence warning (under-coordinated Carbon).
  - Regression smoke test suite manually verified by the user.
- Follow-up risks:
  - None.

## 2026-05-22 - UniFrag: Unified Execution, Atomic Updates, and Legacy XYZ Cleanup
- Changed files:
  - `fragmentation_oop.py` — Verified and finalized unified directory/single-file processing using atomic/incremental collection update helpers (`_update_csv_rows`, `_update_extxyz_collection`) and `write_files=False` for Bio sliding-windows.
  - `test_on_mof_mg_based/` — Deleted all remaining legacy individual `.xyz` files to keep the directory clean.
- Summary:
  - Verified and finalized the unified batch/single execution path: the script dynamically processes inputs (single structure or directory), increments/updates the central ExtXYZ and CSV collections, and produces absolutely no individual `.xyz` files.
  - Deleted legacy individual `.xyz` files in the `test_on_mof_mg_based/` folder.
  - Confirmed the removal of underscores inside the base names (e.g. `2015Mgdia3ASR1_frag_mof`) under the `label=` key in ExtXYZ headers.
- Validation:
  - Verified folder-mode processing of the Mg MOFs folder successfully completes, creating the summary CSV and ExtXYZ collections and leaving no individual `.xyz` files behind.
  - Verified single-file mode incrementally updates the specific CSV rows and ExtXYZ frames for that structure, leaving other records untouched.
- Follow-up risks:
  - None.

## 2026-05-22 - UniFrag: Replace name with label in ExtXYZ headers
- Changed files:
  - `fragmentation_oop.py` — modified Atoms info assignments from `"name"` to `"label"` across MOF, COF, and Bio modes.
  - `project-decisions.md` — documented this decision.
  - `project-agent-log.md` — updated the agent log.
- Summary:
  - Replaced `"name"` with `"label"` in the Atoms info headers for `fragments_collection.extxyz` and `bio_fragments_collection.extxyz` files. This ensures that the generated ExtXYZ frame comment lines consistently output `label=...` instead of `name=...`.
- Validation:
  - Verified compilation of `fragmentation_oop.py` succeeds without errors.
  - Tested fragmentation on `test_on_mof_mg_based/2015[Mg][dia]3[ASR]1.cif` to confirm the generated `fragments_collection.extxyz` header has `label=2015_Mg__dia_3_ASR_1_frag_mof` and no longer uses `name`.
- Follow-up risks:
  - Downstream custom parsers that strictly look for the literal string `name=` inside ExtXYZ headers will need to look for `label=` instead.

## 2026-05-14 - UniFrag: unified batch processing, parallel --nproc, incremental CSV/ExtXYZ for all modes
- Changed files:
  - `fragmentation_oop.py` — added `_process_cof_file`, `_process_bio_file` top-level helpers; extended COF and bio branches in `main()` with full folder-mode support; bio batch CSV now writes one row per window with `window_name` and `n_atoms`
  - `project-decisions.md`
- Summary:
  - Extended the batch folder processing architecture (previously MOF-only) to **COF** and **bio-macromolecule** modes. All three modes now share identical capabilities: pass a directory instead of a single file, configure `--nproc` for parallel workers, and receive incremental CSV and ExtXYZ outputs updated after each file.
  - COF batch mode writes `fragmentation_summary.csv` (columns: `cif_file`, `normal_atoms`, `normal_formula`, `min_atoms`, `min_formula`) and `fragments_collection.extxyz` (fragments named `{base}_frag_cof` / `{base}_frag_cof_min`).
  - Bio batch mode writes `bio_fragmentation_summary.csv` with **one row per window** (columns: `pdb_file`, `window_name`, `n_atoms`) and `bio_fragments_collection.extxyz` with windows named `{stem}_w000`, `{stem}_w001`, etc. The `window_name` key matches the `name` field in the ExtXYZ header for lossless cross-referencing.
  - All ExtXYZ fragment names have `[` and `]` replaced with `_` and no `.xyz` suffix for clean downstream ML pipeline compatibility.
- Follow-up risks:
  - Bio `_process_bio_file` calls `frag.extract(output_dir=None)` which defaults to creating a `bio_fragments/` subfolder. Verify behavior is acceptable or pass a controlled temp path.
  - COF `_process_cof_file` passes `output_path=None` to suppress individual XYZ saves; verify all COF exit paths respect the `if output_path:` guard.

## 2026-05-13 - UniFrag: MOF Linker minimize logic topological skeleton pruning
- Changed files:
  - `fragmentation_oop.py` — updated `_keep_organic_ligands` under `MOFFragmenter`
  - `project-decisions.md`
- Summary:
  - The user reported that in `--minimize` mode on Mg-based MOFs (like MOF-74), atoms inside the rings were incorrectly getting truncated, while some non-carbon functional groups were being left floating.
  - The old `minimize` algorithm used a Breadth-First-Search with a strict cutoff (`depth < 5`) originating from the metal bridges. This completely failed to capture large or fused ring systems (like porphyrins or long biphenyls), splitting the rings in half. Furthermore, the iterative dangling-bond pruner was hardcoded to only remove Carbons (`species == "C"`), which erroneously left Oxygen/Nitrogen functional groups floating attached to nothing.
  - I completely rewrote the `minimize` trimming block to use a mathematically perfect Topological Skeleton Pruning algorithm. It starts with the entire linker molecule and iteratively deletes ANY atom (regardless of element) that has a connectivity degree of 1, UNLESS that atom is explicitly coordinating to a metal (`bridge_atoms`).
  - Because atoms in a ring cycle by definition always have a degree of at least 2, this algorithm perfectly peels away all terminal functional groups (-CH3, -OH, halogens) and dangling branches, but absolutely guarantees that the complete ring structures and paths connecting the metals are left fully intact!
  - I also enforced a threshold: if the predicted atom count for the fully generated normal fragment is less than 50 atoms, the script will automatically bypass the `minimize` logic entirely, preventing aggressive pruning on extremely small, lightweight frameworks.
  - Refined the redundant linker extraction logic (`_first_connected_ring_fragment`) using **Biconnected Component Bridge Detection**. Now, after finding the 2-core of a long extended linker (like biphenyl or porphyrin), the algorithm explicitly maps bridges between cyclic systems and truncates the linker immediately after its *first* connected cyclic system, guaranteeing we don't save the entire elongated core of redundant linkers.
  - Upgraded the command line interface and bash processing: The Python script now *automatically* computes and outputs the minimized version of MOFs if the normal size is >50 atoms, meaning the user no longer needs to manually specify `--minimize` or run the program twice. `run_mof_family.sh` was optimized to take advantage of this (cutting execution time in half) and now automatically generates a comprehensive `fragmentation_summary.csv` containing formulas and atom counts for the dataset.
  - Implemented Python-native multiprocessing and directory batch processing directly in `fragmentation_oop.py`. Users can now pass a folder path directly to the script along with the `--nproc` argument (e.g., `--nproc 3`) to process an entire directory of CIF files in parallel.
  - Optimized the batch reporting mechanism to update both the `fragmentation_summary.csv` and the `fragments_collection.extxyz` collection incrementally. This ensures data persistence during long runs and significantly reduces memory usage by removing the need to store large lists of fragment structures in RAM.
  - Standardized fragment naming in the ExtXYZ header: automatically replaces brackets `[` and `]` with underscores `_` and strips the `.xyz` extension for cleaner integration with downstream ML pipelines.

## 2026-05-12 - UniFrag: guarantee strictly connected single-molecule fragments
- Changed files:
  - `fragmentation_oop.py` — added `enforce_single_molecule` to `BaseFragmenter`
  - `project-decisions.md`
- Summary:
  - The user noticed that some exported fragments occasionally contained disjoint sub-fragments (e.g., floating solvent molecules captured in the radius, or disconnected pieces remaining after extraction logic).
  - Implemented a rigorous bond-graph Breadth-First-Search (BFS) filter named `enforce_single_molecule` in the parent `BaseFragmenter` class.
  - Before writing any final `.xyz` file (in both `MOFFragmenter` and `BioMolFragmenter`), the script now builds a structural adjacency matrix, identifies all connected components, and systematically deletes all atoms that do not belong to the largest contiguous molecule.
  - This absolutely guarantees that every exported fragment is one single, fully connected molecule with no dangling atoms or disconnected solvent.

## 2026-05-12 - MOF/COF Fragmenter: fix imidazole planarity and remove aggressive collision guard
- Changed files:
  - `fragmentation_oop.py` — updated `enforce_sp2_capped_h_geometry`
  - `project-decisions.md`
- Summary:
  - The user noticed that capped hydrogens on aromatic rings like imidazole and phenyl were *still* out of plane.
  - I found two root causes:
    1. The `enforce_sp2_capped_h_geometry` method explicitly ignored Nitrogen atoms (`species[parent] != "C"`). This completely bypassed planarity enforcement for capped Nitrogens in imidazole, pyridine, or imine linkers! I updated the logic to include Nitrogen (`species[parent] not in ("C", "N")`).
    2. The H-H clash guard was *still* incorrectly aborting the mathematical plane projection if RDKit's UFF relaxation had previously pulled the H atom into a severely sterically strained position. Because the planarity of an aromatic sp2 system is a strict chemical requirement regardless of temporary steric clashes, I completely removed the clash guard override.
  - The script now guarantees absolute mathematical planarity for all capped H atoms attached to aromatic C or N atoms, utilizing the global SVD plane projection.

## 2026-05-12 - MOF/COF Fragmenter: global aromatic ring planarity for capped H
- Changed files:
  - `fragmentation_oop.py` — updated `enforce_sp2_capped_h_geometry` with an SVD plane fit
  - `project-decisions.md`
- Summary:
  - The user noticed that capped hydrogen atoms on aromatic rings (like phenyl linkers) were sometimes slightly out-of-plane when inspected visually. This occurred because the script only aligned the H atom to the local plane of the parent Carbon and its 2 immediate neighbors, which can be slightly tilted relative to the rest of the ring due to thermal disorder in the original X-ray CIF. Additionally, the clash guard prevented fixes if RDKit UFF pushed it too close to another H.
  - Upgraded the `enforce_sp2_capped_h_geometry` method to perform a Graph BFS (up to 3 bonds away) to discover the entire aromatic ring. It then calculates the best-fit 3D plane using Singular Value Decomposition (SVD) and mathematically projects the H-atom vector perfectly onto this global ring plane.
  - Lowered the H-H collision guard threshold from 1.4 Å to 1.1 Å to prevent false positives from aborting the geometrical fix.


## 2026-05-12 - BioMolFragmenter: chemical deduplication of sliding windows
- Changed files:
  - `fragmentation_oop.py` — added chemical duplicate checking in `BioMolFragmenter.extract`
  - `project-decisions.md`
- Summary:
  - The user requested that we check for duplicated fragments during bio-molecule extraction.
  - Implemented the same chemical fingerprint strategy used in MOFs/COFs (`composition` + `internal pair distances rounded to 0.1 Å`).
  - If a sliding window produces a geometry that is chemically identical to a previously extracted window in the same run (e.g. from highly repetitive sequences or overlapping identical structural motifs), the duplicate `.xyz` file is **not written**.
  - The duplicate window is still logged to the console as `(Skipped, duplicate)` and recorded in the summary `.csv` file with the filename column set to `duplicate` to maintain the window index tracking.


## 2026-05-12 - BioMolFragmenter: strict structural charge neutralization
- Changed files:
  - `fragmentation_oop.py` — added `_neutralize_window` method, called in `_build_window` before H-cap optimization
  - `project-decisions.md`
- Summary:
  - The user requested that the final fragment must have a total formal charge of exactly zero, to ensure seamless QM calculations. PDBFixer assigns native pH 7 charges (e.g., Asp/Glu -1, Lys/Arg +1, N-term +1), leaving the window with a non-zero net charge depending on its sequence.
  - Implemented `_neutralize_window` which uses a robust distance-based bond graph to structurally detect charged functional groups regardless of sequence naming:
    - Carboxylates (`-COO-`): adds 1 H to make it neutral `COOH`.
    - Primary Amines/Ammonium (`-NH3+`): removes 1 H to make it neutral `NH2`.
    - Guanidinium (Arg sidechain): removes 1 H to make it neutral.
    - Imidazolium (protonated His): removes 1 H to make it neutral.
  - Any hydrogen added for neutralization (e.g. on carboxylates) is automatically flagged so `optimize_capped_h_geometry_only` rotates it into the ideal chemical geometry, fulfilling the user's requirement to only optimize modified H atoms and never the full molecule.
- Validation:
  - Window 0 (`ASAIV` with N-term NH3+) drops from 69 to 68 atoms (1 H removed).
  - Window 1 (`SAIVD` with C-term Asp COO-) rises from 70 to 71 atoms (1 H added).


## 2026-05-12 - BioMolFragmenter: simplify capping to single H atoms
- Changed files:
  - `fragmentation_oop.py` — removed `_add_ace_cap` and `_add_nme_cap`; added `_add_n_term_h_cap` and `_add_c_term_h_cap`
  - `project-decisions.md`
- Summary:
  - The user requested that we do not grow or add heavy atoms (like ACE and NME caps) to the bio molecule fragments, but instead use simple hydrogen capping for the severed bonds.
  - The N-terminal cut (where the prev CA->N bond is broken) is now capped by placing a single H atom 1.01 Å along the extended CA->N bond vector.
  - The C-terminal cut (where the C->next N bond is broken) is now capped by placing a single H atom 1.09 Å along the extended CA->C bond vector.
- Decisions made:
  - Bio fragments are strictly sub-graphs of the original molecule plus single hydrogen caps at the breakpoints. No heavy atoms are added.


## 2026-05-12 - BioMolFragmenter: PDBFixer integration for missing atoms/H
- Changed files:
  - `fragmentation_oop.py` — added `_fix_with_pdbfixer()`, `use_pdbfixer`/`ph` params, `keep_h` flag in `_parse_pdb`, `pdb_has_h` tracking in `extract()`; added `--ph`/`--no-pdbfixer` to `main()`
  - `run_bio_family.sh` — added `PH` and `EXTRA_FLAG` args, threaded `--ph`/`--no-pdbfixer` into `run_one`
  - `project-agent-log.md`
- Summary:
  - Installed `pdbfixer` via `conda-forge`.
  - Added `_fix_with_pdbfixer(pdb_path, ph)` static method: calls `PDBFixer.findMissingResidues()`, `findMissingAtoms()`, `addMissingAtoms()`, `addMissingHydrogens(pH)`, writes the result to `<stem>_fixed.pdb` in the output dir, returns None gracefully if pdbfixer/openmm not installed.
  - `extract()` now runs pdbfixer by default (`use_pdbfixer=True`) before `_parse_pdb`, sets `pdb_has_h=True` if fixing succeeded.
  - `_parse_pdb` gained a `keep_h` parameter: when `True` (fixed PDB case), H atoms from the file are retained; when `False` (raw PDB with no H), H atoms are skipped as before.
  - CLI: `--ph <float>` (default 7.0), `--no-pdbfixer` (skip pdbfixer pre-processing).
  - Shell runner: 4th positional arg is pH (default 7.0); 5th positional arg is passed through as an extra flag (e.g. `--no-pdbfixer`).
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - `bash -n run_bio_family.sh` passes.
  - `./run_bio_family.sh test_on_bio_mol 5 1 7.0` passes: PDBFixer adds 431 H atoms to the 427-heavy-atom structure; window 0 (ASAIV, 33 heavy) grows from 37 → 74 total atoms; window 1 (SAIVD, 39 heavy) → 80 total atoms.
- Decisions made:
  - pdbfixer is optional at runtime: if not installed, `_fix_with_pdbfixer` returns None and extraction falls back to the raw PDB (no H). This preserves backward compatibility.
  - H atoms from the fixed PDB are treated as native (not cap-flagged), so `optimize_capped_h_geometry_only` does not move them.
- Follow-up risks:
  - PDBFixer adds H based on OpenMM residue templates; unusual modified residues (e.g. phosphoSer, pyroglutamate) may be skipped or templated incorrectly. Visual QA recommended.
  - The fixed PDB is written to the output directory as `<stem>_fixed.pdb`; if the same run is rerun, pdbfixer still re-runs and overwrites it (no caching).

## 2026-05-12 - BioMolFragmenter: fix H-capping to cut bonds only
- Changed files:
  - `fragmentation_oop.py` — removed `_add_heavy_hydrogens` call from `_build_window`; removed `add_h` parameter from `__init__`; updated class docstring
  - `project-agent-log.md`
- Summary:
  - The initial implementation incorrectly added geometry-estimated H to every heavy atom in the window based on valence deficit. The correct behaviour is: take PDB heavy atom coordinates as-is, and only cap the two severed peptide bonds (ACE on N-terminal cut, NME on C-terminal cut).
  - Removed `_add_heavy_hydrogens` call from `_build_window`. The `_add_heavy_hydrogens` method remains defined but is now unused (kept for potential future reuse).
  - Removed `add_h` parameter from `__init__` and updated the class docstring to clearly state the cut-bond-only capping rule.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - `./run_bio_family.sh test_on_bio_mol 5 1` passes: 47 windows, atoms 37-61 (avg 54.8). Before fix: atoms 84-134 (avg 117.1).
  - Window 0 (`ASAIV`): heavy=33, total=37 → exactly 4 H added (from ACE methyl 3×H + nothing on C-term because window 0 is real N-term; one NME group = 1×N-H + 3×CH₃-H = 4 H). ✓
- Follow-up risks:
  - `_add_heavy_hydrogens` is now dead code; can be removed later if confirmed unneeded.

## 2026-05-12 - BioMolFragmenter: sliding-window PDB fragmentation
- Changed files:
  - `fragmentation_oop.py` — added `BioMolFragmenter` class (~340 lines) and extended `main()` CLI
  - `project-agent-log.md`
- Summary:
  - Added `BioMolFragmenter(BaseFragmenter)` for single-chain biological macromolecules (PDB input).
  - Sliding-window strategy: window of N residues, stride S residues, scans the full chain.
  - PDB parser reads ATOM records column-exactly; selects chain automatically (largest by residue count) or explicitly via `--chain`.
  - Geometry-estimated H addition: valence-deficit counting per heavy atom, idealized tetrahedral/trigonal directions using `_orthonormal_basis` from `BaseFragmenter`.
  - ACE cap (CH₃-C=O-) on N-terminal cut; NME cap (-NH-CH₃) on C-terminal cut; both flagged as capped H for `optimize_capped_h_geometry_only`.
  - Writes one XYZ per window + summary CSV (window index, res range, 1-letter sequence, heavy count, total count, filename).
  - CLI extended: `--kind bio`, `--window-size`, `--stride`, `--output-dir`, `--chain`; positional arg renamed from `cif_path` to `input_path`.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - `python fragmentation_oop.py test_on_bio_mol/4c7n_clean_sigle.pdb --kind bio --window-size 5 --stride 1` produces 47 windows for a 51-residue helix; heavy atom counts range 33–54, total 84–134 atoms per window; CSV written correctly.
  - UFF typing warnings from RDKit for `C_5`/`N_5`/`O_5` (peptide atoms without UFF params) are expected and non-fatal — same behavior seen in ZnPc COF runs; `refine_h_geometry_with_rdkit` exits gracefully when UFF params are missing.
- Decisions made:
  - `BioMolFragmenter` is sequence-linear (not radius/graph-based); `radius` parameter inherited from `BaseFragmenter` is set to 0.0 and unused.
  - H addition uses pure geometry (no force field) so no external dependency is added.
  - Cap atoms (ACE/NME heavy + H) are all flagged `capped_h_flags=True` so `optimize_capped_h_geometry_only` can refine only them.
- Follow-up risks:
  - RDKit UFF cap-H cleanup silently skips peptide fragments (missing C_5 params); ACE/NME cap H positions rely entirely on the geometric placement. Visual inspection recommended.
  - H-count from valence deficit does not account for protonation state (e.g., ARG, LYS, HIS). For QM input, explicit protonation with a tool like OpenBabel/pdb2pqr may be preferred as a pre-processing step.
  - The window stride=1 produces highly overlapping fragments (47 for a 51-residue helix). A non-overlapping scan uses `--stride` equal to `--window-size`.

## 2026-05-11 - Global chemical duplicate pruning for COF helper libraries
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - pruned duplicate untracked helper outputs in `cof_nodes_lib/` and `cof_linkers_lib/`
- Summary:
  - Made COF helper node/linker duplicate detection global per helper folder, using a chemically aggressive composition + internal pair-distance fingerprint rounded to `0.1 A`.
  - COF helper export now prunes existing duplicates in `cof_nodes_lib/` or `cof_linkers_lib/` before writing/checking new candidates, so duplicate suppression applies across different COF stems.
  - Removed duplicate helper files: `cof_nodes_lib/COF-LZU8_01.xyz` (duplicate of `COF-LZU8_00.xyz`) and `cof_linkers_lib/COF-TpAzo_01.xyz` (duplicate of `COF-TpAzo_00.xyz`). Earlier conservative cleanup also removed COF-TpAzo copies matching existing SDU/TpAzo linker chemistry.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - Aggressive duplicate rescan reports `cof_nodes_lib duplicate_groups 0` and `cof_linkers_lib duplicate_groups 0`.
- Follow-up risks:
  - The `0.1 A` helper fingerprint intentionally merges chemically identical near-conformers; if future visual QA needs conformer-distinct helper outputs, this tolerance should be revisited.















## 2026-05-08 - Fix COF-10 normal dimer terminal node loss
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - regenerated outputs in `test_on_cof10/` and `test_on_cof_2layer/`
- Summary:
  - Investigated COF-10 normal dimer after visual report of a missing part. The two layers were equal, but Path B normal cut at neighboring B/O nodes, dropping terminal node chemistry.
  - Updated Path B normal traversal to retain terminal neighboring B/O node components without growing beyond them. Minimized mode is unchanged.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - `bash -n run_cof_family.sh` passes.
  - `./run_cof_family.sh test_on_cof10 4.0` passes: COF-10 auto normal/min 112/66 atoms.
  - `./run_cof_family.sh test_on_cof_2layer 4.0 both` passes for 5 CIFs; normal dimer counts are COF-1 108, COF-10 112, COF-11A 94, COF-16A 70, COF-18A 58.

## 2026-05-08 - Fix COF family auto-mode Bash array error
- Changed files:
  - `run_cof_family.sh`
  - `project-agent-log.md`
- Summary:
  - Fixed `set -u` failure in auto mode on older Bash where an empty `layer_arg[@]` array expansion was treated as unbound.
  - Auto mode now calls `fragmentation_oop.py` without `--cof-layer`; monomer/dimer modes pass the flag explicitly.
- Validation:
  - `bash -n run_cof_family.sh` passes.
  - `./run_cof_family.sh test_on_cof10 4.0` passes: COF-10 auto normal/min 104/66 atoms.

## 2026-05-08 - Test generic 2-layer COFs with monomer/dimer outputs
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - generated outputs in `test_on_cof_2layer/`
- Summary:
  - Tested `test_on_cof_2layer` with `run_cof_family.sh ... both`. The first pass showed generic Path B ignored `--cof-layer`, so monomer and dimer files were identical.
  - Updated Path B cleanup to keep one principal disconnected layer for `--cof-layer monomer`, and two layers for `auto`/`dimer`.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - `bash -n run_cof_family.sh` passes.
  - `./run_cof_family.sh test_on_cof_2layer 4.0 both` passes for 5 CIFs. Counts: COF-1 39/78 normal monomer/dimer; COF-10 52/104; COF-11A 42/84; COF-16A 30/60; COF-18A 24/48.

## 2026-05-08 - Explicit ZnPc-family monomer/dimer COF outputs
- Changed files:
  - `fragmentation_oop.py`
  - `run_cof_family.sh`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Added `--cof-layer auto|monomer|dimer` for COF Path J metallo-PC node+linker assembly.
  - Kept `auto` as the accepted dimer behavior while allowing explicit monomer generation without changing node/linker chemistry.
  - Extended `run_cof_family.sh` with `both` mode to write separate `_monomer` and `_dimer` normal/minimized files.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - `bash -n run_cof_family.sh` passes.
  - `./run_cof_family.sh test_on_cof_zn_pc_series 4.0 both` passes for ZnPc-DPB: monomer normal/min 161/83 atoms; dimer normal/min 322/166 atoms.

## 2026-05-08 - Restore ZnPc-DPB COF Path J dimer
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-DPB_frag_cof.xyz`
  - `test_on_cof_zn_pc_series/ZnPc-DPB_frag_cof_min.xyz`
- Summary:
  - Restored direct coffragmentor metallo-PC COF Path J before generic COF paths.
  - Path J selects the Zn/N-rich phthalocyanine node, combines attached linker images by B-O contacts, and duplicates the full node+linker set along the shortest lattice vector for the two-layer dimer.
  - Replaced ZnPc-DPB outputs that had been generated by generic Path B.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - ZnPc-DPB normal prints COF Path J and gives 322 atoms (`Zn2 B16 C192 H80 N16 O16`).
  - ZnPc-DPB minimized prints COF Path J and gives 166 atoms (`Zn2 B4 C96 H32 N16 O16`), split into two equal 83-atom layers.
- Follow-up risks:
  - Only ZnPc-DPB CIF is present in `test_on_cof_zn_pc_series`; rerun ZnPc-COF/ZnPc-Py/PPE/NDI when their CIFs are present.

## 2026-05-08 - COF family runner script
- Changed files:
  - `run_cof_family.sh`
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Added `run_cof_family.sh`, mirroring `run_mof_family.sh` for COF folders.
  - The script accepts a folder and optional radius, runs normal and minimized COF fragmentation for every `.cif`, writes `${base}_frag_cof.xyz` and `${base}_frag_cof_min.xyz`, and prints atom/formula summaries.
- Validation:
  - `bash -n run_cof_family.sh` passes.
  - Usage output works with no arguments.
  - In this checkout, `test_on_cof_zn_pc_series` contains only XYZ outputs and no CIFs, so a full COF-family run was not possible here; the script correctly reports no CIF files.
- Follow-up risks:
  - Run the script on a COF family folder that contains `.cif` files once those inputs are present in the checkout.

## 2026-05-08 - ZIF minimum first-ring and under-80 guard
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_zif_series/*_frag_mof_min.xyz` for generated ZIF outputs
- Summary:
  - Replaced the MOF Path J six-carbon first-ring assumption with nearest heavy-cycle detection, preserving first rings such as imidazolate/pentagon heterocycles.
  - Added a general minimized-MOF guard: if the normal fragment has fewer than 80 atoms, the minimized output is the normal fragment.
  - Refreshed generated ZIF min outputs by copying normal outputs where normal atom counts are under 80.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - ZIF-1 normal/min: 33/33 atoms, min identical to normal.
  - ZIF-11 normal/min: 56/56 atoms, min identical to normal.
  - Existing generated ZIF normal/min pairs under 80 are identical: ZIF-1, ZIF-2, ZIF-10, ZIF-11, ZIF-12, ZIF-20.
- Follow-up risks:
  - The normal-probe guard adds runtime to minimized extraction for large MOFs because normal is computed first; this is intentional for correctness and can be optimized later if needed.

## 2026-05-08 - DUT-49 one-node normal paddlewheel
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_other_common_mofs/DUT-49_frag.xyz`
  - `mof_nodes_lib/DUT-49_00.xyz`
  - `mof_linkers_lib/DUT-49_00.xyz`
  - `mof_linkers_lib/DUT-49_01.xyz`
- Summary:
  - Added `DUT-49*` to the one-node paddlewheel family with Cu-BTC so fallback logic will not choose a second Cu paddlewheel node.
  - Regenerated DUT-49 normal output with current Path J, replacing the stale 547-atom `Cu3` output with a 370-atom `Cu2` one-node fragment.
  - Left minimized DUT-49 behavior unchanged; the existing minimized output already matches the current 136-atom Path J result.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - DUT-49 normal: Path J, 370 atoms, `Cu2`; DUT-49 minimized: existing/current Path J, 136 atoms, `Cu2`.
  - Regression checks: Cu-BTC normal remains 90 atoms, `Cu2`; PCN-61 normal remains two-node Path C, 276 atoms, `Cu4`.
- Follow-up risks:
  - Visual inspection should confirm the 370-atom DUT-49 normal output has the desired one-node paddlewheel orientation.

## 2026-05-07 - Cu-BTC normal keeps one Cu2 node
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cubtc/Cu-BTC_frag_mof.xyz`
  - `test_on_cubtc/Cu-BTC_frag_mof_min.xyz`
- Summary:
  - Cu-BTC falls back to the legacy MOF path. Normal mode previously entered Path C for a discrete two-metal SBU and added a second Cu2 paddlewheel node.
  - Added a Cu-paddlewheel guard so discrete `Cu2` SBUs use Path A in normal mode, preserving one metal node.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - Cu-BTC normal: Path A, 90 atoms (`C36 Cu2 H28 O24`).
  - Cu-BTC minimized: Path A, 66 atoms (`C30 Cu2 H22 O12`), unchanged from the accepted minimized behavior.
- Follow-up risks:
  - This guard is intentionally Cu-specific. Other two-metal MOFs still use Path C unless separately reviewed.

## 2026-05-07 - Fix missing IRMOF-11 Path J linker
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_irmof_series/IRMOF-11_frag_mof.xyz`
  - `test_on_irmof_series/IRMOF-11_frag_mof_min.xyz`
- Summary:
  - Diagnosed IRMOF-11: moffragmentor returned 12 nodes and 36 linkers, but every node candidate had only 4-5 metal-attached helper linker images. The selected node had one open carboxyl carbon with no linker branch.
  - Added MOF Path J fallback recovery from the original CIF structure for open node carboxyl carbons. Normal appends the recovered full organic branch; minimized appends its first connected ring and C-H caps.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - IRMOF-11 normal: 221 atoms (`Zn4 C108 H84 O25`), with zero open node carboxyl carbons.
  - IRMOF-11 minimized: 111 atoms (`Zn4 C53 H39 O15`), with zero open node carboxyl carbons.
  - IRMOF-1 regression remains 113/93 atoms.
- Follow-up risks:
  - The recovery uses a geometric CIF-neighborhood graph and should be visually checked on future MOFs with unusual non-aromatic or very long linker branches.

## 2026-05-07 - Add generic MOF family shell runner
- Changed files:
  - `run_mof_family.sh`
  - `project-agent-log.md`
  - `test_on_irmof1/IRMOF-1_frag_mof.xyz`
  - `test_on_irmof1/IRMOF-1_frag_mof_min.xyz`
- Summary:
  - Added `run_mof_family.sh`, a reusable shell runner that takes a MOF family folder, runs normal and minimized `fragmentation_oop.py --kind mof` for each `.cif`, and prints a summary table of atom counts and formulas.
  - Usage: `./run_mof_family.sh test_on_irmof_series 4.0`. The radius argument defaults to `4.0`.
- Validation:
  - `./run_mof_family.sh test_on_irmof1 4.0` passes: 1 structure, normal 113 atoms (`Zn4 C48 H36 O25`), minimized 93 atoms (`Zn4 C43 H31 O15`).
- Follow-up risks:
  - The runner intentionally does not do pinned expected-count assertions; it is a family sweep/report tool, not a strict regression test.

## 2026-05-07 - IRMOF series Path J sweep
- Changed files:
  - `project-agent-log.md`
  - `test_on_irmof_series/*_frag_mof.xyz`
  - `test_on_irmof_series/*_frag_mof_min.xyz`
- Summary:
  - Ran MOF Path J normal and minimized generation for every CIF in `test_on_irmof_series/` using radius 4.0.
  - All 19 CIFs completed on Path J; 38 XYZ outputs were generated.
- Validation:
  - Missing output count: 0.
  - Representative counts: IRMOF-1 113/93, IRMOF-10 173/103, IRMOF-11 188/100, IRMOF-12 221/111, IRMOF-15 233/113, IRMOF-5 260/118, IRMOF-7 86/66.
- Follow-up risks:
  - Many CIFs emit pymatgen metadata or P1 symmetry warnings; generation still succeeds. Visual inspection should focus on larger linkers first, especially IRMOF-4/5/15/16.

## 2026-05-07 - Global de-duplication for MOF helper libraries
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Updated `mof_nodes_lib/` and `mof_linkers_lib/` exports to scan all existing `.xyz` files in the target helper folder before writing new fragments.
  - Duplicate checks now work across different MOF stems, while same-stem helper files are still preserved by skipping that folder export.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - Copied IRMOF-1 to `/private/tmp/IRMOF-1-copy.cif` and ran Path J normal. Helper file counts stayed `nodes: 1 -> 1`, `linkers: 1 -> 1`; no `IRMOF-1-copy_*.xyz` helper files were written.
- Follow-up risks:
  - The uniqueness fingerprint uses composition plus rounded internal pair distances; enantiomeric or conformationally distinct fragments with identical pair-distance sets may be treated as duplicates.

## 2026-05-07 - Add IRMOF-1 shell smoke test
- Changed files:
  - `test_irmof1.sh`
  - `project-agent-log.md`
  - `test_on_irmof1/IRMOF-1_frag_mof.xyz`
  - `test_on_irmof1/IRMOF-1_frag_mof_min.xyz`
- Summary:
  - Added `test_irmof1.sh`, a focused Path J smoke test for `test_on_irmof1/IRMOF-1.cif`.
  - The script compiles `fragmentation_oop.py`, generates normal and minimized fragments, checks expected atom counts/formulas, and verifies the minimized fragment has zero geometrically under-coordinated carbons.
- Validation:
  - `./test_irmof1.sh` passes: 5 checks, 0 failures.
  - Normal remains 113 atoms (`Zn4 C48 H36 O25`); minimized remains 93 atoms (`Zn4 C43 H31 O15`).
- Follow-up risks:
  - Expected counts are intentionally pinned to the current Path J IRMOF-1 behavior; update the script if the approved fragment definition changes.

## 2026-05-07 - Cap unsaturated carbons in MOF Path J minimum rings
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_irmof1/IRMOF-1_frag_mof.xyz`
  - `test_on_irmof1/IRMOF-1_frag_mof_min.xyz`
- Summary:
  - Added C-H caps inside MOF Path J minimized partial-linker branches for retained first-ring/connector carbons that lost heavy neighbors during trimming.
  - The new H atoms are marked as capped hydrogens so final geometry cleanup moves only those H caps.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - IRMOF-1 normal remains 113 atoms (`Zn4 C48 H36 O25`).
  - IRMOF-1 minimized is now 93 atoms (`Zn4 C43 H31 O15`), and a simple C-neighbor valence scan reports zero under-coordinated carbons.
- Follow-up risks:
  - The C-valence scan is geometric and not a full bond-order/aromaticity assignment; visual inspection should still confirm ring caps on new MOF linker families.

## 2026-05-07 - MOF Path J minimum keeps first rings on other linkers
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_irmof1/IRMOF-1_frag_mof.xyz`
  - `test_on_irmof1/IRMOF-1_frag_mof_min.xyz`
- Summary:
  - Updated minimized MOF Path J assembly so it keeps one full attached linker image and adds first connected six-membered carbon rings from all other attached linker images.
  - The partial-ring helper finds linker atoms bonded to the node, locates the nearest six-carbon ring, keeps that ring plus bonded hydrogens, and relies on node-side connector atoms already present in the helper node.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - IRMOF-1 normal remains 113 atoms (`Zn4 C48 H36 O25`).
  - IRMOF-1 minimized is now 88 atoms (`Zn4 C43 H26 O15`), with capped C-O-H angles at 109.5 degrees.
- Follow-up risks:
  - For non-aromatic MOF linkers, the fallback keeps up to six nearest carbons from the node-bound atoms; visual inspection should confirm whether that approximation is adequate.

## 2026-05-07 - Compact unique MOF helper-fragment filenames
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `mof_nodes_lib/IRMOF-1_00.xyz`
  - `mof_linkers_lib/IRMOF-1_00.xyz`
- Summary:
  - Updated MOF helper export naming from verbose composition/smiles names to compact per-folder names such as `IRMOF-1_00.xyz`.
  - Added molecule de-duplication using composition plus rounded internal pair-distance fingerprints. If same-stem helper files already exist, export is skipped so prior visual-check files are preserved.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - Regenerated `test_on_irmof1/IRMOF-1.cif` normal/minimized: Path J remains 113/38 atoms.
  - IRMOF-1 helper exports now list one unique node file and one unique linker file: `mof_nodes_lib/IRMOF-1_00.xyz` and `mof_linkers_lib/IRMOF-1_00.xyz`.
- Follow-up risks:
  - The uniqueness key intentionally ignores absolute position/orientation; this is desired for duplicate helper fragments, but symmetry-distinct conformers with identical internal distances would be treated as duplicates.

## 2026-05-07 - Global capped-H-only final geometry cleanup
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_irmof1/IRMOF-1_frag_mof.xyz`
  - `test_on_irmof1/IRMOF-1_frag_mof_min.xyz`
  - helper exports under `mof_nodes_lib/` and `mof_linkers_lib/` from smoke tests
- Summary:
  - Added shared `optimize_capped_h_geometry_only(...)` in `BaseFragmenter`. It only moves hydrogens tracked as caps, then applies deterministic sp2 C-H and O-H cap geometry cleanup.
  - Routed MOF Path J, legacy MOF finalization, and COF finalization through the shared capped-H-only helper. Extraction finalization no longer calls `refine_h_geometry_with_rdkit(...)`, so full-molecule RDKit/UFF optimization is avoided for every path.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - IRMOF-1 normal/minimized Path J remain 113/38 atoms; capped C-O-H angles are 109.5 degrees.
  - Additional MOF Path J smoke checks: MgMOF74 normal/minimized 64/20 atoms; ZIF-8 normal/minimized 45/12 atoms. No UFF warnings appeared in successful runs.
- Follow-up risks:
  - This checkout did not contain the earlier COF fixture folders, so COF smoke tests could not be rerun here; the COF finalization call site is routed through the same shared helper and should be tested again when those fixtures are available.

## 2026-05-07 - PCN/NU normal fragments keep two metal nodes
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_pcn_series/PCN-61_frag_mof.xyz`
  - `test_on_pcn_series/PCN-68_frag_mof.xyz`
  - `test_on_nu_series/NU-100SP_frag_mof.xyz`
  - `test_on_nu_series/NU-108-Cu_frag_mof.xyz`
- Summary:
  - Narrowed the legacy Cu2 paddlewheel one-node suppression so it only applies to `Cu-BTC*`.
  - Routed normal PCN/NU MOFs around Path J so the legacy Path C candidate-count selector tests possible second nodes and keeps the smallest metal-complete fragment. Minimized fragments are unchanged.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - Normal metal counts: PCN-60 `Zn4` (276 atoms), PCN-61 `Cu4` (276 atoms), PCN-68 `Cu4` (420 atoms), NU-100SP `Cu4` (420 atoms), NU-108-Cu `Cu4` (516 atoms), NU-108-Zn `Zn4` (764 atoms).
  - Cu-BTC normal regression remains `Cu2` and 90 atoms.
  - Minimized smoke checks remain Path J without second-node completion: PCN-61 114 atoms; NU-108-Cu 174 atoms.
- Follow-up risks:
  - Visual inspection should confirm the 516-atom NU-108-Cu candidate is the intended opposite-side two-node fragment.

## 2026-05-06 - MOF moffragmentor Path J for IRMOF-1
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_irmof1/IRMOF-1_frag_mof.xyz`
  - `test_on_irmof1/IRMOF-1_frag_mof_min.xyz`
  - `mof_nodes_lib/*.xyz`
  - `mof_linkers_lib/*.xyz`
- Summary:
  - Added MOF Path J using installed `moffragmentor`: export detected nodes/linkers, select a central node, combine chemically attached linker images, and merge overlapping boundary atoms.
  - Normal mode keeps all attached linker images; minimized mode keeps the best attached linker image. Open terminal linker oxygens are H-capped with UniFrag capping/refinement methods; Path J locally adjusts only capped H geometry and avoids full-molecule RDKit/UFF optimization for Zn-containing fragments.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - `test_on_irmof1/IRMOF-1.cif` normal: Path J, 113 atoms, `Zn4 C48 H36 O25`; min: Path J, 38 atoms, `Zn4 C13 H6 O15`; capped C-O-H angles measure 109.5 degrees.
  - Exported 2 node files and 6 linker files for IRMOF-1 into `mof_nodes_lib/` and `mof_linkers_lib/`.
- Follow-up risks:
  - Visual inspection should confirm whether the minimized one-linker model should keep terminal capping hydrogens or an adjacent node image for MOF cases.

## 2026-05-06 - Apply direct coffragmentor combine to normal/min ZnPc
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof_min.xyz`
  - `test_on_cof_zn_pc_series/ZnPc-DPB_frag_cof_min.xyz`
- Summary:
  - Reverted the failed coffragmentor index/formula helper-selection machinery.
  - Added metallo-PC Path J for both normal and minimized fragments. It directly combines a Zn/N-rich `coffragmentor.py` node molecule with coffragmentor linker molecule image(s), including neighboring-cell images for normal fragments so all four node sides are represented, then duplicates the pair/set along the shortest lattice vector for the ZnPc dimer.
  - Normal ZnPc fragments remain on existing UniFrag Path D.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - ZnPc-DPB normal Path J: 322 atoms, `Zn2 B16 H80 C192 N16 O16`; minimized Path J: 166 atoms, `Zn2 B4 H32 C96 N16 O16`.
  - ZnPc-COF normal Path J: 210 atoms, `Zn2 B16 H48 C112 N16 O16`; minimized Path J: 138 atoms, `Zn2 B4 H24 C76 N16 O16`.
  - Regression minimized checks: COF-202 70 atoms, COF-300 74 atoms, COF-366 107 atoms.
- Follow-up risks:
  - Path J intentionally trusts coffragmentor node/linker molecules. Visual inspection should decide whether the combined node+linker pair needs an additional neighboring node/image or terminal capping in a later pass.

## 2026-05-04 - ZnPc minimized dimer linker balance
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof_min.xyz`
- Summary:
  - Updated metallo-PC minimization from one linker globally to one BDBA linker per retained ZnPc layer in the dimer.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - ZnPc normal: 242 atoms, `Zn2 B16 H64 C112 N16 O32`.
  - ZnPc minimized: 146 atoms, `Zn2 B4 H40 C76 N16 O8`; each ZnPc layer has one BDBA linker (`B-O = 8`, `O-H = 4`, `B-C = 4`).
  - Regression smoke checks: COF-366 Path D 182 atoms; COF-202 Path B 169 atoms; COF-300 Path C 149 atoms.
- Follow-up risks:
  - Visual inspection should confirm the selected one-linker-per-layer orientation is acceptable for stacked ZnPc variants.

## 2026-05-04 - ZnPc metallo-PC dimer support
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof.xyz`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof_min.xyz`
- Summary:
  - Updated ZnPc metallo-PC Path D to keep the nearest stacked ZnPc core along the short lattice axis, producing a two-layer dimer SBU.
  - Updated final component cleanup so metallo-PC mode preserves the two disconnected stacked principal layers.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - ZnPc normal: Path D metallo-PC core dimer, 242 atoms, `Zn2 B16 H64 C112 N16 O32`.
  - ZnPc minimized: Path D metallo-PC core dimer, 146 atoms, `Zn2 B4 H40 C76 N16 O8`.
  - Regression smoke checks: COF-366 Path D 182 atoms; COF-202 Path B 169 atoms; COF-300 Path C 149 atoms.
- Follow-up risks:
  - Dimer selection assumes the stacked partner lies along the shortest lattice axis with about 2.5-4.5 A axial spacing and low perpendicular offset.

## 2026-05-04 - COF fragmentation path tests and agent memory split
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `AGENTS.md`
  - COF test outputs under `test_on_cof_2xx_series/`, `test_on_cof_3xx_series/`, `test_on_cof_Por_series/`, and `test_on_cof_zn_pc_series/`
- Summary:
  - Added/validated COF-202 Path B tied-layer handling.
  - Validated COF-300 and COF-320 Path C behavior.
  - Validated COF-366 Path D porphyrin-core behavior.
  - Added ZnPc metallo-PC Path D priority and Zn radius.
  - Tuned ZnPc normal and minimized fragments: normal keeps all BDBA linkers; minimized keeps one full BDBA linker plus ZnPc fused benzene perimeter.
  - Split project coordination docs into `project-memory.md`, `project-decisions.md`, and `project-agent-log.md`.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - COF-202 normal: Path B layered set, 169 atoms.
  - COF-300 normal: Path C, 149 atoms.
  - COF-366 normal: Path D porphyrin core, 182 atoms.
  - ZnPc normal: Path D metallo-PC core, 121 atoms, `Zn1 B8 H32 C56 N8 O16`.
  - ZnPc minimized: Path D metallo-PC core, 73 atoms, `Zn1 B2 H20 C38 N8 O4`.
- Decisions made:
  - See `project-decisions.md` for Path B tied-neighbor and metallo-PC priority decisions.
- Follow-up risks:
  - RDKit UFF warns about Zn atom typing during H refinement; generation still succeeds.
  - Visual inspection remains important for new COF families because linker/SBU chemistry varies.

## 2026-05-08 - COF helper library export (nodes/linkers)
- Changed files:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Added COF-side helper export mirroring MOF helper libraries.
  - COF extraction now attempts `coffragmentor` and exports helper nodes/linkers to `cof_nodes_lib/` and `cof_linkers_lib/`.
  - Export uses compact filenames (`<stem>_00.xyz`), global duplicate checks per folder, and same-stem skip per folder.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - `ZnPc-DPB.cif` COF run prints helper export availability and writes:
    - `cof_nodes_lib/ZnPc-DPB_00.xyz`
    - `cof_linkers_lib/ZnPc-DPB_00.xyz`
  - Re-running the same CIF keeps file counts stable (`1` node file, `1` linker file for that stem) confirming skip/dupe behavior.
- Follow-up risks:
  - Some COFs (e.g., COF-6 in current heuristics) may return no helper node/linker set from `coffragmentor`; extraction still continues through UniFrag paths.

## 2026-05-09 - Enforce universal node+linker-first path ordering
- Changed files:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Removed MOF normal-mode family bypass before Path J; MOF now always attempts moffragmentor node+linker combine first.
  - Added generic COF Path J attempt (`coffragmentor` node+linker combine) before COF graph/radius fallback paths.
  - COF helper library export remains integrated with this first-pass workflow.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - `IRMOF-1.cif` logs `MOF Path J` first and writes 113-atom fragment.
  - `ZnPc-DPB.cif` logs `COF Path J` first for normal (161 atoms) and minimized (83 atoms).
  - `COF-6.cif` falls back to `COF Path J6` when generic coffragmentor node/linker combine is not usable.
- Follow-up risks:
  - COF Path J currently assembles node+attached linkers as returned/placed by coffragmentor and may differ in size/topology from older family-specialized outputs; visual QA remains required.

## 2026-05-09 - COF-6 / COF-66 layered dimer and topology-routing updates
- Changed files:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Added COF graph node+linker fallback Path J for cases where coffragmentor returns nodes without linkers (e.g., COF-66).
  - Added helper-library fallback exports for COF-6 decomposition to keep `cof_nodes_lib/` and `cof_linkers_lib/` populated.
  - Updated COF-6 minimum branch selection to preserve node/linker edge chemistry targets and maintain expected B/O balance.
  - Refactored COF dimer handling toward crystal-position-based layer selection and introduced topology-based routing gates so COF strict-path selection is driven by node/linker signatures instead of filename matching.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes during these changes.
  - COF-6 and COF-66 were repeatedly regenerated under `test_on_cof6xx_family/` while iterating on dimer placement, wrapping, and edge capping behavior.
- Follow-up risks:
  - Layered dimer geometry for edge cases should still be visually verified; topology-driven routing is in progress and may need one more stabilization pass for universal layered spacing/wrapping robustness.

## 2026-05-11 - Path J layered dimer fallback for face-to-face COFs
- Changed files:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `project-memory.md`
- Summary:
  - Added a global Path J layered-COF dimer fallback for coffragmentor node+linker outputs.
  - If `--cof-layer dimer`/auto requests a non-monomer output and the shortest lattice vector is a plausible face-to-face stacking spacing (2.5-5.0 A), the completed capped monomer fragment is duplicated by that stacking vector.
  - This addresses COF-LZU1 and COF-LZU8 where Path J previously produced dimer files with monomer-only atom counts.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - COF-LZU1 dimer: 102 atoms, centroid layer spacing 3.729 A; min dimer: 62 atoms, spacing 3.729 A.
  - COF-LZU8 dimer: 348 atoms, centroid layer spacing 4.093 A; min dimer: 148 atoms, spacing 4.093 A.
- Follow-up risks:
  - The fallback is conservative and uses the shortest lattice vector for Path J layered COFs; visually verify unusual non-layered COFs with a short lattice axis.

## 2026-05-24 - Automated Hydrogen Saturation and Neutralization Engine
- Changed files:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Integrated `force_qm_readiness` into all remaining finalization paths: the third COF finalization call, the COF fallback finalization call, and the BioMol sliding-window residue finalization call.
  - Implemented the missing `optimize_capped_h_geometry_only` helper method on `BaseFragmenter` to perform both `enforce_sp2_capped_h_geometry` and `enforce_capped_oh_geometry` locally on capped Hydrogen atoms.
  - Upgraded `enforce_sp2_capped_h_geometry` with an inline BFS cycle/ring checker (`is_in_ring`) to strictly gate the global aromatic SVD planarity plane fitting to cyclic/ring systems only. This prevents non-cyclic/peptide backbone cap geometries from being distorted.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes cleanly.
  - MOF folder-mode sweep on `test_on_mof_mg_based` succeeds: both `2015Mgdia3ASR1_frag_mof` and `2015Mgnan3FSR5_frag_mof` are successfully driven to `charge=0 multiplicity=1` and saved in `fragments_collection.extxyz`.
  - BioMol sliding-window sweep on `test_on_bio_mol/4c7n_clean_sigle.pdb` succeeds: all 47 sliding-window protein fragments are successfully neutralized to `charge=0 multiplicity=1` in `bio_fragments_collection.extxyz`. Window 43 has its C-terminal Carbon capping geometry perfectly resolved with no under-coordination warnings.
  - COF folder-mode sweep on `test_on_cof_others` succeeds: all generated normal and minimized COFs are successfully neutralized to `charge=0 multiplicity=1` and saved in `fragments_collection.extxyz`.
- Follow-up risks:
  - The SVD planarity fits now require a minimum cycle size of 3-8 containing the parent atom. Ensure that any future aromatic systems to be flattened are correctly identified by the inline BFS cycle check.

## 2026-05-25 - Revert Automated Hydrogen Saturation and Neutralization Engine
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Completely reverted the codebase changes back to commit `84b53581dbf3dcf52b2426ce98d17a8bdc7fc2ae` to undo the "H Saturation & Neutralization Engine" and related features.
  - This restores the previous performance characteristics of the fragmentation script by removing the expensive iterative RDKit determine bonds, protonation/deprotonation, and local geometry optimization passes.
- Validation:
  - Restored test files to consistent states.
  - Syntax check `/Users/omert/miniconda3/bin/python -c "import fragmentation_oop"` passes cleanly.
- Follow-up risks:
  - Fragments will not be automatically driven to charge=0/multiplicity=1 via RDKit saturation anymore; they will follow the original geometric and chemical capping rules.

## 2026-05-25 - Fast O(N) Electron Parity Check for QM Readiness
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
  - `project-decisions.md`
- Summary:
  - Added a module-level `_ATOMIC_NUMBERS` lookup table covering all elements common in MOFs/COFs/bio molecules.
  - Added `_check_multiplicity(species, label)` — a pure O(N) function that sums atomic numbers (Z_sum), assuming charge=0, and checks whether N_elec = Z_sum is even (singlet, multiplicity=1) or odd (radical doublet, multiplicity=2).
  - Integrated the check into `_write_extxyz`: runs before writing, prints a console status or `[QM WARNING]`, and embeds `multiplicity=N` into the ExtXYZ comment line header.
  - No adjacency matrix, no distance matrix, no RDKit calls. Purely a sum modulo 2.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - `test_on_mof_mg_based` folder sweep completed in 36.9s (identical to reverted baseline of 37.2s — zero overhead added).
  - All 4 fragments correctly classified: `2015Mgdia3ASR1_frag_mof` (Z_sum=154, mult=1 ✅), `2015Mgnan3ASR11_frag_mof` (Z_sum=932, mult=1 ✅), `2015Mgnan3ASR11_frag_mof_min` (Z_sum=434, mult=1 ✅).
  - `[QM WARNING]` correctly fired for `2015Mgnan3FSR5_frag_mof` (Z_sum=293, multiplicity=2), correctly identifying a genuine radical fragment that needs attention before QM calculations.
- Follow-up risks:
  - The check assumes charge=0. If a fragment has a formal nonzero charge, the parity result may not match the true multiplicity. A future charge-estimation step could refine this.




## 2026-05-25 - Odd-Electron Auto-Fix: remove non-C capping H / add H to O/N
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
  - `project-decisions.md`
- Summary:
  - Added `BaseFragmenter.fix_odd_electron_multiplicity(species, coords, capped_h_indices, label)` to `BaseFragmenter`.
  - Logic:
    1. Compute Z_sum parity (O(N)); return immediately if already even.
    2. Scan capping H atoms; find parent heavy atom (nearest non-H within 1.5 Å, O(N) scan).
    3. If parent is NOT Carbon (C) → eligible for removal. Priority: O=3, N=2, B=1, others=0.
    4. Remove highest-priority candidate, remap capped_h_indices, print `QM-Fix: Removed capping H from {sym}[{idx}]...`.
    5. Fallback: if no non-C capping H found, scan for under-protonated O (1 heavy + 0 H) or N (2 heavy + 0 H) → add H geometrically via `place_capping_h`, print note to user.
    6. If neither works, print `[QM WARNING] could not fix...`.
  - Wired into 7 finalization sites across MOFFragmenter, COFFragmenter, BioMolFragmenter, immediately after `optimize_capped_h_geometry_only`.
  - Carbon capping H (sp2 aromatic and sp3 aliphatic) are never removed, per user requirement.
- Validation:
  - Syntax check passes.
  - `test_on_mof_mg_based` sweep: `2015Mgnan3FSR5_frag_mof` auto-fixed:
    `QM-Fix: Removed capping H from O[47] to achieve even electron count for 'mof_fragment'.`
    → Then `_write_extxyz` reports: `electron count OK (Z_sum=292, multiplicity=1)`.
  - Runtime: **12.9 seconds** (down from 36.9s in previous run — cached supercell reuse benefit).
- Follow-up risks:
  - The fallback H-add path has not yet been triggered in real test data. Should be tested with a deliberately radical-only-on-C fragment.
  - The neighbour-count used for the add-H fallback (cutoff 2.2 Å) is heuristic and may misclassify some bridging O atoms in MOFs.


## 2026-05-25 - Functional-Group-Aware Removal Priority (geometric 2-hop walk)
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Added `BaseFragmenter._classify_cap_removal_priority(h_idx, species, coords_arr)` static method.
  - Pure geometric 2-hop neighbour walk: Hop 1 finds parent (nearest non-H ≤ 1.5 Å); Hop 2 finds grandparents (heavy atoms ≤ 2.2 Å of parent). No RDKit, no bond graph object.
  - Priority scores follow pKa ordering:
    - Sulfonate (-SO₃H) pKa ~-1 → 100
    - Phosphonate (-P(O)(OH)) pKa ~2 → 90
    - Carboxylate (-COOH) pKa ~4 → 80
    - Sulfinic (-SO₂H) pKa ~8 → 75
    - Thiol (-SH) pKa ~10 → 70
    - Phenol (Ar-OH) pKa ~10 → 60
    - Alcohol (-C-OH) pKa ~15 → 50
    - Sulfonamide (-SO₂NH-) pKa ~10 → 45
    - Primary amine (-NH₂) pKa ~35 → 20
    - Secondary amine (-NH-) pKa ~35 → 15
    - Amide (-CO-NH-) pKa ~25 → 10 (backbone, avoid)
    - Boron (-BH-) → 5
    - Carbon → 0 (NEVER remove)
  - Updated `fix_odd_electron_multiplicity` to use the new classifier.
  - Console output now includes functional group name and priority score.
  - Fallback add-H site also uses the same classifier for consistent ranking.
- Validation:
  - Syntax check passes.
  - `test_on_mof_mg_based` sweep: FSR5 fixed via [carboxylate] group (priority=80), runtime 13.3s.
  - Group name confirms carboxylate O was selected, not a random O.
- Decision rationale:
  - Pure geometric detection avoids RDKit overhead (~1000x faster).
  - 2-hop walk is sufficient for all practically relevant functional groups.
  - Priority follows pKa chemistry for chemically sound QM fragments.

## 2026-05-25 - Fix geometric flattening bug for capping H
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Fixed a regression in `enforce_sp2_capped_h_geometry` where `is_in_ring` was accidentally removed, causing ALL capping H atoms near 5 heavy atoms to be flattened onto an SVD plane.
  - Restored the `is_in_ring(idx)` BFS check so only actual ring capping H atoms (like phenyls) are projected onto the aromatic plane.
  - Added `overcoordinated-O` (priority=100) to `_classify_cap_removal_priority` to ensure capping Hs mistakenly added to already-saturated water molecules (creating hydronium) are prioritized for removal first.
- Validation:
  - `fragments_collection.extxyz` for `2015Mgnan3ASR10` now has correctly oriented capping Hs on the metal-coordinated water molecules, without the flat degenerate geometries.


## 2026-05-25 - Fix double-counting in stray capping logic
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Addressed a bug where the capping logic double-counted struct hydrogens that were already loaded into the fragment `species` list, causing it to incorrectly skip adding a capping H and leaving OH radicals instead of H2O.
  - Implemented a `capping_hs_added` counter per stray group loop to accurately track explicitly added hydrogens.
  - Lowered `min_hh` and `min_o_contact` to `1.2` for stray group capping to ensure valid capping H placements aren't incorrectly rejected due to proximity to struct atoms.
  - Verified that metal-coordinated oxygen atoms now correctly max out at 2 hydrogens, entirely eliminating the "3 H on O" bug reported by the user while preserving the QM-Fix functionality.
- Known Risks/Follow-ups:
  - None immediately apparent.

## 2026-05-25 - Fix QM-Fix H removal priority: carboxylate vs water-on-metal
- Changed files:
  - `fragmentation_oop.py` (line 624)
  - `project-agent-log.md`
- Summary:
  - Fixed wrong priority in `_classify_cap_removal_priority`: the `overcoordinated-O` check was triggering at `h_neighbors >= 2` (priority 100), which incorrectly flagged normal water molecules (H₂O with 2 H) coordinated to metal atoms as overcoordinated, causing QM-Fix to strip an H from water instead of from a carboxylate group (priority 80).
  - Changed threshold from `>= 2` to `>= 3`. An oxygen with 2 H is normal water; only 3+ H is genuinely overcoordinated (hydronium-like).
  - After the fix, QM-Fix correctly reports: `QM-Fix [carboxylate]: Removed capping H from O[61]` and both water molecules retain their 2 H atoms.
- Validation:
  - Ran `fragmentation_oop.py` on `test_on_mof_mg_based/2015[Mg][nan]3[ASR]10.cif --kind mof`
  - Confirmed QM-Fix now removes from carboxylate (priority 80) instead of water-on-metal (priority 40)
  - Confirmed all metal-coordinated water O atoms have exactly 2 H

## 2026-05-26 - Fix polynuclear SBU merging & O over-capping
- Changed files:
  - `fragmentation_oop.py` (lines ~1354-1397 and ~1522-1536)
  - `project-agent-log.md`
- Summary:
  - **Node merging**: Added logic after moffragmentor node selection to detect and merge nearby nodes that belong to the same polynuclear SBU. moffragmentor may split a dinuclear Mg-O-Mg unit into individual 1-atom nodes. The new code checks metal-metal distances (< 3.5 Å, with periodic images) and merges matching nodes, using a BFS-like while loop to handle chains. Prints "Merged N moffragmentor nodes into one polynuclear SBU."
  - **Smart O capping limit**: When capping an O atom that has a heavy non-metal neighbor (C, P, S), limit total H to 1 (replacing only the missing metal bond). Only pure water-like O (bonded only to metals/H) can have up to 2 H. This prevents over-protonation of coordinating oxygens in phosphonates and carboxylates.
- Validation:
  - Tested on `2015[Mg][nan]3[ASR]1.cif`: Fragment now has 2 Mg (was 1), no suspicious O atoms with 2 H, electron count 312 (even), no QM-Fix needed.
  - Bridging oxygens O[14] and O[30] correctly connect to both Mg atoms with no capping H.
  - Boundary oxygens O[24], O[35], O[40] correctly have 1 H each (replacing the missing 3rd/4th Mg bond).

## 2026-05-26 - Fix --nmetals for infinite SBU chains
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Addressed issue where `--nmetals` flag was ignored for 1D rod MOFs because the `moffragmentor` node extraction path overrode it and simply merged whatever nodes it found in the asymmetric unit.
  - Added a check in `_try_moffragmentor_node_linker_fragment`: if `getattr(result, "has_1d_sbu", False)` is True, we now immediately return `None`.
  - This allows the script to intentionally fall back to the radius-based Path B (Infinite SBU path) which correctly builds the supercell and slices the 1D chain to exactly `--nmetals` length.
  - Also fixed a bug where `structure_stem` was undefined in the fallback export logic.
- Validation:
  - Tested on `2015[Mg][dia]3[ASR]2.cif` with `--nmetals 4`. The fragment correctly extracted a 4-metal segment (Mg4) via Path B, whereas previously it merged the unit cell nodes into an uncontrollable cluster size.

## 2026-05-26 - Format labels to CamelCase
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Updated the generation of fragment labels written to `.extxyz` files (the `label=` property) to use CamelCase instead of underscores.
  - Replaced `_frag_mof_min` with `FragMofMin`, `_frag_cof` with `FragCof`, `_w001` with `W001`, etc.
  - Removed brackets and underscores from the base structure name.
- Validation:
  - Tested on `2015[Mg][ins]3[ASR]1.cif`. The output label successfully formatted as `2015Mgins3ASR1FragMof`.

## 2026-05-28: Fixed Mg MOF fragmentation and Carboxylate Hydrogen capping
- **Changes**:
  - `fragmentation_oop.py`: Updated `_try_moffragmentor_node_linker_fragment` to reject `moffragmentor` output if any generated linker contains a metal atom (e.g. Mg). This forces the tool to use our more robust fallback paths (Path A, B, C).
  - `fragmentation_oop.py`: Increased the topology detection radius cutoff for SBU metals (`get_all_neighbors(r=...)`) from `3.6 Å` to `5.5 Å` to correctly identify dinuclear and infinite metal SBUs linked via longer `Metal-O-Metal` oxygen bridges (like Mg-O-Mg which can be ~5.2 Å).
  - `fragmentation_oop.py`: Fixed the hydrogen capping on carboxylate oxygen atoms. Modified `_cap_path_j_open_oxygens` and the main BFS extraction capping logic to track `capped_central_atoms`. When an `O` atom coordinated to a `C`, `P`, or `S` is capped with an `H`, the central heavy atom is flagged. If a second `O` on the same central atom requires capping, it is skipped. This correctly yields `-COOH` instead of `-C(OH)2`.
- **Validation**:
  - Tested on `2009[Mg][lvt]3[ASR]1.cif` in `test_on_mof_mg_based/`. The script correctly rejected the faulty `moffragmentor` output, identified the SBU as a discrete dinuclear Mg cluster (`SBU size: 2`), and processed it via Path C. The carboxylates are now properly capped as `-COOH`.
- **Risks/Follow-ups**:
  - Increasing the radius to `5.5 Å` could potentially group unrelated metal centers in extremely dense MOFs, although the structural topology checks mitigate this risk.

## 2026-05-28: Fixed min version generation for MOF fragments
- **Changes**:
  - `fragmentation_oop.py`: Updated `_get_first_ring_keep_heavy` to explicitly find the first cyclic structure (ring size <= 8) connected to the bridge atoms and stop traversing further.
  - This logic replaces both the flawed `is_bridge` BFS logic in `_first_connected_ring_fragment` (which failed on infinite unrolled MOF linkers) and the inline `internal_bonds <= 1` core-pruning loop in `_get_fragment` (which kept the entire rigid multi-ring linker).
- **Validation**:
  - Re-ran fragmentation on `2009[Mg][lvt]3[ASR]1.cif`. The `min` version successfully kept exactly 1 complete linker and pruned the remaining 3 linkers down to only their first phenyl ring attached to the node, resulting in exactly the expected formula (`C61H44MgN2O11`).
- **Risks/Follow-ups**:
  - None immediately apparent.

## 2026-06-01 - Fix QM Geometry Flattening and Multiplicity Failures
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Fixed 'Zero distance between atoms' error in QM optimizations. Added a check in `enforce_sp2_capped_h_geometry` to skip SP2 planar enforcement on carbon atoms that already have >= 2 hydrogens (e.g. CH2 groups). Previously, trimming aliphatic linkers caused both capping hydrogens to be artificially collapsed onto the exact same bisector vector.
  - Fixed 'multiplicity is odd' errors by adding progressive fallback relaxation (down to `min_hh=0.0`) in `fix_odd_electron_multiplicity` when strict steric constraints prevented placing the parity-fixing H atom.
- Validation:
  - Multiplicity and zero-distance failures resolved on the reported Mg-based test cases.


## 2026-06-01 - Fix: preserve carboxylate -COOH oxygens in minimized MOF fragments
- Changed files:
  - `fragmentation_oop.py`
  - `project-agent-log.md`
- Summary:
  - Fixed a bug in `_get_first_ring_keep_heavy` where the second oxygen of a carboxylate group (-OH or =O, bonded to the carboxylate C) was being dropped during minimize trimming. The method correctly kept the bridge O (bonded to metal) and the carboxylate C (first-layer neighbor), but the other O was excluded because it is not part of any ring and not on the BFS path to the ring.
  - Added a `species_map` optional parameter to `_get_first_ring_keep_heavy`. After computing the ring-based `keep_atoms`, a single-pass expansion is done: for every kept C/Si/B atom, any directly-bonded O/S/P/N that was not already kept is added to `keep_atoms`. This is intentionally limited to one pass to avoid expanding into the rest of the linker arm.
  - Updated both call sites: `_first_connected_ring_fragment` passes `{i: linker_species[i] for i in heavy}`, and `_get_fragment` passes `{ka: supercell[ka].species_string for ka in comp}`.
- Validation:
  - The full -COOH group (both oxygens) is now always preserved in the minimized fragment of Mg-based MOFs.

## 2026-06-21 - CCDC API activation check, unmodified CIF extraction, and classification
- Changed files:
  - `project-agent-log.md`
  - `runUniFrag/analyze_cifs.py`
- Summary:
  - Verified CCDC Python API (`ccdc` v3.7.1) activation and database connection in Miniconda environment `/Users/omert/miniconda3`.
  - Discovered local SQLite database file at `/Users/omert/CCDC/ccdc-data/csd/as601be_CIP.sqlite`.
  - Resolved `UserWarning` about missing database by setting `CSD_DATA_DIRECTORY=/Users/omert/CCDC/ccdc-data/csd` environment variable.
  - Successfully ran `runUniFrag/fetch_cifs_from_csd.py` using the CCDC API to extract unmodified CIFs for CR (from CSV metadata) and NCR (from folder filename scan) datasets.
  - Wrote and executed `runUniFrag/analyze_cifs.py` to classify the extracted CIF files (renamed to `CSD-unmodified/`) into CR (ASR, FSR, Ion) and NCR subsets.
- Validation:
  - Extracted 11,338 unmodified CIFs out of 11,367 requested REFCODEs into `runUniFrag/CSD-unmodified/` (29 REFCODEs were missing/not found in CSD database).
  - Out of 11,338 extracted: 5,006 are exclusively CR, 5,609 are exclusively NCR, and 723 are present in both subsets.

## 2026-06-22 - Zn-based Unmodified CIF extraction, merging, collection, and auto-continue flag
- Changed files:
  - `project-agent-log.md`
  - `project-memory.md`
  - `fragmentation_oop.py`
  - `runUniFrag/merge_zn_cifs.py`
  - `runUniFrag/collect_all_zn_cifs.py`
- Summary:
  - Wrote and executed `runUniFrag/merge_zn_cifs.py` to extract all Zn-based unmodified CIFs from `runUniFrag/CSD-unmodified/`.
  - Merged these extracted CIFs into parallel subdirectories in `runUniFrag/zn_cifs/` with name prefixes `unmodified_` (e.g. `unmodified_CR_ASR`, `unmodified_CR_FSR`, `unmodified_NCR`, `unmodified_CR_ASR_FSR_merged`).
  - Mapped each Zn-based structure to its corresponding category using CR metadata and NCR folder listings, resolving duplicates by base REFCODE (ensuring unique filename symlinks).
  - Wrote and executed `runUniFrag/collect_all_zn_cifs.py` to compile ALL unique Zn-based MOFs into a single directory `runUniFrag/zn_cifs_noduplicated/`. If a structure exists in both unmodified and modified, it prioritizes the unmodified version and falls back to the modified version if the unmodified file is missing (to maximize coverage). Rename the target files to a clean `[REFCODE].cif` format.
  - Implemented automatic continue/resume behavior by default in `fragmentation_oop.py` for folder mode (MOF, COF, and Bio modes). It skips already completed files listed in the CSV summary and preserves the existing ExtXYZ collection, loading its signatures to properly handle duplicate detection.
  - Added the `--overwrite` flag to allow users to bypass auto-continue and force starting folder-mode runs from scratch.
- Validation:
  - Identified 2,439 Zn-based structures out of the 11,338 unmodified CIFs.
  - Successfully symlinked categories under `zn_cifs/` (1223 ASR, 919 FSR, 0 Ion, 1255 NCR, 1247 merged).
  - Collected exactly **2,443** unique Zn-based structures in `zn_cifs_noduplicated/` (2,439 unmodified structures + 4 modified-only fallback structures).
  - Verified compilation of `fragmentation_oop.py` runs successfully.

## 2026-06-23 - Collection of unique Computationally Ready (CR) MOFs from modified and unmodified subsets
- Changed files:
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Updated `runUniFrag/collect_cr_cifs.py` to include robust metal element parsing and mixed-metal filtering (excluding structures with more than one metal type).
  - Executed the updated script to compile unique single-metal CR MOFs from `CSD-modified` (ASR/FSR/Ion) and `CSD-unmodified` subsets.
  - Used CR metadata `CR_data_CSD_modified_20250227.csv` to map filenames (`coreid`) to base parent REFCODEs and deduplicated by 6-letter base parent CSD REFCODE (ensuring exactly one representative file per framework).
  - Priority hierarchy: ASR -> FSR -> Ion -> Unmodified fallback.
  - Saved output files to `runUniFrag/cr_cifs_noduplicated/` in a clean `[REFCODE].cif` format.
  - Generated summary report at `runUniFrag/cr_collection_summary.txt`.
  - Created and ran `runUniFrag/plot_cr_metals.py` to compile the metal center distribution histogram for the final `cr_cifs_noduplicated/` collection.
- Validation:
  - Confirmed exactly **5,226** single-metal files collected in `runUniFrag/cr_cifs_noduplicated/` (515 mixed-metal structures were successfully filtered out and skipped).
  - Distribution of source files in the final collection:
    - Copied from ASR: 4,810
    - Copied from FSR: 121
    - Copied from Ion: 295
    - Copied from Unmodified fallback: 0
  - Verified that all output filenames are exactly 6 uppercase letters (base REFCODE format).
  - Confirmed no empty (0 bytes) files are present in the output folder.
  - Successfully generated the updated `mof_metals_histogram.png` plot and copied it to the brain artifact directory.
  - Top metal counts in the single-metal CR collection: Zn (1220), Cu (664), Cd (602), Co (539), Mn (209), Ni (204), Eu (170), Ag (147), Tb (138), Gd (113), and Others (1220).
  - Created and ran `runUniFrag/collect_zn_cr_cifs.py` to extract all Zn-based single-metal CR MOFs from `cr_cifs_noduplicated/`.
- Validation:
  - Extracted exactly **1,220** Zn-based structures into `runUniFrag/zn_cr_cifs_noduplicated/`.
  - All filenames are exactly `[REFCODE].cif` and match standard 6-letter base refcodes.
  - Confirmed no empty (0 bytes) files are present in the output folder.

## 2026-06-23 - Atom types coverage analysis and parent-fragment comparison
- Changed files:
  - `runUniFrag/analyze_atom_types.py` [NEW]
  - `runUniFrag/atom_types_analysis.md` [NEW]
  - `runUniFrag/atom_types_distribution.png` [NEW]
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Developed and ran `runUniFrag/analyze_atom_types.py` to compare Sybyl atom types of parent MOFs and fragments.
  - Scanned all 1,220 single-metal Zn parent structures, identifying 27 unique Sybyl atom types.
  - Typed the 1,465 fragments using two approaches:
    - **Method A (Mapped)**: Atoms are mapped back to parent crystal fractional coordinates (modulo 1) to inherit their parent environment types.
    - **Method B (Direct)**: The fragment is typed as an isolated molecule in a temporary CIF file to let CCDC perceive chemistry directly.
  - Achieved **88.89%** coverage (24 out of 27 parent types) for both methods.
  - Identified the three missing types (`S.2`, `S.o`, `Se`) and chemically analyzed their omission:
    - `Se` in `BAFVOX.cif` and `S.o` in `XINFUW.cif` represent pore solvent/guest molecules, which were correctly discarded during fragmentation.
    - `S.2` in `BUCXUT.cif` was omitted during fragment partitioning of the framework.
  - Generated a comparative relative-frequency bar chart of the top 15 atom types and a detailed markdown report.
- Validation:
  - Confirmed 0 mapping errors across all 1,465 fragments.
  - Checked that output markdown and distribution plot are correctly written to the workspace and the brain artifacts directory.

## 2026-06-24 - Zn coordination environment analysis and parent-fragment comparison
- Changed files:
  - `runUniFrag/analyze_zn_coordination.py` [NEW]
  - `runUniFrag/zn_coordination_analysis.md` [NEW]
  - `runUniFrag/zn_coordination_distribution.png` [NEW]
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Developed and ran `runUniFrag/analyze_zn_coordination.py` to compare Zn coordination environments.
  - Identified 68 unique coordination shells in the 1,220 parents (CN=0 to CN=9 with O, N, S, Zn, and Halogens).
  - Achieved **67.65%** coverage (46 / 68 types) for Method A (Mapped) and **70.59%** (48 / 68 types) for Method B (Direct).
  - Explained that the missing environments correspond to lone guest ions, incomplete crystallographic CIF parameters, or rare halogen/metal clusters.
  - Generated a comparative relative-frequency bar chart of the top 10 Zn environments and a detailed markdown report.
- Validation:
  - Confirmed 0 mapping errors across all 3,010 Zn centers in the fragments collection.
  - Checked that output markdown and distribution plot are correctly written to the workspace and the brain artifacts directory.

## 2026-06-24 - Extract parent Zn MOFs with Zn CN=0, 1, or 3
- Changed files:
  - `runUniFrag/extract_low_coordination_mofs.py` [NEW]
  - `runUniFrag/zn_low_coordination_parents.csv` [NEW]
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Wrote and ran `runUniFrag/extract_low_coordination_mofs.py` to filter parent MOFs containing Zn centers with CN in [0, 1, 3].
  - Identified exactly **1,198** parent structures containing at least one low-coordination Zn atom.
  - Exported the results mapping refcodes to their specific CN=0, 1, 3 environments.
- Validation:
  - Checked that output CSV was successfully created in both the workspace and the brain artifacts directory.

## 2026-06-24 - Purify Zn parent MOF CIF collection
- Changed files:
  - `runUniFrag/purify_zn_cifs.py` [NEW]
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Wrote and ran `runUniFrag/purify_zn_cifs.py` to separate CIFs containing heavy/semi-metal elements (`I`, `Si`, `Br`, `B`, `Se`, `As`) into a dedicated folder.
  - Successfully moved exactly **115** parent CIFs to `runUniFrag/zn_cr_cifs_noduplicated/cifs_heavy_elements/`.
  - Left exactly **1,105** purified structures in the main `runUniFrag/zn_cr_cifs_noduplicated/cifs/` directory.
- Validation:
  - Verified that the remaining file count is exactly 1,105 and the moved files are present in the target directory.

## 2026-06-24 - Revert Boron (B) purification exclusion
- Changed files:
  - `runUniFrag/purify_zn_cifs.py`
  - `project-agent-log.md`
- Summary:
  - Removed Boron (`B`) from the exclusion list in `runUniFrag/purify_zn_cifs.py`.
  - Moved the **6** Boron-containing structures (`CUGFAN`, `FAFJIH`, `GATSAY`, `HABREJ`, `HABRIN`, `WORSUT`) from `cifs_heavy_elements/` back to `runUniFrag/zn_cr_cifs_noduplicated/cifs/`.
  - Restored the main `cifs/` directory to exactly **1,111** structures, leaving exactly **109** structures in `cifs_heavy_elements/`.
- Validation:
  - Verified that all 6 Boron files were moved back successfully and file counts match.

## 2026-06-24 - Guest Molecule Detection and Removal
- Changed files:
  - `runUniFrag/remove_guests.py` [NEW]
  - `runUniFrag/guest_removal_report.md` [NEW]
  - Modified 50 parent CIF files under `runUniFrag/zn_cr_cifs_noduplicated/cifs/`
  - Created backup of 50 original structures in `runUniFrag/zn_cr_cifs_noduplicated/cifs_backup_guests/`
  - `project-memory.md`
  - `project-agent-log.md`
- Summary:
  - Developed and executed a python script to scan the 1,111 purified Zn parent MOFs for guest molecules.
  - Used CCDC Python API pre-filtering to identify 50 structures containing guest components (components without Zinc).
  - Isolated the framework and removed the guests from the 50 structures in-place using Pymatgen and NetworkX connected component graph analysis.
  - Implemented JmolNN as a robust fallback bonding strategy for 6 problem structures (5 that failed with Voronoi errors in CrystalNN, and 1 structure, RUGXOI, that over-connected its guests to the framework in CrystalNN).
  - Verified each modified structure with CCDC to confirm complete guest removal while keeping framework structures intact.
  - Generated a comprehensive markdown report listing all cleaned refcodes, original vs cleaned formulas, and removed atom counts.
- Validation:
  - Confirmed all 50 parent CIF structures were successfully purified.
  - CCDC post-cleaning verification shows 0 remaining guest components across all 50 files.
  - Backups of all 50 original structures are safely stored in `cifs_backup_guests/`.






## [2026-09-17] COF stacked-dimer verification and planar amine capping fix
- Files touched:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
- Summary:
  - Verified stacked-dimer construction on a 20-structure HCNO subset (`runUniFrag/CoRECOF/test_dimer20`). Dimer construction itself is correct: 8 structures dimerized, each verified as a pure translation of the monomer with `|T|` matching the parent CIF's shortest lattice vector exactly (3.40-3.57 A).
  - Diagnosed sub-2.4 A interlayer contacts (worst 1.72 A) in 25 of 30 dimer frames. Heavy-atom frameworks are perfectly planar; the clash comes entirely from capping H on terminal nitrogen pyramidalized out of plane by the UFF refinement.
  - Added `COFFragmenter._planarize_conjugated_amine_caps` and wired it into `_cap_severed_double_bond_sites` (see `project-decisions.md` 2026-09-17).
- Validation (`test_dimer20` vs `test_dimer20_planar`, same 20 CIFs, 4 procs, exit 0):
  - Out-of-plane displacement of terminal N-H: mean 0.76 A -> 0.09 A.
  - Minimum interlayer contact across all dimers: 1.72 A -> 2.55 A; median 1.81 A -> 3.45 A.
  - Frames with a sub-2.4 A interlayer contact: 25/30 -> 0/29.
  - Every shared fragment label has **identical composition and identical heavy-atom coordinates** before and after; only hydrogen positions moved.
- Follow-up risks / open items:
  - `411.cif` (parent c = 3.70 A) and `760.cif` (c = 3.56 A) lie inside the 2.5-5.0 A stacking window but still produce monomers: both route through fallback paths (Path A single node, Path D porphyrin) rather than Path J, where the dimer rule lives. Not addressed.
  - `_chemical_identity_key` uses the heavy-atom formula only, with no geometry. `557FragCofOnlyNode` and `637FragCofMin` are both `C18H24N6O6` but are different molecules; they collide on this key, so only one survives and which one depends on multiprocessing arrival order. Pre-existing, unrelated to this change.
  - `1102FragCof` atom 72 remains an under-coordinated carbon (Path D).

## [2026-09-17] COF layered-dimer rule extended to fallback paths (411, 760)
- Files touched:
  - `fragmentation_oop.py`
  - `project-decisions.md`
  - `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_dimer_411_760` (isolated), `runUniFrag/CoRECOF/test_dimer20_pathfix` (regression).
- Summary:
  - Traced the two missing dimers: `411.cif` routes to Path A (single node) and `760.cif` to Path D (porphyrin core); neither reaches Path J, where the layered-dimer rule lived. The Path A/B/C/D branch's own dimer code required an explicit `layer_mode == "dimer"` and never fired under the default `"auto"`.
  - Extended the rule to that branch and gated it on fragment planarity (see `project-decisions.md` 2026-09-17).
- Validation:
  - `411`: now `-> COF layered dimer: second layer added along lattice axis 2 (|T| = 3.70 A)`, verified an exact translation (max atom mismatch 0.000 A) matching the parent c.
  - `760`: now `-> COF layered dimer skipped: ... (flatness = 2.51 A, alignment = 1.00)`, correctly remaining a monomer.
  - 20-structure regression (exit 0, no new warnings): dimer frames 29 -> 33; only `411FragCof`, `411FragCofMin`, `411FragCofOnlyNode` changed; all other labels kept composition and component count.
- Follow-up risks / open items:
  - `411FragCof` and `411FragCofOnlyNode` retain interlayer contacts of 1.84 A and 1.99 A. Cause is *not* the dimer construction: it is the `-CH3` cap that `_cap_severed_double_bond_sites` places on the severed imine carbon. The parent is perfectly flat with no out-of-plane hydrogen anywhere, so that carbon is sp2 with a single in-plane H in the real material; saturating it to a tetrahedral methyl necessarily puts an H ~0.93 A out of the layer. Unlike the amine case this cannot be fixed by projection - a planar `-CH3` would be chemically wrong. The faithful fix is to cap the aldimine terminus as `Ar-CH=NH` instead of `Ar-CH3`, which preserves sp2 and planarity but adds a heavy atom and would affect `_chemical_identity_key` grouping. Deferred pending a decision.
  - `_chemical_identity_key` formula-only collision (`557FragCofOnlyNode` vs `637FragCofMin`, both `C18H24N6O6`) still swaps nondeterministically between runs.

## [2026-09-17] COF building blocks: k+k classification, per-site linker placement, bond-order capping
- Files touched: `coffragmentor.py`, `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_411_blocks` (isolated 411), `runUniFrag/CoRECOF/test_dimer20_final4` (regression).
- Trigger: user review of `411` reported (1) no linker in the extxyz, (2) a node with a redundant carbon past the N on two of three arms, (3) a linker missing parts in both normal and min versions, with reference images of the expected node and linker.
- Root cause: all three were one failure. `coffragmentor` already produced both correct blocks but labelled both as nodes (k+k COF), so Path J bailed and Path A's fallback block growth produced the malformed output. See `project-decisions.md` 2026-09-17 for the three fixes.
- Validation - `411` now matches the user's reference images exactly:
  - `cof_nodes_lib/411_00.xyz` = `H6 C9 N3 O1` (triformylphenol core), `cof_linkers_lib/411_00.xyz` = `H15 C24 N3` (tris(aminophenyl)benzene).
  - Fragments: `411FragCof` C162H126N18O2 (node + 3 full linkers, dimer), `411FragCofMin` C90H70N10O2, `411FragCofOnlyNode` C18H18N6O2, `411FragCofOnlyLinker` C48H42N6.
  - Node nitrogen is now the neutral aldimine `Ar-CH=NH` (1 H at a 1.27 A C=N); linker nitrogen stays aniline `Ar-NH2` (2 H at 1.40 A); the reconnected junction N carries 0 H with both 1.27 and 1.40 A bonds.
- Validation - 20-structure regression (exit 0, no errors):
  - Fragment frames 65 -> 74; `103`, `1102`, `1120` now yield node/linker blocks for the first time, `637` gained a Min.
  - Min coverage 19/20 -> **20/20**.
  - Over-valent terminal sites **74 -> 4** (all 4 in the distorted `760.cif`); correctly capped sites 275 -> 363.
  - Dimers 33, minimum interlayer contact 1.84 -> **2.55 A**, zero sub-2.4 A contacts. The 411 methyl clash reported earlier disappeared on its own: that carbon is now capped as a planar sp2 `=CH-` rather than a tetrahedral `-CH3`.
- Follow-up risks / open items:
  - `760.cif` has a geometrically distorted parent (C-N 1.154-1.338 A, all C-H exactly 1.140 A). Its 4 over-valent sites and 10 bridging-H artefacts trace to that, not to UniFrag. Consider screening the 884 set for non-physical bond lengths before the full run.
  - `_chemical_identity_key` heavy-atom-formula collision still unresolved.
  - `1102FragCof` atom 72 under-coordinated carbon (Path D).

## [2026-09-17] Vinylene COF linkage (704) and sp2 carbon cap planarization
- Files touched: `coffragmentor.py`, `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_704_blocks`, `runUniFrag/CoRECOF/test_dimer20_v2` (regression), `runUniFrag/CoRECOF/test_hcno100` (100-structure run).
- Trigger: user reported "similar problems for 704" after the `411` fixes.
- Root cause: `704` is a vinylene-linked sp2-carbon COF. Its `Ar-CH=CH-Ar` linkage was in none of the cleavage rules, so `COF.fragment()` severed **nothing** and returned whole periodic networks as blocks; Path J bailed and the fallback path produced garbage. Distinct from the `411` failure (which was a classification tie), though both ended in the same fallback.
- Validation:
  - `704`: node lib `H5 C11 N1`, linker lib `H27 C45 N3`; fragments `704FragCof` C140H98N10, `Min` C66H48N4, `OnlyNode` C11H11N, `OnlyLinker` C45H33N3. Previously a single malformed `H10 C14 N1` node block and no linker.
  - 20-structure regression (exit 0, no errors): 75 frames, correctly capped terminal sites 367, Min coverage 20/20, 33 dimers, minimum interlayer contact 2.55 A, **zero** sub-2.4 A contacts.
  - `851` is a mixed imine/vinylene COF; its fragments are now decomposed further, which is correct, and the clashes the new caps briefly introduced were removed by extending planarization to sp2 carbon.
- Follow-up risks / open items:
  - Cleavage-rule coverage is still enumerated by hand (boroxine, imine, oxazole, dioxin, vinylene). A COF whose linkage is absent from that list fails silently through the fallback path rather than reporting that it recognised nothing. Worth a diagnostic that flags "zero severed bonds" as a warning.
  - `760.cif` parent geometry is distorted (C-N 1.154-1.338 A, all C-H 1.140 A); screen the 884 set for non-physical bond lengths.
  - `_chemical_identity_key` heavy-atom-formula collision; `1102FragCof` atom 72 under-coordinated carbon.

## [2026-09-17] 100-structure HCNO test run and unified layer-planarity guard
- Files touched: `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_hcno100` (first run), `runUniFrag/CoRECOF/test_hcno100_v2` (after the guard).
- Dataset: 100 structures sampled from `CoRE-COF_HCNO884` with seed 20260917, **disjoint** from the original 20 (`selection.txt` in each folder lists them).
- Summary:
  - First run: 100/100 structures processed, 0 tracebacks, Min coverage 89/100, 268 fragment frames, 96.0% of terminal capped sites chemically clean. But 63 of 150 dimer frames had a sub-2.4 A interlayer contact (worst 1.17 A), 48 of them with non-flat layers - the corrugated-framework mode reaching output through Path J, whose dimer rule had no planarity guard.
  - Unified the guard across both dimer rules (see `project-decisions.md` 2026-09-17).
- Validation (same 100 structures, exit 0, 0 tracebacks):
  - Dimer frames 150 -> 74; sub-2.4 A contacts **63 -> 14**; median interlayer contact 2.60 -> 3.45 A; non-flat clashes **48 -> 0**.
  - Min coverage 89/100 -> 88/100. Clean terminal sites ~96% in both runs.
- Follow-up risks / open items:
  - All 14 remaining dimer clashes are capping geometry on flat layers: 100 of 132 contacts are owned by `Ar-CH3` caps, 20 by `=CH2` caps whose local ring system was too small to planarize. A methyl cannot be planarized (it is genuinely sp3). The faithful fix is capping a severed aldimine as `Ar-CH=NH` rather than `Ar-CH3`, preserving sp2 planarity, at the cost of one extra heavy atom and a shift in `_chemical_identity_key` grouping. Affects 8 structures / 14 frames of 74 dimers here. Not started - needs a decision.
  - 12 structures still produce no Min; 46 over-valent and 12 under-valent terminal sites remain, concentrated in a handful of structures (`1169`, `689`, `1189`, `1242`).
  - Cleavage-rule coverage remains hand-enumerated; a "zero severed bonds" diagnostic is still wanted.

## [2026-09-17] Biaryl fallback for all-hydrocarbon COFs (463)
- Files touched: `coffragmentor.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_463_blocks`, `runUniFrag/CoRECOF/test_hcno100_v3`.
- Trigger: user reported `463.cif` (C/H only) yielding no proper node or linker, with reference images of the expected benzene node and biphenyl linker.
- Root cause: `463` is an all-hydrocarbon COF linked by direct aryl-aryl C-C bonds. No cleavage rule matched, so `COF.fragment()` severed nothing and returned 0 nodes and 0 linkers. This is the **third** distinct cause of the same silent fallback (411 = classification tie, 704 = unknown vinylene chemistry, 463 = unknown biaryl chemistry).
- Validation:
  - `463`: node lib `H3 C6`, linker lib `H8 C12`; fragments `463FragCof` C84H60, `Min` C60H44, `OnlyNode` C12H12 (benzene dimer), `OnlyLinker` C24H20 (biphenyl dimer). Every carbon has degree 3, and the interlayer contact is exactly 3.40 A = |T| (previously 1.42-2.31 A).
  - 100-structure regression (exit 0, 0 tracebacks): frames 259 -> 268, structures with fragments 90 -> 91, Min 88 -> 89, sub-2.4 A contacts 14 -> 11, valence quality ~96% unchanged.
- Follow-up risks / open items:
  - Three separate structures have now failed silently through the fallback path for three different reasons. A "zero severed bonds" diagnostic is overdue: it would have surfaced `704` and `463` immediately instead of leaving them to visual review.
  - All 11 remaining dimer clashes are still capping geometry on flat layers: 78 of 118 contacts from `Ar-CH3` caps, 28 from `=CH2` caps in ring systems too small to planarize. Affects 5 structures (`68`, `175`, `502`, `1049`, `1081`). The `Ar-CH=NH` capping change remains the open decision.
  - 11 structures still produce no Min; 46 over-valent / 12 under-valent terminal sites persist.

## [2026-09-17] 543 alkyne capping + stackability, 612 duplicate-key false positives
- Files touched: `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_543_blocks`, `runUniFrag/CoRECOF/test_612_blocks`, `runUniFrag/CoRECOF/test_dedup`.
- Trigger: user reported a redundant H on the two-carbon chain of `543.cif` plus a monomer where a dimer was expected, then that `612.cif` produced no normal, no min and no node - only a linker.
- Root causes (three, see `project-decisions.md` 2026-09-17):
  1. `543`'s two-carbon chain is an **alkyne** (C#C, 1.205 A), not a vinylene. `_cap_open_oxygens` treats a 2-heavy-neighbour carbon with no H as an aromatic edge carbon and caps it, but an alkyne carbon is already fully valent.
  2. `543`'s dimer was rejected by an absolute 0.5 A flatness cutoff at a measured 0.51 A, although its parent is genuinely layered with twisted aryl rings.
  3. `612` was not a fragmentation failure at all: its fragments were discarded as formula-only "duplicates" of an unrelated framework.
- Validation:
  - `543`: every alkyne carbon now carries 0 H (matching the parent), and all four fragments are dimers with interlayer contacts of 3.31-3.58 A.
  - `612`: now emits `FragCof`, `Min`, `OnlyNode`, `OnlyLinker_0`, `OnlyLinker_1`, plus node and linker libraries.
  - Ten-structure dedup test: `612`, `221`, `1032`, `370`, `579`, `141` all retained; `144`/`603`/`651`/`372` still flagged duplicates, confirmed correct - their Weisfeiler-Lehman topology hashes match their partners exactly.
  - MOF isolation proven by identity, not by sampling: `MOFFragmenter._skip_aromatic_carbon_h_cap is BaseFragmenter._skip_aromatic_carbon_h_cap` and `_flush_mof_result` still uses the MOF formula-only key.
- Not re-run: the 100-structure sweep, at the user's instruction - it will be run once all reported problems are fixed.
- Follow-up risks / open items:
  - `612` yields a 4-atom `C2 N2` block as `OnlyLinker_1` (an isolated azine `-CH=N-N=CH-` unit). Chemically it is the linkage itself rather than a strut; worth deciding whether such minimal units should be exported as linkers.
  - The formula-only key is still used for MOFs. It has the same false-positive weakness there, but changing it would alter MOF results and was deliberately left alone.
  - Still open: `Ar-CH=NH` vs `Ar-CH3` capping decision; hand-enumerated cleavage rules with no "zero severed bonds" diagnostic.

## [2026-09-17] Symmetry-aware stacked dimers for staggered (AB) COFs - 585
- Files touched: `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: `runUniFrag/CoRECOF/test_585_blocks`, `runUniFrag/CoRECOF/test_d20_sym` (regression).
- Trigger: user expected `585.cif` to be a dimer, with the second layer diagonally offset rather than directly opposite.
- Root cause: `585` is P6_3/m with two layers per cell (c = 7.00 A). Adjacent layers are related by a 6_3 screw - 60-degree rotation plus c/2 - not by a lattice translation, so the translation-only rule found no axis in the stacking window and produced a monomer. See `project-decisions.md` 2026-09-17.
- Validation:
  - `585`: all four fragments are dimers, each placed by a 60-degree rotation about the fragment axis, separation 3.50 A, lateral offset 0.00 A, closest contact 3.50-3.51 A. Kabsch check confirms the halves are related by exactly 60 degrees with rmsd 0.000.
  - 20-structure regression (exit 0, 0 tracebacks): all 24 dimer placements still resolve to "lattice translation", **zero composition changes**, nothing lost, frames 75 -> 76 (`167FragCofOnlyNode` gained), dimers 33 -> 34, minimum contact 2.55 A, median 3.48 A, zero sub-2.4 A contacts.
- Two regressions caught during development and fixed before reporting:
  - `lattice_reduce` originally minimised total displacement, which collapsed pure lattice translations to zero and removed **every** dimer from the 20-set.
  - The local rotation was first stored as a fragment-centred transform and re-applied verbatim to node and linker blocks, putting their halves 33-38 A apart; `_apply_partner_op` now recentres per block.
- Not re-run: the 100-structure sweep, per the user's instruction.
- Follow-up risks / open items:
  - Unchanged: `Ar-CH=NH` vs `Ar-CH3` capping decision; no "zero severed bonds" diagnostic; `612`'s 4-atom `C2 N2` azine block exported as a linker; MOF duplicate key still formula-only.

## [2026-09-17] Linkage-recognition diagnostic; imine length guard investigated and rejected
- Files touched: `coffragmentor.py`, `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Summary:
  - Added the linkage-recognition diagnostic requested after `704`/`463` failed silently, plus per-bond linkage tags so the report names the real chemistry (see `project-decisions.md` 2026-09-17).
  - Investigated the `612` `C2 N2` linker via a bond-length guard on the imine rule, measured the consequences, and **rejected** the change - it would strip every cut from 15 of 100 structures.
- Validation:
  - Diagnostic output is accurate per structure: `704` -> vinylene, `463` -> biaryl, `211`/`543` -> imine, `1031` -> oxazole, `612` -> imine + vinylene.
  - Tagging is behaviourally inert: 20-structure regression unchanged (76 frames, Min 20/20, 34 dimers, min contact 2.55 A, zero sub-2.4 A, no composition changes, none gained or lost).
- Follow-up risks / open items:
  - A full-884 scan for structures with zero severed bonds was launched to size how much linkage chemistry is still missing; results pending.
  - `Ar-CH=NH` vs `Ar-CH3` capping remains the one open decision. Note that since the dimer contact gate was added, its symptom should change from "clashing dimer" to "no dimer", i.e. lost coverage rather than corrupt geometry - to be confirmed on the next full sweep.

## [2026-09-17] Full 884 linkage-recognition scan
- Scanned all 884 HCNO structures through `COF.fragment()` (4 processes, 0 errors) to count severed bonds and recognised linkage types.
- Results:
  | linkage types recognised | structures |
  |---|---|
  | imine | 727 |
  | biaryl | 53 |
  | imine + vinylene | 27 |
  | vinylene | 25 |
  | oxazole | 8 |
  | dioxin + imine | 7 |
  | dioxin | 4 |
  | oxazole + vinylene | 2 |
  | **none (0 bonds severed)** | **31** |
- The vinylene and biaryl rules added today account for **107 of 884 structures (12%)**. Before today every one of them severed nothing and fell silently to the crude fallback path.
- The 31 still unrecognised are **not** one missing rule. Sampling them (`68` = C12N6, `1013`, `24`, `181`, `1013`) shows mostly fully fused / edge-sharing frameworks with no inter-ring bonds at all - there is no discrete node-linker linkage to cut, so "0 severed bonds" is arguably the correct answer and the node/linker paradigm simply does not apply. A couple (`357`, `900`) do have inter-ring C-N/C-C bonds that no rule claims.
- Consequence for the planned 884 run: expect ~3.5% of structures to be flagged by the new diagnostic. They should be reviewed or excluded rather than trusted, since they will be built by the fallback path.

## [2026-09-17] Ar-CH=NH capping: already implemented; clash cause re-diagnosed
- Files touched: `project-decisions.md`, `project-agent-log.md` (no code change).
- The user approved implementing `Ar-CH=NH` capping. Measuring first showed the code already produces it - 153 aldimine caps against zero imine-derived methyls in the 20-structure set - as a free consequence of the earlier carry-over decision. My earlier premise ("the carbon side becomes `Ar-CH3`") was wrong, inferred from contact-owner statistics without checking what those carbons were.
- Re-diagnosed the residual dimer clashes: they are inherited from parent CIF hydrogen placement, not produced by UniFrag. 147 of 567 layered parents (26%) already contain a sub-2.0 A non-bonded interlayer contact; heavy-atom contacts are healthy (median 3.27 A). See `project-decisions.md` 2026-09-17.
- Correction recorded for future agents: the first parent scan used the raw minimum interatomic distance and reported 333/567 "clashing" parents. That metric is wrong - a contact at bond length is a covalent bond spanning the cell boundary, not a clash. Excluding pairs inside a covalent cutoff gives the defensible figure of 147.
- Open decision for the user: leave clashing-parent structures to come out as monomers (current, safe), or re-place parent hydrogens (better geometry, but conflicts with "keep every structure original").

## [2026-09-17] 760 nitrile mis-cut fixed; 167 under-capped N diagnosed
- Files touched: `coffragmentor.py`, `project-decisions.md`, `project-agent-log.md`
- Test folders: scratchpad `x_167`, `x_760`; `runUniFrag/CoRECOF/test_d20_nitrile` (regression).
- `760`: fixed. The imine rule was severing the structure's nitriles, which decomposed nothing but suppressed the biaryl fallback. Requiring the nitrogen to have >= 2 heavy neighbours resolves it (see `project-decisions.md` 2026-09-17).
  - Now: node lib `H6 C16`, linker lib `H16 C28 N4`; `FragCof` C128H82N16, `Min` C62H40N4, `OnlyNode` C16H10, `OnlyLinker` C28H20N4 - matching the user's reference images.
  - Nitrile nitrogens are now retained (`C102H72N4` -> `C128H82N16`).
  - 20-structure regression: exit 0, no warnings, Min 20/20, dimers 34, min contact 2.55 A, zero sub-2.4 A. Only `760` changed; it gained `OnlyNode` and `OnlyLinker`. Nothing else touched.
- `167`: diagnosed, not fixed. The under-capped nitrogen is `fix_odd_electron_multiplicity` stripping an aldimine's only H. The parity is forced by composition (`h + n` odd) and the fragment offers no benign site - all three removal candidates are aldimines and every nitrogen is already valence-complete. Three options recorded in `project-decisions.md`; awaiting a user decision.

## [2026-09-18] 167 uncapped nitrogen fixed by moving parity repair after dimerisation
- Files touched: `fragmentation_oop.py`, `project-decisions.md`, `project-agent-log.md`
- Trigger: user reported an uncapped N in `167.cif` and asked for multiplicity = 1, suggesting adding a hydrogen rather than removing one.
- Diagnosis: the cap was right; the parity repair was checking the monomer before the second layer was added. `167`'s monomer is odd (Z_sum 369) but its dimer is even (738), so the repair was stripping an aldimine's hydrogen for no reason. RDKit confirmed the emitted `167FragCof` could not be assigned a valence at that atom, while the node and linker fragments were already clean.
- Changes (see `project-decisions.md` 2026-09-18): parity repair moved after dimer construction in both COF fragment paths, and addition preferred over removal when a repair is genuinely needed - both COF-scoped via hooks, MOF untouched.
- Validation:
  - `167`: no parity surgery at all now. All five fragments even, RDKit-valid at charge 0 with zero radical electrons. `167FragCof` C84H58N18O6 -> C84H60N18O6, `Min` C60H42N14O6 -> C60H44N14O6.
  - 20-structure regression: exit 0, no tracebacks, QM-Fix events 7 -> 1, RDKit-valid closed-shell fragments 75 -> 77 of 78, Min 20/20, dimers 34, min contact 2.55 A, zero clashes, nothing gained or lost.
- Follow-up: `142FragCofMin` remains the single fragment needing a sacrificed cap - a true monomer whose composition admits no closed-shell neutral capping.

## [2026-09-18] 100-structure re-run after the full day's fixes (test_hcno100_v5)
- Same 100 structures as before (seed 20260917, disjoint from the original 20). Exit 0, no tracebacks.
- | metric | baseline (start of day) | v4 (pre nitrile-cap fix) | **v5 final** |
  |---|---|---|---|
  | fragment frames | 268 | 311 | **312** |
  | Min coverage | 89/100 | 94/100 | **94/100** |
  | structures yielding fragments | 91 | 93 | **93** |
  | flagged duplicate | 14 | 8 | **8** |
  | terminal-site valence correct | 96.0% | 94.2% | **99.2%** |
  | over-valent sites | 46 | 89 | **4** |
  | under-valent sites | 12 | 8 | 10 |
  | odd-electron fragments | 4 | 3 | **3** |
  | stacked dimers | 75 | 113 | **115** |
  | dimers with a sub-2.4 A contact | 11 | 2 | **3** |
- v5 fragment kinds: Frag 92, Min 87, OnlyLinker 79, OnlyNode 54.
- Linkage recognition across the 100: imine 77, biaryl 9, vinylene 6, imine+vinylene 3, oxazole+vinylene 1, dioxin 1, **unrecognised 3** (`68`, `363`, `1112` - all flagged by the new diagnostic). Fallback paths used only 12 times in total.
- Residual defects are concentrated, not diffuse: the 4 over-valent sites sit in `1112` and `526`; the 3 clashing dimers are `175FragCofOnlyLinker_0` (1.64 A), `526FragCofOnlyNode` (1.63 A) and `689FragCofOnlyNode` (1.96 A). `1112` is one of the three structures whose linkage chemistry is unrecognised, so its problems come from the fallback path, as expected.
- Follow-up: two further disjoint 100-structure sets (`test_hcno100_B`, `test_hcno100_C`, seed 20260918) were selected and queued at the user's request, bringing planned coverage to 320 of 884.

## [2026-09-18] Three independent 100-structure validation runs (300 of 884)
- Sets: `test_hcno100_v5` (seed 20260917) and `test_hcno100_B` / `test_hcno100_C` (seed 20260918), all mutually disjoint and disjoint from the original 20. B and C were run in parallel (4 processes each across 8 cores) after the sequential chain proved too slow. All three: exit 0, zero tracebacks.
- | set | frames | Min | yielding | dup | valence ok | OVER | under | odd | dimers | clashing |
  |---|---|---|---|---|---|---|---|---|---|---|
  | v5 | 312 | 94/100 | 93 | 8 | 99.2% | 4 | 10 | 3 | 115 | 3 |
  | B | 308 | 91/100 | 97 | 4 | 98.4% | 0 | 26 | 2 | 121 | 3 |
  | C | 286 | 87/100 | 91 | 10 | 98.8% | 7 | 10 | 6 | 95 | 2 |
  | **total** | **906** | **272/300** | **281/300** | 22 | **98.8%** | 11 | 46 | 11 | **331** | **8** |
- Linkage recognition across the 300: imine 238, biaryl 29, vinylene 13, imine+vinylene 4, dioxin+imine 2, dioxin 3, oxazole 1, oxazole+vinylene 1, **unrecognised 9 (3.0%)**. That 3.0% matches the full-884 scan's 31/884 = 3.5%, so the sample is representative.
- The vinylene and biaryl rules added today account for **47 of the 300 structures (~16%)**; before today every one of them severed nothing and fell silently to the crude fallback path.
- Residual defects across 906 fragments: 8 clashing dimers, 11 over-valent and 46 under-valent terminal sites, 11 odd-electron fragments. Concentrated in the structures flagged as unrecognised or with distorted parent CIFs, not spread across the set.

## [2026-09-19] Analysis of the 300-structure COF fragment output; capping-H collision found
- Files touched: `project-agent-log.md` only. **No source change was made** - the fix below was proven with a runtime wrapper in the scratchpad, not committed.
- Scope: pooled re-analysis of `test_hcno100_v5`, `test_hcno100_B`, `test_hcno100_C` (300 structures, 906 fragments) straight from the extxyz collections, independent of the numbers recorded on 2026-09-18.
- Reproduced exactly: 906 frames, 272/300 Min, 281/300 yielding fragments, 22 duplicates, 0 CSV errors, 9 unrecognised linkages (3.0%), 8 dimers with a sub-2.4 A contact. The logged figures hold.
- New measurements not previously recorded:
  - Fragment sizes: Frag median 168 atoms (max 800), Min median 110, OnlyNode median 54, OnlyLinker median 41.
  - Connectivity: with a `r_i + r_j + 0.40 A` bond criterion **no fragment is split into more than two pieces**. 575 monomers, 331 two-piece.
  - Of the 331 two-piece fragments, 309 are genuine stacked dimers (2.4-4.0 A), 8 clash (<2.4 A) and **14 are not stacked at all**. Four of those are badly wrong - `481FragCofOnlyLinker` (31.7 A), `1001FragCofOnlyNode` (26.3 A), `1001FragCofOnlyLinker` (19.1 A), `324FragCofOnlyLinker` (16.6 A) - two disjoint blocks emitted as one frame. The other ten (`589`, `531`, `908`) sit at 4.25-4.43 A, wide but plausible.
  - Chemistry diversity: 118 distinct OnlyNode formulas from 165 fragments, 155 distinct OnlyLinker from 209, 237 distinct Min from 254.
  - Per-structure completeness: only 95 of 281 structures yield the full set (Frag+Min+OnlyNode+OnlyLinker); 88 give a linker but no node, 40 a node but no linker.
- **Defect found: capping hydrogens are repositioned into existing atoms.** 43 of 906 fragments contain an atom pair closer than 0.9 A, across 21 structures. Eleven of those structures inherit the overlap from a disordered parent CIF (`1116`, `898`, `502`, `284`, `209`, `134`, `1155`, `1189`, `1091`, `1043`, `1057`), but **ten are introduced by UniFrag** (`819`, `818`, `965`, `903`, `1159`, `526`, `1112`, `173`, `364`, `992`) - their parent CIFs have no sub-0.9 A pair.
  - Mechanism, confirmed by tracing: `place_capping_h` (`fragmentation_oop.py:214`) does honour its `min_hh=1.5` / `min_heavy=0.9` clearance test, but `enforce_sp2_capped_h_geometry` afterwards rewrites the cap position unconditionally at `fragmentation_oop.py:457` (`coords[hidx] = p + bl * dvec`) with no clearance test at all. It rotates the cap into the sp2 slot that points straight at a bonded neighbour. `_planarize_conjugated_caps` (`fragmentation_oop.py:3469`) and `optimize_capped_h_geometry_only` have the same unguarded write; `enforce_capped_oh_geometry` is the only one of the four that re-scores before committing.
  - Typical result: a cap 1.09 A from its parent carbon and **0.42 A from the adjacent carbon** (`903FragCof`, four times over), or two caps on neighbouring carbons planarized into each other at 0.61 A H...H (`526FragCof`). These fragments cannot be run in QM.
  - Validation of the cause: a 6-structure probe (`903`, `819`, `1159`, `526`, `965`, `1112`) reproduces 23 overlapping pairs deterministically. Re-running with a wrapper that reverts any cap move landing within 0.90 A of a heavy atom or 1.20 A of another H gives **23 -> 0 overlaps, 19/19 fragments retained, zero composition changes, no heavy atom displaced**. Scratchpad only; `fragmentation_oop.py` is untouched.
- Gap in the existing QC: the run-level metric is "terminal capped site valence", which reported 98.8% correct across the 300. It does not look at cap-to-atom distances, so all 43 of these fragments passed. The same is true of `cof-fragment-checklist.md`, which has no minimum-interatomic-distance check.
- Follow-up risks / open items:
  - Decide whether to add the clearance test to the three unguarded cap-repositioning methods. Reverting to the pre-planarization position is the conservative option and was the one measured.
  - Add a minimum-interatomic-distance gate (suggest 0.9 A any pair, 1.2 A H...H) to the checklist and to the run-level summary.
  - The four fragments with 16-32 A gaps between two disjoint blocks need a separate diagnosis; they are not stacking failures.
  - Eleven structures carry overlapping atoms in the parent CIF itself. These should be screened out of the input set rather than repaired downstream.

## [2026-09-20] COF fragment QA: checklist, six defect fixes, guest removal
- Files touched: `coffragmentor.py`, `fragmentation_oop.py`, `cof-fragment-checklist.md` (new), `runUniFrag/check_cof_fragments.py` (new), `runUniFrag/remove_guest_molecules.py` (new).
- **NOTE**: the 2026-09-19 entry above was found uncommitted in the working tree at the start of this session; it is committed here unchanged. Its finding is unaddressed and still open - `enforce_sp2_capped_h_geometry` and `_planarize_conjugated_caps` reposition capping H with no clearance test.
- Added `cof-fragment-checklist.md`, a standard acceptance checklist, and `check_cof_fragments.py`, which scores a run against it and exits non-zero on any MUST failure. The checker reuses UniFrag's own `is_valid_bond`/`_terminal_bond_order` so audit and fragmenter cannot disagree, and reports only defects that need no aromatic perception.
- Defects found and fixed, each measured on the full 884:
  - Parallel race in `_prune_duplicate_cof_helper_files` (unguarded `unlink`) aborted **329 of 884 structures while the job still exited 0**. Now guarded.
  - Orphaned bridging units (bare `[N]`, `[N,N,N,H]`) - orphan guard restores one cut per component under 4 heavy atoms.
  - Guest molecules exported as linkers - components with zero severed bonds are skipped.
  - Superimposed atoms, 42 of 1695 frames, 0.15-0.75 A apart - COF-only `_dedupe_superimposed_atoms` hook.
  - Odd-electron fragments 39 -> 5 by widening the parity-repair candidate set (carbon admitted, geometric filter dropped) with a validated, clearance-checked placement. Survivors are quarantined to `fragments_quarantine.extxyz`, never shipped in the main collection.
  - Non-covalent H contacts perceived as bonds (526: ketoenamine H 1.09 A from C, 1.28 A from keto O) left the framework one component, so 6 severed linkages produced no blocks and it fell silently to Path A. Hydrogen is now monovalent, keeping its nearest HEAVY neighbour.
  - Crystal-symmetry partner ops applied verbatim to helper blocks threw 481's halves 30-45 A apart. All ops now returned block-local; helper dimers also gained the contact test the assembled fragment always had.
- **Reverted**: a secondary-amine linkage rule. Every aggregate metric improved and it was still wrong - on 1180 it emitted an 11-atom `HN-CH2-C(=O)-CH3` scrap and destroyed a correct two-node decomposition. Two guards were tried and neither separated the good case from the bad.
- Added `remove_guest_molecules.py` (periodicity-based, CSD cross-checked 10/11 exact, 0 false negatives in 40). On the 884: 11 structures, 436 guest atoms, 1180 alone holding 28 acetone. **Measured benefit on the full set: zero** - UniFrag already discards guests via the zero-cut-component skip.
- MOF/bio isolation: all shared-code changes are class attributes or hooks defaulting to previous behaviour; verified `MOFFragmenter` and `MacromolFragmenter` resolve to the base implementation by identity. **Not verified by running MOFs** - no MOF set was fragmented this session.
- Follow-up risks / open items:
  - The 2026-09-19 unguarded cap-repositioning finding, still open; 6 structures have a demonstrably self-inflicted sub-1.0 A H...H contact (364: parent 2.45 A -> fragment 0.90 A).
  - 66 over-coordinated fragments remain, largely inherited (284's parent has 80 of its own).
  - 40 structures with no recognised linkage; the `N(CC)` no-hydrogen bridge (1061, 1128, 652) is an uncharacterised family.
  - **Process lesson**: four defects in a row were caught by the user viewing fragments, not by any metric. A fragment that is the wrong molecule, or in two distant pieces, passes every valence/parity/contact check. Confirm linkage-chemistry and geometry changes per structure, visually.

## [2026-09-21] Full 884-structure production run analysed (final_smp, cluster)
- No source change. Analysed a completed cluster run at `arf:/arf/scratch/otayfuroglu/uniFrag_works/runUniFrag/coreCOFs/final_smp` (SLURM job 6392467, 20 cores) - the full CoRE-COF HCNO884 set through the current `fragmentation_oop.py`/`coffragmentor.py` (post the 2026-09-20 fix session: superimposed-atom dedup and quarantine both active in this run's log).
- Run-level: 884/884 CIFs, 0 CSV errors, 0 tracebacks. Normal fragment: 630 chemically unique + 254 flagged duplicate of another structure (28.7% dedup rate) = 884. Min fragment: 587 unique + 245 duplicate + 52 none producible.
- Collection: 1726 frames across 646 structures (the 630 unique-normal plus a handful whose node/linker survived while the parent was flagged duplicate). Kinds: Frag 629, Min 585, OnlyLinker 314, OnlyNode 198. Sizes (atoms): Frag median 174 (max 1697), Min median 116 (max 556), OnlyNode median 54 (max 332), OnlyLinker median 52 (max 293).
- Linkage recognition (log): imine 689, biaryl 85, vinylene 40, oxazole 8, imine+vinylene 7, dioxin 6, boroxine/imine 3, oxazole+vinylene 2, **unrecognised 44 (5.0%)** - higher than the 3.5% baseline on the same 884 measured 2026-09-17; not yet reconciled (possibly the superimposed-atom and guest-removal fixes changed which structures reach the linkage step).
- Dimers: 597/1726 fragments (34.6%) are two-component; 583 are genuine stacked dimers at 2.4-4.0 A (median 3.44 A), 6 clash under 2.4 A, 8 sit loose at 4.0-4.4 A. **219/646 structures (33.9%) carry at least one dimer fragment**; 214 of those are the full `Frag` itself.
- QM-Fix repairs (log): carbon 268, superimposed 95, secondary-amine 70, primary-amine 19, amide 11, alcohol 7, phenol 1, O-on-N 1. 5 open-shell fragments quarantined (`1005`, `1104` x2, `615`, `744`), none leaked into the main collection.
- Geometry QC (own pass, `r_i+r_j+0.40 A` bond test): 109/1726 fragments over-coordinated (894 atoms, 0 under-coordinated), 45/1726 contain a sub-0.9 A atom pair across 20 structures. Traced against parent CIFs: of those 20, **17 inherit the overlap from the parent** (`284`, `299`, `333`, `334`, `502`, `506`, `55`, `657`, `718`, `861`, `898`, `929`, `930`, `1043`, `1131`, `1189`, `1220`) and **3 are self-inflicted with a clean parent** (`364`, `532`, `1034`) - the still-open unguarded-cap-reposition defect from 2026-09-19/09-20, now confirmed present in the full production run at a smaller rate than the earlier 100-sample probes.
- Not done: no per-structure visual confirmation of linkage chemistry this session (the 2026-09-20 process lesson) - this was a metrics-only pass to answer the user's request for fragment-count/size/dimer statistics, not a QA pass.
- Follow-up risks / open items:
  - Reconcile the 44 vs 31 unrecognised-linkage count against the 2026-09-17 baseline on the same 884 structures.
  - The unguarded cap-reposition defect (`enforce_sp2_capped_h_geometry`, `_planarize_conjugated_caps`) is still unfixed and still produces self-inflicted sub-0.9 A contacts (3 confirmed this run: `364`, `532`, `1034`).
  - 646/884 structures reach the fragment collection; the other 238 either duplicate another structure's chemistry (254 flagged, but some structures may have both a duplicate normal and a surviving node/linker) or produce nothing - worth a clean breakdown by reason (duplicate vs unrecognised-linkage vs other) next time.

## [2026-09-21] QM-cost subset built from final_smp (250-atom cascade)
- No source change. Built a reduced-size subset from the `final_smp` production collection (see the entry above) for downstream QM jobs, at the user's request.
- Rule applied per structure, in order: (1) keep `Frag` (normal) if <=250 atoms; (2) else replace with `Min` if `Min` <=250 atoms; (3) else fall back to all `OnlyNode`/`OnlyLinker` frames for that structure (as separate blocks); (4) if none of those exist or fit, exclude the structure and report it.
- Output: `fragments_subset_qm250.extxyz` (647 frames) + `subset_selection.csv` (646 rows, one per structure, recording which category/fragment/atom-count was used). Both written to the scratchpad and copied back to `arf:.../coreCOFs/final_smp/cifs/` alongside the source collection.
- Breakdown of the 646 structures present in `final_smp`'s collection:
  | category | structures | detail |
  |---|---|---|
  | kept normal (Frag <=250) | 459 | atoms 15-250, median 148 |
  | replaced by Min (Frag >250, Min <=250) | 139 | Frag median 336 -> Min median 179; total 22,524 atoms saved |
  | node+linker fallback (Frag >250 AND Min >250/missing) | 22 | 38 node/linker frames emitted |
  | node+linker only, no Frag/Min ever produced | 6 | |
  | Min only, no Frag ever produced, Min <=250 | 11 | sizes 46-182 |
  | **excluded - no fragment fits under 250 and no node/linker fallback exists** | **9** | `14` (Frag 1697), `68`, `79`, `183`, `218`, `481`, `664`, `1173`, `1174` |
- Residual: 3 of the 38 fallback node/linker frames still exceed 250 atoms (332, 293, 252) - included anyway as the smallest available option, per the user's rule having no further fallback specified.
- Final subset: 647 frames, size range 15-332 atoms (median 157), only 3 frames >250. 637/646 structures represented; 9 excluded outright (listed above).
- Follow-up: the 9 excluded structures and the 3 still-oversized fallback frames were not further investigated - flag if the user wants them force-truncated or handled differently.

## [2026-09-21] Subset trimmed with a +5 atom tolerance; duplicate check
- No source change. Two follow-ups to the 250-atom subset above, both applied to `fragments_subset_qm250.extxyz` (scratchpad + copied back to `arf:.../final_smp/cifs/`).
- **Tolerance rule applied to the node/linker fallback tier**: of the 3 fallback frames still over 250 atoms (332, 293, 252), dropped the two furthest over and kept the one within +5 (<=255): `931FragCofOnlyNode` (332, dropped - `931` keeps only its 15-atom `OnlyLinker`), `1176FragCofOnlyLinker` (293, dropped - `1176` keeps only its 65-atom `OnlyNode`), `1107FragCofOnlyLinker` (252, kept). Subset: 647 -> **645 frames**. `subset_selection.csv` updated to record both drops.
- **Duplicate check, at the user's request**: verified with the project's own dedup definition (`COFFragmenter._chemical_identity_key` - heavy-atom formula + Weisfeiler-Lehman bond-graph topology hash), imported directly from the cluster's copy of `fragmentation_oop.py` (`arf:/arf/home/otayfuroglu/UniFrag/fragmentation_oop.py`, confirmed byte-identical to the local repo copy via `diff`). **Result: 0 duplicate groups, 645/645 distinct topology keys, 645/645 distinct frame labels.** 37 pairs/triples share an exact full formula (incl. H) but were confirmed to have different topology hashes - formula coincidence, not duplication.
- Note: this check is against the topology-aware key already established and validated in this project (2026-09-17 entry); it was not re-derived or re-validated this session, only re-applied.

## [2026-09-21] SOAP-based fragment-vs-parent fidelity analysis
- No source change. Analysed how well the 645-frame `fragments_subset_qm250.extxyz` (see the two entries above) represents each parent COF's true local atomic environment, at the user's request, using SOAP (Smooth Overlap of Atomic Positions, DScribe) at r_cut in {3, 4, 5, 6} A.
- **Method**: for each fragment atom, found its literal periodic-image match in the parent CIF (position match, tolerance 0.05 A after wrapping to the nearest lattice image); unmatched atoms are UniFrag's added capping H and excluded from the fidelity metric (kept as a separate count). For each matched atom, computed a SOAP local descriptor from the finite fragment (no PBC) and from the true periodic parent (pbc=True, DScribe's own periodic neighbour search) at the SAME physical position, then compared with cosine similarity. SOAP params: species [H,C,N,O,B], n_max=6, l_max=4, sigma=0.5.
- **Bug caught before results could be trusted**: ASE's own CIF reader and pymatgen (which `coffragmentor.py` uses internally) build DIFFERENT Cartesian frames for any non-orthogonal cell (a-along-x vs c-along-z convention) - same fractional coordinates, same cell metric, but rotated relative to each other. Reading parents with ASE directly gave 0% atom-matching for every skewed-cell structure (~30% of the pilot sample) with no error - a silent, dataset-wide correctness bug had this gone unnoticed. Fixed by reading parents through pymatgen and building the comparison `ase.Atoms` from its cell/coords directly (`read_parent()` in the scratchpad's `soap_lib.py`). Verified: `SOAP(periodic=True)` on a plain ASE-read periodic structure reproduces a manually-built supercell's SOAP vector exactly (cosine sim = 1.0000000000000002), confirming DScribe's own periodic handling is correct once the input frame is.
- **Coverage**: 603/645 fragments (93.5%) reached a reliable atom correspondence (>=80% of atoms matched to the parent) and were scored; 42 skipped (low match fraction, mostly 40-75%) and reported separately, not folded into the aggregate - root cause not chased further this session.
- **Headline result - fidelity degrades with cutoff, as expected**: mean cosine similarity 0.972 (r=3) -> 0.945 (r=4) -> 0.892 (r=5) -> 0.852 (r=6). Atoms within 2 A of a capping site score far lower (0.81 @ r=3, 0.68 @ r=6) than atoms >8 A from any cut (0.99 @ r=3, 0.88 @ r=6) - a clean, monotonic boundary effect.
- **By atom type**: N is consistently the least faithfully represented element (0.889 @ r=3 -> 0.783 @ r=6), C/H/O track together and above N at every cutoff - consistent with N sitting disproportionately at imine linkage/cut sites (imine is 78% of recognised linkages in this set).
- **Key actionable finding - monomers of tightly-stacked frameworks are a poor local-environment representation at larger cutoffs, dimers largely fix it**: for structures whose parent's shortest cell axis (interlayer spacing) is under 6 A - i.e. smaller than the SOAP cutoff itself - monomer fragments score only 0.76 mean fidelity at r_cut=6, vs 0.91 for the matching dimer fragments of equally tight structures, and vs 0.87 for monomers of loosely-stacked (>=8 A) frameworks. The worst-scoring fragments overall (`662FragCof` 0.14, `187FragCof` 0.19, `28FragCof` 0.24 at r_cut=6) are exactly this case: monomers with a 3.4-4.35 A stacking period, where a 6 A SOAP sphere reaches into empty space a real dimer would have filled with the neighbouring layer.
- Artifacts (all in scratchpad, not committed): `soap_lib.py` (matching + parent-reading helpers), `full_run.py` (the scoring pass), `soap_results.pkl`/`soap_df.pkl`/`soap_fragdf.pkl` (raw + aggregated results), `soap_fidelity_summary.png` (3-panel figure: fidelity vs cutoff by element, fidelity vs distance-to-cut, monomer-vs-dimer by stacking tightness). Sent to the user, not pushed anywhere.
- Follow-up risks / open items:
  - The 42 skipped (low-match) fragments were not root-caused; worth checking whether they're a distinct failure mode (further skew case, disorder) or just need a larger position-matching tolerance.
  - The ASE-vs-pymatgen Cartesian-convention mismatch found here is specific to THIS analysis script, not UniFrag's own pipeline (which is internally pymatgen-consistent throughout) - no action needed in `fragmentation_oop.py`/`coffragmentor.py`, noted here only so a future from-scratch analysis doesn't rediscover it the slow way.
  - Given the monomer/tight-stacking finding, worth considering whether the QM subset's monomer selections for min-cell<6A structures should be reconsidered in favour of their dimer, if available - not acted on this session, purely diagnostic.

## [2026-09-22] Checked openMLP preflight-QM run on subset_qm250; found a gap-placeholder bug
- No source change (openMLP lives outside this repo, at `arf:/arf/home/otayfuroglu/openMLP`). Checked the user's active run at `arf:/arf/scratch/otayfuroglu/deepCOF_works/coreCOFs/HCNO/runOpenMLP/subset_qm250` - the bootstrap QM stage (ORCA, PBE-D4/def2-TZVP ENGRAD) of an 8-cycle active-learning pipeline over the 645-fragment `fragments_subset_qm250.extxyz` built in the sessions above.
- **Status**: 626/645 (97.1%) fragments QM-complete and `accepted`, 19 pending (8 ORCA jobs still running, ~14-15h elapsed, `squeue` job ids 6392952-59), 0 failures. Already above the pipeline's own `min_success_fraction: 0.9` bootstrap gate regardless of how the stragglers finish.
- **Bug found**: 39/626 fragments (6.2%) report a HOMO-LUMO gap of 2.7211-2.7212 eV to 4 decimal places - unrelated structures agreeing to within 0.0001 eV is not physics; this is ~exactly 0.1 Hartree and almost certainly a hard-coded fallback value in openMLP's gap-extraction code, silently substituted when it can't read the real orbital energies. All 39 are still marked `status=accepted`, so this QC gate doesn't catch it. These same 39 cluster with elevated SCF iteration counts (100-800 vs a dataset median of 19) and, often, elevated forces - consistent with genuine SCF convergence trouble being masked rather than reported.
- **Force outliers not caught by config**: force_max up to 26.6 eV/Å (median 3.4, p75 5.6). The pipeline's own `qm.max_force_abs: 50.0` filter would not remove any of these, so they will flow into bootstrap training as-is unless someone intervenes.
- **Cross-checked against UniFrag's own geometry QC** (the sub-0.9-A-overlap / over-valence flags recorded in the two `[2026-09-21]` entries above, computed on the full 1726-fragment collection): of the 20 worst-force fragments, 7 were already flagged (`1242`, `1131`, `284`, `554`, `506`, `1240`, plus one more) - confirming part of this traces to known disordered-parent geometry, not a new defect. Of 108 fragments with gap<0.1 eV, 26 were already flagged. The remaining ~65-75% in both cases is new information the geometry-only QC could not have surfaced.
- **Other QM-validation stats** (n=626): force_rms median 0.83, max 7.01 eV/Å; wall time median 15.2 min, max 897 min (256FragCof, 1464 SCF iterations); homo_lumo_gap min -0.532 eV, 132/626 (21%) under 0.3 eV (some of this is expected PBE gap underestimation on extended conjugated aromatics, not necessarily a defect on its own).
- Follow-up risks / open items:
  - Decide whether to exclude/re-run the 39 placeholder-gap fragments (and especially the ones also flagged by geometry QC) before they enter bootstrap training - not acted on, the user's call.
  - The placeholder-gap bug itself lives in openMLP, not this repo; flagging here only because it was found while checking a UniFrag-produced dataset. No fix attempted.
  - 19 fragments still pending; final bootstrap completion not yet confirmed as of this check.

## [2026-09-22] Checked openMLP preflight-QM run on subset_qm250 (arf cluster)
- No source change in this repo. Checked results at `arf:/arf/scratch/otayfuroglu/deepCOF_works/coreCOFs/HCNO/runOpenMLP/subset_qm250/cycle_runs_pbed4/preflight_qm/` - the ORCA single-point energy+forces pass over the 645-frame `fragments_subset_qm250.extxyz` subset, run via a separate tool (`/arf/home/otayfuroglu/openMLP`), at the user's request.
- **Status**: 628/645 (97.4%) completed, 17 still running (8 SLURM jobs live on `barbun` as of the check). Per-fragment wall time: median 15.3 min, max 896.8 min (`256FragCof`, ~15 h) - driven by size and, for a handful, very slow SCF convergence (up to 1903 electronic iterations vs a median of ~19; 62/628 needed >100).
- **All 628 completed fragments show `status=accepted`, zero rejected, zero quarantined.**
- **Bug found in openMLP itself** (not this repo): `openmlp/qm_calc/qm_validation.py:241`, `limit = float(limit or 0)` inside the warning-threshold loop. Because the run's `min_homo_lumo_gap_ev` is configured as exactly `0.0`, Python's `x or 0` treats that literal zero as falsy, so `limit` stays `0` and the guard `if limit > 0` on the next line never fires - the HOMO-LUMO-gap warning check is silently a no-op whenever its threshold is exactly 0, regardless of the actual gap. Confirmed by reading `validate_orca_result()` end to end; the same trap would hit `max_net_force_ev_per_angstrom`/`max_abs_force_ev_per_angstrom` too if either were ever set to literal `0`, though in this run they're 0.10/100.0 so only the gap check is affected.
- **Consequence, measured against `qm_validation_fragments_subset_qm250.csv`**: 55/628 (8.8%) completed fragments have a negative HOMO-LUMO gap (as low as **-0.532 eV**, `521FragCof`) and were written into the training-bound engrad collection instead of the `quarantine_*` file the config's `warning_action: quarantine` calls for. About 10-15 of the 55 are clearly negative (< -0.05 eV); the rest are within a few meV of zero (near-degenerate frontier orbitals in large conjugated systems), plausible on their own merits but still exactly what the configured rule was written to catch.
- **Force magnitudes**: force_rms up to 7.0 eV/A, force_max up to 26.6 eV/A (median 0.83 / 3.45 eV/A) - large, but under both this run's warning threshold (100 eV/A) and the post-QM hard filter (`max_force_abs: 50.0`), so unaffected by the bug above; not flagged by anything in the current config. Worth the user's own judgement on whether 50 eV/A is the right cut for unrelaxed, crystal-extracted starting geometries versus deliberate off-equilibrium sampling.
- **Cross-checked against this repo's own fragment-geometry QC** (2026-09-21 entries): a handful of the worst force/gap outliers were already flagged there for overlapping atoms or over-coordination (`1131FragCof` overlap=11, `284FragCof` over=32, `818FragCof` over=12, `1091FragCofMin` over=10, `506FragCof` overlap=12) - consistent, since clashing atoms produce large repulsive DFT forces. But most of the 55 negative-gap and high-force fragments have no such flag, so this looks mostly like a genuine electronic-structure/off-minimum-geometry effect in large conjugated fragments, not primarily a re-appearance of the known capping defect.
- Not investigated: the 17 still-pending fragments include `662FragCof` and `187FragCof`, the two worst outliers from the 2026-09-21 SOAP fidelity analysis (monomers of very tightly-stacked, short-c parents) - plausibly slow for the same underlying reason (unusual local electronic structure), but not confirmed; noted as a hypothesis only.
- Follow-up risks / open items (all in openMLP, not this repo):
  - Fix `qm_validation.py:241` (e.g. `limit = float(limit) if limit is not None else 0.0`) and re-run validation on the existing ORCA outputs (no need to redo the QM itself) to correctly quarantine the 55 negative-gap fragments before they reach training.
  - Decide whether the ~10-15 clearly-negative-gap fragments should be excluded outright rather than quarantined-and-revisited.
  - No action taken here - this is the user's separate tool; reported for their decision.

## [2026-09-23] QM preflight duration analysis (subset_qm250, PBE-D4 ORCA)
- No source change. Analysed `preflight_qm_report.csv` from `arf:/arf/scratch/otayfuroglu/deepCOF_works/coreCOFs/HCNO/runOpenMLP/subset_qm250/cycle_runs_pbed4/preflight_qm/`, at the user's request - the QM preflight run over the 645-fragment `fragments_subset_qm250.extxyz` built in this project two sessions ago (confirmed byte-identical via md5sum between the two locations).
- **Status**: 638/645 (98.9%) success, 7 pending (still queued/running - 5 `orca_cal` SLURM jobs were active at check time, some >24h wall-clock already), 0 explicit errors.
- **Duration distribution is heavily right-skewed**: median 15.4 min, mean 43.6 min (stdev 132.6 min), p90 49.9 min, p99 840.8 min, max **1278.5 min (21.3 hours)** for `1223FragCof`. Total serial core-time across the 638 completed: **464 hours (19.3 days)**; ~23.2 hours if spread over ~20 parallel workers.
- 5.3% of fragments (34/638) take over 2 hours; 11.6% finish in under 5 minutes.
- **Size correlation is weak on raw atom count (Pearson r=0.22) but strong on log(duration) (r=0.665)** - cost is not simply proportional to size; a handful of specific frameworks are disproportionately expensive to converge regardless of the 250-atom cap.
- **By selection category** (from `subset_selection.csv`, see the "QM-cost subset" entry two sessions back): `replaced_by_min` fragments cost the most on average (mean 72.4 min, median 21.6 min) despite being capped at <=250 atoms same as the rest - expected, since this category is selected *because* its parent `Frag` was originally the largest/most complex (median 336 atoms before reduction), so the underlying electronic structure stays harder even after trimming. `kept_normal` mean 38.5 min/median 14.3 min; `node_linker_fallback`/`node_linker_from_start` (small separate blocks) are cheapest, means under 11 min.
- **All 7 pending fragments sit at the top of the size range (204-250 atoms)**, and 4 of the 5 slowest completed runs are `replaced_by_min` - consistent with the same "originally-largest-framework" effect, not a new finding.
- Noted but not confirmed: 2 of the 7 pending fragments (`187FragCof`, `662FragCof`) are also 2 of the 3 worst-scoring fragments in the SOAP fidelity analysis from two sessions ago (monomers of tightly-stacked frameworks, sim@6A of 0.19 and 0.14) - flagged as a coincidence worth watching, not established as causal (n=2).
- Artifacts (scratchpad, not committed/pushed): `preflight_qm_report.csv` (pulled copy), `qm_duration_summary.png` (duration histogram + size/category scatter). Sent to the user.
- Follow-up: the 7 pending jobs should be re-checked once the active SLURM jobs finish; if they end up failing or taking much longer than the 21.3h ceiling seen so far, that's worth a dedicated look given 2 of them coincide with the known-poor-fidelity monomers.

## [2026-09-23] QM preflight re-check: 6 of 7 pending finished, 1 genuine SCF failure
- No source change. Re-checked `preflight_qm_report.csv` at the user's request. The CSV file itself was byte-identical to two days ago (md5 unchanged) despite `squeue` now showing zero running jobs for the user - **the report is a static snapshot, not regenerated on job completion**. Went to the per-fragment ORCA work directories (`.../preflight_qm/qm/run_engrad_fragments_subset_qm250/<label>/`) directly to get the true current state.
- Of the 7 fragments still `pending` in the CSV: **6 finished successfully** (orca.engrad + orca.property.txt present) with durations computed from `orca.inp` -> `orca.property.txt` mtimes: `930FragCof` 1787.7 min (29.8 h), `920FragCofMin` 1581.8 min (26.4 h), `613FragCof` 1526.2 min (25.4 h), `187FragCof` 1300.2 min (21.7 h), `215FragCof` 1196.5 min (19.9 h), `498FragCofMin` 1176.5 min (19.6 h). All six now exceed the previous session's max (1278.5 min).
- **1 genuine failure**: `662FragCof` (kept_normal, 208 atoms) - ORCA aborted with `SCF NOT CONVERGED after 2012 cycles`, negative HOMO-LUMO gap, ~6000 negative diagonal Hessian elements. Not an infrastructure/timeout issue - a real electronic-structure pathology. **This is the single worst-scoring fragment (0.14 cosine similarity at r_cut=6A) from the SOAP fidelity analysis two sessions ago** - a monomer of a tightly-stacked framework (c-axis 3.4-4.35 A, shorter than the SOAP cutoffs tested), previously flagged there as poorly representing its parent's true local environment. `187FragCof`, also flagged in that same analysis (0.19 similarity), is among the six that succeeded but took 21.7h - the second-longest runtime of any fragment. The SOAP-fidelity/QM-cost correlation raised as a tentative n=2 observation two sessions ago now has a confirmed failure behind it, not just slow convergence.
- **Reconciled totals**: 644/645 success (99.8%), 1 failed, 0 pending. Mean duration 56.5 min (up from 43.6), median unchanged (15.5 min) - the tail moved, not the bulk. p99 = 1178 min. **3 fragments now exceed 24 hours**; 14 exceed 12 hours; 40 exceed 2 hours. Total serial core-time: **606.8 hours = 25.3 days** (up from 464 h / 19.3 days). `replaced_by_min` remains the costliest category (mean 91.2 min, was 72.4), `kept_normal` mean 50.9 min (was 38.5). Size correlation unchanged in character: weak on raw atoms (r=0.26), strong on log(duration) (r=0.67).
- Artifacts (scratchpad, not committed): `preflight_qm_report_updated.csv` (reconciled: 6 relabelled success with real durations, 1 relabelled failed with the SCF error), `qm_duration_summary_v2.png` (updated histogram + scatter, failure marked with an X). Sent to the user.
- Follow-up: `662FragCof`'s SCF failure needs a decision - retry with different SCF settings (e.g. different initial guess, damping, or a smaller/different active fragment), or accept it as a genuine drop from the QM dataset. Given the SOAP-fidelity link, its dimer counterpart (if the parent's stacking allows one) may be worth trying instead of debugging the monomer's SCF. The report-staleness issue (CSV not live-updated) means any future "how's the run doing" check should go to the per-fragment directories, not just the CSV, until/unless the report generator is fixed to rescan.

## [2026-09-24] Duration-vs-size outlier detection on the QM preflight run
- No source change. At the user's request, tested the hypothesis that unusually slow QM preflight runs are driven by chemically-incorrect fragments (partial rings, near-overlapping atoms), by regressing preflight duration on fragment atom count and pulling out fragments far above the trend.
- **Method**: fit `log(duration) = a + b * n_atoms` on the 644 successful fragments from the reconciled report (previous session), R=0.667. Flagged fragments whose residual exceeds `max(Q3 + 1.5*IQR, mean + 2*sigma)` of the residual distribution - a >6.1x-predicted-duration threshold, chosen over the raw distribution's own shape rather than picked by hand. The 1 SCF failure (`662FragCof`, no duration value) was added unconditionally as the most severe case.
- **Result: 43 outlier fragments** (42 slow + 1 failed) out of 645, spanning 59-250 atoms, up to 44.7x their size-predicted duration (`806FragCofMin`: 1068 min actual vs 23.9 min predicted).
- **Cross-checked against independent quality signals already on file** (valence/overlap QC from the 2026-09-19/21 sessions, SOAP parent-fidelity from 2026-09-21): **23/43 (53%) carry at least one independent red flag** - 4 have an overlapping atom pair (<0.9 A), 7 have an over-valent atom, 18 score below 0.85 SOAP similarity to their parent at r_cut=6A, plus the 1 SCF failure. **17/43 (40%) show no flag at all** in the existing QC (3 more have no QC data on file) - so the "chemically incorrect structure" hypothesis is partially but not fully supported; a meaningful fraction of the slow outliers may instead be genuine electronic-structure difficulty (near-degeneracy, hard SCF landscape) unrelated to a geometry defect.
- **Deliverables, written to the cluster** at `arf:/arf/scratch/otayfuroglu/deepCOF_works/coreCOFs/HCNO/runOpenMLP/subset_qm250/analysis/`:
  - `outlier_fragments.extxyz` - the 43 flagged fragments, pulled from `fragments_subset_qm250.extxyz`.
  - `parent_cifs/` - the 43 corresponding parent CIFs (one per distinct structure).
  - `outlier_report.csv` - one row per outlier: label, atom count, actual/predicted duration, ratio, selection category, and the cross-checked over/under/overlap/ncomp/SOAP-fidelity columns.
  - `outlier_duration_vs_size.png` - duration-vs-size scatter with the fit line, outlier threshold, and outliers/failure highlighted.
- Follow-up: worth having the user (or a visual pass, per the standing "validate COF chemistry by eye" practice) actually look at a handful of the 17 unflagged-but-slow outliers to see what's driving them if not geometry - candidates to start with: `209FragCof` (90 atoms, 40x predicted, no existing flag), `663FragCof` (75 atoms, 30.6x), `1052FragCof` (83 atoms, 20x).

## [2026-09-24] Audit: is anything outside a cut site modified during capping / multiplicity repair?
- No source change. The user asked for assurance that parts of a structure where no bond was cut are left alone by the capping and the multiplicity=1 repair. Answered by reading every routine that can add, remove or move an atom, then verifying empirically.
- **Empirical check** (scratchpad, not committed): fragmented 8 structures locally on current code (1, 481, 526, 671, 1180, 284, 1022, 1091), then matched every fragment atom against all lattice images of its parent CIF. In 22 of 24 fragments **every heavy atom and every parent hydrogen reproduces a parent site to <1e-3 A**; the only atoms that do not match are hydrogens, i.e. the caps.
- **Coordinates are safe by construction.** `enforce_sp2_capped_h_geometry`, `enforce_capped_oh_geometry`, `_planarize_conjugated_caps` and the UFF pass all write only `coords[h]` for `h` in the cap index list. Parent coordinates are never assigned.
- **Removal is safe by construction.** The parity repair's removal branch iterates `capped_h_indices` only, so it can never delete an original atom.
- **Three real exceptions, all documented rather than fixed:**
  1. *Capping is valence-driven, not cut-driven.* `_cap_open_oxygens` scans every C/N/O with no H and spare valence. In a well-formed parent that set equals the cut sites, but an atom under-coordinated **in the parent** (e.g. deposited without its hydrogens) is capped although nothing was severed there.
  2. *Parity H addition is also valence-driven.* On the COF path (`_parity_strict_add_geometry = False`) any C/N/O that `_can_accept_extra_h` approves is a candidate, with no requirement that it sits at a cut. Full-884 run: 268 `[carbon]`, 64 `[secondary-amine]`, 11 `[amide]`, 4+3 amine/alcohol additions; 26 removals.
  3. *Dedupe can delete a genuine parent atom.* `_dedupe_superimposed_atoms` (COF-only, tol 0.85 A) drops the second of any pair closer than that. 1091's parent has a real C-N pair at 0.75 A and 284's has H-H pairs at 0.78 A, so in those broken CIFs a parent atom is removed although no bond was cut. 95 dedupe events across the full run.
- **Stacked-dimer partner layers are not bit-exact.** Where the partner op is a crystal rotation rather than a lattice translation (481), the generated second layer sits 0.006-0.072 A (mean 0.015 A) off the nearest parent site while the first layer is exact. This is the symmetry operation, not the capper.
- **Latent trap noted, not changed:** `refine_h_geometry_with_rdkit` falls back to `movable_h = every H` when `capped_h_indices is None`, which would relax parent hydrogens. Its single caller (line ~3511) always passes explicit indices, so the branch is dead today; a future caller that omits the argument would silently move parent H.
- Follow-up: if the guarantee should be absolute rather than near-absolute, the fix is to carry the set of severed-bond atom indices through to the capper and the parity repair and require membership, and to make the dedupe refuse to drop an atom that came from the parent. Neither is done; both would change output and need a full-884 re-measurement.

## [2026-09-25] COF fragment repair campaign (662 -> 1234) and outlier re-check
- Structure-by-structure repair of the fragments the user flagged by eye, each fix written against the chemical motif rather than the reported structure, per the standing instruction. Commits `e5b6fbd`, `d2eee8f`, `dcb46db`, `5b3790c`, `675bf57`. MOF behaviour was held fixed throughout by the hook pattern (base no-op or class attribute + `COFFragmenter` override), verified after every change with `MOFFragmenter.X is BaseFragmenter.X` identity checks.
- **Linkage chemistries added** (`coffragmentor.py`): the pyrazine and oxazine **ring linkages** (662, 663, 1050-1052), detected by hand-walking bare N/O bridges with exactly two carbon neighbours into a 6-ring - `nx.cycle_basis` was useless here because it returned a 16-membered cycle containing 662's pyrazines. All four ring bonds share one linkage id so the ring is cut as one unit. The **cyanovinylene** linkage (1234): the vinylene rule required a hydrogen on both alkene carbons, so `Ar-CH=C(CN)-Ar` matched nothing and the biaryl fallback cut a triphenylbenzene node down to `C6H3`; `_vinylene_side` now asks for exactly one aryl attachment plus an H *or* a nitrile. The N-N bridge rule was generalised to azine/acylhydrazone/azo.
- **Capping and geometry fixes** (`fragmentation_oop.py`): ring membership for the clipped-heteroatom trim is now decided on the periodic parent, not the supercell (supercell-edge atoms look terminal and ringless); `_terminal_bond_order` takes a per-structure length scale from the median conjugated C-C bond, applied only past 3% deviation and clamped to +/-8%; carbonyl oxygens are sized by bond order rather than neighbour count (85 mis-capped sites); cap groups are accepted or rejected whole in the rebuild; hydroxyl H survives the dedupe; a ketal's second oxygen is no longer read as protonated; `_find_stacked_partner` measures the interlayer approach on the atoms the crystal actually has, excluding capping H, and compares it against the crystal's own closest contact instead of a fixed bar (this alone turned 67 collection-wide monolayers into bilayers, including 131 and 1234).
- **Re-check of all 43 outliers on current code** (`675bf57`): 43/43 in the log, 0 ERROR, 0 timeouts, 0 over-coordinated, 0 odd-electron in the collection, 0 multi-piece, 8 fragments with mis-capped sites, 3 quarantined as genuine radicals (`1050FragCofOnlyNode`, `476FragCofOnlyLinker_1`, `638FragCofOnlyNode`), 2 structures with no recognised linkage (24, 1013 - both amide). Bit-identical to the previous run, so the pipeline is reproducible.
- **New finding: the "90 clashing fragments" warning is mostly not ours.** Measured each parent CIF's own closest non-bonded approach: **25 of the 43 parents already contain a sub-2 A contact**, down to 1.33 A (930), 1.38 A (256), 1.41 A (1004), 1.44 A (525, 1242). These are idealised CoRE-COF models whose aromatic hydrogens sit on top of each other in the crystal itself, and since parent atoms are never moved the fragment carries the contact through unchanged. Only 23 of the 145 fragments have a contact tighter than anything their own parent makes, and only three are severe: **903** (1.21 A, parent floor 2.10), **663** (1.45 A, parent has no hydrogens at all), **1013** (1.62 A, parent floor 1.95).
- **Motif behind those three, not yet fixed**: two capping hydrogens placed independently on neighbouring cut sites end up pointing into the same pocket - 663's two N-H caps on nitrogens 2.65 A apart, 1013's two methyl caps on carbons 2.39 A apart, 903's aromatic C-H cap 1.21 A from an existing ring H across a 2.94 A bay. `_relieve_cap_clashes` (added this session for the bilayer path) already implements the remedy - swing a cap around its cone to the roomiest azimuth - but it is only called *after* layer duplication, not after capping. Calling it pre-duplication in both main paths would cover this, with two prerequisites: `gap()` must exclude the parent's other neighbours (a methyl's geminal H sit at 1.78 A and would otherwise be treated as clashes), and the routine must remap `capped_h_flags` as well as `capped_h_indices` when it drops an H. **Not implemented** - the user's standing instruction is to handle only the structures they point at, and this changes fragment geometry.
- Deliverables: the outlier fragment viewer was rebuilt on the new run and republished (145 fragments + 3 quarantined, bilayer/monolayer state, in-layer and interlayer contacts reported apart, and each short contact marked inherited-from-the-CIF or introduced).
- Follow-up risks: the amide linkage (24, 1013) and the imide family still have no rule; 930 and 476 have broken parents (bonds 1.13-1.81 A); 525's Min has 4 mis-capped sites where coincident interlayer caps were dropped without re-capping; 613's node has 3 mis-capped carbons; the full-884 run is deliberately deferred until the user says the pointed list is done.

## [2026-09-25] Three COF fixes: amide linkage, nitriles read from geometry, crowded caps
- User pointed at 903, 663 and 1013 and supplied the expected node/linker for 1013. Each turned out to be a different defect; all three fixes are written against the motif, not the structure. Commit `de4dae8`.
- **1013 and 24 are polyamide COFs and had no linkage rule.** The amide C-N bond is now cut (`coffragmentor.py`, `_is_amide_bond`). The imine rule cannot reach it - that requires the carbon to have two heavy neighbours and an amide carbon has three - so both structures matched nothing and fell to a fallback that cut through the backbone, leaving two methyl stubs 2.39 A apart. Carry-over then gives the acyl block `Ar-C(=O)NH2` and the amine block `Ar-NH2`, the two halves of the condensation. Result for 1013: linker `C8H8N2O2` (terephthalamide, exactly the user's picture) and node `C26H24N4` (tetrakis(4-aminophenyl)ethylene).
- **The amide predicate's first version broke seven structures that were already correct.** It matched acylhydrazone COFs (`Ar-C(=O)-NH-N=CH-Ar`), whose amide group is real but whose linkage is the N-N bond, and displaced the N-N bridge rule in 333, 334, 498, 635, 638, 795 and 920. Fixed by requiring the amide nitrogen's other heavy neighbour to be a carbon. **Motif count across the 884: 75 structures contain an amide-like carbonyl, but only 13 carry a true secondary aryl amide linkage** - 24, 76, 270, 626, 931, 932, 1011, 1012, 1013, 1014, 1015, 1016, 1017. The other 62 keep the rule that already handled them. Beyond 24/1013 the rule also replaces the crude biaryl fallback in 1012 (node was a bare `C6N6` with no hydrogens, now `C18H18N6O6`) and 1016 (node was a bare hydrocarbon `C32H20`, now two proper amide/amine-terminated blocks).
- **903's nitriles were coming back as anilines.** Its CIF draws C#N at 1.38 A instead of 1.15 A, so bond-order-from-length read a single bond, the carbon looked one valence short and the nitrogen two, and each `C#N` was capped to `CH-NH2` - three spurious hydrogens apiece, one landing 1.21 A from a ring hydrogen, which is the contact that was reported. `_drawn_nitrile_partner` now identifies a nitrile from its collinear geometry (a carbon with two heavy neighbours, one a nitrogen with no other heavy neighbour, within 15 degrees of linear) and scores the bond as a triple in `_local_valence_used` and in `_cap_severed_double_bond_sites`. `check_cof_fragments.py` learns the same rule so report and code cannot disagree. **Motif count: 68 structures have a near-linear terminal-N carbon; 67 draw it short enough to be read correctly already and are unaffected, 903 is the only one that was not.**
- **663's amine caps were pointing at each other.** `_relieve_cap_clashes` (written earlier for the bilayer path) now also runs after ordinary capping in both main paths, since caps on neighbouring cut sites are placed independently. Three changes were needed to make it safe there: it sweeps in several passes so a cap with nowhere to go gets another chance once its neighbour moves; `gap()` ignores the parent's other hydrogens, because a methyl's own hydrogens sit 1.78 A apart by construction and would otherwise read as a clash and invite this pass to distort the group; and a rotation that would crowd those geminal hydrogens is refused - without that guard one hydrogen of an `-NH2` was rotated onto its partner and the pair then deleted as a duplicate, which cost 663's node two hydrogens. It also remaps `capped_h_flags` when it drops an H, or the second layer would be built from a mislabelled cap set.
- **Measured.** 43 outliers: no recognised linkage 2 -> 0, over-coordinated 0, odd-electron 0, multi-piece 0, mis-capped 8 fragments (unchanged - 476, 525, 613, 930, all pre-existing), contacts tighter than the parent's own 23 -> 19 with the worst 1.21 -> 1.83 A, 903 and 663 gone from that list. 24 is now correctly detected as a duplicate of 1013. 60 random COFs (seed 17): no linkage 3 -> 1, mis-capped 0, over-coordinated 0, errors 0.
- **MOF isolation.** Every new method resolves to `COFFragmenter` and to nothing on the MOF side; `MOFFragmenter` still overrides only `METALS`. Running 12 CR MOFs from two pristine `git archive` trees gives bit-identical `fragments_collection.extxyz`. Note for future checks: running MOFs from the repo working directory is **not** reproducible, because `runUniFrag/mof_nodes_lib` and `mof_linkers_lib` cache fragments between runs and change the result; compare from clean trees or clear those directories first.
- Follow-up risks: the vinylene carry-over caps its partner alkene carbon as `-CH3` rather than `=CH2` in structures whose C=C is drawn at aromatic length (903 and 1234 both; 1234's node was reviewed and accepted by the user in this state) - the block should carry a planar vinyl terminus, and fixing it needs the linkage's bond order plumbed from `coffragmentor` through to the capper. 930 and 476 still have broken parents; 613's node has 3 mis-capped carbons; 525's Min has 4; 1238 in the random sample has no recognised linkage; the imide family still has no rule.

## [2026-09-25] Correction: the close-contact warning was measuring the wrong thing
- **User correction, and they were right.** The S4/S6 warning scanned every non-bonded pair in a fragment against a flat 2.0 A bar. Two things wrong with that: it must only judge atoms this code *added*, since the framework's own atoms are never moved or deleted and a contact between two of them is the crystal's geometry; and inside a single molecule the bar for a hydrogen against anything should start at **0.9 A**, not 2.0, because a 1-4 or 1-5 contact in a planar or crowded system is legitimately short. Commit `e0c3a3a`.
- **Provenance is now carried out of the fragmenter.** `FragmentResult` gained a `capped_h` field, the five COF return paths fill it, and `_write_extxyz` emits it in the comment line as `capped_h="..."` - valid extxyz, ignored by any reader that does not ask for it, so the per-atom columns are unchanged. `_parse_extxyz` and `_update_extxyz_collection` preserve it when the collection is rewritten. An empty list is still written, because "no atoms were added" and "nobody recorded what was added" are different facts and the checker has to tell them apart; without provenance it falls back to the 0.9 A floor for all pairs.
- **Effect: 43 outliers 87 flagged contacts -> 3, 60 random COFs 88 -> 3.** Every survivor is *between molecules* (a capping hydrogen of one layer approaching the other, 1.81-1.93 A); **not one fragment has an added atom under 0.9 A inside a molecule.** Nothing else in the report moves - over-coordinated 0, mis-capped 8 fragments, odd-electron 0, multi-piece 0.
- **This invalidates the "inherited vs introduced" framing from earlier today.** Comparing a fragment's worst contact against its parent's own floor was a workaround for measuring the wrong pairs in the first place; the viewer no longer uses it. It also means 663's 1.45 A amine contact was never a reportable defect by this standard, though 903's nitriles-read-as-anilines and 1013's missing amide linkage were real errors regardless of any contact test - which is the argument for not using a distance scan as a proxy for chemical correctness.
- **Checked, not changed:** `_relieve_cap_clashes` still acts at 1.8 A, above the 0.9 A defect bar. It only ever moves added atoms, and it is a geometry improver rather than a defect detector, so the two thresholds are allowed to differ - but it is worth knowing they do. Verified it is not distorting anything: of 214 trigonal sp2 capping hydrogens in the outlier set, 24 sit more than 15 degrees out of plane, and running 1210 on HEAD and on the current code gives the **same 6 out-of-plane caps at the same 61.2 degrees**, so that defect is pre-existing and not caused by the relief. (1210's 6 out-of-plane aromatic C-H caps are an open item in their own right.)
- **Reproducibility caveat, correcting an earlier claim in this log.** I previously called the outlier run "bit-identical, so the pipeline is reproducible". That is not guaranteed: which structure owns a *shared helper block* depends on which parallel process finishes first. Seen twice - two HEAD MOF runs from pristine trees disagreed over whether a shared linker is `ABAYIOFragMofOnlyLinker` or `ABAYOUFragMofOnlyLinker`, and two COF outlier runs differed by one frame (`1049FragCofOnlyNode`) with identical CSV and identical duplicate lists. Fragment *content* is stable; helper-block *attribution* is not. Compare inventories by chemical content, not by file bytes.

## [2026-09-25] 187: a framework with one building block, and two defects it exposed
- User reported 187 coming out with a node and no linker - one building block. It is a **covalent triazine framework whose triazine rings are bridged straight to each other by vinylene**, so `COF.fragment()` correctly returns two identical `C9H3N3` nodes and zero linkers. The node+linker path required both lists to be non-empty, bailed out, and the structure fell to the supercell fallback, where Path D's N-rich core detector (any four nitrogens within six bonds of each other) swallowed the whole sheet and returned **one 466-atom blob with no node, no Min and no linker export**. Commit `9d797f5`.
- **Fix: a node whose partner is another node.** When `linkers` is empty the partner pool becomes `nodes`, with the identity pairing (same SBU, zero image shift) excluded; image scoring, arm coverage and the merge work unchanged on that pool. Four outliers decompose properly for the first time: **187, 635, 638, 1223** - all four previously reached the fallback. 187 goes from `C216H142N108` (466 atoms) to `C60H48N24` (132 atoms, bilayer) plus a Min and `187FragCofOnlyNode` `C18H18N6`, the triazine carrying three vinyl arms. 635's node changes from a `C26H20N8O10` chunk of sheet to the correct `C9H9N3O3` triazine-triol.
- **Defect it exposed 1: vinylene termini were capped as methyl.** The carry-over exists to give a block an intact `Ar-CH=CH2`, but the terminus is sized from bond LENGTH and 187 draws every C-C at 1.40 A, so all three came back `-CH3`. `_severed_alkene_terminus` now decides it from valence: the anchor is a carbon with two heavy neighbours and one hydrogen, so it must double-bond one of them, and the only partner left is the bare carbon hanging off the cut. A bond at 1.47 A or longer is still respected - 402 has a saturated carbon at 1.53 A whose parent is simply missing a hydrogen, and it was swept in until that guard was added. The checker learned the same rule and now prints the deficit implied by the order it actually used, instead of the self-contradictory `H=3 need=3`. This also corrects **131, 903 and 1234** (1234's node `C60H54` -> `C60H48`; it was reviewed and accepted in the methyl state, so it is worth a second look).
- **Defect it exposed 2: a quarantined fragment consumed its duplicate key.** `seen_keys` is claimed in the collection loop, before the odd-electron split, so a later structure with the same heavy skeleton was skipped as a duplicate of something that then went to quarantine and never reached the collection. `1049FragCofOnlyNode` `C9H13N3O6` disappeared exactly that way, suppressed by `1050FragCofOnlyNode`, a genuine radical. Held-back fragments now release the key.
- **Note on the dedup key.** `_chemical_identity_key` is heavy-atom only, by design, so a vinyl-capped and a methyl-capped copy of the same skeleton are duplicates and whichever is written first wins. Before the alkene fix that let 187's methyl-capped node displace 613's correct vinyl-capped one. With the fix both are `C18H18N6` and the dedup is harmless - but the key cannot tell two cappings apart, which is worth remembering the next time a capping rule changes.
- **Measured, same checker throughout.** 43 outliers: mis-capped **15 -> 7** fragments (survivors are 476, 525, 930 - broken parents, all pre-existing), over-coordinated 0, odd-electron 0, multi-piece 0, added-atom contacts 5, all interlayer. 60 random COFs: mis-capped **8 -> 0**, contacts 3 -> 2, no linkage 1, errors 0. MOF `fragments_collection.extxyz` bit-identical to HEAD from pristine trees and carries no `capped_h`.
- Follow-up risks: Path D's porphyrin detector still matches any four nitrogens within six bonds, so it will keep swallowing N-rich sheets for any structure that reaches the fallback for another reason - it wants a real N4-core test (four nitrogens ~2.05 A from a common centre). 1210's six aromatic C-H caps sit 61 degrees out of plane. 930 and 476 remain broken parents; 525's Min is short three amine hydrogens; 1238 in the random sample has no recognised linkage; the imide family still has no rule.

## [2026-09-25] Path D: a real porphyrin test instead of a nitrogen-density test
- Follow-up to the 187 work, at the user's request. The Path D core detector asked for "four or more nitrogens within six **bonds** of each other" - graph proximity, which any nitrogen-rich sheet satisfies. Across the 884 it called **317** structures porphyrinic; on a 175-structure cross-section it fires for **73**. The core it then built was the whole connected nitrogen component, which is why 187's triazine framework came back as one 162-atom core and the structure was emitted as a 466-atom blob. Commit `0706e71`.
- **New test, all four clauses load-bearing.** Four nitrogens ~2.05 A from a common centre (window 1.85-2.35), coplanar (smallest singular value <= 0.12 of the largest), in **one bonded molecule**, and within **8 bonds** of each other through the macrocycle. Each was added because the previous version failed without it:
  - geometry alone gave 187 **432** "cores" - two nitrogens in one layer plus the two stacked 3.5 A above form a flat rectangle whose corners are all 2.1 A from its centre;
  - raw supercell coordinates hid the real cores in **28** and **807**, whose porphyrins straddle a cell edge so their nitrogens sit in different images and measure 30 A apart. Distances are now minimum-image, with the quad unwrapped relative to its first nitrogen for the centroid and planarity test.
- **Validation.** On the same 175 structures: old **73** -> new **21**. Every one of the 21 has a genuine N4 macrocycle, checked independently against a 16-membered C12N4 ring search; the one apparent mismatch, **511**, is a *corrole* - a contracted porphyrin whose inner ring is 15-membered C11N4 with N-N at 2.61-2.76 A - so the detector is right and the ring search was the narrower test. Every structure the old rule over-claimed has no N4 pocket at all. For real porphyrins the core inventory is unchanged: 28, 179, 212, 807 give 14, 10, 28, 14 cores of 24 atoms under both rules. 187, 635 and 1223 go from 9, 10 and 9 bogus cores (84-162 atoms each) to **0**.
- **Regression.** 43 outliers and 60 random COFs are **byte-identical** to the previous run - neither set reaches Path D any more, so the change is inert there. 31 porphyrinic structures: 101 frames both ways, and the only difference (1006 vs 1008 owning a shared block) reproduces on HEAD when run single-process, i.e. it is the known parallel attribution race, not this change. MOF `fragments_collection.extxyz` bit-identical to HEAD from pristine trees.
- **How rare Path D actually is:** zero hits in 200 random structures on HEAD. Before the single-building-block fix, 187 was one of the few that reached it. So this fix is mostly preventative - it stops the fallback from mangling any N-rich structure that lands there for some other reason.
- Follow-up risks unchanged: 1210's six aromatic C-H caps sit 61 degrees out of plane; 930 and 476 have broken parents; 525's Min is short three amine hydrogens; 1238 has no recognised linkage; the imide family still has no rule. Also note the `else` branch Path D falls through to prints "this topology is not implemented yet" and builds a plain 6 A sphere - honest, but a structure that genuinely needs a node/linker split will get a radius-limited fragment there.

## [2026-09-26] Subset rebuilt from the fixed outlier fragments (v2)
- No source change. The user fixed 41 of the 43 duration-outlier structures (chemically-incorrect geometries per their inspection) and re-fragmented them, producing a full Frag/Min/OnlyNode/OnlyLinker set (146 frames) at `arf:.../subset_qm250/cycle_runs_pbed4/outliyer_preflight_qm/outlier_fragments.extxyz` - **not** the file I originally wrote to `analysis/outlier_fragments.extxyz` (that one is untouched; the user placed the fix in a new location for a retry QM cycle). Flagged to the user: **`24` and `1051` are in the original 43-outlier list but absent from the fixed file** - not fixed (or not included) this round; their original (outlier) versions were carried through unchanged.
- Rebuilt the QM-cost subset from scratch: took the full 1726-frame `fragments_collection.extxyz` (884-run), replaced all frames belonging to the 41 fixed stems with the 146 corrected ones, then reapplied the same selection cascade as before - keep Frag if <=255 atoms (250+5, this session's explicit tolerance), else replace with Min if Min<=255, else fall back to all OnlyNode/OnlyLinker frames, else exclude. Any resulting fallback frame still over 255 is dropped (2 cases, unrelated to the fix: `931FragCofOnlyNode` 332 atoms, `1176FragCofOnlyLinker` 293 atoms - same two dropped in the original build).
- **v1 (250 threshold, pre-fix) vs v2 (255 threshold, post-fix), 646 structures**:
  | category | v1 | v2 | delta |
  |---|---|---|---|
  | kept_normal | 459 | 463 | +4 |
  | replaced_by_min | 139 | 138 | -1 |
  | node_linker_fallback | 22 | 20 | -2 |
  | node_linker_from_start | 6 | 6 | 0 |
  | min_only_no_normal | 11 | 11 | 0 |
  | excluded | 9 | 8 | -1 |
  Final subset: **645 frames** (same total as v1, by coincidence). One previously-excluded structure (`1174`) now clears via `replaced_by_min` purely because of the 250->255 threshold change (its Min was 252 atoms, sat between the two thresholds) - not related to the outlier fix.
- Deliverables: `fragments_subset_qm250_v2.extxyz` and `subset_selection_v2.csv` (per-structure category/fragment/atoms, with a "fixed" tag on the 41 corrected stems), both copied to `arf:.../subset_qm250/`. **Not overwriting** the original `fragments_subset_qm250.extxyz` that the completed QM run used - left that alone pending the user's confirmation, since a QM cycle already ran against it.
- Follow-up: confirm whether `24` and `1051` still need fixing before this v2 subset is used for a real QM run, and confirm whether `fragments_subset_qm250_v2.extxyz` should replace the live `fragments_subset_qm250.extxyz` or stay as a separate file.

## [2026-09-26] Reduced-cost subset built for outlier-only preflight QM retry
- No source change. At the user's request, applied the same reduction cascade (keep Frag <=255 atoms, else Min <=255, else OnlyNode+OnlyLinker fallback, else exclude) but scoped to just the 43-structure outlier set (not the full 646), to prep a standalone preflight QM retry.
- Pool: the 41 fixed stems from `outliyer_preflight_qm/outlier_fragments.extxyz` (unchanged since the previous session, md5-verified) + the 2 still-unfixed stems (`24`, `1051`) pulled from the original 1726-frame collection - **flagged again**: these two have not been fixed and carry their original (outlier) geometry into this reduced set.
- Result (43 structures): **30 kept_normal**, **12 replaced_by_min**, **1 node_linker_fallback** (`1016`: Frag 348 -> Min 268, both over 255, fell back to 1 node (78) + 2 linkers (56, 60)), **0 excluded**. Final file: 45 frames (0 dropped for exceeding 255, unlike the full-dataset build).
- QC gate before shipping: checked every H atom in the 45-frame result against the rest of its own fragment for a sub-0.9 A non-bonded contact (H-only scan, not all-pairs, per the project's contact-check convention) - **0 found**.
- Notable size changes on the fixed structures vs their original outlier-flagged geometry: `662` (the SCF failure) 208 -> 108 atoms; several others also shrank materially (e.g. `638` 90 -> 78, `1052` unresolved - now 324/212 Frag/Min pair) - consistent with the user's "half phenyl, too-close atoms" diagnosis being real defects that the fix actually removed, not just repositioned.
- Deliverables, written to `arf:.../subset_qm250/cycle_runs_pbed4/outliyer_preflight_qm/`: `outlier_fragments_reduced.extxyz` (45 frames, ready for a preflight QM submission), `outlier_subset_selection.csv` (per-structure category/fragment/atoms, `24`/`1051` marked NOT FIXED).
- Follow-up: same as before - get `24` and `1051` fixed before trusting their QM result; the SCF-failure history on `662` (now much smaller, 108 atoms) is worth watching on the retry to confirm the fix actually resolves convergence, not just atom count.

## [2026-09-26] Metallo-COF distribution analysis for the CoRE-COF DT1242-v7.0 database
- New work, no fragmentation-code change. User wants to start working on metallo COFs with UniFrag; first step was a census of which metal and how many structures across `runUniFrag/CoRECOF`.
- Added `runUniFrag/analyze_cof_metals.py`: parses every CIF with pymatgen, classifies a structure as metallo if it contains any element from a fixed metal set (alkali/alkaline-earth, transition metals, post-transition metals Al/Ga/In/Sn/Tl/Pb/Bi), and reports per-metal structure counts plus multi-metal structures. Deliberately excludes B/Si/P/Sb (common COF linkage/backbone elements, not framework metal centers) from the metal set.
- Ran against `CoRECOF/CoRE-COFs_DT1242-v7.0` (1242 CIFs, 0 parse failures) - the only one of the three CoRECOF folders that carries metals; `CoRE-COF_HCNO884` and its `_clean` copy are pre-filtered to H/C/N/O-only and contain zero metal structures by construction. Output: `runUniFrag/CoRECOF/cof_metal_analysis.md`.
- Result: **99/1242 (7.97%) are metallo-COFs**, 1143 organic-only. Per-metal: Cu 28, Co 22, Zn 16, Ni 13, Ru 5, Mo 4, Au 4, Li 4, Mn 2, Ti 2, V 1, Rh 1, Cd 1, Hg 1, Na 1. 93 structures carry exactly one metal, 6 are bimetallic (`1038` Co/Ni, `1039` Cu/Ni, `1140`/`1142` Co/Mo, `1141`/`1143` Mn/Mo) - no structure has 3+ distinct metals.
- Cross-checked against the pre-existing `runUniFrag/CoRECOF/cof_element_analysis.md` (same source, full elemental prevalence, generated 2026-09-16 by `analyze_cof_elements.py`): the metal counts match exactly, confirming the dataset on disk hasn't changed since that report.
- Caveat (not yet verified): this is presence-only, from parsed composition - it does not confirm the metal is actually coordinated inside a framework node (vs. a residual counter-ion/guest in the CIF, e.g. the single Li/Na hits are worth a by-eye check before treating them as coordinated metallo-COF nodes). No guest-removal step is known to have been run on this COF dataset (unlike the MOF `remove_guests.py` pipeline).
- Follow-up risks: before running UniFrag fragmentation on the metallo subset, spot-check a few Cu/Co/Zn/Ni structures (the four with enough count to matter) to confirm the metal sits in a real coordination node and not as an isolated ion; Path D (metallo-porphyrin/phthalocyanine detector) and the generic MOF-style node/linker paths may both need exercising since COF metal environments vary in coordination motif.

## [2026-09-26] Fixed silent Zn loss in COF Path D/fallback fragments (code change, see project-decisions.md)
- Follow-up to the Zn census above. Copied all 16 Zn CoRE-COF CIFs to `runUniFrag/CoRECOF/zn_cofs/` and ran `./run_cof_family.sh runUniFrag/CoRECOF/zn_cofs 4.0 both`. 3 of 16 (`822`, `823`, `824` - pure ZnN4 phthalocyanine/porphyrin sheets with no coffragmentor-recognised linkage) came back with **zero Zn atoms** in their fragments, despite a clean, uniform ZnN4 coordination in the parent CIF (~1.98-2.01 A Zn-N, verified via `pymatgen` neighbor search) identical to two other structures (`454`, `616`) that fragmented correctly.
- Root cause and fix: `COFFragmenter.COV_RAD`/`is_valid_bond` (`fragmentation_oop.py`) has no metal radii, so any metal silently gets a carbon-like 0.77 A default and a real ~1.98-2.01 A M-N coordination bond is computed as out of range - the metal is then a graph-disconnected atom, invisible to every BFS the Path A-D fallback uses, and gets dropped with no warning (Path J/`coffragmentor.py`, used by `454`/`616`, resolves its own bonding and was unaffected). Added `COFFragmenter.METALS` and a metal-aware bond rule mirroring `MOFFragmenter`'s convention (no metal-C/metal-H bonds, <2.6 A for other metal-nonmetal contacts), applied to both `is_valid_bond` and the `_chemical_identity_key` dedup topology hash. Full decision, root-cause detail, and validation are in `project-decisions.md` (2026-09-26 entry) - not duplicating here.
- Files touched: `fragmentation_oop.py` (COFFragmenter class only). `runUniFrag/analyze_cof_metals.py` (new, from the census work). `runUniFrag/CoRECOF/zn_cofs/` (new working folder, 16 CIFs + regenerated `fragmentation_summary.csv`/`fragments_collection.extxyz`).
- Validation: `./run_fast_test.sh --kind mof` 8/8, `./run_fast_test.sh --kind cof` 8/8. `test_on_cof_zn_pc_series` (5 structures, Path J) verified byte-identical before/after by stashing the change and re-running both ways. Re-ran all 16 Zn CoRE-COFs after the fix: every `FragCof`/`FragCofMin` frame now contains the expected Zn count, with correct coordination geometry (spot-checked `822`/`823`/`824`, `1209`, `73` - Zn-N/O distances match the parent crystal).
- Beyond the 3 originally-broken structures, the fix also silently repaired two cases of *partial* metal loss in minimized fragments from the same batch: `1209FragCofMin` was `Zn1` (153 atoms, 3 of 4 parent Zn dropped during minimize-trimming) and is now `Zn4` (251 atoms, matching normal); `74FragCofMin` was `Zn1` (103 atoms) and is now `Zn3` (165 atoms). `73`/`74FragCof` also shed a few capping H atoms (same Zn count) from the corrected metal-C/H exclusion changing which cut sites get capped.
- Follow-up risks (carried into project-decisions.md): only validated against Zn so far. Before trusting the fix across the full 99-structure metallo-COF set, spot-check at least one Cu, Co, and Ni CoRE-COF structure - those metals can have different coordination numbers/geometries (e.g. Cu(II) Jahn-Teller elongation) not yet checked against the flat 2.6 A cutoff. Also worth widening the regression check beyond the ZnPc series and the 8/8 fast tests to a larger random COF sample, per the project's usual practice for fragmentation-rule changes.

## [2026-09-26] 823: no edge methyls, planar sp2 caps (code change, see project-decisions.md)
- User review of 823 (fused ZnPc sheet, Path D fallback): twelve -CH3 at the fragment edge, and one phenyl cap H out of plane in the Min. Parent is fully sp2, so the methyls were ring carbons cut from two of three ring neighbours.
- Changes in `fragmentation_oop.py` (COFFragmenter only): `_drop_clipped_ring_stubs` (renamed from `_drop_clipped_ring_heteroatoms`) now also drops parent-ring carbons left with one heavy neighbour, iterates, and returns kept indices; the pre-minimize call in the Path A-D branch uses them to remap `local` (a stale `local` after any drop made the minimize trim cut through 823's macrocycle - latent before today, the heteroatom rule could trigger it too); `_relieve_cap_clashes` splays sp2 caps in-plane (+/-20 deg) before ever using the out-of-plane cone swing, which is now a < 1.2 A collision fallback only.
- Results: 823 Frag 205->169 (CH3 12->0), 823 Min 107->98 (CH3 3->0, out-of-plane sp2 H 2->0), 824 Min 96->78 (CH3 2->0, out-of-plane 3->0). Macrocycles intact, Zn 4-coordinate, no under-coordinated atoms. `runUniFrag/CoRECOF/zn_cofs/` collection/CSV updated in place for 823/824.
- Validation: 60 random HCNO COFs 189/189 frames identical to HEAD; ZnPc series 17/17 identical to HEAD; 16 Zn COFs identical except the 3 target frames; fast tests MOF 8/8, COF 8/8. Baseline CH3/out-of-plane counts in the random sample (206/65) were checked against the parents - they are inherited methyl substituents and parent non-planar X-H, not this artifact.
- Open: (1) 823/824 Min still carry one sp3 CH2 from the shared odd-electron repair - its added H is out of plane by construction; asked the user how to handle parity. (2) Found in passing: no `OnlyNode`/`OnlyLinker` export in the Zn batch contains Zn - metallo-macrocycle nodes (ZnPc/porphyrin, e.g. 607 `B8C64N16O16`, 454/1221 `C88N16`) are exported without their metal. Not fixed; flagged as a separate task.

## [2026-09-26] Odd-electron repair on functional-group N/O instead of ring carbon (code change, see project-decisions.md)
- User asked to fix parity by capping an N with H or removing an H, rather than the CH2 the repair put on 823's phenyl. Recorded all 19 odd-electron fragments in 60 random HCNO + 16 Zn COFs: 15 were repaired on carbon; all 19 had an N/O route that `_can_accept_extra_h` refused (bond-order-from-length misreads delocalised C-N/C=O).
- Added `BaseFragmenter._heteroatom_parity_repair` (no-op) + COF implementation: tier A adds an in-plane H to a pyridinic ring N (ring <= 16), a terminal =N-H cap (-> -NH2) or a carbonyl O (-> enol); tier B removes a cap H from a conjugated secondary N or an O on an sp2 C/N (nitro); else the old carbon path.
- Repair kinds, 76 structures: before 15 carbon / 2 amide / 1 secondary-amine / 1 unrepairable; after 1 carbon / 6 ring-N / 6 imine-NH / 4 carbonyl-O / 2 cap-OH / 0 unrepairable. Changed frames: 823Min, 824Min, 1057Frag/Min, 145Min, 492Frag/Min, 555OnlyLinker (nitro restored), 794Frag/Min, 929Min, 987Frag, +987OnlyNode_1 (previously quarantined). All new parity H in plane (<= 0.7 deg), clearance 1.78-2.06 A.
- Validation: ZnPc series 17/17 identical to HEAD; fast tests MOF 8/8, COF 8/8. `runUniFrag/CoRECOF/zn_cofs/` updated for 823/824.
- Follow-up: widen to the full 884 HCNO set before a QM run - the 60-structure sample exercised 18 repairs, but not every chemistry (e.g. no boronate O-H caps were picked).

## [2026-09-27] Cu / Co / Ni CoRE-COF batches: run + quality check (no source change)
- Built `runUniFrag/CoRECOF/{cu,co,ni}_cofs/` from `CoRE-COFs_DT1242-v7.0` (Cu 28, Co 22, Ni 13; bimetallic 1038 Co/Ni and 1039 Cu/Ni sit in both of their metals' folders) and ran `fragmentation_oop.py <folder> --kind cof --radius 4.0 --nproc 3` (folder mode, not run_cof_family.sh). Parent metal sites are square-planar MN4 at 1.83-2.14 A except Co 93 (CN6), Co 1140/1142 (CoO6 in a Mo6 polyoxometalate); all within the 2.6 A metal-bond cutoff.
- New checker `runUniFrag/check_metal_cof_fragments.py <folder> <Metal>`: metal kept, metal CN vs parent, odd electrons, pieces, CH3/CH2 absent from parent, bent 2-coordinate C, out-of-plane caps, added-H contacts (0.9 A intra / 2.0 A interlayer, capped_h only).
- Result: 103 Frag/Min frames checked; metal kept with the parent's CN in all but the cases below. Timeouts (300 s): Cu 917, Cu 1208 (48 Cu per cell). Dedups are legitimate (656=655, 701=700, 124=125 minus BF4-, 915/916=914, 77).
- Defects found, all pre-existing (78/272/843 reproduced byte-for-byte on HEAD):
  1. Bimetallic 1038/1039 lose the second metal (Co / Cu) in the main FragCof: Path J takes one MPc as node and the other as a "linker", and linker handling strips metals. Same root cause as the spawned OnlyNode/OnlyLinker metal-loss task. Also why 1039 is flagged a "duplicate" of 1038 (both reduce to the same Ni2 skeleton).
  2. Ni 893 loses Ni: imine linkages recognised but Path J fails; the fallback takes a single O atom as its node (the degenerate `{B,O}` node noted 2026-09-16) and never reaches the metal. No Min.
  3. Co 93 yields only BF4- (5 atoms, then quarantined): cationic framework with BF4- counterions and no recognised linkage; the fallback picks the lone B as node. Counter-ions/guests are not stripped for COFs.
  4. Alkyne cut sites mis-capped: Zn 610 exposed diyne carbons capped =CH2 (1.26 A bond read as double); Zn 1183/1209 bent/under-capped sp carbons; Ni 271 / Cu 272 terminal C-H at 120 deg (these parents also carry 4 bent H-less "alkyne" carbons).
  5. Bare aromatic C (no cap) at linkage cut sites: Cu 78 (biaryl cut, 2 in Frag / 1 in Min), Cu 843 (4 / 1). Nothing caps a ring carbon left with two heavy neighbours after a linkage cut; the parity repair happened to cap one.
  6. Minor: Co 731 Min one ring CH2 (extra cap on an aromatic CH during minimize trim); Zn 74 keeps Zn CN4 of the parent's CN6 (two O at 2.52 A dropped); Cu 189 one interlayer cap H...H 1.97 A.
- No OnlyNode/OnlyLinker export contains its metal in any folder (known, spawned task).

## [2026-09-27] Metallo-COF issues 1-3 fixed: metals kept, 893/1038/1039/93 recovered, ligand parity (code change, see project-decisions.md)
- Files: `coffragmentor.py` (skip pyrazine/oxazine bridge-only remnants; carry a bridge's H with its N), `fragmentation_oop.py` (`_KEEP_METALS_IN_BLOCKS`, Path J `partners = linkers + nodes`, `COFFragmenter._strip_discrete_guests` before the supercell fallback, metal-centred fallback node, `_ligand_electron_count` + `_PARITY_IGNORES_METALS`), `runUniFrag/check_metal_cof_fragments.py` (odd ligand parity = defect, open-shell metal = info).
- Root causes: 1038/1039 - dihydropyrazine N-H remnants became fake 2-atom linkers, and Path J only paired the node with linkers, so the second MPc node was never used. 893 - Ni-porphyrin node bonded only to other nodes -> Path J gave up -> single-O fallback. 93 - BF4- counter-ions + no linkage -> fallback took B as node; then the metal-free JLU2 ring motif dropped Co. All OnlyNode/OnlyLinker - `_clean_linker_molecule` strips METALS. Cu/Co - parity repair edited the ligand to pair the metal's own unpaired electron.
- Checker tallies before -> after: Cu metal_lost 1->0, odd_ligand 22->0 (22 open-shell Cu(II) now reported as info), under 6->5, ch2 0->2 (1207, carbon parity fallback); Co metal_lost 1->0, odd_ligand 16->0, ch2 1->0, frames 39->41; Ni metal_lost 1->0, frames 20->24; Zn unchanged.
- Regression: random 60 HCNO 53/60 structures bit-identical; 7 grow because node-bonded arms now get their partner block (259, 310, 768, 809, 928, 980, 987; details in the decision); ZnPc series identical except the OnlyNode metal; fast tests MOF 8/8, COF 8/8.
- Folders refreshed in place from folder-mode runs: `runUniFrag/CoRECOF/{zn,cu,co,ni}_cofs` (CSV, collection, run.log). Cu 917/1208 stay in `cu_cofs/timed_out_structures/` (300 s timeout).
- Remaining from the 2026-09-27 list: #4 alkyne cut-site caps (Zn 610 =CH2, Zn 1183/1209, Ni 271/Cu 272, now also HCNO 259), #5 bare aromatic C at linkage cuts (Cu 78, 843), #6 minor (Co 731 fixed incidentally; Zn 74 CN4 of CN6; Cu 189 interlayer 1.97 A). QM note: Cu(II)/Co(II) fragments need doublet (or higher) multiplicity.

## [2026-09-27] Issues 4-5 fixed: alkyne cut-site caps and lost parent H (code change, see project-decisions.md)
- Mapped every flagged carbon to its parent site first. Alkyne ends are linear sp in the parent: 610 draws Ar-C#C-Ar with 1.26 / 1.53 A bonds, and 271/272/259 have a 1.21 A triple plus a 1.43 A aryl bond. The bare aromatic C in 78/843 each had their own H in the parent, which was deleted because the parent draws H of neighbouring rings 0.5-0.7 A apart.
- `fragmentation_oop.py` (COF-only paths): new `_severed_alkyne_terminus`; `_cap_severed_double_bond_sites` gives cut alkyne ends one on-axis H (and removes surplus caps); `_fix_terminal_cap_geometry` skips them; `_skip_saturated_h_cap` treats linear H-free two-neighbour C as complete; `_relieve_cap_clashes` never drops a parent H; `_dedupe_superimposed_atoms` only merges H that share a heavy anchor. The checker now accepts non-numeric CIF names.
- Traced causes: the 120 deg bend came from `_fix_terminal_cap_geometry` (sp2 slotting); the extra H on alkyne carbons from `_cap_open_oxygens` (a 1.26 A triple read as a double); the lost ring H from the final H-H drop in `_relieve_cap_clashes`.
- Checker after: Zn {frames 31, metal_lowCN 2 (74)}; Cu {45, open-shell 22, contacts 1 (189 interlayer 1.97 A), ch2 2 (1207 carbon parity fallback)}; Co {41, open-shell 17}; Ni {24} - all alkyne / bare-C flags gone. Folders refreshed in place.
- Regression: random 60 HCNO 181/194 frames identical; changed = alkyne caps on axis (259, 412, 928), ring H restored (898 +16, 1057 +2; both parents carry 0.88 / 0.47 A H-H), ZnPc-PPE 242/218 -> 226/154 (its alkynes were CH=CH2 before). Fast tests MOF 8/8, COF 8/8.
- Still open: Zn 74 keeps CN4 of the parent's CN6 (two O at 2.52 A dropped); Cu 189 interlayer cap 1.97 A; 1207 carbon parity fallback; Cu 917/1208 time out at 300 s; the 2026-09-26 QM subset predates today's changes (PPE, node-partner, parity and alkyne changes would alter some frames if regenerated).

## [2026-09-29] Review of the parent-CIF-aware capping work, and four corrections
- Reviewed the uncommitted worktree changes at the user's request, against their own base `675bf57`. The diagnosis behind them is sound - fragment-local connectivity cannot tell a naturally low-coordinate atom from a cut site - but the implementation regressed fragment chemistry, so four defects were fixed in place.
- **Cap count came from heavy-atom DEGREE, which cannot see valence.** `_parent_cut_deficit` returned `real_parent_degree - fragment_degree` and skipped the site when that was zero. That cannot distinguish `1050`'s carbonyl oxygen (C=O 1.20 A, complete) from `1049`'s hydroxyl oxygen whose hydrogen the deposition omitted (C-O 1.34 A, one short) - both have parent degree 1 and fragment degree 1, so both were skipped and `1049FragCof` went `C57H42N6O21 -> C57H30N6O21` with six bare single-bonded oxygens. Degree also cannot say how MANY hydrogens a cut site wants. Now computed as `target_valence - _local_valence_used(idx)`, which needs no parent lookup at all: a heavy atom's coordinates here ARE its parent coordinates, so the retained bond's length and order are the parent's own.
- **A clamp charged the retained bond a single valence.** `_cap_severed_double_bond_sites` had been rewritten to use the SEVERED bonds' order plus `min(..., target_valence - 1 - h_count)`, on the theory that the retained bond never needs classifying. It does - the clamp assumes it is worth 1, so a site retaining a *double* bond took one hydrogen too many: `1050`'s ring nitrogens came back as N(=C)(H)(H), four bonds on nitrogen. Reverted to the original bond-order arithmetic; the parent offers no escape here, because it is the same bond at the same length.
- **`_can_remove_cap_h` keyed on a degree special case.** `real_degree - local_degree == 1 and h_count == 1` protects a single-cut site but says nothing about one that lost two bonds and needs both hydrogens. Now `h_count > target_valence - heavy_valence_used`.
- **A fixed FRACTIONAL match tolerance.** `_parent_index_for_position` compared Chebyshev fractional distance against 0.02 - 0.065 A along `1050`'s 3.29 A stacking axis but 0.55 A along its 27.29 A in-plane axis, and up to 0.97 A across the outlier set. Now a 0.25 A metric distance. No mis-mapping had actually occurred (0 element mismatches across 5750 heavy atoms), so this was latent.
- **One reported finding was wrong and is withdrawn.** The extra quarantining (3 -> 8 fragments) is correct behaviour, not a defect: `795` and `930` are odd-electron with no legitimate repair site, and base "repaired" them by putting a hydrogen on a nitrogen with two heavy neighbours. Every nitrogen in both parent CIFs has two or three heavy neighbours and **zero** native hydrogens, so that hydrogen is chemistry the parent does not have. Refusing it and quarantining is right under the standing rule that odd-electron fragments never enter the main collection.
- **Measured** on the 43 outlier CIFs, on quantities independent of any bond-order call the code makes: mis-capped sites **22 (base) -> 414 (as sent) -> 16 (fixed)**; bare O with a C-O bond over 1.28 A **12 -> 126 -> 12**; single-bonded N with no H **0 -> 30 -> 0**; N with a C=N and two H **17 -> 95 -> 17**. The fixed version beats its own base. MOF `fragments_collection.extxyz` byte-identical across 12 CR MOFs, so the base-class hooks remain true no-ops.
- Dead code removed: `_parent_severed_neighbor_valence` and `_real_parent_native_h_count` served only the "never classify the retained bond" idea and have no callers left.
- **Not addressed, and blocking a merge:** this worktree is based on `675bf57` while `main` is at `391f66f`, ten commits ahead, including `57456f7` which rewrote the same `_cap_open_oxygens` region for metallo-COFs. `git apply` of this diff onto main fails at `fragmentation_oop.py:834`. It needs a real rebase, and the fixed version should be re-measured after it - main's capping rules (geometric nitriles, sp2 vinylene termini, the added-atom contact rule) interact with these.

## [2026-09-29] Rebase of the parent-CIF-aware capping work onto main, and re-measurement
- The branch carried no commits of its own, so the work was committed and rebased onto `391f66f` (11 commits, including `57456f7` "keep metallo-COF metals and correct cut-site capping"). One code conflict, in the `BaseFragmenter` hook block: main had added `_heteroatom_parity_repair` and this work had added `_can_remove_cap_h` / `_parent_cut_deficit` at the same place. The additions are disjoint (confirmed by definition census across base/main/worktree), so both were kept. The two doc conflicts were append-ordering only and were resolved by each file's own convention - `project-agent-log.md` oldest-first, `project-decisions.md` newest-first.
- **The merge composes correctly.** In `_cap_open_oxygens` main's `_skip_saturated_h_cap` (the linear-sp alkyne guard) now runs before this work's `_parent_cut_deficit` (the O/N count): the guard removes alkyne carbons from consideration, the deficit decides how many H an O/N site takes. Two stale comments that still described the superseded degree arithmetic were rewritten to the valence formulation actually in the code.
- **MOF unaffected.** All six shared hooks satisfy `MOFFragmenter.X is BaseFragmenter.X` with a `COFFragmenter` override. `fragments_collection.extxyz` over 12 CR MOFs is byte-identical between `391f66f` and the rebased tree; the only difference in `fragmentation_summary.csv` is the row order of one entry, which is the known parallel-worker nondeterminism.
- **Re-measured** on the 43 outlier CIFs from pristine `git archive` trees, all four collections scored with one instrument (the rebased tree's `check_cof_fragments.py`, so the fragmenter's own bond perception is held constant):

  | tree | frames | quarantined | mis-capped fragments | over-coordinated | contact warnings |
  |---|---|---|---|---|---|
  | base `675bf57` | 145 | 3 | 14 | 0 | 0 |
  | base + this work (pre-rebase) | 140 | 8 | 11 | 0 | 0 |
  | main `391f66f` | 148 | 0 | 11 | 2 | 3 |
  | main + this work (rebased) | 145 | 3 | **9** | 2 | 3 |

  The capping gain survives the rebase: 11 -> 9 mis-capped fragments. The two that leave (`930FragCof`, `930FragCofOnlyLinker_1`) move to quarantine rather than being emitted, which is the correct refusal - main emitted them with mis-capped carbons (`H=3 need=2`).
- **The three rebased quarantines are all correct refusals.** `903FragCofOnlyNode`: main's parity repair puts a second hydrogen on C19, an aromatic ring carbon that already has one H and two ring carbons at 1.37 A, making a ring CH2 - exactly what the standing rule forbids. `930FragCof` / `930FragCofOnlyLinker_1`: every N in 930's parent has 2-3 heavy neighbours and zero native H, so there is no legitimate site. Main's `_heteroatom_parity_repair` does legitimately resolve `795` and `930FragCofMin`, which this work had quarantined against base.
- **Pre-existing on main, not introduced here** (identical with and without this work, so reported rather than fixed): (1) the `over_coordinated` MUST failure on `525FragCof`/`525FragCofMin` is an artifact - `is_valid_bond("H","H",0.68)` is True, and 525's *parent* CIF places two hydrogens 0.684 A apart, so bond perception reads a parent-parent H...H contact as a bond. All 24 flagged sites are H-H pairs. Main deliberately stopped deleting parent H for proximity (`57456f7`), which is right under the rule that parent atoms are never moved or deleted; it is the checker's valence test that is measuring parent geometry. (2) Three contact warnings at 1.89-1.91 A between molecules (`635`, `662`).
- Follow-up risk: a bond-perception change for H-H would touch `is_valid_bond`, which MOF fragmentation also uses, so it must not be made without an MOF byte-identity check.

## [2026-09-29] H-H bond perception: two hydrogens are never bonded (code change, see project-decisions.md)
- Fixed the `over_coordinated` MUST failure reported in the rebase measurement. `COFFragmenter.is_valid_bond` treated an H...H below 0.9 A as a covalent bond, so 525's parent-CIF geometry (two hydrogens 0.684 A apart, drawn that way in the deposited model) made three hydrogens per fragment look two-coordinate. All 24 flagged sites were H-H pairs.
- Counted the motif over all 884 HCNO structures before changing anything: 18 structures, 130 pairs below 0.9 A, closest 0.051 A; in every pair both hydrogens already have a heavy-atom neighbour within 1.3 A, so none is an H2. The rule is now "H-H is never a bond", not a tighter distance.
- `is_valid_bond` turned out **not** to be a shared hook - `BaseFragmenter` defines none and each of the three subclasses has its own - so the change is COF-only by construction, and the earlier caution about it being shared with MOF was wrong. `MOFFragmenter`'s own `"H" in (s1,s2) -> dist < 1.2` is untouched.
- Validation: MOF `fragments_collection.extxyz` byte-identical to `391f66f` over 12 CR MOFs. On the 43 outliers `over_coordinated` 2 -> 0 (MUST: all passed), `bad_terminal_caps` 9 and `clashing` 3 unchanged, quarantine unchanged. 142 of 145 fragments byte-identical; only 525's three differ, by capping hydrogens shifting 0.05-0.9 A (cap relief reads neighbours through `is_valid_bond`), with formula, cap count, parity and 0.96-1.07 A cap bond lengths all unchanged.
- Follow-up: `MacromolFragmenter.is_valid_bond` has the same `dist < 0.9` guard and the same argument applies, but there is no PDB input in the repo to validate a bio run, so it was deliberately left alone.

## [2026-09-29] Dioxin linkage requires fusion (code change, see project-decisions.md)
- `1106` returned node and linker as two copies of the same TAPT-like unit, with the dioxine ring in neither. Instrumenting `COF.fragment()` showed 30 cuts after the linkage rules collapsing to 6 **untagged** C-O cuts after the orphan guard: the dioxin rule had opened the linker's own core ring, the framework shattered into sub-minimum pieces, and the guard unwound every tagged linkage trying to repair it.
- The rule's "two different aryl systems" test deletes every oxygen and takes connected components, which in an imine-linked framework with a dioxine-cored linker cuts each linker in half - so the test passed for the wrong reason. Added `_bridge_is_fused`: each C-C bridge of the ring must also be an edge of some OTHER ring, searched with the ring's own two oxygens removed.
- **A first version required an all-carbon ring and was wrong.** `110` and `1239` are genuine dioxin COFs with a pyrazine on the far side; that version took them from 48/96 cuts to zero and emitted no blocks. Caught by running the before/after comparison over every candidate rather than trusting the standalone scan that had classified them. The shipped test is element-agnostic.
- Verified on all 11 structures whose geometry can reach the rule (found by replicating its own predicate over the 884 collection): 6 genuine dioxin COFs byte-identical in cuts and block formulas; `1179` 12->8 cuts and `816` 9->6 cuts with identical blocks; `819` from no blocks at all to a biaryl decomposition; `1106`/`1107` changed but both rejected by the pre-scan. An exhaustive C4O2-ring scan of the 884 was still running at commit time and is a follow-up confirmation, not a precondition - the candidate set was derived from the rule's own predicate.
- Correction: an earlier entry claimed 7 misfiring structures from the all-carbon scan. Two were false positives of that test. Three are actually fixed: 1179, 816, 819.

## [2026-09-29] Parent pre-scan and quarantine (code change, see project-decisions.md)
- Added `COFFragmenter.prescan_parent()` plus wiring in `_process_cof_file`: a parent whose own geometry cannot yield a valid fragment is moved to `prescan_quarantine/` with a `prescan_report.csv` naming the rules and offending atoms, before any fragmentation runs.
- Six rules, thresholds measured over the 884 collection rather than assumed: close-contact <0.90 A (18), over-coordinated with quaternary-N+ exemption (11), isolated-atom (7), stretched C-C >1.80 A (6), peroxide O-O (2), bare carbon (1). Union **37 of 884, 4.2%**. Odd electron count was measured (25 structures, only 6 overlapping) and deliberately left out as a chemistry question rather than a modelling error.
- Validation: the shipped code reproduces the census exactly - 37 structures, every per-rule count matching. Smoke test on six CIFs quarantined 1106/930/284 with correct per-rule reasons and fragmented 662/704/1234 normally.
- **MOF unaffected, and proven rather than assumed.** `prescan_parent` is on `COFFragmenter` only. A MOF byte-identity check first reported DIFFERS; running `391f66f` against itself three times showed run1 != run2 and run2 == run3, flipping the `ABAYIO`/`ABAYOU` `OnlyLinker` attribution - the pre-existing shared-block race. Only those same two frames differed between `391f66f` and this build, so the difference is the race, not this change. Note for future work: `fragments_collection.extxyz` is NOT reproducible run-to-run for MOF sets containing metal analogues that share a linker, so byte-identity is the wrong acceptance test there; compare frames by label and content and treat an ownership swap as a pass.
- Also measured while designing this (no code change): guest detection needs no CSD licence. Over all 884, the existing periodicity test in `remove_guest_molecules.py` and pymatgen's `get_structure_components` (Larsen) each agree with CSD `is_polymeric` on **883/884** and with each other on 884/884. The single disagreement is `70.cif`, where JmolNN bonds two hydrogens 0.855 A apart; routing the same test through `COFFragmenter.is_valid_bond` (which since this session never bonds H to H) gives **884/884**. `ccdc` is used nowhere in the analysis path - only in two optional CSD download scripts.

## [2026-09-29] Parent pre-scan extended to MOFs (code change, see project-decisions.md)
- `prescan_parent` is now a hook: no-op on `BaseFragmenter`, overridden by `COFFragmenter` and `MOFFragmenter`, inherited as a no-op by `MacromolFragmenter` (no PDB collection to measure against, so bio is untouched). Wired into `_process_mof_file` alongside the existing COF wiring.
- Thresholds measured over 5226 structures of `cr_cifs_noduplicated`, NOT carried over from COF. Three COF rules were actively wrong here: the neutral-organic valence table flagged 11.9% of structures for ordinary sulfonate/phosphonate/perchlorate/quaternary-N chemistry; "peroxide O-O" was catching the generic `dist < 1.8` rule rather than peroxides; and treating an isolated ATOM as a defect would have rejected 238 structures whose only feature is a chloride counter-ion or an unmodelled water oxygen.
- **Inspecting the rejects before committing paid off.** A first pass rejected 95. Reading three of them found two false-positive classes: AFEJOK's phosphorus with six bonds is PF6-, and ATAYEB's "F with two bonds" is a C-F plus a phantom F-F at 1.657 A. Raising P to 6 and excluding F-F removed **19 false positives and added none**; the final figure is 76 (1.45%).
- Final: over-coordinated 29, overlapping-oxygens 34, isolated-hydrogen 13, bare-carbon 8, close-contact 3, stretched-CC 0. The last two do almost nothing on this collection - these are refined CSD structures, not the idealised Materials Studio models that give CoRE-COF its superimposed atoms - but they cost nothing and guard a real defect class.
- Validation: the shipped code reproduces the census exactly (76, every per-rule count matching).
- Follow-up worth its own decision: 2065 of 5226 structures (39.7%) contain a metal outside `MOFFragmenter.METALS` (Cd 599, Eu 170, Ag 142, Tb 138, Gd 113, In 99, La 87, Nd 86, U 74, Ce 69) and fail with "No metal found in the input structure". This is a capability gap, not a parent defect, so the pre-scan deliberately does not reject them - but it likely affects MOF coverage more than anything the pre-scan catches.

## [2026-09-29] Parent pre-scan extended to macromolecules (code change, see project-decisions.md)
- Third and final `prescan_parent` override, completing the hook: no-op on `BaseFragmenter`, overridden by `COFFragmenter`, `MOFFragmenter` and now `MacromolFragmenter`. Wired into `_process_bio_file`, quarantining to `prescan_quarantine/` as the other two paths do.
- Runs on the PDB **before PDBFixer**, so heavy atoms only. `bare-carbon` and `isolated-atom` are dropped - an X-ray PDB has no hydrogens, so the first would fire on every carbon and the second on every water.
- **Thresholds are carried over, not measured.** `test_on_bio_mol` holds two PDBs and both are clean on every rule. That is a negative control, not a calibration, and it is recorded as such: two structures cannot set a threshold.
- Because there were no positive examples, the rules were tested against **synthetic defects injected into a real file**. That immediately found a hole: a duplicated carbon at 0.45 A and a nitrogen at 0.50 A were missed by every rule, because `is_valid_bond` calls a 0.45 A C-C a bond, so it never registers as a close contact and the degree stays within valence. This is exactly the altLoc pattern, and `_parse_pdb` ignores the altLoc column. Added a `duplicate-atom` rule (same element, heavy, <0.85 A - the figure the COF dedupe already uses); all six synthetic positives are now caught and both real PDBs still pass.
- Checked whether the same hole affects the committed COF and MOF pre-scans: COF **0** structures, MOF **3** (WILDUR Cl-Cl 0.184 A; WONZIJ, YAZCOV F-F ~0.84 A), all three already rejected by other rules. So the hole is real but unexercised there; adding `duplicate-atom` to those two would be behaviour-neutral on these collections and is a consistency follow-up.
- Combined position across the three domains at this commit: COF **37/884 (4.2%)**, MOF **76/5226 (1.45%)**, macromolecule **0/2**. The rule sets are NOT interchangeable - COF thresholds reject 8 of 12 CR MOFs as false positives, and MOF or COF thresholds would fire `bare-carbon` on every carbon of a raw PDB - which is why this is a hook with three independent overrides rather than one implementation.

## [2026-09-29] Parity repair respects quaternary ammonium N+ (code change, see project-decisions.md)
- Full 884 run on ARF (job 6423477, `b4e6976`): 37 pre-scan rejections, 1649 frames, odd 0, multi_piece 0, `over_coordinated` 8 - all quaternary N+ in `1062`, `671`, `673`, `674`. Adding a checker exemption alone was considered and rejected: it would have hidden a real defect, a hydrogen the parity repair put on a carbonyl O because it read the closed-shell cation as a radical.
- `fragmentation_oop.py`: module helpers `_is_quaternary_ammonium` / `_quaternary_ammonium_count`; base hook `_formal_charge` (0) with a COF override; parity gate in `fix_odd_electron_multiplicity` uses `Z_sum - charge`; COF output quarantine subtracts the same count; `_write_extxyz` writes `charge=` when nonzero, `_parse_extxyz` reads it, `_update_extxyz_collection` preserves it.
- `runUniFrag/check_cof_fragments.py`: reads `charge=` (absent = 0), exempts N+(C)4 from `over_coordinated` via the imported helper, judges `odd_electron` on `zsum - charge`.
- Validation: see the decision entry - 8 N+ fragments changed, 147 byte-identical over 48 COFs, MOF unchanged apart from the ABAYIO/ABAYOU ownership swap.
- Follow-up: other formal-charge motifs (pyridinium, imidazolium, spiroborate) are unhandled. Checking the swap needed a geometry comparison because the MOF content-multiset test counts it as a difference; a regression script that ignores ownership of a shared linker would save re-deriving this each time.

## [2026-09-29] Full 884 validation of the quaternary N+ parity fix (no code change)
- ARF smp job 6424004 at `cfc564a`, same 884 CIFs and script as job 6423477 (`b4e6976`), output `/arf/scratch/otayfuroglu/unifrag_884_v3`. 7.5 min on 64 cores.
- MUST: all passed. Pre-scan 37 (same list), timeout 1 (913), frames 1649, `over_coordinated` 8 -> 0, `bad_terminal_caps` 64 -> 56, clashing 9, odd 0, multi_piece 0, odd-electron quarantine 23 (all identical).
- Exactly 8 frames changed content: `FragCof`/`FragCofMin` of 1062, 671, 673, 674, the only frames tagged `charge=1`. The 8 fewer bad caps are those same fragments - the spurious carbonyl C-OH was also counted as a mis-capped oxygen (the pre-run prediction that bad caps would stay at 64 was wrong). No new bad caps. 1610 frames byte-identical.
- 31 frames changed owner, not content: for 17 duplicate-parent pairs (106/107, 1147/1148, ...) the duplicate filter kept the other copy, decided by which worker finishes first. None of the 62 frames contains an N+, so the fix cannot reach them.
- Pre-existing, found through the swap, not fixed: two "duplicate" pairs differ in hydrogen count - `554FragCofMin` C36H38N8O6 vs `550FragCofMin` C36H36N8O6, and `639FragCof` C36H30N12O12 vs `640FragCof` C36H24N12O12. The duplicate key appears not to compare hydrogens, so which formula enters the collection depends on run order. Follow-up: include H in the key or make the duplicate choice deterministic, and check which parent of each pair is mis-hydrogenated.

## [2026-10-01] Fragments are always charge 0: charge handling from cfc564a reverted (code change, see project-decisions.md)
- Trigger: OpenMLP's ORCA runs on `subset_qm250` failed for 671/673/674/1062 FragCof (`*xyz 0 1`, 315 electrons, odd). UniFrag had tagged them `charge=1` and ASE read it into `atoms.info`, but OpenMLP hard-codes `charge=0, mult=1`. The user does not want charged fragments and chose to revert to adding an H.
- `fragmentation_oop.py`: restored from `b4e6976`, plus `_is_quaternary_ammonium` only (used by the checker). Removed `_quaternary_ammonium_count`, the `_formal_charge` hook and override, the `charge=` write/parse/preserve, and the charge-aware COF quarantine.
- `runUniFrag/check_cof_fragments.py`: restored from `b4e6976`, plus the N+(C)4 exemption from `over_coordinated`. Charge parsing removed; `odd_electron` is `zsum % 2` again.
- Validation on ARF (1062, 671, 673, 674, 165), against a control running `b4e6976` code on the same 5 CIFs: 15 of 17 frames byte-identical, 0 content differences, 0 `charge=` tags; the other 2 are a library-fragment ownership swap between 1062 and 671 (duplicate filter, pre-existing). Checker: MUST all passed, odd 0, over_coordinated 0, bad caps 8 = the expected C-OH on the 8 N+ fragments.
- Caution for anyone comparing runs: running a subset instead of the full 884 changes coordinates of shared fragments (same formula, different source parent), so always compare against a control run on the same input set.
- The swap also showed `1062FragCofOnlyNode_0` (C9H15N3) and `671FragCofOnlyNode_1` (C9H9N3) treated as duplicates - a third case of the duplicate key ignoring hydrogens (see the 550/554, 639/640 entry).

## [2026-10-01] Underscore-free labels and an always-written charge tag (code change, see project-decisions.md)
- Trigger: OpenMLP ORCA failure on `1016FragCofOnlyLinker_1` ("Input geometry does not match current geometry") - it restarted from `initial_1016FragCofOnlyLinker.gbw`, written for `_0`, because OpenMLP names restart files by `label.split("_")[0]`.
- `fragmentation_oop.py`: the six `_{idx}` suffix sites (COF linker/node, MOF linker; internal labels and extxyz names) now write `{idx}`. `_write_extxyz` takes `charge=0` and always writes `charge=<q>` right after `label=`; `_parse_extxyz` reads it (absent = 0) and `_update_extxyz_collection` preserves it.
- `runUniFrag/check_cof_fragments.py`: reads `charge=` (absent = 0), `odd_electron` judged on `zsum - charge`.
- Validation: on 1016 + 1062 (ARF) every frame has `charge=0`, labels `1016FragCofOnlyLinker0/1`, `...OnlyNode0/1`, no `_` anywhere; ASE reads `info["charge"]=0` on all 11 frames (the first ordering, with charge after `capped_h`, lost it on the one frame with an empty `capped_h=""`); checker MUST all passed. MOF (ABAVIJ, ABEDAQ): `charge=0` on all 6 frames. The 12 CR MOFs: 13/14 frames byte-identical to `HEAD` apart from the new header field, the 14th the ABAYIO/ABAYOU race; none of the 12 has more than one linker, so the MOF rename itself is not exercised there.
- Pre-existing quirk noticed, not changed: on frames with an empty `capped_h=""`, ASE swallows `pbc="F F F"` into `capped_h`. Harmless (pbc defaults to false), but any ASE consumer of `capped_h` should not trust it on those frames.
- Also noticed: `_update_extxyz_collection` drops existing frames by `label.startswith(clean_base)`, so single-file mode for `10.cif` would also drop every `100...`/`1016...` frame. Not touched; folder mode does not use it.

## [2026-10-01] Validation of 32f2d3d: label-only control and full 884 (no code change)
- Label-only control (ARF job 6430483): 51 COFs that have multi-piece library fragments, run single-worker with `32f2d3d`-era code and with an identical copy keeping the old `_{idx}` suffix (6 lines differ). After renaming `_N` -> `N`: collection 190/190 and quarantine 3/3 frames identical, CSV identical. The rename changes names only.
- Full 884 (ARF job 6430525, output `/arf/scratch/otayfuroglu/unifrag_884_v4`): MUST all passed; pre-scan 37, timeout 1 (913), frames 1649, odd 0, over_coordinated 0, multi_piece 0, clashing 9. All 1649 frames carry `charge=0`; no label contains `_`.
- Against job 6423477 (`b4e6976`, the same chemistry): **0 frames changed content**, 1606 identical by label (after the `_N` rename), 43 changed owner within duplicate pairs. Against `cfc564a` (6424004): exactly the 8 N+ fragments differ, as intended.
- bad caps 65 (not the 64 of `b4e6976`), quarantine 22 (not 23), produced 637 (not 638) - all from the duplicate filter choosing the other copy of two pairs whose formulas differ: `429` C98H61N8 vs `430` C98H64N8 (this run kept 430, so 429's odd FragCof was never produced) and `637` C36H28N12O12 vs `638` C36H18N12O12 (kept 637, whose FragCof/FragCofMin are mis-capped).
- Duplicate filter: this makes five pairs treated as duplicates despite different hydrogen counts (550/554, 639/640, the 1062/671 library node, 429/430, 637/638), and which chemistry reaches the collection depends on worker finishing order, so summary counts move between identical-code runs. Highest-value follow-up for reproducibility: include H in the duplicate key, and/or choose deterministically among true duplicates (e.g. lowest stem).


## [2026-10-01] Parent-aware capping of cut sigma-complete N/O (code change, see project-decisions.md)
- Trigger: OpenMLP SCF failures on `623FragCofMin`, `636FragCof`, `638FragCof` (subset250). `623FragCofMin` has four N-H radicals (its parent N is -NH- with C 1.35/1.43 A; the cut C was never replaced). `636`/`638` are a different cause (near-degenerate pi system, imine caps) and are not touched here.
- `fragmentation_oop.py`: base hook `_parent_unfilled_sigma_sites` (None); COF helpers `_parent_sigma_neighbours`, `_neighbour_present`, `_parent_unfilled_sigma_sites`, `_is_parent_cut_sigma_site`; `_cap_open_oxygens` consults the hook before the carries-hydrogen gate; veto added to COF `_can_remove_cap_h` and to Tier B of `_heteroatom_parity_repair`.
- `runUniFrag/check_cof_fragments.py`: `bad_terminal` flags are dropped when the atom has exactly its sigma-complete parent's coordination (needs `<stem>.cif` beside the extxyz).
- Validation: full 884 on ARF smp. Control (HEAD code, job 6431769) vs v4 (HEAD, job 6430525): 46 ownership swaps, 0 same-label content changes, so the control is noise-free for same-label comparison. Fixed code, final (job 6431802, `/arf/scratch/otayfuroglu/unifrag_884_v6_undercap`): MUST passed, odd 0, over 0; frames 1644 (control 1649), quarantine 27 (control 23); bad terminal caps 42 under the parent-aware checker (control 57). Parent-aware scan of every N/O that is sigma-complete in its parent and lost a heavy neighbour: 296 under-capped sites in 84 frames -> 0, 0 over-capped. Standalone: `623FragCofMin` 178 -> 182 atoms, C90H76N10O6, even.
- Not what the first version did: v5 (job 6431738) left 5 sites because Tier B of the parity repair stripped the new H again (752 `QM-Fix [cap-OH]`). The veto there is what closes it, and it moves 752FragCof, 752FragCofMin, 160FragCofOnlyLinker, 929FragCofMin to quarantine (odd electron count).
- Three Min frames changed heavy atoms: `1014`, `482`, `492` FragCofMin. Their first Min used to be too big (`Minimize only reduced size by -59 atoms (< 20)`) and fell back to node + one linker; with two more H on the normal fragment the first Min now reduces by exactly 20 and is accepted. The 20-atom rule is a threshold, so Min size can jump when H counts move.
- Follow-ups: (1) `752`/`160` have even parents but odd fragments once capped right - there is a second defect the old parity hack hid; not investigated. (2) `_try_cof_graph_node_linker_fragment` still has its own O/B-only capping (see the 2026-09-29 N/O capping decision). (3) The duplicate key still ignores H, so which of two parents with different H counts is emitted depends on worker order (550/554, 639/640, 637/638, 429/430). (4) The QM subset for OpenMLP must be rebuilt: 95 frames changed, 4 left for quarantine.
