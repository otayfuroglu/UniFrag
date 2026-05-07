# Project Agent Log

Chronological handoff log for agents working on UniFrag. Add newest entries at the top. Each entry should include changed files, validation, decisions, and follow-up risks.











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
