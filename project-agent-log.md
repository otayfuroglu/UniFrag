# Project Agent Log

Chronological handoff log for agents working on UniFrag. Add newest entries at the top. Each entry should include changed files, validation, decisions, and follow-up risks.

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
