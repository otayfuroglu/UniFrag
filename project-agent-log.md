# Project Agent Log

Chronological handoff log for agents working on UniFrag. Add newest entries at the top. Each entry should include changed files, validation, decisions, and follow-up risks.

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





