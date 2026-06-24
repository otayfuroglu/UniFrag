# Project Memory Template (Reusable)

Use this file as a persistent engineering memory for this project. Keep entries concise, date-stamped, and actionable.

> Coordination note: agents should also read `project-decisions.md` for durable decisions and `project-agent-log.md` for chronological handoffs.

## 0) Metadata

- Project name:
- Repository URL:
- Default branch:
- Primary language(s):
- Last updated: 2026-05-11
- Maintainer(s):

## 1) Project Purpose

### 1.1 Goals
- 
- 

### 1.2 Non-goals
- 
- 

### 1.3 Current scope
- 

## 2) Architecture

### 2.1 Top-level structure
- `path/or/module`: purpose
- `path/or/module`: purpose

### 2.2 Core components
- Component/Class/Service:
  - Responsibility:
  - Inputs:
  - Outputs:
  - Key invariants:

### 2.3 Data flow (high-level)
1. 
2. 
3. 

## 3) Coding Conventions

### 3.1 Style and organization
- 
- 

### 3.2 Error handling
- 

### 3.3 Logging/observability
- 

### 3.4 Performance conventions
- 

## 4) Dependencies and Environment

### 4.1 Runtime dependencies
- Package: why needed
- Package: why needed

### 4.2 Tooling
- Formatter/linter:
- Test framework:
- Build system:

### 4.3 Environment assumptions
- Python/Node/Compiler version: Python 3.13 (Miniconda)
- OS notes: macOS, requires setting CSD_DATA_DIRECTORY environment variable for CCDC Python API
- GPU/CPU notes:

## 5) Commands (Copy/Paste)

### 5.1 Setup
```bash
# Export the CSD database directory for CCDC Python API
export CSD_DATA_DIRECTORY="/Users/omert/CCDC/ccdc-data/csd"
```

### 5.2 Run
```bash
# Run the CIF extraction script
CSD_DATA_DIRECTORY=/Users/omert/CCDC/ccdc-data/csd /Users/omert/miniconda3/bin/python runUniFrag/fetch_cifs_from_csd.py

# Analyze the extracted CIFs and classify them into CR (ASR/FSR) and NCR subsets
/Users/omert/miniconda3/bin/python runUniFrag/analyze_cifs.py

# Extract Zn-based unmodified CIFs and merge them into categories (removing duplicates)
/Users/omert/miniconda3/bin/python runUniFrag/merge_zn_cifs.py

# Collect all unique Zn-based MOFs into a single deduplicated folder (naming files [REFCODE].cif)
/Users/omert/miniconda3/bin/python runUniFrag/collect_all_zn_cifs.py

# Collect all unique CR MOFs into a single deduplicated folder (prioritizing modified files)
/Users/omert/miniconda3/bin/python runUniFrag/collect_cr_cifs.py

# Extract Zn-based CR MOFs from the unique CR collection
/Users/omert/miniconda3/bin/python runUniFrag/collect_zn_cr_cifs.py

# Run atom types coverage analysis and compare parent MOFs vs fragments library
CSD_DATA_DIRECTORY=/Users/omert/CCDC/ccdc-data/csd /Users/omert/miniconda3/bin/python runUniFrag/analyze_atom_types.py

# Run Zn coordination environment analysis and compare parent MOFs vs fragments library
CSD_DATA_DIRECTORY=/Users/omert/CCDC/ccdc-data/csd /Users/omert/miniconda3/bin/python runUniFrag/analyze_zn_coordination.py

# Extract parent Zn MOFs with Zn CN=0, 1, or 3 to a separate CSV file
CSD_DATA_DIRECTORY=/Users/omert/CCDC/ccdc-data/csd /Users/omert/miniconda3/bin/python runUniFrag/extract_low_coordination_mofs.py

# Plot and update comparative metal distribution histogram for CSD-modified and CSD-unmodified subsets
/Users/omert/miniconda3/bin/python runUniFrag/plot_metals.py

# Plot and update metal distribution histogram for the unique CR collection
/Users/omert/miniconda3/bin/python runUniFrag/plot_cr_metals.py

# Collect original unmodified CIFs for standard CSD entries in the 8806 screening list
CSD_DATA_DIRECTORY=/Users/omert/CCDC/ccdc-data/csd /Users/omert/miniconda3/bin/python runUniFrag/collect_screening_mofs.py

# Collect modified CIFs (ASR/FSR/Ion/NCR) matching the 8806 screening list from CSD-modified
/Users/omert/miniconda3/bin/python runUniFrag/collect_modified_screening_mofs.py
```

### 5.3 Test
```bash
# Run MOF family smoke tests
./run_mof_family.sh test_on_irmof_series 4.0

# Run COF family smoke tests
./run_cof_family.sh test_on_cof_zn_pc_series 4.0
```

### 5.4 Lint/format
```bash
# lint/format
```

### 5.5 Release/deploy
```bash
# release/deploy
```

## 6) Development Workflow

### 6.1 Branching and PR flow
1. 
2. 
3. 

### 6.2 Change checklist
- [ ] Reproduce baseline behavior
- [ ] Add/adjust tests
- [ ] Validate key scenarios
- [ ] Update docs/memory

### 6.3 Review focus areas
- 
- 

## 7) Testing Strategy

### 7.1 Test tiers
- Unit:
- Integration:
- End-to-end/manual:

### 7.2 Critical regression cases
- Case name: expected outcome
- Case name: expected outcome

### 7.3 Test datasets/fixtures
- `path/to/fixtures`: purpose

## 8) Known Issues and Risks

For each item:
- ID:
- Symptom:
- Root cause (known/suspected):
- Affected area:
- Workaround:
- Status:

## 9) Decisions Log (ADR-lite)

### Decision 2026-05-28: Heavy-Atom Formula-Only Chemical Identity Key for Robust Conformer Deduplication
- Context: The user noticed that chemically equivalent structures (same formula and connectivity, but with minor conformational/torsional variations or hydrogen-capping differences) were not flagged as duplicates under the collection-level deduplication. Originally, deduplication relied on the strict `_species_coords_unique_key` (which encoded exact coordinates and pairwise distances), meaning conformers were treated as unique.
- Decision: Implemented a robust chemical identity key helper (`_chemical_identity_key`) for collection-level duplicate detection. Instead of using exact coordinates or pairwise distance histograms (which suffer from binning issues due to small lattice parameter shifts across CIF files), the new key filters the fragment down to heavy atoms (excluding Hydrogen atoms entirely) and constructs a sorted count tuple of only heavy element occurrences (e.g. `(("C", 12), ("O", 4), ("Mg", 2))`).
- Consequences: All collection-level duplicate checks across MOF, COF, and Bio modes (including sliding windows and batch directory sweeps) correctly identify and group chemically equivalent structures/conformers as duplicates. This ensures that only unique chemical topologies are saved to the `.extxyz` collections, while duplicates are correctly documented in the summary `.csv` files as duplicates.
- Alternatives considered: Pairwise distance histograms binned at `0.1 A` resolution; rejected because minor lattice shifts across CIFs caused systematic deviations that defeated binning.

### Decision 2026-05-11: COF helper libraries use global chemical duplicate pruning
- Context: COF node/linker helper folders are used for visual QA, and duplicate building blocks can appear under different COF stems. Exact/near-exact coordinate matching was too strict for chemically identical blocks such as COF-LZU8 nodes.
- Decision: COF helper export prunes/checks duplicates across the whole `cof_nodes_lib/` or `cof_linkers_lib/` folder using composition plus internal pair-distance fingerprints rounded to `0.1 A`. This applies globally to every COF helper export before new node/linker files are written.
- Consequences: Helper libraries keep one representative for chemically identical COF building blocks across stems; intentionally similar conformers may be merged.
- Alternatives considered: `0.01 A` duplicate matching; rejected because chemically identical helper blocks remained duplicated.

### Decision YYYY-MM-DD: Title
- Context:
- Decision:
- Consequences:
- Alternatives considered:

### Decision 2026-05-08: Restore ZnPc/metallo-PC COF Path J before generic COF paths
- Context: ZnPc-DPB was using generic COF Path B instead of the accepted direct coffragmentor Path J, causing minimized dimer imbalance.
- Decision: COF extraction should try direct coffragmentor metallo-PC Path J first: combine the Zn/N-rich node with attached linker images and duplicate the full set along the shortest lattice vector for the dimer.
- Consequences: ZnPc-DPB normal/min are restored to 322/166 atoms; minimized has two identical 83-atom layers (`Zn1 B2 C48 H16 N8 O8` each).
- Alternatives considered: Fixing generic Path B for ZnPc; rejected because Path J is the established ZnPc rule.

### Decision 2026-05-08: ZIF/minimum fragments preserve the first linker ring
- Context: ZIF linkers can have first rings that are pentagons/heterocycles, not six-carbon phenyl rings. Minimum trimming must not cut that first connected ring. Also, when the normal fragment is already small, a separate minimum is unnecessary.
- Decision: MOF Path J first-ring trimming now detects the nearest heavy-atom cycle of size 3-8 and preserves the whole ring, including hetero atoms. For any MOF minimized extraction, first generate/probe the normal fragment; if the normal fragment has fewer than 80 atoms, write that normal fragment as the minimized output.
- Consequences: ZIF generated min outputs with normal sizes under 80 atoms are identical to normal outputs (e.g. ZIF-1/ZIF-10/ZIF-2 at 33 atoms, ZIF-11/ZIF-12 at 56 atoms, ZIF-20 at 49 atoms). Larger MOFs still use the first-ring minimum logic.
- Alternatives considered: Keep six-carbon-ring-only trimming; rejected because it can truncate imidazolate/pentagon rings.

### Decision 2026-05-08: DUT-49 normal uses one Cu2 paddlewheel node
- Context: DUT-49 has a Cu paddlewheel structure where the normal fragment should follow Cu-BTC behavior, not PCN/NU behavior. The old saved normal output had `Cu3` and 547 atoms.
- Decision: Treat `DUT-49*` as a one-node paddlewheel family alongside `Cu-BTC*`; normal output should keep one Cu2 node. Minimized Path J behavior remains unchanged.
- Consequences: `test_on_other_common_mofs/DUT-49_frag.xyz` is regenerated as 370 atoms with `Cu2`; minimized remains 136 atoms with `Cu2`. PCN/NU still use two-node atom-count selection.
- Alternatives considered: Apply PCN/NU two-node selection to DUT-49; rejected by visual inspection because one paddlewheel node is enough.

### Decision 2026-05-07: PCN/NU normal fragments keep two metal nodes
- Context: PCN and NU families can have linker structures where the chemically meaningful normal fragment spans two metal nodes, opposite to Cu-BTC and DUT-49 where one Cu2 paddlewheel node is sufficient.
- Decision: Keep the Cu2 one-node suppression only for one-node paddlewheel families such as `Cu-BTC*` and `DUT-49*`. For PCN/NU normal fragments, skip Path J and use the legacy Path C candidate search to choose the smallest metal-complete two-node fragment. Minimized mode remains one-node/one-linker focused and can still use Path J.
- Consequences: PCN/NU normal Cu paddlewheel cases now report `Cu4` or `Zn4` as appropriate, and NU-108-Cu selects 516 atoms instead of 912; `test_on_cubtc/Cu-BTC.cif` normal remains `Cu2` and 90 atoms, and DUT-49 remains one-node `Cu2`.
- Alternatives considered: Apply the Cu2 suppression to every Cu paddlewheel family; rejected because PCN/NU need the opposite normal-fragment topology.

### Decision 2026-05-07: Cu-BTC/DUT-49 normal uses one Cu2 paddlewheel node
- Context: Cu-BTC uses the legacy MOF path. Path C duplicated the discrete Cu2 paddlewheel in normal mode, but one metal node is the desired normal fragment.
- Decision: For `Cu-BTC*` and `DUT-49*` discrete `Cu2` SBUs, skip normal Path C second-SBU expansion and keep one Cu paddlewheel node.
- Consequences: Cu-BTC normal is now `C36 Cu2 H28 O24` (90 atoms); minimized remains `C30 Cu2 H22 O12` (66 atoms).
- Alternatives considered: Generic two-SBU Path C; rejected for Cu-BTC visual correctness.

### Decision 2026-05-07: MOF Path J open-connector recovery
- Context: Some P1/helper-fragmented MOFs can leave an open node carboxyl carbon because a real linker branch exists in the CIF but is not returned as a metal-attached helper linker image.
- Decision: MOF Path J should recover those missing branches from the original CIF structure after helper node+linker assembly. Normal keeps the full recovered branch; minimized keeps its first connected ring and caps unsaturated carbons.
- Consequences: IRMOF-11 is fixed at `Zn4 C108 H84 O25` normal (221 atoms) and `Zn4 C53 H39 O15` minimized (111 atoms); no open node carboxyl carbons remain in either output.
- Alternatives considered: Wider helper image search; rejected because it did not reveal the missing IRMOF-11 linker.

### Decision 2026-05-07: MOF Path J minimum linker context
- Context: Minimum MOF fragments should preserve one full linker and still show first-ring context for the other node attachment directions.
- Decision: Minimized MOF Path J keeps one full attached linker image and adds the first connected six-membered carbon ring, with bonded hydrogens, from every other attached linker image. Node-side connector atoms supplied by the helper node are reused. Retained first-ring/connector carbons that lost heavy neighbors during trimming are H-capped.
- Consequences: IRMOF-1 minimum is now `Zn4 C43 H31 O15` (93 atoms), not the earlier 38-atom one-linker-only model or the uncapped 88-atom first-ring model.
- Alternatives considered: One full linker only; rejected because other node sides lacked first-ring context.

### Decision 2026-05-07: Compact unique MOF helper library names
- Context: MOF helper node/linker exports are manually inspected, so filenames should stay short and duplicate fragments should be suppressed.
- Decision: `mof_nodes_lib/` and `mof_linkers_lib/` use per-folder sequential names like `IRMOF-1_00.xyz`. Export keeps only unique molecules using composition and rounded internal pair-distance fingerprints. If same-stem files already exist in a helper folder, export for that folder is skipped instead of overwriting. For new stems, duplicate checks scan every existing `.xyz` in the target node/linker folder, so duplicate molecules are skipped even across different MOFs.
- Consequences: Do not depend on composition/smiles metadata in helper filenames; inspect the XYZ contents or folder type for details. Existing visual-check helper files are preserved, and helper libraries remain globally unique per folder.
- Alternatives considered: Verbose composition/smiles filenames; rejected as too noisy for visual review.

### Decision 2026-05-07: Global capped-H-only final geometry cleanup
- Context: Final H cleanup must be consistent regardless of which MOF/COF path generates the fragment.
- Decision: Extraction finalization should only move H atoms that were explicitly added as caps. Heavy atoms and non-capped hydrogens are not optimized or moved by the final cleanup step; no extraction path should call full-molecule RDKit/UFF refinement.
- Consequences: Use `optimize_capped_h_geometry_only(...)` for final C-H/O-H cap cleanup. `refine_h_geometry_with_rdkit(...)` is legacy/unused unless a future decision explicitly reintroduces it.
- Alternatives considered: RDKit/UFF minimization with fixed heavy atoms; rejected because the project invariant is capped-H-only movement.

### Decision 2026-05-04: Keep tied Path B layered neighbors
- Context: COF-202 visually missed one linker when Path B selected only one of two equally close B/O node components at 2.85 A.
- Decision: Path B now keeps all same-spacing layered neighbor components within 0.05 A of the best spacing.
- Consequences: COF-202 normal fragment grows from 124 to 169 atoms; minimized fragment grows from 58 to 70 atoms. COF-102 and COF-300 smoke checks remain on Path A/Path C respectively.
- Alternatives considered: Keep only the first nearest layered component; rejected because tie ordering dropped a symmetry-equivalent linker.

### Decision 2026-05-04: Prioritize metallo-PC cores before B/O layered nodes
- Context: ZnPc-COF contains B/O linker nodes, so the generic Path B detector ran before the porphyrin-like N-rich SBU detector and selected a layered B/O fragment.
- Decision: Add Zn bonding radius and detect N-rich metal macrocycles before generic B/O Path B when Zn/Cu/Fe/Co/Ni/Mn are present.
- Consequences: ZnPc-COF now uses Path D (Metallo-PC core dimer). Normal keeps the two stacked ZnPc layers plus all connected benzene-1,4-diboronic acid linkers and cuts only far-side O-C bonds so far B-connected O atoms are H-capped: 242 atoms, Zn2 B16 H64 C112 N16 O32. Minimized keeps the two stacked ZnPc layers plus one full BDBA linker per layer, retains the ZnPc-core fused benzene perimeter, caps discarded linker attachment points with H, and uses the same far-side O-H termination: 146 atoms, Zn2 B4 H40 C76 N16 O8. COF-366 remains Path D porphyrin; COF-202 remains Path B layered set.
- Alternatives considered: Leave ZnPc as Path B; rejected because the intended SBU is the metallo-phthalocyanine core.

### Decision 2026-05-06: Directly combine coffragmentor node and linker for metallo-PC minimum fragments
- Context: The failed helper-selection approach was undone. `coffragmentor.py` should provide the actual node/linker molecular fragments for ZnPc minimum fragments.
- Decision: Metallo-PC COFs now try Path J: combine one Zn/N-rich coffragmentor node molecule with coffragmentor linker molecule image(s), then duplicate the pair/set along the shortest lattice vector for the ZnPc dimer. Minimized mode keeps one nearest attached linker image; normal mode keeps all neighboring-cell linker images that chemically attach to the node.
- Consequences: ZnPc-DPB normal/min are 322/166 atoms; ZnPc-COF normal/min are 210/138 atoms. Non-metallo-PC COFs fall back to existing paths.
- Alternatives considered: coffragmentor-assisted UniFrag linker selection; reverted.

### Decision 2026-05-06: Use moffragmentor node+linker combine as first MOF Path J
- Context: Start generalizing the node+linker approximation to MOFs using `moffragmentor`, analogous to `coffragmentor.py` for COFs.
- Decision: MOF extraction first tries `moffragmentor` Path J, exports helper fragments to `mof_nodes_lib/` and `mof_linkers_lib/`, merges overlapping same-element boundary atoms, H-caps open terminal linker oxygens and locally adjusts only capped H geometry without full-molecule RDKit/UFF optimization, and falls back to legacy MOF paths if helper fragmentation is unavailable.
- Consequences: IRMOF-1 normal is `Zn4 C48 H36 O25` (113 atoms); min is `Zn4 C13 H6 O15` (38 atoms).
- Alternatives considered: Legacy MOF graph path first; kept only as fallback on this branch.

## 10) Operational Notes

### 10.1 Common failure modes
- 

### 10.2 Recovery steps
1. 
2. 

### 10.3 Permissions/secrets notes
- 

## 11) Migration / Reuse Notes

### 11.1 What to copy to a new project
- 

### 11.2 What must be customized
- Paths/directories
- Dependency versions
- Environment-specific commands
- Domain-specific heuristics

## 12) Quick Start for New Contributor

1. Read sections: 1, 2, 5, 7, 8.
2. Run setup and one smoke test.
3. Validate one known critical case.
4. Make small change and run regression checks.

## Appendix A) Project-Specific Filled Example (Optional)

Use this section only if you want this template to include a concrete reference snapshot.



### 2026-05-08 COF metallo-PC layer output
- `fragmentation_oop.py --kind cof` accepts `--cof-layer auto|monomer|dimer` for layered COFs.
- For metallo-PC Path J, `monomer` writes one coffragmentor node+linker layer and `dimer` duplicates along the shortest lattice vector.
- For generic Path B layered COFs, `monomer` keeps one principal layer and `dimer`/`auto` keeps two principal layers.
- `run_cof_family.sh <folder> [radius] both` writes separate `_monomer` and `_dimer` normal/minimized outputs for visual comparison.

### 2026-05-08 COF Path B normal boundary
- Generic layered COF Path B normal mode keeps terminal neighboring B/O node components at linker cut sites. This fixed COF-10 normal dimer missing terminal chemistry.
- Minimized Path B still trims as before.
- COF-10 auto normal/min now: 112/66 atoms; explicit normal monomer/dimer: 56/112 atoms.

### Decision 2026-05-09: COF layered routing by topology similarity (node first, linker second)
- Context: Filename-based routing caused regressions for COF-6/66/8 layered dimers.
- Decision: Route COF strict/fallback paths by node/linker topology signatures rather than structure stem names; keep node+linker-first, then fallback.
- Consequences: Better transfer across COF families, fewer structure-name special cases.


- 2026-05-11: COF Path J layered dimer rule added. For coffragmentor node+linker COFs, non-monomer output now duplicates the completed capped fragment along the shortest lattice vector when that vector is a plausible face-to-face layer spacing (2.5-5.0 A). Validated on COF-LZU1 (3.729 A) and COF-LZU8 (4.093 A), normal and min dimers.
