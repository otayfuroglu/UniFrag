# Project Decisions

Durable implementation and architecture decisions for UniFrag. This file is the source of truth for decisions; keep entries concise, dated, and actionable.

## Decision 2026-05-24: Automated QM-readiness validation, hydrogen saturation, neutralization engine, and ring planarity gating
- Context: For seamless Quantum Chemistry (QM) calculations on the generated fragments, it is essential that all outputs are formally neutral (charge = 0) and closed-shell singlet (spin multiplicity = 1). Relying on manual or external tools to saturate radicals, protonate coordinating anions, or deprotonate acidic groups is error-prone. Additionally, blindly applying aromatic SVD planarity fits to acyclic/peptide backbone carbons distorts local capping geometries and leads to severe under-coordination.
- Decision: Implemented a robust validation helper `QMReadinessChecker` and an automated iterative engine `force_qm_readiness` in `fragmentation_oop.py`. The engine dynamically neutralizes fragments by adding hydrogens to anionic coordinating sites (e.g. carboxylates, phenoxides) or removing hydrogens from acidic groups. It also saturates under-coordinated light non-metals (C, N, O, B) to resolve radical states, iterating up to 10 cycles until exactly `charge = 0` and `multiplicity = 1` are achieved. Furthermore, `enforce_sp2_capped_h_geometry` was upgraded with an inline BFS cycle/ring checker (`is_in_ring`) to strictly gate the global aromatic SVD planarity plane fitting to cyclic/ring systems only, leaving non-cyclic/peptide backbone cap geometries physically realistic.
- Consequences: All exported frames in the ExtXYZ collections (`fragments_collection.extxyz` and `bio_fragments_collection.extxyz`) are guaranteed to be closed-shell singlets with exactly 0 formal charge, completely ready for QM calculations out-of-the-box. Console diagnostics and metadata embedding are preserved, and peptide/acyclic cap geometries are robustly protected against distortions.


## Decision 2026-05-22: Unified Single-File/Folder Execution with Atomic Collection Updates and No Individual XYZ Files
- Context: Previously, running fragmentation on a single file produced separate `.xyz` files and wiped the target collections, whereas directory execution appended frames/rows. The user requested completely eliminating individual `.xyz` files (relying only on `fragments_collection.extxyz` or `bio_fragments_collection.extxyz`) and ensuring that re-running a single file atomically updates its specific entries in the CSV and ExtXYZ collections without affecting other files in the same folder.
- Decision: Unified both single-file and directory modes in `main()` so both write to/update the same central collections in the parent directory. Implemented two atomic, temp-file-safe transactional update helpers: `_update_csv_rows` and `_update_extxyz_collection`. When a structure is processed, any existing entries matching that filename/label are safely deleted and replaced with the new results in both the summary CSV and the ExtXYZ collection. Also, `BioMolFragmenter.extract` was extended with `write_files=False` so that individual window `.xyz` files are completely bypassed.
- Consequences: No individual `.xyz` files are written under any mode (MOF, COF, or Bio). Centralized collections (`fragments_collection.extxyz` and `fragmentation_summary.csv`) serve as the clean, single-point of output, while re-runs remain highly efficient and preserve existing directory progress.

## Decision 2026-05-22: Use label instead of name in ExtXYZ headers
- Context: The user requested replacing `"name"` with `"label"` in the header of the generated `fragments_collection.extxyz` file.
- Decision: Updated the Atoms `info` dictionary key from `"name"` to `"label"` across all three fragmentation modes (MOF, COF, Bio) in `fragmentation_oop.py`. This ensures that the ASE-written ExtXYZ comments consistently serialize as `label=...` rather than `name=...`.
- Consequences: Downstream ExtXYZ outputs automatically write the fragment identifier with a `label` key instead of a `name` key, matching the requested specification.

## Decision 2026-05-14: Unified batch folder processing architecture applies to all three modes (MOF, COF, Bio)
- Context: Batch folder processing, parallel `--nproc`, incremental CSV and ExtXYZ outputs were initially only implemented for the MOF mode. COF and bio modes were still single-file only.
- Decision: Extended the identical architecture to COF (`*.cif` folder) and bio (`*.pdb` folder) modes. Each mode now has a dedicated top-level worker function (`_process_cof_file`, `_process_bio_file`) that is safe to call from a `multiprocessing.Pool`. CSV and ExtXYZ collection files are initialized empty at the start of a batch run and appended to incrementally as each worker returns. Bio CSV specifically writes one row per sliding-window fragment (columns: `pdb_file`, `window_name`, `n_atoms`) so individual windows are identifiable without parsing the ExtXYZ file.
- Consequences: All three modes share a consistent CLI surface (`input_path` can be a file or folder, `--nproc` always applies) and identical output conventions (incremental CSV + ExtXYZ, bracket-free fragment names, no `.xyz` suffix in ExtXYZ names).

## Decision 2026-05-13: MOF Linker Minimization utilizes topological skeleton pruning and bridge detection
- Context: The user noticed that for Mg-based MOFs (like MOF-74) and linkers with complex functional groups, the `--minimize` mode was inappropriately cutting atoms *inside* the aromatic rings. This was due to two flawed algorithms: an arbitrary BFS depth limit (`depth < 5`) in the native UniFrag path, and a hardcoded 6-membered Carbon-only ring finder in the Path J `moffragmentor` fallback which completely obliterated heterocycles like thiophene. Furthermore, redundant linkers in long organic chains were incorrectly saving all connected cycles.
- Decision: Completely rewrote the `minimize` trimming algorithms across all paths. They now use a mathematically perfect Topological Skeleton (2-core) Pruning algorithm. It starts with the entire linker molecule and iteratively deletes ANY atom (regardless of element) that has a connectivity degree of 1, UNLESS that atom is explicitly coordinating to a metal (`bridge_atoms`). For redundant/dangling linkers, it now explicitly maps "bridges" within the 2-core and truncates the structure immediately after the first biconnected ring system is isolated.
- Consequences: This guarantees mathematically perfect topological skeletons. The script effortlessly strips away all dangling functional groups (e.g. uncoordinated carboxylates, -OH, -NH2, -CH3) from the rings, perfectly preserves the entire cyclic ring structure (even non-carbon heterocycles like thiophene), and strictly truncates dangling linkers to just their very first ring system.
- Update: Enforced a new rule that `--minimize` is now handled automatically. The Python script will **automatically abort** and only return the normal un-minimized fragment if the total atom count of the normal fragment is predicted to be strictly less than 50 atoms. If the size is >= 50, it automatically generates both the normal and `_min` versions.
- Architecture Change: Shifted from relying on `run_mof_family.sh` for batch processing to a Python-native multiprocessing approach. `fragmentation_oop.py` now accepts a folder path as input, dynamically processes all enclosed CIF files using a configurable thread pool (`--nproc`), and automatically exports the `fragmentation_summary.csv` file natively. To ensure robustness and memory efficiency, both the CSV and the `fragments_collection.extxyz` collection are updated incrementally. Fragment names in the ExtXYZ header are automatically cleaned by replacing `[`/`]` with `_` and stripping the `.xyz` extension.

## Decision 2026-05-12: All extracted fragments must be a single connected component
- Context: The user noticed that some exported fragments occasionally contained disjoint sub-fragments (e.g., floating solvent molecules captured in the radius, or disconnected pieces remaining after extraction logic).
- Decision: Implemented a rigorous bond-graph Breadth-First-Search (BFS) filter named `enforce_single_molecule` in the parent `BaseFragmenter` class. Before writing any final `.xyz` file (in `MOFFragmenter` and `BioMolFragmenter`), the script systematically deletes all atoms that do not belong to the largest contiguous component.
- Consequences: This completely eliminates the possibility of dangling atoms or disjoint solvent, guaranteeing that every fragment is one single, fully connected molecule.

## Decision 2026-05-12: Aromatic H-caps are enforced globally planar via SVD fit for both Carbon and Nitrogen
- Context: The user noticed that capped hydrogens on aromatic rings (like phenyl linkers) were sometimes visibly out-of-plane. While the previous geometry enforcement calculated the local plane using the parent Carbon and its two immediate neighbors, thermal distortions in the CIF could misalign this local plane relative to the overall aromatic ring. Additionally, the collision guard prevented enforcement if UFF pulled the H into a tight spot, and the code explicitly ignored Nitrogen atoms (leaving imidazole rings broken).
- Decision: Upgraded `enforce_sp2_capped_h_geometry` to perform a Graph-BFS (up to 3 bonds) to find all heavy atoms in the aromatic ring. It computes the best-fit 3D plane using SVD (Singular Value Decomposition) on these atoms, and perfectly projects the local bisector vector onto this global average plane. Crucially, the target parent atom check was expanded to include Nitrogen (`"C"` and `"N"`), and the unreliable H-H collision guard was entirely removed.
- Consequences: All capped sp2 hydrogens (on phenyl, imidazole, pyridine, etc.) are now strictly and aggressively forced to mathematically align with the average 3D plane of their parent aromatic ring, guaranteeing rigorous planarity for QM calculations.

## Decision 2026-05-12: Bio-macromolecule fragments are strictly neutralized, minimally capped, and chemically deduplicated
- Context: The user requested fragmentation of single-chain biological macromolecules (PDBs). Furthermore, fragments must have exactly zero formal charge for QM calculations, and must not add heavy atoms (like ACE/NME) or invoke full-molecule optimization. The user also requested skipping duplicated fragments.
- Decision: Implemented `BioMolFragmenter` with a sequence-based sliding window (`--window-size`, `--stride`). Cut peptide bonds are capped strictly with a single H atom. To enforce zero charge, `_neutralize_window` builds a local bond graph to detect pH-charged functional groups (carboxylates, ammonium, guanidinium, imidazolium) and adds or removes H atoms appropriately. Added hydrogens are automatically geometry-optimized using `optimize_capped_h_geometry_only`; the rest of the molecule remains frozen. Finally, each extracted window is compared against a `seen` set using a chemical fingerprint (`composition` + `internal pairwise distances`). Chemically identical windows (e.g. from highly repetitive motifs) are skipped to prevent redundant `.xyz` files, though they are still documented in the CSV.
- Consequences: Bio fragments are chemically valid for QM without fabricating phantom H atoms on interior atoms. PDB parser reads sequence linearly, independent of radius or graph searches.
- Alternatives considered: Using spatial distance (radius) like MOFs; rejected because biological chains have well-defined sequence topologies that need continuous peptide segments. Using valence-based H addition; rejected because it over-capped interior fully bonded atoms.

## Decision 2026-05-11: COF helper libraries remove chemically duplicate building blocks globally
- Context: Visual QA found duplicate COF helper building blocks across different COF stems, including `COF-LZU8_00.xyz`/`COF-LZU8_01.xyz`. Exact coordinate comparison was too strict because chemically identical helper fragments can differ slightly in coordinates or orientation.
- Decision: COF helper export now treats duplicate nodes/linkers globally per folder using composition plus internal pair-distance fingerprints rounded to `0.1 A`. Before writing any new COF helper node or linker, the full target folder (`cof_nodes_lib/` or `cof_linkers_lib/`) is pruned with the same chemically aggressive fingerprint, then the new candidate is checked against all existing helpers.
- Consequences: COF helper libraries are intended to contain one representative per chemically equivalent building block, even across different COF names. This rule is global for COF helper export and is not filename-specific. A looser `0.1 A` fingerprint may merge near-identical conformers by design because the user wants chemically duplicate building blocks removed.
- Alternatives considered: Keep the safer `0.01 A` duplicate tolerance; rejected because it left visually/chemically identical building blocks such as COF-LZU8 duplicated.

## Decision 2026-05-08: Path B normal fragments retain terminal B/O node components
- Context: COF-10 normal dimer visually missed terminal chemistry because generic layered Path B stopped traversal when it reached the next B/O node and capped the linker-side atom instead.
- Decision: For Path B normal mode only, when growth reaches a neighboring non-core B/O node component, retain that terminal node component but do not continue growing beyond it. Minimized Path B remains unchanged.
- Consequences: COF-10 normal dimer grows from 104 atoms (`B2 C60 H38 O4`) to 112 atoms (`B8 C60 H28 O16`), split as two equal 56-atom layers. COF-10 minimized remains 66 atoms. Other generic 2-layer normal fragments also gain terminal B/O node chemistry; minimized counts stay unchanged.
- Alternatives considered: Increase radius or keep more disconnected layers; rejected because the missing part was at the node/linker boundary, not the layer count.

## Decision 2026-05-08: Make metallo-PC COF monomer and dimer outputs explicit
- Context: ZnPc-family COFs can represent two stacked layers, but visual checking now needs monomer and dimer fragments as separate outputs rather than one implicit dimer file.
- Decision: Add `--cof-layer auto|monomer|dimer` to COF extraction for layered COFs. `auto` remains the accepted dimer behavior for compatibility; metallo-PC Path J monomer skips shortest-axis duplication, while generic Path B monomer keeps one principal disconnected layer instead of two.
- Consequences: ZnPc-DPB and generic 2-layer Path B COFs can be generated as normal/minimized monomer files and normal/minimized dimer files. `run_cof_family.sh` also supports a `both` mode that writes `_monomer` and `_dimer` suffixed files.
- Alternatives considered: Always write both files from the core extractor; rejected because single-output CLI calls and existing test scripts should stay predictable.

## Decision 2026-05-08: Restore ZnPc/metallo-PC COF Path J before generic COF paths
- Context: ZnPc-DPB regenerated through generic COF Path B after the direct coffragmentor Path J machinery was absent from the active COF extraction code, producing a 92-atom minimized fragment where one dimer layer missed linker context.
- Decision: COF extraction again tries direct coffragmentor metallo-PC Path J before generic COF paths. It selects the Zn/N-rich phthalocyanine node, attaches coffragmentor linker images through B-O contacts, and duplicates the whole node+linker set along the shortest lattice vector so both dimer layers receive identical linker context.
- Consequences: ZnPc-DPB returns to Path J sizes: normal 322 atoms (`Zn2 B16 C192 H80 N16 O16`) and minimized 166 atoms (`Zn2 B4 C96 H32 N16 O16`), split as two equal 83-atom layers.
- Alternatives considered: Patch generic Path B minimization; rejected because accepted ZnPc behavior depends on direct coffragmentor node+linker assembly.

## Decision 2026-05-08: ZIF/minimum fragments preserve the first linker ring
- Context: ZIF minimum fragments exposed that the first linker ring may be a pentagon/heterocycle, while the previous helper looked for a six-carbon ring and could cut the actual first ring. The user also clarified that if the normal fragment is already below 80 atoms, no separate minimum is needed.
- Decision: Replace the six-carbon-only partial-linker ring detector with nearest heavy-cycle detection (cycle size 3-8) so the whole first ring is retained. Before any minimized MOF extraction, compute the normal fragment; if it has fewer than 80 atoms, write the normal fragment as the minimized output and skip further minimization.
- Consequences: ZIF min outputs whose normals are under 80 atoms now equal their normal outputs: ZIF-1/ZIF-10/ZIF-2 33 atoms, ZIF-11/ZIF-12 56 atoms, ZIF-20 49 atoms. Larger normal fragments still continue into the minimum Path J/legacy logic with first-ring preservation.
- Alternatives considered: Only special-case ZIF filenames; rejected because first-ring preservation and the <80 rule are general MOF minimum-fragment rules.

## Decision 2026-05-08: DUT-49 normal fragments keep one Cu2 paddlewheel node
- Context: DUT-49 is a Cu paddlewheel MOF, but visually it should behave like Cu-BTC rather than PCN/NU: one paddlewheel node is enough for the normal fragment. The stale saved normal output had `Cu3` and 547 atoms.
- Decision: Add `DUT-49*` to the one-node paddlewheel family used by legacy fallback logic. Keep Path J available for DUT-49; current Path J normal already produces the desired one-node `Cu2` fragment.
- Consequences: `test_on_other_common_mofs/DUT-49_frag.xyz` is regenerated to 370 atoms (`C192 H136 N8 O32 Cu2`). DUT-49 minimized remains 136 atoms (`C69 H49 N2 O14 Cu2`). Cu-BTC remains one-node, while PCN/NU continue using two-node atom-count selection.
- Alternatives considered: Route DUT-49 through PCN/NU two-node candidate selection; rejected because the user visually confirmed it should be one-node like Cu-BTC.

## Decision 2026-05-07: PCN/NU normal fragments keep two metal nodes
- Context: The earlier Cu2 paddlewheel exception fixed Cu-BTC normal fragments but also would suppress the second metal node in PCN/NU families where the linker topology requires two metal nodes for the normal fragment. Some PCN/NU cases use legacy Path C; others use moffragmentor Path J.
- Decision: Restrict the legacy Cu2 one-node suppression to one-node paddlewheel families such as `Cu-BTC*` and `DUT-49*`. For normal-mode `PCN*` and `NU*` structures, skip MOF Path J and use the legacy Path C candidate search, which tests possible second SBUs and selects the smallest metal-complete fragment. Minimized mode still uses Path J when available.
- Consequences: Verified normal outputs now include two metal nodes for PCN/NU (`PCN-61`: `Cu4`, `PCN-68`: `Cu4`, `NU-100SP`: `Cu4`, `NU-108-Cu`: `Cu4`, existing `PCN-60`/`NU-108-Zn`: `Zn4`) while Cu-BTC normal remains 90 atoms with `Cu2`; DUT-49 is intentionally excluded from PCN/NU two-node behavior. NU-108-Cu now selects the 516-atom candidate instead of the wrong 912-atom side.
- Alternatives considered: Keep the broad Cu2 suppression for all paddlewheel-like MOFs; rejected because PCN/NU need the opposite normal topology.

## Decision 2026-05-07: Keep one Cu paddlewheel node for Cu-BTC normal fragments
- Context: Cu-BTC falls back to the legacy MOF path, not Path J. Normal mode previously triggered Path C for any discrete two-metal SBU and added a second Cu2 paddlewheel node, producing `Cu4`; the user visually confirmed one Cu2 metal node is enough.
- Decision: Treat a discrete `Cu2` SBU as a complete single Cu paddlewheel node for one-node paddlewheel families such as `Cu-BTC*` and `DUT-49*`, and skip Path C's second-SBU expansion there. Minimized mode already used Path A and remains unchanged.
- Consequences: `test_on_cubtc/Cu-BTC.cif` normal changes from 156 atoms (`Cu4`) to 90 atoms (`C36 Cu2 H28 O24`); minimized remains 66 atoms (`C30 Cu2 H22 O12`).
- Alternatives considered: Keep generic Path C for all two-metal SBUs; rejected for Cu-BTC because it duplicates the metal node.

## Decision 2026-05-07: Recover missing MOF Path J branches from CIF open connectors
- Context: IRMOF-11 helper node+linker assembly found only five metal-attached linker images, leaving one node carboxyl carbon with no linker branch in both normal and minimized fragments. The missing neighbor exists in the original CIF but is not exposed as a metal-attached moffragmentor linker image.
- Decision: After helper assembly, MOF Path J detects open node carboxyl carbons and recovers missing organic branches directly from the original CIF image. Normal mode appends the recovered full linker branch; minimized mode appends the recovered first connected ring branch with C-H capping. Existing capped-H-only cleanup remains the final geometry step.
- Consequences: IRMOF-11 normal grows from 188 atoms to 221 atoms (`Zn4 C108 H84 O25`); minimized grows from 100 atoms to 111 atoms (`Zn4 C53 H39 O15`). IRMOF-1 remains 113/93 atoms.
- Alternatives considered: Increase periodic image search radius; rejected because IRMOF-11 still found only five attached helper linkers out to image span 3.

## Decision 2026-05-07: MOF Path J minimum keeps first rings on non-primary linkers
- Context: The minimum node+linker fragment should keep one full linker, but the other node attachment directions need chemical context from the first connected linker ring.
- Decision: In minimized MOF Path J, keep the highest-scoring attached linker image fully. For every other chemically attached linker image, keep the first six-membered carbon ring reached from the node-bound linker atoms, plus bonded hydrogens; node-side connector atoms already present in the selected node remain fixed. If trimming leaves retained ring/connector carbons unsaturated, add capped H atoms and pass them through capped-H-only cleanup.
- Consequences: IRMOF-1 minimized Path J grows from 38 atoms to 93 atoms (`Zn4 C43 H31 O15`) while normal remains 113 atoms (`Zn4 C48 H36 O25`). A simple C valence scan finds no under-coordinated carbons in the minimized IRMOF-1 fragment. Capped-H-only cleanup still applies at the end.
- Alternatives considered: Keep only one full linker and drop all other linker branches; rejected because it loses first-ring context around the other node attachment directions.

## Decision 2026-05-07: Use compact unique MOF helper-fragment exports
- Context: Helper node/linker filenames with composition and smiles were too long for visual inspection, and duplicate node/linker fragments cluttered the libraries.
- Decision: `mof_nodes_lib/` and `mof_linkers_lib/` exports now use compact per-folder names, e.g. `IRMOF-1_00.xyz`, `IRMOF-1_01.xyz`, and only write unique molecules based on composition plus internal pair-distance fingerprints. If files for the same structure stem already exist in a helper folder, that folder export is skipped to preserve prior visual-check files. For new stems, the exporter scans all existing `.xyz` files in the target helper folder and skips any molecule whose fingerprint already exists, even if it came from a different MOF.
- Consequences: Visual inspection folders are shorter and contain only unique node/linker geometries across the whole node or linker library. Node and linker files may share the same basename because they live in separate folders. Existing same-stem helper files are not overwritten.
- Alternatives considered: Keep verbose chemical metadata in filenames; rejected because the folder itself separates node/linker and visual checking benefits from compact names.

## Decision 2026-05-07: Final H cleanup only moves capped hydrogens
- Context: The user clarified that capped-H geometry cleanup must apply globally, but full-molecule RDKit/UFF optimization should not move framework/node/linker heavy atoms or pre-existing hydrogens.
- Decision: All MOF and COF extraction finalization paths now route capped hydrogens through `optimize_capped_h_geometry_only(...)`. This local helper adjusts only H atoms recorded as caps, using deterministic C-H sp2 and O-H angle cleanup; `refine_h_geometry_with_rdkit(...)` remains defined but is not called by extraction finalization.
- Consequences: Heavy atoms and non-capped hydrogens remain fixed after fragment assembly. IRMOF-1 Path J remains 113/38 atoms for normal/minimized fragments, with capped C-O-H angles at 109.5 degrees and no RDKit/UFF warnings.
- Alternatives considered: Full-molecule RDKit/UFF minimization with fixed heavy atoms; rejected because metal-containing fragments can lack UFF parameters and the requested invariant is capped-H-only movement.

## Decision 2026-05-04: Keep tied Path B layered neighbors
- Context: COF-202 visually missed one linker when Path B selected only one of two equally close B/O node components at 2.85 A.
- Decision: Path B keeps all same-spacing layered neighbor components within 0.05 A of the best spacing.
- Consequences: COF-202 normal fragment is 169 atoms (`Si2 B3 H75 C83 O6`); minimized fragment is 70 atoms (`Si2 B3 H33 C26 O6`). COF-102 and COF-300 smoke checks remain on Path A/Path C respectively.
- Alternatives considered: Keep only the first nearest layered component; rejected because tie ordering dropped a symmetry-equivalent linker.

## Decision 2026-05-04: Prioritize metallo-PC cores before B/O layered nodes
- Context: ZnPc-COF contains B/O linker nodes, so generic Path B ran before the porphyrin-like N-rich SBU detector and selected a layered B/O fragment.
- Decision: Add Zn bonding radius and detect N-rich metal macrocycles before generic B/O Path B when Zn/Cu/Fe/Co/Ni/Mn are present.
- Consequences: ZnPc-COF uses Path D (Metallo-PC core dimer). Normal keeps the two stacked ZnPc layers plus all connected benzene-1,4-diboronic acid linkers and cuts only far-side O-C bonds so far B-connected O atoms are H-capped: 242 atoms (`Zn2 B16 H64 C112 N16 O32`). Minimized keeps the two stacked ZnPc layers plus one full BDBA linker per layer, retains the ZnPc-core fused benzene perimeter, caps discarded linker attachment points with H, and uses the same far-side O-H termination: 146 atoms (`Zn2 B4 H40 C76 N16 O8`). COF-366 remains Path D porphyrin; COF-202 remains Path B layered set.
- Alternatives considered: Leave ZnPc as Path B; rejected because the intended SBU is the metallo-phthalocyanine core.

## Decision 2026-05-06: Directly combine coffragmentor node and linker for metallo-PC minimum fragments
- Context: Helper-guided selection inside UniFrag still allowed topology choices to truncate or misrepresent the linker. The user asked to undo that approach and chemically combine the node and linker returned by `coffragmentor.py`.
- Decision: For metallo-PC COFs, try a direct Path J that selects a Zn/N-rich node molecule and linker molecule(s) from `coffragmentor.py`, combines their coordinates directly, and duplicates the pair/set along the shortest lattice vector for the ZnPc two-layer dimer. Minimized mode keeps the nearest attached linker image; normal mode keeps all coffragmentor linker images in neighboring cells that chemically attach to the selected node, so all four sides of a ZnPc node are represented. If no metallo-PC node/linker is found, fall back to the existing UniFrag COF paths.
- Consequences: ZnPc-DPB normal/minimized are generated directly from coffragmentor node+linker images as 322/166 atoms. ZnPc-COF normal/minimized are generated as 210/138 atoms. Existing COF-202/300/366 minimized regression sizes remain 70/74/107 atoms and skip the direct path unless metallo-PC metals are present.
- Alternatives considered: Use coffragmentor only as an index/formula helper inside UniFrag linker selection; reverted because it still did not satisfy visual linker completeness.

## Decision 2026-05-06: Use moffragmentor node+linker combine as first MOF Path J
- Context: The project needs a general node+linker approximation for MOFs analogous to the COF `coffragmentor.py` path, with helper-exported node/linker libraries for visual inspection.
- Decision: MOF extraction now first tries `moffragmentor.MOF.from_cif(...).fragment()`. When nodes and linkers are found, Path J exports all helper nodes to `mof_nodes_lib/` and linkers to `mof_linkers_lib/`, selects the node nearest the unit-cell center, combines all chemically attached linker images for normal mode or the nearest attached linker image for minimized mode, and merges overlapping same-element boundary atoms. If moffragmentor is unavailable or no node/linker result is usable, the old MOF paths run unchanged.
- Consequences: `test_on_irmof1/IRMOF-1.cif` normal Path J gives 113 atoms (`Zn4 C48 H36 O25`); minimized gives 38 atoms (`Zn4 C13 H6 O15`). Visual review can inspect generated helper fragments in `mof_nodes_lib/` and `mof_linkers_lib/`.
- Alternatives considered: Keep MOF Path J separate from the main extraction path; rejected for this branch because the goal is to test node+linker as the first approximation.

## Decision 2026-05-08: Export COF helper node/linker libraries for visual QA
- Context: The COF workflow now needs the same eye-check helper libraries used for MOFs, so each COF run should save node/linker fragments for manual inspection.
- Decision: Add coffragmentor-based helper export to COF extraction. On each COF run, attempt `coffragmentor` fragmentation and export helper nodes/linkers to `cof_nodes_lib/` and `cof_linkers_lib/` with compact names (`<stem>_00.xyz`) and duplicate filtering based on composition + internal distance fingerprint. If same-stem files already exist in a given folder, skip writing in that folder.
- Consequences: COF helper libraries are now generated/maintained automatically during COF extraction without changing primary fragment output logic.
- Alternatives considered: Manual one-off helper exports; rejected because recurring visual QA needs automated, consistent exports.

## Decision 2026-05-09: Node+linker-first strategy is mandatory for both MOF and COF
- Context: Project direction requires helper fragmentation output (node + linker) to be the fundamental first attempt for every framework type.
- Decision: MOF extraction always attempts moffragmentor Path J first (no family skip gate). COF extraction now attempts a generic coffragmentor Path J node+linker combine first, then falls back to existing COF-specific paths when helper assembly is unavailable.
- Consequences: Path ordering is now unified across MOF/COF. Legacy family heuristics remain as fallback safety paths.
- Alternatives considered: Keep family-specific bypasses before Path J; rejected because the project’s core strategy is node+linker-first.

## Decision 2026-05-09: COF path selection prioritizes node/linker topology similarity over structure name
- Context: COF-6 and COF-66 regressions showed filename-specific routing is brittle and can misroute layered dimer assembly.
- Decision: Prefer topology-driven routing: detect node/linker character from cut-component signatures (node similarity first, linker similarity second). Use boroxine-like B-rich aromatic node signatures for strict boroxine node+linker handling; otherwise use generic graph/coffragmentor node+linker paths.
- Consequences: COF branching now starts from node/linker similarity rather than structure stem naming. Layered dimer construction should preserve crystal-based placement for selected node/linker components.
- Alternatives considered: Keep explicit per-structure filename gates; rejected as non-general and error-prone.

## Decision 2026-05-11: Path J COFs must produce dimers when face-to-face layers exist
- Context: COF-LZU1 and COF-LZU8 visually had correct monomer fragments, but their dimer outputs were identical to monomer outputs because coffragmentor Path J returned before the later layer-aware COF fallback paths. The user reiterated the global rule: detect face-to-face layers before finishing, and produce dimer versions when layers exist.
- Decision: In COF coffragmentor Path J finalization, when output mode is not `monomer` and the shortest lattice vector is a plausible face-to-face stacking distance (2.5-5.0 A), duplicate the completed capped monomer fragment by that stacking vector. This keeps node+linker-first Path J behavior while enforcing layered dimer generation globally for stacked 2D COFs.
- Consequences: COF-LZU1 dimer/min-dimer now double monomer/min atom counts and use 3.729 A layer spacing. COF-LZU8 dimer/min-dimer now double monomer/min atom counts and use 4.093 A layer spacing.
- Alternatives considered: Let only later Path A/B layer detection handle dimers; rejected because Path J can return early for valid node+linker COFs and must still honor the dimer rule.
