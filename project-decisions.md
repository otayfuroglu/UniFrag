# Project Decisions

Durable implementation and architecture decisions for UniFrag. This file is the source of truth for decisions; keep entries concise, dated, and actionable.

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
