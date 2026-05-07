# Project Memory Template (Reusable)

Use this file as a persistent engineering memory for this project. Keep entries concise, date-stamped, and actionable.

> Coordination note: agents should also read `project-decisions.md` for durable decisions and `project-agent-log.md` for chronological handoffs.

## 0) Metadata

- Project name:
- Repository URL:
- Default branch:
- Primary language(s):
- Last updated: 2026-05-07
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
- Python/Node/Compiler version:
- OS notes:
- GPU/CPU notes:

## 5) Commands (Copy/Paste)

### 5.1 Setup
```bash
# install deps
```

### 5.2 Run
```bash
# run app/script
```

### 5.3 Test
```bash
# run tests
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

### Decision YYYY-MM-DD: Title
- Context:
- Decision:
- Consequences:
- Alternatives considered:

### Decision 2026-05-07: PCN/NU normal fragments keep two metal nodes
- Context: PCN and NU families can have linker structures where the chemically meaningful normal fragment spans two metal nodes, opposite to Cu-BTC where one Cu2 paddlewheel node is sufficient.
- Decision: Keep the Cu2 one-node suppression only for `Cu-BTC*`. For PCN/NU normal fragments, skip Path J and use the legacy Path C candidate search to choose the smallest metal-complete two-node fragment. Minimized mode remains one-node/one-linker focused and can still use Path J.
- Consequences: PCN/NU normal Cu paddlewheel cases now report `Cu4` or `Zn4` as appropriate, and NU-108-Cu selects 516 atoms instead of 912; `test_on_cubtc/Cu-BTC.cif` normal remains `Cu2` and 90 atoms.
- Alternatives considered: Apply the Cu2 suppression to every Cu paddlewheel family; rejected because PCN/NU need the opposite normal-fragment topology.

### Decision 2026-05-07: Cu-BTC normal uses one Cu2 paddlewheel node
- Context: Cu-BTC uses the legacy MOF path. Path C duplicated the discrete Cu2 paddlewheel in normal mode, but one metal node is the desired normal fragment.
- Decision: For `Cu-BTC*` discrete `Cu2` SBUs, skip normal Path C second-SBU expansion and keep Path A with the single Cu paddlewheel node.
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

