# Guest Insertion — Status and Open Issues (parked 2026-11-09)

Branch: `guest-insertion` (not merged into `main`). Base: d5ff4bb.
Spec: `docs/superpowers/specs/2026-11-09-guest-insertion-design.md` · Plan: `docs/superpowers/plans/2026-11-09-guest-insertion.md` · Decision: `project-decisions.md` "Decision 2026-11-09".

## 1. Status

The code works and is validated. It is **parked** because the size of the guest-loaded frames is not settled (section 3). Settle that before any ORCA labelling.

- `python fragmentation_oop.py --config run.yaml` fragments, then adds guests. Every key is documented in `example/run.yaml`.
- `python guest_insertion.py run.yaml` adds guests to a run that is already fragmented.
- `python runUniFrag/check_guest_fragments.py <run_dir>` runs the QA checks G1–G6.
- `python runUniFrag/test_guest_insertion.py` runs 59 unit tests (about 10 s).
- Without `--config`, the output is byte-identical to d5ff4bb.

ARF validation is in `/arf/scratch/otayfuroglu/unifrag_guest_v1`: 20 Zn MOFs and 20 COFs, picked with `random.Random(20261109)` (lists in `picked_{mof,cof}.txt`). Jobs 6935903, 6935970 and 6936032:
- unit tests 59/59;
- the empty-fragment outputs are byte-identical to d5ff4bb;
- the checker passed on 525 MOF and 570 COF guest frames;
- 0 copies were skipped.

Latest outputs: `mof_fix/`, `cof_fix/`. Analysis script: `boxstats.py`.

## 2. What the user decided during design

- **Purpose:** MLP training data. **No yaff dependency.**
- **Sets:**
  - `single` writes one file per gas, e.g. `fragments_guests_single_CO2.extxyz`.
  - `mixture` uses an integer molar ratio, e.g. CO2:CH4 = 5:1. It fills until the box is full while keeping the ratio, and an incomplete last unit is **trimmed to whole units**.
- **Frames:** both normal and Min get guests. 5 copies per frame by default.
- **Config:** one YAML holds all parameters. `guests:` is the last section and has an `enabled` switch.
- **Approved changes vs the old scripts:**
  - vdW radii fall back to Alvarez where ase has NaN;
  - `min_box_extent` was added.

## 3. Open issue: frame size / ORCA cost (DECISION PENDING)

Fill-until-full makes large frames.

| ARF test set | Median host atoms | Median guest atoms |
|---|---|---|
| MOF | 88 | 215 |
| COF | 141 | 347 |

Worst case: WEXWEEFragMof, 554 host atoms in a 31.6 × 29.5 × 43.1 Å box, gets about 1,900–2,000 guest atoms per copy with CO2 alone and **3,745** with the CO2/CH4 mixture.

**Cost:** DFT cost grows roughly as N³.
- A median MOF frame is about (303/88)³ ≈ **40×** the empty fragment.
- WEXWEE is about **460×** and will likely not finish or will run out of memory.

**Why a cap may be wanted even though fill-until-full is physically sensible:**
1. Cost, as above.
2. The insertion box is the fragment's bounding box, not a pore. Corners of the box are vacuum, and some correspond to where the rest of the crystal would be. So most guests in a full box touch other guests, not the framework. The data becomes dominated by guest–guest pairs, while the MLP mainly needs host–guest contacts.
3. Every frame is at saturation. The MLP never sees low loading (one or two molecules near a site), which GCMC or MD at low pressure needs.

**Precedent:** the user's latest `runMixFluidAds.py` used **fixed counts** (`load_mix_fixed_nads(nads1, nads2)`), not fill-until-full.

**Options:**
- **A.** Keep fill-until-full and accept the cost.
- **B.** Random loading per copy: each copy gets a random number of units between 1 and full. This gives low-to-high loading and smaller frames.
- **C.** Cap total atoms per frame (e.g. 300) as a safety net for ORCA.
- **Claude's recommendation: B + C.** Use A + C if the target is saturated adsorption only. The user still has to choose the option and the atom limit for their ORCA setup.

### 3a. What `max_units` means (already implemented, default `null` = no cap)

`max_units` is an upper limit per copy, not a target: if the box fills first, the copy has fewer guests.

| Set | One unit | `max_units: 5` gives |
|---|---|---|
| single CO2 | 1 molecule | ≤ 5 CO2 |
| mixture 5:1 | 5 CO2 + 1 CH4 | ≤ 25 CO2 + 5 CH4 |

So a single value cannot suit both set types. Possible change: separate `max_units_single` and `max_units_mixture`, or replace the setting with B/C above. This is separate from `copies_per_frame` (5 independent random frames per fragment).

## 4. Insertion box: `box_margin` (= ase `vacuum`) and `min_box_extent`

- **What the setting is:** `box = bounding box ± box_margin`. This equals `atoms.center(vacuum=v)` in the old scripts.
- **Default −1.0** comes from the latest `runMixFluidAds.py` (`atoms.center(vacuum=-1.0)`). Older script versions used 0.5 (falling back to 0.0 when nothing fit) and 2.0.
  - −1.0 shrinks the box 1 Å inside the outermost atoms on each side.
  - `vacuum=0` would be `box_margin: 0.0`, a box 2 Å larger per axis.
  - Choosing the value is still open.
- **No retry:** the old scripts retried with a bigger box when nothing fit. The new code skips and logs the copy instead (0 such copies on ARF).
- **`min_box_extent: 8.0`:** a box side shorter than 8 Å is stretched to 8 Å, centred. Without it, flat linkers and COF layers have a box of negative or zero thickness and get **no** guests (every placement would be in the layer's plane).

**The user asked whether the 8 Å minimum causes the cost.** ARF data, single-CO2 files:

| | Frames with ≥1 side stretched | Box volume added by stretching | Guest atoms in stretched frames | Largest frames |
|---|---|---|---|---|
| MOF (175) | 65 | 1.9 % | 10.8 % | unstretched (WEXWEE 1,989) |
| COF (190) | 170 | 51 % | 79.4 % | unstretched (761FragCof 828) |

MOF breakdown by number of stretched sides:
- none: median 113 host / 174 guest atoms;
- 1 side: 76 / 87;
- 2 sides: 52 / 42;
- 3 sides: 46 / 27.

COF breakdown by number of stretched sides:
- none: median 176 host / 484 guest atoms;
- 1 side: 146 / 180;
- 2 sides: 70 / 38.

**Conclusion:**
- The cost comes mainly from large 3D fragments filled completely.
- The 8 Å rule matters for COFs (half of the box volume), but it is what gives flat layers any guests at all. 6 Å is a possible middle value (about one CO2 layer per side).
- Changing `min_box_extent` or `box_margin` does not fix WEXWEE-type frames; only a cap or random loading does.

## 5. Changes made on Claude's own judgement (check)

- **Mixture label tag carries the ratio:** `…-CO2x5CH4x1-3`, instead of the approved `CO2CH4`. Without it, two ratios in one run directory gave identical labels, which collide as OpenMLP keys. Checker G6 now fails duplicate labels across files.
- Duplicate species in a set are rejected.
- The `guests_summary.csv` species columns are the union across runs; the header is rewritten when a gas is added.
- The standalone tool exits 2 when guests are off, for bio, or when `input_path` is missing.
- The error text for a missing `input_path` changed (still exit 2).
- `check()` returns a frame count.

**Fixed after the independent review:**
- a quoted `"false"` in the YAML switched `overwrite` on;
- a single-CIF run deleted or left stale other structures' guests;
- a settings clash was reported only after fragmentation;
- an empty host was filled with pure gas;
- a missing collection exited 0.

**Bug found in the old scripts:** ase `vdw_radii` is NaN for Ti, V, Cr, Mn, Fe, Co, Zr, …, so guests near those metals were never rejected. In addition, yaff `cell.mic` applied a periodic overlap test to a finite cluster.

## 6. Deferred minors

- Guest names only reject `_` and `-`. `/`, `:`, spaces and CSV column names would break paths, tags or the CSV; names should be restricted to `^[A-Za-z0-9]+$`.
- Malformed YAML, a missing YAML file (standalone) and ase read errors give tracebacks instead of clean messages.
- `shlex.split` on a parent comment containing `'` or `\` aborts the guest step.
- There is no `--no-overwrite` / `--no-minimize` CLI switch.
- `guests_config.json` stores absolute guest geometry paths, so moving a run dir breaks the checker and resume for `file:` guests.
- A kill in the middle of a CSV write can leave a partial line.
- `guests_summary.csv` has no `n_atoms` column.

## 7. To resume

1. The user picks a size policy (section 3: A/B/C, atom limit, `max_units` split) and settles `box_margin` / `min_box_extent` (section 4).
2. Implement it test-first in `guest_insertion.py`, and add the frame-size check to `check_guest_fragments.py`.
3. Re-run the ARF set in a new run dir (check that it does not exist first) and compare with `unifrag_guest_v1`.
4. Merge `guest-insertion` into `main` only when the user asks.
