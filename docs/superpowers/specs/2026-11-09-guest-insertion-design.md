# Guest Insertion for UniFrag Fragments - Design

Date: 2026-11-09
Status: approved in conversation, awaiting written-spec review

## 1. Purpose

Produce MLP training data for host-guest interactions. After fragment generation,
UniFrag inserts guest gas molecules (one gas, e.g. CO2, or a binary mixture, e.g.
CO2/CH4) into each fragment, writes many randomly loaded copies, and these are
labelled with ORCA through OpenMLP.

The insertion logic is ported from the user's deepMOF_dev scripts
(`data_gen/gen_geoms/insertionAds.py`, `runFluidAds.py`, `runMixFluidAds.py`), with
no yaff dependency (numpy + ase + yaml only).

Success: one YAML file describes a whole run; with guests switched on, a run
ends with one extxyz file per guest set, every frame closed-shell with `charge=`
written, host atoms untouched, no guest overlapping anything under the vdW rule,
and mixture counts exactly in the requested ratio.

## 2. Out of scope

- Bio (`kind: bio`) fragments. Guests on with `kind: bio` is a config error.
- Periodic insertion. Fragments are finite clusters; distances are plain, never
  minimum-image.
- Charged or open-shell guests. Only neutral closed-shell molecules are expected;
  the electron-parity check in section 8 enforces it.
- Changing anything in the empty-fragment outputs (`fragments_collection.extxyz`,
  quarantine, summary CSV). Without `--config`, or with guests off, a run is
  byte-identical to today.

## 3. Configuration

One YAML file per run. `fragmentation:` first, `guests:` last.

```yaml
fragmentation:
  input_path: runUniFrag/zn_cifs   # CIF/PDB file or folder
  kind: mof                        # mof | cof | bio
  radius: 4.0
  nmetals: 3
  center: -1
  nproc: 8
  timeout: 300                     # s per structure; 0 = off
  overwrite: false
  cof_layer: auto                  # cof only
  # bio only: window_size, stride, chain, ph, no_pdbfixer, output_dir
  # also accepted: output, minimize

guests:
  enabled: true                    # master switch; false = no guests, rest not read

  single:                          # pure-gas loading; one output file per gas
    enabled: true
    species: [CO2, CH4]

  mixture:                         # mixed loading at a fixed molar ratio
    enabled: true
    species:
      - {name: CO2, ratio: 5}
      - {name: CH4, ratio: 1}

  # shared by single and mixture
  frames: [normal, min]            # any of: normal, min, onlylinker, onlynode
  copies_per_frame: 5
  box_margin: -1.0                 # A added to each side of the fragment bounding box
  min_box_extent: 8.0              # A; shorter box sides are stretched to this, centred
  vdw_shift: 0.35                  # A subtracted from every vdW radius
  vdw_scale: 1.0                   # multiplies every (radius - shift)
  vdw_radii: {}                    # optional per-element override in A, e.g. {Zr: 2.2}
  full_after_failures: 1000        # consecutive failed tries for the next guest -> full
  max_units: null                  # optional cap: single = molecules, mixture = ratio units
  seed: 12345
```

### 3.1 Rules

- **Keys.** `fragmentation:` keys are the existing argparse flags with `-` written as
  `_`. A flag typed on the command line overrides the YAML value
  (`--config run.yaml --nproc 1`). Without `--config` nothing changes.
- **Unknown keys** anywhere in the file are an error, and nothing runs.
- **Guest names.** A species is either a name (ase `molecule()` built-in geometry) or
  `{name: X, file: path.xyz}` (the user's own geometry; `name` is then only a tag).
  - Names may not contain `_` or `-`, because they appear in labels and file names.
- **Mixture species.** Needs at least two species. Each `ratio` is a positive integer.
- **Single species.** `single.species` must not be empty when `single.enabled`.
- **Block switches.** `guests.enabled: true` with both `single` and `mixture` off is an
  error, not a silent no-op.
- **`frames`** is a non-empty subset of `[normal, min, onlylinker, onlynode]`.
- **Numeric values.**
  - `copies_per_frame` ≥ 1.
  - `full_after_failures` ≥ 1.
  - `vdw_scale` > 0.
  - `min_box_extent` ≥ 0.
  - `max_units` is null or ≥ 1.

## 4. Components

| Unit | Responsibility |
|---|---|
| `guest_insertion.py` (new, repo root next to `fragmentation_oop.py`) | Config validation for `guests:`, guest geometry loading, radius table, box, single-copy fill, set and frame orchestration, file and CSV writing, resume. Runs on its own: `python guest_insertion.py run.yaml`. |
| `fragmentation_oop.py` `main()` (small edit) | `--config` flag, YAML → argparse defaults, unknown-key check, call into `guest_insertion.run(...)` after the MOF/COF branch has finished writing the collection. |
| `runUniFrag/test_guest_insertion.py` (new) | Plain-assert tests (section 9.1). |
| `runUniFrag/check_guest_fragments.py` (new) | Run-directory checker; every check must pass, exits 1 on any failure (section 8). |
| `example/run.yaml` (new) | Fully commented example of every key. |

Pure functions inside `guest_insertion.py` (testable without files):

- **`load_guest(spec)`**
  - Returns `(symbols, positions)` centred on the centre of mass.
- **`radius_table(overrides, shift, scale)`**
  - Returns a per-Z array of effective radii.
- **`insertion_box(coords, margin, min_extent)`**
  - Returns `(lo, hi)`.
- **`random_rotation(rng)`**
  - Returns a 3x3 matrix. It is a port of `_random_rotation`.
- **`fill_copy(host_species, host_coords, unit, guests, params, rng)`**
  - Returns the placed guests and the number of whole units, or a skip reason.

## 5. Algorithm

For each selected collection frame, for each guest set (each single gas
separately, then the mixture), make `copies_per_frame` independent copies, each
starting from the empty fragment.

1. **Box.** Axis-aligned bounding box of the host atoms, each side moved outward by
   `box_margin` (negative shrinks).
   - Any side shorter than `min_box_extent` is set to `min_box_extent`, centred on the host's box centre.
   - Flat linkers and COF layers therefore get room above and below the plane.
   - 3D MOF fragments wider than the minimum on every axis are unaffected.
2. **Radii.** r(Z) = (R(Z) − `vdw_shift`) × `vdw_scale`. R(Z) is:
   - the user override, if set;
   - else `ase.data.vdw_radii`, if not NaN;
   - else `ase.data.vdw_alvarez.vdw_radii`.

   An element with no radius anywhere is an error naming the element. ase's default table is NaN for Ti, V, Cr, Mn, Fe, Co, Zr and most other transition metals. The original scripts compared `d < nan + r`, which is always False, so guests were never rejected near those metals. This fallback fixes that.
3. **Unit.**
   - **Single gas:** the unit is one molecule.
   - **Mixture:** the unit is the species list expanded by ratio, in list order. For example, 5:1 gives CO2 ×5 then CH4 ×1.
4. **Fill.** Walk the unit order repeatedly. For the next guest, up to `full_after_failures` times:
   - draw a random rotation;
   - draw a centre of mass uniformly in the box;
   - accept if, for every new guest atom i and every atom j already present (host and earlier guests), |r_i − r_j| ≥ r(Z_i) + r(Z_j).

   Only the new guest's atoms are tested. Distances are plain Euclidean.

   The box is **full** when the next guest fails `full_after_failures` times in a row. Filling also stops when `max_units` whole units are placed.
5. **Trim.** Drop guests back to the last complete unit. Zero complete units → the copy is skipped with status `skipped: no whole unit fits`. A single gas never trims.
6. **Randomness.** Each copy has its own generator: `numpy.random.default_rng(SeedSequence([seed, crc32(parent label), crc32(set tag), copy index]))`.
   - Output does not depend on `nproc`, frame order, or whether the run was resumed.
   - `crc32` (zlib), not Python `hash`, which is salted per process.

## 6. Output

### 6.1 Files

All files go in the run directory: the folder for a folder input, or the CIF's directory for a single-file input.

| Set | File |
|---|---|
| single gas G | `fragments_guests_single_G.extxyz` |
| mixture A:B = a:b | `fragments_guests_mixture_A-B_a-b.extxyz` |
| all sets | `guests_summary.csv`, `guests_config.json` |

### 6.2 Frames

Each frame is the parent frame's atoms, unchanged and in the same order, followed by the guest atoms.

```
Properties=species:S:1:pos:R:3 label=IRMOF1FragMofMin-CO2x5CH4x1-3 charge=0 parent=IRMOF1FragMofMin guests="CO2:10 CH4:2" n_host=58 capped_h="93 94 95" pbc="F F F"
```

- **`label`**: `{parent}-{tag}-{copy}`.
  - The tag is the gas name for a single gas, or each species name followed by `x` and its ratio for a mixture (`CO2x5CH4x1`), so two ratios in one run dir never share a label (amended after the final review; was `CO2CH4`).
  - Never contains `_`.
- **`charge`**: always written. Equal to the parent's charge tag (0), because guests are neutral.
- **`parent`**: the parent frame's label.
- **`guests`**: the count of each species in this copy.
- **`n_host`**: the number of host atoms. Atoms `0..n_host-1` are host.
- **Parent tags**: other info tags of the parent (e.g. `capped_h`) are copied unchanged. Host indices do not move, so they stay valid.

### 6.3 Summary and console

- **`guests_summary.csv`**: one row per copy, with columns `parent,set,copy,<one count column per species>,units,status`.
  - `status` is `written` or `skipped: <reason>`.
- **Console**: one total line per set, e.g. `single CO2: 8180 written, 20 skipped`.

### 6.4 Frame selection

- Only `fragments_collection.extxyz` is read. The quarantine file never is.
- Frame class comes from the label ending:
  - `FragMof` / `FragCof` → `normal`
  - `FragMofMin` / `FragCofMin` → `min`
  - `…OnlyLinker[n]` → `onlylinker`
  - `…OnlyNode[n]` → `onlynode`
  - `frames` accepts any of these four. The default `[normal, min]` gives no guests to OnlyLinker or OnlyNode frames.
- For a single-CIF input, only frames for which `_label_belongs_to(label, _clean_base_of(input))` is true are used.

## 7. Wiring, resume, parallel

### 7.1 Wiring

- **`--config`.** `fragmentation_oop.py` gains `--config PATH`. Parsing order:
  1. Pre-parse `--config`.
  2. Load the YAML and validate all keys.
  3. `parser.set_defaults(**fragmentation)`.
  4. Full parse. `input_path` becomes `nargs="?"`, and is an error if given neither on the command line nor in the YAML.
- **After the MOF or COF branch.** When the collection is complete and `guests.enabled` is true, `main()` calls `guest_insertion.run(guests_cfg, run_dir, input_path, nproc, overwrite)`.
- **Bio.** `kind: bio` with guests on is rejected during config validation, before fragmentation starts.

### 7.2 Resume

- **Default (`overwrite: false`).** For each set whose output file exists, parents already listed in `guests_summary.csv` for that set are skipped. New parents are appended.
- **Settings check.** `guests_config.json` records the effective guest settings. A mismatch with the current config and `overwrite: false` stops the run with a message naming the changed keys.
  - A mixture at a different ratio is a different file, so it is not a mismatch.
- **Overwrite (`overwrite: true`).** Guest outputs, CSV and JSON are rewritten from scratch.

### 7.3 Parallel

- **Workers.** A `multiprocessing` pool of `nproc` workers over parent frames.
- **Write order.** The main process writes results in collection order, so file order is the same for any `nproc`.
- **Timeout.** None. `full_after_failures` bounds each copy.
- **Bounding the lookup.** Each trial tests the new guest's few atoms against the atoms already present. The cost is O(n_present) per guest atom, vectorised with numpy.

## 8. Checker: `runUniFrag/check_guest_fragments.py <run_dir>`

MUST (exit 1 on any failure):

- **[G1]** Host atoms (species and coordinates) equal the parent frame's, to 1e-6 Å.
- **[G2]** No guest–host or guest–guest pair is closer than r(Z_i) + r(Z_j). The radii are rebuilt from `guests_config.json`.
- **[G3]** Each guest's internal distances match its template to 1e-4 Å, so the rotation is rigid.
- **[G4]** Mixture frames: each species count is `units × ratio`. Single frames contain only their own gas.
- **[G5]** `charge=` is present on every frame and equals the parent's. The guest atoms' Z sum is even.
  - The host passed its own parity check when it entered the collection, so an even guest sum keeps the frame closed-shell.
  - This check works the same for MOF and COF.
- **[G6]** No `_` in any label. `parent` exists in `fragments_collection.extxyz`.

## 9. Testing and validation

### 9.1 `runUniFrag/test_guest_insertion.py`

Plain asserts, runnable with `python`:

- **Rotation:**
  - The rotation matrix is orthonormal with determinant +1.
  - The rotated guest keeps all pair distances.
- **Overlap rule:** a guest placed at exactly r_i + r_j − 0.01 Å is rejected; at + 0.01 Å it is accepted.
- **Trimming:** a 5:1 fill that stops after 8 CO2 + 1 CH4 is trimmed to 5 + 1; a fill that never completes one unit is skipped.
- **Determinism:** identical output for `nproc` 1 and 2; resume equals an uninterrupted run.
- **Radii:** Cr and Zr get the Alvarez values; an override wins.
- **Box:** a planar host gets a box side of `min_box_extent`; a 3D host keeps its margin box.
- **Labels and tags:** no `_` in labels; `charge=` present on every frame; `n_host` and host atoms equal the parent's.
- **Config errors:** an unknown key, a non-integer ratio, a one-species mixture, both blocks off, `kind: bio` with guests, and `_` or `-` in a guest name.
- **`--config` precedence:** a command-line flag overrides the YAML; without `--config` the existing argparse defaults are unchanged.

### 9.2 Runs (ARF only; no local runs, not the full collection)

All validation runs on ARF. Each run goes in a new run dir under `/arf/scratch/otayfuroglu`; check that the dir does not exist before creating it.

1. **Test set.**
   - 20 MOF structures and 20 COF structures, picked with a fixed seed from the datasets used in the earlier ARF validation runs.
   - The picked file names are saved in the run dir, so the set can be rerun.
2. **Unit tests.** Run `test_guest_insertion.py` on ARF.
3. **Byte-identical without guests.** On the same 40 structures, run current `main` and the new code without `--config`. The empty-fragment outputs (`fragments_collection.extxyz`, quarantine, summary CSV) must be byte-identical.
4. **Guest run.** The same 40 structures with `--config`: single CO2 and CH4 plus a CO2:CH4 = 5:1 mixture, 5 copies per frame, under nohup.
   - Report written and skipped counts per set, for MOF and COF separately.
   - Report the checker result.

## 10. Documentation

- **`example/run.yaml`**: every key, commented.
- **README**: a usage section for `--config` and for running `guest_insertion.py` on its own.
- **`project-decisions.md`**: a new entry covering:
  - the guest-insertion design;
  - the vdW NaN fallback;
  - the `min_box_extent` rule;
  - what each excludes (bio, periodic insertion, charged guests).
- **`project-memory.md`**: commands and outputs.
- **`project-agent-log.md`**: the handoff entry.

## 11. Differences from the deepMOF_dev scripts

| Original | Here | Why |
|---|---|---|
| yaff `System`, `merge`, `cell.mic` | numpy arrays, plain distances | No yaff dependency. A fragment is a finite cluster, so minimum-image wrongly tested guests near one box face against host atoms on the opposite face. |
| `scale_vdw` scales only the elements of the first structure | `vdw_scale` applies to every element | A guest element missing from that structure was left unscaled. |
| NaN vdW radius for many transition metals | Alvarez fallback, optional override | NaN made the overlap test always pass near Cr, Fe, Co, Zr, … |
| `nads1`/`nads2`, `n_trial`, `ratio` via CLI | fill until full at an integer ratio, trimmed to whole units; all settings in YAML | The user's choice. The ratio is exact. |
| `range(200)` copies, labels `{name}_{idx}` | `copies_per_frame` (default 5), labels `{parent}-{tag}-{copy}` | Data size; "_" is not allowed in OpenMLP labels. |
| Unseeded `np.random` | per-copy `SeedSequence` | Reproducible whatever `nproc`, order or resume. |
