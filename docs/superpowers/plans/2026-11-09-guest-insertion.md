# Guest Insertion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** After fragment generation, load guest gases (one gas, or a fixed-ratio mixture) into every selected fragment frame until full, writing reproducible, closed-shell extxyz files for MLP training. One YAML file drives the whole run.

**Architecture:**
- **`guest_insertion.py`** is a new module at the repo root, next to `fragmentation_oop.py`. It holds:
  - the pure geometry functions;
  - the guest-section config validation;
  - the run orchestration: sets, copies, seeds, pool, resume and writing;
  - a standalone CLI.
- **`fragmentation_oop.py`** gains `--config` (YAML → argparse defaults) and one call into `guest_insertion.run(...)` at the end of the MOF and COF branches.
- **`runUniFrag/check_guest_fragments.py`** checks a finished run directory.

**Tech Stack:** Python 3, numpy, ase (≥ 3.23: `ase.build.molecule`, `ase.data.vdw_radii`, `ase.data.vdw_alvarez`, `ase.data.atomic_masses`), PyYAML, `multiprocessing`, `zlib.crc32`. No yaff.

**Spec:** `docs/superpowers/specs/2026-11-09-guest-insertion-design.md`. Read it first; section numbers below refer to it.

## Global Constraints

- No `_` in any frame label. Guest labels are `{parent}-{tag}-{copy}`.
- Guest names may not contain `_` or `-`.
- `charge=` is written on every guest frame and equals the parent frame's charge.
- Host atoms are copied unchanged, in the same order, before the guest atoms. They are never moved or deleted.
- Guest atoms are written in insertion order: the unit order, repeated `units` times.
- Distances are plain Euclidean. No minimum image, no pbc.
- Without `--config`, or with guests off, every existing output is byte-identical to current `main` (d5ff4bb).
- Bio (`kind: bio`) with guests on is a config error.
- Unknown YAML keys anywhere are a config error raised before any work starts.
- Defaults (spec §3):
  - `frames: [normal, min]`
  - `copies_per_frame: 5`
  - `box_margin: -1.0`
  - `min_box_extent: 8.0`
  - `vdw_shift: 0.35`
  - `vdw_scale: 1.0`
  - `vdw_radii: {}`
  - `full_after_failures: 1000`
  - `max_units: null`
  - `seed: 12345`
- Output names:
  - `fragments_guests_single_{G}.extxyz`
  - `fragments_guests_mixture_{A}-{B}_{a}-{b}.extxyz` (more species extend the same pattern)
  - `guests_summary.csv`
  - `guests_config.json`
- **No end-to-end local runs.**
  - Validation fragmentation runs happen only on ARF (Task 8).
  - The unit tests use small synthetic hosts, do no fragmentation, and run in seconds. They are run locally in the TDD loop and again on ARF.
- **No commits** unless the user asks (repo rule). Where a task ends, record the files touched; do not `git commit`.
- AGENTS.md: read `project-memory.md`, `project-decisions.md` and recent `project-agent-log.md` before editing, and update them in Task 7.
- Python locally: `/Users/omert/miniconda3/bin/python`. On ARF: `/arf/home/otayfuroglu/miniconda3/bin/python`.

## Review Focus

1. **A collection frame whose label contains `-`** (e.g. `COF-TpAzoFragCof`). The label class must come from the ending, not from splitting on `-`. The guest label must keep the parent label whole. Test: `test_frame_class_with_dash_in_base` (Task 4).
2. **A host with zero atoms, or one too crowded for any guest.** Expected: every copy is `skipped: no whole unit fits`, nothing crashes, the CSV still has a row per copy. Test: `test_crowded_host_skips_all_copies` (Task 4).
3. **A run interrupted after writing frames but before their CSV rows.** On resume, frames whose parent is not in the CSV for that set are dropped before appending, so no copy appears twice. Test: `test_resume_drops_orphan_frames` (Task 4).
4. **A user geometry file given by a relative path.** It resolves against the YAML file's directory, not the working directory. Test: `test_guest_file_relative_to_yaml` (Task 3).
5. **A YAML value outside argparse `choices`** (e.g. `kind: mfo`). argparse does not check `choices` for defaults, so the wiring must, and must give a clear error. Test: `test_config_bad_choice_rejected` (Task 5).

---

### Task 1: Geometry primitives

**Files:**
- Create: `guest_insertion.py`
- Create: `runUniFrag/test_guest_insertion.py`

**Interfaces:**
- Produces:
  - `class GuestConfigError(ValueError)`
  - `@dataclass(frozen=True) class Guest: name: str; symbols: tuple[str, ...]; positions: np.ndarray`. Positions are in Å, centred on the centre of mass (ase `atomic_masses`).
  - `load_guest(spec: str | dict, base_dir: str) -> Guest`
    - `"CO2"` → `ase.build.molecule("CO2")`.
    - `{"name": "CO2", "file": "geoms/co2.xyz"}` → `ase.io.read(os.path.join(base_dir, file))`. An absolute path is used as is.
    - A name containing `_` or `-`, or an unknown ase name with no `file`, raises `GuestConfigError`.
  - `radius_table(overrides: dict[str, float], shift: float, scale: float) -> np.ndarray`
    - Length 119, indexed by Z.
    - Value: `(R - shift) * scale`, where R is the override, else `ase.data.vdw_radii[Z]` if not NaN, else `ase.data.vdw_alvarez.vdw_radii[Z]`. Remaining NaN stays NaN.
  - `radii_for(symbols, table) -> np.ndarray`
    - Raises `GuestConfigError("no vdW radius for element X")` on NaN.
  - `insertion_box(coords: np.ndarray, margin: float, min_extent: float) -> tuple[np.ndarray, np.ndarray]`
    - Returns `(lo, hi)`.
  - `random_rotation(rng: np.random.Generator) -> np.ndarray`
    - A 3×3 matrix. Port of `insertAds._random_rotation` with `rng.random(3)` replacing `np.random.uniform(size=3)`.
    - The caller applies it as `pos @ M.T`. This equals the original `einsum('ib,ab->ia', pos, M)`.
- Test runner: `runUniFrag/test_guest_insertion.py` runs every `test_*` function in the file, prints `PASS n / FAIL n`, exits 1 on any failure, and needs no pytest. It does `sys.path.insert(0, <repo root>)`.

- [ ] **Step 1: Write failing tests**

```python
def test_rotation_is_proper():
    rng = np.random.default_rng(0)
    for _ in range(200):
        M = gi.random_rotation(rng)
        assert np.allclose(M @ M.T, np.eye(3), atol=1e-12)
        assert abs(np.linalg.det(M) - 1.0) < 1e-12

def test_rotation_keeps_pair_distances():
    g = gi.load_guest("CH4", ".")
    d0 = pdist(g.positions)            # small local helper, all pair distances
    M = gi.random_rotation(np.random.default_rng(1))
    assert np.allclose(pdist(g.positions @ M.T), d0, atol=1e-12)

def test_load_guest_builtin_centred():
    g = gi.load_guest("CO2", ".")
    assert g.symbols == ("C", "O", "O")
    assert np.allclose(com(g), 0.0, atol=1e-10)   # com via ase atomic_masses

def test_load_guest_bad_names():
    for bad in ["CO_2", "CO-2", "NotAMolecule"]:
        assert raises(gi.GuestConfigError, gi.load_guest, bad, ".")

def test_radius_fallback_and_override():
    t = gi.radius_table({}, 0.35, 1.0)
    assert abs(t[24] - (2.45 - 0.35)) < 1e-9      # Cr: Alvarez
    assert abs(t[40] - (2.52 - 0.35)) < 1e-9      # Zr: Alvarez
    assert abs(t[6] - (1.70 - 0.35)) < 1e-9       # C: ase default kept
    t2 = gi.radius_table({"Zr": 2.2}, 0.35, 0.9)
    assert abs(t2[40] - (2.2 - 0.35) * 0.9) < 1e-9

def test_box_margin_and_planar_stretch():
    cube = np.array([[x, y, z] for x in (0, 10) for y in (0, 10) for z in (0, 10)], float)
    lo, hi = gi.insertion_box(cube, -1.0, 8.0)
    assert np.allclose(lo, 1.0) and np.allclose(hi, 9.0)
    flat = np.array([[0, 0, 0], [10, 0, 0], [0, 10, 0], [10, 10, 0]], float)
    lo, hi = gi.insertion_box(flat, -1.0, 8.0)
    assert np.allclose(hi - lo, [8.0, 8.0, 8.0])
    assert np.allclose((hi + lo) / 2, [5.0, 5.0, 0.0])
```

- [ ] **Step 2: Run, expect failure** (`ModuleNotFoundError: guest_insertion`)

Run: `/Users/omert/miniconda3/bin/python runUniFrag/test_guest_insertion.py`

- [ ] **Step 3: Implement the Task 1 interfaces in `guest_insertion.py`**

  `insertion_box`:
  - per axis, `lo = min - margin`, `hi = max + margin`;
  - an axis with `hi - lo < min_extent` becomes `centre ± min_extent/2`, where centre is the unshrunk bounding-box centre.

- [ ] **Step 4: Run, expect `FAIL 0`**

---

### Task 2: Single-copy fill and trim

**Files:**
- Modify: `guest_insertion.py`
- Test: `runUniFrag/test_guest_insertion.py`

**Interfaces:**
- Consumes: `Guest`, `radii_for`, `random_rotation` (Task 1).
- Produces:
  - `@dataclass class FillResult: symbols: list[str]; coords: np.ndarray; counts: dict[str, int]; units: int; status: str`
    - `status` is `"written"` or `"skipped: no whole unit fits"`.
    - `counts` has an entry for every guest name in the unit, including 0.
  - `fill_copy(host_symbols: list[str], host_coords: np.ndarray, unit: list[Guest], table: np.ndarray, box: tuple[np.ndarray, np.ndarray], full_after_failures: int, max_units: int | None, rng: np.random.Generator) -> FillResult`
    - `unit` is the expanded insertion order, e.g. `[CO2]*5 + [CH4]`.
  - `overlaps(new_xyz, new_r, xyz, r) -> bool`
    - True if any `|new_i - x_j| < new_r_i + r_j`.

Algorithm (spec §5, steps 3–5):

```
placed = []; k = 0
loop:
    g = unit[k % len(unit)]
    for attempt in range(full_after_failures):
        M = random_rotation(rng); c = lo + rng.random(3) * (hi - lo)
        xyz = g.positions @ M.T + c
        if not overlaps(xyz, radii(g), host+placed coords, host+placed radii): place, break
    else: stop                        # full
    k += 1
    if k % len(unit) == 0 and max_units and k // len(unit) == max_units: stop
units = k // len(unit); keep the first units*len(unit) guests
```

- [ ] **Step 1: Write failing tests**

```python
def test_overlap_rule_boundary():
    r = np.array([1.0]); xyz = np.zeros((1, 3))
    assert gi.overlaps(np.array([[1.99, 0, 0]]), r, xyz, r)
    assert not gi.overlaps(np.array([[2.01, 0, 0]]), r, xyz, r)

def test_fill_trims_to_whole_units():
    # Fake guests "A" (one Ar atom) and "B" (one Kr atom). The host and box are
    # arranged so exactly 8 A then 1 B fit, but not a second B: a row of 9
    # sites, last one too small for a second B. Expected trim: 5 A + 1 B.
    res = gi.fill_copy(host_sym, host_xyz, [A]*5 + [B], table, box, 2000, None, rng)
    assert res.units == 1 and res.counts == {"A": 5, "B": 1}
    assert len(res.symbols) == 6

def test_fill_skips_when_no_unit_fits():
    res = gi.fill_copy(host_sym, tight_xyz, [A]*5 + [B], table, tight_box, 200, None, rng)
    assert res.status == "skipped: no whole unit fits" and res.units == 0
    assert len(res.symbols) == 0

def test_fill_respects_max_units():
    res = gi.fill_copy([], np.zeros((0, 3)), [CO2], table, big_box, 1000, 3, rng)
    assert res.units == 3 and res.counts == {"CO2": 3}

def test_fill_result_has_no_contacts():
    res = gi.fill_copy(host_sym, host_xyz, [CO2]*5 + [CH4], table, box, 1000, None, rng)
    assert no_pair_closer_than_radius_sum(host + res)   # local helper, all guest-vs-all pairs
```

  The exact host/box for `test_fill_trims_to_whole_units` is the implementer's choice. Use single-atom fake guests (`Guest("A", ("Ar",), np.zeros((1,3)))`) and a 1-D row of cells walled by host atoms. Make the assertion deterministic: every acceptable placement gives the same counts.

- [ ] **Step 2: Run, expect failures** (`AttributeError: fill_copy`)
- [ ] **Step 3: Implement `overlaps`, `FillResult` and `fill_copy`.** Vectorise `overlaps` with numpy broadcasting.
- [ ] **Step 4: Run, expect `FAIL 0`**

---

### Task 3: Guest-section config

**Files:**
- Modify: `guest_insertion.py`
- Test: `runUniFrag/test_guest_insertion.py`

**Interfaces:**
- Consumes: `load_guest`, `GuestConfigError` (Task 1).
- Produces:
  - `@dataclass(frozen=True) class GuestSet: kind: str; guests: tuple[Guest, ...]; ratios: tuple[int, ...]; set_id: str; tag: str; output: str`
    - `kind` is `"single"` or `"mixture"`.
    - Single: `set_id="single:CO2"`, `tag="CO2"`, `output="fragments_guests_single_CO2.extxyz"`, `ratios=(1,)`.
    - Mixture 5:1: `set_id="mixture:CO2-CH4:5-1"`, `tag="CO2x5CH4x1"` (amended after the final review), `output="fragments_guests_mixture_CO2-CH4_5-1.extxyz"`.
    - `GuestSet.unit() -> list[Guest]` expands by ratio, in list order.
  - `@dataclass(frozen=True) class GuestConfig: sets: tuple[GuestSet, ...]; frames: frozenset[str]; copies_per_frame: int; box_margin: float; min_box_extent: float; vdw_shift: float; vdw_scale: float; vdw_radii: dict; full_after_failures: int; max_units: int | None; seed: int`
    - `GuestConfig.shared() -> dict` returns the shared settings: every field except `sets`, with `frames` as a sorted list.
    - `GuestConfig.set_specs() -> dict[str, list]` maps `set_id` to `[[name, ratio, file-or-null], ...]`.
    - Both go into `guests_config.json`.
  - `parse_guests(section: dict | None, base_dir: str) -> GuestConfig | None`
    - Returns `None` when the section is absent or `enabled` is false. Nothing else is read in that case.
    - Raises `GuestConfigError` for every rule in spec §3.1. The message names the key.
  - `load_run_yaml(path: str) -> dict`
    - `yaml.safe_load`.
    - The top level may only contain `fragmentation` and `guests`, otherwise `GuestConfigError`.
    - Returns `{}` for an empty file.

- [ ] **Step 1: Write failing tests**

```python
BASE = {"enabled": True,
        "single": {"enabled": True, "species": ["CO2", "CH4"]},
        "mixture": {"enabled": True, "species": [{"name": "CO2", "ratio": 5},
                                                  {"name": "CH4", "ratio": 1}]}}

def test_defaults_and_sets():
    c = gi.parse_guests(BASE, ".")
    assert [s.set_id for s in c.sets] == ["single:CO2", "single:CH4", "mixture:CO2-CH4:5-1"]
    assert c.sets[2].output == "fragments_guests_mixture_CO2-CH4_5-1.extxyz"
    assert [g.name for g in c.sets[2].unit()] == ["CO2"] * 5 + ["CH4"]
    assert (c.copies_per_frame, c.box_margin, c.min_box_extent, c.full_after_failures,
            c.max_units, c.seed, c.frames) == (5, -1.0, 8.0, 1000, None, 12345,
                                               frozenset({"normal", "min"}))

def test_disabled_or_absent_returns_none():
    assert gi.parse_guests(None, ".") is None
    assert gi.parse_guests({"enabled": False, "bogus": 1}, ".") is None

def test_config_errors():
    bad = [
        {**BASE, "vdw_scal": 1.0},                                             # unknown key
        {**BASE, "mixture": {"enabled": True, "species": [{"name": "CO2", "ratio": 2.5},
                                                          {"name": "CH4", "ratio": 1}]}},
        {**BASE, "mixture": {"enabled": True, "species": [{"name": "CO2", "ratio": 1}]}},
        {"enabled": True, "single": {"enabled": False}, "mixture": {"enabled": False}},
        {**BASE, "single": {"enabled": True, "species": []}},
        {**BASE, "frames": ["normal", "big"]},
        {**BASE, "copies_per_frame": 0},
        {**BASE, "max_units": 0},
    ]
    for b in bad:
        assert raises(gi.GuestConfigError, gi.parse_guests, b, ".")

def test_guest_file_relative_to_yaml(tmp):
    os.makedirs(f"{tmp}/geoms"); write_xyz(f"{tmp}/geoms/co2.xyz", CO2_ATOMS)
    c = gi.parse_guests({"enabled": True,
                         "single": {"enabled": True,
                                    "species": [{"name": "CO2", "file": "geoms/co2.xyz"}]}},
                        base_dir=tmp)
    assert c.sets[0].guests[0].symbols == ("C", "O", "O")

def test_run_yaml_top_level(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {kind: mof}\nguest: {}\n")
    assert raises(gi.GuestConfigError, gi.load_run_yaml, f"{tmp}/r.yaml")
```

- [ ] **Step 2: Run, expect failures**
- [ ] **Step 3: Implement the Task 3 interfaces.**
  - Type checks: `ratio` is an `int` (not `bool`) and ≥ 1; `copies_per_frame`, `full_after_failures`, `seed` are `int`; floats accept int.
  - The two block keys `single` and `mixture` are optional; each accepts only `enabled` and `species`.
- [ ] **Step 4: Run, expect `FAIL 0`**

---

### Task 4: Run orchestration, output, resume, standalone CLI

**Files:**
- Modify: `guest_insertion.py`
- Modify: `fragmentation_oop.py:384-400` (`_parse_extxyz`). Add `"comment": comment` to each frame dict. This is additive; nothing else reads it.
- Test: `runUniFrag/test_guest_insertion.py`

**Interfaces:**
- Consumes:
  - Tasks 1–3.
  - From `fragmentation_oop` (import inside the functions that need it, to keep `import guest_insertion` light): `_parse_extxyz`, `_label_belongs_to`, `_clean_base_of`.
- Produces:
  - `frame_class(label: str) -> str | None`
    - Returns `"normal" | "min" | "onlylinker" | "onlynode"` from the label ending (spec §6.4), otherwise `None`.
    - Use `re.search(r"Frag(Mof|Cof)(Min|OnlyLinker\d*|OnlyNode\d*)?$", label)`.
  - `parse_header_tags(comment: str) -> dict[str, str]`
    - `shlex.split` on the comment line, then `key=value` pairs with quotes removed.
  - `guest_header(parent_comment: str, label: str, charge: int, parent: str, counts: dict[str, int], n_host: int) -> str`
    - Builds `Properties=species:S:1:pos:R:3 label=… charge=… parent=… guests="CO2:10 CH4:2" n_host=…`, then every other parent tag in its original order (re-quoted when it holds a space), then `pbc="F F F"`.
    - The parent's `Properties`, `label`, `charge` and `pbc` are not copied twice.
  - `copy_seed(seed: int, parent: str, set_id: str, copy: int) -> np.random.SeedSequence`
    - Returns `SeedSequence([seed, zlib.crc32(parent.encode()), zlib.crc32(set_id.encode()), copy])`.
  - `run(cfg: GuestConfig, run_dir: str, input_path: str, nproc: int, overwrite: bool) -> dict[str, tuple[int, int]]`
    - Returns `set_id → (written, skipped)`.
    - Prints one line per set: `{kind} {tag}: {w} written, {s} skipped`.
  - `main(argv: list[str] | None = None) -> int`
    - Standalone: `python guest_insertion.py run.yaml [--nproc N] [--overwrite]`.
    - Reads `fragmentation.input_path`, `nproc`, `overwrite` and `kind` from the YAML. The command-line flags override the YAML.
    - `kind: bio` with guests on → error, exit 2.
    - `run_dir` is `input_path` if it is a directory, else its parent directory.

`run` behaviour:
1. **Read and select.**
   - Read `fragments_collection.extxyz` in `run_dir`.
   - Keep frames whose `frame_class` is in `cfg.frames`.
   - If `input_path` is a file, also keep only frames with `_label_belongs_to(label, _clean_base_of(input_path))`.
2. **Resume, or start over.**
   - **`overwrite`**: delete every `cfg.sets[*].output`, `guests_summary.csv` and `guests_config.json` in `run_dir` first.
   - **Otherwise**, if `guests_config.json` exists:
     - Its `shared` must equal `cfg.shared()`. Its `set_specs` must equal the current ones for each `set_id` present in both.
     - On any difference, raise `GuestConfigError`, naming the differing keys and saying to set `overwrite: true`.
     - Then, for each set, collect `done = {parent for CSV rows of that set}`.
     - Rewrite the set's output keeping only frames whose `parent` tag is in `done`. This drops orphans from an interrupted run.
3. **Write the config record.**
   - Write `guests_config.json` as `{"shared": …, "set_specs": {merged old + current}}`.
   - Create `guests_summary.csv` if missing. Header: `parent,set,copy,{one column per species name across all sets, sorted},units,status`.
4. **Work items.**
   - One item per (parent frame) not yet done for at least one set.
   - The worker computes every copy of every pending set for that parent, using `fill_copy` with `copy_seed(...)` and `insertion_box(host, cfg.box_margin, cfg.min_box_extent)`.
   - Pool: `multiprocessing.Pool(nproc).imap` (ordered) when `nproc > 1`, otherwise a plain loop. The config goes to workers as an argument.
5. **Write results in collection order.**
   - For each parent, write the frames of each set to that set's output, then the CSV rows for that parent.
   - Label: `f"{parent}-{set.tag}-{copy}"`.
   - Charge: the parent's `charge`.
   - Atoms: host atoms, then the `FillResult` atoms.
   - Skipped copies get a CSV row only.

- [ ] **Step 1: Write failing tests.** A helper `make_run_dir(tmp, frames)` writes a `fragments_collection.extxyz` with `fragmentation_oop._write_extxyz`. Use small synthetic hosts:
  - a 6-ring of C labelled `XFragMof`;
  - the same as `XFragMofMin`;
  - `XFragMofOnlyLinker`;
  - `COF-TpAzoFragCof`.

```python
def test_frame_class():
    assert gi.frame_class("XFragMof") == "normal"
    assert gi.frame_class("XFragMofMin") == "min"
    assert gi.frame_class("XFragMofOnlyLinker2") == "onlylinker"
    assert gi.frame_class("XFragCofOnlyNode") == "onlynode"
    assert gi.frame_class("XWindow3") is None

def test_frame_class_with_dash_in_base():
    assert gi.frame_class("COF-TpAzoFragCof") == "normal"

def test_run_outputs_and_labels(tmp):
    d = make_run_dir(tmp, STD_FRAMES)            # includes a capped_h tag and charge=0
    out = gi.run(gi.parse_guests(BASE, "."), d, d, nproc=1, overwrite=False)
    fr = parse(f"{d}/fragments_guests_single_CO2.extxyz")
    assert {f["label"] for f in fr} >= {"XFragMof-CO2-0", "XFragMofMin-CO2-4"}
    assert not any(f["label"].startswith("XFragMofOnlyLinker") for f in fr)
    for f in all_guest_frames(d):
        t = gi.parse_header_tags(f["comment"])
        assert "_" not in t["label"] and t["charge"] == "0"
        n = int(t["n_host"]); host = by_label(d, t["parent"])
        assert f["species"][:n] == host["species"]
        assert np.allclose(f["coords"][:n], host["coords"])
        assert t["capped_h"] == gi.parse_header_tags(host["comment"])["capped_h"]
    assert out["mixture:CO2-CH4:5-1"][0] + out["mixture:CO2-CH4:5-1"][1] == 3 * 5

def test_nproc_independent(tmp):
    a = run_copy(tmp, "a", nproc=1); b = run_copy(tmp, "b", nproc=2)
    assert read_bytes_all_outputs(a) == read_bytes_all_outputs(b)

def test_resume_equals_uninterrupted(tmp):
    full = run_copy(tmp, "full", frames=STD_FRAMES)
    part = run_copy(tmp, "part", frames=STD_FRAMES[:1])
    append_frames(part, STD_FRAMES[1:]); gi.run(CFG, part, part, 1, False)
    assert read_bytes_all_outputs(full) == read_bytes_all_outputs(part)

def test_resume_drops_orphan_frames(tmp):
    d = run_copy(tmp, "o", frames=STD_FRAMES)
    drop_csv_rows_for(d, parent="XFragMofMin")          # simulate a crash before CSV write
    gi.run(CFG, d, d, 1, False)
    labels = [f["label"] for f in all_guest_frames(d)]
    assert len(labels) == len(set(labels))

def test_resume_settings_mismatch(tmp):
    d = run_copy(tmp, "m", frames=STD_FRAMES)
    changed = gi.parse_guests({**BASE, "box_margin": 0.0}, ".")
    assert raises(gi.GuestConfigError, gi.run, changed, d, d, 1, False)
    gi.run(changed, d, d, 1, True)                         # overwrite is allowed

def test_new_ratio_is_new_file_not_mismatch(tmp):
    d = run_copy(tmp, "r", frames=STD_FRAMES)
    other = {**BASE, "mixture": {"enabled": True, "species": [{"name": "CO2", "ratio": 3},
                                                              {"name": "CH4", "ratio": 1}]}}
    gi.run(gi.parse_guests(other, "."), d, d, 1, False)
    assert os.path.exists(f"{d}/fragments_guests_mixture_CO2-CH4_3-1.extxyz")

def test_crowded_host_skips_all_copies(tmp):
    d = make_run_dir(tmp, [dense_host("YFragMof")])        # lattice at 1.5 A spacing filling the box
    out = gi.run(gi.parse_guests(BASE, "."), d, d, 1, False)
    assert all(w == 0 for w, s in out.values())
    assert csv_rows(d) == 3 * 5                             # every copy still has a row

def test_single_file_input_only_its_frames(tmp):
    d = make_run_dir(tmp, STD_FRAMES + [ring("ZFragMof")])
    gi.run(gi.parse_guests(BASE, "."), d, f"{d}/Z.cif", 1, False)
    assert {gi.parse_header_tags(f["comment"])["parent"] for f in all_guest_frames(d)} == {"ZFragMof"}
```

- [ ] **Step 2: Run, expect failures**
- [ ] **Step 3: Implement the Task 4 interfaces and add `comment` to `_parse_extxyz`**
- [ ] **Step 4: Run, expect `FAIL 0`**

---

### Task 5: `--config` wiring in `fragmentation_oop.py`

**Files:**
- Modify: `fragmentation_oop.py` `main()` (≈ L10828):
  - Add the `--config` argument.
  - Make `input_path` `nargs="?"`.
  - Replace `args = parser.parse_args()` with `args, guests_cfg = _parse_args_with_config(parser)`.
- Modify: `fragmentation_oop.py`, after the MOF branch's final `print` block (≈ L11073) and after the COF branch's quarantine print (≈ L11268). Add:

```python
if guests_cfg is not None:
    import guest_insertion
    guest_insertion.run(guests_cfg, out_dir, args.input_path, args.nproc, args.overwrite)
```

- Test: `runUniFrag/test_guest_insertion.py`

**Interfaces:**
- Consumes: `load_run_yaml`, `parse_guests`, `GuestConfigError`, `run` (Tasks 3–4).
- Produces: `_parse_args_with_config(parser: argparse.ArgumentParser, argv: list[str] | None = None) -> tuple[argparse.Namespace, GuestConfig | None]`.
  1. `parse_known_args` for `--config` only.
  2. Without it: return `(parser.parse_args(argv), None)`. This is the unchanged path.
  3. With it:
     - `load_run_yaml(path)`.
     - Check the `fragmentation` keys against `{a.dest for a in parser._actions} - {"help", "config"}`; unknown keys → `parser.error`.
     - Check each value against `action.choices` where set → `parser.error`.
     - `parser.set_defaults(**fragmentation)`, then `parser.parse_args(argv)`.
     - `input_path` still `None` → `parser.error("input_path must be given on the command line or in fragmentation:")`.
     - `parse_guests(data.get("guests"), base_dir=dirname(abspath(path)))`. A `GuestConfigError` → `parser.error(str(e))`.
     - Guests on with `args.kind == "bio"` → `parser.error`.
  - `parser.error` exits 2 before any fragmentation starts.

- [ ] **Step 1: Write failing tests.** These call `_parse_args_with_config` on a parser built by a new `_build_parser()`. Extract the existing argparse block of `main()` into it so tests can get the real parser; `main()` calls it.

```python
def test_config_absent_is_unchanged():
    a, g = fo._parse_args_with_config(fo._build_parser(), ["x.cif"])
    assert g is None and a.kind == "mof" and a.radius == 6.0 and a.timeout == 300.0

def test_config_defaults_and_cli_override(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, kind: cof, nproc: 8, cof_layer: dimer}\n")
    a, g = fo._parse_args_with_config(fo._build_parser(), ["--config", f"{tmp}/r.yaml", "--nproc", "1"])
    assert (a.input_path, a.kind, a.nproc, a.cof_layer, g) == ("d", "cof", 1, "dimer", None)

def test_config_unknown_fragmentation_key(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, radious: 4}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])

def test_config_bad_choice_rejected(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, kind: mfo}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])

def test_config_bio_with_guests_rejected(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, kind: bio}\n"
                           "guests: {enabled: true, single: {enabled: true, species: [CO2]}}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])

def test_config_missing_input_path(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {kind: mof}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])
```

- [ ] **Step 2: Run, expect failures**
- [ ] **Step 3: Implement `_build_parser`, `_parse_args_with_config` and the two `run` calls.**
  - `_build_parser` must define the arguments exactly as now: same flags, defaults, help and order, plus `--config`.
  - The no-config path must be byte-identical; this is checked on ARF in Task 8.
- [ ] **Step 4: Run, expect `FAIL 0`**

---

### Task 6: Run-directory checker

**Files:**
- Create: `runUniFrag/check_guest_fragments.py`. Same shape as `check_mof_fragments.py`: a docstring that lists the MUST checks, a JSON report, `MUST: all passed`, exit 0 on pass and 1 on failure.
- Test: `runUniFrag/test_guest_insertion.py`

**Interfaces:**
- Consumes: `radius_table`, `radii_for`, `load_guest`, `parse_header_tags` and `GuestSet.unit()` (rebuilt from `guests_config.json` `set_specs`), plus `fragmentation_oop._parse_extxyz`.
- Produces: `check(run_dir: str) -> dict[str, list]`, mapping check id → failing labels. `main()` prints the report and exits.
  - Checks G1–G6 exactly as spec §8. G3 compares each guest block, segmented by the insertion order, with its template's pair distances to 1e-4 Å.
  - A user `file:` guest is rebuilt from the path recorded in `set_specs`, resolved against `run_dir`. To make that work, `run` records it as an absolute path.

- [ ] **Step 1: Write failing tests**

```python
def test_checker_passes_clean_run(tmp):
    d = run_copy(tmp, "c", frames=STD_FRAMES)
    assert all(not v for v in chk.check(d).values())

def test_checker_catches_moved_host_atom(tmp):
    d = run_copy(tmp, "h", frames=STD_FRAMES)
    nudge_atom(d, "fragments_guests_single_CO2.extxyz", frame=0, atom=0, dx=0.01)
    assert chk.check(d)["G1"]

def test_checker_catches_contact(tmp):
    d = run_copy(tmp, "k", frames=STD_FRAMES)
    move_first_guest_onto_host(d, "fragments_guests_single_CH4.extxyz", frame=0)
    assert chk.check(d)["G2"]

def test_checker_catches_bad_ratio(tmp):
    d = run_copy(tmp, "q", frames=STD_FRAMES)
    drop_last_guest(d, "fragments_guests_mixture_CO2-CH4_5-1.extxyz", frame=0)
    assert chk.check(d)["G4"]
```

- [ ] **Step 2: Run, expect failures**
- [ ] **Step 3: Implement `check_guest_fragments.py`**
- [ ] **Step 4: Run, expect `FAIL 0`**

---

### Task 7: Example config and documentation

**Files:**
- Create: `example/run.yaml`. Every key from spec §3, each with a one-line comment. Guests are on, with single CO2/CH4 and a 5:1 mixture.
- Modify: `README.md`. Add a "Run from a config file / guest insertion" section with three commands:
  - `python fragmentation_oop.py --config example/run.yaml`
  - `python guest_insertion.py example/run.yaml`
  - `python runUniFrag/check_guest_fragments.py <run_dir>`

  Also add `pyyaml ase` to the pip line.
- Modify: `project-decisions.md`. Add a top entry, "Decision 2026-11-09: Guest Insertion for MLP Training Data". It covers:
  - context;
  - decisions: run YAML, fill until full at an integer ratio trimmed to whole units, Alvarez vdW fallback, `min_box_extent`, per-copy seeds, resume rules, OnlyLinker/OnlyNode opt-in;
  - what each excludes: bio, periodic insertion, charged guests;
  - consequences.
- Modify: `project-memory.md`: commands, outputs, the ARF validation dir (filled in Task 8).
- Modify: `project-agent-log.md`: the handoff entry (finished in Task 8).

- [ ] **Step 1: Write the files.**
- [ ] **Step 2: Parse the example config.** `python -c "import guest_insertion as g; c=g.parse_guests(g.load_run_yaml('example/run.yaml')['guests'], 'example'); print([s.set_id for s in c.sets])"` prints the three set ids.

---

### Task 8: ARF validation (20 MOF + 20 COF)

**Files:** none in the repo. Run dirs go under `/arf/scratch/otayfuroglu/unifrag_guest_v1/`. **Check that it does not exist before creating it.**

SSH: `ssh -o BatchMode=yes -o ConnectTimeout=30 arf '…'`, filtering out `**` lines. Long jobs go through Slurm (debug partition) or nohup, because SSH connections drop.

- [ ] **Step 1: Check the path is free, then lay out the run dirs.**
  - **Code.** Two copies: `code_old` from `git show d5ff4bb:<file>` (`fragmentation_oop.py`, `coffragmentor.py`, `mof_linkers_lib/`, `mof_nodes_lib/`, `cof_*_lib/`, `runUniFrag/`), and `code_new` from the working tree, plus `guest_insertion.py`.
  - **Structures.** Pick with `random.Random(20261109).sample(sorted(cifs), 20)`:
    - MOFs from `/arf/scratch/otayfuroglu/unifrag_zn1110_i10b/cifs/*.cif` (1097 CIFs);
    - COFs from `/arf/scratch/otayfuroglu/unifrag_884_v14_dbg/cifs/*.cif` (847 CIFs).
  - **CIF copies.** Copy only `*.cif` into each of `{mof,cof}_{old,new,guest}/`. Save the picks as `picked_mof.txt` and `picked_cof.txt`.
  - **Config.** `guest.yaml` for each kind: `example/run.yaml` with `input_path` pointing at the `_guest` dir, `nproc: 20`, `timeout: 600`, `overwrite: true`.
- [ ] **Step 2: Unit tests on ARF.** `cd code_new && $P runUniFrag/test_guest_insertion.py`. Expected: `FAIL 0`.
- [ ] **Step 3: Submit one Slurm job** (debug partition, 40 CPUs, 3 h):
  - `code_old` on `{mof,cof}_old`, without `--config`;
  - `code_new` on `{mof,cof}_new`, without `--config`;
  - `code_new --config guest_{mof,cof}.yaml`;
  - then `check_guest_fragments.py` on each `_guest` dir, plus `check_mof_fragments.py` and `check_cof_fragments.py` on the `_new` dirs.
- [ ] **Step 4: Byte-identical check.** For each kind, `cmp` the `fragments_collection.extxyz`, `fragments_quarantine.extxyz` and `fragmentation_summary.csv` of `_old` against `_new` and against `_guest`. Expected: no differences.
- [ ] **Step 5: Report.**
  - For each kind and each set: written and skipped counts from the console lines.
  - Skip reasons, with the frames that skipped every copy.
  - The checker JSON for each kind.
  - Wall time of the guest step against the fragmentation step.
- [ ] **Step 6: Finish the docs.** Fill the validation numbers into the `project-decisions.md` entry, `project-memory.md` and the `project-agent-log.md` handoff (files touched, validation, follow-up risks). Offer to commit; do not commit unasked.
