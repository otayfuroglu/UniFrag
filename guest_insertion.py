#!/usr/bin/env python3
"""Guest-molecule insertion into UniFrag fragments (MLP training data).

After fragmentation, each selected frame of fragments_collection.extxyz is
loaded with guest gas molecules - one gas, or a mixture at a fixed integer
molar ratio - until the box around the fragment is full, and written as
several independent copies. Ported from deepMOF_dev data_gen/gen_geoms
(insertionAds.py, runFluidAds.py, runMixFluidAds.py) without yaff.

Design: docs/superpowers/specs/2026-11-09-guest-insertion-design.md

Usage:  python guest_insertion.py run.yaml [--nproc N] [--overwrite]
"""
import csv
import json
import os
import re
import shlex
import zlib
from dataclasses import dataclass

import numpy as np


class GuestConfigError(ValueError):
    """A guest setting in the run YAML is invalid."""


@dataclass(frozen=True)
class Guest:
    """A rigid guest molecule; positions in A, centred on the centre of mass."""
    name: str
    symbols: tuple
    positions: np.ndarray
    file: str = None          # absolute path of a user geometry, None = ase built-in


def _check_name(name):
    if not isinstance(name, str) or not name:
        raise GuestConfigError(f"guest name must be a non-empty string, got {name!r}")
    # Names go into frame labels and file names: "_" breaks OpenMLP GBW keys,
    # "-" is the label field separator.
    if "_" in name or "-" in name:
        raise GuestConfigError(f"guest name {name!r} may not contain '_' or '-'")


def load_guest(spec, base_dir):
    """`"CO2"` -> ase built-in geometry; `{"name": .., "file": ..}` -> that file,
    a relative path resolved against `base_dir` (the YAML file's directory)."""
    from ase.data import atomic_masses

    if isinstance(spec, dict):
        unknown = set(spec) - {"name", "file", "ratio"}
        if unknown:
            raise GuestConfigError(f"unknown guest key(s) {sorted(unknown)} in {spec!r}")
        name, file = spec.get("name"), spec.get("file")
    else:
        name, file = spec, None
    _check_name(name)
    if file is not None:
        from ase.io import read
        path = file if os.path.isabs(file) else os.path.join(base_dir, file)
        path = os.path.abspath(path)
        if not os.path.exists(path):
            raise GuestConfigError(f"guest {name!r}: file not found: {path}")
        atoms = read(path)
    else:
        from ase.build import molecule
        try:
            atoms = molecule(name)
        except KeyError:
            raise GuestConfigError(
                f"guest {name!r} is not an ase built-in molecule; give `file:` "
                f"with its geometry") from None
        path = None
    masses = atomic_masses[atoms.numbers]
    pos = atoms.positions - (atoms.positions * masses[:, None]).sum(0) / masses.sum()
    return Guest(name, tuple(atoms.get_chemical_symbols()), pos, path)


def radius_table(overrides, shift, scale):
    """Effective radius per Z (index 0..118): (R - shift) * scale, R = user
    override, else ase.data.vdw_radii, else the Alvarez 2013 table where the
    default is NaN (most transition metals). Remaining NaN stays NaN."""
    from ase.data import vdw_radii, atomic_numbers
    from ase.data.vdw_alvarez import vdw_radii as alvarez

    base = np.full(119, np.nan)
    n = min(119, len(vdw_radii))
    base[:n] = vdw_radii[:n]
    alv = np.full(119, np.nan)
    m = min(119, len(alvarez))
    alv[:m] = alvarez[:m]
    base = np.where(np.isnan(base), alv, base)
    for sym, r in (overrides or {}).items():
        if sym not in atomic_numbers:
            raise GuestConfigError(f"vdw_radii: unknown element {sym!r}")
        base[atomic_numbers[sym]] = float(r)
    return (base - shift) * scale


def radii_for(symbols, table):
    from ase.data import atomic_numbers
    r = np.array([table[atomic_numbers[s]] for s in symbols], float)
    bad = [s for s, x in zip(symbols, r) if np.isnan(x)]
    if bad:
        raise GuestConfigError(f"no vdW radius for element {bad[0]}; set it in vdw_radii")
    return r


def insertion_box(coords, margin, min_extent):
    """Bounding box moved outward by `margin` per side (negative shrinks); a side
    shorter than `min_extent` becomes min_extent, centred on the bounding box."""
    coords = np.asarray(coords, float).reshape(-1, 3)
    if len(coords):
        mn, mx = coords.min(0), coords.max(0)
    else:
        mn = mx = np.zeros(3)
    lo, hi = mn - margin, mx + margin
    centre = (mn + mx) / 2
    short = (hi - lo) < min_extent
    lo = np.where(short, centre - min_extent / 2, lo)
    hi = np.where(short, centre + min_extent / 2, hi)
    return lo, hi


def random_rotation(rng):
    """Uniform random rotation (Arvo 1992), as insertAds._random_rotation;
    apply as pos @ M.T."""
    theta, phi, z = rng.random(3)
    theta = theta * 2.0 * np.pi
    phi = phi * 2.0 * np.pi
    z = z * 2.0
    r = np.sqrt(z)
    V = np.array([np.sin(phi) * r, np.cos(phi) * r, np.sqrt(2.0 - z)])
    st, ct = np.sin(theta), np.cos(theta)
    R = np.array(((ct, st, 0.0), (-st, ct, 0.0), (0.0, 0.0, 1.0)))
    return (np.outer(V, V) - np.eye(3)).dot(R)


def overlaps(new_xyz, new_r, xyz, r):
    """True when any new atom i and existing atom j are closer than r_i + r_j
    (plain distances: a fragment is a finite cluster, no periodic images)."""
    if len(xyz) == 0:
        return False
    d = np.linalg.norm(np.asarray(new_xyz)[:, None, :] - np.asarray(xyz)[None, :, :], axis=-1)
    return bool((d < np.asarray(new_r)[:, None] + np.asarray(r)[None, :]).any())


@dataclass
class FillResult:
    symbols: list
    coords: np.ndarray
    counts: dict
    units: int
    status: str


def fill_copy(host_symbols, host_coords, unit, table, box, full_after_failures,
              max_units, rng):
    """Fill one copy: insert guests in `unit` order, repeated, until the next
    one fails `full_after_failures` times in a row (box full) or `max_units`
    whole units are placed; then trim back to whole units."""
    lo, hi = (np.asarray(b, float) for b in box)
    guest_r = {id(g): radii_for(g.symbols, table) for g in unit}
    xyz = np.asarray(host_coords, float).reshape(-1, 3)
    rad = radii_for(host_symbols, table) if len(host_symbols) else np.zeros(0)
    placed = []                                      # (guest, coords) in order
    n_unit = len(unit)
    k = 0
    while True:
        g = unit[k % n_unit]
        gr = guest_r[id(g)]
        for _ in range(full_after_failures):
            M = random_rotation(rng)
            c = lo + rng.random(3) * (hi - lo)
            new = g.positions @ M.T + c
            if not overlaps(new, gr, xyz, rad):
                placed.append((g, new))
                xyz = np.vstack([xyz, new])
                rad = np.concatenate([rad, gr])
                break
        else:
            break                                    # full
        k += 1
        if k % n_unit == 0 and max_units and k // n_unit >= max_units:
            break
    units = k // n_unit
    kept = placed[:units * n_unit]
    counts = {g.name: 0 for g in unit}
    for g, _ in kept:
        counts[g.name] += 1
    symbols = [s for g, _ in kept for s in g.symbols]
    coords = np.vstack([c for _, c in kept]) if kept else np.zeros((0, 3))
    status = "written" if units else "skipped: no whole unit fits"
    return FillResult(symbols, coords, counts, units, status)


# ----------------------------------------------------------------- config
FRAME_CLASSES = ("normal", "min", "onlylinker", "onlynode")
DEFAULTS = {
    "frames": ["normal", "min"],
    "copies_per_frame": 5,
    "box_margin": -1.0,
    "min_box_extent": 8.0,
    "vdw_shift": 0.35,
    "vdw_scale": 1.0,
    "vdw_radii": {},
    "full_after_failures": 1000,
    "max_units": None,
    "seed": 12345,
}
_GUEST_KEYS = {"enabled", "single", "mixture"} | set(DEFAULTS)
_BLOCK_KEYS = {"enabled", "species"}


@dataclass(frozen=True)
class GuestSet:
    kind: str                 # "single" | "mixture"
    guests: tuple
    ratios: tuple
    set_id: str
    tag: str
    output: str

    def unit(self):
        """Insertion order of one ratio unit, e.g. CO2 x5 then CH4."""
        return [g for g, r in zip(self.guests, self.ratios) for _ in range(r)]


@dataclass(frozen=True)
class GuestConfig:
    sets: tuple
    frames: frozenset
    copies_per_frame: int
    box_margin: float
    min_box_extent: float
    vdw_shift: float
    vdw_scale: float
    vdw_radii: dict
    full_after_failures: int
    max_units: object
    seed: int

    def shared(self):
        """Settings shared by every set (recorded in guests_config.json)."""
        return {"frames": sorted(self.frames),
                "copies_per_frame": self.copies_per_frame,
                "box_margin": self.box_margin,
                "min_box_extent": self.min_box_extent,
                "vdw_shift": self.vdw_shift,
                "vdw_scale": self.vdw_scale,
                "vdw_radii": dict(sorted(self.vdw_radii.items())),
                "full_after_failures": self.full_after_failures,
                "max_units": self.max_units,
                "seed": self.seed}

    def set_specs(self):
        """set_id -> [[name, ratio, absolute file or None], ...]."""
        return {s.set_id: [[g.name, r, g.file] for g, r in zip(s.guests, s.ratios)]
                for s in self.sets}


def make_set(kind, guests, ratios):
    """The one place set ids, label tags and output file names are built."""
    names = [g.name for g in guests]
    if kind == "single":
        n = names[0]
        return GuestSet("single", tuple(guests), (1,), f"single:{n}", n,
                        f"fragments_guests_single_{n}.extxyz")
    nm, rt = "-".join(names), "-".join(map(str, ratios))
    return GuestSet("mixture", tuple(guests), tuple(ratios), f"mixture:{nm}:{rt}",
                    "".join(f"{n}x{r}" for n, r in zip(names, ratios)),
                    f"fragments_guests_mixture_{nm}_{rt}.extxyz")


def sets_from_specs(set_specs):
    """Rebuild GuestSets from guests_config.json's set_specs (used by the checker)."""
    out = []
    for sid, members in set_specs.items():
        guests = tuple(load_guest({"name": n, "file": f} if f else n, ".")
                       for n, _r, f in members)
        out.append(make_set(sid.split(":", 1)[0], guests, tuple(r for _n, r, _f in members)))
    return out


def _is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def _is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _block(section, key):
    blk = section.get(key)
    if blk is None:
        return None
    if not isinstance(blk, dict):
        raise GuestConfigError(f"guests.{key} must be a mapping")
    unknown = set(blk) - _BLOCK_KEYS
    if unknown:
        raise GuestConfigError(f"unknown key(s) {sorted(unknown)} in guests.{key}")
    en = blk.get("enabled", False)
    if not isinstance(en, bool):
        raise GuestConfigError(f"guests.{key}.enabled must be true or false")
    if not en:
        return None
    sp = blk.get("species")
    if not isinstance(sp, list) or not sp:
        raise GuestConfigError(f"guests.{key}.species must be a non-empty list")
    return sp


def parse_guests(section, base_dir):
    """Validate the `guests:` section. None when it is absent or switched off
    (nothing else is read then); GuestConfigError naming the key otherwise."""
    if section is None:
        return None
    if not isinstance(section, dict):
        raise GuestConfigError("guests must be a mapping")
    en = section.get("enabled", False)
    if not isinstance(en, bool):
        raise GuestConfigError("guests.enabled must be true or false")
    if not en:
        return None
    unknown = set(section) - _GUEST_KEYS
    if unknown:
        raise GuestConfigError(f"unknown key(s) {sorted(unknown)} in guests")
    v = {**DEFAULTS, **{k: section[k] for k in DEFAULTS if k in section}}

    sets = []
    single = _block(section, "single")
    if single is not None:
        seen = set()
        for spec in single:
            if isinstance(spec, dict) and "ratio" in spec:
                raise GuestConfigError("guests.single species take no ratio")
            g = load_guest(spec, base_dir)
            if g.name in seen:
                raise GuestConfigError(f"guests.single lists {g.name} twice")
            seen.add(g.name)
            sets.append(make_set("single", (g,), (1,)))
    mixture = _block(section, "mixture")
    if mixture is not None:
        if len(mixture) < 2:
            raise GuestConfigError("guests.mixture needs at least two species")
        guests, ratios = [], []
        for spec in mixture:
            if not isinstance(spec, dict) or "ratio" not in spec:
                raise GuestConfigError(f"guests.mixture species need name and ratio: {spec!r}")
            r = spec["ratio"]
            if not _is_int(r) or r < 1:
                raise GuestConfigError(f"guests.mixture ratio must be a positive integer, got {r!r}")
            guests.append(load_guest(spec, base_dir))
            ratios.append(r)
        names = [g.name for g in guests]
        if len(set(names)) != len(names):
            raise GuestConfigError(f"guests.mixture lists a species twice: {names}")
        sets.append(make_set("mixture", tuple(guests), tuple(ratios)))
    if not sets:
        raise GuestConfigError("guests.enabled is true but neither single nor mixture is enabled")

    fr = v["frames"]
    if not isinstance(fr, list) or not fr or any(f not in FRAME_CLASSES for f in fr):
        raise GuestConfigError(f"guests.frames must be a non-empty subset of {list(FRAME_CLASSES)}")
    for k in ("copies_per_frame", "full_after_failures"):
        if not _is_int(v[k]) or v[k] < 1:
            raise GuestConfigError(f"guests.{k} must be an integer >= 1")
    if not _is_int(v["seed"]):
        raise GuestConfigError("guests.seed must be an integer")
    if v["max_units"] is not None and (not _is_int(v["max_units"]) or v["max_units"] < 1):
        raise GuestConfigError("guests.max_units must be null or an integer >= 1")
    for k in ("box_margin", "min_box_extent", "vdw_shift", "vdw_scale"):
        if not _is_num(v[k]):
            raise GuestConfigError(f"guests.{k} must be a number")
    if v["vdw_scale"] <= 0:
        raise GuestConfigError("guests.vdw_scale must be > 0")
    if v["min_box_extent"] < 0:
        raise GuestConfigError("guests.min_box_extent must be >= 0")
    vr = v["vdw_radii"] or {}
    if not isinstance(vr, dict) or not all(_is_num(x) for x in vr.values()):
        raise GuestConfigError("guests.vdw_radii must map element symbols to numbers")
    radius_table(vr, 0.0, 1.0)              # unknown element symbols raise here

    return GuestConfig(tuple(sets), frozenset(fr), v["copies_per_frame"],
                       float(v["box_margin"]), float(v["min_box_extent"]),
                       float(v["vdw_shift"]), float(v["vdw_scale"]),
                       {k: float(x) for k, x in vr.items()},
                       v["full_after_failures"], v["max_units"], v["seed"])


def load_run_yaml(path):
    """Read a run YAML; its top level may hold only `fragmentation` and `guests`."""
    import yaml
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise GuestConfigError(f"{path}: the top level must be a mapping")
    unknown = set(data) - {"fragmentation", "guests"}
    if unknown:
        raise GuestConfigError(f"{path}: unknown top-level key(s) {sorted(unknown)}; "
                               f"only fragmentation and guests are allowed")
    return data


# -------------------------------------------------------------- run / output
SUMMARY_CSV = "guests_summary.csv"
CONFIG_JSON = "guests_config.json"
_SUFFIX_CLASS = re.compile(r"Frag(?:Mof|Cof)(Min|OnlyLinker\d*|OnlyNode\d*)?$")


def frame_class(label):
    """normal / min / onlylinker / onlynode from the label ending, else None."""
    m = _SUFFIX_CLASS.search(label)
    if not m:
        return None
    s = m.group(1) or ""
    if not s:
        return "normal"
    if s == "Min":
        return "min"
    return "onlylinker" if s.startswith("OnlyLinker") else "onlynode"


def parse_header_tags(comment):
    """key=value pairs of an extxyz comment line, quotes removed, in order."""
    tags = {}
    for tok in shlex.split(comment):
        if "=" in tok:
            k, v = tok.split("=", 1)
            tags[k] = v
    return tags


def _fmt_tag(k, v):
    return f'{k}="{v}"' if (v == "" or any(c.isspace() for c in v)) else f"{k}={v}"


def guest_header(parent_comment, label, charge, parent, counts, n_host):
    guests = " ".join(f"{n}:{c}" for n, c in counts.items())
    head = (f"Properties=species:S:1:pos:R:3 label={label} charge={int(charge)} "
            f'parent={parent} guests="{guests}" n_host={n_host}')
    rest = [_fmt_tag(k, v) for k, v in parse_header_tags(parent_comment).items()
            if k not in ("Properties", "label", "charge", "pbc")]
    return " ".join([head] + rest + ['pbc="F F F"'])


def copy_seed(seed, parent, set_id, copy):
    return np.random.SeedSequence(
        [seed, zlib.crc32(parent.encode()), zlib.crc32(set_id.encode()), copy])


def _work(item):
    """Every pending copy of one parent frame (runs in a pool worker)."""
    cfg, label, species, coords, set_ids = item
    table = radius_table(cfg.vdw_radii, cfg.vdw_shift, cfg.vdw_scale)
    host = np.asarray(coords, float).reshape(-1, 3)
    box = insertion_box(host, cfg.box_margin, cfg.min_box_extent)
    out = []
    for s in cfg.sets:
        if s.set_id not in set_ids:
            continue
        unit = s.unit()
        for k in range(cfg.copies_per_frame):
            if not len(host):
                out.append((s.set_id, k, FillResult([], np.zeros((0, 3)),
                                                    {g.name: 0 for g in unit}, 0,
                                                    "skipped: empty host")))
                continue
            rng = np.random.default_rng(copy_seed(cfg.seed, label, s.set_id, k))
            out.append((s.set_id, k, fill_copy(species, host, unit, table, box,
                                                cfg.full_after_failures, cfg.max_units, rng)))
    return out


def _write_frame(path, species, coords, header):
    with open(path, "a") as f:
        f.write(f"{len(species)}\n{header}\n")
        for sym, p in zip(species, coords):
            f.write(f"{sym:<8}{p[0]:16.8f}{p[1]:16.8f}{p[2]:16.8f}\n")


def _read_csv(path):
    if not os.path.exists(path):
        return [], []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        return list(r.fieldnames or []), list(r)


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, header, restval="0")
        w.writeheader()
        w.writerows(rows)


def run_dir_of(input_path):
    """Where a run's outputs live: the folder itself, or the input file's folder."""
    return input_path if os.path.isdir(input_path) else \
        os.path.dirname(os.path.abspath(input_path))


def check_resume(cfg, run_dir, input_path, overwrite):
    """Settings recorded by an earlier guest run in run_dir must match, unless
    this call starts the guest outputs again (folder input with overwrite).
    Returns the earlier set_specs to keep; raises GuestConfigError otherwise."""
    json_path = os.path.join(run_dir, CONFIG_JSON)
    if (overwrite and os.path.isdir(input_path)) or not os.path.exists(json_path):
        return {}
    with open(json_path) as f:
        old = json.load(f)
    shared, specs = cfg.shared(), cfg.set_specs()
    diff = sorted(k for k in set(old["shared"]) | set(shared)
                  if old["shared"].get(k) != shared.get(k))
    diff += sorted(sid for sid in specs
                   if sid in old["set_specs"] and old["set_specs"][sid] != specs[sid])
    if diff:
        raise GuestConfigError(
            f"guest settings differ from the earlier guest run in {run_dir} "
            f"({', '.join(diff)}). To redo only the guests, run "
            f"`python guest_insertion.py <run.yaml> --overwrite` on the folder; "
            f"with --config, overwrite: true also redoes the fragmentation.")
    return old["set_specs"]


def run(cfg, run_dir, input_path, nproc, overwrite):
    """Load guests into the selected frames of run_dir/fragments_collection.extxyz.
    Returns set_id -> (written, skipped) for this call."""
    from fragmentation_oop import _parse_extxyz, _label_belongs_to, _clean_base_of

    sets = {s.set_id: s for s in cfg.sets}
    csv_path = os.path.join(run_dir, SUMMARY_CSV)
    json_path = os.path.join(run_dir, CONFIG_JSON)
    outputs = {sid: os.path.join(run_dir, s.output) for sid, s in sets.items()}

    coll = os.path.join(run_dir, "fragments_collection.extxyz")
    if not os.path.exists(coll):
        raise GuestConfigError(f"no fragments_collection.extxyz in {run_dir}; "
                               f"nothing to load guests into")
    # A single-file input touches only its own structure, as fragmentation
    # does: its guests are always redone (its collection frames were just
    # replaced) and the other structures' guests are kept, even with overwrite.
    single = not os.path.isdir(input_path)
    base = _clean_base_of(input_path) if single else None
    frames = [f for f in _parse_extxyz(coll) if frame_class(f["label"]) in cfg.frames]
    if single:
        frames = [f for f in frames if _label_belongs_to(f["label"], base)]

    if overwrite and not single:
        for p in list(outputs.values()) + [csv_path, json_path]:
            if os.path.exists(p):
                os.remove(p)
    old_specs = check_resume(cfg, run_dir, input_path, overwrite)
    with open(json_path, "w") as f:
        json.dump({"shared": cfg.shared(), "set_specs": {**old_specs, **cfg.set_specs()}},
                  f, indent=1)

    # Resume: a parent is done for a set when its rows are in the CSV; frames of
    # parents without rows (an interrupted write) are dropped before appending.
    header, rows = _read_csv(csv_path)
    if single:
        rows = [r for r in rows if not _label_belongs_to(r["parent"], base)]
    species_cols = sorted({g.name for s in cfg.sets for g in s.guests}
                          | (set(header) - {"parent", "set", "copy", "units", "status"}))
    new_header = ["parent", "set", "copy"] + species_cols + ["units", "status"]
    if single or header != new_header:
        _write_csv(csv_path, new_header, rows)
    done = {sid: set() for sid in sets}
    for r in rows:
        if r["set"] in done:
            done[r["set"]].add(r["parent"])
    for sid, p in outputs.items():
        if os.path.exists(p):
            keep = [f for f in _parse_extxyz(p)
                    if parse_header_tags(f["comment"]).get("parent") in done[sid]]
            os.remove(p)
            for f in keep:
                _write_frame(p, f["species"], f["coords"], f["comment"])

    items, parents = [], []
    for f in frames:
        pending = tuple(sid for sid in sets if f["label"] not in done[sid])
        if pending:
            items.append((cfg, f["label"], f["species"], f["coords"], pending))
            parents.append(f)

    totals = {sid: [0, 0] for sid in sets}

    def flush(parent, results):
        n_host = len(parent["species"])
        new_rows = []
        for sid, k, res in results:
            s = sets[sid]
            row = {"parent": parent["label"], "set": sid, "copy": k,
                   "units": res.units, "status": res.status, **res.counts}
            new_rows.append(row)
            if res.units:
                totals[sid][0] += 1
                label = f"{parent['label']}-{s.tag}-{k}"
                header = guest_header(parent["comment"], label, parent["charge"],
                                      parent["label"], res.counts, n_host)
                _write_frame(outputs[sid], list(parent["species"]) + res.symbols,
                             np.vstack([np.asarray(parent["coords"], float), res.coords]),
                             header)
            else:
                totals[sid][1] += 1
        with open(csv_path, "a", newline="") as fh:
            csv.DictWriter(fh, new_header, restval="0").writerows(new_rows)

    if nproc > 1 and len(items) > 1:
        import multiprocessing
        with multiprocessing.Pool(nproc) as pool:
            for parent, results in zip(parents, pool.imap(_work, items)):
                flush(parent, results)
    else:
        for parent, item in zip(parents, items):
            flush(parent, _work(item))

    for sid, s in sets.items():
        w, k = totals[sid]
        print(f"{s.kind} {s.tag}: {w} written, {k} skipped")
    return {sid: tuple(v) for sid, v in totals.items()}


def main(argv=None):
    """python guest_insertion.py run.yaml [--nproc N] [--overwrite]"""
    import argparse
    ap = argparse.ArgumentParser(description="Insert guest molecules into an existing UniFrag run.")
    ap.add_argument("config", help="run YAML (fragmentation: input_path, guests: ...)")
    ap.add_argument("--nproc", type=int, default=None)
    ap.add_argument("--overwrite", action="store_true", default=None)
    a = ap.parse_args(argv)
    try:
        data = load_run_yaml(a.config)
        frag = data.get("fragmentation") or {}
        cfg = parse_guests(data.get("guests"), os.path.dirname(os.path.abspath(a.config)))
        if cfg is None:
            raise GuestConfigError(f"{a.config}: guests are not enabled")
        if frag.get("kind", "mof") == "bio":
            raise GuestConfigError("guest insertion is not supported for kind: bio")
        input_path = frag.get("input_path")
        if not input_path:
            raise GuestConfigError(f"{a.config}: fragmentation.input_path is missing")
        y_nproc, y_over = frag.get("nproc", 1), frag.get("overwrite", False)
        if not _is_int(y_nproc) or y_nproc < 1:
            raise GuestConfigError(f"fragmentation.nproc must be an integer >= 1, got {y_nproc!r}")
        if not isinstance(y_over, bool):
            raise GuestConfigError(f"fragmentation.overwrite must be true or false, got {y_over!r}")
        nproc = a.nproc if a.nproc is not None else y_nproc
        overwrite = a.overwrite if a.overwrite is not None else y_over
        run(cfg, run_dir_of(input_path), input_path, nproc, overwrite)
    except GuestConfigError as e:
        print(f"error: {e}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
