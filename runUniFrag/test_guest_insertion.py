#!/usr/bin/env python3
"""Tests for guest_insertion.py (no pytest needed).

Usage:  python runUniFrag/test_guest_insertion.py [name-substring]

Runs every test_* function, prints PASS/FAIL per test and a total, and exits 1
on any failure. Hosts are small synthetic clusters; nothing is fragmented.
"""
import os
import sys
import shutil
import tempfile
import traceback

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, ROOT)

import guest_insertion as gi  # noqa: E402


# --------------------------------------------------------------------- helpers
def pdist(x):
    x = np.asarray(x, float)
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
    return d[np.triu_indices(len(x), 1)]


def com(g):
    from ase.data import atomic_masses, atomic_numbers
    m = np.array([atomic_masses[atomic_numbers[s]] for s in g.symbols])
    return (g.positions * m[:, None]).sum(0) / m.sum()


def raises(exc, fn, *a, **k):
    try:
        fn(*a, **k)
    except exc:
        return True
    return False


# ------------------------------------------------------------------- Task 1
def test_rotation_is_proper():
    rng = np.random.default_rng(0)
    for _ in range(200):
        M = gi.random_rotation(rng)
        assert np.allclose(M @ M.T, np.eye(3), atol=1e-12)
        assert abs(np.linalg.det(M) - 1.0) < 1e-12


def test_rotation_keeps_pair_distances():
    g = gi.load_guest("CH4", ".")
    d0 = pdist(g.positions)
    M = gi.random_rotation(np.random.default_rng(1))
    assert np.allclose(pdist(g.positions @ M.T), d0, atol=1e-12)


def test_load_guest_builtin_centred():
    g = gi.load_guest("CO2", ".")
    assert g.symbols == ("C", "O", "O")
    assert g.name == "CO2"
    assert np.allclose(com(g), 0.0, atol=1e-10)


def test_load_guest_bad_names():
    for bad in ["CO_2", "CO-2", "NotAMolecule"]:
        assert raises(gi.GuestConfigError, gi.load_guest, bad, "."), bad


def test_radius_fallback_and_override():
    t = gi.radius_table({}, 0.35, 1.0)
    assert len(t) == 119
    assert abs(t[24] - (2.45 - 0.35)) < 1e-9      # Cr: Alvarez
    assert abs(t[40] - (2.52 - 0.35)) < 1e-9      # Zr: Alvarez
    assert abs(t[6] - (1.70 - 0.35)) < 1e-9       # C: ase default kept
    t2 = gi.radius_table({"Zr": 2.2}, 0.35, 0.9)
    assert abs(t2[40] - (2.2 - 0.35) * 0.9) < 1e-9
    assert abs(t2[6] - (1.70 - 0.35) * 0.9) < 1e-9


def test_radii_for_missing_element():
    t = gi.radius_table({}, 0.35, 1.0)
    assert raises(gi.GuestConfigError, gi.radii_for, ["C", "Pm"], t)
    assert np.allclose(gi.radii_for(["C", "Cr"], t), [1.35, 2.10])


def test_box_margin_and_planar_stretch():
    cube = np.array([[x, y, z] for x in (0, 10) for y in (0, 10) for z in (0, 10)], float)
    lo, hi = gi.insertion_box(cube, -1.0, 8.0)
    assert np.allclose(lo, 1.0) and np.allclose(hi, 9.0)
    flat = np.array([[0, 0, 0], [10, 0, 0], [0, 10, 0], [10, 10, 0]], float)
    lo, hi = gi.insertion_box(flat, -1.0, 8.0)
    assert np.allclose(hi - lo, [8.0, 8.0, 8.0])
    assert np.allclose((hi + lo) / 2, [5.0, 5.0, 0.0])


# ------------------------------------------------------------------- Task 2
# Single-atom fake guests with hand-set radii make the fill deterministic:
# r(Ar) = -10, r(Kr) = +10, r(Xe) = +10. A-A and A-B thresholds are <= 0
# (never overlap), B-B and B-host thresholds are 20 A (always overlap inside
# a box whose diagonal is < 20 A).
A = gi.Guest("A", ("Ar",), np.zeros((1, 3)))
B = gi.Guest("B", ("Kr",), np.zeros((1, 3)))


def _fake_table():
    t = gi.radius_table({}, 0.35, 1.0)
    t[18], t[36], t[54] = -10.0, 10.0, 10.0
    return t


def test_overlap_rule_boundary():
    r = np.array([1.0])
    xyz = np.zeros((1, 3))
    assert gi.overlaps(np.array([[1.99, 0, 0]]), r, xyz, r)
    assert not gi.overlaps(np.array([[2.01, 0, 0]]), r, xyz, r)
    assert not gi.overlaps(np.array([[0.5, 0, 0]]), r, np.zeros((0, 3)), np.zeros(0))


def test_fill_trims_to_whole_units():
    # 5 A, 1 B, 5 more A always fit; a second B never does -> 10 A + 1 B placed,
    # trimmed back to one whole unit 5 A + 1 B.
    box = (np.zeros(3), np.full(3, 5.0))
    res = gi.fill_copy([], np.zeros((0, 3)), [A] * 5 + [B], _fake_table(), box,
                       200, None, np.random.default_rng(3))
    assert res.status == "written"
    assert res.units == 1 and res.counts == {"A": 5, "B": 1}
    assert res.symbols == ["Ar"] * 5 + ["Kr"]
    assert res.coords.shape == (6, 3)


def test_fill_skips_when_no_unit_fits():
    # A host Xe atom in the middle blocks every B: 5 A fit, B never -> 0 units.
    box = (np.zeros(3), np.full(3, 5.0))
    res = gi.fill_copy(["Xe"], np.full((1, 3), 2.5), [A] * 5 + [B], _fake_table(),
                       box, 200, None, np.random.default_rng(4))
    assert res.status == "skipped: no whole unit fits"
    assert res.units == 0 and res.counts == {"A": 0, "B": 0}
    assert len(res.symbols) == 0 and res.coords.shape == (0, 3)


def test_fill_respects_max_units():
    co2 = gi.load_guest("CO2", ".")
    t = gi.radius_table({}, 0.35, 1.0)
    box = (np.zeros(3), np.full(3, 40.0))
    res = gi.fill_copy([], np.zeros((0, 3)), [co2], t, box, 1000, 3,
                       np.random.default_rng(5))
    assert res.units == 3 and res.counts == {"CO2": 3}
    assert len(res.symbols) == 9


def test_fill_single_gas_fills_until_full():
    co2 = gi.load_guest("CO2", ".")
    t = gi.radius_table({}, 0.35, 1.0)
    box = (np.zeros(3), np.full(3, 8.0))
    res = gi.fill_copy([], np.zeros((0, 3)), [co2], t, box, 300, None,
                       np.random.default_rng(6))
    assert res.status == "written" and res.units >= 2
    assert res.counts == {"CO2": res.units}


def test_fill_result_has_no_contacts():
    co2, ch4 = gi.load_guest("CO2", "."), gi.load_guest("CH4", ".")
    t = gi.radius_table({}, 0.35, 1.0)
    ang = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    host = np.c_[1.4 * np.cos(ang), 1.4 * np.sin(ang), np.zeros(6)]
    hsym = ["C"] * 6
    box = gi.insertion_box(host, 4.0, 8.0)
    res = gi.fill_copy(hsym, host, [co2] * 5 + [ch4], t, box, 500, None,
                       np.random.default_rng(7))
    assert res.units >= 1
    allsym = hsym + res.symbols
    allxyz = np.vstack([host, res.coords])
    r = gi.radii_for(allsym, t)
    # guest atoms, in insertion order: CO2 x5 then CH4, repeated
    sizes = ([3] * 5 + [5]) * res.units
    starts = np.cumsum([6] + sizes)[:-1]
    for s, n in zip(starts, sizes):
        mine = np.arange(s, s + n)
        others = np.setdiff1d(np.arange(len(allsym)), mine)
        d = np.linalg.norm(allxyz[mine][:, None] - allxyz[others][None], axis=-1)
        assert (d >= r[mine][:, None] + r[others][None] - 1e-9).all()


# ------------------------------------------------------------------- Task 3
BASE = {"enabled": True,
        "single": {"enabled": True, "species": ["CO2", "CH4"]},
        "mixture": {"enabled": True, "species": [{"name": "CO2", "ratio": 5},
                                                  {"name": "CH4", "ratio": 1}]}}

CO2_XYZ = "3\nCO2\nC 0 0 0\nO 0 0 1.16\nO 0 0 -1.16\n"


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def test_defaults_and_sets():
    c = gi.parse_guests(BASE, ".")
    assert [s.set_id for s in c.sets] == ["single:CO2", "single:CH4", "mixture:CO2-CH4:5-1"]
    assert [s.kind for s in c.sets] == ["single", "single", "mixture"]
    assert [s.tag for s in c.sets] == ["CO2", "CH4", "CO2x5CH4x1"]
    assert c.sets[0].output == "fragments_guests_single_CO2.extxyz"
    assert c.sets[2].output == "fragments_guests_mixture_CO2-CH4_5-1.extxyz"
    assert c.sets[2].ratios == (5, 1) and c.sets[0].ratios == (1,)
    assert [g.name for g in c.sets[2].unit()] == ["CO2"] * 5 + ["CH4"]
    assert [g.name for g in c.sets[1].unit()] == ["CH4"]
    assert (c.copies_per_frame, c.box_margin, c.min_box_extent, c.vdw_shift, c.vdw_scale,
            c.vdw_radii, c.full_after_failures, c.max_units, c.seed, c.frames) == \
        (5, -1.0, 8.0, 0.35, 1.0, {}, 1000, None, 12345, frozenset({"normal", "min"}))
    sh = c.shared()
    assert sh["frames"] == ["min", "normal"] and "sets" not in sh
    assert c.set_specs()["mixture:CO2-CH4:5-1"] == [["CO2", 5, None], ["CH4", 1, None]]


def test_disabled_or_absent_returns_none():
    assert gi.parse_guests(None, ".") is None
    assert gi.parse_guests({"enabled": False, "bogus": 1}, ".") is None


def test_only_one_block():
    c = gi.parse_guests({"enabled": True, "mixture": BASE["mixture"]}, ".")
    assert [s.set_id for s in c.sets] == ["mixture:CO2-CH4:5-1"]


def test_config_errors():
    two = lambda r1: {"enabled": True, "species": [{"name": "CO2", "ratio": r1},
                                                   {"name": "CH4", "ratio": 1}]}
    bad = [
        {**BASE, "vdw_scal": 1.0},                                   # unknown key
        {**BASE, "mixture": two(2.5)},                               # non-integer ratio
        {**BASE, "mixture": two(True)},                              # bool is not a ratio
        {**BASE, "mixture": two(0)},
        {**BASE, "mixture": {"enabled": True, "species": [{"name": "CO2", "ratio": 1}]}},
        {"enabled": True, "single": {"enabled": False}, "mixture": {"enabled": False}},
        {"enabled": True},
        {**BASE, "single": {"enabled": True, "species": []}},
        {**BASE, "single": {"enabled": True, "species": ["CO2"], "colour": 1}},
        {**BASE, "frames": ["normal", "big"]},
        {**BASE, "frames": []},
        {**BASE, "copies_per_frame": 0},
        {**BASE, "full_after_failures": 0},
        {**BASE, "max_units": 0},
        {**BASE, "vdw_scale": 0},
        {**BASE, "min_box_extent": -1},
        {**BASE, "seed": 1.5},
        {**BASE, "enabled": "yes"},
    ]
    for b in bad:
        assert raises(gi.GuestConfigError, gi.parse_guests, b, "."), b


def test_frames_accepts_all_four():
    c = gi.parse_guests({**BASE, "frames": ["normal", "min", "onlylinker", "onlynode"]}, ".")
    assert c.frames == frozenset({"normal", "min", "onlylinker", "onlynode"})


def test_guest_file_relative_to_yaml(tmp):
    os.makedirs(f"{tmp}/geoms")
    write(f"{tmp}/geoms/co2.xyz", CO2_XYZ)
    c = gi.parse_guests({"enabled": True,
                         "single": {"enabled": True,
                                    "species": [{"name": "CO2", "file": "geoms/co2.xyz"}]}},
                        base_dir=tmp)
    assert c.sets[0].guests[0].symbols == ("C", "O", "O")
    path = c.set_specs()["single:CO2"][0][2]
    assert os.path.isabs(path) and path.endswith("geoms/co2.xyz")


def test_run_yaml_top_level(tmp):
    write(f"{tmp}/r.yaml", "fragmentation: {kind: mof}\nguest: {}\n")
    assert raises(gi.GuestConfigError, gi.load_run_yaml, f"{tmp}/r.yaml")
    write(f"{tmp}/ok.yaml", "fragmentation: {kind: mof}\n")
    assert gi.load_run_yaml(f"{tmp}/ok.yaml") == {"fragmentation": {"kind": "mof"}}
    write(f"{tmp}/empty.yaml", "")
    assert gi.load_run_yaml(f"{tmp}/empty.yaml") == {}


# ------------------------------------------------------------------- Task 4
def _fo():
    import fragmentation_oop as fo
    return fo


def ring(label, z=0.0, extra="", capped_h=None):
    ang = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    c = np.c_[1.4 * np.cos(ang), 1.4 * np.sin(ang), np.full(6, z)]
    h = np.c_[2.48 * np.cos(ang), 2.48 * np.sin(ang), np.full(6, z)]
    return {"label": label, "species": ["C"] * 6 + ["H"] * 6,
            "coords": np.vstack([c, h]), "extra": extra, "capped_h": capped_h}


def dense_host(label):
    g = np.arange(0.0, 12.01, 1.5)
    xyz = np.array([[x, y, z] for x in g for y in g for z in g])
    return {"label": label, "species": ["C"] * len(xyz), "coords": xyz,
            "extra": "", "capped_h": None}


STD_FRAMES = [ring("XFragMof", extra="parent_charge=unknown", capped_h=[6, 7]),
              ring("XFragMofMin", capped_h=[]),
              ring("XFragMofOnlyLinker"),
              ring("COF-TpAzoFragCof", capped_h=[11])]

CFG_DICT = {**BASE, "full_after_failures": 200}


def cfg():
    return gi.parse_guests(CFG_DICT, ".")


def append_frames(d, frames):
    for f in frames:
        _fo()._write_extxyz(os.path.join(d, "fragments_collection.extxyz"), f["species"],
                            f["coords"], f["label"], extra=f["extra"], capped_h=f["capped_h"])


def make_run_dir(tmp, frames, name="run"):
    d = os.path.join(tmp, name)
    os.makedirs(d)
    append_frames(d, frames)
    return d


def run_copy(tmp, name, frames=STD_FRAMES, nproc=1, c=None):
    d = make_run_dir(tmp, frames, name)
    gi.run(c or cfg(), d, d, nproc, False)
    return d


def guest_files(d):
    return sorted(f for f in os.listdir(d) if f.startswith("fragments_guests_"))


def read_bytes_all_outputs(d):
    names = guest_files(d) + ["guests_summary.csv", "guests_config.json"]
    return {n: open(os.path.join(d, n), "rb").read() for n in names}


def all_guest_frames(d):
    return [f for n in guest_files(d) for f in _fo()._parse_extxyz(os.path.join(d, n))]


def by_label(d, label):
    return next(f for f in _fo()._parse_extxyz(os.path.join(d, "fragments_collection.extxyz"))
                if f["label"] == label)


def csv_rows(d):
    import csv
    with open(os.path.join(d, "guests_summary.csv")) as f:
        return list(csv.DictReader(f))


def drop_csv_rows_for(d, parent):
    import csv
    rows = csv_rows(d)
    with open(os.path.join(d, "guests_summary.csv")) as f:
        header = f.readline().strip().split(",")
    with open(os.path.join(d, "guests_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, header)
        w.writeheader()
        w.writerows(r for r in rows if r["parent"] != parent)


def test_frame_class():
    assert gi.frame_class("XFragMof") == "normal"
    assert gi.frame_class("XFragCof") == "normal"
    assert gi.frame_class("XFragMofMin") == "min"
    assert gi.frame_class("XFragMofOnlyLinker2") == "onlylinker"
    assert gi.frame_class("XFragCofOnlyNode") == "onlynode"
    assert gi.frame_class("XW3") is None
    assert gi.frame_class("XFragMof-CO2-0") is None


def test_frame_class_with_dash_in_base():
    assert gi.frame_class("COF-TpAzoFragCof") == "normal"
    assert gi.frame_class("COF-TpAzoFragCofMin") == "min"


def test_parse_header_tags():
    t = gi.parse_header_tags('Properties=species:S:1:pos:R:3 label=X charge=0 '
                             'capped_h="1 2" pbc="F F F"')
    assert t == {"Properties": "species:S:1:pos:R:3", "label": "X", "charge": "0",
                 "capped_h": "1 2", "pbc": "F F F"}


def test_guest_header_copies_parent_tags():
    h = gi.guest_header('Properties=species:S:1:pos:R:3 label=P charge=0 '
                        'parent_charge=unknown capped_h="" pbc="F F F"',
                        "P-CO2-1", 0, "P", {"CO2": 4}, 12)
    assert h == ('Properties=species:S:1:pos:R:3 label=P-CO2-1 charge=0 parent=P '
                 'guests="CO2:4" n_host=12 parent_charge=unknown capped_h="" pbc="F F F"')


def test_run_outputs_and_labels(tmp):
    d = make_run_dir(tmp, STD_FRAMES)
    out = gi.run(cfg(), d, d, nproc=1, overwrite=False)
    assert guest_files(d) == ["fragments_guests_mixture_CO2-CH4_5-1.extxyz",
                              "fragments_guests_single_CH4.extxyz",
                              "fragments_guests_single_CO2.extxyz"]
    labels = [f["label"] for f in _fo()._parse_extxyz(f"{d}/fragments_guests_single_CO2.extxyz")]
    assert {"XFragMof-CO2-0", "XFragMofMin-CO2-4", "COF-TpAzoFragCof-CO2-2"} <= set(labels)
    assert not any(l.startswith("XFragMofOnlyLinker") for l in labels)
    for f in all_guest_frames(d):
        t = gi.parse_header_tags(f["comment"])
        assert "_" not in t["label"] and t["charge"] == "0"
        n = int(t["n_host"])
        host = by_label(d, t["parent"])
        assert f["species"][:n] == host["species"]
        assert np.allclose(f["coords"][:n], host["coords"])
        assert t.get("capped_h") == gi.parse_header_tags(host["comment"]).get("capped_h")
        if t["parent"] == "XFragMof":
            assert t["parent_charge"] == "unknown"
    for sid, (w, s) in out.items():
        assert w + s == 3 * 5, (sid, w, s)
    assert len(csv_rows(d)) == 3 * 3 * 5
    mix = [f for f in _fo()._parse_extxyz(f"{d}/fragments_guests_mixture_CO2-CH4_5-1.extxyz")]
    for f in mix:
        t = gi.parse_header_tags(f["comment"])
        c = dict(x.split(":") for x in t["guests"].split())
        assert int(c["CO2"]) == 5 * int(c["CH4"]) and int(c["CH4"]) >= 1


def test_nproc_independent(tmp):
    a = run_copy(tmp, "a", nproc=1)
    b = run_copy(tmp, "b", nproc=2)
    assert read_bytes_all_outputs(a) == read_bytes_all_outputs(b)


def test_resume_equals_uninterrupted(tmp):
    full = run_copy(tmp, "full")
    part = run_copy(tmp, "part", frames=STD_FRAMES[:1])
    append_frames(part, STD_FRAMES[1:])
    gi.run(cfg(), part, part, 1, False)
    assert read_bytes_all_outputs(full) == read_bytes_all_outputs(part)


def test_resume_drops_orphan_frames(tmp):
    d = run_copy(tmp, "o")
    drop_csv_rows_for(d, parent="XFragMofMin")          # crash before the CSV write
    gi.run(cfg(), d, d, 1, False)
    labels = [f["label"] for f in all_guest_frames(d)]
    assert len(labels) == len(set(labels)) == 3 * 3 * 5 - sum(
        1 for r in csv_rows(d) if r["status"] != "written")
    assert len(csv_rows(d)) == 3 * 3 * 5


def test_resume_settings_mismatch(tmp):
    d = run_copy(tmp, "m")
    changed = gi.parse_guests({**CFG_DICT, "box_margin": 0.0}, ".")
    assert raises(gi.GuestConfigError, gi.run, changed, d, d, 1, False)
    gi.run(changed, d, d, 1, True)
    assert len(csv_rows(d)) == 3 * 3 * 5


def test_new_ratio_is_new_file_not_mismatch(tmp):
    d = run_copy(tmp, "r")
    other = {**CFG_DICT, "mixture": {"enabled": True,
                                     "species": [{"name": "CO2", "ratio": 3},
                                                 {"name": "CH4", "ratio": 1}]}}
    gi.run(gi.parse_guests(other, "."), d, d, 1, False)
    assert os.path.exists(f"{d}/fragments_guests_mixture_CO2-CH4_3-1.extxyz")
    assert os.path.exists(f"{d}/fragments_guests_mixture_CO2-CH4_5-1.extxyz")
    assert len(csv_rows(d)) == 4 * 3 * 5                   # single x2 (done) + two mixtures


def test_new_species_extends_csv(tmp):
    d = run_copy(tmp, "n")
    other = {**CFG_DICT, "single": {"enabled": True, "species": ["CO2", "N2"]}}
    gi.run(gi.parse_guests(other, "."), d, d, 1, False)
    rows = csv_rows(d)
    assert "N2" in rows[0] and rows[0]["N2"] == "0"
    assert sum(1 for r in rows if r["set"] == "single:N2") == 3 * 5


def test_crowded_host_skips_all_copies(tmp):
    d = make_run_dir(tmp, [dense_host("YFragMof")])
    out = gi.run(cfg(), d, d, 1, False)
    assert all(w == 0 and s == 5 for w, s in out.values())
    rows = csv_rows(d)
    assert len(rows) == 3 * 5
    assert all(r["status"] == "skipped: no whole unit fits" for r in rows)


def test_single_file_input_only_its_frames(tmp):
    d = make_run_dir(tmp, STD_FRAMES + [ring("ZFragMof")])
    gi.run(cfg(), d, os.path.join(d, "Z.cif"), 1, False)
    parents = {gi.parse_header_tags(f["comment"])["parent"] for f in all_guest_frames(d)}
    assert parents == {"ZFragMof"}


def test_standalone_main(tmp):
    d = make_run_dir(tmp, STD_FRAMES[:1])
    write(f"{tmp}/r.yaml", f"fragmentation: {{input_path: {d}, kind: mof, nproc: 1}}\n"
                           "guests: {enabled: true, full_after_failures: 200, "
                           "single: {enabled: true, species: [CO2]}}\n")
    assert gi.main([f"{tmp}/r.yaml"]) == 0
    assert guest_files(d) == ["fragments_guests_single_CO2.extxyz"]
    write(f"{tmp}/b.yaml", f"fragmentation: {{input_path: {d}, kind: bio}}\n"
                           "guests: {enabled: true, single: {enabled: true, species: [CO2]}}\n")
    assert gi.main([f"{tmp}/b.yaml"]) == 2
    write(f"{tmp}/off.yaml", f"fragmentation: {{input_path: {d}}}\n")
    assert gi.main([f"{tmp}/off.yaml"]) == 2


# ------------------------------------------------------------------- Task 5
def exits(code, fn, *a, **k):
    import contextlib
    import io
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            fn(*a, **k)
    except SystemExit as e:
        return e.code == code
    return False


def test_config_absent_is_unchanged():
    fo = _fo()
    a, g = fo._parse_args_with_config(fo._build_parser(), ["x.cif"])
    assert g is None and a.input_path == "x.cif" and a.kind == "mof"
    assert (a.radius, a.timeout, a.nproc, a.center, a.nmetals, a.cof_layer, a.overwrite) == \
        (6.0, 300.0, 1, -1, 3, "auto", False)
    assert a.config is None


def test_config_defaults_and_cli_override(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, kind: cof, nproc: 8, "
                           "cof_layer: dimer, overwrite: true, radius: 4.0}\n")
    a, g = fo._parse_args_with_config(fo._build_parser(),
                                      ["--config", f"{tmp}/r.yaml", "--nproc", "1"])
    assert (a.input_path, a.kind, a.nproc, a.cof_layer, a.overwrite, a.radius, g) == \
        ("d", "cof", 1, "dimer", True, 4.0, None)
    a, _ = fo._parse_args_with_config(fo._build_parser(),
                                      ["other", "--config", f"{tmp}/r.yaml"])
    assert a.input_path == "other"


def test_config_guests_parsed(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d}\n"
                           "guests: {enabled: true, single: {enabled: true, species: [CO2]}}\n")
    a, g = fo._parse_args_with_config(fo._build_parser(), ["--config", f"{tmp}/r.yaml"])
    assert [s.set_id for s in g.sets] == ["single:CO2"]


def test_config_unknown_fragmentation_key(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, radious: 4}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])


def test_config_bad_choice_rejected(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, kind: mfo}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])


def test_config_bad_guest_section(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d}\nguests: {enabled: true}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])


def test_config_bio_with_guests_rejected(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {input_path: d, kind: bio}\n"
                           "guests: {enabled: true, single: {enabled: true, species: [CO2]}}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])


def test_config_missing_input_path(tmp):
    fo = _fo()
    write(f"{tmp}/r.yaml", "fragmentation: {kind: mof}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), [])


# ------------------------------------------------------------------- Task 6
def _chk():
    sys.path.insert(0, os.path.join(ROOT, "runUniFrag"))
    import check_guest_fragments as chk
    return chk


def _edit_frame(d, name, frame, fn):
    """Rewrite file `name` with fn(species, coords, comment) applied to one frame."""
    path = os.path.join(d, name)
    frames = _fo()._parse_extxyz(path)
    os.remove(path)
    for i, f in enumerate(frames):
        sp, co, cm = list(f["species"]), np.asarray(f["coords"], float), f["comment"]
        if i == frame:
            sp, co, cm = fn(sp, co, cm)
        gi._write_frame(path, sp, co, cm)


def _n_host(cm):
    return int(gi.parse_header_tags(cm)["n_host"])


def test_checker_passes_clean_run(tmp):
    d = run_copy(tmp, "c")
    rep = _chk().check(d)
    assert rep["frames"] == 3 * 3 * 5 - sum(1 for r in csv_rows(d) if r["status"] != "written")
    assert all(not v for k, v in rep.items() if k != "frames"), rep


def test_checker_passes_file_guest(tmp):
    write(f"{tmp}/co2.xyz", CO2_XYZ)
    c = gi.parse_guests({"enabled": True, "full_after_failures": 200,
                         "single": {"enabled": True,
                                    "species": [{"name": "CO2", "file": "co2.xyz"}]}}, tmp)
    d = run_copy(tmp, "f", c=c)
    rep = _chk().check(d)
    assert all(not v for k, v in rep.items() if k != "frames"), rep


def test_checker_catches_moved_host_atom(tmp):
    d = run_copy(tmp, "h")

    def nudge(sp, co, cm):
        co = co.copy()
        co[0, 0] += 0.01
        return sp, co, cm
    _edit_frame(d, "fragments_guests_single_CO2.extxyz", 0, nudge)
    assert _chk().check(d)["G1"]


def test_checker_catches_contact(tmp):
    d = run_copy(tmp, "k")

    def onto_host(sp, co, cm):
        n = _n_host(cm)
        co = co.copy()
        co[n:n + 5] += co[0] - co[n]                   # first CH4's C onto host atom 0
        return sp, co, cm
    _edit_frame(d, "fragments_guests_single_CH4.extxyz", 0, onto_host)
    assert _chk().check(d)["G2"]


def test_checker_catches_bent_guest(tmp):
    d = run_copy(tmp, "b")

    def bend(sp, co, cm):
        n = _n_host(cm)
        co = co.copy()
        co[n + 1] += [0.0, 0.0, 0.05]                  # stretch one C=O of the first CO2
        return sp, co, cm
    _edit_frame(d, "fragments_guests_single_CO2.extxyz", 0, bend)
    assert _chk().check(d)["G3"]


def test_checker_catches_bad_ratio(tmp):
    d = run_copy(tmp, "q")

    def drop_last(sp, co, cm):
        return sp[:-5], co[:-5], cm                     # remove the last CH4
    _edit_frame(d, "fragments_guests_mixture_CO2-CH4_5-1.extxyz", 0, drop_last)
    assert _chk().check(d)["G4"]


def test_checker_catches_charge_and_label(tmp):
    d = run_copy(tmp, "z")

    def bad_tags(sp, co, cm):
        return sp, co, cm.replace("charge=0", "charge=1").replace("-CO2-0", "_CO2_0", 1)
    _edit_frame(d, "fragments_guests_single_CO2.extxyz", 0, bad_tags)
    rep = _chk().check(d)
    assert rep["G5"] and rep["G6"]


# --------------------------------------------------------- final-review fixes
def test_yaml_values_are_type_checked(tmp):
    fo = _fo()
    for bad in ['overwrite: "false"', "overwrite: 1", "minimize: \"no\"", "nproc: 8.5",
                "nproc: true", "radius: abc", "chain: 5"]:
        write(f"{tmp}/r.yaml", f"fragmentation: {{input_path: d, {bad}}}\n")
        assert exits(2, fo._parse_args_with_config, fo._build_parser(),
                     ["--config", f"{tmp}/r.yaml"]), bad
    write(f"{tmp}/ok.yaml", "fragmentation: {input_path: d, nproc: 8, radius: 4, chain: null, "
                            "overwrite: false}\n")
    a, _ = fo._parse_args_with_config(fo._build_parser(), ["--config", f"{tmp}/ok.yaml"])
    assert a.nproc == 8 and a.radius == 4.0 and isinstance(a.radius, float)
    assert a.overwrite is False and a.chain is None


def test_standalone_rejects_string_bool(tmp):
    d = make_run_dir(tmp, STD_FRAMES[:1])
    write(f"{tmp}/r.yaml", f"fragmentation: {{input_path: {d}, overwrite: \"false\"}}\n"
                           "guests: {enabled: true, single: {enabled: true, species: [CO2]}}\n")
    assert gi.main([f"{tmp}/r.yaml"]) == 2


def test_single_file_overwrite_keeps_other_structures(tmp):
    d = make_run_dir(tmp, [ring("XFragMof"), ring("ZFragMof")])
    gi.run(cfg(), d, d, 1, False)
    before = [f["label"] for f in all_guest_frames(d) if f["label"].startswith("X")]
    gi.run(cfg(), d, os.path.join(d, "Z.cif"), 1, True)
    after = [f["label"] for f in all_guest_frames(d) if f["label"].startswith("X")]
    assert before and sorted(after) == sorted(before)
    assert sum(1 for r in csv_rows(d) if r["parent"] == "XFragMof") == 3 * 5
    assert sum(1 for r in csv_rows(d) if r["parent"] == "ZFragMof") == 3 * 5


def test_single_file_rerun_regenerates_its_frames(tmp):
    # A single-file re-fragmentation replaces that structure's collection frames;
    # its guest frames must follow the new host, not stay from the old one.
    d = make_run_dir(tmp, [ring("XFragMof"), ring("ZFragMof")])
    gi.run(cfg(), d, d, 1, False)
    coll = os.path.join(d, "fragments_collection.extxyz")
    os.remove(coll)
    append_frames(d, [ring("XFragMof"), ring("ZFragMof", z=0.5)])
    gi.run(cfg(), d, os.path.join(d, "Z.cif"), 1, False)
    rep = _chk().check(d)
    assert not rep["G1"], rep["G1"][:3]
    assert len(csv_rows(d)) == 2 * 3 * 5


def test_settings_mismatch_caught_before_fragmentation(tmp):
    fo = _fo()
    d = run_copy(tmp, "e")
    write(f"{tmp}/r.yaml", f"fragmentation: {{input_path: {d}}}\n"
                           "guests: {enabled: true, box_margin: 0.0, full_after_failures: 200, "
                           "single: {enabled: true, species: [CO2, CH4]}, "
                           "mixture: {enabled: true, species: [{name: CO2, ratio: 5}, "
                           "{name: CH4, ratio: 1}]}}\n")
    assert exits(2, fo._parse_args_with_config, fo._build_parser(), ["--config", f"{tmp}/r.yaml"])
    # with --overwrite the user asked to start again: no error
    fo._parse_args_with_config(fo._build_parser(), ["--config", f"{tmp}/r.yaml", "--overwrite"])


def test_mismatch_message_points_at_guest_only_overwrite(tmp):
    d = run_copy(tmp, "msg")
    changed = gi.parse_guests({**CFG_DICT, "box_margin": 0.0}, ".")
    try:
        gi.run(changed, d, d, 1, False)
        assert False, "no error"
    except gi.GuestConfigError as e:
        assert "box_margin" in str(e) and "guest_insertion.py" in str(e) and "--overwrite" in str(e)


def test_mixture_labels_unique_across_ratios(tmp):
    d = run_copy(tmp, "u")
    other = {**CFG_DICT, "mixture": {"enabled": True,
                                     "species": [{"name": "CO2", "ratio": 3},
                                                 {"name": "CH4", "ratio": 1}]}}
    gi.run(gi.parse_guests(other, "."), d, d, 1, False)
    labels = [f["label"] for f in all_guest_frames(d)]
    assert len(labels) == len(set(labels))
    assert gi.parse_guests(BASE, ".").sets[2].tag == "CO2x5CH4x1"
    assert not _chk().check(d)["G6"]


def test_checker_catches_duplicate_label(tmp):
    d = run_copy(tmp, "dup")

    def relabel(sp, co, cm):
        return sp, co, cm.replace("label=XFragMof-CO2-1", "label=XFragMof-CO2-0")
    _edit_frame(d, "fragments_guests_single_CO2.extxyz", 1, relabel)
    assert _chk().check(d)["G6"]


def test_empty_host_is_skipped(tmp):
    d = make_run_dir(tmp, [{"label": "EFragMof", "species": [], "coords": np.zeros((0, 3)),
                            "extra": "", "capped_h": None}])
    out = gi.run(cfg(), d, d, 1, False)
    assert all(w == 0 and s == 5 for w, s in out.values())
    assert all(r["status"] == "skipped: empty host" for r in csv_rows(d))


def test_missing_collection_is_an_error(tmp):
    os.makedirs(f"{tmp}/empty")
    assert raises(gi.GuestConfigError, gi.run, cfg(), f"{tmp}/empty", f"{tmp}/empty", 1, False)
    assert not os.path.exists(f"{tmp}/empty/guests_config.json")
    write(f"{tmp}/r.yaml", f"fragmentation: {{input_path: {tmp}/nope}}\n"
                           "guests: {enabled: true, single: {enabled: true, species: [CO2]}}\n")
    assert gi.main([f"{tmp}/r.yaml"]) == 2


# ---------------------------------------------------------------------- run
def main():
    sel = sys.argv[1] if len(sys.argv) > 1 else ""
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f) and sel in n]
    fails = 0
    for name, fn in tests:
        tmp = tempfile.mkdtemp(prefix="gi_test_")
        try:
            if fn.__code__.co_argcount:
                fn(tmp)
            else:
                fn()
            print(f"PASS {name}")
        except Exception:
            fails += 1
            print(f"FAIL {name}")
            traceback.print_exc()
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    print(f"PASS {len(tests) - fails} / FAIL {fails}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
