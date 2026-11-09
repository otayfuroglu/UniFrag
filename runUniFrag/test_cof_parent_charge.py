#!/usr/bin/env python3
"""COF parent charge (GitHub #11, follow-up).

    python runUniFrag/test_cof_parent_charge.py [substring]

1. A COF parent drawn with a bare sp3 carbon (two single bonds, no H: a CH2
   whose hydrogens the model omitted, CoRE-COF 141/142/144/145) gets its two
   H before fragmenting. Nothing else is touched: an aromatic C that lost its
   H, an alkyne C, a CH and a site with no room keep their atoms as drawn.
2. The parent's net formal charge is judged (None when it has no plausible
   closed-shell reading), and each fragment's verdict is: odd electron first,
   then formal charge not 0 -> quarantine=charge; parent_charge= on every
   frame of a parent that is not neutral.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentation_oop import COFFragmenter, _cof_frame_verdict  # noqa: E402
from pymatgen.core import Lattice, Structure  # noqa: E402

BOX = 30.0
F = COFFragmenter()


def box(sp, xyz):
    c = np.asarray(xyz, float) + BOX / 2
    return Structure(Lattice.cubic(BOX), sp, c, coords_are_cartesian=True)


def methyl_h(c, away):
    """Three H on carbon c, staggered, pointing away from `away`."""
    c, away = np.asarray(c, float), np.asarray(away, float)
    u = (c - away) / np.linalg.norm(c - away)
    p = np.cross(u, [0, 0, 1.0])
    if np.linalg.norm(p) < 1e-6:
        p = np.cross(u, [0, 1.0, 0])
    p /= np.linalg.norm(p)
    q = np.cross(u, p)
    out = []
    for k in range(3):
        t = np.deg2rad(120 * k)
        d = 0.334 * u + 0.943 * (np.cos(t) * p + np.sin(t) * q)
        out.append(c + 1.09 * d / np.linalg.norm(d))
    return out


def propane(middle_h=0, block=False):
    """CH3-C(Hn)-CH3, C-C 1.54 A, C-C-C 111 deg; the middle C carries
    `middle_h` H. `block` puts an O where an added H would go."""
    a = np.deg2rad(111.0) / 2
    c2 = np.zeros(3)
    c1 = 1.54 * np.array([-np.sin(a), -np.cos(a), 0.0])
    c3 = 1.54 * np.array([np.sin(a), -np.cos(a), 0.0])
    sp = ["C", "C", "C"]
    xyz = [c1, c2, c3]
    for h in methyl_h(c1, c2) + methyl_h(c3, c2):
        sp.append("H")
        xyz.append(h)
    b = np.array([0.0, 1.0, 0.0])
    n = np.array([0.0, 0.0, 1.0])
    tet = [c2 + 1.09 * (0.577 * b + s * 0.816 * n) for s in (1, -1)]
    for h in tet[:middle_h]:
        sp.append("H")
        xyz.append(h)
    if block:
        sp.append("O")
        xyz.append(tet[0] + np.array([0.0, 0.0, 0.3]))
    return sp, xyz


def test_bare_sp3_ch2_gets_two_h():
    sp, xyz = propane()
    s = box(sp, xyz)
    s2, n = F._repair_bare_sp3_ch2(s)
    assert n == 1, n
    assert len(s2) == len(s) + 2
    for i in range(len(s)):                 # parent atoms never move
        assert np.allclose(s2[i].coords, s[i].coords)
    hs = [s2[k].coords for k in range(len(s), len(s2))]
    for h in hs:
        assert abs(np.linalg.norm(h - s2[1].coords) - 1.09) < 0.01
    u = [(h - s2[1].coords) / np.linalg.norm(h - s2[1].coords) for h in hs]
    ang = np.degrees(np.arccos(np.dot(u[0], u[1])))
    assert 104 < ang < 115, ang
    assert F.cof_parent_formal_charge(s2) == 0


def test_ch_is_left_alone():
    sp, xyz = propane(middle_h=1)
    _s2, n = F._repair_bare_sp3_ch2(box(sp, xyz))
    assert n == 0


def test_no_room_leaves_it_alone():
    sp, xyz = propane(block=True)
    s = box(sp, xyz)
    s2, n = F._repair_bare_sp3_ch2(s)
    assert n == 0 and len(s2) == len(s)


def test_aromatic_c_without_h_is_left_alone():
    c = [1.39 * np.array([np.cos(np.deg2rad(60 * k)), np.sin(np.deg2rad(60 * k)), 0.0])
         for k in range(6)]
    sp = ["C"] * 6 + ["H"] * 5
    xyz = c + [c[k] * (1 + 1.08 / 1.39) for k in range(1, 6)]
    _s2, n = F._repair_bare_sp3_ch2(box(sp, xyz))
    assert n == 0


def test_alkyne_c_is_left_alone():
    # H3C-C#C-CH3: the alkyne C have two partners each, linear, 1.46/1.20 A
    c1, c2, c3, c4 = (np.array([x, 0.0, 0.0]) for x in (0.0, 1.46, 2.66, 4.12))
    sp = ["C"] * 4 + ["H"] * 6
    xyz = [c1, c2, c3, c4] + methyl_h(c1, c2) + methyl_h(c4, c3)
    _s2, n = F._repair_bare_sp3_ch2(box(sp, xyz))
    assert n == 0


def vinyl_h(c, partner, away):
    """Two H on a terminal =CH2 carbon c (bonded to `partner`), in the plane
    of c, partner and `away`, 120 deg apart."""
    c, partner, away = (np.asarray(v, float) for v in (c, partner, away))
    u = (c - partner) / np.linalg.norm(c - partner)
    p = away - partner
    p = p - np.dot(p, u) * u
    p /= np.linalg.norm(p)
    return [c + 1.08 * (0.5 * u + s * 0.866 * p) for s in (1, -1)]


def test_ch2_between_sp2_carbons_is_left_alone():
    # CH2=CH-C-CH=CH2 (860: the carbon is drawn sp2, its partners are sp2):
    # the middle C's partners have three partners each, not four
    a = np.deg2rad(111.0) / 2
    c3 = np.zeros(3)
    c2 = 1.50 * np.array([-np.sin(a), -np.cos(a), 0.0])
    c4 = 1.50 * np.array([np.sin(a), -np.cos(a), 0.0])
    c1 = c2 + 1.34 * np.array([-1.0, 0.0, 0.0])
    c5 = c4 + 1.34 * np.array([1.0, 0.0, 0.0])
    h2 = c2 + 1.08 * np.array([0.0, -1.0, 0.0])
    h4 = c4 + 1.08 * np.array([0.0, -1.0, 0.0])
    sp = ["C"] * 5 + ["H"] * 6
    xyz = [c1, c2, c3, c4, c5, h2, h4] + vinyl_h(c1, c2, c3) + vinyl_h(c5, c4, c3)
    _s2, n = F._repair_bare_sp3_ch2(box(sp, xyz))
    assert n == 0, n


def test_wide_angle_bare_c_is_left_alone():
    a = np.deg2rad(124.0) / 2                   # 860's angle: sp2, not sp3
    c1 = 1.54 * np.array([-np.sin(a), -np.cos(a), 0.0])
    c3 = 1.54 * np.array([np.sin(a), -np.cos(a), 0.0])
    sp, xyz = ["C", "C", "C"], [c1, np.zeros(3), c3]
    for h in methyl_h(c1, xyz[1]) + methyl_h(c3, xyz[1]):
        sp.append("H")
        xyz.append(h)
    _s2, n = F._repair_bare_sp3_ch2(box(sp, xyz))
    assert n == 0, n


def test_repair_applies_only_when_it_makes_the_parent_neutral():
    sp, xyz = propane()
    s2, n = F.cof_repair_parent(box(sp, xyz))           # unreadable -> 0: repaired
    assert n == 1 and len(s2) == len(sp) + 2
    sp2, xyz2 = ammonium_methane()                      # unreadable -> +1: left as drawn
    s = box(sp + sp2, list(xyz) + [v + np.array([0.0, 8.0, 0.0]) for v in xyz2])
    s2, n = F.cof_repair_parent(s)
    assert n == 0 and len(s2) == len(s)


def ammonium_methane():
    sp = ["N", "H", "H", "H", "H", "C", "H", "H", "H", "H"]
    t = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
    xyz = [np.zeros(3)] + [1.01 * v for v in t]
    off = np.array([6.0, 0, 0])
    xyz += [off] + [off + 1.09 * v for v in t]
    return sp, xyz


def test_parent_formal_charge_values():
    sp, xyz = ammonium_methane()
    assert F.cof_parent_formal_charge(box(sp, xyz)) == 1
    sp, xyz = propane()
    assert F.cof_parent_formal_charge(box(sp, xyz)) is None      # bare C: implausible
    sp, xyz = propane(middle_h=2)
    assert F.cof_parent_formal_charge(box(sp, xyz)) == 0


def test_fragment_formal_charge_value():
    sp, xyz = ammonium_methane()
    assert F.cof_formal_charge_of(sp, xyz) == 1
    sp, xyz = propane()
    assert F.cof_formal_charge_of(sp, xyz) is None
    sp, xyz = propane(middle_h=2)
    assert F.cof_formal_charge_of(sp, xyz) == 0


def test_frame_verdict():
    even = ["C", "H", "H", "H", "H"]
    odd = ["C", "H", "H", "H"]
    assert _cof_frame_verdict(even, 0, 0) == (None, "")
    assert _cof_frame_verdict(even, 0, 6) == (None, "parent_charge=6")
    assert _cof_frame_verdict(even, 0, None) == (None, "parent_charge=unknown")
    assert _cof_frame_verdict(even, 2, 6) == (
        "quarantine=charge formal_charge=2 parent_charge=6", "parent_charge=6")
    assert _cof_frame_verdict(even, None, 0) == ("quarantine=charge formal_charge=unknown", "")
    r, _ = _cof_frame_verdict(odd, 0, 0)
    assert r.startswith("quarantine=odd_electron"), r


TESTS = [v for k, v in sorted(globals().items()) if k.startswith("test_")]

if __name__ == "__main__":
    sel = sys.argv[1] if len(sys.argv) > 1 else ""
    n_fail = 0
    for t in TESTS:
        if sel not in t.__name__:
            continue
        try:
            t()
            print("PASS", t.__name__)
        except Exception as e:  # noqa: BLE001
            n_fail += 1
            print("FAIL", t.__name__, type(e).__name__, e)
    ran = sum(1 for t in TESTS if sel in t.__name__)
    print(f"PASS {ran - n_fail} / FAIL {n_fail}")
    sys.exit(1 if n_fail else 0)
