#!/usr/bin/env python3
"""MOF capping never puts H on a complete sp atom (GitHub #12).

    python runUniFrag/test_sp_caps.py [substring]

Synthetic molecules, no CIF. An alkyne carbon (two partners, linear) and a
nitrile nitrogen are complete; a bent carbon with two partners and no H is an
aromatic edge that lost its H and keeps its cap, however short its bonds are
drawn (disordered vinylene/imine CIFs draw C=C at 1.11-1.25 A).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentation_oop import MOFFragmenter  # noqa: E402

F = MOFFragmenter()


def methyl(at, away):
    """Three H on carbon `at`, pointing away from direction `away`."""
    away = np.asarray(away, float) / np.linalg.norm(away)
    perp = np.cross(away, [0, 0, 1.0])
    if np.linalg.norm(perp) < 1e-6:
        perp = np.cross(away, [0, 1.0, 0])
    perp /= np.linalg.norm(perp)
    perp2 = np.cross(away, perp)
    out = []
    for k in range(3):
        t = 2 * np.pi * k / 3
        d = -away * 0.36 + (np.cos(t) * perp + np.sin(t) * perp2) * 1.03
        out.append(np.asarray(at, float) + d / np.linalg.norm(d) * 1.09)
    return out


def molecule(heavy):
    """Heavy atoms [(symbol, xyz, cap_with_methyl_H_away_from)] -> species, coords."""
    sp, co = [], []
    for s, p, away in heavy:
        sp.append(s)
        co.append(np.asarray(p, float))
    for s, p, away in heavy:
        if away is not None:
            for h in methyl(p, away):
                sp.append("H")
                co.append(h)
    return sp, co


def dimethylacetylene():
    # H3C-C#C-CH3 along x, the two sp carbons carry no H
    return molecule([("C", [-1.46, 0, 0], [1, 0, 0]), ("C", [0, 0, 0], None),
                     ("C", [1.20, 0, 0], None), ("C", [2.66, 0, 0], [-1, 0, 0])])


def bent_edge(d=1.39, angle=120.0):
    # central C with two methyl partners at `angle`, no H on the centre
    a = np.deg2rad(angle / 2)
    p1 = [-d * np.cos(a), d * np.sin(a), 0]
    p2 = [-d * np.cos(a), -d * np.sin(a), 0]
    return molecule([("C", [0, 0, 0], None), ("C", p1, [-x for x in p1]),
                     ("C", p2, [-x for x in p2])])


def acetonitrile():
    return molecule([("C", [-1.46, 0, 0], [1, 0, 0]), ("C", [0, 0, 0], None),
                     ("N", [1.16, 0, 0], None)])


def bent_amine_n():
    # N with one C partner whose own two bonds are bent: not a nitrile
    sp, co = molecule([("C", [0, 0, 0], None), ("C", [-0.70, 1.20, 0], [0.70, -1.20, 0]),
                       ("N", [1.40, 0, 0], None)])
    sp.append("H")
    co.append(np.array([-0.55, -0.95, 0.0]))     # the centre C keeps its own H
    return sp, co


def caps_on(sp, co, idx):
    sp, co = list(sp), [np.array(c) for c in co]
    flags = [False] * len(sp)
    F._cap_open_oxygens(sp, co, flags)
    return sum(1 for j, s in enumerate(sp) if s == "H" and flags[j]
               and np.linalg.norm(co[j] - co[idx]) < 1.25)


def test_alkyne_carbons_are_complete():
    sp, co = dimethylacetylene()
    assert F._is_complete_sp_atom(1, sp, co) and F._is_complete_sp_atom(2, sp, co)


def test_alkyne_carbons_get_no_cap_in_linker_capping():
    sp, co = dimethylacetylene()
    assert caps_on(sp, co, 1) == 0 and caps_on(sp, co, 2) == 0


def test_alkyne_with_long_drawn_triple_is_still_complete():
    # a model that draws C#C at 1.30 A is linear all the same
    sp, co = molecule([("C", [-1.46, 0, 0], [1, 0, 0]), ("C", [0, 0, 0], None),
                       ("C", [1.30, 0, 0], None), ("C", [2.76, 0, 0], [-1, 0, 0])])
    assert F._is_complete_sp_atom(1, sp, co)


def test_nearly_linear_165_counts_as_sp():
    sp, co = bent_edge(d=1.20, angle=166.0)
    assert F._is_complete_sp_atom(0, sp, co)


def test_bent_edge_carbon_is_not_complete_and_keeps_its_cap():
    sp, co = bent_edge()
    assert not F._is_complete_sp_atom(0, sp, co)
    assert caps_on(sp, co, 0) == 1


def test_short_drawn_bent_carbon_keeps_its_cap():
    # disordered vinylene drawn at 1.22 A, 123 deg: sp2, not sp
    sp, co = bent_edge(d=1.22, angle=123.0)
    assert not F._is_complete_sp_atom(0, sp, co)
    assert caps_on(sp, co, 0) == 1


def test_nitrile_n_is_complete_and_uncapped():
    sp, co = acetonitrile()
    assert F._is_complete_sp_atom(2, sp, co)
    assert F._is_complete_sp_atom(1, sp, co)
    assert caps_on(sp, co, 2) == 0 and caps_on(sp, co, 1) == 0


def test_n_on_a_bent_carbon_is_not_a_nitrile():
    sp, co = bent_amine_n()
    assert not F._is_complete_sp_atom(2, sp, co)


def cut_dicyanamide_arm():
    # one arm of N(C#N)2- whose central N lost its other C in the cut:
    # Ncentral-C 1.30 A, C#N 1.15 A, linear (FALROB, JIQCUK Min frames)
    return molecule([("N", [-1.30, 0, 0], None), ("C", [0, 0, 0], None),
                     ("N", [1.15, 0, 0], None)])


def test_cut_dicyanamide_central_n_is_not_a_nitrile_and_keeps_its_cap():
    sp, co = cut_dicyanamide_arm()
    assert not F._is_complete_sp_atom(0, sp, co)
    assert caps_on(sp, co, 0) == 1


def test_cut_dicyanamide_terminal_n_and_carbon_are_complete():
    sp, co = cut_dicyanamide_arm()
    assert F._is_complete_sp_atom(2, sp, co) and F._is_complete_sp_atom(1, sp, co)
    assert caps_on(sp, co, 2) == 0 and caps_on(sp, co, 1) == 0


def test_metal_partner_does_not_count():
    # Zn<-N#C-CH3: the metal is not a sigma partner of the N
    sp, co = acetonitrile()
    sp.append("Zn")
    co.append(np.array([3.16, 0, 0]))
    assert F._is_complete_sp_atom(2, sp, co)


def test_skip_hook_matches_predicate():
    sp, co = dimethylacetylene()
    assert F._skip_saturated_h_cap(1, sp, co)
    sp, co = bent_edge()
    assert not F._skip_saturated_h_cap(0, sp, co)


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
