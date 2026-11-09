#!/usr/bin/env python3
"""COF: a cut imine-type N gets its H count from the parent, not from a
bond length (GitHub #11).

    python runUniFrag/test_cof_cut_n.py [substring]

Synthetic parents in a 30 A box. Each is drawn with the misleading lengths
the CoRE-COF models use: the C=N long (1.38-1.45 A, reads single) and the
aryl C-N short (1.33-1.35 A, reads double). An aryl C-N is single, so an N
left on its ring carbon needs NH2; an aldimine N left on its CH keeps C=N and
needs one H - unless the parent lost an H (amide, enamine) and the N-C to the
open-chain carbon is single too.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentation_oop import COFFragmenter  # noqa: E402

BOX = 30.0


def ring(center, start_angle, r=1.39):
    """Six ring atoms of a hexagon, the first at `start_angle` (deg)."""
    c = np.asarray(center, float)
    return [c + r * np.array([np.cos(np.deg2rad(start_angle + 60 * k)),
                              np.sin(np.deg2rad(start_angle + 60 * k)), 0.0]) for k in range(6)]


def ring_h(atoms, center, skip):
    c = np.asarray(center, float)
    return [a + 1.08 * (a - c) / np.linalg.norm(a - c) for k, a in enumerate(atoms) if k not in skip]


class Mol:
    def __init__(self):
        self.sp, self.xyz, self.tag = [], [], {}

    def add(self, s, p, tag=None):
        self.sp.append(s)
        self.xyz.append(np.asarray(p, float) + BOX / 2)
        if tag:
            self.tag[tag] = len(self.sp) - 1
        return len(self.sp) - 1


def benzylideneaniline(d_aryl_n=1.33, d_c_n=1.40, n_h=False, ch_partner="ring"):
    """Ar-N=CH-Ar' drawn with misleading lengths. `n_h` adds an H on N (not
    used: parents here have none). `ch_partner` other than 'ring' is set by
    the amide/enamine builders below."""
    m = Mol()
    ca = [0.0, 0.0, 0.0]
    ra = ring(ca, 0.0)                       # aniline ring, ra[0] bonded to N
    for k, p in enumerate(ra):
        m.add("C", p, tag=f"a{k}")
    for p in ring_h(ra, ca, {0}):
        m.add("H", p)
    n = ra[0] + np.array([d_aryl_n, 0, 0])
    m.add("N", n, tag="N")
    u = np.array([np.cos(np.deg2rad(60)), np.sin(np.deg2rad(60)), 0.0])
    ch = n + d_c_n * u
    m.add("C", ch, tag="CH")
    m.add("H", ch + 1.08 * np.array([np.cos(np.deg2rad(120)), np.sin(np.deg2rad(120)), 0.0]))
    cb0 = ch + 1.46 * np.array([1.0, 0, 0])
    cb = cb0 + np.array([1.39, 0, 0])
    rb = ring(cb, 180.0)                     # rb[0] == cb0
    for k, p in enumerate(rb):
        m.add("C", p, tag=f"b{k}")
    for p in ring_h(rb, cb, {0}):
        m.add("H", p)
    return m


def amide_like():
    """Ar-N-C(=O)-Ar' with the amide H omitted (as some models draw it)."""
    m = Mol()
    ca = [0.0, 0.0, 0.0]
    ra = ring(ca, 0.0)
    for k, p in enumerate(ra):
        m.add("C", p, tag=f"a{k}")
    for p in ring_h(ra, ca, {0}):
        m.add("H", p)
    n = ra[0] + np.array([1.35, 0, 0])
    m.add("N", n, tag="N")
    u = np.array([np.cos(np.deg2rad(60)), np.sin(np.deg2rad(60)), 0.0])
    cc = n + 1.35 * u
    m.add("C", cc, tag="CO")
    m.add("O", cc + 1.23 * np.array([np.cos(np.deg2rad(120)), np.sin(np.deg2rad(120)), 0.0]))
    cb0 = cc + 1.48 * np.array([1.0, 0, 0])
    cb = cb0 + np.array([1.39, 0, 0])
    rb = ring(cb, 180.0)
    for k, p in enumerate(rb):
        m.add("C", p, tag=f"b{k}")
    for p in ring_h(rb, cb, {0}):
        m.add("H", p)
    return m


def fragmenter(m):
    from pymatgen.core import Lattice, Structure
    F = COFFragmenter()
    s = Structure(Lattice.cubic(BOX), m.sp, m.xyz, coords_are_cartesian=True)
    F._parent_struct = s
    F.structure = s
    F._length_scale_value = 1.0
    return F


def cut(m, keep):
    """Fragment: the atoms in `keep` plus their parent H (parent coordinates)."""
    keep = set(keep)
    for i, s in enumerate(m.sp):
        if s == "H":
            k = min((j for j in range(len(m.sp)) if m.sp[j] != "H"),
                    key=lambda j: np.linalg.norm(m.xyz[j] - m.xyz[i]))
            if k in keep:
                keep.add(i)
    idx = sorted(keep)
    return [m.sp[i] for i in idx], [m.xyz[i].copy() for i in idx], idx


def n_h_after_capping(m, keep):
    F = fragmenter(m)
    sp, co, idx = cut(m, keep)
    flags = [False] * len(sp)
    F._cap_severed_double_bond_sites(sp, co, flags)
    n = idx.index(m.tag["N"])
    return sum(1 for j, s in enumerate(sp) if s == "H" and np.linalg.norm(co[j] - co[n]) < 1.15), F, sp, co, n


def aniline_side(m):
    return [m.tag[f"a{k}"] for k in range(6)] + [m.tag["N"]]


def benzylidene_side(m, key="CH"):
    return [m.tag[f"b{k}"] for k in range(6)] + [m.tag[key], m.tag["N"]]


# --- the two #11 motifs ------------------------------------------------------

def test_aryl_n_cut_from_its_imine_c_gets_nh2():
    # census group A: retained aryl C-N drawn 1.33 A (old rule: double -> 1 H)
    m = benzylideneaniline(d_aryl_n=1.33, d_c_n=1.35)
    nh, *_ = n_h_after_capping(m, aniline_side(m))
    assert nh == 2, nh


def test_imine_n_cut_from_its_aryl_keeps_c_eq_n_and_gets_one_h():
    # census group B: retained C=N drawn 1.45 A (old rule: single -> 2 H)
    m = benzylideneaniline(d_aryl_n=1.40, d_c_n=1.45)
    nh, *_ = n_h_after_capping(m, benzylidene_side(m))
    assert nh == 1, nh


def test_imine_n_with_c_eq_n_drawn_as_triple_gets_one_h():
    # census group D: retained C=N drawn 1.17 A (old rule: triple -> 0 H)
    m = benzylideneaniline(d_aryl_n=1.40, d_c_n=1.17)
    nh, *_ = n_h_after_capping(m, benzylidene_side(m))
    assert nh == 1, nh


def test_well_drawn_imine_unchanged_both_sides():
    m = benzylideneaniline(d_aryl_n=1.42, d_c_n=1.28)
    assert n_h_after_capping(m, aniline_side(m))[0] == 2
    assert n_h_after_capping(m, benzylidene_side(m))[0] == 1


# --- parents that lost an H: the N-C to the open-chain carbon is single -----

def test_h_deficient_amide_n_on_carbonyl_side_gets_nh2():
    m = amide_like()
    nh, *_ = n_h_after_capping(m, benzylidene_side(m, key="CO") + [m.tag["CO"] + 1])
    assert nh == 2, nh


def test_h_deficient_amide_n_on_aryl_side_gets_nh2():
    m = amide_like()
    nh, *_ = n_h_after_capping(m, aniline_side(m))
    assert nh == 2, nh


# --- the parent rule itself ----------------------------------------------------

def test_parent_rule_reports_orders():
    m = benzylideneaniline(d_aryl_n=1.33, d_c_n=1.45)
    F = fragmenter(m)
    sp, co, idx = cut(m, aniline_side(m))
    n, a0 = idx.index(m.tag["N"]), idx.index(m.tag["a0"])
    assert F._parent_cut_n_retained_order(n, a0, sp, co) == 1
    sp, co, idx = cut(m, benzylidene_side(m))
    n, ch = idx.index(m.tag["N"]), idx.index(m.tag["CH"])
    assert F._parent_cut_n_retained_order(n, ch, sp, co) == 2


def test_parent_rule_silent_without_parent():
    m = benzylideneaniline()
    F = COFFragmenter()
    sp, co, idx = cut(m, aniline_side(m))
    assert F._parent_cut_n_retained_order(idx.index(m.tag["N"]), idx.index(m.tag["a0"]), sp, co) is None


def test_parent_rule_silent_for_n_n_linkage():
    # hydrazone/azine C=N-N: a different motif (635, 637 read it as azo)
    m = benzylideneaniline(d_aryl_n=1.40, d_c_n=1.29)
    m.sp[m.tag["a0"]] = "N"
    F = fragmenter(m)
    sp, co, idx = cut(m, benzylidene_side(m))
    assert F._parent_cut_n_retained_order(idx.index(m.tag["N"]), idx.index(m.tag["CH"]), sp, co) is None


def test_parent_rule_silent_for_n_with_parent_h():
    # an N that already has an H in the parent is not an imine-type N
    m = benzylideneaniline()
    m.add("H", m.xyz[m.tag["N"]] - BOX / 2 + np.array([0.3, -1.0, 0.0]))
    F = fragmenter(m)
    sp, co, idx = cut(m, aniline_side(m))
    assert F._parent_cut_n_retained_order(idx.index(m.tag["N"]), idx.index(m.tag["a0"]), sp, co) is None


# --- fragment-level choice at cut imine N when the fragment is charged -------

def capped_fragment(m, keep, n_caps):
    """Fragment with the cut N given exactly `n_caps` H, placed trigonal."""
    F = fragmenter(m)
    sp, co, idx = cut(m, keep)
    n = idx.index(m.tag["N"])
    c = next(j for j in range(len(sp)) if sp[j] != "H" and j != n
             and np.linalg.norm(co[j] - co[n]) < 1.6)
    u = (co[n] - co[c]) / np.linalg.norm(co[n] - co[c])
    w = np.array([-u[1], u[0], 0.0])
    caps = []
    for k, ang in enumerate((60.0, -60.0)[:n_caps]):
        d = np.cos(np.deg2rad(ang)) * u + np.sin(np.deg2rad(ang)) * w
        sp.append("H")
        co.append(co[n] + 1.01 * d / np.linalg.norm(d))
        caps.append(len(sp) - 1)
    return F, sp, co, caps, n


def h_count(sp, co, n):
    return sum(1 for j, s in enumerate(sp) if s == "H" and np.linalg.norm(co[j] - co[n]) < 1.15)


def test_balance_adds_h_to_undercapped_aryl_n():
    m = benzylideneaniline()
    F, sp, co, caps, n = capped_fragment(m, aniline_side(m), 1)      # Ar-NH-  (-1)
    sp2, co2, caps2, balanced = F._balance_formal_charge(sp, co, caps, label="t")
    assert balanced is False
    assert h_count(sp2, co2, n) == 2 and len(caps2) == 2
    assert F._fragment_formal_charge(sp2, co2) == (0, 0)


def test_balance_removes_h_from_overcapped_imine_n():
    m = benzylideneaniline()
    F, sp, co, caps, n = capped_fragment(m, benzylidene_side(m), 2)  # CH=NH2+ (+1)
    sp2, co2, caps2, _ = F._balance_formal_charge(sp, co, caps, label="t")
    assert h_count(sp2, co2, n) == 1 and len(caps2) == 1
    assert F._fragment_formal_charge(sp2, co2) == (0, 0)


def test_balance_leaves_neutral_fragment_alone():
    m = benzylideneaniline()
    F, sp, co, caps, n = capped_fragment(m, aniline_side(m), 2)      # Ar-NH2, neutral
    sp2, co2, caps2, _ = F._balance_formal_charge(list(sp), list(co), list(caps), label="t")
    assert sp2 == sp and len(caps2) == 2 and all(np.allclose(a, b) for a, b in zip(co, co2))


def test_balance_is_all_or_nothing():
    # a charge the cut N cannot cancel (a ring C that lost its H): the cut N
    # must not be half-moved, which would only make the electron count odd
    # Ar-NH- (-1, fixable) plus a phenolate O- on the ring (-1, not a cut N):
    # one move reaches -1, never 0, so nothing may move.
    m = benzylideneaniline()
    F, sp, co, caps, n = capped_fragment(m, aniline_side(m), 1)
    ring_h = [j for j in range(len(sp)) if sp[j] == "H" and j not in caps]
    h = ring_h[1]
    c = min((j for j in range(len(sp)) if sp[j] == "C"), key=lambda j: np.linalg.norm(co[j] - co[h]))
    u = (co[h] - co[c]) / np.linalg.norm(co[h] - co[c])
    sp[h], co[h] = "O", co[c] + 1.36 * u
    assert F._fragment_formal_charge(sp, co)[0] == -2
    sp2, co2, caps2, _ = F._balance_formal_charge(list(sp), list(co), list(caps), label="t")
    assert sp2 == sp and caps2 == caps and h_count(sp2, co2, n) == 1


def test_balance_never_touches_n_without_parent():
    m = benzylideneaniline()
    F, sp, co, caps, n = capped_fragment(m, aniline_side(m), 1)
    F._parent_struct = None
    sp2, co2, caps2, _ = F._balance_formal_charge(list(sp), list(co), list(caps), label="t")
    assert sp2 == sp and len(caps2) == 1


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
