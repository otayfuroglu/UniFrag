#!/usr/bin/env python3
"""COF minimized-fragment route (GitHub #13).

    python runUniFrag/test_cof_min_route.py [substring]

The Min route - first-ring trim, node + one linker, or no Min - is decided on
how many HEAVY atoms the minimized fragment removes from the normal one
(>= MIN_HEAVY_ATOM_REDUCTION_FOR_MINIMIZE_COF = 15). Cap H do not count: a
trim also adds a cap H per new cut, so counting H let cap-only changes switch
the route (1043) and judged some trims by their H (1079: 8 heavy + 16 H).
"""
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fragmentation_oop import (  # noqa: E402
    MIN_HEAVY_ATOM_REDUCTION_FOR_MINIMIZE_COF, _choose_cof_min)


def frag(heavy, h):
    return SimpleNamespace(species=["C"] * heavy + ["H"] * h)


class Once:
    """one-linker builder that records whether it was called."""

    def __init__(self, res):
        self.res, self.calls = res, 0

    def __call__(self):
        self.calls += 1
        return self.res


def test_threshold_is_15_heavy_atoms():
    assert MIN_HEAVY_ATOM_REDUCTION_FOR_MINIMIZE_COF == 15


def test_trim_that_adds_cap_h_is_kept():
    # 245: first-ring removes 16 heavy atoms but adds 8 cap H (all-atom count: 8)
    n, ring = frag(108, 18), frag(92, 26)
    alt = Once(frag(44, 14))
    chosen, red, route = _choose_cof_min(n, ring, alt)
    assert chosen is ring and red == 16 and route == "first_ring"
    assert alt.calls == 0


def test_h_rich_trim_falls_back_to_one_linker():
    # 1079: first-ring removes 8 heavy + 16 H (all-atom count: 24)
    n, ring, link = frag(48, 54), frag(40, 38), frag(26, 28)
    chosen, red, route = _choose_cof_min(n, ring, Once(link))
    assert chosen is link and red == 22 and route == "one_linker"


def test_exactly_15_is_enough():
    n, ring = frag(60, 30), frag(45, 40)
    chosen, red, route = _choose_cof_min(n, ring, Once(None))
    assert chosen is ring and red == 15 and route == "first_ring"


def test_neither_reaches_threshold_gives_no_min():
    n, ring, link = frag(60, 30), frag(50, 20), frag(48, 20)
    chosen, red, route = _choose_cof_min(n, ring, Once(link))
    assert chosen is None and red == 12 and route is None


def test_one_linker_only_when_it_removes_more():
    n, ring, link = frag(60, 30), frag(50, 20), frag(52, 10)
    chosen, red, route = _choose_cof_min(n, ring, Once(link))
    assert chosen is None and red == 10 and route is None


def test_cap_h_change_cannot_switch_route():
    # same skeletons, the Min carries 3 more or fewer cap H
    n = frag(70, 40)
    for dh in (-3, -1, 0, 1, 3):
        chosen, red, route = _choose_cof_min(n, frag(55, 30 + dh), Once(frag(40, 20)))
        assert route == "first_ring" and red == 15, (dh, route, red)


def test_no_first_ring_candidate_gives_no_min():
    alt = Once(frag(10, 5))
    chosen, red, route = _choose_cof_min(frag(60, 30), None, alt)
    assert chosen is None and route is None and alt.calls == 0


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
