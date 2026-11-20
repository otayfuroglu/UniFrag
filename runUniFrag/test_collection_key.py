#!/usr/bin/env python3
"""Cross-structure duplicate key of the MOF and COF collections (GitHub #14).

    python runUniFrag/test_collection_key.py [substring]

Two frames from different structures are duplicates only if they are the same
molecule: same heavy-atom bond graph (metal-N/O bonds up to 2.6 A) AND the same
H on the same atoms. The MOF collection used the heavy-atom formula alone, so
terephthalic acid (C8H6O4) was dropped as a "duplicate" of cyclohexane-1,4-
dicarboxylic acid (C8H12O4), and MAQBEO's linker (C18H18O8) of KEFBEE's
(C18H10O8). Conformers and lattice noise (same molecule, slightly different
coordinates) must still collapse.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import unifrag  # noqa: E402
from unifrag import _collection_key, _load_seen_keys_from_extxyz, _write_extxyz  # noqa: E402

RNG = np.random.default_rng(0)


def ring(n_h_per_c, r=1.39, z_pucker=0.0):
    """C6 ring with n H per carbon (1 = benzene, 2 = cyclohexane)."""
    sp, xyz = [], []
    for k in range(6):
        t = np.deg2rad(60 * k)
        c = np.array([r * np.cos(t), r * np.sin(t), z_pucker * (-1) ** k])
        sp.append("C")
        xyz.append(c)
        out = c / np.linalg.norm(c[:2]) if np.linalg.norm(c[:2]) else c
        out = np.array([out[0], out[1], 0.0])
        if n_h_per_c == 1:
            sp.append("H")
            xyz.append(c + 1.08 * out)
        else:
            for s in (1, -1):
                sp.append("H")
                xyz.append(c + 1.09 * (0.5 * out + s * 0.866 * np.array([0, 0, 1.0])))
    return sp, np.array(xyz)


def rotate(xyz, seed):
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return xyz @ q.T + np.array([3.0, -2.0, 7.0])


def test_benzene_and_cyclohexane_differ():
    b, c = ring(1), ring(2, r=1.53, z_pucker=0.25)
    assert _collection_key(*b) != _collection_key(*c)


def test_same_molecule_moved_and_jittered_is_a_duplicate():
    sp, xyz = ring(1)
    moved = rotate(xyz, 1) + RNG.normal(scale=0.02, size=xyz.shape)
    assert _collection_key(sp, xyz) == _collection_key(sp, moved)


def test_h_on_different_atoms_differ():
    # CH3-CH2-OH vs CH3-O-CH3 analogue on one skeleton: same H count,
    # H moved from O to C (C-C-O chain with 6 H)
    sp = ["C", "C", "O"] + ["H"] * 6
    c1, c2, o = np.array([0, 0, 0.0]), np.array([1.52, 0, 0]), np.array([2.0, 1.35, 0])
    h_on_c1 = [c1 + v for v in ([-0.5, 0.9, 0], [-0.5, -0.9, 0], [-0.4, 0, 1.0])]
    ethanol = h_on_c1 + [c2 + [0.3, -0.9, 0.5], c2 + [0.3, -0.9, -0.5], o + [0.9, 0.3, 0]]
    other = h_on_c1 + [c2 + [0.3, -0.9, 0.5], c2 + [0.3, -0.9, -0.5], c2 + [0.2, 0.3, -1.0]]
    k1 = _collection_key(sp, np.array([c1, c2, o] + ethanol))
    k2 = _collection_key(sp, np.array([c1, c2, o] + other))
    assert k1 != k2


def test_metal_coordination_counts():
    # Zn with 4 O at 1.98 A vs 3 O at 1.98 A and one at 3.0 A (not bonded)
    d = [np.array(v, float) for v in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
    d = [v / np.linalg.norm(v) for v in d]
    sp = ["Zn"] + ["O"] * 4 + ["H"] * 8
    def build(far):
        xyz = [np.zeros(3)]
        os_ = [v * (3.0 if (far and k == 3) else 1.98) for k, v in enumerate(d)]
        xyz += os_
        for v, p in zip(d, os_):
            perp = np.cross(v, [0.3, 0.5, 0.8])
            perp /= np.linalg.norm(perp)
            xyz += [p + 0.96 * (0.33 * v + 0.94 * perp), p + 0.96 * (0.33 * v - 0.94 * perp)]
        return np.array(xyz)
    assert _collection_key(sp, build(False)) != _collection_key(sp, build(True))


def test_resume_loader_uses_the_collection_key():
    sp, xyz = ring(1)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "fragments_collection.extxyz")
        _write_extxyz(p, sp, xyz, "XFragMof")
        keys = _load_seen_keys_from_extxyz(p)
    assert keys == {_collection_key(sp, xyz)}
    assert _collection_key(*ring(2, r=1.53, z_pucker=0.25)) not in keys


def test_mof_writer_uses_the_collection_key():
    # the MOF and COF writers live inside main(); no heavy-formula-only key may
    # be left in them (the within-structure dedupe keeps the heavy key)
    src = open(unifrag.__file__).read()
    mof = src[src.index("def _flush_mof_result"):src.index("def _flush_cof_result")]
    assert "_chemical_identity_key" not in mof, "MOF writer still uses the heavy-formula key"
    assert mof.count("_collection_key(") >= 5


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
