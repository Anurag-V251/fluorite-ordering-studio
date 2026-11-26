#!/usr/bin/env python3
"""
Regression test for the shell-1 Warren-Cowley alpha (no NEP, numpy only, ~10 s).

  python tests/test_wc_shell1.py

1. enumerate_substitutions.py labels (title + filename) == wc_sro.py == an independent count on the
   integer fcc lattice, for several supercell shapes (incl. repeats of 1) and compositions.
2. Exact sum rule: the mean alpha over ALL configurations of a cell is -1/(N-1) (not 0).
3. Known orderings (2x2x2, 50 %): L1_0 (001) layers -1/3, (111) planes 0, Ce slab +1/3.
4. alpha is unchanged when the same decoration is put on the m-ZrO2 cation pattern (cell and cation
   positions of refs/POSCAR_m, embedded below; first-shell cation distances 3.33-4.55 A):
   wc_sro.py --shell-edges 0 4.8 gives 12 neighbours and the parent alpha; 0 4.3 misses one
   neighbour (M = 11, wrong alpha) and must warn.
"""
import itertools
import os
import re
import subprocess
import sys
import tempfile
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, S)
import enumerate_substitutions as ES   # noqa: E402
import wc_sro                          # noqa: E402

FCC = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
TET = np.array([[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)])
# refs/POSCAR_m (m-ZrO2, 12-atom cell, axes along the cube directions): lattice rows and the 4 cations
M_LAT = np.array([[0.0, 5.1333199981058977, 0.0084648081920438],
                  [5.2103974283849110, 0.0, 0.0],
                  [0.0, -0.8789926184819942, -5.2304577856144707]])
M_CAT = np.array([[0.2771200240411886, 0.0431790825065548, 0.2096288825968884],
                  [0.2771200240411886, 0.4568209174934452, 0.7096288825968884],
                  [0.7228799759588114, 0.5431790825065548, 0.2903711174031116],
                  [0.7228799759588114, 0.9568209174934452, 0.7903711174031116]])
V1 = np.array([v for v in itertools.product((-1, 0, 1), repeat=3) if sorted(map(abs, v)) == [0, 1, 1]])
FAIL = []


def check(ok, msg):
    print(("  ok    " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def build(n, a=5.10):
    cells = np.array(list(itertools.product(*[range(k) for k in n])))
    cat = ((cells[:, None] + FCC[None]) / n).reshape(-1, 3)
    ox = ((cells[:, None] + TET[None]) / n).reshape(-1, 3)
    return np.diag(np.asarray(n) * a), np.vstack([cat, ox]), np.array(["Zr"] * len(cat) + ["O"] * len(ox))


def alpha_indep(frac_cat, is_ce, n):
    """1 - P(Zr | Ce)/x_Zr from integer lattice coordinates (units of a/2, periodic modulo 2n)."""
    n2 = 2 * np.asarray(n)
    g = np.rint(frac_cat * n2).astype(int) % n2
    lut = {tuple(k): c for k, c in zip(g, is_ce)}
    slots = [lut[tuple((k + v) % n2)] for k, c in zip(g, is_ce) if c for v in V1]
    return 1 - (1 - np.mean(slots)) / (1 - np.mean(is_ce))


def wc(path, edges):
    lat, frac, sym = wc_sro.read_poscar(path)
    names, _, res, _, _ = wc_sro.warren_cowley(lat, frac, sym, species=["Zr", "Ce"], shell_edges=edges)
    return res[0]["alpha"][names.index("Ce"), names.index("Zr")], res[0]["M"]


def main():
    tmp = tempfile.mkdtemp(prefix="wc_test_")
    print("== 1. labels vs wc_sro.py vs independent count")
    for n in [(1, 1, 1), (1, 1, 2), (2, 2, 2), (1, 2, 3), (2, 2, 3), (3, 3, 3)]:
        lat, frac, sp = build(n)
        d = os.path.join(tmp, "x".join(map(str, n)))
        os.makedirs(d)
        ES.write_poscar(os.path.join(d, "POSCAR"), "test", lat, frac, sp)
        N = int((sp == "Zr").sum())
        for nd in sorted({1, N // 4, N // 2} - {0, N}):
            r = subprocess.run([sys.executable, os.path.join(S, "enumerate_substitutions.py"), os.path.join(d, "POSCAR"),
                                "--host", "Zr", "--dopant", "Ce", "--n", str(nd), "--random", "15", "--seed", "1",
                                "--outdir", d], capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
            folder = [f for f in os.listdir(d) if f.startswith("Ce_")]
            worst = 0.0
            files = [os.path.join(d, fo, f) for fo in folder for f in os.listdir(os.path.join(d, fo))
                     if f.startswith("POSCAR_") and f" {nd}/{N} " in open(os.path.join(d, fo, f)).readline()]
            for f in files:
                title = float(re.search(r"WC_shell1\(Ce-Zr\)=([+-][\d.]+)", open(f).readline()).group(1))
                fname = float(re.search(r"_WC([+-][\d.]+)$", f).group(1))
                lat2, fr2, sy2 = wc_sro.read_poscar(f)
                c = np.isin(sy2, ["Zr", "Ce"])
                ind = alpha_indep(fr2[c], sy2[c] == "Ce", n)
                w, M = wc(f, [0, 4.8])
                worst = max(worst, abs(title - ind), abs(fname - ind), abs(w - ind), abs(M - 12))
            check(worst < 1e-4, f"{'x'.join(map(str, n)):6s} Ce {nd:3d}/{N:<3d}: {len(files)} structures, max deviation {worst:.1e}")

    print("== 2. mean over all configurations = -1/(N-1)")
    lat, frac, sp = build((2, 2, 2))
    hidx = np.nonzero(sp == "Zr")[0]
    A, _, _ = ES.shell_adjacency(lat, frac[hidx], 1)
    for nd in (1, 2, 3):
        cb = np.array(list(itertools.combinations(range(32), nd)), dtype=np.int32)
        m = ES.compute_alphas(cb, A, 32, nd)[:, 0].mean()
        check(len(cb) == comb(32, nd) and abs(m + 1 / 31) < 1e-6, f"Ce {nd}/32: mean {m:+.6f}, expected {-1 / 31:+.6f}")

    print("== 3. known orderings")
    fc = frac[hidx]
    ijk = np.rint(fc * 4).astype(int)
    for label, mask, expect in [("L1_0 (001) layers", ijk[:, 2] % 2 == 1, -1 / 3),
                                ("(111) planes", ijk.sum(1) % 4 == 2, 0.0),
                                ("Ce slab z < 1/2", fc[:, 2] < 0.49, 1 / 3)]:
        e = ES.compute_alphas(np.nonzero(mask)[0][None].astype(np.int32), A, 32, int(mask.sum()))[0, 0]
        check(abs(e - expect) < 1e-6 and abs(alpha_indep(fc, mask, (2, 2, 2)) - expect) < 1e-9,
              f"{label:18s}: {e:+.4f} (expected {expect:+.4f})")

    print("== 4. monoclinic cell: same decoration, wc_sro.py shell edges")
    f0 = sorted(f for f in os.listdir(os.path.join(tmp, "2x2x2", "Ce_25")) if f.startswith("POSCAR_"))[3]
    f0 = os.path.join(tmp, "2x2x2", "Ce_25", f0)
    lat, frac, sym = wc_sro.read_poscar(f0)
    title = float(re.search(r"WC_shell1\(Ce-Zr\)=([+-][\d.]+)", open(f0).readline()).group(1))
    # displacement of each m cation from its ideal fcc site, laid on the same lattice points of the 2x2x2 cell
    sh = (M_CAT - M_CAT[0]) % 1
    k_of = [int(np.argmin(np.linalg.norm((sh[i] - FCC + .5) % 1 - .5, axis=1))) for i in range(4)]
    u = {k: (sh[i] - FCC[k] + .5) % 1 - .5 for i, k in enumerate(k_of)}
    assert sorted(k_of) == [0, 1, 2, 3]
    cat = np.isin(sym, ["Zr", "Ce"])
    new = frac.copy()
    for i in np.nonzero(cat)[0]:
        x = frac[i] * 2
        cell = np.floor(x + 1e-6)
        k = int(np.argmin(np.linalg.norm((x - cell - FCC + .5) % 1 - .5, axis=1)))
        new[i] = (cell + FCC[k] + u[k]) / 2
    fm = os.path.join(tmp, "POSCAR_m_pattern")
    ES.write_poscar(fm, "m-ZrO2 cation pattern on " + os.path.basename(f0), M_LAT * 2, new % 1, sym, ["Zr", "Ce", "O"])
    w48, M48 = wc(fm, [0, 4.8])
    w43, M43 = wc(fm, [0, 4.3])
    out = subprocess.run([sys.executable, os.path.join(S, "wc_sro.py"), fm, "--shell-edges", "0", "4.3",
                          "--species", "Zr", "Ce"], capture_output=True, text=True).stdout
    check(abs(w48 - title) < 1e-4 and abs(M48 - 12) < 1e-9, f"0-4.8 A: alpha {w48:+.4f} = label {title:+.4f}, M = {M48:g}")
    check(abs(M43 - 11) < 1e-9 and "not 12" in out, f"0-4.3 A: M = {M43:g}, alpha {w43:+.4f}, warning printed")

    print("\nALL PASSED" if not FAIL else f"\n{len(FAIL)} FAILED")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
