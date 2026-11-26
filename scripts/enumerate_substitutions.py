#!/usr/bin/env python3
"""
Generate cation substitutions in a POSCAR (all of them, or N random ones) and label
each with its Warren-Cowley short-range-order parameter.

Example: ZrO2 cell with 12 Zr, substitute 4 Zr -> Ce
    -> folder  Ce_33.33/   (Ce fraction of the substituted sublattice, in %)
    -> each POSCAR's first (comment) line holds alpha(Ce-Zr), e.g.
         Zr8Ce4O24 | Ce 4/12 = 33.33% | WC_shell1(Ce-Zr)=-0.2500 | id=0001
    -> filename also carries it:  POSCAR_0001_WC-0.2500
    -> summary.csv lists every structure in the folder with its alpha and substituted sites

Random mode (--random N): draws N distinct configurations that are NOT already in the
folder (read from summary.csv), so you can run it again to add N more new ones.

Warren-Cowley (same formula as sqsgenerator / icet), on the substituted sublattice:
    alpha_DH = 1 - N_DH / (N * M * x_D * x_H)
    alpha < 0 : D-H neighbours preferred (mixing / ordering)
    alpha = 0 : random
    alpha > 0 : D-H neighbours avoided (clustering / segregation)

Usage
    python enumerate_substitutions.py POSCAR                                   # asks interactively
    python enumerate_substitutions.py POSCAR --host Zr --dopant Ce --n 4        # all combinations
    python enumerate_substitutions.py POSCAR --host Zr --dopant Ce --n 8 --random 50 --seed 1
    python enumerate_substitutions.py POSCAR --host Zr --dopant Ce --n 4 --unique   # needs pymatgen
    python enumerate_substitutions.py POSCAR --host Zr --dopant Ce --n 0            # end members (n = 0 or all):
                                                       # one structure, alpha written as nan

Only numpy is required (pymatgen only for --unique).
"""
import argparse
import csv
import os
import sys
from itertools import combinations
from math import comb

import numpy as np


# --------------------------------------------------------------------------- POSCAR I/O
def read_poscar(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    scale = float(lines[1].split()[0])
    lat = np.array([[float(v) for v in lines[k].split()[:3]] for k in (2, 3, 4)])
    lat = lat * ((abs(scale) / abs(np.linalg.det(lat))) ** (1 / 3) if scale < 0 else scale)

    tok = lines[5].split()
    if all(t.lstrip("-").isdigit() for t in tok):
        sys.exit("POSCAR has no species line (VASP 4 format). Add the element names on line 6.")
    names = [t.split("/")[0].split("_")[0] for t in tok]
    counts = [int(v) for v in lines[6].split()]
    k = 7
    if lines[k].strip()[0] in "sS":                   # Selective dynamics -> dropped
        k += 1
    cart = lines[k].strip()[0] in "cCkK"
    n = sum(counts)
    pos = np.array([[float(v) for v in lines[k + 1 + i].split()[:3]] for i in range(n)])
    if cart:
        pos = np.linalg.solve(lat.T, (pos * (scale if scale > 0 else 1.0)).T).T
    return lat, pos % 1.0, np.repeat(names, counts)


def write_poscar(path, comment, lat, frac, symbols, order=None):
    order = [o for o in (order or list(dict.fromkeys(symbols))) if np.any(symbols == o)]   # drop absent species
    with open(path, "w") as f:
        f.write(comment + "\n1.0\n")
        for v in lat:
            f.write(f"  {v[0]:20.12f} {v[1]:20.12f} {v[2]:20.12f}\n")
        f.write("  " + "  ".join(order) + "\n")
        f.write("  " + "  ".join(str(int(np.sum(symbols == s))) for s in order) + "\n")
        f.write("Direct\n")
        for s in order:
            for p in frac[symbols == s]:
                f.write(f"  {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}  {s}\n")


# --------------------------------------------------------------------------- neighbours
def pair_distances(lat, frac, cutoff):
    """All directed pairs (i, j, d) with 0 < d <= cutoff, periodic images included."""
    widths = 1.0 / np.linalg.norm(np.linalg.inv(lat), axis=0)       # perpendicular cell widths
    nrep = np.ceil(cutoff / widths).astype(int)
    grid = np.array([[a, b, c] for a in range(-nrep[0], nrep[0] + 1)
                     for b in range(-nrep[1], nrep[1] + 1)
                     for c in range(-nrep[2], nrep[2] + 1)]) @ lat
    cart = frac @ lat
    I, J, D = [], [], []
    for i in range(len(cart)):
        d = np.linalg.norm(cart[None, :, :] + grid[:, None, :] - cart[i], axis=2)
        s, j = np.nonzero((d > 1e-6) & (d <= cutoff))
        I.append(np.full(len(j), i)); J.append(j); D.append(d[s, j])
    return np.concatenate(I), np.concatenate(J), np.concatenate(D), widths


def shell_adjacency(lat, frac, n_shells, tol=0.1):
    """
    Adjacency matrices A[k] (n x n) of the sublattice for shells k = 0..n_shells-1.
    A[k][i, j] = number of times j (incl. periodic images) is a shell-k neighbour of i.
    """
    r = 1.6 * (abs(np.linalg.det(lat)) / len(frac)) ** (1 / 3)
    while True:
        i, j, d, widths = pair_distances(lat, frac, r)
        ds = np.sort(d)
        if len(ds):
            br = np.nonzero(np.diff(ds) > tol)[0]
            edges = np.concatenate(([0.0], (ds[br] + ds[br + 1]) / 2))
            if len(edges) > n_shells:                                 # n_shells complete shells
                break
        r *= 1.3
        if r > 50:
            sys.exit("Could not find enough coordination shells.")
    edges = edges[: n_shells + 1]
    shell = np.searchsorted(edges, d, side="left") - 1
    n = len(frac)
    A, info = [], []
    for k in range(n_shells):
        sel = shell == k
        M = np.zeros((n, n))
        np.add.at(M, (i[sel], j[sel]), 1)
        A.append(M)
        info.append((d[sel].mean(), sel.sum() / n))                   # (radius, coordination)
    return A, info, widths


def compute_alphas(combos, A, n_host, n_sub, chunk=50000):
    """alpha(dopant-host) per shell for an array of combinations (rows = sorted site indices)."""
    xD = n_sub / n_host
    xH = 1 - xD
    out = np.empty((len(combos), len(A)))
    for start in range(0, len(combos), chunk):
        cb = combos[start:start + chunk]
        S = np.zeros((len(cb), n_host))
        S[np.arange(len(cb))[:, None], cb] = 1.0                      # 1 = dopant
        for k, Ak in enumerate(A):
            M = Ak.sum() / n_host
            N_DH = np.einsum("ci,ij,cj->c", S, Ak, 1 - S)             # dopant -> host pairs
            if xD * xH == 0:                                          # pure host / pure dopant: undefined
                out[start:start + chunk, k] = np.nan
            else:
                out[start:start + chunk, k] = 1 - N_DH / (n_host * M * xD * xH)
    return np.round(out, 6) + 0.0                                     # kill -0.0


# --------------------------------------------------------------------------- selection
def random_combos(n_host, n_sub, k, exclude, rng, max_enum):
    """k distinct sorted combinations not in `exclude` (set of tuples)."""
    total = comb(n_host, n_sub)
    remaining = total - len(exclude)
    if k > remaining:
        print(f"only {remaining} new configurations left; generating all of them")
        k = remaining
    if k == 0:
        return np.empty((0, n_sub), dtype=np.int32)
    if total <= max_enum and k > 0.3 * remaining:                    # small space: enumerate + pick
        pool = [c for c in combinations(range(n_host), n_sub) if c not in exclude]
        pick = rng.choice(len(pool), size=k, replace=False)
        return np.array([pool[i] for i in sorted(pick)], dtype=np.int32)
    seen, out = set(exclude), []                                      # large space: rejection sampling
    while len(out) < k:
        c = tuple(sorted(rng.choice(n_host, size=n_sub, replace=False).tolist()))
        if c not in seen:
            seen.add(c); out.append(c)
    return np.array(out, dtype=np.int32)


def read_existing(folder, hidx):
    """Configurations (as host-sublattice indices) and max id already in folder/summary.csv."""
    path = os.path.join(folder, "summary.csv")
    if not os.path.exists(path):
        return set(), 0
    pos = {int(h): k for k, h in enumerate(hidx)}
    done, last = set(), 0
    with open(path) as f:
        for row in csv.DictReader(f):
            last = max(last, int(row["id"]))
            done.add(tuple(sorted(pos[int(s)] for s in row["substituted_atom_indices"].split())))
    return done, last


# --------------------------------------------------------------------------- main
def ask(prompt, cast=str, valid=None, default=None):
    while True:
        v = input(prompt).strip()
        if v == "" and default is not None:
            return default
        try:
            v = cast(v)
        except ValueError:
            print("  invalid input"); continue
        if valid is None or valid(v):
            return v
        print("  invalid choice")


def fmt_pct(x):
    return f"{100 * x:.2f}".rstrip("0").rstrip(".")


def main():
    p = argparse.ArgumentParser(description="Generate substitutions and label them with Warren-Cowley SRO.")
    p.add_argument("poscar")
    p.add_argument("--host", help="element to substitute (e.g. Zr)")
    p.add_argument("--dopant", help="new element (e.g. Ce)")
    p.add_argument("--n", type=int, help="number of host atoms to substitute")
    p.add_argument("--random", type=int, metavar="N",
                   help="generate N random distinct configurations not already in the folder")
    p.add_argument("--all", action="store_true", help="generate all combinations (no question asked)")
    p.add_argument("--seed", type=int, help="random seed (for reproducible --random)")
    p.add_argument("--shells", type=int, default=1, help="shells to report (default 1; filename uses shell 1)")
    p.add_argument("--tol", type=float, default=0.1, help="distance tolerance (Å) for grouping shells")
    p.add_argument("--unique", action="store_true",
                   help="all combinations, keep only symmetry-inequivalent ones (needs pymatgen)")
    p.add_argument("--outdir", default=".", help="where to create the Ce_xx folder")
    p.add_argument("--max", type=int, default=200000,
                   help="max combinations to enumerate for 'all' (default 200000)")
    a = p.parse_args()
    if a.random is not None and (a.all or a.unique):
        p.error("--random cannot be combined with --all or --unique")

    lat, frac, symbols = read_poscar(a.poscar)
    species = [str(s) for s in dict.fromkeys(symbols)]
    print(f"{a.poscar}: " + ", ".join(f"{s} = {np.sum(symbols == s)}" for s in species))

    host = a.host or ask(f"Which element to substitute {species}? ", str, lambda v: v in species)
    if host not in species:
        sys.exit(f"{host} not in structure")
    n_host = int(np.sum(symbols == host))
    dopant = a.dopant or ask("Substitute with which element? ", str, lambda v: v.isalpha())
    if dopant in species:
        sys.exit(f"{dopant} is already in the structure; choose a new element.")
    n_sub = a.n if a.n is not None else ask(f"How many of the {n_host} {host} to substitute? ",
                                            int, lambda v: 0 <= v <= n_host)
    if not 0 <= n_sub <= n_host:
        sys.exit(f"--n must be between 0 and {n_host}")
    if n_sub in (0, n_host):
        print(f"--n {n_sub}: end member (one structure); Warren-Cowley alpha is undefined there and written as nan")

    total = comb(n_host, n_sub)
    xD = n_sub / n_host
    folder = os.path.join(a.outdir, f"{dopant}_{fmt_pct(xD)}")
    hidx = np.nonzero(symbols == host)[0]
    done, last_id = read_existing(folder, hidx)
    print(f"{total:,} possible configurations (C({n_host},{n_sub}))"
          + (f", {len(done)} already in {folder}/" if done else ""))

    # ---- how many?
    n_random = a.random
    if n_random is None and not (a.all or a.unique):
        if a.host and a.dopant and a.n is not None and total <= a.max:
            n_random = None                                           # fully scripted call: all
        else:
            hint = f"[Enter = all {total - len(done)}]" if total <= a.max else "(too many for all)"
            n_random = ask(f"How many structures to generate? {hint} ", int, lambda v: v > 0,
                           default=-1 if total <= a.max else None)
            n_random = None if n_random == -1 else n_random
    if n_random is None and total > a.max:
        sys.exit(f"{total:,} combinations > --max {a.max:,}. Use --random N.")

    # ---- shells on the host sublattice
    A, info, widths = shell_adjacency(lat, frac[hidx], a.shells, a.tol)
    print(f"\nshells of the {host} sublattice:")
    for k, (r, Z) in enumerate(info, 1):
        print(f"  shell {k}: r = {r:.3f} Å, Z = {Z:.2f}")
    half = widths.min() / 2
    if info[-1][0] > half - 1e-6:
        print(f"WARNING: shell radius >= half the smallest cell width ({half:.3f} Å) -> "
              "neighbours include periodic images; alpha describes this periodic cell only.")

    # ---- choose configurations
    rng = np.random.default_rng(a.seed)
    if n_random is not None:
        combos = random_combos(n_host, n_sub, n_random, done, rng, a.max)
        print(f"\ndrew {len(combos)} new random configurations"
              + (f" (seed {a.seed})" if a.seed is not None else ""))
    else:
        combos = np.array([c for c in combinations(range(n_host), n_sub) if c not in done],
                          dtype=np.int32).reshape(-1, n_sub)
    if len(combos) == 0:
        print("nothing new to write."); return
    alphas = compute_alphas(combos, A, n_host, n_sub)

    # ---- optional symmetry reduction
    mult = np.ones(len(combos), dtype=int)
    keep = np.arange(len(combos))
    if a.unique:
        try:
            from pymatgen.core import Lattice, Structure
            from pymatgen.analysis.structure_matcher import StructureMatcher
        except ImportError:
            sys.exit("--unique needs pymatgen (pip install pymatgen)")
        sm = StructureMatcher(primitive_cell=False, scale=False, attempt_supercell=False)
        pl = Lattice(lat)
        keep, mult = [], []
        for key in np.unique(alphas, axis=0):                         # equivalent => same alpha
            grp = np.nonzero(np.all(alphas == key, axis=1))[0]
            structs = []
            for c in grp:
                sy = symbols.copy(); sy[hidx[combos[c]]] = dopant
                st = Structure(pl, list(sy), frac)
                st.properties = {"id": int(c)}
                structs.append(st)
            for g in sm.group_structures(structs):
                keep.append(g[0].properties["id"]); mult.append(len(g))
        keep, mult = np.array(keep), np.array(mult)
        print(f"\n{len(keep)} symmetry-inequivalent structures (from {len(combos)})")

    # ---- write: this batch numbered from most negative alpha to most positive
    srt = np.lexsort(alphas[keep].T[::-1])
    keep, mult = keep[srt], mult[srt]
    os.makedirs(folder, exist_ok=True)
    pair = f"{dopant}-{host}"
    order = species[:species.index(host) + 1] + [dopant] + species[species.index(host) + 1:]
    width = max(4, len(str(last_id + len(keep))))
    summary = os.path.join(folder, "summary.csv")
    new_file = not os.path.exists(summary)
    with open(summary, "a") as fs:
        if new_file:
            fs.write("id,file,multiplicity,"
                     + ",".join(f"alpha_shell{k + 1}_{pair}" for k in range(a.shells))
                     + ",substituted_atom_indices\n")
        for num, (c, m) in enumerate(zip(keep, mult), last_id + 1):
            sy = symbols.copy()
            sites = hidx[combos[c]]
            sy[sites] = dopant
            al = alphas[c]
            sro_txt = " ".join(f"WC_shell{k + 1}({pair})={al[k]:+.4f}" for k in range(a.shells))
            comp = "".join(f"{s}{np.sum(sy == s)}" for s in order if np.any(sy == s))
            comment = (f"{comp} | {dopant} {n_sub}/{n_host} = {fmt_pct(xD)}% | {sro_txt}"
                       + (f" | multiplicity={m}" if a.unique else "") + f" | id={num:0{width}d}")
            fname = f"POSCAR_{num:0{width}d}_WC{al[0]:+.4f}"
            write_poscar(os.path.join(folder, fname), comment, lat, frac, sy, order)
            fs.write(f"{num},{fname},{m}," + ",".join(f"{v:.6f}" for v in al)
                     + "," + " ".join(map(str, sites)) + "\n")

    # ---- overview
    print(f"\nwrote {len(keep)} POSCARs to {folder}/  (ids {last_id + 1}-{last_id + len(keep)}, "
          "summary.csv updated)")
    print(f"\nalpha_shell1({pair}) in this batch:")
    vals, cnt = np.unique(alphas[keep, 0], return_counts=True, equal_nan=True)
    for v, c in zip(vals, cnt):
        print(f"  {v:+.4f} : {c:6d}")
    print(f"\n0 = random, negative = {dopant}-{host} mixing/ordering, positive = clustering")


if __name__ == "__main__":
    main()
