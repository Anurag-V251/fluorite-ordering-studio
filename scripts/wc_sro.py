#!/usr/bin/env python3
"""
Warren-Cowley short-range order parameters from a POSCAR, sqsgenerator convention.

    alpha^i_ab = 1 - N^i_ab / (N * M^i * x_a * x_b)        (= 1 - P_i(b|a) / x_b)

    N^i_ab : number of a->b neighbour pairs in shell i (each bond counted from both ends)
    M^i    : mean coordination number of shell i
    N, x   : number of atoms and mole fractions of the analysed (sub)lattice

alpha < 0 : a-b pairs preferred (ordering)
alpha = 0 : random
alpha > 0 : a-b pairs avoided (clustering of like atoms)

Only numpy is needed.

Examples
    python wc_sro.py POSCAR --cutoff 5.7
    python wc_sro.py POSCAR --cutoff 5.7 --species Zr Y          # cation sublattice only
    python wc_sro.py POSCAR --shell-edges 0 4.8 5.6 --species Zr Y    # explicit shells (relaxed/MD cells)
    python wc_sro.py POSCAR --cutoff 5.7 --species Zr Y --csv sro.csv

Fluorite cation sublattice (fcc, 12 first neighbours): use --shell-edges 0 4.8 [5.6] for every phase.
In monoclinic cells the 12 first-shell cation-cation distances spread over ~3.3-4.6 Å (2nd shell from
~5.1 Å), and in oIII over ~3.4-3.9 Å, so a smaller first edge (e.g. 4.3 Å) or automatic shell detection
(--cutoff) splits or truncates shell 1. The script warns when shell 1 does not have 12 neighbours.
"""
import argparse
import sys

import numpy as np


# --------------------------------------------------------------------------- POSCAR
def read_poscar(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]

    scale = float(lines[1].split()[0])
    lat = np.array([[float(v) for v in lines[k].split()[:3]] for k in (2, 3, 4)])
    if scale < 0:                                    # negative scale = target volume
        lat *= (abs(scale) / abs(np.linalg.det(lat))) ** (1 / 3)
    else:
        lat *= scale

    tok = lines[5].split()
    if all(t.lstrip("-").isdigit() for t in tok):    # VASP 4: no species line
        sys.exit("POSCAR has no species line (VASP 4 format). Add the element names on line 6.")
    names = [t.split("/")[0].split("_")[0] for t in tok]
    counts = [int(v) for v in lines[6].split()]
    k = 7
    if lines[k].strip()[0] in "sS":                  # Selective dynamics
        k += 1
    cartesian = lines[k].strip()[0] in "cCkK"
    n = sum(counts)
    pos = np.array([[float(v) for v in lines[k + 1 + i].split()[:3]] for i in range(n)])
    frac = np.linalg.solve(lat.T, (pos * (scale if scale > 0 else 1)).T).T if cartesian else pos
    if cartesian and scale < 0:
        sys.exit("Negative scale factor with Cartesian coordinates is not supported.")
    symbols = np.repeat(names, counts)
    return lat, frac % 1.0, symbols


# --------------------------------------------------------------------------- neighbours
def pair_distances(lat, frac, cutoff):
    """All (i, j, d) with 0 < d <= cutoff, periodic images included, directed pairs."""
    vol = abs(np.linalg.det(lat))
    widths = np.array([vol / np.linalg.norm(np.cross(lat[(a + 1) % 3], lat[(a + 2) % 3]))
                       for a in range(3)])
    nrep = np.ceil(cutoff / widths).astype(int)
    shifts = np.array([[a, b, c] for a in range(-nrep[0], nrep[0] + 1)
                       for b in range(-nrep[1], nrep[1] + 1)
                       for c in range(-nrep[2], nrep[2] + 1)]) @ lat      # (S, 3)
    cart = frac @ lat
    I, J, D = [], [], []
    for i in range(len(cart)):
        vec = cart[None, :, :] + shifts[:, None, :] - cart[i]               # (S, N, 3)
        d = np.linalg.norm(vec, axis=2)
        s, j = np.nonzero((d > 1e-6) & (d <= cutoff))
        I.append(np.full(len(j), i)); J.append(j); D.append(d[s, j])
    return np.concatenate(I), np.concatenate(J), np.concatenate(D), widths


def find_shell_edges(d, cutoff, tol):
    """Split sorted distances wherever the gap exceeds tol; edges at gap midpoints."""
    ds = np.sort(d)
    gap = np.nonzero(np.diff(ds) > tol)[0]
    return np.concatenate(([0.0], (ds[gap] + ds[gap + 1]) / 2, [cutoff]))


# --------------------------------------------------------------------------- SRO
def warren_cowley(lat, frac, symbols, cutoff=None, species=None, shell_edges=None, tol=0.1):
    if species:
        missing = set(species) - set(symbols)
        if missing:
            sys.exit(f"Species not in structure: {', '.join(sorted(missing))}")
        keep = np.isin(symbols, species)
        frac, symbols = frac[keep], symbols[keep]

    if shell_edges is not None:
        cutoff = shell_edges[-1]
    i, j, d, widths = pair_distances(lat, frac, cutoff)
    if len(d) == 0:
        sys.exit("No neighbours within the cutoff.")
    edges = np.asarray(shell_edges if shell_edges is not None else find_shell_edges(d, cutoff, tol))

    names = list(dict.fromkeys(symbols))             # POSCAR order
    idx = np.array([names.index(s) for s in symbols])
    N = len(symbols)
    x = np.bincount(idx, minlength=len(names)) / N

    shell = np.searchsorted(edges, d, side="left") - 1
    results = []
    for k in range(len(edges) - 1):
        sel = shell == k
        if not sel.any():
            continue
        counts = np.zeros((len(names), len(names)))
        np.add.at(counts, (idx[i[sel]], idx[j[sel]]), 1)
        M = sel.sum() / N
        alpha = 1.0 - counts / (N * M * np.outer(x, x))
        results.append(dict(shell=len(results) + 1, rmin=d[sel].min(), rmax=d[sel].max(),
                            rmean=d[sel].mean(), M=M, alpha=alpha))
    return names, x, results, widths, edges


# --------------------------------------------------------------------------- CLI
def main():
    p = argparse.ArgumentParser(description="Warren-Cowley SRO parameters from a POSCAR.")
    p.add_argument("poscar")
    p.add_argument("--cutoff", type=float, help="max pair distance in Å (shells found automatically)")
    p.add_argument("--shell-edges", nargs="+", type=float,
                   help="explicit shell boundaries in Å, e.g. 0 4.8 5.6 for fluorite cations (overrides --cutoff)")
    p.add_argument("--species", nargs="+", help="analyse only these species (sublattice), e.g. Zr Y")
    p.add_argument("--tol", type=float, default=0.1,
                   help="min gap (Å) between distances to start a new shell (default 0.1)")
    p.add_argument("--csv", help="also write results to this CSV file")
    a = p.parse_args()
    if a.cutoff is None and a.shell_edges is None:
        p.error("give --cutoff or --shell-edges")

    lat, frac, symbols = read_poscar(a.poscar)
    names, x, res, widths, edges = warren_cowley(lat, frac, symbols, a.cutoff, a.species,
                                                 a.shell_edges, a.tol)
    cutoff = edges[-1]

    print(f"structure : {a.poscar}  ({len(symbols)} atoms)")
    n_sub = np.isin(symbols, a.species).sum() if a.species else len(symbols)
    print(f"sublattice: {n_sub} sites,  x: {'  '.join(f'{n}={c:.4f}' for n, c in zip(names, x))}")
    print(f"cutoff    : {cutoff:.3f} Å   shell edges: {' '.join(f'{e:.3f}' for e in edges)}")
    if cutoff > widths.min() / 2:
        print(f"WARNING   : cutoff > half the smallest cell width ({widths.min() / 2:.3f} Å); "
              "outer shells include periodic images of the same sites.")
    if res and abs(res[0]["M"] - 12) > 1e-6:
        print(f"WARNING   : shell 1 has {res[0]['M']:.2f} neighbours per site, not 12 as on the fcc cation "
              "sublattice of fluorite; alpha is then NOT the lattice shell-1 value written by "
              "enumerate_substitutions.py. Use --shell-edges 0 4.8 (monoclinic shell 1 reaches ~4.6 Å).")
    for r in res:
        if r["rmax"] - r["rmin"] > 0.3:
            print(f"WARNING   : shell {r['shell']} spans {r['rmin']:.2f}-{r['rmax']:.2f} Å; "
                  "consider --shell-edges or a smaller --tol.")
    print()

    rows = []
    print(f"{'shell':>5} {'r_mean':>7} {'r_range':>13} {'M':>6}   {'pair':<8} {'alpha':>9}")
    for r in res:
        first = True
        for ai in range(len(names)):
            for bi in range(ai, len(names)):
                al = r["alpha"][ai, bi]
                lead = (f"{r['shell']:>5} {r['rmean']:7.3f} {r['rmin']:6.3f}-{r['rmax']:<6.3f} {r['M']:6.2f}"
                        if first else " " * 34)
                print(f"{lead}   {names[ai] + '-' + names[bi]:<8} {al:+9.5f}")
                rows.append((r["shell"], r["rmean"], r["rmin"], r["rmax"], r["M"], names[ai], names[bi], al))
                first = False

    if a.csv:
        with open(a.csv, "w") as f:
            f.write("shell,r_mean,r_min,r_max,M,a,b,alpha\n")
            for row in rows:
                f.write(f"{row[0]},{row[1]:.5f},{row[2]:.5f},{row[3]:.5f},{row[4]:.4f},"
                        f"{row[5]},{row[6]},{row[7]:.6f}\n")
        print(f"\nwritten {a.csv}")


if __name__ == "__main__":
    main()
