#!/usr/bin/env python3
"""
Warren-Cowley SRO parameters from a POSCAR using icet's built-in ShortRangeOrderObserver.
Just pass how many coordination shells you want; the cutoff is found automatically (ASE).

    alpha_AB^m = 1 - N_AB^m / <N_AB^m>_random       (per shell m, per sublattice)

Requires:  pip install icet       (pulls in ase, numpy, spglib)

Examples
    # first 2 cation shells, Zr/Y on the cation sublattice, O fixed
    python icet_sro.py POSCAR --shells 2 --sublattice Zr Y

    # O vacancies stored as 'X' atoms in the POSCAR
    python icet_sro.py POSCAR --shells 2 --sublattice Zr Y --sublattice O X

    # relaxed / MD structure: map onto the ideal lattice first (shells taken from the ideal one)
    python icet_sro.py CONTCAR --shells 3 --sublattice Zr Y --reference POSCAR_ideal

    # set the cutoff yourself instead of --shells
    python icet_sro.py POSCAR --cutoff 5.7 --sublattice Zr Y --csv sro.csv
"""
import argparse
import re
import sys
from itertools import combinations, combinations_with_replacement

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
from icet import ClusterSpace
from mchammer.observers import ShortRangeOrderObserver


# --------------------------------------------------------------------------- shells (ASE)
def find_shells(atoms, species, n_shells, tol=0.1, r_max=40.0):
    """
    Coordination shells of the sublattice made of `species`.
    Returns list of (r_mean, r_min, r_max, coordination) for at least n_shells + 1 shells.
    """
    idx = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s in species]
    if not idx:
        sys.exit(f"No atoms of {species} in the structure.")
    sub = atoms[idx]
    r = 1.6 * (sub.get_volume() / len(sub)) ** (1 / 3)        # ~ a bit beyond the 1st shell
    while True:
        d = np.sort(neighbor_list("d", sub, r))
        if len(d):
            breaks = np.nonzero(np.diff(d) > tol)[0] + 1
            groups = np.split(d, breaks)
            # last group may be cut by r; only trust groups well inside r
            complete = [g for g in groups if g.max() < r - tol]
            if len(complete) > n_shells:
                return [(g.mean(), g.min(), g.max(), len(g) / len(sub)) for g in complete]
        r *= 1.3
        if r > r_max:
            sys.exit(f"Could not find {n_shells + 1} shells below {r_max} Å.")


def half_width(atoms):
    """Half of the smallest perpendicular cell width (minimum-image limit)."""
    return 0.5 / np.linalg.norm(atoms.cell.reciprocal(), axis=1).max()


# --------------------------------------------------------------------------- icet setup
def allowed_species(atoms, sublattices):
    """Per-site allowed species: a site whose species is in a sublattice gets that whole list."""
    allowed = []
    for sym in atoms.get_chemical_symbols():
        for sub in sublattices:
            if sym in sub:
                allowed.append(sorted(sub))
                break
        else:
            allowed.append([sym])          # not listed -> fixed, ignored by the observer
    return allowed


def main():
    p = argparse.ArgumentParser(description="Warren-Cowley SRO from a POSCAR with icet.")
    p.add_argument("poscar", help="structure on the ideal lattice (or use --reference)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--shells", type=int, help="number of coordination shells (cutoff found automatically)")
    g.add_argument("--cutoff", type=float, help="max pair distance in Å (manual)")
    p.add_argument("--sublattice", nargs="+", action="append", required=True, metavar="EL",
                   help="species sharing one sublattice, e.g. --sublattice Zr Y (repeatable); "
                        "shells are determined on the FIRST one")
    p.add_argument("--tol", type=float, default=0.1,
                   help="distances closer than this (Å) belong to the same shell (default 0.1)")
    p.add_argument("--reference", help="ideal-lattice structure to map a relaxed POSCAR onto")
    p.add_argument("--no-like-pairs", action="store_true", help="only report unlike pairs (A-B)")
    p.add_argument("--csv", help="also write results to this CSV file")
    a = p.parse_args()

    atoms = read(a.poscar, format="vasp")
    if a.reference:
        from icet.tools import map_structure_to_reference
        ref = read(a.reference, format="vasp")
        atoms, info = map_structure_to_reference(atoms, ref)
        print(f"mapped onto reference: max displacement {info['drmax']:.3f} Å, "
              f"avg {info['dravg']:.3f} Å")

    # ---- composition check
    syms = atoms.get_chemical_symbols()
    print(f"species in file: " + ", ".join(f"{s}={syms.count(s)}" for s in dict.fromkeys(syms)))
    for sub in a.sublattice:
        missing = [s for s in sub if s not in syms]
        if missing:
            print(f"WARNING  : {', '.join(missing)} not in the structure -> alpha for "
                  f"{'/'.join(missing)} pairs will be nan. Check the species line of the POSCAR.")

    # ---- shells and cutoff
    shells = find_shells(atoms, a.sublattice[0], (a.shells or 1), a.tol)
    if a.shells:
        n = a.shells
        cutoff = 0.5 * (shells[n - 1][2] + shells[n][1])       # midpoint of the gap after shell n
    else:
        cutoff = a.cutoff
        n = sum(s[2] <= cutoff for s in shells)

    print(f"structure: {a.poscar}  ({len(atoms)} atoms)")
    print(f"shells of the {'/'.join(a.sublattice[0])} sublattice:")
    print(f"  {'shell':>5} {'r (Å)':>8} {'range':>15} {'Z':>6}")
    for k, (rm, r0, r1, Z) in enumerate(shells[: n + 1], 1):
        mark = "" if k <= n else "   <- excluded"
        print(f"  {k:>5} {rm:8.3f} {r0:7.3f}-{r1:<7.3f} {Z:6.2f}{mark}")
    print(f"cutoff   : {cutoff:.3f} Å  ({n} shell{'s' if n != 1 else ''})")
    if cutoff > half_width(atoms):
        print(f"WARNING  : cutoff > half the smallest cell width ({half_width(atoms):.3f} Å); "
              "outer shells contain periodic images of the same sites.")

    # ---- icet
    cs = ClusterSpace(atoms, cutoffs=[cutoff], chemical_symbols=allowed_species(atoms, a.sublattice))
    pairs = []
    for sub in a.sublattice:
        sub = sorted(sub)
        combo = combinations(sub, 2) if a.no_like_pairs else combinations_with_replacement(sub, 2)
        pairs += [tuple(pr) for pr in combo]
    obs = ShortRangeOrderObserver(cs, atoms, radius=cutoff, pairs=pairs)
    sro = obs.get_observable(atoms)

    rows = []
    for key, val in sro.items():
        m = re.match(r"sro_(.+)_(.+)_(\d+)$", key)
        rows.append((int(m.group(3)), m.group(1), m.group(2), val) if m else (0, key, "", val))
    rows.sort()

    print(f"\n{'shell':>5}  {'pair':<10} {'alpha':>10}")
    for shell, A, B, val in rows:
        print(f"{shell:>5}  {A + '-' + B:<10} {val:+10.5f}")

    if a.csv:
        with open(a.csv, "w") as f:
            f.write("shell,a,b,alpha\n")
            for shell, A, B, val in rows:
                f.write(f"{shell},{A},{B},{val:.6f}\n")
        print(f"\nwritten {a.csv}")


if __name__ == "__main__":
    main()
