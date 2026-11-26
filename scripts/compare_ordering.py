#!/usr/bin/env python3
"""
compare_ordering.py A_POSCAR B_POSCAR  - do two fluorite-derived supercells have the same cation ordering?

Both structures (any phase: cubic, tetragonal, monoclinic, oIII, ...) are put on the ideal fluorite lattice
(every cation tied to its nearest ideal fcc site, as in map_fluorite.py) and the species on those lattice points
are compared.  Phase and distortion therefore do not matter, only which species sits on which lattice point.

Verdicts
  IDENTICAL        same species on the same lattice points up to a global translation (the frame of the two
                   files agrees, e.g. both made from one input by map_fluorite.py)
  EQUIVALENT       identical only after a rotation/reflection of the lattice axes (relabelled axes)
  DIFFERENT        no lattice symmetry operation maps one ordering onto the other
Also printed: composition, pair counts per neighbour shell, and (if the atom order of the two files is the same)
whether atom i of A sits on the same lattice point as atom i of B.

  python3 compare_ordering.py A B [--n NX NY NZ]
"""
import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import map_fluorite as M


def lattice_decoration(fn, n_arg):
    cell, frac, sp = M.read_poscar(fn)
    n = tuple(n_arg) if n_arg else M.detect_n(cell, frac, sp)
    nv = np.array(n)
    iscat = sp != 'O'
    if int(iscat.sum()) != 4 * int(np.prod(n)) or int((~iscat).sum()) != 8 * int(np.prod(n)):
        M.die('%s: atom counts do not match a %dx%dx%d repeat of the 12-atom cell' % (fn, n[0], n[1], n[2]))
    u, species = M.decoration_on_lattice(cell, frac, sp, n)
    site_dist = M.match(frac[iscat], M.ideal_sites(M.find_frame(frac[iscat], cell, n)[0], n)[0], cell)[1].max()
    if site_dist > M.MAX_CAT_DISP:
        M.die('%s is not a small distortion of fluorite (cation %.2f A from ideal sites)' % (fn, site_dist))
    return dict(fn=fn, n=nv, u=u, sp=species, cell=cell, cat_frac=frac[iscat])


def pair_counts(d):
    """neighbour pairs per fcc shell (1: 12 neighbours, 2: 6, 3: 24) by species, on the periodic lattice"""
    w = d['u'] * d['n']
    nn = d['n'].astype(float)
    diff = w[:, None, :] - w[None, :, :]
    diff = diff - nn * np.round(diff / nn)
    r2 = (diff ** 2).sum(-1)
    out = {}
    N = len(w)
    for shell, target in ((1, 0.5), (2, 1.0), (3, 1.5)):
        ii, jj = np.where((np.abs(r2 - target) < 1e-6) & (np.arange(N)[:, None] < np.arange(N)[None, :]))
        c = {}
        for i, j in zip(ii, jj):
            k = '-'.join(sorted((d['sp'][i], d['sp'][j])))
            c[k] = c.get(k, 0) + 1
        out[shell] = c
    return out


def transformed(d, g, n_target):
    """lattice points of d after the lattice symmetry g (signed axis permutation), in the supercell of n_target"""
    w = (d['u'] * d['n']) @ g
    n_new = np.abs(d['n'] @ g).round().astype(int)
    if not np.array_equal(n_new, n_target):
        return None
    return (w / n_target) % 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('a')
    ap.add_argument('b')
    ap.add_argument('--n', type=int, nargs=3, metavar=('NX', 'NY', 'NZ'), help='supercell repeats (default: auto, same for both)')
    args = ap.parse_args()
    A = lattice_decoration(args.a, args.n)
    B = lattice_decoration(args.b, args.n)
    ca = {str(s): int((A['sp'] == s).sum()) for s in sorted(set(A['sp']))}
    cb = {str(s): int((B['sp'] == s).sum()) for s in sorted(set(B['sp']))}
    print('A: %s  supercell %s  cations %s' % (args.a, 'x'.join(map(str, A['n'])), ca))
    print('B: %s  supercell %s  cations %s' % (args.b, 'x'.join(map(str, B['n'])), cb))
    if ca != cb:
        print('\nVERDICT: DIFFERENT - the compositions differ, so the orderings cannot be the same.')
        return 1

    verdict, op = None, None
    if np.array_equal(A['n'], B['n']) and M.same_decoration(A['u'], A['sp'], B['u'], B['sp']):
        verdict = 'IDENTICAL'
    else:
        for g in M.signed_perms():
            ub = transformed(B, g, A['n'])
            if ub is not None and M.same_decoration(A['u'], A['sp'], ub, B['sp']):
                verdict, op = 'EQUIVALENT', g
                break
    pa, pb = pair_counts(A), pair_counts(B)
    print('\npair counts per fcc neighbour shell (A | B):')
    for sh in (1, 2, 3):
        keys = sorted(set(pa[sh]) | set(pb[sh]))
        print('  shell %d: %s' % (sh, '   '.join('%s %d|%d' % (k, pa[sh].get(k, 0), pb[sh].get(k, 0)) for k in keys)))

    same_order = len(A['sp']) == len(B['sp']) and np.array_equal(A['sp'], B['sp'])
    if same_order and np.array_equal(A['n'], B['n']):
        d = A['cat_frac'] - B['cat_frac']
        d -= np.round(d)
        d -= d.mean(0)
        shift = np.linalg.norm(d @ A['cell'], axis=1)
        du = A['u'] - B['u']                        # lattice points of atom i in A and in B (same frame)
        du -= np.round(du)
        e = du - du[0]                              # all equal up to ONE global translation?
        e -= np.round(e)
        same_sites = bool(np.abs(e).max() < 1e-6)
        print('\natom-by-atom: atom i of A and atom i of B sit on the same lattice point for %s'
              % ('all %d cations' % len(du) if same_sites else 'only some cations'))
        print('              (positions differ by up to %.2f A after removing one global shift)' % shift.max())

    if verdict == 'IDENTICAL':
        print('\nVERDICT: IDENTICAL - the same species sit on the same lattice points (up to a global translation).')
        return 0
    if verdict == 'EQUIVALENT':
        print('\nVERDICT: EQUIVALENT - the same ordering after relabelling the lattice axes of B with\n%s' % g.astype(int))
        return 0
    print('\nVERDICT: DIFFERENT - no lattice translation, rotation or reflection maps one ordering onto the other.')
    same_pairs = all(pa[s] == pb[s] for s in (1, 2, 3))
    print('         (pair counts in the first three shells %s)' % ('are equal, so they are only distinguishable by longer-range correlations' if same_pairs else 'differ as well'))
    return 1


if __name__ == '__main__':
    sys.exit(main())
