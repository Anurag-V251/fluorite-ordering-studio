#!/usr/bin/env python3
"""
map_fluorite.py - move a fluorite-derived ZrO2/CeO2-type supercell into another phase
                  while keeping the cation ordering (Ce/Zr/... decoration) exactly the same.

Idea
----
Every phase of this family (cubic fluorite, tetragonal, monoclinic, orthorhombic Pca2_1, ...)
is a small distortion of the same fluorite lattice: the cations sit on a distorted fcc lattice
and the oxygens on its tetrahedral holes.  The ideal fluorite lattice is therefore a common
"bond network" for all phases.  The script

  1. finds that ideal fluorite frame for your input structure (nearest-site matching; every
     atom is tied to the lattice point it came from),
  2. reads the displacement pattern of the target phase from a 12-atom reference cell
     (4 cations + 8 O, axes along the fluorite directions), measured from ITS ideal sites,
  3. puts the target pattern on the same lattice points.  Every cation keeps its lattice
     point and therefore its species, so the cation ordering is identical by construction.

Among the 192 symmetry-equivalent ways of laying the target pattern on the lattice (which cubic axis
becomes c_t, which direction is polar, ...) the script takes the one that needs the smallest distortion
from the input, as a displacive transformation does: smallest rms atomic shift plus smallest cell strain.
Choices with exactly the same distortion (for example every cubic input) are reported as such.
Atom i of the input is atom i of the output (file order is kept).

Usage
-----
  python3 map_fluorite.py INPUT_POSCAR --template REFERENCE_POSCAR [--n NX NY NZ] [--out FILE]
                          [--identity] [--keep-volume]

  INPUT_POSCAR       supercell, nx*ny*nz repeats of the 12-atom conventional cell (any fluorite-derived
                     phase, any species on the cation sites, Direct/Cartesian, selective dynamics ok)
  REFERENCE_POSCAR   12-atom cell of the TARGET phase (cubic, tetragonal pseudo-cubic setting, oIII, m, ...)
  --n                supercell repeats; detected automatically from the cell if omitted
  --identity         do not search orientations; lay the template pattern on the input axes as they are
  --keep-volume      rescale the output isotropically to the input's volume per formula unit
                     (default: keep the reference's lattice parameters; relax afterwards in either case)

Exit status is non-zero if anything is inconsistent (the script verifies its own result:
cation ordering, handedness, collisions along the distortion path).
"""
import argparse
import itertools
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment

S_FCC = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
Q_O = np.array([[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)])
MAX_CAT_DISP = 0.8      # A, largest cation displacement from an ideal fluorite site that we accept
MAX_O_DISP = 1.4        # A, same for oxygen


def die(msg):
    sys.stderr.write('ERROR: %s\n' % msg)
    sys.exit(2)


# ----------------------------------------------------------------------------- I/O
def read_poscar(fn):
    try:
        L = open(fn).read().split('\n')
    except OSError as e:
        die('cannot read %s (%s)' % (fn, e))
    try:
        scale = float(L[1].split()[0])
        cell = abs(scale) * np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])
        names = L[5].split()
        cnt = [int(x) for x in L[6].split()]
        n = sum(cnt)
        idx = 7
        if L[idx].strip()[0] in 'Ss':                 # Selective dynamics
            idx += 1
        mode = L[idx].strip()[0]
        xyz = np.array([[float(x) for x in L[idx + 1 + i].split()[:3]] for i in range(n)])
    except (ValueError, IndexError):
        die('%s is not a VASP5 POSCAR (needs a species-name line, counts, Direct/Cartesian coordinates)' % fn)
    frac = (abs(scale) * xyz) @ np.linalg.inv(cell) if mode in 'CcKk' else xyz
    sp = np.array(sum([[s] * c for s, c in zip(names, cnt)], []))
    return cell, frac, sp


def write_poscar(fn, cell, species, frac, title):
    order = list(dict.fromkeys(species))
    idx = np.concatenate([np.where(species == s)[0] for s in order])
    out = [title, '1.0'] + ['  %.12f  %.12f  %.12f' % tuple(r) for r in cell]
    out += [' '.join(order), ' '.join(str(int((species == s).sum())) for s in order), 'Direct']
    out += ['  %.12f  %.12f  %.12f' % tuple(x % 1.0) for x in frac[idx]]
    open(fn, 'w').write('\n'.join(out) + '\n')
    return idx


# ----------------------------------------------------------------------------- ideal fluorite frame
def pdiff(f1, f2):
    d = f1[:, None, :] - f2[None, :, :]
    return d - np.round(d)


def ideal_sites(shift, n):
    """Ideal fluorite sites of an n-repeat supercell (supercell fractional coords).
    `shift` (conventional-cell units) places the cation fcc lattice; oxygens sit at (1/4|3/4)^3."""
    nv = np.array(n)
    cells = np.array(list(itertools.product(*[range(m) for m in n])))
    cat = ((cells[:, None, :] + S_FCC[None, :, :] + shift) / nv).reshape(-1, 3) % 1.0
    ox = ((cells[:, None, :] + Q_O[None, :, :] + shift) / nv).reshape(-1, 3) % 1.0
    return cat, ox


def match(frac, ideal, cell):
    """one-to-one nearest-site assignment (Hungarian); returns site index per atom and distances in A"""
    dist = np.linalg.norm(pdiff(frac, ideal) @ cell, axis=2)
    r, c = linear_sum_assignment(dist ** 2)
    site = np.empty(len(frac), int)
    site[r] = c
    return site, dist[r, c]


def find_frame(cat_frac, cell, n):
    """Position of the ideal fluorite frame (conventional-cell units): coarse grid search, then refined so
    that the cation displacement field has zero mean (a global translation of the input is irrelevant).
    Returns (shift, rms cation displacement from the ideal sites)."""
    nv = np.array(n)
    half = np.arange(0, .5, 1 / 8)                     # shifts differing by an fcc translation are equivalent
    best = (np.inf, None)
    for sx, sy, sz in itertools.product(half, half, np.arange(0, 1, 1 / 8)):
        sh = np.array([sx, sy, sz])
        cat, _ = ideal_sites(sh, n)
        _, dist = match(cat_frac, cat, cell)
        tot = (dist ** 2).sum()
        if tot < best[0] - 1e-12:
            best = (tot, sh)
    sh = best[1]
    for _ in range(3):
        cat, _ = ideal_sites(sh, n)
        site, dist = match(cat_frac, cat, cell)
        d = cat_frac - cat[site]
        d -= np.round(d)
        sh = sh + d.mean(0) * nv
    cat, _ = ideal_sites(sh, n)
    _, dist = match(cat_frac, cat, cell)
    return sh, float(np.sqrt((dist ** 2).mean()))


def best_shift(cat_frac, cell, n):                      # kept for the test-suite
    return find_frame(cat_frac, cell, n)[0]


def codes(q):
    """integer code of site positions q (conventional-cell units, relative to the frame): quarters mod 1"""
    k = np.round(np.asarray(q) * 4).astype(int) % 4
    return k[..., 0] * 16 + k[..., 1] * 4 + k[..., 2]


def signed_perms():
    mats = []
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            g = np.zeros((3, 3))
            for i in range(3):
                g[i, perm[i]] = sg[i]
            mats.append(g)
    return mats                                          # identity is first


def detect_n(cell, frac, sp):
    ncat = int((sp != 'O').sum())
    nO = int((sp == 'O').sum())
    if nO != 2 * ncat or ncat % 4:
        die('input has %d cations and %d O: expected a fluorite-derived cell with 4 cations + 8 O per '
            'conventional cell (O = 2 x cations, cations a multiple of 4)' % (ncat, nO))
    M = ncat // 4
    lens = np.linalg.norm(cell, axis=1)
    cands = []
    for nx in range(1, M + 1):
        if M % nx:
            continue
        for ny in range(1, M // nx + 1):
            if (M // nx) % ny:
                continue
            nz = M // (nx * ny)
            conv = lens / np.array([nx, ny, nz])
            if np.all((conv > 4.3) & (conv < 6.2)):
                cands.append((nx, ny, nz))
    if not cands:
        die('cannot find a supercell nx*ny*nz = %d whose conventional-cell lengths are 4.3-6.2 A '
            '(cell lengths %s); pass --n explicitly' % (M, np.round(lens, 2)))
    if len(cands) == 1:
        return cands[0]
    cat = frac[sp != 'O']
    scored = sorted((find_frame(cat, cell, c)[1], c) for c in cands)
    if len(scored) > 1 and scored[1][0] < 3 * scored[0][0] + 0.05:
        die('supercell repeats are ambiguous (%s fit about equally); pass --n explicitly' %
            ', '.join(str(c) for _, c in scored[:3]))
    return scored[0][1]


# ----------------------------------------------------------------------------- template (target phase)
def template_fields(tfile):
    """displacement field of the 12-atom reference on its own ideal sites (conventional-cell units),
    as lookup arrays indexed by site code"""
    cell, frac, sp = read_poscar(tfile)
    iscat = sp != 'O'
    if iscat.sum() != 4 or (~iscat).sum() != 8:
        die('template must be the 12-atom conventional cell (4 cations + 8 O) with axes along the fluorite '
            'directions (pseudo-cubic setting for tetragonal); %s has %d cations and %d O' %
            (tfile, iscat.sum(), (~iscat).sum()))
    sh, _ = find_frame(frac[iscat], cell, (1, 1, 1))
    ic, io = ideal_sites(sh, (1, 1, 1))
    sC, dC = match(frac[iscat], ic, cell)
    sO, dO = match(frac[~iscat], io, cell)
    if dC.max() > MAX_CAT_DISP or dO.max() > MAX_O_DISP:
        die('template does not look like a distorted fluorite cell in the right setting (max cation displacement '
            '%.2f A, max O displacement %.2f A from the ideal sites). Use the 12-atom cell with axes along '
            'the cubic fluorite directions.' % (dC.max(), dO.max()))
    Tc = np.zeros((64, 3))
    To = np.zeros((64, 3))
    used = set()
    for arr, fr, ids, site, tag in ((Tc, frac[iscat], ic, sC, 'cation'), (To, frac[~iscat], io, sO, 'O')):
        for a, s in enumerate(site):
            d = fr[a] - ids[s]
            d -= np.round(d)
            k = int(codes(ids[s] - sh))
            if (tag, k) in used:
                die('template %s sites are not distinct on the fluorite lattice' % tag)
            used.add((tag, k))
            arr[k] = d
    return cell, Tc, To


def strain_measure(G_in, G_out):
    """rotation-free lattice strain between two cells given their metrics G = A A^T.
    Principal stretches s_k are the singular values of the deformation; the common volume change is divided out
    (it is identical for every correspondence), S = |s_k/(prod s)^(1/3) - 1|."""
    lam = np.linalg.eigvals(np.linalg.solve(G_in, G_out)).real
    s = np.sqrt(np.clip(lam, 1e-12, None))
    s_dev = s / np.prod(s) ** (1 / 3)
    return float(np.sqrt(((s_dev - 1) ** 2).sum())), np.sort(s_dev)


def distortion(cell_in, f_in, A_out, f_out, nv):
    """distortion needed to take the input into the output while every atom keeps its lattice point:
    non-affine atomic shifts (A, global translation removed) and lattice strain S.
    Score J (A^2) = mean square shift + (S * a)^2, a = conventional cell edge, i.e. the strain counted as the
    shift of an atom one cell away."""
    G_in = cell_in @ cell_in.T
    G_out = A_out @ A_out.T
    S, stretches = strain_measure(G_in, G_out)
    d = f_out - f_in
    d -= np.round(d)
    d -= d.mean(0)
    u2 = np.einsum('ij,jk,ik->i', d, 0.5 * (G_in + G_out), d)
    a_ref = float(np.mean(np.linalg.norm(cell_in, axis=1) / nv))
    return dict(drms=float(np.sqrt(u2.mean())), dmax=float(np.sqrt(u2.max())), S=S, stretches=stretches,
                J=float(u2.mean() + (S * a_ref) ** 2))


def rank_correspondences(cell, f_in, tcell, Tc, To, idealc, idealo, cat_site, o_site, qC, qO, nv, identity_only):
    """try all 192 correspondences (signed axis permutation g x fcc origin t) and keep the minimum-distortion ones.
    Returns (variants within numerical noise of the minimum, best score, next-best score)."""
    ops = [(g, t) for g in signed_perms() for t in S_FCC]          # identity first
    if identity_only:
        ops = ops[:1]
    scored = []
    for g, t in ops:
        v = build_variant(g, t, tcell, Tc, To, idealc, idealo, cat_site, o_site, qC, qO, nv)
        v.update(distortion(cell, f_in, v['Anew'] * nv[:, None], v['pos'], nv))
        scored.append(v)
    Jmin = min(v['J'] for v in scored)
    tol = 5e-4 + 5e-3 * Jmin                                        # relaxed references carry ~1e-3 A noise
    tied = [v for v in scored if v['J'] <= Jmin + tol]
    rest = [v['J'] for v in scored if v['J'] > Jmin + tol]
    return tied, Jmin, (min(rest) if rest else float('nan'))


def mesh_check(ideal_cat, new_cat, A_out, nv):
    """bond mesh = the 12 nearest cation neighbours of every cation on the ideal fcc lattice.
    Returns (number of cations whose 12 nearest neighbours in the final structure are NOT the mesh neighbours,
    smallest ratio d13/d12 over all cations = margin before a mesh bond would be exchanged)."""
    N = len(ideal_cat)
    dlt = new_cat - ideal_cat                                        # unwrap: every atom next to its own lattice point
    dlt -= np.round(dlt)
    new_cat = ideal_cat + dlt
    u = ideal_cat * nv
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    self_idx = 13                                                    # shift (0,0,0)
    bad, margin = 0, np.inf
    for i in range(N):
        vi = u[None, :, :] - u[i] + shifts[:, None, :] * nv
        mesh = set(np.nonzero((np.abs((vi ** 2).sum(-1) - 0.5) < 1e-6).ravel())[0].tolist())
        if len(mesh) != 12:
            return -1, float('nan')
        vf = (new_cat[None, :, :] - new_cat[i] + shifts[:, None, :]) @ A_out
        d = np.linalg.norm(vf, axis=-1)
        d[self_idx, i] = np.inf
        flat = d.ravel()
        order = np.argsort(flat)
        if set(order[:12].tolist()) != mesh:
            bad += 1
        margin = min(margin, flat[order[12]] / flat[order[11]])
    return bad, float(margin)


def build_variant(g, t, tcell, Tc, To, idealc, idealo, cat_site, o_site, qC, qO, nv):
    perm = [int(np.nonzero(g[i])[0][0]) for i in range(3)]
    sgn = [int(g[i, perm[i]]) for i in range(3)]
    Anew = np.array([sgn[i] * tcell[perm[i]] for i in range(3)])
    mirrored = bool(np.linalg.det(Anew) < 0)
    if mirrored:
        # right-handed cell by mirroring the CELL in Cartesian space: fractional coordinates, and hence the
        # cation ordering on the supercell lattice, stay untouched (the target is the template's mirror image,
        # the same crystal for these achiral phases)
        Anew = Anew @ np.diag([1.0, 1.0, -1.0])
    new_cat = idealc[cat_site] + (Tc[codes((qC - t) @ g)] @ g.T) / nv
    new_O = idealo[o_site] + (To[codes((qO - t) @ g)] @ g.T) / nv
    return dict(perm=perm, sgn=sgn, mirrored=mirrored, Anew=Anew, pos=np.vstack([new_cat, new_O]))


def distinct_variants(variants, nv):
    """drop variants congruent to an earlier one (global translation, rotation or mirror)"""
    out = []
    for v in variants:
        for w in out:
            # congruent structures (same metric + same fractional positions up to a translation, e.g. mirror images)
            # have identical energies and count as one
            if not np.allclose(v['Anew'] @ v['Anew'].T, w['Anew'] @ w['Anew'].T, atol=0.05):      # A^2; references carry ~0.01 skew noise
                continue
            d = v['pos'] - w['pos']
            d -= np.round(d)
            e = d - d[0]
            e -= np.round(e)
            if np.linalg.norm(e @ (v['Anew'] * nv[:, None]), axis=1).max() < 0.02:     # A
                break
        else:
            out.append(v)
    return out


# ----------------------------------------------------------------------------- checks
def min_distance(cell, frac):
    best = 99.0
    for sh in itertools.product((-1, 0, 1), repeat=3):
        d = (frac[:, None, :] - frac[None, :, :] + np.array(sh)) @ cell
        r = np.linalg.norm(d, axis=2)
        if sh == (0, 0, 0):
            r = r + np.eye(len(frac)) * 99
        best = min(best, r.min())
    return best


def path_check(ideal, final, cell, nv, steps=6):
    """distort the ideal lattice into the final structure in `steps` equal steps; return the closest
    approach of any two atoms along the way (collision test; atoms keep their lattice identity)"""
    disp = final - ideal
    disp -= np.round(disp)
    A_end = np.linalg.cholesky(cell @ cell.T)             # same metric, standard orientation, det > 0
    a0 = (abs(np.linalg.det(cell)) / np.prod(nv)) ** (1 / 3)
    A0 = np.diag(nv * a0)
    return [min_distance((1 - l) * A0 + l * A_end, ideal + l * disp) for l in np.linspace(0, 1, steps)]


def same_decoration(u1, sp1, u2, sp2, tol=1e-5):
    """species decoration of two sets of lattice points identical up to a global translation?"""
    if sorted(sp1) != sorted(sp2):
        return False
    for j in range(len(u2)):
        if sp2[j] != sp1[0]:
            continue
        sh = u1 + (u2[j] - u1[0])
        d = sh[:, None, :] - u2[None, :, :]
        d -= np.round(d)
        m = np.abs(d).max(axis=2) < tol
        if m.sum(1).min() < 1:
            continue
        idx = m.argmax(1)
        if len(set(idx)) == len(idx) and np.all(sp2[idx] == sp1):
            return True
    return False


def decoration_on_lattice(cell, frac, sp, n):
    iscat = sp != 'O'
    sh, _ = find_frame(frac[iscat], cell, n)
    ic, _ = ideal_sites(sh, n)
    site, _ = match(frac[iscat], ic, cell)
    return ic[site], sp[iscat]


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('poscar')
    ap.add_argument('--template', required=True, help='12-atom POSCAR of the target phase')
    ap.add_argument('--n', type=int, nargs=3, metavar=('NX', 'NY', 'NZ'), help='supercell repeats (default: auto)')
    ap.add_argument('--out', default=None)
    ap.add_argument('--identity', action='store_true', help='skip the orientation search')
    ap.add_argument('--keep-volume', action='store_true', help='rescale output to the input volume per f.u.')
    ap.add_argument('--variant', type=int, default=1, help='which of the equally good variants to write (default 1)')
    ap.add_argument('--all-variants', action='store_true', help='also write every distinct equally good variant as OUT_v1, OUT_v2, ...')
    ap.add_argument('--label', default=None,
                    help='phase label written into the title line (default: template name without POSCAR_)')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    say = (lambda *a: None) if args.quiet else print

    cell, frac, sp = read_poscar(args.poscar)
    n = tuple(args.n) if args.n else detect_n(cell, frac, sp)
    nv = np.array(n)
    iscat = sp != 'O'
    ncat, nO = int(iscat.sum()), int((~iscat).sum())
    if ncat != 4 * int(np.prod(n)) or nO != 8 * int(np.prod(n)):
        die('atom counts (%d cations, %d O) do not match a %dx%dx%d repeat of the 12-atom cell '
            '(needs %d + %d)' % (ncat, nO, n[0], n[1], n[2], 4 * np.prod(n), 8 * np.prod(n)))
    comp = {str(s): int((sp == s).sum()) for s in dict.fromkeys(sp)}
    say('input   : %s  %s  supercell %dx%dx%d%s' % (args.poscar, comp, n[0], n[1], n[2], '' if args.n else ' (auto)'))

    # 1. ideal fluorite frame of the input
    cat, O = frac[iscat], frac[~iscat]
    shift, rms = find_frame(cat, cell, n)
    idealc, idealo = ideal_sites(shift, n)
    cat_site, dcat = match(cat, idealc, cell)
    o_site, dO_ = match(O, idealo, cell)
    say('frame   : max cation displacement %.3f A, max O displacement %.3f A (from ideal fluorite sites)' %
        (dcat.max(), dO_.max()))
    if dcat.max() > MAX_CAT_DISP or dO_.max() > MAX_O_DISP:
        die('input is not a small distortion of fluorite for a %dx%dx%d repeat (cation %.2f A, O %.2f A from ideal '
            'sites). Check --n and that the cell axes follow the conventional 4-cation cell.' %
            (n[0], n[1], n[2], dcat.max(), dO_.max()))

    # 2. target pattern and orientation
    tcell, Tc, To = template_fields(args.template)
    dCf = (cat - idealc[cat_site])
    dCf -= np.round(dCf)
    dOf = (O - idealo[o_site])
    dOf -= np.round(dOf)
    qC = idealc[cat_site] * nv - shift
    qO = idealo[o_site] * nv - shift
    f_in_all = np.empty_like(frac)
    f_in_all[:ncat] = cat
    f_in_all[ncat:] = O
    tied, Jmin, Jnext = rank_correspondences(cell, f_in_all, tcell, Tc, To, idealc, idealo, cat_site, o_site,
                                             qC, qO, nv, args.identity)
    variants = distinct_variants(tied, nv)
    if not 1 <= args.variant <= len(variants):
        die('--variant %d requested but only %d distinct variants exist' % (args.variant, len(variants)))
    if args.identity:
        say('orient  : identity correspondence (search skipped)')
    elif np.isnan(Jnext):
        say('orient  : all 192 correspondences need the same distortion')
    else:
        say('orient  : smallest-distortion correspondence of 192 (score %.3f A rms-equivalent; next best %.3f A)' %
            (np.sqrt(Jmin), np.sqrt(Jnext)))
    if len(variants) > 1:
        say('variants: %d structures need exactly the same distortion; they differ only through the cation ordering, '
            'so distortion cannot rank them (--all-variants writes them all for an energy comparison). Using #%d.' %
            (len(variants), args.variant))
    else:
        say('variants: 1 (unique minimum)')

    out = args.out or 'POSCAR_%s' % args.template.split('/')[-1].replace('POSCAR_', '')
    # keep the input's first line (e.g. composition + Warren-Cowley alpha): the cation ordering on the lattice is
    # identical, so lattice-based SRO parameters are unchanged by the mapping
    label = args.label or args.template.split('/')[-1].replace('POSCAR_', '')
    first = open(args.poscar).readline().strip()
    title = '%s | phase=%s | mapped from %s onto %s (cation ordering kept)' % (
        first, label, args.poscar.split('/')[-1], args.template.split('/')[-1])
    todo = [(args.variant, False)]
    if args.all_variants:
        todo += [(k + 1, True) for k in range(len(variants))]
    all_problems = []
    for k, to_variant_file in todo:
        v = variants[k - 1]
        verbose = not to_variant_file
        Anew = v['Anew']
        Atot = Anew * nv[:, None]
        vin = abs(np.linalg.det(cell)) / ncat
        if args.keep_volume:
            Atot = Atot * (vin / (abs(np.linalg.det(Atot)) / ncat)) ** (1 / 3)
        vout = abs(np.linalg.det(Atot)) / ncat
        newfrac = np.empty_like(frac)
        newfrac[iscat] = v['pos'][:ncat]
        newfrac[~iscat] = v['pos'][ncat:]
        if verbose:
            lens = np.linalg.norm(Anew, axis=1)
            say('axes    : input a,b,c <- template axes %s (signs %s)%s, conventional-cell lengths %s A' %
                (['abc'[p] for p in v['perm']], v['sgn'], ', mirrored' if v['mirrored'] else '', np.round(lens, 3)))
            say('volume  : %.2f A^3 per f.u. (input %.2f)%s' % (vout, vin, ', rescaled to the input volume' if args.keep_volume else ''))
        # self-verification
        problems = []
        if np.linalg.det(Atot) <= 0:
            problems.append('left-handed output cell')
        mind = min_distance(Atot, newfrac)
        path = path_check(np.vstack([idealc[cat_site], idealo[o_site]]), v['pos'], Atot, nv)
        if mind < 1.5 or min(path) < 0.7 * path[0]:
            problems.append('atoms collide (closest approach %.2f A)' % mind)
        u_in, s_in = decoration_on_lattice(cell, frac, sp, n)
        u_out, s_out = decoration_on_lattice(Atot, newfrac, sp, n)
        order_ok = same_decoration(u_in, s_in, u_out, s_out)
        if not order_ok:
            problems.append('cation ordering not preserved')
        if verbose:
            dis = distortion(cell, f_in_all, Atot, np.vstack([newfrac[iscat], newfrac[~iscat]]), nv)
            say('distort : rms atomic shift %.3f A (max %.3f A); cell strain S = %.4f, principal stretches %s' %
                (dis['drms'], dis['dmax'], dis['S'], np.round(dis['stretches'], 4)))
            nbad, margin = mesh_check(idealc[cat_site], newfrac[iscat], Atot, nv)
            if nbad < 0:
                say('mesh    : bond-mesh check skipped (supercell too small to identify the 12 neighbours)')
            elif nbad == 0:
                say('mesh    : the 12 cation-cation bonds of every cation are the same atoms as in the cubic mesh '
                    '(no bond exchanged; margin d13/d12 = %.2f)' % margin)
            else:
                say('mesh    : WARNING %d of %d cations changed nearest neighbours relative to the cubic mesh' % (nbad, ncat))
            say('check   : closest approach %.2f A; along the distortion path %.2f A (start %.2f A)' % (mind, min(path), path[0]))
            say('check   : cation ordering on the lattice identical to the input: %s' % ('YES' if order_ok else 'NO'))
        if problems:
            all_problems.append('variant %d: %s' % (k, '; '.join(problems)))
            continue
        fn = '%s_v%d' % (out, k) if to_variant_file else out
        write_poscar(fn, Atot, sp, newfrac, title + (' variant %d' % k if to_variant_file else ''))
        say('wrote   : %s%s' % (fn, ' (atoms keep the order of the input file)' if verbose else ''))
    if all_problems:
        die('; '.join(all_problems))


if __name__ == '__main__':
    main()
