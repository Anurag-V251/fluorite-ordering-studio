#!/usr/bin/env python3
"""
transform_set.py - the physical phase-transformation partners of a cation-ordered fluorite supercell.

  cubic -> tetragonal -> orthorhombic oIII (Pca2_1)         and         tetragonal -> monoclinic

The tetragonal phase is the hub.  Any input (c, t, oIII or m, any supercell, any cation ordering) is first
mapped onto t with the smallest-distortion correspondence of map_fluorite.py; oIII, m (and c) are then made FROM
THAT t.  Because every step keeps the cation ordering on the fluorite lattice, the ordering is identical in all
outputs, and because every step is the minimum-distortion correspondence the orientation relation is the real one:

    t main axis (unique axis)  ->  LONGEST axis of oIII / of m (c_m)
    oIII polar axis            ->  one of the two in-plane t axes (never the t main axis)
    m unique axis (b_m)        ->  one of the two in-plane t axes

The equally good variants (which in-plane axis, which sign / shear sense, which origin) cannot be ranked by
distortion; use the options below to pick one, or --all-variants to write all of them for an energy comparison.

  python3 transform_set.py INPUT [INPUT ...] [--refs REF_DIR] [--targets t oIII m c] [--out-dir DIR]
          [--main-axis a|b|c] [--t-variant K]                       how the hub t is oriented (inputs that are not t)
          [--polar-axis a|b|c] [--polar-sign +|-]                    oIII variant
          [--unique-axis a|b|c] [--shear obtuse|acute]               monoclinic variant
          [--all-variants] [--jobs N]

* INPUT: POSCAR files and/or folders (all regular files of the folder, not recursive).
* REF_DIR: folder with the four 12-atom references POSCAR_c, POSCAR_t, POSCAR_oIII, POSCAR_m (default: refs/ next to
  this script if it exists, else the script's own folder).  map_fluorite.py must be next to this script too.
* Output: OUT_DIR/<phase>/<input name>  (with --all-variants: <input name>__<tag>), OUT_DIR/summary.csv.
  Axes a, b, c of every output are the lattice vectors a, b, c of the hub t, i.e. of the input.
* --main-axis (default b, like your references) only matters when several t orientations tie (cubic input).  An input
  that already is tetragonal / oIII / m has one smallest-distortion t; its main axis is kept (a note is printed).
* Variant tags: oIII  polar+c_r1 (polar axis, sign, origin r1/r2);  m  uniq-a_ang99.63_r1 (unique axis, angle between
  the other two axes, origin).  t:  main-b_r1.
* Exit status 1 if any input failed.
"""
import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import map_fluorite as M                      # noqa: E402

SCRIPT = os.path.join(HERE, 'map_fluorite.py')
DEFAULT_REFS = os.path.join(HERE, 'refs') if os.path.isfile(os.path.join(HERE, 'refs', 'POSCAR_t')) else HERE
AX = 'abc'
ALL_TARGETS = ('t', 'oIII', 'm', 'c')


class Fail(Exception):
    pass


# --------------------------------------------------------------------------------------------- helpers
def run_map(src, tpl, out):
    """map_fluorite.py with --all-variants; returns the variant files (in its order) and its log"""
    r = subprocess.run([sys.executable, SCRIPT, src, '--template', tpl, '--out', out, '--all-variants'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        err = [l for l in r.stderr.strip().split('\n') if l.strip()]
        raise Fail((err[-1] if err else 'map_fluorite.py exit status %d' % r.returncode).replace('ERROR: ', ''))
    files = sorted(glob.glob(out + '_v[0-9]*'), key=lambda f: int(f.rsplit('_v', 1)[1]))
    if not files:
        raise Fail('map_fluorite.py wrote no output')
    m = re.search(r'rms atomic shift ([\d.]+) A.*?S = ([\d.]+)', r.stdout)
    return files, (m.group(1), m.group(2)) if m else ('', '')


def describe(fn):
    """lattice lengths per 12-atom cell, angles, and the formal-charge dipole of the displacements from the
    ideal fluorite sites, projected on the lattice vectors a, b, c (per 12-atom cell)"""
    cell, frac, sp = M.read_poscar(fn)
    n = M.detect_n(cell, frac, sp)
    iscat = sp != 'O'
    shift, _ = M.find_frame(frac[iscat], cell, n)
    ic, io = M.ideal_sites(shift, n)
    cs, _ = M.match(frac[iscat], ic, cell)
    os_, _ = M.match(frac[~iscat], io, cell)
    dc = frac[iscat] - ic[cs]
    dc -= np.round(dc)
    do = frac[~iscat] - io[os_]
    do -= np.round(do)
    p = 4.0 * (dc @ cell).sum(0) - 2.0 * (do @ cell).sum(0)
    ln = np.linalg.norm(cell, axis=1)
    comp = np.array([p @ cell[i] / ln[i] for i in range(3)]) / float(np.prod(n))

    def ang(u, v):
        return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
    angles = (ang(cell[1], cell[2]), ang(cell[0], cell[2]), ang(cell[0], cell[1]))      # alpha, beta, gamma
    return dict(L=ln / np.array(n), dip=comp, ang=angles)


def tet_main_axis(L):
    """index of the axis whose length differs from the two equal ones (tetragonal), else None"""
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        if abs(L[j] - L[k]) < 2e-3 and abs(L[i] - L[j]) > 0.02:
            return i
    return None


def add_registry(items, keys):
    """r1, r2, ... = rank among the variants of the same (keys) group, in map_fluorite's order"""
    seen = {}
    for it in items:
        g = tuple(it[k] for k in keys)
        seen[g] = seen.get(g, 0) + 1
        it['reg'] = seen[g]
        it['tag'] = it['tag'] + '_r%d' % it['reg']


def classify(phase, files):
    items = []
    for f in files:
        d = describe(f)
        it = dict(file=f, L=d['L'])
        if phase == 't':
            i = tet_main_axis(d['L'])
            if i is None:
                raise Fail('the t reference does not have one unique axis: cannot define the main axis')
            it.update(main=AX[i], tag='main-' + AX[i])
        elif phase == 'oIII':
            i = int(np.argmax(np.abs(d['dip'])))
            if abs(d['dip'][i]) < 0.5:
                raise Fail('the oIII reference is not polar: cannot define the polar axis')
            sign = '+' if d['dip'][i] > 0 else '-'
            it.update(polar=AX[i], sign=sign, tag='polar%s%s' % (sign, AX[i]))
        elif phase == 'm':
            dev = [abs(a - 90.0) for a in d['ang']]
            i = int(np.argmax(dev))
            if dev[i] < 0.5:
                raise Fail('the monoclinic reference has no angle different from 90 degrees')
            it.update(uniq=AX[i], shear='obtuse' if d['ang'][i] > 90 else 'acute', angle=d['ang'][i],
                      tag='uniq-%s_ang%.2f' % (AX[i], d['ang'][i]))
        else:
            it.update(tag='c')
        items.append(it)
    keyset = {'t': ('main',), 'oIII': ('polar', 'sign'), 'm': ('uniq', 'shear'), 'c': ()}[phase]
    add_registry(items, keyset)
    return items


def pick(items, want, label):
    sel = [it for it in items if all(v is None or it[k] == v for k, v in want.items())]
    if not sel:
        have = sorted({', '.join('%s=%s' % (k, it[k]) for k in want) for it in items})
        asked = ', '.join('%s=%s' % (k, v) for k, v in want.items() if v is not None)
        raise Fail('no %s variant has %s; the equally good variants are: %s' % (label, asked, '; '.join(have)))
    return sel


def decoration(fn):
    cell, frac, sp = M.read_poscar(fn)
    n = M.detect_n(cell, frac, sp)
    u, species = M.decoration_on_lattice(cell, frac, sp, n)
    return np.array(n), u, species


def ordering_identical(a, b):
    na, ua, sa = decoration(a)
    nb, ub, sb = decoration(b)
    return bool(np.array_equal(na, nb) and M.same_decoration(ua, sa, ub, sb))


# --------------------------------------------------------------------------------------------- one input
def process(src, args, refs, outdir):
    name = os.path.basename(src)
    rows, notes = [], []
    tmp = tempfile.mkdtemp(prefix='tset_')
    try:
        def row(target, it, dest, extra, dist):
            rows.append(dict(file=name, target=target, tag=it['tag'], status='ok', output=os.path.relpath(dest, outdir),
                             main_axis=hub_main, detail=extra,
                             cell_abc_A='%.3f %.3f %.3f' % tuple(it['L']), step_rms_shift_A=dist[0], step_strain_S=dist[1],
                             ordering_vs_input='IDENTICAL', note=''))

        def write(target, it, multi, dist, extra):
            dest = os.path.join(outdir, target, name + ('__' + it['tag'] if multi else ''))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copyfile(it['file'], dest)
            if not ordering_identical(src, dest):
                raise Fail('internal check failed: cation ordering of %s differs from the input' % dest)
            row(target, it, dest, extra, dist)

        # ---- hub: tetragonal
        tfiles, tdist = run_map(src, refs['t'], os.path.join(tmp, 't'))
        titems = classify('t', tfiles)
        cands = [it for it in titems if it['main'] == args.main_axis]
        if not cands:
            mains = sorted({it['main'] for it in titems})
            if len(mains) == 1:
                cands = titems
                notes.append('input fixes the tetragonal main axis to %s (asked %s): kept' % (mains[0], args.main_axis))
            else:
                raise Fail('no tetragonal orientation with main axis %s; available: %s' % (args.main_axis, mains))
        if args.t_variant > len(cands):
            raise Fail('--t-variant %d: only %d tetragonal variant(s) with this main axis' % (args.t_variant, len(cands)))
        hub = cands[args.t_variant - 1]
        hub_main = hub['main']
        if len(cands) > 1:
            notes.append('%d equally good tetragonal variants with main axis %s; used #%d (--t-variant to change)'
                         % (len(cands), hub_main, args.t_variant))
        if 't' in args.targets:
            write('t', hub, False, tdist, 'main axis %s' % hub_main)

        # ---- oIII, m, c from the hub
        for target in args.targets:
            if target == 't':
                continue
            files, dist = run_map(hub['file'], refs[target], os.path.join(tmp, target))
            items = classify(target, files)
            if target == 'oIII':
                want = dict(polar=args.polar_axis, sign=args.polar_sign)
            elif target == 'm':
                want = dict(uniq=args.unique_axis, shear=args.shear)
            else:
                want = {}
            sel = pick(items, want, target) if want else items
            if not args.all_variants:
                sel = sel[:1]
            if len(items) > 1 and not args.all_variants and not any(v is not None for v in want.values()):
                notes.append('%s: %d equally good variants, wrote #1 (%s); see --all-variants / selection options'
                             % (target, len(items), items[0]['tag']))
            for it in sel:
                longest = AX[int(np.argmax(it['L']))] if target in ('oIII', 'm') else ''
                if target == 'oIII':
                    extra = 'polar axis %s%s, longest %s' % (it['sign'], it['polar'], longest)
                    ok = longest == hub_main and it['polar'] != hub_main
                elif target == 'm':
                    extra = 'unique axis %s, angle %.2f, longest %s' % (it['uniq'], it['angle'], longest)
                    ok = longest == hub_main and it['uniq'] != hub_main
                else:
                    extra, ok = '', True
                write(target, it, args.all_variants, dist, extra)
                if not ok:
                    rows[-1]['note'] = 'WARNING: orientation differs from the usual relation (longest axis %s, main axis %s)' % (longest, hub_main)
        if notes:
            rows[0]['note'] = (rows[0]['note'] + ' ' if rows[0]['note'] else '') + '; '.join(notes)
        return rows
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='+', help='POSCAR files and/or folders')
    ap.add_argument('--refs', default=DEFAULT_REFS, help='folder with POSCAR_c, POSCAR_t, POSCAR_oIII, POSCAR_m '
                    '(default: refs/ next to this script if it exists, else the folder of this script)')
    ap.add_argument('--targets', nargs='+', choices=ALL_TARGETS, default=['t', 'oIII', 'm'])
    ap.add_argument('--out-dir', default='transformed')
    ap.add_argument('--main-axis', choices=list(AX), default='b', help='tetragonal main axis when several tie (default b)')
    ap.add_argument('--t-variant', type=int, default=1, help='which tetragonal variant with that main axis (default 1)')
    ap.add_argument('--polar-axis', choices=list(AX), default=None)
    ap.add_argument('--polar-sign', choices=['+', '-'], default=None)
    ap.add_argument('--unique-axis', choices=list(AX), default=None)
    ap.add_argument('--shear', choices=['obtuse', 'acute'], default=None, help='angle between the two non-unique axes of m: 99.63 (obtuse) or 80.37 (acute)')
    ap.add_argument('--all-variants', action='store_true', help='write every equally good oIII / m variant (matching the selection)')
    ap.add_argument('--jobs', type=int, default=min(4, os.cpu_count() or 1))
    args = ap.parse_args()

    if not os.path.isfile(SCRIPT):
        sys.exit('ERROR: map_fluorite.py not found next to transform_set.py (%s)' % SCRIPT)
    refs = {}
    for ph in ('c', 't', 'oIII', 'm'):
        refs[ph] = os.path.join(args.refs, 'POSCAR_' + ph)
        if ph in set(args.targets) | {'t'} and not os.path.isfile(refs[ph]):
            sys.exit('ERROR: reference %s not found (use --refs DIR with POSCAR_c, POSCAR_t, POSCAR_oIII, POSCAR_m)' % refs[ph])
    if args.polar_axis or args.polar_sign:
        if 'oIII' not in args.targets:
            sys.exit('ERROR: --polar-axis / --polar-sign need oIII in --targets')
    if args.unique_axis or args.shear:
        if 'm' not in args.targets:
            sys.exit('ERROR: --unique-axis / --shear need m in --targets')
    outdir = os.path.abspath(args.out_dir)

    skip = {os.path.realpath(p) for p in refs.values()}
    files = []
    for item in args.inputs:
        if os.path.isdir(item):
            cand = sorted(glob.glob(os.path.join(item, '*')))
        elif os.path.isfile(item):
            cand = [item]
        else:
            sys.exit('ERROR: %s is not a file or folder' % item)
        for f in cand:
            rp = os.path.realpath(f)
            if (not os.path.isfile(f)) or os.path.basename(f).startswith('.') or rp in skip \
                    or rp.startswith(outdir + os.sep) or os.path.basename(f) == 'summary.csv':
                continue
            files.append(f)
    if not files:
        sys.exit('ERROR: no input files')
    names = [os.path.basename(f) for f in files]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        sys.exit('ERROR: several inputs have the same file name %s; outputs would overwrite each other' % dup)

    print('%d input file(s), targets %s, %d in parallel; output in %s' % (len(files), ' '.join(args.targets), args.jobs, outdir))
    os.makedirs(outdir, exist_ok=True)
    allrows, failed = [], []

    def job(f):
        try:
            return f, process(f, args, refs, outdir), None
        except Fail as e:
            return f, [], str(e)
        except Exception as e:                                   # unreadable file etc.
            return f, [], '%s: %s' % (type(e).__name__, e)

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        for i, (f, rows, err) in enumerate(ex.map(job, files), 1):
            name = os.path.basename(f)
            if err:
                failed.append((name, err))
                allrows.append(dict(file=name, target='', tag='', status='FAILED', output='', main_axis='', detail='',
                                    cell_abc_A='', step_rms_shift_A='', step_strain_S='', ordering_vs_input='', note=err[:300]))
                print('[%d/%d] %-30s FAILED: %s' % (i, len(files), name, err))
                continue
            allrows += rows
            print('[%d/%d] %s' % (i, len(files), name))
            for r in rows:
                print('        %-5s %-26s %-44s step rms shift %s A, S %s  ordering %s' % (
                    r['target'], r['tag'], r['detail'], r['step_rms_shift_A'], r['step_strain_S'], r['ordering_vs_input']))
                if r['note']:
                    print('        note: %s' % r['note'])

    cols = ['file', 'target', 'tag', 'status', 'output', 'main_axis', 'detail', 'cell_abc_A', 'step_rms_shift_A',
            'step_strain_S', 'ordering_vs_input', 'note']
    with open(os.path.join(outdir, 'summary.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(allrows)
    print('\ndone: %d ok, %d failed. Summary: %s' % (len(files) - len(failed), len(failed), os.path.join(outdir, 'summary.csv')))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
