#!/usr/bin/env python3
"""
batch_map.py - run map_fluorite.py on every POSCAR in a folder, for one or several target phases.

  python3 batch_map.py INPUT_DIR --templates POSCAR_t POSCAR_oIII [--out-dir DIR] [--glob PATTERN]
                       [--jobs N] [--args "--keep-volume"]

* Every file matched in INPUT_DIR (not recursive; default pattern 'POSCAR*') is mapped onto every template.
* Results go to one sub-folder per target phase, with the phase label in the file name:
      OUT_DIR/t/POSCAR_t_0001_WC-0.2500    OUT_DIR/oIII/POSCAR_oIII_0001_WC-0.2500   ...
  (default OUT_DIR = INPUT_DIR/mapped; a name not starting with POSCAR_ gets _<label> appended)
* The first line of every output POSCAR is the input's first line (composition, Warren-Cowley alpha, id)
  followed by '| phase=<label> | mapped from ...'.  The cation ordering on the lattice is kept, so the
  lattice-based Warren-Cowley values carry over unchanged.
* The phase label is the template name without the leading 'POSCAR_' (POSCAR_oIII -> oIII);
  use label=path (e.g. tetra=refs/POSCAR_t) to choose your own label.
* A file that cannot be mapped (not a POSCAR, not fluorite-like, wrong atom counts, ...) is reported and skipped;
  the others carry on.  OUT_DIR/summary.csv lists every file/target with status, supercell, rms shift, strain and the
  number of equally good variants.
* --args passes extra options to map_fluorite.py for every run, e.g. "--keep-volume" or "--n 2 2 2".
* --skip-existing (restart): a file/target pair whose output exists AND is listed as ok in the existing
  OUT_DIR/summary.csv is not mapped again; its summary row is kept.  Everything else is (re)mapped.
* map_fluorite.py must sit in the same folder as this script.  Exit status 1 if any run failed.
"""
import argparse
import csv
import glob
import os
import re
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, 'map_fluorite.py')


def parse_template(spec):
    if '=' in spec:
        label, path = spec.split('=', 1)
    else:
        path = spec
        base = os.path.basename(path)
        label = re.sub(r'^POSCAR[_.\-]?', '', base) or base
    return label, path


def out_name(src, label):
    """POSCAR_0001_WC-0.2500 -> POSCAR_t_0001_WC-0.2500 ; any other name -> <name>_t"""
    base = os.path.basename(src)
    if base.startswith('POSCAR_'):
        return 'POSCAR_%s_%s' % (label, base[len('POSCAR_'):])
    if base == 'POSCAR':
        return 'POSCAR_%s' % label
    return '%s_%s' % (base, label)


def read_wc(src):
    """Warren-Cowley entries (e.g. 'WC_shell1(Ce-Zr)=-0.2500') from the first line of a POSCAR"""
    try:
        first = open(src).readline()
    except OSError:
        return ''
    return ' '.join(re.findall(r'WC_shell\d+\([^)]*\)=[+-]?[\d.]+', first))


def run_one(job):
    src, label, tpl, outdir, extra = job
    outfile = os.path.join(outdir, label, out_name(src, label))
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    cmd = [sys.executable, SCRIPT, src, '--template', tpl, '--out', outfile, '--label', label] + extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    row = dict(file=os.path.basename(src), target=label, out_file=os.path.basename(outfile), wc=read_wc(src),
               status='ok', supercell='', rms_shift_A='', strain_S='', tied_variants='', note='')
    out = r.stdout
    m = re.search(r'supercell (\d+x\d+x\d+)', out)
    if m:
        row['supercell'] = m.group(1)
    m = re.search(r'rms atomic shift ([\d.]+) A.*S = ([\d.]+)', out)
    if m:
        row['rms_shift_A'], row['strain_S'] = m.group(1), m.group(2)
    m = re.search(r'^variants: (\d+)', out, re.M)
    if m:
        row['tied_variants'] = m.group(1)
    if r.returncode != 0:
        row['status'] = 'FAILED'
        err = [l for l in r.stderr.strip().split('\n') if l.strip()]
        row['note'] = (err[-1] if err else 'exit status %d' % r.returncode).replace('ERROR: ', '')[:300]
        if os.path.exists(outfile):
            os.remove(outfile)
    elif row['tied_variants'] not in ('', '1'):
        row['note'] = '%s equally good variants: the written one is #1; see README / --all-variants' % row['tied_variants']
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input_dir')
    ap.add_argument('--templates', nargs='+', required=True, help='target-phase 12-atom POSCARs (or label=path)')
    ap.add_argument('--out-dir', default=None)
    ap.add_argument('--glob', default='POSCAR*',
                    help="which files in INPUT_DIR to use (default 'POSCAR*', so summary.csv etc. are skipped)")
    ap.add_argument('--jobs', type=int, default=min(4, os.cpu_count() or 1), help='parallel runs')
    ap.add_argument('--args', default='', help='extra options for map_fluorite.py, as one quoted string')
    ap.add_argument('--skip-existing', action='store_true',
                    help='restart: keep outputs that exist and are ok in OUT_DIR/summary.csv')
    argv, merged, i = sys.argv[1:], [], 0
    while i < len(argv):                       # allow:  --args "--keep-volume"   (argparse would take it for an option)
        if argv[i] == '--args' and i + 1 < len(argv):
            merged.append('--args=' + argv[i + 1])
            i += 2
        else:
            merged.append(argv[i])
            i += 1
    args = ap.parse_args(merged)

    if not os.path.isfile(SCRIPT):
        sys.exit('ERROR: map_fluorite.py not found next to batch_map.py (%s)' % SCRIPT)
    if not os.path.isdir(args.input_dir):
        sys.exit('ERROR: %s is not a folder' % args.input_dir)
    templates = [parse_template(t) for t in args.templates]
    for label, path in templates:
        if not os.path.isfile(path):
            sys.exit('ERROR: template %s not found' % path)
    labels = [l for l, _ in templates]
    if len(set(labels)) != len(labels):
        sys.exit('ERROR: two templates have the same label %s; use label=path' % labels)
    outdir = args.out_dir or os.path.join(args.input_dir, 'mapped')
    extra = shlex.split(args.args)
    if any(a in extra for a in ('--out', '--template')):
        sys.exit('ERROR: do not pass --out/--template in --args (batch_map sets them)')

    skip = {os.path.realpath(p) for _, p in templates}
    files = []
    for f in sorted(glob.glob(os.path.join(args.input_dir, args.glob))):
        if (not os.path.isfile(f)) or os.path.basename(f).startswith('.'):
            continue
        if os.path.realpath(f) in skip or os.path.realpath(f).startswith(os.path.realpath(outdir) + os.sep):
            continue
        files.append(f)
    if not files:
        sys.exit('ERROR: no files matching %r in %s' % (args.glob, args.input_dir))

    jobs = [(f, label, path, outdir, extra) for f in files for label, path in templates]
    rows = []
    if args.skip_existing:
        old = {}
        sp = os.path.join(outdir, 'summary.csv')
        if os.path.isfile(sp):
            with open(sp) as fh:
                old = {(r['file'], r['target']): r for r in csv.DictReader(fh) if r.get('status') == 'ok'}
        keep = []
        for j in jobs:
            r = old.get((os.path.basename(j[0]), j[1]))
            if r is not None and os.path.isfile(os.path.join(outdir, j[1], r['out_file'])):
                rows.append(r)
            else:
                keep.append(j)
        if rows:
            print('%d of %d runs already done (--skip-existing)' % (len(rows), len(jobs)))
        jobs = keep
    print('%d input files x %d target phases = %d runs, %d in parallel; output in %s' %
          (len(files), len(templates), len(jobs), args.jobs, outdir))
    os.makedirs(outdir, exist_ok=True)
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        for i, row in enumerate(ex.map(run_one, jobs), 1):
            rows.append(row)
            detail = ('rms shift %s A, strain %s' % (row['rms_shift_A'], row['strain_S'])) if row['status'] == 'ok' else row['note']
            print('[%d/%d] %-30s -> %-6s %s  %s' % (i, len(jobs), row['file'], row['target'], row['status'], detail))

    cols = ['file', 'target', 'out_file', 'wc', 'status', 'supercell', 'rms_shift_A', 'strain_S', 'tied_variants', 'note']
    with open(os.path.join(outdir, 'summary.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r['file'], r['target'])))
    failed = [r for r in rows if r['status'] != 'ok']
    tied = [r for r in rows if r['status'] == 'ok' and r['tied_variants'] not in ('', '1')]
    print('\ndone: %d ok, %d failed. Summary: %s' % (len(rows) - len(failed), len(failed), os.path.join(outdir, 'summary.csv')))
    if tied:
        print('%d runs have several equally good variants (listed in the summary); only variant 1 was written.' % len(tied))
    for r in failed:
        print('  FAILED %s -> %s: %s' % (r['file'], r['target'], r['note']))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
