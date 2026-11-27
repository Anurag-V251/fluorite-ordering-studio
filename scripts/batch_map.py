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
* Speed: the --jobs workers are Python processes that import map_fluorite.py ONCE and then map file after file
  (starting a fresh Python + numpy/scipy for every mapping cost 2/3 of the time, much more on a network filesystem).
  Each worker is limited to --threads BLAS/OpenMP threads (default 1), so N workers use N cores.
  --subprocess restores the old one-process-per-mapping mode.
* Progress: one line per run (default), or with --progress S one line every S seconds with the count, rate and
  ETA (failed runs are always printed at once).  Output is flushed, so it shows up live in logs and pipes.
  summary.csv is rewritten every 60 s while running, so a killed run restarts with --skip-existing where it stopped.
* map_fluorite.py must sit in the same folder as this script.  Exit status 1 if any run failed.
"""
import argparse
import contextlib
import csv
import glob
import io
import os
import re
import shlex
import signal
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, 'map_fluorite.py')
THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
SAVE_EVERY = float(os.environ.get('BATCH_MAP_SAVE_EVERY', 60))   # s between summary.csv saves (tests use 0)
COLS = ['file', 'target', 'out_file', 'wc', 'status', 'supercell', 'rms_shift_A', 'strain_S', 'tied_variants', 'note']


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


def map_in_process(cmd_args):
    """run map_fluorite.main() in this (worker) process; returns (returncode, stdout, stderr)"""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import map_fluorite                       # imported once per worker, then cached
    out, err = io.StringIO(), io.StringIO()
    argv, code = sys.argv, 0
    sys.argv = ['map_fluorite.py'] + cmd_args
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            map_fluorite.main()
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        if e.code is not None and not isinstance(e.code, int):
            err.write(str(e.code) + '\n')
    except Exception:                         # a bug or odd input in one file must not stop the batch
        code = 1
        err.write(traceback.format_exc())
    finally:
        sys.argv = argv
    return code, out.getvalue(), err.getvalue()


def worker_init():
    """pool worker: die with the parent (kill / scancel of batch_map must not leave orphans behind; Linux)"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)          # Ctrl+C is handled by the parent
    try:
        import ctypes
        ctypes.CDLL('libc.so.6', use_errno=True).prctl(1, signal.SIGKILL)   # PR_SET_PDEATHSIG
    except (OSError, AttributeError):
        pass
    if os.getppid() == 1:                                 # parent already gone
        os._exit(1)


def run_one(job):
    src, label, tpl, outdir, extra, use_subprocess = job
    outfile = os.path.join(outdir, label, out_name(src, label))
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    cmd_args = [src, '--template', tpl, '--out', outfile, '--label', label] + extra
    if use_subprocess:
        r = subprocess.run([sys.executable, SCRIPT] + cmd_args, capture_output=True, text=True)
        code, out, stderr = r.returncode, r.stdout, r.stderr
    else:
        code, out, stderr = map_in_process(cmd_args)
    row = dict(file=os.path.basename(src), target=label, out_file=os.path.basename(outfile), wc=read_wc(src),
               status='ok', supercell='', rms_shift_A='', strain_S='', tied_variants='', note='')
    m = re.search(r'supercell (\d+x\d+x\d+)', out)
    if m:
        row['supercell'] = m.group(1)
    m = re.search(r'rms atomic shift ([\d.]+) A.*S = ([\d.]+)', out)
    if m:
        row['rms_shift_A'], row['strain_S'] = m.group(1), m.group(2)
    m = re.search(r'^variants: (\d+)', out, re.M)
    if m:
        row['tied_variants'] = m.group(1)
    if code != 0:
        row['status'] = 'FAILED'
        err = [l for l in stderr.strip().split('\n') if l.strip()]
        row['note'] = (err[-1] if err else 'exit status %d' % code).replace('ERROR: ', '')[:300]
        if os.path.exists(outfile):
            os.remove(outfile)
    elif row['tied_variants'] not in ('', '1'):
        row['note'] = '%s equally good variants: the written one is #1; see README / --all-variants' % row['tied_variants']
    return row


def write_summary(path, rows):
    tmp = path + '.tmp'
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction='ignore')
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r['file'], r['target'])))
    os.replace(tmp, path)


def hms(s):
    s = int(round(s))
    return '%d:%02d:%02d' % (s // 3600, s % 3600 // 60, s % 60)


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
    ap.add_argument('--threads', type=int, default=1, help='BLAS/OpenMP threads per worker (default 1)')
    ap.add_argument('--progress', type=float, default=0, metavar='S',
                    help='print a progress line (count, rate, ETA) every S seconds instead of one line per run')
    ap.add_argument('--subprocess', action='store_true',
                    help='old mode: a fresh Python for every mapping (slower; isolates crashes)')
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
            print('%d of %d runs already done (--skip-existing)' % (len(rows), len(jobs)), flush=True)
        jobs = keep
    print('%d input files x %d target phases = %d runs, %d in parallel; output in %s' %
          (len(files), len(templates), len(jobs), args.jobs, outdir), flush=True)
    os.makedirs(outdir, exist_ok=True)
    summary = os.path.join(outdir, 'summary.csv')
    for var in THREAD_VARS:                   # read by numpy/scipy when the workers import them
        os.environ[var] = str(max(1, args.threads))
    jobs = [j + (args.subprocess,) for j in jobs]
    if args.subprocess:
        ex = ThreadPoolExecutor(max_workers=max(1, args.jobs))
    else:
        ex = ProcessPoolExecutor(max_workers=max(1, args.jobs), initializer=worker_init)

    def stop(signum, frame):                  # kill / scancel / Ctrl+C: keep what is done, end the workers
        write_summary(summary, rows)
        for proc in list((getattr(ex, '_processes', None) or {}).values()):
            try:
                proc.kill()
            except Exception:
                pass
        print('\nstopped (signal %d) after %d of %d runs; summary.csv saved, restart with --skip-existing'
              % (signum, len(rows) - nold, len(jobs)), flush=True)
        os._exit(128 + signum)

    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, stop)
    nold = len(rows)
    t0 = last_print = last_save = time.time()
    nfail = 0
    with ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            row = fut.result()
            rows.append(row)
            now = time.time()
            rate = i / max(now - t0, 1e-9)
            eta = hms((len(jobs) - i) / rate)
            ok = row['status'] == 'ok'
            nfail += not ok
            if args.progress <= 0 or not ok:
                detail = ('rms shift %s A, strain %s' % (row['rms_shift_A'], row['strain_S'])) if ok else row['note']
                print('[%d/%d] %-30s -> %-6s %s  %s  (ETA %s)' % (i, len(jobs), row['file'], row['target'],
                                                                 row['status'], detail, eta), flush=True)
            if args.progress > 0 and (now - last_print >= args.progress or i == len(jobs)):
                print('mapped %d/%d (%.1f %%), %d failed, %.2f runs/s, elapsed %s, ETA %s' %
                      (i, len(jobs), 100.0 * i / len(jobs), nfail, rate, hms(now - t0), eta), flush=True)
                last_print = now
            if now - last_save >= SAVE_EVERY:
                write_summary(summary, rows)
                last_save = now

    write_summary(summary, rows)
    failed = [r for r in rows if r['status'] != 'ok']
    tied = [r for r in rows if r['status'] == 'ok' and r['tied_variants'] not in ('', '1')]
    print('\ndone: %d ok, %d failed in %s. Summary: %s' % (len(rows) - len(failed), len(failed), hms(time.time() - t0), summary))
    if tied:
        print('%d runs have several equally good variants (listed in the summary); only variant 1 was written.' % len(tied))
    for r in failed:
        print('  FAILED %s -> %s: %s' % (r['file'], r['target'], r['note']))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
