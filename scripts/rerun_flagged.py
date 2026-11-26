#!/usr/bin/env python3
"""
rerun_flagged.py - recompute the structures of a campaign whose free energy is flagged (imaginary phonon modes
and/or unconverged relaxation), with stricter settings (default) or the same ones, on SLURM or locally.

Each call is one "round" R = 1, 2, ...:
  1. reads every result in <workdir>/results/ (latest row per structure, as campaign.py merge does) and picks
     the flagged ones (--what, --freq-window)
  2. start geometry: the newest relaxed structure of that file (from the previous round, else from the
     campaign run, else the unrelaxed input; --start original forces the input), copied to <workdir>/rerun_RRR/
     under the same relative path, so the file name - and with it the merge key - stays the same
  3. writes rerun_RRR/list_K.txt (split over --nodes array tasks) and rerun_RRR/rerun.sh (+ merge.sh), which run
     batch_free_energy.py --files-from list_K.txt --out <workdir>/results/chunk_rerunRRR_K.csv --relax ...
  4. --submit: sbatch the array and a merge job that starts after it (afterany).  --local: run here instead.
After the round, `campaign.py CONFIG merge` (done by the merge job) uses the newest row of every structure,
so rerun rows replace the old ones; a structure whose rerun FAILED keeps its old row.  Restartable: submitting
the same round's rerun.sh again skips structures already done in it.  Run the next round only after merging.

Settings (meV and Å as in the config):
  default    fmax = config fmax / 10,  steps = config steps x 10,  min_length = config min_length + 5 Å
  --same     exactly the config values (useful after a crash / walltime, or with --start original)
  --fmax/--steps/--min-length/--supercell/--mesh override single values (also together with --same)
Physics: imaginary modes that survive tighter settings and a larger phonon cell mean the structure is a saddle
point in that phase; its harmonic F is not meaningful and more reruns will not change that.  The rerun rows use
other settings than the rest of the campaign (usually < 1 meV/atom effect on F; the settings are in rerun.json).

Usage
  python scripts/rerun_flagged.py campaign.toml --dry-run                 # what would be rerun
  python scripts/rerun_flagged.py campaign.toml --submit                  # round with stricter settings
  python scripts/rerun_flagged.py campaign.toml --what unconverged --same --steps 20000 --submit
  python scripts/rerun_flagged.py campaign.toml --what imaginary --freq-window -0.5 0 --submit   # weak ones only
  python scripts/rerun_flagged.py campaign.toml --nodes 4 --walltime 06:00:00 --submit
  python scripts/rerun_flagged.py campaign.toml --local --jobs 8          # no SLURM
"""
import argparse
import csv
import json
import math
import os
import shlex
import shutil
import signal
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import campaign as C          # noqa: E402  (same scripts/ folder)


def latest_rows(cfg):
    """newest row per structure, exactly as campaign.py merge chooses it (ok beats FAILED, later beats earlier)"""
    best = {}
    for r in C.all_results(cfg):
        f = r["file"]
        if f not in best or r["status"] == "ok" or best[f]["status"] != "ok":
            best[f] = r
    return best


def is_imag(r):
    return int(float(r.get("n_imaginary_modes") or 0)) > 0


def is_unconv(r):
    return r.get("relax_converged") == "False"


def select(best, what, window):
    out = []
    for f, r in sorted(best.items()):
        if r["status"] != "ok":
            continue
        imag, unc = is_imag(r), is_unconv(r)
        if imag and window:
            mf = float(r["min_freq_THz"])
            imag = window[0] <= mf <= window[1]
        if (imag and what in ("both", "imaginary")) or (unc and what in ("both", "unconverged")):
            out.append((f, imag, unc))
    return out


def rounds(cfg):
    return sorted(int(p.name.split("_")[1]) for p in cfg.work.glob("rerun_[0-9][0-9][0-9]") if p.is_dir())


def start_geometry(cfg, rel, start):
    """newest relaxed geometry of rel: last rerun round, then the campaign run, then the input structure"""
    if start == "relaxed":
        for R in reversed(rounds(cfg)):
            p = cfg.work / f"rerun_{R:03d}" / "relaxed" / rel
            if p.is_file():
                return p, f"round {R}"
        p = cfg.work / "relaxed" / rel
        if p.is_file():
            return p, "campaign"
    return cfg.work / rel, "input"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("config", help="the campaign.toml of the campaign")
    ap.add_argument("--what", choices=["both", "imaginary", "unconverged"], default="both")
    ap.add_argument("--freq-window", type=float, nargs=2, metavar=("LO", "HI"),
                    help="only imaginary rows whose min_freq_THz is in [LO, HI], e.g. -0.5 0 (weak ones)")
    ap.add_argument("--same", action="store_true", help="the config's settings instead of stricter ones")
    ap.add_argument("--fmax", type=float, help="eV/Å (default config/10, or config with --same)")
    ap.add_argument("--steps", type=int, help="max relaxation steps (default config x 10, or config with --same)")
    ap.add_argument("--min-length", type=float, help="phonon supercell widths >= this, Å (default config + 5)")
    ap.add_argument("--supercell", type=int, nargs=3, metavar=("NX", "NY", "NZ"), help="fixed phonon supercell")
    ap.add_argument("--mesh", type=int, help="q mesh density (default config)")
    ap.add_argument("--start", choices=["relaxed", "original"], default="relaxed",
                    help="start from the newest relaxed geometry (default) or from the unrelaxed input")
    ap.add_argument("--nodes", type=int, default=1, help="split over this many array tasks / nodes (default 1)")
    ap.add_argument("--walltime", help="per task (default [cluster] walltime)")
    ap.add_argument("--dry-run", action="store_true", help="only print what would be rerun")
    ap.add_argument("--submit", action="store_true", help="sbatch the rerun array + a merge job after it")
    ap.add_argument("--no-merge", dest="merge", action="store_false", help="with --submit: no merge job")
    ap.add_argument("--local", action="store_true", help="run here (no SLURM), then merge")
    ap.add_argument("--jobs", type=int, help="--local: parallel structures (default: all cores)")
    a = ap.parse_args()
    if hasattr(signal, "SIGPIPE"):                       # `rerun_flagged.py ... | head` must not crash
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    cfg = C.Cfg(a.config)
    fe, cl = cfg["free_energy"], cfg["cluster"]
    if not cfg.nep or not cfg.nep.is_file():
        sys.exit(f"NEP model {cfg.nep} not found ([free_energy] nep)")
    best = latest_rows(cfg)
    if not best:
        sys.exit(f"no results in {cfg.work / 'results'}")
    sel = select(best, a.what, a.freq_window)
    n_imag = sum(1 for r in best.values() if r["status"] == "ok" and is_imag(r))
    n_unc = sum(1 for r in best.values() if r["status"] == "ok" and is_unconv(r))
    print(f"{len(best)} structures with results: {n_imag} with imaginary modes, {n_unc} unconverged "
          f"(newest row of each)")
    print(f"selected ({a.what}{f', min freq in {a.freq_window} THz' if a.freq_window else ''}): {len(sel)}  "
          f"[{sum(1 for s in sel if s[1])} imaginary, {sum(1 for s in sel if s[2])} unconverged]")
    if not sel:
        print("nothing to rerun"); return

    fmax = a.fmax if a.fmax is not None else float(fe["fmax"]) / (1 if a.same else 10)
    steps = a.steps if a.steps is not None else int(fe["steps"]) * (1 if a.same else 10)
    min_len = a.min_length if a.min_length is not None else float(fe["min_length"]) + (0 if a.same else 5)
    supercell = a.supercell or fe["supercell"] or []
    mesh = a.mesh if a.mesh is not None else int(fe["mesh"])
    settings = dict(fmax=fmax, steps=steps, min_length=min_len, supercell=list(supercell), mesh=mesh,
                    T=float(fe["T"]), fixed_cell=bool(fe["fixed_cell"]), threads=int(fe["threads"]),
                    start=a.start, what=a.what, freq_window=a.freq_window)
    starts = [start_geometry(cfg, f, a.start) for f, _, _ in sel]
    src_count = {}
    for _, s in starts:
        src_count[s] = src_count.get(s, 0) + 1
    print("settings: fmax {fmax:g} eV/Å, steps {steps}, phonon cell {sc}, mesh {mesh}, T {T:g} K".format(
        sc=("x".join(map(str, supercell)) if supercell else f"auto, widths >= {min_len:g} Å"), **settings))
    print("start geometries: " + ", ".join(f"{v} from {k}" for k, v in src_count.items()))
    if a.dry_run:
        return

    R = (rounds(cfg) or [0])[-1] + 1
    rd = cfg.work / f"rerun_{R:03d}"
    rd.mkdir(parents=True)
    for (f, _, _), (src, _) in zip(sel, starts):
        dst = rd / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    k = max(1, min(a.nodes, len(sel)))
    per = math.ceil(len(sel) / k)
    files = [f for f, _, _ in sel]
    for t in range(k):
        (rd / f"list_{t}.txt").write_text("\n".join(files[t * per:(t + 1) * per]) + "\n")
    (rd / "rerun.json").write_text(json.dumps(dict(round=R, n=len(sel), tasks=k, settings=settings,
                                                   files=files), indent=1))

    def bfe_cmd(list_file, out_csv, jobs):
        cmd = [cfg.py(), str(HERE / "batch_free_energy.py"), str(rd), "--files-from", str(list_file),
               "--out", str(out_csv), "--nep", str(cfg.nep), "--T", f"{settings['T']:g}", "--relax",
               "--fmax", f"{fmax:g}", "--steps", str(steps), "--min-length", f"{min_len:g}", "--mesh", str(mesh),
               "--threads", str(settings["threads"])]
        if settings["fixed_cell"]:
            cmd.append("--fixed-cell")
        if supercell:
            cmd += ["--supercell", *map(str, supercell)]
        return cmd + ["--jobs", jobs]

    res_dir = cfg.work / "results"
    res_dir.mkdir(exist_ok=True)
    if a.local:
        jobs = str(a.jobs or max(1, (os.cpu_count() or 1) // settings["threads"]))
        for t in range(k):
            subprocess.run(bfe_cmd(rd / f"list_{t}.txt", res_dir / f"chunk_rerun{R:03d}_{t}.csv", jobs), check=True)
        subprocess.run([cfg.py(), str(HERE / "campaign.py"), str(cfg.path), "merge"], check=True)
        return

    (cfg.work / "logs").mkdir(exist_ok=True)
    thr = settings["threads"]
    raw = (f'"{rd}/list_${{SLURM_ARRAY_TASK_ID}}.txt"',                       # expanded by bash at run time
           f'"{res_dir}/chunk_rerun{R:03d}_${{SLURM_ARRAY_TASK_ID}}.csv"',
           f'"$(( ${{SLURM_CPUS_ON_NODE:-{cl["cores_per_node"]}}} / {thr} ))"')
    body = " ".join(c if c in raw else shlex.quote(c) for c in bfe_cmd(*raw))
    sh = rd / "rerun.sh"
    sh.write_text(C.sbatch_header(cfg, f"fo_rerun{R}", a.walltime or cl["walltime"], cl["cores_per_node"],
                                  cfg.work / "logs" / f"rerun{R:03d}_%A_%a.out")
                  + "export RAYON_NUM_THREADS=1\n" + body + "\n")
    mg = rd / "merge.sh"
    me = f"{cfg.py()} {shlex.quote(str(HERE / 'campaign.py'))} {shlex.quote(str(cfg.path))}"
    mg.write_text(C.sbatch_header(cfg, f"fo_rerun{R}_merge", cl["merge_walltime"], 1,
                                  cfg.work / "logs" / f"rerun{R:03d}_merge_%j.out") + f"{me} merge\n")
    arr = ["sbatch", "--parsable", f"--array=0-{k - 1}", str(sh)]
    print(f"round {R}: {len(sel)} structures in {k} task(s) -> {rd}")
    print(f"wrote {sh}\n      {mg}")
    if not a.submit:
        print("submit with:\n  " + " ".join(arr) + "\n  sbatch --dependency=afterany:<jobid> " + str(mg))
        return
    jid = C.run(arr).stdout.strip().split(";")[0]
    print(f"  submitted rerun array {jid}")
    if a.merge:
        mid = C.run(["sbatch", "--parsable", f"--dependency=afterany:{jid}", str(mg)]).stdout.strip().split(";")[0]
        print(f"  submitted merge job {mid} (after {jid})")


if __name__ == "__main__":
    main()
