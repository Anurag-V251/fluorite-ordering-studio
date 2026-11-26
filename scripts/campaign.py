#!/usr/bin/env python3
"""
campaign.py - whole-composition-range free-energy campaign, driven by one config file.

For a host supercell (e.g. POSCAR_96: 32 Zr + 64 O) and every dopant count n in the config:
  1. enumerate_substitutions.py --n n --random MAX --seed SEED+n
       MAX distinct random orderings (uniform draw of distinct site sets, no symmetry reduction);
       all C(N, n) of them when fewer exist; n = 0 and n = N give the single end member
  2. batch_map.py --skip-existing     every ordering mapped onto each target phase (t, oIII, m, ...);
       variant 1 of equally good orientations unless map_args = "--all-variants"
  3. compare_ordering.py              every mapped structure against its parent: must be IDENTICAL
  4. plan                             all structures (parents + mapped) split into chunks that fit one node
                                      for one walltime (cost model below)
  5. run (one SLURM array task per chunk)   batch_free_energy.py --files-from chunk --relax --T ...
  6. merge                            one CSV for everything + one per composition (plot scripts work on it)

Every step is restart-safe: run it again and it continues where it stopped.

Commands (CONFIG = campaign.toml)
  campaign.py CONFIG dry-run                 counts per composition/phase, chunks, time estimate (no files written)
  campaign.py CONFIG prepare [--comp 0 6 32] [--task I]   steps 1-3 (all compositions, or some, or the I-th)
  campaign.py CONFIG plan [--append|--force]             step 4
  campaign.py CONFIG run --chunk K           step 5 for one chunk (or --task I --offset O on SLURM)
  campaign.py CONFIG status [--json]         progress, failures, flags, measured time per structure (recalibrate)
  campaign.py CONFIG merge                   step 6
  campaign.py CONFIG slurm prepare|compute|merge [--submit] [--unfinished]   write (and submit) SLURM scripts

Cost model for chunking:  seconds(structure) = sec_per_structure * (n_atoms * n_phonon_supercell_atoms) / N0^2,
N0 = atoms in the host POSCAR; i.e. sec_per_structure is the wall time of ONE host-sized structure whose
phonon supercell is 1x1x1, on `threads` cores (relaxation + all displacements).  `status` prints the value
measured from finished structures; put it into the config and re-plan the remaining work (plan --force: finished structures are kept and merged).

Physics caveats: harmonic F at fixed (relaxed) cell, not QHA; rows with imaginary modes or unconverged
relaxations are flagged in the CSV and the merge report.  Equally good mapping variants can only be ranked by
energy (map_args = "--all-variants" computes them all).
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from math import comb
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:                                   # Python < 3.11
    import tomli as tomllib                                   # type: ignore

import numpy as np

HERE = Path(__file__).resolve().parent

DEFAULTS = {
    "structure": {"poscar": None, "host": "Zr", "dopant": "Ce", "compositions": "all",
                  "max_orderings": 1000, "seed": 1},
    "paths": {"scripts": str(HERE), "workdir": "campaign"},
    "mapping": {"templates": None, "map_args": "", "include_parent": True, "verify": True, "jobs": 4},
    "free_energy": {"nep": None, "T": 300.0, "relax": True, "fmax": 1e-3, "steps": 2000, "fixed_cell": False,
                    "min_length": 10.0, "mesh": 100, "threads": 1, "supercell": []},
    "cluster": {"cores_per_node": 64, "walltime": "24:00:00", "fill": 0.8, "sec_per_structure": 300.0,
                "max_array": 1000, "max_concurrent": 100, "prepare_walltime": "04:00:00",
                "merge_walltime": "01:00:00", "partition": "", "account": "", "qos": "", "mem": "",
                "sbatch_extra": [], "setup": "", "python": "python3"},
}


# ============================================================================ config
class Cfg:
    def __init__(self, path):
        self.path = Path(path).resolve()
        with open(self.path, "rb") as f:
            raw = tomllib.load(f)
        unknown = set(raw) - set(DEFAULTS)
        if unknown:
            sys.exit(f"{path}: unknown section(s) {sorted(unknown)}")
        self.c = {}
        for sec, d in DEFAULTS.items():
            got = raw.get(sec, {})
            bad = set(got) - set(d)
            if bad:
                sys.exit(f"{path}: unknown key(s) in [{sec}]: {sorted(bad)}")
            self.c[sec] = {**d, **got}
        base = self.path.parent
        rp = lambda p: (base / p).resolve() if p else None
        s, m, fe = self.c["structure"], self.c["mapping"], self.c["free_energy"]
        if not s["poscar"] or not m["templates"]:
            sys.exit(f"{path}: [structure] poscar and [mapping] templates are required")
        self.poscar = rp(s["poscar"])
        self.scripts = rp(self.c["paths"]["scripts"])
        self.work = rp(self.c["paths"]["workdir"])
        self.nep = rp(fe["nep"]) if fe["nep"] else None
        self.templates = []
        for spec in m["templates"]:
            label, p = spec.split("=", 1) if "=" in spec else (re.sub(r"^POSCAR[_.\-]?", "", Path(spec).name), spec)
            self.templates.append((label, rp(p)))
        for label, p in self.templates:
            if not p.is_file():
                sys.exit(f"template {p} ({label}) not found")
        if not self.poscar.is_file():
            sys.exit(f"host POSCAR {self.poscar} not found")
        for scr in ("enumerate_substitutions.py", "batch_map.py", "map_fluorite.py", "compare_ordering.py",
                    "batch_free_energy.py"):
            if not (self.scripts / scr).is_file():
                sys.exit(f"{scr} not found in scripts folder {self.scripts}")
        sys.path.insert(0, str(self.scripts))
        import enumerate_substitutions as ES                     # noqa: E402  (same folder as the package)
        self.ES = ES
        lat, frac, sym = ES.read_poscar(str(self.poscar))
        self.host_lat, self.n_atoms0 = lat, len(sym)
        self.n_host = int(np.sum(sym == s["host"]))
        if self.n_host == 0:
            sys.exit(f"{self.poscar} has no {s['host']}")
        self.comps = self._parse_comps(s["compositions"])
        self.T_tag = f"{float(fe['T']):g}"

    def __getitem__(self, k):
        return self.c[k]

    def _parse_comps(self, spec):
        N = self.n_host
        if spec == "all":
            out = list(range(N + 1))
        elif isinstance(spec, str):
            m = re.fullmatch(r"\s*(\d+)\s*:\s*(\d+)\s*(?::\s*(\d+))?\s*", spec)
            if not m:
                sys.exit(f"compositions = {spec!r}: use 'all', 'start:stop[:step]' (inclusive) or a list")
            a, b, st = int(m[1]), int(m[2]), int(m[3] or 1)
            out = list(range(a, b + 1, st))
        else:
            out = sorted({int(v) for v in spec})
        bad = [n for n in out if not 0 <= n <= N]
        if bad:
            sys.exit(f"compositions {bad} outside 0..{N} ({self.c['structure']['host']} sites)")
        return out

    # ---- naming helpers (the same rules as enumerate_substitutions.py)
    def comp_dir(self, n):
        return self.work / f"{self.c['structure']['dopant']}_{self.ES.fmt_pct(n / self.n_host)}"

    def target(self, n):
        return min(int(self.c["structure"]["max_orderings"]), comb(self.n_host, n))

    def seed(self, n):
        return int(self.c["structure"]["seed"]) * 1000 + n

    def workers(self):
        cpn = int(os.environ.get("SLURM_CPUS_ON_NODE", 0)) or int(self.c["cluster"]["cores_per_node"])
        return max(1, cpn // int(self.c["free_energy"]["threads"]))

    def py(self):
        return self.c["cluster"]["python"]


def run(cmd, log=None, check=True):
    show = [Path(str(c)).name if os.sep in str(c) else str(c) for c in cmd]     # short echo; full one in the log
    print("  $ " + " ".join(shlex.quote(c) for c in show), flush=True)
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if log:
        with open(log, "a") as f:
            f.write("$ " + " ".join(map(str, cmd)) + "\n" + r.stdout + r.stderr + "\n")
    if check and r.returncode != 0:
        tail = "\n".join((r.stdout + r.stderr).strip().split("\n")[-15:])
        sys.exit(f"command failed (exit {r.returncode}):\n{tail}")
    return r


def hms(s):
    s = int(round(s))
    return f"{s // 3600:d}:{s % 3600 // 60:02d}:{s % 60:02d}"


def parse_walltime(w):
    parts = [int(x) for x in re.split(r"[-:]", w)]
    if "-" in w:                                               # D-HH:MM:SS
        d, rest = parts[0], parts[1:] + [0] * (3 - len(parts[1:]))
        return d * 86400 + rest[0] * 3600 + rest[1] * 60 + rest[2]
    while len(parts) < 3:
        parts = [0] + parts
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


# ============================================================================ cost model
def read_cell(path):
    lines = Path(path).read_text().split("\n", 8)
    sc = float(lines[1].split()[0])
    lat = np.array([[float(v) for v in lines[k].split()[:3]] for k in (2, 3, 4)])
    lat = lat * ((abs(sc) / abs(np.linalg.det(lat))) ** (1 / 3) if sc < 0 else sc)
    n = sum(int(v) for v in lines[6].split())
    return lat, n


def phonon_sc(cfg, lat):
    fe = cfg["free_energy"]
    if fe["supercell"]:
        return [int(v) for v in fe["supercell"]]
    widths = 1.0 / np.linalg.norm(np.linalg.inv(lat), axis=0)
    return [max(1, int(np.ceil(float(fe["min_length"]) / w - 1e-6))) for w in widths]


def cost_units(cfg, lat, n_atoms):
    """(n_atoms * phonon-supercell atoms) / N0^2  ->  seconds = units * sec_per_structure"""
    sc = phonon_sc(cfg, lat)
    return n_atoms * n_atoms * int(np.prod(sc)) / cfg.n_atoms0 ** 2, sc


def predicted_cells(cfg):
    """cell of each phase for the host's repeat, before anything is mapped (template axes x repeats)"""
    import map_fluorite as M
    with contextlib.redirect_stderr(io.StringIO()):
        cell, frac, sp = M.read_poscar(str(cfg.poscar))
        try:
            nrep = np.array(M.detect_n(cell, frac, sp))
        except SystemExit:
            nrep = None
    out = {"parent": cfg.host_lat}
    for label, p in cfg.templates:
        tcell = M.read_poscar(str(p))[0]
        out[label] = tcell * (nrep[:, None] if nrep is not None else 1)
    return out, nrep


# ============================================================================ prepare
def prepare_one(cfg, n):
    s, m = cfg["structure"], cfg["mapping"]
    d = cfg.comp_dir(n)
    cfg.work.mkdir(parents=True, exist_ok=True)
    log = cfg.work / "logs" / f"prepare_{d.name}.log"
    log.parent.mkdir(exist_ok=True)
    want = cfg.target(n)
    old = d / "prepared.json"
    if old.is_file():
        j = json.loads(old.read_text())
        tpl_now = [f"{label}={p}" for label, p in cfg.templates]
        if j.get("templates") != tpl_now or j.get("map_args") != m["map_args"]:
            sys.exit(f"{d} was prepared with templates {j.get('templates')} / map_args {j.get('map_args')!r}; "
                     "the config now differs. Use a new workdir (or delete the folder) to change them.")
    print(f"\n=== {d.name}: n = {n} of {cfg.n_host} {s['host']} -> {s['dopant']}; {want} orderings "
          f"(of C({cfg.n_host},{n}) = {comb(cfg.n_host, n):,}), seed {cfg.seed(n)}")

    # 1. orderings (enumerate --random draws only the missing ones and never repeats)
    have = count_parents(d)
    if have < want:
        r = run([cfg.py(), cfg.scripts / "enumerate_substitutions.py", cfg.poscar, "--host", s["host"],
                 "--dopant", s["dopant"], "--n", n, "--random", want - have, "--seed", cfg.seed(n),
                 "--outdir", cfg.work], log)
        print("  " + "\n  ".join(l for l in r.stdout.split("\n") if l.startswith(("wrote", "WARNING", "only"))))
    have = count_parents(d)
    if have != want:
        sys.exit(f"{d}: {have} orderings, expected {want}")
    print(f"  orderings: {have}")

    # 2. mapping (restart-safe)
    tpl = [f"{label}={p}" for label, p in cfg.templates]
    cmd = [cfg.py(), cfg.scripts / "batch_map.py", d, "--templates", *tpl, "--glob", "POSCAR_*",
           "--jobs", m["jobs"], "--skip-existing"]
    if m["map_args"]:
        cmd += ["--args", m["map_args"]]
    r = run(cmd, log, check=False)
    rows = read_csv(d / "mapped" / "summary.csv")
    bad = [x for x in rows if x["status"] != "ok"]
    if r.returncode != 0 or bad or len(rows) != have * len(cfg.templates):
        sys.exit(f"{d}: mapping failed for {len(bad)} of {len(rows)} runs (see {log})\n"
                 + "\n".join(f"  {x['file']} -> {x['target']}: {x['note']}" for x in bad[:10]))
    report_mapping(rows)

    # 3. verification with compare_ordering.py
    if m["verify"]:
        verify_comp(cfg, d, rows)
    (d / "prepared.json").write_text(json.dumps(
        {"n": n, "orderings": have, "seed": cfg.seed(n), "templates": tpl, "map_args": m["map_args"],
         "time": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
    print(f"  {d.name} prepared")


def count_parents(d):
    p = d / "summary.csv"
    if not p.is_file():
        return 0
    rows = read_csv(p)
    return sum(1 for r in rows if (d / r["file"]).is_file())


def read_csv(p):
    p = Path(p)
    if not p.is_file():
        return []
    with p.open() as f:
        return list(csv.DictReader(f))


def report_mapping(rows):
    by = {}
    for r in rows:
        by.setdefault(r["target"], []).append(r)
    for ph, rr in by.items():
        sh = np.array([float(r["rms_shift_A"] or "nan") for r in rr])
        st = np.array([float(r["strain_S"] or "nan") for r in rr])
        tv = np.array([int(r["tied_variants"] or 1) for r in rr])
        print(f"  mapped {ph:<6}: {len(rr):5d} ok   rms shift {np.nanmin(sh):.3f}-{np.nanmax(sh):.3f} A,  "
              f"strain S {np.nanmin(st):.4f}-{np.nanmax(st):.4f},  equally good variants {tv.min()}-{tv.max()} "
              f"(variant 1 written; rank them only by energy)")


def verify_one(args):
    """compare_ordering.py A B, in-process; returns the verdict word"""
    scripts, a, b = args
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    import compare_ordering as CO
    out, err = io.StringIO(), io.StringIO()
    argv = sys.argv
    sys.argv = ["compare_ordering.py", a, b]
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            CO.main()
    except SystemExit:
        pass
    finally:
        sys.argv = argv
    mm = re.search(r"VERDICT: (\w+)", out.getvalue())
    return mm.group(1) if mm else ("ERROR " + err.getvalue().strip()[-200:])


def verify_comp(cfg, d, rows):
    from concurrent.futures import ProcessPoolExecutor
    vpath = d / "mapped" / "verify.csv"
    done = {(r["file"], r["out_file"]): r["verdict"] for r in read_csv(vpath)}
    todo = []
    for r in rows:
        outs = [r["out_file"]]
        outs += sorted(p.name for p in (d / "mapped" / r["target"]).glob(r["out_file"] + "_v*"))   # --all-variants
        for o in outs:
            if (r["file"], o) not in done:
                todo.append((r["file"], r["target"], o))
    if todo:
        jobs = [(str(cfg.scripts), str(d / f), str(d / "mapped" / t / o)) for f, t, o in todo]
        with ProcessPoolExecutor(max_workers=max(1, int(cfg["mapping"]["jobs"]))) as ex:
            verdicts = list(ex.map(verify_one, jobs, chunksize=16))
        new = not vpath.exists()
        with vpath.open("a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["file", "target", "out_file", "verdict"])
            for (f, t, o), v in zip(todo, verdicts):
                w.writerow([f, t, o, v])
                done[(f, o)] = v
    bad = {k: v for k, v in done.items() if v != "IDENTICAL"}
    print(f"  compare_ordering: {len(done) - len(bad)} of {len(done)} mapped structures IDENTICAL to their parent")
    if bad:
        sys.exit(f"{d}: {len(bad)} mapped structures are NOT identical to their parent ordering, e.g. "
                 + "; ".join(f"{a} -> {b}: {v}" for (a, b), v in list(bad.items())[:5]))


# ============================================================================ structures / plan
def structures(cfg, n):
    """[(relative path to workdir, phase)] for one prepared composition, in a stable order"""
    d = cfg.comp_dir(n)
    out = []
    if cfg["mapping"]["include_parent"]:
        for r in read_csv(d / "summary.csv"):
            out.append((f"{d.name}/{r['file']}", "parent"))
    all_var = "--all-variants" in cfg["mapping"]["map_args"]
    for label, _ in cfg.templates:
        for r in read_csv(d / "mapped" / "summary.csv"):
            if r["target"] != label or r["status"] != "ok":
                continue
            vs = sorted((d / "mapped" / label).glob(r["out_file"] + "_v*")) if all_var else []
            if vs:                                             # OUT == OUT_v1: compute the _v files only
                out += [(f"{d.name}/mapped/{label}/{v.name}", label) for v in vs]
            else:
                out.append((f"{d.name}/mapped/{label}/{r['out_file']}", label))
    return out


def prepared(cfg, n):
    """composition fully prepared for the CURRENT config (orderings, templates, map options)"""
    p = cfg.comp_dir(n) / "prepared.json"
    if not p.is_file():
        return False
    j = json.loads(p.read_text())
    return (j.get("orderings") == cfg.target(n) and j.get("map_args") == cfg["mapping"]["map_args"]
            and j.get("templates") == [f"{label}={p}" for label, p in cfg.templates])


def cmd_plan(cfg, a):
    cdir = cfg.work / "chunks"
    plan_p = cdir / "plan.json"
    missing = [n for n in cfg.comps if not prepared(cfg, n)]
    if missing and not a.partial:
        sys.exit(f"compositions not prepared yet: {missing}  (run prepare, or plan --partial)")
    old = json.loads(plan_p.read_text()) if plan_p.exists() else None
    if old and not (a.append or a.force):
        sys.exit(f"{plan_p} exists. --append adds structures not yet planned as new chunks; --force re-plans "
                 "everything (only if no compute job is running).")
    if old and a.force:
        for f in cdir.glob("chunk_*.txt"):
            f.unlink()
        old = None
    planned = set()
    chunks = old["chunks"] if old else []
    for ch in chunks:
        planned.update(Path(cdir / ch["file"]).read_text().split())
    done_ok = {r["file"] for r in all_results(cfg) if r["status"] == "ok"}

    todo = []
    for n in cfg.comps:
        if not prepared(cfg, n):
            continue
        for rel, ph in structures(cfg, n):
            if rel in planned or rel in done_ok:
                continue
            lat, nat = read_cell(cfg.work / rel)
            u, _ = cost_units(cfg, lat, nat)
            todo.append((rel, u))
    cl = cfg["cluster"]
    workers = int(cl["cores_per_node"]) // int(cfg["free_energy"]["threads"])
    cap = workers * parse_walltime(cl["walltime"]) * float(cl["fill"]) / float(cl["sec_per_structure"])  # units
    new, cur, cur_u = [], [], 0.0
    for rel, u in todo:
        if cur and cur_u + u > cap:
            new.append((cur, cur_u)); cur, cur_u = [], 0.0
        cur.append(rel); cur_u += u
    if cur:
        new.append((cur, cur_u))
    cdir.mkdir(parents=True, exist_ok=True)
    k0 = len(chunks)
    for i, (files, u) in enumerate(new, k0):
        name = f"chunk_{i:05d}.txt"
        (cdir / name).write_text("\n".join(files) + "\n")
        chunks.append({"id": i, "file": name, "n": len(files), "units": round(u, 3),
                       "est_s": round(u * float(cl["sec_per_structure"]) / workers, 1)})
    plan = {"workers_per_node": workers, "walltime": cl["walltime"], "sec_per_structure": cl["sec_per_structure"],
            "chunks": chunks}
    plan_p.write_text(json.dumps(plan, indent=1))
    biggest = max((c["est_s"] for c in chunks[k0:]), default=0)
    print(f"{len(todo)} structures -> {len(new)} new chunks (ids {k0}-{len(chunks) - 1}), "
          f"{workers} parallel structures per node, largest chunk ~{hms(biggest)} of walltime {cl['walltime']}")
    if any(u > cap for _, u in todo):
        print("WARNING: single structures predicted to take longer than the walltime; raise walltime or threads")
    print(f"plan: {plan_p}")


def load_plan(cfg):
    p = cfg.work / "chunks" / "plan.json"
    if not p.exists():
        sys.exit("no plan yet: run  campaign.py CONFIG plan")
    return json.loads(p.read_text())


def chunk_csv(cfg, k):
    return cfg.work / "results" / f"chunk_{k:05d}.csv"


def all_results(cfg):
    rows = []
    for p in sorted((cfg.work / "results").glob("chunk_*.csv")):
        rows += read_csv(p)
    return rows


# ============================================================================ run
def cmd_run(cfg, a):
    plan = load_plan(cfg)
    k = a.chunk
    if k is None:
        tid = a.task if a.task is not None else os.environ.get("SLURM_ARRAY_TASK_ID")
        if tid is None:
            sys.exit("give --chunk K (or --task / SLURM_ARRAY_TASK_ID)")
        k = int(tid) + int(a.offset if a.offset is not None else os.environ.get("CHUNK_OFFSET", 0))
    if not 0 <= k < len(plan["chunks"]):
        sys.exit(f"chunk {k} not in plan (0-{len(plan['chunks']) - 1})")
    if not cfg.nep or not cfg.nep.is_file():
        sys.exit(f"NEP model {cfg.nep} not found ([free_energy] nep)")
    fe = cfg["free_energy"]
    out = chunk_csv(cfg, k)
    out.parent.mkdir(exist_ok=True)
    cmd = [sys.executable, cfg.scripts / "batch_free_energy.py", cfg.work,
           "--files-from", cfg.work / "chunks" / plan["chunks"][k]["file"], "--out", out,
           "--nep", cfg.nep, "--T", fe["T"], "--fmax", fe["fmax"], "--steps", fe["steps"],
           "--min-length", fe["min_length"], "--mesh", fe["mesh"],
           "--jobs", cfg.workers(), "--threads", fe["threads"]]
    if fe["relax"]:
        cmd.append("--relax")
    if fe["fixed_cell"]:
        cmd.append("--fixed-cell")
    if fe["supercell"]:
        cmd += ["--supercell", *fe["supercell"]]
    print(f"chunk {k}: {plan['chunks'][k]['n']} structures, {cfg.workers()} in parallel", flush=True)
    print("$ " + " ".join(map(str, cmd)), flush=True)
    sys.exit(subprocess.call([str(c) for c in cmd]))


# ============================================================================ status / merge
def chunk_progress(cfg, plan):
    out = []
    for ch in plan["chunks"]:
        files = (cfg.work / "chunks" / ch["file"]).read_text().split()
        rows = {r["file"]: r for r in read_csv(chunk_csv(cfg, ch["id"]))}
        ok = sum(1 for f in files if rows.get(f, {}).get("status") == "ok")
        failed = sum(1 for f in files if rows.get(f, {}).get("status") == "FAILED")
        out.append((ch["id"], len(files), ok, failed))
    return out


def _num(v, default=0.0):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return default
    return x if np.isfinite(x) else default


def row_flags(r):
    """(imaginary modes, unconverged relaxation) of an ok row; tolerant of empty/odd values in old CSVs"""
    return _num(r.get("n_imaginary_modes")) > 0, str(r.get("relax_converged", "")).strip() == "False"


def status_summary(cfg):
    """everything `status` reports, as a dict (also printed by status --json).  Counts use the same rule as merge
    (latest row per structure, ok beats FAILED, rerun rows after campaign rows), so status, merge and the
    Studio's Results page agree."""
    plan = load_plan(cfg)
    best, n_rows, temps, campaign_ok = {}, 0, set(), []
    for p in sorted((cfg.work / "results").glob("chunk_*.csv")):
        is_rerun = p.name.startswith("chunk_rerun")
        for r in read_csv(p):
            if not r.get("file"):
                continue
            n_rows += 1
            if r.get("T_K"):
                temps.add(f"{_num(r['T_K'], float('nan')):g}")
            f = r["file"]
            if f not in best or r.get("status") == "ok" or best[f].get("status") != "ok":
                best[f] = r
            if r.get("status") == "ok" and not is_rerun:
                campaign_ok.append(r)
    # ---- chunks (progress from each chunk's own CSV: that is what `run` / `slurm compute --unfinished` redo)
    planned, chunks, missing_lists = set(), [], []
    for ch in plan["chunks"]:
        lp = cfg.work / "chunks" / ch["file"]
        if not lp.is_file():
            missing_lists.append(ch["file"])
            continue
        files = lp.read_text().split()
        planned.update(files)
        own = {}
        for r in read_csv(chunk_csv(cfg, ch["id"])):
            own[r.get("file")] = r.get("status")
        ok = sum(1 for f in files if own.get(f) == "ok")
        failed = sum(1 for f in files if own.get(f) == "FAILED")
        chunks.append({"id": ch["id"], "n": len(files), "ok": ok, "failed": failed, "not_run": len(files) - ok - failed})
    unf = [c for c in chunks if c["ok"] < c["n"]]
    kinds = {"not_started": [c["id"] for c in unf if c["ok"] + c["failed"] == 0],
             "partly_run": [c["id"] for c in unf if c["not_run"] and c["ok"] + c["failed"]],
             "only_failed_left": [c["id"] for c in unf if not c["not_run"]]}
    # ---- structures
    allf = planned | set(best)
    ok_rows = [best[f] for f in allf if f in best and best[f].get("status") == "ok"]
    failed_f = sorted(f for f in allf if f in best and best[f].get("status") != "ok")
    not_run = sorted(f for f in planned if f not in best)
    earlier = sorted(f for f in best if f not in planned)
    unplanned = []                                         # prepared, but in no plan and without a result
    for n in cfg.comps:
        if prepared(cfg, n):
            unplanned += [rel for rel, _ in structures(cfg, n) if rel not in allf]
    notes = {}
    for f in failed_f:
        k = (best[f].get("note") or "no note").strip()[:120]
        notes[k] = notes.get(k, 0) + 1
    imag = sum(1 for r in ok_rows if row_flags(r)[0])
    unc = sum(1 for r in ok_rows if row_flags(r)[1])
    both = sum(1 for r in ok_rows if all(row_flags(r)))
    rounds = sorted(int(p.name.split("_")[1]) for p in cfg.work.glob("rerun_[0-9][0-9][0-9]") if p.is_dir())
    # ---- measured cost (campaign rows only: rerun rows use stricter settings; one row per structure)
    meas, seen = [], set()
    for r in reversed(campaign_ok):                        # newest row of a structure first
        if r["file"] in seen or _num(r.get("time_s"), -1) <= 0:
            continue
        seen.add(r["file"])
        try:
            lat, nat = read_cell(cfg.work / r["file"])
        except (OSError, ValueError, IndexError):
            continue
        meas.append(_num(r["time_s"]) / cost_units(cfg, lat, nat)[0])
    cal = None
    if meas:
        m = np.array(meas)
        cfg_s = float(cfg["cluster"]["sec_per_structure"])
        med = float(np.median(m))
        cal = {"median": round(med, 2), "p10": round(float(np.percentile(m, 10)), 2),
               "p90": round(float(np.percentile(m, 90)), 2), "n": len(m), "config": cfg_s,
               "recalibrate": bool((not_run or unplanned) and not 2 / 3 <= med / cfg_s <= 1.5)}
    return {"structures": len(allf) + len(unplanned), "ok": len(ok_rows), "failed": len(failed_f),
            "not_run": len(not_run), "not_planned": len(unplanned), "planned": len(planned),
            "earlier_plan": len(earlier), "result_rows": n_rows,
            "chunks": {"n": len(chunks), "complete": len(chunks) - len(unf), "unfinished": [c["id"] for c in unf],
                       **kinds, "missing_lists": missing_lists, "list": chunks},
            "flags": {"imaginary": imag, "unconverged": unc, "both": both, "clean": len(ok_rows) - imag - unc + both},
            "failed_notes": sorted(notes.items(), key=lambda kv: -kv[1]), "rerun_rounds": rounds,
            "temperatures": sorted(temps, key=float), "sec_per_structure": cal,
            "workers_per_node_planned": plan.get("workers_per_node"), "T_config": cfg.T_tag}


def cmd_status(cfg, a):
    s = status_summary(cfg)
    if a.json:
        print(json.dumps(s, indent=1))
        return
    c = s["chunks"]
    print(f"{s['structures']} structures: {s['ok']} ok, {s['failed']} failed, {s['not_run']} not run"
          + (f", {s['not_planned']} not planned yet (plan --append)" if s["not_planned"] else "")
          + "   (same counts as merge)")
    print(f"plan: {c['n']} chunk(s), {c['complete']} complete, {len(c['unfinished'])} unfinished; "
          f"{s['planned']} structures planned"
          + (f", {s['earlier_plan']} more finished under an earlier plan (kept, merged)" if s["earlier_plan"] else ""))
    if c["unfinished"]:
        parts = [f"{lab} {compress(c[k])}" for k, lab in (("not_started", "not started"), ("partly_run", "partly run"),
                                                          ("only_failed_left", "only FAILED left")) if c[k]]
        print("unfinished chunk ids: " + "; ".join(parts) + "   (slurm compute --unfinished reruns them, "
              "FAILED ones are retried)")
    if c["missing_lists"]:
        print(f"WARNING: {len(c['missing_lists'])} chunk lists in plan.json are missing (e.g. {c['missing_lists'][0]}); "
              "re-plan with plan --force")
    for note, k in s["failed_notes"][:3]:
        print(f"FAILED {k}x: {note}")
    if len(s["failed_notes"]) > 3:
        print(f"FAILED: {len(s['failed_notes']) - 3} more kinds of error (see the note column)")
    fl = s["flags"]
    if s["ok"]:
        rr = s["rerun_rounds"]
        print(f"flags (latest row per structure{', after rerun round ' + compress(rr) if rr else ''}): "
              f"{fl['imaginary']} with imaginary modes, {fl['unconverged']} with unconverged relaxation, "
              f"{fl['clean']} clean" + ("   (rerun_flagged.py recomputes flagged ones)"
                                        if fl["imaginary"] + fl["unconverged"] and not rr else ""))
    if len(s["temperatures"]) > 1:
        print(f"WARNING: results at several temperatures {s['temperatures']} K; merge mixes them "
              f"(config T = {s['T_config']} K)")
    cal = s["sec_per_structure"]
    if cal:
        print(f"measured sec_per_structure: {cal['median']:.3g} s (median of {cal['n']} campaign structures, "
              f"10-90 %: {cal['p10']:.3g}-{cal['p90']:.3g}); config has {cal['config']:g}")
        if cal["recalibrate"]:
            print(f"  -> set sec_per_structure = {cal['median']:.3g} in [cluster] and run plan --force before "
                  "computing the rest (finished results are kept)")


def compress(ids):
    ids = sorted(ids)
    parts, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        parts.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ",".join(parts)


def cmd_merge(cfg, a):
    import batch_free_energy as BFE
    plan = load_plan(cfg)
    planned = []
    for ch in plan["chunks"]:
        planned += (cfg.work / "chunks" / ch["file"]).read_text().split()
    best = {}
    for r in all_results(cfg):                                 # ok beats FAILED; later beats earlier
        f = r["file"]
        if f not in best or r["status"] == "ok" or best[f]["status"] != "ok":
            best[f] = r
    planned_set = set(planned)
    # structures finished under an earlier plan (e.g. before `plan --force`) are not in the current chunks: keep them
    earlier = sorted(f for f in best if f not in planned_set)
    rows = [best[f] for f in planned if f in best] + [best[f] for f in earlier]
    missing = [f for f in planned if f not in best]
    if earlier:
        print(f"{len(earlier)} structures from earlier plans included "
              f"({sum(1 for f in earlier if best[f]['status'] == 'ok')} ok)")
    lowest = {}
    for r in rows:
        if r["status"] == "ok":
            F = float(r["F_total_eV_atom"])
            lowest[r["composition"]] = min(lowest.get(r["composition"], F), F)
    for r in rows:
        r["dF_total_meV_atom"] = (f"{1000 * (float(r['F_total_eV_atom']) - lowest[r['composition']]):.4f}"
                                  if r["status"] == "ok" else "")
        r["comp_dir"] = r["file"].split("/")[0]
    by_dir = {cfg.comp_dir(n).name: n for n in cfg.comps}
    for r in rows:
        n = by_dir.get(r["comp_dir"])
        r["n_dopant"] = n if n is not None else ""
        r["x_dopant"] = f"{n / cfg.n_host:.6f}" if n is not None else ""
    key = lambda r: (float(r["x_dopant"] or 0), r["status"] != "ok",
                     float(r["F_total_eV_atom"]) if r["status"] == "ok" else 0.0)
    rows.sort(key=key)
    cols = ["n_dopant", "x_dopant"] + BFE.CSV_COLS
    allp = cfg.work / f"free_energy_{cfg.T_tag}K_all.csv"
    with allp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows({k: r.get(k, "") or "" for k in cols} for r in rows)
    # one CSV per composition, file paths relative to the composition folder (the plot scripts read these)
    per = {}
    for r in rows:
        per.setdefault(r["comp_dir"], []).append(r)
    for dname, rr in per.items():
        with (cfg.work / dname / f"free_energy_{cfg.T_tag}K.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=BFE.CSV_COLS, extrasaction="ignore")
            w.writeheader()
            for r in rr:
                w.writerow({**{k: r.get(k, "") or "" for k in BFE.CSV_COLS}, "file": r["file"].split("/", 1)[1]})
    # report
    rep = io.StringIO()
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"{len(rows) + len(missing)} structures ({len(planned)} in the current plan), {len(ok)} ok, {len(rows) - len(ok)} failed, {len(missing)} not run", file=rep)
    print(f"{'comp':<10} {'phase':<7} {'ok':>5} {'fail':>5} {'imag':>5} {'unconv':>6}  {'lowest dF':>9}  "
          f"{'lowest dF, unflagged':>20}   (meV/atom vs lowest F of the composition)", file=rep)
    phases = ["parent"] * bool(cfg["mapping"]["include_parent"]) + [l for l, _ in cfg.templates]
    for n in cfg.comps:
        dname = cfg.comp_dir(n).name
        for ph in phases:
            rr = [r for r in per.get(dname, []) if r["phase"] == ph]
            if not rr:
                continue
            o = [r for r in rr if r["status"] == "ok"]
            im = sum(1 for r in o if row_flags(r)[0])
            un = sum(1 for r in o if row_flags(r)[1])
            lo = min((float(r["dF_total_meV_atom"]) for r in o), default=float("nan"))
            lc = min((float(r["dF_total_meV_atom"]) for r in o if not any(row_flags(r))), default=float("nan"))
            print(f"{dname:<10} {ph:<7} {len(o):5d} {len(rr) - len(o):5d} {im:5d} {un:6d}  {lo:9.2f}  {lc:20.2f}",
                  file=rep)
    if missing:
        print(f"\nnot run yet ({len(missing)}), e.g. {missing[:3]}", file=rep)
    print(f"\nharmonic F at {cfg.T_tag} K, fixed relaxed cell (not QHA). Rows with imaginary modes or "
          "unconverged relaxation are kept but flagged (n_imaginary_modes, relax_converged, note).", file=rep)
    (cfg.work / f"merge_report_{cfg.T_tag}K.txt").write_text(rep.getvalue())
    print(rep.getvalue())
    print(f"wrote {allp}\n      {len(per)} per-composition CSVs <comp>/free_energy_{cfg.T_tag}K.csv\n"
          f"      {cfg.work / f'merge_report_{cfg.T_tag}K.txt'}")


# ============================================================================ dry run
def cmd_dry(cfg, a):
    cl, fe = cfg["cluster"], cfg["free_energy"]
    cells, nrep = predicted_cells(cfg)
    phases = (["parent"] if cfg["mapping"]["include_parent"] else []) + [l for l, _ in cfg.templates]
    units, scs = {}, {}
    for ph in phases:
        units[ph], scs[ph] = cost_units(cfg, cells[ph], cfg.n_atoms0)
    print(f"host {cfg.poscar.name}: {cfg.n_atoms0} atoms, {cfg.n_host} {cfg['structure']['host']} sites"
          + (f", repeat {'x'.join(map(str, nrep))} of the 12-atom cell" if nrep is not None else ""))
    print(f"orderings: up to {cfg['structure']['max_orderings']} distinct random per composition "
          f"(all of them where fewer exist), seed {cfg['structure']['seed']}*1000+n")
    print("phonon supercell / relative cost per phase: "
          + ", ".join(f"{ph} {'x'.join(map(str, scs[ph]))} ({units[ph]:.2g})" for ph in phases))
    nv = "--all-variants" in cfg["mapping"]["map_args"]
    if nv:
        print("NOTE: --all-variants: mapped counts below are per variant 1 only; each ordering adds all its "
              "equally good variants (often 6-48 per phase) -> real counts are known after prepare")
    print(f"\n{'n':>3} {'folder':<10} {'C(N,n)':>14} {'orderings':>9} " + " ".join(f"{p:>7}" for p in phases)
          + f" {'total':>7} {'status':>9}")
    tot_s, tot_u, grand = {p: 0 for p in phases}, 0.0, 0
    for n in cfg.comps:
        d = cfg.comp_dir(n)
        if prepared(cfg, n):
            st = structures(cfg, n)
            cnt = {p: sum(1 for _, ph in st if ph == p) for p in phases}
            u = sum(cost_units(cfg, *read_cell(cfg.work / rel))[0] for rel, _ in st)
            status = "prepared"
        else:
            k = cfg.target(n)
            cnt = {p: k for p in phases}
            u = sum(k * units[p] for p in phases)
            status = "-"
        for p in phases:
            tot_s[p] += cnt[p]
        tot_u += u
        grand += sum(cnt.values())
        print(f"{n:3d} {d.name:<10} {comb(cfg.n_host, n):>14,} {cfg.target(n):>9} "
              + " ".join(f"{cnt[p]:7d}" for p in phases) + f" {sum(cnt.values()):7d} {status:>9}")
    print(f"{'':3} {'TOTAL':<10} {'':>14} {'':>9} " + " ".join(f"{tot_s[p]:7d}" for p in phases) + f" {grand:7d}")
    sec = float(cl["sec_per_structure"])
    workers = int(cl["cores_per_node"]) // int(fe["threads"])
    wall = parse_walltime(cl["walltime"])
    cap = workers * wall * float(cl["fill"]) / sec
    core_h = tot_u * sec * int(fe["threads"]) / 3600
    n_chunks = int(np.ceil(tot_u / cap)) if cap > 0 else 0
    print(f"\nestimate with sec_per_structure = {sec:g} s (host-sized structure, 1x1x1 phonon cell, "
          f"{fe['threads']} thread(s); measure it: run one chunk, then `status`):")
    print(f"  {core_h:,.1f} core-hours = {core_h / int(cl['cores_per_node']):,.0f} node-hours "
          f"({cl['cores_per_node']} cores/node)")
    print(f"  ~{n_chunks} chunks of <= {cl['walltime']} on one node each ({workers} structures in parallel, "
          f"{float(cl['fill']):.0%} of the walltime planned)")
    if n_chunks:
        per_chunk = tot_u * sec / workers / n_chunks
        rounds = int(np.ceil(n_chunks / max(1, int(cl["max_concurrent"]))))
        print(f"  each chunk ~{hms(per_chunk)}; with up to {cl['max_concurrent']} nodes at a time: ~{hms(rounds * per_chunk)} "
              "wall (+ queueing). A shorter walltime gives more, smaller chunks.")
    print(f"  mapping + verification: {grand - tot_s.get('parent', 0)} map_fluorite runs, ~1 core-s each")


# ============================================================================ SLURM
def sbatch_header(cfg, name, time_, cpus, log):
    cl = cfg["cluster"]
    h = ["#!/bin/bash", f"#SBATCH --job-name={name}", "#SBATCH --nodes=1", "#SBATCH --ntasks=1",
         f"#SBATCH --cpus-per-task={cpus}", f"#SBATCH --time={time_}", f"#SBATCH --output={log}"]
    for k in ("partition", "account", "qos", "mem"):
        if cl[k] != "":
            h.append(f"#SBATCH --{k}={cl[k]}")
    h += [f"#SBATCH {x}" for x in cl["sbatch_extra"]]
    h += ["set -euo pipefail", "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"]
    if cl["setup"]:
        h += ["# --- environment (config [cluster] setup)", cl["setup"].strip()]
    return "\n".join(h) + "\n"


def array_blocks(ids, max_array):
    """split chunk ids into SLURM arrays whose indices stay below MaxArraySize: [(offset, 'i,j-k')]"""
    out = {}
    for i in sorted(ids):
        out.setdefault(i // max_array * max_array, []).append(i)
    return [(off, compress([i - off for i in v])) for off, v in sorted(out.items())]


def cmd_slurm(cfg, a):
    cl = cfg["cluster"]
    sd = cfg.work / "slurm"
    (cfg.work / "logs").mkdir(parents=True, exist_ok=True)
    sd.mkdir(parents=True, exist_ok=True)
    me = f"{cfg.py()} {shlex.quote(str(Path(__file__).resolve()))} {shlex.quote(str(cfg.path))}"
    cmds = []
    if a.stage == "prepare":
        comps = [n for n in cfg.comps if not (a.unfinished and prepared(cfg, n))]
        if not comps:
            print("all compositions prepared"); return
        idx = [cfg.comps.index(n) for n in comps]
        sh = sd / "prepare.sh"
        sh.write_text(sbatch_header(cfg, "fo_prepare", cl["prepare_walltime"], cfg["mapping"]["jobs"],
                                    cfg.work / "logs" / "prepare_%a.out")
                      + f'{me} prepare --task "$SLURM_ARRAY_TASK_ID"\n')
        cmds.append(["sbatch", "--parsable", f"--array={compress(idx)}%{cl['max_concurrent']}", sh])
    elif a.stage == "compute":
        plan = load_plan(cfg)
        ids = [c["id"] for c in plan["chunks"]]
        if a.unfinished:
            ids = [p[0] for p in chunk_progress(cfg, plan) if p[2] < p[1]]
        if not ids:
            print("nothing to compute"); return
        sh = sd / "compute.sh"
        sh.write_text(sbatch_header(cfg, "fo_free_energy", cl["walltime"], cl["cores_per_node"],
                                    cfg.work / "logs" / "chunk_%A_%a.out")
                      + f'{me} run --task "$SLURM_ARRAY_TASK_ID" --offset "${{CHUNK_OFFSET:-0}}"\n')
        for off, spec in array_blocks(ids, int(cl["max_array"])):
            cmds.append(["sbatch", "--parsable", f"--array={spec}%{cl['max_concurrent']}",
                         f"--export=ALL,CHUNK_OFFSET={off}", sh])
    else:
        sh = sd / "merge.sh"
        sh.write_text(sbatch_header(cfg, "fo_merge", cl["merge_walltime"], 1, cfg.work / "logs" / "merge_%j.out")
                      + f"{me} status || true\n{me} merge\n")
        dep = [f"--dependency=afterany:{':'.join(a.after)}"] if a.after else []
        cmds.append(["sbatch", "--parsable", *dep, sh])
    print(f"wrote {sh}")
    ids = []
    for c in cmds:
        if a.submit:
            r = run(c)
            ids.append(r.stdout.strip().split(";")[0])
            print(f"  submitted job {ids[-1]}")
        else:
            print("  " + " ".join(map(str, c)))
    if a.submit and a.stage == "compute" and a.with_merge:
        dep = f"--dependency=afterany:{':'.join(ids)}"
        sh = sd / "merge.sh"
        sh.write_text(sbatch_header(cfg, "fo_merge", cl["merge_walltime"], 1, cfg.work / "logs" / "merge_%j.out")
                      + f"{me} status || true\n{me} merge\n")
        r = run(["sbatch", "--parsable", dep, sh])
        print(f"  submitted merge job {r.stdout.strip()} (after {', '.join(ids)})")
    if not a.submit:
        print("(not submitted: add --submit)")


# ============================================================================ main
def cmd_prepare(cfg, a):
    if a.task is not None:
        if not 0 <= a.task < len(cfg.comps):
            sys.exit(f"--task {a.task}: the config has {len(cfg.comps)} compositions (0-{len(cfg.comps) - 1})")
        comps = [cfg.comps[a.task]]
    elif a.comp is not None:
        bad = [n for n in a.comp if n not in cfg.comps]
        if bad:
            sys.exit(f"--comp {bad} not in the config's compositions")
        comps = a.comp
    else:
        comps = cfg.comps
    cfg.work.mkdir(parents=True, exist_ok=True)
    (cfg.work / "campaign_config.toml").write_bytes(cfg.path.read_bytes())    # record what was used
    for n in comps:
        prepare_one(cfg, n)


def main():
    if hasattr(signal, "SIGPIPE"):                         # `campaign.py ... status | head` must not crash
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("config")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("dry-run")
    p = sub.add_parser("prepare")
    p.add_argument("--comp", type=int, nargs="+", help="dopant counts (default: all in the config)")
    p.add_argument("--task", type=int, help="index into the config's composition list (SLURM array)")
    p = sub.add_parser("plan")
    p.add_argument("--append", action="store_true", help="add structures that are not planned/done yet")
    p.add_argument("--force", action="store_true", help="discard the old plan (results are kept)")
    p.add_argument("--partial", action="store_true", help="plan only the compositions prepared so far")
    p = sub.add_parser("run")
    p.add_argument("--chunk", type=int)
    p.add_argument("--task", type=int)
    p.add_argument("--offset", type=int)
    p = sub.add_parser("status")
    p.add_argument("--json", action="store_true", help="print the status as JSON (for scripts and the Studio)")
    sub.add_parser("merge")
    p = sub.add_parser("slurm")
    p.add_argument("stage", choices=["prepare", "compute", "merge"])
    p.add_argument("--submit", action="store_true", help="call sbatch (default: only write and print)")
    p.add_argument("--unfinished", action="store_true", help="only compositions/chunks not finished yet")
    p.add_argument("--with-merge", action="store_true", help="compute: also submit merge after the arrays")
    p.add_argument("--after", nargs="*", help="merge: job ids to wait for")
    a = ap.parse_args()
    cfg = Cfg(a.config)
    {"dry-run": cmd_dry, "prepare": cmd_prepare, "plan": cmd_plan, "run": cmd_run, "status": cmd_status,
     "merge": cmd_merge, "slurm": cmd_slurm}[a.cmd](cfg, a)


if __name__ == "__main__":
    main()
