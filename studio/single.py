"""
single.py - the Single runs page's backend: one run folder per submission, run locally (job) or on SLURM.

A single run is a folder  <base>/<YYYYmmdd-HHMMSS>_<name>/  with
    request.json     what the page sent (to show it again / copy it into a new run)
    inputs/          copies of the POSCARs and of the NEP model (md5 in request.json)
    scripts/         frozen copies of the scripts used (single_runs.py and what it calls)
    run.sh           the command (bash run.sh  reproduces the run)
    job.sbatch       SLURM version (when written / submitted)
    results/         written by scripts/single_runs.py (status.json, summary.csv, <label>/, plots/)
The calculation itself is scripts/single_runs.py; plots are scripts/plot_harmonic_compare.py.  Nothing here
knows the elements: species come from the POSCARs.

The endpoints are a mixin of the server's request handler (server.H), so server.py only lists the class.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parent
SCRIPTS = ("single_runs.py", "batch_free_energy.py", "xrd_average.py", "plot_harmonic_compare.py")
TASKS = ("relax", "harmonic", "xrd", "plot")
FORMATS = ("png", "pdf", "svg", "eps")
NEP_COPY_LIMIT = 200 * 1024 * 1024       # bigger NEP files are used in place (path + md5 recorded)


def _md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _rjson(p, default=None):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return default


def _num(v, what, cast=float, lo=None):
    try:
        x = cast(v)
    except (TypeError, ValueError):
        raise ValueError(f"{what}: {v!r} is not {'a whole number' if cast is int else 'a number'}")
    if x != x or (lo is not None and x < lo):
        raise ValueError(f"{what} must be >= {lo}")
    return x


def safe_name(s, default="run"):
    s = re.sub(r"[^\w.+-]+", "_", str(s or "").strip()).strip("._")
    return s[:60] or default


def with_names(spec, poscars):
    """labels = the file names when the page gave none (not the numbered copies in inputs/)"""
    names = spec.get("names") or []
    if any(str(n).strip() for n in names):
        return spec
    labs = [safe_name(Path(p).name, f"s{i + 1}") for i, p in enumerate(poscars)]
    if len(set(labs)) != len(labs):                       # same file name from different folders
        labs = [safe_name(f"{Path(p).parent.name}_{Path(p).name}", f"s{i + 1}") for i, p in enumerate(poscars)]
        if len(set(labs)) != len(labs):
            labs = [f"{n}_{i + 1}" for i, n in enumerate(labs)]
    return {**spec, "names": labs}


def build_args(spec, inputs_rel, nep_rel):
    """single_runs.py arguments from the page's request (relative to the run folder); raises ValueError"""
    tasks = [t for t in (spec.get("tasks") or []) if t in TASKS]
    if not tasks:
        raise ValueError("choose at least one task")
    if not inputs_rel:
        raise ValueError("add at least one POSCAR")
    args = [*inputs_rel, "--tasks", *tasks, "--out-dir", "results"]
    if {"relax", "harmonic"} & set(tasks):
        if not nep_rel:
            raise ValueError("relax / harmonic need a NEP model (nep.txt)")
        args += ["--nep", nep_rel]
    names = spec.get("names") or []
    if names and any(str(n).strip() for n in names):
        if len(names) != len(inputs_rel):
            raise ValueError("one name per POSCAR")
        labs = [safe_name(n, f"s{i + 1}") for i, n in enumerate(names)]
        if len(set(labs)) != len(labs):
            raise ValueError("the structure names must be different")
        args += ["--names", *labs]
    if "harmonic" in tasks:
        T = [_num(t, "T", lo=0) for t in (spec.get("T") or []) if str(t).strip() != ""]
        rng = spec.get("T_range")
        if rng:
            lo, hi, st = (_num(x, "T range", lo=0) for x in rng)
            if hi < lo or st <= 0:
                raise ValueError("T range: TMIN <= TMAX and STEP > 0")
            if (hi - lo) / st > 20000:
                raise ValueError("T range: more than 20000 temperatures")
            args += ["--T-range", f"{lo:g}", f"{hi:g}", f"{st:g}"]
        if T:
            args += ["--T", *[f"{t:g}" for t in T]]
        if not T and not rng:
            raise ValueError("harmonic: give a temperature or a temperature range")
        sc = spec.get("supercell")
        if sc:
            args += ["--supercell", *[str(_num(x, "supercell", int, 1)) for x in sc]]
        else:
            args += ["--min-length", f"{_num(spec.get('min_length', 10), 'min length', lo=0.1):g}"]
        args += ["--mesh", f"{_num(spec.get('mesh', 100), 'mesh', lo=1):g}"]
    if "relax" in tasks:
        args += ["--fmax", f"{_num(spec.get('fmax', 1e-4), 'fmax', lo=1e-12):g}",
                 "--steps", str(_num(spec.get("steps", 2000), "steps", int, 1))]
        if spec.get("fixed_cell"):
            args.append("--fixed-cell")
    if spec.get("threads"):
        args += ["--threads", str(_num(spec["threads"], "threads", int, 1))]
    if "xrd" in tasks:
        refs = [str(r).strip() for r in (spec.get("xrd_ref") or []) if str(r).strip()]
        if refs:
            args += ["--xrd-ref", *refs]
        if spec.get("refs_dir"):
            args += ["--refs", str(Path(spec["refs_dir"]).expanduser())]
        if str(spec.get("xrd_args") or "").strip():
            shlex.split(spec["xrd_args"])            # raises ValueError on bad quoting
            args += ["--xrd-args", spec["xrd_args"]]
    if "plot" in tasks:
        ref = str(spec.get("reference") or "lowest")
        args += ["--reference", safe_name(ref, "lowest") if ref not in ("lowest", "none") else ref]
        at = [_num(t, "dF at T", lo=0) for t in (spec.get("at_T") or []) if str(t).strip() != ""]
        if at:
            args += ["--at-T", *[f"{t:g}" for t in at]]
        fm = [f for f in (spec.get("formats") or ["png", "pdf", "svg"]) if f in FORMATS] or ["png"]
        args += ["--formats", *fm]
    return args


def sbatch_text(spec, cmd_line, run_dir):
    sl = spec.get("slurm") or {}
    name = safe_name(spec.get("name"), "single")
    lines = ["#!/bin/bash", f"#SBATCH --job-name=sr_{name}", "#SBATCH --nodes=1", "#SBATCH --ntasks=1",
             f"#SBATCH --cpus-per-task={_num(sl.get('cpus') or 1, 'cpus', int, 1)}",
             f"#SBATCH --time={sl.get('time') or '04:00:00'}", f"#SBATCH --output={run_dir}/slurm-%j.out"]
    for k, flag in (("partition", "--partition"), ("account", "--account"), ("mem", "--mem"), ("qos", "--qos")):
        if str(sl.get(k) or "").strip():
            v = str(sl[k]).strip()
            if not re.fullmatch(r"[\w.,:=+-]+", v):
                raise ValueError(f"SLURM {k}: {v!r} has characters that are not allowed")
            lines.append(f"#SBATCH {flag}={v}")
    if not re.fullmatch(r"[\d:-]+", str(sl.get("time") or "04:00:00")):
        raise ValueError("SLURM time: use D-HH:MM:SS or HH:MM:SS")
    for x in str(sl.get("extra") or "").splitlines():
        x = x.strip()
        if x:
            if not x.startswith("#SBATCH "):
                raise ValueError(f"extra SLURM lines must start with '#SBATCH ': {x}")
            lines.append(x)
    lines += ["", "set -e", f"cd {shlex.quote(str(run_dir))}"]
    if str(sl.get("setup") or "").strip():
        lines += ["# setup (modules, conda ...)", *str(sl["setup"]).splitlines()]
    cpus = _num(sl.get("cpus") or 1, "cpus", int, 1)
    lines += [f"export OMP_NUM_THREADS={cpus} OPENBLAS_NUM_THREADS={cpus} MKL_NUM_THREADS={cpus}", cmd_line, ""]
    return "\n".join(lines)


def create_run(spec, base):
    """make the run folder, copy inputs and scripts, write run.sh (and job.sbatch); returns (dir, cmd_line)"""
    base = Path(base).expanduser()
    base.mkdir(parents=True, exist_ok=True)
    poscars = [Path(p).expanduser() for p in (spec.get("poscars") or [])]
    for p in poscars:
        if not p.is_file():
            raise ValueError(f"{p} not found")
    tasks = set(spec.get("tasks") or [])
    nep = Path(spec["nep"]).expanduser() if spec.get("nep") else None
    if {"relax", "harmonic"} & tasks and (not nep or not nep.is_file()):
        raise ValueError(f"NEP model {nep or '(none)'} not found")
    if not {"relax", "harmonic"} & tasks:
        nep = None
    rd = base / (time.strftime("%Y%m%d-%H%M%S") + "_" + safe_name(spec.get("name"), "single"))
    k = 2
    while rd.exists():
        rd = rd.with_name(rd.name.rsplit("~", 1)[0] + f"~{k}")
        k += 1
    # check the arguments before anything is written
    inp_rel = [f"inputs/{i + 1:02d}_{safe_name(p.name, 'POSCAR')}" for i, p in enumerate(poscars)]
    spec = with_names(spec, poscars)
    nep_rel = None
    if nep:
        nep_rel = "inputs/" + safe_name(nep.name, "nep.txt") if nep.stat().st_size <= NEP_COPY_LIMIT else str(nep.resolve())
    args = build_args(spec, inp_rel, nep_rel)
    mode = spec.get("mode", "local")
    py = str(spec.get("python") or "").strip() or sys.executable
    if mode in ("slurm_write", "slurm_submit") and spec.get("slurm", {}).get("python"):
        py = str(spec["slurm"]["python"]).strip()
    cmd_line = shlex.join([py, "scripts/single_runs.py", *args])
    sb = sbatch_text(spec, cmd_line, rd.resolve()) if mode in ("slurm_write", "slurm_submit") else None

    (rd / "inputs").mkdir(parents=True)
    (rd / "scripts").mkdir()
    for p, r in zip(poscars, inp_rel):
        shutil.copyfile(p, rd / r)
    if nep and not nep_rel.startswith("/"):
        shutil.copyfile(nep, rd / nep_rel)
    for s in SCRIPTS:
        shutil.copyfile(APP / "scripts" / s, rd / "scripts" / s)
    (rd / "run.sh").write_text("#!/bin/bash\n# single run written by Fluorite Ordering Studio; run it again with: bash run.sh\n"
                               f"set -e\ncd \"$(dirname \"$0\")\"\n{cmd_line}\n")
    if sb:
        (rd / "job.sbatch").write_text(sb)
    from . import runs as R
    rec = {"app_version": R.APP_VERSION, "created": time.strftime("%Y-%m-%d %H:%M:%S"), "mode": mode,
           "spec": spec, "command": cmd_line,
           "inputs": [{"source": str(p.resolve()), "copy": r, "md5": _md5(p)} for p, r in zip(poscars, inp_rel)],
           "nep": {"source": str(nep.resolve()), "used": nep_rel, "md5": _md5(nep)} if nep else None,
           "scripts": {s: _md5(rd / "scripts" / s) for s in SCRIPTS}}
    (rd / "request.json").write_text(json.dumps(rec, indent=1))
    return rd, cmd_line


def list_runs(base):
    base = Path(base).expanduser()
    out = []
    if not base.is_dir():
        return out
    for d in sorted(base.iterdir(), reverse=True)[:300]:
        rq = _rjson(d / "request.json") if d.is_dir() else None
        if not isinstance(rq, dict):
            continue
        st = _rjson(d / "results" / "status.json", {}) or {}
        sp = rq.get("spec") or {}
        out.append({"folder": str(d), "name": d.name, "created": rq.get("created"), "mode": rq.get("mode"),
                    "tasks": sp.get("tasks"), "n": len(rq.get("inputs") or []),
                    "state": st.get("state") or ("written" if rq.get("mode") == "slurm_write" else "waiting"),
                    "done": st.get("done"), "failed": st.get("failed")})
    return out


def _csv(p):
    try:
        with open(p, newline="") as f:
            return list(csv.DictReader(f))
    except OSError:
        return []


def result(folder):
    d = Path(folder).expanduser()
    rq = _rjson(d / "request.json")
    if not isinstance(rq, dict):
        raise ValueError(f"{d} is not a single run (no request.json)")
    res = d / "results"
    st = _rjson(res / "status.json", {}) or {}
    structs = []
    for lab in (st.get("structures") or {}):
        r = _rjson(res / lab / "result.json", {}) or {}
        files = sorted(str(p) for p in (res / lab).rglob("*") if p.is_file()) if (res / lab).is_dir() else []
        structs.append({"label": lab, "state": st["structures"][lab], "flags": r.get("flags", []),
                        "notes": r.get("notes", []), "files": files[:200]})
    plots = sorted(str(p) for p in (res / "plots").iterdir() if p.is_file()) if (res / "plots").is_dir() else []
    slurm_out = sorted(str(p) for p in d.glob("slurm-*.out"))
    return {"ok": True, "folder": str(d), "request": rq, "status": st, "summary": _csv(res / "summary.csv"),
            "structures": structs, "plots": plots, "summary_csv": str(res / "summary.csv"),
            "comparison_csv": str(res / "plots" / "comparison.csv"), "slurm_out": slurm_out,
            "has_sbatch": (d / "job.sbatch").is_file(), "results": str(res)}


def replot_args(spec, res_dir):
    args = [str(res_dir)]
    ref = str(spec.get("reference") or "lowest")
    args += ["--reference", ref if ref in ("lowest", "none") else safe_name(ref)]
    inc = [safe_name(x) for x in (spec.get("include") or []) if str(x).strip()]
    if inc:
        args += ["--include", *inc]
    w = spec.get("T_window")
    if w and all(str(x).strip() != "" for x in w):
        args += ["--T-window", *[f"{_num(x, 'T window', lo=0):g}" for x in w]]
    at = [_num(t, "dF at T", lo=0) for t in (spec.get("at_T") or []) if str(t).strip() != ""]
    if at:
        args += ["--at-T", *[f"{t:g}" for t in at]]
    fm = [f for f in (spec.get("formats") or ["png", "pdf", "svg"]) if f in FORMATS] or ["png"]
    args += ["--formats", *fm]
    for k, flag in (("width", "--width"), ("height", "--height"), ("font_size", "--font-size")):
        if spec.get(k) not in (None, ""):
            args += [flag, f"{_num(spec[k], k, lo=0.5):g}"]
    return args


# ============================================================================ HTTP endpoints (mixin)
class SingleRunsAPI:
    """mixed into server.H: api_<method>_single_<name>  ->  /api/single/<name>"""

    def _sr_state(self):
        from . import server
        return server.ST

    def _sr_base(self, q=None, b=None):
        ST = self._sr_state()
        v = (b or {}).get("base") or ((q or {}).get("base") or [None])[0]
        return Path(v).expanduser() if v else (ST.prefs().get("single_base") and Path(ST.prefs()["single_base"])) \
            or ST.start_dir / "single_runs"

    def api_get_single_defaults(self, q):
        from . import jobs as J
        refs = APP / "refs"
        self._json({"ok": True, "base": str(self._sr_base(q)), "refs_dir": str(refs),
                    "refs": sorted(p.name[len("POSCAR_"):] for p in refs.glob("POSCAR_*")),
                    "slurm": J.slurm_available(), "python": sys.executable,
                    "scripts_ok": {s: (APP / "scripts" / s).is_file() for s in SCRIPTS},
                    "needs": {m: _has(m) for m in ("numpy", "ase", "calorine", "phonopy", "pymatgen", "matplotlib")}})

    def api_get_single_refs(self, q):
        d = Path((q.get("dir") or [str(APP / "refs")])[0]).expanduser()
        if not d.is_dir():
            return self._err(f"{d} is not a folder")
        self._json({"ok": True, "dir": str(d), "refs": sorted(p.name[len("POSCAR_"):] for p in d.glob("POSCAR_*"))})

    def api_post_single_start(self, q):
        ST = self._sr_state()
        spec = json.loads(self._body() or b"{}")
        base = self._sr_base(b=spec)
        mode = spec.get("mode", "local")
        if mode not in ("local", "slurm_write", "slurm_submit", "preview"):
            return self._err("mode must be local, slurm_write, slurm_submit or preview")
        if mode == "preview":
            try:
                ps = spec.get("poscars") or []
                inp = [f"inputs/{i + 1:02d}_{safe_name(Path(p).name)}" for i, p in enumerate(ps)]
                args = build_args(with_names(spec, ps), inp, "inputs/" + safe_name(Path(spec["nep"]).name) if spec.get("nep") else None)
            except ValueError as e:
                return self._err(str(e))
            return self._json({"ok": True, "cmd": shlex.join(["python3", "scripts/single_runs.py", *args])})
        try:
            rd, cmd_line = create_run(spec, base)
        except (ValueError, OSError) as e:
            return self._err(str(e))
        ST.save_prefs({"single_base": str(base)})
        out = {"ok": True, "folder": str(rd), "cmd": cmd_line}
        title = f"Single run - {rd.name}"
        watch = [str(rd), "R:" + str(rd / "results")]
        if mode == "local":
            out["job"] = ST.jobs.start(["bash", str(rd / "run.sh")], rd, title=title, module="single_runs",
                                       watch=watch, meta={"kind": "single", "folder": str(rd)})
        elif mode == "slurm_submit":
            out["job"] = ST.jobs.start(["sbatch", str(rd / "job.sbatch")], rd, title=title + " [SLURM]",
                                       module="single_runs", watch=watch,
                                       meta={"kind": "single", "mode": "slurm", "folder": str(rd)})
        else:
            out["sbatch"] = str(rd / "job.sbatch")
        self._json(out)

    def api_get_single_list(self, q):
        base = self._sr_base(q)
        self._json({"ok": True, "base": str(base), "runs": list_runs(base)})

    def api_get_single_result(self, q):
        try:
            self._json(result((q.get("folder") or [""])[0]))
        except ValueError as e:
            self._err(str(e))

    def api_post_single_replot(self, q):
        ST = self._sr_state()
        b = json.loads(self._body() or b"{}")
        rd = Path(str(b.get("folder") or "")).expanduser()
        if not (rd / "request.json").is_file() or not (rd / "results").is_dir():
            return self._err("not a single run with results")
        script = rd / "scripts" / "plot_harmonic_compare.py"
        from .server import file_untrusted
        if not script.is_file() or file_untrusted(script) or file_untrusted(script.parent):
            script = APP / "scripts" / "plot_harmonic_compare.py"     # never run someone else's copy
        try:
            args = replot_args(b, rd / "results")
        except ValueError as e:
            return self._err(str(e))
        cmd = [sys.executable, str(script), *args]
        jid = ST.jobs.start(cmd, rd, title=f"Single run plots - {rd.name}", module="single_runs_plot",
                            watch=["R:" + str(rd / "results" / "plots")], meta={"kind": "single", "folder": str(rd)})
        self._json({"ok": True, "job": jid, "cmd": shlex.join(cmd)})


def _has(m):
    import importlib.util
    return importlib.util.find_spec(m) is not None
