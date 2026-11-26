"""
state.py - the state of one run for the run dashboard: counts, the 7 pipeline stages, SLURM tasks, the
"Attention" list (recovery actions in the order to do them) and the next step.  Standard library only.

Everything is derived from the files campaign.py / rerun_flagged.py write (the scan), the run's plan.json and
result CSVs, the app's job records (which hold the sbatch output) and squeue / sacct.  Nothing is written.

    state = run_state(config_path, scan_result, jobs, slurm)      # slurm = {"squeue": [...], "sacct": [...]}
"""
from __future__ import annotations

import csv
import json
import math
import re
import statistics
import time
from pathlib import Path

from . import runs as R

STAGES = [  # key, module id, title, sub-label
    ("check", "dry_run", "Check", "dry-run"),
    ("prepare", "prepare", "Prepare", "enumerate · map"),
    ("plan", "plan", "Plan", "chunks"),
    ("compute", "compute", "Compute", "relax + F"),
    ("status", "status", "Status", "timing"),
    ("merge", "merge", "Merge", "CSV"),
    ("rerun", "rerun_flagged", "Rerun flagged", "recovery"),
]
ACTIVE = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED", "RESIZING", "SUSPENDED", "SUBMITTED"}
SPS_TOLERANCE = 0.25           # measured sec_per_structure differs from the plan by more -> re-plan suggested


# ---------------------------------------------------------------------- cost model (same as campaign.py)
def _read_cell(path):
    L = Path(path).read_text(errors="replace").split("\n", 8)
    sc = float(L[1].split()[0])
    lat = [[float(v) for v in L[k].split()[:3]] for k in (2, 3, 4)]
    det = _det(lat)
    f = (abs(sc) / abs(det)) ** (1 / 3) if sc < 0 else sc
    lat = [[x * f for x in r] for r in lat]
    names = L[5].split()
    cnt = L[6].split() if names and not names[0].isdigit() else names
    return lat, sum(int(x) for x in cnt)


def _det(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def _inv(m):
    d = _det(m)
    c = [[(m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3] -
           m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3]) / d for j in range(3)] for i in range(3)]
    return c


def phonon_sc(fe, lat):
    if fe.get("supercell"):
        return [int(v) for v in fe["supercell"]]
    inv = _inv(lat)
    out = []
    for j in range(3):                        # np.linalg.norm(inv(lat), axis=0): column norms of the inverse
        w = 1.0 / math.sqrt(sum(inv[i][j] ** 2 for i in range(3)))
        out.append(max(1, int(math.ceil(float(fe.get("min_length", 10.0)) / w - 1e-6))))
    return out


_UNITS = {}


def cost_units(path, fe, n0):
    p = Path(path)
    try:
        key = (str(p), p.stat().st_mtime, json.dumps(fe.get("supercell")), fe.get("min_length"), n0)
    except OSError:
        return None
    if key not in _UNITS:
        try:
            lat, n = _read_cell(p)
            sc = phonon_sc(fe, lat)
            _UNITS[key] = n * n * sc[0] * sc[1] * sc[2] / float(n0) ** 2
        except (OSError, ValueError, IndexError, ZeroDivisionError):
            _UNITS[key] = None
        if len(_UNITS) > 400000:
            _UNITS.clear()
    return _UNITS[key]


# ---------------------------------------------------------------------- results: timing and core-hours
_TIMING = {}


def timing(work: Path, fe, n0, threads, T, max_rows=3000):
    """measured sec_per_structure (median, 10-90 %) and core-hours used, from the ok rows of results/chunk_*.csv"""
    res = sorted((work / "results").glob("chunk_*.csv")) if (work / "results").is_dir() else []
    sig = tuple((p.name, p.stat().st_mtime, p.stat().st_size) for p in res) + (json.dumps(fe, sort_keys=True, default=str), n0, T)
    hit = _TIMING.get(str(work))
    if hit and hit[0] == sig:
        return hit[1]
    best, camp = {}, {}
    for p in res:
        rerun = p.name.startswith("chunk_rerun")
        try:
            with p.open(newline="") as f:
                for r in csv.DictReader(f):
                    fn = r.get("file")
                    if not fn:
                        continue
                    t = _f(r.get("time_s"))
                    if T is None or r.get("T_K") in (None, "") or abs(_f(r["T_K"]) - T) < 1e-6:
                        if fn not in best or r.get("status") == "ok" or best[fn][0] != "ok":
                            best[fn] = (r.get("status"), t)
                    if r.get("status") == "ok" and not rerun and math.isfinite(t) and t > 0:
                        camp[fn] = t               # newest campaign row per structure (as campaign.py status)
        except (OSError, csv.Error, UnicodeDecodeError):
            continue
    secs = [t for st, t in best.values() if math.isfinite(t)]
    used_core_h = sum(secs) * max(1, int(threads)) / 3600.0
    # measured cost: campaign rows only (rerun rows use stricter settings), one per structure
    ok = list(camp.items())
    step = max(1, len(ok) // max_rows)
    meas = []
    for fn, t in ok[::step]:
        u = cost_units(work / fn, fe, n0)
        if u:
            meas.append(t / u)
    out = {"used_core_h": used_core_h, "measured": None}
    if meas:
        meas.sort()
        q = lambda x: meas[min(len(meas) - 1, max(0, int(round(x * (len(meas) - 1)))))]
        out["measured"] = {"n": len(ok), "sampled": len(meas), "median": statistics.median(meas),
                           "p10": q(0.10), "p90": q(0.90)}
    _TIMING[str(work)] = (sig, out)
    return out


def _f(v, d=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def parse_walltime(w):
    try:
        parts = [int(x) for x in re.split(r"[-:]", str(w))]
    except ValueError:
        return None
    if "-" in str(w):
        d, rest = parts[0], parts[1:] + [0] * (3 - len(parts[1:]))
        return d * 86400 + rest[0] * 3600 + rest[1] * 60 + rest[2]
    while len(parts) < 3:
        parts = [0] + parts
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


def hms(s):
    s = int(round(max(0, s)))
    return f"{s // 3600:d}:{s % 3600 // 60:02d}:{s % 60:02d}"


# ---------------------------------------------------------------------- SLURM: what the app submitted
def array_ids(spec):
    """'0-4,7%50' -> [0, 1, 2, 3, 4, 7]"""
    out = []
    for part in str(spec).split("%")[0].strip("[]").split(","):
        m = re.fullmatch(r"(\d+)(?:-(\d+))?(?::(\d+))?", part.strip())
        if m:
            a, b, st = int(m[1]), int(m[2] or m[1]), int(m[3] or 1)
            out += list(range(a, b + 1, st))
    return out


def submissions(log_text, job_id=None, t=None):
    """sbatch calls printed by campaign.py / rerun_flagged.py ('$ sbatch ... X.sh' then 'submitted ... ID')"""
    out, pending = [], None
    for line in log_text.splitlines():
        s = line.strip()
        if s.startswith("$ sbatch") or s.startswith("sbatch --parsable"):
            arr = re.search(r"--array=(\S+)", s)
            off = re.search(r"CHUNK_OFFSET=(\d+)", s)
            script = s.split()[-1]
            name = Path(script).name
            kind = ("merge" if "merge" in name else "rerun" if "rerun" in name else
                    "compute" if "compute" in name else "prepare" if "prepare" in name else "other")
            pending = {"script": name, "kind": kind, "array": arr[1] if arr else None,
                       "offset": int(off[1]) if off else 0}
            continue
        m = re.search(r"submitted (?:merge job|job|\w+ array) (\d+)", s) or re.search(r"Submitted batch job (\d+)", s)
        if m:
            sub = pending or {"script": "", "kind": "merge" if "merge" in s else "other", "array": None, "offset": 0}
            if "merge job" in s:
                sub = {"script": "merge.sh", "kind": "merge", "array": None, "offset": 0}
            sub = dict(sub, id=m[1], app_job=job_id, time=t)
            sub["tasks"] = array_ids(sub["array"]) if sub["array"] else None
            out.append(sub)
            pending = None
    return out


_SUBS = {}


def slurm_ids(jobs, config):
    """SLURM job ids submitted by the app's jobs of this run"""
    out = []
    for j in jobs.list(400):
        if (j.get("meta") or {}).get("config") == str(config):
            out += j.get("slurm_ids") or []
    return sorted(set(out), key=int)


def _slurm_rows(sq, sa):
    """{base id: {task index or None: {state, time, limit, where}}} from squeue (live) and sacct (history)"""
    out = {}

    def put(jid, st):
        m = re.fullmatch(r"(\d+)(?:_(\[[^\]]*\]|\d+))?", str(jid))
        if not m:
            return
        base, task = m[1], m[2]
        d = out.setdefault(base, {})
        if task is None:
            d.setdefault(None, {}).update(st)
        elif task.startswith("["):
            for i in array_ids(task):
                d.setdefault(i, {}).update(st)
        else:
            d.setdefault(int(task), {}).update(st)
    for r in sa or []:
        put(r.get("id"), {"state": (r.get("state") or "").split()[0], "time": r.get("elapsed", ""),
                          "limit": r.get("limit", ""), "where": r.get("nodes", ""), "exit": r.get("exit", "")})
    for r in sq or []:                         # live view wins
        put(r.get("id"), {"state": r.get("state", ""), "time": r.get("time", ""), "limit": r.get("limit", ""),
                          "where": r.get("reason", "")})
    return out


def _tsec(t):
    """squeue/sacct time 'D-HH:MM:SS' | 'HH:MM:SS' | 'MM:SS' -> seconds"""
    return parse_walltime(t) if t and re.fullmatch(r"(\d+-)?\d+(:\d+){0,2}", str(t)) else None


# ---------------------------------------------------------------------- the state
def comp_list(spec, N):
    if spec == "all" or spec is None:
        return list(range(N + 1))
    if isinstance(spec, str):
        m = re.fullmatch(r"\s*(\d+)\s*:\s*(\d+)\s*(?::\s*(\d+))?\s*", spec)
        if not m:
            return []
        return list(range(int(m[1]), int(m[2]) + 1, int(m[3] or 1)))
    return sorted({int(v) for v in spec})


def _merged_rows(p: Path, cache={}):          # noqa: B006 - module-level cache on purpose
    try:
        key = (str(p), p.stat().st_mtime, p.stat().st_size)
    except OSError:
        return None
    if key not in cache:
        with p.open("rb") as f:
            cache[key] = max(0, sum(1 for _ in f) - 1)
    return cache[key]


def run_state(config_path, scan, jobs, slurm=None, now=None):
    now = now or time.time()
    cp = Path(config_path).expanduser().resolve()
    info = R.read_run(cp)
    cfg = info.get("config") or {}
    rec = info.get("record") or {}
    sec = lambda k: {**R.DEFAULTS.get(k, {}), **(cfg.get(k) or {})}
    st_, fe, cl, mp = sec("structure"), sec("free_energy"), sec("cluster"), sec("mapping")
    work = Path(info.get("workdir") or cp.parent / "campaign")
    T = float(fe.get("T", 300.0))
    Tk = f"{T:g}"

    host = R.peek_poscar((cp.parent / st_["poscar"]) if st_.get("poscar") else cp.parent / "_none_")
    n0 = host.get("n_atoms") if host.get("ok") else None
    n_host = (host.get("species") or {}).get(st_.get("host")) if host.get("ok") else None
    comps_cfg = comp_list(st_.get("compositions"), n_host) if n_host else []

    # ---- counts from the scan (at the run's temperature)
    tot = (scan or {}).get("totals", {}).get(Tk) or {}
    counts = {k: int(tot.get(k, 0)) for k in ("structures", "planned", "relaxed", "ok", "failed", "imaginary",
                                               "unconverged", "clean", "not_run", "unplanned")}
    counts["flagged"] = counts["ok"] - counts["clean"]
    comps = (scan or {}).get("compositions", [])
    prepared = [c for c in comps if c.get("prepared") and (c.get("n") in comps_cfg or not comps_cfg)]
    n_comp = len(comps_cfg) or len(comps)

    # ---- plan
    plan = None
    pj = work / "chunks" / "plan.json"
    if pj.is_file():
        try:
            plan = json.loads(pj.read_text())
        except (OSError, ValueError):
            plan = None
    chunks = (scan or {}).get("chunks") or None
    chunk_by_id = {c["id"]: c for c in (chunks or {}).get("list", [])}
    unfinished = [c["id"] for c in (chunks or {}).get("list", []) if c["ok"] < c["n"]]
    threads = int(fe.get("threads", 1) or 1)
    planned_core_h = None
    sps_plan = None
    if plan:
        sps_plan = _f(plan.get("sec_per_structure"), None)
        u = sum(_f(c.get("units"), 0) for c in plan.get("chunks", []))
        if sps_plan:
            planned_core_h = u * sps_plan * threads / 3600.0
    sps_cfg = _f(cl.get("sec_per_structure"), None)
    tm = timing(work, fe, n0, threads, T) if n0 and work.is_dir() else {"used_core_h": 0.0, "measured": None}
    meas = tm["measured"]

    # ---- the app's jobs for this run, and what they submitted to SLURM
    mine = [j for j in jobs.list(400) if (j.get("meta") or {}).get("config") == str(cp)]
    last = {}
    for j in mine:                                  # newest first
        last.setdefault(j.get("module"), j)
    running_mods = {j.get("module") for j in mine if j["state"] in ("running", "starting")}
    subs = []
    for j in mine[:60]:
        if not j.get("slurm_ids"):
            continue
        hit = _SUBS.get(j["id"])
        if hit is None:
            g = jobs.get(j["id"], log_tail=2_000_000, with_outputs=False) or {}
            hit = submissions(g.get("log", ""), j["id"], j.get("start"))
            if j["state"] not in ("running", "starting"):
                _SUBS[j["id"]] = hit             # an ended job's log does not change any more
        subs += hit
    sl = _slurm_rows((slurm or {}).get("squeue"), (slurm or {}).get("sacct"))
    tasks = []
    for s in subs:
        rows = sl.get(s["id"], {})
        idx = s["tasks"] if s["tasks"] is not None else [None]
        recent = s.get("time") and now - s["time"] < 120
        for i in idx:
            r = rows.get(i) or {}
            state = r.get("state") or ("SUBMITTED" if recent else "ENDED" if (slurm or {}).get("available")
                                        else "UNKNOWN")   # gone from squeue and no sacct: finished
            chunk = (i + s["offset"]) if (s["kind"] == "compute" and i is not None) else None
            c = chunk_by_id.get(chunk) if chunk is not None else None
            tasks.append({"job": s["id"], "task": i, "kind": s["kind"], "state": state, "time": r.get("time", ""),
                          "limit": r.get("limit", ""), "where": r.get("where", ""), "chunk": chunk,
                          "n": c["n"] if c else None, "done": (c["ok"] + c["failed"]) if c else None,
                          "elapsed_s": _tsec(r.get("time")), "limit_s": _tsec(r.get("limit")),
                          "submitted": s.get("time"), "app_job": s.get("app_job")})
    for t in tasks:                                  # projection of a running chunk: elapsed * n / done
        if t["state"] == "RUNNING" and t["elapsed_s"] and t["n"]:
            t["projected_s"] = t["elapsed_s"] * t["n"] / t["done"] if t["done"] else None
    active = lambda kind: [t for t in tasks if t["kind"] == kind and t["state"] in ACTIVE]
    compute_running = "compute" in running_mods or bool(active("compute"))
    rerun_running = "rerun_flagged" in running_mods or bool(active("rerun"))
    prepare_running = "prepare" in running_mods or bool(active("prepare"))
    merge_running = "merge" in running_mods or bool(active("merge"))
    timeouts = [t for t in tasks if t["state"] == "TIMEOUT"]
    slurm_sum = {}
    for t in tasks:
        slurm_sum[t["state"]] = slurm_sum.get(t["state"], 0) + 1

    # ---- merged CSV
    merged = work / f"free_energy_{Tk}K_all.csv"
    res_m = max((p.stat().st_mtime for p in (work / "results").glob("chunk_*.csv")), default=0) \
        if (work / "results").is_dir() else 0
    merge = {"exists": merged.is_file(), "path": str(merged), "rows": _merged_rows(merged) if merged.is_file() else None,
             "mtime": merged.stat().st_mtime if merged.is_file() else None}
    merge["stale"] = bool(res_m and (not merge["exists"] or res_m > merge["mtime"] + 1))
    reruns = (scan or {}).get("reruns", [])

    def jstate(mod):
        j = last.get(mod)
        return None if not j else {"id": j["id"], "state": j["state"], "start": j["start"], "end": j.get("end"),
                                   "mode": (j.get("meta") or {}).get("mode"), "exit": j.get("exit")}

    # ---- the 7 stages
    stages = []

    def add(key, state, pill, lines):
        k = [s for s in STAGES if s[0] == key][0]
        stages.append({"key": key, "module": k[1], "title": k[2], "sub": k[3], "state": state,
                       "pill": pill, "lines": lines, "job": jstate(k[1])})
    j = jstate("dry_run")
    if j and j["state"] in ("running", "starting"):
        add("check", "running", ["run", "Running"], ["dry-run in progress"])
    elif j and j["state"] == "done":
        add("check", "done", ["ok", "Passed"], [time.strftime("%b %d, %H:%M", time.localtime(j["end"] or j["start"]))])
    elif j and j["state"] in ("failed", "lost"):
        add("check", "bad", ["bad", "Failed"], ["see the log"])
    else:
        add("check", "done" if prepared else "idle", ["ok", "Skipped"] if prepared else ["idle", "Not run"],
            ["checks inputs, counts, estimate"])
    np_, nc = len(prepared), n_comp
    if prepare_running:
        add("prepare", "running", ["run", f"{np_}/{nc}"], [f"{np_} of {nc} compositions prepared"])
    elif nc and np_ >= nc:
        add("prepare", "done", ["ok", "Done"], [f"{counts['structures']:,} structures", f"{nc} compositions"])
    elif np_:
        add("prepare", "warn", ["warn", f"{np_}/{nc}"], [f"{nc - np_} compositions not prepared"])
    else:
        add("prepare", "idle", ["idle", "Not run"], ["enumerate, map, verify"])
    sps_off = bool(meas and sps_plan and meas["n"] >= 3 and abs(meas["median"] / sps_plan - 1) > SPS_TOLERANCE)
    work_left = bool(unfinished or counts["unplanned"])
    if "plan" in running_mods:
        add("plan", "running", ["run", "Running"], ["packing chunks"])
    elif not plan:
        add("plan", "idle", ["idle", "Not planned"], [f"sec_per_structure {sps_cfg:g} s" if sps_cfg else ""])
    else:
        lines = [f"{len(plan.get('chunks', []))} chunks · {plan.get('walltime')}", f"{sps_plan:g} s/structure planned" if sps_plan else ""]
        if sps_off and work_left:
            add("plan", "warn", ["warn", "Stale"], lines)
        elif counts["unplanned"]:
            add("plan", "warn", ["warn", f"{counts['unplanned']:,} new"], lines)
        elif sps_cfg and sps_plan and abs(sps_cfg - sps_plan) > 1e-9 and work_left:
            add("plan", "warn", ["warn", "Re-pack"], lines + [f"config now {sps_cfg:g} s"])
        else:
            add("plan", "done", ["ok", "Planned"], lines)
    nch = (chunks or {}).get("n", 0)
    ncomp = (chunks or {}).get("complete", 0)
    seg = {"done": ncomp, "unfinished": len(unfinished), "n": nch,
           "running": len([t for t in tasks if t["kind"] == "compute" and t["state"] == "RUNNING"]),
           "pending": len([t for t in tasks if t["kind"] == "compute" and t["state"] in ("PENDING", "SUBMITTED")]),
           "timeout": len([t for t in timeouts if t["kind"] == "compute"])}
    if compute_running:
        add("compute", "running", ["run", f"{ncomp}/{nch}"], [f"{seg['running']} running · {seg['pending']} pending"
                                                               if seg["running"] + seg["pending"] else "running here"])
    elif nch and not unfinished:
        add("compute", "done", ["ok", f"{ncomp}/{nch}"], ["all chunks complete"])
    elif nch and (chunks or {}).get("started"):
        add("compute", "warn", ["warn", f"{ncomp}/{nch}"], [f"{len(unfinished)} unfinished" +
                                                             (f" · {seg['timeout']} timeout" if seg["timeout"] else "")])
    else:
        add("compute", "idle", ["idle", "Not started" if nch else "No plan"], [f"{nch} chunks" if nch else ""])
    stages[-1]["segments"] = seg
    if meas:
        rel = (meas["median"] / sps_plan - 1) if sps_plan else None
        add("status", "warn" if sps_off else "done", ["warn" if sps_off else "ok",
                                                    f"{rel:+.0%}" if rel is not None else "Measured"],
            [f"measured {meas['median']:.0f} s/structure", f"over {meas['n']:,} structures"])
    else:
        add("status", "idle", ["idle", "No timing"], ["after the first finished structures"])
    if merge_running:
        add("merge", "running", ["run", "Running"], ["merging"])
    elif merge["exists"]:
        add("merge", "warn" if merge["stale"] else "done", ["warn", "Stale"] if merge["stale"] else ["ok", "Merged"],
            [f"{merge['rows'] or 0:,} rows", "newer results exist" if merge["stale"] else "up to date"])
    else:
        add("merge", "idle", ["idle", "Not merged"], ["one CSV for the plots"])
    fl = counts["flagged"]
    if rerun_running:
        add("rerun", "running", ["run", "Running"], [f"round {len(reruns)}"])
    elif fl:
        add("rerun", "warn", ["warn", f"{fl:,} flagged"], [f"{counts['imaginary']:,} imaginary · {counts['unconverged']:,} unconv.",
                                                          f"{len(reruns)} round(s) so far"])
    elif reruns:
        add("rerun", "done", ["ok", f"{len(reruns)} round(s)"], ["nothing flagged"])
    else:
        add("rerun", "idle", ["idle", "Nothing flagged" if counts["ok"] else "Later"], ["imaginary / unconverged rows"])

    # ---- attention: recovery actions in the recommended order
    att = []
    wall = parse_walltime(cl.get("walltime"))
    if sps_off and work_left:
        biggest = max((_f(c.get("est_s"), 0) for c in (plan or {}).get("chunks", [])), default=0)
        need = biggest * meas["median"] / sps_plan if sps_plan else 0
        att.append({"id": "replan", "sev": "warn", "button": "Re-plan", "primary": True,
                    "title": f"sec_per_structure is off by {meas['median'] / sps_plan - 1:+.0%}",
                    "desc": f"Planned {sps_plan:g} s, measured median {meas['median']:.0f} s over {meas['n']:,} structures "
                            f"(10–90 %: {meas['p10']:.0f}–{meas['p90']:.0f} s). The largest chunk now needs ≈ {hms(need)}"
                            f" of the {cl.get('walltime')} walltime." + (" Chunks will time out." if wall and need > wall else ""),
                    "action": {"panel": "plan", "preset": {"apply_measured": round(meas["median"], 1), "choice": "--force"}},
                    "blocked": "a compute job of this run is still running: wait for it (or cancel it) before re-packing"
                    if compute_running else None})
    elif plan and sps_cfg and sps_plan and abs(sps_cfg - sps_plan) > 1e-9 and work_left:
        att.append({"id": "repack", "sev": "warn", "button": "Re-pack", "primary": True,
                    "title": "Chunks were packed with an older sec_per_structure",
                    "desc": f"The plan used {sps_plan:g} s, the config now has {sps_cfg:g} s. Re-pack the remaining structures.",
                    "action": {"panel": "plan", "preset": {"choice": "--force"}},
                    "blocked": "a compute job of this run is still running" if compute_running else None})
    if nc and np_ < nc and not prepare_running:
        att.append({"id": "prepare", "sev": "info", "button": "Prepare", "title": f"{nc - np_} composition(s) not prepared",
                    "desc": "Enumerate, map and verify the missing compositions.", "action": {"panel": "prepare"}})
    if counts["unplanned"] and plan and "plan" not in running_mods:
        att.append({"id": "append", "sev": "info", "button": "Add to plan",
                    "title": f"{counts['unplanned']:,} prepared structures are in no chunk",
                    "desc": "Add them as new chunks (plan --append); finished chunks are not touched.",
                    "action": {"panel": "plan", "preset": {"choice": "--append"}}})
    if unfinished and not compute_running:
        to = [t["chunk"] for t in timeouts if t["kind"] == "compute"]
        att.append({"id": "resubmit", "sev": "bad" if to else "warn", "button": "Resubmit unfinished",
                    "title": f"{len(unfinished)} chunk(s) unfinished" + (f" ({len(to)} TIMEOUT)" if to else ""),
                    "desc": f"Chunks {_compress(unfinished)[:120]} have structures that are not done. Resubmitting runs "
                            "only those chunks; finished structures are skipped.",
                    "action": {"panel": "compute", "preset": {"--unfinished": True}}})
    done_enough = nch and not unfinished
    if fl and not rerun_running and not compute_running:
        att.append({"id": "rerun", "sev": "warn" if done_enough else "info", "button": "Rerun",
                    "title": f"{fl:,} flagged structures",
                    "desc": f"{counts['imaginary']:,} with imaginary modes, {counts['unconverged']:,} with an unconverged "
                            f"relaxation (fmax {fe.get('fmax')} eV/Å in {fe.get('steps')} steps).",
                    "action": {"panel": "rerun"}})
    if counts["failed"]:
        att.append({"id": "failed", "sev": "info", "button": "Inspect", "title": f"{counts['failed']:,} failed structures",
                    "desc": "Look at the failures (log, NEP species) before computing them again.",
                    "action": {"results": {"status": "failed"}}})
    if merge["stale"] and not merge_running and not compute_running:
        att.append({"id": "merge", "sev": "info", "button": "Merge", "title": "Merged CSV is out of date" if merge["exists"] else "No merged CSV yet",
                    "desc": "Merge the results so Results, Plots and XRD use the newest rows.", "action": {"panel": "merge"}})

    # ---- next step
    first = next((a for a in att if not a.get("blocked")), None)
    if first:
        nxt = {"what": first["title"], "why": first["desc"], "button": first["button"], "action": first["action"],
               "att": first["id"]}
    elif prepare_running or compute_running or rerun_running or merge_running:
        nxt = {"what": "Wait for the running jobs", "why": "Nothing needs you right now. "
               + (f"{seg['running']} task(s) running, {seg['pending']} pending." if tasks else "A job of this run is running here."),
               "button": "Watch jobs", "action": {"go": "jobs"}}
    elif not stages[0]["state"] == "done":
        nxt = {"what": "Run the check (dry-run)", "why": "Reads the config and inputs and estimates the cost. Writes nothing.",
               "button": "Run check", "action": {"panel": "check"}}
    elif not prepared:
        nxt = {"what": "Prepare the structures", "why": "Enumerate orderings, map them onto every phase, verify.",
               "button": "Prepare", "action": {"panel": "prepare"}}
    elif not plan:
        nxt = {"what": "Plan the chunks", "why": "Pack the structures into chunks of one node × one walltime.",
               "button": "Plan", "action": {"panel": "plan", "preset": {"choice": ""}}}
    elif nch and not (chunks or {}).get("started"):
        nxt = {"what": "Start the computation", "why": f"{nch} chunk(s) are planned and none has started.",
               "button": "Compute", "action": {"panel": "compute"}}
    else:
        nxt = {"what": "Done: look at the results", "why": "Everything is computed and merged.", "button": "Open Results",
               "action": {"go": "results"}}

    started = rec.get("created")
    return {
        "ok": True, "time": now, "config": str(cp), "folder": str(cp.parent), "workdir": str(work),
        "run": {"name": rec.get("name") or cp.parent.name, "starting_phase": rec.get("starting_phase") or "",
                "created": started, "user": rec.get("user"), "host": rec.get("host"), "adopted": bool(rec.get("adopted")),
                "has_record": bool(rec), "notes": rec.get("notes", ""), "app_version": rec.get("app_version")},
        "material": {"host": st_.get("host"), "dopant": st_.get("dopant"), "n_atoms": n0, "n_host": n_host,
                     "species": host.get("species"), "compositions": comps_cfg, "phases": (scan or {}).get("phases", []),
                     "T": T, "nep": Path(str(fe.get("nep") or "")).name},
        "cluster": {k: cl.get(k) for k in ("partition", "account", "qos", "cores_per_node", "walltime", "max_concurrent",
                                           "sec_per_structure", "fill", "max_array", "mem", "python", "setup",
                                           "sbatch_extra", "prepare_walltime", "merge_walltime")},
        "free_energy": {k: fe.get(k) for k in ("T", "fmax", "steps", "threads", "min_length", "mesh", "relax", "fixed_cell")},
        "counts": counts, "prepared": np_, "n_comp": n_comp,
        "plan": None if not plan else {"chunks": len(plan.get("chunks", [])), "walltime": plan.get("walltime"),
                                       "sec_per_structure": sps_plan, "workers_per_node": plan.get("workers_per_node")},
        "chunks": None if not chunks else {k: chunks[k] for k in ("n", "complete", "started")} | {"unfinished": unfinished},
        "timing": {"measured": meas, "sps_config": sps_cfg, "sps_plan": sps_plan, "off": sps_off,
                   "used_core_h": tm["used_core_h"], "planned_core_h": planned_core_h},
        "merge": merge, "reruns": len(reruns), "stages": stages, "attention": att, "next": nxt,
        "jobs": [{k: j.get(k) for k in ("id", "title", "module", "state", "start", "end", "exit", "slurm_ids")}
                 | {"mode": (j.get("meta") or {}).get("mode")} for j in mine[:40]],
        "running": sorted(m for m in running_mods if m),
        "compute_running": compute_running,
        "slurm": {"available": (slurm or {}).get("available", False), "tasks": tasks, "summary": slurm_sum,
                  "submissions": [{k: s[k] for k in ("id", "kind", "array", "offset", "app_job", "time")}
                                  | {"n_tasks": len(s["tasks"]) if s["tasks"] is not None else 1} for s in subs]},
        "warnings": (scan or {}).get("warnings", []),
    }


def _compress(ids):
    ids = sorted(ids)
    parts, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        parts.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ",".join(parts)


# ---------------------------------------------------------------------- New run: live estimate (no files written)
def estimate(cfg_in):
    """what campaign.py dry-run would estimate for a config that is not saved yet: structures per phase,
    core-hours, chunks and the walltime of one chunk.  Problems are returned as a list, not raised."""
    cfg = R.coerce(cfg_in)
    st_, fe, cl, mp = cfg["structure"], cfg["free_energy"], cfg["cluster"], cfg["mapping"]
    problems = []
    try:
        R.validate(cfg)
    except ValueError as e:
        problems.append(str(e))
    host = R.peek_poscar(st_["poscar"]) if st_.get("poscar") else {"ok": False, "error": "no host POSCAR yet"}
    out = {"problems": problems, "host": None, "structures": None, "phases": [], "per_phase": {}, "compositions": [],
           "core_h": None, "node_h": None, "chunks": None, "chunk_s": None, "walltime_s": parse_walltime(cl.get("walltime")),
           "all_variants": "--all-variants" in str(mp.get("map_args", "")), "rounds": None, "wall_total_s": None}
    nep = str(fe.get("nep") or "").strip()
    if nep:
        np_ = Path(nep).expanduser()
        info = {"found": np_.is_file(), "species": None}
        if np_.is_file():
            try:
                with np_.open(errors="replace") as f:
                    w = f.readline().split()
                if w and w[0].lower().startswith("nep") and len(w) > 2 and w[1].isdigit():
                    info["species"] = w[2:2 + int(w[1])]
            except OSError:
                pass
        out["nep"] = info
    if not host.get("ok"):
        return out
    n0 = host["n_atoms"]
    sp = (out.get("nep") or {}).get("species")
    if sp:
        need = [e for e in list(host["species"]) + [st_.get("dopant")] if e]
        miss = [e for e in dict.fromkeys(need) if e not in sp]
        if miss:
            problems.append(f"the NEP model has no {', '.join(miss)} (it covers {', '.join(sp)})")
    n_host = host["species"].get(st_["host"])
    out["host"] = {"n_atoms": n0, "n_host": n_host, "species": host["species"], "repeat": host.get("repeat_guess")}
    if not n_host:
        return out
    try:
        comps = comp_list(st_["compositions"], n_host)
    except (ValueError, TypeError):
        comps = []
    comps = [n for n in comps if 0 <= n <= n_host]
    lat_host, _ = _read_cell(Path(st_["poscar"]).expanduser())
    rep = host.get("repeat_guess")
    cells = {}
    if mp.get("include_parent", True):
        cells["parent"] = lat_host
    for t in mp.get("templates") or []:
        label, p = R.split_template(t)
        try:
            tl, n = _read_cell(Path(p).expanduser())
        except (OSError, ValueError, IndexError):
            continue
        cells[label] = [[x * (rep[i] if rep else 1) for x in tl[i]] for i in range(3)] if rep else lat_host
    units = {}
    for ph, lat in cells.items():
        sc = phonon_sc(fe, lat)
        units[ph] = sc[0] * sc[1] * sc[2]          # n_atoms = n0 for every phase (same cations, same repeat)
    per_phase = {ph: 0 for ph in cells}
    tot_u = 0.0
    rows = []
    for n in comps:
        k = min(int(st_["max_orderings"]), math.comb(n_host, n))
        rows.append({"n": n, "pct": round(100 * n / n_host, 2), "orderings": k, "possible": math.comb(n_host, n)})
        for ph in cells:
            per_phase[ph] += k
            tot_u += k * units[ph]
    sec = float(cl["sec_per_structure"])
    threads = max(1, int(fe["threads"]))
    workers = max(1, int(cl["cores_per_node"]) // threads)
    wall = parse_walltime(cl["walltime"]) or 0
    cap = workers * wall * float(cl["fill"]) / sec if sec > 0 else 0
    n_chunks = int(math.ceil(tot_u / cap)) if cap > 0 and tot_u else 0
    out.update({"structures": sum(per_phase.values()), "per_phase": per_phase, "phases": list(cells),
                "units": {k: v for k, v in units.items()}, "compositions": rows,
                "core_h": tot_u * sec * threads / 3600.0, "node_h": tot_u * sec * threads / 3600.0 / max(1, int(cl["cores_per_node"])),
                "chunks": n_chunks, "workers": workers})
    if n_chunks:
        per = tot_u * sec / workers / n_chunks
        rounds = int(math.ceil(n_chunks / max(1, int(cl["max_concurrent"]))))
        biggest = max(units.values()) * sec / 1.0       # one structure of the most expensive phase (s, one worker)
        out.update({"chunk_s": per, "rounds": rounds, "wall_total_s": rounds * per, "largest_structure_s": biggest,
                    "array_too_big": n_chunks > int(cl["max_array"])})
    return out
