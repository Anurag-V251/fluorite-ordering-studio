"""
scan.py - read a results folder and report what is there, without changing anything.

Layouts recognised (derived from what campaign.py, rerun_flagged.py and batch_free_energy.py write):

  campaign workdir                        (campaign.py; the folder given as [paths] workdir)
    <Dop>_<pct>/                          one per composition, e.g. Ce_18.75 (name rule of enumerate_substitutions.py)
        summary.csv                       orderings (enumerate_substitutions.py)
        POSCAR_<id>_WC<alpha>             parent orderings      -> phase "parent"
        prepared.json                     written when enumerate + map + verify finished
        mapped/summary.csv                batch_map.py: one row per ordering x target (status, rms shift, strain, variants)
        mapped/verify.csv                 compare_ordering.py verdict of every mapped file
        mapped/<phase>/POSCAR_<phase>_... mapped structures    -> phase <phase>
        free_energy_<T>K.csv              per-composition CSV written by merge
    chunks/plan.json, chunk_*.txt         plan (which structure is in which chunk)
    results/chunk_*.csv                   batch_free_energy.py rows (chunk_rerunRRR_K.csv from rerun_flagged.py)
    relaxed/<rel path>                    relaxed structures; rerun_RRR/relaxed/<rel path> for rerun rounds
    free_energy_<T>K_all.csv              merged CSV (plot scripts read this)
    merge_report_<T>K.txt, campaign_config.toml, logs/, slurm/

  campaign start folder                   a folder with campaign.toml (-> its workdir is scanned)

  plain folder                            Ce_xx/ folders processed by hand with batch_map.py / batch_free_energy.py:
                                          Ce_xx/free_energy_<T>K.csv and Ce_xx/relaxed/...

Rows are combined like campaign.py merge does: per structure the newest row, an ok row beats a FAILED one.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import statistics
import time
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:          # Python < 3.11
    try:
        import tomli as tomllib      # type: ignore
    except ModuleNotFoundError:
        tomllib = None

COMP_RE = re.compile(r"^(?P<el>[A-Z][a-z]?)_(?P<pct>\d+(?:\.\d+)?)$")
FE_RE = re.compile(r"^free_energy_(?P<T>[\d.]+)K\.csv$")
ALL_RE = re.compile(r"^free_energy_(?P<T>[\d.]+)K_all\.csv$")
TITLE_RE = re.compile(r"\|\s*([A-Z][a-z]?)\s+(\d+)\s*/\s*(\d+)\s*=")


def read_csv(p):
    p = Path(p)
    if not p.is_file():
        return []
    try:
        with p.open(newline="") as f:
            return list(csv.DictReader(f))
    except (OSError, csv.Error, UnicodeDecodeError):
        return []


RESULT_KEYS = ("file", "status", "T_K", "phase", "n_imaginary_modes", "relax_converged", "F_total_eV_atom")


def iter_result_rows(p):
    """rows of a result CSV, streamed, keeping only the columns the scan uses (130k structures stay small)"""
    p = Path(p)
    try:
        with p.open(newline="") as f:
            for r in csv.DictReader(f):
                yield {k: r[k] for k in RESULT_KEYS if r.get(k) not in (None, "")}
    except (OSError, csv.Error, UnicodeDecodeError):
        return


def add_best(best, r):
    """fold one row into {file: row}: latest row per structure, ok beats FAILED (campaign.py merge rule)"""
    f = r.get("file")
    if f and (f not in best or r.get("status") == "ok" or best[f].get("status") != "ok"):
        best[f] = r


def fnum(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def load_toml(p):
    if tomllib is None:
        return None
    try:
        with open(p, "rb") as f:
            return tomllib.load(f)
    except Exception:            # noqa: BLE001 - a broken file is reported, not fatal
        return None


def find_config(root: Path):
    """campaign config that belongs to this folder: (config dict, config path, workdir)"""
    # 1. root is a workdir: campaign.py prepare copies the config there
    rec = root / "campaign_config.toml"
    if rec.is_file():
        return load_toml(rec), rec, root
    # 2. root holds campaign config(s) whose workdir exists (several -> the caller lists them)
    found = []
    for p in sorted(root.glob("*.toml")):
        c = load_toml(p)
        if not c or "structure" not in c or "mapping" not in c:
            continue
        wd = (root / c.get("paths", {}).get("workdir", "campaign")).resolve()
        if wd.is_dir():
            found.append((c, p, wd))
    if len(found) == 1:
        return found[0]
    if len(found) > 1:
        return None, [str(f[2]) for f in found], root
    # 3. config present but workdir not created yet
    for p in sorted(root.glob("*.toml")):
        c = load_toml(p)
        if c and "structure" in c and "mapping" in c:
            return c, p, (root / c.get("paths", {}).get("workdir", "campaign")).resolve()
    return None, None, root


def list_files(d: Path, prefix="POSCAR"):
    try:
        return sorted(e.name for e in os.scandir(d) if e.is_file() and e.name.startswith(prefix))
    except OSError:
        return []


def walk_rel(base: Path):
    """set of file paths (relative to base, '/'-separated) below base"""
    out = set()
    if not base.is_dir():
        return out
    for dp, _dn, fn in os.walk(base):
        rel = os.path.relpath(dp, base)
        for f in fn:
            out.add(f if rel == "." else f"{rel}/{f}".replace(os.sep, "/"))
    return out


def composition_dirs(work: Path):
    out = []
    try:
        entries = sorted(os.scandir(work), key=lambda e: e.name)
    except OSError:
        return out
    for e in entries:
        m = COMP_RE.match(e.name)
        if not (m and e.is_dir()):
            continue
        d = Path(e.path)
        if (d / "summary.csv").is_file() or list_files(d) or (d / "mapped").is_dir() or any(
                FE_RE.match(x) for x in os.listdir(d)):
            out.append((d, m["el"], float(m["pct"])))
    out.sort(key=lambda x: (x[1], x[2]))
    return out


def comp_n(d: Path, parents):
    """(n_dopant, n_sites) from prepared.json or the POSCAR title line 'Ce 6/32 = ...'"""
    pj = d / "prepared.json"
    n = N = None
    if pj.is_file():
        try:
            n = json.loads(pj.read_text()).get("n")
        except (OSError, ValueError):
            pass
    for f in parents[:1]:
        try:
            first = (d / f).open().readline()
        except OSError:
            break
        m = TITLE_RE.search(first)
        if m:
            n = int(m[2]) if n is None else n
            N = int(m[3])
    return n, N


def structure_universe(d: Path, cname: str, include_parent=True, labels=None):
    """{rel path to workdir: phase} of the structures that belong to one composition, the way campaign.py
    lists them: parents from summary.csv, mapped files from mapped/summary.csv rows with status ok and a target
    of the config (labels); directory listings only where those CSVs are missing (plain layouts)"""
    uni = {}
    summ = read_csv(d / "summary.csv")
    if summ and "file" in summ[0]:
        parents = [r["file"] for r in summ if (d / r["file"]).is_file()]
    else:
        parents = list_files(d)
    if include_parent:
        for f in parents:
            uni[f"{cname}/{f}"] = "parent"
    md = d / "mapped"
    msum = read_csv(md / "summary.csv")
    if msum and "out_file" in msum[0]:
        variants = {}                              # per phase dir: OUT -> [OUT_v1, OUT_v2, ...] (one listing)
        for r in msum:
            t = r.get("target", "")
            if r.get("status") != "ok" or (labels and t not in labels):
                continue
            if t not in variants:
                vv = {}
                for x in list_files(md / t, ""):
                    m = re.fullmatch(r"(.*)_v(\d+)", x)
                    if m:
                        vv.setdefault(m[1], []).append((int(m[2]), x))
                variants[t] = {k: [x for _, x in sorted(v)] for k, v in vv.items()}
            vs = variants[t].get(r["out_file"], [])
            for n in (vs or [r["out_file"]]):     # --all-variants: OUT == OUT_v1, only the _v files are computed
                uni[f"{cname}/mapped/{t}/{n}"] = t
        return uni, parents
    if md.is_dir():
        for e in sorted(os.scandir(md), key=lambda e: e.name):
            if not e.is_dir():
                continue
            names = list_files(Path(e.path))
            have_v = {re.sub(r"_v\d+$", "", n) for n in names if re.search(r"_v\d+$", n)}
            for n in names:
                if n in have_v:             # --all-variants: OUT == OUT_v1, only the _v files are computed
                    continue
                uni[f"{cname}/mapped/{e.name}/{n}"] = e.name
    return uni, parents


def combine_rows(rows):
    """latest row per structure; ok beats FAILED (same rule as campaign.py merge)"""
    best = {}
    for r in rows:
        add_best(best, r)
    return best


def row_flags(r):
    v = fnum(r.get("n_imaginary_modes"), 0)
    imag = math.isfinite(v) and v > 0
    unconv = str(r.get("relax_converged", "")).strip() == "False"
    return imag, unconv


def _stats(vals):
    v = [x for x in vals if math.isfinite(x)]
    if not v:
        return None
    return {"min": min(v), "median": statistics.median(v), "max": max(v)}


def scan(path: str, max_rows=None, relaxed_walk=True) -> dict:
    t0 = time.time()
    root = Path(path).expanduser().resolve()
    if root.is_file() and root.suffix == ".toml":         # a campaign config: scan its workdir
        c = load_toml(root)
        if not c or "structure" not in c:
            return {"ok": False, "error": f"{root} is not a campaign config"}
        cfg, cfg_path = c, root
        work = (root.parent / c.get("paths", {}).get("workdir", "campaign")).resolve()
        root = root.parent
    elif not root.is_dir():
        return {"ok": False, "error": f"{root} is not a folder"}
    else:
        cfg, cfg_path, work = find_config(root)
    warnings = []
    if isinstance(cfg_path, list):
        return {"ok": True, "root": str(root), "layout": "none", "candidates": cfg_path,
                "warnings": ["several campaign configs here; pick one of their workdirs"]}
    if cfg_path is not None and not work.is_dir():
        warnings.append(f"config {cfg_path.name} found, but its workdir {work} does not exist yet "
                        "(run prepare first)")
    include_parent = True
    templates = []
    T_cfg = None
    if cfg:
        include_parent = bool(cfg.get("mapping", {}).get("include_parent", True))
        for spec in cfg.get("mapping", {}).get("templates", []) or []:
            templates.append(spec.split("=", 1)[0] if "=" in spec else re.sub(r"^POSCAR[_.\-]?", "", Path(spec).name))
        T_cfg = cfg.get("free_energy", {}).get("T", 300.0)

    comps = composition_dirs(work) if work.is_dir() else []
    layout = "empty"
    if (work / "results").is_dir() or (work / "chunks").is_dir() or cfg:
        layout = "campaign"
    elif comps:
        layout = "plain"
    if not comps and layout == "empty":
        # maybe one level down (user pointed at a parent of several workdirs)
        subs = [str(p) for p in sorted(root.iterdir()) if p.is_dir() and (
            (p / "campaign_config.toml").is_file() or (p / "results").is_dir() or composition_dirs(p))]
        return {"ok": True, "root": str(root), "layout": "none", "candidates": subs,
                "warnings": ["no composition folders (e.g. Ce_18.75/) found here" +
                             ("; candidates below" if subs else "")]}

    # ---------------- structures that should exist
    universe, comp_info = {}, []
    for d, el, pct in comps:
        uni, parents = structure_universe(d, d.name, include_parent, set(templates) or None)
        universe.update(uni)
        n, N = comp_n(d, parents)
        prep = None
        if (d / "prepared.json").is_file():
            try:
                prep = json.loads((d / "prepared.json").read_text())
            except (OSError, ValueError):
                prep = {}
        mapped = read_csv(d / "mapped" / "summary.csv")
        verify = read_csv(d / "mapped" / "verify.csv")
        comp_info.append({"dir": d.name, "dopant": el, "pct": pct, "n": n, "N": N,
                          "orderings": len(parents), "prepared": prep is not None, "prepared_info": prep,
                          "mapped_rows": mapped, "verify_rows": verify})

    # ---------------- results
    best_by_T = {}                 # T -> {file: best row}, folded while streaming (no full row list in memory)
    sources = []
    chunk_best = {}                # chunk csv name -> {file: best row} (any T), for the chunk progress
    if layout == "campaign" and (work / "results").is_dir():
        res = sorted((work / "results").glob("chunk_*.csv"))
        sources = [str(p.relative_to(work)) for p in res]
        for p in res:
            cb = chunk_best[p.name] = {}
            for r in iter_result_rows(p):
                add_best(cb, r)
                add_best(best_by_T.setdefault(f"{fnum(r.get('T_K')):g}", {}), r)
    else:
        for d, _, _ in comps:
            for f in sorted(os.listdir(d)):
                m = FE_RE.match(f)
                if not m:
                    continue
                sources.append(f"{d.name}/{f}")
                for r in iter_result_rows(d / f):
                    if not r.get("file"):
                        continue
                    r["file"] = f"{d.name}/{r['file']}"
                    add_best(best_by_T.setdefault(f"{fnum(m['T']):g}", {}), r)
        if not best_by_T:          # only a merged CSV
            for p in sorted(work.glob("free_energy_*K_all.csv")):
                sources.append(p.name)
                for r in iter_result_rows(p):
                    add_best(best_by_T.setdefault(f"{fnum(r.get('T_K')):g}", {}), r)
    rows_by_T = best_by_T
    temps = sorted(rows_by_T, key=fnum)
    if T_cfg is not None and f"{float(T_cfg):g}" not in temps:
        temps.append(f"{float(T_cfg):g}")

    # ---------------- relaxed structures on disk
    relaxed = set()
    if relaxed_walk:
        relaxed = walk_rel(work / "relaxed")
        for rd in sorted(work.glob("rerun_[0-9][0-9][0-9]")):
            relaxed |= walk_rel(rd / "relaxed")
        for d, _, _ in comps:                              # plain layout: Ce_xx/relaxed/...
            if (d / "relaxed").is_dir():
                relaxed |= {f"{d.name}/{x}" for x in walk_rel(d / "relaxed")}

    # ---------------- plan
    planned, chunks = set(), None
    pj = work / "chunks" / "plan.json"
    if pj.is_file():
        try:
            plan = json.loads(pj.read_text())
            chunk_list = []
            for ch in plan.get("chunks", []):
                files = (work / "chunks" / ch["file"]).read_text().split() if (work / "chunks" / ch["file"]).is_file() else []
                planned.update(files)
                rr = chunk_best.get(f"chunk_{int(ch['id']):05d}.csv", {})
                ok = sum(1 for f in files if rr.get(f, {}).get("status") == "ok")
                fl = sum(1 for f in files if rr.get(f, {}).get("status") == "FAILED")
                chunk_list.append({"id": ch["id"], "n": len(files), "ok": ok, "failed": fl,
                                   "est_s": ch.get("est_s")})
            chunks = {"n": len(chunk_list), "complete": sum(1 for c in chunk_list if c["ok"] == c["n"]),
                      "started": sum(1 for c in chunk_list if c["ok"] + c["failed"] > 0),
                      "workers_per_node": plan.get("workers_per_node"), "walltime": plan.get("walltime"),
                      "sec_per_structure": plan.get("sec_per_structure"), "list": chunk_list}
        except (OSError, ValueError, KeyError) as e:
            warnings.append(f"chunks/plan.json unreadable: {e}")

    # ---------------- per composition x phase table, per temperature
    def phase_of(rel, row=None):
        if rel in universe:
            return universe[rel]
        if row and row.get("phase"):
            return row["phase"]
        parts = rel.split("/")
        return parts[2] if len(parts) > 3 and parts[1] == "mapped" else "parent"

    phases_seen = []
    if include_parent and any(c["orderings"] for c in comp_info):
        phases_seen.append("parent")
    for t in templates:
        if t not in phases_seen:
            phases_seen.append(t)
    for ph in universe.values():
        if ph not in phases_seen:
            phases_seen.append(ph)

    tables = {}
    for T in temps or ["-"]:
        best = best_by_T.get(T, {})
        cells = {}
        keys = set(universe) | set(best) | planned
        for rel in keys:
            cdir = rel.split("/", 1)[0]
            ph = phase_of(rel, best.get(rel))
            if ph not in phases_seen:
                phases_seen.append(ph)
            c = cells.setdefault((cdir, ph), {"structures": 0, "planned": 0, "relaxed": 0, "ok": 0, "failed": 0,
                                              "imaginary": 0, "unconverged": 0, "clean": 0, "not_run": 0,
                                              "F": [], "F_clean": []})
            c["structures"] += 1
            c["planned"] += rel in planned
            c["relaxed"] += rel in relaxed
            r = best.get(rel)
            if r is None:
                c["not_run"] += 1
            elif r.get("status") == "ok":
                c["ok"] += 1
                im, un = row_flags(r)
                c["imaginary"] += im
                c["unconverged"] += un
                F = fnum(r.get("F_total_eV_atom"))
                c["F"].append(F)
                if not (im or un):
                    c["clean"] += 1
                    c["F_clean"].append(F)
            else:
                c["failed"] += 1
        # lowest F per composition (any phase, unflagged) -> dF of each phase's lowest
        low = {}
        for (cdir, ph), c in cells.items():
            v = [x for x in c["F_clean"] if math.isfinite(x)]
            if v:
                low[cdir] = min(low.get(cdir, v[0]), min(v))
        out = []
        for (cdir, ph), c in cells.items():
            fc = [x for x in c["F_clean"] if math.isfinite(x)]
            item = {k: v for k, v in c.items() if k not in ("F", "F_clean")}
            item.update({"comp": cdir, "phase": ph,
                         "F_stats": _stats(c["F"]),
                         "dF_lowest_meV": (1000 * (min(fc) - low[cdir])) if fc and cdir in low else None,
                         "dF_median_meV": (1000 * (statistics.median(fc) - low[cdir])) if fc and cdir in low else None})
            out.append(item)
        tables[T] = out

    comps_out = []
    for c in comp_info:
        mr = c.pop("mapped_rows")
        vr = c.pop("verify_rows")
        mp = {}
        for r in mr:
            m = mp.setdefault(r.get("target", "?"), {"ok": 0, "failed": 0, "rms": [], "S": [], "tied": []})
            if r.get("status") == "ok":
                m["ok"] += 1
                m["rms"].append(fnum(r.get("rms_shift_A")))
                m["S"].append(fnum(r.get("strain_S")))
                m["tied"].append(fnum(r.get("tied_variants"), 1))
            else:
                m["failed"] += 1
        mapping = {}
        for t, m in mp.items():
            mapping[t] = {"ok": m["ok"], "failed": m["failed"], "rms_shift_A": _stats(m["rms"]),
                          "strain_S": _stats(m["S"]), "tied_variants": _stats(m["tied"])}
        ver = {}
        for r in vr:
            v = ver.setdefault(r.get("target", "?"), {"identical": 0, "other": 0})
            if r.get("verdict") == "IDENTICAL":
                v["identical"] += 1
            else:
                v["other"] += 1
        c["mapping"] = mapping
        c["verify"] = ver
        comps_out.append(c)

    merged = [p.name for p in sorted(work.glob("free_energy_*K_all.csv"))] if work.is_dir() else []
    per_comp_csv = {}
    for d, _, _ in comps:
        per_comp_csv[d.name] = sorted(f for f in os.listdir(d) if FE_RE.match(f))
    reports = {}
    for p in sorted(work.glob("merge_report_*K.txt")) if work.is_dir() else []:
        try:
            reports[p.name] = p.read_text()[:20000]
        except OSError:
            pass
    reruns = []
    for rd in sorted(work.glob("rerun_[0-9][0-9][0-9]")) if work.is_dir() else []:
        info = {"dir": rd.name}
        if (rd / "rerun.json").is_file():
            try:
                info.update(json.loads((rd / "rerun.json").read_text()))
            except (OSError, ValueError):
                pass
        info["listed"] = sum(len((p.read_text().split())) for p in rd.glob("list_*.txt"))
        reruns.append({k: v for k, v in info.items() if isinstance(v, (str, int, float, bool, type(None)))})

    totals = {}
    for T, tab in tables.items():
        tt = {k: sum(x[k] for x in tab) for k in ("structures", "planned", "relaxed", "ok", "failed",
                                                   "imaginary", "unconverged", "clean", "not_run")}
        totals[T] = tt
    dopants = sorted({c["dopant"] for c in comps_out})
    if dopants and len(dopants) > 1:
        warnings.append(f"several dopant prefixes in one folder: {dopants}")
    cdop = (cfg or {}).get("structure", {}).get("dopant", "Ce")
    if cfg and dopants and cdop not in dopants:
        warnings.append(f"config dopant {cdop!r} but folders are {dopants}")
    bad_ver = sum(v["other"] for c in comps_out for v in c["verify"].values())
    if bad_ver:
        warnings.append(f"{bad_ver} mapped structures are NOT identical to their parent ordering (verify.csv)")
    run_name = None
    rr = (Path(cfg_path).parent if cfg_path and not isinstance(cfg_path, list) else work.parent) / "run_record.json"
    if rr.is_file():
        try:
            run_name = json.loads(rr.read_text()).get("name") or rr.parent.name
        except (OSError, ValueError):
            pass
    return {
        "ok": True, "root": str(root), "workdir": str(work), "layout": layout, "run_name": run_name,
        "config_path": str(cfg_path) if cfg_path else None, "config": cfg,
        "phases": phases_seen, "temperatures": temps, "compositions": comps_out,
        "tables": tables, "totals": totals, "chunks": chunks, "result_sources": sources[:200],
        "n_result_sources": len(sources), "merged_csv": merged, "per_comp_csv": per_comp_csv,
        "merge_reports": reports, "reruns": reruns, "warnings": warnings,
        "has_logs": (work / "logs").is_dir(), "has_slurm": (work / "slurm").is_dir(),
        "scan_seconds": round(time.time() - t0, 2), "relaxed_counted": relaxed_walk,
    }


if __name__ == "__main__":                    # python -m studio.scan FOLDER  -> JSON summary
    import sys
    out = scan(sys.argv[1] if len(sys.argv) > 1 else ".")
    out.pop("config", None)
    print(json.dumps(out, indent=1, default=str)[:20000])
