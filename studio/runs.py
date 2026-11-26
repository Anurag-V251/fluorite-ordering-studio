"""
runs.py - create a reproducible run folder from the settings on the "New run" page, and read it back.

A run folder:
    <run>/campaign.toml        the config campaign.py reads (all paths relative to this folder)
    <run>/run_record.json      everything needed to reproduce: settings, md5 of every input and script,
                               package versions, host, user, and the history of every command started from the app
    <run>/inputs/              copies of the host POSCAR, templates and NEP model (unless "copy inputs" is off)
    <run>/scripts/             frozen copy of the scripts used (unless "freeze scripts" is off)
    <run>/campaign/            workdir written by campaign.py
"""
from __future__ import annotations

import getpass
import hashlib
import json
import platform
import re
import shutil
import socket
import sys
import time
import threading
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:                     # Python < 3.11
    try:
        import tomli as tomllib  # type: ignore
    except ModuleNotFoundError:
        sys.exit("Python < 3.11 needs the package tomli:  pip install tomli")

APP_VERSION = "1.2.0"
HERE = Path(__file__).resolve().parent

DEFAULTS = {   # same keys as campaign.py DEFAULTS (values are the app's start values); only these are written
    "structure": {"poscar": "", "host": "Zr", "dopant": "Ce", "compositions": "all", "max_orderings": 1000, "seed": 1},
    "paths": {"scripts": "scripts", "workdir": "campaign"},
    "mapping": {"templates": [], "map_args": "", "include_parent": True, "verify": True, "jobs": 8},
    "free_energy": {"nep": "", "T": 300.0, "relax": True, "fmax": 1e-3, "steps": 2000, "fixed_cell": False,
                    "min_length": 10.0, "supercell": [], "mesh": 100, "threads": 1},
    "cluster": {"cores_per_node": 64, "walltime": "24:00:00", "fill": 0.8, "sec_per_structure": 300.0,
                "max_array": 1000, "max_concurrent": 50, "prepare_walltime": "04:00:00",
                "merge_walltime": "01:00:00", "partition": "", "account": "", "qos": "", "mem": "",
                "sbatch_extra": [], "python": "python3", "setup": ""},
}


def _user():
    try:
        return getpass.getuser()
    except Exception:
        import os
        return str(os.getuid())


RECORD_LOCK = threading.Lock()
CREATE_LOCK = threading.Lock()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def versions():
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    try:
        from importlib.metadata import version, PackageNotFoundError
    except ImportError:
        return out
    for pkg in ("numpy", "scipy", "matplotlib", "ase", "calorine", "phonopy", "pymatgen", "icet"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


# ---------------------------------------------------------------------- POSCAR peek (validation on the page)
def peek_poscar(path):
    p = Path(path).expanduser()
    try:
        L = p.read_text().split("\n")
        sc = float(L[1].split()[0])
        lat = [[float(v) for v in L[k].split()[:3]] for k in (2, 3, 4)]
        names = L[5].split()
        if names and all(re.fullmatch(r"\d+", x) for x in names):
            return {"ok": False, "error": "VASP4 POSCAR (no element line): add the element names on line 6"}
        counts = [int(x) for x in L[6].split()]
    except (OSError, ValueError, IndexError) as e:
        return {"ok": False, "error": f"not a readable VASP5 POSCAR: {e}"}
    import math
    lens = [math.sqrt(sum(x * x for x in v)) * (sc if sc > 0 else 1) for v in lat]
    off = max(abs(lat[i][j]) for i in range(3) for j in range(3) if i != j)
    comp = {}
    for n, c in zip(names, counts):            # repeated element blocks ("Zr O Zr O") add up
        comp[n] = comp.get(n, 0) + c
    n_cat = sum(c for e, c in comp.items() if e not in ("O", "F", "S", "N"))
    rep = None
    if n_cat and n_cat % 4 == 0:
        guess = [max(1, round(x / 5.1)) for x in lens]
        if guess[0] * guess[1] * guess[2] * 4 == n_cat:
            rep = guess
    return {"ok": True, "title": L[0].strip(), "species": comp, "n_atoms": sum(counts),
            "lengths": [round(x, 4) for x in lens], "off_diagonal_max": round(off, 4),
            "repeat_guess": rep, "md5": md5(p)}


# ---------------------------------------------------------------------- TOML writer (only what campaign.toml needs)
def _tv(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_tv(x) for x in v) + "]"
    s = str(v).encode("utf-8", "replace").decode("utf-8")          # no lone surrogates
    s = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]", "", s)                  # TOML forbids these control characters
    if "\n" in s:
        return '"""\n' + s.replace("\r", "").replace("\\", "\\\\").replace('"""', '\\"""').strip("\n") + '\n"""'
    return json.dumps(s, ensure_ascii=False)


def split_template(t):
    """'label=path' -> (label, path); a bare path -> (label from the file name, path).  Only the first '=' counts,
    and only when the part before it is a valid label (so /data/a=b/POSCAR_x is a bare path)."""
    t = str(t).strip()
    if "=" in t:
        label, p = t.split("=", 1)
        if re.fullmatch(r"[\w.+\-]+", label):
            return label, p
    return re.sub(r"^POSCAR[_.\-]?", "", Path(t).name), t


def to_toml(cfg, header=""):
    out = [header.rstrip()] if header else []
    for sec, d in cfg.items():
        out.append(f"\n[{sec}]")
        for k, v in d.items():
            out.append(f"{k:<17} = {_tv(v)}")
    return "\n".join(out) + "\n"


def coerce(cfg_in):
    """settings from the page -> campaign config dict with the right types (unknown keys dropped)"""
    out = {}
    for sec, d in DEFAULTS.items():
        got = (cfg_in or {}).get(sec, {})
        o = {}
        for k, dv in d.items():
            v = got.get(k, dv)
            if isinstance(dv, bool):
                v = v if isinstance(v, bool) else str(v).lower() in ("1", "true", "yes", "on")
            elif isinstance(dv, int) and not isinstance(dv, bool):
                sv = str(v).strip()
                try:
                    v = int(sv) if re.fullmatch(r"-?\d+", sv) else int(float(sv)) if sv else dv
                except (ValueError, OverflowError):
                    raise ValueError(f"[{sec}] {k} = {v!r} is not a whole number")
            elif isinstance(dv, float):
                try:
                    v = float(v) if str(v).strip() != "" else dv
                except ValueError:
                    raise ValueError(f"[{sec}] {k} = {v!r} is not a number")
            elif isinstance(dv, list):
                if isinstance(v, str):
                    v = [x for x in re.split(r"[\s,]+" if k != "sbatch_extra" else r"\n+", v.strip()) if x]
                if k == "supercell":
                    v = [int(x) for x in v]
            o[k] = v
        out[sec] = o
    c = out["structure"]["compositions"]
    if isinstance(c, (list, tuple)):           # a list from the API: same rules as typed text
        c = " ".join(str(x) for x in c)
    if isinstance(c, str):
        cs = c.strip()
        if cs == "all" or re.fullmatch(r"\d+\s*:\s*\d+(\s*:\s*\d+)?", cs):
            out["structure"]["compositions"] = cs.replace(" ", "")
        elif re.fullmatch(r"[\d\s,\[\]]+", cs) and re.search(r"\d", cs):
            out["structure"]["compositions"] = sorted({int(x) for x in re.split(r"[\s,\[\]]+", cs) if x})
        else:
            raise ValueError(f"compositions = {c!r}: use 'all', 'start:stop[:step]' (inclusive) or a list of "
                             "dopant counts such as '0, 6, 16, 32' or '0 6 16 32'")
    return out


def validate(cfg):
    """the checks campaign.py dry-run would do, at create time, with plain messages"""
    import math
    s, fe = cfg["structure"], cfg["free_energy"]
    for sec, d in cfg.items():
        for k, v in d.items():
            if isinstance(v, float) and not math.isfinite(v):
                raise ValueError(f"[{sec}] {k} = {v}: not a finite number")
    if not s["poscar"]:
        raise ValueError("host POSCAR (starting structure) is required")
    pk = peek_poscar(s["poscar"])
    if not pk["ok"]:
        raise ValueError(f"host POSCAR {s['poscar']}: {pk['error']}")
    if s["host"] not in pk["species"]:
        raise ValueError(f"host element {s['host']} is not in the host POSCAR (it has {', '.join(pk['species'])})")
    if not re.fullmatch(r"[A-Z][a-z]?", str(s["dopant"])):
        raise ValueError(f"dopant {s['dopant']!r} is not an element symbol")
    N = pk["species"][s["host"]]
    c = s["compositions"]
    if isinstance(c, list):
        bad = [n for n in c if not 0 <= n <= N]
    elif c == "all":
        bad = []
    else:
        a, b, *st = [int(x) for x in c.split(":")]
        if a > b:
            raise ValueError(f"compositions {c}: start {a} is larger than stop {b}")
        if st and st[0] < 1:
            raise ValueError(f"compositions {c}: the step must be at least 1")
        bad = [n for n in (a, b) if not 0 <= n <= N]
    if bad:
        raise ValueError(f"compositions {bad} outside 0..{N} (number of {s['host']} sites in the host POSCAR)")
    if s["max_orderings"] < 1:
        raise ValueError("max orderings must be at least 1")
    for k in ("T", "fmax", "min_length"):
        if fe[k] <= 0:
            raise ValueError(f"{k} must be positive")
    if fe["supercell"] and (len(fe["supercell"]) != 3 or min(fe["supercell"]) < 1):
        raise ValueError("phonon supercell: three whole numbers >= 1, or empty for automatic")
    if cfg["cluster"]["cores_per_node"] < 1 or fe["threads"] < 1 or fe["steps"] < 1:
        raise ValueError("cores_per_node, threads and steps must be >= 1")
    for k in ("walltime", "prepare_walltime", "merge_walltime"):
        w = cfg["cluster"].get(k)
        if k == "walltime" or w not in (None, ""):
            if not re.fullmatch(r"(\d+-)?\d+(:\d+){0,2}", str(w)):
                raise ValueError(f"{k} {w!r}: use HH:MM:SS or D-HH:MM:SS")
    wd = str(cfg["paths"].get("workdir", "campaign")).strip()
    if not wd or Path(wd).is_absolute() or ".." in Path(wd).parts or wd in (".", "./") or \
            Path(wd).parts[0] in ("inputs", "scripts", "slurm"):
        raise ValueError(f"workdir {wd!r}: a plain sub-folder name inside the run folder, such as 'campaign'")
    tpl = cfg["mapping"]["templates"]
    if not tpl:
        raise ValueError("at least one target phase template is needed")
    seen = set()
    for t in tpl:
        label, p = split_template(t)
        if not re.fullmatch(r"[\w.+\-]+", label):
            raise ValueError(f"phase label {label!r}: letters, digits, _ . + - only")
        if label in seen:
            raise ValueError(f"phase label {label!r} used twice")
        seen.add(label)
        pk = peek_poscar(Path(p).expanduser())
        if not pk["ok"]:
            raise ValueError(f"template {label}: {p}: {pk['error']}")
        if pk["n_atoms"] != 12:
            raise ValueError(f"template {label}: {p} has {pk['n_atoms']} atoms; map_fluorite.py needs the 12-atom "
                             "conventional cell (4 cations + 8 anions) of the target phase")
    if fe["nep"] and not Path(fe["nep"]).expanduser().is_file():
        raise ValueError(f"NEP model {fe['nep']} not found")


# ---------------------------------------------------------------------- create / read runs
def create_run(spec, app_scripts: Path):
    """spec = {"folder", "name", "starting_phase", "notes", "copy_inputs", "freeze_scripts", "config": {...}}"""
    folder = Path(spec["folder"]).expanduser().resolve()
    if not isinstance(spec.get("config", {}), dict) or not all(isinstance(v, dict) for v in spec.get("config", {}).values()):
        raise ValueError("config must be a table of sections")
    cfg = coerce(spec.get("config", {}))
    validate(cfg)
    with CREATE_LOCK:                           # two quick clicks: the second waits and then finds the run
        made = False
        if not folder.exists():
            folder.mkdir(parents=True)
            made = True
        try:
            return _create_run(spec, folder, cfg, app_scripts)
        except Exception:
            if made:                            # nothing half-written left behind: the same request can be retried
                shutil.rmtree(folder, ignore_errors=True)
            raise


def _create_run(spec, folder, cfg, app_scripts):
    if (folder / "campaign.toml").exists() and not spec.get("overwrite"):
        raise ValueError(f"{folder}/campaign.toml exists: pick a new folder, or tick 'overwrite settings'")
    is_run = (folder / "run_record.json").is_file()
    wds = [folder / cfg["paths"]["workdir"]]
    if (folder / "campaign.toml").is_file():    # the workdir of the settings being replaced counts too
        try:
            with (folder / "campaign.toml").open("rb") as f:
                wds.append(folder / tomllib.load(f).get("paths", {}).get("workdir", "campaign"))
        except Exception:                       # noqa: BLE001
            pass
    if spec.get("overwrite"):
        for wd in wds:
            if any((wd / "results").glob("chunk_*.csv")):
                raise ValueError(f"{wd} already has results: rewriting the settings would make the record wrong "
                                 "for them. Use 'copy settings to a new run' instead.")
            if (wd / "chunks" / "plan.json").is_file():
                raise ValueError(f"{wd} is already planned: rewriting the settings now would not match the planned "
                                 "chunks. Use 'copy settings to a new run' instead.")
    if spec.get("freeze_scripts", True) and (folder / "scripts").exists() and not is_run:
        raise ValueError(f"{folder}/scripts exists and was not made by this app: choose an empty or new folder")
    if folder.is_dir() and any(folder.iterdir()) and not is_run and not spec.get("overwrite"):
        raise ValueError(f"{folder} is not empty: choose a new folder (or an empty one)")
    folder.mkdir(parents=True, exist_ok=True)
    old_rec = folder / "run_record.json"
    if old_rec.is_file():                       # keep the earlier record next to the new one
        shutil.copy2(old_rec, folder / f"run_record.{time.strftime('%Y%m%d-%H%M%S')}.json")
    inputs = {}

    def take(src, sub=""):
        src = Path(src).expanduser().resolve()
        if not src.is_file():
            raise ValueError(f"file {src} not found")
        if spec.get("copy_inputs", True):
            dst = folder / "inputs" / sub / src.name
            if str(dst.relative_to(folder)) in inputs and inputs[str(dst.relative_to(folder))]["source"] != str(src):
                raise ValueError(f"two inputs would be copied to {dst.relative_to(folder)}")
            dst.parent.mkdir(parents=True, exist_ok=True)
            if src != dst:
                shutil.copy2(src, dst)
            rel = str(dst.relative_to(folder))
        else:
            rel = str(src)
        inputs[rel] = {"source": str(src), "md5": md5(src)}
        return rel

    s, m, fe = cfg["structure"], cfg["mapping"], cfg["free_energy"]
    if not s["poscar"]:
        raise ValueError("host POSCAR (starting structure) is required")
    s["poscar"] = take(s["poscar"])
    tpl = []
    for t in m["templates"]:
        label, p = split_template(t)
        if not re.fullmatch(r"[\w.+\-]+", label):
            raise ValueError(f"phase label {label!r}: letters, digits, _ . + - only")
        if label in [x.split("=", 1)[0] for x in tpl]:
            raise ValueError(f"phase label {label!r} used twice")
        tpl.append(f"{label}={take(p, 'refs/' + label)}")
    if not tpl:
        raise ValueError("at least one target phase template is needed")
    m["templates"] = tpl
    if fe["nep"]:
        fe["nep"] = take(fe["nep"])
    else:
        fe["nep"] = "nep.txt"
    if spec.get("freeze_scripts", True):
        dst = folder / "scripts"
        if dst.resolve() != app_scripts.resolve():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(app_scripts, dst, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            import os
            for p in [dst, *dst.rglob("*")]:     # frozen copy: only you may change it (umask 002 sources)
                os.chmod(p, p.stat().st_mode & ~0o022)
        cfg["paths"]["scripts"] = "scripts"
        script_dir = dst
    else:
        cfg["paths"]["scripts"] = str(app_scripts)
        script_dir = app_scripts
    one = lambda x: re.sub(r"[\x00-\x1f\x7f]+", " ", str(x)).encode("utf-8", "replace").decode("utf-8")
    spec = dict(spec, name=one(spec.get("name", "")), starting_phase=one(spec.get("starting_phase", "")))
    header = (f"# campaign.py config written by Fluorite Ordering Studio {APP_VERSION} on "
              f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n# run: {spec.get('name', '')}  starting phase: "
              f"{spec.get('starting_phase', '')}\n# Paths are relative to this file.  Reproduce: copy this folder "
              "(or use run_record.json) and run\n#   python3 scripts/campaign.py campaign.toml dry-run")
    (folder / "campaign.toml").write_text(to_toml(cfg, header), encoding="utf-8")
    rec = {
        "app": "fluorite-ordering-studio", "app_version": APP_VERSION,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"), "user": _user(), "host": socket.gethostname(),
        "name": spec.get("name", ""), "starting_phase": spec.get("starting_phase", ""),
        "notes": spec.get("notes", ""), "config": cfg, "inputs": inputs,
        "scripts": {p.name: md5(p) for p in sorted(script_dir.glob("*.py"))},
        "scripts_folder": "scripts" if spec.get("freeze_scripts", True) else str(script_dir),
        "campaign_toml_md5": md5(folder / "campaign.toml"),
        "versions": versions(),
        "versions_note": "packages of the Python that ran the app; SLURM jobs use [cluster] python + setup",
        "history": [],
    }
    old = folder / "run_record.json"
    with RECORD_LOCK:
        if old.is_file():
            try:
                rec["history"] = json.loads(old.read_text(encoding="utf-8")).get("history", [])
                rec["history"].append({"time": rec["created"], "event": "settings rewritten"})
            except ValueError:
                pass
        _write_record(old, rec)
    return {"folder": str(folder), "config": str(folder / "campaign.toml"), "record": rec}


def find_runs(base: Path, depth=2):
    """run folders (campaign.toml present) below base"""
    out = []
    base = Path(base).expanduser()
    if not base.is_dir():
        return out
    pats = ["campaign.toml"] + ["*/campaign.toml", "*/*/campaign.toml"][:depth]
    for pat in pats:
        for p in sorted(base.glob(pat)):
            out.append(read_run(p.parent))
    return out


def read_run(folder):
    folder = Path(folder).expanduser().resolve()
    info = {"folder": str(folder), "config_path": None, "record": None, "config": None}
    cp = folder / "campaign.toml"
    if not cp.is_file():
        tom = sorted(folder.glob("*.toml"))
        cp = tom[0] if tom else None
    if cp:
        info["config_path"] = str(cp)
        try:
            with cp.open("rb") as f:
                info["config"] = tomllib.load(f)
        except Exception as e:   # noqa: BLE001
            info["error"] = f"{cp.name}: {e}"
    rp = folder / "run_record.json"
    if rp.is_file():
        try:
            info["record"] = json.loads(rp.read_text(encoding="utf-8"))
        except ValueError as e:
            info["error"] = f"run_record.json: {e}"
    if info["config"]:
        wd = (folder / info["config"].get("paths", {}).get("workdir", "campaign")).resolve()
        info["workdir"] = str(wd)
        sp = info["config"].get("paths", {}).get("scripts")
        info["scripts"] = str((folder / sp).resolve()) if sp else None
    return info


def check_reproducible(folder):
    """compare the md5 of inputs and scripts with the record"""
    r = read_run(folder)
    rec = r.get("record")
    if not rec:
        return {"ok": False, "error": "no run_record.json in this folder"}
    folder = Path(r["folder"])
    diffs = []
    for rel, d in rec.get("inputs", {}).items():
        p = (folder / rel) if not Path(rel).is_absolute() else Path(rel)
        now = md5(p) if p.is_file() else None
        want = d.get("md5") if isinstance(d, dict) else None
        if now != want:
            diffs.append({"file": rel, "recorded": want, "now": now})
    sf = rec.get("scripts_folder") or "scripts"
    sdir = folder / sf if not Path(sf).is_absolute() else Path(sf)   # relative = frozen copy, travels with the folder
    if rec.get("campaign_toml_md5") and (folder / "campaign.toml").is_file():
        now = md5(folder / "campaign.toml")
        if now != rec["campaign_toml_md5"]:
            diffs.append({"file": "campaign.toml", "recorded": rec["campaign_toml_md5"], "now": now})
    for name, h in rec.get("scripts", {}).items():
        p = sdir / name
        now = md5(p) if p.is_file() else None
        if now != h:
            diffs.append({"file": f"scripts/{name}", "recorded": h, "now": now})
    v_now = versions()
    vdiff = {k: [v, v_now.get(k)] for k, v in rec.get("versions", {}).items() if v != v_now.get(k) and k != "platform"}
    return {"ok": True, "identical": not diffs, "differences": diffs, "version_differences": vdiff}


def _write_record(rp, rec):
    tmp = rp.with_name(rp.name + f".tmp{time.time_ns()}")
    tmp.write_text(json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
    import os
    os.replace(tmp, rp)


def add_history(folder, entry):
    rp = Path(folder) / "run_record.json"
    with RECORD_LOCK:
        if not rp.is_file():
            return
        try:
            rec = json.loads(rp.read_text(encoding="utf-8"))
        except ValueError:
            return
        rec.setdefault("history", []).append({"time": time.strftime("%Y-%m-%d %H:%M:%S"), **entry})
        _write_record(rp, rec)
