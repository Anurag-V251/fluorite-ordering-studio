#!/usr/bin/env python3
"""
xrd_average.py - ONE averaged powder XRD pattern of many relaxed structures (ASE + pymatgen)

Selects the structures of one composition (and phase) whose free energy lies in a chosen F window, takes
each one's relaxed geometry, computes its powder pattern with pymatgen's XRDCalculator, broadens every
pattern on a common 2theta grid and averages them (mean, standard deviation, min, max).  Optionally the
12-atom reference phases (refs/POSCAR_{c,t,oIII,m}) are computed with the same settings for comparison.

The relaxed POSCARs are only READ.  --replace / --replace-cations change the species in memory, for the
XRD only (e.g. every cation -> Zr shows the pattern of the distortion alone, without the Ce/Zr contrast).

Selection
    SOURCE  = a run folder (campaign.toml), a campaign workdir, a plain folder with Ce_xx/, or a CSV written
              by batch_free_energy.py / campaign.py merge (free_energy_<T>K_all.csv or Ce_xx/free_energy_<T>K.csv)
    --composition Ce_25.00 (folder name, or the percentage: 25)   --phase oIII [m ...]   --T 300
    --F-range LO HI   meV/atom above the lowest unflagged F of the composition (any phase)
                      = the axis of the Studio's XRD page and of the Results page's dF
    --F-abs LO HI     absolute F_total in eV/atom instead
    flagged rows (imaginary modes, unconverged relaxation) are left out unless --include-flagged
    --max N --seed S  random subset of N of the selected structures (reproducible)
    or give the structures directly: --structures POSCAR ... / --files LIST (one path per line)

Geometry of a selected row: the relaxed file that belongs to the row (campaign: <workdir>/relaxed/<file>,
rerun rounds: <workdir>/rerun_RRR/relaxed/<file> for rows of results/chunk_rerunRRR_*.csv; plain folders:
Ce_xx/relaxed/<file>).  Rows computed without relaxation (relaxed = False) use the input structure.  Rows
whose relaxed file is missing are skipped and reported (--allow-unrelaxed uses the input structure).

Pattern
    --wavelength CuKa (pymatgen names: CuKa CuKa1 MoKa CoKa ...) or a number in Angstrom
    --two-theta 10 90 --step 0.02     --fwhm 0.10 (deg, constant) or --caglioti U V W (FWHM^2 = U tan^2 + V tan + W)
    --profile gaussian | lorentzian | pseudo-voigt (--eta 0.5 = Lorentzian fraction); profiles have peak
              height = stick intensity
    --norm max       every structure scaled to max 100 before averaging (default)
           absolute  intensities / n_atoms^2 (independent of the supercell size), then averaged
    --symprec 0 (pymatgen refinement off: the cell is used as it is)  --debye-waller Zr=0.3 O=0.5 (B, A^2)

Output (--out-dir, default ./xrd_average):
    <stem>.csv              two_theta, mean, std, min, max, n, then one column per reference (max 100)
    <stem>_structures.csv   the selected structures: key, phase, F, dF, geometry used
    <stem>_ref_peaks.csv    reference sticks with hkl and d
    <stem>.png/.pdf         figure (matplotlib)
    <stem>.json             every setting, md5 of this script, versions, the command line
    --sticks-json FILE      instead: the unbroadened sticks of every structure (used by the Studio's XRD page)

Examples
    python xrd_average.py run_c/ --composition Ce_25.00 --phase oIII --F-range 0 2.5
    python xrd_average.py free_energy_300K_all.csv --composition 25 --phase m --F-range 1 4 --replace-cations Zr
    python xrd_average.py --structures relaxed/Ce_25.00/mapped/t/POSCAR_* --ref t c --fwhm 0.08
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import sys
import time
from pathlib import Path

__version__ = "1.0"
HERE = Path(__file__).resolve().parent
COMP_RE = re.compile(r"^(?P<el>[A-Z][a-z]?)_(?P<pct>\d+(?:\.\d+)?)$")
FE_RE = re.compile(r"^free_energy_(?P<T>[\d.]+)K\.csv$")
RERUN_CSV_RE = re.compile(r"^chunk_rerun(?P<R>\d{3})_")
NUM = r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
WINDOW = {"gaussian": 6.0, "lorentzian": 20.0, "pseudo-voigt": 20.0}     # profile cut-off (+-), in FWHM
CACHE_VERSION = 1


# =========================================================================== rows (standard library only:
# the Studio imports this part, so it must not need numpy / ASE / pymatgen)
def fnum(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def load_toml(p):
    try:
        import tomllib
    except ModuleNotFoundError:
        try:
            import tomli as tomllib            # type: ignore
        except ModuleNotFoundError:
            return None
    try:
        with open(p, "rb") as f:
            return tomllib.load(f)
    except Exception:                          # noqa: BLE001
        return None


def _csv_rows(p):
    try:
        with open(p, newline="") as f:
            yield from csv.DictReader(f)
    except (OSError, csv.Error, UnicodeDecodeError):
        return


def _fold(best, r, extra):
    """latest row per structure, an ok row beats a FAILED one (the rule of campaign.py merge)"""
    f = r.get("file")
    if not f:
        return
    if f not in best or r.get("status") == "ok" or best[f]["status"] != "ok":
        best[f] = {**extra, "file": f, "status": r.get("status", ""), "phase": r.get("phase") or "",
                   "F": fnum(r.get("F_total_eV_atom")), "T": f"{fnum(r.get('T_K')):g}",
                   "imag": fnum(r.get("n_imaginary_modes"), 0) > 0,
                   "unconv": str(r.get("relax_converged", "")).strip() == "False",
                   "relaxed": str(r.get("relaxed", "")).strip() != "False",
                   "n_atoms": r.get("n_atoms", ""), "id": r.get("id", "")}


def find_workdir(src: Path):
    """(workdir, layout) of a folder: run folder (campaign.toml), campaign workdir or plain folder"""
    if src.is_file() and src.suffix == ".toml":
        c = load_toml(src) or {}
        return (src.parent / c.get("paths", {}).get("workdir", "campaign")).resolve(), "campaign"
    if (src / "results").is_dir() or (src / "campaign_config.toml").is_file() or (src / "chunks").is_dir():
        return src, "campaign"
    for p in sorted(src.glob("*.toml")):
        c = load_toml(p)
        if c and "structure" in c and "mapping" in c:
            wd = (src / c.get("paths", {}).get("workdir", "campaign")).resolve()
            if wd.is_dir():
                return wd, "campaign"
    return src, "plain"


def load_rows(source, T=None, first_T=False):
    """every structure with a result: {key: row}, key = path relative to the workdir (Ce_xx/...).
    Returns (rows, info) with info = {workdir, layout, temperatures, source}"""
    src = Path(source).expanduser().resolve()
    by_T = {}
    if src.is_file() and src.suffix.lower() == ".csv":
        work, prefix = src.parent, ""
        if COMP_RE.match(src.parent.name):                       # Ce_xx/free_energy_<T>K.csv
            first = next(iter(_csv_rows(src)), {}).get("file", "")
            if first and not COMP_RE.match(first.split("/")[0]):
                work, prefix = src.parent.parent, src.parent.name + "/"
        layout = "campaign" if (work / "results").is_dir() else "plain"
        src_round = _rounds_of_results(work)                    # which round's row merge kept (geometry)
        for r in _csv_rows(src):
            if r.get("file"):
                r = {**r, "file": prefix + r["file"]}
                rd = src_round.get(r["file"], "unknown") if src_round else "unknown"
                _fold(by_T.setdefault(f"{fnum(r.get('T_K')):g}", {}), r, {"round": rd})
    elif src.is_dir() or src.suffix == ".toml":
        work, layout = find_workdir(src)
        if layout == "campaign" and (work / "results").is_dir():
            for p in sorted((work / "results").glob("chunk_*.csv")):
                m = RERUN_CSV_RE.match(p.name)
                for r in _csv_rows(p):
                    _fold(by_T.setdefault(f"{fnum(r.get('T_K')):g}", {}), r, {"round": int(m["R"]) if m else None})
        else:
            if layout == "campaign":
                layout = "plain"
            for d in sorted(work.iterdir()) if work.is_dir() else []:
                if not (d.is_dir() and COMP_RE.match(d.name)):
                    continue
                for f in sorted(os.listdir(d)):
                    m = FE_RE.match(f)
                    if m:
                        for r in _csv_rows(d / f):
                            if r.get("file"):
                                r = {**r, "file": f"{d.name}/{r['file']}"}
                                _fold(by_T.setdefault(f"{fnum(r.get('T_K')):g}", {}), r, {"round": "unknown"})
            if not by_T:                                           # only a merged CSV
                for p in sorted(work.glob("free_energy_*K_all.csv")):
                    for r in _csv_rows(p):
                        _fold(by_T.setdefault(f"{fnum(r.get('T_K')):g}", {}), r, {"round": "unknown"})
    else:
        raise ValueError(f"{src} is not a folder or a CSV")
    temps = sorted(by_T, key=fnum)
    if T is None:
        if len(temps) > 1 and not first_T:
            raise ValueError(f"several temperatures in {src}: {', '.join(temps)} K - choose one with --T")
        T = temps[0] if temps else None
    else:
        T = f"{float(T):g}"
        if T not in by_T:
            raise ValueError(f"no rows at T = {T} K in {src} (found: {', '.join(temps) or 'none'})")
    rows = by_T.get(T, {}) if T is not None else {}
    for k, r in rows.items():
        r["comp"] = k.split("/", 1)[0]
        if not r["phase"]:
            parts = k.split("/")
            r["phase"] = parts[2] if len(parts) > 3 and parts[1] == "mapped" else "parent"
    return rows, {"workdir": str(work), "layout": layout, "temperatures": temps, "T": T, "source": str(src)}


def _rounds_of_results(work: Path):
    """{file: rerun round or None} of the row that campaign.py merge keeps, from results/chunk_*.csv"""
    best = {}
    if not (work / "results").is_dir():
        return {}
    for p in sorted((work / "results").glob("chunk_*.csv")):
        m = RERUN_CSV_RE.match(p.name)
        for r in _csv_rows(p):
            f = r.get("file")
            if f and (f not in best or r.get("status") == "ok" or best[f][0] != "ok"):
                best[f] = (r.get("status"), int(m["R"]) if m else None)
    return {f: v[1] for f, v in best.items()}


def compositions(rows):
    """[(dir, element, pct, n_rows)] sorted by element and percentage"""
    cnt = {}
    for r in rows.values():
        cnt[r["comp"]] = cnt.get(r["comp"], 0) + 1
    out = []
    for c, n in cnt.items():
        m = COMP_RE.match(c)
        out.append((c, m["el"] if m else "", float(m["pct"]) if m else float("nan"), n))
    out.sort(key=lambda x: (x[1], x[2] if math.isfinite(x[2]) else 1e9, x[0]))
    return out


def match_comp(rows, comp):
    """composition folder name from 'Ce_25.00', '25' or '25.0'"""
    cs = compositions(rows)
    names = [c[0] for c in cs]
    if comp in names:
        return comp
    v = fnum(str(comp).split("_")[-1])
    hits = [c[0] for c in cs if math.isfinite(v) and abs(c[2] - v) < 1e-6]
    if len(hits) == 1:
        return hits[0]
    raise ValueError(f"composition {comp!r} not found (have: {', '.join(names[:40])}{' ...' if len(names) > 40 else ''})")


def reference_F(rows, comp):
    """lowest unflagged F of the composition, any phase (lowest ok F if every row is flagged)"""
    ok = [r for r in rows.values() if r["comp"] == comp and r["status"] == "ok" and math.isfinite(r["F"])]
    clean = [r["F"] for r in ok if not (r["imag"] or r["unconv"])]
    return min(clean) if clean else (min(r["F"] for r in ok) if ok else float("nan"))


def select(rows, comp, phases=None, F_range=None, F_abs=None, include_flagged=False, max_n=None, seed=0):
    """the rows in the window, sorted by F; (selected, counts).  Bounds are inclusive."""
    Fref = reference_F(rows, comp)
    lo, hi = -math.inf, math.inf
    if F_range is not None:
        lo, hi = Fref + F_range[0] / 1000.0, Fref + F_range[1] / 1000.0
    if F_abs is not None:
        lo, hi = F_abs
    eps = 1e-12
    cnt = {"in_window": 0, "flagged_excluded": 0, "not_ok": 0, "subsampled_from": None}
    sel = []
    for k, r in rows.items():
        if r["comp"] != comp or (phases and r["phase"] not in phases):
            continue
        if r["status"] != "ok" or not math.isfinite(r["F"]):
            cnt["not_ok"] += 1
            continue
        if not (lo - eps <= r["F"] <= hi + eps):
            continue
        cnt["in_window"] += 1
        if (r["imag"] or r["unconv"]) and not include_flagged:
            cnt["flagged_excluded"] += 1
            continue
        sel.append({**r, "key": k, "dF_meV": 1000.0 * (r["F"] - Fref)})
    sel.sort(key=lambda r: (r["F"], r["key"]))
    if max_n and len(sel) > max_n:
        cnt["subsampled_from"] = len(sel)
        keep = set(random.Random(seed).sample(range(len(sel)), max_n))
        sel = [r for i, r in enumerate(sel) if i in keep]
    cnt["selected"] = len(sel)
    cnt["F_ref_eV_atom"] = Fref
    return sel, cnt


def geometry(work, row, allow_unrelaxed=False):
    """(path, how) of the structure to use for a row; path None if no relaxed file is there"""
    work = Path(work)
    k = row["key"]
    inp = work / k
    if not row.get("relaxed", True):
        return (inp if inp.is_file() else None), "input (computed without relaxation)"
    rd = row.get("round", "unknown")       # int: rerun round of the row; None: campaign row; "unknown": newest
    if isinstance(rd, int):
        cands = [(work / f"rerun_{rd:03d}" / "relaxed" / k, f"rerun round {rd}")]
    else:
        cands = [(work / "relaxed" / k, "relaxed")]
        if rd == "unknown":
            rounds = sorted(int(p.name.split("_")[1]) for p in work.glob("rerun_[0-9][0-9][0-9]") if p.is_dir())
            cands = [(work / f"rerun_{R:03d}" / "relaxed" / k, f"rerun round {R}") for R in reversed(rounds)] + cands
    comp, _, rest = k.partition("/")
    cands.append((work / comp / "relaxed" / rest, "relaxed"))                   # plain layout
    for p, how in cands:
        if p.is_file():
            return p, how
    if allow_unrelaxed and inp.is_file():
        return inp, "INPUT (no relaxed file found)"
    return None, "no relaxed file"


# =========================================================================== XRD (numpy, ASE, pymatgen)
def parse_replace(specs, cations_to=None, anions=("O",)):
    """--replace Ce:Zr ... and --replace-cations X -> (mapping dict or None, description)"""
    mp = {}
    for s in specs or []:
        for part in re.split(r"[,\s]+", s.strip()):
            if not part:
                continue
            m = re.fullmatch(r"([A-Z][a-z]?)\s*[:=>]+\s*([A-Z][a-z]?)", part)
            if not m:
                raise ValueError(f"--replace {part!r}: write OLD:NEW, e.g. Ce:Zr")
            mp[m[1]] = m[2]
    rule = {"map": mp, "cations_to": cations_to or None, "anions": sorted(set(anions or []))}
    if not mp and not cations_to:
        return None
    return rule


def replace_species(structure, rule):
    """copy of a pymatgen Structure with the species replaced (the file is not touched)"""
    if not rule:
        return structure
    from pymatgen.core import Element
    s = structure.copy()
    mp = {}
    for el in s.composition.elements:
        sym = el.symbol
        new = rule["map"].get(sym)
        if new is None and rule.get("cations_to") and sym not in rule.get("anions", []):
            new = rule["cations_to"]
        if new and new != sym:
            mp[el] = Element(new)
    if mp:
        s.replace_species(mp)
    return s


def wavelength_A(w):
    from pymatgen.analysis.diffraction.xrd import WAVELENGTHS
    if isinstance(w, (int, float)):
        return float(w), f"{float(w):g} A"
    if str(w) in WAVELENGTHS:
        return float(WAVELENGTHS[str(w)]), str(w)
    v = fnum(w)
    if math.isfinite(v) and v > 0:
        return v, f"{v:g} A"
    raise ValueError(f"wavelength {w!r}: a pymatgen name ({', '.join(WAVELENGTHS)}) or a number in A")


def file_md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stick_params(a):
    lam, lam_name = wavelength_A(a.wavelength)
    lo, hi = a.two_theta
    pad = getattr(a, "pad", 0.0) or 0.0
    return {"wavelength_A": lam, "wavelength": lam_name, "range": [max(0.0, lo - pad), min(180.0, hi + pad)],
            "symprec": a.symprec, "debye_waller": dict(sorted(dw_dict(a.debye_waller).items())),
            "replace": parse_replace(a.replace, a.replace_cations, a.anions), "min_rel": a.min_rel}


def dw_dict(specs):
    out = {}
    for s in specs or []:
        m = re.fullmatch(r"([A-Z][a-z]?)=(" + NUM + ")", s.strip())
        if not m:
            raise ValueError(f"--debye-waller {s!r}: write El=B, e.g. Zr=0.3")
        out[m[1]] = float(m[2])
    return out


def sticks_of(path, P, scale_volume_to=None, replace=True):
    """sticks of one structure: dict(x, y (per n_atoms^2), hkl, d, n_atoms, formula, vol_atom)"""
    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    atoms = read(str(path), format="vasp")
    s = AseAtomsAdaptor.get_structure(atoms)
    if scale_volume_to:
        s = s.copy()
        s.scale_lattice(scale_volume_to * len(s))
    if replace:
        s = replace_species(s, P["replace"])
    calc = XRDCalculator(wavelength=P["wavelength_A"], symprec=P["symprec"],
                         debye_waller_factors=P["debye_waller"] or None)
    pat = calc.get_pattern(s, scaled=False, two_theta_range=tuple(P["range"]))
    n = len(s)
    x = [float(v) for v in pat.x]
    y = [float(v) / n ** 2 for v in pat.y]
    keep = [i for i, v in enumerate(y) if v >= P["min_rel"] * (max(y) if y else 0)]
    return {"x": [round(x[i], 5) for i in keep], "y": [float(f"{y[i]:.5g}") for i in keep],
            "hkl": [_hkl_text(pat.hkls[i]) for i in keep], "d": [round(float(pat.d_hkls[i]), 5) for i in keep],
            "n_atoms": n, "formula": s.composition.formula.replace(" ", ""),
            "formula_file": atoms.get_chemical_formula(mode="metal"), "vol_atom": s.volume / n}


def _hkl_text(fams):
    """'(1 1 1)x8' of the first family of a pymatgen peak"""
    if not fams:
        return ""
    f = fams[0]
    return "(" + " ".join(str(int(v)) for v in f["hkl"]) + ")" + (f"x{f['multiplicity']}" if f.get("multiplicity") else "")


def _cache_key(md5, P):
    return hashlib.sha1(json.dumps([CACHE_VERSION, md5, P], sort_keys=True).encode()).hexdigest()


def _work(args):
    path, P, cache_dir = args
    try:
        md5 = file_md5(path)
        key = _cache_key(md5, P)
        cp = Path(cache_dir) / key[:2] / f"{key}.json" if cache_dir else None
        if cp and cp.is_file():
            try:
                d = json.loads(cp.read_text())
                d["cached"] = True
                return path, d, None
            except (OSError, ValueError):
                pass
        d = sticks_of(path, P)
        d.pop("hkl", None)                                 # hkl and d are kept for the references only
        d.pop("d", None)
        d["md5"] = md5
        if cp:
            cp.parent.mkdir(parents=True, exist_ok=True)
            tmp = cp.with_suffix(f".tmp{os.getpid()}")
            tmp.write_text(json.dumps(d))
            os.replace(tmp, cp)
        d["cached"] = False
        return path, d, None
    except Exception as e:                                 # noqa: BLE001 - reported per structure
        return path, None, f"{type(e).__name__}: {e}"


def compute_sticks(paths, P, jobs=1, cache_dir=None, progress=True):
    out, errs = {}, {}
    tasks = [(str(p), P, cache_dir) for p in paths]
    t0 = time.time()
    done = 0
    every = max(1, len(tasks) // 20)

    def note(path, d, e):
        nonlocal done
        done += 1
        if e:
            errs[path] = e
        else:
            out[path] = d
        if progress and (done % every == 0 or done == len(tasks)):
            print(f"XRD {done}/{len(tasks)} ({time.time() - t0:.1f} s)", flush=True)

    if jobs > 1 and len(tasks) > 4:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for path, d, e in ex.map(_work, tasks, chunksize=max(1, len(tasks) // (jobs * 8))):
                note(path, d, e)
    else:
        for t in tasks:
            note(*_work(t))
    return out, errs


def fwhm_at(x_deg, a):
    """FWHM (deg 2theta) at the given 2theta values"""
    import numpy as np
    x = np.asarray(x_deg, float)
    if a.caglioti:
        U, V, W = a.caglioti
        t = np.tan(np.radians(x / 2))
        return np.sqrt(np.maximum(U * t * t + V * t + W, 1e-8))
    return np.full_like(x, float(a.fwhm))


def broaden(x, y, grid, a):
    """sum of profiles (height = stick intensity) on the grid; the same algorithm as the Studio page"""
    import numpy as np
    out = np.zeros_like(grid)
    if not len(x):
        return out
    g0, step = grid[0], grid[1] - grid[0] if len(grid) > 1 else 1.0
    fw = fwhm_at(x, a)
    win = WINDOW[a.profile]
    for x0, I, w in zip(x, y, fw):
        i0 = max(0, int(math.floor((x0 - win * w - g0) / step)))
        i1 = min(len(grid), int(math.ceil((x0 + win * w - g0) / step)) + 1)
        if i1 <= i0:
            continue
        u = (grid[i0:i1] - x0) / w
        if a.profile == "gaussian":
            p = np.exp(-4 * math.log(2) * u * u)
        elif a.profile == "lorentzian":
            p = 1.0 / (1.0 + 4 * u * u)
        else:
            p = a.eta / (1.0 + 4 * u * u) + (1 - a.eta) * np.exp(-4 * math.log(2) * u * u)
        out[i0:i1] += I * p
    return out


def make_grid(a):
    import numpy as np
    lo, hi = a.two_theta
    n = int(round((hi - lo) / a.step)) + 1
    return lo + a.step * np.arange(n)


def average(patterns, grid, a):
    """mean, std, min, max over structures of the broadened patterns (--norm)"""
    import numpy as np
    Y = []
    for d in patterns:
        y = broaden(np.asarray(d["x"]), np.asarray(d["y"]), grid, a)
        if a.norm == "max":
            m = y.max()
            y = 100.0 * y / m if m > 0 else y
        Y.append(y)
    Y = np.array(Y) if Y else np.zeros((0, len(grid)))
    if not len(Y):
        z = np.zeros_like(grid)
        return z, z, z, z
    return Y.mean(0), Y.std(0), Y.min(0), Y.max(0)


# =========================================================================== command line
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].split(" - ", 1)[-1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("source", nargs="?", help="run folder, campaign workdir, plain folder or result CSV")
    g = ap.add_argument_group("selection")
    g.add_argument("--composition", help="composition folder (Ce_25.00) or its percentage (25)")
    g.add_argument("--phase", nargs="+", help="phase label(s), e.g. oIII; default: every phase")
    g.add_argument("--F-range", dest="F_range", type=float, nargs=2, metavar=("LO", "HI"),
                   help="meV/atom above the lowest unflagged F of the composition (inclusive)")
    g.add_argument("--F-abs", dest="F_abs", type=float, nargs=2, metavar=("LO", "HI"),
                   help="absolute F_total window in eV/atom (instead of --F-range)")
    g.add_argument("--T", help="temperature of the rows (needed when the source has several)")
    g.add_argument("--include-flagged", action="store_true",
                   help="also rows with imaginary modes / unconverged relaxation (left out by default)")
    g.add_argument("--allow-unrelaxed", action="store_true",
                   help="use the input structure when the relaxed file is missing (default: skip the row)")
    g.add_argument("--max", dest="max_n", type=int, help="random subset of this many selected structures")
    g.add_argument("--seed", type=int, default=0, help="seed of the random subset (default 0)")
    g.add_argument("--structures", nargs="+", metavar="POSCAR", help="structure files instead of a selection")
    g.add_argument("--files", help="text file with one structure path per line (relative to its folder)")
    g.add_argument("--list", action="store_true", help="only print the selected structures and stop")
    g = ap.add_argument_group("XRD (pymatgen XRDCalculator)")
    g.add_argument("--wavelength", default="CuKa", help="pymatgen name (CuKa, MoKa, CuKa1 ...) or Angstrom")
    g.add_argument("--two-theta", dest="two_theta", type=float, nargs=2, default=[10.0, 90.0], metavar=("MIN", "MAX"),
                   help="2theta range in degrees (default 10 90)")
    g.add_argument("--symprec", type=float, default=0.0,
                   help="pymatgen symmetry refinement before the pattern (default 0 = off: the cell as it is)")
    g.add_argument("--debye-waller", nargs="+", metavar="El=B", help="Debye-Waller B factors in A^2, e.g. Zr=0.3 O=0.5")
    g.add_argument("--min-rel", type=float, default=1e-4,
                   help="drop sticks weaker than this fraction of the structure's strongest (default 1e-4)")
    g = ap.add_argument_group("species, for the XRD only (the files are never changed)")
    g.add_argument("--replace", nargs="+", metavar="OLD:NEW", help="e.g. Ce:Zr (several allowed)")
    g.add_argument("--replace-cations", metavar="EL", help="every species that is not an anion -> EL (e.g. Zr)")
    g.add_argument("--anions", nargs="+", default=["O"], help="species --replace-cations leaves alone (default O)")
    g = ap.add_argument_group("broadening and averaging")
    g.add_argument("--step", type=float, default=0.02, help="2theta grid step in degrees (default 0.02)")
    g.add_argument("--fwhm", type=float, default=0.10, help="peak FWHM in degrees 2theta (default 0.10)")
    g.add_argument("--caglioti", type=float, nargs=3, metavar=("U", "V", "W"),
                   help="2theta-dependent FWHM^2 = U tan^2(theta) + V tan(theta) + W (deg^2) instead of --fwhm")
    g.add_argument("--profile", choices=["gaussian", "lorentzian", "pseudo-voigt"], default="gaussian")
    g.add_argument("--eta", type=float, default=0.5, help="Lorentzian fraction of the pseudo-Voigt (default 0.5)")
    g.add_argument("--norm", choices=["max", "absolute"], default="max",
                   help="max: each pattern scaled to max 100 before averaging; absolute: per n_atoms^2")
    g = ap.add_argument_group("references")
    g.add_argument("--ref", nargs="+", metavar="NAME[=POSCAR]",
                   help="reference phases: names in --refs (c t oIII m) or NAME=path")
    g.add_argument("--refs", default=str(HERE.parent / "refs"), help="folder with POSCAR_<name> (default ../refs)")
    g.add_argument("--ref-replace", action=argparse.BooleanOptionalAction, default=True,
                   help="apply --replace / --replace-cations to the references too (default yes)")
    g.add_argument("--ref-scale-volume", action="store_true",
                   help="scale each reference cell isotropically to the mean volume/atom of the selection")
    g = ap.add_argument_group("output")
    g.add_argument("--out-dir", default="xrd_average", help="folder for the results (default ./xrd_average)")
    g.add_argument("--stem", help="file name stem (default from the selection)")
    g.add_argument("--formats", nargs="+", default=["png"], help="figure formats (default png; e.g. png pdf svg)")
    g.add_argument("--no-plot", action="store_true", help="no figure (CSV only)")
    g.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1), help="parallel processes (default <= 4)")
    g.add_argument("--cache-dir", help="keep the sticks of every structure here (re-runs with the same settings are instant)")
    g.add_argument("--sticks-json", help="write the unbroadened sticks of every structure to this JSON and stop")
    g.add_argument("--pad", type=float, default=0.0, help=argparse.SUPPRESS)
    return ap


def check_args(ap, a):
    if a.two_theta[0] >= a.two_theta[1] or a.two_theta[0] < 0 or a.two_theta[1] > 180:
        ap.error("--two-theta MIN MAX needs 0 <= MIN < MAX <= 180")
    if not a.step > 0 or (a.two_theta[1] - a.two_theta[0]) / a.step > 2_000_000:
        ap.error("--step must be > 0 (and give at most 2 million points)")
    if not a.caglioti and not a.fwhm > 0:
        ap.error("--fwhm must be > 0")
    if not 0 <= a.eta <= 1:
        ap.error("--eta must be in [0, 1]")
    if a.F_range and a.F_abs:
        ap.error("give --F-range or --F-abs, not both")
    for w in (a.F_range, a.F_abs):
        if w and w[0] > w[1]:
            ap.error("F window: LO must not be above HI")
    if a.max_n is not None and a.max_n < 1:
        ap.error("--max must be >= 1")
    try:
        parse_replace(a.replace, a.replace_cations, a.anions)
        dw_dict(a.debye_waller)
    except ValueError as e:
        ap.error(str(e))
    for el in [a.replace_cations] if a.replace_cations else []:
        if not re.fullmatch(r"[A-Z][a-z]?", el):
            ap.error(f"--replace-cations {el!r} is not an element symbol")


def gather(ap, a):
    """[(path, info dict)], selection report"""
    if a.structures or a.files:
        paths = list(a.structures or [])
        if a.files:
            base = Path(a.files).resolve().parent
            for line in Path(a.files).read_text().split("\n"):
                line = line.strip()
                if line and not line.startswith("#"):
                    p = Path(line).expanduser()
                    paths.append(str(p if p.is_absolute() else base / p))
        out = [(str(Path(p).resolve()), {"key": str(p), "phase": "", "F": float("nan"), "dF_meV": float("nan"),
                                         "how": "given"}) for p in paths]
        missing = [p for p, _ in out if not Path(p).is_file()]
        if missing:
            ap.error(f"{len(missing)} structure file(s) not found, e.g. {missing[0]}")
        return out, {"mode": "files", "selected": len(out)}
    if not a.source:
        ap.error("give a SOURCE (run folder, workdir or CSV) or --structures / --files")
    try:
        rows, info = load_rows(a.source, a.T)
    except ValueError as e:
        ap.error(str(e))
    if not rows:
        ap.error(f"no result rows in {a.source}")
    comps = compositions(rows)
    if not a.composition:
        if len(comps) == 1:
            a.composition = comps[0][0]
        else:
            ap.error("choose --composition: " + ", ".join(c[0] for c in comps[:60]))
    try:
        comp = match_comp(rows, a.composition)
    except ValueError as e:
        ap.error(str(e))
    sel, cnt = select(rows, comp, a.phase, a.F_range, a.F_abs, a.include_flagged, a.max_n, a.seed)
    out, skipped = [], []
    for r in sel:
        p, how = geometry(info["workdir"], r, a.allow_unrelaxed)
        if p is None:
            skipped.append(r["key"])
            continue
        out.append((str(p), {**r, "how": how}))
    rep = {"mode": "selection", **info, "composition": comp, "phases": a.phase or "all", **cnt,
           "no_geometry": len(skipped), "no_geometry_examples": skipped[:5], "computed": len(out)}
    return out, rep


def default_stem(a, rep):
    if rep.get("mode") != "selection":
        return "xrd_structures"
    ph = "-".join(a.phase) if a.phase else "all"
    if a.F_range:
        w = f"_{a.F_range[0]:g}-{a.F_range[1]:g}meV"
    elif a.F_abs:
        w = f"_{a.F_abs[0]:g}-{a.F_abs[1]:g}eV"
    else:
        w = ""
    return f"xrd_{rep['composition']}_{ph}{w}"


def ref_list(a):
    out = []
    for r in a.ref or []:
        name, _, p = r.partition("=")
        path = Path(p).expanduser() if p else Path(a.refs) / f"POSCAR_{name}"
        if not path.is_file():
            raise SystemExit(f"reference {name}: {path} not found")
        out.append((name, path))
    return out


def ref_sticks(a, P, mean_vol):
    out = []
    for name, path in ref_list(a):
        P2 = P if a.ref_replace else {**P, "replace": None}
        d = sticks_of(path, P2, scale_volume_to=mean_vol if a.ref_scale_volume else None)
        d.update(name=name, path=str(path))
        out.append(d)
    return out


def versions():
    v = {"python": sys.version.split()[0]}
    for m in ("numpy", "ase", "pymatgen.core"):
        try:
            v[m.split(".")[0]] = __import__(m, fromlist=["__version__"]).__version__
        except Exception:                       # noqa: BLE001
            v[m.split(".")[0]] = None
    return v


def main(argv=None):
    ap = build_parser()
    a = ap.parse_args(argv)
    check_args(ap, a)
    items, rep = gather(ap, a)
    if rep.get("mode") == "selection":
        print(f"{rep['composition']} | phases {rep['phases']} | T {rep['T']} K | F_ref {rep['F_ref_eV_atom']:.6f} eV/atom | "
              f"{rep['in_window']} in the window, {rep['flagged_excluded']} flagged left out, "
              f"{rep['no_geometry']} without relaxed file, {len(items)} to compute"
              + (f" (random {a.max_n} of {rep['subsampled_from']}, seed {a.seed})" if rep.get("subsampled_from") else ""),
              flush=True)
        if rep["no_geometry"]:
            print(f"  no relaxed file for e.g. {', '.join(rep['no_geometry_examples'])} (--allow-unrelaxed uses the input)")
    if a.list:
        for p, r in items:
            print(f"{r['key']}\t{r['phase']}\t{r['dF_meV']:.4f}\t{r['how']}\t{p}")
        return 0
    if not items and not a.ref:
        print("nothing selected", file=sys.stderr)
        return 2
    try:
        P = stick_params(a)
    except ValueError as e:
        ap.error(str(e))
    sticks, errs = compute_sticks([p for p, _ in items], P, a.jobs, a.cache_dir)
    for p, e in list(errs.items())[:10]:
        print(f"  FAILED {p}: {e}", file=sys.stderr)
    ok = [(p, r) for p, r in items if p in sticks]
    vols = [sticks[p]["vol_atom"] for p, _ in ok]
    mean_vol = sum(vols) / len(vols) if vols else None
    try:
        refs = ref_sticks(a, P, mean_vol)
    except Exception as e:                      # noqa: BLE001
        raise SystemExit(f"reference failed: {type(e).__name__}: {e}")
    record = {"script": "xrd_average.py", "version": __version__, "script_md5": file_md5(__file__),
              "command": [sys.executable, *sys.argv] if argv is None else ["xrd_average.py", *argv],
              "stick_params": P, "selection": rep, "versions": versions(), "failed": errs,
              "created": time.strftime("%Y-%m-%d %H:%M:%S"),
              "note": "relaxed POSCARs are only read; species replacement is applied in memory for the XRD only"}
    if a.sticks_json:
        out = {**record, "structures": [{"path": p, "key": r["key"], "phase": r["phase"], "F": r["F"],
                                         "dF_meV": r["dF_meV"], "how": r["how"], **compact(sticks[p])} for p, r in ok],
               "refs": refs, "encoding": "x = xi / 1e4 deg, y = yi * ymax / 1e6 (per n_atoms^2)"}
        Path(a.sticks_json).parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(f"{a.sticks_json}.tmp{os.getpid()}")
        tmp.write_text(json.dumps(out, allow_nan=False, default=_nan_none))
        os.replace(tmp, a.sticks_json)
        print(f"wrote {a.sticks_json} ({len(ok)} structures, {len(refs)} references, "
              f"{sum(1 for p, _ in ok if sticks[p].get('cached'))} from the cache, {len(errs)} failed)")
        return 0 if ok or refs else 2
    import numpy as np
    grid = make_grid(a)
    mean, std, mn, mx = average([sticks[p] for p, _ in ok], grid, a)
    ref_y = {}
    for d in refs:
        y = broaden(np.asarray(d["x"]), np.asarray(d["y"]), grid, a)
        ref_y[d["name"]] = 100 * y / y.max() if y.max() > 0 else y
    od = Path(a.out_dir)
    od.mkdir(parents=True, exist_ok=True)
    stem = a.stem or default_stem(a, rep)
    with open(od / f"{stem}.csv", "w", newline="") as f:
        w = csv.writer(f)
        unit = "rel_max100" if a.norm == "max" else "per_atom2"
        w.writerow(["two_theta_deg", f"mean_{unit}", "std", "min", "max", "n"] + [f"ref_{k}_max100" for k in ref_y])
        for i, x in enumerate(grid):
            w.writerow([f"{x:.5f}", f"{mean[i]:.6g}", f"{std[i]:.6g}", f"{mn[i]:.6g}", f"{mx[i]:.6g}", len(ok)] +
                       [f"{ref_y[k][i]:.6g}" for k in ref_y])
    with open(od / f"{stem}_structures.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["key", "phase", "F_total_eV_atom", "dF_meV_atom", "geometry", "path", "n_atoms", "formula_used", "md5"])
        for p, r in ok:
            s = sticks[p]
            w.writerow([r["key"], r["phase"], _fmt(r["F"], 8), _fmt(r["dF_meV"], 4), r["how"], p, s["n_atoms"],
                        s["formula"], s.get("md5", "")])
    if refs:
        with open(od / f"{stem}_ref_peaks.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ref", "two_theta_deg", "d_A", "I_rel_max100", "hkl", "formula_used"])
            for d in refs:
                m = max(d["y"]) if d["y"] else 1
                for x, y, hk, dd in zip(d["x"], d["y"], d["hkl"], d["d"]):
                    if 100 * y / m >= 0.5:
                        w.writerow([d["name"], f"{x:.4f}", f"{dd:.5f}", f"{100 * y / m:.3f}", hk, d["formula"]])
    record.update(grid={"two_theta": list(a.two_theta), "step": a.step}, broadening={
        "fwhm": a.fwhm, "caglioti": a.caglioti, "profile": a.profile, "eta": a.eta}, norm=a.norm,
        refs=[{k: d[k] for k in ("name", "path", "formula", "formula_file", "vol_atom", "n_atoms")} for d in refs],
        n_structures=len(ok))
    (od / f"{stem}.json").write_text(json.dumps(record, indent=1, default=_nan_none))
    if not a.no_plot:
        plot(od / stem, grid, mean, std, mn, mx, ref_y, refs, a, rep, len(ok), P)
    print(f"wrote {od / stem}.csv (+ _structures.csv, .json{', _ref_peaks.csv' if refs else ''}"
          f"{', ' + ', '.join(a.formats) if not a.no_plot else ''}) from {len(ok)} structures"
          + (f", {len(errs)} failed" if errs else ""))
    return 0


def compact(d):
    """sticks as integers (small JSON): xi = 2theta * 1e4, yi = y / ymax * 1e6"""
    ym = max(d["y"]) if d["y"] else 0.0
    out = {k: v for k, v in d.items() if k not in ("x", "y")}
    out["xi"] = [int(round(x * 1e4)) for x in d["x"]]
    out["yi"] = [int(round(y / ym * 1e6)) if ym > 0 else 0 for y in d["y"]]
    out["ymax"] = ym
    return out


def _nan_none(o):
    return None


def _fmt(v, d):
    return f"{v:.{d}f}" if isinstance(v, float) and math.isfinite(v) else ""


def plot(base, grid, mean, std, mn, mx, ref_y, refs, a, rep, n, P):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed: no figure")
        return
    colors = {"c": "#52514e", "parent": "#52514e", "t": "#2a78d6", "m": "#eb6834", "oIII": "#1baf7a"}
    extra = ["#6250d6", "#e34948", "#eda100", "#e87ba4"]
    scale = 100.0 / mean.max() if a.norm == "absolute" and mean.max() > 0 else 1.0
    nref = len(ref_y)
    fig, ax = plt.subplots(figsize=(7.0, 3.2 + 0.55 * nref), dpi=150)
    sel_col = colors.get((a.phase or [""])[0], "#1f1f1f") if a.phase and len(a.phase) == 1 else "#1f1f1f"
    ax.fill_between(grid, scale * (mean - std), scale * (mean + std), color=sel_col, alpha=0.18, lw=0, label="±1 std")
    ax.plot(grid, scale * mean, color=sel_col, lw=1.0, label=f"mean of {n}")
    off = 0.0
    for i, (k, y) in enumerate(ref_y.items()):
        off -= 60.0
        c = colors.get(k, extra[i % len(extra)])
        ax.plot(grid, off + 0.5 * y, color=c, lw=0.8)
        ax.text(grid[-1], off + 4, f"{k} ({refs[i]['formula']})", ha="right", va="bottom", fontsize=7, color=c)
    ax.set_xlim(grid[0], grid[-1])
    ax.set_ylim(off - 5 if nref else -3, 110)
    ax.set_yticks([0, 50, 100])
    ax.set_xlabel(r"$2\theta$ (deg), " + P["wavelength"] + f" ({P['wavelength_A']:.5f} Å)")
    ax.set_ylabel("intensity (mean = 100)" if a.norm == "absolute" else "intensity (each max = 100)")
    t = (f"{rep.get('composition', '')} {('/'.join(a.phase) if a.phase else 'all phases')}"
         + (f", F {a.F_range[0]:g}..{a.F_range[1]:g} meV/atom" if a.F_range else "")
         + (f", F {a.F_abs[0]:g}..{a.F_abs[1]:g} eV/atom" if a.F_abs else "") + f", n = {n}")
    rp = P.get("replace")
    if rp:
        t += " | XRD with " + ", ".join([f"{k}->{v}" for k, v in rp["map"].items()] +
                                         ([f"cations->{rp['cations_to']}"] if rp.get("cations_to") else []))
    ax.set_title(t, fontsize=8)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    for fm in a.formats:
        fig.savefig(f"{base}.{fm}")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
