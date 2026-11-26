"""
xrd.py - server side of the XRD page (standard library only).

The rows (which structures exist, their F and flags) and the selection rule come from the app's own
scripts/xrd_average.py (its first part needs no numpy), so the page, the job and the command line all select
the same structures.  The pymatgen work runs as a job: xrd_average.py --sticks-json writes the unbroadened
sticks of every selected structure (cached per file md5 + settings), the browser broadens and averages them,
so the FWHM, step, profile and normalisation change the plot at once.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
SCRIPT = APP / "scripts" / "xrd_average.py"
_LIB = {}
EL = re.compile(r"^[A-Z][a-z]?$")
NAME = re.compile(r"^[A-Za-z0-9_.+-]{1,40}$")


def lib():
    """scripts/xrd_average.py as a module (reloaded when the file changes)"""
    mt = SCRIPT.stat().st_mtime
    if _LIB.get("mtime") != mt:
        spec = importlib.util.spec_from_file_location("xrd_average_lib", SCRIPT)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _LIB.update(mtime=mt, mod=m)
    return _LIB["mod"]


def rows(source, T=None, comp=None):
    L = lib()
    rws, info = L.load_rows(source, T or None, first_T=True)
    comps = [{"dir": c, "el": el, "pct": (pct if math.isfinite(pct) else None), "n": n}
             for c, el, pct, n in L.compositions(rws)]
    if not comps:
        return {"ok": False, "error": f"no result rows in {source}", **info}
    names = [c["dir"] for c in comps]
    if comp not in names:
        comp = names[0]
    Fref = L.reference_F(rws, comp)
    out = []
    for k, r in rws.items():
        if r["comp"] != comp:
            continue
        out.append({"k": k, "p": r["phase"], "st": r["status"],
                    "F": r["F"] if math.isfinite(r["F"]) else None,
                    "im": bool(r["imag"]), "un": bool(r["unconv"])})
    out.sort(key=lambda x: (x["F"] is None, x["F"] or 0.0, x["k"]))
    return {"ok": True, **info, "compositions": comps, "comp": comp,
            "F_ref": Fref if math.isfinite(Fref) else None, "rows": out}


def _num(v, what, lo=-math.inf, hi=math.inf):
    try:
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{what}: {v!r} is not a number")
    if not (math.isfinite(x) and lo <= x <= hi):
        raise ValueError(f"{what}: {v!r} is outside [{lo}, {hi}]")
    return x


def build_args(b, mode, state_dir: Path):
    """argument list of xrd_average.py from the page's settings.  mode: sticks (page) | export (files)"""
    src = str(Path(str(b.get("source") or "")).expanduser())
    if not src or not Path(src).is_absolute() or not Path(src).exists():
        raise ValueError("choose the results folder or CSV first")
    a = [src]
    if b.get("T") not in (None, ""):
        a += ["--T", f"{_num(b['T'], 'T', 0):g}"]
    comp = str(b.get("comp") or "")
    if not re.fullmatch(r"[A-Za-z0-9_.+-]{1,60}", comp):
        raise ValueError("choose a composition")
    a += ["--composition", comp]
    ph = [str(p) for p in b.get("phases") or []]
    if any(not NAME.match(p) for p in ph):
        raise ValueError("bad phase name")
    if not ph:
        raise ValueError("choose at least one phase")
    a += ["--phase", *ph]
    if b.get("lo") is not None and b.get("hi") is not None:
        lo, hi = _num(b["lo"], "F from"), _num(b["hi"], "F to")
        if lo > hi:
            lo, hi = hi, lo
        a += ["--F-range", f"{lo:.6g}", f"{hi:.6g}"]
    if b.get("include_flagged"):
        a.append("--include-flagged")
    if b.get("allow_unrelaxed"):
        a.append("--allow-unrelaxed")
    if b.get("max") not in (None, "", 0):
        a += ["--max", str(int(_num(b["max"], "max structures", 1))), "--seed", str(int(_num(b.get("seed") or 0, "seed", 0)))]
    p = b.get("params") or {}
    w = str(p.get("wavelength") or "CuKa")
    if not re.fullmatch(r"[A-Za-z]{2}K[ab]\d?|\d*\.?\d+", w):
        raise ValueError(f"wavelength {w!r}")
    a += ["--wavelength", w]
    t0, t1 = _num(p.get("tt_lo", 10), "2theta min", 0, 180), _num(p.get("tt_hi", 90), "2theta max", 0, 180)
    if t0 >= t1:
        raise ValueError("2theta min must be below 2theta max")
    a += ["--two-theta", f"{t0:g}", f"{t1:g}"]
    a += ["--symprec", f"{_num(p.get('symprec', 0), 'symprec', 0, 1):g}",
          "--min-rel", f"{_num(p.get('min_rel', 1e-4), 'min-rel', 0, 0.5):g}"]
    dw = [x for x in re.split(r"[\s,]+", str(p.get("debye_waller") or "").strip()) if x]
    for x in dw:
        if not re.fullmatch(r"[A-Z][a-z]?=\d*\.?\d+", x):
            raise ValueError(f"Debye-Waller {x!r}: write El=B, e.g. Zr=0.3")
    if dw:
        a += ["--debye-waller", *dw]
    sp = p.get("species") or {}
    mode_sp = sp.get("mode", "none")
    if mode_sp == "cations":
        el = str(sp.get("to") or "")
        if not EL.match(el):
            raise ValueError(f"replace cations by {el!r}: give an element symbol")
        an = [x for x in re.split(r"[\s,]+", str(sp.get("anions") or "O").strip()) if x]
        if any(not EL.match(x) for x in an):
            raise ValueError("anions: element symbols, e.g. O")
        a += ["--replace-cations", el, "--anions", *an]
    elif mode_sp == "map":
        mp = [x for x in re.split(r"[\s,]+", str(sp.get("map") or "").strip()) if x]
        if not mp or any(not re.fullmatch(r"[A-Z][a-z]?[:=][A-Z][a-z]?", x) for x in mp):
            raise ValueError("species map: write OLD:NEW pairs, e.g. Ce:Zr")
        a += ["--replace", *mp]
    refs = [str(r) for r in p.get("refs") or []]
    if any(not NAME.match(r) for r in refs):
        raise ValueError("bad reference name")
    if refs:
        rd = str(Path(str(p.get("refs_dir") or APP / "refs")).expanduser())
        a += ["--ref", *refs, "--refs", rd]
        if not p.get("ref_replace", True):
            a.append("--no-ref-replace")
        if p.get("ref_scale_volume"):
            a.append("--ref-scale-volume")
    jobs = int(_num(b.get("jobs") or 4, "jobs", 1, 256))
    a += ["--jobs", str(jobs)]
    if mode == "sticks":
        key = hashlib.sha1(json.dumps(a).encode()).hexdigest()[:16]
        d = state_dir / "xrd"
        a += ["--pad", "5", "--sticks-json", str(d / f"sticks_{key}.json"), "--cache-dir", str(state_dir / "xrd_cache")]
        return a, key
    # export: broadening + output files
    br = b.get("broadening") or {}
    a += ["--step", f"{_num(br.get('step', 0.02), 'step', 1e-5, 10):g}"]
    if br.get("caglioti"):
        U, V, W = (_num(x, "Caglioti U V W") for x in br["caglioti"])
        a += ["--caglioti", f"{U:g}", f"{V:g}", f"{W:g}"]
    else:
        a += ["--fwhm", f"{_num(br.get('fwhm', 0.1), 'FWHM', 1e-5, 30):g}"]
    prof = br.get("profile", "gaussian")
    if prof not in ("gaussian", "lorentzian", "pseudo-voigt"):
        raise ValueError("profile")
    a += ["--profile", prof, "--eta", f"{_num(br.get('eta', 0.5), 'eta', 0, 1):g}"]
    norm = br.get("norm", "max")
    if norm not in ("max", "absolute"):
        raise ValueError("norm")
    a += ["--norm", norm]
    od = str(Path(str(b.get("out_dir") or "")).expanduser())
    if not od or not Path(od).is_absolute():
        raise ValueError("choose an output folder")
    a += ["--out-dir", od, "--formats", "png", "pdf"]
    return a, None


def sticks_path(state_dir: Path, key: str):
    if not re.fullmatch(r"[0-9a-f]{16}", key or ""):
        return None
    return state_dir / "xrd" / f"sticks_{key}.json"
