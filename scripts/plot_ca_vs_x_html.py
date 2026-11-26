#!/usr/bin/env python3
"""
plot_ca_vs_x_html.py - one self-contained HTML page of c/a vs composition with P(c/a) at every composition,
live controls, no server.

    python plot_ca_vs_x_html.py campaign/free_energy_300K_all.csv
    python plot_ca_vs_x_html.py all.csv --xrange 0 40 --quantities ca theta a_n c_n vol
    python plot_ca_vs_x_html.py all.csv --repeats 2 2 2 --out ca_vs_x.html

Reads the merged CSV of campaign.py (free_energy_<T>K_all.csv: x_dopant, phase, a_A b_A c_A, alpha_deg
beta_deg gamma_deg, ...), works out the cell shape of every structure, fits the peaks of P(c/a) of every phase
at every composition, and writes <stem>_ca_vs_x.html next to the CSV (or --out).  Data and fits are stored in
that one file: copy it to your own computer and double click it - no Python, no server, no internet needed.

Same page as plot_F_vs_x_html.py: sideways histograms P(c/a) at each composition (one colour per phase) with
Gaussian-mixture fits and their peaks, peaks joined across compositions, one panel or one row per phase.
The controls redraw at once; Save SVG / PNG; the page shows the command that writes it with the same settings.

c/a definition (identical to plot_cell_distributions.py / plot_cell_annotations.py, imported from there)
    Every lattice-vector length is divided by its supercell repeat (repeats found from the cation count and
    the lengths, or --repeats), giving conventional fluorite-cell parameters a_n, b_n, c_n whatever the
    supercell shape.  c = the longest normalised axis, decided per structure (it may lie along x, y or z),
    a_n = mean of the other two, c/a = c_n / a_n >= 1.  For t that is the tetragonal main axis; for oIII and m
    the longest axis over the mean of the two shorter ones.  --c-axis unique takes the axis most different
    from the other two instead (allows c/a < 1; splits near-cubic cells into two bumps around 1).
    Which Cartesian axis was c is counted per composition and phase and shown in the tooltips.
    Other quantities (--quantities): theta = cell angle farthest from 90 deg (folded >= 90, the monoclinic
    beta), a_n, c_n (Angstrom), vol = volume per atom (Angstrom^3).

Peaks: plot_free_energy.fit_peaks (1..3 Gaussians by BIC, maxima of the fit, bootstrap check); fit options are
fixed when the page is written.  Rows with imaginary modes or an unconverged relaxation are left out
(--include-flagged keeps them).  Needs plot_cell_annotations.py and plot_free_energy.py in the same folder;
no other script is changed.
"""
import argparse
import csv
import hashlib
import json
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from plot_free_energy import fit_peaks                     # noqa: E402  same peak method as the other scripts
    from plot_cell_annotations import cell_shape, detect_repeats, n_cations   # noqa: E402  same c/a definition
except ImportError as e:
    sys.exit(f"plot_free_energy.py and plot_cell_annotations.py must be in the same folder as this script ({e})")

PHASE_STYLE = {"parent": dict(color="#52514e", marker="o", label="c (initial)"),
               "t": dict(color="#2a78d6", marker="s", label="t"),
               "m": dict(color="#eb6834", marker="^", label="m"),
               "oIII": dict(color="#1baf7a", marker="D", label="oIII")}
EXTRA_STYLES = [dict(color="#6250d6", marker="v"), dict(color="#e34948", marker="P"),
                dict(color="#eda100", marker="X"), dict(color="#e87ba4", marker="h")]
QUANTS = {"ca": dict(name="c/a", unit="", dec=4, ref=1.0),
          "theta": dict(name="θ", unit="deg", dec=3, ref=90.0),
          "a_n": dict(name="a_n", unit="Å", dec=4, ref=None),
          "c_n": dict(name="c_n", unit="Å", dec=4, ref=None),
          "vol": dict(name="V/atom", unit="Å³", dec=4, ref=None)}
CELL = ["a_A", "b_A", "c_A", "alpha_deg", "beta_deg", "gamma_deg"]


def flagged(r):
    imag = int(float(r.get("n_imaginary_modes") or 0)) > 0
    return imag or r.get("relax_converged") == "False"


def load(path, include_flagged, repeats, c_mode):
    if not Path(path).is_file():
        sys.exit(f"{path}: no such file")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or any(c not in rows[0] for c in CELL):
        sys.exit(f"{path}: no cell columns ({', '.join(CELL)}) - not a batch_free_energy / campaign CSV")
    if "x_dopant" not in rows[0]:
        sys.exit(f"{path}: no x_dopant column - use the merged free_energy_<T>K_all.csv from campaign.py merge")
    data, nflag, nbad, reps = {}, {}, 0, {}
    for r in rows:
        if r.get("status", "ok") != "ok" or r["x_dopant"] == "":
            continue
        try:
            L = np.array([float(r[c]) for c in CELL[:3]])
            A = np.array([float(r[c]) for c in CELL[3:]])
        except ValueError:
            nbad += 1
            continue
        ncat = n_cations(r)
        rep = tuple(repeats) if repeats else detect_repeats(L, ncat)
        if rep is None or (ncat and 4 * int(np.prod(rep)) != ncat):     # --repeats must fit the cation count
            nbad += 1
            continue
        x, ph = float(r["x_dopant"]), r["phase"] or "parent"
        if flagged(r):
            nflag[(x, ph)] = nflag.get((x, ph), 0) + 1
            if not include_flagged:
                continue
        sh = cell_shape(L, A, rep, c_mode)
        try:
            F = float(r.get("F_total_eV_atom") or "nan")
        except ValueError:
            F = float("nan")
        try:
            vol = float(r.get("volume_A3_atom") or "nan")
        except ValueError:
            vol = float("nan")
        if not np.isfinite(vol):
            vol = float(abs(np.linalg.det(cell_matrix(L, A)))) / float(r["n_atoms"]) if r.get("n_atoms") else np.nan
        data.setdefault((x, ph), []).append(dict(ca=sh["ca"], theta=sh["theta"], a_n=sh["a_n"], c_n=sh["c_n"],
                                                 vol=vol, F=F, cax=sh["c_axis"]))
        reps.setdefault((x, ph), set()).add("×".join(map(str, rep)))
    if not data:
        hint = f" ({sum(nflag.values())} flagged rows left out: try --include-flagged)" if nflag else ""
        hint += f" ({nbad} rows without cell data or repeats: try --repeats)" if nbad else ""
        sys.exit(f"{path}: no usable rows{hint}")
    T = next((r["T_K"] for r in rows if r.get("T_K")), "")
    return data, nflag, nbad, reps, T


def cell_matrix(L, A):
    """lattice vectors from lengths and angles (only for the volume when the CSV has none)"""
    a, b, c = L
    al, be, ga = np.radians(A)
    v1 = [a, 0, 0]
    v2 = [b * np.cos(ga), b * np.sin(ga), 0]
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    return np.array([v1, v2, [cx, cy, np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))]])


def _fit_job(job):
    key, y, seed, kmax, kfix, min_peak, n_boot = job
    return key, fit_peaks(y, np.random.default_rng(seed), kmax=kmax, k_fixed=kfix, min_peak=min_peak, n_boot=n_boot)


def default_jobs():
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    return max(1, min(n, 8))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="merged campaign CSV (free_energy_<T>K_all.csv)")
    ap.add_argument("--out", help="HTML file to write (default: <stem>_ca_vs_x.html next to the CSV)")
    g = ap.add_argument_group("cell shape and peak fits (fixed in the page)")
    g.add_argument("--quantities", nargs="+", choices=list(QUANTS), default=["ca"],
                   help="quantities stored in the page (default: ca); e.g. ca theta a_n c_n vol")
    g.add_argument("--c-axis", choices=["longest", "unique"], default="longest",
                   help="c = longest normalised axis (default) or the axis most different from the other two")
    g.add_argument("--repeats", type=int, nargs=3, metavar=("N1", "N2", "N3"),
                   help="supercell repeats for ALL structures (default: automatic, per structure)")
    g.add_argument("--max-components", type=int, default=3, help="most Gaussians per distribution (default 3)")
    g.add_argument("--components", action="append", default=[], metavar="[PHASE=]K",
                   help="force K Gaussians (all phases, or one phase)")
    g.add_argument("--min-peak", type=float, default=0.15,
                   help="a peak must reach this fraction of the highest one (default 0.15)")
    g.add_argument("--bootstrap", type=int, default=20,
                   help="bootstrap re-fits per distribution (default 20; 0 = fast, no peak errors)")
    g.add_argument("--seed", type=int, default=0)
    g.add_argument("--include-flagged", action="store_true",
                   help="keep rows with imaginary modes / unconverged relaxation")
    g.add_argument("--jobs", type=int, default=default_jobs(), help="parallel processes for the fits")
    g.add_argument("--no-cache", dest="cache", action="store_false",
                   help="do not reuse / save the fits (<stem>_ca_vs_x_html_fitcache.pkl)")
    p = ap.add_argument_group("start values of the page controls")
    p.add_argument("--quantity", choices=list(QUANTS), help="quantity shown first (default: the first stored)")
    p.add_argument("--xrange", type=float, nargs=2, metavar=("LO", "HI"), help="composition range, %% of cations")
    p.add_argument("--phases", nargs="+", help="phases ticked at the start (default: all)")
    p.add_argument("--label", action="append", default=[], metavar="PHASE=TEXT", help="legend text of a phase")
    p.add_argument("--style", choices=["hist", "band"], default="hist")
    p.add_argument("--layout", choices=["overlay", "rows"], default="overlay")
    p.add_argument("--norm", choices=["each", "common"], default="each")
    p.add_argument("--side", choices=["left", "right"], default="left")
    p.add_argument("--width", type=float, default=0.85)
    p.add_argument("--dodge", type=float, default=0.0)
    p.add_argument("--tail", type=float, default=0.02)
    p.add_argument("--bin", type=float, default=0.0, help="bin width of the shown quantity (default: automatic)")
    p.add_argument("--ylim", type=float, nargs=2, metavar=("LO", "HI"))
    p.add_argument("--no-lines", action="store_true")
    p.add_argument("--no-fit", action="store_true")
    p.add_argument("--no-bars", action="store_true")
    p.add_argument("--show-lowest-f", action="store_true", help="mark the structure with the lowest F")
    p.add_argument("--no-ref-line", action="store_true", help="no dotted line at c/a = 1 (θ = 90°)")
    p.add_argument("--same-y", action="store_true")
    p.add_argument("--fraction", action="store_true")
    p.add_argument("--dopant", default="Ce")
    p.add_argument("--font", choices=["serif", "sans"], default="serif")
    a = ap.parse_args()

    quants = list(dict.fromkeys(a.quantities))
    if a.quantity and a.quantity not in quants:
        quants.append(a.quantity)
    labels, k_all, k_phase = {}, None, {}
    for item in a.label:
        if "=" not in item:
            sys.exit(f"--label needs PHASE=TEXT, got {item!r}")
        k, v = item.split("=", 1); labels[k] = v
    for item in a.components:
        if "=" in item:
            q, v = item.split("=", 1); k_phase[q] = int(v)
        else:
            k_all = int(item)

    data, nflag, nbad, reps, T = load(a.csv, a.include_flagged, a.repeats, a.c_axis)
    present = {q for (x, q) in data}
    phases = [q for q in PHASE_STYLE if q in present] + sorted(present - set(PHASE_STYLE))
    if a.phases and any(q not in phases for q in a.phases):
        sys.exit(f"phases {[q for q in a.phases if q not in phases]} not in the CSV (have {phases})")
    style, k = {}, 0
    for q in phases:
        if q in PHASE_STYLE:
            style[q] = dict(PHASE_STYLE[q])
        else:
            style[q] = dict(EXTRA_STYLES[k % len(EXTRA_STYLES)], label=q); k += 1

    # ------------------------------------------------------------ fits (on value - lowest value)
    csvp = Path(a.csv).resolve()
    outp = Path(a.out) if a.out else csvp.parent / f"{csvp.stem}_ca_vs_x.html"
    cache_p = outp.parent / f"{csvp.stem}_ca_vs_x_html_fitcache.pkl"
    cache = {}
    if a.cache and cache_p.exists():
        try:
            cache = pickle.loads(cache_p.read_bytes())
        except Exception:
            cache = {}
    keys = sorted(data, key=lambda t: (t[0], phases.index(t[1])))
    dist, jobs, queued = {}, [], set()
    for (x, ph) in keys:
        kfix = k_phase.get(ph, k_all)
        opts = f"{a.max_components}|{kfix}|{a.min_peak}|{a.bootstrap}|{a.seed}".encode()
        for qn in quants:
            vals = np.array([s[qn] for s in data[(x, ph)]], float)
            vals = np.sort(vals[np.isfinite(vals)])
            if not len(vals):
                continue
            low = float(vals[0]); y0 = np.round(vals - low, 7)
            h = hashlib.sha1(b"cahtml1" + np.ascontiguousarray(y0).tobytes() + opts).hexdigest()
            dist[(x, ph, qn)] = (low, y0, h)
            if h not in cache and h not in queued:
                queued.add(h)
                jobs.append((h, y0, (int(h[:12], 16) + a.seed) % 2 ** 63, a.max_components, kfix, a.min_peak,
                             a.bootstrap))
    nuniq = len({v[2] for v in dist.values()})
    print(f"peak fits: {len(dist)} distributions ({nuniq} different), {nuniq - len(jobs)} from cache, "
          f"{len(jobs)} to fit", flush=True)
    t0 = time.time(); done = 0
    def got(rr):
        nonlocal done
        cache[rr[0]] = rr[1]; done += 1
        if done == len(jobs) or done % max(1, len(jobs) // 10) == 0:
            print(f"  {done}/{len(jobs)} fitted, {time.time() - t0:.0f} s", flush=True)
    if len(jobs) > 1 and a.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(min(a.jobs, len(jobs))) as pool:
            for rr in pool.imap_unordered(_fit_job, jobs):
                got(rr)
    else:
        for j in jobs:
            got(_fit_job(j))
    if a.cache and jobs:
        cache_p.write_bytes(pickle.dumps(cache))

    dists = []
    for (x, ph) in keys:
        S = data[(x, ph)]
        Fs = np.array([s["F"] for s in S])
        ilow = int(np.nanargmin(Fs)) if np.isfinite(Fs).any() else None
        cax = [sum(1 for s in S if s["cax"] == i) for i in range(3)]
        Q = {}
        for qn in quants:
            if (x, ph, qn) not in dist:
                continue
            low, y0, h = dist[(x, ph, qn)]
            f = cache[h]
            peaks = sorted(f["peaks"], key=lambda t: t[0])
            atF = S[ilow][qn] if ilow is not None else None
            Q[qn] = dict(low=round(low, 7), v=[round(float(t), 7) for t in y0], k=int(f["k"]),
                         deg=bool(f.get("degenerate")), w=[float(t) for t in f["w"]],
                         mu=[float(t) for t in f["mu"]], sg=[float(t) for t in f["sigma"]],
                         peaks=[[float(t[0]), float(t[1]), float(t[2])] for t in peaks],
                         main=int(np.argmax([t[2] for t in peaks])),
                         atF=(round(float(atF), 7) if atF is not None and np.isfinite(atF) else None))
        dists.append(dict(x=round(100 * x, 6), p=phases.index(ph), n=len(S), nflag=nflag.get((x, ph), 0),
                          cax=cax, reps=sorted(reps[(x, ph)]), Q=Q))
    page = dict(
        csv=csvp.name, T=(f"{float(T):g}" if T else ""), nf=sum(nflag.values()), flagged_in=a.include_flagged,
        nbad=nbad, c_axis=a.c_axis, repeats=a.repeats,
        quants=[dict(id=q, **QUANTS[q]) for q in quants],
        phases=[dict(id=q, label=labels.get(q, style[q]["label"]), deflabel=style[q]["label"],
                     color=style[q]["color"], marker=style[q]["marker"]) for q in phases],
        fit=dict(max_components=a.max_components, components=a.components, min_peak=a.min_peak,
                 bootstrap=a.bootstrap, seed=a.seed),
        init=dict(quantity=a.quantity or quants[0], layout=a.layout, style=a.style, norm=a.norm, side=a.side,
                  font=a.font, width=a.width, dodge=a.dodge, tail=a.tail, bin=a.bin, xrange=a.xrange,
                  phases=a.phases, ylim=a.ylim, no_fit=a.no_fit, no_bars=a.no_bars, no_lines=a.no_lines,
                  show_lowest_f=a.show_lowest_f, no_ref_line=a.no_ref_line, same_y=a.same_y,
                  fraction=a.fraction, dopant=a.dopant),
        dists=dists)
    blob = json.dumps(page, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(PAGE.replace("__DATA__", blob), encoding="utf-8")
    nrow = sum(len(v) for v in data.values())
    print(f"{nrow} rows, {len({x for x, q in data})} compositions, phases {phases}; flagged rows "
          f"{'included' if a.include_flagged else 'left out'}: {sum(nflag.values())}"
          + (f"; {nbad} rows skipped (no cell data / repeats not found: use --repeats)" if nbad else ""))
    print(f"c/a = c_n/a_n, c = {a.c_axis} normalised axis per structure, a_n = mean of the other two; "
          f"repeats {'×'.join(map(str, a.repeats)) if a.repeats else 'automatic'}")
    print(f"wrote {outp} ({outp.stat().st_size / 1e6:.1f} MB): open it in any browser, no server needed")


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>c/a vs composition</title>
<style>
:root{--bg:#f6f5f2;--panel:#ffffff;--ink:#21201d;--ink2:#6b6964;--line:#dcd9d2;--acc:#2a78d6;--bad:#c0392b;--code:#f0eee9}
@media (prefers-color-scheme: dark){:root{--bg:#1b1b1a;--panel:#252523;--ink:#eceae5;--ink2:#a3a09a;--line:#3a3936;--acc:#6aa7ee;--bad:#ff8a7a;--code:#2e2e2b}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif}
header{padding:10px 16px;border-bottom:1px solid var(--line);display:flex;gap:12px;align-items:baseline;flex-wrap:wrap}
header h1{font-size:16px;margin:0;font-weight:600}
header span{color:var(--ink2)}
main{display:grid;grid-template-columns:330px 1fr;gap:16px;padding:16px;align-items:start}
@media (max-width:860px){main{grid-template-columns:1fr}}
#controls{display:flex;flex-direction:column;gap:10px;max-height:calc(100vh - 80px);overflow:auto;padding-right:4px}
@media (max-width:860px){#controls{max-height:none}}
fieldset{background:var(--panel);border:1px solid var(--line);border-radius:8px;margin:0;padding:10px 12px 12px;min-width:0}
legend{font-weight:600;padding:0 6px;font-size:13px}
.row{display:grid;grid-template-columns:112px 1fr;gap:8px;align-items:center;margin-top:7px}
.row>label{color:var(--ink2);font-size:13px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:6px}
input[type=number],input[type=text],select{width:100%;padding:4px 6px;border:1px solid var(--line);border-radius:5px;
  background:var(--bg);color:var(--ink);font:inherit}
input[type=range]{width:100%;accent-color:var(--acc)}
.sl{display:grid;grid-template-columns:1fr 44px;gap:6px;align-items:center}
.sl output{font-variant-numeric:tabular-nums;text-align:right;color:var(--ink2);font-size:13px}
.checks{display:flex;flex-wrap:wrap;gap:4px 14px;margin-top:7px}
.checks label{display:flex;gap:5px;align-items:center;font-size:13px}
input[type=checkbox]{accent-color:var(--acc)}
.note{color:var(--ink2);font-size:12px;margin:6px 0 0}
#view{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px;min-width:0}
#bar{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:10px}
#status{color:var(--ink2);flex:1;min-width:160px}
#status.bad{color:var(--bad)}
button{font:inherit;padding:5px 12px;border-radius:6px;border:1px solid var(--line);background:var(--bg);color:var(--ink);cursor:pointer}
button.main{background:var(--acc);border-color:var(--acc);color:#fff}
#plotbox{background:#fff;border-radius:6px;overflow-x:auto}
#plotbox svg{display:block;margin:0 auto;max-width:100%;height:auto}
pre{background:var(--code);border-radius:6px;padding:8px 10px;margin:10px 0 0;white-space:pre-wrap;word-break:break-word;
  font:12px/1.45 ui-monospace,Menlo,Consolas,monospace}
h2{font-size:12px;margin:14px 0 0;color:var(--ink2);font-weight:600;text-transform:uppercase;letter-spacing:.04em}
</style></head><body>
<header><h1>c/a vs composition</h1><span id="file"></span></header>
<main>
<form id="controls" onsubmit="return false">
<fieldset><legend>Data</legend>
  <div class="row"><label for="quantity">y axis</label><select id="quantity"></select></div>
  <div class="row"><label>x range (%)</label><div class="two"><input type="number" id="xlo" step="any" placeholder="from"><input type="number" id="xhi" step="any" placeholder="to"></div></div>
  <div class="row"><label>phases</label><div class="checks" id="phases" style="margin-top:0"></div></div>
</fieldset>
<fieldset><legend>Look</legend>
  <div class="row"><label for="style">style</label><select id="style"><option value="hist">histograms</option><option value="band">bands</option></select></div>
  <div class="row"><label for="layout">layout</label><select id="layout"><option value="overlay">one panel</option><option value="rows">one row per phase</option></select></div>
  <div class="row"><label for="norm">histogram length</label><select id="norm"><option value="each">each the same</option><option value="common">one P scale</option></select></div>
  <div class="row"><label for="side">side</label><select id="side"><option value="left">left of x</option><option value="right">right of x</option></select></div>
  <div class="row"><label for="width">width</label><div class="sl"><input type="range" id="width" min="0.2" max="1.5" step="0.05"><output></output></div></div>
  <div class="row"><label for="dodge">dodge</label><div class="sl"><input type="range" id="dodge" min="0" max="0.5" step="0.02"><output></output></div></div>
  <div class="row"><label for="tail">hide tail below</label><div class="sl"><input type="range" id="tail" min="0" max="0.2" step="0.005"><output></output></div></div>
  <div class="row"><label for="bin" id="binlab">bin</label><input type="number" id="bin" min="0" step="any" placeholder="automatic"></div>
  <div class="row"><label>y range</label><div class="two"><input type="number" id="ylo" step="any" placeholder="from"><input type="number" id="yhi" step="any" placeholder="to"></div></div>
  <div class="checks">
    <label><input type="checkbox" id="no_lines">no connecting lines</label>
    <label><input type="checkbox" id="no_fit">no fitted curves</label>
    <label><input type="checkbox" id="no_bars">no bars</label>
    <label><input type="checkbox" id="show_lowest_f">mark lowest-F structure</label>
    <label><input type="checkbox" id="no_ref_line">no line at c/a = 1</label>
    <label><input type="checkbox" id="same_y">same y in rows</label>
  </div>
</fieldset>
<fieldset><legend>Labels</legend>
  <div id="labels"></div>
  <div class="row"><label for="dopant">dopant name</label><input type="text" id="dopant" maxlength="20"></div>
  <div class="row"><label for="font">font</label><select id="font"><option value="serif">serif</option><option value="sans">sans</option></select></div>
  <div class="checks"><label><input type="checkbox" id="fraction">x as fraction</label></div>
</fieldset>
<fieldset><legend>Cell shape and peak fits (fixed in this file)</legend>
  <p class="note" id="fitnote" style="margin-top:2px"></p>
</fieldset>
</form>
<section id="view">
  <div id="bar"><span id="status"></span>
    <button id="reset" type="button">Reset</button>
    <button id="svgb" type="button">Save SVG</button>
    <button id="pngb" type="button" class="main">Save PNG</button></div>
  <div id="plotbox"></div>
  <h2>Command that writes this page with these settings</h2><pre id="cmd"></pre>
</section>
</main>
<script id="data" type="application/json">__DATA__</script>
<script>
"use strict";
const D = JSON.parse(document.getElementById("data").textContent);
const $ = id => document.getElementById(id);
const INK = "#1f1f1f", INK2 = "#6b6b6b";
const SL = ["width", "dodge", "tail"], SE = ["quantity", "style", "layout", "norm", "side", "font"];
const CH = ["no_lines", "no_fit", "no_bars", "show_lowest_f", "no_ref_line", "same_y", "fraction"];
const DEF = {quantity: D.quants[0].id, layout: "overlay", style: "hist", norm: "each", side: "left", font: "serif",
             width: 0.85, dodge: 0, tail: 0.02, dopant: "Ce"};
const fmt = (v, n) => (+v).toFixed(n === undefined ? 2 : n);
const g6 = v => String(+(+v).toPrecision(6));

// ---------------------------------------------------------------- controls
function setControls(I){
  SE.forEach(k => $(k).value = I[k]);
  SL.forEach(k => { $(k).value = I[k]; $(k).nextElementSibling.textContent = fmt(I[k], k === "tail" ? 3 : 2); });
  CH.forEach(k => $(k).checked = !!I[k]);
  $("bin").value = I.bin > 0 ? I.bin : "";
  const pair = (a, b, v) => { $(a).value = v ? v[0] : ""; $(b).value = v ? v[1] : ""; };
  pair("xlo", "xhi", I.xrange); pair("ylo", "yhi", I.ylim);
  $("dopant").value = I.dopant;
  D.phases.forEach((q, i) => { $("ph" + i).checked = !I.phases || I.phases.includes(q.id); $("lb" + i).value = q.label; });
}
function state(){
  const s = {};
  SE.forEach(k => s[k] = $(k).value);
  SL.forEach(k => s[k] = +$(k).value);
  CH.forEach(k => s[k] = $(k).checked);
  const num = id => $(id).value === "" ? null : +$(id).value;
  s.bin = num("bin"); s.xlo = num("xlo"); s.xhi = num("xhi");
  s.ylo = num("ylo"); s.yhi = num("yhi");
  s.phases = D.phases.map((q, i) => i).filter(i => $("ph" + i).checked);
  s.labels = D.phases.map((q, i) => $("lb" + i).value.trim() || q.label);
  s.dopant = $("dopant").value.trim() || "Ce";
  return s;
}

// ---------------------------------------------------------------- numerics (same rules as the Python script)
function pct(sorted, q){                       // numpy percentile, linear interpolation
  const h = (sorted.length - 1) * q / 100, i = Math.floor(h), j = Math.min(i + 1, sorted.length - 1);
  return sorted[i] + (sorted[j] - sorted[i]) * (h - i);
}
function commonBin(ds){
  const ws = [];
  ds.forEach(d => { if (d.v.length >= 5){ const iqr = pct(d.v, 75) - pct(d.v, 25); if (iqr > 0) ws.push(2 * iqr * Math.pow(d.v.length, -1 / 3)); } });
  let w = 0.5;
  if (ws.length){ ws.sort((a, b) => a - b); w = pct(ws, 50); }
  const k = Math.floor(Math.log10(w)), m = w / Math.pow(10, k);
  let best = 1, bd = Infinity;
  [1, 2, 2.5, 5, 10].forEach(c => { const dd = Math.abs(Math.log(c / m)); if (dd < bd){ bd = dd; best = c; } });
  return +(best * Math.pow(10, k)).toPrecision(12);
}
function halfMax(g, pdf, F0){
  let k = 0, bd = Infinity;
  for (let i = 0; i < g.length; i++){ const dd = Math.abs(g[i] - F0); if (dd < bd){ bd = dd; k = i; } }
  const half = 0.5 * pdf[k]; let i = k, j = k;
  while (i > 0 && pdf[i - 1] >= half && pdf[i - 1] <= pdf[i] + 1e-12) i--;
  while (j < g.length - 1 && pdf[j + 1] >= half && pdf[j + 1] <= pdf[j] + 1e-12) j++;
  return [g[i], g[j]];
}
function link(prev, cur){
  if (!prev.length || !cur.length) return [];
  const amin = (v, arr) => { let b = 0; arr.forEach((t, i) => { if (Math.abs(v - t) < Math.abs(v - arr[b])) b = i; }); return b; };
  const seg = new Set();
  cur.forEach((c, j) => seg.add(amin(c, prev) + "," + j));
  prev.forEach((p, i) => seg.add(i + "," + amin(p, cur)));
  return [...seg].map(t => t.split(",").map(Number));
}
function prepare(d, off, bw, tail){
  // d: one distribution; off: shift of its lowest value on the chosen y axis
  const n = d.v.length, r = {d: d, n: n, low: off, peaks: d.peaks.map(q => [q[0] + off, q[1], q[2]]), main: d.main};
  const ymin = off, ymax = d.v[n - 1] + off;
  const lo = (Math.floor(ymin / bw + 0.5) - 0.5) * bw, hi = (Math.floor(ymax / bw + 0.5) + 0.5) * bw;
  const nb = Math.max(1, Math.round((hi - lo) / bw));
  const edges = Array.from({length: nb + 1}, (_, i) => lo + i * bw), prob = new Array(nb).fill(0);
  d.v.forEach(t => { let b = Math.floor((t + off - lo) / bw); if (b < 0) b = 0; if (b >= nb) b = nb - 1; prob[b] += 1 / n; });
  const flo = bw / Math.sqrt(2 * Math.PI), sg = d.sg.map(t => d.deg ? flo : Math.max(t, flo));
  const span = 3 * Math.max(...sg), g0 = edges[0] - span, g1 = edges[nb] + span, N = 300;
  let g = [], pdf = [];
  for (let i = 0; i < N; i++){
    const t = g0 + (g1 - g0) * i / (N - 1); let s = 0;
    for (let c = 0; c < sg.length; c++){ const z = (t - d.mu[c] - off) / sg[c]; s += d.w[c] * Math.exp(-0.5 * z * z) / (sg[c] * Math.sqrt(2 * Math.PI)); }
    g.push(t); pdf.push(s * bw);
  }
  const pm = Math.max(...pdf); let j0 = 0, j1 = N - 1;
  while (j0 < N - 1 && pdf[j0] < 0.02 * pm) j0++;
  while (j1 > 0 && pdf[j1] < 0.02 * pm) j1--;
  g = g.slice(j0, j1 + 1); pdf = pdf.slice(j0, j1 + 1);
  r.fwhm = r.peaks.map(q => halfMax(g, pdf, q[0]));
  r.nclip = 0;
  let P = prob.slice(), C = pdf.slice();
  if (n > 1 && tail > 0){
    const thr = tail * Math.max(...prob); let offb = prob.map(t => t < thr);
    for (let k = 0; k < nb;){
      if (offb[k]){ k++; continue; }
      let j = k, s = 0; while (j < nb && !offb[j]){ s += prob[j]; j++; }
      if (s < 0.01) for (let t = k; t < j; t++) offb[t] = true;
      k = j;
    }
    if (offb.every(t => t)) offb = prob.map(t => t < thr);
    r.nclip = Math.round(prob.reduce((s, t, i) => s + (offb[i] ? t : 0), 0) * n);
    P = prob.map((t, i) => offb[i] ? null : t);
    const cm = Math.max(...pdf);
    C = pdf.map((t, i) => {
      if (t < tail * cm) return null;
      for (let b = 0; b < nb; b++) if (!offb[b] && g[i] >= edges[b] - bw && g[i] <= edges[b + 1] + bw) return t;
      return null;
    });
    if (C.every(t => t === null)) C = pdf.slice();
  }
  const on = P.map((t, i) => t === null ? -1 : i).filter(i => i >= 0), gon = g.filter((t, i) => C[i] !== null);
  r.span = [Math.min(edges[on[0]], gon[0]), Math.max(edges[on[on.length - 1] + 1], gon[gon.length - 1])];
  r.edges = edges; r.prob = P; r.g = g; r.pdf = C;
  r.top = Math.max(...P.filter(t => t !== null), ...C.filter(t => t !== null));
  return r;
}

// ---------------------------------------------------------------- svg helpers
function niceTicks(lo, hi, target){
  const span = hi - lo, raw = span / target, k = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / k;
  const step = (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * k, out = [];
  for (let t = Math.ceil(lo / step - 1e-9) * step; t <= hi + 1e-9 * span; t += step) out.push(Math.abs(t) < step * 1e-9 ? 0 : t);
  const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return {ticks: out, dec: dec};
}
const esc = t => String(t).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
function marker(kind, x, y, s, fill, stroke, sw, title){
  let el; const a = `fill="${fill}" stroke="${stroke}" stroke-width="${sw}"`;
  const poly = pts => `<polygon points="${pts.map(p => (x + p[0] * s).toFixed(2) + "," + (y + p[1] * s).toFixed(2)).join(" ")}" ${a}>`;
  if (kind === "s") el = `<rect x="${(x - 0.82 * s).toFixed(2)}" y="${(y - 0.82 * s).toFixed(2)}" width="${(1.64 * s).toFixed(2)}" height="${(1.64 * s).toFixed(2)}" ${a}>`;
  else if (kind === "^") el = poly([[0, -1.1], [1, 0.75], [-1, 0.75]]);
  else if (kind === "v") el = poly([[0, 1.1], [1, -0.75], [-1, -0.75]]);
  else if (kind === "D") el = poly([[0, -1.15], [1.15, 0], [0, 1.15], [-1.15, 0]]);
  else if (kind === "P") el = poly([[-.4, -1], [.4, -1], [.4, -.4], [1, -.4], [1, .4], [.4, .4], [.4, 1], [-.4, 1], [-.4, .4], [-1, .4], [-1, -.4], [-.4, -.4]]);
  else if (kind === "X") el = poly([[0, -.5], [.6, -1.1], [1.1, -.6], [.5, 0], [1.1, .6], [.6, 1.1], [0, .5], [-.6, 1.1], [-1.1, .6], [-.5, 0], [-1.1, -.6], [-.6, -1.1]]);
  else if (kind === "h") el = poly([[0, -1.1], [.95, -.55], [.95, .55], [0, 1.1], [-.95, .55], [-.95, -.55]]);
  else el = `<circle cx="${x.toFixed(2)}" cy="${y.toFixed(2)}" r="${(0.95 * s).toFixed(2)}" ${a}>`;
  return el + (title ? `<title>${esc(title)}</title>` : "") + el.replace(/^<(\w+).*$/, "</$1>");
}

// ---------------------------------------------------------------- the plot
const qname = q => q.name + (q.unit ? ` (${q.unit})` : "");
const AX = "xyz";
function qlabel(q){                              // y-axis label as svg tspans
  const it = t => `<tspan font-style="italic">${t}</tspan>`, sub = t => `<tspan font-size="11" dy="4">${t}</tspan>`;
  const u = q.unit ? `<tspan dy="-4"> (${q.unit})</tspan>` : "";
  if (q.id === "ca") return it("c") + sub("n") + `<tspan dy="-4">/</tspan>` + it("a") + sub("n");
  if (q.id === "theta") return it("θ") + (q.unit ? ` (${q.unit})` : "");
  if (q.id === "a_n") return it("a") + sub("n") + u;
  if (q.id === "c_n") return it("c") + sub("n") + u;
  if (q.id === "vol") return it("V") + `/atom (${q.unit})`;
  return esc(qname(q));
}
function build(s){
  const ph = s.phases; if (!ph.length) throw "choose at least one phase";
  const Q = D.quants.find(q => q.id === s.quantity), dec = Q.dec, un = Q.unit ? " " + Q.unit : "";
  const all = D.dists.filter(d => ph.includes(d.p) && d.Q[Q.id]).map(d => Object.assign({}, d.Q[Q.id], {x: d.x, p: d.p, n: d.n, nflag: d.nflag, cax: d.cax, reps: d.reps}));
  if (!all.length) throw "no data for these phases";
  if ((s.xlo === null) !== (s.xhi === null)) throw "x range: give both values or none";
  if (s.xlo !== null && !(s.xlo < s.xhi)) throw "x range: the first value must be below the second";
  const use = s.xlo === null ? all : all.filter(d => d.x >= s.xlo - 1e-9 && d.x <= s.xhi + 1e-9);
  if (!use.length) throw "no compositions inside the x range";
  const phases = ph.filter(p => use.some(d => d.p === p));
  const comps = [...new Set(use.map(d => d.x))].sort((a, b) => a - b);
  const xs = s.fraction ? 0.01 : 1;
  let dx = 1; if (comps.length > 1){ dx = Infinity; for (let i = 1; i < comps.length; i++) dx = Math.min(dx, comps[i] - comps[i - 1]); dx *= xs; }
  const bw = s.bin > 0 ? s.bin : commonBin(use);
  const res = {}, key = (x, p) => x + "|" + p;
  phases.forEach(p => comps.forEach(x => { const d = use.find(t => t.x === x && t.p === p); if (d) res[key(x, p)] = prepare(d, d.low, bw, s.tail); }));
  const R = Object.values(res), rows = s.layout === "rows", band = s.style === "band", sgn = s.side === "left" ? -1 : 1;
  const slot = s.width * dx;
  if (s.norm === "each") R.forEach(r => r.sc = slot / r.top);
  else { const tops = R.map(r => r.top).sort((a, b) => a - b), pref = pct(tops, 90); R.forEach(r => r.sc = slot / pref); }
  const shift = {}; phases.forEach((p, i) => shift[p] = rows ? 0 : (i - (phases.length - 1) / 2) * s.dodge * dx);
  const col = p => D.phases[p].color, mk = p => D.phases[p].marker;

  // geometry
  const W = 900, ML = 84, MR = 44, PW = W - ML - MR, GAP = 14, serif = s.font === "serif";
  const FF = serif ? "'STIX Two Text','STIXGeneral','Times New Roman',Times,serif" : "Helvetica,Arial,sans-serif";
  const nA = rows ? phases.length : 1, hA = rows ? Math.max(120, Math.min(190, 620 / nA)) : 440;
  const top = 46; let y0 = top; const panels = [];
  for (let i = 0; i < nA; i++){ panels.push({y: y0, h: hA, p: rows ? phases[i] : null}); y0 += hA + GAP; }
  const bottom = y0 - GAP;
  const dd = rows ? 0 : s.dodge * dx * phases.length / 2, far = (band ? 0.4 * dx : 1.15 * slot) + dd, near = 0.4 * dx + dd;
  const xmin = comps[0] * xs - (sgn < 0 ? far : near), xmax = comps[comps.length - 1] * xs + (sgn > 0 ? far : near);
  const X = v => ML + (v - xmin) / (xmax - xmin) * PW;
  const yr = rr => {
    let lo, hi;
    if (band){ lo = Math.min(...rr.flatMap(r => r.fwhm.map(f => f[0]))); hi = Math.max(...rr.flatMap(r => r.fwhm.map(f => f[1]))); }
    else { lo = Math.min(...rr.map(r => r.n > 1 ? r.span[0] : r.low)); hi = Math.max(...rr.map(r => r.n > 1 ? r.span[1] : r.low)); }
    if (s.show_lowest_f) rr.forEach(r => { if (r.d.atF !== null){ lo = Math.min(lo, r.d.atF); hi = Math.max(hi, r.d.atF); } });
    const pd = 0.06 * ((hi - lo) || Math.abs(hi) * 1e-3 || 1); return [lo - pd, hi + pd];
  };
  if ((s.ylo === null) !== (s.yhi === null) || (s.ylo !== null && !(s.ylo < s.yhi))) throw "y range: give two values, the first below the second";
  let o = "", defs = "", cid = 0;
  const line = (x1, y1, x2, y2, c, w, dash) => `<line x1="${x1.toFixed(2)}" y1="${y1.toFixed(2)}" x2="${x2.toFixed(2)}" y2="${y2.toFixed(2)}" stroke="${c}" stroke-width="${w}"${dash ? ` stroke-dasharray="${dash}"` : ""}/>`;
  const frame = (pn, yl, last) => {
    const Y = v => pn.y + pn.h - (v - yl[0]) / (yl[1] - yl[0]) * pn.h;
    let t = "";
    const xt = niceTicks(xmin, xmax, 7), yt = niceTicks(yl[0], yl[1], Math.max(3, Math.round(pn.h / 70)));
    xt.ticks.forEach(v => { const px = X(v); t += line(px, pn.y + pn.h, px, pn.y + pn.h - 5, INK, 0.9) + line(px, pn.y, px, pn.y + 5, INK, 0.9);
      if (last) t += `<text x="${px.toFixed(1)}" y="${pn.y + pn.h + 18}" text-anchor="middle" font-size="14">${v.toFixed(xt.dec)}</text>`; });
    yt.ticks.forEach(v => { const py = Y(v); t += line(ML, py, ML + 5, py, INK, 0.9) + line(ML + PW, py, ML + PW - 5, py, INK, 0.9) +
      `<text x="${ML - 8}" y="${(py + 4.8).toFixed(1)}" text-anchor="end" font-size="14">${(v < 0 ? "−" : "") + Math.abs(v).toFixed(yt.dec)}</text>`; });
    t += `<rect x="${ML}" y="${pn.y}" width="${PW}" height="${pn.h}" fill="none" stroke="${INK}" stroke-width="1.1"/>`;
    return {Y: Y, deco: t};
  };
  const clip = pn => { const id = "c" + (cid++); defs += `<clipPath id="${id}"><rect x="${ML}" y="${pn.y}" width="${PW}" height="${pn.h}"/></clipPath>`; return id; };
  const DASH = "4.2 2.8";
  const axes = r => { const c = r.d.cax, n = c[0] + c[1] + c[2]; return AX.split("").map((a, i) => c[i] ? `${a} ${c[i]}` : "").filter(t => t).join(", ") + (n ? "" : "-"); };
  const peakMarks = (r, x, p, px, Y) => {
    let t = "";
    r.peaks.forEach((q, i) => {
      const main = i === r.main, v = q[0], py = Y(v);
      const tip = `${s.labels[p]}, ${g6(x)} % ${s.dopant}: peak ${i + 1} of ${r.peaks.length} at ${fmt(v, dec)} ± ${fmt(q[1], dec)}${un} (n = ${r.n}${r.d.nflag ? ", flagged " + r.d.nflag : ""}, ${r.d.k} Gaussian${r.d.k > 1 ? "s" : ""}; c axis along ${axes(r)}; repeats ${r.d.reps.join(", ")})`;
      t += marker(mk(p), px, py, main ? 4.4 : 4.0, main ? col(p) : "#fff", main ? "#fff" : col(p), main ? 0.6 : 1.1, tip);
    });
    return t;
  };
  const branches = (Y, pset) => {
    let t = "";
    pset.forEach(p => {
      const xl = comps.filter(x => res[key(x, p)]);
      for (let i = 1; i < xl.length; i++){
        const P0 = res[key(xl[i - 1], p)].peaks.map(q => q[0]), P1 = res[key(xl[i], p)].peaks.map(q => q[0]);
        const X0 = X(xl[i - 1] * xs + shift[p]), X1 = X(xl[i] * xs + shift[p]);
        link(P0, P1).forEach(([a, b]) => { t += line(X0, Y(P0[a]), X1, Y(P1[b]), col(p), 1.2, DASH); });
      }
    });
    return t;
  };
  let nMulti = 0, nClip = 0; R.forEach(r => { nMulti += r.peaks.length > 1; nClip += r.nclip; });
  const letters = "abcdefghijklmnop";
  panels.forEach((pn, idx) => {
    const last = idx === panels.length - 1, id = clip(pn);
    let body = "";
    const pset = rows ? [pn.p] : phases, own = R.filter(r => !rows || r.d.p === pn.p);
    const yl = s.ylo !== null ? [s.ylo, s.yhi] : yr((s.same_y || !rows) ? R : own);
    const fr = frame(pn, yl, last), Y = fr.Y;
    if (Q.ref !== null && !s.no_ref_line && Q.ref >= yl[0] && Q.ref <= yl[1]) body += line(ML, Y(Q.ref), ML + PW, Y(Q.ref), INK2, 0.9, "1.4 2.4");
    pset.forEach(p => {
      const xl = comps.filter(x => res[key(x, p)]);
      if (band) for (let i = 1; i < xl.length; i++){
        const r0 = res[key(xl[i - 1], p)], r1 = res[key(xl[i], p)], X0 = X(xl[i - 1] * xs + shift[p]), X1 = X(xl[i] * xs + shift[p]);
        link(r0.peaks.map(q => q[0]), r1.peaks.map(q => q[0])).forEach(([a, b]) => {
          body += `<polygon points="${X0.toFixed(2)},${Y(r0.fwhm[a][0]).toFixed(2)} ${X1.toFixed(2)},${Y(r1.fwhm[b][0]).toFixed(2)} ${X1.toFixed(2)},${Y(r1.fwhm[b][1]).toFixed(2)} ${X0.toFixed(2)},${Y(r0.fwhm[a][1]).toFixed(2)}" fill="${col(p)}" fill-opacity="0.28"/>`; });
      }
      if (!band) xl.forEach(x => {
        const r = res[key(x, p)]; if (r.n < 2) return;
        const X0 = x * xs + shift[p], cap = (rows ? 1.15 : 1.0) * slot, px0 = X(X0);
        if (!s.no_bars){
          for (let k = 0; k < r.prob.length;){
            if (r.prob[k] === null){ k++; continue; }
            let j = k, pts = [[px0, Y(r.edges[k])]];
            while (j < r.prob.length && r.prob[j] !== null){ const lx = X(X0 + sgn * Math.min(r.prob[j] * r.sc, cap)); pts.push([lx, Y(r.edges[j])], [lx, Y(r.edges[j + 1])]); j++; }
            pts.push([px0, Y(r.edges[j])]);
            const dstr = pts.map((q, i) => (i ? "L" : "M") + q[0].toFixed(2) + " " + q[1].toFixed(2)).join("");
            body += `<path d="${dstr}Z" fill="${col(p)}" fill-opacity="${rows ? 0.35 : 0.22}" stroke="none"/>`;
            if (!rows) body += `<path d="${dstr}" fill="none" stroke="${col(p)}" stroke-opacity="0.9" stroke-width="0.75"/>`;
            k = j;
          }
        }
        if (!s.no_fit || s.no_bars){
          for (let k = 0; k < r.pdf.length;){
            if (r.pdf[k] === null){ k++; continue; }
            let j = k; const pts = [];
            while (j < r.pdf.length && r.pdf[j] !== null){ pts.push([X(X0 + sgn * Math.min(r.pdf[j] * r.sc, cap)), Y(r.g[j])]); j++; }
            const dstr = pts.map((q, i) => (i ? "L" : "M") + q[0].toFixed(2) + " " + q[1].toFixed(2)).join("");
            if (s.no_bars) body += `<path d="M${px0.toFixed(2)} ${pts[0][1].toFixed(2)}${dstr.replace(/^M/, "L")}L${px0.toFixed(2)} ${pts[pts.length - 1][1].toFixed(2)}Z" fill="${col(p)}" fill-opacity="0.25"/>`;
            body += `<path d="${dstr}" fill="none" stroke="${col(p)}" stroke-width="1.2" stroke-linejoin="round"/>`;
            k = j;
          }
        }
      });
    });
    if (!s.no_lines) body += branches(Y, pset);
    if (s.show_lowest_f) pset.forEach(p => {
      const xl = comps.filter(x => res[key(x, p)] && res[key(x, p)].d.atF !== null);
      const pts = xl.map(x => [X(x * xs + shift[p]), Y(res[key(x, p)].d.atF)]);
      if (pts.length > 1) body += `<polyline points="${pts.map(q => q[0].toFixed(2) + "," + q[1].toFixed(2)).join(" ")}" fill="none" stroke="${col(p)}" stroke-width="1.1" stroke-dasharray="1.4 2.2"/>`;
      xl.forEach((x, i) => { const v = res[key(x, p)].d.atF;
        body += marker("X", pts[i][0], pts[i][1], 3.6, col(p), "#fff", 0.5, `${s.labels[p]}, ${g6(x)} % ${s.dopant}: lowest-F structure, ${Q.name} = ${fmt(v, dec)}${un}`); });
    });
    pset.forEach(p => comps.forEach(x => { const r = res[key(x, p)]; if (r) body += peakMarks(r, x, p, X(x * xs + shift[p]), Y); }));
    if (rows || panels.length > 1){
      const tag = `(${letters[idx]})` + (rows ? "  " + esc(s.labels[pn.p]) : "");
      fr.deco += `<text transform="translate(${ML + PW + 16},${pn.y + pn.h / 2}) rotate(90)" text-anchor="middle" font-size="13.5">${tag}</text>`;
    }
    o += `<g clip-path="url(#${id})">${body}</g>` + fr.deco;
  });
  const a0 = panels[0], a1 = panels[panels.length - 1];
  o += `<text transform="translate(${ML - 62},${(a0.y + a1.y + a1.h) / 2}) rotate(-90)" text-anchor="middle" font-size="15.5">${qlabel(Q)}</text>`;
  // legend
  const items = phases.map(p => ({p: p, text: s.labels[p]})); items.push({p: -1, text: "second peak"});
  if (s.show_lowest_f) items.push({p: -2, text: "lowest F"});
  const wd = items.map(it => 44 + 8 + it.text.length * 7.6 + 22), tot = wd.reduce((a, b) => a + b, 0) - 22;
  let lx = ML + PW / 2 - tot / 2; const ly = 24;
  items.forEach((it, i) => {
    if (it.p >= 0){ o += `<rect x="${lx.toFixed(1)}" y="${ly - 6}" width="44" height="12" fill="${col(it.p)}" fill-opacity="0.35"/>` + line(lx, ly, lx + 44, ly, col(it.p), 1.2, DASH) +
      marker(mk(it.p), lx + 22, ly, 4.4, col(it.p), "#fff", 0.6); }
    else if (it.p === -1) o += marker("o", lx + 22, ly, 4.0, "#fff", INK2, 1.1);
    else o += marker("X", lx + 22, ly, 3.6, INK2, "#fff", 0.5);
    o += `<text x="${(lx + 52).toFixed(1)}" y="${ly + 4.8}" font-size="13.5">${esc(it.text)}</text>`; lx += wd[i];
  });
  // x label + notes
  o += `<text x="${ML + PW / 2}" y="${bottom + 42}" text-anchor="middle" font-size="15.5">${esc(s.dopant)} ${s.fraction ? "fraction of cations" : "concentration (% of cations)"}</text>`;
  const info = [];
  const def = {ca: `c_n/a_n: cell lengths / supercell repeats, c = ${D.c_axis === "unique" ? "most different" : "longest"} axis of each structure, a_n = mean of the other two`,
               theta: "θ: cell angle farthest from 90°, folded to ≥ 90°", a_n: "a_n: mean of the two axes other than c (normalised to one fluorite cell)",
               c_n: `c_n: ${D.c_axis === "unique" ? "most different" : "longest"} normalised axis`, vol: "V/atom: cell volume per atom"}[Q.id] || "";
  info.push((D.T ? `T = ${D.T} K relaxed cells; ` : "") + def);
  info.push(band ? "bands: width at half maximum of every peak of the Gaussian fit; markers: the peaks"
    : `histograms: P(${Q.name}), bin ${g6(bw)}${un}, ${s.norm === "each" ? "each scaled to the same length" : "one scale"}; line: Gaussian fit; markers: its peaks`);
  if (nClip && !band) info.push(`${nClip} outlying values (bins below ${g6(100 * s.tail)} % of the highest bar, or lone groups under 1 %) not drawn`);
  const nf = D.dists.filter(d => phases.includes(d.p)).reduce((a, d) => a + d.nflag, 0);
  if (nf) info.push(`${nf} flagged rows (imaginary modes / not converged) ${D.flagged_in ? "included" : "left out"}`);
  info.forEach((t, i) => { o += `<text x="${ML + PW}" y="${bottom + 64 + 15 * i}" text-anchor="end" font-size="11" fill="${INK2}">${esc(t)}</text>`; });
  const H = bottom + 64 + 15 * info.length;
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="${FF}" fill="${INK}"><rect width="${W}" height="${H}" fill="#fff"/><defs>${defs}</defs>${o}</svg>`;
  return {svg: svg, W: W, H: H, text: `${R.reduce((a, r) => a + r.n, 0)} rows, ${comps.length} compositions, ${phases.length} phases; ${nMulti} of ${R.length} distributions have more than one peak; bin ${g6(bw)}${un}`};
}

function command(s){
  const a = [];
  const f = D.fit, qs = D.quants.map(q => q.id);
  if (qs.length !== 1 || qs[0] !== "ca") a.push("--quantities", ...qs);
  if (D.c_axis !== "longest") a.push("--c-axis", D.c_axis);
  if (D.repeats) a.push("--repeats", ...D.repeats);
  SE.forEach(k => { if (s[k] !== DEF[k]) a.push("--" + k, s[k]); });
  if (s.phases.length !== D.phases.length) a.push("--phases", ...s.phases.map(i => D.phases[i].id));
  if (s.xlo !== null && s.xhi !== null) a.push("--xrange", g6(s.xlo), g6(s.xhi));
  if (s.ylo !== null && s.yhi !== null) a.push("--ylim", g6(s.ylo), g6(s.yhi));
  SL.forEach(k => { if (Math.abs(s[k] - DEF[k]) > 1e-9) a.push("--" + k, g6(s[k])); });
  if (s.bin > 0) a.push("--bin", g6(s.bin));
  CH.forEach(k => { if (s[k]) a.push("--" + k.replace(/_/g, "-")); });
  D.phases.forEach((q, i) => { if (s.labels[i] !== q.deflabel) a.push("--label", `${q.id}="${s.labels[i]}"`); });
  if (s.dopant !== DEF.dopant) a.push("--dopant", s.dopant);
  if (f.max_components !== 3) a.push("--max-components", f.max_components);
  f.components.forEach(c => a.push("--components", c));
  if (f.min_peak !== 0.15) a.push("--min-peak", f.min_peak);
  if (f.bootstrap !== 20) a.push("--bootstrap", f.bootstrap);
  if (f.seed !== 0) a.push("--seed", f.seed);
  if (D.flagged_in) a.push("--include-flagged");
  return "python scripts/plot_ca_vs_x_html.py " + D.csv + (a.length ? " " + a.join(" ") : "");
}

let last = null;
function draw(){
  const s = state();
  try {
    last = build(s);
    $("plotbox").innerHTML = last.svg; $("status").textContent = last.text; $("status").className = "";
    $("cmd").textContent = command(s);
  } catch (e){ $("status").textContent = typeof e === "string" ? e : "could not draw: " + e; $("status").className = "bad"; if (typeof e !== "string") console.error(e); }
}
function download(name, blob){ const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }
const stem = () => D.csv.replace(/\.csv$/i, "") + "_" + state().quantity + "_vs_x";
function start(){
  $("quantity").innerHTML = D.quants.map(q => `<option value="${q.id}">${esc(qname(q))}</option>`).join("");
  $("file").textContent = `${D.csv} · ${D.dists.reduce((a, d) => a + d.n, 0)} rows · ${new Set(D.dists.map(d => d.x)).size} compositions · phases ${D.phases.map(q => q.id).join(", ")}`;
  $("phases").innerHTML = D.phases.map((q, i) => `<label><input type="checkbox" id="ph${i}">${esc(q.id)}</label>`).join("");
  $("labels").innerHTML = D.phases.map((q, i) => `<div class="row"><label for="lb${i}">label of ${esc(q.id)}</label><input type="text" id="lb${i}"></div>`).join("");
  const f = D.fit;
  $("fitnote").textContent = `c/a = c_n/a_n: lattice lengths divided by the supercell repeats (${D.repeats ? D.repeats.join("×") + " for all" : "found per structure"}); c = the ${D.c_axis === "unique" ? "axis most different from the other two" : "longest normalised axis"} of each structure (x, y or z), a_n = mean of the other two. At most ${f.max_components} Gaussians` + (f.components.length ? `, forced: ${f.components.join(" ")}` : "") + `; min peak height ${f.min_peak}; bootstrap ${f.bootstrap}; seed ${f.seed}; flagged rows ${D.flagged_in ? "included" : "left out"}. ${D.nbad ? " " + D.nbad + " rows skipped (no cell data or repeats not found)." : ""} To change these, write the page again: python scripts/plot_ca_vs_x_html.py ${D.csv} --c-axis unique --quantities ca theta ...`;
  setControls(D.init);
  $("controls").addEventListener("change", e => { if (e.target.id === "quantity"){ $("bin").value = ""; $("ylo").value = ""; $("yhi").value = ""; draw(); } });
  $("controls").addEventListener("input", e => { if (e.target.id === "quantity") return; if (e.target.type === "range") e.target.nextElementSibling.textContent = fmt(e.target.value, e.target.id === "tail" ? 3 : 2); draw(); });
  $("reset").onclick = () => { setControls(D.init); draw(); };
  $("svgb").onclick = () => { if (last) download(stem() + ".svg", new Blob([last.svg], {type: "image/svg+xml"})); };
  $("pngb").onclick = () => {
    if (!last) return;
    const img = new Image(), k = 3;
    img.onload = () => { const c = document.createElement("canvas"); c.width = last.W * k; c.height = last.H * k; c.getContext("2d").drawImage(img, 0, 0, c.width, c.height); c.toBlob(b => download(stem() + ".png", b), "image/png"); };
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(last.svg);
  };
  draw();
}
start();
</script></body></html>
"""


if __name__ == "__main__":
    main()
