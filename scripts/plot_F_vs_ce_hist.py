#!/usr/bin/env python3
"""
plot_F_vs_ce_hist.py - free energy vs dopant (Ce) concentration, with the full distribution of F at every
concentration drawn as a horizontal histogram ("violin made of bars"), one colour per starting phase.

    x  = Ce concentration in % of the cation sites  (100 * n_Ce / n_cation)
    y  = F_total (harmonic, fixed relaxed cell, not QHA) of every structure at that concentration

Input: the merged CSV of campaign.py  (free_energy_<T>K_all.csv: columns n_dopant, x_dopant + batch_free_energy.py
columns), or one or more per-composition CSVs of batch_free_energy.py (concentration then taken from the
'composition' formula, e.g. Zr26Ce6O64 -> 6/32 = 18.75 %).

y reference (--ref):
    none          absolute F_total in eV/atom (default)
    global        F - lowest F of all plotted rows                                 [meV/atom]
    composition   F - lowest F of the same concentration (= dF_total of the CSV)   [meV/atom]
    mixing        F - [(1-x) F(0 %) + x F(100 %)], end members = lowest F at Ce 0 % and 100 %
                  (or of --mix-ref-phase), i.e. the free energy of mixing           [meV/atom]
--per fu multiplies everything by 3 (per MO2 formula unit = per cation, as all structures are MO2).
The minimum is taken AFTER the row filters (phases, exclusions), so it is the minimum of what is drawn.

Histogram layout (--layout); with --anchor line (default) every histogram starts at x = its concentration
(thin vertical baseline there) and grows to the right; --anchor center centres the slot on the concentration:
    overlay   all phases on the same baseline, semi-transparent (default)
    dodge     one histogram per phase side by side inside each concentration slot
    stack     one histogram per concentration, bins stacked by phase
    combine   one histogram per concentration, all phases pooled, one colour
Bin width is common to all concentrations (--bin-width in y units, or auto: median Freedman-Diaconis width);
the bin edges are a common grid, so bins line up between concentrations.
Bar length (--pf / --norm pf): P(F) as in plot_free_energy.py: fraction of the structures of that concentration in
each bin (all phases of one concentration sum to 1), the same scale for all concentrations, with a P(F) scale bar
under the plot. Otherwise (--norm): composition = fraction of the structures at that concentration, scaled per concentration
(default: phases at one concentration are comparable, the longest bar fills the slot), phase = every histogram
scaled to its own maximum (shape only), global = raw counts on one scale for the whole plot (populations comparable
across concentrations; end members with few orderings then look tiny or huge).

Gaussian fits (default, --no-fit to hide): every histogram gets the mixture-of-Gaussians fit of
plot_free_energy.py (1..--max-components Gaussians chosen by BIC, maximum likelihood on the F values, peaks =
maxima of the fitted curve, bootstrap for robustness). Solid line = fitted curve on the histogram's scale,
dotted = individual Gaussians (dodge/combine), dashed = peaks; stack draws cumulative curves. Fewer than 10
structures: one Gaussian from mean/std (fit_ok = False in the peak table). Drawn Gaussians are never narrower than
bin/sqrt(2 pi) (curve top <= a one-bin bar holding the same structures); the peak table uses the real fit.
Peaks go to <stem>_F_vs_Ce_peaks.csv (pct, phase, N, n_gaussians, peak, uncertainty, fit_ok).

Flags: rows with imaginary modes (n_imaginary_modes > 0) or an unconverged relaxation (relax_converged = False)
are kept and drawn hatched by default; --exclude-imaginary / --exclude-unconverged / --exclude-flagged drop them.
FAILED rows are always dropped. A thin line per phase through the lowest F of each concentration shows the lower
envelope (--no-min-line hides it); --mark-min adds a short black tick at every lowest F.

Output (next to the CSV or --out-dir):  <stem>_F_vs_Ce.pdf/.png   and   <stem>_F_vs_Ce_summary.csv
(per concentration and phase: N, flagged, min, median, mean, std of the plotted y quantity).

Usage
    python plot_F_vs_ce_hist.py free_energy_300K_all.csv
    python plot_F_vs_ce_hist.py free_energy_300K_all.csv --ref global
    python plot_F_vs_ce_hist.py free_energy_300K_all.csv --ref mixing --phases t m --exclude-flagged
    python plot_F_vs_ce_hist.py free_energy_300K_all.csv --ref composition --layout stack --bin-width 0.5
    python plot_F_vs_ce_hist.py Ce_*/free_energy_300K.csv --ref global --xrange 10 40 --formats pdf
"""
import argparse
import csv
import hashlib
import os
import pickle
import re
import sys
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_free_energy as PFE  # noqa: E402  (same scripts/ folder: style, phase colours)

plt = PFE.plt
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

PHASE_ORDER = list(PFE.PHASE_STYLE)           # parent, t, m, oIII  (fixed order, colour follows the phase)
COMBINED_COLOR = "#52514e"
NORM_TEXT = {"pf": "P(F) per concentration (one scale)",
             "composition": "fraction per concentration (scaled per concentration)",
             "phase": "scaled per histogram", "global": "count (one scale)"}


# --------------------------------------------------------------------------- data
def formula_counts(formula):
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or "")}


def truthy_false(v):
    return str(v).strip().lower() in ("false", "0", "no")


def load(paths, dopant, anion):
    rows = []
    for p in paths:
        with open(p, newline="") as f:
            for r in csv.DictReader(f):
                r["_src"] = str(p)
                rows.append(r)
    if not rows:
        sys.exit("no rows in the CSV file(s)")
    out, n_failed = [], 0
    for r in rows:
        if r.get("status", "ok") != "ok":
            n_failed += 1
            continue
        try:
            F = float(r["F_total_eV_atom"])
        except (KeyError, ValueError):
            n_failed += 1
            continue
        x = None
        if r.get("x_dopant", "") not in ("", None):
            x = float(r["x_dopant"])
        else:
            c = formula_counts(r.get("composition", ""))
            ncat = sum(v for k, v in c.items() if k != anion)
            if ncat == 0:
                sys.exit(f"cannot get the concentration of row {r.get('file')} in {r['_src']}: "
                         f"no x_dopant column and composition {r.get('composition')!r} has no cations")
            x = c.get(dopant, 0) / ncat
        out.append(dict(x=x, pct=round(100 * x, 6), phase=r.get("phase") or "parent", F=F,
                        imag=int(float(r.get("n_imaginary_modes") or 0)) > 0,
                        unconv=truthy_false(r.get("relax_converged", "")),
                        T=r.get("T_K", ""), file=r.get("file", "")))
    return out, n_failed


def apply_reference(data, ref, mix_phase, per_factor):
    """returns (y array, axis label, reference description)"""
    F = np.array([d["F"] for d in data])
    if ref == "none":
        return F * per_factor, None
    if ref == "global":
        i = int(np.argmin(F))
        y = F - F[i]
        desc = f"global min: {data[i]['phase']} at {data[i]['pct']:g} %, F = {F[i]:.6f} eV/atom"
    elif ref == "composition":
        lo = {}
        for d in data:
            lo[d["pct"]] = min(lo.get(d["pct"], d["F"]), d["F"])
        y = F - np.array([lo[d["pct"]] for d in data])
        desc = "lowest F of each concentration"
    else:  # mixing
        xs = np.array([d["x"] for d in data])
        ends = {}
        for e in (0.0, 1.0):
            sel = [d["F"] for d in data if abs(d["x"] - e) < 1e-9 and (mix_phase is None or d["phase"] == mix_phase)]
            if not sel:
                sys.exit(f"--ref mixing needs rows at Ce {100 * e:g} %"
                         + (f" with phase {mix_phase}" if mix_phase else "") + " after the filters")
            ends[e] = min(sel)
        y = F - ((1 - xs) * ends[0.0] + xs * ends[1.0])
        desc = (f"end members{' (' + mix_phase + ')' if mix_phase else ''}: F(0 %) = {ends[0.0]:.6f}, "
                f"F(100 %) = {ends[1.0]:.6f} eV/atom")
    return 1000 * y * per_factor, desc


def bin_width_auto(groups):
    """median Freedman-Diaconis width over all (concentration, phase) groups with >= 5 values"""
    ws = []
    for v in groups:
        if len(v) >= 5:
            q75, q25 = np.percentile(v, [75, 25])
            w = 2 * (q75 - q25) / len(v) ** (1 / 3)
            if w > 0:
                ws.append(w)
    if ws:
        return float(np.median(ws))
    allv = np.concatenate([np.asarray(v) for v in groups])
    span = np.ptp(allv)
    return span / 20 if span > 0 else 1.0


def nice(w):
    """round a bin width to 1, 2, 2.5 or 5 x 10^k"""
    k = np.floor(np.log10(w))
    m = w / 10 ** k
    m = min((1, 2, 2.5, 5, 10), key=lambda c: abs(np.log(c / m)))
    return m * 10 ** k


# --------------------------------------------------------------------------- drawing
def hatched_bar(ax, x0, y0, w, h, color, flagged, hatch_lw, alpha=1.0):
    if w <= 0:
        return
    if flagged:
        ax.add_patch(Rectangle((x0, y0), w, h, facecolor="white", edgecolor=color, hatch="//////",
                               linewidth=hatch_lw, zorder=2))
    else:
        ax.add_patch(Rectangle((x0, y0), w, h, facecolor=color, edgecolor="none", alpha=alpha, zorder=2))


def _fit_one(job):
    """worker: Gaussian-mixture fit of one histogram (top level so multiprocessing can pickle it)"""
    key, y, seed, kmax, min_peak, n_boot = job
    f = PFE.fit_peaks(y, np.random.default_rng(seed), kmax=kmax, min_peak=min_peak, n_boot=n_boot)
    return key, f


def default_jobs():
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    return max(1, min(n, 8))


def main():
    ap = argparse.ArgumentParser(description="F vs Ce concentration with a histogram of F at every concentration")
    ap.add_argument("csv", nargs="+", help="free_energy_<T>K_all.csv from campaign.py merge, or per-composition CSVs")
    ap.add_argument("--ref", choices=["none", "global", "composition", "mixing"], default="none",
                    help="y reference (see header); 'global' = relative to the global minimum")
    ap.add_argument("--relative-to-global-min", action="store_true", help="same as --ref global")
    ap.add_argument("--mix-ref-phase", help="--ref mixing: take the end members from this phase only")
    ap.add_argument("--per", choices=["atom", "fu"], default="atom", help="per atom (default) or per MO2 formula unit")
    ap.add_argument("--layout", choices=["overlay", "dodge", "stack", "combine"], default="overlay")
    ap.add_argument("--anchor", choices=["line", "center"], default="line",
                    help="line: every histogram starts at its concentration (x = Ce %%) and grows to the right "
                         "(default); center: slot centred on the concentration")
    ap.add_argument("--mark-min", action="store_true", help="short black tick at the lowest F of every histogram")
    ap.add_argument("--no-peak-lines", dest="peak_lines", action="store_false",
                    help="no dashed horizontal line at the fitted peaks")
    ap.add_argument("--no-zero-line", dest="zero_line", action="store_false", help="no dashed line at y = 0")
    ap.add_argument("--pf", action="store_true",
                    help="bar length = P(F) as in plot_free_energy.py (fraction of the structures of that "
                         "concentration per bin, all phases sum to 1), one scale for the whole plot + scale bar")
    ap.add_argument("--norm", choices=["composition", "phase", "global", "pf"], default="composition",
                    help="bar length scaling: composition (default), phase (shape only), global (raw counts)")
    ap.add_argument("--bin-width", type=float, help="bin width in y units (meV or eV, as on the axis); default auto")
    ap.add_argument("--bins", type=int, help="instead of --bin-width: this many bins over the whole y range")
    ap.add_argument("--slot", type=float, default=0.85,
                    help="fraction of the smallest concentration step a histogram slot may fill (default 0.85)")
    ap.add_argument("--phases", nargs="+", help="phases to plot, e.g. t m oIII parent (default all in the CSV)")
    ap.add_argument("--exclude-imaginary", action="store_true", help="drop rows with imaginary modes")
    ap.add_argument("--exclude-unconverged", action="store_true", help="drop rows with relax_converged = False")
    ap.add_argument("--exclude-flagged", action="store_true", help="both of the above")
    ap.add_argument("--xrange", nargs=2, type=float, metavar=("MIN", "MAX"), help="Ce %% range to plot")
    ap.add_argument("--ylim", nargs=2, type=float, metavar=("MIN", "MAX"))
    ap.add_argument("--no-min-line", dest="min_line", action="store_false",
                    help="no line through the lowest F of each phase")
    ap.add_argument("--no-fit", dest="fit", action="store_false",
                    help="no Gaussian-mixture curves (default: fit and draw them, as in plot_free_energy.py)")
    ap.add_argument("--max-components", type=int, default=3, help="max Gaussians per histogram (default 3, BIC)")
    ap.add_argument("--min-peak", type=float, default=0.15,
                    help="drop peaks lower than this fraction of the highest one (default 0.15)")
    ap.add_argument("--bootstrap", type=int, default=30,
                    help="bootstrap re-fits per histogram for peak robustness/uncertainty (default 30; 0 = off)")
    ap.add_argument("--no-components", dest="show_components", action="store_false",
                    help="do not draw the individual Gaussians (dotted) of multi-Gaussian fits")
    ap.add_argument("--seed", type=int, default=0, help="random seed for fits/bootstrap")
    ap.add_argument("--jobs", type=int, default=default_jobs(),
                    help="parallel processes for the fits (default: available cores, at most 8)")
    ap.add_argument("--no-cache", dest="cache", action="store_false",
                    help="do not reuse/save fits in <stem>_F_vs_Ce_fitcache.pkl (fits are cached by their data and "
                         "fit options, so re-plotting with other layout/label/bin options is fast)")
    ap.add_argument("--dopant", default="Ce", help="dopant element (default Ce), used if there is no x_dopant column")
    ap.add_argument("--anion", default="O", help="anion element excluded from the cation count (default O)")
    ap.add_argument("--label", action="append", default=[], metavar="PHASE=TEXT", help="legend text for a phase")
    ap.add_argument("--width-cm", type=float, default=17.8, help="figure width (default 17.8 = double column)")
    ap.add_argument("--height-cm", type=float, default=9.0)
    ap.add_argument("--font", choices=["serif", "sans"], default="serif")
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"], choices=["pdf", "png", "svg", "eps"])
    ap.add_argument("--out-dir", help="output folder (default: next to the first CSV)")
    ap.add_argument("--stem", help="output file stem (default: stem of the first CSV)")
    a = ap.parse_args()
    if a.relative_to_global_min:
        if a.ref not in ("none", "global"):
            sys.exit("--relative-to-global-min conflicts with --ref " + a.ref)
        a.ref = "global"
    if a.pf:
        a.norm = "pf"
    if a.exclude_flagged:
        a.exclude_imaginary = a.exclude_unconverged = True
    if a.bins is not None and a.bin_width is not None:
        sys.exit("give --bins or --bin-width, not both")
    labels = dict(s.split("=", 1) for s in a.label)

    data, n_failed = load(a.csv, a.dopant, a.anion)
    n_read = len(data)
    present = list(OrderedDict.fromkeys(d["phase"] for d in data))
    if a.phases:
        miss = [p for p in a.phases if p not in present]
        if miss:
            sys.exit(f"phase(s) {miss} not in the CSV (has {present})")
        data = [d for d in data if d["phase"] in a.phases]
    if a.xrange:
        data = [d for d in data if a.xrange[0] - 1e-9 <= d["pct"] <= a.xrange[1] + 1e-9]
    n_imag = sum(d["imag"] for d in data)
    n_unc = sum(d["unconv"] for d in data)
    if a.exclude_imaginary:
        data = [d for d in data if not d["imag"]]
    if a.exclude_unconverged:
        data = [d for d in data if not d["unconv"]]
    if not data:
        sys.exit("no rows left after the filters")
    phases = [p for p in PHASE_ORDER if p in {d["phase"] for d in data}]
    phases += sorted({d["phase"] for d in data} - set(phases))
    extra = {}
    styles = {p: PFE.style_for(p, extra, labels) for p in phases}

    per_factor = 3.0 if a.per == "fu" else 1.0
    y, ref_desc = apply_reference(data, a.ref, a.mix_ref_phase, per_factor)
    for d, v in zip(data, y):
        d["y"] = float(v)
    unit = ("meV" if a.ref != "none" else "eV") + ("/f.u." if a.per == "fu" else "/atom")
    ylab = {"none": r"$F$", "global": r"$F - F_{\min}^{\mathrm{global}}$",
            "composition": r"$F - F_{\min}(x)$", "mixing": r"$\Delta F_{\mathrm{mix}}$"}[a.ref] + f" ({unit})"

    pcts = sorted({d["pct"] for d in data})
    groups = {(c, p): [d for d in data if d["pct"] == c and d["phase"] == p] for c in pcts for p in phases}
    groups = {k: v for k, v in groups.items() if v}

    # ---- bins: one global grid
    yall = np.array([d["y"] for d in data])
    if a.bins:
        span = np.ptp(yall) or 1.0
        bw = span / a.bins
    elif a.bin_width:
        bw = a.bin_width
    else:
        bw = nice(bin_width_auto([[d["y"] for d in v] for v in groups.values()]))
    if bw <= 0:
        sys.exit("bin width must be > 0")
    y0 = np.floor(yall.min() / bw) * bw

    def hist(vals):
        idx = np.floor((np.asarray(vals) - y0) / bw + 1e-9).astype(int)
        return idx

    # ---- geometry
    steps = np.diff(pcts)
    step = float(steps.min()) if len(steps) else 10.0
    slot = a.slot * step
    sub = {"overlay": 1, "dodge": len(phases), "stack": 1, "combine": 1}[a.layout]
    sub_w = slot / sub

    # counts per (pct, phase): dict bin -> [clean, flagged]
    counts = {}
    for k, v in groups.items():
        cc = {}
        for d, b in zip(v, hist([d["y"] for d in v])):
            cc.setdefault(int(b), [0, 0])[1 if (d["imag"] or d["unconv"]) else 0] += 1
        counts[k] = cc
    n_at = {c: sum(len(v) for (cc, _), v in groups.items() if cc == c) for c in pcts}

    def denom(c, p):
        if a.norm in ("composition", "pf"):
            return n_at[c]
        if a.norm == "phase":
            return len(groups[(c, p)])
        return 1

    # ---- Gaussian-mixture fits (same method as plot_free_energy.py: 1..K Gaussians by BIC, peaks = maxima)
    fits = {}
    if a.fit:
        fit_groups = ({(c, "__all__"): [d for d in data if d["pct"] == c] for c in pcts} if a.layout == "combine"
                      else groups)
        # cache: key = hash of the fitted values + fit options (independent of layout, bins, labels, job count)
        od0 = Path(a.out_dir) if a.out_dir else Path(a.csv[0]).parent
        cache_p = od0 / f"{a.stem or Path(a.csv[0]).stem}_F_vs_{a.dopant}_fitcache.pkl"
        cache = {}
        if a.cache and cache_p.exists():
            try:
                cache = pickle.loads(cache_p.read_bytes())
            except Exception:
                cache = {}
        opts = f"{a.max_components}|{a.min_peak}|{a.bootstrap}|{a.seed}".encode()
        hkey, jobs = {}, []
        for k, v in fit_groups.items():
            y_ = np.array([d["y"] for d in v], float)
            h = hashlib.sha1(np.ascontiguousarray(np.sort(y_)).tobytes() + opts).hexdigest()
            hkey[k] = h
            if h not in cache:                     # seed from the data hash: same result for any order / --jobs
                jobs.append((h, y_, (int(h[:12], 16) + a.seed) % 2**63, a.max_components, a.min_peak, a.bootstrap))
        t0 = time.time()
        print(f"Gaussian fits: {len(fit_groups)} histograms, {len(fit_groups) - len(jobs)} from cache, "
              f"{len(jobs)} to fit on {min(a.jobs, max(len(jobs), 1))} process(es)", flush=True)
        done = 0
        def note(r):
            nonlocal done
            done += 1
            cache[r[0]] = r[1]
            if done == len(jobs) or done % max(1, len(jobs) // 10) == 0:
                print(f"  {done}/{len(jobs)} fitted, {time.time() - t0:.0f} s", flush=True)
        if jobs and a.jobs > 1 and len(jobs) > 1:
            import multiprocessing as mp
            with mp.Pool(min(a.jobs, len(jobs))) as pool:
                for r in pool.imap_unordered(_fit_one, jobs):
                    note(r)
        else:
            for j in jobs:
                note(_fit_one(j))
        if a.cache and jobs:
            od0.mkdir(parents=True, exist_ok=True)
            cache_p.write_bytes(pickle.dumps(cache))
        for k, v in fit_groups.items():
            f = dict(cache[hkey[k]])
            w_, mu_, sg_ = PFE.display_params(f, bw)
            sg_ = np.maximum(sg_, bw / np.sqrt(2 * np.pi))   # drawn width >= one bin: no spikes past the bars
            lo_, hi_ = min(d["y"] for d in v), max(d["y"] for d in v)
            span = max(3 * float(np.max(sg_)), bw)
            yy = np.linspace(lo_ - span, hi_ + span, 400)
            f["yy"], f["dens"] = yy, len(v) * bw * PFE.mix_pdf(yy, w_, mu_, sg_)   # expected count per bin
            f["comp"] = [len(v) * bw * PFE.mix_pdf(yy, [wj], [mj], [sj]) for wj, mj, sj in zip(w_, mu_, sg_)]
            fits[k] = f

    def curve(c, p):
        """fitted curve of one unit as bar length (before scale): (y grid, length, components) or None"""
        if a.layout == "combine":
            f = fits.get((c, "__all__"))
            if f is None or p != phases[0]:
                return None
            dn = 1 if a.norm == "global" else n_at[c]
            return f["yy"], f["dens"] / dn, [g / dn for g in f["comp"]]
        f = fits.get((c, p))
        if f is None:
            return None
        if a.layout in ("dodge", "overlay"):
            return f["yy"], f["dens"] / denom(c, p), [g / denom(c, p) for g in f["comp"]]
        # stack: cumulative over the phases up to p, on p's grid
        tot = np.zeros_like(f["yy"])
        for q in phases[:phases.index(p) + 1]:
            fq = fits.get((c, q))
            if fq is not None:
                tot += np.interp(f["yy"], fq["yy"], fq["dens"], left=0, right=0) / denom(c, q)
        return f["yy"], tot, []

    # bar-length scale.  composition: per concentration (longest bar of that concentration fills its unit),
    # phase: per histogram (shape only), global: one scale for the whole plot (raw counts comparable)
    def longest_of(c, p=None):
        if a.layout in ("dodge", "overlay"):
            ps = [p] if p else phases
            L = max((max(sum(v) for v in counts[(c, q)].values()) / denom(c, q)
                     for q in ps if (c, q) in counts), default=0.0)
        else:
            tot = {}
            for q in phases:
                for b, v in counts.get((c, q), {}).items():
                    tot[b] = tot.get(b, 0) + sum(v) / denom(c, q)
            L = max(tot.values(), default=0.0)
            ps = phases
        for q in ps:                                   # the fitted curve must fit in the slot too
            cv = curve(c, q)
            if cv is not None:
                L = max(L, float(cv[1].max()))
        return L

    glob_long = max(longest_of(c) for c in pcts)
    def scale_of(c, p):
        if a.norm in ("global", "pf"):
            L = glob_long
        elif a.norm == "phase" and a.layout in ("dodge", "overlay"):
            L = longest_of(c, p)
        else:
            L = longest_of(c)
        return sub_w * 0.95 / L if L > 0 else 0.0

    PFE.set_style(a.font)
    fig, ax = plt.subplots(figsize=(a.width_cm * PFE.CM, a.height_cm * PFE.CM))
    gap = 0.12 * bw                                # thin surface gap between neighbouring bins
    for c in pcts:
        left = c if a.anchor == "line" else c - slot / 2
        if a.anchor == "line":                         # the baseline the histograms grow from
            yc = [d["y"] for d in data if d["pct"] == c]
            ax.vlines(c, min(yc) - bw, max(yc) + bw, color=PFE.INK2, lw=0.5, zorder=1)
        for i, p in enumerate(phases):
            cc = counts.get((c, p))
            if cc is None:
                continue
            col = COMBINED_COLOR if a.layout == "combine" else styles[p]["color"]
            for b, (nc, nf) in cc.items():
                yb = y0 + b * bw + gap / 2
                if a.layout == "dodge":
                    x0 = left + i * sub_w
                elif a.layout == "overlay":
                    x0 = left
                else:                                   # stack/combine: start after the phases drawn before
                    x0 = left + sum(sum(counts.get((c, q), {}).get(b, [0, 0])) / denom(c, q)
                                    for q in phases[:i]) * scale_of(c, p)
                sc = scale_of(c, p)
                wc, wf = nc / denom(c, p) * sc, nf / denom(c, p) * sc
                hatched_bar(ax, x0, yb, wc, bw - gap, col, False, 0.4,
                            alpha=(0.3 if a.fit else 0.5) if a.layout == "overlay" else (0.35 if a.fit else 1.0))
                hatched_bar(ax, x0 + wc, yb, wf, bw - gap, col, True, 0.4)
    # fitted curves (solid), components (dotted, k > 1, dodge/combine), dashed line at every peak
    peak_rows = []
    for c in pcts:
        left = c if a.anchor == "line" else c - slot / 2
        for i, p in enumerate(phases):
            cv = curve(c, p)
            if cv is None:
                continue
            yy, ln, comps = cv
            sc = scale_of(c, p)
            base_x = left + (i * sub_w if a.layout == "dodge" else 0.0)
            col = COMBINED_COLOR if a.layout == "combine" else styles[p]["color"]
            keep = ln >= 0.05 * ln.max()                           # no long zero tails along the baseline
            if keep.sum() > 1:
                j0, j1 = np.nonzero(keep)[0][[0, -1]]
                ax.plot(base_x + ln[j0:j1 + 1] * sc, yy[j0:j1 + 1], color=col, lw=1.0, zorder=5)
            if a.show_components and len(comps) > 1:
                for g in comps:
                    m = g >= 0.05 * ln.max()
                    if m.sum() > 1:
                        j0, j1 = np.nonzero(m)[0][[0, -1]]
                        ax.plot(base_x + g[j0:j1 + 1] * sc, yy[j0:j1 + 1], color=col, lw=0.5,
                                ls=(0, (1, 1.2)), zorder=5)
            f = fits[(c, "__all__") if a.layout == "combine" else (c, p)]
            for F0, err, hgt in f["peaks"]:
                xl = float(np.interp(F0, yy, ln)) * sc
                if a.peak_lines:
                    ax.plot([base_x, base_x + xl], [F0, F0], color=col, lw=0.6, ls=(0, (3, 2)), zorder=5)
                peak_rows.append((c, "all" if a.layout == "combine" else p, f["n"], f["k"], F0, err,
                                  f["ok"]))

    # lowest-F markers + envelope line per phase
    for i, p in enumerate(phases):
        xs, ys = [], []
        for c in pcts:
            v = groups.get((c, p))
            if not v:
                continue
            left = c if a.anchor == "line" else c - slot / 2
            xm = left + (i + 0.5) * sub_w if a.layout == "dodge" else (c if a.anchor == "center" else left)
            xs.append(xm)
            ys.append(min(d["y"] for d in v))
        col = COMBINED_COLOR if a.layout == "combine" else styles[p]["color"]
        half = (sub_w if a.layout == "dodge" else slot) * 0.5
        for xm, ym in zip(xs, ys):
            if not a.mark_min:
                break
            if a.anchor == "line" and a.layout != "dodge":
                xm = xm + half
            ax.plot([xm - half, xm + half], [ym, ym], color=PFE.INK, lw=0.9, zorder=4, solid_capstyle="butt")
        if a.min_line and len(xs) > 1 and a.layout != "combine":
            ax.plot(xs, ys, color=col, lw=1.0, ls=(0, (1, 2)), zorder=3, alpha=0.9)   # dotted envelope
    if a.min_line and a.layout == "combine" and len(pcts) > 1:
        lo = [min(d["y"] for d in data if d["pct"] == c) for c in pcts]
        ax.plot(pcts, lo, color=PFE.INK2, lw=1.0, ls=(0, (1, 2)), zorder=3)
    if a.zero_line and a.ref in ("global", "composition", "mixing"):
        ax.axhline(0, color=PFE.INK2, lw=0.5, ls=(0, (4, 3)), zorder=1)

    if a.anchor == "line":
        ax.set_xlim(pcts[0] - 0.4 * step, pcts[-1] + slot + 0.4 * step)
    else:
        ax.set_xlim(pcts[0] - slot / 2 - 0.4 * step, pcts[-1] + slot / 2 + 0.4 * step)
    if len(pcts) > 1:                              # ticks only inside the concentration range (no -20 %, 120 %)
        tk = MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10]).tick_values(pcts[0], pcts[-1])
        ax.set_xticks([t for t in tk if pcts[0] - 1e-9 <= t <= pcts[-1] + 1e-9])
    ymin, ymax = y0, y0 + (max(max(cc) for cc in counts.values()) + 1) * bw
    pad = 0.04 * (ymax - ymin or 1)
    ax.set_ylim(*(a.ylim or (ymin - pad, ymax + pad)))
    ax.set_xlabel(f"{a.dopant} concentration (% of cations)")
    ax.ticklabel_format(axis="y", useOffset=False)       # no '-9.498 +' offset text under the legend
    ax.set_ylabel(ylab)
    PFE.finish_axes(ax)

    # legend above the axes (data fill the whole box at every concentration)
    hs = ([Patch(facecolor=COMBINED_COLOR, label="all phases")] if a.layout == "combine"
          else [Patch(facecolor=styles[p]["color"], label=styles[p]["label"]) for p in phases])
    if any(d["imag"] or d["unconv"] for d in data):
        hs.append(Patch(facecolor="white", edgecolor=PFE.INK2, hatch="//////", linewidth=0.4,
                        label="imaginary modes / not converged"))
    if a.mark_min:
        hs.append(Line2D([], [], color=PFE.INK, lw=0.9, label="lowest F"))
    if a.min_line and len(pcts) > 1:
        hs.append(Line2D([], [], color=PFE.INK2, lw=1.0, ls=(0, (1, 2)), label="lowest F per phase"))
    if a.fit:
        hs.append(Line2D([], [], color=PFE.INK2, lw=1.0, label="Gaussian fit"))
        if a.peak_lines:
            hs.append(Line2D([], [], color=PFE.INK2, lw=0.6, ls=(0, (3, 2)), label="peak"))
    ax.legend(handles=hs, loc="lower left", bbox_to_anchor=(0, 1.01), ncol=len(hs), borderaxespad=0,
              columnspacing=1.2, handlelength=1.4)
    if a.norm in ("pf", "global"):                  # scale bar: how long P(F) = v (or N = v) is on the x axis
        sc = scale_of(pcts[0], phases[0])
        v = nice(0.5 * glob_long)
        x0 = ax.get_xlim()[0] + 0.01 * np.ptp(ax.get_xlim())
        yb = -1.7                                     # own line under the info text
        ax.annotate("", xy=(x0 + v * sc, yb), xytext=(x0, yb), xycoords=("data", ax.xaxis.label),
                    textcoords=("data", ax.xaxis.label), annotation_clip=False,
                    arrowprops=dict(arrowstyle="|-|,widthA=0.25,widthB=0.25", lw=0.8, color=PFE.INK,
                                    shrinkA=0, shrinkB=0))
        ax.annotate((f"$P(F)$ = {v:g}" if a.norm == "pf" else f"N = {v:g}"), xy=(x0 + v * sc, yb),
                    xycoords=("data", ax.xaxis.label), xytext=(4, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=7, color=PFE.INK, annotation_clip=False)
    T = sorted({d["T"] for d in data if d["T"]})
    info = (f"T = {', '.join(T)} K, harmonic, fixed cell" if T else "") + \
           f"   bin {bw:g} {unit}, bar length = {NORM_TEXT[a.norm]}"
    ax.annotate(info.strip(), xy=(1, 0), xycoords=("axes fraction", ax.xaxis.label), xytext=(0, -2),
                textcoords="offset points", ha="right", va="top", fontsize=6.5, color=PFE.INK2)

    first = Path(a.csv[0])
    od = Path(a.out_dir) if a.out_dir else first.parent
    od.mkdir(parents=True, exist_ok=True)
    stem = a.stem or first.stem
    base = od / f"{stem}_F_vs_{a.dopant}"
    files = PFE.save(fig, str(base), a.formats)

    # summary table
    sp = od / f"{stem}_F_vs_{a.dopant}_summary.csv"
    with sp.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pct", "phase", "N", "n_imaginary", "n_unconverged", f"min_{unit}", f"median_{unit}",
                    f"mean_{unit}", f"std_{unit}", "file_of_min"])
        for (c, p), v in sorted(groups.items(), key=lambda kv: (kv[0][0], phases.index(kv[0][1]))):
            ys = np.array([d["y"] for d in v])
            j = int(np.argmin(ys))
            w.writerow([f"{c:g}", p, len(v), sum(d["imag"] for d in v), sum(d["unconv"] for d in v),
                        f"{ys.min():.6g}", f"{np.median(ys):.6g}", f"{ys.mean():.6g}", f"{ys.std():.6g}", v[j]["file"]])

    if a.fit:
        pp = od / f"{stem}_F_vs_{a.dopant}_peaks.csv"
        with pp.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pct", "phase", "N", "n_gaussians", f"peak_{unit}", f"peak_err_{unit}", "fit_ok"])
            for c, p, n, k, F0, err, ok in peak_rows:
                w.writerow([f"{c:g}", p, n, k, f"{F0:.6g}", f"{err:.3g}", ok])
        files.append(str(pp))
        print(f"Gaussian fits: {len(fits)} histograms, {len(peak_rows)} peaks "
              f"({sum(1 for f in fits.values() if not f['ok'])} with < 10 structures: one Gaussian from mean/std)")
    print(f"{n_read} ok rows read ({n_failed} failed/unreadable dropped); {len(data)} plotted, "
          f"{len(pcts)} concentrations, phases {phases}")
    print(f"flagged among selected rows: {n_imag} with imaginary modes, {n_unc} unconverged relaxations "
          + ("(excluded)" if (a.exclude_imaginary and a.exclude_unconverged) else
             "(imaginary excluded)" if a.exclude_imaginary else "(unconverged excluded)" if a.exclude_unconverged
             else "(kept, hatched)"))
    if ref_desc:
        print(f"y reference: {ref_desc}")
    print(f"bin width {bw:g} {unit}; harmonic F at fixed cell (not QHA)")
    print("wrote", ", ".join(files + [str(sp)]))


if __name__ == "__main__":
    main()
