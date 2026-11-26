#!/usr/bin/env python3
"""
plot_F_vs_composition.py - free energy of every phase across the whole composition range, from the merged CSV
written by campaign.py (free_energy_<T>K_all.csv).

  (a) at every composition and for every phase: the distribution of F over all orderings, drawn sideways from
      the composition as a filled histogram P(F) (one bin width for every histogram) with the fitted mixture of
      Gaussians as a thin line on top.  Default layout: all phases in one panel, each histogram outlined in its phase
      colour so overlapping phases stay readable; --layout rows gives one row per phase (no overlap).
      Peaks = maxima of the fit (same method as plot_free_energy.py: 1..3 Gaussians by BIC, bootstrap-checked),
      so a phase can have two or more peaks.  The peaks of a phase are joined across compositions into branches:
      where a phase goes from one peak to two the line splits, where two merge it joins.  Main (highest) peak
      filled marker, other peaks hollow.
  (b) every peak minus the lowest peak (any phase) at the same composition, all phases in one panel, with
      bootstrap error bars: which phase / which branch is most stable where.

y axis (--ref):
  formation    F_form = F - [(1-x) F(x=0) + x F(x=1)], the formation (mixing) free energy per atom from the two
               pure end members (x = Ce fraction of the cations; MO2 at every x, so per atom is per 1/3 formula
               unit and the reference is linear in x).  End members (--ends, one value or X0,X1):
                 same    (default) each phase relative to its own pure end members (phase-consistent)
                 lowest  the lowest-F phase of the CSV at each end (e.g. m-ZrO2 and c-CeO2), the same for every
                         phase: the thermodynamic formation free energy
                 PHASE   that phase, e.g. --ends m,parent
               Needs x = 0 and x = 1 rows (status ok, not flagged unless --include-flagged); stops with a clear
               message if one is missing.  Panel (b) and dF in the peaks CSV stay the plain F differences.
  mixing       F - [(1-x) F_A + x F_B]   F_A, F_B = lowest F of the two end compositions (any phase).
               Removes the large linear trend between the end members, so phase differences become visible.
  composition  (default) F - lowest F of the same composition (any phase)
  global       F - lowest F of the whole data set
Rows with imaginary modes or an unconverged relaxation are left out (counted on the figure and in the peaks CSV);
--include-flagged keeps them.  Harmonic, fixed cell: not QHA.  Needs plot_free_energy.py in the same folder.

Look (--style):
  hist    (default) sideways histogram + fitted curve at every composition
  band    no histograms: every peak branch of a phase is a shaded ribbon across the compositions, as wide as the
          peak at half maximum.  Easiest to read when the distributions are narrow compared with the y range.
When the phases differ by much less than the y range, use --ref composition (removes the common curve) and/or
--xrange.  Bins lower than --tail (2 %) of the highest bar, and lone groups of bins with < 1 % of the orderings,
are not drawn, so a few outlying orderings do not stretch the y axis (--tail 0 draws everything).

Histogram length (--norm):
  each    (default) every histogram scaled so its highest bar fills the slot: shapes and peaks of all
          distributions equally readable (P(F) values are not comparable between compositions)
  common  P(F) on one scale for all histograms (fraction of that phase's orderings per bin); the 10 % spikiest
          histograms (few orderings, e.g. near the end members) are clipped so they do not shrink the rest

Output (next to the CSV, or --out-dir):
    <stem>_F_vs_x.pdf/.png     (a) + (b), double column (17.8 cm); --no-panel-b for (a) only
    <stem>_F_vs_x_peaks.csv    every peak: composition, phase, branch, F (chosen ref), error, height, dF to lowest
                               peak, n, n_flagged, number of Gaussians
    <stem>_F_vs_x_fitcache.pkl fits, reused when only the look changes (--no-cache to switch off)

Usage
    python plot_F_vs_composition.py campaign/free_energy_300K_all.csv
    python plot_F_vs_composition.py all.csv --xrange 0 50 --phases t m oIII
    python plot_F_vs_composition.py all.csv --dodge 0.2          # phases side by side at each x
    python plot_F_vs_composition.py all.csv --style band
    python plot_F_vs_composition.py all.csv --layout rows        # one row per phase
    python plot_F_vs_composition.py all.csv --ref mixing --bin 0.5 --norm common
    python plot_F_vs_composition.py all.csv --components m=2 --label parent="c (initial)" --font sans --jobs 8
"""
import argparse
import csv
import hashlib
import os
import pickle
import sys
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402
from matplotlib.patches import Patch                 # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_free_energy import fit_peaks, mix_pdf, display_params   # noqa: E402

PHASE_STYLE = OrderedDict([
    ("parent", dict(color="#52514e", marker="o", label="c (initial)")),
    ("t",      dict(color="#2a78d6", marker="s", label="t")),
    ("m",      dict(color="#eb6834", marker="^", label="m")),
    ("oIII",   dict(color="#1baf7a", marker="D", label="oIII")),
])
EXTRA_STYLES = [dict(color="#6250d6", marker="v"), dict(color="#e34948", marker="P"),
                dict(color="#eda100", marker="X"), dict(color="#e87ba4", marker="h")]
INK, INK2, GRID = "#1f1f1f", "#6b6b6b", "#e6e6e6"


def flagged(r):
    imag = int(float(r.get("n_imaginary_modes") or 0)) > 0
    return imag or r.get("relax_converged") == "False"


def load(path, include_flagged):
    rows = list(csv.DictReader(open(path)))
    if not rows or "F_total_eV_atom" not in rows[0]:
        sys.exit(f"{path}: not a batch_free_energy / campaign CSV")
    if "x_dopant" not in rows[0]:
        sys.exit(f"{path}: no x_dopant column - use the merged free_energy_<T>K_all.csv from campaign.py merge")
    data, nflag = [], {}
    for r in rows:
        if r.get("status") != "ok" or r["x_dopant"] == "" or r["F_total_eV_atom"] == "":
            continue
        x, ph = float(r["x_dopant"]), r["phase"] or "parent"
        if flagged(r):
            nflag[(x, ph)] = nflag.get((x, ph), 0) + 1
            if not include_flagged:
                continue
        data.append((x, ph, float(r["F_total_eV_atom"])))
    if not data:
        hint = f" ({sum(nflag.values())} flagged rows left out: try --include-flagged)" if nflag else ""
        sys.exit(f"{path}: no usable rows{hint}")
    T = next((r["T_K"] for r in rows if r.get("T_K")), "")
    return data, nflag, T


def parse_ends(text):
    parts = [t.strip() for t in text.split(",")]
    if len(parts) == 1:
        parts *= 2
    if len(parts) != 2 or not all(parts):
        sys.exit(f"--ends needs same, lowest or a phase, or two of them as X0,X1 (got {text!r})")
    return parts


def formation_reference(data, nflag, ends, phases, include_flagged, name=str):
    """F_ref(x, phase) = (1-x) F0 + x F1 (eV/atom).  data: (x, phase, F) of every phase; ends: [choice x=0, x=1]."""
    have = sorted({p for x, p, F in data} | {p for x, p in nflag})
    for c in ends:
        if c not in ("same", "lowest") and c not in have:
            sys.exit(f"--ends: phase {c!r} not in the CSV (have {have}; or use same / lowest)")
    comps = sorted({x for x, p, F in data})
    def at(x0, p=None):
        return [F for x, q, F in data if abs(x - x0) < 1e-9 and (p is None or q == p)]
    def missing(x0, p):
        nf = sum(v for (x, q), v in nflag.items() if abs(x - x0) < 1e-9 and (p is None or q == p))
        what = f"phase {p}" if p else "any phase"
        hint = f" ({nf} flagged rows left out: --include-flagged keeps them)" if nf and not include_flagged else ""
        alt = ("choose other end members (--ends lowest or --ends PHASE0,PHASE1)" if ends != ["lowest", "lowest"]
               else "use --ref mixing")
        span = f"; compositions in the CSV: {100 * comps[0]:g} .. {100 * comps[-1]:g} %" if comps else ""
        sys.exit(f"--ref formation: no usable x = {x0:g} row for {what}{hint}{span}.  Pure end members are needed; "
                 f"or {alt}.")
    def end(x0, choice, p):
        if choice == "lowest":
            v = at(x0)
            if not v:
                missing(x0, None)
            best = min(v)
            return best, next(q for x, q, F in data if abs(x - x0) < 1e-9 and F == best)
        q = p if choice == "same" else choice
        v = at(x0, q)
        if not v:
            missing(x0, q)
        return min(v), q
    ref = {}
    for p in phases:
        (f0, q0), (f1, q1) = end(0.0, ends[0], p), end(1.0, ends[1], p)
        ref[p] = (f0, f1, q0, q1)
    def word(i, choice):
        if choice == "same":
            return "same phase"
        if choice == "lowest":
            return "lowest F = " + "/".join(sorted({name(r[2 + i]) for r in ref.values()}))
        return name(choice)
    note = ("F_form = F - [(1-x) F(x=0) + x F(x=1)]; end members: "
            + (f"same phase at both ends" if ends == ["same", "same"]
               else f"x = 0: {word(0, ends[0])}, x = 1: {word(1, ends[1])}"))
    detail = "; ".join(f"{p}: F(0) = {r[0]:.6f} ({r[2]}), F(1) = {r[1]:.6f} ({r[3]}) eV/atom" for p, r in ref.items())
    return (lambda x, p: ref[p][0] + x * (ref[p][1] - ref[p][0])), note, detail


def reference(data, ref):
    xs = np.array([d[0] for d in data]); F = np.array([d[2] for d in data])
    if ref == "global":
        return lambda x: F.min(), f"global minimum, F = {F.min():.6f} eV/atom"
    if ref == "composition":
        low = {x: F[xs == x].min() for x in np.unique(xs)}
        return lambda x: low[x], "lowest F of each composition"
    xa, xb = xs.min(), xs.max()
    if xa == xb:
        sys.exit("--ref mixing needs at least two compositions (use --ref composition)")
    fa, fb = F[xs == xa].min(), F[xs == xb].min()
    note = "" if (xa == 0 and xb == 1) else f" (end compositions present: {100 * xa:g} % and {100 * xb:g} %)"
    return (lambda x: fa + (x - xa) / (xb - xa) * (fb - fa),
            f"line between F({100 * xa:g} %) = {fa:.6f} and F({100 * xb:g} %) = {fb:.6f} eV/atom{note}")


def common_bin(groups):
    """one bin width for all histograms: median Freedman-Diaconis width, rounded to 1/2/2.5/5 x 10^k"""
    ws = []
    for y in groups:
        if len(y) >= 5:
            iqr = np.subtract(*np.percentile(y, [75, 25]))
            if iqr > 0:
                ws.append(2 * iqr * len(y) ** (-1 / 3))
    w = float(np.median(ws)) if ws else 0.5
    k = np.floor(np.log10(w)); m = w / 10 ** k
    return float(min((1, 2, 2.5, 5, 10), key=lambda c: abs(np.log(c / m))) * 10 ** k)


def link(prev, cur):
    """segments joining the peaks of a phase at two neighbouring compositions: every peak is joined to the
    nearest peak on the other side, so 1 -> 2 peaks branches and 2 -> 1 merges"""
    if not prev or not cur:
        return []
    seg = set()
    for j, c in enumerate(cur):
        seg.add((int(np.argmin([abs(c - p) for p in prev])), j))
    for i, p in enumerate(prev):
        seg.add((i, int(np.argmin([abs(p - c) for c in cur]))))
    return sorted(seg)


def style_for(phases, labels):
    out, k = {}, 0
    for p in phases:
        if p in PHASE_STYLE:
            st = dict(PHASE_STYLE[p])
        else:
            st = dict(EXTRA_STYLES[k % len(EXTRA_STYLES)], label=p); k += 1
        if p in labels:
            st["label"] = labels[p]
        out[p] = st
    return out


def _half_max(g, pdf, F0):
    """F range where the curve stays above half of the peak at F0 (stops at a valley between two peaks)."""
    k = int(np.argmin(np.abs(g - F0))); half = 0.5 * pdf[k]; i = j = k
    while i > 0 and pdf[i - 1] >= half and pdf[i - 1] <= pdf[i] + 1e-12:
        i -= 1
    while j < len(g) - 1 and pdf[j + 1] >= half and pdf[j + 1] <= pdf[j] + 1e-12:
        j += 1
    return float(g[i]), float(g[j])


def _fit_job(job):
    """worker (top level for multiprocessing): one Gaussian-mixture fit"""
    key, y, seed, kmax, kfix, min_peak, n_boot = job
    return key, fit_peaks(y, np.random.default_rng(seed), kmax=kmax, k_fixed=kfix, min_peak=min_peak, n_boot=n_boot)


def default_jobs():
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    return max(1, min(n, 8))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("csv")
    ap.add_argument("--ref", choices=["composition", "formation", "mixing", "global"], default="composition",
                    help="y axis: composition = F - lowest F at that x (default); formation = formation free energy "
                         "from the pure end members (see --ends); mixing = F - line between the lowest F of the end "
                         "compositions; global = F - lowest F overall")
    ap.add_argument("--ends", default="same", metavar="X0[,X1]",
                    help="--ref formation: end members at x = 0 and x = 1: same (each phase its own, default), "
                         "lowest (lowest-F phase of the CSV at that end) or a phase name; one value or two as X0,X1")
    ap.add_argument("--phases", nargs="+", help="phases to show, in this order (default: all in the CSV)")
    ap.add_argument("--xrange", type=float, nargs=2, metavar=("LO", "HI"), help="composition range, %% of cations")
    ap.add_argument("--label", action="append", default=[], metavar="PHASE=TEXT", help="legend text of a phase")
    ap.add_argument("--layout", choices=["overlay", "rows"], default="overlay",
                    help="overlay: all phases in one panel (default); rows: one panel per phase, no overlap")
    ap.add_argument("--style", choices=["hist", "band"], default="hist",
                    help="hist: sideways histogram at every composition (default); band: one shaded band per peak "
                         "branch across the compositions (band = width of the peak at half maximum)")
    ap.add_argument("--tail", type=float, default=0.02, metavar="FRAC",
                    help="bins and curve parts lower than FRAC of the highest bar of a distribution are not drawn "
                         "(a few outlying orderings otherwise stretch the y axis); default 0.02, 0 draws everything."
                         "  Fits and the peaks CSV always use all rows")
    ap.add_argument("--bin", type=float, default=0.0,
                    help="histogram bin width, meV/atom, the same for every histogram (default: automatic)")
    ap.add_argument("--width", type=float, default=0.85,
                    help="longest histogram bar as a fraction of the composition step (default 0.85)")
    ap.add_argument("--norm", choices=["each", "common"], default="each",
                    help="each: every histogram scaled to the same length (default); common: one P(F) scale")
    ap.add_argument("--side", choices=["left", "right"], default="left", help="histograms drawn to this side")
    ap.add_argument("--dodge", type=float, default=0.0, metavar="FRAC",
                    help="overlay layout: shift the phases apart by FRAC of the composition step (e.g. 0.2)")
    ap.add_argument("--no-fit", dest="show_fit", action="store_false", help="no fitted curves (peaks still fitted)")
    ap.add_argument("--no-bars", dest="show_bars", action="store_false", help="fitted curves only, no histogram")
    ap.add_argument("--max-components", type=int, default=3, help="most Gaussians per distribution (default 3)")
    ap.add_argument("--components", action="append", default=[], metavar="[PHASE=]K",
                    help="force K Gaussians (all phases, or one phase)")
    ap.add_argument("--min-peak", type=float, default=0.15,
                    help="a peak must reach this fraction of the highest one (default 0.15)")
    ap.add_argument("--bootstrap", type=int, default=20,
                    help="bootstrap re-fits per distribution, for peak errors; peaks that do not reappear in 70 %% of "
                         "them are dropped (default 20; 0 = fast preview without that check)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=default_jobs(),
                    help="parallel processes for the fits (default: available cores, at most 8)")
    ap.add_argument("--no-cache", dest="cache", action="store_false", help="do not reuse / save the fits")
    ap.add_argument("--include-flagged", action="store_true",
                    help="keep rows with imaginary modes / unconverged relaxation")
    ap.add_argument("--no-lines", dest="lines", action="store_false",
                    help="panel (a): no dashed lines joining the peaks across compositions (markers stay)")
    ap.add_argument("--show-min", action="store_true", help="also join the lowest F of every phase (dotted)")
    ap.add_argument("--no-panel-b", dest="panel_b", action="store_false", help="only panel (a)")
    ap.add_argument("--fraction", action="store_true", help="x axis as fraction instead of %% of cations")
    ap.add_argument("--dopant", default="Ce", help="name on the x axis (default Ce)")
    ap.add_argument("--same-y", action="store_true",
                    help="rows layout: the same y range in every phase row (default: each row zoomed to its phase)")
    ap.add_argument("--ylim", type=float, nargs=2, metavar=("LO", "HI"), help="y range of panel (a), meV/atom")
    ap.add_argument("--ylim-b", type=float, nargs=2, metavar=("LO", "HI"), help="y range of panel (b), meV/atom")
    ap.add_argument("--font", choices=["serif", "sans"], default="serif")
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"], choices=["pdf", "png", "svg", "eps"])
    ap.add_argument("--out-dir", help="where to write (default: next to the CSV)")
    a = ap.parse_args()

    labels = {}
    for item in a.label:
        if "=" not in item:
            sys.exit(f"--label needs PHASE=TEXT, got {item!r}")
        k, v = item.split("=", 1); labels[k] = v
    k_all, k_phase = None, {}
    for item in a.components:
        if "=" in item:
            p, v = item.split("=", 1); k_phase[p] = int(v)
        else:
            k_all = int(item)

    data, nflag, T = load(a.csv, a.include_flagged)
    present = list(OrderedDict.fromkeys(p for p in list(PHASE_STYLE) + sorted({d[1] for d in data})
                                        if any(d[1] == p for d in data)))
    phases = a.phases or present
    missing = [p for p in phases if p not in present]
    if missing:
        sys.exit(f"phases {missing} not in the CSV (have {present})")
    ends = parse_ends(a.ends)
    if a.ref == "formation":     # end members from every phase of the CSV (lowest / named phase need not be shown)
        refF, ref_note, ref_detail = formation_reference(data, nflag, ends, phases, a.include_flagged,
                                                         lambda q: labels.get(q, PHASE_STYLE.get(q, {}).get("label", q)))
    data = [d for d in data if d[1] in phases]
    if a.ref != "formation":
        f1, ref_note = reference(data, a.ref)     # reference from all compositions, before --xrange
        refF, ref_detail = (lambda x, p: f1(x)), ""
    if a.xrange:
        data = [d for d in data if a.xrange[0] - 1e-9 <= 100 * d[0] <= a.xrange[1] + 1e-9]
        if not data:
            sys.exit("no rows inside --xrange")
    phases = [p for p in phases if any(d[1] == p for d in data)]
    st = style_for(phases, labels)
    xscale = 1.0 if a.fraction else 100.0

    # ------------------------------------------------------------ distributions
    comps = sorted({d[0] for d in data})
    dx = np.diff(comps).min() * xscale if len(comps) > 1 else 1.0
    res = OrderedDict()
    for p in phases:
        for x in comps:
            y = np.array([1000 * (F - refF(x, p)) for xx, pp, F in data if xx == x and pp == p])
            if len(y):      # base: the reference, meV/atom; panel (b) uses F itself, not F - reference
                res[(x, p)] = dict(y=y, low=float(y.min()), med=float(np.median(y)), base=1000 * refF(x, p))
    if not res:
        sys.exit("nothing to plot")
    bw = a.bin if a.bin > 0 else common_bin([r["y"] for r in res.values()])
    for r in res.values():
        y = r["y"]
        lo = (np.floor(y.min() / bw + 0.5) - 0.5) * bw      # bins centred on multiples of bw (0 is a bin centre)
        hi = (np.floor(y.max() / bw + 0.5) + 0.5) * bw
        r["edges"] = np.arange(lo, hi + 0.5 * bw, bw)
        r["prob"] = np.histogram(y, r["edges"])[0] / len(y)

    # ------------------------------------------------------------ peak fits (parallel, cached)
    out = Path(a.out_dir) if a.out_dir else Path(a.csv).parent
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(a.csv).stem + "_F_vs_x"
    cache_p = out / f"{stem}_fitcache.pkl"
    cache = {}
    if a.cache and cache_p.exists():
        try:
            cache = pickle.loads(cache_p.read_bytes())
        except Exception:
            cache = {}
    jobs, hkey = [], {}
    for (x, p), r in res.items():
        kfix = k_phase.get(p, k_all)
        opts = f"{a.max_components}|{kfix}|{a.min_peak}|{a.bootstrap}|{a.seed}".encode()
        h = hashlib.sha1(np.ascontiguousarray(np.sort(r["y"])).tobytes() + opts).hexdigest()
        hkey[(x, p)] = h
        if h not in cache:
            jobs.append((h, r["y"], (int(h[:12], 16) + a.seed) % 2 ** 63, a.max_components, kfix, a.min_peak,
                         a.bootstrap))
    t0 = time.time()
    print(f"peak fits: {len(res)} distributions, {len(res) - len(jobs)} from cache, {len(jobs)} to fit"
          + (f" on {min(a.jobs, len(jobs))} process(es)" if jobs else ""), flush=True)
    done = 0

    def got(rr):
        nonlocal done
        cache[rr[0]] = rr[1]
        done += 1
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
    for key, r in res.items():
        f = cache[hkey[key]]
        peaks = sorted(f["peaks"], key=lambda q: q[0])                  # (F, err, height), ascending F
        r.update(fit=f, peaks=peaks, main=int(np.argmax([q[2] for q in peaks])))
        # drawn curve: fitted density, never narrower than one bin (no spikes), as probability per bin
        w_, mu_, sg_ = display_params(f, bw)
        sg_ = np.maximum(np.asarray(sg_, float), bw / np.sqrt(2 * np.pi))
        span = 3 * float(sg_.max())
        g = np.linspace(r["edges"][0] - span, r["edges"][-1] + span, 300)
        pdf = mix_pdf(g, w_, mu_, sg_) * bw
        keep = pdf >= 0.02 * pdf.max()
        j0, j1 = np.nonzero(keep)[0][[0, -1]]
        r["g"], r["pdf"] = g[j0:j1 + 1], pdf[j0:j1 + 1]
        r["fwhm"] = [_half_max(g, pdf, q[0]) for q in peaks]            # (low, high) F at half maximum
        # shown: only bins / curve parts above --tail of the highest bar, so a few outlying orderings do not
        # draw long thin lines or stretch the y axis (fit and peaks CSV keep every row)
        if len(r["y"]) > 1 and a.tail > 0:
            thr = a.tail * r["prob"].max()
            off = r["prob"] < thr
            k, nb = 0, len(off)                  # a lone group of bins holding < 1 % of the orderings: also off
            while k < nb:
                if off[k]:
                    k += 1; continue
                j = k
                while j < nb and not off[j]:
                    j += 1
                if r["prob"][k:j].sum() < 0.01:
                    off[k:j] = True
                k = j
            if off.all():
                off[:] = r["prob"] < thr
            r["nclip"] = int(round(float(r["prob"][off].sum()) * len(r["y"])))
            r["prob"] = np.where(off, np.nan, r["prob"])
            on_ = np.nonzero(~off)[0]; e = r["edges"]
            shown = np.zeros(len(r["g"]), bool)
            for b in on_:
                shown |= (r["g"] >= e[b] - bw) & (r["g"] <= e[b + 1] + bw)
            r["pdf"] = np.where((r["pdf"] < a.tail * r["pdf"].max()) | ~shown, np.nan, r["pdf"])
            if np.isnan(r["pdf"]).all():
                r["pdf"] = pdf[j0:j1 + 1]
        else:
            r["nclip"] = 0
        on = np.nonzero(~np.isnan(r["prob"]))[0]
        gon = r["g"][~np.isnan(r["pdf"])]
        r["span"] = (min(r["edges"][on[0]], gon[0]), max(r["edges"][on[-1] + 1], gon[-1]))

    # ------------------------------------------------------------ scales
    slot = a.width * dx
    tops = [max(np.nanmax(r["prob"]), np.nanmax(r["pdf"])) for r in res.values()]
    if a.norm == "each":
        for r in res.values():
            r["sc"] = slot / max(np.nanmax(r["prob"]), np.nanmax(r["pdf"]))
    else:     # one scale; the 10 % spikiest (few orderings) must not shrink all others: they are clipped
        pref = float(np.percentile(tops, 90))
        for r in res.values():
            r["sc"] = slot / pref
    sgn = -1 if a.side == "left" else 1
    rows = a.layout == "rows"
    shift = {p: 0.0 if rows else (i - (len(phases) - 1) / 2) * a.dodge * dx for i, p in enumerate(phases)}
    # panel (b) / dF: plain F differences at the same x (reference added back, so --ref formation with a
    # phase-dependent reference does not change them); stored relative to the lowest peak
    lowest_peak = {x: min(q[0] + r["base"] for (xx, p), r in res.items() if xx == x for q in r["peaks"])
                   for x in comps}
    for (x, p), r in res.items():
        r["db"] = r["base"] - lowest_peak[x]        # dF of a peak = F0 + db

    # ------------------------------------------------------------ figure
    plt.rcParams.update({
        "font.family": "serif" if a.font == "serif" else "sans-serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix" if a.font == "serif" else "dejavusans",
        "font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
        "ytick.right": True, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.edgecolor": INK,
        "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
    W = 17.8 / 2.54
    na = len(phases) if rows else 1
    ha = [1.0] * na if rows else [2.2]
    if a.panel_b:
        ha = ha + [1.1 if rows else 1.1]
    H = (2.6 * na / 2.54 + (3.6 / 2.54 if a.panel_b else 0)) if rows else (W * (0.82 if a.panel_b else 0.52))
    fig, axes = plt.subplots(len(ha), 1, figsize=(W, max(H, 2.2)), sharex=True, sharey=False,
                             gridspec_kw=dict(height_ratios=ha, hspace=0.08))
    axes = list(np.atleast_1d(axes))
    axA = axes[:na]
    bx = axes[na] if a.panel_b else None
    ax_of = {p: (axA[i] if rows else axA[0]) for i, p in enumerate(phases)}

    def draw_dist(ax, X, r, col):
        if len(r["y"]) < 2:                      # a single ordering: marker only, no distribution
            return
        sc = r["sc"]; cap = (1.15 if rows else 1.0) * slot
        if a.show_bars:
            ln = np.minimum(r["prob"] * sc, cap)
            ax.stairs(X + sgn * ln, r["edges"], baseline=X, orientation="horizontal", fill=True,
                      color=col, alpha=0.22 if not rows else 0.35, lw=0, zorder=2)
            if not rows:                         # outline keeps overlapping phases apart in one panel
                ax.stairs(X + sgn * ln, r["edges"], baseline=X, orientation="horizontal", fill=False,
                          color=col, alpha=0.9, lw=0.5, zorder=3)
        if a.show_fit or not a.show_bars:
            ax.plot(X + sgn * np.minimum(r["pdf"] * sc, cap), r["g"], color=col, lw=0.8, zorder=4)
        if not a.show_bars:
            ax.fill_betweenx(r["g"], X, X + sgn * np.minimum(r["pdf"] * sc, cap), color=col, alpha=0.25, lw=0,
                             zorder=2)

    band = a.style == "band"
    if band:       # one ribbon per peak branch: full width at half maximum of the peak, joined like the lines
        for p in phases:
            xs = [x for x in comps if (x, p) in res]
            for x0, x1 in zip(xs[:-1], xs[1:]):
                r0, r1 = res[(x0, p)], res[(x1, p)]
                X0, X1 = x0 * xscale + shift[p], x1 * xscale + shift[p]
                for i, j in link([q[0] for q in r0["peaks"]], [q[0] for q in r1["peaks"]]):
                    (l0, h0), (l1, h1) = r0["fwhm"][i], r1["fwhm"][j]
                    ax_of[p].fill([X0, X1, X1, X0], [l0, l1, h1, h0], color=st[p]["color"], alpha=0.28, lw=0,
                                  zorder=2)
    n_multi = 0
    for (x, p), r in res.items():
        s, X, ax = st[p], x * xscale + shift[p], ax_of[p]
        if not band:
            draw_dist(ax, X, r, s["color"])
        n_multi += len(r["peaks"]) > 1
        for i, (F0, err, h) in enumerate(r["peaks"]):
            main = i == r["main"]
            ax.plot([X], [F0], ls="", marker=s["marker"], ms=3.6 if main else 3.2, zorder=6,
                    mew=0.4 if main else 0.7, mfc=s["color"] if main else "white", mec="white" if main else s["color"])

    def branches(axis_of, value, shifted):
        for p in phases:
            xs = [x for x in comps if (x, p) in res]
            for x0, x1 in zip(xs[:-1], xs[1:]):
                P0 = [value(x0, p, q) for q in res[(x0, p)]["peaks"]]
                P1 = [value(x1, p, q) for q in res[(x1, p)]["peaks"]]
                X0 = x0 * xscale + (shift[p] if shifted else 0)
                X1 = x1 * xscale + (shift[p] if shifted else 0)
                for i, j in link(P0, P1):
                    axis_of(p).plot([X0, X1], [P0[i], P1[j]], color=st[p]["color"], lw=0.8, ls=(0, (3, 2)),
                                    zorder=5)

    if a.lines:
        branches(lambda p: ax_of[p], lambda x, p, q: q[0], True)
    if a.show_min:
        for p in phases:
            pts = [(x * xscale + shift[p], res[(x, p)]["low"]) for x in comps if (x, p) in res]
            ax_of[p].plot(*zip(*pts), color=st[p]["color"], lw=0.8, ls=(0, (1, 1.5)), zorder=4)
    ylab = {"formation": r"$F_\mathrm{form}$ (meV/atom)",
            "mixing": r"$F - F_\mathrm{lin}$ (meV/atom)", "composition": r"$F - F_\mathrm{min}(x)$ (meV/atom)",
            "global": r"$F - F_\mathrm{min}^\mathrm{global}$ (meV/atom)"}[a.ref]
    letters = iter("abcdefghijklmnop")
    def yrange(rr):
        rr = list(rr)
        if band:
            lo_ = min(min(l for l, h in r["fwhm"]) for r in rr)
            hi_ = max(max(h for l, h in r["fwhm"]) for r in rr)
        else:
            lo_ = min(r["span"][0] if len(r["y"]) > 1 else r["low"] for r in rr)
            hi_ = max(r["span"][1] if len(r["y"]) > 1 else r["low"] for r in rr)
        pd_ = 0.06 * (hi_ - lo_ or 1)
        return lo_ - pd_, hi_ + pd_
    for ax in axes:
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)   # no offset text on top of the legend
    for i, ax in enumerate(axA):
        ax.axhline(0, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=1)
        own = [r for (x, p), r in res.items() if not rows or p == phases[i]]
        ax.set_ylim(*(a.ylim or yrange(res.values() if (a.same_y or not rows) else own)))
        # panel letter + phase name outside the right edge: never on top of data
        tag = f"({next(letters)})" + (f"  {st[phases[i]]['label']}" if rows else "")
        ax.text(1.012, 0.5, tag, transform=ax.transAxes, ha="left", va="center", rotation=270, fontsize=9,
                color=INK)
    if not rows:
        axA[0].set_ylabel(ylab)
    handles = [(Patch(facecolor=st[p]["color"], alpha=0.35, lw=0),
                Line2D([], [], color=st[p]["color"], marker=st[p]["marker"], ms=3.6, mec="white", mew=0.4,
                       lw=0.8, ls=(0, (3, 2)))) for p in phases]
    names = [st[p]["label"] for p in phases]
    handles.append(Line2D([], [], ls="", marker="o", ms=3.4, mfc="white", mec=INK2, mew=0.7))
    names.append("second peak")
    if a.show_min:
        handles.append(Line2D([], [], color=INK2, lw=0.8, ls=(0, (1, 1.5)))); names.append("lowest F")
    axA[0].legend(handles, names, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=len(names), frameon=False,
                  handlelength=2.6, columnspacing=1.2, borderaxespad=0)

    nf = sum(v for (x, p), v in nflag.items() if p in phases)
    info = (f"T = {float(T):g} K, harmonic, fixed cell; " if T else "") + \
        (f"histograms: P(F), bin {bw:g} meV/atom, each scaled to the same length" if a.norm == "each" else
         f"histograms: P(F), bin {bw:g} meV/atom, one scale") + "; line: Gaussian fit; markers: its peaks"
    if band:
        info = (f"T = {float(T):g} K, harmonic, fixed cell; " if T else "") + \
            "bands: width at half maximum of every peak of the Gaussian fit; markers: the peaks"
    ncl = sum(r["nclip"] for r in res.values())
    if ncl and not band:
        info += f"\n{ncl} outlying values (bins below {100 * a.tail:g} % of the highest bar, or lone groups under 1 %) not drawn"
    if nf:
        info += f"\n{nf} flagged rows (imaginary modes / not converged) {'included' if a.include_flagged else 'left out'}"
    if a.ref == "formation":
        info += ("; " if nf else "\n") + ref_note.replace("F_form", r"$F_\mathrm{form}$")
    if a.ref == "mixing":
        info += ("; " if nf else "\n") + r"$F_\mathrm{lin}$: line between the lowest F of the end compositions"

    if bx is not None:
        branches(lambda p: bx, lambda x, p, q: q[0] + res[(x, p)]["db"], False)
        for (x, p), r in res.items():
            s = st[p]
            for i, (F0, err, h) in enumerate(r["peaks"]):
                main = i == r["main"]
                bx.errorbar([x * xscale], [F0 + r["db"]], yerr=[err], ls="", marker=s["marker"],
                            ms=3.4, mfc=s["color"] if main else "white", mec="white" if main else s["color"],
                            mew=0.4 if main else 0.7, ecolor=s["color"], elinewidth=0.6, capsize=0, zorder=6)
        bx.axhline(0, color=INK2, lw=0.6, ls=(0, (3, 2)))
        bx.set_ylabel(r"$\Delta F_\mathrm{peak}$ (meV/atom)")
        if a.ylim_b:
            bx.set_ylim(*a.ylim_b)
        bx.text(1.012, 0.5, f"({next(letters)})", transform=bx.transAxes, ha="left", va="center", rotation=270,
                fontsize=9)
        info += f"\n({'abcdefghijklmnop'[na]}) every peak minus the lowest peak (any phase) at the same composition; bars: bootstrap error"
    last = axes[-1]
    last.set_xlabel(f"{a.dopant} " + ("fraction of cations" if a.fraction else "concentration (% of cations)"))
    last.annotate(info, xy=(1, 0), xycoords=("axes fraction", last.xaxis.label), xytext=(0, -4),
                  textcoords="offset points", ha="right", va="top", fontsize=7, color=INK2)
    lo, hi = comps[0] * xscale, comps[-1] * xscale
    dd = 0 if rows else a.dodge * dx * len(phases) / 2
    far = (0.4 * dx if band else 1.15 * slot) + dd
    axes[0].set_xlim(lo - (far if sgn < 0 else 0.4 * dx + dd), hi + (far if sgn > 0 else 0.4 * dx + dd))

    if rows:                                     # one y label centred on the phase rows
        fig.canvas.draw()
        y0, y1 = axA[-1].get_position().y0, axA[0].get_position().y1
        lab_x = min(t.get_window_extent().x0 for ax in axA for t in ax.get_yticklabels() if t.get_text())
        lab_x = fig.transFigure.inverted().transform((lab_x, 0))[0]
        fig.text(lab_x - 0.012, (y0 + y1) / 2, ylab, rotation=90, ha="right", va="center", fontsize=10)
        if bx is not None:
            fig.align_ylabels([bx])
    for fmt in a.formats:
        fig.savefig(out / f"{stem}.{fmt}", dpi=300, bbox_inches="tight")
    with (out / f"{stem}_peaks.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_dopant", "phase", "peak_no", "main_peak", f"F_peak_meV_atom_ref_{a.ref}", "err_meV_atom",
                    "height_rel", "dF_to_lowest_peak_meV_atom", "n", "n_flagged", "n_gaussians",
                    "lowest_meV_atom", "median_meV_atom"])
        for (x, p), r in res.items():
            hmax = max(q[2] for q in r["peaks"])
            for i, (F0, err, h) in enumerate(r["peaks"]):
                w.writerow([f"{x:.6f}", p, i + 1, i == r["main"], f"{F0:.3f}", f"{err:.3f}", f"{h / hmax:.3f}",
                            f"{round(F0 + r['db'], 3) + 0.0:.3f}", len(r["y"]), nflag.get((x, p), 0), r["fit"]["k"],
                            f"{r['low']:.3f}", f"{r['med']:.3f}"])
    print(f"{len(data)} rows, {len(comps)} compositions, phases {phases}; flagged rows "
          f"{'included' if a.include_flagged else 'left out'}: {nf}; distributions with more than one peak: "
          f"{n_multi} of {len(res)}; bin {bw:g} meV/atom")
    print(f"y reference: {ref_note}" + (f"\n  {ref_detail}" if ref_detail else ""))
    print(f"wrote {out / stem}.{{{','.join(a.formats)}}} and {stem}_peaks.csv")


if __name__ == "__main__":
    main()
