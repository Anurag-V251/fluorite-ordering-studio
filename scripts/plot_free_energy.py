#!/usr/bin/env python3
"""
plot_free_energy.py - publication-ready figures from the CSV written by batch_free_energy.py

  (a) P(F)      one histogram per phase (all phases together sum to 1), each fitted with a mixture of
                Gaussians (maximum likelihood on the F values; 1..3 Gaussians, number chosen by BIC).
                Every maximum of a fitted distribution is a peak, so a phase can have several.
                The lowest-F peak of all phases is the reference; every peak gets
                    dF_peak = F_peak - F_peak(reference)   [meV/atom, bootstrap uncertainty in brackets]
                (one Gaussian that is merely skewed does not create a second peak: peaks are the maxima of the
                 fitted curve, and maxima without a clear dip between them are merged)
  (b) WC vs F   F on the x axis, Warren-Cowley alpha (shell 1) on the y axis, one colour + marker per phase;
                structures with imaginary phonon modes are hollow.
Legends are placed automatically where they cover no data (extra headroom is added if needed).

Both F axes:  F - F_peak(reference)  in meV/atom  (or absolute F in eV/atom with --absolute).

Output (next to the CSV, or --out-dir), Physical Review sizes, Times-like fonts (STIX), editable PDF text:
    <stem>_PF.pdf/.png          single column  (8.6 cm)
    <stem>_WC_vs_F.pdf/.png     single column  (8.6 cm)
    <stem>_combined.pdf/.png    two panels (a)(b), double column (17.8 cm)
    <stem>_peak_fit.csv         every peak: F, dF_peak, uncertainty, fitted Gaussians

Usage
    python plot_free_energy.py Ce_33.33/free_energy_300K.csv
    python plot_free_energy.py free_energy_300K.csv --bins 30 --exclude-imaginary
    python plot_free_energy.py free_energy_300K.csv --label parent="c (parent)" --label oIII="o-III"
    python plot_free_energy.py free_energy_300K.csv --font sans --formats pdf svg
    python plot_free_energy.py free_energy_300K.csv --components m=2      # force 2 Gaussians for m
"""
import argparse
import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import AutoMinorLocator  # noqa: E402

CM = 1 / 2.54
WIDTH_SINGLE, WIDTH_DOUBLE = 8.6 * CM, 17.8 * CM          # Physical Review column widths

# fixed colour + marker per phase (colour follows the phase, never its rank); blue/orange/aqua validated
# colour-blind safe as a set, grey = the initial (un-mapped) reference structures
PHASE_STYLE = OrderedDict([
    ("parent", dict(color="#52514e", marker="o", label="c (initial)")),
    ("t",      dict(color="#2a78d6", marker="s", label="t")),
    ("m",      dict(color="#eb6834", marker="^", label="m")),
    ("oIII",   dict(color="#1baf7a", marker="D", label="oIII")),
])
EXTRA_STYLES = [dict(color="#6250d6", marker="v"), dict(color="#e34948", marker="P"),
                dict(color="#eda100", marker="X"), dict(color="#e87ba4", marker="h")]
INK, INK2 = "#000000", "#444444"


# --------------------------------------------------------------------------- style
def set_style(font):
    serif = font == "serif"
    plt.rcParams.update({
        "font.family": "STIXGeneral" if serif else "DejaVu Sans",
        "mathtext.fontset": "stix" if serif else "dejavusans",
        "font.size": 9, "axes.labelsize": 9.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "legend.fontsize": 7.5, "legend.title_fontsize": 7.5, "legend.frameon": False,
        "legend.handlelength": 1.6, "legend.handletextpad": 0.5, "legend.borderaxespad": 0.4,
        "axes.linewidth": 0.6, "axes.labelpad": 3,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5, "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "lines.linewidth": 1.2, "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",       # editable text in PDF/SVG
    })


def finish_axes(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())


def chem(formula):
    """Zr8Ce4O24 -> Zr$_{8}$Ce$_{4}$O$_{24}$"""
    return re.sub(r"(\d+)", r"$_{\1}$", formula)


def unc(v, e, min_dec=1, max_dec=3):
    """+6.4(8), +2.50(3): value with one significant digit of uncertainty in brackets"""
    if not np.isfinite(e) or e <= 0:
        return f"{v:+.{min_dec}f}"
    dec = int(np.clip(-np.floor(np.log10(e)), min_dec, max_dec))
    d = int(round(e * 10 ** dec))
    if d >= 10 and dec > min_dec:                     # e.g. 0.096 -> 0.10 -> keep one digit
        dec -= 1
        d = int(round(e * 10 ** dec))
    return f"{v:+.{dec}f}({max(d, 1)})"


def place_legend(ax, make, pts, filled, loc="auto", step=0.06, max_steps=30):
    """
    Put the legend where it covers no data. pts = (N, 2) data-coordinate points of everything drawn.
    filled = True for histograms: everything below a bar top counts as data, so only 'upper' positions are
    tried. For each candidate (corner, columns) the y-range is extended (top for 'upper', bottom for 'lower')
    until the legend box is clear; the candidate needing the least extra room wins.
    """
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    y0, y1 = ax.get_ylim()
    if loc != "auto":
        return make(loc, 1)
    corners = ["upper right", "upper left"] + ([] if filled else ["lower right", "lower left"])
    best = None
    for ncol in (1, 2):
        for corner in corners:
            ax.set_ylim(y0, y1)
            leg = make(corner, ncol)
            for k in range(max_steps + 1):
                bb = leg.get_window_extent(renderer).transformed(ax.transData.inverted())
                lo, hi = ax.get_ylim()
                pad = 0.02 * (hi - lo)
                inx = (pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1)
                if filled:
                    clash = inx.any() and pts[inx, 1].max() > bb.y0 - pad
                else:
                    clash = (inx & (pts[:, 1] >= bb.y0 - pad) & (pts[:, 1] <= bb.y1 + pad)).any()
                if not clash:
                    break
                span = hi - lo
                ax.set_ylim(lo, hi + step * span) if corner.startswith("upper") else ax.set_ylim(lo - step * span, hi)
            leg.remove()
            if not clash and (best is None or k < best[0]):
                best = (k, corner, ncol)
    ax.set_ylim(y0, y1)
    if best is None:                                   # nothing fits: fall back to matplotlib's choice
        return make("best", 1)
    k, corner, ncol = best
    span = y1 - y0
    for _ in range(k):
        lo, hi = ax.get_ylim()
        span = hi - lo
        ax.set_ylim(lo, hi + step * span) if corner.startswith("upper") else ax.set_ylim(lo - step * span, hi)
    return make(corner, ncol)


# --------------------------------------------------------------------------- data
def style_for(phase, extra_used, labels):
    st = dict(PHASE_STYLE[phase]) if phase in PHASE_STYLE else None
    if st is None:
        if phase not in extra_used:
            extra_used[phase] = dict(EXTRA_STYLES[len(extra_used) % len(EXTRA_STYLES)], label=phase)
        st = dict(extra_used[phase])
    if phase in labels:
        st["label"] = labels[phase]
    return st


def load(csv_path, shell, exclude_imag):
    rows = [r for r in csv.DictReader(open(csv_path)) if r.get("status", "ok") == "ok"]
    if not rows:
        sys.exit("no successful rows in the CSV")
    out = []
    for r in rows:
        m = re.findall(r"WC_shell%d\(([^)]*)\)=([+-]?[\d.]+)" % shell, r.get("WC", ""))
        imag = int(float(r.get("n_imaginary_modes") or 0)) > 0
        if exclude_imag and imag:
            continue
        out.append(dict(phase=r.get("phase") or "parent", comp=r.get("composition", ""),
                        F=float(r["F_total_eV_atom"]), wc=float(m[0][1]) if m else np.nan,
                        pair=m[0][0] if m else "", imag=imag, T=r.get("T_K", "")))
    return out


LOG2PI = np.log(2 * np.pi)


def mix_pdf(x, w, mu, sig):
    x = np.asarray(x, float)[..., None]
    w, mu, sig = (np.asarray(v, float) for v in (w, mu, sig))
    return (w * np.exp(-0.5 * ((x - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))).sum(-1)


def em_fit(x, k, rng, n_init=6, iters=2000, init=None):
    """maximum-likelihood mixture of k Gaussians (1D EM, several starts). Returns (logL, w, mu, sigma)."""
    n, s = len(x), float(np.std(x)) or 1.0
    floor = 0.02 * s + 1e-9
    best = None
    for t in range(n_init):
        if t == 0 and init is not None:
            w, mu, sig = (np.array(v, float) for v in init)
        else:
            mu = np.quantile(x, (np.arange(k) + 0.5) / k) if t <= 1 else rng.choice(x, k, replace=False).astype(float)
            sig, w = np.full(k, max(s / k, floor)), np.full(k, 1.0 / k)
        prev = -np.inf
        for _ in range(iters):
            logp = np.log(w) - np.log(sig) - 0.5 * ((x[:, None] - mu) / sig) ** 2 - 0.5 * LOG2PI
            m = logp.max(1, keepdims=True)
            lse = m + np.log(np.exp(logp - m).sum(1, keepdims=True))
            ll = float(lse.sum())
            r = np.exp(logp - lse)
            nk = r.sum(0) + 1e-12
            w = nk / n
            mu = (r * x[:, None]).sum(0) / nk
            sig = np.maximum(np.sqrt((r * (x[:, None] - mu) ** 2).sum(0) / nk), floor)
            if ll - prev < 1e-10 * abs(ll) + 1e-12:
                break
            prev = ll
        if best is None or ll > best[0]:
            best = (ll, w, mu, sig)
    return best


def display_params(f, bw):
    """mixture parameters for drawing: a phase whose values are all identical is drawn as a single-bin peak
    (sigma = bw / sqrt(2 pi), so the curve top equals the bar height) instead of an invisible spike"""
    sig = np.asarray(f["sigma"], float)
    if f.get("degenerate"):
        sig = np.full_like(sig, bw / np.sqrt(2 * np.pi))
    return f["w"], f["mu"], sig


def find_modes(w, mu, sig, lo, hi, min_peak, min_dip=0.1, npts=4000):
    """local maxima of the mixture density; peaks lower than min_peak x the highest one, or not separated
    from a higher neighbour by a dip of at least min_dip (fraction of the lower peak), are dropped"""
    g = np.linspace(lo, hi, npts)
    d = mix_pdf(g, w, mu, sig)
    idx = list(np.nonzero((d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:]))[0] + 1)
    idx = [i for i in idx if d[i] >= min_peak * d.max()]
    changed = True
    while changed and len(idx) > 1:
        changed = False
        for a_, b_ in zip(idx[:-1], idx[1:]):
            if d[a_:b_ + 1].min() > (1 - min_dip) * min(d[a_], d[b_]):
                idx.remove(a_ if d[a_] < d[b_] else b_)
                changed = True
                break
    return g[idx], d[idx]


def fit_peaks(xp, rng, kmax=3, k_fixed=None, min_peak=0.15, n_boot=100, min_presence=0.7, dbic=10.0):
    """
    Mixture-of-Gaussians fit to one phase's F values (meV/atom). Number of components: k_fixed, or the
    smallest k = 1..kmax whose BIC is within dbic of every larger k (an extra Gaussian must lower the BIC by
    more than dbic = strong evidence). Peaks = maxima of the fitted density. Bootstrap re-fits (same k) give
    the uncertainties; a peak that reappears (within one Gaussian width) in fewer than min_presence of the
    re-fits is not robust and is dropped (the highest peak is always kept).
    Returns dict(n, k, w, mu, sigma, peaks=[(F, err, height)], dropped, ok).
    """
    xp = np.asarray(xp, float)
    n = len(xp)
    s = float(np.std(xp)) or 1.0
    if n and np.ptp(xp) <= 1e-9 * max(1.0, abs(float(np.mean(xp)))):   # all values identical (e.g. fixed cells)
        m0 = float(np.mean(xp))
        return dict(n=n, k=1, w=np.array([1.0]), mu=np.array([m0]), sigma=np.array([1e-9 * max(1.0, abs(m0))]),
                    bic={}, ok=True, dropped=[], degenerate=True, peaks=[(m0, 0.0, 1.0)])
    if n < 10:                                                       # too few points: one Gaussian from moments
        w, mu, sig = np.array([1.0]), np.array([xp.mean()]), np.array([s])
        return dict(n=n, k=1, w=w, mu=mu, sigma=sig, bic={}, ok=False, dropped=[],
                    peaks=[(float(xp.mean()), s / np.sqrt(max(n, 1)), 1.0)])
    ks = [k_fixed] if k_fixed else [k for k in range(1, kmax + 1) if n >= 10 * k]
    fits, bic = {}, {}
    for k in ks:
        fits[k] = em_fit(xp, k, rng)
        bic[k] = -2 * fits[k][0] + (3 * k - 1) * np.log(n)
    k = ks[0]
    for kk in ks[1:]:                                                # accept more Gaussians only for strong evidence
        if bic[kk] < bic[k] - dbic:
            k = kk
    _, w, mu, sig = fits[k]
    lo, hi = xp.min() - 3 * s, xp.max() + 3 * s
    modes, heights = find_modes(w, mu, sig, lo, hi, min_peak)
    if not len(modes):                                                # grid too coarse for a very narrow fit
        j = int(np.argmax(w / sig))
        modes, heights = np.array([mu[j]]), mix_pdf([mu[j]], w, mu, sig)
    width = [float(sig[int(np.argmin(np.abs(mu - m0)))]) for m0 in modes]   # match radius per peak
    boots, n_ok = [[] for _ in modes], 0
    for _ in range(n_boot):
        xb = rng.choice(xp, n, replace=True)
        _, wb, mub, sb = em_fit(xb, k, rng, n_init=1, iters=300, init=(w, mu, sig))   # warm start
        mb, _ = find_modes(wb, mub, sb, lo, hi, min_peak, npts=1500)
        if not len(mb):
            continue
        n_ok += 1
        for i, m0 in enumerate(modes):
            near = mb[np.argmin(np.abs(mb - m0))]
            if abs(near - m0) <= width[i]:
                boots[i].append(near)
    peaks, dropped = [], []
    top = int(np.argmax(heights))
    for i, m0 in enumerate(modes):
        presence = len(boots[i]) / n_ok if n_ok else 1.0
        if n_ok and presence < min_presence and i != top:
            dropped.append((float(m0), presence))
            continue
        if len(boots[i]) > 5:
            err = float(np.std(boots[i], ddof=1))
        else:                                                         # no bootstrap: nearest component
            j = int(np.argmin(np.abs(mu - m0)))
            err = float(sig[j] / np.sqrt(max(w[j] * n, 1)))
        peaks.append((float(m0), err, float(heights[i])))
    return dict(n=n, k=k, w=w, mu=mu, sigma=sig, bic=bic, ok=True, peaks=peaks, dropped=dropped)


# --------------------------------------------------------------------------- panels
def draw_pf(ax, P):
    """P(F) histograms + Gaussian-mixture fits, dashed line at every peak"""
    edges, x0, ph_of, fits, styles, to_axis = (P[k] for k in ("edges", "x0", "ph_of", "fits", "styles", "to_axis"))
    ex = to_axis(edges)
    bw = edges[1] - edges[0]
    fine = np.linspace(edges[0] - 2 * bw, edges[-1] + 2 * bw, 1200)
    ymax, handles, labels, pts = 0.0, [], [], []
    for p in sorted(P["phases"], key=lambda q: fits[q]["peaks"][0][0]):
        st, f = styles[p], fits[p]
        scale = f["n"] / len(x0) * bw                                   # density -> histogram height
        h = np.histogram(x0[ph_of == p], bins=edges)[0] / len(x0)
        ax.stairs(h, ex, fill=True, color=st["color"], alpha=0.18, lw=0, zorder=1)
        ax.stairs(h, ex, color=st["color"], lw=0.5, alpha=0.9, zorder=2)
        g = scale * mix_pdf(fine, *display_params(f, bw))
        if P["show_components"] and f["k"] > 1:
            for wj, mj, sj in zip(*display_params(f, bw)):
                ax.plot(to_axis(fine), scale * mix_pdf(fine, [wj], [mj], [sj]), color=st["color"], lw=0.6,
                        ls=(0, (1, 1.2)), zorder=3)
        ax.plot(to_axis(fine), g, color=st["color"], lw=1.3, zorder=4)
        for F0, _, _ in f["peaks"]:                                     # peak line up to the curve only
            ax.vlines(to_axis(F0), 0, scale * mix_pdf([F0], *display_params(f, bw))[0], color=st["color"],
                      lw=0.7, ls=(0, (3, 2)), zorder=3)
        ymax = max(ymax, h.max(), g.max())
        for k in range(len(h)):                                         # bar tops, sampled across each bin
            xs = np.linspace(ex[k], ex[k + 1], 5)
            pts.append(np.column_stack([xs, np.full(5, h[k])]))
        pts.append(np.column_stack([to_axis(fine), g]))
        handles.append((Patch(facecolor=st["color"], alpha=0.25, lw=0), Line2D([], [], color=st["color"], lw=1.3)))
        vals = ["0 (ref.)" if pk["is_ref"] else unc(pk["dF"], pk["dF_err"]) for pk in f["peak_info"]]
        labels.append(f"{st['label']}: " + ", ".join(vals) + ("" if f["ok"] else r"$^{*}$"))
    ax.set_ylim(0, ymax * 1.05)
    if not P["absolute"]:
        ax.set_xlim(ex[0] - bw, ex[-1] + bw)
    ax.set_xlabel(P["flabel"])
    ax.set_ylabel(r"$P(F)$")
    finish_axes(ax)
    pts = np.vstack(pts)
    P["_handles"] = handles                                             # for scripts that build their own legend
    if not P.get("legend", True):
        return pts
    place_legend(ax, lambda loc, nc: ax.legend(handles, labels, loc=loc, ncol=nc, alignment="left",
                                               title=r"$\Delta F_\mathrm{peak}$ (meV/atom)", columnspacing=1.0),
                 pts, filled=True, loc=P["legend_loc"])
    return pts                                                          # drawn data, for label placement


def draw_wc(ax, P):
    """F on x, Warren-Cowley alpha on y"""
    sub, x, styles, phases = P["sub"], P["x"], P["styles"], P["phases"]
    uw = np.unique([d["wc"] for d in sub if not np.isnan(d["wc"])])
    dodge = {p: 0.0 for p in phases}
    if P["dodge"] and len(phases) > 1 and 1 < len(uw) <= 40:
        step = np.diff(uw).min() * 0.6 / len(phases)
        dodge = {p: (k - (len(phases) - 1) / 2) * step for k, p in enumerate(phases)}
    ax.axhline(0, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=1)
    n_imag, allpts = 0, []
    msize = P["msize"] * float(np.clip((800.0 / max(len(sub), 1)) ** 0.35, 0.45, 1.0))   # smaller when crowded
    for p in phases:
        st = styles[p]
        for hollow in (False, True):
            pts = [(xi, d["wc"]) for xi, d in zip(x, sub) if d["phase"] == p and d["imag"] == hollow
                   and not np.isnan(d["wc"])]
            if not pts:
                continue
            ff, wc = np.array(pts).T
            wc = wc + dodge[p]
            allpts.append(np.column_stack([ff, wc]))
            n_imag += len(pts) if hollow else 0
            ax.scatter(ff, wc, s=msize, marker=st["marker"],
                       facecolors="none" if hollow else st["color"],
                       edgecolors=st["color"] if hollow else "white",
                       linewidths=0.6 if hollow else 0.25, alpha=0.9, zorder=3, rasterized=P["raster"])
    handles = [Line2D([], [], ls="", marker=styles[p]["marker"], ms=4.2, mfc=styles[p]["color"], mec="white",
                      mew=0.3) for p in phases]
    labels = [styles[p]["label"] for p in phases]
    if n_imag:
        handles.append(Line2D([], [], ls="", marker="o", ms=4.2, mfc="none", mec=INK2, mew=0.6))
        labels.append("imag. modes")
    ax.set_ylabel(r"Warren–Cowley $\alpha^{(%d)}_{\mathrm{%s}}$" % (P["shell"], P["pair"].replace("-", "{-}"))
                  if P["pair"] else f"Warren–Cowley $\\alpha^{{({P['shell']})}}$")
    ax.set_xlabel(P["flabel"])
    finish_axes(ax)
    if allpts:
        allpts = np.vstack(allpts)
        span = np.ptp(allpts[:, 1]) or 1.0
        ax.set_ylim(allpts[:, 1].min() - 0.06 * span, allpts[:, 1].max() + 0.06 * span)
        xl = ax.get_xlim()                                              # keep the legend off the alpha = 0 line
        allpts = np.vstack([allpts, np.column_stack([np.linspace(xl[0], xl[1], 200), np.zeros(200)])])
        place_legend(ax, lambda loc, nc: ax.legend(handles, labels, loc=loc, ncol=nc, handletextpad=0.2,
                                                   borderaxespad=0.3, columnspacing=0.8),
                     allpts, filled=False, loc=P["legend_loc_wc"])
    return any(dodge.values())


def info_text(ax, P):
    if P["info"]:
        ax.text(0.0, 1.015, P["info"], transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=INK2)


def save(fig, base, formats):
    out = []
    for fmt in formats:
        fn = f"{base}.{fmt}"
        fig.savefig(fn)
        out.append(fn)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="Publication figures: P(F) with Gaussian fits, and WC vs F")
    ap.add_argument("csv")
    ap.add_argument("--bins", default="auto", help="number of bins or 'auto' (default)")
    ap.add_argument("--shell", type=int, default=1, help="Warren-Cowley shell (default 1)")
    ap.add_argument("--absolute", action="store_true", help="absolute F in eV/atom on the axes")
    ap.add_argument("--exclude-imaginary", action="store_true", help="drop structures with imaginary modes")
    ap.add_argument("--no-dodge", dest="dodge", action="store_false",
                    help="do not offset the phases vertically at each WC value")
    ap.add_argument("--label", action="append", default=[], metavar="PHASE=TEXT",
                    help='rename a phase in the figures, e.g. --label parent="c (parent)" (repeatable)')
    ap.add_argument("--font", choices=["serif", "sans"], default="serif",
                    help="serif = Times-like STIX (Phys. Rev. style, default); sans = DejaVu Sans")
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"], choices=["pdf", "png", "svg", "eps"])
    ap.add_argument("--legend-loc", default="auto",
                    help="legend position in (a): auto (default, clear of the data) or e.g. 'upper right'")
    ap.add_argument("--legend-loc-wc", default="auto", help="legend position in (b), as --legend-loc")
    ap.add_argument("--no-info", action="store_true", help="omit the composition/temperature line above the plots")
    ap.add_argument("--out-dir", help="where to write the figures (default: next to the CSV)")
    ap.add_argument("--max-components", type=int, default=3,
                    help="most Gaussians per phase; the number is chosen by BIC (default 3)")
    ap.add_argument("--components", action="append", default=[], metavar="[PHASE=]K",
                    help="force K Gaussians, for all phases (--components 2) or one (--components m=2); repeatable")
    ap.add_argument("--min-peak", type=float, default=0.15,
                    help="ignore maxima lower than this fraction of the phase's highest peak (default 0.15)")
    ap.add_argument("--bootstrap", type=int, default=100, help="bootstrap re-fits for peak uncertainties (default 100)")
    ap.add_argument("--no-components", dest="show_components", action="store_false",
                    help="do not draw the individual Gaussians (thin dotted lines)")
    ap.add_argument("--reference", choices=["lowest", "main"], default="lowest",
                    help="0 of the F axis: the lowest-F peak of all phases (default), or 'main' = the lowest of "
                         "the phases' tallest peaks (ignores small low-F side peaks)")
    ap.add_argument("--seed", type=int, default=0, help="random seed for the fits/bootstrap (default 0)")
    a = ap.parse_args()
    forced = {}
    for item in a.components:
        k_, _, v_ = item.rpartition("=")
        forced[k_ or "*"] = int(v_)

    labels = {}
    for item in a.label:
        if "=" not in item:
            sys.exit(f"--label needs PHASE=TEXT, got {item!r}")
        k, v = item.split("=", 1)
        labels[k] = v

    set_style(a.font)
    csv_path = Path(a.csv)
    data = load(csv_path, a.shell, a.exclude_imaginary)
    out_dir = Path(a.out_dir) if a.out_dir else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    T = data[0]["T"]

    comps = list(OrderedDict.fromkeys(d["comp"] for d in data))
    extra = OrderedDict()
    for comp in comps:
        sub = [d for d in data if d["comp"] == comp]
        Fmin = min(d["F"] for d in sub)
        x0 = np.array([1000 * (d["F"] - Fmin) for d in sub])                  # meV/atom above F_min
        ph_of = np.array([d["phase"] for d in sub])
        phases = [p for p in list(PHASE_STYLE) + sorted(set(ph_of) - set(PHASE_STYLE)) if p in set(ph_of)]
        styles = {p: style_for(p, extra, labels) for p in phases}
        tag = csv_path.stem + (f"_{comp}" if len(comps) > 1 and comp else "")

        bins = int(a.bins) if str(a.bins).isdigit() else "auto"
        edges = np.histogram_bin_edges(x0, bins=bins)
        if len(edges) < 2 or edges[-1] == edges[0]:
            edges = np.linspace(x0.min() - 0.5, x0.max() + 0.5, 11)
        rng = np.random.default_rng(a.seed)
        fits = {p: fit_peaks(x0[ph_of == p], rng, kmax=a.max_components, k_fixed=forced.get(p, forced.get("*")),
                             min_peak=a.min_peak, n_boot=a.bootstrap) for p in phases}
        # reference = lowest-F peak of all phases (or, --reference main, the lowest of the phases' tallest peaks)
        cands = [(p, pk) for p in phases for pk in fits[p]["peaks"]]
        if a.reference == "main":
            cands = [(p, max(fits[p]["peaks"], key=lambda q: q[2])) for p in phases]
        ref, (mu_ref, err_ref, _) = min(cands, key=lambda t: t[1][0])
        for p in phases:
            fits[p]["peak_info"] = [dict(F=F0, err=e, height=hgt, is_ref=(p == ref and F0 == mu_ref),
                                         dF=F0 - mu_ref,
                                         dF_err=0.0 if (p == ref and F0 == mu_ref) else float(np.hypot(e, err_ref)))
                                    for F0, e, hgt in fits[p]["peaks"]]

        ref_lbl = styles[ref]["label"].split(" (")[0]
        if a.absolute:
            to_axis = lambda v: Fmin + v / 1000.0                            # noqa: E731
            flabel = r"$F$ (eV/atom)"
        else:
            to_axis = lambda v: v - mu_ref                                   # noqa: E731
            flabel = r"$F - F_\mathrm{peak}^{\,\mathrm{%s}}$ (meV/atom)" % ref_lbl.replace(" ", r"\ ")
        pair = next((d["pair"] for d in sub if d["pair"]), "")
        info = "" if a.no_info else ", ".join(b for b in (chem(comp) if comp else "",
                                                           f"$T$ = {float(T):g} K" if T else "",
                                                           f"$N$ = {len(sub)}") if b)
        P = dict(edges=edges, x0=x0, ph_of=ph_of, fits=fits, styles=styles, ref=ref, to_axis=to_axis,
                 show_components=a.show_components,
                 phases=phases, absolute=a.absolute, flabel=flabel, sub=sub, x=to_axis(x0), shell=a.shell,
                 pair=pair, dodge=a.dodge, legend_loc=a.legend_loc, legend_loc_wc=a.legend_loc_wc, info=info,
                 raster=len(sub) > 3000)

        written = []
        # (a) alone, single column
        fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 0.78 * WIDTH_SINGLE))
        P["msize"] = 9
        draw_pf(ax, P)
        info_text(ax, P)
        written += save(fig, out_dir / f"{tag}_PF", a.formats)

        # (b) alone, single column
        fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 0.82 * WIDTH_SINGLE))
        dodged = draw_wc(ax, P)
        info_text(ax, P)
        written += save(fig, out_dir / f"{tag}_WC_vs_F", a.formats)

        # (a)+(b), double column
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(WIDTH_DOUBLE, 0.40 * WIDTH_DOUBLE),
                                       gridspec_kw=dict(wspace=0.28))
        P["msize"] = 11
        draw_pf(axa, P)
        draw_wc(axb, P)
        info_text(axa, P)
        for lab, axx in (("(a)", axa), ("(b)", axb)):
            axx.text(-0.16, 1.0, lab, transform=axx.transAxes, ha="left", va="bottom", fontsize=10,
                     fontweight="bold" if a.font == "sans" else "normal")
        written += save(fig, out_dir / f"{tag}_combined", a.formats)

        # fit table: one row per peak
        f_csv = out_dir / f"{tag}_peak_fit.csv"
        order = sorted(phases, key=lambda q: fits[q]["peaks"][0][0])
        with open(f_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["phase", "peak", "n_structures_phase", "n_gaussians", "peak_F_eV_atom",
                        "peak_minus_ref_meV_atom", "err_meV_atom", "peak_height_rel", "reference",
                        "gaussians (weight, mean meV/atom above F_min, sigma)", "method"])
            for p in order:
                f = fits[p]
                hmax = max(pk["height"] for pk in f["peak_info"])
                comps = "; ".join(f"{wj:.3f}, {mj:.3f}, {sj:.3f}" for wj, mj, sj in zip(f["w"], f["mu"], f["sigma"]))
                for i, pk in enumerate(f["peak_info"], 1):
                    w.writerow([styles[p]["label"], i, f["n"], f["k"], f"{Fmin + pk['F'] / 1000:.8f}",
                                f"{pk['dF']:.4f}", f"{pk['dF_err']:.4f}", f"{pk['height'] / hmax:.3f}",
                                "yes" if pk["is_ref"] else "", comps,
                                f"Gaussian mixture (BIC), bootstrap x{a.bootstrap}" if f["ok"] else "mean/std (n<10)"])
        written.append(str(f_csv))

        print(f"{comp or 'all'}: {len(sub)} structures; peaks relative to the reference peak ({styles[ref]['label']}, {a.reference}):")
        for p in order:
            f = fits[p]
            bic = "  ".join(f"k={k}: {b:.0f}" for k, b in f["bic"].items())
            print(f"  {styles[p]['label']:<14} n={f['n']:<6} {f['k']} Gaussian(s)" + (f"   [BIC {bic}]" if bic else ""))
            for i, pk in enumerate(f["peak_info"], 1):
                print(f"      peak {i}: F = {Fmin + pk['F'] / 1000:.6f} eV/atom   "
                      f"ΔF = {pk['dF']:+8.2f} ± {pk['dF_err']:.2f} meV/atom" + ("   (reference)" if pk["is_ref"] else ""))
            for F0, pres in f["dropped"]:
                print(f"      (maximum at {F0 - mu_ref:+.2f} meV/atom ignored: found in only {100 * pres:.0f}% of "
                      f"bootstrap re-fits)")
        for fn in written:
            print(f"  wrote {fn}")
        if dodged:
            print("  caption note: in (b) the phases are offset slightly in alpha at each alpha value for visibility.")
        if any(d["imag"] for d in sub):
            print("  caption note: hollow symbols = structures with imaginary phonon modes.")
        if any(not fits[p]["ok"] for p in phases):
            print("  caption note: * = fewer than 10 structures, single Gaussian from mean and standard deviation.")
        print("  caption note: curves = Gaussian-mixture fits (number of Gaussians by BIC); dashed lines = peaks of "
              "the fitted distributions; uncertainties from bootstrap resampling.")


if __name__ == "__main__":
    main()
