#!/usr/bin/env python3
"""
plot_cell_annotations.py - P(F) histograms (same fits and peaks as plot_free_energy.py) with every peak
annotated by the cell shape of the structures in that peak:

  (a) c/a   ratio of the normalised lattice parameters, c_n / a_n
  (b) θ     the cell angle that deviates most from 90° (folded to >= 90°, e.g. the monoclinic β)

Normalised lattice parameters
    The supercell is n1 x n2 x n3 conventional fluorite cells. Each lattice vector length is divided by its
    repeat count, L_i / n_i, giving the conventional-cell parameters (a_n, b_n, c_n) whatever the supercell
    shape. Repeats are found from the number of cations (4 per conventional cell) and the cell lengths
    (each L_i / n_i must be a fluorite-like 4.0-6.5 Å; the most cubic-like choice wins), or set with --repeats.

Which axis is c - decided per structure, because it is along x, y or z in different POSCARs
    c = the longest normalised axis (default; in t-, m- and oIII-ZrO2 the distinct axis is the long one),
    a_n = mean of the other two,   c/a = c_n / a_n  (>= 1).
    --c-axis unique instead takes the axis that differs most from the other two (allows c/a < 1); for nearly
    cubic cells that choice flips between the longest and the shortest axis, which splits the c/a distribution
    of cubic structures into two artificial bumps around 1, hence the default.

Which structures belong to a peak
    The fitted distribution of each phase is split at its minima between peaks; every structure of that phase
    in a peak's basin counts for that peak. The label gives the mean over those structures, with the standard
    deviation across them in brackets (a spread, not an error bar): e.g. 1.018(4) = 1.018 ± 0.004.

Cell parameters come from the CSV written by batch_free_energy.py (a_A ... gamma_deg: after relaxation when
--relax was used). Needs plot_free_energy.py in the same folder.

Usage
    python plot_cell_annotations.py Ce_33.33/free_energy_300K.csv
    python plot_cell_annotations.py free_energy_300K.csv --reference main --exclude-imaginary
    python plot_cell_annotations.py free_energy_300K.csv --repeats 2 2 2

Output (next to the CSV, or --out-dir)
    <stem>_PF_ca.pdf/.png            single column, peaks labelled with c/a
    <stem>_PF_angle.pdf/.png         single column, peaks labelled with θ
    <stem>_PF_cell_combined.pdf/.png (a) c/a and (b) θ, double column
    <stem>_peak_cell.csv             per peak: F, ΔF, structures in the peak, c/a, a_n, c_n, θ (mean, std)
    <stem>_cell_params.csv           per structure: repeats, normalised a/b/c along x/y/z, c axis, c/a, θ
"""
import argparse
import csv
import itertools
import re
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import plot_free_energy as pfe
except ImportError:
    sys.exit("plot_free_energy.py must be in the same folder as this script")
_needed = ["fit_peaks", "mix_pdf", "set_style", "finish_axes", "chem", "style_for", "info_text", "save",
           "PHASE_STYLE", "WIDTH_SINGLE", "WIDTH_DOUBLE", "INK"]
_missing = [n for n in _needed if not hasattr(pfe, n)]
if _missing:
    sys.exit("the plot_free_energy.py next to this script is an older version (missing: %s).\n"
             "Replace it with the latest plot_free_energy.py (the one with multi-peak fits)." % ", ".join(_missing))

_dp = getattr(pfe, "display_params", lambda f, bw: (f["w"], f["mu"], f["sigma"]))  # drawing parameters
import matplotlib.pyplot as plt  # noqa: E402  (backend already set by plot_free_energy)
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402


def draw_hist(ax, P):
    """histograms + mixture fits + dashed peak lines (no legend); returns drawn points and legend handles"""
    edges, x0, ph_of, fits, styles, to_axis = (P[k] for k in ("edges", "x0", "ph_of", "fits", "styles", "to_axis"))
    ex = to_axis(edges)
    bw = edges[1] - edges[0]
    fine = np.linspace(edges[0] - 2 * bw, edges[-1] + 2 * bw, 1200)
    ymax, handles, pts = 0.0, [], []
    for p in sorted(P["phases"], key=lambda q: fits[q]["peaks"][0][0]):
        st, f = styles[p], fits[p]
        scale = f["n"] / len(x0) * bw
        h = np.histogram(x0[ph_of == p], bins=edges)[0] / len(x0)
        ax.stairs(h, ex, fill=True, color=st["color"], alpha=0.18, lw=0, zorder=1)
        ax.stairs(h, ex, color=st["color"], lw=0.5, alpha=0.9, zorder=2)
        g = scale * pfe.mix_pdf(fine, *_dp(f, bw))
        if P["show_components"] and f["k"] > 1:
            for wj, mj, sj in zip(*_dp(f, bw)):
                ax.plot(to_axis(fine), scale * pfe.mix_pdf(fine, [wj], [mj], [sj]), color=st["color"], lw=0.6,
                        ls=(0, (1, 1.2)), zorder=3)
        ax.plot(to_axis(fine), g, color=st["color"], lw=1.3, zorder=4)
        for F0, _, _ in f["peaks"]:
            ax.vlines(to_axis(F0), 0, scale * pfe.mix_pdf([F0], *_dp(f, bw))[0], color=st["color"],
                      lw=0.7, ls=(0, (3, 2)), zorder=3)
        ymax = max(ymax, h.max(), g.max())
        for k in range(len(h)):
            pts.append(np.column_stack([np.linspace(ex[k], ex[k + 1], 5), np.full(5, h[k])]))
        pts.append(np.column_stack([to_axis(fine), g]))
        handles.append((Patch(facecolor=st["color"], alpha=0.25, lw=0), Line2D([], [], color=st["color"], lw=1.3)))
    ax.set_ylim(0, ymax * 1.05)
    ax.set_xlim(ex[0] - bw, ex[-1] + bw)
    ax.set_xlabel(P["flabel"])
    ax.set_ylabel(r"$P(F)$")
    pfe.finish_axes(ax)
    return np.vstack(pts), handles

AXES = "abc"
ANGLE_NAMES = {0: (r"\alpha", "b,c"), 1: (r"\beta", "a,c"), 2: (r"\gamma", "a,b")}


# --------------------------------------------------------------------------- cell analysis
def n_cations(row):
    comp = row.get("composition", "")
    parts = re.findall(r"([A-Z][a-z]?)(\d*)", comp)
    if parts and any(el for el, _ in parts):
        n = sum(int(c or 1) for el, c in parts if el and el != "O")
        if n:
            return n
    na = float(row.get("n_atoms") or 0)
    return int(round(na / 3)) if na else None


def detect_repeats(lengths, ncat, lo=4.0, hi=6.5):
    """(n1, n2, n3) with n1*n2*n3 = ncat/4 and every L_i/n_i in [lo, hi] Å; most cubic-like if several"""
    if not ncat or ncat % 4:
        return None
    M = ncat // 4
    best = None
    for n1 in range(1, M + 1):
        if M % n1:
            continue
        for n2 in range(1, M // n1 + 1):
            if (M // n1) % n2:
                continue
            n = np.array([n1, n2, M // (n1 * n2)])
            conv = lengths / n
            if np.all((conv >= lo) & (conv <= hi)):
                spread = float(np.std(conv) / np.mean(conv))
                if best is None or spread < best[0]:
                    best = (spread, tuple(int(v) for v in n))
    return best[1] if best else None


def cell_shape(lengths, angles, repeats, c_mode="longest"):
    """normalised lengths, c axis, c/a, and the angle deviating most from 90°"""
    Ln = np.asarray(lengths, float) / np.asarray(repeats, float)
    if c_mode == "unique":
        dev = [abs(Ln[i] - np.mean(np.delete(Ln, i))) for i in range(3)]
        ic = int(np.argmax(dev))
    else:
        ic = int(np.argmax(Ln))
    a_n = float(np.mean(np.delete(Ln, ic)))
    ang = np.asarray(angles, float)
    k = int(np.argmax(np.abs(ang - 90.0)))
    theta = 90.0 + abs(ang[k] - 90.0)                                 # folded: 81° and 99° are the same cell
    return dict(Ln=Ln, c_axis=ic, a_n=a_n, c_n=float(Ln[ic]), ca=float(Ln[ic] / a_n), theta=float(theta),
                theta_idx=k)


def load_rows(csv_path, exclude_imag, repeats_override, c_mode="longest"):
    rows = [r for r in csv.DictReader(open(csv_path)) if r.get("status", "ok") == "ok"]
    if not rows:
        sys.exit("no successful rows in the CSV")
    need = ["a_A", "b_A", "c_A", "alpha_deg", "beta_deg", "gamma_deg"]
    if any(c not in rows[0] for c in need):
        sys.exit(f"the CSV has no cell columns ({', '.join(need)}); it must come from batch_free_energy.py")
    out, bad = [], 0
    for r in rows:
        imag = int(float(r.get("n_imaginary_modes") or 0)) > 0
        if exclude_imag and imag:
            continue
        try:
            L = np.array([float(r[c]) for c in need[:3]])
            A = np.array([float(r[c]) for c in need[3:]])
        except ValueError:
            bad += 1
            continue
        rep = tuple(repeats_override) if repeats_override else detect_repeats(L, n_cations(r))
        if rep is None:
            bad += 1
            continue
        sh = cell_shape(L, A, rep, c_mode)
        out.append(dict(file=r.get("file", ""), phase=r.get("phase") or "parent", comp=r.get("composition", ""),
                        F=float(r["F_total_eV_atom"]), T=r.get("T_K", ""), imag=imag, rep=rep, **sh))
    if bad:
        print(f"note: {bad} structure(s) skipped (missing cell data, or supercell repeats not found: use --repeats)")
    if not out:
        sys.exit("nothing to plot")
    return out


# --------------------------------------------------------------------------- peaks -> structures
def basins(f, lo, hi, npts=4000):
    """split points between consecutive peaks = minimum of the fitted density there"""
    peaks = sorted(pk["F"] for pk in f["peak_info"])
    g = np.linspace(lo, hi, npts)
    d = pfe.mix_pdf(g, f["w"], f["mu"], f["sigma"])
    cuts = []
    for p1, p2 in zip(peaks[:-1], peaks[1:]):
        m = (g > p1) & (g < p2)
        cuts.append(float(g[m][np.argmin(d[m])]) if m.any() else 0.5 * (p1 + p2))
    return peaks, [-np.inf] + cuts + [np.inf]


def summarise(vals):
    vals = np.asarray(vals, float)
    if not len(vals):
        return np.nan, np.nan
    return float(vals.mean()), float(vals.std(ddof=1)) if len(vals) > 1 else 0.0


def spread_txt(m, s, dec):
    """1.018(4): mean with the standard deviation in the last digit(s)"""
    if not np.isfinite(m):
        return "–"
    if not np.isfinite(s) or s <= 0:
        return f"{m:.{dec}f}"
    d = int(round(s * 10 ** dec))
    while d >= 100 and dec > 0:
        dec -= 1
        d = int(round(s * 10 ** dec))
    return f"{m:.{dec}f}({max(d, 1)})"


# --------------------------------------------------------------------------- label placement
def spread_1d(d, w, lo, hi, gap, iters=1000):
    """centres as close as possible to d (sorted), no overlaps, inside [lo, hi] (display units)"""
    c = np.array(d, float)
    w = np.asarray(w, float)
    for _ in range(iters):
        moved = False
        for i in range(len(c) - 1):
            ov = (w[i] + w[i + 1]) / 2 + gap - (c[i + 1] - c[i])
            if ov > 1e-6:
                c[i] -= ov / 2
                c[i + 1] += ov / 2
                moved = True
        c = np.clip(c, lo + w / 2, hi - w / 2)
        if not moved:
            break
    return c


def annotated_panel(ax, P, items, legend_title, fontsize=7.0):
    """
    Histograms + fits, a one-row legend of phase names at the top, and below it a band of labels (one per
    peak) spread sideways so they never overlap, each joined to its peak by a thin line in the phase colour.
    The y range is set so that data, labels and legend each get their own space.
    """
    pts, handles = draw_hist(ax, P)
    fig = ax.figure
    r = fig.canvas.get_renderer()
    order = sorted(P["phases"], key=lambda q: P["fits"][q]["peaks"][0][0])
    names = [P["styles"][p]["label"] for p in order]
    ncol = len(names)
    while True:
        leg = ax.legend(handles, names, loc="upper left", ncol=ncol, title=legend_title, alignment="left",
                        handlelength=1.3, columnspacing=0.9, borderaxespad=0.35)
        if leg.get_window_extent(r).width <= ax.get_window_extent(r).width * 0.98 or ncol == 1:
            break
        leg.remove()
        ncol = max(1, (ncol + 1) // 2)
    axbb = ax.get_window_extent(r)
    H, pad = axbb.height, 3.0
    hl = axbb.y1 - leg.get_window_extent(r).y0                        # legend height from the top of the axes

    texts = []
    for it in items:
        t = ax.text(0, 0, it["text"], ha="center", va="bottom", fontsize=fontsize, color=pfe.INK, zorder=7,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.95))
        texts.append(t)
    wd = [t.get_window_extent(r).width for t in texts]
    hb = max(t.get_window_extent(r).height for t in texts)
    gap = 4.0
    idx = sorted(range(len(items)), key=lambda i: items[i]["x"])
    dx = {i: ax.transData.transform((items[i]["x"], 0))[0] for i in idx}  # x mapping does not depend on ylim
    best = None
    for nrow in range(1, min(3, len(items)) + 1):                     # fewest rows with small label shifts
        lay, disp, ok = {}, 0.0, True
        for rrow in range(nrow):
            mem = idx[rrow::nrow]
            if sum(wd[i] for i in mem) + gap * len(mem) > axbb.width - 2 * pad:
                ok = False
                break
            cs = spread_1d([dx[i] for i in mem], [wd[i] for i in mem], axbb.x0 + pad, axbb.x1 - pad, gap)
            lay.update(dict(zip(mem, cs)))
            disp += float(np.sum(np.abs(cs - np.array([dx[i] for i in mem]))))
        if not ok:
            continue
        cost = disp + (nrow - 1) * 0.12 * axbb.width
        if best is None or cost < best[0]:
            best = (cost, nrow, lay)
    if best is None:                                                  # labels wider than the axes: one per row
        best = (0, len(items), {i: dx[i] for i in idx})
    _, rows, lay = best

    ymax = float(pts[:, 1].max())
    row_h = hb * 1.25
    free = H - hl - rows * row_h - 0.06 * H
    ytop = ymax * H / max(free, 0.35 * H)
    ax.set_ylim(0, ytop)
    to_disp = ax.transData.transform
    to_data = ax.transData.inverted().transform
    base = to_disp((0, ymax))[1] + 0.03 * H

    for rrow in range(rows):
        yb = base + rrow * row_h
        for i in idx[rrow::rows]:
            xd, yd = to_data((lay[i], yb))
            texts[i].set_position((xd, yd))
            y_from = items[i]["y"] + 0.012 * ytop
            ax.plot([items[i]["x"], xd], [y_from, yd], color=items[i]["color"], lw=0.55, zorder=5,
                    solid_capstyle="round")
    return texts


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="P(F) histograms with peaks annotated by c/a and cell angle")
    ap.add_argument("csv")
    ap.add_argument("--repeats", type=int, nargs=3, metavar=("N1", "N2", "N3"),
                    help="supercell repeats along the three lattice vectors for ALL structures (default: automatic)")
    ap.add_argument("--c-axis", choices=["longest", "unique"], default="longest",
                    help="c = longest normalised axis (default) or the axis most different from the other two")
    ap.add_argument("--bins", default="auto")
    ap.add_argument("--exclude-imaginary", action="store_true")
    ap.add_argument("--max-components", type=int, default=3)
    ap.add_argument("--components", action="append", default=[], metavar="[PHASE=]K")
    ap.add_argument("--min-peak", type=float, default=0.15)
    ap.add_argument("--bootstrap", type=int, default=100)
    ap.add_argument("--reference", choices=["lowest", "main"], default="lowest")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--label", action="append", default=[], metavar="PHASE=TEXT")
    ap.add_argument("--font", choices=["serif", "sans"], default="serif")
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"], choices=["pdf", "png", "svg", "eps"])
    ap.add_argument("--legend-loc", default="auto")
    ap.add_argument("--no-components", dest="show_components", action="store_false")
    ap.add_argument("--no-info", action="store_true")
    ap.add_argument("--out-dir")
    a = ap.parse_args()

    forced = {}
    for item in a.components:
        k_, _, v_ = item.rpartition("=")
        forced[k_ or "*"] = int(v_)
    labels = dict(item.split("=", 1) for item in a.label if "=" in item)

    pfe.set_style(a.font)
    csv_path = Path(a.csv)
    data = load_rows(csv_path, a.exclude_imaginary, a.repeats, a.c_axis)
    out_dir = Path(a.out_dir) if a.out_dir else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    T = data[0]["T"]

    # per-structure table
    f_struct = out_dir / f"{csv_path.stem}_cell_params.csv"
    with open(f_struct, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "phase", "F_eV_atom", "repeats", "L1_n_A", "L2_n_A", "L3_n_A", "c_axis",
                    "a_n_A", "c_n_A", "c_over_a", "theta_deg", "theta_between"])
        for d in data:
            w.writerow([d["file"], d["phase"], f"{d['F']:.8f}", "x".join(map(str, d["rep"])),
                        *[f"{v:.5f}" for v in d["Ln"]], "L" + str(d["c_axis"] + 1), f"{d['a_n']:.5f}",
                        f"{d['c_n']:.5f}", f"{d['ca']:.5f}", f"{d['theta']:.3f}",
                        "L" + ANGLE_NAMES[d["theta_idx"]][1].replace("a", "1").replace("b", "2").replace("c", "3")
                        .replace(",", "-L")])

    comps = list(OrderedDict.fromkeys(d["comp"] for d in data))
    extra = OrderedDict()
    for comp in comps:
        sub = [d for d in data if d["comp"] == comp]
        Fmin = min(d["F"] for d in sub)
        x0 = np.array([1000 * (d["F"] - Fmin) for d in sub])
        ph_of = np.array([d["phase"] for d in sub])
        phases = [p for p in list(pfe.PHASE_STYLE) + sorted(set(ph_of) - set(pfe.PHASE_STYLE)) if p in set(ph_of)]
        styles = {p: pfe.style_for(p, extra, labels) for p in phases}
        tag = csv_path.stem + (f"_{comp}" if len(comps) > 1 and comp else "")

        bins = int(a.bins) if str(a.bins).isdigit() else "auto"
        edges = np.histogram_bin_edges(x0, bins=bins)
        if len(edges) < 2 or edges[-1] == edges[0]:
            edges = np.linspace(x0.min() - 0.5, x0.max() + 0.5, 11)
        rng = np.random.default_rng(a.seed)
        fits = {p: pfe.fit_peaks(x0[ph_of == p], rng, kmax=a.max_components, k_fixed=forced.get(p, forced.get("*")),
                                 min_peak=a.min_peak, n_boot=a.bootstrap) for p in phases}
        cands = [(p, pk) for p in phases for pk in fits[p]["peaks"]]
        if a.reference == "main":
            cands = [(p, max(fits[p]["peaks"], key=lambda q: q[2])) for p in phases]
        ref, (mu_ref, err_ref, _) = min(cands, key=lambda t: t[1][0])
        for p in phases:
            fits[p]["peak_info"] = [dict(F=F0, err=e, height=hgt, is_ref=(p == ref and F0 == mu_ref), dF=F0 - mu_ref,
                                         dF_err=0.0 if (p == ref and F0 == mu_ref) else float(np.hypot(e, err_ref)))
                                    for F0, e, hgt in fits[p]["peaks"]]

        # structures in each peak -> c/a and angle statistics
        bw = edges[1] - edges[0]
        lo, hi = edges[0] - 5 * bw, edges[-1] + 5 * bw
        peak_rows = []
        for p in phases:
            idx = np.nonzero(ph_of == p)[0]
            peaks_sorted, cut = basins(fits[p], lo, hi)
            for pk in fits[p]["peak_info"]:
                j = peaks_sorted.index(pk["F"])
                members = [sub[i] for i in idx if cut[j] < x0[i] <= cut[j + 1]]
                pk["n_in"] = len(members)
                pk["ca"] = summarise([m["ca"] for m in members])
                pk["a_n"] = summarise([m["a_n"] for m in members])
                pk["c_n"] = summarise([m["c_n"] for m in members])
                pk["theta"] = summarise([m["theta"] for m in members])
                names = [ANGLE_NAMES[m["theta_idx"]][0] for m in members]
                pk["theta_name"] = max(set(names), key=names.count) if names else r"\theta"
                peak_rows.append((p, pk))

        ref_lbl = styles[ref]["label"].split(" (")[0]
        to_axis = lambda v: v - mu_ref                                    # noqa: E731
        flabel = r"$F - F_\mathrm{peak}^{\,\mathrm{%s}}$ (meV/atom)" % ref_lbl.replace(" ", r"\ ")
        info = "" if a.no_info else ", ".join(b for b in (pfe.chem(comp) if comp else "",
                                                          f"$T$ = {float(T):g} K" if T else "",
                                                          f"$N$ = {len(sub)}") if b)
        P = dict(edges=edges, x0=x0, ph_of=ph_of, fits=fits, styles=styles, ref=ref, to_axis=to_axis,
                 show_components=a.show_components, phases=phases, absolute=False, flabel=flabel,
                 legend_loc=a.legend_loc, info=info)

        def items_for(kind):
            out = []
            for p, pk in peak_rows:
                f = fits[p]
                scale = f["n"] / len(x0) * bw
                ytop = scale * pfe.mix_pdf([pk["F"]], *_dp(f, bw))[0]
                short = styles[p]["label"].split(" (")[0]
                if kind == "ca":
                    txt = f"{short}: {spread_txt(*pk['ca'], 3)}"
                else:
                    txt = f"{short}: {spread_txt(*pk['theta'], 1)}°"
                out.append(dict(x=to_axis(pk["F"]), y=ytop, text=txt, color=styles[p]["color"]))
            return out

        def panel(ax, kind):
            title = (r"labels: $c_n/a_n$ per peak" if kind == "ca"
                     else r"labels: cell angle $\theta$ farthest from 90° per peak")
            annotated_panel(ax, P, items_for(kind), title)

        written = []
        for kind, name in (("ca", "PF_ca"), ("angle", "PF_angle")):
            fig, ax = plt.subplots(figsize=(pfe.WIDTH_SINGLE, 0.85 * pfe.WIDTH_SINGLE))
            fig.canvas.draw()
            panel(ax, kind)
            pfe.info_text(ax, P)
            written += pfe.save(fig, out_dir / f"{tag}_{name}", a.formats)

        fig, (axa, axb) = plt.subplots(1, 2, figsize=(pfe.WIDTH_DOUBLE, 0.42 * pfe.WIDTH_DOUBLE),
                                       gridspec_kw=dict(wspace=0.25))
        fig.canvas.draw()
        panel(axa, "ca")
        panel(axb, "angle")
        pfe.info_text(axa, P)
        for lab, axx in (("(a)", axa), ("(b)", axb)):
            axx.text(-0.15, 1.0, lab, transform=axx.transAxes, ha="left", va="bottom", fontsize=10,
                     fontweight="bold" if a.font == "sans" else "normal")
        written += pfe.save(fig, out_dir / f"{tag}_PF_cell_combined", a.formats)

        # per-peak table
        f_peak = out_dir / f"{tag}_peak_cell.csv"
        with open(f_peak, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["phase", "peak_F_eV_atom", "dF_meV_atom", "dF_err_meV_atom", "reference",
                        "structures_in_peak", "c_over_a_mean", "c_over_a_std", "a_n_mean_A", "a_n_std_A",
                        "c_n_mean_A", "c_n_std_A", "theta_mean_deg", "theta_std_deg", "theta_is"])
            for p, pk in sorted(peak_rows, key=lambda t: t[1]["F"]):
                w.writerow([styles[p]["label"], f"{Fmin + pk['F'] / 1000:.8f}", f"{pk['dF']:.4f}",
                            f"{pk['dF_err']:.4f}", "yes" if pk["is_ref"] else "", pk["n_in"],
                            *[f"{v:.5f}" for v in (*pk["ca"], *pk["a_n"], *pk["c_n"])],
                            *[f"{v:.3f}" for v in pk["theta"]], pk["theta_name"].strip("\\")])
        written += [str(f_peak), str(f_struct)]

        print(f"{comp or 'all'}: {len(sub)} structures; reference peak: {styles[ref]['label']} ({a.reference})")
        print(f"  {'phase':<13} {'ΔF (meV/atom)':>15} {'n':>6} {'c/a':>16} {'a_n (Å)':>15} {'c_n (Å)':>15} {'θ (°)':>13}")
        for p, pk in sorted(peak_rows, key=lambda t: t[1]["F"]):
            print(f"  {styles[p]['label']:<13} {pk['dF']:+9.2f} ± {pk['dF_err']:<4.2f} {pk['n_in']:>6} "
                  f"{pk['ca'][0]:>9.4f} ± {pk['ca'][1]:<5.4f} {pk['a_n'][0]:>7.4f} ± {pk['a_n'][1]:<5.4f} "
                  f"{pk['c_n'][0]:>7.4f} ± {pk['c_n'][1]:<5.4f} {pk['theta'][0]:>6.2f} ± {pk['theta'][1]:<4.2f}")
        reps = OrderedDict()
        for d in sub:
            reps[d["rep"]] = reps.get(d["rep"], 0) + 1
        print("  supercell repeats used: " + ", ".join(f"{'x'.join(map(str, k))} ({v})" for k, v in reps.items()))
        for fn in written:
            print(f"  wrote {fn}")
        print("  caption note: labels = mean over the structures in each peak, standard deviation across them in "
              "brackets; c/a = c_n/a_n with c the unique axis of each cell; θ = cell angle farthest from 90°.")


if __name__ == "__main__":
    main()
