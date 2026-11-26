#!/usr/bin/env python3
"""
plot_cell_distributions.py - distributions of the cell shape over all structures, one figure each:

  P(c_n/a_n) vs c_n/a_n     (tetragonality)
  P(θ)       vs θ           (cell angle farthest from 90°, folded to >= 90°; e.g. the monoclinic β)

Bins: one common width, set by the narrowest phase (Freedman-Diaconis), so every peak is resolved. If the
values fall into groups far apart (e.g. θ ~ 90° and ~ 99°), the x axis is broken into one panel per group.

One histogram per phase (all phases together sum to 1), each fitted with a Gaussian mixture (number of
Gaussians by BIC, as in plot_free_energy.py). Every maximum of a fitted curve is a peak; the legend gives
the peak positions with bootstrap uncertainties in brackets, e.g. t: 1.0257(1).

c_n/a_n and θ are computed per structure exactly as in plot_cell_annotations.py: lattice vector lengths are
divided by the supercell repeats; c is the longest normalised axis (x, y or z, decided per structure; or the
axis most different from the other two with --c-axis unique), a_n the mean of the other two. Cell data come from the CSV of batch_free_energy.py.
Needs plot_free_energy.py and plot_cell_annotations.py in the same folder.

Usage
    python plot_cell_distributions.py Ce_33.33/free_energy_300K.csv
    python plot_cell_distributions.py free_energy_300K.csv --angle-range 89.5 100 --no-break
    python plot_cell_distributions.py free_energy_300K.csv --repeats 2 2 2 --exclude-imaginary

Output (next to the CSV, or --out-dir)
    <stem>_P_ca.pdf/.png        single column
    <stem>_P_angle.pdf/.png     single column
    <stem>_cell_peaks.csv       every peak of both distributions per phase, plus mean and std per phase
"""
import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import plot_cell_annotations as pca        # also checks that plot_free_energy.py is recent enough
except ImportError as e:
    sys.exit(f"plot_cell_annotations.py and plot_free_energy.py must be in the same folder as this script ({e})")
pfe = pca.pfe

import matplotlib  # noqa: E402
import matplotlib.ticker  # noqa: E402
_dp = getattr(pfe, "display_params", lambda f, bw: (f["w"], f["mu"], f["sigma"]))  # drawing parameters
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402


def fmt_val(v, e, max_dec=6):
    """1.0257(1): value with one significant digit of uncertainty in brackets"""
    if not np.isfinite(e) or e <= 0:
        return f"{v:.4f}"
    dec = int(np.clip(-np.floor(np.log10(e)), 0, max_dec))
    d = int(round(e * 10 ** dec))
    if d >= 10 and dec > 0:
        dec -= 1
        d = int(round(e * 10 ** dec))
    return f"{v:.{dec}f}({max(d, 1)})"


def fd_width(x):
    """Freedman-Diaconis bin width of one phase"""
    x = np.asarray(x, float)
    if len(x) < 5:
        return np.inf
    q75, q25 = np.percentile(x, [75, 25])
    spread = (q75 - q25) or np.std(x)
    return 2 * spread / len(x) ** (1 / 3) if spread > 0 else np.inf


def plan_axes(values, lo, hi, bin_width=None, max_bins=120, break_ok=True):
    """
    Common bin width (narrowest phase, Freedman-Diaconis) and x segments. If the values form groups so far
    apart that one axis cannot resolve the narrow peaks, the axis is broken into one segment per group.
    Returns (bin width, [(seg_lo, seg_hi), ...]).
    """
    allv = np.sort(np.concatenate([v[(v >= lo) & (v <= hi)] for v in values.values()]))
    span = hi - lo
    bw_fd = min([fd_width(v[(v >= lo) & (v <= hi)]) for v in values.values()] + [np.inf])
    if not np.isfinite(bw_fd) or bw_fd <= 0:
        bw_fd = span / 40 if span > 0 else 1.0
    bw = bin_width or bw_fd
    segs = [(lo, hi)]
    if break_ok and not bin_width and span / bw_fd > 300:            # one axis would hide the narrow peaks
        gaps = np.diff(allv)
        cut = np.nonzero(gaps > max(30 * bw_fd, 0.08 * span))[0]
        if len(cut):
            bounds = np.r_[0, cut + 1, len(allv)]
            segs = [(allv[b0], allv[b1 - 1]) for b0, b1 in zip(bounds[:-1], bounds[1:])]
    if not bin_width:
        widest = max(b - a for a, b in segs) or span or 1.0
        bw = max(bw_fd, widest / max_bins)
    out = []
    for a_, b_ in segs:                                              # margins for the fitted tails
        m = max(6 * bw, 0.12 * (b_ - a_))
        out.append((a_ - m, b_ + m))
    return bw, out


def draw_distribution(fig, values, phases, styles, fits, bw, segs, xlabel, ylabel, ref_line, legend_title,
                      show_components, legend_loc, info):
    """one axis, or a broken x axis with one panel per segment (shared y, break marks between panels)"""
    n_all = sum(len(v) for v in values.values())
    widths = [max(b - a, 1e-12) for a, b in segs]
    ratios = [max(w / sum(widths), 0.3) for w in widths]
    axes = fig.subplots(1, len(segs), sharey=True, gridspec_kw=dict(width_ratios=ratios, wspace=0.06))
    axes = np.atleast_1d(axes)
    order = sorted(phases, key=lambda p: fits[p]["peaks"][0][0])
    handles, labels = [], []
    seg_pts = [[] for _ in segs]
    allv = np.concatenate(list(values.values()))
    seg_edges = []                                                   # bins centred on the data (median)
    for a_, b_ in segs:
        inseg = allv[(allv >= a_) & (allv <= b_)]
        anchor = float(np.median(inseg)) if len(inseg) else 0.5 * (a_ + b_)
        k0 = int(np.floor((a_ - anchor + bw / 2) / bw))
        k1 = int(np.ceil((b_ - anchor + bw / 2) / bw))
        seg_edges.append(anchor - bw / 2 + bw * np.arange(k0, k1 + 1))
    ymax = 0.0
    for p in order:
        st, f, x = styles[p], fits[p], values[p]
        scale = len(x) / n_all * bw
        for si, (ax, (a_, b_)) in enumerate(zip(axes, segs)):
            edges = seg_edges[si]
            if len(edges) < 2:
                continue
            h = np.histogram(x, bins=edges)[0] / n_all
            fine = np.linspace(a_, b_, 1500)
            g = scale * pfe.mix_pdf(fine, *_dp(f, bw))
            if h.max() == 0 and g.max() < 1e-6 * max(ymax, 1e-12):
                continue
            ax.stairs(h, edges, fill=True, color=st["color"], alpha=0.18, lw=0, zorder=1)
            ax.stairs(h, edges, color=st["color"], lw=0.5, alpha=0.9, zorder=2)
            if show_components and f["k"] > 1:
                for wj, mj, sj in zip(*_dp(f, bw)):
                    ax.plot(fine, scale * pfe.mix_pdf(fine, [wj], [mj], [sj]), color=st["color"], lw=0.6,
                            ls=(0, (1, 1.2)), zorder=3)
            ax.plot(fine, g, color=st["color"], lw=1.3, zorder=4)
            for x0, _, _ in f["peaks"]:
                if a_ <= x0 <= b_:
                    ax.vlines(x0, 0, scale * pfe.mix_pdf([x0], *_dp(f, bw))[0], color=st["color"],
                              lw=0.7, ls=(0, (3, 2)), zorder=3)
            ymax = max(ymax, h.max(), g.max())
            for k in range(len(h)):
                seg_pts[si].append(np.column_stack([np.linspace(edges[k], edges[k + 1], 5), np.full(5, h[k])]))
            seg_pts[si].append(np.column_stack([fine, g]))
        handles.append((Patch(facecolor=st["color"], alpha=0.25, lw=0), Line2D([], [], color=st["color"], lw=1.3)))
        labels.append(f"{st['label']}: " + ", ".join(fmt_val(x0, e) for x0, e, _ in sorted(f["peaks"])))

    for si, (ax, (a_, b_)) in enumerate(zip(axes, segs)):
        ax.set_xlim(a_, b_)
        if ref_line is not None and a_ <= ref_line <= b_:
            ax.axvline(ref_line, color=pfe.INK2, lw=0.6, ls=(0, (1, 1.5)), zorder=0)
        pfe.finish_axes(ax)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(
            nbins=max(2, int(5 * ratios[si] + 1)), prune=("both" if len(segs) > 1 else None)))
        ax.ticklabel_format(axis="x", useOffset=False)               # no '+1' offset labels
    axes[0].set_ylim(0, ymax * 1.05)
    axes[0].set_ylabel(ylabel)
    if len(axes) > 1:                                                # broken-axis cosmetics
        for i, ax in enumerate(axes):
            if i > 0:
                ax.spines["left"].set_visible(False)
                ax.tick_params(which="both", left=False, labelleft=False)
            if i < len(axes) - 1:
                ax.spines["right"].set_visible(False)
                ax.tick_params(which="both", right=False)
            d = 0.012
            kw = dict(transform=ax.transAxes, color="k", lw=0.6, clip_on=False)
            if i < len(axes) - 1:
                for yy in (0, 1):
                    ax.plot((1 - d * 2, 1 + d * 2), (yy - d * 4, yy + d * 4), **kw)
            if i > 0:
                for yy in (0, 1):
                    ax.plot((-d * 2, d * 2), (yy - d * 4, yy + d * 4), **kw)
        fig.supxlabel(xlabel, fontsize=plt.rcParams["axes.labelsize"], y=0.0)
    else:
        axes[0].set_xlabel(xlabel)

    # legend in the panel whose data leave the most room at the top
    tallest = [(np.vstack(sp)[:, 1].max() if sp else 0.0) for sp in seg_pts]
    li = int(np.argmin(tallest)) if len(axes) > 1 else 0
    lax = axes[li]
    pts = np.vstack(seg_pts[li]) if seg_pts[li] else np.zeros((1, 2))
    pfe.place_legend(lax, lambda loc, nc: lax.legend(handles, labels, loc=loc, ncol=1, alignment="left",
                                                     title=legend_title, columnspacing=1.0),
                     pts, filled=True, loc=legend_loc)
    if info:
        axes[0].text(0.0, 1.07 if len(axes) > 1 else 1.015, info, transform=axes[0].transAxes, ha="left",
                     va="bottom", fontsize=8, color=pfe.INK2)                 # clear of the break marks
    return axes


def main():
    ap = argparse.ArgumentParser(description="P(c/a) and P(angle) distributions per phase")
    ap.add_argument("csv")
    ap.add_argument("--repeats", type=int, nargs=3, metavar=("N1", "N2", "N3"),
                    help="supercell repeats for ALL structures (default: automatic)")
    ap.add_argument("--c-axis", choices=["longest", "unique"], default="longest",
                    help="c = longest normalised axis (default) or the axis most different from the other two")
    ap.add_argument("--bin-width-ca", type=float, help="bin width for c/a (default: from the narrowest phase)")
    ap.add_argument("--bin-width-angle", type=float, help="bin width for θ in degrees (default: automatic)")
    ap.add_argument("--no-break", action="store_true", help="never break the x axis")
    ap.add_argument("--ca-range", type=float, nargs=2, metavar=("MIN", "MAX"), help="x range of P(c/a)")
    ap.add_argument("--angle-range", type=float, nargs=2, metavar=("MIN", "MAX"), help="x range of P(θ) (degrees)")
    ap.add_argument("--exclude-imaginary", action="store_true")
    ap.add_argument("--max-components", type=int, default=3)
    ap.add_argument("--min-peak", type=float, default=0.15)
    ap.add_argument("--bootstrap", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--label", action="append", default=[], metavar="PHASE=TEXT")
    ap.add_argument("--font", choices=["serif", "sans"], default="serif")
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"], choices=["pdf", "png", "svg", "eps"])
    ap.add_argument("--legend-loc", default="auto")
    ap.add_argument("--no-components", dest="show_components", action="store_false")
    ap.add_argument("--no-info", action="store_true")
    ap.add_argument("--out-dir")
    a = ap.parse_args()
    labels = dict(item.split("=", 1) for item in a.label if "=" in item)

    pfe.set_style(a.font)
    csv_path = Path(a.csv)
    data = pca.load_rows(csv_path, a.exclude_imaginary, a.repeats)   # works with old and new plot_cell_annotations
    for d in data:                                                    # c-axis rule applied here, independent of it
        Ln = np.asarray(d["Ln"], float)
        if a.c_axis == "unique":
            ic = int(np.argmax([abs(Ln[i] - np.mean(np.delete(Ln, i))) for i in range(3)]))
        else:
            ic = int(np.argmax(Ln))
        d["c_axis"], d["c_n"], d["a_n"] = ic, float(Ln[ic]), float(np.mean(np.delete(Ln, ic)))
        d["ca"] = d["c_n"] / d["a_n"]
    out_dir = Path(a.out_dir) if a.out_dir else csv_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    T = data[0]["T"]

    comps = list(OrderedDict.fromkeys(d["comp"] for d in data))
    extra = OrderedDict()
    for comp in comps:
        sub = [d for d in data if d["comp"] == comp]
        tag = csv_path.stem + (f"_{comp}" if len(comps) > 1 and comp else "")
        phases_all = [d["phase"] for d in sub]
        phases = [p for p in list(pfe.PHASE_STYLE) + sorted(set(phases_all) - set(pfe.PHASE_STYLE))
                  if p in set(phases_all)]
        styles = {p: pfe.style_for(p, extra, labels) for p in phases}
        info = "" if a.no_info else ", ".join(b for b in (pfe.chem(comp) if comp else "",
                                                          f"$T$ = {float(T):g} K" if T else "",
                                                          f"$N$ = {len(sub)}") if b)
        rng = np.random.default_rng(a.seed)
        written, table = [], []

        for key, name, xlabel, ylabel, ref_line, title, rng_opt, bwidth in (
                ("ca", "P_ca", r"$c_n/a_n$", r"$P(c_n/a_n)$", 1.0, r"$c_n/a_n$ peaks", a.ca_range, a.bin_width_ca),
                ("theta", "P_angle", r"$\theta$ (deg)", r"$P(\theta)$", 90.0,
                 r"$\theta$ peaks (deg)", a.angle_range, a.bin_width_angle)):
            values = {p: np.array([d[key] for d in sub if d["phase"] == p]) for p in phases}
            allv = np.concatenate(list(values.values()))
            lo, hi = (rng_opt if rng_opt else (allv.min(), allv.max()))
            if hi <= lo:
                lo, hi = lo - 0.5, hi + 0.5
            bw, segs = plan_axes(values, lo, hi, bwidth, break_ok=not a.no_break)
            fits = {p: pfe.fit_peaks(values[p], rng, kmax=a.max_components, min_peak=a.min_peak,
                                     n_boot=a.bootstrap) for p in phases}

            fig = plt.figure(figsize=(pfe.WIDTH_SINGLE, 0.78 * pfe.WIDTH_SINGLE))
            draw_distribution(fig, values, phases, styles, fits, bw, segs, xlabel, ylabel, ref_line, title,
                              a.show_components, a.legend_loc, info)
            written += pfe.save(fig, out_dir / f"{tag}_{name}", a.formats)
            if len(segs) > 1:
                print(f"  {name}: x axis broken into {len(segs)} panels "
                      + ", ".join(f"[{s0:.4g}, {s1:.4g}]" for s0, s1 in segs) + " (--no-break to avoid)")

            for p in phases:
                v = values[p]
                for i, (x0, e, hgt) in enumerate(sorted(fits[p]["peaks"]), 1):
                    table.append(dict(quantity="c_n/a_n" if key == "ca" else "theta_deg", phase=styles[p]["label"],
                                      peak=i, peak_value=f"{x0:.6f}", peak_err=f"{e:.6f}", n_structures=len(v),
                                      n_gaussians=fits[p]["k"], mean=f"{v.mean():.6f}",
                                      std=f"{v.std(ddof=1) if len(v) > 1 else 0:.6f}",
                                      min=f"{v.min():.6f}", max=f"{v.max():.6f}"))

        f_csv = out_dir / f"{tag}_cell_peaks.csv"
        with open(f_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(table[0]))
            w.writeheader()
            w.writerows(table)
        written.append(str(f_csv))

        print(f"{comp or 'all'}: {len(sub)} structures")
        for q in ("c_n/a_n", "theta_deg"):
            print(f"  {q}")
            for row in [t for t in table if t["quantity"] == q]:
                print(f"    {row['phase']:<13} peak {row['peak']}: {float(row['peak_value']):.5f} ± "
                      f"{float(row['peak_err']):.5f}   (phase mean {float(row['mean']):.5f}, std {float(row['std']):.5f}, "
                      f"n = {row['n_structures']})")
        for fn in written:
            print(f"  wrote {fn}")
        print("  caption note: histograms per phase sum to 1 over all phases; curves = Gaussian-mixture fits; "
              "peak positions with bootstrap uncertainties in brackets; dotted line = c/a = 1 or θ = 90°.")


if __name__ == "__main__":
    main()
