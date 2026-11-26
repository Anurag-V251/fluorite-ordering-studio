#!/usr/bin/env python3
"""
plot_harmonic_compare.py - compare harmonic free energies of several structures: F(T), dF(T) against a chosen
reference, S / Cv / U_vib / F_vib panels, dF at chosen temperatures, phonon DOS.  Publication figures
(PNG / PDF / SVG, matplotlib), a comparison CSV and one self-contained interactive HTML.

Input: results of single_runs.py (its --out-dir, or several), or thermo.csv files directly (label = folder name)
    python plot_harmonic_compare.py single_runs/
    python plot_harmonic_compare.py single_runs/ --reference POSCAR_m --T-window 0 1200 --formats pdf svg
    python plot_harmonic_compare.py runA/ runB/ --include m_relaxed t_relaxed --at-T 300 900

Reference of dF (--reference):
    lowest   the lowest F_total at each temperature among the unflagged structures (default)
    LABEL    one structure: dF = F(structure) - F(LABEL) at every T
    none     no dF panel
Energies are per atom (eV/atom; dF in meV/atom).  Structures with imaginary modes or an unconverged relaxation
are drawn dashed and marked in the legend and the CSV: their harmonic F is not reliable.

Method note printed on every figure: fixed-cell harmonic approximation (not quasi-harmonic), so no thermal
expansion; F_total(T) = E_static + F_vib(T).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

PALETTE = ["#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"]
PALETTE_DARK = ["#338ced", "#e76624", "#00a888", "#d28e05", "#e24d8b", "#0aac23", "#9887ef", "#eb6263"]
OVERFLOW = "#8a8a8a"
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]
METHOD = "Fixed-cell harmonic approximation (not quasi-harmonic): no thermal expansion"
QUANT = {   # key: (column, label, unit, scale)
    "F": ("F_total_eV_atom", "F$_\\mathrm{total}$", "eV/atom", 1.0),
    "F_vib": ("F_vib_eV_atom", "F$_\\mathrm{vib}$", "meV/atom", 1000.0),
    "U_vib": ("U_vib_eV_atom", "U$_\\mathrm{vib}$ (incl. ZPE)", "meV/atom", 1000.0),
    "S": ("entropy_J_molK_atom", "S$_\\mathrm{vib}$", "J mol$^{-1}$ K$^{-1}$ per atom", 1.0),
    "Cv": ("Cv_J_molK_atom", "C$_V$", "J mol$^{-1}$ K$^{-1}$ per atom", 1.0),
}


def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


# --------------------------------------------------------------------------- loading
def load(sources):
    """[{label, thermo: {col: [..]}, flags, dos, meta}] from single_runs folders or thermo.csv files"""
    out = []
    for s in sources:
        p = Path(s)
        items = []
        if p.is_file():
            items.append((p.parent.name, p))
        elif (p / "thermo.csv").is_file():
            items.append((p.name, p / "thermo.csv"))
        elif p.is_dir():
            summ = {r["label"]: r for r in read_csv(p / "summary.csv")} if (p / "summary.csv").is_file() else {}
            order = list(summ) or sorted(d.name for d in p.iterdir() if (d / "thermo.csv").is_file())
            for lab in order:
                if (p / lab / "thermo.csv").is_file():
                    items.append((lab, p / lab / "thermo.csv"))
        else:
            sys.exit(f"{s}: not a single_runs folder or thermo.csv")
        for lab, tp in items:
            rows = [r for r in read_csv(tp) if fnum(r.get("T_K")) is not None]
            if not rows:
                continue
            rows.sort(key=lambda r: float(r["T_K"]))
            th = {k: [fnum(r.get(k)) for r in rows] for k in rows[0]}
            res = {}
            rj = tp.parent / "result.json"
            if rj.is_file():
                try:
                    res = json.loads(rj.read_text())
                except ValueError:
                    res = {}
            dos = None
            dp = tp.parent / "dos.csv"
            if dp.is_file():
                dr = read_csv(dp)
                dos = ([fnum(r["frequency_THz"]) for r in dr], [fnum(r["dos_states_per_THz_atom"]) for r in dr])
            name = lab
            k = 2
            while name in {o["label"] for o in out}:                 # same label from two sources
                name = f"{lab}_{k}"
                k += 1
            out.append({"label": name, "thermo": th, "flags": list(res.get("flags") or []), "dos": dos,
                        "phase": res.get("phase") or "", "formula": res.get("formula") or "",
                        "notes": res.get("notes") or [], "source": str(tp)})
    return out


def common_T(series, window=None):
    """temperatures present in every structure (within 1e-6 K), optionally limited to a window"""
    sets = [set(round(t, 6) for t in s["thermo"]["T_K"]) for s in series]
    T = sorted(set.intersection(*sets)) if sets else []
    if window:
        T = [t for t in T if window[0] - 1e-9 <= t <= window[1] + 1e-9]
    return T


def values_at(s, col, T):
    idx = {round(t, 6): i for i, t in enumerate(s["thermo"]["T_K"])}
    return [s["thermo"][col][idx[t]] if t in idx and col in s["thermo"] else None for t in T]


def reference_curve(series, T, ref):
    if ref == "none":
        return None, "none"
    F = {s["label"]: values_at(s, "F_total_eV_atom", T) for s in series}
    if ref == "lowest":
        pool = [s["label"] for s in series if not s["flags"]] or [s["label"] for s in series]
        curve = [min((F[l][i] for l in pool if F[l][i] is not None), default=None) for i in range(len(T))]
        return curve, "lowest F at each T" + ("" if len(pool) == len(series) else " (unflagged structures)")
    if ref not in F:
        sys.exit(f"--reference {ref}: no such structure; choose one of: lowest, none, " + ", ".join(F))
    return F[ref], ref


# --------------------------------------------------------------------------- figures
def style(plt, fs):
    plt.rcParams.update({
        "font.size": fs, "axes.labelsize": fs, "axes.titlesize": fs, "legend.fontsize": fs - 1,
        "xtick.labelsize": fs - 1, "ytick.labelsize": fs - 1, "axes.linewidth": 0.8,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True, "lines.linewidth": 1.6,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular"})


def colour(i):
    return PALETTE[i] if i < len(PALETTE) else OVERFLOW


def legend_label(s):
    tags = []
    if "imaginary_modes" in s["flags"]:
        tags.append("imag. modes")
    if "relax_not_converged" in s["flags"]:
        tags.append("relax not conv.")
    if "not_relaxed" in s["flags"]:
        tags.append("not relaxed")
    return s["label"] + (f"  ({', '.join(tags)})" if tags else "")


def footnote(fig, series, extra=""):
    """method note under everything already drawn (axis labels included), so it never overlaps them"""
    txt = METHOD + "." + (" Dashed: flagged structure, harmonic F not reliable." if any(s["flags"] for s in series) else "")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = [ax.get_tightbbox(r) for ax in fig.axes]
    y0 = min(b.y0 for b in boxes) / fig.bbox.height
    x0 = min(b.x0 for b in boxes) / fig.bbox.width
    fig.text(x0, y0 - 0.012, txt + extra, ha="left", va="top", fontsize=plt_fs(fig) - 2, color="#555555", wrap=True)


def plt_fs(_fig=None):
    import matplotlib.pyplot as plt
    return plt.rcParams["font.size"]


def direct_labels(ax, ends):
    """labels at the line ends, pushed apart vertically so equal end values stay readable"""
    if not ends:
        return
    fig = ax.figure
    fig.canvas.draw()
    to_px = ax.transData.transform
    gap = plt_fs() * 1.25 * fig.dpi / 72.0
    pts = sorted(((to_px((x, y))[1], x, y, lab) for x, y, lab in ends))
    clusters = []                                  # runs of labels closer than one line: spread around their mean
    for p in pts:
        if clusters and p[0] - clusters[-1][-1][0] < gap:
            clusters[-1].append(p)
        else:
            clusters.append([p])
    placed = []
    for c in clusters:
        mid = sum(p[0] for p in c) / len(c)
        for k, (py, x, y, lab) in enumerate(c):
            placed.append((mid + (k - (len(c) - 1) / 2) * gap, x, y, lab))
    for py, x, y, lab in placed:
        dy = (py - to_px((x, y))[1]) * 72.0 / fig.dpi
        ax.annotate(lab, (x, y), xytext=(5, dy), textcoords="offset points", va="center",
                    fontsize=plt_fs() - 2, color="#333333")


def draw_lines(ax, series, T, ys, direct):
    ends = []
    for i, (s, y) in enumerate(zip(series, ys)):
        pts = [(t, v) for t, v in zip(T, y) if v is not None]
        if not pts:
            continue
        x, v = zip(*pts)
        ls = "--" if s["flags"] else "-"
        mk = MARKERS[i % len(MARKERS)]
        every = max(1, len(x) // 8)
        ax.plot(x, v, ls, color=colour(i), marker=mk, markersize=4, markevery=every, label=legend_label(s),
                markeredgewidth=0.6, markeredgecolor="white")
        ends.append((x[-1], v[-1], s["label"]))
    return ends


def save(fig, out, stem, formats, dpi):
    paths = []
    for f in formats:
        p = out / f"{stem}.{f}"
        fig.savefig(p, dpi=dpi)
        paths.append(p)
    return paths


def fig_F_vs_T(plt, series, T, ref_curve, ref_name, a, out):
    has_ref = ref_curve is not None
    fig, axes = plt.subplots(2 if has_ref else 1, 1, sharex=True, figsize=(a.width, a.height * (1.55 if has_ref else 1)),
                             gridspec_kw={"height_ratios": [1, 1]} if has_ref else None, squeeze=False)
    axes = axes[:, 0]
    direct = len(series) <= 4
    draw_lines(axes[0], series, T, [values_at(s, "F_total_eV_atom", T) for s in series], False)
    axes[0].set_ylabel("F$_\\mathrm{total}$ (eV/atom)")
    if has_ref:
        dF = [[None if (v is None or r is None) else 1000 * (v - r) for v, r in zip(values_at(s, "F_total_eV_atom", T), ref_curve)]
              for s in series]
        ends = draw_lines(axes[1], series, T, dF, direct)
        axes[1].axhline(0, color="#999999", lw=0.8, zorder=0)
        axes[1].set_ylabel(f"$\\Delta$F (meV/atom)\nvs {ref_name}")
    axes[-1].set_xlabel("T (K)")
    if direct and has_ref:
        axes[-1].set_xlim(T[0], T[-1] + 0.18 * (T[-1] - T[0] or 1))   # room for the direct labels
    if len(series) > 1:
        axes[0].legend(loc="best", handlelength=2.4)
    fig.align_ylabels(axes)
    if direct and has_ref:
        direct_labels(axes[1], ends)
    footnote(fig, series)
    paths = save(fig, out, "F_vs_T", a.formats, a.dpi)
    plt.close(fig)
    return paths


def fig_panels(plt, series, T, a, out):
    keys = ["F_vib", "U_vib", "S", "Cv"]
    fig, axes = plt.subplots(2, 2, sharex=True, figsize=(a.width * 1.9, a.height * 1.6))
    for ax, k in zip(axes.flat, keys):
        col, lab, unit, sc = QUANT[k]
        draw_lines(ax, series, T, [[None if v is None else sc * v for v in values_at(s, col, T)] for s in series], False)
        ax.set_ylabel(f"{lab} ({unit})")
    for ax in axes[-1]:
        ax.set_xlabel("T (K)")
    if len(series) > 1:
        axes[0, 0].legend(loc="best", handlelength=2.4)
    fig.align_ylabels(axes[:, 0])
    footnote(fig, series)
    fig.tight_layout()
    paths = save(fig, out, "thermo_panels", a.formats, a.dpi)
    plt.close(fig)
    return paths


def fig_dF_at_T(plt, series, temps, ref, a, out):
    """dot chart: one row per structure, one marker per temperature (dF vs the reference at that T)"""
    data = []
    for t in temps:
        T = [round(t, 6)]
        rc, rname = reference_curve(series, T, "lowest" if ref == "none" else ref)
        data.append((t, rname, [None if (v is None or rc[0] is None) else 1000 * (v - rc[0])
                                for v in (values_at(s, "F_total_eV_atom", T)[0] for s in series)]))
    n = len(series)
    fig, ax = plt.subplots(figsize=(a.width, max(1.6, 0.32 * n + 0.9)))
    for j, (t, rname, d) in enumerate(data):
        idx = [i for i, v in enumerate(d) if v is not None]
        x = [d[i] for i in idx]
        dodge = (j - (len(data) - 1) / 2) * min(0.22, 0.7 / max(1, len(data)))   # one row per T inside a structure
        y = [i + dodge for i in idx]
        ax.plot(x, y, MARKERS[j % len(MARKERS)], color=PALETTE[j % len(PALETTE)] if len(data) > 1 else "#257de0",
                markersize=7, markeredgecolor="white", markeredgewidth=0.8, label=f"{t:g} K", linestyle="none")
        if len(data) == 1:
            for xi, yi in zip(x, y):
                ax.annotate(f"{xi + 0.0:.2f}".replace("-0.00", "0.00"), (xi, yi), xytext=(6, 0), textcoords="offset points", va="center",
                            fontsize=plt_fs(None) - 2, color="#333333")
    ax.set_yticks(range(n))
    ax.set_yticklabels([legend_label(s) for s in series])
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.invert_yaxis()
    ax.axvline(0, color="#999999", lw=0.8, zorder=0)
    rn = data[0][1] if data else ""
    ax.set_xlabel(f"$\\Delta$F (meV/atom) vs {rn}")
    if len(data) > 1:
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=min(5, len(data)), handletextpad=0.2,
                  columnspacing=1.0, borderaxespad=0.2)
    lo, hi = ax.get_xlim()
    ax.set_xlim(lo, hi + 0.15 * (hi - lo or 1))
    footnote(fig, series)
    paths = save(fig, out, "dF_at_T", a.formats, a.dpi)
    plt.close(fig)
    return paths, data


def fig_dos(plt, series, a, out):
    with_dos = [(i, s) for i, s in enumerate(series) if s["dos"]]
    if not with_dos:
        return []
    fig, ax = plt.subplots(figsize=(a.width, a.height))
    for i, s in with_dos:
        x, y = s["dos"]
        ax.plot(x, y, "--" if s["flags"] else "-", color=colour(i), lw=1.2, label=legend_label(s))
    ax.axvline(0, color="#999999", lw=0.8, zorder=0)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Phonon DOS (states THz$^{-1}$ atom$^{-1}$)")
    ax.set_ylim(bottom=0)
    if len(with_dos) > 1:
        ax.legend(loc="best")
    footnote(fig, series, " Negative frequencies = imaginary modes.")
    paths = save(fig, out, "phonon_dos", a.formats, a.dpi)
    plt.close(fig)
    return paths


# --------------------------------------------------------------------------- CSV / HTML
def write_csv(series, T, ref_curve, ref_name, out):
    p = out / "comparison.csv"
    with p.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "label", "F_total_eV_atom", "dF_meV_atom", "reference", "E_static_eV_atom", "F_vib_eV_atom",
                    "U_vib_eV_atom", "entropy_J_molK_atom", "Cv_J_molK_atom", "flags"])
        for s in series:
            cols = {c: values_at(s, c, T) for c in ("F_total_eV_atom", "E_static_eV_atom", "F_vib_eV_atom",
                                                     "U_vib_eV_atom", "entropy_J_molK_atom", "Cv_J_molK_atom")}
            for i, t in enumerate(T):
                F = cols["F_total_eV_atom"][i]
                dF = "" if ref_curve is None or F is None or ref_curve[i] is None else f"{1000 * (F - ref_curve[i]):.10g}"
                w.writerow([f"{t:g}", s["label"], *(("" if cols[c][i] is None else f"{cols[c][i]:.10g}")
                                                     for c in ("F_total_eV_atom",)), dF, ref_name,
                            *(("" if cols[c][i] is None else f"{cols[c][i]:.10g}") for c in
                              ("E_static_eV_atom", "F_vib_eV_atom", "U_vib_eV_atom", "entropy_J_molK_atom", "Cv_J_molK_atom")),
                            " ".join(s["flags"])])
    return p


HTML = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Harmonic free energies</title><style>
:root{--bg:#fafafa;--panel:#fff;--ink:#1d2330;--muted:#5f6b7a;--line:#dde2e8;--grid:#eceff3}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#14171c;--panel:#191919;--ink:#efefef;--muted:#9c9c9c;--line:#2e343e;--grid:#262a31}}
:root[data-theme=dark]{--bg:#14171c;--panel:#191919;--ink:#efefef;--muted:#9c9c9c;--line:#2e343e;--grid:#262a31}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1100px;margin:0 auto;padding:16px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin-bottom:14px}
h1{font-size:17px;margin:0 0 4px}.muted{color:var(--muted)}.small{font-size:12px}
.row{display:flex;gap:12px;flex-wrap:wrap;align-items:center}
select,button{background:var(--bg);color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:4px 8px;font:inherit}
svg{width:100%;height:auto;display:block}svg text{fill:var(--ink);font-size:12px}
.leg{display:flex;gap:6px 16px;flex-wrap:wrap;margin-top:6px;font-size:13px}.leg label{cursor:pointer;display:flex;gap:6px;align-items:center}
.sw{display:inline-block;width:22px;height:0;border-top:2px solid}
.tip{position:fixed;pointer-events:none;background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:6px 8px;font-size:12px;display:none;box-shadow:0 2px 8px rgba(0,0,0,.15);z-index:5}
table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right}th:first-child,td:first-child{text-align:left}
.tw{overflow:auto;max-height:340px}.warn{color:#9a6200}
</style></head><body><main>
<div class="card"><h1>Harmonic free energies</h1><div class="muted small" id="method"></div></div>
<div class="card"><div class="row">
<label>Quantity <select id="q"></select></label><label>Reference <select id="ref"></select></label>
<span class="muted small" id="refnote"></span></div>
<div id="plot"></div><div class="leg" id="leg"></div><div class="small warn" id="flags"></div></div>
<div class="card"><details><summary>Table</summary><div class="tw" id="table"></div></details></div>
<div class="tip" id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
(function(){
const D=JSON.parse(document.getElementById("data").textContent);
const dark=()=>document.documentElement.dataset.theme==="dark"||(document.documentElement.dataset.theme!=="light"&&matchMedia("(prefers-color-scheme: dark)").matches);
const col=i=>i<8?(dark()?D.palette_dark:D.palette)[i]:"#8a8a8a";
const Q={F:["F_total_eV_atom","F total","eV/atom",1],dF:["dF","ΔF vs reference","meV/atom",1],F_vib:["F_vib_eV_atom","F vib","meV/atom",1000],
 U_vib:["U_vib_eV_atom","U vib (incl. ZPE)","meV/atom",1000],S:["entropy_J_molK_atom","S vib","J/(mol K) per atom",1],Cv:["Cv_J_molK_atom","Cv","J/(mol K) per atom",1]};
const qs=document.getElementById("q"),rs=document.getElementById("ref");
for(const k in Q)qs.append(new Option(Q[k][1],k));qs.value=D.T.length>1?"dF":"F";
rs.append(new Option("lowest F at each T","lowest"));for(const s of D.series)rs.append(new Option(s.label,s.label));rs.value=D.reference==="none"?"lowest":D.reference;
document.getElementById("method").textContent=D.method+". F_total(T) = E_static + F_vib(T).";
const hidden=new Set();
function refCurve(){const r=rs.value;if(r!=="lowest")return D.series.find(s=>s.label===r).F_total_eV_atom;
 const pool=D.series.filter(s=>!s.flags.length);const P=pool.length?pool:D.series;
 return D.T.map((t,i)=>Math.min(...P.map(s=>s.F_total_eV_atom[i]).filter(v=>v!=null)));}
function ys(s){const k=qs.value;if(k==="dF"){const rc=refCurve();return s.F_total_eV_atom.map((v,i)=>v==null||rc[i]==null||!isFinite(rc[i])?null:1000*(v-rc[i]));}
 return s[Q[k][0]].map(v=>v==null?null:v*Q[k][3]);}
function ticks(lo,hi,n){const sp=(hi-lo)||1,st=Math.pow(10,Math.floor(Math.log10(sp/n))),e=sp/n/st,m=e>=5?10:e>=2?5:e>=1?2:1,s=m*st;const out=[];for(let v=Math.ceil(lo/s)*s;v<=hi+s*1e-9;v+=s)out.push(+v.toFixed(10));return out;}
const NS="http://www.w3.org/2000/svg";const E=(t,a,p)=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);if(p)p.append(e);return e;};
function draw(){
 document.getElementById("refnote").textContent=qs.value==="dF"?"":"(reference used for ΔF only)";
 const W=900,H=440,L=78,R=110,T0=16,B=46,pw=W-L-R,ph=H-T0-B;const svg=E("svg",{viewBox:`0 0 ${W} ${H}`,role:"img"});
 const vis=D.series.map((s,i)=>({s,i,y:ys(s)})).filter(o=>!hidden.has(o.s.label));
 const X=D.T;let lo=Infinity,hi=-Infinity;for(const o of vis)for(const v of o.y)if(v!=null&&isFinite(v)){lo=Math.min(lo,v);hi=Math.max(hi,v);}
 if(!isFinite(lo)){lo=0;hi=1;}if(hi-lo<1e-12){lo-=0.5;hi+=0.5;}const pad=(hi-lo)*0.06;lo-=pad;hi+=pad;
 const x0=X[0],x1=X.length>1?X[X.length-1]:X[0]+1;const sx=v=>L+(X.length>1?(v-x0)/(x1-x0):0.5)*pw,sy=v=>T0+(hi-v)/(hi-lo)*ph;
 for(const t of ticks(lo,hi,6)){E("line",{x1:L,x2:L+pw,y1:sy(t),y2:sy(t),stroke:"var(--grid)"},svg);const tx=E("text",{x:L-6,y:sy(t)+4,"text-anchor":"end"},svg);tx.textContent=t;}
 if(X.length>1)for(const t of ticks(x0,x1,8)){const tx=E("text",{x:sx(t),y:T0+ph+18,"text-anchor":"middle"},svg);tx.textContent=t;E("line",{x1:sx(t),x2:sx(t),y1:T0+ph,y2:T0+ph+4,stroke:"var(--muted)"},svg);}
 E("rect",{x:L,y:T0,width:pw,height:ph,fill:"none",stroke:"var(--line)"},svg);
 const xl=E("text",{x:L+pw/2,y:H-6,"text-anchor":"middle"},svg);xl.textContent="T (K)";
 const q=Q[qs.value];const yl=E("text",{x:16,y:T0+ph/2,"text-anchor":"middle",transform:`rotate(-90 16 ${T0+ph/2})`},svg);yl.textContent=`${q[1]} (${q[2]})`+(qs.value==="dF"?` vs ${rs.value==="lowest"?"lowest":rs.value}`:"");
 if(qs.value==="dF"&&lo<0&&hi>0)E("line",{x1:L,x2:L+pw,y1:sy(0),y2:sy(0),stroke:"var(--muted)","stroke-width":1},svg);
 for(const o of vis){let d="";o.y.forEach((v,k)=>{if(v==null||!isFinite(v))return;d+=(d?"L":"M")+sx(X[k]).toFixed(1)+","+sy(v).toFixed(1);});
  if(X.length===1){E("circle",{cx:sx(X[0]),cy:sy(o.y[0]),r:5,fill:col(o.i),stroke:"var(--panel)","stroke-width":2},svg);}
  else E("path",{d,fill:"none",stroke:col(o.i),"stroke-width":2,"stroke-dasharray":o.s.flags.length?"6 4":"none"},svg);
  if(vis.length<=4){const k=o.y.length-1;if(o.y[k]!=null){const tx=E("text",{x:sx(X[k])+6,y:sy(o.y[k])+4},svg);tx.textContent=o.s.label;}}}
 const cross=E("line",{y1:T0,y2:T0+ph,stroke:"var(--muted)","stroke-dasharray":"3 3",visibility:"hidden"},svg);
 const hit=E("rect",{x:L,y:T0,width:pw,height:ph,fill:"transparent"},svg);const tip=document.getElementById("tip");
 hit.addEventListener("mousemove",ev=>{const r=svg.getBoundingClientRect(),px=(ev.clientX-r.left)*W/r.width;let k=0,b=1e99;X.forEach((t,j)=>{const d=Math.abs(sx(t)-px);if(d<b){b=d;k=j;}});
  cross.setAttribute("x1",sx(X[k]));cross.setAttribute("x2",sx(X[k]));cross.setAttribute("visibility","visible");
  const dec=qs.value==="F"?6:3;tip.innerHTML=`<b>${X[k]} K</b><br>`+vis.slice().sort((a,b)=>(b.y[k]??-1e99)-(a.y[k]??-1e99)).map(o=>`<span class="sw" style="border-color:${col(o.i)}"></span> ${o.s.label}: ${o.y[k]==null?"–":o.y[k].toFixed(dec)}`).join("<br>");
  tip.style.display="block";tip.style.left=Math.min(ev.clientX+14,innerWidth-260)+"px";tip.style.top=(ev.clientY+10)+"px";});
 hit.addEventListener("mouseleave",()=>{tip.style.display="none";cross.setAttribute("visibility","hidden");});
 const p=document.getElementById("plot");p.innerHTML="";p.append(svg);
 const lg=document.getElementById("leg");lg.innerHTML="";D.series.forEach((s,i)=>{const l=document.createElement("label");const c=document.createElement("input");c.type="checkbox";c.checked=!hidden.has(s.label);
  c.onchange=()=>{c.checked?hidden.delete(s.label):hidden.add(s.label);draw();};l.append(c);const sw=document.createElement("span");sw.className="sw";sw.style.borderColor=col(i);if(s.flags.length)sw.style.borderTopStyle="dashed";l.append(sw,s.label+(s.flags.length?" ("+s.flags.join(", ")+")":""));lg.append(l);});
 const fl=D.series.filter(s=>s.flags.length);document.getElementById("flags").textContent=fl.length?"Dashed = flagged, harmonic F not reliable: "+fl.map(s=>s.label+": "+s.notes.join("; ")).join(" | "):"";
 let t="<table><thead><tr><th>T (K)</th>"+vis.map(o=>`<th>${o.s.label}</th>`).join("")+"</tr></thead><tbody>";
 X.forEach((x,k)=>{t+=`<tr><td>${x}</td>`+vis.map(o=>`<td>${o.y[k]==null?"–":o.y[k].toFixed(qs.value==="F"?6:3)}</td>`).join("")+"</tr>";});
 document.getElementById("table").innerHTML=t+"</tbody></table>";}
qs.onchange=draw;rs.onchange=draw;matchMedia("(prefers-color-scheme: dark)").addEventListener("change",draw);draw();
})();
</script></main></body></html>
"""


def write_html(series, T, ref_name_key, out):
    cols = ["F_total_eV_atom", "F_vib_eV_atom", "U_vib_eV_atom", "entropy_J_molK_atom", "Cv_J_molK_atom"]
    data = {"T": T, "method": METHOD, "reference": ref_name_key, "palette": PALETTE, "palette_dark": PALETTE_DARK,
            "series": [{"label": s["label"], "flags": s["flags"], "notes": s["notes"],
                        **{c: values_at(s, c, T) for c in cols}} for s in series]}
    p = out / "harmonic_compare.html"
    p.write_text(HTML.replace("__DATA__", json.dumps(data).replace("</", "<\\/")))
    return p


# --------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="+", help="single_runs.py result folder(s) or thermo.csv file(s)")
    ap.add_argument("--reference", default="lowest", help="lowest (default), a structure label, or none")
    ap.add_argument("--include", nargs="+", metavar="LABEL", help="only these structures (in this order)")
    ap.add_argument("--exclude", nargs="+", metavar="LABEL", help="leave these structures out")
    ap.add_argument("--T-window", dest="T_window", type=float, nargs=2, metavar=("TMIN", "TMAX"),
                    help="temperature range of the plots")
    ap.add_argument("--at-T", type=float, nargs="+", help="temperatures of the dF dot chart (default: a few of the grid)")
    ap.add_argument("--out-dir", help="where to write (default <first source>/plots)")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"], help="figure formats (default png pdf svg)")
    ap.add_argument("--width", type=float, default=3.4, help="figure width in inches (default 3.4 = one column)")
    ap.add_argument("--height", type=float, default=2.6, help="height of one panel in inches (default 2.6)")
    ap.add_argument("--font-size", type=float, default=9, help="font size in pt (default 9)")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    ap.add_argument("--no-html", action="store_true", help="no interactive HTML")
    a = ap.parse_args(argv)

    series = load(a.sources)
    if a.include:
        by = {s["label"]: s for s in series}
        miss = [l for l in a.include if l not in by]
        if miss:
            sys.exit(f"--include: unknown label(s) {' '.join(miss)}; known: {' '.join(by)}")
        series = [by[l] for l in a.include]
    if a.exclude:
        series = [s for s in series if s["label"] not in set(a.exclude)]
    if not series:
        sys.exit("no structures with thermo.csv found")
    T = common_T(series, a.T_window)
    if not T:
        sys.exit("the structures have no temperature in common (inside --T-window)")
    if a.reference not in ("lowest", "none") and a.reference not in {s["label"] for s in series}:
        sys.exit(f"--reference {a.reference}: no such structure; choose lowest, none or one of: "
                 + " ".join(s["label"] for s in series))
    out = Path(a.out_dir) if a.out_dir else (Path(a.sources[0]) if Path(a.sources[0]).is_dir()
                                             else Path(a.sources[0]).parent) / "plots"
    out.mkdir(parents=True, exist_ok=True)
    if len(series) > len(PALETTE):
        print(f"note: {len(series)} structures, only {len(PALETTE)} colours: structures beyond the 8th are grey "
              "(use --include to choose)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(plt, a.font_size)
    ref_curve, ref_name = reference_curve(series, T, a.reference)
    written = []
    if len(T) > 1:
        written += fig_F_vs_T(plt, series, T, ref_curve, ref_name, a, out)
        written += fig_panels(plt, series, T, a, out)
    at = a.at_T or ([T[0]] if len(T) == 1 else sorted({T[0], *[t for t in (300.0, 600.0, 900.0, 1200.0) if T[0] <= t <= T[-1]], T[-1]}))
    at_ok = [t for t in at if any(abs(t - x) < 1e-6 for x in T)]
    if len(at_ok) < len(at):
        print("note: not on the common grid, skipped: " + " ".join(f"{t:g}" for t in at if t not in at_ok) + " K")
    dF_data = []
    if at_ok and len(series) > 1:
        p, dF_data = fig_dF_at_T(plt, series, at_ok, a.reference, a, out)
        written += p
    written += fig_dos(plt, series, a, out)
    written.append(write_csv(series, T, ref_curve, ref_name, out))
    if not a.no_html:
        written.append(write_html(series, T, a.reference, out))

    print(f"{len(series)} structure(s), {len(T)} temperature(s) {T[0]:g}..{T[-1]:g} K, reference: {ref_name}")
    print(METHOD)
    for t, rname, d in dF_data:
        print(f"  dF at {t:g} K vs {rname} (meV/atom): " +
              ", ".join(f"{s['label']} {v:+.3f}" for s, v in zip(series, d) if v is not None))
    for s in series:
        if s["flags"]:
            print(f"  FLAGGED {s['label']}: {'; '.join(s['notes']) or ' '.join(s['flags'])}")
    for p in written:
        print(f"  wrote {p}")


if __name__ == "__main__":
    main()
