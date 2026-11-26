#!/usr/bin/env python3
"""
single_runs.py - single runs on a few POSCARs with one NEP model: relaxation, harmonic free energy at one
temperature or over a temperature ramp, powder XRD with optional references, and comparison plots.

    NEP -> calorine CPUNEP -> (optional relaxation) -> phonopy finite displacements -> thermal properties
    F_total(T) = E_static + F_vib(T)   [eV/atom]
FIXED-CELL HARMONIC approximation (not quasi-harmonic): one phonon calculation at the (relaxed) cell is used
for every temperature, so there is no thermal expansion.  A ramp costs the same as one temperature.

Usage
    python single_runs.py POSCAR_a POSCAR_b --nep nep.txt --T 300
    python single_runs.py POSCAR_* --nep nep.txt --relax --T-range 0 1500 10 --reference POSCAR_a
    python single_runs.py POSCAR_m --nep nep.txt --relax --tasks relax harmonic xrd plot --xrd-ref c t oIII m
    python single_runs.py POSCAR_m --tasks xrd --xrd-ref m=refs/POSCAR_m            (XRD only: no NEP needed)

Tasks (--tasks, default: harmonic plot; --relax adds relax)
    relax     NEP relaxation (cell + positions, or --fixed-cell), same settings and convergence check as
              batch_free_energy.py: relax_converged = False if the max force (incl. the cell filter) > --fmax
    harmonic  phonons once, then F_vib, F_total, U_vib (vibrational internal energy, includes the zero-point
              energy), ZPE, S, Cv on the temperature grid (--T and/or --T-range), phonon DOS
    xrd       powder XRD of the structure used (relaxed if relaxed) with scripts/xrd_average.py, optional
              references (--xrd-ref NAME or NAME=POSCAR, from --refs); species are never hard-coded
    plot      plot_harmonic_compare.py on the results: F(T), dF(T) vs --reference, S / Cv / U_vib panels,
              dF at chosen temperatures, CSV + interactive HTML

Output (--out-dir, default ./single_runs):
    <label>/input.vasp, <label>/relaxed.vasp, <label>/thermo.csv (one row per T), <label>/dos.csv,
    <label>/result.json (flags, settings, md5s), <label>/xrd/...
    summary.csv           one row per structure (at --T-report, default the first T of the grid)
    status.json           progress, read by the Studio's Single runs page
    single_runs.json      every setting, script md5, versions, the command line
    plots/                figures (PNG / PDF / SVG), comparison CSV, HTML
Re-running the same command skips structures already finished with the same settings (safe to restart);
--force recomputes.  Imaginary modes and unconverged relaxations are flagged in every output: harmonic F is
not reliable for such structures.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from batch_free_energy import auto_supercell, parse_title  # noqa: E402  (same folder, same conventions)

TASKS = ("relax", "harmonic", "xrd", "plot")
THERMO_COLS = ["T_K", "E_static_eV_atom", "F_vib_eV_atom", "F_total_eV_atom", "U_vib_eV_atom",
               "entropy_J_molK_atom", "Cv_J_molK_atom", "TS_eV_atom"]
SUMMARY_COLS = ["label", "file", "formula", "phase", "id", "composition", "n_atoms", "status",
                "T_K", "E_static_eV_atom", "F_vib_eV_atom", "F_total_eV_atom", "ZPE_eV_atom",
                "entropy_J_molK_atom", "Cv_J_molK_atom",
                "volume_A3_atom", "a_A", "b_A", "c_A", "alpha_deg", "beta_deg", "gamma_deg",
                "relaxed", "relax_converged", "max_force_eV_A", "rms_relax_shift_A",
                "phonon_supercell", "min_freq_THz", "n_imaginary_modes", "flags", "note", "time_s"]
IMAG_THz = -0.05                      # same threshold as batch_free_energy.py
EV_PER_J_MOL = 1.0 / 96485.33212      # J/mol -> eV per particle


# --------------------------------------------------------------------------- small helpers
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fmt(v):
    return f"{v:.10g}" if isinstance(v, float) else ("" if v is None else v)


def wjson(p, obj):
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=str))
    os.replace(tmp, p)


def temperature_grid(T, T_range):
    temps = []
    if T_range:
        lo, hi, step = T_range
        if step <= 0 or hi < lo:
            raise ValueError("--T-range needs TMIN <= TMAX and STEP > 0")
        n = int(np.floor((hi - lo) / step + 1e-9)) + 1
        temps += [lo + i * step for i in range(n)]
        if abs(temps[-1] - hi) > 1e-9:
            temps.append(hi)
    temps += list(T or [])
    temps = sorted({round(float(t), 6) for t in temps})
    if any(t < 0 for t in temps):
        raise ValueError("temperatures must be >= 0 K")
    return temps


def unique_labels(paths, names=None):
    if names:
        if len(names) != len(paths):
            raise ValueError("--names needs one name per POSCAR")
        labels = [re.sub(r"[^\w.+-]+", "_", n) for n in names]
        if len(set(labels)) != len(labels):
            raise ValueError("--names must be different from each other")
        return labels
    stems = [re.sub(r"[^\w.+-]+", "_", Path(p).name) for p in paths]
    out = []
    for p, s in zip(paths, stems):
        lab = s
        if stems.count(s) > 1:                                   # same file name in different folders
            lab = re.sub(r"[^\w.+-]+", "_", Path(p).parent.name) + "_" + s
        k, base = 2, lab
        while lab in out:
            lab = f"{base}_{k}"
            k += 1
        out.append(lab)
    return out


def settings_key(a, temps):
    d = {"relax": a.relax, "fmax": a.fmax, "steps": a.steps, "fixed_cell": a.fixed_cell,
         "supercell": a.supercell, "min_length": a.min_length, "mesh": a.mesh, "temps": temps,
         "nep_md5": md5(a.nep) if a.nep else None}
    return hashlib.md5(json.dumps(d, sort_keys=True).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------- one structure
def relax(atoms, a, src):
    from calorine.tools import relax_structure
    start = atoms.get_positions().copy()
    relax_structure(atoms, fmax=a.fmax, steps=a.steps, minimizer="lbfgs", constant_cell=a.fixed_cell, logfile=None)
    d = atoms.get_positions() - start
    d -= d.mean(axis=0)
    rms = float(np.sqrt((d ** 2).sum(axis=1).mean()))
    if a.fixed_cell:
        mf = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
    else:                                            # what the optimiser checks: atoms + cell "forces"
        from ase.filters import FrechetCellFilter
        mf = float(np.linalg.norm(FrechetCellFilter(atoms).get_forces(), axis=1).max())
    return rms, mf, bool(mf <= 1.0001 * a.fmax)


def write_vasp(path, atoms, title):
    from ase.io import write
    write(str(path), atoms, format="vasp", direct=True, sort=False)
    lines = path.read_text().split("\n")
    lines[0] = title
    path.write_text("\n".join(lines))


def harmonic(atoms, a, temps):
    from ase.units import kJ, mol
    from calorine.tools import get_force_constants
    sc = a.supercell or auto_supercell(np.array(atoms.cell), a.min_length)
    phonon = get_force_constants(atoms, atoms.calc, sc)
    mesh = phonon.run_mesh(a.mesh)
    freqs = np.asarray(mesh.frequencies)
    grid = sorted(set([0.0] + list(temps)))                      # T = 0 gives the zero-point energy
    tp = phonon.run_thermal_properties(temperatures=grid)
    n_prim = len(phonon.primitive)
    T = np.asarray(tp.temperatures, dtype=float)
    F = np.asarray(tp.free_energy, dtype=float) * kJ / mol / n_prim          # eV/atom
    S = np.asarray(tp.entropy, dtype=float) / n_prim                         # J/(mol K) per atom
    Cv = np.asarray(tp.heat_capacity, dtype=float) / n_prim
    TS = T * S * EV_PER_J_MOL                                                 # eV/atom
    zpe = float(F[np.argmin(np.abs(T))])
    dos = None
    try:
        phonon.run_total_dos()
        td = phonon.total_dos
        dos = (np.asarray(td.frequency_points), np.asarray(td.dos) / n_prim)
    except Exception:                                # noqa: BLE001  DOS is only a figure: never fail for it
        pass
    keep = [i for i, t in enumerate(T) if any(abs(t - x) < 1e-6 for x in temps)]
    return {"sc": sc, "freqs": freqs, "T": T[keep], "F_vib": F[keep], "S": S[keep], "Cv": Cv[keep],
            "TS": TS[keep], "U_vib": F[keep] + TS[keep], "zpe": zpe, "dos": dos}


def run_structure(path, label, a, temps, tasks, out_dir, key):
    d = out_dir / label
    d.mkdir(parents=True, exist_ok=True)
    res_p = d / "result.json"
    old = json.loads(res_p.read_text()) if res_p.is_file() else {}
    need_calc = bool({"relax", "harmonic"} & tasks)
    if need_calc and not a.force and old.get("status") == "ok" and old.get("key") == key \
            and old.get("input_md5") == md5(path):
        print(f"  {label}: already done with the same settings (use --force to recompute)", flush=True)
        return old
    from ase.io import read
    shutil.copyfile(path, d / "input.vasp")
    info = parse_title(Path(path))
    atoms = read(str(path), format="vasp")
    res = {"label": label, "file": str(Path(path).resolve()), "input_md5": md5(path), "key": key,
           "formula": atoms.get_chemical_formula(mode="metal"), "n_atoms": len(atoms), "status": "ok",
           "phase": info.get("phase") if info.get("phase") != "parent" else "", "id": info.get("id", ""),
           "composition": info.get("composition", ""), "relaxed": False, "flags": [], "notes": []}
    t0 = time.time()
    try:
        if need_calc:
            from calorine.calculators import CPUNEP
            atoms.calc = CPUNEP(a.nep)
            if "relax" in tasks:
                print(f"  {label}: relaxing ({'positions only' if a.fixed_cell else 'cell + positions'}, "
                      f"fmax {a.fmax:g} eV/Å, max {a.steps} steps)", flush=True)
                rms, mf, conv = relax(atoms, a, path)
                res.update(relaxed=True, rms_relax_shift_A=rms, max_force_eV_A=mf, relax_converged=conv)
                title = Path(path).open().readline().strip() + f" | relaxed with {Path(a.nep).name}"
                write_vasp(d / "relaxed.vasp", atoms, title)
                if not conv:
                    res["flags"].append("relax_not_converged")
                    res["notes"].append(f"relaxation NOT converged (max force {mf:.2e} > fmax {a.fmax:g} after "
                                        f"{a.steps} steps)")
            res["E_static_eV_atom"] = atoms.get_potential_energy() / len(atoms)
            if "relax" not in tasks:
                res["max_force_eV_A"] = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
            L, A = atoms.cell.lengths(), atoms.cell.angles()
            res.update(volume_A3_atom=atoms.get_volume() / len(atoms), a_A=L[0], b_A=L[1], c_A=L[2],
                       alpha_deg=A[0], beta_deg=A[1], gamma_deg=A[2])
        if "harmonic" in tasks:
            print(f"  {label}: phonons ({len(atoms)} atoms) ...", flush=True)
            h = harmonic(atoms, a, temps)
            e0 = res["E_static_eV_atom"]
            with (d / "thermo.csv").open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(THERMO_COLS)
                for i, T in enumerate(h["T"]):
                    w.writerow([fmt(float(x)) for x in (T, e0, h["F_vib"][i], e0 + h["F_vib"][i], h["U_vib"][i],
                                                       h["S"][i], h["Cv"][i], h["TS"][i])])
            if h["dos"] is not None:
                with (d / "dos.csv").open("w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["frequency_THz", "dos_states_per_THz_atom"])
                    w.writerows([[fmt(float(x)), fmt(float(y))] for x, y in zip(*h["dos"])])
            nimag = int(np.count_nonzero(h["freqs"] < IMAG_THz))
            res.update(phonon_supercell="x".join(map(str, h["sc"])), min_freq_THz=float(h["freqs"].min()),
                       n_imaginary_modes=nimag, ZPE_eV_atom=h["zpe"], harmonic=True)
            if nimag:
                res["flags"].append("imaginary_modes")
                res["notes"].append(f"{nimag} imaginary modes (min {h['freqs'].min():.2f} THz): harmonic F not "
                                    "reliable (structure may be unstable in this phase)")
            if "relax" not in tasks and res.get("max_force_eV_A", 0) > 1e-2:
                res["flags"].append("not_relaxed")
                res["notes"].append(f"input not relaxed (max force {res['max_force_eV_A']:.2e} eV/Å): phonons "
                                    "of a structure away from its minimum")
    except Exception as exc:                                     # keep going with the next structure
        res.update(status="FAILED", notes=res["notes"] + [f"{type(exc).__name__}: {exc}".replace("\n", " ")[:400]])
        print(f"  {label}: FAILED {res['notes'][-1]}", flush=True)
    res["time_s"] = round(time.time() - t0, 1)
    wjson(res_p, res)
    return res


def run_xrd(label, a, out_dir, res):
    d = out_dir / label
    geo = d / "relaxed.vasp" if (d / "relaxed.vasp").is_file() else d / "input.vasp"
    if not geo.is_file():
        geo = Path(res["file"])
    cmd = [sys.executable, str(HERE / "xrd_average.py"), "--structures", str(geo), "--out-dir", str(d / "xrd"),
           "--stem", "xrd", "--formats", *a.formats, "--jobs", "1", "--refs", a.refs]
    if a.xrd_ref:
        cmd += ["--ref", *a.xrd_ref]
    cmd += shlex.split(a.xrd_args or "")
    print(f"  {label}: XRD of {geo.name}" + (f" + references {' '.join(a.xrd_ref)}" if a.xrd_ref else ""), flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(r.stdout[-1500:] + r.stderr[-2500:], flush=True)
        return f"XRD failed (exit {r.returncode})"
    return None


# --------------------------------------------------------------------------- outputs
def plot_xrd_compare(out_dir, labels, formats):
    """every structure's pattern stacked (same colour per structure as in the free-energy figures), the
    references once below; intensities each max = 100"""
    pats, refs, x = [], {}, None
    for i, lab in enumerate(labels):
        p = out_dir / lab / "xrd" / "xrd.csv"
        if not p.is_file():
            continue
        with p.open() as f:
            rows = list(csv.DictReader(f))
        if not rows:
            continue
        x = np.array([float(r["two_theta_deg"]) for r in rows])
        pats.append((i, lab, np.array([float(r["mean_rel_max100"]) for r in rows])))
        for k in rows[0]:
            if k.startswith("ref_") and k.endswith("_max100") and k not in refs:
                refs[k] = np.array([float(r[k]) for r in rows])
    if not pats:
        return []
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE))
    from plot_harmonic_compare import METHOD, PALETTE, OVERFLOW, style  # noqa: F401
    style(plt, 9)
    n = len(pats) + len(refs)
    fig, ax = plt.subplots(figsize=(7.0, 1.0 + 0.62 * n))
    off, step = 0.0, 110.0
    for i, lab, y in pats:
        c = PALETTE[i] if i < len(PALETTE) else OVERFLOW
        ax.plot(x, off + y, color=c, lw=0.9)
        ax.text(x[0] + 0.01 * (x[-1] - x[0]), off + 12, lab, ha="left", va="bottom", fontsize=8, color="#222222")
        off -= step
    for k, y in refs.items():
        name = k[4:-7]
        ax.plot(x, off + y, color="#6b6b6b", lw=0.8)
        ax.text(x[0] + 0.01 * (x[-1] - x[0]), off + 12, f"reference {name}", ha="left", va="bottom", fontsize=8, color="#444444")
        off -= step
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(off + step - 6, 108)
    ax.set_yticks([])
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel(r"2$\theta$ (deg)")
    ax.set_ylabel("Intensity (each pattern max = 100, offset)")
    paths = []
    (out_dir / "plots").mkdir(exist_ok=True)
    for f in formats:
        pth = out_dir / "plots" / f"xrd_compare.{f}"
        fig.savefig(pth, dpi=300)
        paths.append(pth)
    plt.close(fig)
    return paths


def write_summary(out_dir, results, T_report):
    rows = []
    for r in results:
        row = {k: r.get(k, "") for k in SUMMARY_COLS}
        row["flags"] = " ".join(r.get("flags", []))
        row["note"] = "; ".join(r.get("notes", []))
        tp = out_dir / r["label"] / "thermo.csv"
        if r.get("status") == "ok" and r.get("harmonic") and tp.is_file():
            with tp.open() as f:
                th = list(csv.DictReader(f))
            if th:
                best = min(th, key=lambda x: abs(float(x["T_K"]) - T_report))
                for k in ("T_K", "F_vib_eV_atom", "F_total_eV_atom", "entropy_J_molK_atom", "Cv_J_molK_atom"):
                    row[k] = best[k]
        rows.append(row)
    with (out_dir / "summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLS)
        w.writeheader()
        w.writerows({k: fmt(v) for k, v in row.items()} for row in rows)
    return rows


def versions():
    out = {"python": platform.python_version()}
    for m in ("numpy", "ase", "phonopy", "calorine", "pymatgen", "matplotlib"):
        try:
            out[m] = __import__(m).__version__
        except Exception:                            # noqa: BLE001
            pass
    return out


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("poscars", nargs="+", metavar="POSCAR", help="structure files (VASP format)")
    ap.add_argument("--nep", help="NEP model file (nep.txt); needed for relax / harmonic")
    ap.add_argument("--tasks", nargs="+", choices=TASKS, help="what to do (default: harmonic plot; --relax adds relax)")
    ap.add_argument("--out-dir", default="single_runs", help="results folder (default ./single_runs)")
    ap.add_argument("--names", nargs="+", help="labels of the structures (default: file names)")
    ap.add_argument("--force", action="store_true", help="recompute structures already done")
    g = ap.add_argument_group("temperatures (fixed-cell harmonic: one phonon run serves every T)")
    g.add_argument("--T", type=float, nargs="+", help="temperature(s) in K, e.g. 300 or 300 600 900")
    g.add_argument("--T-range", dest="T_range", type=float, nargs=3, metavar=("TMIN", "TMAX", "STEP"),
                   help="temperature ramp, e.g. 0 1500 10")
    g.add_argument("--T-report", type=float, help="temperature of summary.csv (default: the first --T, else TMIN)")
    g = ap.add_argument_group("relaxation (same settings as batch_free_energy.py)")
    g.add_argument("--relax", action="store_true", help="relax with the NEP before the phonons (= task relax)")
    g.add_argument("--fmax", type=float, default=1e-4, help="force threshold, eV/Å (default 1e-4)")
    g.add_argument("--steps", type=int, default=2000, help="max relaxation steps (default 2000)")
    g.add_argument("--fixed-cell", action="store_true", help="relax positions only, keep the cell")
    g = ap.add_argument_group("phonons")
    g.add_argument("--supercell", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                   help="phonon supercell (default: automatic from --min-length)")
    g.add_argument("--min-length", type=float, default=10.0, help="automatic supercell: every width >= this (Å)")
    g.add_argument("--mesh", type=float, default=100, help="phonopy q-point mesh density (default 100)")
    g.add_argument("--threads", type=int, default=None, help="BLAS/OpenMP threads (default: as set in the environment)")
    g = ap.add_argument_group("XRD (scripts/xrd_average.py)")
    g.add_argument("--xrd-ref", nargs="+", metavar="NAME[=POSCAR]",
                   help="reference structures: names of POSCAR_<name> in --refs, or NAME=path")
    g.add_argument("--refs", default=str(HERE.parent / "refs"), help="folder with POSCAR_<name> (default ../refs)")
    g.add_argument("--xrd-args", default="", help="more xrd_average.py options, e.g. \"--fwhm 0.08 --wavelength MoKa\"")
    g = ap.add_argument_group("plots (scripts/plot_harmonic_compare.py)")
    g.add_argument("--reference", default="lowest",
                   help="reference of dF: a label, 'lowest' (lowest F at each T, default) or 'none'")
    g.add_argument("--at-T", type=float, nargs="+", help="temperatures of the dF bar chart (default: --T)")
    g.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"], help="figure formats (default png pdf svg)")
    return ap


def main(argv=None):
    ap = build_parser()
    a = ap.parse_args(argv)
    tasks = set(a.tasks or ["harmonic", "plot"])
    if a.relax:
        tasks.add("relax")
    a.relax = "relax" in tasks
    if {"relax", "harmonic"} & tasks:
        if not a.nep:
            ap.error("--nep is needed for relax / harmonic")
        if not Path(a.nep).is_file():
            sys.exit(f"NEP model {a.nep} not found")
    else:
        a.nep = None
    temps = []
    if "harmonic" in tasks:
        if not a.T and not a.T_range:
            ap.error("harmonic needs --T and/or --T-range")
        try:
            temps = temperature_grid(a.T, a.T_range)
        except ValueError as e:
            ap.error(str(e))
    for p in a.poscars:
        if not Path(p).is_file():
            sys.exit(f"{p} not found")
    if a.threads:
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[v] = str(a.threads)
    try:
        labels = unique_labels(a.poscars, a.names)
    except ValueError as e:
        ap.error(str(e))
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    T_report = a.T_report if a.T_report is not None else (a.T[0] if a.T else (temps[0] if temps else 0.0))
    key = settings_key(a, temps)
    wjson(out_dir / "single_runs.json", {
        "script": str(Path(__file__).resolve()), "script_md5": md5(__file__),
        "command": shlex.join([sys.executable, *sys.argv]) if argv is None else shlex.join(["single_runs.py", *argv]),
        "cwd": os.getcwd(), "started": time.strftime("%Y-%m-%d %H:%M:%S"), "host": platform.node(),
        "tasks": sorted(tasks, key=TASKS.index), "temperatures_K": temps, "T_report": T_report,
        "structures": [{"label": l, "file": str(Path(p).resolve()), "md5": md5(p)} for l, p in zip(labels, a.poscars)],
        "nep": str(Path(a.nep).resolve()) if a.nep else None, "nep_md5": md5(a.nep) if a.nep else None,
        "settings": {k: v for k, v in vars(a).items() if k not in ("poscars",)}, "versions": versions(),
        "method": "fixed-cell harmonic approximation (not quasi-harmonic): no thermal expansion"})
    status = {"state": "running", "start": time.time(), "n": len(labels), "done": 0, "current": None,
              "structures": {l: "waiting" for l in labels}}
    wjson(out_dir / "status.json", status)
    print(f"{len(labels)} structure(s), tasks: {' '.join(sorted(tasks, key=TASKS.index))}"
          + (f", {len(temps)} temperature(s) {temps[0]:g}..{temps[-1]:g} K" if temps else "")
          + f" -> {out_dir}", flush=True)
    if "harmonic" in tasks:
        print("method: fixed-cell harmonic approximation (not quasi-harmonic): no thermal expansion", flush=True)
    results, n_fail = [], 0
    for i, (p, lab) in enumerate(zip(a.poscars, labels), 1):
        status.update(current=lab)
        status["structures"][lab] = "running"
        wjson(out_dir / "status.json", status)
        print(f"[{i}/{len(labels)}] {lab}  ({p})", flush=True)
        res = run_structure(Path(p), lab, a, temps, tasks, out_dir, key)
        if "xrd" in tasks and res.get("status") == "ok":
            err = run_xrd(lab, a, out_dir, res)
            res["xrd"] = err or "ok"
            if err:
                res["notes"].append(err)
            wjson(out_dir / lab / "result.json", res)
        n_fail += res.get("status") != "ok"
        results.append(res)
        status["done"] = i
        status["structures"][lab] = res["status"] + (" (" + ", ".join(res["flags"]) + ")" if res.get("flags") else "")
        wjson(out_dir / "status.json", status)
        if res.get("status") == "ok" and "E_static_eV_atom" in res:
            print(f"  {lab}: E_static {res['E_static_eV_atom']:.6f} eV/atom" +
                  (f", ZPE {res['ZPE_eV_atom'] * 1000:.2f} meV/atom, sc {res['phonon_supercell']}, "
                   f"imag {res['n_imaginary_modes']}" if res.get("harmonic") else "") +
                  ("  FLAGS: " + ", ".join(res["flags"]) if res["flags"] else ""), flush=True)
            for n in res.get("notes", []):
                print(f"    ! {n}", flush=True)
    write_summary(out_dir, results, T_report)
    plot_err = None
    if "plot" in tasks and any(r.get("harmonic") and r.get("status") == "ok" for r in results):
        cmd = [sys.executable, str(HERE / "plot_harmonic_compare.py"), str(out_dir), "--reference", a.reference,
               "--formats", *a.formats]
        at = a.at_T or a.T
        if at:
            cmd += ["--at-T", *[f"{t:g}" for t in at]]
        print("plots: " + shlex.join(cmd[1:]), flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout[-3000:], flush=True)
        if r.returncode:
            print(r.stderr[-3000:], flush=True)
            plot_err = f"plot failed (exit {r.returncode})"
    if "xrd" in tasks and "plot" in tasks:
        try:
            for pth in plot_xrd_compare(out_dir, labels, a.formats):
                print(f"  wrote {pth}", flush=True)
        except Exception as exc:                     # noqa: BLE001  a figure must not lose the results
            print(f"XRD comparison figure failed: {type(exc).__name__}: {exc}", flush=True)
    status.update(state="done" if not n_fail and not plot_err else "done_with_errors", end=time.time(), current=None,
                  failed=n_fail, plot_error=plot_err)
    wjson(out_dir / "status.json", status)
    print(f"\n{len(results) - n_fail} ok, {n_fail} failed -> {out_dir / 'summary.csv'}")
    flagged = [r["label"] for r in results if r.get("flags")]
    if flagged:
        print("FLAGGED (do not trust their harmonic F): " + ", ".join(flagged))
    if n_fail or plot_err:
        sys.exit(1)


if __name__ == "__main__":
    main()
