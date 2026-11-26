#!/usr/bin/env python3
"""
batch_free_energy.py - harmonic free energy at ONE temperature for every POSCAR in a folder
(including sub-folders such as mapped/t, mapped/m, ...), written to one CSV.

Same method and units as thermal_harmonic_nep.py:
    NEP -> calorine CPUNEP -> (optional relaxation) -> phonopy finite displacements
    F_total(T) = E_static + F_vib(T)        [eV/atom], fixed-cell harmonic (not quasi-harmonic)

Usage
    python batch_free_energy.py Ce_33.33 --nep nep.txt --T 300
    python batch_free_energy.py Ce_33.33 --nep nep.txt --T 300 --relax --fmax 1e-3
    python batch_free_energy.py Ce_33.33 --nep nep.txt --T 300 --supercell 2 2 2 --mesh 50

Output  (default  <folder>/free_energy_<T>K.csv), one row per structure:
    file, phase, id, composition, WC values (from the POSCAR's first line),
    E_static, F_vib, F_total, dF_total (meV/atom, relative to the lowest F_total of the same composition),
    entropy, Cv, volume/atom, cell lengths and angles (after relaxation), rms relaxation shift,
    phonon supercell, lowest frequency, number of imaginary modes, status
* Rows are written as soon as each structure finishes; re-running the same command skips
  structures already done (safe to restart after an interruption) and retries failed ones.
* With --relax the relaxed structures are saved to <folder>/relaxed/<same sub-path>.
* Parallel: --jobs N runs N structures at the same time (separate processes, default = all cores),
  each limited to --threads BLAS/OpenMP threads (default 1) so the workers do not oversubscribe the CPU.
* --files-from LIST: only the structures listed in LIST (one path per line, relative to the folder),
  e.g. one chunk of a large campaign (campaign.py); use a separate --out per list.
* Relaxation convergence is checked: max_force_eV_A after relaxing, relax_converged = False (and a note)
  if it is above --fmax after --steps steps.  Such rows are kept but should not be trusted.
* Phonon supercell: by default the smallest repeat with every cell width >= --min-length
  (10 Å), since a fixed 5x5x5 of a 96-atom cell would be 12000 atoms.
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import multiprocessing as mp
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

CSV_COLS = [
    "file", "phase", "id", "composition", "WC",
    "T_K", "E_static_eV_atom", "F_vib_eV_atom", "F_total_eV_atom", "dF_total_meV_atom",
    "entropy_J_molK_atom", "Cv_J_molK_atom",
    "n_atoms", "volume_A3_atom", "a_A", "b_A", "c_A", "alpha_deg", "beta_deg", "gamma_deg",
    "relaxed", "rms_relax_shift_A", "phonon_supercell", "min_freq_THz", "n_imaginary_modes",
    "status", "note", "time_s",
    "max_force_eV_A", "relax_converged",
]


# --------------------------------------------------------------------------- files / labels
def find_structures(folder: Path, pattern: str, skip_dirs: set[str]) -> list[Path]:
    out = []
    for root, dirs, files in os.walk(folder):
        dirs[:] = sorted(d for d in dirs if d not in skip_dirs and not d.startswith("."))
        for f in sorted(files):
            if fnmatch.fnmatch(f, pattern):
                out.append(Path(root) / f)
    return out


def parse_title(path: Path) -> dict:
    """composition, WC values, id and phase from the first line written by the other scripts"""
    first = path.open().readline().strip()
    info = {"composition": "", "WC": "", "id": "", "phase": "parent"}
    if "|" in first:
        info["composition"] = first.split("|")[0].strip()
    info["WC"] = " ".join(re.findall(r"WC_shell\d+\([^)]*\)=[+-]?[\d.]+", first))
    m = re.search(r"id=(\d+)", first)
    if m:
        info["id"] = m.group(1)
    m = re.search(r"phase=([^\s|]+)", first)
    if m:
        info["phase"] = m.group(1)
    return info


def auto_supercell(cell: np.ndarray, min_length: float) -> list[int]:
    widths = 1.0 / np.linalg.norm(np.linalg.inv(cell), axis=0)        # perpendicular widths
    return [max(1, int(np.ceil(min_length / w - 1e-6))) for w in widths]


# --------------------------------------------------------------------------- calculation
def free_energy_at_T(path, nep_file, T, supercell, min_length, mesh, relax, fmax, fixed_cell,
                     relaxed_out, steps=2000):
    from ase.io import read, write
    from ase.units import kJ, mol
    from calorine.calculators import CPUNEP
    from calorine.tools import get_force_constants, relax_structure

    atoms = read(str(path), format="vasp")
    atoms.calc = CPUNEP(nep_file)
    start = atoms.get_positions().copy()
    rms_shift = 0.0
    converged = ""
    if relax:
        relax_structure(atoms, fmax=fmax, steps=steps, minimizer="lbfgs", constant_cell=fixed_cell,
                        logfile=None)
        d = atoms.get_positions() - start
        d -= d.mean(axis=0)
        rms_shift = float(np.sqrt((d ** 2).sum(axis=1).mean()))
        if relaxed_out is not None:
            relaxed_out.parent.mkdir(parents=True, exist_ok=True)
            write(str(relaxed_out), atoms, format="vasp", direct=True, sort=False)
            # put the original first line (composition, WC, id, phase) back; works for every ASE version
            title = path.open().readline().strip() + f" | relaxed with {Path(nep_file).name}"
            lines = relaxed_out.read_text().split("\n")
            lines[0] = title
            relaxed_out.write_text("\n".join(lines))

    e_static = atoms.get_potential_energy() / len(atoms)
    max_force = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
    if relax:
        if not fixed_cell:                       # what the optimiser checks: atom forces + cell "forces" of the filter
            from ase.filters import FrechetCellFilter
            g = FrechetCellFilter(atoms).get_forces()
            max_force = float(np.linalg.norm(g, axis=1).max())
        converged = bool(max_force <= 1.0001 * fmax)

    sc = supercell or auto_supercell(np.array(atoms.cell), min_length)
    phonon = get_force_constants(atoms, atoms.calc, sc)
    mesh_result = phonon.run_mesh(mesh)
    # phonopy includes t_max; a 1 K window around T is enough
    tp = phonon.run_thermal_properties(t_min=T, t_max=T + 1.0, t_step=1.0)
    temps = np.asarray(tp.temperatures, dtype=float)
    k = int(np.argmin(np.abs(temps - T)))
    n_prim = len(phonon.primitive)
    f_vib = float(np.asarray(tp.free_energy)[k]) * kJ / mol / n_prim
    s_atom = float(np.asarray(tp.entropy)[k]) / n_prim
    cv_atom = float(np.asarray(tp.heat_capacity)[k]) / n_prim

    freqs = np.asarray(mesh_result.frequencies)
    lengths = atoms.cell.lengths()
    angles = atoms.cell.angles()
    return {
        "T_K": float(temps[k]),
        "E_static_eV_atom": e_static,
        "F_vib_eV_atom": f_vib,
        "F_total_eV_atom": e_static + f_vib,
        "entropy_J_molK_atom": s_atom,
        "Cv_J_molK_atom": cv_atom,
        "n_atoms": len(atoms),
        "volume_A3_atom": atoms.get_volume() / len(atoms),
        "a_A": lengths[0], "b_A": lengths[1], "c_A": lengths[2],
        "alpha_deg": angles[0], "beta_deg": angles[1], "gamma_deg": angles[2],
        "relaxed": relax,
        "rms_relax_shift_A": rms_shift,
        "phonon_supercell": "x".join(map(str, sc)),
        "min_freq_THz": float(freqs.min()),
        "n_imaginary_modes": int(np.count_nonzero(freqs < -0.05)),
        "max_force_eV_A": max_force,
        "relax_converged": converged,
    }


def run_one(job: dict) -> dict:
    """one structure, executed in a worker process; never raises"""
    path, folder = Path(job["path"]), Path(job["folder"])
    rel = str(path.relative_to(folder))
    info = parse_title(path)
    row = {"file": rel, **info, "status": "ok", "note": ""}
    t0 = time.time()
    try:
        res = free_energy_at_T(path, job["nep"], job["T"], job["supercell"], job["min_length"], job["mesh"],
                               job["relax"], job["fmax"], job["fixed_cell"],
                               folder / "relaxed" / rel if job["relax"] else None, job.get("steps", 2000))
        row.update(res)
        notes = []
        if res["relax_converged"] is False:
            notes.append(f"relaxation NOT converged (max force {res['max_force_eV_A']:.2e} > fmax "
                         f"{job['fmax']:g} after {job.get('steps', 2000)} steps)")
        if res["n_imaginary_modes"] > 0:
            notes.append("imaginary modes: harmonic F not reliable (structure may be unstable in this phase)")
        row["note"] = "; ".join(notes)
    except Exception as exc:                                         # keep going with the next file
        row.update(status="FAILED", note=f"{type(exc).__name__}: {exc}".replace("\n", " ")[:300], T_K=job["T"])
    row["time_s"] = f"{time.time() - t0:.1f}"
    return row


def describe(row: dict) -> str:
    if row["status"] != "ok":
        return f"FAILED: {row['note']}"
    return (f"F = {row['F_total_eV_atom']:.6f} eV/atom  (E0 {row['E_static_eV_atom']:.6f}, "
            f"Fvib {row['F_vib_eV_atom']:+.6f}), sc {row['phonon_supercell']}, imag {row['n_imaginary_modes']}"
            + ("  RELAX NOT CONVERGED" if row.get("relax_converged") is False else ""))


# --------------------------------------------------------------------------- CSV
def read_done(csv_path: Path) -> list[dict]:
    if not csv_path.exists():
        return []
    with csv_path.open() as f:
        return list(csv.DictReader(f))


def append_row(csv_path: Path, row: dict) -> None:
    new = not csv_path.exists()
    if not new:                                   # CSV from an older version (fewer columns): rewrite its header
        with csv_path.open() as f:
            header = f.readline().strip().split(",")
        if header != CSV_COLS:
            old = read_done(csv_path)
            with csv_path.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=CSV_COLS)
                w.writeheader()
                w.writerows({k: r.get(k, "") or "" for k in CSV_COLS} for r in old)
    with csv_path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLS)
        if new:
            w.writeheader()
        w.writerow({k: fmt(row.get(k, "")) for k in CSV_COLS})


def fmt(v):
    if isinstance(v, float):
        return f"{v:.10g}"
    return v


def finalize(csv_path: Path) -> list[dict]:
    """fill dF_total (meV/atom vs the lowest F_total of the same composition) and sort"""
    rows = list({r["file"]: r for r in read_done(csv_path)}.values())   # latest row per file (retries)
    best = {}
    for r in rows:
        if r["status"] == "ok":
            f = float(r["F_total_eV_atom"])
            best[r["composition"]] = min(best.get(r["composition"], f), f)
    for r in rows:
        if r["status"] == "ok":
            r["dF_total_meV_atom"] = f"{1000 * (float(r['F_total_eV_atom']) - best[r['composition']]):.4f}"
    rows.sort(key=lambda r: (r["composition"], r["status"] != "ok",
                             float(r["F_total_eV_atom"]) if r["status"] == "ok" else 0.0))
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLS)
        w.writeheader()
        w.writerows({k: r.get(k, "") or "" for k in CSV_COLS} for r in rows)
    return rows


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", help="e.g. Ce_33.33 (sub-folders such as mapped/* are included)")
    ap.add_argument("--nep", required=True, help="NEP model file (nep.txt)")
    ap.add_argument("--T", type=float, required=True, help="temperature in K")
    ap.add_argument("--glob", default="POSCAR*", help="file-name pattern (default POSCAR*)")
    ap.add_argument("--supercell", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                    help="phonon supercell (default: automatic from --min-length)")
    ap.add_argument("--min-length", type=float, default=10.0,
                    help="automatic supercell: every width >= this many Å (default 10)")
    ap.add_argument("--mesh", type=int, default=100, help="phonopy q-point mesh density (default 100)")
    ap.add_argument("--relax", action="store_true", help="relax with the NEP before phonons")
    ap.add_argument("--fmax", type=float, default=1e-4, help="relaxation force threshold, eV/Å")
    ap.add_argument("--fixed-cell", action="store_true", help="relax positions only, keep the cell")
    ap.add_argument("--steps", type=int, default=2000, help="max relaxation steps (default 2000)")
    ap.add_argument("--files-from", metavar="LIST",
                    help="only the structures listed in LIST (one path per line, relative to FOLDER)")
    ap.add_argument("--out", help="CSV file (default <folder>/free_energy_<T>K.csv)")
    ap.add_argument("--skip-dirs", nargs="*", default=["relaxed"], help="sub-folders to ignore")
    ap.add_argument("--jobs", type=int, default=None,
                    help="structures computed in parallel (default: number of cores // --threads)")
    ap.add_argument("--threads", type=int, default=1, help="BLAS/OpenMP threads per job (default 1)")
    a = ap.parse_args()
    jobs = a.jobs or max(1, (os.cpu_count() or 1) // a.threads)

    folder = Path(a.folder)
    if not folder.is_dir():
        sys.exit(f"{folder} is not a folder")
    if not Path(a.nep).is_file():
        sys.exit(f"NEP file {a.nep} not found")
    T_tag = f"{a.T:g}"
    csv_path = Path(a.out) if a.out else folder / f"free_energy_{T_tag}K.csv"

    if a.files_from:
        files = []
        for line in Path(a.files_from).read_text().split("\n"):
            line = line.strip()
            if line and not line.startswith("#"):
                p = folder / line
                if not p.is_file():
                    sys.exit(f"{a.files_from}: {p} not found")
                files.append(p)
    else:
        files = find_structures(folder, a.glob, set(a.skip_dirs))
    done = {r["file"] for r in read_done(csv_path) if r["status"] == "ok"}     # failed ones are retried
    todo = [p for p in files if str(p.relative_to(folder)) not in done]
    jobs = min(jobs, max(1, len(todo)))
    print(f"{len(files)} structures in {folder} (incl. sub-folders), {len(files) - len(todo)} already done in "
          f"{csv_path.name}, {len(todo)} to do;  T = {T_tag} K,  relax = {a.relax},  "
          f"{jobs} parallel jobs x {a.threads} thread(s)")

    common = dict(folder=str(folder), nep=a.nep, T=a.T, supercell=a.supercell, min_length=a.min_length,
                  mesh=a.mesh, relax=a.relax, fmax=a.fmax, fixed_cell=a.fixed_cell, steps=a.steps)
    work = [dict(common, path=str(p)) for p in todo]

    def report(n, row):
        append_row(csv_path, row)                                   # only the main process writes the CSV
        print(f"[{n}/{len(todo)}] {row['file']:<45} {row['phase']:<7} {row['WC']:<26} {describe(row)}  "
              f"({row['time_s']} s)", flush=True)

    t_start = time.time()
    if jobs == 1:
        for n, job in enumerate(work, 1):
            report(n, run_one(job))
    else:
        # limit threads in the workers BEFORE they start (they read these when importing numpy / NEP)
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[var] = str(a.threads)
        ctx = mp.get_context("spawn")                                # fresh interpreters: safe with OpenMP
        with ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
            futures = {ex.submit(run_one, job): job for job in work}
            for n, fut in enumerate(as_completed(futures), 1):
                job = futures[fut]
                try:
                    row = fut.result()
                except Exception as exc:                             # worker died (e.g. out of memory)
                    p = Path(job["path"])
                    row = {"file": str(p.relative_to(folder)), **parse_title(p), "status": "FAILED",
                           "note": f"worker crashed: {type(exc).__name__}: {exc}"[:300], "T_K": a.T,
                           "time_s": ""}
                report(n, row)
    if todo:
        print(f"\nwall time {time.time() - t_start:.0f} s for {len(todo)} structures")

    rows = finalize(csv_path)
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"\n{len(ok)} ok, {len(rows) - len(ok)} failed -> {csv_path}")
    if ok:
        print("\nlowest F_total per phase (meV/atom above the overall lowest of the same composition):")
        seen = {}
        for r in ok:
            key = (r["composition"], r["phase"])
            if key not in seen:
                seen[key] = r
        for (comp, ph), r in sorted(seen.items(), key=lambda kv: float(kv[1]["dF_total_meV_atom"])):
            print(f"  {comp:<14} {ph:<8} {float(r['dF_total_meV_atom']):8.2f}   {r['WC']:<26} {r['file']}")


if __name__ == "__main__":
    main()
