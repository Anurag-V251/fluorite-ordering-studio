#!/usr/bin/env python3
"""
make_xrd_mock.py OUTDIR [--n 8] [--orderings 120] [--rep 2 2 2] - MOCK campaign for the XRD tool (no NEP needed)

  * the orderings and the t / oIII / m structures are made by the real pipeline (campaign.py prepare =
    enumerate_substitutions.py + batch_map.py + compare_ordering.py) from a cubic host and refs/POSCAR_*
  * the "relaxed" structures and the free energies are PLACEHOLDERS: every mapped structure gets random atomic
    displacements (0.03 A) and a small random strain, and F is drawn from made-up Gaussian mixtures per phase
    (oIII has two peaks).  Some rows are flagged (imaginary modes / unconverged); a rerun round fixes a few.
Only for showing and testing the XRD page.  Nothing here is a physical result.
"""
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = HERE.parent
sys.path.insert(0, str(APP / "scripts"))

FCC = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
TET = np.array([[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)])
# made-up F distributions, meV/atom above an arbitrary zero: phase -> [(weight, mu, sigma)]
MIX = {"parent": [(1, 9.0, 1.2)], "t": [(1, 5.0, 0.9)], "oIII": [(0.55, 2.2, 0.55), (0.45, 4.6, 0.7)],
       "m": [(1, 1.0, 0.8)]}


def write_host(path, rep, a=5.15):
    rep = np.array(rep)
    cells = np.array([[i, j, k] for i in range(rep[0]) for j in range(rep[1]) for k in range(rep[2])])
    cat = ((cells[:, None] + FCC[None]) / rep).reshape(-1, 3)
    ox = ((cells[:, None] + TET[None]) / rep).reshape(-1, 3)
    with open(path, "w") as f:
        f.write(f"ZrO2 cubic fluorite {'x'.join(map(str, rep))} (mock host)\n1.0\n")
        for v in np.diag(rep * a):
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write(f"  Zr O\n  {len(cat)} {len(ox)}\nDirect\n")
        for p in np.vstack([cat, ox]):
            f.write("  %.10f %.10f %.10f\n" % tuple(p))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out")
    ap.add_argument("--n", type=int, default=8, help="Ce atoms (default 8 of 32 = 25 %%)")
    ap.add_argument("--orderings", type=int, default=120)
    ap.add_argument("--rep", type=int, nargs=3, default=[2, 2, 2])
    ap.add_argument("--T", type=float, default=300)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    from ase.io import read, write
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_host(out / "POSCAR_host", a.rep)
    (out / "campaign.toml").write_text(f'''# MOCK campaign for the XRD tool (placeholder F, placeholder relaxation)
[structure]
poscar = "POSCAR_host"
host = "Zr"
dopant = "Ce"
compositions = [{a.n}]
max_orderings = {a.orderings}
seed = 3
[paths]
scripts = "{APP / 'scripts'}"
workdir = "campaign"
[mapping]
templates = ["t={APP / 'refs/POSCAR_t'}", "oIII={APP / 'refs/POSCAR_oIII'}", "m={APP / 'refs/POSCAR_m'}"]
include_parent = true
verify = true
jobs = 4
[free_energy]
nep = "nep.txt"
T = {a.T:g}
''')
    (out / "nep.txt").write_text("mock\n")
    subprocess.run([sys.executable, str(APP / "scripts" / "campaign.py"), "campaign.toml", "prepare"], cwd=out, check=True,
                   stdout=subprocess.DEVNULL)
    work = out / "campaign"
    import batch_free_energy as bfe
    rng = np.random.default_rng(a.seed)
    rows, flagged = [], []
    comp = next(p for p in work.iterdir() if p.is_dir() and p.name.startswith("Ce_"))
    files = sorted(str(p.relative_to(work)) for p in comp.rglob("POSCAR_*") if "relaxed" not in p.parts)
    files = [f for f in files if "/mapped/" in f or f.count("/") == 1]
    F0 = -9.21
    for rel in files:
        info = bfe.parse_title(work / rel)
        ph = info["phase"]
        mix = MIX.get(ph, MIX["t"])
        w = np.array([m[0] for m in mix], float)
        j = rng.choice(len(mix), p=w / w.sum())
        F = F0 + 1e-3 * rng.normal(mix[j][1], mix[j][2])
        at = read(str(work / rel), format="vasp")
        strain = np.eye(3) + rng.normal(0, 0.002, (3, 3))
        at.set_cell(at.cell[:] @ strain, scale_atoms=True)
        at.positions += rng.normal(0, 0.03, at.positions.shape)
        dst = work / "relaxed" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        write(str(dst), at, format="vasp", direct=True, sort=False)
        title = (work / rel).open().readline().strip() + " | relaxed with mock (PLACEHOLDER)"
        lines = dst.read_text().split("\n")
        lines[0] = title
        dst.write_text("\n".join(lines))
        u = rng.random()
        imag, unc = (2 if u < 0.03 else 0), ("False" if 0.03 <= u < 0.05 else "True")
        L = at.cell.lengths()
        ang = at.cell.angles()
        r = dict(file=rel, **info, T_K=a.T, E_static_eV_atom=F + 0.05, F_vib_eV_atom=-0.05, F_total_eV_atom=f"{F:.8f}",
                 n_atoms=len(at), volume_A3_atom=at.get_volume() / len(at), a_A=L[0], b_A=L[1], c_A=L[2],
                 alpha_deg=ang[0], beta_deg=ang[1], gamma_deg=ang[2], relaxed=True, n_imaginary_modes=imag,
                 min_freq_THz=-0.3 if imag else 0.5, status="ok", relax_converged=unc,
                 note="PLACEHOLDER F (mock for the XRD tool)")
        rows.append(r)
        if imag or unc == "False":
            flagged.append((r, at))
    (work / "results").mkdir(exist_ok=True)
    with open(work / "results" / "chunk_00000.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=bfe.CSV_COLS, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    # rerun round 1: half of the flagged rows become clean, with a new (slightly different) geometry
    rr = []
    for r, at in flagged[::2]:
        at2 = at.copy()
        at2.positions += rng.normal(0, 0.01, at2.positions.shape)
        dst = work / "rerun_001" / "relaxed" / r["file"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        write(str(dst), at2, format="vasp", direct=True, sort=False)
        rr.append({**r, "n_imaginary_modes": 0, "relax_converged": "True", "min_freq_THz": 0.4,
                   "F_total_eV_atom": f"{float(r['F_total_eV_atom']) - 2e-4:.8f}", "note": "rerun PLACEHOLDER"})
    with open(work / "results" / "chunk_rerun001_00000.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=bfe.CSV_COLS, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rr)
    print(f"mock campaign in {out}: {comp.name}, {len(rows)} structures, {len(flagged)} flagged, "
          f"{len(rr)} fixed in rerun_001 (PLACEHOLDER F and relaxation)")


if __name__ == "__main__":
    main()
