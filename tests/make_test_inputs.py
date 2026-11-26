#!/usr/bin/env python3
"""
Smoke-test inputs (synthetic, for checking the installation only - not physical results):

  make_test_inputs.py structures OUTDIR
      OUTDIR/POSCAR        2x2x2 cubic fluorite ZrO2 (32 Zr, 64 O), a = 5.10 Å
      OUTDIR/POSCAR_t      12-atom tetragonal template (pseudo-cubic setting, a = 5.07, c = 5.18, O shifted along c)
      OUTDIR/POSCAR_c      12-atom cubic template (a = 5.10)

  make_test_inputs.py fake-csv FOLDER T
      Writes FOLDER/free_energy_<T>K.csv in the batch_free_energy.py format for every POSCAR* under FOLDER
      (incl. mapped/*). Cell lengths/angles are read from the real POSCARs; the free energies are PLACEHOLDERS
      (a smooth function of phase and Warren-Cowley alpha + noise), because the NEP step needs your model.
"""
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

FCC = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
TET = np.array([[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)])


def write(path, title, lat, names, counts, frac):
    with open(path, "w") as f:
        f.write(title + "\n1.0\n")
        for v in lat:
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write("  " + " ".join(names) + "\n  " + " ".join(map(str, counts)) + "\nDirect\n")
        for p in frac:
            f.write("  %.10f %.10f %.10f\n" % tuple(p))


def structures(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    n = np.array([2, 2, 2])
    cells = np.array([[i, j, k] for i in range(2) for j in range(2) for k in range(2)])
    cat = ((cells[:, None] + FCC[None]) / n).reshape(-1, 3)
    ox = ((cells[:, None] + TET[None]) / n).reshape(-1, 3)
    write(out / "POSCAR", "ZrO2 cubic fluorite 2x2x2", np.diag(n * 5.10), ["Zr", "O"], [32, 64], np.vstack([cat, ox]))
    o_t = []
    for q in TET:
        s = 1 if (round(q[0] * 4) + round(q[1] * 4)) % 4 == 0 else -1          # t-ZrO2 checkerboard O shift
        o_t.append([q[0], q[1], q[2] + s * 0.045])
    write(out / "POSCAR_t", "synthetic t-ZrO2 template (pseudo-cubic setting)", np.diag([5.07, 5.07, 5.18]),
          ["Zr", "O"], [4, 8], np.vstack([FCC, o_t]))
    write(out / "POSCAR_c", "synthetic c-ZrO2 template", np.diag([5.10] * 3), ["Zr", "O"], [4, 8],
          np.vstack([FCC, TET]))
    print(f"wrote {out}/POSCAR, POSCAR_t, POSCAR_c")


def fake_csv(folder, T):
    import batch_free_energy as bfe
    folder = Path(folder)
    rng = np.random.default_rng(0)
    rows = []
    for root, dirs, files in os.walk(folder):
        dirs[:] = [d for d in dirs if d != "relaxed"]
        for fn in sorted(files):
            if not fn.startswith("POSCAR"):
                continue
            p = Path(root) / fn
            lines = p.read_text().split("\n")
            lat = np.array([[float(x) for x in lines[k].split()[:3]] for k in (2, 3, 4)]) * float(lines[1].split()[0])
            L = np.linalg.norm(lat, axis=1)
            ang = [np.degrees(np.arccos(np.dot(lat[i], lat[j]) / (L[i] * L[j]))) for i, j in ((1, 2), (0, 2), (0, 1))]
            info = bfe.parse_title(p)
            m = re.search(r"=([+-][\d.]+)", info["WC"])
            alpha = float(m.group(1)) if m else 0.0
            off = {"parent": 0.0, "c": 0.3, "t": -1.0}.get(info["phase"], 0.0)
            F = -9.20 + 1e-3 * (off - 6.0 * alpha + rng.normal(0, 0.4))   # PLACEHOLDER energy
            n_at = sum(int(x) for x in lines[6].split())
            rows.append(dict(file=str(p.relative_to(folder)), **info, T_K=T, E_static_eV_atom=F + 0.05,
                             F_vib_eV_atom=-0.05, F_total_eV_atom=F, n_atoms=n_at,
                             volume_A3_atom=abs(np.linalg.det(lat)) / n_at, a_A=L[0], b_A=L[1], c_A=L[2],
                             alpha_deg=ang[0], beta_deg=ang[1], gamma_deg=ang[2], relaxed=False,
                             n_imaginary_modes=0, status="ok", note="PLACEHOLDER energies (smoke test)"))
    out = folder / f"free_energy_{T:g}K.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=bfe.CSV_COLS)
        w.writeheader()
        w.writerows({k: r.get(k, "") for k in bfe.CSV_COLS} for r in rows)
    print(f"wrote {out} ({len(rows)} structures, placeholder energies)")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "structures":
        structures(sys.argv[2])
    elif len(sys.argv) >= 4 and sys.argv[1] == "fake-csv":
        fake_csv(sys.argv[2], float(sys.argv[3]))
    else:
        sys.exit(__doc__)
