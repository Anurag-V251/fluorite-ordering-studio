#!/usr/bin/env python3
"""
test_single_runs.py - checks of scripts/single_runs.py, scripts/plot_harmonic_compare.py and the Single runs
endpoints (mock calorine = Einstein crystal, no NEP; needs ase, phonopy, matplotlib, pymatgen).
    python3 tests/test_single_runs.py            (~3 min; work files in tests/work_single/)

  1. physics: the harmonic F_vib, U_vib, ZPE, S, Cv of the mock (Einstein crystal, omega = sqrt(k/m) per atom)
     agree with the analytic formulas at every T of a ramp; U = F + TS; one T == the same T of a ramp ==
     batch_free_energy.py (the campaign's numbers)
  2. phases and shapes: refs c, t, oIII, m and their 2x2x2, 1x2x4, 3x2x2 repeats (with Ce on some cation sites)
     give the same F per atom for the same composition; relaxation converges and keeps the cell shape type
  3. flags: imaginary modes (MOCK_NEP_SOFT), unconverged relaxation (--steps 1), unrelaxed input; restart skips
     finished structures, --force / changed settings recompute
  4. plots: dF vs a structure and vs the lowest at each T (unflagged), --include/--exclude, --T-window, single T,
     every figure format, interactive HTML; XRD with references, XRD only without NEP for another element (Hf)
  5. app: /api/single/{defaults,refs,start (preview, local, SLURM write, SLURM submit with a fake sbatch),list,
     result,replot}; bad input refused; the run folder reproduces the run (bash run.sh)
"""
import csv
import json
import math
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SR = ROOT / "scripts" / "single_runs.py"
PC = ROOT / "scripts" / "plot_harmonic_compare.py"
W = ROOT / "tests" / "work_single"
PY = sys.executable
MOCK = str(ROOT / "tests" / "mock_calorine")
ENV = dict(os.environ, PYTHONPATH=MOCK, OMP_NUM_THREADS="1", PYTHONWARNINGS="ignore")
FAILS = []
K = {"Zr": 10.0, "Ce": 8.0, "O": 5.0}           # stiffness of the mock (tests/mock_calorine/calorine/calculators.py)


def check(cond, msg):
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def run(script, args, env=None, cwd=None):
    r = subprocess.run([PY, str(script), *map(str, args)], capture_output=True, text=True, env=env or ENV, cwd=cwd or W)
    if r.returncode:
        print(r.stdout[-2500:], r.stderr[-2500:])
    return r


def rcsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def make_structures():
    from ase.io import read, write
    d = W / "inp"
    d.mkdir(parents=True, exist_ok=True)
    out = {}
    for ph in ("c", "t", "oIII", "m"):
        a = read(str(ROOT / "refs" / f"POSCAR_{ph}"), format="vasp")
        for rep in [(1, 1, 1), (2, 2, 2), (1, 2, 4), (3, 2, 2)]:
            b = a.repeat(rep)
            s = b.get_chemical_symbols()
            cat = [i for i, x in enumerate(s) if x != "O"]
            for i in cat[::4]:                      # a quarter of the cations -> Ce (same fraction in every cell)
                s[i] = "Ce"
            b.set_chemical_symbols(s)
            f = d / f"POSCAR_{ph}_{'x'.join(map(str, rep))}"
            write(str(f), b, format="vasp", direct=True, sort=False)
            lines = f.read_text().split("\n")
            lines[0] = f"Ce_25.00 | id={len(out)} phase={ph}"
            f.write_text("\n".join(lines))
            out[(ph, rep)] = f
    (W / "nep.txt").write_text("mock model\n")
    return out


def einstein(symbols, masses, T):
    """F_vib, U_vib, S (J/mol K), Cv per atom of independent 3D oscillators, omega = sqrt(k/m)"""
    hbar, kB = 6.582119569e-16, 8.617333262e-5
    w = np.array([math.sqrt(K[s] * 1.602176634e-19 / 1e-20 / (m * 1.66053906660e-27)) for s, m in zip(symbols, masses)])
    e = hbar * w                                                    # eV
    F = 3 * (e / 2).mean()
    U = F
    S = Cv = 0.0
    if T > 0:
        x = e / (kB * T)
        F = 3 * (e / 2 + kB * T * np.log1p(-np.exp(-x))).mean()
        U = 3 * (e / 2 + e / np.expm1(x)).mean()
        S = 3 * (kB * (x / np.expm1(x) - np.log1p(-np.exp(-x)))).mean() * 96485.33212
        Cv = 3 * (kB * x ** 2 * np.exp(x) / np.expm1(x) ** 2).mean() * 96485.33212
    return F, U, S, Cv


# ============================================================================ 1-3 script
def physics(st):
    from ase.io import read
    print("== 1. harmonic thermodynamics = analytic Einstein crystal")
    f = st[("m", (1, 1, 1))]
    r = run(SR, [f, "--nep", W / "nep.txt", "--relax", "--T-range", 0, 1500, 100, "--T", 300, 333, "--out-dir", W / "ramp",
                 "--min-length", 6, "--mesh", 10, "--tasks", "relax", "harmonic"])
    check(r.returncode == 0, "ramp run finished")
    lab = f.name
    th = rcsv(W / "ramp" / lab / "thermo.csv")
    a = read(str(W / "ramp" / lab / "relaxed.vasp"), format="vasp")
    Ts = [float(x["T_K"]) for x in th]
    check(Ts[:3] == [0.0, 100.0, 200.0] and 300.0 in Ts and 333.0 in Ts and Ts[-1] == 1500.0 and len(Ts) == 17,
          f"grid = ramp + extra T ({len(Ts)} temperatures)")
    worst = 0.0
    for x in th:
        F, U, S, Cv = einstein(a.get_chemical_symbols(), a.get_masses(), float(x["T_K"]))
        worst = max(worst, abs(float(x["F_vib_eV_atom"]) - F), abs(float(x["U_vib_eV_atom"]) - U),
                    abs(float(x["entropy_J_molK_atom"]) - S) / 96485.0, abs(float(x["Cv_J_molK_atom"]) - Cv) / 96485.0)
    check(worst < 2e-6, f"F_vib, U_vib, S, Cv match the analytic values at every T (worst {worst:.1e} eV/atom)")
    res = json.loads((W / "ramp" / lab / "result.json").read_text())
    F0 = einstein(a.get_chemical_symbols(), a.get_masses(), 0)[0]
    check(abs(res["ZPE_eV_atom"] - F0) < 1e-7 and abs(float(th[0]["F_vib_eV_atom"]) - F0) < 1e-7, "ZPE = F_vib(0 K)")
    tsd = max(abs(float(x["F_total_eV_atom"]) + float(x["TS_eV_atom"]) - float(x["E_static_eV_atom"]) - float(x["U_vib_eV_atom"])) for x in th)
    check(tsd < 1e-9, "F_total + TS = E_static + U_vib")
    check(res["relax_converged"] is True and res["max_force_eV_A"] <= 1e-4 and not res["flags"], "relaxation converged, no flags")

    r = run(SR, [f, "--nep", W / "nep.txt", "--relax", "--T", 333, "--out-dir", W / "one", "--min-length", 6, "--mesh", 10,
                 "--tasks", "relax", "harmonic"])
    one = rcsv(W / "one" / lab / "thermo.csv")
    ramp333 = [x for x in th if float(x["T_K"]) == 333.0][0]
    check(r.returncode == 0 and len(one) == 1 and abs(float(one[0]["F_total_eV_atom"]) - float(ramp333["F_total_eV_atom"])) < 1e-10,
          "one temperature == the same T of the ramp")
    bf = W / "bfe"
    bf.mkdir(exist_ok=True)
    shutil.copy(f, bf / "POSCAR_x")
    r = run(ROOT / "scripts" / "batch_free_energy.py", [bf, "--nep", W / "nep.txt", "--T", 333, "--relax", "--min-length", 6,
                                                      "--mesh", 10, "--jobs", 1])
    b = rcsv(bf / "free_energy_333K.csv")[0]
    check(abs(float(b["F_total_eV_atom"]) - float(one[0]["F_total_eV_atom"])) < 1e-9 and
          abs(float(b["entropy_J_molK_atom"]) - float(one[0]["entropy_J_molK_atom"])) < 1e-6,
          "same F_total and S as batch_free_energy.py (campaign numbers)")
    summ = rcsv(W / "one" / "summary.csv")[0]
    check(summ["T_K"] == "333" and summ["label"] == lab and summ["formula"] == "CeZr3O8", "summary.csv row (T, label, formula)")

    print("== 2. four phases x four cell shapes")
    files = [st[k] for k in sorted(st)]
    r = run(SR, [*files, "--nep", W / "nep.txt", "--relax", "--T", 300, 900, "--out-dir", W / "phases", "--min-length", 6,
                 "--mesh", 8, "--tasks", "relax", "harmonic", "plot", "--formats", "png", "--reference", "POSCAR_c_1x1x1"])
    check(r.returncode == 0, f"{len(files)} structures (c/t/oIII/m x 1x1x1, 2x2x2, 1x2x4, 3x2x2) finished")
    rows = {x["label"]: x for x in rcsv(W / "phases" / "summary.csv")}
    # refs/POSCAR_c is CeO2 (Ce4O8), t/oIII/m are ZrO2: a quarter of the cations -> Ce gives CeO2 and Ce1Zr3O8
    groups = {"CeO2": [x for x in rows.values() if x["phase"] == "c"], "CeZr3O8": [x for x in rows.values() if x["phase"] != "c"]}
    for g, xs in groups.items():
        Fs = [float(x["F_total_eV_atom"]) for x in xs]
        check(len(xs) == (4 if g == "CeO2" else 12) and max(Fs) - min(Fs) < 1e-8,
              f"{g}: same F per atom in every phase and cell shape ({len(xs)} structures, spread {1000 * (max(Fs) - min(Fs)):.1e} meV/atom)")
    check(all(x["relax_converged"] == "True" and x["n_imaginary_modes"] == "0" for x in rows.values()), "all relaxed, no imaginary modes")
    a0 = read(str(st[("t", (1, 2, 4))]), format="vasp")
    a1 = read(str(W / "phases" / "POSCAR_t_1x2x4" / "relaxed.vasp"), format="vasp")
    check(np.allclose(a0.cell.angles(), a1.cell.angles(), atol=1e-3) and a0.get_chemical_symbols() == a1.get_chemical_symbols(),
          "relaxed t 1x2x4 keeps its cell angles and atom order")
    cmp_ = rcsv(W / "phases" / "plots" / "comparison.csv")
    check(len(cmp_) == 32 and all(abs(float(x["dF_meV_atom"])) < 1e-3 for x in cmp_ if x["label"].startswith("POSCAR_c"))
          and cmp_[0]["reference"] == "POSCAR_c_1x1x1",
          "comparison.csv: 16 structures x 2 T, dF vs the chosen structure")
    check(rows["POSCAR_oIII_1x1x1"]["phase"] == "oIII" and rows["POSCAR_m_3x2x2"]["composition"] == "Ce_25.00",
          "phase and composition read from the POSCAR title")

    print("== 3. flags and restart")
    soft = dict(ENV, MOCK_NEP_SOFT="1")
    r = run(SR, [st[("c", (1, 1, 1))], "--nep", W / "nep.txt", "--T", 300, "--out-dir", W / "soft", "--min-length", 6,
                 "--mesh", 8, "--tasks", "harmonic", "--T", 333], env=soft)
    res = json.loads((W / "soft" / "POSCAR_c_1x1x1" / "result.json").read_text())
    check(r.returncode == 0 and "imaginary_modes" in res["flags"] and res["n_imaginary_modes"] > 0, "imaginary modes flagged")
    check("not_relaxed" in res["flags"], "unrelaxed input (forces) flagged")
    r = run(SR, [st[("m", (1, 1, 1))], "--nep", W / "nep.txt", "--relax", "--steps", 1, "--T", 300, "--out-dir", W / "unconv",
                 "--min-length", 6, "--mesh", 8, "--tasks", "relax", "harmonic"])
    res = json.loads((W / "unconv" / "POSCAR_m_1x1x1" / "result.json").read_text())
    check(res["relax_converged"] is False and "relax_not_converged" in res["flags"] and "NOT converged" in r.stdout,
          "unconverged relaxation flagged")
    t0 = (W / "one" / lab / "thermo.csv").stat().st_mtime
    r = run(SR, [f, "--nep", W / "nep.txt", "--relax", "--T", 333, "--out-dir", W / "one", "--min-length", 6, "--mesh", 10,
                 "--tasks", "relax", "harmonic"])
    check("already done" in r.stdout and (W / "one" / lab / "thermo.csv").stat().st_mtime == t0, "restart skips a finished structure")
    r = run(SR, [f, "--nep", W / "nep.txt", "--relax", "--T", 333, 400, "--out-dir", W / "one", "--min-length", 6, "--mesh", 10,
                 "--tasks", "relax", "harmonic"])
    check("already done" not in r.stdout and len(rcsv(W / "one" / lab / "thermo.csv")) == 2, "changed temperatures recompute")
    r = run(SR, [f, "--nep", W / "nep.txt", "--relax", "--T", 333, 400, "--out-dir", W / "one", "--min-length", 6, "--mesh", 10,
                 "--tasks", "relax", "harmonic", "--force"])
    check("already done" not in r.stdout, "--force recomputes")
    r = run(SR, [f, "--T", 300, "--out-dir", W / "x"])
    check(r.returncode != 0 and "--nep" in r.stderr, "harmonic without a NEP model refused")
    r = run(SR, [f, "--nep", W / "nep.txt", "--out-dir", W / "x"])
    check(r.returncode != 0 and "--T" in r.stderr, "harmonic without a temperature refused")


# ============================================================================ 4 plots / XRD
def plots(st):
    from ase.io import read, write
    print("== 4. comparison plots and XRD")
    # three structures with different F: change the Ce count (Einstein F depends on the species)
    d = W / "cmp_in"
    d.mkdir(exist_ok=True)
    fs = []
    for n_ce, ph in [(0, "t"), (1, "oIII"), (2, "m")]:
        a = read(str(ROOT / "refs" / f"POSCAR_{ph}"), format="vasp")
        s = a.get_chemical_symbols()
        for i in range(n_ce):
            s[i] = "Ce"
        a.set_chemical_symbols(s)
        p = d / f"POSCAR_{ph}"
        write(str(p), a, format="vasp", direct=True, sort=False)
        fs.append(p)
    r = run(SR, [*fs, "--nep", W / "nep.txt", "--relax", "--T-range", 0, 1200, 100, "--out-dir", W / "cmp", "--min-length", 6,
                 "--mesh", 8, "--tasks", "relax", "harmonic", "xrd", "plot", "--xrd-ref", "c", "t", "oIII", "m",
                 "--names", "zro2_t", "ce1_oIII", "ce2_m", "--reference", "zro2_t", "--at-T", 300, 900])
    check(r.returncode == 0, "3 structures: relax + ramp + XRD + plots")
    P = W / "cmp" / "plots"
    want = [f"{s}.{e}" for s in ("F_vs_T", "thermo_panels", "dF_at_T", "phonon_dos", "xrd_compare") for e in ("png", "pdf", "svg")]
    miss = [w for w in want if not (P / w).is_file()]
    check(not miss and (P / "harmonic_compare.html").is_file() and (P / "comparison.csv").is_file(),
          "figures PNG/PDF/SVG, comparison.csv and HTML written" + (f" (missing {miss})" if miss else ""))
    cmp_ = rcsv(P / "comparison.csv")
    by = {(x["label"], x["T_K"]): float(x["F_total_eV_atom"]) for x in cmp_}
    ok = all(abs(float(x["dF_meV_atom"]) - 1000 * (by[(x["label"], x["T_K"])] - by[("zro2_t", x["T_K"])])) < 1e-6 for x in cmp_)
    check(ok and all(abs(float(x["dF_meV_atom"])) < 1e-9 for x in cmp_ if x["label"] == "zro2_t"), "dF = F - F(reference) at every T")
    check(any(float(x["dF_meV_atom"]) < -1 for x in cmp_ if x["label"] == "ce2_m" and x["T_K"] == "1200"),
          "softer Ce-rich structure is lower in F at high T (Einstein: lower frequencies)")
    r = run(PC, [W / "cmp", "--reference", "lowest", "--out-dir", W / "cmp" / "p2", "--formats", "png"])
    c2 = rcsv(W / "cmp" / "p2" / "comparison.csv")
    for T in ("0", "600", "1200"):
        m = min(float(x["F_total_eV_atom"]) for x in c2 if x["T_K"] == T)
        low = [x for x in c2 if x["T_K"] == T and abs(float(x["F_total_eV_atom"]) - m) < 1e-12][0]
        check(abs(float(low["dF_meV_atom"])) < 1e-9 and all(float(x["dF_meV_atom"]) >= -1e-9 for x in c2 if x["T_K"] == T),
              f"reference 'lowest' at {T} K: lowest is 0, the others >= 0")
    r = run(PC, [W / "cmp", "--include", "ce1_oIII", "ce2_m", "--reference", "ce1_oIII", "--T-window", 200, 800,
                 "--out-dir", W / "cmp" / "p3", "--formats", "svg"])
    c3 = rcsv(W / "cmp" / "p3" / "comparison.csv")
    check(r.returncode == 0 and {x["label"] for x in c3} == {"ce1_oIII", "ce2_m"} and {x["T_K"] for x in c3} == {"200", "300", "400", "500", "600", "700", "800"},
          "--include, --T-window")
    r = run(PC, [W / "cmp", "--reference", "nope", "--out-dir", W / "cmp" / "p4"])
    check(r.returncode != 0 and "no such structure" in r.stderr, "unknown reference refused with the list of labels")
    r = run(PC, [W / "one", "--out-dir", W / "one" / "p", "--formats", "png"])   # one structure, two T
    check(r.returncode == 0 and (W / "one" / "p" / "F_vs_T.png").is_file(), "one structure works")
    # flagged series: lowest uses the unflagged ones only
    shutil.copytree(W / "soft" / "POSCAR_c_1x1x1", W / "cmp_flag" / "soft_c", dirs_exist_ok=True)
    shutil.copytree(W / "one" / lab_of(W / "one"), W / "cmp_flag" / "good", dirs_exist_ok=True)
    r = run(PC, [W / "cmp_flag" / "soft_c" / "thermo.csv", W / "cmp_flag" / "good" / "thermo.csv", "--out-dir", W / "cmp_flag" / "p",
                 "--formats", "png", "--at-T", 300])
    c5 = rcsv(W / "cmp_flag" / "p" / "comparison.csv")
    g = [x for x in c5 if x["label"] == "good" and x["T_K"] == "333"]
    check(r.returncode == 0 and "FLAGGED soft_c" in r.stdout and any(x["flags"] for x in c5) and g and abs(float(g[0]["dF_meV_atom"])) < 1e-9,
          "flagged structure marked; 'lowest' reference ignores it")
    html = (W / "cmp" / "plots" / "harmonic_compare.html").read_text()
    check("Fixed-cell harmonic" in html and '"zro2_t"' in html and "<script src" not in html, "HTML self-contained, method note")
    xr = rcsv(W / "cmp" / "ce2_m" / "xrd" / "xrd.csv")
    check({"ref_c_max100", "ref_t_max100", "ref_oIII_max100", "ref_m_max100"} <= set(xr[0]) and
          json.loads((W / "cmp" / "ce2_m" / "xrd" / "xrd.json").read_text())["structures"][0]["path"].endswith("relaxed.vasp")
          if "structures" in json.loads((W / "cmp" / "ce2_m" / "xrd" / "xrd.json").read_text()) else
          {"ref_c_max100", "ref_m_max100"} <= set(xr[0]), "XRD of the relaxed structure with the four references")
    # XRD only, no NEP, another element
    a = read(str(ROOT / "refs" / "POSCAR_m"), format="vasp")
    a.set_chemical_symbols(["Hf" if x == "Zr" else x for x in a.get_chemical_symbols()])
    write(str(d / "POSCAR_HfO2_m"), a, format="vasp", direct=True, sort=False)
    r = run(SR, [d / "POSCAR_HfO2_m", "--tasks", "xrd", "plot", "--xrd-ref", "m", "--out-dir", W / "hf", "--formats", "png"])
    check(r.returncode == 0 and (W / "hf" / "POSCAR_HfO2_m" / "xrd" / "xrd.png").is_file() and (W / "hf" / "plots" / "xrd_compare.png").is_file(),
          "XRD only (no NEP) for HfO2: nothing hard-coded to Zr/Ce")


def lab_of(d):
    return [p.name for p in d.iterdir() if (p / "thermo.csv").is_file()][0]


# ============================================================================ 5 app
def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def app(st):
    print("== 5. app endpoints")
    sd = W / "state"
    (W / "bin").mkdir(exist_ok=True)
    (W / "bin" / "sbatch").write_text("#!/bin/sh\necho Submitted batch job 4242\n")
    os.chmod(W / "bin" / "sbatch", 0o755)
    env = dict(ENV, PATH=f"{W / 'bin'}:{os.environ['PATH']}")
    port = free_port()
    srv = subprocess.Popen([PY, str(ROOT / "start_studio.py"), "--port", str(port), "--state", str(sd), "--dir", str(W)],
                           env=env, stdout=open(W / "server.log", "w"), stderr=subprocess.STDOUT)
    try:
        for _ in range(60):
            if (sd / "token").is_file():
                break
            time.sleep(0.2)
        time.sleep(1)
        tok = (sd / "token").read_text().strip()

        def call(path, body=None):
            req = urllib.request.Request(f"http://127.0.0.1:{port}" + path, data=None if body is None else json.dumps(body).encode(),
                                         headers={"X-Token": tok, "X-Requested-With": "studio", "Content-Type": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=120) as f:
                    return json.loads(f.read())
            except urllib.error.HTTPError as e:
                return json.loads(e.read())

        def wait(jid, t=240):
            for _ in range(t * 2):
                j = call(f"/api/jobs?id={jid}&tail=20000")["job"]
                if j["state"] in ("done", "failed", "cancelled", "lost"):
                    return j
                time.sleep(0.5)
            return j

        info = call("/api/info")
        check("22_single.js" in info["pages"], "page registered")
        mods = {m["id"] for m in call("/api/modules")["modules"]}
        check({"harmonic_compare", "single_runs_cli"} <= mods, "plot + tool modules registered")
        for mid in ("harmonic_compare", "single_runs_cli"):
            a = call(f"/api/argspec?module={mid}")
            check(a.get("ok") and len(a.get("actions", [])) > 5, f"form of {mid} built from argparse")
        d = call("/api/single/defaults")
        check(d["ok"] and set(d["refs"]) >= {"c", "t", "oIII", "m"} and d["base"].endswith("single_runs"), "defaults: refs and base")
        check(call("/api/single/refs?dir=" + str(ROOT / "refs"))["refs"] == sorted(["c", "m", "oIII", "t"]), "refs of a folder")
        base = str(W / "app_runs")
        ps = [str(st[("c", (2, 2, 2))]), str(st[("oIII", (1, 1, 1))])]
        spec = {"base": base, "name": "two phases", "poscars": ps, "nep": str(W / "nep.txt"),
                "tasks": ["relax", "harmonic", "xrd", "plot"], "T_range": [0, 1000, 100], "T": ["300"], "min_length": 6, "mesh": 8,
                "fmax": 1e-4, "steps": 500, "xrd_ref": ["c", "oIII"], "refs_dir": str(ROOT / "refs"), "reference": "POSCAR_oIII_1x1x1",
                "at_T": ["300"], "formats": ["png", "svg"]}
        r = call("/api/single/start", {**spec, "mode": "preview"})
        check(r["ok"] and "--T-range 0 1000 100" in r["cmd"] and "--xrd-ref c oIII" in r["cmd"], "preview command")
        for bad, msg in [({"tasks": []}, "task"), ({"poscars": []}, "POSCAR"), ({"nep": str(W / "missing.txt")}, "not found"),
                         ({"T_range": [500, 100, 10], "T": []}, "TMIN"), ({"T_range": None, "T": []}, "temperature"),
                         ({"T": ["abc"]}, "number"), ({"names": ["a", "a"]}, "different")]:
            r = call("/api/single/start", {**spec, **bad, "mode": "local"})
            check(not r["ok"] and msg in r["error"], f"refused: {msg} ({r.get('error', '')[:60]})")
        check(not (W / "app_runs").exists() or not any((W / "app_runs").iterdir()), "nothing written for refused requests")
        r = call("/api/single/start", {**spec, "mode": "local"})
        check(r["ok"] and r.get("job"), "local run started")
        j = wait(r["job"])
        check(j["state"] == "done", f"local run done ({j['state']})")
        rd = Path(r["folder"])
        res = call("/api/single/result?folder=" + str(rd))
        check(res["ok"] and res["status"]["state"] == "done" and len(res["summary"]) == 2, "result: state done, 2 summary rows")
        check(any(p.endswith("F_vs_T.png") for p in res["plots"]) and any(p.endswith("xrd_compare.svg") for p in res["plots"]),
              "result lists the figures")
        rq = json.loads((rd / "request.json").read_text())
        check(len(rq["inputs"]) == 2 and rq["nep"]["md5"] and set(rq["scripts"]) == {"single_runs.py", "batch_free_energy.py",
              "xrd_average.py", "plot_harmonic_compare.py"}, "request.json: inputs, NEP md5, frozen scripts")
        lst = call("/api/single/list?base=" + base)
        check(lst["ok"] and lst["runs"][0]["folder"] == str(rd) and lst["runs"][0]["state"] == "done", "list shows the run")
        r = call("/api/single/replot", {"folder": str(rd), "reference": "lowest", "formats": ["pdf"], "at_T": ["300", "700"],
                                        "T_window": ["100", "900"], "width": 7, "font_size": 11})
        j = wait(r["job"])
        cc = rcsv(rd / "results" / "plots" / "comparison.csv")
        check(j["state"] == "done" and {x["T_K"] for x in cc} == {str(t) for t in range(100, 1000, 100)} and cc[0]["reference"].startswith("lowest"),
              "re-plot with reference 'lowest' and a T window")
        before = {x["label"]: x["F_total_eV_atom"] for x in rcsv(rd / "results" / "summary.csv")}
        shutil.rmtree(rd / "results")
        rr = subprocess.run(["bash", str(rd / "run.sh")], capture_output=True, text=True, env=ENV)
        after = {x["label"]: x["F_total_eV_atom"] for x in rcsv(rd / "results" / "summary.csv")}
        check(rr.returncode == 0 and before == after, "bash run.sh reproduces the run (same F)")
        r = call("/api/single/start", {**spec, "mode": "slurm_write", "slurm": {"partition": "cpu", "time": "02:00:00", "cpus": 4,
                                                                                "setup": "module load Python", "extra": "#SBATCH --constraint=x"}})
        sb = Path(r["sbatch"]).read_text()
        check(r["ok"] and "#SBATCH --partition=cpu" in sb and "#SBATCH --cpus-per-task=4" in sb and "module load Python" in sb
              and "scripts/single_runs.py" in sb and "#SBATCH --constraint=x" in sb, "SLURM script written")
        r = call("/api/single/start", {**spec, "mode": "slurm_write", "slurm": {"partition": "cpu; rm -rf ~"}})
        check(not r["ok"] and "not allowed" in r["error"], "SLURM option with shell characters refused")
        r = call("/api/single/start", {**spec, "mode": "slurm_submit", "slurm": {"cpus": 2}})
        j = wait(r["job"])
        check(j["state"] == "done" and "4242" in j.get("slurm_ids", []), "submitted with sbatch, SLURM id picked up")
        r = call("/api/single/start", {**spec, "tasks": ["xrd"], "nep": "", "mode": "local", "name": "xrd only"})
        j = wait(r["job"])
        rq = json.loads((Path(r["folder"]) / "request.json").read_text())
        check(j["state"] == "done" and rq["nep"] is None, "XRD-only run without a NEP model")
    finally:
        srv.terminate()


def main():
    shutil.rmtree(W, ignore_errors=True)
    W.mkdir(parents=True)
    st = make_structures()
    t0 = time.time()
    physics(st)
    plots(st)
    app(st)
    print(f"\n{time.time() - t0:.0f} s")
    if FAILS:
        print(f"{len(FAILS)} FAILED:\n  " + "\n  ".join(FAILS))
        sys.exit(1)
    print("SINGLE RUNS TESTS PASSED")


if __name__ == "__main__":
    main()
