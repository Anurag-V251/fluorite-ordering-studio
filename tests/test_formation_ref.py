#!/usr/bin/env python3
"""
test_formation_ref.py - checks --ref formation (--ends, --end-x) of plot_F_vs_composition.py and plot_F_vs_x_html.py (~1 min, no NEP).

Mock merged CSVs F = (1-x) F0[phase] + x F1[phase] + W x(1-x) + noise for the four phases (parent = c, t, oIII, m)
and three cell sizes (8, 32, 48 cations).  Checks: formation F of the end members is 0 (--ends same) or the end
member difference (--ends lowest / named phases); main peaks follow the expected curve; panel (b) dF and the old
references are unchanged by the new option; missing or flagged end members stop with a clear message.

    python tests/test_formation_ref.py          (from the package folder; scripts/ next to tests/)
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
F0 = {"parent": -9.30, "t": -9.36, "oIII": -9.37, "m": -9.38}      # x = 0: m lowest
F1 = {"parent": -8.20, "t": -8.195, "oIII": -8.19, "m": -8.17}     # x = 1: parent (c) lowest
W = {"parent": 0.060, "t": 0.050, "oIII": 0.055, "m": 0.070}
ok = True


def check(cond, msg):
    global ok
    print(("PASS " if cond else "FAIL ") + msg)
    ok &= bool(cond)


def mock(path, nhost, nord, phases, drop=(), flag_t0=False, seed=1):
    rng = np.random.default_rng(seed)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n_dopant", "x_dopant", "phase", "status", "T_K", "F_total_eV_atom", "n_imaginary_modes",
                    "relax_converged"])
        for n in range(nhost + 1):
            x = n / nhost
            if x in drop:
                continue
            for p in phases:
                for i in range(1 if n in (0, nhost) else nord):
                    F = (1 - x) * F0[p] + x * F1[p] + W[p] * x * (1 - x)
                    if 0 < n < nhost:
                        F += rng.normal(0, 0.0005)
                    bad = flag_t0 and p == "t" and n == 0
                    w.writerow([n, f"{x:.6f}", "" if p == "parent" else p, "ok", 300, f"{F:.7f}", 2 if bad else 0,
                                "True"])


def run(script, csvp, *args):
    return subprocess.run([sys.executable, str(SCRIPTS / script), str(csvp), "--bootstrap", "0", "--jobs", "2",
                           "--no-cache", *args], capture_output=True, text=True)


def peaks(d, stem):
    return {(r["x_dopant"], r["phase"], r["peak_no"]): r for r in csv.DictReader(open(Path(d) / f"{stem}_F_vs_x_peaks.csv"))}


def expected(x, p, e0, e1, phases):
    r0 = F0[p] if e0 == "same" else (min(F0[q] for q in phases) if e0 == "lowest" else F0[e0])
    r1 = F1[p] if e1 == "same" else (min(F1[q] for q in phases) if e1 == "lowest" else F1[e1])
    return 1000 * ((1 - x) * (F0[p] - r0) + x * (F1[p] - r1) + W[p] * x * (1 - x))


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    for nhost, phases in [(8, ["t", "oIII", "m"]), (32, ["parent", "t", "oIII", "m"]), (48, ["parent", "t", "oIII", "m"])]:
        c = tmp / f"mock{nhost}.csv"; mock(c, nhost, 12, phases)
        base = tmp / f"comp{nhost}"
        r = run("plot_F_vs_composition.py", c, "--ref", "composition", "--formats", "png", "--out-dir", str(base))
        check(r.returncode == 0, f"{nhost} cations: --ref composition runs")
        ref_b = peaks(base, c.stem)
        for e0, e1 in [("same", "same"), ("lowest", "lowest"), ("m", "parent" if "parent" in phases else "t"),
                       ("lowest", "same")]:
            od = tmp / f"f{nhost}_{e0}_{e1}"
            r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--ends", f"{e0},{e1}", "--formats", "png",
                    "--out-dir", str(od))
            check(r.returncode == 0, f"{nhost} cations, --ends {e0},{e1}: runs" + (f" ({r.stderr[-300:]})" if r.returncode else ""))
            if r.returncode:
                continue
            P = peaks(od, c.stem)
            dev = max(abs(float(v["F_peak_meV_atom_ref_formation"]) - expected(float(k[0]), k[1], e0, e1, phases))
                      for k, v in P.items() if v["main_peak"] == "True")
            check(dev < 0.6, f"   main peaks follow the expected formation curve (worst {dev:.2f} meV/atom)")
            ends = [float(v["F_peak_meV_atom_ref_formation"]) for k, v in P.items() if k[0] in ("0.000000", "1.000000")]
            if (e0, e1) == ("same", "same"):
                check(max(map(abs, ends)) < 1e-6, "   end members are exactly 0")
            db = max(abs(float(v["dF_to_lowest_peak_meV_atom"]) - float(ref_b[k]["dF_to_lowest_peak_meV_atom"]))
                     for k, v in P.items())
            check(db < 1e-9, "   panel (b) dF identical to --ref composition")
            check("end points" in r.stdout, "   reference printed")
        h = tmp / f"page{nhost}.html"
        r = run("plot_F_vs_x_html.py", c, "--ref", "formation", "--ends", "lowest", "--out", str(h))
        page = h.read_text() if h.exists() else ""
        check(r.returncode == 0 and '"ends":["lowest","lowest"]' in page and 'value="formation"' in page,
              f"{nhost} cations: HTML page written with formation preset")

    # chosen end points: other end compositions, fixed F values
    c = tmp / "mock32.csv"
    od = tmp / "ex"
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--ends", "lowest", "--end-x", "0", "50",
            "--formats", "png", "--out-dir", str(od))
    check(r.returncode == 0 and "line through F(x=0%) and F(x=50%)" in r.stdout, "--end-x 0 50: runs, says which line")
    if r.returncode == 0:
        P = peaks(od, c.stem)
        low50 = min(float(v["lowest_meV_atom"]) for k, v in P.items() if k[0] == "0.500000")
        low0 = min(float(v["lowest_meV_atom"]) for k, v in P.items() if k[0] == "0.000000")
        check(abs(low50) < 1e-6 and abs(low0) < 1e-6, "   lowest ordering at both end compositions is exactly 0")
    r1 = run("plot_F_vs_composition.py", c, "--ref", "formation", "--ends", "m,parent", "--formats", "png",
             "--out-dir", str(tmp / "named"))
    r2 = run("plot_F_vs_composition.py", c, "--ref", "formation", f"--ends={F0['m']},{F1['parent']}", "--formats",
             "png", "--out-dir", str(tmp / "fixed"))
    check(r1.returncode == 0 and r2.returncode == 0 and "fixed F" in r2.stdout and
          peaks(tmp / "named", c.stem) == peaks(tmp / "fixed", c.stem),
          "fixed end F values give the same result as the phases with those F")
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--ends=same,-8.2", "--end-x", "0", "50",
            "--formats", "png", "--out-dir", str(tmp / "mixfix"))
    check(r.returncode == 0, "same phase at one end, fixed F at the other, --end-x 0 50")
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--end-x", "0", "41", "--formats", "png",
            "--out-dir", str(tmp / "e"))
    check(r.returncode != 0 and "nearest" in r.stderr, "--end-x at a composition not in the CSV: error lists nearest")
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--end-x", "50", "20", "--formats", "png",
            "--out-dir", str(tmp / "e"))
    check(r.returncode != 0 and "below" in r.stderr, "--end-x in the wrong order: error")
    h = tmp / "page_ex.html"
    r = run("plot_F_vs_x_html.py", c, "--ref", "formation", "--ends=lowest,-8.2", "--end-x", "0", "50", "--out", str(h))
    page = h.read_text() if h.exists() else ""
    check(r.returncode == 0 and '"ends":["lowest","value"]' in page and '"end_val":[null,-8.2]' in page
          and '"end_x":[0.0,50.0]' in page, "HTML page: end x and fixed F preset")

    # end members missing / flagged
    c = tmp / "no1.csv"; mock(c, 16, 10, ["parent", "t", "oIII", "m"], drop=(1.0,))
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--formats", "png", "--out-dir", str(tmp / "e"))
    check(r.returncode != 0 and "no composition at 100 %" in r.stderr, "missing x = 1: clear error")
    c = tmp / "flag.csv"; mock(c, 16, 10, ["parent", "t", "oIII", "m"], flag_t0=True)
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--formats", "png", "--out-dir", str(tmp / "e"))
    check(r.returncode != 0 and "phase t" in r.stderr and "flagged" in r.stderr, "flagged t end member: error names it")
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--include-flagged", "--formats", "png",
            "--out-dir", str(tmp / "e"))
    check(r.returncode == 0, "   --include-flagged uses it")
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--ends", "lowest", "--formats", "png",
            "--out-dir", str(tmp / "e"))
    check(r.returncode == 0 and "left out flagged end points: t" in r.stdout, "   --ends lowest does not need it, warns")
    r = run("plot_F_vs_composition.py", c, "--ref", "formation", "--phases", "m", "oIII", "--formats", "png",
            "--out-dir", str(tmp / "e"))
    check(r.returncode == 0, "   not needed when t is not shown")
    for s in ("plot_F_vs_composition.py", "plot_F_vs_x_html.py"):
        r = run(s, c, "--ref", "formation", "--ends", "m,foo", "--out", str(tmp / "x.html")) if "html" in s else \
            run(s, c, "--ref", "formation", "--ends", "m,foo", "--out-dir", str(tmp / "e"))
        check(r.returncode != 0 and "foo" in r.stderr, f"{s}: unknown --ends phase rejected")

print("ALL PASSED" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
