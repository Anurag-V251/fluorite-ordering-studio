#!/usr/bin/env python3
"""
test_xrd.py - checks of scripts/xrd_average.py and the XRD endpoints of the app (needs ase, pymatgen; no NEP).
    python3 tests/test_xrd.py            (~2 min; work files in tests/work_xrd/)

  1. physics: refs c, t, oIII, m and their 2x2x2, 1x2x4, 3x2x2 repeats give the same pattern (--norm absolute,
     per n_atoms^2) and the same sticks; CeO2 (111) is where Bragg's law puts it, for CuKa and MoKa
  2. species: --replace-cations Zr == --replace Ce:Zr == the same geometry written with Zr; references too;
     the relaxed files are never changed (md5 + mtime)
  3. averaging: N copies of one structure -> std 0 and the single pattern; profile heights = stick intensities
  4. selection on mock campaigns (2x2x2 and 1x2x4 hosts): window, phases, flags, rerun geometry, the
     per-composition CSV, a plain folder, missing relaxed files, --max subsets; same rows as the app's endpoint
  5. app: /api/xrd/rows, /api/xrd/compute (job -> sticks JSON, cache), /api/xrd/export (files)
"""
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "xrd_average.py"
W = ROOT / "tests" / "work_xrd"
PY = sys.executable
sys.path.insert(0, str(ROOT / "scripts"))
import xrd_average as XA  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def run(args, **kw):
    r = subprocess.run([PY, str(SCRIPT), *map(str, args)], capture_output=True, text=True, cwd=kw.get("cwd", W))
    if r.returncode != 0 and not kw.get("may_fail"):
        print(r.stdout[-2000:], r.stderr[-3000:])
    return r


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def pattern(out_dir, stem):
    rows = read_csv(Path(out_dir) / f"{stem}.csv")
    keys = list(rows[0].keys())
    return np.array([[float(r[k]) for k in keys[:5]] for r in rows]), rows


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def physics():
    from ase.io import read, write
    print("== 1. the four phases and supercells of different shapes give the same pattern")
    cells = W / "cells"
    cells.mkdir(parents=True, exist_ok=True)
    for ph in ("c", "t", "oIII", "m"):
        a = read(str(ROOT / "refs" / f"POSCAR_{ph}"), format="vasp")
        ref = None
        for rep in [(1, 1, 1), (2, 2, 2), (1, 2, 4), (3, 2, 2)]:
            f = cells / f"POSCAR_{ph}_{'x'.join(map(str, rep))}"
            write(str(f), a.repeat(rep), format="vasp", direct=True, sort=False)
            r = run(["--structures", f, "--norm", "absolute", "--out-dir", cells / "out", "--stem", f.name, "--no-plot",
                     "--jobs", "1", "--min-rel", "0"])
            if r.returncode:
                check(False, f"{f.name} computed")
                continue
            y, _ = pattern(cells / "out", f.name)
            if ref is None:
                ref = y
                npk = len(XA.sticks_of(f, XA.stick_params(_ns(two_theta=[10, 90])))["x"])
                continue
            dev = np.max(np.abs(y[:, 1] - ref[:, 1])) / np.max(ref[:, 1])
            s = XA.sticks_of(f, XA.stick_params(_ns(two_theta=[10, 90])))
            check(dev < 1e-6 and len(s["x"]) == npk,
                  f"{ph} {'x'.join(map(str, rep))} ({len(a) * np.prod(rep)} atoms): pattern = 12-atom cell (max dev {dev:.1e}), {len(s['x'])} sticks")
    lam = {"CuKa": 1.54184, "MoKa": 0.71073}
    a0 = np.linalg.norm(read(str(ROOT / "refs" / "POSCAR_c"), format="vasp").cell[0])
    for w, L in lam.items():
        s = XA.sticks_of(ROOT / "refs" / "POSCAR_c", XA.stick_params(_ns(wavelength=w, two_theta=[5, 60])))
        tt = 2 * math.degrees(math.asin(L / (2 * a0 / math.sqrt(3))))
        i = int(np.argmax(s["y"]))
        check(abs(s["x"][i] - tt) < 2e-3 and s["hkl"][i].startswith("(1 1 1)"),
              f"CeO2 strongest peak (111) at {s['x'][i]:.3f} deg with {w} (Bragg: {tt:.3f})")


def _ns(**kw):
    ap = XA.build_parser()
    a = ap.parse_args([])
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def species(mock):
    from ase.io import read, write
    print("== 2. species replacement happens in memory only")
    work = mock / "campaign"
    rel = sorted((work / "relaxed").rglob("POSCAR_*"))
    files = [p for p in rel if "/mapped/m/" in str(p)][:3] + [p for p in rel if "/mapped/" not in str(p)][:2]
    before = {p: (md5(p), p.stat().st_mtime_ns) for p in files}
    outs = {}
    for tag, extra in [("cat", ["--replace-cations", "Zr"]), ("map", ["--replace", "Ce:Zr"])]:
        r = run(["--structures", *files, *extra, "--ref", "c", "--out-dir", W / "sp", "--stem", tag, "--no-plot", "--jobs", "2"])
        outs[tag] = pattern(W / "sp", tag)[0]
        st = read_csv(W / "sp" / f"{tag}_structures.csv")
        check(r.returncode == 0 and all(x["formula_used"].startswith("Zr") and "Ce" not in x["formula_used"] for x in st),
              f"{tag}: every structure computed as {st[0]['formula_used']}")
        rp = read_csv(W / "sp" / f"{tag}_ref_peaks.csv")
        check(all(x["formula_used"] == "Zr4O8" for x in rp), f"{tag}: reference c computed as Zr4O8 (species change applied to refs)")
    pure = W / "sp" / "pure"
    pure.mkdir(exist_ok=True)
    copies = []
    for p in files:
        a = read(str(p), format="vasp")
        a.set_chemical_symbols(["Zr" if s == "Ce" else s for s in a.get_chemical_symbols()])
        q = pure / (p.parent.name + "_" + p.name)
        write(str(q), a, format="vasp", direct=True, sort=False)
        copies.append(q)
    run(["--structures", *copies, "--out-dir", W / "sp", "--stem", "pure", "--no-plot", "--jobs", "2"])
    yp = pattern(W / "sp", "pure")[0]
    d1 = np.max(np.abs(outs["cat"][:, 1] - yp[:, 1]))
    d2 = np.max(np.abs(outs["map"][:, 1] - yp[:, 1]))
    check(d1 < 1e-6 and d2 < 1e-6, f"--replace-cations Zr = --replace Ce:Zr = files written with Zr (max dev {max(d1, d2):.1e})")
    r = run(["--structures", *files, "--replace-cations", "Zr", "--ref", "c", "--no-ref-replace", "--out-dir", W / "sp",
             "--stem", "norefrep", "--no-plot"])
    check(all(x["formula_used"] == "Ce4O8" for x in read_csv(W / "sp" / "norefrep_ref_peaks.csv")), "--no-ref-replace keeps CeO2 for the reference")
    run(["--structures", *files, "--out-dir", W / "sp", "--stem", "asis", "--no-plot"])
    ya = pattern(W / "sp", "asis")[0]
    check(np.max(np.abs(ya[:, 1] - yp[:, 1])) > 1, "without replacement the pattern differs (Ce/Zr contrast)")
    after = {p: (md5(p), p.stat().st_mtime_ns) for p in files}
    check(before == after, f"the {len(files)} relaxed POSCARs are unchanged (md5 and mtime)")


def averaging():
    print("== 3. averaging and broadening")
    f = ROOT / "refs" / "POSCAR_oIII"
    run(["--structures", f, f, f, f, "--out-dir", W / "avg", "--stem", "four", "--no-plot", "--jobs", "1"])
    run(["--structures", f, "--out-dir", W / "avg", "--stem", "one", "--no-plot", "--jobs", "1"])
    y4, y1 = pattern(W / "avg", "four")[0], pattern(W / "avg", "one")[0]
    check(np.max(np.abs(y4[:, 1] - y1[:, 1])) < 1e-9 and np.max(y4[:, 2]) < 1e-6, "4 copies: mean = single pattern, std = 0")
    check(abs(np.max(y1[:, 1]) - 100) < 1e-6, "--norm max: the single pattern peaks at 100")
    for prof in ("gaussian", "lorentzian", "pseudo-voigt"):
        a = _ns(profile=prof, fwhm=0.2, step=0.001, two_theta=[20.0, 40.0], eta=0.4)
        g = XA.make_grid(a)
        y = XA.broaden(np.array([30.0005]), np.array([7.0]), g, a)
        half = g[y >= 3.5]
        check(abs(y.max() - 7.0) < 1e-3 and abs((half[-1] - half[0]) - 0.2) < 3e-3,
              f"{prof}: height = stick intensity ({y.max():.4f}), FWHM {half[-1] - half[0]:.3f}")
    a = _ns(caglioti=[0.01, -0.005, 0.004])
    w = XA.fwhm_at(np.array([20.0, 80.0]), a)
    check(w[1] > w[0] > 0, f"Caglioti FWHM grows with 2theta ({w[0]:.3f} -> {w[1]:.3f})")
    r = run(["--structures", f, "--profile", "pseudo-voigt", "--eta", "0.3", "--caglioti", "0.01", "-0.005", "0.004",
             "--wavelength", "1.0", "--two-theta", "5", "60", "--step", "0.01", "--debye-waller", "Zr=0.3", "O=0.5",
             "--ref", "c", "t", "oIII", "m", "--out-dir", W / "avg", "--stem", "opts", "--formats", "png", "pdf"])
    check(r.returncode == 0 and (W / "avg" / "opts.png").is_file() and (W / "avg" / "opts.pdf").is_file(),
          "pseudo-Voigt + Caglioti + 1.0 A + Debye-Waller + 4 references, PNG + PDF")
    rec = json.loads((W / "avg" / "opts.json").read_text())
    check(rec["stick_params"]["debye_waller"] == {"O": 0.5, "Zr": 0.3} and rec["script_md5"] == md5(SCRIPT), "JSON record: settings + script md5")
    r = run(["--structures", f, "--replace", "Ce-Zr"], may_fail=True)
    check(r.returncode != 0 and "OLD:NEW" in r.stderr, "bad --replace is refused with a hint")


def manual_selection(work, comp, phases, lo, hi, include_flagged):
    best = {}
    for p in sorted((work / "results").glob("chunk_*.csv")):
        rd = int(p.name[11:14]) if p.name.startswith("chunk_rerun") else None
        for r in read_csv(p):
            f = r["file"]
            if f not in best or r["status"] == "ok" or best[f][0]["status"] != "ok":
                best[f] = (r, rd)
    rows = {f: v for f, v in best.items() if f.startswith(comp + "/")}
    clean = [float(r["F_total_eV_atom"]) for r, _ in rows.values() if r["n_imaginary_modes"] in ("0", "") and r["relax_converged"] != "False"]
    Fr = min(clean)
    out = {}
    for f, (r, rd) in rows.items():
        F = float(r["F_total_eV_atom"])
        flag = r["n_imaginary_modes"] not in ("0", "") or r["relax_converged"] == "False"
        if r["phase"] in phases and Fr + lo / 1000 <= F <= Fr + hi / 1000 and (include_flagged or not flag):
            out[f] = rd
    return out, Fr


def listing(args):
    r = run([*args, "--list"])
    rows = [l.split("\t") for l in r.stdout.strip().split("\n") if "\t" in l]
    return {x[0]: (x[3], x[4]) for x in rows}, r


def selection(mock, tag):
    print(f"== 4. selection on the mock campaign {tag}")
    work = mock / "campaign"
    comp = next(p.name for p in work.iterdir() if p.is_dir() and p.name.startswith("Ce_"))
    for phases, lo, hi, fl in [(["oIII"], 1.5, 4.0, False), (["m", "t"], 0, 3, True), (["parent", "t", "oIII", "m"], -1, 99, False)]:
        exp, Fr = manual_selection(work, comp, phases, lo, hi, fl)
        got, r = listing([mock, "--composition", comp, "--phase", *phases, "--F-range", lo, hi] + (["--include-flagged"] if fl else []))
        check(set(got) == set(exp) and len(got) > 0, f"{'/'.join(phases)} {lo}..{hi} meV{' +flagged' if fl else ''}: {len(got)} structures, same as from the CSVs")
        wrong = [k for k, (how, p) in got.items() if (exp[k] is not None and f"rerun_{exp[k]:03d}/relaxed/{k}" not in p) or
                 (exp[k] is None and not p.endswith(f"campaign/relaxed/{k}"))]
        check(not wrong, f"geometry: the relaxed file of the row that counts ({sum(1 for k in got if exp[k] is not None)} from rerun_001)")
    # the app's endpoint gives the same rows
    import importlib
    sys.path.insert(0, str(ROOT))
    X = importlib.import_module("studio.xrd")
    d = X.rows(str(mock), None, comp)
    app_sel = [x["k"] for x in d["rows"] if x["p"] == "oIII" and x["st"] == "ok" and not (x["im"] or x["un"]) and
               d["F_ref"] + 1.5e-3 - 1e-12 <= x["F"] <= d["F_ref"] + 4e-3 + 1e-12]
    exp, _ = manual_selection(work, comp, ["oIII"], 1.5, 4.0, False)
    check(set(app_sel) == set(exp), f"app rows endpoint: same {len(app_sel)} oIII structures")
    # per-composition CSV (file relative to Ce_xx) inside the workdir -> same geometry
    best = {}
    for p in sorted((work / "results").glob("chunk_*.csv")):
        for r in read_csv(p):
            f = r["file"]
            if f not in best or r["status"] == "ok" or best[f]["status"] != "ok":
                best[f] = r
    pc = work / comp / "free_energy_300K.csv"
    with open(pc, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(next(iter(best.values())).keys()))
        wr.writeheader()
        for r in best.values():
            wr.writerow({**r, "file": r["file"].split("/", 1)[1]})
    got2, _ = listing([pc, "--phase", "oIII", "--F-range", 1.5, 4.0])
    got1, _ = listing([mock, "--composition", comp, "--phase", "oIII", "--F-range", 1.5, 4.0])
    check(got2 == got1, f"per-composition CSV: same {len(got2)} structures and the same relaxed files")
    # plain folder (Ce_xx/free_energy_300K.csv + Ce_xx/relaxed/...)
    plain = W / f"plain_{tag}"
    shutil.rmtree(plain, ignore_errors=True)
    (plain / comp).mkdir(parents=True)
    shutil.copy(pc, plain / comp / pc.name)
    shutil.copytree(work / "relaxed" / comp, plain / comp / "relaxed")
    got3, _ = listing([plain, "--phase", "oIII", "--F-range", 1.5, 4.0])
    check(set(got3) == set(got1) and all(f"/{comp}/relaxed/" in p for _, p in got3.values()),
          f"plain folder: same {len(got3)} structures, geometry from {comp}/relaxed/")
    # missing relaxed file
    k = sorted(got3)[0]
    victim = plain / comp / "relaxed" / k.split("/", 1)[1]
    victim.unlink()
    got4, r = listing([plain, "--phase", "oIII", "--F-range", 1.5, 4.0])
    check(k not in got4 and "1 without relaxed file" in r.stdout, "a missing relaxed file is skipped and reported")
    got5, _ = listing([plain, "--phase", "oIII", "--F-range", 1.5, 4.0, "--allow-unrelaxed"])
    check(k not in got5 or "INPUT" in got5[k][0], "--allow-unrelaxed: input only if it exists (not in a plain copy: still skipped)")
    # --max subset reproducible
    a1, _ = listing([mock, "--composition", comp, "--F-range", -1, 99, "--max", 7, "--seed", 3])
    a2, _ = listing([mock, "--composition", comp, "--F-range", -1, 99, "--max", 7, "--seed", 3])
    a3, _ = listing([mock, "--composition", comp, "--F-range", -1, 99, "--max", 7, "--seed", 4])
    check(len(a1) == 7 and a1 == a2 and a1 != a3, "--max 7 --seed: reproducible random subset")
    # a selection computed end to end
    r = run([mock, "--composition", comp, "--phase", "oIII", "--F-range", 1.5, 4.0, "--ref", "c", "t", "oIII", "m",
             "--out-dir", W / f"out_{tag}", "--cache-dir", W / "cache", "--formats", "png"])
    check(r.returncode == 0 and len(read_csv(W / f"out_{tag}" / f"xrd_{comp}_oIII_1.5-4meV_structures.csv")) == len(got1),
          f"computed: {len(got1)} structures averaged, PNG written")


def app(mock):
    print("== 5. app endpoints")
    st = W / "state"
    shutil.rmtree(st, ignore_errors=True)
    port = 8700 + os.getpid() % 200
    srv = subprocess.Popen([PY, str(ROOT / "start_studio.py"), "--port", str(port), "--state", str(st), "--dir", str(W)],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        for _ in range(50):
            if (st / "token").is_file():
                break
            time.sleep(0.2)
        time.sleep(0.8)
        line = ""
        while "open" not in line:
            line = srv.stdout.readline()
        base = line.split("open    : ")[1].split("/?token=")[0]
        tok = (st / "token").read_text().strip()

        def call(path, body=None):
            req = urllib.request.Request(base + path, data=None if body is None else json.dumps(body).encode(),
                                         headers={"X-Token": tok, "X-Requested-With": "studio", "Content-Type": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=60) as f:
                    return json.loads(f.read())
            except urllib.error.HTTPError as e:
                return json.loads(e.read())

        d = call(f"/api/xrd/rows?source={mock}")
        comp = d["comp"]
        check(d["ok"] and len(d["rows"]) > 0 and d["F_ref"] is not None, f"rows: {len(d['rows'])} structures of {comp}")
        b = {"source": str(mock), "comp": comp, "phases": ["oIII"], "lo": 1.5, "hi": 4.0, "max": 500, "seed": 0, "jobs": 2,
             "params": {"wavelength": "CuKa", "tt_lo": 10, "tt_hi": 90, "refs": ["c", "t", "oIII", "m"],
                        "species": {"mode": "cations", "to": "Zr", "anions": "O"}},
             "broadening": {"fwhm": 0.1, "step": 0.02, "profile": "gaussian", "norm": "max"}}
        bad = call("/api/xrd/compute", {**b, "params": {**b["params"], "species": {"mode": "cations", "to": "Zr;rm"}}})
        check(not bad["ok"] and "element" in bad["error"], "bad element refused")
        times = []
        for k in range(2):
            t0 = time.time()
            r = call("/api/xrd/compute", b)
            for _ in range(600):
                j = call(f"/api/jobs?id={r['job']}")["job"]
                if j["state"] not in ("running", "starting"):
                    break
                time.sleep(0.3)
            times.append(time.time() - t0)
            check(j["state"] == "done", f"compute job {k + 1}: {j['state']} in {times[-1]:.1f} s")
        req = urllib.request.Request(base + f"/api/xrd/sticks?key={r['key']}", headers={"X-Token": tok, "Accept-Encoding": "gzip"})
        import gzip
        with urllib.request.urlopen(req) as f:
            raw = f.read()
            enc = f.headers.get("Content-Encoding")
        s = json.loads(gzip.decompress(raw) if enc == "gzip" else raw)
        exp, _ = manual_selection(mock / "campaign", comp, ["oIII"], 1.5, 4.0, False)
        check(len(s["structures"]) == len(exp) and all(x["formula"].startswith("Zr") for x in s["structures"]),
              f"sticks JSON ({enc}, {len(raw) // 1024} kB): {len(exp)} structures, computed as {s['structures'][0]['formula']}")
        check(len(s["refs"]) == 4 and all(x["formula"] == "Zr4O8" for x in s["refs"]), "4 references, species change applied")
        check("from the cache" in j["log"] and f"{len(exp)} from the cache" in j["log"], "second run: every structure from the cache")
        od = W / "export"
        r = call("/api/xrd/export", {**b, "out_dir": str(od)})
        for _ in range(600):
            j = call(f"/api/jobs?id={r['job']}")["job"]
            if j["state"] not in ("running", "starting"):
                break
            time.sleep(0.3)
        files = sorted(p.name for p in od.glob("*"))
        check(j["state"] == "done" and any(f.endswith(".png") for f in files) and any(f.endswith(".pdf") for f in files),
              f"export job writes {', '.join(files)}")
        check(call("/api/xrd/sticks?key=../../etc")["ok"] is False, "sticks key is checked")
    finally:
        srv.terminate()


def main():
    W.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    physics()
    averaging()
    mocks = {}
    for tag, rep in [("2x2x2", (2, 2, 2)), ("1x2x4", (1, 2, 4))]:
        m = W / f"mock_{tag}"
        if not (m / "campaign" / "results" / "chunk_00000.csv").is_file():
            shutil.rmtree(m, ignore_errors=True)
            subprocess.run([PY, str(ROOT / "tests" / "make_xrd_mock.py"), str(m), "--orderings", "30", "--rep", *map(str, rep)],
                           check=True)
        mocks[tag] = m
    species(mocks["2x2x2"])
    for tag, m in mocks.items():
        selection(m, tag)
    app(mocks["2x2x2"])
    print(f"\n{'ALL OK' if not FAILS else str(len(FAILS)) + ' FAILED'} ({time.time() - t0:.0f} s)")
    for f in FAILS:
        print("  FAIL:", f)
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
