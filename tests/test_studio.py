#!/usr/bin/env python3
"""
test_studio.py - end-to-end test of the web app through its HTTP API (no browser, no NEP, no cluster), ~3-6 min.

    python3 tests/test_studio.py            # from the package root; keeps its files in tests/work_studio/

What it does, all through the same endpoints the page uses:
  * security: no token -> 403, POST without the page header -> 403
  * every registered module: the option form can be read from its script (argparse)
  * New run from a cubic host (1x1x2 and 2x2x2) and from a TETRAGONAL host (starting phase t), real 12-atom
    references refs/POSCAR_{t,oIII,m}: run_record.json + campaign.toml written, inputs and scripts frozen
  * stages run locally as jobs: dry-run, prepare, plan, compute (mock calorine = Einstein crystal), status, merge
  * the folder scan agrees with the files (orderings, mapped, IDENTICAL, relaxed, F done, flags)
  * every plot module runs on the merged / per-composition CSV and writes its files
  * map_check (map_fluorite and transform_set) on c, t, oIII and m supercells of shapes 2x2x2, 1x2x4, 3x2x2:
    every output IDENTICAL to the input ordering
  * SLURM: scripts written and "submitted" with a fake sbatch, job ids picked up
  * reproducibility check: identical right after creation, a changed input is reported
"""
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
W = ROOT / "tests" / "work_studio"
PY = sys.executable
FCC = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]
TET = [[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)]
FAILS = []


def check(cond, msg):
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def write_cubic(fn, n, a=5.10):
    cells = [[i, j, k] for i in range(n[0]) for j in range(n[1]) for k in range(n[2])]
    cat = [[(c[d] + f[d]) / n[d] for d in range(3)] for c in cells for f in FCC]
    ox = [[(c[d] + f[d]) / n[d] for d in range(3)] for c in cells for f in TET]
    with open(fn, "w") as f:
        f.write(f"ZrO2 cubic {n[0]}x{n[1]}x{n[2]}\n1.0\n")
        for d in range(3):
            v = [0.0, 0.0, 0.0]
            v[d] = a * n[d]
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write(f"  Zr O\n  {len(cat)} {len(ox)}\nDirect\n")
        for p in cat + ox:
            f.write("  %.10f %.10f %.10f\n" % tuple(p))


class Api:
    def __init__(self, port, token):
        self.base, self.token = f"http://127.0.0.1:{port}", token

    def call(self, path, body=None, token=True, header=True):
        req = urllib.request.Request(self.base + path, data=None if body is None else json.dumps(body).encode())
        if token:
            req.add_header("X-Token", self.token)
        if header:
            req.add_header("X-Requested-With", "studio")
        if body is not None:
            req.add_header("Content-Type", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.status, json.loads(r.read() or b"{}")
        except urllib.error.HTTPError as e:
            return e.code, json.loads(e.read() or b"{}")

    def get(self, path):
        return self.call(path)[1]

    def post(self, path, body):
        return self.call(path, body)[1]

    def wait(self, jid, limit=900):
        t0 = time.time()
        while time.time() - t0 < limit:
            j = self.get(f"/api/jobs?id={jid}&tail=4000")["job"]
            if j["state"] in ("done", "failed", "lost", "cancelled"):
                return j
            time.sleep(1)
        raise TimeoutError(jid)


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def stage(api, cfg, module, mode="local", options=None, expect_ok=True):
    r = api.post("/api/stage", {"module": module, "mode": mode, "config": cfg, "options": options or {}})
    assert r.get("ok"), r
    j = api.wait(r["job"])
    if expect_ok:
        check(j["state"] == "done", f"stage {module} [{mode}] exit {j['exit']}")
        if j["state"] != "done":
            print(j.get("log", "")[-3000:])
    return j


def main():
    if W.exists():
        shutil.rmtree(W)
    (W / "bin").mkdir(parents=True)
    os.chdir(W)
    # fake sbatch (prints a job id, runs nothing) for the SLURM part
    (W / "bin" / "sbatch").write_text("#!/bin/sh\necho $$\n")
    os.chmod(W / "bin" / "sbatch", 0o755)
    env = dict(os.environ, PYTHONPATH=str(ROOT / "tests" / "mock_calorine"), OMP_NUM_THREADS="1",
               PATH=f"{W / 'bin'}:{os.environ['PATH']}")
    port = free_port()
    srv = subprocess.Popen([PY, str(ROOT / "start_studio.py"), "--port", str(port), "--state", str(W / "state"),
                            "--dir", str(W)], env=env, stdout=open(W / "server.log", "w"), stderr=subprocess.STDOUT)
    try:
        for _ in range(50):
            if (W / "state" / "token").is_file():
                break
            time.sleep(0.2)
        time.sleep(1)
        api = Api(port, (W / "state" / "token").read_text().strip())
        run_tests(api, env)
    finally:
        srv.terminate()
    print("\nSTUDIO TEST " + ("PASSED" if not FAILS else f"FAILED ({len(FAILS)}): " + "; ".join(FAILS)))
    sys.exit(1 if FAILS else 0)


def run_tests(api, env):
    print("== security")
    check(api.call("/api/info", token=False)[0] == 403, "request without token refused")
    check(api.call("/api/mkdir", {"path": str(W / "x")}, header=False)[0] == 403, "POST without page header refused")
    info = api.get("/api/info")
    check(info["ok"] and len(info["pages"]) >= 5, f"info: pages {info['pages']}")

    print("== every module: option form from its script")
    mods = api.get("/api/modules")["modules"]
    check(not api.get("/api/modules")["errors"], "module files parse")
    for m in mods:
        if m.get("script"):
            sp = api.get(f"/api/argspec?module={m['id']}")
            check(sp.get("ok") and sp.get("actions"), f"argspec {m['id']} ({len(sp.get('actions', []))} options)")

    print("== inputs: hosts and supercells")
    refs = ROOT / "refs"
    write_cubic(W / "POSCAR_c_1x1x2", (1, 1, 2))
    write_cubic(W / "POSCAR_c_2x2x2", (2, 2, 2))
    (W / "nep.txt").write_text("mock nep\n")
    # tetragonal host = the cubic 2x2x2 mapped onto t (map_fluorite.py, not by hand)
    subprocess.run([PY, str(ROOT / "scripts" / "map_fluorite.py"), "POSCAR_c_2x2x2", "--template", str(refs / "POSCAR_t"),
                    "--out", "POSCAR_t_2x2x2", "--quiet"], check=True)
    pk = api.get("/api/poscar?path=" + str(W / "POSCAR_c_2x2x2"))
    check(pk["ok"] and pk["species"] == {"Zr": 32, "O": 64} and pk["repeat_guess"] == [2, 2, 2], f"POSCAR peek {pk.get('species')} {pk.get('repeat_guess')}")

    runs = {}
    for name, host, start, comps, maxo in [("run_c112", "POSCAR_c_1x1x2", "c", "all", 6),
                                           ("run_c222", "POSCAR_c_2x2x2", "c", "0, 1, 16, 31, 32", 3),
                                           ("run_t222", "POSCAR_t_2x2x2", "t", "0, 8, 16", 3)]:
        print(f"== new run {name} (starting phase {start}, {host})")
        tpl = [f"{p}={refs / ('POSCAR_' + p)}" for p in ("t", "oIII", "m") if p != start] + ([f"c={refs / 'POSCAR_c'}"] if start != "c" else [])
        r = api.post("/api/runs/create", {
            "folder": str(W / name), "name": name, "starting_phase": start, "notes": "test",
            "copy_inputs": True, "freeze_scripts": True,
            "config": {"structure": {"poscar": str(W / host), "host": "Zr", "dopant": "Ce", "compositions": comps,
                                     "max_orderings": maxo, "seed": 2},
                       "mapping": {"templates": tpl, "jobs": 4},
                       "free_energy": {"nep": str(W / "nep.txt"), "T": 300, "supercell": "1 1 1", "mesh": 4, "fmax": 1e-3,
                                       "steps": 200},
                       "cluster": {"cores_per_node": 4, "walltime": "00:10:00", "sec_per_structure": 2, "max_concurrent": 8,
                                   "setup": "export RAYON_NUM_THREADS=1"}}})
        check(r.get("ok"), f"run created {r.get('error', '')}")
        cfg = r["config"]
        runs[name] = cfg
        rec = json.loads((W / name / "run_record.json").read_text())
        check((W / name / "scripts" / "campaign.py").is_file() and len(rec["scripts"]) >= 18, "scripts frozen + md5 recorded")
        check(all((W / name / k).is_file() for k in rec["inputs"]), f"{len(rec['inputs'])} inputs copied with md5")
        for st in ("dry_run", "prepare", "plan", "compute", "status", "merge"):
            stage(api, cfg, st)
        sc = api.get("/api/scan?path=" + str(W / name))
        T = "300"
        tot = sc["totals"][T]
        n_ord = sum(c["orderings"] for c in sc["compositions"])
        n_phase = 1 + len(tpl)
        check(sc["layout"] == "campaign" and tot["structures"] == n_ord * n_phase,
              f"scan: {len(sc['compositions'])} compositions, {n_ord} orderings x {n_phase} phases = {tot['structures']} structures")
        check(tot["ok"] == tot["structures"] and tot["not_run"] == 0 and tot["failed"] == 0,
              f"scan: F done {tot['ok']}, failed {tot['failed']}, not run {tot['not_run']}")
        check(tot["relaxed"] == tot["structures"], f"scan: relaxed files {tot['relaxed']}")
        ident = sum(v["identical"] for c in sc["compositions"] for v in c["verify"].values())
        other = sum(v["other"] for c in sc["compositions"] for v in c["verify"].values())
        check(ident == n_ord * len(tpl) and other == 0, f"scan: {ident} mapped structures IDENTICAL, {other} not")
        mapped = {t: c["mapping"][t] for c in sc["compositions"] for t in c["mapping"]}
        print("        mapping (last composition): " + "; ".join(
            f"{t}: rms {m['rms_shift_A']['max']:.3f} A, S {m['strain_S']['max']:.4f}, variants {m['tied_variants']['max']:.0f}"
            for t, m in mapped.items()))
        check(sc["merged_csv"] == ["free_energy_300K_all.csv"] and sc["chunks"]["complete"] == sc["chunks"]["n"],
              f"merged CSV + {sc['chunks']['n']} chunks complete")
        # same scan when pointed at the workdir itself
        sc2 = api.get("/api/scan?path=" + sc["workdir"])
        check(sc2["totals"][T] == tot, "scan of the workdir = scan of the run folder")
        runs[name + "_scan"] = sc

    print("== flags: a row with imaginary modes / unconverged relaxation is counted")
    sc = runs["run_c112_scan"]
    wd = Path(sc["workdir"])
    res = sorted((wd / "results").glob("chunk_*.csv"))[0]
    lines = res.read_text().split("\n")
    head = lines[0].split(",")
    row = lines[1].split(",")
    row[head.index("n_imaginary_modes")] = "3"
    row[head.index("relax_converged")] = "False"
    lines[1] = ",".join(row)
    res.write_text("\n".join(lines))
    t2 = api.get("/api/scan?fresh=1&path=" + str(wd))["totals"]["300"]
    check(t2["imaginary"] == 1 and t2["unconverged"] == 1 and t2["clean"] == t2["ok"] - 1, f"flags counted {t2['imaginary']} / {t2['unconverged']}")

    print("== plots on the merged and per-composition CSV (cubic 1x1x2 run, all compositions)")
    stage(api, runs["run_c112"], "merge")
    merged = str(wd / "free_energy_300K_all.csv")
    comp_csv = str(wd / "Ce_50" / "free_energy_300K.csv")
    for mid, inp, args in [("f_vs_x_html", merged, ["--bootstrap", "5", "--jobs", "2"]),
                           ("ca_vs_x_html", merged, ["--bootstrap", "5", "--jobs", "2", "--quantities", "ca", "theta"]),
                           ("f_vs_composition", merged, ["--bootstrap", "5", "--jobs", "2", "--formats", "png"]),
                           ("f_vs_ce_hist", merged, ["--ref", "composition", "--bootstrap", "5", "--formats", "png"]),
                           ("free_energy_pf", comp_csv, ["--bootstrap", "5"]),
                           ("cell_annotations", comp_csv, ["--bootstrap", "5"]),
                           ("cell_distributions", comp_csv, [])]:
        r = api.post("/api/run", {"module": mid, "input": inp, "args": args})
        j = api.wait(r["job"])
        outs = [o["name"] for o in j.get("outputs", [])]
        good = [o for o in outs if o.endswith((".html", ".png", ".pdf", ".svg"))]
        check(j["state"] == "done" and good, f"plot {mid}: {', '.join(good[:4])}")
        if j["state"] != "done":
            print(j.get("log", "")[-2500:])

    print("== map_check: c, t, oIII, m supercells of several shapes, map_fluorite and transform_set")
    shapes = [(2, 2, 2), (1, 2, 4), (3, 2, 2)]
    for n in shapes:
        tag = "x".join(map(str, n))
        write_cubic(W / f"P_c_{tag}", n)
        for ph in ("t", "oIII", "m"):
            subprocess.run([PY, str(ROOT / "scripts" / "map_fluorite.py"), f"P_c_{tag}", "--template", str(refs / f"POSCAR_{ph}"),
                            "--out", f"P_{ph}_{tag}", "--quiet"], check=True)
        # decorate: put Ce on a few cation sites of every input (same sites), so the ordering is not trivial
        for ph in ("c", "t", "oIII", "m"):
            p = W / f"P_{ph}_{tag}"
            L = p.read_text().split("\n")
            ncat = int(L[6].split()[0])
            sites = [0, 3, 5] if ncat > 6 else [0]
            coords = L[8:8 + int(L[6].split()[0]) + int(L[6].split()[1])]
            cat = [coords[i] for i in range(ncat)]
            ce = [cat[i] for i in sites]
            zr = [cat[i] for i in range(ncat) if i not in sites]
            L[5] = "  Zr Ce O"
            L[6] = f"  {len(zr)} {len(ce)} {L[6].split()[1]}"
            L[8:8 + ncat] = zr + ce
            p.write_text("\n".join(L))
        for ph in ("c", "t", "oIII", "m"):
            inp = str(W / f"P_{ph}_{tag}")
            for method, args in [("map", ["--method", "map", "--templates"] + [f"{q}={refs / ('POSCAR_' + q)}" for q in ("c", "t", "oIII", "m") if q != ph]),
                                 ("transform", ["--method", "transform", "--refs", str(refs), "--targets", "t", "oIII", "m", "c"])]:
                out = W / "mapcheck" / f"{ph}_{tag}_{method}"
                r = api.post("/api/run", {"module": "map_check", "input": inp, "args": args + ["--out-dir", str(out)]})
                j = api.wait(r["job"])
                rows = []
                if (out / "map_check.csv").is_file():
                    import csv
                    rows = list(csv.DictReader(open(out / "map_check.csv")))
                ok = j["state"] == "done" and rows and all(x["verdict"] == "IDENTICAL" for x in rows)
                summ = ", ".join(f"{x['target']} rms {x['rms_shift_A']} S {x['strain_S']}" + (f" var {x['equally_good_variants']}" if x['equally_good_variants'] else "") for x in rows)
                check(ok, f"{ph} {tag} {method}: {len(rows)} outputs IDENTICAL [{summ}]")
                if not ok:
                    print(j.get("log", "")[-2500:])

    print("== SLURM: write and submit with a fake sbatch")
    cfg = runs["run_c222"]
    j = stage(api, cfg, "compute", "slurm_write", {"--unfinished": False})
    check((Path(cfg).parent / "campaign" / "slurm" / "compute.sh").is_file(), "compute.sh written")
    pv = api.post("/api/stage", {"module": "compute", "mode": "slurm", "config": cfg, "options": {}, "preview": True})
    check(pv["ok"] and "--unfinished" in pv["cmd"] and pv["confirm"], f"SLURM compute: --unfinished by default, asks first: {pv['cmd'][-40:]}")
    j = stage(api, cfg, "compute", "slurm", {"--unfinished": False, "--with-merge": True})
    check(len(j["slurm_ids"]) >= 2, f"submitted job ids picked up: {j['slurm_ids']}")
    sh = (Path(cfg).parent / "campaign" / "slurm" / "compute.sh").read_text()
    check("RAYON_NUM_THREADS=1" in sh and "/scripts/campaign.py" in sh, "SLURM script uses the setup lines and the run's frozen scripts")
    code, rr = api.call("/api/stage", {"module": "plan", "mode": "local", "config": cfg, "options": {"plan_mode": "--bogus"}, "preview": True})
    check(code == 400, "plan: a choice outside the list is refused")
    code, rr = api.call("/api/stage", {"module": "rerun_flagged", "mode": "dry", "config": cfg, "options": {"--what": "everything"}, "preview": True})
    check(code == 400, "rerun: --what outside both|imaginary|unconverged refused")
    pv = api.post("/api/stage", {"module": "rerun_flagged", "mode": "slurm", "config": cfg, "preview": True,
                                 "options": {"--what": "imaginary", "--freq-window": "-0.5 0", "--fmax": "0.001", "--steps": "4000",
                                             "--start": "original", "--nodes": "4", "--walltime": "02:00:00"}})
    check(pv["ok"] and all(x in pv["cmd"] for x in ("--what imaginary", "--freq-window -0.5 0", "--fmax 0.001", "--steps 4000",
                                                    "--start original", "--nodes 4", "--walltime 02:00:00")) and pv["confirm"],
          "rerun options reach rerun_flagged.py: " + pv["cmd"].split("rerun_flagged.py", 1)[-1])
    pv = api.post("/api/stage", {"module": "rerun_flagged", "mode": "local", "config": cfg, "preview": True, "options": {"--nodes": "4"}})
    check("--nodes" not in pv["cmd"] and "--what" not in pv["cmd"], "SLURM-only rerun options left out when running here; defaults not repeated")
    j = stage(api, cfg, "rerun_flagged", "dry")
    check("rerun" in j.get("log", "").lower() or "flagged" in j.get("log", "").lower(), "rerun_flagged dry-run ran")

    print("== run dashboard state (/api/run/state)")
    st = api.get("/api/run/state?config=" + cfg)
    check(st.get("ok") and [x["key"] for x in st["stages"]] == ["check", "prepare", "plan", "compute", "status", "merge", "rerun"],
          "7 stages in order")
    sc = runs["run_c222_scan"]
    check(st["counts"]["structures"] == sc["totals"]["300"]["structures"] and st["counts"]["clean"] == sc["totals"]["300"]["clean"],
          f"state counts = scan counts ({st['counts']['structures']} structures, {st['counts']['clean']} clean)")
    subs = st["slurm"]["submissions"]
    check({x["kind"] for x in subs} >= {"compute", "merge"} and any(x["kind"] == "compute" and x["n_tasks"] >= 1 for x in subs),
          f"submitted SLURM jobs read from the job log: {[(x['id'], x['kind'], x['n_tasks']) for x in subs]}")
    js = subprocess.run([PY, str(Path(cfg).parent / "scripts" / "campaign.py"), cfg, "status", "--json"], capture_output=True, text=True, env=env)
    sj = json.loads(js.stdout)
    m1, m2 = sj["sec_per_structure"]["median"], st["timing"]["measured"]["median"]
    check(abs(m1 - m2) <= 0.006, f"measured sec_per_structure: app {m2:.4f} = campaign.py status --json {m1:.2f}")
    check(sj["planned"] == st["counts"]["planned"] and sj["chunks"]["unfinished"] == st["chunks"]["unfinished"], "status --json agrees on chunks")
    # live SLURM view with a fake squeue / sacct: the compute task runs, the merge job waits
    cid = [x["id"] for x in subs if x["kind"] == "compute"][-1]
    mid = [x["id"] for x in subs if x["kind"] == "merge"][-1]
    (W / "bin" / "squeue").write_text(f"#!/bin/sh\necho '{cid}_0|RUNNING|0:04:10|0:10:00|1|normal|node7|fo_free_energy'\n"
                                      f"echo '{mid}|PENDING|0:00|1:00:00|1|normal|(Dependency)|fo_merge'\n")
    (W / "bin" / "sacct").write_text("#!/bin/sh\nexit 0\n")
    for f in ("squeue", "sacct"):
        os.chmod(W / "bin" / f, 0o755)
    st = api.get("/api/run/state?fresh=1&config=" + cfg)
    t0 = [t for t in st["slurm"]["tasks"] if t["job"] == cid and t["task"] == 0]
    check(st["slurm"]["available"] and t0 and t0[0]["state"] == "RUNNING" and t0[0]["chunk"] == 0 and t0[0]["elapsed_s"] == 250,
          f"squeue task -> chunk 0 RUNNING 250 s: {t0[:1]}")
    check(st["compute_running"] and [x for x in st["stages"] if x["key"] == "compute"][0]["state"] == "running",
          "compute shown as running while its SLURM task runs")
    check(st["next"]["action"].get("go") == "jobs", f"next step while running: {st['next']['what']}")
    pv = api.post("/api/stage", {"module": "plan", "mode": "local", "config": cfg, "options": {"plan_mode": "--force"}, "preview": True})
    check(pv["ok"] and pv["blocked"] and "--force" in pv["cmd"] and pv["confirm"], f"re-pack refused while compute runs: {pv.get('blocked')}")
    code, rr = api.call("/api/stage", {"module": "plan", "mode": "local", "config": cfg, "options": {"plan_mode": "--force"}})
    check(code == 400 and "compute" in rr.get("error", ""), f"POST plan --force refused while compute runs ({code})")
    pv = api.post("/api/stage", {"module": "plan", "mode": "local", "config": cfg, "options": {"plan_mode": "--append"}, "preview": True})
    check(pv["ok"] and not pv["blocked"] and pv["cmd"].endswith("plan --append"), "plan --append allowed while compute runs")
    for f in ("squeue", "sacct"):
        (W / "bin" / f).unlink()
    # an unfinished chunk + a measured time far from the plan -> re-plan first, then resubmit
    cfg_t = runs["run_t222"]
    wd_t = Path(cfg_t).parent / "campaign"
    res = sorted((wd_t / "results").glob("chunk_*.csv"))[0]
    keep = res.read_text()
    lines = keep.strip().split("\n")
    res.write_text("\n".join(lines[:1 + (len(lines) - 1) // 2]) + "\n")
    st = api.get("/api/run/state?fresh=1&config=" + cfg_t)
    ids = [a["id"] for a in st["attention"]]
    check(ids[:1] == ["replan"] and "resubmit" in ids, f"attention order: {ids}")
    check(st["next"].get("att") == "replan" and st["attention"][0]["action"]["preset"]["choice"] == "--force",
          f"next step = re-plan with the measured time: {st['next']['what']}")
    check([x for x in st["stages"] if x["key"] == "plan"][0]["pill"][1] == "Stale", "plan stage marked stale")
    res.write_text(keep)
    st = api.get("/api/run/state?fresh=1&config=" + cfg_t)
    check(not [a for a in st["attention"] if a["id"] in ("replan", "resubmit")], "attention clears when the chunk is complete again")
    # a run made without the app: open its workdir
    r = api.get("/api/runs/read?folder=" + str(wd_t))
    check(r.get("ok") and r["config_path"] == cfg_t, f"workdir opens its run: {r.get('config_path')}")
    st = api.get("/api/run/state?config=" + str(wd_t))
    check(st.get("ok") and st["config"] == cfg_t, "state of a workdir = state of its run")

    print("== reproducibility check")
    c = api.get("/api/runs/check?folder=" + str(W / "run_t222"))
    check(c["ok"] and c["identical"], "fresh run: inputs and scripts match the record")
    p = W / "run_t222" / "inputs" / "nep.txt"
    p.write_text("changed\n")
    c = api.get("/api/runs/check?folder=" + str(W / "run_t222"))
    check(c["ok"] and not c["identical"] and any(d["file"].endswith("nep.txt") for d in c["differences"]), "changed NEP file reported")
    rd = api.get("/api/runs/read?folder=" + str(W / "run_t222"))
    check(rd["ok"] and len(rd["record"]["history"]) >= 6, f"history of {len(rd['record']['history'])} commands in run_record.json")

    print("== cluster settings after creation, adopting a run made without the app")
    tom = Path(cfg_t)
    r = api.post("/api/runs/settings", {"config": cfg_t, "changes": {"cluster": {"walltime": "00:20:00", "sec_per_structure": "0.6"}},
                                        "note": "test"})
    txt = tom.read_text()
    check(r.get("ok") and set(r["changed"]) == {"cluster.walltime", "cluster.sec_per_structure"} and "00:20:00" in txt,
          f"settings changed: {r.get('changed')}")
    check(r.get("backup") and (tom.parent / r["backup"]).is_file() and "00:10:00" in (tom.parent / r["backup"]).read_text(), "old campaign.toml kept as .bak")
    rd = api.get("/api/runs/read?folder=" + str(tom.parent))
    h = rd["record"]["history"][-1]
    check(h.get("event") == "settings changed" and h["changes"]["cluster.walltime"] == ["00:10:00", "00:20:00"] and
          rd["record"]["config"]["cluster"]["sec_per_structure"] == 0.6, f"history entry: {h}")
    c = api.get("/api/runs/check?folder=" + str(tom.parent))
    check(not any(d["file"] == "campaign.toml" for d in c["differences"]), "record md5 follows the changed campaign.toml")
    code, rr = api.call("/api/runs/settings", {"config": cfg_t, "changes": {"free_energy": {"T": 500}}})
    check(code == 400 and "cannot be changed" in rr.get("error", ""), "physics settings cannot change after creation")
    code, rr = api.call("/api/runs/settings", {"config": cfg_t, "changes": {"cluster": {"walltime": "ten hours"}}})
    check(code == 400 and "walltime" in rr.get("error", ""), f"bad walltime refused: {rr.get('error')}")
    st = api.get("/api/run/state?fresh=1&config=" + cfg_t)
    check(st["cluster"]["walltime"] == "00:20:00", "state shows the new walltime")
    ad = W / "adopt_me"
    shutil.rmtree(ad, ignore_errors=True); ad.mkdir()
    (ad / "campaign.toml").write_text(txt)
    r = api.post("/api/runs/adopt", {"config": str(ad / "campaign.toml")})
    check(r.get("ok") and (ad / "run_record.json").is_file() and r["record"]["history"][0]["event"] == "adopted by the app",
          "a run made without the app gets run_record.json (adopted)")

    print("== results: one row per structure, New run estimate")
    rr = api.get("/api/results/rows?path=" + str(W / "run_t222") + "&sort=dF_meV&limit=500")
    st = api.get("/api/run/state?config=" + cfg_t)
    check(rr["ok"] and rr["all"] == st["counts"]["structures"] and len(rr["rows"]) == rr["total"],
          f"results rows: {rr.get('all')} structures = state count {st['counts']['structures']}")
    dfs = [x["dF_meV"] for x in rr["rows"] if x["dF_meV"] is not None]
    check(dfs == sorted(dfs) and min(dfs) == 0, "sorted by ΔF, the most stable structure of a composition has ΔF 0")
    lows = {}
    for x in rr["rows"]:
        if x["status"] == "clean":
            lows[x["comp"]] = min(lows.get(x["comp"], 9e9), x["F"])
    x0 = max((x for x in rr["rows"] if x["status"] == "clean"), key=lambda x: x["dF_meV"])
    check(abs(x0["dF_meV"] - 1000 * (x0["F"] - lows[x0["comp"]])) < 1e-3, "ΔF = F minus the lowest clean F of its composition")
    one = api.get("/api/results/rows?path=" + str(W / "run_t222") + "&comp=" + x0["comp"] + "&phase=" + x0["phase"])
    check(one["total"] == sum(1 for x in rr["rows"] if x["comp"] == x0["comp"] and x["phase"] == x0["phase"]) and
          all(x["comp"] == x0["comp"] and x["phase"] == x0["phase"] for x in one["rows"]), f"filter by composition and phase: {one['total']} rows")
    rp = x0.get("relaxed_path")
    body = b""
    if rp:
        rq = urllib.request.Request(api.base + "/api/file?download=1&path=" + urllib.request.quote(rp), headers={"X-Token": api.token})
        with urllib.request.urlopen(rq, timeout=60) as resp:
            body = resp.read()
    check(rp and "relaxed" in rp and body == Path(rp).read_bytes(), f"relaxed POSCAR of a row downloads: {rp}")
    nf = api.get("/api/results/rows?path=" + str(W / "run_t222") + "&status=failed")
    check(nf["ok"] and nf["total"] == 0 and nf["counts"].get("clean") == rr["all"], "status filter (failed: none) with counts of every status")
    rd = api.get("/api/runs/read?folder=" + str(W / "run_c222"))
    ec = json.loads(json.dumps(rd["config"]))
    base_c = W / "run_c222"
    ab = lambda p: p if not p or p.startswith("/") else str(base_c / p)
    ec["structure"]["poscar"] = ab(ec["structure"]["poscar"])
    ec["free_energy"]["nep"] = ab(ec["free_energy"]["nep"])
    ec["mapping"]["templates"] = [(t.split("=", 1)[0] + "=" + ab(t.split("=", 1)[1])) if "=" in t else ab(t) for t in ec["mapping"]["templates"]]
    est = api.post("/api/runs/estimate", {"config": ec})
    dry = subprocess.run([PY, str(base_c / "scripts" / "campaign.py"), str(base_c / "campaign.toml"), "dry-run"], capture_output=True, text=True, env=env).stdout
    import re as _re
    m_ch = _re.search(r"~(\d+) chunks", dry)
    m_ct = _re.search(r"TOTAL.*?(\d+)\s*$", dry, _re.M)
    m_cs = _re.search(r"([\d,.]+) core-hours", dry)
    check(est["ok"] and not est["problems"] and m_ct and est["structures"] == int(m_ct[1]) and m_ch and est["chunks"] == int(m_ch[1]) and
          m_cs and abs(est["core_h"] - float(m_cs[1].replace(",", ""))) <= 0.051 + 0.01 * est["core_h"],
          f"New run estimate = campaign.py dry-run: {est.get('structures')} structures, {est.get('chunks')} chunks, {est.get('core_h')} core-h "
          f"(dry-run {m_ct and m_ct[1]}, {m_ch and m_ch[1]}, {m_cs and m_cs[1]})")
    # the same with a big sec_per_structure and no prepared folders: pure prediction, compare the numbers closely
    sys.path.insert(0, str(ROOT))
    from studio import runs as R_
    ec2 = json.loads(json.dumps(ec))
    ec2["cluster"].update({"sec_per_structure": 3600.0, "walltime": "02:00:00", "cores_per_node": 8})
    ec2["structure"]["compositions"] = "0:32:4"
    ec2["paths"] = {"scripts": str(base_c / "scripts"), "workdir": "est_wd_none"}
    (W / "est_check.toml").write_text(R_.to_toml(R_.coerce(ec2)))
    est2 = api.post("/api/runs/estimate", {"config": ec2})
    dry2 = subprocess.run([PY, str(base_c / "scripts" / "campaign.py"), str(W / "est_check.toml"), "dry-run"], capture_output=True, text=True, env=env).stdout
    m_ch, m_ct, m_cs = _re.search(r"~(\d+) chunks", dry2), _re.search(r"TOTAL.*?(\d+)\s*$", dry2, _re.M), _re.search(r"([\d,.]+) core-hours", dry2)
    check(m_ct and m_ch and m_cs and est2["structures"] == int(m_ct[1]) and est2["chunks"] == int(m_ch[1]) and
          abs(est2["core_h"] - float(m_cs[1].replace(",", ""))) <= 0.051,
          f"estimate before prepare = dry-run: {est2['structures']} structures, {est2['chunks']} chunks, {est2['core_h']:.1f} core-h "
          f"(dry-run {m_ct and m_ct[1]}, {m_ch and m_ch[1]}, {m_cs and m_cs[1]})")
    e2 = api.post("/api/runs/estimate", {"config": {**ec, "free_energy": {**ec["free_energy"], "nep": ""}}})
    check(any("NEP model" in p for p in e2["problems"]), "estimate: missing NEP model reported")
    code, rr2 = api.call("/api/runs/create", {"folder": str(W / "run_nonep"), "config": {**ec, "free_energy": {**ec["free_energy"], "nep": ""}}})
    check(code == 400 and "NEP" in rr2.get("error", ""), "create refuses a run without a NEP model")
    xr = api.get("/api/xrd/rows?source=" + str(W / "run_t222"))
    check(xr.get("host") == "Zr" and xr.get("dopant") == "Ce", f"XRD page gets host/dopant from the run: {xr.get('host')} {xr.get('dopant')}")

    print("== input validation, double launch, stop")
    base = {"structure": {"poscar": str(W / "POSCAR_c_1x1x2"), "compositions": "0 4 8", "max_orderings": 2},
            "mapping": {"templates": [f"t={refs / 'POSCAR_t'}"]}, "free_energy": {"nep": str(W / "nep.txt")}}
    r = api.post("/api/runs/create", {"folder": str(W / "run_val"), "config": base})
    check(r.get("ok") and "compositions        = [0, 4, 8]" in (W / "run_val" / "campaign.toml").read_text().replace("  ", " ") or
          "[0, 4, 8]" in (W / "run_val" / "campaign.toml").read_text(), "space-separated compositions accepted")
    bad = json.loads(json.dumps(base)); bad["structure"]["compositions"] = "0-8"
    r = api.post("/api/runs/create", {"folder": str(W / "run_bad"), "config": bad})
    check(not r.get("ok") and "compositions" in r.get("error", ""), f"bad compositions refused: {r.get('error', '')[:60]}")
    r = api.post("/api/runs/create", {"folder": str(W / "run_c112"), "config": base})
    check(not r.get("ok"), f"existing run not overwritten: {r.get('error', '')[:60]}")
    (W / "notarun" / "scripts").mkdir(parents=True)
    r = api.post("/api/runs/create", {"folder": str(W / "notarun"), "config": base, "overwrite": True})
    check(not r.get("ok") and (W / "notarun" / "scripts").is_dir(), "foreign scripts/ folder left alone")
    two = json.loads(json.dumps(base)); two["mapping"]["templates"] = [f"t={refs / 'POSCAR_t'}", f"m={refs / 'POSCAR_t'}"]
    r = api.post("/api/runs/create", {"folder": str(W / "run_two"), "config": two})
    rt = (W / "run_two" / "campaign.toml").read_text() if r.get("ok") else ""
    check(r.get("ok") and "inputs/refs/t/POSCAR_t" in rt and "inputs/refs/m/POSCAR_t" in rt, "same file name for two phases kept apart")
    cfgv = str(W / "run_val" / "campaign.toml")
    stage(api, cfgv, "prepare")
    stage(api, cfgv, "plan")
    r1 = api.post("/api/stage", {"module": "compute", "mode": "local", "config": cfgv, "options": {}})
    r2 = api.post("/api/stage", {"module": "compute", "mode": "local", "config": cfgv, "options": {}})
    check(r1.get("ok") and not r2.get("ok"), f"second launch of the same stage refused: {r2.get('error', '')[:50]}")
    time.sleep(2)
    c = api.post("/api/jobs/cancel", {"id": r1["job"]})
    j = api.wait(r1["job"], 60)
    check(c.get("ok") and j["state"] == "cancelled", f"stop: job state {j['state']}")

    print("== jobs list")
    jl = api.get("/api/jobs")["jobs"]
    check(len(jl) > 20 and all(j["state"] in ("done", "failed", "cancelled") for j in jl), f"{len(jl)} jobs listed, none left running")


if __name__ == "__main__":
    main()
