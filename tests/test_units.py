#!/usr/bin/env python3
"""test_units.py - fast checks of app internals (no server, ~15 s):  python3 tests/test_units.py"""
import json
import os
import sys
import tempfile
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ["FO_STUDIO_KILL_GRACE"] = "2"
from studio import jobs as J, runs as R, scan as S  # noqa: E402

fails = 0


def check(c, msg):
    global fails
    print(("   ok   " if c else "   FAIL ") + msg)
    fails += not c


def err(fn):
    try:
        fn()
    except ValueError as e:
        return str(e)
    return None


T = Path(tempfile.mkdtemp(prefix="fo_units_"))
host = ROOT / "tests" / "work_studio" / "POSCAR_c_2x2x2"
if not host.is_file():                          # make a 2x2x2 cubic host from the 12-atom reference
    import subprocess
    host = T / "POSCAR_c_2x2x2"
    subprocess.run([sys.executable, "-c", f"""
from ase.io import read, write
write('{host}', read('{ROOT}/refs/POSCAR_c').repeat((2,2,2)), format='vasp', direct=True)"""], check=True)
refs = ROOT / "refs"
nep = T / "nep.txt"
nep.write_text("dummy")
from studio.runs import peek_poscar  # noqa: E402
_sp = peek_poscar(host)["species"]
_h = "Zr" if "Zr" in _sp else [e for e in _sp if e != "O"][0]   # refs/ may hold CeO2 cells
base = {"structure": {"poscar": str(host), "compositions": "0 8", "host": _h, "dopant": "Ce" if _h != "Ce" else "Zr"},
        "mapping": {"templates": [f"t={refs}/POSCAR_t"]}, "free_energy": {"nep": str(nep)}}

print("== create: validation before anything is written")
f1 = T / "run_missing"
e = err(lambda: R.create_run({"folder": str(f1), "config": json.loads(json.dumps(base).replace("POSCAR_t", "POSCAR_nope"))}, ROOT / "scripts"))
check(e and "POSCAR_nope" in e and not f1.exists(), f"missing template refused, no folder left: {e}")
e = err(lambda: R.create_run({"folder": str(f1), "config": dict(base, mapping={"templates": [f"x={host}"]})}, ROOT / "scripts"))
check(e and "12-atom" in e, f"supercell as template refused: {e}")
e = err(lambda: R.create_run({"folder": str(f1), "config": dict(base, structure=dict(base["structure"], compositions="8:0"))}, ROOT / "scripts"))
check(e and "larger than stop" in e, f"empty range refused: {e}")
e = err(lambda: R.create_run({"folder": str(f1), "config": dict(base, free_energy={"nep": str(nep), "T": "inf"})}, ROOT / "scripts"))
check(e and "finite" in e, f"T = inf refused: {e}")
e = err(lambda: R.create_run({"folder": str(f1), "config": dict(base, paths={"workdir": "../x"})}, ROOT / "scripts"))
check(e and "workdir" in e, f"workdir outside the run refused: {e}")
out = R.create_run({"folder": str(f1), "name": "a\x7fb \"q\" \\ ü", "config": dict(base, cluster={"setup": "module load x\x7f\nsource y"})}, ROOT / "scripts")
check(bool(tomllib.loads((f1 / "campaign.toml").read_text(encoding="utf-8"))), "campaign.toml with odd characters parses")
check(R.split_template("/data/a=b/POSCAR_m") == ("m", "/data/a=b/POSCAR_m") and R.split_template("t=/x=y/P") == ("t", "/x=y/P"),
      "template 'label=path' split on the first '=' only before a label")
check(not any(p.stat().st_mode & 0o022 for p in [f1 / "scripts", *(f1 / "scripts").rglob("*")]),
      "frozen scripts are not group/other-writable (umask 002 sources)")
from studio import server as SV  # noqa: E402
check(SV.scripts_untrusted(f1 / "scripts") is None, "own frozen scripts are trusted")
os.chmod(f1 / "scripts" / "campaign.py", 0o666)
why = SV.scripts_untrusted(f1 / "scripts")
check(why and "chmod -R go-w" in why, f"world-writable script refused with the fix in the message: {why}")
os.chmod(f1 / "scripts" / "campaign.py", 0o644)
R.create_run({"folder": str(T / "run_sp"), "starting_phase": "t", "config": base}, ROOT / "scripts")
info = SV.run_info_of(T / "run_sp" / "campaign" / "free_energy_300K_all.csv")
check(info.get("starting_phase") == "t" and info.get("dopant") == base["structure"]["dopant"] and info.get("T") == "300",
      f"plots get dopant, T and starting phase from the run: {info}")
pz = T / "POSCAR_split"
pz.write_text("x\n1\n5 0 0\n0 5 0\n0 0 5\nZr O Zr O\n1 2 1 2\nDirect\n" + "0 0 0\n" * 6)
check(R.peek_poscar(pz)["species"] == {"Zr": 2, "O": 4}, "repeated element blocks add up")
for i in range(5):
    R.add_history(f1, {"event": f"e{i}"})
check(len(json.loads((f1 / "run_record.json").read_text())["history"]) == 5, "history entries kept")

print("== jobs: stop escalates to SIGKILL, ids, outputs frozen")
jobs = J.Jobs(T / "state")
jid = jobs.start(["bash", "-c", "trap '' TERM; echo submitted rerun array 1234; echo Submitted batch job 999; "
                  "echo submitted merge job 1000; sleep 60"], T, title="stubborn")
for _ in range(50):
    if jobs.get(jid)["state"] == "running":
        break
    time.sleep(0.1)
time.sleep(0.5)
check(jobs.cancel(jid), "stop signal sent")
for _ in range(80):
    g = jobs.get(jid)
    if g["state"] != "running":
        break
    time.sleep(0.1)
check(g["state"] == "cancelled", f"job that ignores SIGTERM ends as cancelled after the grace time: {g['state']}")
check(g["slurm_ids"] == ["999", "1000", "1234"], f"SLURM ids incl. rerun array, sorted as numbers: {g['slurm_ids']}")
check((jobs.dir / jid / "outputs.json").is_file(), "outputs list frozen at the end")
(jobs.dir / jid / "slurm_ids.json").write_text("")       # half-written cache: no crash
check(jobs.get(jid)["slurm_ids"] == ["999", "1000", "1234"], "empty cache file re-read from the log")
check(isinstance(jobs.list(), list) and "outputs" not in jobs.list()[0], "list() does not walk outputs")

print("== scan: streamed rows")
wd = T / "w"
(wd / "results").mkdir(parents=True)
(wd / "results" / "chunk_00000.csv").write_text("file,status,T_K,F_total_eV_atom,n_imaginary_modes,relax_converged,extra\n"
                                                "Ce_00/a,FAILED,300,,,,x\nCe_00/a,ok,300,-1.0,0,True,x\nCe_00/b,ok,300,-1.1,2,True,x\n")
b = {}
for r in S.iter_result_rows(wd / "results" / "chunk_00000.csv"):
    S.add_best(b, r)
check(b["Ce_00/a"]["status"] == "ok" and "extra" not in b["Ce_00/a"], "ok beats FAILED, unused columns dropped")

print("== state: sbatch output, array ids, squeue/sacct rows, cost model")
from studio import state as STATE  # noqa: E402
log = """wrote /w/slurm/compute.sh
  $ sbatch --parsable --array=0-4,9%50 --export=ALL,CHUNK_OFFSET=1000 compute.sh
  submitted job 777
  $ sbatch --parsable --dependency=afterany:777 merge.sh
  submitted merge job 778 (after 777)
round 2: 3 structures in 2 task(s)
  $ sbatch --parsable --array=0-1 rerun_002.sh
  submitted rerun array 779
"""
sub = STATE.submissions(log, "J", 0)
check([(x["id"], x["kind"], x["offset"], x["tasks"]) for x in sub] ==
      [("777", "compute", 1000, [0, 1, 2, 3, 4, 9]), ("778", "merge", 0, None), ("779", "rerun", 0, [0, 1])],
      f"submissions parsed: {[(x['id'], x['kind'], x['tasks']) for x in sub]}")
rows = STATE._slurm_rows([{"id": "777_[5-9%50]", "state": "PENDING"}, {"id": "777_3", "state": "RUNNING", "time": "1-00:00:05"}],
                         [{"id": "777_0", "state": "TIMEOUT", "elapsed": "02:00:00"}, {"id": "777_3", "state": "RUNNING"}])
check(rows["777"][0]["state"] == "TIMEOUT" and rows["777"][9]["state"] == "PENDING" and rows["777"][3]["state"] == "RUNNING"
      and STATE._tsec(rows["777"][3]["time"]) == 86405, "squeue wins over sacct, pending ranges expanded, D-HH:MM:SS read")
check(STATE.array_ids("[0-6:2]") == [0, 2, 4, 6] and STATE.comp_list("0:8:4", 32) == [0, 4, 8] and STATE.comp_list("all", 2) == [0, 1, 2],
      "array steps and composition specs")
import numpy as np  # noqa: E402
lat = [[5.1, 0.2, 0], [0, 5.3, 0.1], [0.3, 0, 10.4]]
inv = np.linalg.inv(np.array(lat))
want = [max(1, int(np.ceil(10.0 / w - 1e-6))) for w in 1.0 / np.linalg.norm(inv, axis=0)]
check(STATE.phonon_sc({"min_length": 10.0}, lat) == want, f"phonon supercell = campaign.py's numpy version {want}")

print("UNIT TESTS FAILED" if fails else "UNIT TESTS PASSED")
sys.exit(1 if fails else 0)
