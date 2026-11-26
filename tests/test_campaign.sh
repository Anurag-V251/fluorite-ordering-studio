#!/usr/bin/env bash
# End-to-end test of campaign.py (no NEP model, no cluster needed), ~2-4 min:
#   * synthetic host cells of several shapes: 2x2x2 (96 atoms), 1x2x4 (96), 3x2x2 (144), 1x1x2 (24),
#     and a tetragonal (non-cubic) 2x2x2 parent
#   * synthetic t, oIII (Pca2_1) and m (P2_1/c) templates
#   * end members (Ce = 0 and Ce = all), compositions with fewer orderings than requested
#   * every mapped structure checked with compare_ordering.py (must be IDENTICAL)
#   * free energies with a MOCK calorine (Einstein crystal: F_vib checked against the analytic value)
#   * kill + restart of a chunk, plan --append after raising max_orderings, SLURM scripts run by a fake sbatch
#     (array split at MaxArraySize, --unfinished, merge dependency), flags for imaginary modes / unconverged relax
#
#   bash tests/test_campaign.sh            # from the package root
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
S="$ROOT/scripts"
W="$ROOT/tests/work_campaign"
rm -rf "$W" && mkdir -p "$W/bin" && cd "$W"
export PYTHONPATH="$ROOT/tests/mock_calorine"
export OMP_NUM_THREADS=1
PY=python3

# ---------------------------------------------------------------- fake sbatch (runs array tasks right away)
cat > bin/sbatch <<'EOF'
#!/usr/bin/env python3
import os, random, subprocess, sys
args = sys.argv[1:]; script = args[-1]; env = dict(os.environ); arr = None
for a in args[:-1]:
    if a.startswith("--array="):
        arr = a.split("=", 1)[1].split("%")[0]
    if a.startswith("--export="):
        for kv in a.split("=", 1)[1].split(","):
            if "=" in kv:
                k, v = kv.split("=", 1); env[k] = v
idx = []
for part in (arr.split(",") if arr else []):
    a, _, b = part.partition("-"); idx += range(int(a), int(b or a) + 1)
jid = str(random.randint(1000, 9999))
with open(os.environ["FAKE_SBATCH_LOG"], "a") as log:
    log.write(" ".join(sys.argv[1:]) + "\n")
    for i in idx or [None]:
        e = dict(env, SLURM_JOB_ID=jid, **({} if i is None else {"SLURM_ARRAY_TASK_ID": str(i)}))
        r = subprocess.run(["bash", script], env=e, capture_output=True, text=True)
        log.write(f"  task {i}: exit {r.returncode}\n")
        if r.returncode:
            log.write(r.stdout[-2000:] + r.stderr[-2000:])
print(jid)
EOF
chmod +x bin/sbatch
export PATH="$W/bin:$PATH" FAKE_SBATCH_LOG="$W/sbatch.log"

# ---------------------------------------------------------------- inputs
$PY - <<'EOF'
import os
import numpy as np
FCC = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
TET = np.array([[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)])

def write(fn, title, lat, names, counts, frac):
    with open(fn, "w") as f:
        f.write(title + "\n1.0\n")
        for v in lat:
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write("  " + " ".join(names) + "\n  " + " ".join(map(str, counts)) + "\nDirect\n")
        for p in frac:
            f.write("  %.10f %.10f %.10f\n" % tuple(np.array(p) % 1))

os.makedirs("refs", exist_ok=True)
o_t = [[q[0], q[1], q[2] + (0.045 if (round(q[0] * 4) + round(q[1] * 4)) % 4 == 0 else -0.045)] for q in TET]
write("refs/POSCAR_t", "synthetic t template", np.diag([5.07, 5.07, 5.18]), ["Zr", "O"], [4, 8], np.vstack([FCC, o_t]))
a, b, c, be = 5.1505, 5.2116, 5.3173, np.radians(99.23)                      # m-ZrO2-like, P2_1/c
ops = lambda p: [(p[0], p[1], p[2]), (-p[0], p[1] + .5, -p[2] + .5), (-p[0], -p[1], -p[2]), (p[0], -p[1] + .5, p[2] + .5)]
write("refs/POSCAR_m", "synthetic m template", np.array([[a, 0, 0], [0, b, 0], [c * np.cos(be), 0, c * np.sin(be)]]),
      ["Zr", "O"], [4, 8], ops((0.2758, 0.0404, 0.2089)) + ops((0.0703, 0.3359, 0.3406)) + ops((0.4423, 0.7549, 0.4789)))
ops = lambda p: [(p[0], p[1], p[2]), (-p[0], -p[1], p[2] + .5), (p[0] + .5, -p[1], p[2]), (-p[0] + .5, p[1], p[2] + .5)]
write("refs/POSCAR_oIII", "synthetic oIII template (Pca2_1, polar c)", np.diag([5.08, 5.26, 5.07]), ["Zr", "O"], [4, 8],
      ops((0.03, 0.27, 0.25)) + ops((0.28, 0.52, 0.53)) + ops((0.23, 0.03, 0.47)))
for n in [(2, 2, 2), (1, 2, 4), (3, 2, 2), (1, 1, 2)]:
    nv = np.array(n)
    cells = np.array([[i, j, k] for i in range(n[0]) for j in range(n[1]) for k in range(n[2])])
    cat = ((cells[:, None] + FCC[None]) / nv).reshape(-1, 3)
    ox = ((cells[:, None] + TET[None]) / nv).reshape(-1, 3)
    write("POSCAR_%dx%dx%d" % n, "ZrO2 cubic %dx%dx%d" % n, np.diag(nv * 5.10), ["Zr", "O"], [len(cat), len(ox)],
          np.vstack([cat, ox]))
open("nep.txt", "w").write("mock\n")
EOF
$PY "$S/map_fluorite.py" POSCAR_2x2x2 --template refs/POSCAR_t --out POSCAR_t2x2x2 --quiet      # tetragonal parent

config() {   # name host compositions max_orderings [extra cluster lines]
cat > "$1.toml" <<EOF
[structure]
poscar = "$2"
compositions = $3
max_orderings = $4
seed = 3
[paths]
scripts = "$S"
workdir = "$1"
[mapping]
templates = ["t=refs/POSCAR_t", "oIII=refs/POSCAR_oIII", "m=refs/POSCAR_m"]
jobs = 4
[free_energy]
nep = "nep.txt"
mesh = 8
fmax = 1e-3
[cluster]
cores_per_node = 4
walltime = "00:02:00"
sec_per_structure = 10
max_concurrent = 8
setup = "export PYTHONPATH=$PYTHONPATH"
EOF
if [ -n "${5:-}" ]; then                       # an extra [cluster] line overrides the default above
  sed -i "/^${5%% =*} = /d" "$1.toml"; echo "$5" >> "$1.toml"
fi
}

check_identical() {   # every verify.csv row IDENTICAL, count = orderings x 3
  $PY - "$1" <<'EOF'
import csv, glob, sys
n = bad = 0
for f in glob.glob(sys.argv[1] + "/*/mapped/verify.csv"):
    for r in csv.DictReader(open(f)):
        n += 1; bad += r["verdict"] != "IDENTICAL"
assert n > 0 and bad == 0, (n, bad)
print(f"   compare_ordering: {n}/{n} IDENTICAL")
EOF
}

check_fvib() {        # mock = Einstein crystal: F_vib known analytically
  $PY - "$1" <<'EOF'
import csv, re, sys
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.units import _amu, _e, _hbar, _k
K = {"Zr": 10.0, "Ce": 8.0, "O": 5.0}
rows = list(csv.DictReader(open(sys.argv[1])))
err = 0
for r in rows:
    T = float(r["T_K"]); comp = re.findall(r"([A-Z][a-z]?)(\d+)", r["composition"])
    F = sum(int(c) * 3 * _k * T * np.log(2 * np.sinh(_hbar * np.sqrt(K[s] * _e * 1e20 / (atomic_masses[atomic_numbers[s]] * _amu)) / (2 * _k * T)))
            for s, c in comp) / _e / sum(int(c) for _, c in comp)
    err = max(err, abs(float(r["F_vib_eV_atom"]) - F))
assert err < 1e-5, err
print(f"   F_vib vs analytic Einstein value: max error {err:.1e} eV/atom over {len(rows)} rows")
EOF
}

# ---------------------------------------------------------------- 1. shapes: full local chain
for shape in ${SHAPES:-2x2x2 1x2x4 3x2x2 1x1x2 t2x2x2}; do
  host=POSCAR_$shape; [ "$shape" = t2x2x2 ] && host=POSCAR_t2x2x2
  N=$(awk 'NR==7{print $1}' $host)
  echo "== shape $shape ($N cation sites): compositions 0, 1, 2, N/2, N-1, N; 3 orderings each"
  config "c_$shape" "$host" "[0, 1, 2, $((N / 2)), $((N - 1)), $N]" 3
  $PY "$S/campaign.py" "c_$shape.toml" dry-run > "dry_$shape.log"; sed -n 3p "dry_$shape.log"
  $PY "$S/campaign.py" "c_$shape.toml" prepare > "prepare_$shape.log"
  grep -h "mapped\|compare_ordering" "prepare_$shape.log" | sort | uniq -c | sort -rn | head -4 | sed 's/^/  /'
  check_identical "c_$shape"
  $PY "$S/campaign.py" "c_$shape.toml" plan > "plan_$shape.log"; head -1 "plan_$shape.log"
  n_chunks=$($PY -c "import json;print(len(json.load(open('c_$shape/chunks/plan.json'))['chunks']))")
  for ((k = 0; k < n_chunks; k++)); do $PY "$S/campaign.py" "c_$shape.toml" run --chunk $k > "run_${shape}_$k.log" 2>&1; done
  $PY "$S/campaign.py" "c_$shape.toml" merge > "merge_$shape.log"; head -1 "merge_$shape.log"
  check_fvib "c_$shape/free_energy_300K_all.csv"
done
# end members exist and have no Zr / no Ce
for d in c_*/Ce_0; do [ "$(sed -n 6p $d/POSCAR_0001_WC+nan | xargs)" = "Zr O" ] || { echo "$d species wrong"; exit 1; }; done
for d in c_*/Ce_100; do [ "$(sed -n 6p $d/POSCAR_0001_WC+nan | xargs)" = "Ce O" ] || { echo "$d species wrong"; exit 1; }; done
echo "   end members ok (Ce_0 = pure ZrO2, Ce_100 = pure CeO2, alpha = nan)"

# ---------------------------------------------------------------- 2. random draws: distinct, seeded, capped
echo "== random orderings: distinct, reproducible, capped at C(N,n)"
config c_rand POSCAR_2x2x2 "[1, 2, 6]" 40
$PY "$S/campaign.py" c_rand.toml prepare > prepare_rand.log
config c_rand2 POSCAR_2x2x2 "[6]" 40
$PY "$S/campaign.py" c_rand2.toml prepare > /dev/null
$PY - <<'EOF'
import csv
r = lambda f: [x["substituted_atom_indices"] for x in csv.DictReader(open(f))]
a, b = r("c_rand/Ce_18.75/summary.csv"), r("c_rand2/Ce_18.75/summary.csv")
assert len(a) == 40 and len(set(a)) == 40, "not 40 distinct"
assert sorted(a) == sorted(b), "same seed gave different orderings"
assert len(r("c_rand/Ce_3.12/summary.csv")) == 32, "n=1 must give all 32"
print("   n=6: 40 distinct, identical with the same seed; n=1: all 32 (fewer than requested)")
EOF

# ---------------------------------------------------------------- 3. kill + restart of a chunk
echo "== restart after an interrupted chunk"
config c_rs POSCAR_2x2x2 "[6]" 6 'walltime = "01:00:00"'
$PY "$S/campaign.py" c_rs.toml prepare > /dev/null && $PY "$S/campaign.py" c_rs.toml plan > /dev/null
( $PY "$S/campaign.py" c_rs.toml run --chunk 0 > rs1.log 2>&1 & )
for i in $(seq 1 60); do [ -f c_rs/results/chunk_00000.csv ] && [ "$(wc -l < c_rs/results/chunk_00000.csv)" -ge 5 ] && break; sleep 0.5; done
pkill -f "batch_free_energy.py.*c_rs" || true; sleep 1
before=$(($(wc -l < c_rs/results/chunk_00000.csv) - 1))
$PY "$S/campaign.py" c_rs.toml run --chunk 0 > rs2.log 2>&1
grep -o "[0-9]* already done" rs2.log | sed "s/^/   killed after $before rows; rerun: /"
st=$($PY "$S/campaign.py" c_rs.toml status); echo "   ${st%%$'\n'*}"
[[ "$st" == *"24 ok, 0 failed, 0 not run"* ]] || { echo "restart failed"; exit 1; }

# ---------------------------------------------------------------- 3b. cancel, re-plan with --force, merge keeps everything
echo "== cancel mid-chunk, plan --force with a new sec_per_structure, merge keeps earlier results"
config c_rp POSCAR_2x2x2 "[6, 16]" 6 'walltime = "01:00:00"'
$PY "$S/campaign.py" c_rp.toml prepare > /dev/null && $PY "$S/campaign.py" c_rp.toml plan > /dev/null
total=$($PY -c "import json;print(sum(c['n'] for c in json.load(open('c_rp/chunks/plan.json'))['chunks']))")
( $PY "$S/campaign.py" c_rp.toml run --chunk 0 > rp1.log 2>&1 & )
for i in $(seq 1 60); do [ -f c_rp/results/chunk_00000.csv ] && [ "$(wc -l < c_rp/results/chunk_00000.csv)" -ge 8 ] && break; sleep 0.5; done
pkill -f "batch_free_energy.py.*c_rp" || true; sleep 1
sed -i 's/sec_per_structure = .*/sec_per_structure = 3/' c_rp.toml
grep -q "sec_per_structure = 3" c_rp.toml || sed -i 's/^\[cluster\]/[cluster]\nsec_per_structure = 3/' c_rp.toml
$PY "$S/campaign.py" c_rp.toml plan --force > rp_plan.log; head -1 rp_plan.log | sed 's/^/   /'
n_chunks=$($PY -c "import json;print(len(json.load(open('c_rp/chunks/plan.json'))['chunks']))")
for ((k = 0; k < n_chunks; k++)); do $PY "$S/campaign.py" c_rp.toml run --chunk $k > "rp_run_$k.log" 2>&1; done
$PY "$S/campaign.py" c_rp.toml status | grep "earlier" | sed 's/^/   /' || true
$PY "$S/campaign.py" c_rp.toml merge > rp_merge.log; grep "earlier" rp_merge.log | sed 's/^/   /' || true
$PY - "$total" <<'PYEOF'
import csv, sys
rows = list(csv.DictReader(open("c_rp/free_energy_300K_all.csv")))
files = [r["file"] for r in rows]
assert len(files) == len(set(files)), "duplicate rows in merge"
assert len(rows) == int(sys.argv[1]) and all(r["status"] == "ok" for r in rows), (len(rows), sys.argv[1])
print(f"   merged {len(rows)}/{sys.argv[1]} structures, all ok, no duplicates")
PYEOF
check_fvib c_rp/free_energy_300K_all.csv

# ---------------------------------------------------------------- 4. SLURM: arrays, MaxArraySize split, append, --unfinished
echo "== SLURM scripts with a fake sbatch"
config c_sl POSCAR_2x2x2 '"0:32:8"' 2 'max_array = 1'
$PY "$S/campaign.py" c_sl.toml slurm prepare --submit | grep sbatch
$PY "$S/campaign.py" c_sl.toml plan > sl_plan.log; head -1 sl_plan.log
$PY "$S/campaign.py" c_sl.toml slurm compute --submit --with-merge | grep "sbatch\|merge"
sed -i 's/max_orderings = 2/max_orderings = 3/' c_sl.toml
$PY "$S/campaign.py" c_sl.toml slurm prepare --unfinished --submit | grep sbatch
$PY "$S/campaign.py" c_sl.toml plan --append > sl_plan2.log; head -1 sl_plan2.log
$PY "$S/campaign.py" c_sl.toml slurm compute --unfinished --submit --with-merge | grep "sbatch\|merge"
grep -q "exit [1-9]" sbatch.log && { cat sbatch.log; echo "an array task failed"; exit 1; }
st=$($PY "$S/campaign.py" c_sl.toml status); echo "   ${st%%$'\n'*}"
[[ "$st" == *"0 failed, 0 not run"* ]] || { echo "slurm flow incomplete"; exit 1; }
check_fvib c_sl/free_energy_300K_all.csv

# ---------------------------------------------------------------- 5. flags
echo "== flags: imaginary modes and unconverged relaxation"
config c_fl POSCAR_2x2x2 "[6]" 1
sed -i 's/fmax = 1e-3/fmax = 1e-7\nsteps = 2/' c_fl.toml
$PY "$S/campaign.py" c_fl.toml prepare > /dev/null && $PY "$S/campaign.py" c_fl.toml plan > /dev/null
MOCK_NEP_SOFT=1 $PY "$S/campaign.py" c_fl.toml run --chunk 0 > fl.log 2>&1
$PY "$S/campaign.py" c_fl.toml merge > fl_merge.log; sed -n '2,3p' fl_merge.log | sed 's/^/   /'
$PY - <<'EOF'
import csv
rows = list(csv.DictReader(open("c_fl/free_energy_300K_all.csv")))
assert all(r["relax_converged"] == "False" and int(r["n_imaginary_modes"]) > 0 and "NOT converged" in r["note"] for r in rows)
print(f"   {len(rows)}/{len(rows)} rows flagged (relax_converged=False, n_imaginary_modes>0, note)")
EOF

# ---------------------------------------------------------------- 6. rerun_flagged.py
echo "== rerun_flagged.py: flagged rows recomputed in rounds, merge keeps the newest row"
config c_rr POSCAR_2x2x2 "[6, 16]" 3
sed -i 's/fmax = 1e-3/fmax = 1e-7\nsteps = 2/' c_rr.toml
$PY "$S/campaign.py" c_rr.toml prepare > /dev/null && $PY "$S/campaign.py" c_rr.toml plan > /dev/null
n_chunks=$($PY -c "import json;print(len(json.load(open('c_rr/chunks/plan.json'))['chunks']))")
for ((k = 0; k < n_chunks; k++)); do MOCK_NEP_SOFT=1 $PY "$S/campaign.py" c_rr.toml run --chunk $k > "rr_run_$k.log" 2>&1; done
$PY "$S/campaign.py" c_rr.toml merge > /dev/null
total=$($PY -c "import csv;print(len(list(csv.DictReader(open('c_rr/free_energy_300K_all.csv')))))")
$PY "$S/rerun_flagged.py" c_rr.toml --dry-run > rr_dry.log; sed -n 2p rr_dry.log | sed 's/^/   /'
grep -q "selected (both): $total " rr_dry.log || { cat rr_dry.log; echo "selection wrong"; exit 1; }
$PY "$S/rerun_flagged.py" c_rr.toml --what imaginary --freq-window -0.5 0 --dry-run > rr_win.log
grep -q "selected (imaginary, min freq in \[-0.5, 0.0\] THz): 0 " rr_win.log || { cat rr_win.log; echo "window"; exit 1; }
# round 1: same settings, still soft -> still flagged
MOCK_NEP_SOFT=1 $PY "$S/rerun_flagged.py" c_rr.toml --same --local --jobs 4 > rr1.log 2>&1
grep "start geometries" rr1.log | sed 's/^/   round 1 /'
# round 2: through (fake) SLURM, 2 array tasks, normal force field + converging settings -> clean
$PY "$S/rerun_flagged.py" c_rr.toml --same --fmax 1e-3 --steps 2000 --nodes 2 --submit > rr2.log
grep "start geometries\|submitted" rr2.log | sed 's/^/   round 2 /'
grep -q "from round 1" rr2.log || { cat rr2.log; echo "round 2 did not start from round 1 geometries"; exit 1; }
grep -q "exit [1-9]" sbatch.log && { cat sbatch.log; echo "a rerun task failed"; exit 1; }
$PY - "$total" <<'PYEOF'
import csv, sys
rows = list(csv.DictReader(open("c_rr/free_energy_300K_all.csv")))
files = [r["file"] for r in rows]
assert len(files) == len(set(files)) == int(sys.argv[1]), (len(files), len(set(files)), sys.argv[1])
bad = [r for r in rows if r["status"] != "ok" or int(r["n_imaginary_modes"]) > 0 or r["relax_converged"] != "True"]
assert not bad, bad[0]
print(f"   after round 2: {len(rows)}/{sys.argv[1]} structures, all ok, 0 imaginary, 0 unconverged, no duplicates")
PYEOF
check_fvib c_rr/free_energy_300K_all.csv
$PY "$S/rerun_flagged.py" c_rr.toml --dry-run > rr3.log; grep -q "nothing to rerun" rr3.log || { cat rr3.log; echo "round 3 should find nothing"; exit 1; }
echo "   round 3: nothing left to rerun"

echo
echo "CAMPAIGN TEST PASSED - outputs in $W"
