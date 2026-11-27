#!/usr/bin/env bash
# batch_map.py test (~1-2 min, numpy + scipy only):
#   * host cells of several shapes: 2x2x2 (96 atoms), 1x2x4 (96), 3x2x2 (144), 1x1x2 (24)
#   * all four phases: refs/POSCAR_c, _t, _oIII, _m
#   * new in-process workers give byte-identical POSCARs and summary.csv to the old --subprocess mode
#   * every mapped structure IDENTICAL to its parent (compare_ordering.py)
#   * a bad input file is reported as FAILED and the others carry on (exit status 1)
#   * --progress lines, summary.csv saved while running, kill (no workers left behind) + restart with --skip-existing
#
#   bash tests/test_batch_map.sh            # from the package root
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
S="$ROOT/scripts"
W="$ROOT/tests/work_batch_map"
rm -rf "$W" && mkdir -p "$W" && cd "$W"
PY=python3
TPL="$ROOT/refs/POSCAR_c $ROOT/refs/POSCAR_t $ROOT/refs/POSCAR_oIII $ROOT/refs/POSCAR_m"

$PY - <<'EOF'
import numpy as np
FCC = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
TET = np.array([[a, b, c] for a in (.25, .75) for b in (.25, .75) for c in (.25, .75)])
for n in [(2, 2, 2), (1, 2, 4), (3, 2, 2), (1, 1, 2)]:
    nv = np.array(n)
    cells = np.array([[i, j, k] for i in range(n[0]) for j in range(n[1]) for k in range(n[2])])
    cat = ((cells[:, None] + FCC[None]) / nv).reshape(-1, 3)
    ox = ((cells[:, None] + TET[None]) / nv).reshape(-1, 3)
    with open("POSCAR_%dx%dx%d" % n, "w") as f:
        f.write("ZrO2 cubic %dx%dx%d\n1.0\n" % n)
        for v in np.diag(nv * 5.10):
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write("  Zr O\n  %d %d\nDirect\n" % (len(cat), len(ox)))
        for p in np.vstack([cat, ox]):
            f.write("  %.10f %.10f %.10f\n" % tuple(p % 1))
EOF

fail() { echo "FAIL: $*"; exit 1; }

for shape in ${SHAPES:-2x2x2 1x2x4 3x2x2 1x1x2}; do
  N=$(awk 'NR==7{print $1}' POSCAR_$shape)
  $PY "$S/enumerate_substitutions.py" POSCAR_$shape --host Zr --dopant Ce --n $((N / 4)) --random 3 --seed 7 \
      --outdir "e_$shape" > /dev/null
  D=$(ls -d e_$shape/*/)
  # new mode (in-process workers) vs old mode (one Python per mapping)
  $PY "$S/batch_map.py" "$D" --templates $TPL --glob 'POSCAR_*' --jobs 4 --progress 1 > "new_$shape.log"
  $PY "$S/batch_map.py" "$D" --templates $TPL --glob 'POSCAR_*' --jobs 4 --subprocess --out-dir "old_$shape" \
      > "old_$shape.log"
  diff -r "${D}mapped" "old_$shape" > /dev/null || fail "$shape: in-process and --subprocess outputs differ"
  grep -q "^mapped 12/12 (100.0 %), 0 failed" "new_$shape.log" || fail "$shape: no final progress line"
  nid=0
  for f in "${D}"mapped/*/POSCAR_*; do
    b=$(basename "$f"); parent="${D}POSCAR_${b#POSCAR_*_}"
    $PY "$S/compare_ordering.py" "$parent" "$f" | grep -q "VERDICT: IDENTICAL" || fail "$f not IDENTICAL to $parent"
    nid=$((nid + 1))
  done
  [ "$nid" = 12 ] || fail "$shape: $nid mapped files, expected 12"
  echo "  ok  $shape ($N cations, Ce $((N / 4))): c/t/oIII/m x 3 orderings, same as --subprocess, 12/12 IDENTICAL"
  $PY - "${D}mapped/summary.csv" <<'EOF'
import csv, sys
for ph in ("c", "t", "oIII", "m"):
    rr = [r for r in csv.DictReader(open(sys.argv[1])) if r["target"] == ph]
    print("      %-5s rms shift %s A, strain S %s, equally good variants %s" % (
        ph, "/".join(r["rms_shift_A"] for r in rr), "/".join(r["strain_S"] for r in rr), rr[0]["tied_variants"]))
EOF
done

# bad input: reported, others carry on, exit 1
D=$(ls -d e_*/*/ | head -1)
echo "not a poscar" > "${D}POSCAR_9999_bad"
set +e
$PY "$S/batch_map.py" "$D" --templates "$ROOT/refs/POSCAR_t" --glob 'POSCAR_*' --jobs 2 --out-dir bad > bad.log
rc=$?
set -e
rm "${D}POSCAR_9999_bad"
[ $rc = 1 ] || fail "bad input: exit status $rc, expected 1"
grep -q "^done: 3 ok, 1 failed" bad.log || fail "bad input: wrong summary line"
grep -q "POSCAR_9999_bad.*FAILED.*not a VASP5 POSCAR" bad.log || fail "bad input: failure not reported"
echo "  ok  bad input reported as FAILED, the other 3 mapped, exit status 1"

# kill + restart: summary.csv is saved while running, --skip-existing continues from it
$PY "$S/enumerate_substitutions.py" POSCAR_3x2x2 --host Zr --dopant Ce --n 12 --random 40 --seed 9 \
    --outdir kill > /dev/null
D=$(ls -d kill/*/)
BATCH_MAP_SAVE_EVERY=0 $PY "$S/batch_map.py" "$D" --templates $TPL --glob 'POSCAR_*' --jobs 1 > kill1.log &
pid=$!
for _ in $(seq 100); do
  [ -f "${D}mapped/summary.csv" ] && [ "$(wc -l < "${D}mapped/summary.csv")" -ge 6 ] && break
  sleep 0.2
done
kill $pid 2>/dev/null || fail "kill test: batch_map finished before it could be killed"
wait $pid 2>/dev/null || true
grep -q "^stopped (signal 15) after" kill1.log || fail "kill test: no 'stopped' message"
sleep 1
if pgrep -f "batch_map.py ${D}" > /dev/null; then pkill -f "batch_map.py ${D}"; fail "kill test: workers left running"; fi
done1=$(($(wc -l < "${D}mapped/summary.csv") - 1))
[ "$done1" -ge 5 ] && [ "$done1" -lt 160 ] || fail "kill test: $done1 rows saved"
$PY "$S/batch_map.py" "$D" --templates $TPL --glob 'POSCAR_*' --jobs 4 --skip-existing --progress 5 > kill2.log
grep -q "^$done1 of 160 runs already done" kill2.log || fail "restart did not reuse the $done1 saved runs"
grep -q "^done: 160 ok, 0 failed" kill2.log || fail "restart did not finish all 160 runs"
[ "$(wc -l < "${D}mapped/summary.csv")" = 161 ] || fail "restart: summary.csv rows"
echo "  ok  killed after $done1 of 160 runs; --skip-existing redid only the other $((160 - done1))"

echo "test_batch_map.sh: all checks passed"
