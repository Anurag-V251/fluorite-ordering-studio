#!/usr/bin/env bash
# Smoke test: runs the whole pipeline on a small synthetic ZrO2 cell.
# Checks the installation and that the scripts work together. The free energies are PLACEHOLDERS
# (the NEP/phonon step needs your model); with calorine + phonopy installed and NEP=path/to/nep.txt set,
# the real batch_free_energy.py step is run on two structures as well.
#
#   bash tests/smoke_test.sh            # from the project root
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
S="$ROOT/scripts"
W="$ROOT/tests/work"
rm -rf "$W" && mkdir -p "$W" && cd "$W"

echo "== 1. test structures"
python "$ROOT/tests/make_test_inputs.py" structures .

echo "== 2. enumerate 20 random Ce orderings (6 of 32 Zr)"
python "$S/enumerate_substitutions.py" POSCAR --host Zr --dopant Ce --n 6 --random 20 --seed 1 | tail -3

echo "== 3. Warren-Cowley check: first line vs independent recalculation"
f=$(ls Ce_18.75/POSCAR_* | head -1)
head -1 "$f"
python "$S/wc_sro.py" "$f" --shell-edges 0 4.8 --species Zr Ce | tail -3

echo "== 3b. Warren-Cowley shell-1 regression test (shapes, sum rule, known orderings, m cell)"
python "$ROOT/tests/test_wc_shell1.py" | tail -1

echo "== 4. map onto t and c templates"
python "$S/batch_map.py" Ce_18.75 --templates POSCAR_t POSCAR_c --jobs 4 | tail -3
head -1 "$(ls Ce_18.75/mapped/t/POSCAR_t_* | head -1)"

echo "== 5. free energies"
if [ -n "${NEP:-}" ] && python -c "import calorine, phonopy" 2>/dev/null; then
  python "$S/batch_free_energy.py" Ce_18.75 --nep "$NEP" --T 300 --glob 'POSCAR_000[12]*' --jobs 2 --out real_test.csv
fi
python "$ROOT/tests/make_test_inputs.py" fake-csv Ce_18.75 300

echo "== 6. plots"
python "$S/plot_free_energy.py" Ce_18.75/free_energy_300K.csv --bootstrap 10 --formats png | grep -E "peak|wrote" | head -8
python "$S/plot_cell_annotations.py" Ce_18.75/free_energy_300K.csv --bootstrap 10 --formats png | head -6
python "$S/plot_cell_distributions.py" Ce_18.75/free_energy_300K.csv --bootstrap 10 --formats png | head -6

echo
echo "SMOKE TEST PASSED - outputs in $W (energies are placeholders)"
