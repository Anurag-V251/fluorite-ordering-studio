#!/usr/bin/env bash
# Full workflow for one composition. Edit the variables, then run from the folder that holds
# POSCAR (host supercell), the phase templates and the NEP model.
set -euo pipefail

SCRIPTS=/path/to/fluorite-ordering-thermo/scripts   # <- edit
HOST=Zr; DOPANT=Ce; N_SUB=6                          # substitute 6 Zr by Ce
N_RANDOM=500                                         # random distinct orderings (omit --random for all)
TEMPLATES="POSCAR_t POSCAR_m POSCAR_oIII"            # 12-atom target-phase cells
NEP=nep.txt
T=300
JOBS=16                                              # e.g. $SLURM_CPUS_PER_TASK on a cluster

# 1. orderings + Warren-Cowley alpha  -> e.g. Ce_18.75/
python "$SCRIPTS/enumerate_substitutions.py" POSCAR --host $HOST --dopant $DOPANT --n $N_SUB \
       --random $N_RANDOM --seed 1
DIR=$(ls -d ${DOPANT}_* | head -1)                   # folder just created
echo "working in $DIR"

# 2. same orderings in the other phases  -> $DIR/mapped/<phase>/
python "$SCRIPTS/batch_map.py" "$DIR" --templates $TEMPLATES --jobs $JOBS

# 3. harmonic free energies of everything (restart-safe: just run again if interrupted)
python "$SCRIPTS/batch_free_energy.py" "$DIR" --nep $NEP --T $T --relax --fmax 1e-3 --jobs $JOBS

# 4. figures
CSV="$DIR/free_energy_${T}K.csv"
python "$SCRIPTS/plot_free_energy.py"        "$CSV"
python "$SCRIPTS/plot_cell_annotations.py"   "$CSV"
python "$SCRIPTS/plot_cell_distributions.py" "$CSV"
