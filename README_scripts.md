# fluorite-ordering-thermo

Cation-ordering → phase → free-energy workflow for fluorite-derived oxides (ZrO₂–CeO₂ and similar),
built around a NEP machine-learned potential.

For one composition it:

1. **enumerates** cation orderings of a host supercell (all of them, or N random distinct ones) and labels each
   with its **Warren–Cowley short-range-order parameter** α;
2. **maps** every ordering onto other fluorite-derived phases (t, m, oIII, …) while keeping the cation ordering
   exactly the same;
3. computes the **harmonic free energy** F(T) = E_static + F_vib(T) of every structure with NEP + calorine +
   phonopy, in parallel, into one CSV;
4. makes **publication-ready figures**: P(F) per phase with multi-peak fits and ΔF between peaks, α vs F,
   P(F) peaks annotated with c/a and cell angle, and the distributions P(c/a) and P(θ).

---

## Pipeline

```
POSCAR (host, e.g. ZrO2 2x2x2)
   │  enumerate_substitutions.py            → Ce_18.75/POSCAR_0001_WC-0.1282 …  (+ summary.csv)
   ▼
Ce_xx/   ── batch_map.py (map_fluorite.py)  → Ce_xx/mapped/t/POSCAR_t_0001_WC-0.1282 …  (one folder per phase)
   │
   │  batch_free_energy.py --relax --jobs N → Ce_xx/free_energy_300K.csv  (+ Ce_xx/relaxed/…)
   ▼
free_energy_300K.csv
   ├─ plot_free_energy.py         → P(F) with Gaussian-mixture peaks, ΔF table; α vs F
   ├─ plot_cell_annotations.py    → P(F) peaks labelled with c/a and θ; per-peak cell table
   └─ plot_cell_distributions.py  → P(c/a) and P(θ) per phase
```

Every POSCAR carries its labels in the **first line**, and they survive mapping and relaxation:

```
Zr26Ce6O64 | Ce 6/32 = 18.75% | WC_shell1(Ce-Zr)=-0.1282 | id=0001 | phase=t | mapped from … | relaxed with nep.txt
```

---

## Installation

```bash
pip install -r requirements.txt          # numpy scipy matplotlib ase calorine phonopy
bash tests/smoke_test.sh                 # whole chain on a synthetic cell, ~20 s, no NEP model needed
```

Keep **all scripts in one folder** (`scripts/`): `batch_map.py` calls `map_fluorite.py`, and the two
`plot_cell_*` scripts import `plot_free_energy.py`. Mixing versions from different folders causes import errors.

---

## Step by step

### 1. Enumerate orderings + Warren–Cowley α

```bash
python scripts/enumerate_substitutions.py POSCAR                                   # interactive
python scripts/enumerate_substitutions.py POSCAR --host Zr --dopant Ce --n 6 --random 500 --seed 1
python scripts/enumerate_substitutions.py POSCAR --host Zr --dopant Ce --n 4      # all C(n_host, n)
```

* Output folder `Ce_<% of host sites>` (e.g. `Ce_18.75`), files `POSCAR_<id>_WC<α>`, `summary.csv`.
* `--random N` draws N new distinct orderings; re-running adds N more and never repeats (reads `summary.csv`).
* `--shells 2` also writes the 2nd-shell α; `--unique` (all combinations, needs pymatgen) keeps only
  symmetry-inequivalent ones.

Independent check of α for any structure: `python scripts/wc_sro.py POSCAR --shell-edges 0 4.8 --species Zr Ce`
(`scripts/icet_sro.py` does the same with icet).

### 2. Map onto other phases

```bash
python scripts/batch_map.py Ce_18.75 --templates POSCAR_t POSCAR_m POSCAR_oIII --jobs 8
```

* Templates: 12-atom cells (4 cations + 8 O) of the target phases, axes along the cubic fluorite directions
  (pseudo-cubic setting for t). The phase label is the template name without `POSCAR_`.
* Output `Ce_18.75/mapped/<phase>/POSCAR_<phase>_<id>_WC<α>` and `mapped/summary.csv`.
* `--jobs N` workers, each with `--threads` BLAS threads (default 1); `--progress 30` prints a line with count and
  ETA every 30 s; `summary.csv` is saved every minute, so `--skip-existing` restarts a stopped run where it stopped.
* The cation on every lattice site is unchanged, so α carries over. Several equally good orientations →
  variant 1 is written (`--args "--all-variants"` writes all).

### 3. Free energies (parallel)

```bash
python scripts/batch_free_energy.py Ce_18.75 --nep nep.txt --T 300 --relax --fmax 1e-3 --jobs 16
```

* All `POSCAR*` in the folder **and sub-folders** (parents + mapped phases); `relaxed/` is skipped.
* One CSV row per structure: E_static, F_vib, F_total (eV/atom), S, Cv, relaxed cell (a, b, c, α, β, γ),
  phonon supercell, lowest frequency, number of imaginary modes, phase, WC, status.
* Restart-safe: rows are written as each structure finishes; re-running skips finished ones and retries failed ones.
* `--jobs N` processes in parallel, each limited to `--threads` (default 1) BLAS/OpenMP threads; results are
  identical to a serial run. Phonon supercell: smallest repeat with every width ≥ `--min-length` (10 Å) unless
  `--supercell` is given. `--fixed-cell` relaxes positions only.
* `scripts/thermal_harmonic_nep.py` is the original single-structure, full-F(T) script this is based on.

### 4. Figures

```bash
python scripts/plot_free_energy.py         Ce_18.75/free_energy_300K.csv
python scripts/plot_F_vs_composition.py    campaign/free_energy_300K_all.csv      # whole composition range
python scripts/plot_F_vs_x_html.py         campaign/free_energy_300K_all.csv      # one HTML file: same plot with live controls, no server
python scripts/plot_F_vs_composition.py    all.csv --ref formation               # formation free energy per atom, F - [(1-x)F(x=0) + x F(x=1)]
                                           # --ends same (own phase, default) | lowest | PHASE_A,PHASE_B | fixed F, e.g. --ends=same,-8.2034
                                           # --end-x A B: other end compositions in % (default 0 100), e.g. --end-x 0 50
python scripts/plot_cell_annotations.py    Ce_18.75/free_energy_300K.csv
python scripts/plot_cell_distributions.py  Ce_18.75/free_energy_300K.csv
```

| Script | Figures (PDF + 600 dpi PNG) | Tables |
|---|---|---|
| `plot_free_energy.py` | `_PF` (P(F) + fits), `_WC_vs_F` (F on x, α on y), `_combined` | `_peak_fit.csv` |
| `plot_cell_annotations.py` | `_PF_ca`, `_PF_angle`, `_PF_cell_combined` | `_peak_cell.csv`, `_cell_params.csv` |
| `plot_cell_distributions.py` | `_P_ca`, `_P_angle` | `_cell_peaks.csv` |

Common options: `--reference main`, `--components m=2`, `--exclude-imaginary`, `--label parent="c (parent)"`,
`--font sans`, `--formats pdf svg`, `--legend-loc "upper right"`, `--bootstrap 30` (faster).
Each script prints the sentences to put in the figure caption.

---

## Whole composition range on a cluster: `campaign.py`

One config file (`examples/campaign.toml`) drives everything for many compositions at once, e.g. Ce = 0 … 32 in a
96-atom ZrO₂ cell, 1000 random orderings each, mapped onto t, oIII and m, harmonic F at 300 K:

```bash
python scripts/campaign.py campaign.toml dry-run                 # counts per composition/phase + cost estimate
python scripts/campaign.py campaign.toml slurm prepare --submit  # array: one task per composition
python scripts/campaign.py campaign.toml plan                    # split all structures into one-node chunks
python scripts/campaign.py campaign.toml slurm compute --submit --with-merge   # array of chunks + merge job
python scripts/campaign.py campaign.toml status                  # progress, failures, flags, measured s/structure (--json)
python scripts/rerun_flagged.py campaign.toml --submit          # afterwards: recompute flagged rows (stricter settings) + merge
```

Run `plan` once all prepare tasks have finished (`dry-run` then shows every composition as `prepared`).
(Without SLURM: `prepare`, `plan`, `run --chunk K`, `merge` do the same steps locally.)

* **Orderings**: for each dopant count n, `enumerate_substitutions.py --random MAX --seed SEED*1000+n` — MAX distinct
  site sets drawn uniformly at random; **all** C(N, n) when fewer exist (n = 1 → 32, n = 2 → 496 for 32 sites).
  n = 0 and n = N give the single end member (α undefined, written `nan`). No symmetry reduction: in a cubic parent
  some draws are symmetry-equivalent (many at n ≤ 7 and n ≥ 25), which costs compute but keeps the statistics of a
  uniform random sample. The seed is stored in `<comp>/prepared.json`.
* **Mapping**: `batch_map.py --skip-existing` onto every template; variant 1 of equally good orientations
  (`map_args = "--all-variants"` computes all of them: 6 (t) to 48 (oIII, m) per ordering for a cubic parent —
  they can only be ranked by energy). Every mapped file is checked with `compare_ordering.py` against its parent
  (`<comp>/mapped/verify.csv`); `prepare` stops unless all are IDENTICAL. Mapping statistics (rms shift, strain S,
  number of equally good variants) are printed per composition and phase.
* **Chunks**: structures are packed into chunks that fill `fill` × `walltime` on one node with
  `cores_per_node / threads` structures in parallel, using
  seconds = `sec_per_structure` × (n_atoms × phonon-supercell atoms) / N₀². `sec_per_structure` is a guess until
  measured: run one chunk, read `status`, put the value in the config, then `plan --force` (no compute job running) to re-pack the unfinished structures; finished ones are kept and `merge` includes them.
* **Status**: counts use the merge rule (newest row per structure, ok beats FAILED, rerun rows replace campaign
  rows), so `status`, `merge` and the Results page agree; results from an earlier plan are included. Unfinished
  chunks are split into *not started*, *partly run* and *only FAILED left* (the most common FAILED notes are
  printed). Flags count only the newest row of each structure (fixed by a rerun = not flagged). The measured
  `sec_per_structure` uses campaign rows only (not rerun rows) and says when to re-plan. `status --json` gives
  the same as JSON.
* **Restart**: every step skips finished work. Re-submit unfinished chunks with
  `slurm compute --unfinished --submit --with-merge`. Raising `max_orderings` later: `slurm prepare --unfinished`,
  `plan --append`, `slurm compute --unfinished`. Arrays larger than `max_array` (MaxArraySize) are split.
* **Output**: `free_energy_<T>K_all.csv` (all compositions, with `n_dopant`, `x_dopant`), one
  `<comp>/free_energy_<T>K.csv` per composition for the plot scripts, `merge_report_<T>K.txt` (ok / failed /
  imaginary / unconverged per composition and phase, lowest ΔF with and without flagged rows). Relaxed structures in
  `relaxed/<comp>/…`.
* `batch_free_energy.py` now records `max_force_eV_A` and `relax_converged` (False if the force is still above
  `--fmax` after `--steps`); such rows are kept and flagged in the note.

Test without NEP or a cluster: `bash tests/test_campaign.sh` (≈ 3 min; mock calorine = Einstein crystal whose F_vib is
checked against the analytic value; fake `sbatch`; cell shapes 2×2×2, 1×2×4, 3×2×2, 1×1×2 and a tetragonal parent).

---

## Definitions and conventions

**Warren–Cowley α** (same as sqsgenerator / icet), on the substituted cation sublattice:
α⁽ⁱ⁾_AB = 1 − N⁽ⁱ⁾_AB / (N · M⁽ⁱ⁾ · x_A · x_B). α < 0: A–B neighbours preferred (mixing/ordering);
α = 0: random; α > 0: avoided (clustering). For a binary sublattice one α per shell holds all the information.
α is lattice-based: it describes which site holds which cation, so it is the same before and after mapping and
relaxation.

**Phases.** `phase=parent` = the un-mapped enumerated structures (plotted as "c (initial)"); other labels come
from the template names.

**Free-energy peaks.** Per phase, a mixture of 1–3 Gaussians is fitted to the F values (maximum likelihood; an
extra Gaussian only if the BIC drops by > 10). Peaks = maxima of the fitted curve, kept only if they reappear in
≥ 70 % of 100 bootstrap re-fits (which also give the uncertainties). ΔF_peak is measured from the lowest-F peak
of all phases (`--reference main`: the lowest of the phases' tallest peaks). Values like +2.50(3) mean
2.50 ± 0.03 meV/atom.

**Histogram normalisation.** P(bin) = count of that phase in the bin / total number of structures of that
composition (all phases together sum to 1); "probability per bin", not a density. Peak heights therefore scale
with the number of structures per phase; peak positions do not depend on the normalisation.

**Normalised lattice parameters.** a_n, b_n, c_n = supercell lengths / repeats (repeats from the number of cations,
4 per conventional cell, and fluorite-like lengths 4.0–6.5 Å; `--repeats` overrides). c = the longest normalised
axis, decided per structure (x, y or z), a_n = mean of the other two, c/a = c_n/a_n. This is the fluorite
(pseudo-cubic) setting, i.e. the experimental c/(√2·a_t) of the primitive P4₂/nmc cell — do not multiply by √2
again. θ = cell angle farthest from 90°, folded to ≥ 90° (monoclinic β). Labels like 1.026(2) in the annotation
figures are the mean with the **standard deviation across structures** in brackets.

---

## Caveats

* **Cell size.** α in a shell is only a property of the structure (not of the periodic cell) if the shell radius
  is below half the smallest cell width; the scripts warn otherwise. A 2×2×2 conventional cell resolves the 1st
  cation shell only; 3×3×3 resolves two.
* **Imaginary modes** make the harmonic F unreliable; they are flagged (hollow symbols, `--exclude-imaginary`).
* **Relaxation can change the phase.** Check the relaxed c/a and θ (`_cell_params.csv`) — a mapped m structure
  that relaxes to θ ≈ 90° is no longer monoclinic. `--fixed-cell` keeps the template cell (c/a and θ are then the
  template values).
* **Harmonic, fixed volume** — no thermal expansion (not quasi-harmonic). Compare trends with experiment, not
  absolute lattice parameters.
* **Accuracy.** ΔF of a few tenths of meV/atom can be statistically resolved with many structures but is below
  the accuracy of a NEP potential; compare against undoped references computed with the same potential.
* **c/a ≈ 1 does not prove cubic:** the t″ form has c/a = 1 with the oxygen sublattice still displaced along c.
* `GPUNEP` (instead of `CPUNEP`) writes files in a working directory and is not safe with `--jobs > 1`.

---

## Files

```
scripts/
  enumerate_substitutions.py   orderings + α, random or exhaustive, no repeats across runs
  wc_sro.py                    α for any POSCAR (numpy only; independent check)
  icet_sro.py                  α with icet's ShortRangeOrderObserver (optional)
  map_fluorite.py              map one structure onto a target phase, keeping the cation ordering
  batch_map.py                 map a folder onto several phases (parallel)
  batch_free_energy.py         harmonic F at one temperature for a folder tree (parallel, restartable)
  campaign.py                  many compositions: enumerate -> map -> verify -> chunk -> SLURM arrays -> merge
  rerun_flagged.py             recompute flagged rows (imaginary modes / unconverged) of a campaign in rounds, SLURM or local
  compare_ordering.py          same cation ordering in two structures? (IDENTICAL / EQUIVALENT / DIFFERENT)
  thermal_harmonic_nep.py      original single-structure F(T) script (reference)
  single_runs.py               a few POSCARs: relax, harmonic F over a T grid (F, U_vib, ZPE, S, Cv, DOS), XRD, plots
  plot_harmonic_compare.py     F(T), dF(T) vs a chosen reference, S/Cv panels, dF at T, DOS; PNG/PDF/SVG + HTML
  plot_free_energy.py          P(F) with multi-peak fits + ΔF; α vs F
  plot_F_vs_composition.py     all compositions: P(F) histograms per phase, multi-peak fits, peak branches; ΔF of every peak
  plot_F_vs_x_html.py          one self-contained HTML page of F vs composition with live controls (no server)
  plot_cell_annotations.py     P(F) peaks annotated with c/a and θ
  plot_cell_distributions.py   P(c/a) and P(θ)
tests/
  smoke_test.sh                end-to-end run on a synthetic cell (placeholder energies)
  make_test_inputs.py          synthetic POSCAR / templates / CSV for the smoke test
  test_campaign.sh             campaign.py end to end (shapes, restart, SLURM, flags), mock NEP
  test_batch_map.sh            batch_map.py: four phases x four cell shapes, same output as --subprocess, restart
  mock_calorine/               TEST-ONLY stand-in for calorine (Einstein crystal)
examples/
  run_pipeline.sh              the full workflow for one composition
  campaign.toml                config for campaign.py (all compositions, t/oIII/m, 300 K)
requirements.txt
```
