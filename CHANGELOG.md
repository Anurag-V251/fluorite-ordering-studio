# Changelog

Every release is a git tag (`v1.5.0`, `v1.4.0`, `v1.3.0`, `v1.2.1`, `v1.2.0`, ...) on `main`. Each change arrives as a pull request whose
commits are kept, so `git log` shows what changed and when.

## How to go back to an earlier version

```bash
cd ~/fluorite-ordering-studio
git fetch --tags
git checkout v1.2.0          # run the old version (read-only "detached" checkout); restart the app afterwards
git checkout main            # back to the newest version
```

To undo a release on `main` itself (keeps the history, adds a commit that reverses it):

```bash
git log --oneline            # find the merge commit of the release, e.g. "Merge pull request #2 ..."
git revert -m 1 <merge-commit>
git push
```

Run folders keep their own frozen copy of the scripts, so going back never changes a finished run.
To use a new `campaign.py` in an old run folder, copy it into `<run>/scripts/` (the old copy can be put back the
same way). Single runs freeze their scripts the same way (`<single run>/scripts/`).

## 1.5.0 (2026-11-26) - run dashboard: pipeline, side panel, attention, jobs, settings (GUI redesign, phase 2 of 3)

* **Pipeline**: the 7 stages as one row with a live state each (pill with icon + word, short facts, chunk bar on
  Compute). Clicking a stage opens a **side panel**: facts for that stage, where to run it (here / submit to SLURM /
  write the script only), its options, the **equivalent command** (updates while you change options) and, after
  Run, its output. Escape or the scrim closes it. The old stage cards are gone.
* **Plan** has three explicit choices instead of check boxes: *First plan*, *Add new structures* (`--append`) and
  *Re-pack unfinished* (`--force`). Once a plan exists the default is *Add new structures*. Re-pack asks first and
  is **refused while a compute job of this run runs** (here or on SLURM), in the page and in the server.
* **Compute on SLURM** submits only unfinished chunks by default (`--unfinished` ticked; also the server default
  when a caller leaves it out).
* **Rerun flagged** options on the page: which flags (`--what`), imaginary window (`--freq-window`), `--same`,
  `--fmax`, `--steps`, `--start`, and for SLURM `--nodes` and `--walltime`.
* **Confirmations** before every SLURM submission, re-pack, Stop of a job and `scancel` (also on the Jobs page).
* **Attention** list: recovery actions in the order to do them (re-plan with the measured time before resubmitting,
  rerun flagged after compute, inspect failures, merge last), one button each, blocked ones say why.
* **Cluster settings** after creation: edit `[cluster]` (walltime, cores, sec_per_structure, partition, setup ...)
  from the page, or apply the measured sec_per_structure with one button. The old file is kept as
  `campaign.<time>.bak`, campaign.toml gets a comment line, run_record.json gets the new values, the new md5 and a
  history entry with every old → new value. Physics settings still need *Copy to new run*.
  New endpoints `POST /api/runs/settings`, `POST /api/runs/adopt`.
* **Runs made without the app**: *Keep a history* writes run_record.json next to their config (nothing else
  changes), so settings changes and stages are recorded from then on.
* **Jobs** panel of the run: the app's jobs (a SLURM submission shows "submitted N tasks · M active" with the live
  squeue state, not "done") and every SLURM task the app submitted: chunk, structures done, walltime used with the
  projected end (red when it will not fit), node; `scancel` per submission.
* **Log drawer** at the bottom: the newest job of the run (or one you pick), level filter All / Info / Warn / Error,
  remembers open / closed.
* Module files: new option keys `default`, `help`, `choices`, `kind = "choice"`, `confirm` (see `studio/modules.py`).
  Old module files keep working.
* Tests: test_studio.py (SLURM compute defaults to --unfinished and asks first, re-pack refused while compute runs,
  bad choices refused, rerun options reach rerun_flagged.py, settings change + backup + history + md5, physics and
  bad walltime refused, adopting a run), test_ui.js (every stage from its side panel, re-pack confirmation,
  rerun options, log drawer, editing settings, applying the measured time).
* Revert: `git checkout v1.4.0` (or `git revert -m 1 <merge commit of the 1.5.0 PR>`).

## 1.4.0 (2026-11-26) - new look, app shell and run dashboard (GUI redesign, phase 1 of 3)

* New design system (the chosen "mission control" design): one set of colour tokens, light and dark theme
  with a button in the top bar (remembered; first visit follows the system), visible focus rings, tabular
  numbers, status pills with an icon and a word (never colour alone). Every existing page uses it.
* App shell: top bar with icon navigation (Pipeline, New run, Jobs with a running-jobs badge, Results, Plots,
  XRD, Tools, Modules) and where the app runs (host, SLURM found or not, version).
* Run dashboard (Pipeline → a run): run header with a run switcher, the material, NEP model, cluster settings
  and age; a **state strip** with the real counts (structures, clean, flagged, failed, not run, chunks done,
  core-hours used / planned) as KPIs and a stacked bar, and a **Next step** box that says what to do now.
  Pipeline opens the run you used last; *Runs* lists all of them.
* New endpoint `GET /api/run/state?config=…` (studio/state.py): counts, the 7 stages, the SLURM tasks the app
  submitted (read from the sbatch output, live state from squeue/sacct, chunk of each array task), measured
  sec_per_structure, and the recovery actions in the recommended order. It also accepts a run folder or a
  campaign **workdir** (runs made without the app: the config next to the workdir is found).
* Tests: new checks in test_units.py (sbatch output parsing, squeue/sacct merge, cost model), test_studio.py
  (state = scan, measured time = `status --json` of 1.2.1, fake squeue task → running chunk, re-plan before resubmit,
  workdir opens its run) and test_ui.js (state strip, next step, theme button).
* Revert: `git checkout v1.3.0` (or `git revert -m 1 <merge commit of the 1.4.0 PR>`).

## 1.3.0 (2026-11-26)

New **Single runs** page: a few POSCARs + one NEP model, run here or on SLURM.

* `scripts/single_runs.py` (new): NEP relaxation (same settings and convergence check as
  `batch_free_energy.py`), harmonic free energy at one temperature or over a temperature ramp from ONE phonon
  calculation (F_total, F_vib, U_vib including the zero-point energy, ZPE, S, Cv per atom, phonon DOS), powder
  XRD of the structure used with optional references (`xrd_average.py`), comparison plots. Fixed-cell
  harmonic approximation, not quasi-harmonic. Flags imaginary modes, unconverged relaxations and unrelaxed
  inputs. Restartable. Species are read from the POSCARs (nothing assumes Zr/Ce).
* `scripts/plot_harmonic_compare.py` (new): F(T) and ΔF(T) against a chosen reference (a structure, or the
  lowest unflagged F at each T), F_vib / U_vib / S / Cv panels, ΔF at chosen temperatures, phonon DOS,
  comparison CSV, self-contained interactive HTML; PNG / PDF / SVG / EPS for publication.
* Studio: `studio/single.py` (endpoints `/api/single/*`), page `static/pages/22_single.js`, modules
  `plot_harmonic_compare.toml` (Plots page) and `tool_single_runs.toml` (Tools page). `server.py` only gained
  the import and the mixin class. Run folders: `request.json`, `inputs/`, frozen `scripts/`, `run.sh`,
  `job.sbatch`, `results/`.
* Tests: `tests/test_single_runs.py` (analytic Einstein crystal of the mock calorine, c/t/oIII/m x 1x1x1,
  2x2x2, 1x2x4, 3x2x2, flags, restart, plots, XRD only for HfO2, every endpoint incl. SLURM),
  `tests/test_single_ui.js` (the page in Chromium). The existing suites pass unchanged.
* No existing script changed.

## 1.2.1 (2026-11-26)

`campaign.py status` fixed for old and rerun campaigns (`scripts/campaign.py`, md5 in `MD5SUMS.txt`):

* Flags (imaginary modes, unconverged relaxation) now count only the newest row of each structure, as merge
  does. Before, structures fixed by `rerun_flagged.py` were still reported as flagged, and flagged reruns
  were counted twice.
* Measured `sec_per_structure` uses campaign rows only, one per structure (rerun rows use 10x steps and a
  larger phonon cell and biased it), and says when to put it into the config and `plan --force`.
* Counts now match `merge` and the Results page: one headline over all structures (current plan, earlier plan,
  FAILED, not run, prepared but not planned yet).
* Unfinished chunks listed as "not started / partly run / only FAILED left" (the old line
  `unfinished chunks: 0` read like "none"); the most common FAILED notes are printed.
* No crash on old CSVs with missing or odd columns (`n_imaginary_modes = nan`, empty `time_s`, no
  `relax_converged`), on structure files that were moved, or on missing chunk lists; warns when results at
  several temperatures are mixed.
* New `status --json` (same content, for scripts and the Studio).
* Tests: `tests/test_campaign.sh` section 7 (old messy run: `plan --force`, FAILED rows, rerun round,
  unplanned structures, old-style CSV); mock NEP can fail on purpose (`MOCK_NEP_FAIL=<Ce count>`).

Revert just this change: `git revert <merge commit of the 1.2.1 PR>`, or `git checkout v1.2.0 -- scripts/campaign.py`.

## 1.2.0 (2026-11-26)

First version in git: XRD page and `scripts/xrd_average.py` (see README "Changes in 1.2.0" and earlier).
