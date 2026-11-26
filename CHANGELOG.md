# Changelog

Every release is a git tag (`v1.3.0`, `v1.2.1`, `v1.2.0`, ...) on `main`. Each change arrives as a pull request whose
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
