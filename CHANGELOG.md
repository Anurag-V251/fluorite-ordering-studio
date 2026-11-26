# Changelog

Every release is a git tag (`v1.3.0`, `v1.2.0`, ...) on `main`. Each change arrives as a pull request whose
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

## 1.2.0 (2026-11-22)

XRD page and `scripts/xrd_average.py`. Earlier versions: see "Changes in ..." in README.md.
