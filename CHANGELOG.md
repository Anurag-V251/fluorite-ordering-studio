# Changelog

Every release is a tag (`vX.Y.Z`) on `main`, merged from a pull request with its commits kept.
To go back if something fails:

```bash
git fetch --tags
git checkout v1.2.0            # run an older release (detached); `git checkout main` to come back
# or undo one release on main for everyone (keeps the history):
git revert -m 1 <merge commit of that release>    # the merge commit is listed in `git log --merges`
```

Run folders keep their own frozen scripts, so going back to an older app version never changes a run that
is already prepared. `~/.fluorite_studio` (token, job records) works with every version.

## 1.3.0 - new look, app shell and run dashboard (GUI redesign, phase 1 of 3)

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
* `campaign.py CONFIG status --json` prints the status as one JSON object (same numbers as before).
* Tests: new checks in test_units.py (sbatch output parsing, squeue/sacct merge, cost model), test_studio.py
  (state = scan, measured time = `status --json`, fake squeue task → running chunk, re-plan before resubmit,
  workdir opens its run) and test_ui.js (state strip, next step, theme button).
* Revert: `git checkout v1.2.0`.
