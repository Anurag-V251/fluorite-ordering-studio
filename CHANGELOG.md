# Changelog

Every release is a git tag (`v1.2.1`, ...). To go back to an earlier release:

```bash
git fetch --tags
git checkout v1.2.0            # look at / run the old version (detached); `git checkout main` to come back
git revert <commit>            # or undo one change on main, keeping history (commit ids below / `git log`)
```

Run folders keep their own frozen copy of the scripts (`<run>/scripts/`), so changing the version of the
repository never changes a run that is already prepared. To use a new `campaign.py` in an old run folder,
copy it into `<run>/scripts/` (the old copy can be put back the same way).

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
