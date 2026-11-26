# Fluorite Ordering Studio

A browser interface for the cation-ordering → phase mapping → NEP relaxation → harmonic free-energy workflow
(`campaign.py` and the scripts of *fluorite-ordering-thermo*). It runs on your machine or on the cluster login
node. It only uses the Python standard library and has no internet or CDN dependency. The scripts in
`scripts/` are the same files as before, unchanged; the app only calls them.

```bash
python3 start_studio.py
```

It prints a link such as `http://127.0.0.1:8765/?token=…`. On a cluster, start it on the login node and open an
SSH tunnel from your laptop. The command is printed too:

```bash
ssh -N -L 8765:127.0.0.1:8765 you@login-node        # then open the printed link in your browser
```

### Your files are on the cluster: `connect_cluster.sh` (run on your laptop)

Get the app on the cluster once (`git clone`, see *Getting and updating the code* below, or
`scp -r fluorite-ordering-studio you@cluster:~/`), then on your laptop:

```bash
./connect_cluster.sh -d /path/to/your/runs -s 'module load foss Python' you@login.cluster.org
```

It logs in once (password / 2FA asked once), starts the app on the login node it landed on (or reuses the one
already running there), forwards a free local port through the same connection and opens a `http://127.0.0.1`
link. Everything stays on the cluster, so Results, Plots and SLURM Submit all work on your real run folders.
Ctrl+C disconnects and stops the app; `-k` keeps it running, `--stop` stops it later. Jobs started from the
app keep running either way. Save your settings once in `~/.fluorite_studio/cluster.conf` and then just run
`./connect_cluster.sh`:

```bash
HOST=you@login.cluster.org
APP=~/fluorite-ordering-studio          # where start_studio.py is on the cluster
DIR=/path/to/your/runs                  # folder the app opens in
SETUP='module load foss Python'         # whatever gives python3 with numpy/scipy/matplotlib
```

Other options: `-P python` (Python on the cluster), `-p 8765` (first local port tried), `-o 'IdentityFile=…'`
(any ssh option), `-n` (no browser). Needs `ssh` with ControlMaster support (Linux, macOS, WSL on Windows).

Start it from the Python environment that has numpy / scipy / matplotlib (plots) and ase / calorine / phonopy
(free energies). Options: `--port`, `--dir START_FOLDER`, `--scripts DIR`, `--state DIR` (default
`~/.fluorite_studio`), `--new-token`.

**Security.** The server listens on 127.0.0.1 only. Every request needs the random token, which is kept in a
cookie after the first visit, so other users of a shared login node cannot use your server. The app runs
commands as you, and only the commands of registered modules.

### Getting and updating the code (git)

The code lives in a private GitHub repository. Clone it once on every machine (laptop, cluster login node):

```bash
git clone git@github.com:Anurag-V251/fluorite-ordering-studio.git     # SSH key on GitHub
# or: git clone https://github.com/Anurag-V251/fluorite-ordering-studio.git   (asks for a token)
cd fluorite-ordering-studio
```

Update to the latest version (stop a running app first, or restart it after the pull):

```bash
cd ~/fluorite-ordering-studio && git pull
git checkout v1.3.0        # optional: go to a tagged release; `git checkout main` to come back
```

Run folders keep their own frozen copy of the scripts (`<run>/scripts/`), so a `git pull` never changes a
run that is already prepared. Keep run folders, results, NEP models and `~/.fluorite_studio` outside the
repository (they are also listed in `.gitignore`). Personal modules go in `modules.d/` or `$FO_STUDIO_MODULES`.
`MD5SUMS.txt` describes the release it was made for; `git status` / `git diff` are the check for local edits.
What changed in each release, and how to go back to an older one: `CHANGELOG.md`.

---

## Pages

| Page | What it does |
|---|---|
| **Results** | Point it at a folder and it scans it. It reports compositions, x, orderings, whether each composition is prepared, structures per phase, relaxed files, F done, flagged (imaginary modes / unconverged relaxation), failed, not run, chunk progress, the lowest-F phase per composition, mapping statistics (rms shift, strain S, equally good variants, IDENTICAL count), the merge report and the logs. Temperature selector, auto-refresh. |
| **Plots** | Every plot script (F vs x HTML, c/a vs x HTML, F vs composition, F vs Ce hist, P(F), cell annotations, cell distributions). The option form is built from the script's own `argparse`, so new options appear by themselves. The result is shown in the page (HTML in a frame, PNG inline), with the equivalent command line. Options are remembered per plot. |
| **Pipeline** | The run dashboard. Header with a run switcher; a state strip with the counts (clean / flagged / failed / not run, chunks, core-hours used vs planned) and the **Next step**; the stages: dry-run → prepare → plan → compute → status → merge → rerun flagged, each run here, written as a SLURM script, or submitted. Opens the run you used last; *Runs* lists all runs and opens any run folder or campaign workdir (also runs made without the app). |
| **New run** | Host POSCAR (= starting phase), host/dopant species, compositions, target phase templates ("phases for transformation", same cation ordering), all-variants, NEP model, T, relaxation, phonon and cluster settings. Writes a reproducible run folder. |
| **Jobs** | Every command started from the app, with a live log and the files it wrote. Jobs keep running if you close the page or stop the server. Your SLURM queue (`squeue`), `sacct` of jobs submitted from the app, and `scancel`. |
| **XRD** | One composition at a time: P(F) per phase (same colours as plot_free_energy.py, Gaussian-mixture fits). Drag on it, type the window (meV/atom above the lowest unflagged F) or pick a Gaussian (μ ± kσ), choose the phases with the chips. Every relaxed structure in the window gets its powder pattern from pymatgen, and the page shows ONE averaged pattern (mean, ±1 std or min–max, single patterns), with the references c, t, oIII, m (12-atom cells as they are, stacked, overlaid or as Bragg ticks with hkl). Details below. |
| **Tools** | `map_check` (map_fluorite.py or transform_set.py on one structure, then compare_ordering.py on every output: IDENTICAL, rms shift, S, equally good variants), compare_ordering, enumerate, batch_map, wc_sro, transform_set, batch_free_energy. |
| **Modules** | What is registered, where each module comes from, and how to add or replace one. |

## XRD page (`scripts/xrd_average.py`)

* **Selection**: composition + phases + F window; flagged rows (imaginary modes, unconverged) are left out
  unless ticked; at most N structures (random subset, seed). The geometry is the relaxed file of the row that
  counts (`relaxed/<file>`, or `rerun_RRR/relaxed/<file>` when the rerun row won; `Ce_xx/relaxed/` in plain folders).
  Rows without a relaxed file are skipped and reported.
* **Recomputed (a job, sticks cached per file md5 + settings, so only new structures cost time)**: selection,
  radiation (pymatgen names or Å), 2θ range beyond what was computed (5° margin), species change, references,
  symprec, Debye–Waller B, stick cut-off.
* **Instant (in the browser)**: step, FWHM (or Caglioti U V W), Gaussian / Lorentzian / pseudo-Voigt (η),
  normalisation (each pattern max = 100, or absolute I/n_atoms², which is the same for any supercell of a cell),
  spread, sqrt axis, zoom (drag; double-click resets), hover readout (2θ, d, mean ± std, nearest reference hkl).
* **Species for the XRD**: every cation → Zr / Ce / any element, or pairs such as `Ce:Zr`; applied in memory
  only (the relaxed POSCARs are only read). Applied to the references too unless unticked. References can be
  scaled to the selection's volume/atom.
* **Output**: CSV / SVG / PNG of what is shown, the structure list, the equivalent command line, and
  "Write files here" (a job running xrd_average.py: CSV, PNG, PDF, `_structures.csv`, `_ref_peaks.csv`, JSON record).
* Needs `ase` and `pymatgen` where the app runs. 100 structures of 96 atoms take about 5 s on 4 cores.

```bash
python3 scripts/xrd_average.py RUN_FOLDER --composition Ce_25 --phase oIII --F-range 3.8 4.8 --ref c t oIII m
python3 scripts/xrd_average.py free_energy_300K_all.csv --composition 25 --phase m --F-range 0 2 --replace-cations Zr
```

## Folder layouts the scan understands

Derived from what `campaign.py` / `rerun_flagged.py` / `batch_free_energy.py` write:

```
<run>/campaign.toml                     (pointing at the run folder also works)
<workdir>/Ce_<pct>/summary.csv, POSCAR_<id>_WC<α>         orderings  (phase "parent")
<workdir>/Ce_<pct>/mapped/<phase>/…, mapped/summary.csv, mapped/verify.csv
<workdir>/Ce_<pct>/prepared.json
<workdir>/chunks/plan.json, chunk_*.txt
<workdir>/results/chunk_*.csv            (chunk_rerunRRR_K.csv from rerun_flagged.py)
<workdir>/relaxed/…, rerun_RRR/relaxed/…
<workdir>/free_energy_<T>K_all.csv, merge_report_<T>K.txt, logs/, slurm/
```

It also reads plain folders processed by hand (`Ce_xx/free_energy_<T>K.csv`, `Ce_xx/relaxed/`). Rows are
combined the way `campaign.py merge` combines them: the newest row per structure, and an ok row beats a FAILED one.

## Reproducible runs

*New run* writes a run folder:

```
<run>/campaign.toml      what campaign.py reads (paths relative to the folder)
<run>/run_record.json    settings, md5 of every input and script, package versions, user, host, and the
                         history of every command started from the app
<run>/inputs/            copies of the host POSCAR, the templates and the NEP model
<run>/scripts/           frozen copy of the scripts; stages of this run use these, not the app's
<run>/campaign/          workdir
```

*Check reproducibility* compares the md5 values now with the record. *Copy settings to a new run* starts a new
run from the old settings, for example with another material, NEP model or starting phase. A different starting
phase means a separate run, because the orderings are drawn on its host cell. The folder alone is enough to
re-run without the app: `python3 scripts/campaign.py campaign.toml dry-run`.

## Changing or adding parts (modules)

Every pipeline stage, plot and tool is one small TOML file in `studio/modules/`. To add a module or replace a
shipped one, put a file in `modules.d/` (or a folder listed in `$FO_STUDIO_MODULES`). A file with the same `id`
replaces the shipped module, and `enabled = false` hides it. Then click *Reload modules*. You do not need to change
any code.

```toml
# modules.d/my_plot.toml: a new plot; its argparse options become the form
id     = "my_plot"
kind   = "plot"                      # plot | stage | tool
title  = "My new plot"
script = "plot_my_thing.py"          # in scripts/, or an absolute path
input  = "merged_csv"                # merged_csv | comp_csv | any_csv | config | poscar | folder | none
watch  = ["input_dir", "--out-dir"]  # where its files appear

# modules.d/stage_40_compute.toml: replace a stage, e.g. by your own GPUMD driver
id = "compute"
kind = "stage"
title = "Relax + F (GPUMD)"
order = 40
[modes.local]
cmd = ["{python}", "/path/my_driver.py", "{config}"]
[modes.slurm]
cmd = ["{python}", "/path/my_driver.py", "{config}", "--submit"]
[[options]]
flag = "--unfinished"
label = "only unfinished"
value = false
```

Placeholders: `{python} {scripts} {config} {workdir} {rundir} {studio}`. Pages are separate files in
`studio/static/pages/`. A new `NN_name.js` that calls `Studio.registerPage({...})` adds a page.

```
start_studio.py              launcher
connect_cluster.sh           laptop -> cluster connection (ssh, port forward, start/stop)
studio/server.py             HTTP API (stdlib)        studio/scan.py      folder scanner
studio/jobs.py, runner.py    detached jobs, SLURM     studio/runs.py      run folders, campaign.toml, md5 record
studio/modules.py            module registry          studio/argspec.py   forms from argparse
studio/helpers/map_check.py  map/transform + verify   studio/helpers/run_chunks.py  local compute of all chunks
studio/modules/*.toml        shipped modules          studio/static/      page (core.js + pages/*.js, no libraries)
scripts/                     the workflow scripts (unchanged) + plot_ca_vs_x_html.py + transform_set.py
refs/POSCAR_{c,t,oIII,m}     12-atom references       tests/test_studio.py  end-to-end test of the app
```

## Tests

```bash
python3 tests/test_units.py     # ~15 s: create-time checks, stop/kill, SLURM ids, scan rows
python3 tests/test_studio.py    # ~5 min, no NEP / cluster: mock calorine (Einstein crystal), fake sbatch
node tests/test_ui.js           # ~3 min, after test_studio.py: clicks through the page in Chromium (Playwright)
bash tests/test_campaign.sh     # the campaign.py test, as before
python3 tests/test_xrd.py       # ~2 min: XRD script + endpoints (needs ase, pymatgen); mock campaigns from make_xrd_mock.py
node tests/test_xrd_ui.js       # ~1 min, after test_xrd.py: the XRD page in Chromium
```

`test_studio.py` drives the app through its HTTP API. It creates three runs: a cubic 1×1×2 host, a cubic
2×2×2 host and a tetragonal 2×2×2 host (starting phase t), with the real references. It then runs every stage
locally and checks that the scan agrees with the files. It runs all seven plot modules, and runs map_check
(map_fluorite and transform_set) on c, t, oIII and m supercells of shapes 2×2×2, 1×2×4 and 3×2×2; every
output must be IDENTICAL. It also checks SLURM writing and submission, the reproducibility check, and that
requests without the token are refused.

`test_ui.js` uses a real browser: it fills the New-run form (a wrong composition first, which must be refused),
clicks every stage, checks the Results page, runs a plot and map_check, checks dark mode, and fails on any page
error. It needs Node and Playwright (`npm i playwright`); it is optional.

## Notes

* Harmonic F at the fixed relaxed cell (not QHA). Rows with imaginary modes or unconverged relaxation are
  flagged everywhere and are not reliable.
* Equally good mapping variants cannot be ranked by distortion, only by energy. Tick *all variants* to compute
  them all.
* `scripts/map_fluorite.py` is the package version (md5 afb8ff11…). It differs from the validated one
  (7a6b3022…) only by the `--label` option and the title line, which `batch_map.py` needs. The geometry code
  is identical.

## Changes in 1.3.0

* New look and run dashboard (state strip, next step, light/dark). Details and how to revert: `CHANGELOG.md`.

## Changes in 1.2.1

* `campaign.py status` fixed for old and rerun campaigns (see `CHANGELOG.md`): flags now count only the newest
  row per structure (a structure fixed by `rerun_flagged.py` is no longer reported as flagged), counts match
  `merge`, measured `sec_per_structure` ignores rerun rows, clearer unfinished-chunk and FAILED lines,
  old CSVs with missing or odd columns no longer crash. New `status --json`.

## Changes in 1.2.0

* New **XRD** page and `scripts/xrd_average.py` (also on the Tools page): averaged powder XRD of the structures
  in an F window, species replacement in memory, reference phases. `tests/make_xrd_mock.py` makes a mock
  campaign (real enumerate/map/verify, PLACEHOLDER relaxation and F) to try it without NEP.

## Changes in 1.1.2

* F vs x plots (`plot_F_vs_composition.py`, `plot_F_vs_x_html.py`, Plots page): new y axis `--ref formation`,
  the formation free energy per atom F_form = F − [(1−x) F(x=0) + x F(x=1)] from the pure end members.
  `--ends same` (default: each phase relative to its own end members), `lowest` (lowest-F phase of the CSV at
  each end) or phase names (`--ends m,parent`). Missing or flagged end members stop with a clear message.
  Panel (b) ΔF is unchanged. Test: `python tests/test_formation_ref.py` (~30 s).

## Changes in 1.1.1

* `connect_cluster.sh`: one command on the laptop to use the app on the cluster files (see above). The app
  itself is unchanged.

## Changes in 1.1.0 (three rounds of user testing and code review)

* Create run checks every input first (host element, 12-atom templates, NEP file, compositions, numbers,
  walltimes, workdir) and leaves nothing behind when it fails.
* Plots take the dopant, the temperature and the starting phase from the run that holds the CSV
  (`[from_run]` in a module file), so another material gets the right axis labels and legends.
* Running stages come back after a page reload; failed stages show the next step; "done, N failed" when some
  structures failed; Stop in the Jobs list; a job that ignores Stop is killed after 10 s.
* The app runs only scripts that belong to you and that nobody else can change; frozen scripts are made go-w.
* Large campaigns: result CSVs are streamed (only the needed columns are kept); the Jobs page no longer walks
  output folders; big files are streamed; log views say when they show only the end.
* Results and Pipeline list the runs found below the start folder; host element set from the POSCAR;
  file picker can show hidden files and make folders; the Modules page lists disabled modules.
