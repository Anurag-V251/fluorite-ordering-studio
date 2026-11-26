/* Single runs page: a few POSCARs + one NEP model -> relaxation, harmonic free energy at one T or over a
   T ramp (fixed-cell harmonic, not quasi-harmonic), XRD with optional references, comparison plots.
   Every submission is a run folder (inputs, frozen scripts, run.sh, job.sbatch, results/); the work is done
   by scripts/single_runs.py and scripts/plot_harmonic_compare.py (backend: studio/single.py). */
(function () {
  const S = Studio, h = S.h;
  const DEF = {name: "", base: "", poscars: [], names: [], nep: "", tasks: ["relax", "harmonic", "xrd", "plot"],
    Tmode: "ramp", T: "300", Tmin: 0, Tmax: 1500, Tstep: 10, Textra: "300", fmax: 1e-4, steps: 2000, fixed_cell: false,
    scMode: "auto", min_length: 10, sc: "2 2 2", mesh: 100, threads: "", xrd_ref: [], refs_dir: "", xrd_args: "",
    reference: "lowest", at_T: "300", formats: ["png", "pdf", "svg"], mode: "local",
    slurm: {partition: "", time: "04:00:00", cpus: 8, mem: "", account: "", qos: "", extra: "", setup: "", python: ""}};
  const TASKS = [["relax", "Relax (NEP)"], ["harmonic", "Harmonic free energy"], ["xrd", "XRD"], ["plot", "Plots"]];

  S.registerPage({id: "single", title: "Single runs", order: 22, render});

  function render(main, params) {
    let st = Object.assign({}, DEF, S.store.get("single.state", {}));
    st.slurm = Object.assign({}, DEF.slurm, st.slurm || {});
    const save = () => S.store.set("single.state", st);
    let D = null;
    const newCard = h("div", {class: "card"});
    const runsCard = h("div", {class: "card"});
    const resCard = h("div", {});
    main.append(newCard, runsCard, resCard);

    S.api("/api/single/defaults").then((r) => {
      if (!r.ok) return;
      D = r;
      if (!st.base) st.base = r.base;
      if (!st.refs_dir) st.refs_dir = r.refs_dir;
      drawNew(); drawRuns();
      if (params.folder) showResult(params.folder);
    });

    // ------------------------------------------------------------------ new run
    function drawNew() {
      newCard.innerHTML = "";
      const missing = Object.entries(D.needs).filter(([, ok]) => !ok).map(([m]) => m);
      newCard.append(h("h2", {}, "New single run"),
        h("div", {class: "small muted"}, "A few structures and one NEP model: relaxation, harmonic free energy at one temperature or a ramp " +
          "(one phonon calculation serves every T), XRD with references, and comparison plots. ",
          h("b", {}, "Fixed-cell harmonic approximation, not quasi-harmonic: no thermal expansion.")),
        missing.length ? h("div", {class: "notice warn small"}, "Not installed in the app's Python: " + missing.join(", ") +
          ". Local runs that need them will fail; SLURM runs use the Python set below.") : null);

      // structures
      const list = h("div", {});
      const drawList = () => {
        list.innerHTML = "";
        if (!st.poscars.length) list.append(h("div", {class: "muted small"}, "No structures yet."));
        st.poscars.forEach((p, i) => {
          const nm = h("input", {type: "text", value: st.names[i] || "", placeholder: S.basename(p), style: {width: "170px"},
            title: "label in tables and plots (default: file name)", onchange: () => { st.names[i] = nm.value.trim(); save(); drawRefSel(); }});
          list.append(h("div", {class: "row", style: {marginBottom: "4px"}}, nm, h("span", {class: "mono small grow"}, p),
            h("button", {title: "remove", onclick: () => { st.poscars.splice(i, 1); st.names.splice(i, 1); save(); drawList(); drawRefSel(); }}, "✕")));
        });
      };
      const addPick = async () => {
        const p = await S.pick({title: "Choose a POSCAR (upload from your computer with 'Upload here…')"});
        if (p) { st.poscars.push(p); st.names.push(""); save(); drawList(); drawRefSel(); }
      };
      const up = h("input", {type: "file", multiple: true, style: {display: "none"}, onchange: async (e) => {
        const dir = await uploadDir();
        if (!dir) return;
        for (const f of e.target.files) {
          const r = await S.api("/api/upload?overwrite=1&dir=" + encodeURIComponent(dir) + "&name=" + encodeURIComponent(f.name), f);
          if (r.ok) { st.poscars.push(r.path); st.names.push(""); }
        }
        up.value = ""; save(); drawList(); drawRefSel();
      }});
      const addFolder = async () => {
        const d = await S.pick({dir: true, title: "Folder: every POSCAR* / *.vasp file in it is added"});
        if (!d) return;
        const r = await S.api("/api/ls?path=" + encodeURIComponent(d));
        if (!r.ok) return;
        const fs = r.entries.filter((e) => !e.dir && /^(POSCAR|CONTCAR)|\.vasp$/i.test(e.name));
        if (!fs.length) return S.toast("no POSCAR* / CONTCAR* / *.vasp files in " + d);
        for (const e of fs) { st.poscars.push(r.path + "/" + e.name); st.names.push(""); }
        save(); drawList(); drawRefSel();
      };
      drawList();

      // NEP
      const nep = S.pathField(st.nep, {title: "NEP model (nep.txt)", placeholder: "nep.txt"});
      nep.input.addEventListener("change", () => { st.nep = nep.input.value.trim(); save(); });
      const nepUp = h("input", {type: "file", style: {display: "none"}, onchange: async (e) => {
        const f = e.target.files[0]; if (!f) return;
        const dir = await uploadDir(); if (!dir) return;
        const r = await S.api("/api/upload?overwrite=1&dir=" + encodeURIComponent(dir) + "&name=" + encodeURIComponent(f.name), f);
        if (r.ok) { st.nep = r.path; nep.input.value = r.path; save(); S.toast("uploaded " + f.name); }
        nepUp.value = "";
      }});

      // tasks
      const tasks = h("div", {class: "row"}, TASKS.map(([k, lab]) => {
        const c = h("input", {type: "checkbox", onchange: () => {
          st.tasks = TASKS.map(([x]) => x).filter((x) => x === k ? c.checked : st.tasks.includes(x)); save(); vis();
        }});
        c.checked = st.tasks.includes(k);
        return h("label", {}, c, " " + lab);
      }));

      // settings
      const num = (key, w = "90px", obj = st) => {
        const i = h("input", {type: "text", value: obj[key], style: {width: w}, onchange: () => { obj[key] = i.value.trim(); save(); }});
        return i;
      };
      const sel = (key, opts, obj = st, on) => {
        const s = h("select", {onchange: () => { obj[key] = s.value; save(); if (on) on(); }}, opts.map(([v, l]) => h("option", {value: v}, l)));
        s.value = obj[key];
        return s;
      };
      const Tbox = h("div", {});
      const drawT = () => {
        Tbox.innerHTML = "";
        Tbox.append(h("div", {class: "row"}, sel("Tmode", [["ramp", "temperature ramp"], ["const", "constant temperature(s)"]], st, drawT),
          st.Tmode === "ramp" ? [h("label", {}, "from ", num("Tmin", "70px")), h("label", {}, "to ", num("Tmax", "70px")),
            h("label", {}, "step ", num("Tstep", "60px")), "K", h("label", {class: "small"}, "also at ", num("Textra", "110px"), " K")]
            : [h("label", {}, "T ", num("T", "160px"), " K"), h("span", {class: "small muted"}, "several: 300 600 900")]));
      };
      drawT();
      const scBox = h("span", {});
      const drawSc = () => { scBox.innerHTML = ""; scBox.append(st.scMode === "auto" ? h("label", {}, "every width ≥ ", num("min_length", "60px"), " Å") : h("label", {}, "NX NY NZ ", num("sc", "90px"))); };
      drawSc();
      const refsBox = h("div", {class: "row"});
      const drawRefs = async () => {
        refsBox.innerHTML = "";
        const r = await S.api("/api/single/refs?dir=" + encodeURIComponent(st.refs_dir), undefined, {quiet: true});
        const names = r.ok ? r.refs : [];
        if (!names.length) refsBox.append(h("span", {class: "muted small"}, "no POSCAR_<name> files in this folder"));
        for (const n of names) {
          const c = h("input", {type: "checkbox", onchange: () => {
            st.xrd_ref = c.checked ? [...new Set([...st.xrd_ref, n])] : st.xrd_ref.filter((x) => x !== n); save();
          }});
          c.checked = st.xrd_ref.includes(n);
          refsBox.append(h("label", {}, c, " " + n));
        }
      };
      const refsDir = S.pathField(st.refs_dir, {dir: true, title: "Folder with reference structures POSCAR_<name>"});
      refsDir.input.addEventListener("change", () => { st.refs_dir = refsDir.input.value.trim(); st.xrd_ref = []; save(); drawRefs(); });
      drawRefs();
      const refSel = h("span", {});
      const drawRefSel = () => {
        refSel.innerHTML = "";
        const labs = labels();
        if (st.reference !== "lowest" && st.reference !== "none" && !labs.includes(st.reference)) st.reference = "lowest";
        refSel.append(sel("reference", [["lowest", "lowest F at each T"], ["none", "none (no ΔF)"], ...labs.map((l) => [l, l])]));
      };
      drawRefSel();
      const fmts = h("span", {class: "row"}, ["png", "pdf", "svg", "eps"].map((f) => {
        const c = h("input", {type: "checkbox", onchange: () => { st.formats = ["png", "pdf", "svg", "eps"].filter((x) => x === f ? c.checked : st.formats.includes(x)); save(); }});
        c.checked = st.formats.includes(f);
        return h("label", {}, c, " " + f.toUpperCase());
      }));

      const secRelax = h("div", {}, h("h3", {}, "Relaxation"), h("div", {class: "row"},
        h("label", {}, "fmax ", num("fmax", "80px"), " eV/Å"), h("label", {}, "max steps ", num("steps", "70px")),
        (() => { const c = h("input", {type: "checkbox", onchange: () => { st.fixed_cell = c.checked; save(); }}); c.checked = st.fixed_cell; return h("label", {}, c, " positions only (fixed cell)"); })()),
        h("div", {class: "small muted"}, "Same LBFGS settings and convergence check as the campaign (batch_free_energy.py); unconverged relaxations are flagged."));
      const secHarm = h("div", {}, h("h3", {}, "Harmonic free energy"), Tbox,
        h("div", {class: "row", style: {marginTop: "6px"}}, h("label", {}, "phonon supercell ", sel("scMode", [["auto", "automatic"], ["fixed", "fixed"]], st, drawSc)), scBox,
          h("label", {}, "q-mesh density ", num("mesh", "60px")), h("label", {}, "threads ", num("threads", "50px"))),
        h("div", {class: "small muted"}, "Gives F_total = E_static + F_vib, U_vib (with zero-point energy), ZPE, S and Cv per atom at every T, and the phonon DOS. Imaginary modes are flagged."));
      const secXrd = h("div", {}, h("h3", {}, "XRD (pymatgen, of the relaxed structure if relaxed)"),
        h("div", {class: "row"}, h("span", {class: "small"}, "references folder"), refsDir), refsBox,
        h("div", {class: "row", style: {marginTop: "4px"}}, h("label", {class: "grow"}, "more xrd_average.py options ", num("xrd_args", "320px")),
          h("span", {class: "small muted"}, "e.g. --wavelength MoKa --fwhm 0.08 --replace-cations Zr")));
      const secPlot = h("div", {}, h("h3", {}, "Plots"), h("div", {class: "row"},
        h("label", {}, "ΔF reference ", refSel), h("label", {}, "ΔF at ", num("at_T", "120px"), " K"), fmts),
        h("div", {class: "small muted"}, "F(T) and ΔF(T), S/Cv/U_vib/F_vib, ΔF at chosen T, phonon DOS, XRD stack, comparison CSV and an interactive HTML. You can re-plot with another reference afterwards."));
      const vis = () => {
        secRelax.style.display = st.tasks.includes("relax") ? "" : "none";
        secHarm.style.display = st.tasks.includes("harmonic") ? "" : "none";
        secXrd.style.display = st.tasks.includes("xrd") ? "" : "none";
        secPlot.style.display = st.tasks.includes("plot") ? "" : "none";
        nepRow.style.opacity = st.tasks.includes("relax") || st.tasks.includes("harmonic") ? "1" : ".45";
      };

      // where / how
      const slurmBox = h("div", {});
      const drawWhere = () => {
        slurmBox.innerHTML = "";
        if (st.mode === "local") return;
        const sl = st.slurm;
        slurmBox.append(h("div", {class: "row", style: {marginTop: "6px"}},
          h("label", {}, "partition ", num("partition", "110px", sl)), h("label", {}, "time ", num("time", "90px", sl)),
          h("label", {}, "cpus ", num("cpus", "50px", sl)), h("label", {}, "mem ", num("mem", "70px", sl)),
          h("label", {}, "account ", num("account", "100px", sl)), h("label", {}, "qos ", num("qos", "80px", sl))),
          h("div", {class: "form", style: {marginTop: "6px"}},
            h("label", {class: "k"}, "setup lines"), (() => { const t = h("textarea", {rows: 2, placeholder: "module load Python\nsource ~/venvs/nep/bin/activate", onchange: () => { sl.setup = t.value; save(); }}); t.value = sl.setup; return t; })(),
            h("label", {class: "k"}, "python on the node"), num("python", "260px", sl),
            h("label", {class: "k"}, "more #SBATCH lines"), (() => { const t = h("textarea", {rows: 1, placeholder: "#SBATCH --constraint=...", onchange: () => { sl.extra = t.value; save(); }}); t.value = sl.extra; return t; })()),
          !D.slurm && st.mode === "slurm_submit" ? h("div", {class: "notice warn small"}, "sbatch is not available on this machine: use 'write SLURM script' and submit it on the cluster.") : null);
      };
      drawWhere();
      const base = S.pathField(st.base, {dir: true, title: "Folder that holds the single-run folders"});
      base.input.addEventListener("change", () => { st.base = base.input.value.trim(); save(); drawRuns(); });
      const name = h("input", {type: "text", value: st.name, placeholder: "name (optional)", style: {width: "180px"}, onchange: () => { st.name = name.value.trim(); save(); }});
      const out = h("div", {});
      const go = async (mode) => {
        const spec = buildSpec(mode);
        if (!spec) return;
        const r = await S.api("/api/single/start", spec);
        if (!r.ok) return;
        out.innerHTML = "";
        if (mode === "preview") { out.append(h("pre", {class: "log"}, r.cmd)); return; }
        out.append(h("div", {class: "small"}, "Run folder ", h("a", {href: "javascript:void 0", onclick: () => showResult(r.folder)}, r.folder)));
        if (r.job) out.append(S.jobView(r.job, {title: "Single run", onDone: (j) => {
          out.innerHTML = "";            // the results card below shows the files: keep only a one-line note here
          out.append(h("div", {class: "small"}, h("span", {class: "pill " + (j.state === "done" ? "ok" : "bad")}, j.state), " ",
            h("a", {href: "javascript:void 0", onclick: () => showResult(r.folder)}, S.basename(r.folder)), " · ",
            h("a", {href: "#/jobs?id=" + r.job}, "log in Jobs")));
          drawRuns(); showResult(r.folder);
        }}));
        if (r.sbatch) out.append(h("div", {class: "notice ok small"}, "SLURM script written: ", h("span", {class: "mono"}, r.sbatch),
          " · submit on the cluster with ", h("span", {class: "mono"}, "sbatch " + r.sbatch), " · ", h("a", {href: "javascript:void 0", onclick: () => newCard.append(S.preview(r.sbatch))}, "show")));
        drawRuns();
      };
      const goLabel = () => (st.mode === "local" ? "Run" : st.mode === "slurm_submit" ? "Submit" : "Write script");
      const goBtn = h("button", {class: "primary", onclick: () => go(st.mode)}, goLabel());
      const nepRow = h("div", {class: "row"}, h("span", {style: {width: "120px"}}, "NEP model"), nep,
        h("button", {onclick: () => nepUp.click(), title: "upload nep.txt from your computer"}, "Upload…"), nepUp);
      newCard.append(
        h("h3", {}, "Structures"), list,
        h("div", {class: "row"}, h("button", {onclick: addPick}, "Add POSCAR…"), h("button", {onclick: addFolder}, "Add a folder…"),
          h("button", {onclick: () => up.click(), title: "upload POSCARs from your computer (stored in the run base folder, _uploads/)"}, "Upload from computer…"), up,
          h("button", {onclick: () => { st.poscars = []; st.names = []; save(); drawList(); drawRefSel(); }}, "Clear"),
          h("span", {class: "small muted"}, "VASP POSCARs; species are read from the files.")),
        h("h3", {}, "Model"), nepRow,
        h("h3", {}, "Tasks"), tasks,
        secRelax, secHarm, secXrd, secPlot,
        h("h3", {}, "Where"),
        h("div", {class: "row"}, h("span", {class: "small"}, "runs folder"), base, name),
        h("div", {class: "row", style: {marginTop: "6px"}}, h("label", {}, "run ", sel("mode", [["local", "here (this machine)"], ["slurm_submit", "submit to SLURM"], ["slurm_write", "write SLURM script only"]], st, () => { drawWhere(); goBtn.textContent = goLabel(); }))),
        slurmBox,
        h("div", {class: "row", style: {marginTop: "10px"}}, h("button", {onclick: () => go("preview")}, "Show command"),
          goBtn),
        out);
      vis();
    }

    function labels() {
      return st.poscars.map((p, i) => clean((st.names[i] || "").trim() || S.basename(p)));
    }
    function clean(s) { return String(s).replace(/[^\w.+-]+/g, "_").replace(/^[._]+|[._]+$/g, "").slice(0, 60); }
    async function uploadDir() {
      const d = (st.base || D.base) + "/_uploads/" + new Date().toISOString().slice(0, 10);
      const r = await S.api("/api/mkdir", {path: d});
      return r.ok ? r.path : null;
    }
    function nums(s) { return String(s || "").split(/[\s,;]+/).filter(Boolean); }
    function buildSpec(mode) {
      if (!st.poscars.length) { S.toast("add at least one POSCAR"); return null; }
      const spec = {mode, base: st.base, name: st.name, poscars: st.poscars, tasks: st.tasks, nep: st.nep,
        names: st.names.some((x) => x) ? labels() : [], fmax: st.fmax, steps: st.steps, fixed_cell: st.fixed_cell, mesh: st.mesh,
        threads: st.threads || null, xrd_ref: st.xrd_ref, refs_dir: st.refs_dir, xrd_args: st.xrd_args,
        reference: st.reference, at_T: nums(st.at_T), formats: st.formats, slurm: st.slurm};
      if (st.Tmode === "ramp") { spec.T_range = [st.Tmin, st.Tmax, st.Tstep]; spec.T = nums(st.Textra); }
      else spec.T = nums(st.T);
      if (st.scMode === "fixed") spec.supercell = nums(st.sc); else spec.min_length = st.min_length;
      return spec;
    }

    // ------------------------------------------------------------------ runs list
    async function drawRuns() {
      const r = await S.api("/api/single/list?base=" + encodeURIComponent(st.base || ""), undefined, {quiet: true});
      runsCard.innerHTML = "";
      runsCard.append(h("h2", {}, "Single runs in ", h("span", {class: "mono small"}, r.base || st.base)));
      if (!r.ok || !r.runs.length) { runsCard.append(h("div", {class: "muted small"}, "None yet.")); return; }
      const pill = (s) => h("span", {class: "pill " + ({done: "ok", done_with_errors: "warn", running: "run", written: "warn"}[s] || "run")}, s.replace(/_/g, " "));
      runsCard.append(h("div", {class: "tablewrap", style: {maxHeight: "300px"}}, h("table", {},
        h("thead", {}, h("tr", {}, ["run", "created", "structures", "tasks", "state"].map((x) => h("th", {class: "l"}, x)))),
        h("tbody", {}, r.runs.map((x) => h("tr", {},
          h("td", {class: "l"}, h("a", {href: "javascript:void 0", onclick: () => showResult(x.folder)}, x.name)),
          h("td", {class: "l small"}, x.created || ""), h("td", {}, x.n),
          h("td", {class: "l small"}, (x.tasks || []).join(" ")), h("td", {class: "l"}, pill(x.state), x.failed ? h("span", {class: "small muted"}, ` ${x.failed} failed`) : null)))))));
    }

    // ------------------------------------------------------------------ results
    async function showResult(folder) {
      const r = await S.api("/api/single/result?folder=" + encodeURIComponent(folder));
      if (!r.ok) return;
      resCard.innerHTML = "";
      const card = h("div", {class: "card"});
      resCard.append(card);
      const stt = r.status || {};
      const rq = r.request || {};
      card.append(h("div", {class: "row"}, h("h2", {class: "grow", style: {margin: 0}}, S.basename(r.folder)),
        h("button", {onclick: () => showResult(folder)}, "Refresh"),
        h("button", {onclick: () => copySettings(rq.spec || {})}, "Copy settings to the form")),
        h("div", {class: "small muted mono"}, r.folder),
        h("div", {class: "small"}, (stt.state ? `state: ${stt.state.replace(/_/g, " ")} · ${stt.done || 0}/${stt.n || "?"} structures` : "not started yet (or still in the SLURM queue)") +
          (stt.current ? ` · now: ${stt.current}` : "")),
        h("div", {class: "small muted"}, "Fixed-cell harmonic approximation (not quasi-harmonic): no thermal expansion."));
      if (r.slurm_out.length) card.append(h("div", {class: "small"}, "SLURM output: ", r.slurm_out.map((p) => h("a", {href: "javascript:void 0", onclick: () => card.append(S.preview(p))}, S.basename(p) + " "))));
      const flagged = r.structures.filter((s) => s.flags.length || /FAILED/.test(s.state));
      if (flagged.length) card.append(h("div", {class: "notice warn small"}, flagged.map((s) => h("div", {}, h("b", {}, s.label + ": "), s.notes.join("; ") || s.state))));
      if (r.summary.length) {
        const cols = [["label", "structure", 0], ["formula", "formula", 0], ["T_K", "T (K)", 0], ["E_static_eV_atom", "E static (eV/atom)", 6],
          ["F_total_eV_atom", "F total (eV/atom)", 6], ["ZPE_eV_atom", "ZPE (meV/atom)", -3], ["entropy_J_molK_atom", "S (J/mol K atom)", 3],
          ["volume_A3_atom", "V (Å³/atom)", 3], ["n_imaginary_modes", "imag. modes", 0], ["relax_converged", "relax conv.", 0], ["flags", "flags", 0]];
        const cell = (row, [k, , d]) => {
          const v = row[k];
          if (v === "" || v == null) return "–";
          if (d === -3) return S.fmt(1000 * Number(v), 2);
          return d ? S.fmt(Number(v), d) : v;
        };
        card.append(h("h3", {}, "Summary"), h("div", {class: "tablewrap"}, h("table", {},
          h("thead", {}, h("tr", {}, cols.map(([, l], i) => h("th", {class: i < 2 || i === 10 ? "l" : ""}, l)))),
          h("tbody", {}, r.summary.map((row) => h("tr", {}, cols.map((c, i) => h("td", {class: i < 2 || i === 10 ? "l" : ""}, cell(row, c)))))))),
          h("div", {class: "row small"}, h("a", {href: S.fileUrl(r.summary_csv, true)}, "summary.csv"),
            r.plots.some((p) => p.endsWith("comparison.csv")) ? h("a", {href: S.fileUrl(r.comparison_csv, true)}, "comparison.csv (every T)") : null));
      }
      // figures
      if (r.plots.length) {
        const prev = h("div", {});
        const figs = r.plots.filter((p) => /\.(png|html)$/.test(p));
        const pick = (p) => { prev.innerHTML = ""; prev.append(S.preview(p)); };
        card.append(h("h3", {}, "Figures"), h("div", {class: "row"}, figs.map((p) => h("a", {href: "javascript:void 0", onclick: () => pick(p)}, S.basename(p)))),
          h("div", {class: "row small muted"}, "download: ", r.plots.filter((p) => !/\.(png|html)$/.test(p)).map((p) => h("a", {href: S.fileUrl(p, true)}, S.basename(p)))),
          prev);
        const first = figs.find((p) => /harmonic_compare\.html$/.test(p)) || figs.find((p) => /F_vs_T\.png$/.test(p)) || figs[0];
        if (first) pick(first);
      }
      // per structure files
      if (r.structures.length) {
        card.append(h("details", {}, h("summary", {}, "Files of each structure (relaxed POSCAR, thermo.csv, DOS, XRD)"),
          r.structures.map((s) => h("div", {class: "small", style: {margin: "4px 0"}}, h("b", {}, s.label + " "), h("span", {class: "muted"}, s.state + " · "),
            s.files.map((p) => h("span", {}, h("a", {href: "javascript:void 0", onclick: () => { const x = S.preview(p); card.append(x); x.scrollIntoView(); }}, p.slice(p.indexOf("/results/") + 9 + s.label.length + 1)), " ",
              h("a", {href: S.fileUrl(p, true), class: "muted"}, "⤓"), "  "))))));
      }
      // re-plot
      if (r.summary.some((x) => x.F_total_eV_atom)) {
        const labs = r.summary.filter((x) => x.F_total_eV_atom).map((x) => x.label);
        const rp = Object.assign({reference: "lowest", include: [], tlo: "", thi: "", at_T: (rq.spec && rq.spec.at_T || []).join(" "), formats: ["png", "pdf", "svg"], width: 3.4, height: 2.6, font_size: 9}, S.store.get("single.replot", {}));
        if (!["lowest", "none", ...labs].includes(rp.reference)) rp.reference = "lowest";
        const ref = h("select", {}, h("option", {value: "lowest"}, "lowest F at each T"), h("option", {value: "none"}, "none"), labs.map((l) => h("option", {value: l}, l)));
        ref.value = rp.reference;
        const inc = h("div", {class: "row"}, labs.map((l) => { const c = h("input", {type: "checkbox"}); c.checked = !rp.include.length || rp.include.includes(l); c.dataset.l = l; return h("label", {}, c, " " + l); }));
        const f = (k, w) => { const i = h("input", {type: "text", value: rp[k], style: {width: w}}); i.dataset.k = k; return i; };
        const inputs = [f("tlo", "60px"), f("thi", "60px"), f("at_T", "110px"), f("width", "50px"), f("height", "50px"), f("font_size", "40px")];
        const jobBox = h("div", {});
        card.append(h("h3", {}, "Re-plot (other reference, structures, T window, size)"),
          h("div", {class: "row"}, h("label", {}, "ΔF reference ", ref), h("label", {}, "T from ", inputs[0]), h("label", {}, "to ", inputs[1]),
            h("label", {}, "ΔF at ", inputs[2], " K"), h("label", {}, "width ", inputs[3], " in"), h("label", {}, "panel height ", inputs[4], " in"), h("label", {}, "font ", inputs[5], " pt")),
          h("div", {class: "row small"}, "structures: ", inc),
          h("div", {class: "row", style: {marginTop: "6px"}}, h("button", {class: "primary", onclick: async () => {
            inputs.forEach((i) => { rp[i.dataset.k] = i.value.trim(); });
            rp.reference = ref.value;
            const chosen = [...inc.querySelectorAll("input")].filter((c) => c.checked).map((c) => c.dataset.l);
            rp.include = chosen.length === labs.length ? [] : chosen;
            if (rp.reference !== "lowest" && rp.reference !== "none" && rp.include.length && !rp.include.includes(rp.reference)) return S.toast("the reference must be one of the chosen structures");
            S.store.set("single.replot", rp);
            const body = {folder: r.folder, reference: rp.reference, include: rp.include, at_T: nums(rp.at_T), formats: rp.formats,
              width: rp.width, height: rp.height, font_size: rp.font_size, T_window: rp.tlo !== "" && rp.thi !== "" ? [rp.tlo, rp.thi] : null};
            const j = await S.api("/api/single/replot", body);
            if (j.ok) { jobBox.innerHTML = ""; jobBox.append(S.jobView(j.job, {title: "Re-plot", onDone: () => showResult(r.folder)})); }
          }}, "Re-plot")), jobBox);
      }
      card.scrollIntoView({behavior: "smooth", block: "start"});
    }

    function copySettings(spec) {
      const keep = ["tasks", "fmax", "steps", "fixed_cell", "mesh", "threads", "xrd_ref", "refs_dir", "xrd_args", "reference", "formats"];
      for (const k of keep) if (spec[k] != null) st[k] = spec[k];
      st.poscars = []; st.names = [];
      if (spec.T_range) { st.Tmode = "ramp"; [st.Tmin, st.Tmax, st.Tstep] = spec.T_range; st.Textra = (spec.T || []).join(" "); }
      else if (spec.T) { st.Tmode = "const"; st.T = spec.T.join(" "); }
      if (spec.supercell) { st.scMode = "fixed"; st.sc = spec.supercell.join(" "); } else if (spec.min_length) { st.scMode = "auto"; st.min_length = spec.min_length; }
      st.at_T = (spec.at_T || []).join(" ");
      if (spec.nep) st.nep = spec.nep;
      if (spec.slurm) st.slurm = Object.assign({}, DEF.slurm, spec.slurm);
      if (spec.mode && spec.mode !== "preview") st.mode = spec.mode;
      save(); drawNew();
      newCard.scrollIntoView({behavior: "smooth"});
      S.toast("settings copied; add the structures");
    }
  }
})();
