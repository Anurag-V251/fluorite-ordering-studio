/* Pipeline page: create a reproducible run (campaign.toml + run_record.json) and drive its stages. */
(function () {
  // "label=path" -> [label, path]; only the first "=" counts, and only before a valid label (else a bare path)
  const splitTpl = (t) => { const i = t.indexOf("="); const l = i > 0 ? t.slice(0, i) : "";
    return /^[\w.+\-]+$/.test(l) ? [l, t.slice(i + 1)] : ["", t]; };
  const S = Studio, h = S.h;
  S.registerPage({id: "pipeline", title: "Pipeline", order: 30, render});

  async function render(main, params) {
    if (params.config) return showRun(main, params.config);
    if (params.new) return newRun(main, params.from);
    const recent = (S.info.prefs || {}).recent_runs || [];
    const pf = S.pathField("", {dir: true, title: "Run folder (contains campaign.toml)"});
    main.append(h("div", {class: "card"}, h("h2", {}, "Runs"),
      h("div", {class: "row"}, h("button", {class: "primary", onclick: () => S.go("pipeline", {new: 1})}, "New run…"),
        h("span", {class: "muted"}, "or open an existing run folder:"), pf,
        h("button", {onclick: async () => {
          const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(pf.input.value));
          if (r.ok) S.go("pipeline", {config: r.config_path});
        }}, "Open")),
      h("div", {class: "notice ok small", style: {marginTop: "10px"}},
        "Workflow: New run (writes the folder) → 1 dry-run → 2 prepare (enumerate, map, verify) → 3 plan → 4 compute → 5 status → 6 merge → (7 rerun flagged → 6 merge) → Results and Plots. ",
        "Compositions are dopant counts n (0 … number of host sites), not %."),
      recent.length ? h("div", {}, h("h3", {}, "Recent"), recent.map((c) => h("div", {}, h("a", {href: "#/pipeline?config=" + encodeURIComponent(c)}, c)))) : null,
      S.runsList((x) => x.config && S.go("pipeline", {config: x.config})),
      h("div", {class: "small muted", style: {marginTop: "10px"}},
        "A run folder holds campaign.toml (what campaign.py reads), run_record.json (settings, md5 of every input and script, package versions, history of every command) and, by default, copies of the inputs and a frozen copy of the scripts, so the run can be reproduced from the folder alone. Several starting phases = several runs.")));
  }

  // ------------------------------------------------------------------ new run form
  async function newRun(main, from) {
    const d = (await S.api("/api/runs/defaults")).defaults;
    let pre = null;
    if (from) {
      const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(from));
      if (r.ok) pre = r;
    }
    const cfg = JSON.parse(JSON.stringify(d));
    const preCfg = pre ? (pre.record && pre.record.config) || pre.config : null;
    const base = pre ? pre.folder : "";
    if (preCfg) for (const sec of Object.keys(cfg)) Object.assign(cfg[sec], preCfg[sec] || {});
    const abs = (p) => (!p ? "" : p.startsWith("/") ? p : base + "/" + p);
    if (pre) {
      cfg.structure.poscar = abs(cfg.structure.poscar);
      cfg.free_energy.nep = abs(cfg.free_energy.nep);
      cfg.mapping.templates = (cfg.mapping.templates || []).map((t) => { const [l, p] = splitTpl(t); return (l ? l + "=" : "") + abs(p); });
    }
    const F = {};
    const txt = (sec, k, w) => (F[sec + "." + k] = h("input", {type: "text", value: fmtv(cfg[sec][k], k), style: {width: w || "140px"}}));
    const chk = (sec, k) => { const c = h("input", {type: "checkbox"}); c.checked = !!cfg[sec][k]; F[sec + "." + k] = c; return c; };
    const rowf = (label, el, help) => [h("label", {class: "k"}, label), h("div", {}, el, help ? h("div", {class: "help"}, help) : null)];

    const folder = S.pathField(pre ? pre.folder + "_copy" : S.info.start_dir + "/runs/run_" + new Date().toISOString().slice(0, 16).replace(/[-:]/g, "").replace("T", "_"), {dir: true, title: "Run folder (created if missing)"});
    const name = h("input", {type: "text", class: "grow", value: pre && pre.record ? pre.record.name : ""});
    const notes = h("textarea", {rows: 2}, pre && pre.record ? pre.record.notes || "" : "");
    const startPhase = h("input", {type: "text", value: pre && pre.record ? pre.record.starting_phase || "c" : "c", style: {width: "80px"}});
    const poscar = S.pathField(cfg.structure.poscar, {title: "Host supercell POSCAR (starting phase)"});
    const peek = h("div", {class: "small"});
    const hostSel = h("input", {type: "text", value: cfg.structure.host, style: {width: "80px"}});
    F["structure.host"] = hostSel;
    poscar.input.addEventListener("change", doPeek);
    async function doPeek() {
      peek.innerHTML = "";
      if (!poscar.input.value) return;
      const r = await S.api("/api/poscar?path=" + encodeURIComponent(poscar.input.value), undefined, {quiet: true});
      if (!r.ok) { peek.append(h("span", {style: {color: "var(--bad)"}}, r.error)); return; }
      const sp = Object.entries(r.species).map(([e, n]) => e + n).join(" ");
      const cations = Object.keys(r.species).filter((e) => !["O", "F", "S", "N", "Se", "Te", "Cl"].includes(e));
      let note = null;
      if (!(hostSel.value in r.species) && cations.length === 1) {   // another material: take its cation as the host
        note = h("div", {class: "small", style: {color: "var(--ok)"}}, `host element set to ${cations[0]} (the cation of this POSCAR); check the dopant below`);
        hostSel.value = cations[0];
      }
      peek.append(h("span", {class: "mono"}, `${sp} · ${r.n_atoms} atoms · lengths ${r.lengths.join(" × ")} Å` +
        (r.repeat_guess ? ` · looks like ${r.repeat_guess.join("×")} of the 12-atom cell` : "")),
        note, !(hostSel.value in r.species) ? h("div", {style: {color: "var(--warn)"}}, `host element ${hostSel.value} is not in this POSCAR (${Object.keys(r.species).join(", ")})`) : null);
    }
    hostSel.addEventListener("change", doPeek);

    const tplBox = h("div", {});
    const addTpl = (label, path) => {
      const l = h("input", {type: "text", value: label || "", placeholder: "label (t, oIII, m…)", style: {width: "110px"}});
      const p = S.pathField(path || "", {title: "12-atom template POSCAR of the target phase"});
      const r = h("div", {class: "row tplrow", style: {marginBottom: "4px"}}, l, p, h("button", {onclick: () => r.remove()}, "✕"));
      r.get = () => (p.input.value ? (l.value || S.basename(p.input.value).replace(/^POSCAR[_.-]?/, "")) + "=" + p.input.value : null);
      tplBox.append(r);
    };
    for (const t of cfg.mapping.templates || []) { const [l, p] = splitTpl(t); addTpl(l, p); }
    if (!(cfg.mapping.templates || []).length) addTpl("t", "");
    const fromRefs = async () => {
      const dir = await S.pick({dir: true, title: "Folder with POSCAR_<phase> templates"});
      if (!dir) return;
      const r = await S.api("/api/ls?path=" + encodeURIComponent(dir));
      tplBox.querySelectorAll(".tplrow").forEach((x) => { if (!x.get()) x.remove(); });
      for (const e of r.entries || []) if (!e.dir && /^POSCAR_/.test(e.name)) {
        const lab = e.name.replace(/^POSCAR_/, "");
        if (lab === startPhase.value) continue;            // the starting phase itself is the parent
        addTpl(lab, dir + "/" + e.name);
      }
    };
    const allVar = h("input", {type: "checkbox"});
    allVar.checked = /--all-variants/.test(cfg.mapping.map_args);
    const nep = S.pathField(cfg.free_energy.nep, {title: "NEP model (nep.txt)"});
    const setup = h("textarea", {rows: 4}, cfg.cluster.setup || "");
    const sbx = h("textarea", {rows: 2, placeholder: "one per line, e.g. --exclusive"}, (cfg.cluster.sbatch_extra || []).join("\n"));
    const copyIn = h("input", {type: "checkbox", checked: true});
    const freeze = h("input", {type: "checkbox", checked: true});
    const overwrite = h("input", {type: "checkbox"});
    const msg = h("div", {});

    const form = h("div", {class: "form"},
      h("div", {class: "grp"}, "Run"),
      rowf("folder", folder, "created if missing; holds campaign.toml, run_record.json, inputs/, scripts/ and the workdir"),
      rowf("name", name), rowf("notes", notes),
      h("div", {class: "grp"}, "Material and starting phase"),
      rowf("host POSCAR", h("div", {}, poscar, peek), "starting phase: VASP5 supercell, diagonal repeat of the 12-atom cell, axes along the cube directions"),
      rowf("starting phase label", startPhase, "name of the phase of the host POSCAR (for the records and plot legends; campaign.py calls it 'parent')"),
      rowf("host element", hostSel, "element on the substituted sites"),
      rowf("dopant element", txt("structure", "dopant", "80px")),
      rowf("compositions", txt("structure", "compositions", "220px"), "dopant counts: all | start:stop[:step] (inclusive) | list, e.g. 0, 6, 16, 32"),
      rowf("max orderings", txt("structure", "max_orderings"), "distinct random orderings per composition (all of them where fewer exist)"),
      rowf("seed", txt("structure", "seed"), "composition n uses seed*1000+n"),
      h("div", {class: "grp"}, "Phases for transformation (same cation ordering)"),
      rowf("target templates", h("div", {}, tplBox, h("div", {class: "row"}, h("button", {onclick: () => addTpl()}, "+ phase"), h("button", {onclick: fromRefs}, "Add all from a folder…"))),
        "12-atom references (4 cations + 8 O) of each target phase; every ordering is mapped onto each with map_fluorite.py and verified with compare_ordering.py (must be IDENTICAL)"),
      rowf("all variants", allVar, "compute every equally good orientation (6 for t, 48 for oIII / m from a cubic parent); they can only be ranked by energy"),
      rowf("other map_fluorite options", txt("mapping", "map_args", "260px")),
      rowf("include parent", chk("mapping", "include_parent"), "also compute the un-mapped orderings"),
      rowf("verify", chk("mapping", "verify")),
      rowf("mapping jobs", txt("mapping", "jobs")),
      h("div", {class: "grp"}, "NEP relaxation and harmonic free energy"),
      rowf("NEP model", nep),
      rowf("T (K)", txt("free_energy", "T")),
      rowf("relax", chk("free_energy", "relax")), rowf("fmax (eV/Å)", txt("free_energy", "fmax"), "rows above it after `steps` are flagged relax_converged = False"),
      rowf("steps", txt("free_energy", "steps")), rowf("fixed cell", chk("free_energy", "fixed_cell"), "relax positions only (c/a and θ then stay the template values)"),
      rowf("min_length (Å)", txt("free_energy", "min_length"), "automatic phonon supercell: every width ≥ this"),
      rowf("phonon supercell", txt("free_energy", "supercell"), "empty = automatic, or e.g. 1 1 1"),
      rowf("mesh", txt("free_energy", "mesh")), rowf("threads", txt("free_energy", "threads"), "BLAS/OpenMP threads per structure"),
      h("div", {class: "grp"}, "Cluster (SLURM) — also used for the chunk plan when running locally"),
      ...["cores_per_node", "walltime", "fill", "sec_per_structure", "max_array", "max_concurrent", "prepare_walltime", "merge_walltime",
        "partition", "account", "qos", "mem", "python"].flatMap((k) => rowf(k, txt("cluster", k), {
        cores_per_node: "set it to what one job really uses (e.g. 256 on nodes with 128 cores × 2 threads, whole node allocated)",
        sec_per_structure: "a guess until measured: run one chunk, read Status, put the value here, re-plan with --force",
        max_array: "MaxArraySize of your cluster", python: "python used inside the SLURM jobs"}[k])),
      rowf("setup", setup, "shell lines run at the start of every job (module load …, source …/activate, export RAYON_NUM_THREADS=1)"),
      rowf("sbatch extra", sbx),
      h("div", {class: "grp"}, "Reproducibility"),
      rowf("copy inputs", copyIn, "copy POSCAR, templates and NEP into <run>/inputs (md5 recorded either way)"),
      rowf("freeze scripts", freeze, "copy the scripts into <run>/scripts so later changes of the app do not change this run"),
      rowf("overwrite settings", overwrite, "allow rewriting campaign.toml of an existing run (campaign.py refuses changed templates / map_args in a used workdir)"));

    main.append(h("div", {class: "card"}, h("h2", {}, pre ? "New run from " + pre.folder : "New run"), form,
      msg, h("div", {class: "row", style: {marginTop: "12px"}}, h("button", {class: "primary", onclick: create}, "Create run"),
        h("span", {class: "muted small"}, "then: Pipeline → dry-run to check"))));
    doPeek();

    async function create() {
      const c = JSON.parse(JSON.stringify(d));
      for (const [k, el] of Object.entries(F)) { const [sec, key] = k.split("."); c[sec][key] = el.type === "checkbox" ? el.checked : el.value; }
      c.structure.poscar = poscar.input.value;
      c.free_energy.nep = nep.input.value;
      c.mapping.templates = [...tplBox.querySelectorAll(".tplrow")].map((r) => r.get()).filter(Boolean);
      let ma = String(c.mapping.map_args || "").replace(/--all-variants/g, "").trim();
      c.mapping.map_args = allVar.checked ? ("--all-variants " + ma).trim() : ma;
      c.cluster.setup = setup.value;
      c.cluster.sbatch_extra = sbx.value;
      if (c.structure.compositions.trim() === "all" || /^\d+:\d+(:\d+)?$/.test(c.structure.compositions.trim())) c.structure.compositions = c.structure.compositions.trim();
      const btn = document.activeElement && document.activeElement.tagName === "BUTTON" ? document.activeElement : null;
      if (create.busy) return;
      create.busy = true; if (btn) btn.disabled = true;     // no double submit
      let r;
      try {
        r = await S.api("/api/runs/create", {folder: folder.input.value, name: name.value, notes: notes.value,
          starting_phase: startPhase.value, copy_inputs: copyIn.checked, freeze_scripts: freeze.checked, overwrite: overwrite.checked, config: c});
      } finally { create.busy = false; if (btn) btn.disabled = false; }
      msg.innerHTML = "";
      if (!r.ok) { msg.append(h("div", {class: "notice bad"}, r.error)); msg.scrollIntoView({block: "center"}); return; }
      S.info.prefs.recent_runs = [r.config, ...((S.info.prefs || {}).recent_runs || [])];
      S.go("pipeline", {config: r.config});
    }
  }
  function fmtv(v, k) { return Array.isArray(v) ? v.join(k === "supercell" ? " " : ", ") : v == null ? "" : String(v); }

  // ------------------------------------------------------------------ one run: stages
  async function showRun(main, config) {
    const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(config));
    if (!r.ok) { main.append(h("div", {class: "notice bad"}, r.error)); return; }
    const rec = r.record;
    const repro = h("div", {});
    const head = h("div", {class: "card"},
      h("div", {class: "row"}, h("h2", {style: {margin: 0}}, rec && rec.name ? rec.name : S.basename(r.folder)),
        rec && rec.starting_phase ? h("span", {class: "pill run"}, "starting phase " + rec.starting_phase) : null,
        h("span", {class: "mono small muted grow"}, r.config_path),
        h("a", {href: "#/results?path=" + encodeURIComponent(r.workdir || r.folder)}, "results"),
        h("a", {href: "#/pipeline?new=1&from=" + encodeURIComponent(r.folder)}, "copy settings to a new run"),
        rec ? h("button", {onclick: async () => {
          const c = await S.api("/api/runs/check?folder=" + encodeURIComponent(r.folder));
          repro.innerHTML = "";
          if (!c.ok) return repro.append(h("div", {class: "notice bad"}, c.error));
          repro.append(h("div", {class: "notice " + (c.identical ? "ok" : "warn")}, c.identical ? "All inputs and scripts have the md5 recorded when the run was created." :
            [c.differences.length + " file(s) differ from the record:",
              h("ul", {style: {margin: "4px 0"}}, c.differences.map((x) => h("li", {class: "mono small"},
                `${x.file}: md5 ${(x.recorded || "-").slice(0, 10)}… → ${x.now ? x.now.slice(0, 10) + "…" : "missing"}`))),
              c.differences.some((x) => x.file === "campaign.toml") ? h("div", {class: "small"}, "campaign.toml changed: compare it with the settings in run_record.json (the 'config' entry).") : null]),
            Object.keys(c.version_differences).length ? h("div", {class: "notice warn small"}, "package versions differ: " +
              Object.entries(c.version_differences).map(([k, v]) => `${k} ${v[0]} → ${v[1]}`).join(", ")) : null);
        }}, "Check reproducibility") : null),
      rec ? h("div", {class: "small muted"}, `created ${rec.created} by ${rec.user}@${rec.host} · app ${rec.app_version} · python ${rec.versions.python}` +
        ["numpy", "ase", "calorine", "phonopy"].map((k) => rec.versions[k] ? ` · ${k} ${rec.versions[k]}` : "").join("")) :
        h("div", {class: "small muted"}, "no run_record.json (folder not created by the app): stages still work"),
      repro,
      h("details", {}, h("summary", {}, "campaign.toml"), (() => { const p = h("pre", {class: "log"}); S.api("/api/text?path=" + encodeURIComponent(r.config_path)).then((t) => { p.textContent = t.text; }); return p; })()),
      rec && rec.notes ? h("div", {class: "small"}, rec.notes) : null);
    main.append(head);
    const stages = S.mods("stage");
    const jobs = h("div", {});
    const grid = h("div", {style: {display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(330px, 1fr))", gap: "0 14px"}});
    stages.forEach((m, i) => grid.append(stageCard(m, r, jobs, i + 1)));
    const clear = h("button", {class: "small", onclick: () => { jobs.innerHTML = ""; }}, "Clear finished cards");
    main.append(grid, h("div", {class: "row", style: {margin: "0 0 8px"}}, h("b", {}, "Output of the stages started here (newest first)"), clear),
      h("div", {}, h("div", {}, jobs,
      rec && rec.history && rec.history.length ? h("details", {class: "card"}, h("summary", {}, `History (${rec.history.length})`),
        h("pre", {class: "log"}, rec.history.slice().reverse().map((x) => `${x.time}  ${x.stage || x.event || ""} [${x.mode || ""}]  ${x.cmd || ""}`).join("\n"))) : null)));
    // stages of this run still running (started before a reload, or in another tab): show them again
    const jl = await S.api("/api/jobs", undefined, {quiet: true});
    if (jl.ok) for (const j of jl.jobs.filter((j) => ["running", "starting"].includes(j.state) && (j.meta || {}).config === r.config_path).reverse())
      jobs.prepend(S.jobView(j.id, {title: j.title}));
    if (!jobs.childElementCount) jobs.append(h("div", {class: "card muted"}, "Start a stage above. Its log and output appear here. Order: dry-run → prepare → plan → compute → status → merge → (rerun flagged → merge)."));
  }

  function stageCard(m, run, jobs, num) {
    const card = h("div", {class: "card", style: {padding: "10px 12px"}});
    const optEls = {};
    const opts = h("div", {});
    for (const o of m.options || []) {
      const key = o.flag || "__free__";
      const el = o.value ? h("input", {type: "text", class: "grow mono", placeholder: o.flag || "arguments"}) : h("input", {type: "checkbox"});
      optEls[key] = el;
      const lbl = h("span", {}, o.label || o.flag, o.modes ? h("span", {class: "muted"}, " (" + o.modes.join(", ") + ")") : null);
      opts.append(o.value ? h("div", {class: "small", style: {margin: "4px 0"}}, h("div", {class: "muted"}, lbl), el)
        : h("label", {class: "small", style: {display: "flex", gap: "6px", alignItems: "flex-start", margin: "3px 0"}}, el, lbl));
    }
    const getOpts = () => Object.fromEntries(Object.entries(optEls).map(([k, el]) => [k, el.type === "checkbox" ? el.checked : el.value]));
    const btns = h("div", {class: "row"});
    const order = ["dry", "local", "slurm", "slurm_write"];
    const rank = (k) => (order.indexOf(k) < 0 ? 99 : order.indexOf(k));
    const modes = Object.entries(m.modes || {}).sort((a, b) => rank(a[0]) - rank(b[0]));
    for (const [mode, md] of modes) {
      if (mode.startsWith("slurm") && !S.info.slurm && mode === "slurm") {
        btns.append(h("button", {disabled: true, "data-off": "1", title: "SLURM not found on this machine (squeue is not on PATH): use 'write script' and submit on the cluster"}, md.label || mode));
        continue;
      }
      btns.append(h("button", {class: mode === "local" && modes.length === 1 ? "primary" : "", onclick: () => launch(mode, md)}, md.label || mode));
    }
    btns.append(h("a", {href: "javascript:void 0", class: "small", onclick: async () => {
      const lines = [];
      for (const [mode] of modes) {
        const p = await S.api("/api/stage", {module: m.id, mode, config: run.config_path, options: getOpts(), preview: true});
        if (p.ok) lines.push(`# ${mode}\n${p.cmd}`);
      }
      pre.style.display = ""; pre.textContent = lines.join("\n\n");
    }}, "commands"));
    const pre = h("pre", {class: "log", style: {display: "none"}});
    card.append(h("b", {}, num + ". " + m.title), m.description ? h("div", {class: "small muted"}, m.description) : null, opts, btns, pre);
    async function launch(mode, md) {
      const bs = [...btns.querySelectorAll("button")];
      bs.forEach((b) => (b.disabled = true));             // no double submit
      let r;
      try { r = await S.api("/api/stage", {module: m.id, mode, config: run.config_path, options: getOpts()}); }
      finally { setTimeout(() => bs.forEach((b) => (b.disabled = b.dataset.off === "1")), 1500); }
      if (!r.ok) return;
      const first = jobs.querySelector(".muted.card");
      if (first) first.remove();
      jobs.prepend(S.jobView(r.job, {title: `${m.title} · ${md.label || mode}`}));
    }
    return card;
  }
})();
