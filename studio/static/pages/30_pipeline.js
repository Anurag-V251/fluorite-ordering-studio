/* Pipeline page: create a reproducible run (campaign.toml + run_record.json) and drive its stages. */
(function () {
  // "label=path" -> [label, path]; only the first "=" counts, and only before a valid label (else a bare path)
  const splitTpl = (t) => { const i = t.indexOf("="); const l = i > 0 ? t.slice(0, i) : "";
    return /^[\w.+\-]+$/.test(l) ? [l, t.slice(i + 1)] : ["", t]; };
  const S = Studio, h = S.h;
  S.registerPage({id: "pipeline", title: "Pipeline", order: 10, icon: "flow", render,
    navAs: (p) => (p.new ? "new" : null)});
  S.registerPage({id: "new", title: "New run", order: 15, icon: "plus", render: (main, p) => newRun(main, p.from)});

  async function render(main, params) {
    if (params.config) return showRun(main, params.config);
    if (params.new) return newRun(main, params.from);
    const recent = (S.info.prefs || {}).recent_runs || [];
    if (!params.list && recent.length) {        // straight to the run you had open last (if it still exists)
      const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(recent[0]), undefined, {quiet: true});
      if (r.ok && r.config_path) return S.go("pipeline", {config: r.config_path});
    }
    return runsPage(main);
  }

  function runsPage(main) {
    const recent = (S.info.prefs || {}).recent_runs || [];
    const pf = S.pathField("", {dir: true, title: "Run folder (campaign.toml), or a campaign workdir"});
    const open = async () => {
      const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(pf.input.value));
      if (r.ok) S.go("pipeline", {config: r.config_path});
    };
    pf.input.addEventListener("keydown", (e) => { if (e.key === "Enter") open(); });
    main.append(
      h("div", {class: "runhead"}, h("div", {}, h("div", {class: "crumb"}, "Runs"), h("h1", {}, "Runs"),
        h("div", {class: "runmeta"}, h("span", {}, "One run = one starting phase of one material. Open a run to see its pipeline, jobs and results."))),
        h("div", {class: "actions"}, h("a", {class: "btn primary", href: "#/new"}, S.icon("plus"), "New run"))),
      h("div", {class: "card"}, h("h2", {}, "Open a run folder"),
        h("div", {class: "row"}, pf, h("button", {class: "btn primary", onclick: open}, "Open")),
        h("div", {class: "small muted", style: {marginTop: "8px"}},
          "A folder made by New run, a folder with your own campaign.toml, or a campaign workdir (it has campaign_config.toml) all work. ",
          "Runs made without the app keep working as they are.")),
      recent.length ? h("div", {class: "card"}, h("h2", {}, "Recent"), recent.slice(0, 8).map((c) =>
        h("div", {style: {padding: "3px 0"}}, h("a", {href: "#/pipeline?config=" + encodeURIComponent(c), class: "mono"}, c)))) : null,
      h("div", {class: "card"}, S.runsList((x) => x.config && S.go("pipeline", {config: x.config}))),
      h("div", {class: "small muted"},
        "A run folder holds campaign.toml (what campaign.py reads), run_record.json (settings, md5 of every input and script, package versions, history of every command) and, by default, copies of the inputs and a frozen copy of the scripts, so the run can be reproduced from the folder alone."));
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

  // ------------------------------------------------------------------ one run: dashboard
  async function showRun(main, config) {
    const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(config));
    if (!r.ok) { main.append(h("div", {class: "notice bad"}, r.error), h("a", {href: "#/pipeline?list=1"}, "all runs")); return; }
    const rec = r.record;
    S.currentRun = r;
    const head = h("div", {class: "runhead"});
    const strip = h("div", {class: "card strip"}, h("div", {class: "left"}, "reading the run…"), h("div", {class: "next"}));
    const warnBox = h("div", {});
    const repro = h("div", {});
    main.append(head, warnBox, strip);
    let st = null, polled = 0;
    const pollAgo = h("span", {class: "muted num"});
    const draw = () => {
      drawHead(head, r, st, pollAgo, poll, repro);
      if (st) { drawStrip(strip, st, r); warnBox.innerHTML = ""; (st.warnings || []).forEach((w) => warnBox.append(h("div", {class: "notice warn"}, w))); }
      S.emit("run-state", st);
    };
    async function poll(fresh) {
      const s2 = await S.api("/api/run/state?config=" + encodeURIComponent(r.config_path) + (fresh ? "&fresh=1" : ""), undefined, {quiet: !fresh});
      if (!s2.ok) return;
      st = s2; polled = Date.now(); S.runState = st; draw();
    }
    S.every(1000, () => { pollAgo.textContent = polled ? Math.round((Date.now() - polled) / 1000) + " s" : ""; });
    S.every(15000, () => poll(false));
    draw();
    await poll(false);
    main.append(repro);

    // the stages (cards with their options; the next version shows them as a pipeline with a side panel)
    const stages = S.mods("stage");
    const jobs = h("div", {});
    const grid = h("div", {style: {display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(330px, 1fr))", gap: "0 14px"}});
    stages.forEach((m, i) => grid.append(stageCard(m, r, jobs, i + 1)));
    const clear = h("button", {class: "btn sm", onclick: () => { jobs.innerHTML = ""; }}, "Clear finished cards");
    main.append(h("div", {class: "card-h", style: {padding: "4px 0 10px"}}, h("h2", {}, "Stages"),
        h("span", {class: "sub"}, "Run here, write the SLURM scripts, or submit them. Every stage keeps the equivalent command.")),
      grid, h("div", {class: "row", style: {margin: "0 0 8px"}}, h("b", {}, "Output of the stages started here (newest first)"), clear),
      h("div", {}, h("div", {}, jobs,
      rec && rec.history && rec.history.length ? h("details", {class: "card"}, h("summary", {}, `History (${rec.history.length})`),
        h("pre", {class: "log"}, rec.history.slice().reverse().map((x) => `${x.time}  ${x.stage || x.event || ""} [${x.mode || ""}]  ${x.cmd || ""}`).join("\n"))) : null,
      h("details", {class: "card"}, h("summary", {}, "campaign.toml"), (() => { const p = h("pre", {class: "log"}); S.api("/api/text?path=" + encodeURIComponent(r.config_path)).then((t) => { p.textContent = t.text; }); return p; })()))));
    S.on("stage-started", () => poll(true));
    // stages of this run still running (started before a reload, or in another tab): show them again
    const jl = await S.api("/api/jobs", undefined, {quiet: true});
    if (jl.ok) for (const j of jl.jobs.filter((j) => ["running", "starting"].includes(j.state) && (j.meta || {}).config === r.config_path).reverse())
      jobs.prepend(S.jobView(j.id, {title: j.title, onDone: () => poll(true)}));
    if (!jobs.childElementCount) jobs.append(h("div", {class: "card muted"}, "Start a stage above. Its log and output appear here. Order: dry-run → prepare → plan → compute → status → merge → (rerun flagged → merge)."));
  }

  function runStatus(st) {
    if (!st) return S.pill("idle", "…");
    const run = st.running || [];
    if (run.includes("compute") || st.compute_running) return S.pill("run", "Computing");
    if (run.includes("prepare")) return S.pill("run", "Preparing");
    if (run.includes("rerun_flagged")) return S.pill("run", "Rerunning");
    if (run.length) return S.pill("run", "Running");
    if ((st.attention || []).some((a) => a.sev !== "info")) return S.pill("warn", "Needs attention");
    if (st.next && st.next.action && st.next.action.go === "results") return S.pill("ok", "Complete");
    if (!st.counts.structures) return S.pill("idle", "Not started");
    return S.pill("acc", "In progress");
  }

  function drawHead(head, r, st, pollAgo, poll, repro) {
    const rec = r.record || {};
    const m = st ? st.material : {};
    const cl = st ? st.cluster : {};
    const name = rec.name || S.basename(r.folder);
    const sw = h("button", {class: "select-run", title: "Switch run", "aria-haspopup": "listbox",
      onclick: (e) => runSwitcher(e.currentTarget, r)}, h("span", {class: "mono", style: {overflow: "hidden", textOverflow: "ellipsis"}}, S.basename(r.folder)), S.icon("chevdown", "sm"));
    const created = rec.created ? new Date(rec.created.replace(" ", "T")) : null;
    const el = created && !isNaN(created) ? (Date.now() - created) / 1000 : null;
    const elTxt = el == null ? "" : el > 86400 ? `${Math.floor(el / 86400)} d ${S.hms(el % 86400).slice(0, -3)}` : S.hms(el).slice(0, -3);
    head.innerHTML = "";
    head.append(h("div", {style: {minWidth: 0}},
      h("div", {class: "crumb"}, h("a", {href: "#/pipeline?list=1"}, "Runs"), S.icon("chev", "sm"), sw),
      h("h1", {}, name, rec.starting_phase ? h("span", {class: "pill p-idle", title: "starting phase (the phase of the host POSCAR)"}, "from " + rec.starting_phase) : null, runStatus(st)),
      h("div", {class: "runmeta"},
        m.n_atoms ? h("span", {}, S.icon("layers", "sm"), `${m.n_atoms}-atom host · ${m.n_host ?? "?"} ${m.host || ""} sites → ${m.dopant || "?"}`) : null,
        m.compositions ? h("span", {}, S.icon("grid", "sm"), `${m.compositions.length} compositions · ${(m.phases || []).length} phases · T ${m.T} K`) : null,
        m.nep ? h("span", {}, S.icon("flask", "sm"), "NEP " + m.nep) : null,
        cl.cores_per_node ? h("span", {}, S.icon("cpu", "sm"), (cl.partition ? `partition ${cl.partition} · ` : "") + `${cl.cores_per_node} cores/node · ${cl.walltime}`) : null,
        rec.created ? h("span", {}, S.icon("clock", "sm"), `created ${rec.created.slice(0, 16)}` + (elTxt ? ` · ${elTxt} ago` : "")) : null,
        !rec.created ? h("span", {class: "t2"}, S.icon("alert", "sm"), "no run_record.json (made without the app)") : null)),
      h("div", {class: "actions"},
        h("button", {class: "btn", onclick: () => poll(true), title: "Read the files, squeue and sacct again"}, S.icon("refresh"), "Refresh ", pollAgo),
        h("a", {class: "btn", href: "#/results?path=" + encodeURIComponent(r.workdir || r.folder)}, S.icon("table"), "Results"),
        h("a", {class: "btn", href: "#/new?from=" + encodeURIComponent(r.folder), title: "New run with these settings (other material, NEP model or starting phase)"}, S.icon("plus"), "Copy to new run"),
        r.record ? h("button", {class: "btn ghost", onclick: () => checkRepro(r, repro)}, S.icon("check"), "Check reproducibility") : null));
  }

  async function checkRepro(r, repro) {
    const c = await S.api("/api/runs/check?folder=" + encodeURIComponent(r.folder));
    repro.innerHTML = "";
    if (!c.ok) return repro.append(h("div", {class: "notice bad"}, c.error));
    repro.append(h("div", {class: "notice " + (c.identical ? "ok" : "warn")}, c.identical ? "All inputs and scripts have the md5 recorded when the run was created." :
      [c.differences.length + " file(s) differ from the record:",
        h("ul", {style: {margin: "4px 0"}}, c.differences.map((x) => h("li", {class: "mono small"},
          `${x.file}: md5 ${(x.recorded || "-").slice(0, 10)}… → ${x.now ? x.now.slice(0, 10) + "…" : "missing"}`))),
        c.differences.some((x) => x.file === "campaign.toml") ? h("div", {class: "small"}, "campaign.toml changed: compare it with the settings in run_record.json (the 'config' entry, and the history of settings changes).") : null]),
      Object.keys(c.version_differences).length ? h("div", {class: "notice warn small"}, "package versions differ: " +
        Object.entries(c.version_differences).map(([k, v]) => `${k} ${v[0]} → ${v[1]}`).join(", ")) : null);
    repro.scrollIntoView({block: "nearest"});
  }

  async function runSwitcher(btn, cur) {
    document.querySelectorAll(".runpop").forEach((x) => x.remove());
    const pop = h("div", {class: "card runpop", role: "listbox", style: {position: "absolute", zIndex: 35, padding: "6px", minWidth: "320px", maxHeight: "60vh", overflow: "auto"}});
    const rc = btn.getBoundingClientRect();
    Object.assign(pop.style, {left: rc.left + window.scrollX + "px", top: rc.bottom + window.scrollY + 4 + "px"});
    document.body.append(pop);
    const close = (e) => { if (!pop.contains(e.target) && e.target !== btn) { pop.remove(); document.removeEventListener("click", close, true); } };
    setTimeout(() => document.addEventListener("click", close, true), 0);
    const recent = (S.info.prefs || {}).recent_runs || [];
    const l = await S.api("/api/runs/list", undefined, {quiet: true});
    const seen = new Set();
    const items = [];
    for (const c of recent) if (!seen.has(c)) { seen.add(c); items.push({config: c, name: S.basename(S.dirname(c)), folder: S.dirname(c)}); }
    for (const x of (l.runs || [])) if (x.config && !seen.has(x.config)) { seen.add(x.config); items.push(x); }
    pop.append(...items.slice(0, 30).map((x) => h("a", {href: "#/pipeline?config=" + encodeURIComponent(x.config), role: "option",
      class: "fl", style: {textDecoration: "none", color: "var(--text)", fontWeight: x.config === cur.config_path ? 700 : 500}},
      S.icon(x.config === cur.config_path ? "check" : "folder", "sm"), h("span", {}, x.name || S.basename(x.folder)), h("span", {class: "sz mono"}, x.starting_phase || ""))),
      h("div", {style: {borderTop: "1px solid var(--border)", marginTop: "4px", paddingTop: "4px"}},
        h("a", {class: "fl", href: "#/pipeline?list=1", style: {textDecoration: "none"}}, S.icon("folder", "sm"), "All runs / open a folder…"),
        h("a", {class: "fl", href: "#/new", style: {textDecoration: "none"}}, S.icon("plus", "sm"), "New run")));
    const first = pop.querySelector("a"); if (first) first.focus();
  }

  const COL = {clean: "var(--ok)", flagged: "var(--warn)", failed: "var(--bad)", not_run: "var(--idle)"};
  const LAB = {clean: "Clean", flagged: "Flagged", failed: "Failed", not_run: "Not run"};
  const TIP = {clean: "free energy done, no imaginary modes, relaxation converged",
    flagged: "free energy done but imaginary modes and/or relaxation not converged: not reliable until rerun",
    failed: "relaxation or phonons failed (status FAILED in the result CSV)", not_run: "prepared, no result yet"};
  function drawStrip(strip, st, r) {
    const c = st.counts, n = c.structures || 0;
    const pct = (v) => (n ? (100 * v / n).toFixed(1) : "0.0") + "%";
    const left = strip.querySelector(".left"), next = strip.querySelector(".next");
    const tm = st.timing || {};
    left.innerHTML = "";
    left.append(h("div", {class: "st-kpis"},
      h("div", {class: "st-kpi total"}, h("div", {class: "kpi-l"}, "Structures"), h("div", {class: "kpi-v num"}, S.int(n))),
      h("div", {class: "kpi-sep"}),
      ["clean", "flagged", "failed", "not_run"].map((k) => h("div", {class: "st-kpi", title: TIP[k]},
        h("div", {class: "kpi-l"}, h("i", {style: {background: COL[k]}}), LAB[k]),
        h("div", {class: "kpi-v num"}, S.int(c[k]), h("small", {}, pct(c[k]))))),
      h("div", {class: "kpi-sep"}),
      h("div", {class: "st-kpi", title: "chunks whose structures are all done (ok)"}, h("div", {class: "kpi-l"}, "Chunks done"),
        h("div", {class: "kpi-v num"}, st.chunks ? S.int(st.chunks.complete) : "–", h("small", {}, st.chunks ? "/ " + S.int(st.chunks.n) : "no plan"))),
      h("div", {class: "st-kpi", title: "sum of time_s × threads of the computed rows; planned = chunk cost model × sec_per_structure of the plan"},
        h("div", {class: "kpi-l"}, "Core-hours used"),
        h("div", {class: "kpi-v num"}, S.int(tm.used_core_h), h("small", {}, tm.planned_core_h != null ? "/ " + S.int(tm.planned_core_h) + " planned" : "")))));
    left.append(h("div", {class: "stack", role: "img", "aria-label": `Outcome of ${n} structures: ${["clean", "flagged", "failed", "not_run"].map((k) => c[k] + " " + LAB[k].toLowerCase()).join(", ")}`},
      ["clean", "flagged", "failed", "not_run"].filter((k) => c[k]).map((k) => h("div", {style: {width: (100 * c[k] / Math.max(1, n)) + "%", background: COL[k]}, title: `${LAB[k]}: ${c[k]}`}))),
      h("div", {class: "stack-axis"}, h("span", {}, "0"), h("span", {}, n ? `${(100 * (c.ok + c.failed) / n).toFixed(1)} % computed` : "nothing prepared yet"),
        h("span", {}, S.int(n) + " structures")));
    next.innerHTML = "";
    const nx = st.next || {};
    const nOpen = (st.attention || []).length;
    next.append(h("div", {class: "eyebrow"}, S.icon(nx.att ? "target" : "check", "sm"), "Next step"),
      h("div", {class: "what"}, nx.what || ""), h("div", {class: "why"}, nx.why || ""),
      h("div", {class: "btns"}, nx.button ? h("button", {class: "btn primary sm", "data-next": "1", onclick: () => S.runAction(nx.action, st)}, nx.button, S.icon("arrow", "sm")) : null,
        nOpen > 1 ? h("button", {class: "btn ghost sm", onclick: () => S.emit("show-attention")}, `See all ${nOpen}`) : null));
  }

  // what a "Next step" / attention button does: open a stage, go to a page, or filter Results
  S.runAction = function (act, st) {
    if (!act) return;
    if (act.go) return S.go(act.go, act.go === "results" ? {path: st.workdir} : undefined);
    if (act.results) return S.go("results", {path: st.workdir, ...act.results});
    if (act.panel) return S.emit("open-stage", act);
  };

  const PANEL_MOD = {check: "dry_run", prepare: "prepare", plan: "plan", compute: "compute", status: "status", merge: "merge", rerun: "rerun_flagged"};
  function stageCard(m, run, jobs, num) {
    const card = h("div", {class: "card", style: {padding: "10px 12px"}, "data-module": m.id});
    S.on("open-stage", (act) => {
      if (PANEL_MOD[act.panel] !== m.id) return;
      card.scrollIntoView({block: "center", behavior: "smooth"});
      card.style.boxShadow = "0 0 0 3px var(--accent-weak)"; card.style.borderColor = "var(--accent)";
      setTimeout(() => { card.style.boxShadow = ""; card.style.borderColor = ""; }, 2500);
      for (const [k, v] of Object.entries(act.preset || {})) {
        const el = card.querySelector(`[data-flag="${k === "choice" ? v : k}"]`);
        if (el && el.type === "checkbox") el.checked = true;
      }
    });
    const optEls = {};
    const opts = h("div", {});
    for (const o of m.options || []) {
      const key = o.flag || "__free__";
      const el = o.value ? h("input", {type: "text", class: "grow mono", placeholder: o.flag || "arguments", "data-flag": o.flag || "__free__"}) : h("input", {type: "checkbox", "data-flag": o.flag});
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
      jobs.prepend(S.jobView(r.job, {title: `${m.title} · ${md.label || mode}`, onDone: () => S.emit("stage-started")}));
      S.emit("stage-started");
    }
    return card;
  }
})();
