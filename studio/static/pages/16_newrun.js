/* New run: one page in sections with a live estimate on the right (structures, chunks, core-hours, walltime of a
   chunk).  Writes a reproducible run folder (campaign.toml + run_record.json) through /api/runs/create. */
(function () {
  const S = Studio, h = S.h;
  const splitTpl = (t) => { const i = t.indexOf("="); const l = i > 0 ? t.slice(0, i) : "";
    return /^[\w.+\-]+$/.test(l) ? [l, t.slice(i + 1)] : ["", t]; };
  const fmtv = (v, k) => (Array.isArray(v) ? v.join(k === "supercell" ? " " : ", ") : v == null ? "" : String(v));
  const term = (word, tip) => h("span", {class: "term", title: tip, tabindex: "0"}, word);
  const GLOSS = {
    sps: "sec_per_structure: core-seconds one host-sized structure needs (relaxation + phonons) with one thread. A guess until the first chunk is measured; the dashboard offers the measured value later.",
    fill: "fill: fraction of the walltime a chunk is planned to use (the rest is a safety margin).",
    orderings: "Distinct random arrangements of dopant and host cations per composition (all of them where fewer exist).",
    variants: "A cubic ordering can map onto a lower-symmetry phase in several equally good orientations (6 for t, up to 48 for oIII / m). They differ only in energy.",
    minlen: "Phonon supercell: repeated until every width is at least this length.",
    threads: "BLAS/OpenMP threads per structure; structures in parallel per node = cores per node / threads.",
    parent: "Also compute the orderings in the starting phase itself (un-mapped).",
  };

  S.newRun = async function (main, from) {
    const d = (await S.api("/api/runs/defaults")).defaults;
    let pre = null;
    if (from) { const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(from)); if (r.ok) pre = r; }
    const cfg = JSON.parse(JSON.stringify(d));
    cfg.structure.dopant = "";                       // never assume an element: the user names it
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
    const onCh = () => schedule();
    const inp = (sec, k, o = {}) => {
      const e = h("input", {type: "text", class: "inp" + (o.mono ? " mono" : ""), value: fmtv(cfg[sec][k], k), "data-k": sec + "." + k,
        placeholder: o.ph || "", "aria-required": o.req ? "true" : null});
      e.addEventListener("input", onCh); F[sec + "." + k] = e; return e;
    };
    const fld = (label, el, o = {}) => h("label", {class: "field"}, h("span", {class: "flabel"}, label, o.em ? h("em", {}, o.em) : null),
      o.unit ? h("div", {class: "inp-u"}, el, h("span", {class: "u"}, o.unit)) : el, o.hint ? h("span", {class: "fhint"}, o.hint) : null);
    const chk = (sec, k, label, hint) => { const c = h("input", {type: "checkbox", "data-k": sec + "." + k}); c.checked = !!cfg[sec][k];
      c.addEventListener("change", onCh); F[sec + "." + k] = c;
      return h("label", {class: "chk"}, c, h("span", {}, label, hint ? h("span", {class: "d"}, hint) : null)); };
    const adv = (n, ...kids) => h("details", {class: "adv"}, h("summary", {}, S.icon("chev", "sm"), `Advanced (${n})`), ...kids);

    // ---- 1 run
    const folder = S.pathField(pre ? pre.folder + "_copy" : S.info.start_dir + "/runs/run_" + new Date().toISOString().slice(0, 16).replace(/[-:]/g, "").replace("T", "_"), {dir: true, title: "Run folder (created if missing)"});
    folder.input.dataset.k = "run.folder";
    const name = h("input", {type: "text", class: "inp", value: pre && pre.record ? pre.record.name : "", "data-k": "run.name", placeholder: "e.g. ZrCeO2 2x2x2 from cubic"});
    const notes = h("textarea", {rows: 2, "data-k": "run.notes"}, pre && pre.record ? pre.record.notes || "" : "");
    const startPhase = h("input", {type: "text", class: "inp", value: pre && pre.record ? pre.record.starting_phase || "c" : "c", "data-k": "run.starting_phase"});
    // ---- 2 material
    const poscar = S.pathField(cfg.structure.poscar, {title: "Host supercell POSCAR (starting phase)"});
    poscar.input.dataset.k = "structure.poscar";
    const peek = h("div", {class: "fhint"});
    const host = inp("structure", "host", {req: true}), dop = inp("structure", "dopant", {req: true, ph: "e.g. Ce"});
    const comps = inp("structure", "compositions", {req: true, mono: true});
    const compList = h("div", {class: "fhint"});
    // ---- 3 phases
    const tplBox = h("div", {class: "tpls"});
    const addTpl = (label, path) => {
      const l = h("input", {type: "text", class: "inp", value: label || "", placeholder: "label (t, oIII, m…)", style: {width: "120px", flex: "none"}, "data-tpl": "label"});
      const p = S.pathField(path || "", {title: "12-atom template POSCAR of the target phase"});
      p.input.dataset.tpl = "path";
      [l, p.input].forEach((e) => { e.addEventListener("input", onCh); e.addEventListener("change", onCh); });
      const r = h("div", {class: "row tplrow"}, l, p, h("button", {class: "btn sm ghost", "aria-label": "Remove phase", onclick: () => { r.remove(); onCh(); }}, S.icon("close", "sm")));
      r.get = () => (p.input.value ? (l.value || S.basename(p.input.value).replace(/^POSCAR[_.-]?/, "")) + "=" + p.input.value : null);
      tplBox.append(r);
    };
    for (const t of cfg.mapping.templates || []) { const [l, p] = splitTpl(t); addTpl(l, p); }
    if (!(cfg.mapping.templates || []).length) addTpl("t", "");
    const fromRefs = async () => {
      const dir = await S.pick({dir: true, title: "Folder with POSCAR_<phase> templates", start: S.info.app_dir ? S.info.app_dir + "/refs" : undefined});
      if (!dir) return;
      const r = await S.api("/api/ls?path=" + encodeURIComponent(dir));
      tplBox.querySelectorAll(".tplrow").forEach((x) => { if (!x.get()) x.remove(); });
      for (const e of r.entries || []) if (!e.dir && /^POSCAR_/.test(e.name)) {
        const lab = e.name.replace(/^POSCAR_/, "");
        if (lab === startPhase.value) continue;           // the starting phase itself is the parent
        addTpl(lab, dir + "/" + e.name);
      }
      onCh();
    };
    const allVar = h("input", {type: "checkbox", "data-k": "mapping.all_variants"});
    allVar.checked = /--all-variants/.test(cfg.mapping.map_args);
    allVar.addEventListener("change", onCh);
    // ---- 4 relaxation
    const nep = S.pathField(cfg.free_energy.nep, {title: "NEP model (nep.txt)"});
    nep.input.dataset.k = "free_energy.nep";
    const nepInfo = h("span", {class: "fhint"});
    // ---- 5 cluster
    const setup = h("textarea", {rows: 4, "data-k": "cluster.setup", placeholder: "module load …\nsource …/activate\nexport RAYON_NUM_THREADS=1"}, cfg.cluster.setup || "");
    const sbx = h("textarea", {rows: 2, "data-k": "cluster.sbatch_extra", placeholder: "one per line, e.g. --exclusive"}, (cfg.cluster.sbatch_extra || []).join("\n"));
    // ---- 6 review
    const copyIn = h("input", {type: "checkbox", checked: true}), freeze = h("input", {type: "checkbox", checked: true}), overwrite = h("input", {type: "checkbox"});
    const checks = h("ul", {class: "checks"});
    const msg = h("div", {});
    [poscar.input, nep.input, folder.input].forEach((e) => { e.addEventListener("change", onCh); e.addEventListener("input", onCh); });
    [setup, sbx, notes, name, startPhase].forEach((e) => e.addEventListener("input", onCh));

    const SECS = [
      ["s-run", "Run", "Where the run folder goes and what to call it."],
      ["s-mat", "Material", "The host supercell (= starting phase), the elements and the compositions."],
      ["s-ph", "Phases", "Target phases for the transformation: the same cation ordering, mapped onto each."],
      ["s-fe", "Relaxation & phonons", "NEP relaxation, then the harmonic free energy at one temperature (fixed cell, not QHA)."],
      ["s-cl", "Cluster", "How the work is packed into SLURM jobs. Also used when running here."],
      ["s-rev", "Review & create", "Nothing runs until you start a stage on the run's page."],
    ];
    const sec = (i, ...body) => h("div", {class: "card sec", id: SECS[i][0]},
      h("div", {class: "sec-h"}, h("span", {class: "n"}, String(i + 1)), h("div", {}, h("h2", {}, SECS[i][1]), h("p", {}, SECS[i][2]))),
      h("div", {class: "sec-b"}, ...body));
    const createBtn = h("button", {class: "btn primary", "data-create": "1", onclick: create}, S.icon("check"), "Create run");
    const form = h("div", {class: "nrform"},
      sec(0, fld("Run folder", folder, {hint: "created if missing; holds campaign.toml, run_record.json, inputs/, scripts/ and the workdir"}),
        h("div", {class: "row2"}, fld("Name", name), fld("Starting phase label", startPhase, {hint: "phase of the host POSCAR (records and plot legends)"})),
        fld("Notes", notes)),
      sec(1, fld("Host POSCAR", h("div", {}, poscar, peek), {hint: "VASP5 supercell, a diagonal repeat of the 12-atom cell, axes along the cube directions"}),
        h("div", {class: "row3"}, fld("Host element", host, {hint: "element on the substituted sites"}), fld("Dopant element", dop, {hint: "required; any element"}),
          fld(h("span", {}, "Max ", term("orderings", GLOSS.orderings)), inp("structure", "max_orderings"), {em: "per composition"})),
        fld("Compositions", comps, {em: "dopant counts n", hint: "all | start:stop[:step] (inclusive) | a list, e.g. 0, 6, 16, 32"}), compList,
        adv(1, h("div", {class: "row3"}, fld("Seed", inp("structure", "seed"), {hint: "composition n uses seed×1000+n"})))),
      sec(2, h("div", {class: "field"}, h("span", {class: "flabel"}, "Target templates ", h("em", {}, "12 atoms: 4 cations + 8 anions")), tplBox,
        h("div", {class: "row"}, h("button", {class: "btn sm", onclick: () => { addTpl(); onCh(); }}, S.icon("plus", "sm"), "Phase"),
          h("button", {class: "btn sm", onclick: fromRefs}, S.icon("folder", "sm"), "Add all from a folder…")),
        h("span", {class: "fhint"}, "Every ordering is mapped onto each phase with map_fluorite.py and verified with compare_ordering.py (must be IDENTICAL).")),
        h("label", {class: "chk"}, allVar, h("span", {}, "Compute ", term("all variants", GLOSS.variants), h("span", {class: "d"}, "more structures; they can only be ranked by energy"))),
        chk("mapping", "include_parent", h("span", {}, "Include the ", term("starting phase", GLOSS.parent)), "the un-mapped orderings are computed too"),
        adv(3, h("div", {class: "row2"}, fld("Other map_fluorite options", inp("mapping", "map_args", {mono: true})), fld("Mapping jobs", inp("mapping", "jobs"))),
          chk("mapping", "verify", "Verify every mapped structure (compare_ordering.py)"))),
      sec(3, fld("NEP model", h("div", {}, nep, nepInfo)),
        h("div", {class: "row3"}, fld("Temperature", inp("free_energy", "T"), {unit: "K"}), fld("fmax", inp("free_energy", "fmax"), {unit: "eV/Å", hint: "rows above it after `steps` are flagged"}),
          fld("Max steps", inp("free_energy", "steps"))),
        chk("free_energy", "relax", "Relax with NEP before the phonons"),
        adv(5, h("div", {class: "row3"}, fld(term("min_length", GLOSS.minlen), inp("free_energy", "min_length"), {unit: "Å"}),
          fld("Phonon supercell", inp("free_energy", "supercell"), {hint: "empty = automatic, or e.g. 1 1 1"}),
          fld("q mesh", inp("free_energy", "mesh"))),
          h("div", {class: "row3"}, fld(term("Threads", GLOSS.threads), inp("free_energy", "threads"), {em: "per structure"})),
          chk("free_energy", "fixed_cell", "Fixed cell", "relax positions only (c/a and angles stay the template values)"))),
      sec(4, h("div", {class: "row3"}, fld("Cores per node", inp("cluster", "cores_per_node"), {hint: "what one job really uses"}),
        fld("Walltime", inp("cluster", "walltime", {mono: true}), {hint: "per compute task, HH:MM:SS"}),
        fld(term("sec_per_structure", GLOSS.sps), inp("cluster", "sec_per_structure"), {unit: "s", hint: "a guess until measured"})),
        h("div", {class: "row3"}, fld("Partition", inp("cluster", "partition")), fld("Max concurrent tasks", inp("cluster", "max_concurrent")),
          fld("Python in the jobs", inp("cluster", "python", {mono: true}))),
        fld("Setup lines", setup, {hint: "run at the start of every job (module load, source …/activate, export …)"}),
        adv(8, h("div", {class: "row3"}, fld(term("fill", GLOSS.fill), inp("cluster", "fill")), fld("max_array", inp("cluster", "max_array"), {hint: "MaxArraySize of the cluster"}),
          fld("mem", inp("cluster", "mem"))),
          h("div", {class: "row3"}, fld("Account", inp("cluster", "account")), fld("QOS", inp("cluster", "qos")), fld("Prepare walltime", inp("cluster", "prepare_walltime", {mono: true}))),
          h("div", {class: "row3"}, fld("Merge walltime", inp("cluster", "merge_walltime", {mono: true}))),
          fld("sbatch extra", sbx))),
      sec(5, checks,
        h("label", {class: "chk"}, copyIn, h("span", {}, "Copy inputs into <run>/inputs", h("span", {class: "d"}, "POSCAR, templates and NEP (md5 recorded either way)"))),
        h("label", {class: "chk"}, freeze, h("span", {}, "Freeze the scripts into <run>/scripts", h("span", {class: "d"}, "later app updates do not change this run"))),
        h("label", {class: "chk"}, overwrite, h("span", {}, "Overwrite the settings of an existing run", h("span", {class: "d"}, "campaign.py refuses changed templates / map_args in a used workdir"))),
        msg, h("div", {class: "row"}, createBtn, h("span", {class: "fhint"}, "then: the run's page → Check (dry-run)"))));

    // ---- summary (live estimate)
    const eS = h("span", {class: "v num"}, "–"), eC = h("span", {class: "v num"}, "–"), eH = h("span", {class: "v num"}, "–");
    const eCu = h("span", {class: "u"}), eHu = h("span", {class: "u"});
    const formula = h("div", {class: "formula"});
    const gaugeB = h("i", {}), gaugeT = h("b", {class: "num"}, "–"), gaugeH = h("span", {class: "fhint"});
    const kv = h("dl", {class: "kv"});
    const sumProblems = h("ul", {class: "checks"});
    const summary = h("aside", {class: "sum", "aria-label": "Estimate"}, h("div", {class: "card flush"},
      h("div", {class: "card-h"}, h("h2", {}, "Estimate"), h("span", {class: "sub"}, "updates as you type")),
      h("div", {class: "card-b"},
        h("div", {class: "big3"}, h("div", {}, h("span", {class: "l"}, "Structures"), eS),
          h("div", {}, h("span", {class: "l"}, "Chunks"), eC, eCu), h("div", {}, h("span", {class: "l"}, "Core-hours"), eH, eHu)),
        formula,
        h("div", {class: "gauge"}, h("div", {class: "lbl"}, h("span", {}, "One chunk of the walltime"), gaugeT),
          h("span", {class: "pbar"}, gaugeB), gaugeH),
        kv, sumProblems,
        h("button", {class: "btn primary", onclick: create}, S.icon("check"), "Create run"))));
    const toc = h("nav", {class: "toc", "aria-label": "Form sections"}, SECS.map(([id, t], i) =>
      h("a", {href: "javascript:void 0", "data-sec": id, onclick: () => document.getElementById(id).scrollIntoView({behavior: "smooth", block: "start"})},
        h("span", {class: "n"}, String(i + 1)), t)));

    main.append(h("div", {class: "runhead"}, h("div", {},
      h("div", {class: "crumb"}, h("a", {href: "#/pipeline?list=1"}, "Runs"), S.icon("chev", "sm"), "New"),
      h("h1", {}, pre ? "New run from " + S.basename(pre.folder) : "New run"),
      h("div", {class: "runmeta"}, h("span", {}, "One run = one starting phase of one material. Six short sections; the estimate on the right updates as you type."))),
      h("div", {class: "actions"}, h("a", {class: "btn ghost", href: "#/pipeline?list=1"}, "Cancel"))),
      h("div", {class: "nr"}, toc, form, summary));

    function collect() {
      const c = JSON.parse(JSON.stringify(d));
      for (const [k, el] of Object.entries(F)) { const [s, key] = k.split("."); c[s][key] = el.type === "checkbox" ? el.checked : el.value; }
      c.structure.poscar = poscar.input.value.trim();
      c.free_energy.nep = nep.input.value.trim();
      c.mapping.templates = [...tplBox.querySelectorAll(".tplrow")].map((r) => r.get()).filter(Boolean);
      const ma = String(c.mapping.map_args || "").replace(/--all-variants/g, "").trim();
      c.mapping.map_args = allVar.checked ? ("--all-variants " + ma).trim() : ma;
      c.cluster.setup = setup.value;
      c.cluster.sbatch_extra = sbx.value;
      c.structure.compositions = String(c.structure.compositions).trim();
      return c;
    }

    let tmr = null, seq = 0, peekFor = null;
    function schedule() { clearTimeout(tmr); tmr = setTimeout(update, 300); }
    async function update() {
      const my = ++seq;
      const c = collect();
      if (peekFor !== c.structure.poscar) { peekFor = c.structure.poscar; await doPeek(); }
      const r = await S.api("/api/runs/estimate", {config: collect()}, {quiet: true});
      if (my !== seq || !r) return;
      draw(r, collect());
    }
    async function doPeek() {
      peek.innerHTML = "";
      const p = poscar.input.value.trim();
      if (!p) return;
      const r = await S.api("/api/poscar?path=" + encodeURIComponent(p), undefined, {quiet: true});
      if (!r.ok) { peek.className = "fhint bad"; peek.append(r.error); return; }
      peek.className = "fhint";
      const sp = Object.entries(r.species).map(([e, n]) => e + n).join(" ");
      const cations = Object.keys(r.species).filter((e) => !["O", "F", "S", "N", "Se", "Te", "Cl"].includes(e));
      if (!(host.value in r.species) && cations.length === 1) { host.value = cations[0]; S.toast(`host element set to ${cations[0]} (the cation of this POSCAR)`); }
      peek.append(h("span", {class: "mono"}, `${sp} · ${r.n_atoms} atoms · ${r.lengths.join(" × ")} Å` + (r.repeat_guess ? ` · ${r.repeat_guess.join("×")} of the 12-atom cell` : "")));
    }
    function draw(r, c) {
      const ok = r.ok !== false;
      const est = ok ? r : {};
      eS.textContent = est.structures != null ? S.int(est.structures) + (est.all_variants ? "+" : "") : "–";
      eC.textContent = est.chunks != null ? S.int(est.chunks) : "–";
      eCu.textContent = est.workers ? `× ${est.workers} in parallel` : "";
      eH.textContent = est.core_h != null ? S.int(est.core_h) : "–";
      eHu.textContent = est.node_h != null ? `${S.fmt(est.node_h, est.node_h < 10 ? 1 : 0)} node-h` : "";
      formula.textContent = est.structures != null ? `${S.int(est.structures)} × ${c.cluster.sec_per_structure} s × ${c.free_energy.threads} thread(s) / 3600` +
        (Object.values(est.units || {}).some((u) => u !== 1) ? " × phonon cell size" : "") : "";
      const wall = est.walltime_s || 0, ch = est.chunk_s || 0;
      const frac = wall ? ch / wall : 0;
      gaugeB.style.width = Math.min(100, 100 * frac) + "%";
      gaugeB.parentNode.className = "pbar" + (frac > 0.95 ? " bad" : frac > Number(c.cluster.fill || 0.8) + 0.01 ? " warn" : " ok");
      gaugeT.textContent = ch ? `${S.hms(ch)} / ${c.cluster.walltime}` : "–";
      gaugeH.textContent = !ch ? "" : est.largest_structure_s > wall ? "one structure alone needs more than the walltime: raise the walltime" :
        `planned to use ${Math.round(100 * frac)} % (fill ${c.cluster.fill}); a shorter walltime gives more, smaller chunks`;
      kv.innerHTML = "";
      const kvs = [];
      if (est.host) kvs.push(["host", `${est.host.n_atoms} atoms · ${est.host.n_host ?? "?"} ${c.structure.host} sites`]);
      if ((est.compositions || []).length) kvs.push(["compositions", String(est.compositions.length)]);
      for (const [ph, n] of Object.entries(est.per_phase || {})) kvs.push([ph === "parent" ? `${startPhase.value || "parent"} (starting)` : ph, S.int(n) + (est.all_variants && ph !== "parent" ? "+" : "")]);
      if (est.rounds) kvs.push(["wall time", `≈ ${S.hms(est.wall_total_s)} with ${c.cluster.max_concurrent} nodes at a time (+ queue)`]);
      kvs.forEach(([k, v]) => kv.append(h("dt", {}, k), h("dd", {class: "num"}, v)));
      compList.textContent = (est.compositions || []).length ? "n = " + est.compositions.slice(0, 40).map((x) => `${x.n} (${S.fmt(x.pct, 1)} %${x.orderings < Number(c.structure.max_orderings) ? `, all ${x.orderings}` : ""})`).join(", ") +
        (est.compositions.length > 40 ? " …" : "") : "";
      nepInfo.className = "fhint" + (est.nep ? (est.nep.found ? " ok" : " bad") : "");
      nepInfo.textContent = !c.free_energy.nep ? "required" : !est.nep ? "" : !est.nep.found ? "not found" :
        "✓ found" + (est.nep.species ? ` · species ${est.nep.species.join(", ")}` : "");
      // checks: what blocks Create, and the numbers to look at
      const probs = ok ? (est.problems || []) : [r.error];
      const items = [];
      for (const p of probs) items.push(["bad", p]);
      if (est.array_too_big) items.push(["warn", `${est.chunks} chunks > max_array ${c.cluster.max_array}: campaign.py splits them into several arrays`]);
      if (est.all_variants) items.push(["warn", "all variants: the real count is known after prepare"]);
      if (!probs.length) items.push(["ok", "settings pass the checks campaign.py dry-run makes"]);
      for (const box of [checks, sumProblems]) {
        box.innerHTML = "";
        items.forEach(([k, t]) => box.append(h("li", {}, S.icon(k === "ok" ? "check" : k === "warn" ? "alert" : "xcircle", k), h("span", {}, t))));
      }
      // TOC: a section with a problem gets a warning mark
      const secOf = (p) => /POSCAR|host|dopant|composition|orderings/i.test(p) ? "s-mat" : /template|phase label/i.test(p) ? "s-ph" :
        /NEP|fmax|min_length|T must|supercell|steps|threads/i.test(p) ? "s-fe" : /walltime|cores|workdir/i.test(p) ? "s-cl" : "s-rev";
      const bad = new Set(probs.map(secOf));
      toc.querySelectorAll("a").forEach((a) => { a.classList.toggle("wn", bad.has(a.dataset.sec)); a.classList.toggle("ok", !bad.has(a.dataset.sec) && a.dataset.sec !== "s-rev"); });
    }

    async function create(e) {
      const btn = e && e.currentTarget;
      if (create.busy) return;
      create.busy = true; if (btn) btn.disabled = true;     // no double submit
      let r;
      try {
        r = await S.api("/api/runs/create", {folder: folder.input.value, name: name.value, notes: notes.value,
          starting_phase: startPhase.value, copy_inputs: copyIn.checked, freeze_scripts: freeze.checked, overwrite: overwrite.checked, config: collect()});
      } finally { create.busy = false; if (btn) btn.disabled = false; }
      msg.innerHTML = "";
      if (!r.ok) { msg.append(h("div", {class: "callout bad", role: "alert"}, S.icon("xcircle"), h("div", {}, r.error)));
        msg.scrollIntoView({block: "center"}); return; }
      S.info.prefs.recent_runs = [r.config, ...((S.info.prefs || {}).recent_runs || [])];
      S.go("pipeline", {config: r.config});
    }
    // section highlight in the TOC while scrolling
    const io = new IntersectionObserver((es) => es.forEach((x) => { if (x.isIntersecting) toc.querySelectorAll("a").forEach((a) => a.classList.toggle("cur", a.dataset.sec === x.target.id)); }),
      {rootMargin: "-40% 0px -55% 0px"});
    form.querySelectorAll(".sec").forEach((s) => io.observe(s));
    S.onLeave(() => io.disconnect());
    update();
  };
})();
