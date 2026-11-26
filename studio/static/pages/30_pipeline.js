/* Pipeline page: create a reproducible run (campaign.toml + run_record.json) and drive its stages. */
(function () {
  const S = Studio, h = S.h;
  S.registerPage({id: "pipeline", title: "Pipeline", order: 10, icon: "flow", render,
    navAs: (p) => (p.new ? "new" : null)});
  S.registerPage({id: "new", title: "New run", order: 15, icon: "plus", render: (main, p) => S.newRun(main, p.from)});

  async function render(main, params) {
    if (params.config) return showRun(main, params.config);
    if (params.new) return S.newRun(main, params.from);
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

  // ------------------------------------------------------------------ one run: dashboard
  async function showRun(main, config) {
    const r = await S.api("/api/runs/read?folder=" + encodeURIComponent(config));
    if (!r.ok) { main.append(h("div", {class: "notice bad"}, r.error), h("a", {href: "#/pipeline?list=1"}, "all runs")); return; }
    const rec = r.record;
    S.currentRun = r;
    const D = S.dash;
    const head = h("div", {class: "runhead"});
    const strip = h("div", {class: "card strip"}, h("div", {class: "left"}, "reading the run…"), h("div", {class: "next"}));
    const warnBox = h("div", {});
    const repro = h("div", {});
    const pipe = h("div", {});
    const attBox = h("div", {}), setBox = h("div", {});
    const jobsBox = h("div", {});
    const hmBox = h("div", {});
    main.append(head, warnBox, strip, repro, pipe, h("div", {class: "grid g-main"}, hmBox, h("div", {}, attBox, setBox)), jobsBox);
    let hmAt = 0;
    const drawHm = async () => {                 // the scan is cached on the server; redraw at most every 30 s
      if (Date.now() - hmAt < 30000) return;
      hmAt = Date.now();
      const sc = await S.api("/api/scan?path=" + encodeURIComponent(r.config_path), undefined, {quiet: true});
      if (!sc.ok || sc.layout === "none") return;
      const T = st ? String(st.material.T) : (sc.temperatures || [])[0];
      const Tk = (sc.temperatures || []).find((t) => Number(t) === Number(T)) || (sc.temperatures || [])[0];
      hmBox.innerHTML = "";
      hmBox.append(S.heatmap(sc, Tk, {onPick: (comp, phase) => S.go("results", {path: r.folder, comp, phase})}));
    };
    let st = null, polled = 0, active = null;
    const pollAgo = h("span", {class: "muted num"});
    const openStage = (key, preset) => {
      if (!st) return;
      active = key; drawPipe();
      D.openStage(key, st, r, preset, {refresh: () => poll(true), onclose: () => { active = null; drawPipe(); }});
    };
    const drawPipe = () => { if (st) { pipe.innerHTML = ""; pipe.append(D.pipeline(st, active, (k) => openStage(k))); } };
    const draw = () => {
      drawHead(head, r, st, pollAgo, poll, repro);
      if (st) {
        drawStrip(strip, st, r); warnBox.innerHTML = "";
        (st.warnings || []).forEach((w) => warnBox.append(h("div", {class: "notice warn"}, w)));
        drawPipe();
        attBox.innerHTML = ""; attBox.append(D.attention(st));
        setBox.innerHTML = ""; setBox.append(D.settings(st, r, () => poll(true)));
        jobsBox.innerHTML = ""; jobsBox.append(D.jobs(st, () => poll(true)));
        D.drawer.setJobs(st.jobs || []);
        drawHm();
      }
      S.emit("run-state", st);
    };
    async function poll(fresh) {
      const s2 = await S.api("/api/run/state?config=" + encodeURIComponent(r.config_path) + (fresh ? "&fresh=1" : ""), undefined, {quiet: !fresh});
      if (!s2.ok) return;
      st = s2; polled = Date.now(); S.runState = st; draw();
    }
    S.every(1000, () => { pollAgo.textContent = polled ? Math.round((Date.now() - polled) / 1000) + " s" : ""; });
    S.every(5000, () => {                          // every 5 s while something runs, else every 15 s
      const busy = st && ((st.running || []).length || st.compute_running);
      if (busy || Date.now() - polled > 15000) poll(false);
    });
    S.on("open-stage", (act) => openStage(Object.keys(D.STAGE_MOD).find((k) => k === act.panel) || act.panel, act.preset));
    S.on("show-attention", () => { const a = document.getElementById("attention"); if (a) { a.scrollIntoView({block: "center", behavior: "smooth"}); const b = a.querySelector("button"); if (b) b.focus(); } });
    S.on("stage-started", () => { hmAt = 0; setTimeout(() => poll(true), 600); });
    draw();
    await poll(false);
    main.append(h("div", {}, rec && rec.history && rec.history.length ? h("details", {class: "card"}, h("summary", {}, `History (${rec.history.length})`),
        h("pre", {class: "log"}, rec.history.slice().reverse().map((x) => `${x.time}  ${x.stage || x.event || ""}` + (x.mode ? ` [${x.mode}]` : "") + `  ${x.cmd || ""}` +
          (x.changes ? Object.entries(x.changes).map(([k, v]) => `${k}: ${v[0]} -> ${v[1]}`).join(", ") + (x.note ? `  (${x.note})` : "") : "")).join("\n"))) : null,
      h("details", {class: "card"}, h("summary", {}, S.basename(r.config_path)), (() => { const p = h("pre", {class: "log"}); S.api("/api/text?path=" + encodeURIComponent(r.config_path)).then((t) => { p.textContent = t.text; }); return p; })())));
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
        !rec.created ? h("span", {class: "t2"}, S.icon("alert", "sm"), "no run_record.json (made without the app) ",
          h("button", {class: "btn sm", title: "Write run_record.json next to the config so settings changes and stages get a history (nothing else changes)",
            onclick: async () => { const a = await S.api("/api/runs/adopt", {config: r.config_path}); if (a.ok) { S.toast("run adopted: history kept from now on"); S.route(); } }}, "Keep a history")) : null)),
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

})();
