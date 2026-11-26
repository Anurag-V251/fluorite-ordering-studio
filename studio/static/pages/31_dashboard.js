/* Run dashboard pieces (used by the Pipeline page): the 7-stage pipeline, the Attention list, the jobs panel
   (this run's jobs here + its SLURM tasks), cluster settings, the stage side panel and the log drawer. */
(function () {
  const S = Studio, h = S.h;
  const D = S.dash = {};
  const STAGE_MOD = {check: "dry_run", prepare: "prepare", plan: "plan", compute: "compute", status: "status", merge: "merge", rerun: "rerun_flagged"};
  D.STAGE_MOD = STAGE_MOD;
  const MODE_ORDER = ["dry", "local", "slurm", "slurm_write"];
  const TIPS = {
    check: "campaign.py dry-run: reads the config and the inputs, counts the structures and estimates the core-hours. Writes nothing.",
    prepare: "Random orderings per composition (enumerate_substitutions.py), mapped onto every target phase (map_fluorite.py) and verified (compare_ordering.py must say IDENTICAL).",
    plan: "Packs the structures into chunks; one chunk = one node for one walltime. Uses sec_per_structure.",
    compute: "NEP relaxation + harmonic phonons (batch_free_energy.py) for every chunk. Finished structures are skipped on a restart.",
    status: "Chunk progress and the measured seconds per structure (compare with sec_per_structure).",
    merge: "One CSV of the newest row per structure: what Results, Plots and XRD read.",
    rerun: "Rows with imaginary modes or an unconverged relaxation, computed again with stricter settings.",
  };

  const sec = (v) => (v == null || !isFinite(v) ? "–" : (v < 10 ? Math.round(v * 10) / 10 : Math.round(v)) + " s");
  D.sec = sec;

  // ------------------------------------------------------------------ pipeline
  D.pipeline = function (st, active, onOpen) {
    const card = h("div", {class: "card flush"});
    card.append(h("div", {class: "card-h"}, h("h2", {}, "Pipeline"),
      h("span", {class: "sub"}, "Click a stage for its options, where to run it and the equivalent command")));
    const row = h("div", {class: "pipeline", role: "list"});
    st.stages.forEach((s, i) => {
      const seg = s.segments;
      const segbar = seg && seg.n ? h("div", {class: "segbar", role: "img",
        "aria-label": `${seg.done} done, ${seg.running} running, ${seg.unfinished - seg.running} not done of ${seg.n} chunks`},
        h("div", {style: {width: (100 * seg.done / seg.n) + "%", background: "var(--ok)"}}),
        h("div", {style: {width: (100 * seg.running / seg.n) + "%", background: "var(--run)"}}),
        seg.timeout ? h("div", {style: {width: (100 * seg.timeout / seg.n) + "%", background: "var(--bad)"}}) : null) : null;
      const btn = h("button", {class: "sbtn", "data-stage": s.key, title: TIPS[s.key] || s.title,
        "aria-label": `Stage ${i + 1}: ${s.title}, ${s.pill[1]}`, onclick: () => onOpen(s.key)},
        h("div", {class: "sn"}, h("b", {}, String(i + 1)), s.title, h("span", {class: "sub-l"}, s.sub)),
        S.pill(s.pill[0], s.pill[1]),
        h("div", {class: "sm"}, (s.lines || []).filter(Boolean).map((l) => h("div", {}, l))),
        segbar, h("span", {class: "open"}, S.icon("chev", "sm")));
      row.append(h("div", {class: "stage " + s.state + (active === s.key ? " is-active" : ""), role: "listitem"},
        i ? h("span", {class: "conn-line"}) : null, btn));
    });
    card.append(row);
    return card;
  };

  // ------------------------------------------------------------------ attention
  D.attention = function (st) {
    const att = st.attention || [];
    const card = h("div", {class: "card", id: "attention"});
    card.append(h("div", {class: "row", style: {marginBottom: "10px"}}, h("h2", {style: {margin: 0}}, "Attention"),
      att.length ? S.pill(att.some((a) => a.sev === "bad") ? "bad" : att.some((a) => a.sev === "warn") ? "warn" : "acc", `${att.length} open`) : S.pill("ok", "all clear")));
    if (!att.length) {
      card.append(h("div", {class: "allclear"}, S.icon("check"), h("div", {}, h("b", {}, "Nothing needs you. "),
        st.running && st.running.length ? "Jobs are running; this list updates by itself." : "Use the pipeline above for the next stage.")));
      return card;
    }
    card.append(h("ol", {class: "att"}, att.map((a, i) => h("li", {class: a.sev, "data-att": a.id},
      h("span", {class: "ord", "aria-hidden": "true"}, String(i + 1)),
      h("div", {}, h("div", {class: "ttl"}, a.title), h("div", {class: "dsc"}, a.desc),
        a.blocked ? h("div", {class: "blk"}, S.icon("alert", "sm"), " ", a.blocked) : null),
      h("div", {class: "act"}, h("button", {class: "btn sm" + (a.primary && !a.blocked ? " primary" : ""), disabled: a.blocked ? true : null,
        title: a.blocked || null, onclick: () => S.runAction(a.action, st)}, a.button))))));
    card.append(h("div", {class: "att-foot"}, S.icon("help", "sm"),
      "In the order to do them: fix the time estimate before resubmitting, rerun flagged rows after the compute ends, merge last."));
    return card;
  };

  // ------------------------------------------------------------------ cluster settings (editable after creation)
  D.settings = function (st, run, refresh) {
    const cl = st.cluster || {}, tm = st.timing || {};
    const card = h("div", {class: "card"});
    const meas = tm.measured;
    const fact = (l, v, cls, tip) => h("div", {class: "fact", title: tip || null}, h("div", {class: "l"}, l), h("div", {class: "v num " + (cls || "")}, v));
    card.append(h("div", {class: "row", style: {marginBottom: "10px"}}, h("h2", {style: {margin: 0}}, "Cluster settings"),
      h("span", {class: "grow"}),
      h("button", {class: "btn sm", onclick: () => D.editSettings(st, run, refresh)}, S.icon("edit", "sm"), "Edit")));
    card.append(h("div", {class: "facts"},
      fact("cores / node", S.int(cl.cores_per_node)), fact("walltime", cl.walltime || "–"),
      fact("partition", cl.partition || "default"),
      fact("sec_per_structure", (cl.sec_per_structure ?? "–") + " s", "", "what the next plan uses"),
      fact("planned with", tm.sps_plan != null ? tm.sps_plan + " s" : "no plan", "", "what the current chunks were packed with"),
      fact("measured", meas ? sec(meas.median) : "–", tm.off ? "warn" : meas ? "ok" : "",
        meas ? `median of ${meas.n} structures, 10–90 %: ${Math.round(meas.p10)}–${Math.round(meas.p90)} s` : "after the first finished structures")));
    if (meas && cl.sec_per_structure && Math.abs(meas.median / cl.sec_per_structure - 1) > 0.1) {
      const v = meas.median < 10 ? Math.round(meas.median * 100) / 100 : Math.round(meas.median * 10) / 10;
      card.append(h("div", {class: "callout info", style: {marginTop: "10px"}}, S.icon("pulse"),
        h("div", {class: "grow"}, `Measured ${sec(meas.median)} per structure, the config says ${cl.sec_per_structure} s.`),
        h("button", {class: "btn sm primary", "data-apply": "1", onclick: () => D.applyMeasured(st, v, refresh)}, `Use ${v} s`)));
    }
    return card;
  };

  D.applyMeasured = async function (st, v, refresh) {
    const ok = await S.confirm({title: `Set sec_per_structure to ${v} s?`, icon: "pulse", ok: "Change setting",
      body: `campaign.toml is changed (the old file is kept as a .bak copy and the change goes into the run's history). ` +
        `Chunks already planned do not change: re-pack them afterwards (Plan → Re-pack unfinished).`});
    if (!ok) return false;
    const r = await S.api("/api/runs/settings", {config: st.config, changes: {cluster: {sec_per_structure: v}}, note: "measured value applied"});
    if (!r.ok) return false;
    S.toast(Object.keys(r.changed || {}).length ? `sec_per_structure ${v} s saved (backup ${r.backup})` : "already set");
    if (refresh) refresh();
    return true;
  };

  const EDIT = [["cores_per_node", "cores per node", "what one job really uses"], ["walltime", "walltime", "HH:MM:SS or D-HH:MM:SS per compute task"],
    ["fill", "fill", "fraction of the walltime a chunk may use (0–1]"], ["sec_per_structure", "sec_per_structure", "seconds × threads per structure of the reference size"],
    ["max_array", "max_array", "MaxArraySize of the cluster"], ["max_concurrent", "max_concurrent", "array tasks running at once (%N)"],
    ["prepare_walltime", "prepare walltime", ""], ["merge_walltime", "merge walltime", ""],
    ["partition", "partition", ""], ["account", "account", ""], ["qos", "qos", ""], ["mem", "mem", "e.g. 200G; empty = cluster default"],
    ["python", "python", "python used inside the SLURM jobs"]];
  D.editSettings = function (st, run, refresh) {
    const cl = st.cluster || {};
    const F = {};
    const field = ([k, l, hint]) => { F[k] = h("input", {type: "text", class: "inp mono", value: cl[k] ?? "", "data-k": k});
      return h("label", {class: "field"}, h("span", {class: "flabel"}, l), F[k], hint ? h("span", {class: "fhint"}, hint) : null); };
    F.setup = h("textarea", {rows: 4, "data-k": "setup"}, cl.setup || "");
    F.sbatch_extra = h("textarea", {rows: 2, "data-k": "sbatch_extra", placeholder: "one per line"}, (cl.sbatch_extra || []).join("\n"));
    const note = h("input", {type: "text", class: "inp", placeholder: "why (goes into the history)"});
    const msg = h("div", {});
    const body = h("div", {class: "slide-b"},
      h("div", {class: "callout info"}, S.icon("help"), h("div", {}, "Only the cluster settings can change after a run is created. ",
        "The old campaign.toml is kept as a .bak copy and every change goes into the run's history. ",
        "For other settings (NEP model, phases, T) use ", h("a", {href: "#/new?from=" + encodeURIComponent(st.folder)}, "Copy to new run"), ".")),
      h("div", {class: "grp"}, h("div", {class: "grp-t"}, "Time and size"), h("div", {class: "row2"}, ...EDIT.slice(0, 8).map(field))),
      h("div", {class: "grp"}, h("div", {class: "grp-t"}, "Queue"), h("div", {class: "row2"}, ...EDIT.slice(8).map(field))),
      h("div", {class: "grp"}, h("div", {class: "grp-t"}, "Job script"),
        h("label", {class: "field"}, h("span", {class: "flabel"}, "setup ", h("em", {}, "shell lines at the start of every job")), F.setup),
        h("label", {class: "field"}, h("span", {class: "flabel"}, "sbatch extra"), F.sbatch_extra)),
      h("label", {class: "field"}, h("span", {class: "flabel"}, "note"), note), msg);
    const save = h("button", {class: "btn primary", onclick: async () => {
      const ch = {};
      for (const [k, el] of Object.entries(F)) {
        const old = k === "sbatch_extra" ? (cl.sbatch_extra || []).join("\n") : String(cl[k] ?? "");
        if (el.value !== old) ch[k] = k === "sbatch_extra" ? el.value.split("\n").map((x) => x.trim()).filter(Boolean) : el.value;
      }
      msg.innerHTML = "";
      if (!Object.keys(ch).length) { msg.append(h("div", {class: "notice warn"}, "Nothing changed.")); return; }
      save.disabled = true;
      const r = await S.api("/api/runs/settings", {config: st.config, changes: {cluster: ch}, note: note.value}, {quiet: true});
      save.disabled = false;
      if (!r.ok) { msg.append(h("div", {class: "notice bad"}, r.error)); return; }
      S.toast(`${Object.keys(r.changed).length} setting(s) changed · backup ${r.backup || "–"}`);
      D.panel.close(); if (refresh) refresh();
    }}, S.icon("check"), "Save settings");
    D.panel.open({eyebrow: "Run settings", title: "Cluster settings", desc: "Used by the next plan, SLURM script and submission.",
      body, foot: [h("span", {class: "sp"}), h("button", {class: "btn ghost", onclick: () => D.panel.close()}, "Cancel"), save]});
  };

  // ------------------------------------------------------------------ jobs of this run (here + SLURM)
  const SL_PILL = {RUNNING: "run", PENDING: "clock", SUBMITTED: "clock", CONFIGURING: "run", COMPLETING: "run", COMPLETED: "ok",
    ENDED: "ok", TIMEOUT: "bad", FAILED: "bad", CANCELLED: "warn", OUT_OF_MEMORY: "bad", NODE_FAIL: "bad", UNKNOWN: "idle"};
  const ACTIVE = new Set(["PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED", "RESIZING", "SUSPENDED", "SUBMITTED"]);
  D.jobs = function (st, refresh) {
    const card = h("div", {class: "card flush", id: "runjobs"});
    const sl = st.slurm || {}, tasks = sl.tasks || [];
    const sum = sl.summary || {};
    card.append(h("div", {class: "card-h"}, h("h2", {}, "Jobs"),
      h("span", {class: "sub"}, sl.available ? "squeue / sacct, refreshed with the page" : "SLURM not found here: states of submitted jobs are unknown"),
      h("div", {class: "right jobs-sum"}, Object.entries(sum).sort().map(([k, n]) => S.pill(SL_PILL[k] || "idle", `${n} ${k.toLowerCase()}`)))));
    const body = h("div", {class: "card-b"});
    card.append(body);
    // jobs started from the app for this run
    const mine = (st.jobs || []).slice(0, 8);
    if (!mine.length && !tasks.length) { body.append(h("div", {class: "muted small"}, "No jobs yet. Start a stage from the pipeline.")); return card; }
    if (mine.length) {
      body.append(h("div", {class: "tbl-wrap", style: {maxHeight: "240px", borderTop: 0}}, h("table", {class: "tbl"},
        h("thead", {}, h("tr", {}, ["started", "stage", "where", "state", "took", ""].map((x) => h("th", {}, x)))),
        h("tbody", {}, mine.map((j) => {
          const subs = (sl.submissions || []).filter((s) => s.app_job === j.id);
          let pill;
          if (subs.length && j.state === "done") {           // a submission: what matters is the SLURM state, not 'done'
            const tk = tasks.filter((t) => t.app_job === j.id);
            const act = tk.filter((t) => ACTIVE.has(t.state)).length;
            const n = subs.reduce((a, s) => a + s.n_tasks, 0);
            pill = S.pill(act ? "run" : "ok", `submitted ${n} task${n === 1 ? "" : "s"}` + (act ? ` · ${act} active` : " · ended"),
              "SLURM job " + subs.map((s) => s.id).join(", "));
          } else pill = S.pill({done: "ok", failed: "bad", lost: "warn", cancelled: "warn", running: "run", starting: "run"}[j.state] || "idle",
            j.state + (j.exit != null && j.exit !== 0 ? ` (exit ${j.exit})` : ""));
          return h("tr", {"data-job": j.id},
            h("td", {class: "num"}, new Date(j.start * 1000).toLocaleString([], {month: "short", day: "numeric", hour: "2-digit", minute: "2-digit"})),
            h("td", {}, (j.title || j.module || "").split(" [")[0]),
            h("td", {}, j.mode === "slurm" ? "SLURM" : j.mode === "slurm_write" ? "script" : "here"),
            h("td", {}, pill),
            h("td", {class: "num"}, j.end ? S.dur(j.end - j.start) : S.dur(Date.now() / 1000 - j.start)),
            h("td", {}, h("div", {class: "row", style: {gap: "4px", flexWrap: "nowrap"}},
              h("button", {class: "btn sm ghost", title: "Show its log in the drawer", onclick: () => D.drawer.follow(j.id)}, S.icon("term", "sm"), "Log"),
              ["running", "starting"].includes(j.state) ? h("button", {class: "btn sm danger", onclick: () => D.stopJob(j, refresh)}, S.icon("stop", "sm"), "Stop") : null)));
        })))));
    }
    // SLURM tasks the app submitted, grouped by submission
    for (const s of (sl.submissions || []).slice(0, 6)) {
      const tk = tasks.filter((t) => t.job === s.id);
      const act = tk.filter((t) => ACTIVE.has(t.state));
      body.append(h("div", {class: "row", style: {margin: "14px 0 6px"}},
        h("b", {}, `${s.kind} · job ${s.id}`), h("span", {class: "muted small"}, `${s.n_tasks} task${s.n_tasks === 1 ? "" : "s"}` +
          (s.array ? ` · array ${s.array}` : "") + (s.offset ? ` · chunks from ${s.offset}` : "") + (s.time ? " · " + S.ago(s.time) : "")),
        h("span", {class: "grow"}),
        act.length ? h("button", {class: "btn sm danger", onclick: () => D.scancel(s, act.length, refresh)}, S.icon("xcircle", "sm"), "scancel") : null));
      const LIM = 12;
      const rows = tk.slice(0, LIM).map((t) => {
        const prog = t.n ? h("div", {class: "bcell"}, h("span", {class: "pbar" + (t.done >= t.n ? " ok" : ""), title: `${t.done}/${t.n} structures`},
          h("i", {style: {width: (100 * Math.min(1, t.done / t.n)) + "%"}})), h("span", {class: "num"}, `${t.done}/${t.n}`)) : h("span", {class: "muted"}, "–");
        const over = t.projected_s && t.limit_s && t.projected_s > t.limit_s;
        const wall = t.limit_s ? h("div", {class: "bcell"}, h("span", {class: "pbar" + (over ? " bad" : ""), title: `${t.time || "0"} of ${t.limit}`},
          h("i", {style: {width: (100 * Math.min(1, (t.elapsed_s || 0) / t.limit_s)) + "%"}}),
          t.projected_s ? h("span", {class: "mark", style: {left: Math.min(99, 100 * t.projected_s / t.limit_s) + "%"}, title: "projected end"}) : null),
          h("span", {class: "num proj" + (over ? " bad" : "")}, over ? "won't fit" : t.time || "")) : h("span", {class: "muted small"}, t.time || "");
        return h("tr", {}, h("td", {class: "mono"}, t.task == null ? s.id : `${s.id}_${t.task}`), h("td", {class: "num"}, t.chunk ?? "–"),
          h("td", {}, S.pill(SL_PILL[t.state] || "idle", t.state.toLowerCase())), h("td", {}, prog), h("td", {}, wall), h("td", {class: "muted small"}, t.where || ""));
      });
      if (tk.length > LIM) rows.push(h("tr", {class: "fold-row"}, h("td", {colspan: 6}, `… ${tk.length - LIM} more tasks (` +
        Object.entries(tk.slice(LIM).reduce((a, t) => ((a[t.state] = (a[t.state] || 0) + 1), a), {})).map(([k, n]) => `${n} ${k.toLowerCase()}`).join(", ") + ")")));
      body.append(h("div", {class: "tbl-wrap"}, h("table", {class: "tbl"},
        h("thead", {}, h("tr", {}, ["task", "chunk", "state", "structures", "walltime", "node"].map((x) => h("th", {}, x)))), h("tbody", {}, rows))));
    }
    return card;
  };

  D.stopJob = async function (j, refresh) {
    const ok = await S.confirm({title: "Stop this job?", danger: true, ok: "Stop job", icon: "stop",
      body: `${j.title || j.module}: the process and its children get a stop signal. Finished structures are kept; a restart skips them.` +
        (j.mode === "slurm" ? " This does not cancel jobs it already submitted to SLURM (use scancel for those)." : "")});
    if (!ok) return;
    const r = await S.api("/api/jobs/cancel", {id: j.id});
    S.toast(r.ok ? "stop signal sent" : "could not stop (already finished?)");
    if (refresh) setTimeout(refresh, 1200);
  };

  D.scancel = async function (s, n, refresh) {
    const ok = await S.confirm({title: `scancel job ${s.id}?`, danger: true, ok: "scancel", icon: "xcircle",
      body: `${n} pending/running task(s) of this ${s.kind} submission are cancelled on the cluster. Structures finished so far are kept; ` +
        "running ones are lost and are computed again on the next submission."});
    if (!ok) return;
    const r = await S.api("/api/slurm/cancel", {ids: [String(s.id)]});
    S.toast("scancel: " + (r.output || "sent"));
    if (refresh) setTimeout(refresh, 1500);
  };

  // ------------------------------------------------------------------ slide-over panel (one at a time)
  D.panel = (function () {
    let scrim = null, slide = null, prev = null, onClose = null;
    const key = (e) => { if (e.key === "Escape" && slide && !document.querySelector(".dlg")) { e.preventDefault(); close(); } };
    function close() {
      if (!slide) return;
      scrim.classList.remove("on"); slide.classList.remove("on");
      const s = slide, c = scrim; setTimeout(() => { s.remove(); c.remove(); }, 230);
      slide = scrim = null;
      document.removeEventListener("keydown", key, true);
      if (onClose) { const f = onClose; onClose = null; f(); }
      if (prev && prev.focus && document.contains(prev)) prev.focus();
    }
    function open({eyebrow, title, desc, body, foot, onclose}) {
      if (slide) { const f = onClose; onClose = null; slide.remove(); scrim.remove(); slide = scrim = null; if (f) f(); }
      else prev = document.activeElement;
      onClose = onclose || null;
      scrim = h("div", {class: "scrim", onclick: close});
      const x = h("button", {class: "iconbtn", "aria-label": "Close panel", onclick: close, style: {marginLeft: "auto"}}, S.icon("close"));
      slide = h("aside", {class: "slide", role: "dialog", "aria-modal": "true", "aria-labelledby": "slideT"},
        h("div", {class: "slide-h"}, h("div", {style: {minWidth: 0}}, h("div", {class: "eyebrow"}, eyebrow || ""), h("h2", {id: "slideT"}, title),
          desc ? h("p", {}, desc) : null), x),
        body, h("div", {class: "slide-f"}, ...(foot || [])));
      document.body.append(scrim, slide);
      document.addEventListener("keydown", key, true);
      S.onLeave(() => { if (slide) close(); });
      requestAnimationFrame(() => { scrim.classList.add("on"); slide.classList.add("on");
        const f = slide.querySelector(".slide-b input:not([type=hidden]), .slide-b select, .slide-b button, .slide-f button.primary"); (f || x).focus(); });
      return slide;
    }
    return {open, close, isOpen: () => !!slide};
  })();

  // ------------------------------------------------------------------ one stage in the side panel
  D.openStage = function (key, st, run, preset, hooks) {
    const m = S.mods("stage").find((x) => x.id === (STAGE_MOD[key] || key));
    if (!m) return S.toast("no module for stage " + key, 4000, true);
    const sidx = st.stages.findIndex((s) => s.key === key);
    const stg = st.stages[sidx] || {};
    preset = Object.assign({}, preset || {});
    if (key === "plan" && !("choice" in preset)) preset.choice = st.plan ? "--append" : "";   // a plan exists: adding is the safe default
    const modes = Object.keys(m.modes || {}).sort((a, b) => (MODE_ORDER.indexOf(a) + 99) % 99 - (MODE_ORDER.indexOf(b) + 99) % 99);
    const saved = S.store.get("stage.mode." + m.id, null);
    let mode = preset.mode && modes.includes(preset.mode) ? preset.mode
      : saved && modes.includes(saved) && (saved !== "slurm" || S.info.slurm) ? saved
      : (preset["--unfinished"] && modes.includes(S.info.slurm ? "slurm" : "slurm_write")) ? (S.info.slurm ? "slurm" : "slurm_write")
      : modes.includes("local") ? "local" : modes[0];

    // facts at the top: what this stage is about, from the run state
    const c = st.counts || {}, tm = st.timing || {}, ch = st.chunks || {};
    const fact = (l, v, cls) => h("div", {class: "fact"}, h("div", {class: "l"}, l), h("div", {class: "v num " + (cls || "")}, v));
    const FACTS = {
      check: [fact("compositions", S.int(st.n_comp)), fact("phases", S.int((st.material.phases || []).length) || "–"), fact("T", st.material.T + " K")],
      prepare: [fact("prepared", `${st.prepared}/${st.n_comp}`, st.prepared >= st.n_comp ? "ok" : "warn"), fact("structures", S.int(c.structures)), fact("not in a chunk", S.int(c.unplanned), c.unplanned ? "warn" : "")],
      plan: [fact("chunks", st.plan ? S.int(st.plan.chunks) : "–"), fact("planned with", tm.sps_plan != null ? tm.sps_plan + " s" : "–"),
        fact("measured", sec(tm.measured && tm.measured.median), tm.off ? "warn" : "")],
      compute: [fact("chunks done", ch.n ? `${ch.complete}/${ch.n}` : "–", ch.n && ch.complete === ch.n ? "ok" : ""), fact("unfinished", S.int((ch.unfinished || []).length), (ch.unfinished || []).length ? "warn" : ""),
        fact("not run", S.int(c.not_run))],
      status: [fact("free energies", S.int(c.ok)), fact("failed", S.int(c.failed), c.failed ? "bad" : ""), fact("measured", sec(tm.measured && tm.measured.median), tm.off ? "warn" : "")],
      merge: [fact("merged rows", S.int(st.merge.rows)), fact("state", st.merge.exists ? (st.merge.stale ? "stale" : "up to date") : "not merged", st.merge.stale ? "warn" : st.merge.exists ? "ok" : ""), fact("rerun rounds", S.int(st.reruns))],
      rerun: [fact("flagged", S.int(c.flagged), c.flagged ? "warn" : "ok"), fact("imaginary", S.int(c.imaginary)), fact("unconverged", S.int(c.unconverged))],
    };
    const callouts = h("div", {class: "grp"});
    if (key === "plan" && preset.apply_measured) {
      callouts.append(h("div", {class: "callout"}, S.icon("pulse"), h("div", {class: "grow"},
        h("b", {}, `Measured ${preset.apply_measured} s per structure`), ` (the plan used ${tm.sps_plan} s, the config has ${st.cluster.sec_per_structure} s). `,
        "Put the measured value into the config first, then re-pack."),
        h("button", {class: "btn sm primary", "data-apply": "1", onclick: async () => {
          if (await D.applyMeasured(st, preset.apply_measured, hooks.refresh)) { await hooks.refresh(); D.openStage(key, S.runState || st, run, {choice: "--force"}, hooks); }
        }}, `Use ${preset.apply_measured} s`)));
    }
    if (key === "compute" && !S.info.slurm) callouts.append(h("div", {class: "callout info"}, S.icon("server"),
      h("div", {}, "SLURM is not found on this machine: run here, or write the script and submit it on the cluster with sbatch.")));
    if (!callouts.childElementCount) callouts.style.display = "none";
    const blockBox = h("div", {});

    // where to run
    const seg = h("div", {class: "seg", role: "radiogroup", "aria-label": "Where to run"});
    for (const md of modes) {
      const off = md === "slurm" && !S.info.slurm;
      const inp = h("input", {type: "radio", name: "stmode", value: md, disabled: off ? true : null});
      inp.checked = md === mode;
      inp.addEventListener("change", () => { mode = md; S.store.set("stage.mode." + m.id, md); drawOpts(); });
      seg.append(h("label", {title: off ? "SLURM not found on this machine (squeue is not on PATH)" : null}, inp,
        h("span", {}, S.icon({local: "play", slurm: "server", slurm_write: "download", dry: "search"}[md] || "play", "sm"), m.modes[md].label || md)));
    }

    // options of the chosen mode
    const vals = {};
    for (const o of m.options || []) {
      const k = o.key || o.flag || "__free__";
      const pk = o.kind === "choice" ? "choice" : k;
      vals[k] = pk in preset ? preset[pk] : o.default !== undefined ? o.default : o.kind === "choice" ? "" : o.value ? "" : false;
    }
    const optBox = h("div", {class: "grp"});
    const cli = h("div", {class: "cli", "aria-live": "polite"}, "…");
    const runBtn = h("button", {class: "btn primary", "data-run": "1"}, S.icon("play"), "Run");
    let last = null, tmr = null;
    function drawOpts() {
      optBox.innerHTML = "";
      const shown = (m.options || []).filter((o) => !o.modes || o.modes.includes(mode));
      if (shown.length) optBox.append(h("div", {class: "grp-t"}, "Options"));
      for (const o of shown) {
        const k = o.key || o.flag || "__free__";
        if (o.kind === "choice") {
          const box = h("div", {class: "field", role: "radiogroup", "aria-label": o.label || k, "data-opt": k},
            h("span", {class: "flabel"}, o.label || k));
          for (const chc of o.choices || []) {
            const r = h("input", {type: "radio", name: "opt_" + k, value: chc.flag || "", "data-flag": chc.flag || ""});
            r.checked = String(vals[k] || "") === String(chc.flag || "");
            r.addEventListener("change", () => { vals[k] = chc.flag || ""; preview(); });
            box.append(h("label", {class: "chk"}, r, h("span", {}, h("b", {}, chc.label || chc.flag || "default"),
              chc.flag ? h("span", {class: "mono muted"}, " " + chc.flag) : null, chc.help ? h("span", {class: "d"}, chc.help) : null)));
          }
          optBox.append(box);
        } else if (o.value && o.choices) {
          const sel = h("select", {class: "inp", "data-flag": o.flag}, o.choices.map((x) => h("option", {value: x}, x)));
          sel.value = vals[k] || o.default || o.choices[0];
          sel.addEventListener("change", () => { vals[k] = sel.value; preview(); });
          optBox.append(h("label", {class: "field"}, h("span", {class: "flabel"}, o.label || o.flag, h("em", {class: "mono"}, o.flag)), sel,
            o.help ? h("span", {class: "fhint"}, o.help) : null));
        } else if (o.value) {
          const inp = h("input", {type: "text", class: "inp mono", value: vals[k] || "", placeholder: o.flag || "arguments", "data-flag": o.flag || "__free__"});
          inp.addEventListener("input", () => { vals[k] = inp.value; preview(); });
          optBox.append(h("label", {class: "field"}, h("span", {class: "flabel"}, o.label || o.flag, o.flag ? h("em", {class: "mono"}, o.flag) : null), inp,
            o.help ? h("span", {class: "fhint"}, o.help) : null));
        } else {
          const cb = h("input", {type: "checkbox", "data-flag": o.flag});
          cb.checked = !!vals[k];
          cb.addEventListener("change", () => { vals[k] = cb.checked; preview(); });
          optBox.append(h("label", {class: "chk"}, cb, h("span", {}, o.label || o.flag, h("span", {class: "mono muted"}, " " + o.flag),
            o.help ? h("span", {class: "d"}, o.help) : null)));
        }
      }
      if (!shown.length) optBox.append(h("div", {class: "fhint"}, "No options for this mode."));
      runBtn.lastChild.textContent = m.modes[mode].label || "Run";
      preview();
    }
    function preview() {
      clearTimeout(tmr);
      tmr = setTimeout(async () => {
        const p = await S.api("/api/stage", {module: m.id, mode, config: st.config, options: vals, preview: true}, {quiet: true});
        last = p;
        cli.innerHTML = "";
        blockBox.innerHTML = "";
        if (!p.ok) { cli.append(h("span", {style: {color: "var(--bad)"}}, p.error)); runBtn.disabled = true; return; }
        cli.append(h("span", {class: "muted"}, "$ "), h("b", {}, p.cmd));
        runBtn.disabled = !!p.blocked;
        if (p.blocked) blockBox.append(h("div", {class: "callout bad"}, S.icon("stop"), h("div", {}, p.blocked)));
      }, 120);
    }
    const out = h("div", {class: "grp", "data-out": "1"});
    const lastJob = stg.job ? h("div", {class: "fhint"}, `Last run: ${stg.job.state}` + (stg.job.mode ? ` (${stg.job.mode})` : "") +
      ` · ${S.ago(stg.job.end || stg.job.start)} · `, h("a", {href: "javascript:void 0", onclick: () => D.drawer.follow(stg.job.id)}, "log")) : null;
    const body = h("div", {class: "slide-b"},
      h("div", {class: "facts"}, ...(FACTS[key] || [])), callouts, blockBox,
      h("div", {class: "grp"}, h("div", {class: "grp-t"}, "Where to run"), seg, lastJob),
      optBox, h("div", {class: "grp"}, h("div", {class: "grp-t"}, "Equivalent command"), cli,
        h("div", {class: "fhint"}, "Runs in the run folder; you can paste it into a terminal on the cluster.")), out);
    runBtn.onclick = async () => {
      if (runBtn.disabled) return;
      const p = last && last.ok ? last : await S.api("/api/stage", {module: m.id, mode, config: st.config, options: vals, preview: true});
      if (!p.ok) return;
      if (p.blocked) return S.toast(p.blocked, 6000, true);
      if ((p.confirm || []).length) {
        const ok = await S.confirm({title: m.modes[mode].label || "Run this stage?", ok: m.modes[mode].label || "Run",
          icon: mode === "slurm" ? "server" : "alert", body: h("div", {}, p.confirm.map((x) => h("p", {style: {margin: "0 0 8px"}}, x)),
            h("div", {class: "cli", style: {marginTop: "6px"}}, p.cmd))});
        if (!ok) return;
      }
      runBtn.disabled = true;
      let r;
      try { r = await S.api("/api/stage", {module: m.id, mode, config: st.config, options: vals}); }
      finally { setTimeout(() => { runBtn.disabled = false; }, 1500); }
      if (!r.ok) return;
      S.store.set("stage.mode." + m.id, mode);
      out.innerHTML = "";
      out.append(h("div", {class: "grp-t"}, "Output"), S.jobView(r.job, {title: `${m.title} · ${m.modes[mode].label || mode}`, onDone: () => hooks.refresh && hooks.refresh()}));
      out.scrollIntoView({block: "nearest", behavior: "smooth"});
      D.drawer.follow(r.job);
      S.emit("stage-started", {module: m.id, job: r.job});
    };
    D.panel.open({eyebrow: `Stage ${sidx + 1} of ${st.stages.length}` + (stg.pill ? " · " + stg.pill[1] : ""), title: m.title, desc: m.description,
      body, foot: [h("a", {class: "btn ghost sm", href: "#/modules", title: "Where this stage is defined (" + (m.source || "") + ")"}, S.icon("layers", "sm"), "Module"),
        h("span", {class: "sp"}), h("button", {class: "btn ghost", onclick: () => D.panel.close()}, "Close"), runBtn],
      onclose: () => { out.querySelectorAll(".jobcard").forEach((b) => b.stop && b.stop()); if (hooks.onclose) hooks.onclose(); }});
    drawOpts();
  };

  // ------------------------------------------------------------------ log drawer (bottom): the newest job of this run
  const LVL = (s0) => { const s = s0.replace(/\b0 (failed|errors?|FAILED)\b/g, ""); return LVL1(s); };
  const LVL1 = (s) => (/traceback|error|exception|\bfailed\b|FAILED|cannot|refus|not found/i.test(s) ? "ERROR"
    : /warn|flagged|imaginary|unconverged|timeout|stale|skipp/i.test(s) ? "WARN"
    : /\bok\b|\bdone\b|passed|written|submitted|complete|IDENTICAL/i.test(s) ? "OK" : "INFO");
  D.drawer = (function () {
    let el = null, jobs = [], jid = null, auto = true, lvl = "all", timer = null, lines = [], title = "";
    const stored = () => { try { return S.store.get("drawer.open", false); } catch (e) { return false; } };
    function mount() {
      if (el && document.body.contains(el)) return;
      const sel = h("select", {class: "inp", style: {height: "26px", width: "auto", maxWidth: "260px", fontSize: "12px"}, "aria-label": "Job shown in the log",
        onchange: () => { auto = sel.value === "auto"; jid = auto ? pickAuto() : sel.value; load(); }});
      const lv = h("div", {class: "lvl", role: "group", "aria-label": "Log level"}, ["all", "INFO", "WARN", "ERROR"].map((l) =>
        h("button", {"aria-pressed": String(l === lvl), "data-l": l, onclick: (e) => { lvl = l; el.querySelectorAll(".lvl button").forEach((b) => b.setAttribute("aria-pressed", String(b === e.currentTarget))); render(); }},
          l === "all" ? "All" : l[0] + l.slice(1).toLowerCase())));
      const tog = h("button", {class: "tog", "aria-expanded": "false", onclick: () => setOpen(!el.classList.contains("open"))},
        S.icon("term"), "Log", S.icon("chevup", "sm chev"));
      el = h("section", {class: "drawer", "aria-label": "Log of the newest job"},
        h("div", {class: "drawer-h"}, tog, h("span", {class: "last"}, ""), sel, lv,
          h("a", {class: "btn sm ghost", href: "#/jobs", title: "All jobs"}, "Jobs")),
        h("div", {class: "log", tabindex: "0"}));
      document.body.append(el);
      document.body.classList.add("has-drawer");
      setOpen(stored());
      S.onLeave(unmount);
    }
    function unmount() { if (el) el.remove(); el = null; clearInterval(timer); timer = null; document.body.classList.remove("has-drawer"); }
    function setOpen(o) {
      if (!el) return;
      el.classList.toggle("open", o);
      el.querySelector(".tog").setAttribute("aria-expanded", String(o));
      try { S.store.set("drawer.open", o); } catch (e) { /* private mode */ }
      if (o) { load(); const lg = el.querySelector(".log"); lg.scrollTop = lg.scrollHeight; }
    }
    function pickAuto() { const r = jobs.find((j) => ["running", "starting"].includes(j.state)); return (r || jobs[0] || {}).id || null; }
    function setJobs(list) {
      mount();
      jobs = list || [];
      if (auto) { const a = pickAuto(); if (a !== jid) { jid = a; lines = []; } }
      const sel = el.querySelector("select");
      sel.innerHTML = "";
      sel.append(h("option", {value: "auto"}, "newest job (follows)"), ...jobs.slice(0, 20).map((j) =>
        h("option", {value: j.id}, `${(j.title || j.module || j.id).split(" [")[0]} · ${j.state}`)));
      sel.value = auto ? "auto" : jid;
      if (!timer) timer = setInterval(load, 2500);
      load();
    }
    function follow(id) { mount(); auto = false; jid = id; lines = [];
      const sel = el.querySelector("select");
      if (![...sel.options].some((o) => o.value === id)) sel.append(h("option", {value: id}, id));
      sel.value = id; setOpen(true); }
    async function load() {
      if (!el || !jid) { if (el) el.querySelector(".last").textContent = "no job yet"; return; }
      const open = el.classList.contains("open");
      const r = await S.api("/api/jobs?id=" + jid + "&tail=" + (open ? 40000 : 2000), undefined, {quiet: true});
      if (!r.ok || !el) return;
      const j = r.job;
      title = (j.title || "").split(" [")[0];
      lines = (j.log || "").split("\n").filter((x) => x.trim());
      const lastL = lines[lines.length - 1] || "";
      el.querySelector(".last").textContent = `${title} · ${j.state} · ${lastL}`.slice(0, 300);
      if (open) render();
      if (["done", "failed", "lost", "cancelled"].includes(j.state) && !auto && !open) { /* keep showing it */ }
    }
    function render() {
      const lg = el.querySelector(".log");
      const atBottom = lg.scrollTop + lg.clientHeight >= lg.scrollHeight - 30;
      lg.innerHTML = "";
      const show = lines.map((t, i) => [i + 1, LVL(t), t]).filter((x) => lvl === "all" || x[1] === lvl).slice(-600);
      if (!show.length) lg.append(h("div", {class: "muted"}, h("span", {}), h("span", {}), h("span", {}), lines.length ? "no lines at this level" : "empty log"));
      for (const [n, L, t] of show) lg.append(h("div", {}, h("span", {class: "ts num"}, String(n)), h("span", {class: L}, L), h("span", {class: "src"}, title), h("span", {}, t)));
      if (atBottom) lg.scrollTop = lg.scrollHeight;
    }
    return {setJobs, follow, setOpen, unmount};
  })();
})();
