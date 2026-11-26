/* Plots and Tools pages: one generic runner for any module with a script (form built from its argparse). */
(function () {
  const S = Studio, h = S.h;

  const INPUT_LABEL = {merged_csv: "merged CSV (free_energy_<T>K_all.csv)", comp_csv: "per-composition CSV (Ce_xx/free_energy_<T>K.csv)",
    any_csv: "CSV", config: "campaign.toml", poscar: "structure (POSCAR)", folder: "folder", none: ""};
  const INPUT_FILTER = {merged_csv: "\\.csv$", comp_csv: "\\.csv$", any_csv: "\\.csv$", config: "\\.toml$"};

  S.moduleRunner = function (main, kind, params, title) {
    const mods = S.mods(kind);
    const left = h("div", {class: "list"});
    const right = h("div", {});
    main.append(h("div", {class: "cols"}, h("div", {}, h("div", {class: "card"}, h("h2", {}, title), left,
      h("div", {class: "small muted"}, "Add your own: drop a .toml in modules.d/ (see the Modules page).")), ), right));
    let cur = params.module && mods.find((m) => m.id === params.module) ? params.module : S.store.get(kind + ".module", mods[0] && mods[0].id);
    if (!mods.find((m) => m.id === cur) && mods.length) cur = mods[0].id;
    for (const m of mods) {
      const it = h("div", {class: "item" + (m.id === cur ? " on" : ""), onclick: () => { cur = m.id; S.store.set(kind + ".module", cur);
        left.querySelectorAll(".item").forEach((x) => x.classList.remove("on")); it.classList.add("on"); show(m); }},
        h("b", {}, m.title), m.script_found ? null : h("div", {class: "small", style: {color: "var(--bad)"}}, "script not found"),
        h("div", {class: "small muted"}, m.script ? S.basename(m.script) : ""));
      left.append(it);
    }
    const m0 = mods.find((m) => m.id === cur);
    if (m0) show(m0);

    async function show(m) {
      right.innerHTML = "";
      const card = h("div", {class: "card"}, h("h2", {}, m.title), m.description ? h("div", {class: "muted"}, m.description) : null);
      right.append(card);
      const missing = Object.entries(m.needs_ok || {}).filter(([, ok]) => !ok).map(([n]) => n);
      if (missing.length) card.append(h("div", {class: "notice warn small"}, "Python packages missing in the server's environment: " + missing.join(", ") + " (start the app from the environment that has them)."));
      if (!m.script_found) { card.append(h("div", {class: "notice bad"}, "Script not found: " + m.script_path)); return; }
      let inField = null;
      if (m.input && m.input !== "none") {
        let start = params.csv || params.input || S.store.get("input." + m.input, "");
        if (!start && S.lastScan && /csv/.test(m.input)) {
          const sc = S.lastScan;
          const T = S.store.get("results.T", null);
          const pick = sc.merged_csv.find((f) => T && f === `free_energy_${T}K_all.csv`) || sc.merged_csv[0];
          if (m.input === "merged_csv" && pick) start = sc.workdir + "/" + pick;
        }
        inField = S.pathField(start, {dir: m.input === "folder", filter: INPUT_FILTER[m.input], title: INPUT_LABEL[m.input]});
        card.append(h("h3", {}, "Input: " + INPUT_LABEL[m.input]), inField);
      }
      const specBox = h("div", {class: "muted"}, "reading the script's options…");
      const optsBox = h("details", {open: true}, h("summary", {}, "Options (from the script's --help)"), specBox);
      card.append(optsBox);
      const cmdLine = h("pre", {class: "log", style: {display: "none"}});
      const runBtn = h("button", {class: "primary", disabled: true, title: "reading the script's options…"}, "Run");
      const prevBtn = h("button", {disabled: true}, "Show command");
      card.append(h("div", {class: "row", style: {marginTop: "10px"}}, runBtn, prevBtn,
        h("button", {onclick: () => { S.store.set("args." + m.id, null); show(m); }}, "Reset options")), cmdLine);
      const jobs = h("div", {});
      card.append(jobs);
      const spec = await S.api("/api/argspec?module=" + m.id, undefined, {quiet: true});
      specBox.innerHTML = "";
      if (!spec.ok) {
        specBox.append(h("div", {class: "notice bad"}, "Could not read the options: " + spec.error + (spec.hint ? " — " + spec.hint : "")),
          spec.output ? h("pre", {class: "log"}, spec.output) : null);
      }
      const saved = S.store.get("args." + m.id, null);
      const form = S.argForm(spec.ok ? spec : {actions: []}, {skipFirstPositional: m.input && m.input !== "none",
        saved: saved ? saved.args : [], savedExtra: saved ? saved.extra : ""});
      if (spec.ok && spec.description) specBox.append(h("details", {}, h("summary", {}, "Script description"), h("pre", {class: "log"}, spec.description + "\n\n" + (spec.epilog || ""))));
      specBox.append(form);
      const argv = () => form.getArgs();
      const command = () => ["python3", m.script_path, ...(inField && inField.input.value ? [inField.input.value] : []), ...(m.args || []), ...argv()].map(S.quote).join(" ");
      runBtn.disabled = prevBtn.disabled = false; runBtn.title = "";
      prevBtn.onclick = () => { cmdLine.style.display = ""; cmdLine.textContent = command(); };
      runBtn.onclick = async () => {
        const input = inField ? inField.input.value.trim() : null;
        if (inField && !input) return S.toast("choose the input first");
        const bad = form.check();
        if (bad.length) { optsBox.open = true; return S.toast(bad.join("; ")); }
        if (inField) S.store.set("input." + m.input, input);
        S.store.set("args." + m.id, {args: form.getArgs(true), extra: form.getExtra()});
        const r = await S.api("/api/run", {module: m.id, input, args: argv()});
        if (!r.ok) return;
        cmdLine.style.display = ""; cmdLine.textContent = r.cmd;
        const v = S.jobView(r.job, {title: m.title + (input ? " · " + S.basename(input) : "")});
        jobs.prepend(v);
        v.scrollIntoView({block: "start", behavior: "smooth"});
      };
    }
  };

  S.registerPage({id: "plots", title: "Plots", order: 40, icon: "chart", render: (main, p) => S.moduleRunner(main, "plot", p, "Plots")});
  S.registerPage({id: "tools", title: "Tools", order: 50, icon: "tool", render: (main, p) => S.moduleRunner(main, "tool", p, "Tools")});
})();
