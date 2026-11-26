/* Results page: scan a folder and show what has been done. */
(function () {
  const S = Studio, h = S.h;
  S.registerPage({id: "results", title: "Results", order: 30, icon: "table", render});

  function render(main, params) {
    const start = params.path || S.store.get("results.path", (S.info.prefs || {}).last_results || "");
    const pf = S.pathField(start, {dir: true, title: "Results folder (campaign workdir, or the folder with campaign.toml)"});
    const auto = h("input", {type: "checkbox"});
    const out = h("div", {});
    const recent = ((S.info.prefs || {}).recent_results || []).slice(0, 6);
    const filt = {comp: params.comp || "", phase: params.phase || "", status: params.status || ""};
    const go = () => { S.store.set("results.path", pf.input.value); scan(pf.input.value, out, false, filt); };
    pf.input.addEventListener("keydown", (e) => { if (e.key === "Enter") go(); });
    main.append(h("div", {class: "card"},
      h("h2", {}, "Results folder"),
      h("div", {class: "row"}, pf, h("button", {class: "primary", onclick: go}, "Scan"),
        h("label", {class: "small"}, auto, " refresh every 30 s")),
      recent.length ? h("div", {class: "small muted", style: {marginTop: "6px"}}, "recent: ",
        recent.map((r) => h("a", {href: "javascript:void 0", title: r, style: {marginRight: "10px"}, onclick: () => { pf.input.value = r; go(); }},
          r.split("/").slice(-2).join("/")))) : null),
      out);
    if (!start) out.append(h("div", {class: "card"}, h("div", {class: "muted"},
      "Pick a folder above: a run folder (with campaign.toml), its workdir, or any folder with Ce_xx/ composition folders. ",
      "No run yet? ", h("a", {href: "#/pipeline?new=1"}, "Pipeline → New run"), "."),
      S.runsList((x) => { pf.input.value = x.folder; go(); })));
    S.every(30000, () => { if (auto.checked && pf.input.value) scan(pf.input.value, out, true, filt); });
    if (start) go();
  }

  async function scan(path, out, quiet, filt) {
    if (!quiet) out.innerHTML = '<div class="card muted">scanning…</div>';
    const r = await S.api("/api/scan?path=" + encodeURIComponent(path));
    out.innerHTML = "";
    if (!r.ok) { out.append(h("div", {class: "notice bad"}, r.error)); return; }
    for (const w of r.warnings || []) out.append(h("div", {class: "notice warn"}, w));
    if (r.layout === "none") {
      if ((r.candidates || []).length) out.append(h("div", {class: "card"}, h("h2", {}, "Pick one"),
        r.candidates.map((c) => h("div", {}, h("a", {href: "#/results?path=" + encodeURIComponent(c)}, c)))));
      return;
    }
    S.lastScan = r;
    const temps = r.temperatures.length ? r.temperatures : ["-"];
    let T = S.store.get("results.T", temps[0]);
    if (!temps.includes(T)) T = temps[temps.length - 1];
    const body = h("div", {});
    const tsel = h("select", {onchange: () => { S.store.set("results.T", tsel.value); draw(tsel.value); }}, temps.map((t) => h("option", {value: t}, t + " K")));
    tsel.value = T;
    out.append(h("div", {class: "card"},
      h("div", {class: "row"}, h("h2", {style: {margin: 0}}, r.run_name || (S.basename(r.workdir) === "campaign" ? S.basename(S.dirname(r.workdir)) : S.basename(r.workdir))),
        h("span", {class: "pill run"}, r.layout === "campaign" ? "campaign workdir" : "plain folders"),
        h("span", {class: "muted small mono grow"}, r.workdir),
        temps.length > 1 ? h("label", {}, "T ", tsel) : h("span", {class: "muted"}, temps[0] + " K"),
        r.config_path ? h("a", {href: "#/pipeline?config=" + encodeURIComponent(r.config_path)}, "open run") : null),
      h("div", {class: "small muted"}, `layout derived from campaign.py output · ${r.n_result_sources} result file(s) · scanned in ${r.scan_seconds} s`)),
      body);
    draw(T);

    function draw(T) {
      body.innerHTML = "";
      const tab = r.tables[T] || [];
      const tot = r.totals[T] || {};
      const flagged = tab.reduce((s, x) => s + x.ok - x.clean, 0);
      const kp = (v, l, cls) => h("div", {class: "kpi"}, h("b", {class: cls || ""}, v), h("span", {}, l));
      body.append(h("div", {class: "card"}, h("div", {class: "kpis"},
        kp(r.compositions.length, "compositions"),
        kp(r.phases.join(" · "), "phases"),
        kp(tot.structures ?? 0, "structures"),
        kp(tot.relaxed ?? 0, "relaxed (files on disk)"),
        kp(tot.ok ?? 0, "free energies done"),
        kp(tot.imaginary ?? 0, "with imaginary modes"),
        kp(tot.unconverged ?? 0, "unconverged relaxation"),
        kp(tot.failed ?? 0, "failed"),
        kp(tot.not_run ?? 0, "not run yet"),
        r.chunks ? kp(`${r.chunks.complete}/${r.chunks.n}`, "chunks complete") : null),
        h("div", {class: "small muted", style: {marginTop: "6px"}},
          "Harmonic F at fixed relaxed cell (not QHA). Flagged = imaginary modes and/or unconverged relaxation; those rows are kept in the CSV but are not reliable.")));

      // heatmap (composition x phase) and one row per structure
      const tblBox = h("div", {});
      const tbl = S.structTable(r.root, T, filt || {}, tblBox);
      body.append(S.heatmap(r, T, {onPick: (comp, phase) => tbl.setFilter({comp, phase}), selected: filt || {}}), tblBox);

      // per composition overview
      const comps = r.compositions;
      const by = {};
      for (const x of tab) (by[x.comp] = by[x.comp] || {})[x.phase] = x;
      const phases = r.phases.filter((p) => comps.some((c) => by[c.dir] && by[c.dir][p]));
      const t1 = h("table", {});
      t1.append(h("thead", {}, h("tr", {}, h("th", {class: "l"}, "composition"), h("th", {}, "x %"), h("th", {}, "n"),
        h("th", {}, "orderings"), h("th", {class: "l"}, "prepared"),
        phases.map((p) => h("th", {class: "l"}, p)), h("th", {class: "l"}, "lowest F (unflagged)"))));
      const tb = h("tbody", {});
      for (const c of comps) {
        const row = by[c.dir] || {};
        const best = Object.values(row).filter((x) => x.dF_lowest_meV != null).sort((a, b) => a.dF_lowest_meV - b.dF_lowest_meV);
        tb.append(h("tr", {},
          h("td", {class: "l mono"}, c.dir), h("td", {}, S.fmt(c.pct, 2)), h("td", {}, c.n ?? "–"),
          h("td", {}, c.orderings), h("td", {class: "l"}, c.prepared ? h("span", {class: "pill ok"}, "yes") : h("span", {class: "pill warn"}, "no")),
          phases.map((p) => h("td", {class: "l"}, cellBar(row[p]))),
          h("td", {class: "l small"}, !best.length ? "–" :
            best[1] && best[1].dF_lowest_meV < 0.005 ? "tie: " + best.filter((x) => x.dF_lowest_meV < 0.005).map((x) => x.phase).join(" = ") + " (< 0.005 meV/atom apart)" :
            [h("b", {}, best[0].phase), best[1] ? ` (next ${best[1].phase} +${S.fmt(best[1].dF_lowest_meV, 2)} meV/atom)` : ""])));
      }
      t1.append(tb);
      body.append(h("details", {class: "card"}, h("summary", {}, "Progress per composition"),
        h("div", {class: "small muted"}, "bar: ", h("span", {class: "pill ok"}, "F done, clean"), " ", h("span", {class: "pill warn"}, "flagged"), " ",
          h("span", {class: "pill bad"}, "failed"), " grey = not run. Numbers: done / structures. Hover for details."),
        h("div", {class: "tablewrap"}, t1)));

      // detail table
      const t2 = h("table", {});
      t2.append(h("thead", {}, h("tr", {}, ["composition", "phase"].map((x) => h("th", {class: "l"}, x)),
        ["structures", "planned", "relaxed", "F ok", "clean", "imag.", "unconv.", "failed", "not run", "F min eV/atom", "F median", "lowest ΔF meV/atom", "median ΔF"].map((x) => h("th", {}, x)))));
      const tb2 = h("tbody", {});
      for (const c of comps) for (const p of phases) {
        const x = (by[c.dir] || {})[p];
        if (!x) continue;
        tb2.append(h("tr", {}, h("td", {class: "l mono"}, c.dir), h("td", {class: "l"}, p),
          [x.structures, x.planned, x.relaxed, x.ok, x.clean].map((v) => h("td", {}, v)),
          h("td", {class: x.imaginary ? "" : "muted"}, x.imaginary), h("td", {class: x.unconverged ? "" : "muted"}, x.unconverged),
          h("td", {class: x.failed ? "" : "muted"}, x.failed), h("td", {class: x.not_run ? "" : "muted"}, x.not_run),
          h("td", {}, x.F_stats ? S.fmt(x.F_stats.min, 5) : "–"), h("td", {}, x.F_stats ? S.fmt(x.F_stats.median, 5) : "–"),
          h("td", {}, S.fmt(x.dF_lowest_meV, 2)), h("td", {}, S.fmt(x.dF_median_meV, 2))));
      }
      t2.append(tb2);
      body.append(h("details", {class: "card"}, h("summary", {}, "Table per composition and phase"),
        h("div", {class: "small muted"}, "ΔF relative to the lowest unflagged F of the same composition (any phase)."),
        h("div", {class: "tablewrap"}, t2)));

      // mapping
      const t3 = h("table", {});
      t3.append(h("thead", {}, h("tr", {}, h("th", {class: "l"}, "composition"), h("th", {class: "l"}, "target"),
        ["mapped ok", "failed", "IDENTICAL", "not identical", "rms shift Å", "strain S", "equally good variants"].map((x) => h("th", {}, x)))));
      const tb3 = h("tbody", {});
      let anyTied = false;
      for (const c of comps) for (const [t, m] of Object.entries(c.mapping || {})) {
        const v = (c.verify || {})[t] || {identical: 0, other: 0};
        const rng = (s, d) => (s ? (s.min === s.max ? S.fmt(s.min, d) : S.fmt(s.min, d) + "–" + S.fmt(s.max, d)) : "–");
        if (m.tied_variants && m.tied_variants.max > 1) anyTied = true;
        tb3.append(h("tr", {}, h("td", {class: "l mono"}, c.dir), h("td", {class: "l"}, t), h("td", {}, m.ok), h("td", {}, m.failed),
          h("td", {}, v.identical), h("td", {}, v.other ? h("span", {class: "pill bad"}, v.other) : 0),
          h("td", {}, rng(m.rms_shift_A, 3)), h("td", {}, rng(m.strain_S, 4)), h("td", {}, rng(m.tied_variants, 0))));
      }
      t3.append(tb3);
      body.append(h("details", {class: "card"}, h("summary", {}, "Mapping and ordering check (batch_map + compare_ordering)"),
        anyTied ? h("div", {class: "notice warn small"}, "Several equally good variants: they cannot be ranked by distortion, only by energy. Variant 1 was computed unless map_args has --all-variants.") : null,
        h("div", {class: "tablewrap"}, t3)));

      // CSVs, reports, reruns, logs
      const files = h("div", {class: "card"}, h("h2", {}, "Data files"));
      if (r.merged_csv.length) files.append(h("h3", {}, "Merged CSV (whole composition range)"),
        r.merged_csv.map((f) => { const p = r.workdir + "/" + f; return h("div", {class: "row"}, h("span", {class: "mono"}, f),
          h("a", {href: "#/plots?csv=" + encodeURIComponent(p)}, "plot"), h("a", {href: S.fileUrl(p, true)}, "download")); }));
      else files.append(h("div", {class: "notice warn small"}, "No merged CSV yet: run Merge on the Pipeline page" +
        (r.config_path ? "" : " (or campaign.py CONFIG merge)") + ". The plots of the whole composition range need it."),
        r.config_path ? h("a", {href: "#/pipeline?config=" + encodeURIComponent(r.config_path)}, "open the run's pipeline") : null);
      const per = Object.entries(r.per_comp_csv).filter(([, v]) => v.length);
      if (per.length) files.append(h("h3", {}, "Per composition (written by merge)"), h("div", {class: "row small"},
        per.map(([d, v]) => v.map((f) => h("a", {href: "#/plots?csv=" + encodeURIComponent(r.workdir + "/" + d + "/" + f)}, d + "/" + f)))));
      for (const [n, txt] of Object.entries(r.merge_reports)) files.append(h("details", {}, h("summary", {}, n), h("pre", {class: "log"}, txt)));
      if (r.reruns.length) files.append(h("h3", {}, "Rerun rounds"), r.reruns.map((x) => h("div", {class: "small mono"}, JSON.stringify(x))));
      if (r.has_logs || r.has_slurm) {
        const lg = h("div", {});
        files.append(h("details", {ontoggle: async (e) => {
          if (!e.target.open || lg.childElementCount) return;
          const l = await S.api("/api/logs?workdir=" + encodeURIComponent(r.workdir));
          const pv = h("div", {});
          lg.append(h("div", {class: "small"}, (l.logs || []).slice(0, 60).map((x) => h("div", {}, h("a", {href: "javascript:void 0", onclick: () => { pv.innerHTML = ""; pv.append(S.preview(x.path)); }}, x.name), " ", h("span", {class: "muted"}, S.ago(x.mtime) + " · " + S.size(x.size))))), pv);
        }}, h("summary", {}, "Logs and SLURM scripts"), lg));
      }
      if (r.chunks) {
        const unf = r.chunks.list.filter((c) => c.ok < c.n);
        files.append(h("h3", {}, "Chunks"), h("div", {class: "small"},
          `${r.chunks.n} chunks, ${r.chunks.complete} complete, ${r.chunks.started} started; ${r.chunks.workers_per_node} structures in parallel per node, walltime ${r.chunks.walltime}, sec_per_structure ${r.chunks.sec_per_structure}`),
          unf.length ? h("div", {class: "small muted"}, "unfinished: " + unf.slice(0, 80).map((c) => `${c.id} (${c.ok}/${c.n})`).join(", ") + (unf.length > 80 ? " …" : "")) : null);
      }
      body.append(files);
    }
  }

  function cellBar(x) {
    if (!x) return h("span", {class: "muted"}, "–");
    const n = x.structures || 1;
    const w = (v) => (100 * v / n).toFixed(1) + "%";
    const flagged = x.ok - x.clean;
    const tip = `${x.structures} structures, ${x.relaxed} relaxed, ${x.ok} F done (${x.clean} clean, ${x.imaginary} imaginary, ${x.unconverged} unconverged), ${x.failed} failed, ${x.not_run} not run` +
      (x.dF_lowest_meV != null ? `; lowest ΔF ${S.fmt(x.dF_lowest_meV, 2)} meV/atom` : "");
    return h("div", {class: "row", title: tip, style: {gap: "6px", flexWrap: "nowrap"}},
      h("div", {class: "bar"}, h("i", {class: "c-ok", style: {width: w(x.clean)}}), h("i", {class: "c-flag", style: {width: w(flagged)}}),
        h("i", {class: "c-bad", style: {width: w(x.failed)}})),
      h("span", {class: "small"}, `${x.ok}/${x.structures}`));
  }
})();
