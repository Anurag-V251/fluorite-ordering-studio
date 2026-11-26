/* Results pieces: the composition x phase heatmap and the table with one row per structure (sortable, filter chips,
   sticky first column, relaxed POSCAR download).  Used by the Results page and the run dashboard. */
(function () {
  const S = Studio, h = S.h;

  // ------------------------------------------------------------------ heatmap
  const METRICS = {
    done: {label: "Progress", tip: "share of the structures with a free energy (clean or flagged)"},
    dF: {label: "Lowest ΔF", tip: "lowest unflagged F of the phase minus the lowest of the composition, meV/atom (0 = most stable phase)"},
    flagged: {label: "Flagged", tip: "structures with imaginary modes or an unconverged relaxation"},
  };
  S.heatmap = function (scan, T, opts = {}) {
    const tab = (scan.tables || {})[T] || [];
    const by = {};
    for (const x of tab) (by[x.comp] = by[x.comp] || {})[x.phase] = x;
    const comps = (scan.compositions || []).filter((c) => by[c.dir]);
    const phases = (scan.phases || []).filter((p) => comps.some((c) => by[c.dir][p]));
    let metric = S.store.get("hm.metric", "done");
    if (!METRICS[metric]) metric = "done";
    let sel = {comp: opts.selected && opts.selected.comp || "", phase: opts.selected && opts.selected.phase || ""};
    const card = h("div", {class: "card flush", id: "heatmap"});
    const tools = h("div", {class: "hm-tools", role: "group", "aria-label": "Colour by"});
    const grid = h("div", {class: "hm", role: "grid", "aria-label": "Compositions × phases"});
    const legend = h("div", {class: "hm-legend"});
    card.append(h("div", {class: "card-h"}, h("h2", {}, "Compositions × phases"),
      h("span", {class: "sub"}, opts.onPick ? "click a cell to filter the structures" : ""), h("div", {class: "right"}, tools)),
      h("div", {class: "card-b"}, comps.length ? grid : h("div", {class: "muted small"}, "Nothing prepared yet."), comps.length ? legend : null));
    for (const [k, m] of Object.entries(METRICS)) tools.append(h("button", {"aria-pressed": String(k === metric), title: m.tip, "data-m": k,
      onclick: () => { metric = k; S.store.set("hm.metric", k); tools.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.m === k))); draw(); }}, m.label));
    function val(x) {
      if (!x) return null;
      if (metric === "done") return x.structures ? (x.ok) / x.structures : 0;
      if (metric === "dF") return x.dF_lowest_meV;
      return x.ok - x.clean;
    }
    function draw() {
      grid.innerHTML = "";
      grid.style.gridTemplateColumns = `64px repeat(${comps.length}, minmax(34px, 1fr))`;
      const vals = comps.flatMap((c) => phases.map((p) => val(by[c.dir][p]))).filter((v) => v != null && isFinite(v));
      const vmax = Math.max(1e-9, ...vals);
      for (const p of phases) {
        grid.append(h("button", {class: "rl" + (sel.phase === p && !sel.comp ? " sel" : ""), title: "all compositions of " + p,
          onclick: () => pick("", sel.phase === p && !sel.comp ? "" : p)}, p));
        for (const c of comps) {
          const x = by[c.dir][p];
          if (!x) { grid.append(h("div", {class: "cell none", "aria-hidden": "true"})); continue; }
          const v = val(x);
          let pct, txt;
          if (metric === "done") { pct = 15 + 85 * v; txt = Math.round(100 * v) + ""; }
          else if (metric === "dF") { pct = v == null ? 0 : 100 - 85 * Math.min(1, v / Math.max(5, vmax)); txt = v == null ? "–" : v < 0.005 ? "0" : v < 10 ? S.fmt(v, 1) : S.fmt(v, 0); }
          else { pct = v ? 25 + 75 * v / vmax : 0; txt = String(v); }
          const col = metric === "flagged" ? "var(--hm-warn)" : "var(--hm)";
          const fl = x.ok - x.clean;
          const isSel = sel.comp === c.dir && (!sel.phase || sel.phase === p);
          grid.append(h("button", {class: "cell" + (pct > 55 ? " hi" : "") + (isSel ? " sel" : "") + (sel.comp && !isSel || sel.phase && sel.phase !== p ? " dim" : ""),
            style: {background: `color-mix(in srgb, ${col} ${pct.toFixed(0)}%, var(--hm-0))`}, "data-comp": c.dir, "data-phase": p,
            title: `${c.dir} · ${p}: ${x.structures} structures, ${x.ok} done (${x.clean} clean, ${fl} flagged), ${x.failed} failed, ${x.not_run} not run` +
              (x.dF_lowest_meV != null ? `; lowest ΔF ${S.fmt(x.dF_lowest_meV, 2)} meV/atom` : ""),
            onclick: () => pick(isSel && sel.phase === p ? "" : c.dir, isSel && sel.phase === p ? "" : p)},
            txt, fl || x.failed ? h("span", {class: "fl"}, fl ? h("i", {}) : null, x.failed ? h("i", {class: "f"}) : null) : null));
        }
      }
      grid.append(h("span", {}));
      for (const c of comps) grid.append(h("button", {class: "cl" + (sel.comp === c.dir && !sel.phase ? " sel" : ""), title: `${c.dir}: ${S.fmt(c.pct, 2)} % (n = ${c.n ?? "?"})`,
        onclick: () => pick(sel.comp === c.dir && !sel.phase ? "" : c.dir, "")}, S.fmt(c.pct, c.pct % 1 ? 1 : 0)));
      legend.innerHTML = "";
      legend.append(h("span", {}, metric === "done" ? "% done" : metric === "dF" ? "lowest ΔF, meV/atom (darker = more stable)" : "flagged structures"),
        h("span", {class: "ramp", style: {background: `linear-gradient(90deg, var(--hm-0), ${metric === "flagged" ? "var(--hm-warn)" : "var(--hm)"})`}}),
        h("span", {class: "lg"}, h("i", {}), "has flagged"), h("span", {class: "lg"}, h("i", {class: "f"}), "has failed"),
        h("span", {class: "muted"}, "columns: dopant %"));
    }
    function pick(comp, phase) {
      sel = {comp, phase};
      draw();
      if (opts.onPick) opts.onPick(comp, phase);
    }
    S.on("hm-clear", () => { sel = {comp: "", phase: ""}; draw(); });
    draw();
    return card;
  };

  // ------------------------------------------------------------------ one row per structure
  const STATUS = [["clean", "ok", "Clean"], ["flagged", "warn", "Flagged"], ["imaginary", "warn", "Imaginary"], ["unconverged", "warn", "Unconverged"],
    ["failed", "bad", "Failed"], ["not_run", "idle", "Not run"]];
  const COLS = [["file", "structure", "l"], ["comp", "composition", "l"], ["phase", "phase", "l"], ["status", "status", "l"],
    ["dF_meV", "ΔF meV/atom", "r", "F minus the lowest unflagged F of the composition"], ["F", "F eV/atom", "r"],
    ["min_freq", "min freq THz", "r", "lowest phonon frequency (negative = imaginary)"], ["fmax", "max force eV/Å", "r"],
    ["a", "a Å", "r"], ["b", "b Å", "r"], ["c", "c Å", "r"], ["time_s", "time s", "r"], ["WC", "Warren-Cowley", "l", "short-range order of the ordering"]];
  S.structTable = function (path, T, filt, box) {
    let st = {comp: filt.comp || "", phase: filt.phase || "", status: filt.status ? filt.status.split(",") : [], q: "",
      sort: S.store.get("rt.sort", "dF_meV"), desc: S.store.get("rt.desc", false), offset: 0, limit: 200};
    const card = h("div", {class: "card flush", id: "structures"});
    box.append(card);
    let seq = 0;
    async function load() {
      const my = ++seq;
      const q = new URLSearchParams({path, T: T || "", comp: st.comp, phase: st.phase, status: st.status.join(","), q: st.q,
        sort: st.sort, desc: st.desc ? "1" : "0", offset: String(st.offset), limit: String(st.limit)});
      const r = await S.api("/api/results/rows?" + q.toString(), undefined, {quiet: true});
      if (my !== seq) return;
      draw(r);
    }
    const search = h("input", {type: "search", placeholder: "file or WC…", "aria-label": "Search structures", value: st.q});
    let stmr = null;
    search.addEventListener("input", () => { clearTimeout(stmr); stmr = setTimeout(() => { st.q = search.value; st.offset = 0; load(); }, 250); });
    function draw(r) {
      card.innerHTML = "";
      if (!r.ok) { card.append(h("div", {class: "card-b"}, h("div", {class: "notice warn"}, r.error || "no rows"))); return; }
      const chips = h("div", {class: "chips", role: "group", "aria-label": "Status"}, STATUS.map(([k, kind, l]) =>
        h("button", {class: "chip", "aria-pressed": String(st.status.includes(k)), "data-status": k,
          onclick: () => { st.status = st.status.includes(k) ? st.status.filter((x) => x !== k) : [...st.status, k]; st.offset = 0; load(); }},
          S.icon({ok: "check", warn: "alert", bad: "xcircle", idle: "circle"}[kind], "sm " + kind), l, h("span", {class: "c num"}, S.int(r.counts[k] || 0)))));
      const phaseSel = h("select", {class: "inp", style: {width: "auto"}, "aria-label": "Phase", onchange: () => { st.phase = phaseSel.value; st.offset = 0; load(); }},
        h("option", {value: ""}, "all phases"), r.phases.map((p) => h("option", {value: p}, p)));
      phaseSel.value = st.phase;
      const compSel = h("select", {class: "inp", style: {width: "auto"}, "aria-label": "Composition", onchange: () => { st.comp = compSel.value; st.offset = 0; load(); }},
        h("option", {value: ""}, "all compositions"), r.compositions.map((c) => h("option", {value: c.dir}, `${c.dir} (${S.fmt(c.pct, 2)} %)`)));
      compSel.value = st.comp;
      const anyF = st.comp || st.phase || st.status.length || st.q;
      card.append(h("div", {class: "card-h"}, h("h2", {}, "Structures"),
        h("span", {class: "sub"}, `${S.int(r.total)} of ${S.int(r.all)} at ${r.T} K · newest result per structure`)),
        h("div", {class: "rbar"}, chips, h("span", {class: "vsep"}), compSel, phaseSel,
          h("label", {class: "search"}, S.icon("search", "sm"), search),
          anyF ? h("button", {class: "btn sm ghost", onclick: () => { st = {...st, comp: "", phase: "", status: [], q: "", offset: 0}; search.value = ""; S.emit("hm-clear"); load(); }}, S.icon("close", "sm"), "Clear filters") : null));
      const th = COLS.map(([k, l, al, tip]) => {
        const on = st.sort === k;
        return h("th", {class: al === "r" ? "r" : "", "aria-sort": on ? (st.desc ? "descending" : "ascending") : null, title: tip || null},
          h("button", {onclick: () => { if (st.sort === k) st.desc = !st.desc; else { st.sort = k; st.desc = false; }
            S.store.set("rt.sort", st.sort); S.store.set("rt.desc", st.desc); st.offset = 0; load(); }},
          l, h("span", {class: "arr"}, on ? (st.desc ? "▼" : "▲") : "↕")));
      });
      th.push(h("th", {}, "relaxed"));
      const pill = (x) => x.status === "clean" ? S.pill("ok", "clean") : x.status === "failed" ? S.pill("bad", "failed", x.note) :
        x.status === "not_run" ? S.pill("idle", x.planned ? "not run" : "not planned") :
        S.pill("warn", [x.imaginary ? `${x.n_imag || ""} imag.`.trim() : "", x.unconverged ? "unconv." : ""].filter(Boolean).join(" · "));
      const num = (v, d) => (v == null ? h("td", {class: "r na"}, "–") : h("td", {class: "r num"}, S.fmt(v, d)));
      const rows = r.rows.map((x) => h("tr", {"data-file": x.file},
        h("td", {class: "mono small", title: x.file}, x.file.split("/").slice(-1)[0]),
        h("td", {class: "mono small"}, x.comp), h("td", {}, x.phase), h("td", {}, pill(x)),
        num(x.dF_meV, 2), num(x.F, 5), x.min_freq != null && x.min_freq < 0 ? h("td", {class: "r num neg"}, S.fmt(x.min_freq, 3)) : num(x.min_freq, 3),
        num(x.fmax, 4), num(x.a, 3), num(x.b, 3), num(x.c, 3), num(x.time_s, 0), h("td", {class: "small mono"}, x.WC || ""),
        h("td", {}, x.relaxed_path ? h("a", {class: "btn sm ghost", href: S.fileUrl(x.relaxed_path, true), title: "download the relaxed POSCAR",
          download: S.basename(x.file) + "_relaxed"}, S.icon("download", "sm"), "POSCAR") : h("span", {class: "muted small"}, "–"))));
      card.append(h("div", {class: "rwrap"}, h("table", {class: "tbl rt"}, h("thead", {}, h("tr", {}, th)),
        h("tbody", {}, rows.length ? rows : h("tr", {}, h("td", {colspan: th.length, class: "muted"}, "No structure matches."))))));
      const last = Math.min(r.total, st.offset + st.limit);
      card.append(h("div", {class: "pager"}, `${r.total ? S.int(st.offset + 1) : 0}–${S.int(last)} of ${S.int(r.total)}`,
        h("button", {class: "btn sm", disabled: st.offset ? null : true, onclick: () => { st.offset = Math.max(0, st.offset - st.limit); load(); }}, "Previous"),
        h("button", {class: "btn sm", disabled: last < r.total ? null : true, onclick: () => { st.offset += st.limit; load(); }}, "Next"),
        h("span", {class: "grow"}), h("span", {class: "muted"}, "ΔF and flags as in the merged CSV: the newest row of a structure wins (reruns included).")));
    }
    load();
    return {setFilter: (f) => { st = {...st, ...f, offset: 0}; load(); card.scrollIntoView({block: "start", behavior: "smooth"}); }, reload: load};
  };
})();
