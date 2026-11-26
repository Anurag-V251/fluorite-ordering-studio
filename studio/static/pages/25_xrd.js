/* XRD page: pick an F window on P(F) of one composition, compute the powder XRD of every structure in it
   (scripts/xrd_average.py: ASE + pymatgen, run as a job, sticks cached per file), and show ONE averaged
   pattern.  The browser broadens and averages the sticks, so FWHM, step, profile, normalisation and the
   display change the plot at once; wavelength, 2theta range (beyond what was computed), species replacement,
   references and the selection recompute the sticks (only the new ones: the cache keeps the rest).
   The relaxed POSCARs are only read: species replacement happens in memory, for the XRD only. */
(function () {
  const S = Studio, h = S.h;
  const NS = "http://www.w3.org/2000/svg";
  const sv = (tag, attrs, ...kids) => {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs || {})) if (v != null && v !== false) e.setAttribute(k, v);
    for (const c of kids.flat(9)) if (c != null && c !== false) e.append(c.nodeType ? c : document.createTextNode(String(c)));
    return e;
  };
  const PHASE = {parent: "#52514e", c: "#52514e", t: "#2a78d6", m: "#eb6834", oIII: "#1baf7a"};
  const EXTRA = ["#6250d6", "#e34948", "#eda100", "#e87ba4"];
  const WL = {CuKa: 1.54184, CuKa1: 1.54056, CuKa2: 1.54439, CuKb1: 1.39222, MoKa: 0.71073, MoKa1: 0.7093, MoKa2: 0.71359,
    MoKb1: 0.63229, CoKa: 1.79026, CoKa1: 1.78896, CoKa2: 1.79285, CoKb1: 1.63079, FeKa: 1.93735, FeKa1: 1.93604,
    FeKa2: 1.93998, FeKb1: 1.75661, CrKa: 2.291, CrKa1: 2.2897, CrKa2: 2.29361, CrKb1: 2.08487, AgKa: 0.560885,
    AgKa1: 0.559421, AgKa2: 0.563813, AgKb1: 0.497082};
  const WINDOW = {gaussian: 6, lorentzian: 20, "pseudo-voigt": 20};    // same cut-offs as xrd_average.py
  const LN2x4 = 4 * Math.log(2);
  const DEF = {source: "", T: "", comp: "", phases: null, lo: null, hi: null, flagged: false, unrelaxed: false, max: 500, seed: 0,
    bins: 40, fits: true, ksig: 1, xlo: "", xhi: "",
    wl: "CuKa", wlc: "1.5406", tt0: 10, tt1: 90, step: 0.02, fwhm: 0.1, profile: "gaussian", eta: 0.5, cag: false, U: 0.004, V: -0.002, W: 0.006,
    norm: "max", spread: "std", indiv: false, yscale: "linear", species: "none", spTo: "", spMap: "", spAuto: true, anions: "O",
    dw: "", symprec: 0, minrel: 1e-4, refs: ["c", "t", "oIII", "m"], refsDir: "", refReplace: true, refScale: false,
    refMode: "stack", ticks: true, refLabels: true, peakLabels: true, auto: true, jobs: 4, outDir: ""};
  const RECOMPUTE = ["source", "T", "comp", "phases", "lo", "hi", "flagged", "unrelaxed", "max", "seed", "wl", "wlc", "species", "spTo",
    "spMap", "anions", "dw", "symprec", "minrel", "refs", "refsDir", "refReplace", "refScale", "jobs"];

  S.registerPage({id: "xrd", title: "XRD", order: 45, icon: "xrd", render});

  function render(main, params) {
    let st = Object.assign({}, DEF, S.store.get("xrd.state", {}));
    if (params.source) st.source = params.source;
    if (!st.source) st.source = (S.info.prefs || {}).last_results || "";
    if (!st.refsDir) st.refsDir = S.info.app + "/refs";
    const save = () => S.store.set("xrd.state", st);
    let D = null;            // rows of the composition
    let X = null;            // computed sticks {key, sig, range, structures, refs, ...}
    let job = null, pollT = null, debT = null, zoom = null, ycache = null;
    S.onLeave(() => { clearTimeout(pollT); clearTimeout(debT); });
    const mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
    const reTheme = () => { drawPF(); drawXRD(); };            // plot colours follow the light / dark theme
    if (mq && mq.addEventListener) { mq.addEventListener("change", reTheme); S.onLeave(() => mq.removeEventListener("change", reTheme)); }
    window.addEventListener("studio-theme", reTheme);                  // the theme button in the top bar
    S.onLeave(() => window.removeEventListener("studio-theme", reTheme));

    // ---------------------------------------------------------------- data card
    const pf = S.pathField(st.source, {title: "Results: run folder, campaign workdir, folder with Ce_xx/, or a result CSV"});
    const tsel = h("select", {onchange: () => { st.T = tsel.value; load(); }});
    const csel = h("select", {onchange: () => { st.comp = csel.value; st.lo = st.hi = null; load(); }});
    const dinfo = h("span", {class: "muted small grow"});
    main.append(h("div", {class: "card"},
      h("h2", {}, "XRD averaged over an F window"),
      h("div", {class: "row"}, pf, h("label", {}, "T ", tsel), h("label", {}, "composition ", csel),
        h("button", {class: "primary", onclick: () => { st.source = pf.input.value.trim(); st.lo = st.hi = null; load(); }}, "Load")),
      h("div", {class: "row", style: {marginTop: "6px"}}, dinfo),
      h("div", {class: "small muted", style: {marginTop: "4px"}},
        "Drag on P(F) to choose a window (or type it), pick the phases, and every relaxed structure in the window is computed with ",
        "pymatgen; the plot below is their average. The relaxed POSCARs are only read: a species change is applied in memory, for the XRD only.")));
    pf.input.addEventListener("keydown", (e) => { if (e.key === "Enter") { st.source = pf.input.value.trim(); load(); } });

    // ---------------------------------------------------------------- P(F) card
    const pfBox = h("div", {class: "xrd-plot"});
    const chips = h("div", {class: "row xrd-chips"});
    const loIn = h("input", {type: "number", step: "any", style: {width: "90px"}});
    const hiIn = h("input", {type: "number", step: "any", style: {width: "90px"}});
    const absTxt = h("span", {class: "muted small"});
    const gSel = h("select", {style: {maxWidth: "330px"}});
    const kIn = h("input", {type: "number", step: "0.1", min: "0.1", value: st.ksig, style: {width: "60px"}});
    const binsIn = h("input", {type: "number", min: "5", max: "400", value: st.bins, style: {width: "64px"}});
    const xloIn = h("input", {type: "number", step: "any", value: st.xlo, placeholder: "auto", style: {width: "74px"}});
    const xhiIn = h("input", {type: "number", step: "any", value: st.xhi, placeholder: "auto", style: {width: "74px"}});
    const fitsCb = h("input", {type: "checkbox"}); fitsCb.checked = st.fits;
    const flagCb = h("input", {type: "checkbox"}); flagCb.checked = st.flagged;
    const unrCb = h("input", {type: "checkbox"}); unrCb.checked = st.unrelaxed;
    const maxIn = h("input", {type: "number", min: "1", value: st.max, style: {width: "80px"}});
    const seedIn = h("input", {type: "number", min: "0", value: st.seed, style: {width: "60px"}});
    const selTxt = h("div", {class: "small", style: {marginTop: "6px"}});
    const autoCb = h("input", {type: "checkbox"}); autoCb.checked = st.auto;
    const goBtn = h("button", {class: "primary", onclick: () => compute(true)}, "Compute XRD");
    main.append(h("div", {class: "card"},
      h("div", {class: "row"}, h("h2", {style: {margin: 0}}, "P(F) and the window"), h("span", {class: "grow"}), chips),
      pfBox,
      h("div", {class: "row small", style: {marginTop: "6px"}},
        h("b", {}, "window"), loIn, "to", hiIn, "meV/atom above the lowest unflagged F", absTxt,
        h("button", {onclick: () => { st.lo = st.hi = null; save(); drawSel(); sched(); }}, "whole range")),
      h("div", {class: "row small", style: {marginTop: "6px"}},
        h("span", {}, "Gaussian:"), gSel, h("label", {}, "±", kIn, " σ"),
        h("label", {}, fitsCb, " show fits"), h("label", {}, "bins ", binsIn), h("label", {}, "F axis ", xloIn, " – ", xhiIn)),
      h("div", {class: "row small", style: {marginTop: "6px"}},
        h("label", {title: "rows with imaginary phonon modes or an unconverged relaxation"}, flagCb, " include flagged rows"),
        h("label", {title: "use the unrelaxed input when the relaxed file is missing"}, unrCb, " allow unrelaxed"),
        h("label", {}, "at most ", maxIn, " structures (random, seed ", seedIn, ")")),
      selTxt,
      h("div", {class: "row", style: {marginTop: "8px"}}, goBtn, h("label", {class: "small"}, autoCb, " recompute automatically"))));

    // ---------------------------------------------------------------- XRD card
    const num = (k, w, attrs) => { const e = h("input", {type: "number", step: "any", value: st[k], style: {width: (w || 70) + "px"}, ...(attrs || {})});
      e.addEventListener("input", () => { if (e.value === "" || !isFinite(Number(e.value))) return; st[k] = Number(e.value); changed(k); }); return e; };
    const txt = (k, w, ph) => { const e = h("input", {type: "text", value: st[k], placeholder: ph || "", style: {width: (w || 90) + "px"}});
      e.addEventListener("change", () => { st[k] = e.value.trim(); changed(k); }); return e; };
    const sel = (k, opts) => { const e = h("select", {}, opts.map(([v, l]) => h("option", {value: v}, l)));
      e.value = String(st[k]); e.addEventListener("change", () => { st[k] = e.value; changed(k); }); return e; };
    const cb = (k) => { const e = h("input", {type: "checkbox"}); e.checked = !!st[k];
      e.addEventListener("change", () => { st[k] = e.checked; changed(k); }); return e; };
    const wlSel = sel("wl", [...Object.keys(WL).map((k) => [k, `${k} (${WL[k]} Å)`]), ["custom", "custom (Å)"]]);
    const wlc = txt("wlc", 70);
    const spSel = sel("species", [["none", "as in the files"], ["cations", "every cation →"], ["map", "replace pairs"]]);
    const spTo = txt("spTo", 44, "host"), spMap = txt("spMap", 110, "dopant:host"), anions = txt("anions", 50, "O");
    [spTo, spMap].forEach((e) => e.addEventListener("change", () => { st.spAuto = false; save(); }));
    const refBoxes = h("span", {class: "row", style: {gap: "8px"}});
    const refsDir = S.pathField(st.refsDir, {dir: true});
    refsDir.style.flexWrap = "nowrap"; refsDir.input.style.minWidth = "0"; refsDir.style.minWidth = "0";
    refsDir.input.addEventListener("change", () => { st.refsDir = refsDir.input.value.trim(); listRefs(); changed("refsDir"); });
    const cagRow = h("span", {}, "U ", num("U", 64), " V ", num("V", 64), " W ", num("W", 64));
    const etaIn = num("eta", 56, {min: 0, max: 1});
    const P = (label, ...el) => h("div", {class: "xrd-p"}, h("span", {class: "k"}, label), h("span", {class: "row", style: {gap: "6px"}}, el));
    const parBox = h("div", {class: "xrd-params"},
      h("div", {class: "grp"}, "Pattern (pymatgen)"),
      P("radiation", wlSel, wlc),
      P("2θ range (°)", num("tt0", 60), "–", num("tt1", 60)),
      P("species for the XRD", spSel, spTo, spMap),
      P("anions (kept)", anions),
      h("details", {}, h("summary", {class: "small"}, "more pymatgen options"),
        P("Debye–Waller B (Å²)", txt("dw", 130, "Zr=0.3 O=0.5")),
        P("symprec", num("symprec", 70), h("span", {class: "muted small"}, "0 = cell as it is")),
        P("drop sticks below", num("minrel", 80), h("span", {class: "muted small"}, "× strongest")),
        P("parallel processes", num("jobs", 56, {min: 1}))),
      h("div", {class: "grp"}, "Broadening and averaging (instant)"),
      P("step (°)", num("step", 70, {min: 0.0005})),
      P("profile", sel("profile", [["gaussian", "Gaussian"], ["lorentzian", "Lorentzian"], ["pseudo-voigt", "pseudo-Voigt"]]), h("span", {class: "small"}, "η"), etaIn),
      P("FWHM (°)", num("fwhm", 70, {min: 0.001}), h("label", {class: "small"}, cb("cag"), " Caglioti")),
      P("", cagRow),
      P("normalise", sel("norm", [["max", "each pattern max = 100"], ["absolute", "absolute (÷ n_atoms²)"]])),
      h("div", {class: "grp"}, "References (12-atom cells, as they are)"),
      P("phases", refBoxes),
      P("folder", refsDir),
      P("", h("label", {class: "small"}, cb("refReplace"), " same species change"), h("label", {class: "small"}, cb("refScale"), " scale to the selection's volume/atom")),
      P("show", sel("refMode", [["stack", "stacked below"], ["overlay", "overlaid"], ["none", "hidden"]])),
      P("", h("label", {class: "small"}, cb("ticks"), " Bragg ticks"), h("label", {class: "small"}, cb("refLabels"), " hkl labels")),
      h("div", {class: "grp"}, "Display"),
      P("spread", sel("spread", [["std", "± 1 std"], ["minmax", "min – max"], ["none", "none"]]), h("label", {class: "small"}, cb("indiv"), " single patterns")),
      P("y axis", sel("yscale", [["linear", "linear"], ["sqrt", "square root"]]), h("label", {class: "small"}, cb("peakLabels"), " peak 2θ")));
    const xrdBox = h("div", {class: "xrd-plot"});
    const readout = h("div", {class: "small mono xrd-readout"}, " ");
    const status = h("div", {class: "small"});
    const cliPre = h("pre", {class: "log", style: {maxHeight: "120px"}}, "");
    const outDir = S.pathField(st.outDir, {dir: true, placeholder: "output folder"});
    const exportBox = h("div", {});
    const listBox = h("div", {});
    main.append(h("div", {class: "card"},
      h("div", {class: "row"}, h("h2", {style: {margin: 0}}, "Averaged XRD"), status, h("span", {class: "grow"}),
        h("button", {onclick: () => { zoom = null; drawXRD(); }}, "full range"),
        h("button", {onclick: dlCSV}, "CSV"), h("button", {onclick: () => dlSVG(false)}, "SVG"), h("button", {onclick: () => dlSVG(true)}, "PNG")),
      h("div", {class: "xrd-grid"}, parBox, h("div", {}, xrdBox, readout,
        h("div", {class: "small muted"}, "Drag across the plot to zoom, double-click for the whole range. Hover for 2θ, d and the reference peaks."))),
      h("details", {}, h("summary", {}, "Structures used"), listBox),
      h("details", {}, h("summary", {}, "Command line and files"),
        h("div", {class: "small muted"}, "The same average with the standalone script (writes CSV, PNG, PDF and a JSON record):"), cliPre,
        h("div", {class: "row"}, outDir, h("button", {onclick: exportRun}, "Write files here (job)")), exportBox)));

    // ---------------------------------------------------------------- wiring
    const syncVis = () => {
      wlc.style.display = st.wl === "custom" ? "" : "none";
      spTo.style.display = st.species === "cations" ? "" : "none";
      spMap.style.display = st.species === "map" ? "" : "none";
      cagRow.style.display = st.cag ? "" : "none";
      etaIn.disabled = st.profile !== "pseudo-voigt";
    };
    syncVis();
    function changed(k) {
      save(); syncVis();
      if (RECOMPUTE.includes(k)) return sched();
      if (k === "tt0" || k === "tt1") { zoom = null; if (X && (st.tt0 < X.range[0] - 1e-9 || st.tt1 > X.range[1] + 1e-9)) return sched(); }
      ycache = null;
      drawXRD();
    }
    const winIn = () => {
      const a = loIn.value === "" ? null : Number(loIn.value), b = hiIn.value === "" ? null : Number(hiIn.value);
      if (a == null || b == null || !isFinite(a) || !isFinite(b)) return;
      st.lo = Math.min(a, b); st.hi = Math.max(a, b); save(); drawSel(); sched();
    };
    loIn.addEventListener("change", winIn); hiIn.addEventListener("change", winIn);
    binsIn.addEventListener("change", () => { st.bins = Math.max(5, Math.min(400, Number(binsIn.value) || 40)); save(); drawPF(); });
    xloIn.addEventListener("change", () => { st.xlo = xloIn.value; save(); drawPF(); });
    xhiIn.addEventListener("change", () => { st.xhi = xhiIn.value; save(); drawPF(); });
    fitsCb.addEventListener("change", () => { st.fits = fitsCb.checked; save(); drawPF(); });
    flagCb.addEventListener("change", () => { st.flagged = flagCb.checked; save(); drawPF(); sched(); });
    unrCb.addEventListener("change", () => { st.unrelaxed = unrCb.checked; save(); sched(); });
    maxIn.addEventListener("change", () => { st.max = Math.max(1, Number(maxIn.value) || 500); save(); drawSel(); sched(); });
    seedIn.addEventListener("change", () => { st.seed = Math.max(0, Number(seedIn.value) || 0); save(); sched(); });
    kIn.addEventListener("change", () => { st.ksig = Math.max(0.05, Number(kIn.value) || 1); save(); if (gSel.value) pickGauss(gSel.value); });
    autoCb.addEventListener("change", () => { st.auto = autoCb.checked; save(); });
    gSel.addEventListener("change", () => { if (gSel.value) pickGauss(gSel.value); });
    outDir.input.addEventListener("change", () => { st.outDir = outDir.input.value.trim(); save(); });

    listRefs();
    if (st.source) load(); else dinfo.textContent = "choose a results folder";

    async function listRefs() {
      const r = await S.api("/api/ls?path=" + encodeURIComponent(st.refsDir), undefined, {quiet: true});
      const names = r.ok ? r.entries.filter((e) => !e.dir && /^POSCAR_/.test(e.name)).map((e) => e.name.slice(7)) : [];
      const order = ["c", "t", "oIII", "m"];
      names.sort((a, b) => ((order.indexOf(a) + 1 || 9) - (order.indexOf(b) + 1 || 9)) || a.localeCompare(b));
      refBoxes.innerHTML = "";
      if (!names.length) refBoxes.append(h("span", {class: "muted small"}, "no POSCAR_* in this folder"));
      for (const n of names) {
        const c = h("input", {type: "checkbox"}); c.checked = st.refs.includes(n);
        c.addEventListener("change", () => { st.refs = c.checked ? [...new Set([...st.refs, n])] : st.refs.filter((x) => x !== n); changed("refs"); });
        refBoxes.append(h("label", {class: "small", style: {color: PHASE[n] || "inherit", fontWeight: 600}}, c, " " + n));
      }
    }

    // ---------------------------------------------------------------- rows
    async function load() {
      save();
      if (!st.source) return;
      dinfo.textContent = "reading…";
      const q = "/api/xrd/rows?source=" + encodeURIComponent(st.source) + (st.T ? "&T=" + encodeURIComponent(st.T) : "") +
        (st.comp ? "&comp=" + encodeURIComponent(st.comp) : "");
      const r = await S.api(q);
      if (!r.ok) { dinfo.textContent = r.error || "failed"; return; }
      D = r;
      st.comp = r.comp; st.T = r.T || "";
      if (st.spAuto && r.host) {             // species defaults from the run (no element is assumed)
        st.spTo = r.host; st.spMap = (r.dopant ? r.dopant + ":" : "") + r.host;
        spTo.value = st.spTo; spMap.value = st.spMap; save();
      }
      tsel.innerHTML = ""; (r.temperatures.length ? r.temperatures : [""]).forEach((t) => tsel.append(h("option", {value: t}, t ? t + " K" : "–")));
      tsel.value = st.T;
      csel.innerHTML = ""; r.compositions.forEach((c) => csel.append(h("option", {value: c.dir}, `${c.dir} (${c.n})`)));
      csel.value = st.comp;
      if (!st.outDir) { st.outDir = r.workdir + "/xrd_average"; outDir.input.value = st.outDir; }
      D.phases = [...new Set(r.rows.map((x) => x.p))];
      const pOrder = ["parent", "c", "t", "oIII", "m"];
      D.phases.sort((a, b) => ((pOrder.indexOf(a) + 1 || 9) - (pOrder.indexOf(b) + 1 || 9)) || a.localeCompare(b));
      if (!st.phases || !st.phases.some((p) => D.phases.includes(p))) st.phases = D.phases.slice();
      st.phases = st.phases.filter((p) => D.phases.includes(p));
      for (const x of r.rows) x.dF = x.F == null || r.F_ref == null ? null : 1000 * (x.F - r.F_ref);
      const nok = r.rows.filter((x) => x.st === "ok").length, nfl = r.rows.filter((x) => x.st === "ok" && (x.im || x.un)).length;
      dinfo.textContent = `${r.layout} · ${r.workdir} · ${r.rows.length} structures in ${r.comp}, ${nok} with F, ${nfl} flagged · F_ref = ${S.fmt(r.F_ref, 6)} eV/atom (lowest unflagged)`;
      save();
      drawPF();
      sched(true);
    }
    const color = (p) => PHASE[p] || EXTRA[Math.max(0, (D ? D.phases.indexOf(p) : 0)) % EXTRA.length];
    const usable = (x) => x.st === "ok" && x.dF != null && (st.flagged || !(x.im || x.un));

    // ---------------------------------------------------------------- P(F) plot
    const PW = 900, PH = 280, PM = {l: 58, r: 14, t: 12, b: 42};
    let pfx = null;                   // x scale of the P(F) plot
    function drawPF() {
      pfBox.innerHTML = "";
      chips.innerHTML = "";
      gSel.innerHTML = "";
      if (!D) return;
      const ok = D.rows.filter((x) => x.st === "ok" && x.dF != null);
      if (!ok.length) { pfBox.append(h("div", {class: "muted"}, "no free energies in this composition")); return; }
      for (const p of D.phases) {
        const n = ok.filter((x) => x.p === p).length;
        const on = st.phases.includes(p);
        chips.append(h("button", {class: "xrd-chip" + (on ? " on" : ""), style: `--c: ${color(p)}`, title: "include / leave out this phase",
          onclick: () => { st.phases = on ? st.phases.filter((x) => x !== p) : [...st.phases, p]; save(); drawPF(); sched(); }},
          h("i", {}), `${p} (${n})`));
      }
      let lo = Math.min(...ok.map((x) => x.dF)), hi = Math.max(...ok.map((x) => x.dF));
      const pad = 0.04 * (hi - lo || 1); lo -= pad; hi += pad;
      if (st.xlo !== "" && isFinite(Number(st.xlo))) lo = Number(st.xlo);
      if (st.xhi !== "" && isFinite(Number(st.xhi))) hi = Number(st.xhi);
      if (hi <= lo) hi = lo + 1;
      const nb = st.bins, bw = (hi - lo) / nb;
      const inr = ok.filter((x) => x.dF >= lo && x.dF <= hi);
      const Ntot = inr.length || 1;
      const hist = {}, fits = {};
      let ymax = 0;
      for (const p of D.phases) {
        const c = new Array(nb).fill(0);
        for (const x of inr) if (x.p === p) c[Math.min(nb - 1, Math.floor((x.dF - lo) / bw))]++;
        hist[p] = c.map((v) => v / (Ntot * bw));
        ymax = Math.max(ymax, ...hist[p]);
        const vals = ok.filter((x) => x.p === p && (st.flagged || !(x.im || x.un))).map((x) => x.dF);
        fits[p] = vals.length >= 2 ? gmm(vals) : null;
        if (fits[p]) fits[p].scale = vals.length / Ntot;
      }
      const X0 = PM.l, X1 = PW - PM.r, Y0 = PH - PM.b, Y1 = PM.t;
      const fx = (v) => X0 + (v - lo) / (hi - lo) * (X1 - X0);
      const fy = (v) => Y0 - v / (ymax * 1.08 || 1) * (Y0 - Y1 - 10);
      pfx = {fx, inv: (px) => lo + (px - X0) / (X1 - X0) * (hi - lo), X0, X1, lo, hi};
      const col = cssCols();
      const svg = sv("svg", {viewBox: `0 0 ${PW} ${PH}`, class: "xrd-svg", "font-family": "system-ui, sans-serif", "font-size": "11"});
      svg.append(sv("rect", {x: 0, y: 0, width: PW, height: PH, fill: col.panel}));
      const gAxes = axes(svg, {X0, X1, Y0, Y1, lo, hi, ymin: 0, ymax: ymax * 1.08, fyTicks: fy, col,
        xlabel: "F − F_ref  (meV/atom)", ylabel: "P(F)  (meV⁻¹)"});
      svg.append(gAxes);
      for (const p of D.phases) {
        const on = st.phases.includes(p), c = color(p);
        let d = `M${fx(lo)},${Y0}`;
        hist[p].forEach((v, i) => { d += `L${fx(lo + i * bw)},${fy(v)}L${fx(lo + (i + 1) * bw)},${fy(v)}`; });
        d += `L${fx(hi)},${Y0}`;
        svg.append(sv("path", {d, fill: c, "fill-opacity": on ? 0.2 : 0.06, stroke: c, "stroke-opacity": on ? 0.9 : 0.3, "stroke-width": 0.8}));
        const f = fits[p];
        if (f && st.fits) {
          let pd = "";
          for (let i = 0; i <= 300; i++) { const xv = lo + (hi - lo) * i / 300; pd += (i ? "L" : "M") + fx(xv) + "," + fy(f.scale * mix(xv, f)); }
          svg.append(sv("path", {d: pd, fill: "none", stroke: c, "stroke-width": on ? 1.6 : 0.8, "stroke-opacity": on ? 1 : 0.35}));
          if (f.k > 1) f.mu.forEach((m, j) => {
            let cd = "";
            for (let i = 0; i <= 200; i++) { const xv = lo + (hi - lo) * i / 200; cd += (i ? "L" : "M") + fx(xv) + "," + fy(f.scale * f.w[j] * gpdf(xv, m, f.sig[j])); }
            svg.append(sv("path", {d: cd, fill: "none", stroke: c, "stroke-width": 0.7, "stroke-dasharray": "3 2", "stroke-opacity": on ? 0.9 : 0.3}));
          });
        }
        if (f) f.mu.forEach((m, j) => gSel.append(h("option", {value: `${p}|${j}`},
          `${p} · ${j + 1}/${f.k}: μ ${m.toFixed(2)}, σ ${f.sig[j].toFixed(2)}, w ${f.w[j].toFixed(2)}`)));
      }
      gSel.prepend(h("option", {value: ""}, "pick a Gaussian → window μ ± kσ"));
      gSel.value = "";
      // rug: one tick per structure, flagged ones red
      for (const x of ok) {
        if (x.dF < lo || x.dF > hi) continue;
        const fl = x.im || x.un;
        svg.append(sv("line", {x1: fx(x.dF), x2: fx(x.dF), y1: Y0 - (fl ? 9 : 6), y2: Y0, stroke: fl ? "#d93025" : color(x.p),
          "stroke-opacity": st.phases.includes(x.p) ? 0.8 : 0.25, "stroke-width": fl ? 1.2 : 0.8}));
      }
      const gS = sv("g", {class: "xrd-sel"});
      svg.append(gS);
      svg.append(sv("rect", {x: X0, y: Y1, width: X1 - X0, height: Y0 - Y1, fill: "transparent", class: "xrd-hit", style: "cursor: crosshair"}));
      pfBox.append(svg);
      pfx.svg = svg; pfx.g = gS; pfx.Y0 = Y0; pfx.Y1 = Y1;
      pfx.fits = fits;
      dragPF(svg);
      drawSel();
    }
    function pickGauss(v) {
      const [p, j] = v.split("|");
      const f = pfx && pfx.fits[p];
      if (!f) return;
      st.phases = [p];
      st.lo = f.mu[+j] - st.ksig * f.sig[+j]; st.hi = f.mu[+j] + st.ksig * f.sig[+j];
      save(); drawPF(); gSel.value = v; sched();
    }
    function drawSel() {
      if (!D) return;
      loIn.value = st.lo == null ? "" : st.lo.toFixed(3);
      hiIn.value = st.hi == null ? "" : st.hi.toFixed(3);
      absTxt.textContent = st.lo == null || D.F_ref == null ? "(no window: the whole range)" :
        `= ${(D.F_ref + st.lo / 1000).toFixed(6)} … ${(D.F_ref + st.hi / 1000).toFixed(6)} eV/atom`;
      if (pfx) {
        pfx.g.innerHTML = "";
        if (st.lo != null) {
          const a = Math.max(pfx.X0, Math.min(pfx.X1, pfx.fx(st.lo))), b = Math.max(pfx.X0, Math.min(pfx.X1, pfx.fx(st.hi)));
          const ac = cssCols().acc;
          pfx.g.append(sv("rect", {x: a, y: pfx.Y1, width: Math.max(0, b - a), height: pfx.Y0 - pfx.Y1, fill: ac, "fill-opacity": 0.12}),
            sv("line", {x1: a, x2: a, y1: pfx.Y1, y2: pfx.Y0, stroke: ac, "stroke-width": 1.5}),
            sv("line", {x1: b, x2: b, y1: pfx.Y1, y2: pfx.Y0, stroke: ac, "stroke-width": 1.5}),
            sv("rect", {x: a - 4, y: pfx.Y1 + (pfx.Y0 - pfx.Y1) / 2 - 12, width: 8, height: 24, rx: 2, fill: ac}),
            sv("rect", {x: b - 4, y: pfx.Y1 + (pfx.Y0 - pfx.Y1) / 2 - 12, width: 8, height: 24, rx: 2, fill: ac}));
        }
      }
      const s = selection();
      const by = {};
      s.sel.forEach((x) => { by[x.p] = (by[x.p] || 0) + 1; });
      selTxt.innerHTML = "";
      selTxt.append(h("b", {}, `${s.sel.length} structures selected`),
        s.sel.length ? " (" + Object.entries(by).map(([p, n]) => `${p} ${n}`).join(", ") + ")" : "",
        s.flagged ? ` · ${s.flagged} flagged left out` : "",
        s.sel.length > st.max ? h("span", {class: "pill warn", style: {marginLeft: "6px"}}, `random ${st.max} of ${s.sel.length} will be used`) : "",
        !st.phases.length ? h("span", {class: "pill warn", style: {marginLeft: "6px"}}, "no phase chosen") : "");
    }
    function selection() {
      const out = {sel: [], flagged: 0};
      if (!D) return out;
      const Fl = st.lo == null ? -Infinity : D.F_ref + st.lo / 1000, Fh = st.hi == null ? Infinity : D.F_ref + st.hi / 1000;
      for (const x of D.rows) {
        if (x.st !== "ok" || x.F == null || !st.phases.includes(x.p)) continue;
        if (!(x.F >= Fl - 1e-12 && x.F <= Fh + 1e-12)) continue;
        if ((x.im || x.un) && !st.flagged) { out.flagged++; continue; }
        out.sel.push(x);
      }
      return out;
    }
    function dragPF(svg) {
      let mode = null, x0 = 0, l0 = 0, h0 = 0;
      const toX = (e) => { const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY; return pt.matrixTransform(svg.getScreenCTM().inverse()).x; };
      svg.addEventListener("pointerdown", (e) => {
        const px = toX(e);
        if (px < pfx.X0 - 6 || px > pfx.X1 + 6) return;
        const v = pfx.inv(px);
        const near = (b) => b != null && Math.abs(pfx.fx(b) - px) < 7;
        if (near(st.lo)) mode = "lo"; else if (near(st.hi)) mode = "hi";
        else if (st.lo != null && v > st.lo && v < st.hi) { mode = "move"; x0 = v; l0 = st.lo; h0 = st.hi; }
        else { mode = "new"; x0 = v; st.lo = st.hi = v; }
        svg.setPointerCapture(e.pointerId);
        e.preventDefault();
      });
      svg.addEventListener("pointermove", (e) => {
        if (!mode) return;
        const v = Math.max(pfx.lo, Math.min(pfx.hi, pfx.inv(toX(e))));
        if (mode === "lo") st.lo = v; else if (mode === "hi") st.hi = v;
        else if (mode === "move") { st.lo = l0 + v - x0; st.hi = h0 + v - x0; }
        else { st.lo = Math.min(x0, v); st.hi = Math.max(x0, v); }
        if (st.lo > st.hi) { const t = st.lo; st.lo = st.hi; st.hi = t; mode = mode === "lo" ? "hi" : mode === "hi" ? "lo" : mode; }
        drawSel();
      });
      const up = () => {
        if (!mode) return;
        if (mode === "new" && Math.abs(pfx.fx(st.hi) - pfx.fx(st.lo)) < 3) { st.lo = st.hi = null; }
        mode = null; save(); drawSel(); sched();
      };
      svg.addEventListener("pointerup", up); svg.addEventListener("pointercancel", up);
    }

    // ---------------------------------------------------------------- compute (job) and result
    function body() {
      return {source: st.source, T: st.T || null, comp: st.comp, phases: st.phases, lo: st.lo, hi: st.hi,
        include_flagged: st.flagged, allow_unrelaxed: st.unrelaxed, max: st.max, seed: st.seed, jobs: st.jobs,
        params: {wavelength: st.wl === "custom" ? String(st.wlc) : st.wl, tt_lo: st.tt0, tt_hi: st.tt1, symprec: st.symprec,
          min_rel: st.minrel, debye_waller: st.dw,
          species: {mode: st.species, to: st.spTo, map: st.spMap, anions: st.anions},
          refs: st.refs, refs_dir: st.refsDir, ref_replace: st.refReplace, ref_scale_volume: st.refScale},
        broadening: {step: st.step, fwhm: st.fwhm, profile: st.profile, eta: st.eta, norm: st.norm,
          caglioti: st.cag ? [st.U, st.V, st.W] : null},
        out_dir: st.outDir};
    }
    const sigOf = (b) => JSON.stringify([b.source, b.T, b.comp, b.phases, b.lo, b.hi, b.include_flagged, b.allow_unrelaxed, b.max, b.seed,
      {...b.params, tt_lo: 0, tt_hi: 0}]);
    function sched(force) {
      clearTimeout(debT);
      if (!st.auto && !force) { status.innerHTML = ""; status.append(h("span", {class: "pill warn"}, "settings changed: press Compute XRD")); return; }
      debT = setTimeout(() => compute(false), 700);
    }
    async function compute(force) {
      if (!D || !st.phases.length) return;
      if (!force && st.lo == null) { setStatus("warn", "drag on P(F) to choose a window (or press Compute XRD for the whole range)"); return; }
      const b = body();
      const sig = sigOf(b);
      if (!force && X && X.sig === sig && st.tt0 >= X.range[0] - 1e-9 && st.tt1 <= X.range[1] + 1e-9) { drawXRD(); return; }
      if (!selection().sel.length && !st.refs.length) { X = null; drawXRD(); return; }
      clearTimeout(pollT);
      b.cancel = job && job.running ? job.id : null;
      setStatus("run", "starting…");
      const r = await S.api("/api/xrd/compute", b);
      if (!r.ok) { setStatus("bad", r.error); return; }
      job = {id: r.job, key: r.key, sig, running: true};
      if (r.cli) cliPre.textContent = r.cli;
      poll(job);
    }
    async function poll(j) {
      if (j !== job) return;
      const r = await S.api("/api/jobs?id=" + j.id + "&tail=6000", undefined, {quiet: true});
      if (j !== job) return;
      if (!r.ok) { pollT = setTimeout(() => poll(j), 1000); return; }
      const jj = r.job, log = jj.log || "";
      const m = [...log.matchAll(/XRD (\d+)\/(\d+)/g)].pop();
      const head = (log.match(/^.*\| T .*$/m) || [""])[0];
      if (jj.state === "running" || jj.state === "starting") {
        setStatus("run", m ? `computing ${m[1]} / ${m[2]}` : "starting (loading pymatgen)…", m ? Number(m[1]) / Number(m[2]) : 0);
        pollT = setTimeout(() => poll(j), 600);
        return;
      }
      j.running = false;
      if (jj.state !== "done") {
        setStatus("bad", `XRD job ${jj.state}`);
        status.append(h("a", {href: "#/jobs?id=" + j.id, style: {marginLeft: "8px"}}, "log"));
        status.append(h("pre", {class: "log", style: {maxHeight: "160px"}}, log.split("\n").slice(-12).join("\n")));
        return;
      }
      const res = await fetch("/api/xrd/sticks?key=" + j.key, {headers: {"X-Requested-With": "studio"}});
      const d = await res.json().catch(() => null);
      if (!d || d.ok === false) { setStatus("bad", (d && d.error) || "result unreadable"); return; }
      X = decode(d, j);
      ycache = null;
      const sel = d.selection || {};
      setStatus("ok", `${X.structures.length} structures` + (sel.no_geometry ? `, ${sel.no_geometry} without relaxed file (skipped)` : "") +
        (Object.keys(d.failed || {}).length ? `, ${Object.keys(d.failed).length} failed` : "") + ` · ${X.wavelength}`);
      if (head) status.title = head;
      drawList();
      drawXRD();
    }
    function setStatus(cls, msg, frac) {
      status.innerHTML = "";
      status.append(h("span", {class: "pill " + cls}, msg));
      if (frac != null) status.append(h("span", {class: "xrd-prog"}, h("i", {style: {width: (100 * frac).toFixed(0) + "%"}})));
    }
    function decode(d, j) {
      const dec = (s) => ({x: Float64Array.from(s.xi, (v) => v / 1e4), y: Float64Array.from(s.yi, (v) => v * s.ymax / 1e6)});
      return {key: j.key, sig: j.sig, range: d.stick_params.range, wavelength: d.stick_params.wavelength, lambda: d.stick_params.wavelength_A,
        replace: d.stick_params.replace, selection: d.selection,
        structures: d.structures.map((s) => ({...dec(s), key: s.key, phase: s.phase, F: s.F, dF: s.dF_meV, how: s.how, path: s.path,
          n: s.n_atoms, formula: s.formula})),
        refs: (d.refs || []).map((r) => ({name: r.name, x: Float64Array.from(r.x), y: Float64Array.from(r.y), hkl: r.hkl, d: r.d,
          formula: r.formula, n: r.n_atoms}))};
    }
    function drawList() {
      listBox.innerHTML = "";
      if (!X) return;
      const t = h("table", {}, h("thead", {}, h("tr", {}, ["structure", "phase", "F − F_ref meV/atom", "species used", "geometry"].map((x) => h("th", {class: "l"}, x)))),
        h("tbody", {}, X.structures.slice(0, 2000).map((s) => h("tr", {}, h("td", {class: "l mono"}, s.key), h("td", {class: "l"}, s.phase),
          h("td", {}, S.fmt(s.dF, 3)), h("td", {class: "l mono"}, s.formula), h("td", {class: "l small"}, s.how)))));
      listBox.append(h("div", {class: "row small"}, h("a", {href: "javascript:void 0", onclick: dlList}, "download list (CSV)")),
        h("div", {class: "tablewrap", style: {maxHeight: "300px"}}, t));
    }

    // ---------------------------------------------------------------- broadening and averaging (same as xrd_average.py)
    function grid() {
      const n = Math.round((st.tt1 - st.tt0) / st.step) + 1;
      const g = new Float64Array(Math.max(2, Math.min(n, 2000000)));
      for (let i = 0; i < g.length; i++) g[i] = st.tt0 + st.step * i;
      return g;
    }
    function fwhmAt(x) {
      if (!st.cag) return st.fwhm;
      const t = Math.tan(x * Math.PI / 360);
      return Math.sqrt(Math.max(st.U * t * t + st.V * t + st.W, 1e-8));
    }
    function broaden(xs, ys, g) {
      const out = new Float64Array(g.length), g0 = g[0], step = g[1] - g[0], win = WINDOW[st.profile], eta = st.eta, prof = st.profile;
      for (let k = 0; k < xs.length; k++) {
        const x0 = xs[k], I = ys[k], w = fwhmAt(x0);
        const i0 = Math.max(0, Math.floor((x0 - win * w - g0) / step)), i1 = Math.min(g.length, Math.ceil((x0 + win * w - g0) / step) + 1);
        for (let i = i0; i < i1; i++) {
          const u = (g[i] - x0) / w, u2 = u * u;
          out[i] += I * (prof === "gaussian" ? Math.exp(-LN2x4 * u2) : prof === "lorentzian" ? 1 / (1 + 4 * u2) :
            eta / (1 + 4 * u2) + (1 - eta) * Math.exp(-LN2x4 * u2));
        }
      }
      return out;
    }
    function averaged() {
      const sig = JSON.stringify([X.key, st.tt0, st.tt1, st.step, st.fwhm, st.profile, st.eta, st.cag, st.U, st.V, st.W, st.norm]);
      if (ycache && ycache.sig === sig) return ycache;
      const g = grid(), n = g.length, N = X.structures.length;
      const mean = new Float64Array(n), m2 = new Float64Array(n), mn = new Float64Array(n).fill(N ? Infinity : 0), mx = new Float64Array(n).fill(N ? -Infinity : 0);
      const indiv = [];
      X.structures.forEach((s, j) => {
        let y = broaden(s.x, s.y, g);
        if (st.norm === "max") { let m = 0; for (const v of y) if (v > m) m = v; if (m > 0) for (let i = 0; i < n; i++) y[i] *= 100 / m; }
        for (let i = 0; i < n; i++) { mean[i] += y[i]; m2[i] += y[i] * y[i]; if (y[i] < mn[i]) mn[i] = y[i]; if (y[i] > mx[i]) mx[i] = y[i]; }
        if (j < 30) indiv.push(y);
      });
      const std = new Float64Array(n);
      for (let i = 0; i < n; i++) { if (N) { mean[i] /= N; std[i] = Math.sqrt(Math.max(0, m2[i] / N - mean[i] * mean[i])); } }
      const refs = X.refs.map((r) => { const y = broaden(r.x, r.y, g); let m = 0; for (const v of y) if (v > m) m = v;
        if (m > 0) for (let i = 0; i < n; i++) y[i] *= 100 / m; return {...r, prof: y}; });
      ycache = {sig, g, mean, std, mn, mx, indiv, refs, N};
      return ycache;
    }

    // ---------------------------------------------------------------- XRD plot
    function drawXRD() {
      xrdBox.innerHTML = "";
      if (!X) { xrdBox.append(h("div", {class: "muted xrd-empty"}, D ? "choose a window and phases on P(F)" : "load a results folder first")); return; }
      const A = averaged(), g = A.g, col = cssCols();
      const [lo, hi] = zoom || [st.tt0, st.tt1];
      const showRefs = st.refMode === "stack" ? A.refs : [];
      const W = 1000, MH = 330, RH = 74, M = {l: 58, r: 14, t: 14, b: 44};
      const ticksH = st.ticks && A.refs.length && st.refMode !== "none" ? 10 * A.refs.length + 6 : 0;
      const H = M.t + MH + ticksH + showRefs.length * RH + M.b;
      const X0 = M.l, X1 = W - M.r;
      const fx = (v) => X0 + (v - lo) / (hi - lo) * (X1 - X0);
      const i0 = Math.max(0, Math.floor((lo - g[0]) / st.step)), i1 = Math.min(g.length - 1, Math.ceil((hi - g[0]) / st.step));
      const tr = st.yscale === "sqrt" ? (v) => Math.sqrt(Math.max(0, v)) : (v) => v;
      let top = 0;
      for (let i = i0; i <= i1; i++) top = Math.max(top, st.spread === "minmax" ? A.mx[i] : st.spread === "std" ? A.mean[i] + A.std[i] : A.mean[i]);
      if (st.indiv) A.indiv.forEach((y) => { for (let i = i0; i <= i1; i++) top = Math.max(top, y[i]); });
      if (st.refMode === "overlay") top = Math.max(top, 100);
      top = top || 1;
      const mY0 = M.t + MH, mY1 = M.t + 8;
      const fy = (v) => mY0 - tr(v) / tr(top * 1.06) * (mY0 - mY1);
      const svg = sv("svg", {viewBox: `0 0 ${W} ${H}`, class: "xrd-svg", "font-family": "system-ui, sans-serif", "font-size": "11"});
      svg.append(sv("rect", {x: 0, y: 0, width: W, height: H, fill: col.panel}));
      const bottom = H - M.b;
      svg.append(axes(svg, {X0, X1, Y0: mY0, Y1: mY1, lo, hi, ymin: 0, ymax: top * 1.06, fyTicks: fy, col, sqrt: st.yscale === "sqrt",
        xAxisY: bottom, xlabel: `2θ (°)   ${X.wavelength}, λ = ${X.lambda.toFixed(5)} Å`,
        ylabel: st.norm === "max" ? "intensity (each pattern max = 100)" : "intensity / n_atoms²"}));
      const clipId = "xc" + Math.random().toString(36).slice(2, 8);
      svg.append(sv("defs", {}, sv("clipPath", {id: clipId}, sv("rect", {x: X0, y: 0, width: X1 - X0, height: H}))));
      const gp = sv("g", {"clip-path": `url(#${clipId})`});
      svg.append(gp);
      const line = (y, f) => { let d = ""; const stp = Math.max(1, Math.floor((i1 - i0) / 4000)); for (let i = i0; i <= i1; i += stp) d += (d ? "L" : "M") + fx(g[i]).toFixed(2) + "," + f(y[i]).toFixed(2); return d; };
      const selCol = st.phases.length === 1 ? color(st.phases[0]) : col.ink;
      if (st.spread !== "none" && A.N > 1) {
        const up = st.spread === "std" ? (i) => A.mean[i] + A.std[i] : (i) => A.mx[i];
        const dn = st.spread === "std" ? (i) => Math.max(0, A.mean[i] - A.std[i]) : (i) => A.mn[i];
        let d = "";
        const stp = Math.max(1, Math.floor((i1 - i0) / 3000));
        for (let i = i0; i <= i1; i += stp) d += (d ? "L" : "M") + fx(g[i]).toFixed(2) + "," + fy(up(i)).toFixed(2);
        for (let i = i1; i >= i0; i -= stp) d += "L" + fx(g[i]).toFixed(2) + "," + fy(dn(i)).toFixed(2);
        gp.append(sv("path", {d: d + "Z", fill: selCol, "fill-opacity": 0.18, stroke: "none"}));
      }
      if (st.indiv) A.indiv.forEach((y) => gp.append(sv("path", {d: line(y, fy), fill: "none", stroke: col.muted, "stroke-width": 0.4, "stroke-opacity": 0.5})));
      if (st.refMode === "overlay") A.refs.forEach((r) => gp.append(sv("path", {d: line(r.prof, fy), fill: "none", stroke: color(r.name), "stroke-width": 0.9, "stroke-dasharray": "4 2"})));
      gp.append(sv("path", {d: line(A.mean, fy), fill: "none", stroke: selCol, "stroke-width": 1.3}));
      // peak labels of the mean
      if (st.peakLabels && A.N) {
        const pk = [];
        for (let i = Math.max(i0, 2); i <= Math.min(i1, g.length - 3); i++) {
          const v = A.mean[i];
          if (v > A.mean[i - 1] && v >= A.mean[i + 1] && v > A.mean[i - 2] && v >= A.mean[i + 2] && v > 0.08 * top) pk.push(i);
        }
        pk.sort((a, b) => A.mean[b] - A.mean[a]);
        const used = [];
        for (const i of pk.slice(0, 14)) {
          const px = fx(g[i]);
          if (used.some((u) => Math.abs(u - px) < 26)) continue;
          used.push(px);
          gp.append(sv("text", {x: px, y: fy(A.mean[i]) - 4, "text-anchor": "middle", fill: col.muted, "font-size": 9}, g[i].toFixed(2)));
        }
      }
      // Bragg ticks of the references under the main panel
      let yT = mY0 + 6;
      if (ticksH) A.refs.forEach((r) => {
        const m = Math.max(...r.y, 0);
        for (let k = 0; k < r.x.length; k++) if (r.x[k] >= lo && r.x[k] <= hi && r.y[k] >= 0.005 * m)
          gp.append(sv("line", {x1: fx(r.x[k]), x2: fx(r.x[k]), y1: yT, y2: yT + 7, stroke: color(r.name), "stroke-width": 1}));
        svg.append(sv("text", {x: X0 - 4, y: yT + 7, "text-anchor": "end", fill: color(r.name), "font-size": 9}, r.name));
        yT += 10;
      });
      // stacked reference profiles
      let yR = mY0 + ticksH;
      showRefs.forEach((r) => {
        const rY0 = yR + RH - 6, rY1 = yR + 10;
        const fr = (v) => rY0 - tr(v) / tr(105) * (rY0 - rY1);
        svg.append(sv("line", {x1: X0, x2: X1, y1: rY0, y2: rY0, stroke: col.line}));
        gp.append(sv("path", {d: line(r.prof, fr), fill: "none", stroke: color(r.name), "stroke-width": 1}));
        svg.append(sv("text", {x: X1 - 4, y: yR + 20, "text-anchor": "end", fill: color(r.name), "font-weight": 600},
          `${r.name}  (${r.formula}${X.replace && st.refReplace ? ", species changed" : ""})`));
        if (st.refLabels) {
          const m = Math.max(...r.y, 0), idx = [];
          for (let k = 0; k < r.x.length; k++) if (r.x[k] >= lo && r.x[k] <= hi && r.y[k] >= 0.05 * m) idx.push(k);
          idx.sort((a, b) => r.y[b] - r.y[a]);
          const used = [];
          for (const k of idx.slice(0, 10)) {
            const px = fx(r.x[k]);
            if (used.some((u) => Math.abs(u - px) < 34)) continue;
            used.push(px);
            gp.append(sv("text", {x: px, y: Math.max(rY1 + 2, fr(100 * r.y[k] / m) - 3), "text-anchor": "middle", fill: col.muted, "font-size": 8.5},
              r.hkl[k].replace(/x\d+$/, "")));
          }
        }
        yR += RH;
      });
      // legend
      const leg = sv("g", {});
      const lg = [[selCol, `mean of ${A.N}` + (st.phases.length ? ` (${st.phases.join(", ")}, ${winText()})` : ""), 1.3]];
      if (st.spread !== "none" && A.N > 1) lg.push([selCol, st.spread === "std" ? "± 1 std" : "min – max", 0]);
      if (st.refMode === "overlay") A.refs.forEach((r) => lg.push([color(r.name), r.name + " (" + r.formula + ")", -1]));
      lg.forEach(([c, t, w], k) => {
        const y = M.t + 10 + 14 * k;
        leg.append(w > 0 ? sv("line", {x1: X1 - 250, x2: X1 - 228, y1: y, y2: y, stroke: c, "stroke-width": 1.6}) :
          w === 0 ? sv("rect", {x: X1 - 250, y: y - 5, width: 22, height: 10, fill: c, "fill-opacity": 0.2}) :
            sv("line", {x1: X1 - 250, x2: X1 - 228, y1: y, y2: y, stroke: c, "stroke-width": 1, "stroke-dasharray": "4 2"}),
        sv("text", {x: X1 - 222, y: y + 4, fill: col.ink}, t));
      });
      if (X.replace) {
        const rp = X.replace;
        leg.append(sv("text", {x: X1 - 250, y: M.t + 14 + 14 * lg.length, fill: col.muted, "font-size": 10},
          "XRD species: " + [...Object.entries(rp.map || {}).map(([a, b]) => `${a}→${b}`), ...(rp.cations_to ? [`cations→${rp.cations_to}`] : [])].join(", ") + " (in memory)"));
      }
      svg.append(leg);
      // hover + zoom
      const hl = sv("line", {y1: M.t, y2: bottom, stroke: col.muted, "stroke-dasharray": "2 2", visibility: "hidden"});
      const zr = sv("rect", {y: M.t, height: bottom - M.t, fill: col.acc, "fill-opacity": 0.12, visibility: "hidden"});
      svg.append(hl, zr, sv("rect", {x: X0, y: M.t, width: X1 - X0, height: bottom - M.t, fill: "transparent", style: "cursor: crosshair"}));
      const toX = (e) => { const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY; return pt.matrixTransform(svg.getScreenCTM().inverse()).x; };
      const inv = (px) => lo + (px - X0) / (X1 - X0) * (hi - lo);
      let z0 = null;
      svg.addEventListener("pointerdown", (e) => { const px = toX(e); if (px < X0 || px > X1) return; z0 = px; svg.setPointerCapture(e.pointerId); });
      svg.addEventListener("pointermove", (e) => {
        const px = Math.max(X0, Math.min(X1, toX(e)));
        hl.setAttribute("x1", px); hl.setAttribute("x2", px); hl.setAttribute("visibility", "visible");
        if (z0 != null) { zr.setAttribute("x", Math.min(z0, px)); zr.setAttribute("width", Math.abs(px - z0)); zr.setAttribute("visibility", "visible"); }
        const tt = inv(px), i = Math.max(0, Math.min(g.length - 1, Math.round((tt - g[0]) / st.step)));
        const dsp = X.lambda / (2 * Math.sin(tt * Math.PI / 360));
        let near = "";
        A.refs.forEach((r) => {
          let best = -1;
          const m = Math.max(...r.y, 0);
          for (let k = 0; k < r.x.length; k++) if (r.y[k] >= 0.005 * m && Math.abs(r.x[k] - tt) < 0.6 && (best < 0 || Math.abs(r.x[k] - tt) < Math.abs(r.x[best] - tt))) best = k;
          if (best >= 0) near += `  ${r.name} ${r.hkl[best]} ${r.x[best].toFixed(2)}° (${(100 * r.y[best] / m).toFixed(1)})`;
        });
        readout.textContent = `2θ ${tt.toFixed(3)}°  d ${dsp.toFixed(4)} Å  mean ${S.fmt(A.mean[i], 2)}` +
          (A.N > 1 ? ` ± ${S.fmt(A.std[i], 2)} [${S.fmt(A.mn[i], 1)}…${S.fmt(A.mx[i], 1)}]` : "") + (near ? "  |" + near : "");
      });
      svg.addEventListener("pointerleave", () => hl.setAttribute("visibility", "hidden"));
      svg.addEventListener("pointerup", (e) => {
        if (z0 == null) return;
        const px = Math.max(X0, Math.min(X1, toX(e)));
        if (Math.abs(px - z0) > 6) zoom = [inv(Math.min(z0, px)), inv(Math.max(z0, px))];
        z0 = null; zr.setAttribute("visibility", "hidden");
        if (zoom) drawXRD();
      });
      svg.addEventListener("dblclick", () => { zoom = null; drawXRD(); });
      xrdBox.append(svg);
    }
    const winText = () => st.lo == null ? "whole F range" : `${st.lo.toFixed(2)}…${st.hi.toFixed(2)} meV/atom`;

    // ---------------------------------------------------------------- shared drawing helpers
    function cssCols() {
      const cs = getComputedStyle(document.documentElement);
      const v = (n, d) => (cs.getPropertyValue(n) || d).trim() || d;
      return {ink: v("--ink", "#1d2330"), muted: v("--muted", "#5f6b7a"), line: v("--line", "#dde2e8"), panel: v("--panel", "#fff"), acc: v("--acc", "#2f5fb3")};
    }
    function niceTicks(a, b, n) {
      const span = b - a, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
      const st = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
      const out = [];
      for (let v = Math.ceil(a / st) * st; v <= b + 1e-9 * span; v += st) out.push(Math.abs(v) < 1e-12 * span ? 0 : v);
      return out;
    }
    function axes(svg, o) {
      const g = sv("g", {});
      const xY = o.xAxisY || o.Y0;
      g.append(sv("line", {x1: o.X0, x2: o.X1, y1: xY, y2: xY, stroke: o.col.ink, "stroke-width": 0.8}),
        sv("line", {x1: o.X0, x2: o.X0, y1: o.Y1 - 6, y2: o.Y0, stroke: o.col.ink, "stroke-width": 0.8}));
      const fx = (v) => o.X0 + (v - o.lo) / (o.hi - o.lo) * (o.X1 - o.X0);
      for (const t of niceTicks(o.lo, o.hi, 10)) {
        g.append(sv("line", {x1: fx(t), x2: fx(t), y1: xY, y2: xY + 4, stroke: o.col.ink}),
          sv("text", {x: fx(t), y: xY + 16, "text-anchor": "middle", fill: o.col.ink}, +t.toFixed(6)));
      }
      const yt = o.sqrt ? [0, 0.01, 0.05, 0.2, 0.5, 1].map((f) => f * o.ymax).filter((v) => v <= o.ymax) : niceTicks(0, o.ymax, 4);
      for (const t of yt) {
        const y = o.fyTicks(t);
        g.append(sv("line", {x1: o.X0 - 4, x2: o.X0, y1: y, y2: y, stroke: o.col.ink}),
          sv("text", {x: o.X0 - 7, y: y + 4, "text-anchor": "end", fill: o.col.ink}, t >= 100 || t === 0 ? Math.round(t) : +t.toPrecision(2)));
      }
      g.append(sv("text", {x: (o.X0 + o.X1) / 2, y: xY + 34, "text-anchor": "middle", fill: o.col.ink}, o.xlabel),
        sv("text", {x: 14, y: (o.Y0 + o.Y1) / 2, transform: `rotate(-90 14 ${(o.Y0 + o.Y1) / 2})`, "text-anchor": "middle", fill: o.col.ink}, o.ylabel));
      return g;
    }

    // ---------------------------------------------------------------- downloads and export job
    function download(name, blob) {
      const a = h("a", {href: URL.createObjectURL(blob), download: name});
      document.body.append(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    }
    const stem = () => `xrd_${st.comp}_${st.phases.join("-") || "none"}` + (st.lo == null ? "" : `_${st.lo.toFixed(2)}-${st.hi.toFixed(2)}meV`);
    function dlCSV() {
      if (!X) return S.toast("nothing computed yet");
      const A = averaged();
      const unit = st.norm === "max" ? "rel_max100" : "per_atom2";
      const rows = [["two_theta_deg", `mean_${unit}`, "std", "min", "max", "n", ...A.refs.map((r) => `ref_${r.name}_max100`)].join(",")];
      for (let i = 0; i < A.g.length; i++)
        rows.push([A.g[i].toFixed(5), A.mean[i].toPrecision(6), A.std[i].toPrecision(6), A.mn[i].toPrecision(6), A.mx[i].toPrecision(6), A.N,
          ...A.refs.map((r) => r.prof[i].toPrecision(6))].join(","));
      download(stem() + ".csv", new Blob([rows.join("\n") + "\n"], {type: "text/csv"}));
    }
    function dlList() {
      const rows = [["key", "phase", "F_total_eV_atom", "dF_meV_atom", "geometry", "path", "formula_used"].join(",")];
      X.structures.forEach((s) => rows.push([s.key, s.phase, s.F, S.fmt(s.dF, 4), JSON.stringify(s.how), JSON.stringify(s.path), s.formula].join(",")));
      download(stem() + "_structures.csv", new Blob([rows.join("\n") + "\n"], {type: "text/csv"}));
    }
    function dlSVG(png) {
      const svg = xrdBox.querySelector("svg");
      if (!svg) return S.toast("nothing to save yet");
      const c = svg.cloneNode(true);
      c.setAttribute("xmlns", NS);
      const vb = svg.viewBox.baseVal;
      c.setAttribute("width", vb.width); c.setAttribute("height", vb.height);
      const txt = new XMLSerializer().serializeToString(c);
      if (!png) return download(stem() + ".svg", new Blob([txt], {type: "image/svg+xml"}));
      const img = new Image();
      img.onload = () => {
        const cv = document.createElement("canvas");
        cv.width = vb.width * 2; cv.height = vb.height * 2;
        const ctx = cv.getContext("2d"); ctx.scale(2, 2); ctx.drawImage(img, 0, 0);
        cv.toBlob((b) => download(stem() + ".png", b), "image/png");
      };
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(txt);
    }
    async function exportRun() {
      st.outDir = outDir.input.value.trim(); save();
      const r = await S.api("/api/xrd/export", body());
      if (!r.ok) return;
      exportBox.prepend(S.jobView(r.job, {title: "XRD files → " + st.outDir}));
    }
  }

  // ---------------------------------------------------------------- Gaussian mixture (EM, 1-3 components by BIC)
  const gpdf = (x, m, s) => Math.exp(-0.5 * ((x - m) / s) ** 2) / (s * Math.sqrt(2 * Math.PI));
  const mix = (x, f) => f.w.reduce((a, w, j) => a + w * gpdf(x, f.mu[j], f.sig[j]), 0);
  function em(x, k, mu0, sd) {
    const n = x.length;
    let w = new Array(k).fill(1 / k), mu = mu0.slice(), sig = new Array(k).fill(sd / k);
    const floor = 1e-3 * sd + 1e-12;
    let ll = -Infinity;
    const r = Array.from({length: k}, () => new Float64Array(n));
    for (let it = 0; it < 500; it++) {
      let L = 0;
      for (let i = 0; i < n; i++) {
        let mxl = -Infinity;
        const lp = new Array(k);
        for (let j = 0; j < k; j++) { lp[j] = Math.log(w[j]) - Math.log(sig[j]) - 0.5 * ((x[i] - mu[j]) / sig[j]) ** 2; if (lp[j] > mxl) mxl = lp[j]; }
        let s = 0;
        for (let j = 0; j < k; j++) s += Math.exp(lp[j] - mxl);
        L += mxl + Math.log(s) - 0.5 * Math.log(2 * Math.PI);
        for (let j = 0; j < k; j++) r[j][i] = Math.exp(lp[j] - mxl) / s;
      }
      for (let j = 0; j < k; j++) {
        let nj = 0, m = 0;
        for (let i = 0; i < n; i++) { nj += r[j][i]; m += r[j][i] * x[i]; }
        nj = Math.max(nj, 1e-9); m /= nj;
        let v = 0;
        for (let i = 0; i < n; i++) v += r[j][i] * (x[i] - m) ** 2;
        w[j] = nj / n; mu[j] = m; sig[j] = Math.max(Math.sqrt(v / nj), floor);
      }
      if (Math.abs(L - ll) < 1e-9 * Math.max(1, Math.abs(L))) { ll = L; break; }
      ll = L;
    }
    return {ll, w, mu, sig};
  }
  function gmm(vals) {
    const x = vals.slice().sort((a, b) => a - b), n = x.length;
    const mean = x.reduce((a, b) => a + b, 0) / n;
    const sd = Math.sqrt(x.reduce((a, b) => a + (b - mean) ** 2, 0) / n);
    if (n < 10 || sd < 1e-9) return {k: 1, w: [1], mu: [mean], sig: [Math.max(sd, 1e-3)]};
    let best = null;
    for (let k = 1; k <= 3; k++) {
      if (n < 10 * k) break;
      const inits = [Array.from({length: k}, (_, j) => x[Math.floor((j + 0.5) / k * n)]),
        Array.from({length: k}, (_, j) => x[0] + (j + 1) / (k + 1) * (x[n - 1] - x[0])),
        Array.from({length: k}, (_, j) => x[Math.floor((j + 0.2) / k * n)])];
      let fit = null;
      for (const i0 of inits) { const f = em(x, k, i0, sd); if (!fit || f.ll > fit.ll) fit = f; }
      fit.k = k;
      fit.bic = -2 * fit.ll + (3 * k - 1) * Math.log(n);
      if (!best || fit.bic < best.bic - 10) best = fit;       // a further Gaussian only for strong evidence
    }
    const o = best.mu.map((m, j) => j).sort((a, b) => best.mu[a] - best.mu[b]);
    return {k: best.k, w: o.map((j) => best.w[j]), mu: o.map((j) => best.mu[j]), sig: o.map((j) => best.sig[j])};
  }
})();
