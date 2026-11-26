/* Fluorite Ordering Studio - core: API, router, file picker, option forms built from argparse, job viewer.
   Pages live in static/pages/*.js and call Studio.registerPage({id, title, order, render(el, params)}).
   Drop a new file there to add a page; nothing else needs to change. */
(function () {
  "use strict";
  const S = (window.Studio = {pages: {}, info: null, modules: null});
  // pages may pass arrays, null or false to append(): flatten and skip them (only inside this app)
  const nativeAppend = Element.prototype.append;
  Element.prototype.append = function (...kids) {
    return nativeAppend.apply(this, kids.flat(9).filter((c) => c != null && c !== false));
  };

  // ------------------------------------------------------------------ helpers
  S.h = function (tag, attrs, ...kids) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v == null || v === false || (k === "class" && v === "")) continue;
      if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else if (k === "html") e.innerHTML = v;
      else if (k === "style" && typeof v === "object") Object.assign(e.style, v);
      else e.setAttribute(k, v === true ? "" : v);
    }
    for (const c of kids.flat(9)) if (c != null && c !== false) e.append(c.nodeType ? c : document.createTextNode(String(c)));
    return e;
  };
  const h = S.h;
  S.esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c]));
  S.fmt = (v, d = 2) => (v == null || !isFinite(v) ? "–" : Number(v).toFixed(d));
  S.size = (b) => (b == null ? "" : b < 1024 ? b + " B" : b < 1048576 ? (b / 1024).toFixed(1) + " kB" : (b / 1048576).toFixed(1) + " MB");
  S.ago = (t) => { const s = Date.now() / 1000 - t; return s < 60 ? Math.round(s) + " s ago" : s < 3600 ? Math.round(s / 60) + " min ago" : s < 86400 ? (s / 3600).toFixed(1) + " h ago" : new Date(t * 1000).toLocaleString(); };
  S.dur = (s) => { s = Math.round(s); return s < 60 ? s + " s" : s < 3600 ? Math.floor(s / 60) + " min " + (s % 60) + " s" : Math.floor(s / 3600) + " h " + Math.floor((s % 3600) / 60) + " min"; };
  S.basename = (p) => String(p).replace(/\/+$/, "").split("/").pop();
  S.dirname = (p) => String(p).replace(/\/[^/]*$/, "") || "/";
  // icons from the sprite in index.html; pills = icon + text (never colour alone)
  S.ic = (n, c = "") => `<svg class="ic ${c}" aria-hidden="true"><use href="#i-${n}"/></svg>`;
  S.icon = (n, c = "") => { const t = document.createElement("span"); t.innerHTML = S.ic(n, c); return t.firstChild; };
  const PILLS = {ok: ["p-ok", "check"], warn: ["p-warn", "alert"], bad: ["p-bad", "xcircle"], run: ["p-run", "spin"],
    idle: ["p-idle", "circle"], clock: ["p-bad", "clock"], acc: ["p-acc", "pulse"]};
  S.pill = (kind, text, title) => { const [c, i] = PILLS[kind] || PILLS.idle;
    return h("span", {class: "pill " + c, title: title || null}, S.icon(i), text); };
  S.int = (v) => (v == null || !isFinite(v) ? "–" : Math.round(v).toLocaleString("en-US"));
  S.hms = (s) => { if (s == null || !isFinite(s)) return "–"; s = Math.max(0, Math.round(s));
    return `${Math.floor(s / 3600)}:${String(Math.floor(s % 3600 / 60)).padStart(2, "0")}:${String(s % 60).padStart(2, "0")}`; };
  S.toast = function (msg, ms = 3800, kind) {
    const box = document.getElementById("toasts") || document.body;
    const bad = kind === "bad" || (kind == null && /fail|error|not |refus|missing|wrong|cannot|could not/i.test(String(msg)));
    const t = h("div", {class: "toast" + (bad ? " bad" : "")}, S.icon(bad ? "alert" : "check"), h("div", {}, msg));
    box.append(t);
    setTimeout(() => t.remove(), ms);
  };
  // a confirm dialog for steps that submit, cancel or replace something.  Resolves true / false.
  S.confirm = function ({title, body, ok = "Continue", danger = false, icon = "alert"}) {
    return new Promise((resolve) => {
      const prev = document.activeElement;
      const done = (v) => { d.remove(); document.removeEventListener("keydown", key, true); if (prev && prev.focus) prev.focus(); resolve(v); };
      const okB = h("button", {class: "btn " + (danger ? "danger solid" : "primary"), onclick: () => done(true)}, ok);
      const noB = h("button", {class: "btn ghost", onclick: () => done(false)}, "Cancel");
      const d = h("div", {class: "dlg", role: "alertdialog", "aria-modal": "true", "aria-labelledby": "dlgT",
        onclick: (e) => { if (e.target === d) done(false); }},
        h("div", {class: "box"}, h("h2", {id: "dlgT"}, S.icon(icon, danger ? "bad" : "warn"), title),
          h("div", {class: "body"}, body || ""), h("div", {class: "btns"}, noB, okB)));
      const key = (e) => {
        if (e.key === "Escape") { e.preventDefault(); done(false); }
        if (e.key === "Tab") { e.preventDefault(); (document.activeElement === okB ? noB : okB).focus(); }
      };
      document.addEventListener("keydown", key, true);
      document.body.append(d);
      setTimeout(() => (danger ? noB : okB).focus(), 30);
    });
  };
  S.setTheme = function (t) {
    document.documentElement.dataset.theme = t;
    S.store.set("theme", t);
    const b = document.getElementById("themeBtn");
    if (b) { b.innerHTML = S.ic(t === "dark" ? "sun" : "moon"); b.setAttribute("aria-label", t === "dark" ? "Switch to light theme" : "Switch to dark theme"); }
    window.dispatchEvent(new Event("studio-theme"));
  };
  S.store = {
    get(k, d) { try { const v = localStorage.getItem("fo." + k); return v == null ? d : JSON.parse(v); } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem("fo." + k, JSON.stringify(v)); } catch (e) { /* private mode */ } },
  };
  S.split = function (s) {             // shell-like split: quotes keep spaces
    const out = []; let cur = "", q = null, any = false;
    for (const c of String(s || "")) {
      if (q) { if (c === q) q = null; else cur += c; }
      else if (c === '"' || c === "'") { q = c; any = true; }
      else if (/\s/.test(c)) { if (cur || any) out.push(cur); cur = ""; any = false; }
      else cur += c;
    }
    if (cur || any) out.push(cur);
    return out;
  };
  S.quote = (a) => (/^[\w@%+=:,./-]+$/.test(a) ? a : !String(a).includes("'") ? "'" + a + "'" :
    '"' + String(a).replace(/"/g, "") + '"');

  S.api = async function (path, body, opts = {}) {
    const init = {headers: {"X-Requested-With": "studio"}};
    if (body !== undefined) {
      init.method = "POST";
      if (body instanceof Blob || body instanceof ArrayBuffer) init.body = body;
      else { init.body = JSON.stringify(body); init.headers["Content-Type"] = "application/json"; }
    }
    const r = await fetch(path, init);
    let j;
    try { j = await r.json(); } catch (e) { j = {ok: false, error: "bad response " + r.status}; }
    if (!j.ok && !opts.quiet) S.toast(j.error || "request failed", 6000, "bad");
    return j;
  };
  S.fileUrl = (p, dl) => "/api/file?path=" + encodeURIComponent(p) + (dl ? "&download=1" : "");

  // ------------------------------------------------------------------ pages / router
  S.registerPage = function (p) { S.pages[p.id] = p; };
  S.go = function (id, params) {
    const q = params ? "?" + new URLSearchParams(params).toString() : "";
    location.hash = "#/" + id + q;
  };
  S.route = function () {
    const m = location.hash.match(/^#\/([\w-]+)(?:\?(.*))?/);
    const ids = Object.values(S.pages).sort((a, b) => (a.order || 99) - (b.order || 99)).map((p) => p.id);
    const id = m && S.pages[m[1]] ? m[1] : ids[0];
    const params = Object.fromEntries(new URLSearchParams((m && m[2]) || ""));
    const page = S.pages[id];
    const navId = page.navAs ? page.navAs(params) || id : id;
    const nav = document.getElementById("nav");
    nav.innerHTML = "";
    for (const pid of ids) {
      const p = S.pages[pid];
      if (p.nav === false) continue;
      const a = h("a", {href: p.href || "#/" + pid, class: pid === navId ? "on" : "", "aria-current": pid === navId ? "page" : null},
        p.icon ? S.icon(p.icon) : null, p.title);
      if (pid === "jobs") a.append(S._jobBadge = h("span", {class: "count", style: {display: "none"}}));
      nav.append(a);
    }
    S.updateJobBadge();
    if (S._stop) { S._stop.forEach((f) => f()); }
    S._stop = [];
    const main = document.getElementById("main");
    main.innerHTML = "";
    // a fresh container per visit: late async output of a page you already left lands in a detached node
    const box = document.createElement("div");
    main.append(box);
    page.render(box, params);
  };
  S.updateJobBadge = async function () {
    const r = await S.api("/api/jobs/running", undefined, {quiet: true});
    const b = S._jobBadge;
    if (!b || !r.ok) return;
    b.style.display = r.n ? "" : "none";
    b.textContent = r.n + " running";
  };
  S.onLeave = (f) => S._stop.push(f);          // pages register timers to stop when leaving
  // tiny event bus between the parts of a page (cleared when the page is left)
  S.on = (ev, f) => { const g = (e) => f(e.detail); window.addEventListener("studio:" + ev, g); S.onLeave(() => window.removeEventListener("studio:" + ev, g)); };
  S.emit = (ev, d) => window.dispatchEvent(new CustomEvent("studio:" + ev, {detail: d}));
  S.every = function (ms, f) { const t = setInterval(f, ms); S.onLeave(() => clearInterval(t)); return t; };

  S.boot = async function () {
    const info = await S.api("/api/info");
    if (!info.ok) { document.getElementById("main").textContent = info.error; return; }
    S.info = info;
    const hi = document.getElementById("hostinfo");
    hi.innerHTML = "";
    hi.append(h("span", {class: "dot" + (info.slurm ? "" : " off")}), info.host + " · " + (info.slurm ? "SLURM" : "no SLURM here") + " · v" + info.version);
    hi.title = info.slurm ? "SLURM commands found: submit, squeue and scancel work from here" :
      "squeue/sbatch not found on this machine: SLURM scripts can be written, not submitted";
    const tb = document.getElementById("themeBtn");
    if (tb) { tb.onclick = () => S.setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark"); S.setTheme(document.documentElement.dataset.theme || "light"); }
    setInterval(S.updateJobBadge, 15000);
    await S.loadModules();
    for (const p of info.pages) {
      await new Promise((res) => {
        const s = document.createElement("script");
        s.src = "/static/pages/" + p; s.onload = res;
        s.onerror = () => { S.toast("page " + p + " failed to load"); res(); };
        document.body.append(s);
      });
    }
    window.addEventListener("hashchange", S.route);
    S.route();
  };
  S.loadModules = async function () {
    const m = await S.api("/api/modules");
    S.modules = m.ok ? m.modules : [];
    S.moduleErrors = m.errors || [];
    S.moduleDisabled = m.disabled || {};
    S.moduleDirs = m.module_dirs || [];
    return S.modules;
  };
  S.mods = (kind) => S.modules.filter((m) => m.kind === kind);

  // ------------------------------------------------------------------ path picker (server-side files)
  S.pick = function (opts = {}) {
    // opts: {dir: bool, start, title, filter: regex string, upload: bool}
    return new Promise(async (resolve) => {
      let cur = opts.start || S.store.get("lastdir", S.info.start_dir);
      let sel = null;
      const list = h("div", {class: "bd"});
      const pathIn = h("input", {type: "text", class: "path grow", value: cur});
      const selLbl = h("span", {class: "muted small grow"});
      const hid = h("input", {type: "checkbox", onchange: () => load(cur)});
      const newName = h("input", {type: "text", placeholder: "folder name", style: {width: "130px", display: "none"}});
      const close = (v) => { m.remove(); resolve(v); };
      const upload = h("input", {type: "file", multiple: true, style: {display: "none"}, onchange: async (e) => {
        for (const f of e.target.files) {
          const r = await S.api("/api/upload?dir=" + encodeURIComponent(cur) + "&name=" + encodeURIComponent(f.name), f);
          if (r.ok) S.toast("uploaded " + f.name);
        }
        load(cur);
      }});
      const m = h("div", {class: "modal", onclick: (e) => { if (e.target === m) close(null); }},
        h("div", {class: "box"},
          h("div", {class: "hd"}, h("b", {}, opts.title || (opts.dir ? "Choose a folder" : "Choose a file"))),
          h("div", {class: "hd"},
            h("button", {onclick: () => load(S.dirname(cur))}, "↑"),
            h("button", {onclick: () => load(S.info.home)}, "~"),
            pathIn, h("button", {onclick: () => load(pathIn.value)}, "Go")),
          list,
          h("div", {class: "ft"}, selLbl,
            h("label", {class: "small"}, hid, " hidden"),
            h("button", {onclick: async () => {
              const n = (newName.value || "").trim();
              if (!n) { newName.style.display = ""; newName.focus(); return S.toast("type the new folder's name, then press New folder again"); }
              const r = await S.api("/api/mkdir", {path: cur + "/" + n});
              if (r.ok) { newName.value = ""; newName.style.display = "none"; load(cur); }
            }}, "New folder"), newName,
            opts.upload !== false ? h("button", {onclick: () => upload.click(), title: "copy files from your computer into this folder"}, "Upload here…") : null,
            upload,
            h("button", {onclick: () => close(null)}, "Cancel"),
            h("button", {class: "primary", onclick: () => {
              if (opts.dir) { S.store.set("lastdir", cur); close(sel && sel.dir ? sel.path : cur); }
              else if (sel && !sel.dir) { S.store.set("lastdir", cur); close(sel.path); }
              else S.toast("select a file");
            }}, opts.dir ? "Use this folder" : "Use this file"))));
      document.body.append(m);
      pathIn.addEventListener("keydown", (e) => { if (e.key === "Enter") load(pathIn.value); });
      const re = opts.filter ? new RegExp(opts.filter) : null;
      async function load(p) {
        const r = await S.api("/api/ls?path=" + encodeURIComponent(p) + (hid.checked ? "&hidden=1" : ""));
        if (!r.ok) return;
        cur = r.path; pathIn.value = cur; sel = null;
        selLbl.textContent = opts.dir ? cur : "";
        list.innerHTML = "";
        for (const e of r.entries) {
          if (!e.dir && opts.dir) continue;
          if (!e.dir && re && !re.test(e.name)) continue;
          const full = (cur === "/" ? "" : cur) + "/" + e.name;
          const row = h("div", {class: "fl"}, e.dir ? "📁" : "📄", e.name, h("span", {class: "sz"}, e.dir ? "" : S.size(e.size)));
          row.onclick = () => { list.querySelectorAll(".sel").forEach((x) => x.classList.remove("sel")); row.classList.add("sel"); sel = {path: full, dir: e.dir}; selLbl.textContent = full; };
          row.ondblclick = () => { if (e.dir) load(full); else { S.store.set("lastdir", cur); close(full); } };
          list.append(row);
        }
        if (!list.children.length) list.append(h("div", {class: "muted small", style: {padding: "8px 14px"}}, "(empty)"));
      }
      load(cur);
    });
  };
  // text field + "…" button
  S.pathField = function (value, opts = {}) {
    const inp = h("input", {type: "text", class: "path grow", value: value || "", placeholder: opts.placeholder || ""});
    const b = h("button", {onclick: async () => {
      const p = await S.pick({...opts, start: inp.value ? (opts.dir ? inp.value : S.dirname(inp.value)) : undefined});
      if (p) { inp.value = p; inp.dispatchEvent(new Event("change")); }
    }}, "…");
    const wrap = h("div", {class: "row grow"}, inp, b);
    wrap.input = inp;
    return wrap;
  };

  // ------------------------------------------------------------------ option form from an argparse description
  S.argForm = function (spec, opts = {}) {
    // spec: from /api/argspec; opts.skipFirstPositional, opts.saved (array of args to preset)
    const box = h("div", {class: "form"});
    const fields = [];
    const saved = opts.saved || [];
    let firstPosSkipped = !opts.skipFirstPositional;
    let lastGroup = null;
    const acts = (spec.actions || []).filter((a) => a.kind !== "subcommand");
    for (const a of acts) {
      if (a.positional && !firstPosSkipped) { firstPosSkipped = true; continue; }
      if (opts.hide && opts.hide.includes(a.dest)) continue;
      const g = a.group && !/^(options|optional arguments|positional arguments)$/.test(a.group) ? a.group : null;
      if (g && g !== lastGroup) { box.append(h("div", {class: "grp"}, g)); lastGroup = g; }
      const flag = a.positional ? null : a.flags.find((f) => f.startsWith("--")) || a.flags[0];
      const pre = presetFor(a, flag, saved);
      let input, get;
      if (a.kind === "bool") {                     // BooleanOptionalAction: --x / --no-x
        const neg = a.flags.find((f) => f.startsWith("--no-"));
        input = h("select", {}, h("option", {value: ""}, `(default: ${a.default})`), h("option", {value: "y"}, flag), h("option", {value: "n"}, neg));
        input.value = saved.includes(neg) ? "n" : saved.includes(flag) ? "y" : "";
        get = () => (input.value === "y" ? [flag] : input.value === "n" ? [neg] : []);
      } else if (a.kind === "true" || a.kind === "false" || a.kind === "const") {
        input = h("input", {type: "checkbox"});
        input.checked = pre === true;
        get = () => (input.checked ? [flag] : []);
      } else if (a.kind === "count") {
        input = h("input", {type: "number", min: 0, value: pre || 0, style: {width: "80px"}});
        get = () => Array(Number(input.value) || 0).fill(flag);
      } else if (a.choices && (a.nargs == null || a.nargs === "?")) {
        input = h("select", {}, h("option", {value: ""}, a.default != null ? `(default: ${a.default})` : "(not set)"),
          a.choices.map((c) => h("option", {value: c}, c)));
        input.value = pre != null ? pre : "";
        get = () => (input.value === "" ? [] : flag ? [flag, input.value] : [input.value]);
      } else {
        const multi = (a.nargs != null && a.nargs !== "?") || a.kind === "extend";
        const ph = a.default != null && a.default !== "" && !(Array.isArray(a.default) && !a.default.length) ?
          "default: " + (Array.isArray(a.default) ? a.default.join(" ") : a.default) :
          a.choices ? a.choices.join(" | ") : a.metavar ? [].concat(a.metavar).join(" ") : "";
        const isPath = /file|path|dir|poscar|template|out|nep|csv|refs|folder|config|from/i.test(a.dest + " " + (a.metavar || ""));
        const pf = isPath ? S.pathField(pre || "", {dir: /dir|folder|refs/i.test(a.dest), placeholder: ph}) : null;
        input = pf || h("input", {type: "text", class: "grow", value: pre || "", placeholder: ph});
        const ti = pf ? pf.input : input;
        if (a.kind === "append") ti.title = "several values: separate with ;";
        get = () => {
          const v = ti.value.trim();
          if (!v) return [];
          if (a.kind === "append") return v.split(";").map((x) => x.trim()).filter(Boolean).flatMap((x) => [flag, ...(multi ? S.split(x) : [x])]);
          const vals = multi ? S.split(v) : [v];
          return flag ? [flag, ...vals] : vals;
        };
      }
      const req = a.required && !a.positional ? " *" : a.positional ? " *" : "";
      box.append(h("label", {class: "k", title: a.flags.join(", ")}, (flag || a.dest) + req),
        h("div", {}, input, a.help ? h("div", {class: "help"}, a.help + (a.kind === "append" ? " (several: separate with ;)" : "")) : null));
      const numRx = a.type === "int" ? /^[+-]?\d+$/ : a.type === "float" ? /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$|^[+-]?(inf|nan)$/i : null;
      const check = () => {                       // a number option with a non-number: say so before running
        if (!numRx) return null;
        const vals = get().filter((x) => x !== flag);
        const bad = vals.find((x) => !numRx.test(x));
        return bad != null ? `${flag || a.dest}: "${bad}" is not ${a.type === "int" ? "a whole number" : "a number"}` : null;
      };
      fields.push({a, flag, get, check});
    }
    const extra = h("input", {type: "text", class: "grow mono", placeholder: "more arguments, typed as on the command line",
      value: opts.savedExtra || ""});
    box.append(h("label", {class: "k"}, "extra"), h("div", {}, extra));
    box.getArgs = (noExtra) => {
      const pos = [], optv = [];
      for (const f of fields) (f.a.positional ? pos : optv).push(...f.get());
      return [...pos, ...optv, ...(noExtra ? [] : S.split(extra.value))];
    };
    box.getExtra = () => extra.value;
    box.check = () => fields.map((f) => f.check()).filter(Boolean);
    return box;
  };
  function presetFor(a, flag, saved) {
    if (!saved || !saved.length) return a.kind === "true" ? false : null;
    if (a.positional) return null;
    const i = saved.indexOf(flag);
    if (a.kind === "true" || a.kind === "false" || a.kind === "const") return i >= 0;
    if (i < 0) return null;
    if (a.kind === "count") return saved.filter((x) => x === flag).length;
    const vals = [];
    for (let k = i + 1; k < saved.length && !/^--?[a-zA-Z]/.test(saved[k]); k++) {
      vals.push(saved[k]);
      if (a.nargs == null || a.nargs === "?") break;
    }
    if (a.kind === "append") {
      const all = [];
      saved.forEach((x, k) => { if (x === flag && k + 1 < saved.length) all.push(saved[k + 1]); });
      return all.join("; ");
    }
    return vals.map(S.quote).join(" ");
  }

  // ------------------------------------------------------------------ job viewer
  // what to do next when a job failed with a known message (pages and modules may push more)
  // runs found below the start folder, as a small clickable table (used by Results and Pipeline)
  S.runsList = function (onPick, opts = {}) {
    const box = h("div", {class: "small muted"}, "looking for runs below the start folder…");
    S.api("/api/runs/list", undefined, {quiet: true}).then((r) => {
      box.innerHTML = "";
      if (!r.ok || !r.runs.length) { box.append(opts.empty || `No run folders (campaign.toml) below ${r.base || "the start folder"} yet.`); return; }
      box.className = "";
      const t = h("table", {}, h("thead", {}, h("tr", {}, ["run", "starting phase", "created", "folder"].map((x) => h("th", {class: "l"}, x)))),
        h("tbody", {}, r.runs.slice(0, 30).map((x) => h("tr", {},
          h("td", {class: "l"}, h("a", {href: "javascript:void 0", onclick: () => onPick(x)}, x.name || S.basename(x.folder))),
          h("td", {class: "l"}, x.starting_phase || ""), h("td", {class: "l small"}, x.created || ""),
          h("td", {class: "l mono small"}, x.folder)))));
      box.append(h("h3", {}, `Runs below ${r.base}`), h("div", {class: "tablewrap", style: {maxHeight: "260px"}}, t));
    });
    return box;
  };
  S.hints = [
    [/compositions not prepared yet/, "Run stage 2 (prepare) first, or tick “plan only the compositions prepared so far” in stage 3."],
    [/no plan yet/, "Run stage 3 (plan) first."],
    [/exists\. --append adds structures/, "Already planned: tick “add structures not yet planned” (append) or “re-plan” (force) in stage 3."],
    [/mapping failed for \d+ of/, "Check the target templates: 12-atom cells of each phase with axes along the fluorite cube axes (see the log named above)."],
    [/NOT identical to their parent ordering/, "A mapped structure changed the cation ordering: check the templates and the map_args of the run."],
    [/was prepared with templates/, "This workdir was prepared with other templates: use “copy settings to a new run” instead of changing them."],
    [/NEP model .* not found/, "The NEP model file is missing: fix [free_energy] nep (copy settings to a new run with the right file)."],
    [/No module named '?(calorine|ase|phonopy|numpy|scipy|matplotlib)/, "A Python package is missing: start the app from the environment that has it (or set [cluster] setup for SLURM jobs)."],
    [/sbatch: (command )?not found|No such file or directory: 'sbatch'/, "sbatch is not available on this machine: use the 'write script' mode and submit on the cluster."],
  ];
  S.jobView = function (jid, opts = {}) {
    const st = h("span", {class: "pill run"}, "starting");
    const meta = h("span", {class: "muted small"});
    const log = h("pre", {class: "log"}, "");
    const outs = h("div", {});
    const cancel = h("button", {class: "danger", onclick: async () => {
      const r = await S.api("/api/jobs/cancel", {id: jid});
      S.toast(r.ok ? "stop signal sent" : "could not stop (already finished?)");
    }}, "Stop");
    const hint = h("div", {});
    const trunc = h("div", {});
    const box = h("div", {class: "card jobcard"},
      h("div", {class: "row"}, h("b", {}, opts.title || "Job"), st, meta, h("span", {class: "grow"}), cancel,
        opts.inJobs ? null : h("a", {href: "#/jobs?id=" + jid}, "open in Jobs")),
      hint, trunc, log, outs);
    let done = false, timer = null;
    async function poll() {
      const r = await S.api("/api/jobs?id=" + jid + "&tail=60000", undefined, {quiet: true});
      if (!r.ok) return;
      const j = r.job;
      st.className = "pill " + ({done: "ok", failed: "bad", lost: "warn", cancelled: "warn"}[j.state] || "run");
      st.textContent = j.state + (j.exit != null && j.exit !== 0 ? " (exit " + j.exit + ")" : "");
      meta.textContent = (j.end ? "took " + S.dur(j.end - j.start) : "running " + S.dur(Date.now() / 1000 - j.start)) +
        (j.slurm_ids && j.slurm_ids.length ? " · SLURM " + j.slurm_ids.join(", ") : "");
      const atBottom = log.scrollTop + log.clientHeight >= log.scrollHeight - 30;
      const cut = j.log_size && j.log && j.log_size > j.log.length + 10;
      log.textContent = (cut ? j.log.slice(j.log.indexOf("\n") + 1) : j.log) || "";
      trunc.innerHTML = "";
      if (cut) trunc.append(h("div", {class: "small muted"}, `log is ${(j.log_size / 1e6).toFixed(1)} MB: showing its end · `,
        h("a", {href: S.fileUrl(j.log_path), target: "_blank"}, "full log"), " · ", h("a", {href: S.fileUrl(j.log_path, true)}, "download")));
      if (atBottom) log.scrollTop = log.scrollHeight;
      cancel.style.display = j.state === "running" ? "" : "none";
      if (["done", "failed", "lost", "cancelled"].includes(j.state)) {
        done = true; clearInterval(timer);
        outs.innerHTML = "";
        outs.append(S.outputsView(j.outputs || []));
        hint.innerHTML = "";
        const tail = (j.log || "").slice(-20000);
        if (j.state === "failed") {
          const hit = S.hints.find(([rx]) => rx.test(tail));
          if (hit) hint.append(h("div", {class: "notice warn"}, "Next step: " + hit[1]));
        } else if (j.state === "done") {        // exit 0 but some structures failed (batch_free_energy / merge)
          let nf = 0;
          for (const m of tail.matchAll(/(\d+) ok, (\d+) failed/g)) nf += Number(m[2]);
          if (nf > 0) { st.className = "pill warn"; st.textContent = `done, ${nf} failed`;
            hint.append(h("div", {class: "notice warn"}, `${nf} structure(s) FAILED (see the log). Check the NEP model covers every element; run Status, then Rerun flagged or compute again (failed ones are retried).`)); }
        }
        if (opts.onDone) opts.onDone(j);
      }
    }
    timer = setInterval(() => { if (!done) poll(); }, 1500);
    S.onLeave(() => clearInterval(timer));
    box.stop = () => clearInterval(timer);
    poll();
    return box;
  };
  S.outputsView = function (outs) {
    if (!outs.length) return h("div", {class: "muted small"}, "no new files in the watched folders");
    const box = h("div", {});
    box.append(h("h3", {}, "Files written"));
    const prev = h("div", {});
    const lst = h("div", {class: "row"});
    let firstShown = false;
    const order = (o) => [".html", ".png", ".svg"].indexOf(o.name.slice(o.name.lastIndexOf(".")).toLowerCase());
    const sorted = [...outs].sort((a, b) => ((order(a) + 1 || 9) - (order(b) + 1 || 9)));
    for (const o of sorted) {
      const ext = o.name.slice(o.name.lastIndexOf(".")).toLowerCase();
      const show = () => { prev.innerHTML = ""; prev.append(S.preview(o.path)); };
      lst.append(h("span", {}, h("a", {href: "javascript:void 0", onclick: show}, o.name), " ",
        h("a", {href: S.fileUrl(o.path, true), class: "small muted"}, "⤓")));
      if (!firstShown && [".html", ".png", ".svg"].includes(ext)) { firstShown = true; show(); }
    }
    box.append(lst, prev);
    return box;
  };
  S.preview = function (path) {
    const ext = path.slice(path.lastIndexOf(".")).toLowerCase();
    const url = S.fileUrl(path);
    const head = h("div", {class: "row small", style: {margin: "8px 0"}}, h("span", {class: "mono grow"}, path),
      h("a", {href: url, target: "_blank"}, "open in new tab"), h("a", {href: S.fileUrl(path, true)}, "download"));
    if (ext === ".html") return h("div", {}, head, h("iframe", {class: "view", src: url, sandbox: "allow-scripts allow-downloads allow-popups"}));
    if ([".png", ".svg", ".jpg", ".jpeg", ".gif"].includes(ext)) return h("div", {}, head, h("img", {class: "view", src: url}));
    if (ext === ".pdf") return h("div", {}, head, h("iframe", {class: "view", src: url}));
    const pre = h("pre", {class: "log"}, "…");
    S.api("/api/text?path=" + encodeURIComponent(path) + "&tail=300000").then((r) => { pre.textContent = r.ok ? r.text : r.error; });
    return h("div", {}, head, pre);
  };
})();
