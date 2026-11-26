// test_xrd_ui.js - the XRD page in a real browser (Playwright + Chromium), on the mock campaign of test_xrd.py.
//   python3 tests/test_xrd.py && node tests/test_xrd_ui.js        (~1 min)
// Window by typing and by dragging, phase chips, Gaussian pick, live broadening (no new job), the page's
// average = xrd_average.py's average (CSV of the page vs files of the export job), species change -> new job,
// hover readout, zoom, dark mode, no page errors.  Screenshots: tests/work_xrd/ui_xrd_*.png
const {chromium} = require("playwright");
const {spawn} = require("child_process");
const fs = require("fs"), path = require("path");
const ROOT = path.resolve(__dirname, ".."), W = path.join(ROOT, "tests", "work_xrd");
const MOCK = process.env.XRD_MOCK || path.join(W, "mock_2x2x2");
const PORT = 8950 + Math.floor(Math.random() * 40);
let fails = 0;
const check = (c, m) => { console.log((c ? "   ok   " : "   FAIL ") + m); if (!c) fails++; };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const readCsv = (t) => { const L = t.trim().split("\n"); const hd = L[0].split(","); return L.slice(1).map((l) => Object.fromEntries(l.split(",").map((v, i) => [hd[i], v]))); };

(async () => {
  if (!fs.existsSync(path.join(MOCK, "campaign", "results"))) { console.log("run python3 tests/test_xrd.py first"); process.exit(2); }
  const state = path.join(W, "ui_state");
  fs.rmSync(state, {recursive: true, force: true});
  const srv = spawn("python3", [path.join(ROOT, "start_studio.py"), "--port", String(PORT), "--state", state, "--dir", W], {stdio: "ignore"});
  try {
    for (let i = 0; i < 50 && !fs.existsSync(path.join(state, "token")); i++) await sleep(200);
    await sleep(800);
    const tok = fs.readFileSync(path.join(state, "token"), "utf8").trim();
    const b = await chromium.launch();
    const ctx = await b.newContext({viewport: {width: 1400, height: 1000}, acceptDownloads: true});
    const p = await ctx.newPage();
    const errs = [];
    p.on("pageerror", (e) => errs.push(e.message));
    p.on("console", (m) => { if (m.type() === "error" && !/Failed to load resource/.test(m.text())) errs.push(m.text()); });
    p.on("response", (r) => { if (r.status() >= 400) errs.push(r.status() + " " + r.url().replace(/token=[^&]+/, "token=…")); });
    const base = `http://127.0.0.1:${PORT}`;
    await p.goto(`${base}/?token=${tok}`);
    await p.waitForSelector("nav a");
    const nJobs = async () => (await p.evaluate(async () => (await (await fetch("/api/jobs")).json()).jobs.length));

    console.log("== load, phases, window");
    await p.goto(`${base}/#/xrd?source=${encodeURIComponent(MOCK)}`);
    await p.waitForSelector(".xrd-chip");
    check((await p.locator(".xrd-chip").count()) === 4, "4 phase chips (parent t oIII m)");
    check(/drag on P\(F\)/.test(await p.locator("text=/drag on P/").first().textContent()), "no window yet: nothing computed, hint shown");
    for (const c of await p.locator(".xrd-chip").all()) if (!/oIII/.test(await c.textContent())) await c.click();
    const win = p.locator("input[type=number]");
    await win.nth(0).fill("1.5"); await win.nth(0).dispatchEvent("change");
    await win.nth(1).fill("4"); await win.nth(1).dispatchEvent("change");
    await p.waitForSelector("text=/\\d+ structures · CuKa/", {timeout: 60000});
    const st1 = await p.locator("text=/\\d+ structures · CuKa/").textContent();
    check(/^16 structures/.test(st1), `window 1.5..4 meV/atom, oIII only: "${st1}" (16 expected, as test_xrd.py)`);
    check(/16 structures selected \(oIII 16\)/.test(await p.locator("b", {hasText: "structures selected"}).locator("..").textContent()), "selection count shown before computing");

    console.log("== drag on P(F)");
    const pf = await p.locator(".xrd-plot svg").first().boundingBox();
    await p.mouse.move(pf.x + pf.width * 0.15, pf.y + pf.height * 0.5);
    await p.mouse.down(); await p.mouse.move(pf.x + pf.width * 0.3, pf.y + pf.height * 0.5, {steps: 6}); await p.mouse.up();
    const lo = Number(await win.nth(0).inputValue()), hi = Number(await win.nth(1).inputValue());
    check(hi > lo && !(lo === 1.5 && hi === 4), `drag sets the window (${lo}..${hi})`);
    const opt = await p.locator("select option", {hasText: "oIII · 1/"}).first().getAttribute("value");
    await p.locator(`select:has(option[value='${opt}'])`).selectOption(opt);
    check(Number(await win.nth(0).inputValue()) < Number(await win.nth(1).inputValue()), "pick a Gaussian: window = μ ± σ");
    await win.nth(0).fill("1.5"); await win.nth(0).dispatchEvent("change");
    await win.nth(1).fill("4"); await win.nth(1).dispatchEvent("change");
    await sleep(1500);
    await p.waitForSelector("text=/^16 structures · CuKa/", {timeout: 60000});

    console.log("== live broadening (no new job)");
    const meanD = async () => (await p.locator(".xrd-plot svg").nth(1).locator("path[stroke-width='1.3']").getAttribute("d"));
    const d0 = await meanD(), j0 = await nJobs();
    const fw = p.locator(".xrd-p:has(span.k:text-is('FWHM (°)')) input[type=number]");
    await fw.fill("0.3");
    await sleep(400);
    check((await meanD()) !== d0 && (await nJobs()) === j0, "FWHM 0.1 -> 0.3 redraws at once, no new job");
    await fw.fill("0.1"); await sleep(300);

    console.log("== the page's average = xrd_average.py's average");
    const [dl] = await Promise.all([p.waitForEvent("download"), p.click("button:text-is('CSV')")]);
    const pageCsv = readCsv(fs.readFileSync(await dl.path(), "utf8"));
    const od = path.join(W, "ui_export");
    fs.rmSync(od, {recursive: true, force: true});
    await p.click("summary:text-is('Command line and files')");
    const odIn = p.locator("input[placeholder='output folder']");
    await odIn.fill(od); await odIn.dispatchEvent("change");
    await p.click("text=Write files here (job)");
    await p.waitForSelector(".jobcard .pill.ok", {timeout: 120000});
    const f = fs.readdirSync(od).find((x) => /^xrd_.*meV\.csv$/.test(x));
    const pyCsv = readCsv(fs.readFileSync(path.join(od, f), "utf8"));
    let dmax = 0;
    for (let i = 0; i < pyCsv.length; i++) dmax = Math.max(dmax, Math.abs(Number(pyCsv[i].mean_rel_max100) - Number(pageCsv[i].mean_rel_max100)));
    check(pageCsv.length === pyCsv.length && dmax < 0.05, `same grid (${pyCsv.length} points), max |page - script| = ${dmax.toExponential(2)} (scale 0-100)`);
    const cli = await p.locator("pre.log").first().textContent();
    check(/xrd_average\.py .* --composition Ce_25 --phase oIII --F-range 1\.5 4 /.test(cli), "equivalent command line shown");
    await p.screenshot({path: path.join(W, "ui_xrd_light.png"), fullPage: true});

    console.log("== species change -> new sticks, references follow");
    const j1 = await nJobs();
    await p.locator(".xrd-p:has(span.k:text-is('species for the XRD')) select").selectOption("cations");
    await p.waitForSelector("text=/XRD species: cations→Zr/", {timeout: 60000});
    check((await nJobs()) > j1, "species change starts a new job");
    await p.click("summary:text-is('Structures used')");
    check((await p.locator("td:text-is('Zr32O64')").count()) === 16, "all 16 structures computed as Zr32O64 (files unchanged)");
    check((await p.locator("text=/\\(Zr4O8, species changed\\)/").count()) >= 1, "reference c shown as Zr4O8");

    console.log("== hover, zoom, display options");
    const xb = await p.locator(".xrd-plot svg").nth(1).boundingBox();
    await p.mouse.move(xb.x + xb.width * 0.30, xb.y + xb.height * 0.2);
    const ro = await p.locator(".xrd-readout").textContent();
    check(/2θ [\d.]+°  d [\d.]+ Å  mean/.test(ro), `hover readout: ${ro.slice(0, 70)}…`);
    await p.mouse.move(xb.x + xb.width * 0.22, xb.y + xb.height * 0.3);
    await p.mouse.down(); await p.mouse.move(xb.x + xb.width * 0.40, xb.y + xb.height * 0.3, {steps: 5}); await p.mouse.up();
    await sleep(200);
    const ticks = await p.locator(".xrd-plot svg").nth(1).locator("text").allTextContents();
    check(ticks.includes("26") && ticks.includes("36") && !ticks.includes("50"), "drag zooms the 2θ axis (about 24-40°)");
    await p.locator(".xrd-p:has(span.k:text-is('show')) select").selectOption("overlay");
    await p.locator(".xrd-p:has(span.k:text-is('spread')) select").selectOption("minmax");
    await p.locator(".xrd-p:has(span.k:text-is('y axis')) select").selectOption("sqrt");
    await sleep(300);
    await p.screenshot({path: path.join(W, "ui_xrd_zoom_overlay.png"), fullPage: true});
    await p.locator(".xrd-plot svg").nth(1).dblclick();
    await p.emulateMedia({colorScheme: "dark"});
    await p.locator(".xrd-p:has(span.k:text-is('show')) select").selectOption("stack");
    await sleep(400);
    await p.screenshot({path: path.join(W, "ui_xrd_dark.png"), fullPage: true});
    check(errs.length === 0, "no page errors" + (errs.length ? ": " + errs.slice(0, 5).join(" | ") : ""));
    await b.close();
  } finally {
    srv.kill();
  }
  console.log(fails ? `\n${fails} FAILED` : "\nALL OK");
  process.exit(fails ? 1 : 0);
})();
