// test_single_ui.js - the Single runs page in a real browser (Playwright + Chromium), mock calorine (no NEP).
//   python3 tests/test_single_runs.py && node tests/test_single_ui.js        (~1 min)
// Upload POSCARs and nep.txt from "your computer", add one from the server by path, tasks and settings,
// show the command, run locally, results (summary table, figures, flags), re-plot with another reference,
// write a SLURM script, copy settings back, dark mode and phone width, no page errors.
// Screenshots: tests/work_single/ui_single_*.png
const {chromium} = require("playwright");
const {spawn} = require("child_process");
const fs = require("fs"), path = require("path");
const ROOT = path.resolve(__dirname, ".."), W = path.join(ROOT, "tests", "work_single");
const PORT = 9050 + Math.floor(Math.random() * 40);
let fails = 0;
const check = (c, m) => { console.log((c ? "   ok   " : "   FAIL ") + m); if (!c) fails++; };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

(async () => {
  const inp = path.join(W, "inp");
  if (!fs.existsSync(path.join(inp, "POSCAR_m_1x1x1"))) { console.log("run python3 tests/test_single_runs.py first"); process.exit(2); }
  const state = path.join(W, "ui_state"), base = path.join(W, "ui_runs");
  fs.rmSync(state, {recursive: true, force: true});
  fs.rmSync(base, {recursive: true, force: true});
  const env = {...process.env, PYTHONPATH: path.join(ROOT, "tests", "mock_calorine"), PYTHONWARNINGS: "ignore", OMP_NUM_THREADS: "1"};
  const srv = spawn("python3", [path.join(ROOT, "start_studio.py"), "--port", String(PORT), "--state", state, "--dir", W], {stdio: "ignore", env});
  try {
    for (let i = 0; i < 50 && !fs.existsSync(path.join(state, "token")); i++) await sleep(200);
    await sleep(800);
    const tok = fs.readFileSync(path.join(state, "token"), "utf8").trim();
    const b = await chromium.launch();
    const ctx = await b.newContext({viewport: {width: 1400, height: 1000}});
    const p = await ctx.newPage();
    const errs = [];
    p.on("pageerror", (e) => errs.push(e.message));
    p.on("console", (m) => { if (m.type() === "error" && !/Failed to load resource/.test(m.text())) errs.push(m.text()); });
    p.on("response", (r) => { if (r.status() >= 400) errs.push(r.status() + " " + r.url().replace(/token=[^&]+/, "token=…")); });
    const url = `http://127.0.0.1:${PORT}`;
    await p.goto(`${url}/?token=${tok}`);
    await p.waitForSelector("nav a");
    check(await p.locator("nav a", {hasText: "Single runs"}).count() === 1, "Single runs in the navigation");
    // runs folder of this test, a short ramp and a small phonon mesh (fast)
    await p.evaluate((o) => localStorage.setItem("fo.single.state", JSON.stringify(o)),
      {base, Tmode: "ramp", Tmin: 0, Tmax: 1000, Tstep: 100, Textra: "300", min_length: 6, mesh: 8, formats: ["png", "svg"],
       xrd_ref: ["c", "m"], at_T: "300 900"});
    await p.goto(`${url}/#/single`);
    await p.waitForSelector("text=New single run");

    console.log("== structures and model");
    await p.locator("input[type=file][multiple]").setInputFiles([path.join(inp, "POSCAR_t_1x2x4"), path.join(inp, "POSCAR_oIII_1x1x1")]);
    for (let i = 0; i < 50 && await p.locator(".mono.small.grow", {hasText: "_uploads"}).count() < 2; i++) await sleep(200);
    check(await p.locator(".mono.small.grow", {hasText: "_uploads"}).count() === 2, "two POSCARs uploaded into <runs>/_uploads/");
    await p.locator("button", {hasText: "Add POSCAR…"}).click();
    await p.locator(".modal input.path").fill(inp);
    await p.locator(".modal button", {hasText: "Go"}).click();
    await p.locator(".modal .fl", {hasText: "POSCAR_m_3x2x2"}).dblclick();
    check(await p.locator(".mono.small.grow", {hasText: "POSCAR_m_3x2x2"}).count() === 1, "third POSCAR picked on the server");
    await p.locator("input[placeholder='POSCAR_m_3x2x2']").fill("mono");
    await p.locator("input[placeholder='POSCAR_m_3x2x2']").dispatchEvent("change");
    const nepUp = p.locator("input[type=file]:not([multiple])").first();
    await nepUp.setInputFiles(path.join(W, "nep.txt"));
    await sleep(500);
    check(/_uploads\/.*nep\.txt$/.test(await p.locator("input[placeholder='nep.txt']").inputValue()), "nep.txt uploaded");
    check(await p.locator("text=references folder").count() === 1 && await p.locator("label", {hasText: /^ oIII$/}).count() >= 1, "XRD reference choices from refs/");
    // reference choice lists the labels
    const refOpts = await p.locator("select:has(option[value=lowest])").first().locator("option").allTextContents();
    check(refOpts.includes("mono") && refOpts.includes("POSCAR_t_1x2x4"), "ΔF reference choices = structure labels: " + refOpts.join(", "));
    await p.locator("select:has(option[value=lowest])").first().selectOption("mono");

    console.log("== command, run, results");
    await p.locator("button", {hasText: "Show command"}).click();
    const cmd = await p.locator("pre.log").first().textContent();
    check(/--T-range 0 1000 100/.test(cmd) && /--reference mono/.test(cmd) && /--xrd-ref c m/.test(cmd) && /--names POSCAR_t_1x2x4 POSCAR_oIII_1x1x1 mono/.test(cmd),
      "command shows ramp, reference, XRD refs, names");
    await p.screenshot({path: path.join(W, "ui_single_form.png"), fullPage: true});
    await p.locator("button.primary", {hasText: /^Run$/}).click();
    await p.waitForSelector("text=Summary", {timeout: 180000});
    // the summary table fills in after the heading (slow under load): wait for it, then count
    await p.waitForSelector("h3:has-text('Summary') + .tablewrap tbody tr", {timeout: 30000}).catch(() => {});
    const rows = await p.locator("h3:has-text('Summary') + .tablewrap tbody tr").count();
    check(rows === 3, `summary table: ${rows} rows`);
    await p.waitForSelector("iframe.view, img.view", {timeout: 10000});
    const figs = await p.locator("h3:has-text('Figures') + .row a").allTextContents();
    check(figs.includes("harmonic_compare.html") && figs.includes("F_vs_T.png") && figs.includes("xrd_compare.png"), "figures listed: " + figs.join(" "));
    const fr = p.frameLocator("iframe.view").first();
    await fr.locator("select#ref").waitFor({timeout: 10000});
    check((await fr.locator("select#ref").inputValue()) === "mono", "interactive HTML opens with the chosen reference");
    await fr.locator("select#q").selectOption("S");
    check(/S vib/.test(await fr.locator("svg").first().textContent()), "HTML: switch quantity to S");
    await fr.locator("svg rect[fill=transparent]").hover({position: {x: 300, y: 150}});
    await sleep(200);
    check(/K/.test(await fr.locator("#tip").textContent()), "HTML: hover tooltip");
    await p.screenshot({path: path.join(W, "ui_single_results.png"), fullPage: true});

    console.log("== re-plot, SLURM script, copy settings");
    await p.locator("h3:has-text('Re-plot') + .row select").selectOption("lowest");
    await p.locator("button.primary", {hasText: "Re-plot"}).click();
    const runs = fs.readdirSync(base).filter((x) => !x.startsWith("_"));
    const cmpPath = path.join(base, runs[0], "results", "plots", "comparison.csv");
    let cmp = "";
    for (let i = 0; i < 120 && !/lowest F at each T/.test(cmp); i++) { await sleep(500); cmp = fs.readFileSync(cmpPath, "utf8"); }
    await p.waitForSelector("text=Summary", {timeout: 20000});
    check(/lowest F at each T/.test(cmp), "re-plot wrote comparison.csv vs lowest");
    await p.locator("button", {hasText: "Copy settings to the form"}).click();
    await sleep(300);
    check(await p.locator("text=No structures yet.").count() === 1, "copy settings: form filled, structures cleared");
    await p.locator("select:has(option[value=slurm_write])").selectOption("slurm_write");
    await p.locator("input[type=file][multiple]").setInputFiles([path.join(inp, "POSCAR_c_1x1x1")]);
    await sleep(500);
    await p.locator("button.primary", {hasText: "Write script"}).click();
    await p.waitForSelector("text=SLURM script written", {timeout: 10000});
    check(await p.locator("text=SLURM script written").count() === 1, "SLURM script written from the page");
    await p.waitForSelector("td >> text=written");
    check(await p.locator("td .pill", {hasText: "written"}).count() === 1, "runs list: the written run");

    console.log("== dark mode, phone width");
    await p.emulateMedia({colorScheme: "dark"});
    await p.screenshot({path: path.join(W, "ui_single_dark.png"), fullPage: true});
    await p.setViewportSize({width: 390, height: 900});
    const over = await p.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
    check(over < 40, `phone width: horizontal overflow ${over}px`);
    await p.screenshot({path: path.join(W, "ui_single_phone.png"), fullPage: false});
    check(errs.length === 0, "no page errors" + (errs.length ? ": " + errs.slice(0, 5).join(" | ") : ""));
    await b.close();
  } finally {
    srv.kill();
  }
  console.log(fails ? `${fails} FAILED` : "SINGLE RUNS UI TESTS PASSED");
  process.exit(fails ? 1 : 0);
})();
