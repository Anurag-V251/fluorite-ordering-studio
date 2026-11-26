// test_ui.js - click-through test of the page in a real browser (Playwright + Chromium), no NEP needed.
//   node tests/test_ui.js            (after python3 tests/test_studio.py, which makes tests/work_studio/)
// Starts its own server with a RELATIVE --state path, creates a run through the New-run form, runs dry-run,
// prepare, plan, compute, merge by clicking, scans it on Results, runs a plot and map_check, checks dark mode,
// and fails on any page error.  Screenshots: tests/work_studio/ui_*.png
const {chromium} = require("playwright");
const {spawn} = require("child_process");
const fs = require("fs"), path = require("path");
const ROOT = path.resolve(__dirname, ".."), W = path.join(ROOT, "tests", "work_studio");
const PORT = 8899 + Math.floor(Math.random() * 50);
let fails = 0;
const check = (c, m) => { console.log((c ? "   ok   " : "   FAIL ") + m); if (!c) fails++; };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

(async () => {
  if (!fs.existsSync(path.join(W, "POSCAR_c_2x2x2"))) { console.log("run python3 tests/test_studio.py first"); process.exit(2); }
  fs.rmSync(path.join(W, "ui_state"), {recursive: true, force: true});
  fs.rmSync(path.join(W, "runs_ui"), {recursive: true, force: true});
  const srv = spawn("python3", [path.join(ROOT, "start_studio.py"), "--port", String(PORT), "--state", "ui_state", "--dir", W],
    {cwd: W, env: {...process.env, PYTHONPATH: path.join(ROOT, "tests", "mock_calorine")}, stdio: "ignore"});
  try {
    for (let i = 0; i < 50 && !fs.existsSync(path.join(W, "ui_state", "token")); i++) await sleep(200);
    await sleep(800);
    const tok = fs.readFileSync(path.join(W, "ui_state", "token"), "utf8").trim();
    const b = await chromium.launch();
    const p = await b.newPage({viewport: {width: 1300, height: 950}});
    const errs = [];
    const expect400 = ["/api/runs/create"];     // the deliberately wrong composition
    p.on("pageerror", (e) => errs.push(e.message));
    p.on("console", (m) => { if (m.type() === "error" && !/Failed to load resource/.test(m.text())) errs.push(m.text()); });
    p.on("response", (r) => { if (r.status() >= 400 && !expect400.some((x) => r.url().includes(x))) errs.push(r.status() + " " + r.url().replace(/token=[^&]+/, "token=…")); });
    const base = `http://localhost:${PORT}`;
    await p.goto(`${base}/?token=${tok}`);
    await p.waitForSelector("nav a");

    console.log("== new run through the form");
    await p.goto(`${base}/#/pipeline?new=1`);
    await p.waitForSelector("text=Create run");
    const form = p.locator(".form");
    const fieldAfter = (label) => form.locator(`label.k:text-is("${label}") + div input`).first();
    await fieldAfter("folder").fill(path.join(W, "runs_ui", "run1"));
    await fieldAfter("name").fill("ui test run");
    await fieldAfter("host POSCAR").fill(path.join(W, "POSCAR_c_2x2x2"));
    await fieldAfter("host POSCAR").dispatchEvent("change");
    await fieldAfter("compositions").fill("0 8 16");
    await fieldAfter("max orderings").fill("2");
    const tl = form.locator(".tplrow input");
    await tl.nth(0).fill("t"); await tl.nth(1).fill(path.join(ROOT, "refs", "POSCAR_t"));
    await p.click("text=+ phase");
    await tl.nth(2).fill("m"); await tl.nth(3).fill(path.join(ROOT, "refs", "POSCAR_m"));
    await fieldAfter("NEP model").fill(path.join(W, "nep.txt"));
    await fieldAfter("phonon supercell").fill("1 1 1");
    await fieldAfter("mesh").fill("4");
    await fieldAfter("cores_per_node").fill("4");
    await fieldAfter("walltime").fill("00:10:00");
    await fieldAfter("sec_per_structure").fill("2");
    // a wrong composition first: the error must be shown next to the button
    await fieldAfter("compositions").fill("0, 40");
    await p.click("text=Create run");
    await p.waitForSelector(".notice.bad", {timeout: 10000});
    check((await p.textContent(".notice.bad")).includes("outside 0..32"), "bad composition refused with a clear message");
    await fieldAfter("compositions").fill("0 8 16");
    await p.click("text=Create run");
    await p.waitForURL(/pipeline\?config=/, {timeout: 20000});
    check(fs.existsSync(path.join(W, "runs_ui", "run1", "campaign.toml")), "Create run made campaign.toml");
    await p.screenshot({path: path.join(W, "ui_pipeline.png"), fullPage: true});

    console.log("== stages by clicking");
    for (const [num, btn] of [["1.", "Run"], ["2.", "Run here"], ["3.", "Run"], ["4.", "Run all unfinished chunks here"], ["6.", "Run here"]]) {
      const card = p.locator(".card", {has: p.locator(`b:text-matches("^${num.replace(".", "\\.")} ")`)}).first();
      await card.locator(`button:text-is("${btn}")`).click();
      let st = "";
      for (let i = 0; i < 240; i++) {
        const t = (await p.locator(".jobcard .pill").first().textContent().catch(() => "")) || "";
        if (/done|failed|lost|cancelled/.test(t)) { st = t; break; }
        await sleep(1000);
      }
      check(st.startsWith("done"), `stage ${num} ${btn}: ${st}`);
      await p.click("text=Clear finished cards");
    }

    console.log("== run dashboard: state strip, next step, theme");
    await p.click('button:has-text("Refresh")');
    await p.waitForFunction(() => { const e = document.querySelector(".st-kpi.total .kpi-v"); return e && e.textContent.trim() === "15"; }, null, {timeout: 20000});
    const kp2 = await p.locator(".st-kpi").allTextContents();
    check(kp2.some((t) => /Clean\s*15/.test(t)) && kp2.some((t) => /Chunks done\s*\d+\s*\/\s*\d+/.test(t)), "state strip: 15 structures, 15 clean, chunks " + kp2.join(" | ").replace(/\s+/g, " ").slice(0, 160));
    check((await p.locator(".next .what").textContent()).length > 3, "next step shown: " + (await p.locator(".next .what").textContent()));
    check(await p.locator(".nav a.on", {hasText: "Pipeline"}).count() === 1, "Pipeline highlighted in the nav");
    const th0 = await p.evaluate(() => document.documentElement.dataset.theme);
    await p.click("#themeBtn");
    const th1 = await p.evaluate(() => document.documentElement.dataset.theme);
    check(th0 !== th1 && ["light", "dark"].includes(th1), `theme button switches ${th0} -> ${th1}`);
    await p.screenshot({path: path.join(W, "ui_dashboard_" + th1 + ".png"), fullPage: false});
    await p.click("#themeBtn");
    await p.reload(); await p.waitForSelector(".strip .kpi-v");
    check(await p.evaluate(() => document.documentElement.dataset.theme) === th0, "theme choice kept after reload");

    console.log("== results of the new run");
    await p.goto(`${base}/#/results?path=${encodeURIComponent(path.join(W, "runs_ui", "run1"))}`);
    await p.waitForSelector("text=Progress per composition", {timeout: 30000});
    const kp = await p.locator(".kpi").allTextContents();
    check(kp.some((t) => /^\s*15\s*free energies done/.test(t)), "results: 15 free energies done (1+2+2 orderings x 3 phases) " + kp.slice(2, 5).join(" | "));
    check(await p.locator("h2", {hasText: "ui test run"}).count() === 1, "results header shows the run name");
    await p.screenshot({path: path.join(W, "ui_results.png"), fullPage: true});

    console.log("== a plot and map_check");
    await p.click('a:text-is("plot")');
    await p.waitForSelector("text=Run");
    await p.locator(".list .item", {hasText: "P(c/a) and P(angle)"}).click().catch(() => {});
    await p.locator(".list .item", {hasText: "F vs dopant %"}).click();
    await p.locator('button:text-is("Run")').click();
    await p.waitForSelector("img.view", {timeout: 120000});
    check(true, "F vs dopant % plot shown");
    await p.goto(`${base}/#/tools?module=map_check`);
    await p.waitForSelector('button:text-is("Run")');
    await p.locator("input.path").first().fill(path.join(W, "P_m_1x2x4"));
    await p.locator('button:text-is("Run")').click();
    await p.waitForSelector("text=Files written", {timeout: 120000});
    const files = await p.locator("main a").allTextContents();
    check(files.some((t) => t.includes("map_check.csv")), "map_check outputs listed (default folder in the app state)");
    await p.screenshot({path: path.join(W, "ui_tools.png"), fullPage: false});

    console.log("== dark mode links");
    await p.emulateMedia({colorScheme: "dark"});
    await p.goto(`${base}/#/results?path=${encodeURIComponent(path.join(W, "runs_ui", "run1"))}`);
    await p.waitForSelector("text=Progress per composition");
    const col = await p.locator("a", {hasText: "open run"}).evaluate((e) => getComputedStyle(e).color);
    check(col !== "rgb(0, 0, 238)" && col !== "rgb(85, 26, 139)", "dark-mode link colour " + col);
    await p.screenshot({path: path.join(W, "ui_dark.png"), fullPage: false});
    check(errs.length === 0, "no page errors" + (errs.length ? ": " + errs.slice(0, 3).join(" | ") : ""));
    await b.close();
  } finally {
    srv.kill();
  }
  console.log(fails ? `UI TEST FAILED (${fails})` : "UI TEST PASSED");
  process.exit(fails ? 1 : 0);
})();
