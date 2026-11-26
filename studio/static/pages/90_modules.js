/* Modules page: what is registered, where it comes from, and how to add or replace a part. */
(function () {
  const S = Studio, h = S.h;
  S.registerPage({id: "modules", title: "Modules", order: 90, render});

  async function render(main) {
    await S.loadModules();
    const card = h("div", {class: "card"}, h("h2", {}, "Modules"),
      h("div", {class: "small"}, "Every pipeline stage, plot and tool is one small TOML file. To add one, or to replace a shipped one, put a file with the same ",
        h("span", {class: "mono"}, "id"), " in one of these folders (later folders win; ", h("span", {class: "mono"}, "enabled = false"), " hides a module):"),
      h("pre", {class: "log"}, S.info.module_dirs.join("\n")),
      h("div", {class: "small"}, "Pages are files in ", h("span", {class: "mono"}, "studio/static/pages/"), " (" + S.info.pages.join(", ") + "). Scripts are taken from ",
        h("span", {class: "mono"}, S.info.scripts), "; a run with a frozen scripts/ folder uses its own copy."),
      h("button", {style: {marginTop: "8px"}, onclick: () => S.loadModules().then(() => S.route())}, "Reload modules"));
    if (S.moduleErrors.length) card.append(h("div", {class: "notice bad"}, S.moduleErrors.join("\n")));
    const dis = Object.entries(S.moduleDisabled || {});
    if (dis.length) card.append(h("div", {class: "notice warn small"}, "Disabled (enabled = false): " + dis.map(([id, f]) => `${id} by ${f}`).join("; ")));
    if ((S.moduleDirs || []).length) card.append(h("div", {class: "small muted"}, "Module folders, later ones override earlier: ",
      S.moduleDirs.map((d) => h("span", {class: "mono", style: {marginRight: "10px"}}, d)),
      " (set FO_STUDIO_MODULES=dir1:dir2 to keep your own modules outside the app folder)."));
    main.append(card);
    for (const kind of ["stage", "plot", "tool"]) {
      const t = h("table", {});
      t.append(h("thead", {}, h("tr", {}, ["id", "title", "script / commands", "input", "file"].map((x) => h("th", {class: "l"}, x)))));
      t.append(h("tbody", {}, S.mods(kind).map((m) => h("tr", {},
        h("td", {class: "l mono"}, m.id), h("td", {class: "l"}, m.title),
        h("td", {class: "l mono small"}, m.script ? [m.script, m.script_found ? "" : h("span", {style: {color: "var(--bad)"}}, " (not found)")] :
          Object.entries(m.modes || {}).map(([k, v]) => h("div", {}, k + ": " + v.cmd.join(" ")))),
        h("td", {class: "l"}, m.input), h("td", {class: "l mono small muted"}, m.source)))));
      main.append(h("div", {class: "card"}, h("h2", {}, {stage: "Pipeline stages", plot: "Plots", tool: "Tools"}[kind]), h("div", {class: "tablewrap"}, t)));
    }
    main.append(h("div", {class: "card"}, h("h2", {}, "Example: add a plot script"),
      h("pre", {class: "log"}, `# modules.d/my_plot.toml
id          = "my_plot"
kind        = "plot"                       # plot | stage | tool
title       = "My new plot"
script      = "plot_my_thing.py"           # in scripts/, or an absolute path; its argparse options become the form
input       = "merged_csv"                 # merged_csv | comp_csv | any_csv | config | poscar | folder | none
watch       = ["input_dir", "--out-dir"]   # where its output files appear
needs       = ["matplotlib"]

# a pipeline stage instead (e.g. swap the relaxation for your own GPUMD driver):
# id = "compute"   kind = "stage"   title = "Relax + F (GPUMD)"
# [modes.local]  cmd = ["{python}", "/path/my_driver.py", "{config}"]
# [modes.slurm]  cmd = ["{python}", "/path/my_driver.py", "{config}", "--submit"]
# [[options]]    flag = "--unfinished"  label = "only unfinished"  value = false`)));
  }
})();
