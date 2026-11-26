/* Jobs page: commands started from the app (local, detached) and your SLURM jobs. */
(function () {
  const S = Studio, h = S.h;
  S.registerPage({id: "jobs", title: "Jobs", order: 20, icon: "server", render});

  async function render(main, params) {
    const top = h("div", {});
    const list = h("div", {class: "card"});
    const slurm = h("div", {class: "card"});
    main.append(top, list, slurm);
    let view = null;
    const showJob = (id, title) => { if (view && view.stop) view.stop(); top.innerHTML = ""; view = S.jobView(id, {title, inJobs: true}); top.append(view); window.scrollTo(0, 0); };
    if (params.id) showJob(params.id, "Job " + params.id);
    await drawList();
    await drawSlurm();
    S.every(5000, drawList);
    S.every(30000, drawSlurm);

    async function drawList() {
      const r = await S.api("/api/jobs", undefined, {quiet: true});
      if (!r.ok) return;
      list.innerHTML = "";
      list.append(h("h2", {}, "Started from the app"), h("div", {class: "small muted"}, "Jobs run detached: they keep running if you close the page or stop the server."));
      const t = h("table", {});
      t.append(h("thead", {}, h("tr", {}, ["state", "job", "started", "duration"].map((x) => h("th", {class: "l"}, x)))));
      const tb = h("tbody", {});
      for (const j of r.jobs) {
        const pill = h("span", {class: "pill " + ({done: "ok", failed: "bad", lost: "warn", cancelled: "warn"}[j.state] || "run")}, j.state);
        tb.append(h("tr", {}, h("td", {class: "l"}, pill),
          h("td", {class: "l"}, h("a", {href: "javascript:void 0", onclick: () => showJob(j.id, j.title)}, j.title),
            j.slurm_ids.length ? h("span", {class: "small muted"}, " · SLURM " + j.slurm_ids.join(", ")) : null),
          h("td", {class: "l small"}, S.ago(j.start)),
          h("td", {class: "l small"}, j.end ? S.dur(j.end - j.start) : j.state === "running" ? S.dur(Date.now() / 1000 - j.start) : ""),
          h("td", {}, !["running", "starting"].includes(j.state)
            ? h("button", {class: "small", title: "delete this job's record and log", onclick: async () => { await S.api("/api/jobs/delete", {id: j.id}); drawList(); }}, "✕")
            : j.state === "running" ? h("button", {class: "small danger", onclick: async () => {
              const c = await S.api("/api/jobs/cancel", {id: j.id}); S.toast(c.ok ? "stop signal sent" : "could not stop"); setTimeout(drawList, 1200); }}, "Stop") : null)));
      }
      t.append(tb);
      list.append(r.jobs.length ? h("div", {class: "tablewrap"}, t) : h("div", {class: "muted"}, "none yet"));
    }
    async function drawSlurm() {
      const r = await S.api("/api/slurm", undefined, {quiet: true});
      slurm.innerHTML = "";
      slurm.append(h("h2", {}, "SLURM"));
      if (!r.ok || !r.available) { slurm.append(h("div", {class: "muted"}, "squeue not found on this machine (SLURM view is available on the cluster login node).")); return; }
      if (r.error) slurm.append(h("div", {class: "notice warn small"}, r.error));
      const sel = new Set();
      const t = h("table", {});
      t.append(h("thead", {}, h("tr", {}, ["", "id", "name", "state", "time", "limit", "nodes", "reason"].map((x) => h("th", {class: "l"}, x)))));
      const tb = h("tbody", {});
      for (const j of r.jobs) {
        const c = h("input", {type: "checkbox", onchange: () => (c.checked ? sel.add(j.id) : sel.delete(j.id))});
        const known = r.known_ids.includes(j.id.split("_")[0]);
        tb.append(h("tr", {}, h("td", {class: "l"}, c), h("td", {class: "l mono"}, j.id, known ? " ★" : ""), h("td", {class: "l"}, j.name),
          h("td", {class: "l"}, h("span", {class: "pill " + (j.state === "RUNNING" ? "run" : j.state === "PENDING" ? "warn" : "")}, j.state)),
          h("td", {class: "l"}, j.time), h("td", {class: "l"}, j.limit), h("td", {class: "l"}, j.nodes), h("td", {class: "l small"}, j.reason)));
      }
      t.append(tb);
      slurm.append(h("div", {class: "small muted"}, `${r.jobs.length} queued / running jobs of yours (★ = submitted from the app). Refreshed every 30 s.`),
        h("div", {class: "tablewrap"}, t),
        h("div", {class: "row", style: {marginTop: "6px"}}, h("button", {class: "danger", onclick: async () => {
          if (!sel.size) return S.toast("tick jobs first");
          const o = await S.api("/api/slurm/cancel", {ids: [...sel]});
          S.toast("scancel: " + o.output); drawSlurm();
        }}, "scancel selected")));
      if ((r.history || []).length) {
        const t2 = h("table", {});
        t2.append(h("thead", {}, h("tr", {}, ["id", "name", "state", "elapsed", "exit"].map((x) => h("th", {class: "l"}, x)))));
        t2.append(h("tbody", {}, r.history.slice(-100).reverse().map((j) => h("tr", {}, ["id", "name", "state", "elapsed", "exit"].map((k) => h("td", {class: "l small"}, j[k]))))));
        slurm.append(h("h3", {}, "Submitted from the app (sacct)"), h("div", {class: "tablewrap"}, t2));
      }
    }
  }
})();
