"""
modules.py - the plug-in registry.  Everything the app can run is a module described by one small TOML file.

Where modules are read from (later folders override earlier ones with the same id):
  1. <app>/studio/modules/*.toml          shipped with the app
  2. <app>/modules.d/*.toml                your own additions / replacements (kept when you update the app)
  3. $FO_STUDIO_MODULES (folders separated by ':')
A file with `enabled = false` removes a module of the same id.

Module file (all keys except id / kind / title optional):

  id          = "f_vs_x_html"             unique name
  kind        = "plot"                    plot | stage | tool
  title       = "F vs x (interactive)"
  description = "..."
  order       = 10                        position in its page
  script      = "plot_F_vs_x_html.py"     file in the scripts folder (or an absolute path); its argparse options
                                          become the form on the page (no other code needed)
  input       = "merged_csv"              what the first positional argument is chosen from:
                                          merged_csv | comp_csv | any_csv | config | poscar | folder | none
  args        = []                        fixed arguments placed after the input
  watch       = ["input_dir", "--out", "--out-dir"]   where the module's output files appear
                                          ("input_dir", "workdir", "rundir" or an option whose value is a path)
  needs       = ["matplotlib"]            python modules checked on the page (warning only)
  [from_run]                              arguments filled in from the run (campaign.toml + run_record.json) that
  "--dopant" = "{dopant}"                 holds the input; keys: host, dopant, T, starting_phase, name.  Skipped when
  "--label" = "parent={starting_phase} (initial)"   unknown or when you give the flag yourself

  Stages (kind = "stage", the pipeline page) run commands instead of a form:
  [modes.local]  cmd = ["{python}", "{scripts}/campaign.py", "{config}", "prepare"]   label = "Run here"
  [modes.slurm]  cmd = [...]                                                          label = "Submit to SLURM"
  [[options]]    flag = "--unfinished"  label = "..."  value = false  modes = ["slurm"]
                 (value = false: a check box; value = true: a text field whose words follow the flag;
                  flag = "" with value = true: free arguments)
                 default = true / "text"   start value (also used when a caller leaves the option out)
                 help = "..."              one line under the option
                 choices = ["a", "b"]      with value = true: a drop-down, the word must be one of them
                 kind = "choice", key = "name", [[options.choices]] flag = "--x" label = "..." help = "..."
                     radio buttons: exactly one flag of the list ("" = none); a choice may have
                     confirm = "..." (asked before it runs) and block_if_compute_running = true
  [modes.X] confirm = "..."     asked before this mode runs (e.g. SLURM submissions)
  Placeholders: {python} {scripts} {config} {workdir} {rundir} {studio}
"""
from __future__ import annotations

import os
import re
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:                     # Python < 3.11
    try:
        import tomli as tomllib  # type: ignore
    except ModuleNotFoundError:
        raise SystemExit("Python < 3.11 needs the package tomli:  pip install tomli")

HERE = Path(__file__).resolve().parent
APP = HERE.parent
KINDS = ("plot", "stage", "tool")


def module_dirs():
    dirs = [HERE / "modules", APP / "modules.d"]
    for d in os.environ.get("FO_STUDIO_MODULES", "").split(":"):
        if d:
            dirs.append(Path(d).expanduser())
    return dirs


DISABLED = {}      # id -> file that disabled it (last load)


def load_modules():
    mods, errors = {}, []
    DISABLED.clear()
    for d in module_dirs():
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*.toml")):
            try:
                with p.open("rb") as f:
                    m = tomllib.load(f)
            except Exception as e:      # noqa: BLE001
                errors.append(f"{p}: {e}")
                continue
            mid = m.get("id") or p.stem
            if m.get("enabled", True) is False:
                mods.pop(mid, None)
                DISABLED[mid] = str(p)
                continue
            DISABLED.pop(mid, None)
            if m.get("kind") not in KINDS:
                errors.append(f"{p}: kind must be one of {KINDS}")
                continue
            if not re.fullmatch(r"[\w\-]+", mid):
                errors.append(f"{p}: id {mid!r} may only use letters, digits, _ and -")
                continue
            try:
                m["order"] = float(m.get("order", 100))
            except (TypeError, ValueError):
                errors.append(f"{p}: order must be a number")
                continue
            bad = [k for k in ("args", "watch", "needs") if k in m and not isinstance(m[k], list)]
            bad += [k for k in ("script", "title", "input", "description") if k in m and not isinstance(m[k], str)]
            bad += [k for k in ("from_run", "watch_defaults") if k in m and not isinstance(m[k], dict)]
            if bad or ("modes" in m and not isinstance(m["modes"], dict)) or ("options" in m and not isinstance(m["options"], list)):
                errors.append(f"{p}: wrong type for {bad or ['modes/options']}")
                continue
            if any(not isinstance(md, dict) or not isinstance(md.get("cmd"), list) for md in (m.get("modes") or {}).values()):
                errors.append(f"{p}: every [modes.X] needs cmd = [...]")
                continue
            m["id"] = mid
            m["source"] = str(p)
            m.setdefault("title", mid)
            m.setdefault("order", 100)
            m.setdefault("input", "none")
            m.setdefault("args", [])
            m.setdefault("watch", ["input_dir"])
            mods[mid] = m
    return dict(sorted(mods.items(), key=lambda kv: (kv[1]["kind"], kv[1]["order"], kv[0]))), errors


def option_key(o):
    """the name a stage option has in the page's options dict"""
    return o.get("key") or o.get("flag") or "__free__"


def script_path(mod, scripts: Path):
    s = mod.get("script")
    if not s:
        return None
    p = Path(s).expanduser()
    return p if p.is_absolute() else scripts / p


def expand(tokens, ctx):
    """replace {python} {scripts} ... in a command template"""
    out = []
    for t in tokens:
        t = str(t)
        for k, v in ctx.items():
            t = t.replace("{" + k + "}", str(v))
        out.append(t)
    return out
