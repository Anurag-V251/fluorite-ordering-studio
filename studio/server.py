"""
server.py - the web server (Python standard library only).  Start it with  python3 start_studio.py

Security: binds to 127.0.0.1 by default and every request needs the random token printed at start
(kept in a cookie after the first visit), so other users of a shared login node cannot use it.
It runs commands as you, on this machine, only the scripts of registered modules.
"""
from __future__ import annotations

import json
import mimetypes
import os
import re
import secrets
import shlex
import subprocess
import sys
import threading
import time
import shutil
import traceback
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, quote, urlparse

from . import jobs as J
from . import modules as MOD
from . import runs as R
from . import scan as S
from . import single as SR
from . import xrd as X


class BadRequest(Exception):
    pass

HERE = Path(__file__).resolve().parent
APP = HERE.parent
STATIC = HERE / "static"


class State:
    def __init__(self, scripts: Path, state_dir: Path, token: str, start_dir: Path):
        self.scripts = scripts
        self.state_dir = state_dir
        self.token = token
        self.start_dir = start_dir
        self.jobs = J.Jobs(state_dir)
        self.argspec_cache = {}
        self.lock = threading.Lock()
        self.prefs_path = state_dir / "prefs.json"
        self.prefs_lock = threading.Lock()
        self.launch_lock = threading.Lock()
        self.scan_locks = {}
        self.scan_cache = {}
        self.port = 0

    def prefs(self):
        try:
            return json.loads(self.prefs_path.read_text())
        except (OSError, ValueError):
            return {}

    def save_prefs(self, d):
        with self.prefs_lock:
            p = self.prefs()
            p.update(d)
            tmp = self.prefs_path.with_suffix(".tmp")
            tmp.write_text(json.dumps(p, indent=1))
            os.replace(tmp, self.prefs_path)
            return p


ST: State | None = None


def argspec(script: Path):
    key = (str(script), script.stat().st_mtime if script.is_file() else 0)
    with ST.lock:
        if key in ST.argspec_cache:
            return ST.argspec_cache[key]
    if not script.is_file():
        return {"ok": False, "error": f"{script} not found"}
    try:
        r = subprocess.run([sys.executable, str(HERE / "argspec.py"), str(script)], capture_output=True, text=True,
                           timeout=120, cwd=str(script.parent))
        out = json.loads(r.stdout.strip().split("\n")[-1]) if r.stdout.strip() else \
            {"ok": False, "error": r.stderr.strip()[-2000:] or "no output"}
    except (subprocess.TimeoutExpired, ValueError) as e:
        out = {"ok": False, "error": str(e)}
    if not out.get("ok") and "No module named" in out.get("error", "") + out.get("output", ""):
        out["hint"] = "a Python package the script imports is missing in this environment"
    if out.get("ok"):                           # failures (missing package, timeout) are retried next time
        with ST.lock:
            ST.argspec_cache[key] = out
    return out


def module_ctx(config=None):
    ctx = {"python": sys.executable, "scripts": str(ST.scripts), "studio": str(HERE), "app": str(APP),
           "output": str(ST.state_dir / "outputs")}
    if config:
        cp = Path(config).expanduser().resolve()
        ctx["config"] = str(cp)
        ctx["rundir"] = str(cp.parent)
        info = R.read_run(cp.parent)
        if info.get("workdir"):
            ctx["workdir"] = info["workdir"]
        # a run uses ITS scripts folder (frozen copy = reproducibility).  Code is only run from the app's own
        # scripts/ or from a folder whose scripts belong to you and nobody else can change (opening someone
        # else's run folder must not run their code as you)
        sp = info.get("scripts")
        if sp and Path(sp, "campaign.py").is_file():
            spp = Path(sp).resolve()
            why = None if spp == ST.scripts else scripts_untrusted(spp)
            if why:
                ctx["untrusted_scripts"] = str(spp)
                ctx["untrusted_why"] = why
            else:
                ctx["scripts"] = str(spp)
    return ctx


def run_info_of(path):
    """{starting_phase, host, dopant, T, name} of the run folder (campaign.toml) that holds `path`, 4 levels up at
    most; keys that are not known are left out"""
    p = Path(path).expanduser().resolve()
    for d in list(p.parents)[:4]:
        if (d / "campaign.toml").is_file() or (d / "run_record.json").is_file():
            out = {}
            r = R.read_run(d)
            st = ((r.get("config") or {}).get("structure") or {})
            for k in ("host", "dopant"):
                if st.get(k):
                    out[k] = str(st[k])
            T = ((r.get("config") or {}).get("free_energy") or {}).get("T")
            if T is not None:
                out["T"] = f"{float(T):g}"
            rec = r.get("record") or {}
            if str(rec.get("starting_phase") or "").strip():
                out["starting_phase"] = str(rec["starting_phase"]).strip()
            if rec.get("name"):
                out["name"] = str(rec["name"])
            return out
    return {}


def _private_group(gid):
    """True if the group is the user's own private group (umask 002 clusters): nobody else can be in it"""
    try:
        import grp
        import pwd
        g = grp.getgrgid(gid)
        return g.gr_name == pwd.getpwuid(os.getuid()).pw_name and not [m for m in g.gr_mem if m != g.gr_name]
    except (ImportError, KeyError):
        return False


def file_untrusted(p: Path):
    """None if p is yours and nobody else can change it, else the reason"""
    try:
        st = p.stat()
    except OSError as e:
        return str(e)
    if st.st_uid != os.getuid():
        return f"{p} belongs to another user"
    if st.st_mode & 0o002 or (st.st_mode & 0o020 and not _private_group(st.st_gid)):
        return f"{p} can be changed by other users (fix: chmod -R go-w {p.parent if p.is_file() else p})"
    return None


def scripts_untrusted(d: Path):
    """None if d and every .py in it are yours and nobody else can change them"""
    for p in [d, *sorted(d.glob("*.py"))]:
        why = file_untrusted(p)
        if why:
            return why
    return None


def resolve_script(mod, ctx=None):
    s = mod.get("script")
    if not s:
        return None
    s = MOD.expand([s], ctx or module_ctx())[0]
    p = Path(s).expanduser()
    return p if p.is_absolute() else Path((ctx or module_ctx())["scripts"]) / p


def needs_report(mod):
    out = {}
    import importlib.util
    for n in mod.get("needs", []):
        out[n] = importlib.util.find_spec(n) is not None
    return out


# ============================================================================ handler
class H(SR.SingleRunsAPI, BaseHTTPRequestHandler):
    timeout = 120                                # idle connections do not hold a thread forever
    server_version = "FluoriteStudio/" + R.APP_VERSION

    def log_message(self, fmt, *args):            # quiet; errors still go to stderr
        if os.environ.get("FO_STUDIO_DEBUG"):
            super().log_message(fmt, *args)

    # ---------------------------------------------------------------- helpers
    def _authed(self, q):
        tok = (q.get("token") or [None])[0] or self.headers.get("X-Token")
        if not tok and self.headers.get("Cookie"):
            c = SimpleCookie(self.headers.get("Cookie"))
            if f"fo_token_{ST.port}" in c:
                tok = c[f"fo_token_{ST.port}"].value
        return tok is not None and secrets.compare_digest(tok, ST.token)

    def _send(self, code, body: bytes, ctype="application/json", extra=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, default=str).encode())

    def _err(self, msg, code=400):
        self._json({"ok": False, "error": msg}, code)

    def _body(self, limit=10 * 1024 * 1024):
        try:
            n = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            n = -1
        if n < 0 or n > limit:
            raise BadRequest(f"request body length {self.headers.get('Content-Length')} not accepted")
        return self.rfile.read(n) if n else b""

    # ---------------------------------------------------------------- routing
    def do_GET(self):
        self._route("GET")

    def do_POST(self):
        self._route("POST")

    def _route(self, method):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        try:
            if not self._authed(q):
                if u.path == "/" or u.path.startswith("/static"):
                    return self._send(403, b"<h3>Fluorite Ordering Studio</h3><p>Open the link with the token that "
                                           b"was printed when the server started.</p>", "text/html")
                return self._err("missing or wrong token", 403)
            if u.path == "/" and "token" in q:            # first visit: keep the token in a cookie
                return self._send(303, b"", "text/plain", {
                    "Location": "/", "Set-Cookie": f"fo_token_{ST.port}={ST.token}; Path=/; HttpOnly; SameSite=Strict"})
            if u.path == "/":
                return self._static("index.html")
            if u.path.startswith("/static/"):
                return self._static(u.path[len("/static/"):])
            fn = getattr(self, f"api_{method.lower()}_" + u.path.strip("/").replace("api/", "", 1).replace("/", "_"), None)
            if fn is None:
                return self._err(f"no such endpoint {method} {u.path}", 404)
            if method == "POST":
                if self.headers.get("X-Requested-With") != "studio":      # simple CSRF guard
                    return self._err("bad request origin", 403)
            return fn(q)
        except BadRequest as e:
            return self._err(str(e), 400)
        except (BrokenPipeError, ConnectionResetError):
            return None
        except Exception as e:      # noqa: BLE001
            traceback.print_exc()
            return self._err(f"{type(e).__name__}: {e}", 500)

    def _static(self, rel):
        p = (STATIC / rel).resolve()
        if STATIC not in p.parents and p != STATIC or not p.is_file():
            return self._err("not found", 404)
        ctype = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        if p.suffix == ".js":
            ctype = "text/javascript"
        self._send(200, p.read_bytes(), ctype + ("; charset=utf-8" if ctype.startswith("text") else ""))

    # ---------------------------------------------------------------- info / modules
    def api_get_info(self, q):
        mods, errs = MOD.load_modules()
        self._json({"ok": True, "version": R.APP_VERSION, "scripts": str(ST.scripts), "state_dir": str(ST.state_dir),
                    "app": str(APP), "home": str(Path.home()), "start_dir": str(ST.start_dir),
                    "python": sys.executable, "slurm": J.slurm_available(), "module_errors": errs,
                    "module_dirs": [str(d) for d in MOD.module_dirs()],
                    "pages": sorted(p.name for p in (STATIC / "pages").glob("*.js")),
                    "prefs": ST.prefs(), "host": os.uname().nodename})

    def api_post_prefs(self, q):
        self._json({"ok": True, "prefs": ST.save_prefs(json.loads(self._body() or b"{}"))})

    def api_get_modules(self, q):
        mods, errs = MOD.load_modules()
        out = []
        for m in mods.values():
            d = {k: v for k, v in m.items()}
            d["needs_ok"] = needs_report(m)
            sp = resolve_script(m)
            d["script_path"] = str(sp) if sp else None
            d["script_found"] = bool(sp and sp.is_file())
            out.append(d)
        self._json({"ok": True, "modules": out, "errors": errs, "disabled": dict(MOD.DISABLED),
                    "module_dirs": [str(d) for d in MOD.module_dirs()]})

    def api_get_argspec(self, q):
        mods, _ = MOD.load_modules()
        m = mods.get((q.get("module") or [""])[0])
        if not m:
            return self._err("unknown module")
        ctx = module_ctx((q.get("config") or [None])[0])
        sp = resolve_script(m, ctx)
        if not sp:
            return self._json({"ok": True, "actions": [], "description": m.get("description", "")})
        self._json(argspec(sp))

    # ---------------------------------------------------------------- files
    def api_get_ls(self, q):
        p = Path((q.get("path") or [str(ST.start_dir)])[0]).expanduser()
        try:
            p = p.resolve()
            ents = []
            with os.scandir(p) as it:
                for e in it:
                    if e.name.startswith(".") and (q.get("hidden") or ["0"])[0] != "1":
                        continue
                    try:
                        isdir = e.is_dir()
                        st = e.stat()
                    except OSError:
                        continue
                    ents.append({"name": e.name, "dir": isdir, "size": None if isdir else st.st_size,
                                 "mtime": st.st_mtime})
        except OSError as e:
            return self._err(str(e))
        ents.sort(key=lambda x: (not x["dir"], x["name"].lower()))
        self._json({"ok": True, "path": str(p), "parent": str(p.parent), "entries": ents[:5000],
                    "truncated": len(ents) > 5000})

    def api_get_file(self, q):
        p = Path((q.get("path") or [""])[0]).expanduser()
        if not p.is_file():
            return self._err("not a file", 404)
        ctype = mimetypes.guess_type(p.name)[0] or "text/plain"
        if p.suffix.lower() in (".csv", ".txt", ".log", ".out", ".toml", ".sh", ".json", "") or p.name.startswith("POSCAR"):
            ctype = "text/plain"
        extra = {}
        if (q.get("download") or ["0"])[0] == "1":
            extra["Content-Disposition"] = f"attachment; filename*=UTF-8''{quote(p.name)}"
        if not (ctype in ("text/plain", "text/csv", "image/png", "image/jpeg", "image/gif", "application/pdf")):
            # anything that a browser could run (html, svg, xhtml, xml, ...) runs sandboxed, never with the app's rights
            extra["Content-Security-Policy"] = "sandbox allow-scripts allow-downloads allow-popups"
        size = p.stat().st_size
        if size > 512 * 1024 * 1024:
            return self._err("file larger than 512 MB: open it on the machine itself", 413)
        try:
            f = p.open("rb")
        except OSError as e:
            return self._err(str(e), 403)
        with f:                                 # streamed: a big file never sits in memory
            self.send_response(200)
            self.send_header("Content-Type", ctype + ("; charset=utf-8" if ctype.startswith("text") else ""))
            self.send_header("Content-Length", str(size))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            for k, v in extra.items():
                self.send_header(k, v)
            self.end_headers()
            left = size                         # exactly the announced length, even if the file grows
            while left > 0:
                buf = f.read(min(left, 1024 * 1024))
                if not buf:
                    break
                self.wfile.write(buf)
                left -= len(buf)

    def api_get_text(self, q):
        p = Path((q.get("path") or [""])[0]).expanduser()
        try:
            n = min(max(int((q.get("tail") or ["200000"])[0]), 0), 5_000_000)
        except ValueError:
            return self._err("tail must be a whole number")
        if not p.is_file():
            return self._err("not a file", 404)
        size = p.stat().st_size
        with p.open("rb") as f:
            f.seek(max(0, size - n))
            txt = f.read().decode(errors="replace")
        self._json({"ok": True, "path": str(p), "size": size, "text": txt, "truncated": size > n})

    def api_post_upload(self, q):
        d = Path((q.get("dir") or [""])[0]).expanduser()
        name = Path((q.get("name") or [""])[0]).name
        if not d.is_dir() or not name:
            return self._err("upload needs an existing folder and a file name")
        dst = d / name
        if dst.exists() and (q.get("overwrite") or ["0"])[0] != "1":
            return self._err(f"{dst} exists")
        try:
            n = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            return self._err("bad Content-Length")
        if n > 1024 * 1024 * 1024:
            return self._err("upload larger than 1 GB: copy it with scp instead", 413)
        if n < 0:
            return self._err("bad Content-Length")
        tmp = dst.with_name(f".upload-{secrets.token_hex(6)}")
        try:
            with tmp.open("wb") as f:          # streamed to disk in pieces
                left = n
                while left > 0:
                    buf = self.rfile.read(min(left, 1024 * 1024))
                    if not buf:
                        break
                    f.write(buf)
                    left -= len(buf)
            if left > 0:
                return self._err("upload cut off")
            os.replace(tmp, dst)
        finally:
            if tmp.exists():
                tmp.unlink()
        self._json({"ok": True, "path": str(dst)})

    def api_post_mkdir(self, q):
        d = Path(json.loads(self._body())["path"]).expanduser()
        d.mkdir(parents=True, exist_ok=True)
        self._json({"ok": True, "path": str(d.resolve())})

    def api_get_poscar(self, q):
        self._json(R.peek_poscar((q.get("path") or [""])[0]))

    # ---------------------------------------------------------------- scan
    def api_get_scan(self, q):
        path = (q.get("path") or [""])[0]
        key = str(Path(path).expanduser())
        with ST.lock:
            lk = ST.scan_locks.setdefault(key, threading.Lock())
        with lk:                                 # one scan per folder at a time; repeated calls within 5 s reuse it
            hit = ST.scan_cache.get(key)
            if hit and time.time() - hit[0] < 5 and (q.get("fresh") or ["0"])[0] != "1":
                res = hit[1]
            else:
                res = S.scan(path)
                with ST.lock:                   # a short-lived cache: keep only fresh results, at most 6
                    now = time.time()
                    ST.scan_cache[key] = (now, res)
                    for k in sorted(ST.scan_cache, key=lambda k: ST.scan_cache[k][0])[:-6]:
                        ST.scan_cache.pop(k, None)
                    for k in [k for k, v in ST.scan_cache.items() if now - v[0] > 30]:
                        ST.scan_cache.pop(k, None)
                    if len(ST.scan_locks) > 200:
                        ST.scan_locks = {k: v for k, v in ST.scan_locks.items() if k in ST.scan_cache or k == key}
        if res.get("ok") and res.get("layout") != "none":
            ST.save_prefs({"last_results": res["root"]})
            rec = ST.prefs().get("recent_results", [])
            rec = [res["root"]] + [r for r in rec if r != res["root"]]
            ST.save_prefs({"recent_results": rec[:12]})
        self._json(res)

    def api_get_logs(self, q):
        wd = Path((q.get("workdir") or [""])[0]).expanduser()
        out = []
        for sub in ("logs", "slurm"):
            d = wd / sub
            if d.is_dir():
                for e in os.scandir(d):
                    if e.is_file():
                        st = e.stat()
                        out.append({"path": e.path, "name": f"{sub}/{e.name}", "size": st.st_size, "mtime": st.st_mtime})
        out.sort(key=lambda x: -x["mtime"])
        self._json({"ok": True, "logs": out[:500]})

    # ---------------------------------------------------------------- run a module with a form
    def api_post_run(self, q):
        b = json.loads(self._body())
        mods, _ = MOD.load_modules()
        m = mods.get(b.get("module"))
        if not m:
            return self._err("unknown module")
        ctx = module_ctx(b.get("config"))
        sp = resolve_script(m, ctx)
        if not sp or not sp.is_file():
            return self._err(f"script of module {m['id']} not found: {sp}")
        args = [str(a) for a in b.get("args", [])]
        inp = b.get("input")
        auto = []
        if m.get("from_run") and inp:           # values taken from the run that holds the input (dopant, phase ...)
            info = run_info_of(inp)
            for flag, tpl in m["from_run"].items():
                if flag in args and not (flag == "--label" and not any(x.startswith("parent=") for x in args)):
                    continue                    # the user's own value wins
                try:
                    val = tpl.format(**info)
                except (KeyError, IndexError, ValueError):
                    continue                    # not a run folder, or the value is not known
                auto += [flag, val]
        cmd = [sys.executable, str(sp)] + ([inp] if inp else []) + MOD.expand(m.get("args", []), ctx) + auto + args
        cwd = Path(b.get("cwd") or (Path(inp).parent if inp else ST.start_dir)).expanduser()
        if not cwd.is_dir():
            cwd = Path(inp).parent if inp and Path(inp).parent.is_dir() else ST.start_dir
        watch = []
        for w in m.get("watch", []):
            if w in ("input_dir", "R:input_dir") and inp:     # R: = also sub-folders (3 levels)
                watch.append(("R:" if w.startswith("R:") else "") +
                             str(Path(inp).parent if Path(inp).is_file() else Path(inp)))
            elif w in ("workdir", "rundir") and w in ctx:
                watch.append(Path(ctx[w]))
            elif w.startswith("-"):
                val = None
                if w in args and args.index(w) + 1 < len(args):
                    val = args[args.index(w) + 1]
                elif (m.get("watch_defaults") or {}).get(w):
                    val = MOD.expand([m["watch_defaults"][w]], ctx)[0]
                if val:
                    v = Path(val).expanduser()
                    v = v if v.is_absolute() else cwd / v
                    watch.append("R:" + str(v if (v.is_dir() or not v.suffix) else v.parent))
        watch.append(cwd)
        jid = ST.jobs.start(cmd, cwd, title=f"{m['title']}" + (f" - {Path(inp).name}" if inp else ""),
                            module=m["id"], watch=sorted({str(w) for w in watch}),
                            env={"FO_STUDIO_OUT": str(ST.state_dir / "outputs")},
                            meta={"kind": m["kind"], "input": inp})
        self._json({"ok": True, "job": jid, "cmd": shlex.join(cmd)})

    # ---------------------------------------------------------------- pipeline stage
    def api_post_stage(self, q):
        b = json.loads(self._body())
        mods, _ = MOD.load_modules()
        m = mods.get(b.get("module"))
        if not m or m["kind"] != "stage":
            return self._err("unknown stage")
        mode = b.get("mode", "local")
        md = (m.get("modes") or {}).get(mode)
        if not md:
            return self._err(f"stage {m['id']} has no mode {mode}")
        if not b.get("config") or not Path(b["config"]).expanduser().is_file():
            return self._err("choose the run (campaign.toml) first")
        ctx = module_ctx(b["config"])
        if ctx.get("untrusted_scripts"):
            other = "another user" in ctx["untrusted_why"]
            return self._err(f"not running this run's scripts ({ctx['untrusted_scripts']}): {ctx['untrusted_why']}." +
                             (" To use the run, 'copy settings to a new run' (it freezes the app's own scripts)."
                              if other else ""))
        why = file_untrusted(Path(ctx["config"]))   # its [cluster] setup / python lines go into SLURM scripts
        if why:
            return self._err(f"not running stages of {ctx['config']}: {why}." +
                             (" Use 'copy settings to a new run' and check the settings (setup, python) first."
                              if "another user" in why else ""))
        cmd = MOD.expand(md["cmd"], ctx)
        for o in m.get("options", []):
            if o.get("modes") and mode not in o["modes"]:
                continue
            v = (b.get("options") or {}).get(o.get("flag") or "__free__")
            if o.get("value"):
                if v and str(v).strip():
                    try:
                        words = shlex.split(str(v))
                    except ValueError as e:
                        return self._err(f"option {o.get('flag') or 'arguments'}: {e}")
                    cmd += ([o["flag"]] if o.get("flag") else []) + words
            elif v:
                cmd.append(o["flag"])
        if b.get("preview"):
            return self._json({"ok": True, "cmd": shlex.join(cmd)})
        rundir = Path(ctx["rundir"])
        wd = Path(ctx.get("workdir", rundir))
        watch = [str(rundir), str(wd)] + ["R:" + str(wd / d) for d in ("slurm", "results", "logs", "chunks")] + \
            ["R:" + str(d) for d in sorted(wd.glob("rerun_*"))[-3:]]
        with ST.launch_lock:                    # check + start together: two quick clicks launch once
            busy = ST.jobs.running(module=m["id"], config=ctx["config"])
            if busy:
                return self._err(f"'{m['title']}' is already running for this run (job {busy[0]['id']}); stop it first")
            jid = ST.jobs.start(cmd, rundir, title=f"{m['title']} [{md.get('label', mode)}] - {rundir.name}",
                                module=m["id"], watch=watch, meta={"kind": "stage", "mode": mode,
                                                                   "config": ctx["config"]})
        R.add_history(rundir, {"stage": m["id"], "mode": mode, "cmd": shlex.join(cmd), "job": jid})
        self._json({"ok": True, "job": jid, "cmd": shlex.join(cmd)})

    # ---------------------------------------------------------------- jobs
    def api_get_jobs(self, q):
        jid = (q.get("id") or [None])[0]
        if jid:
            try:
                tail = min(max(int((q.get("tail") or ["60000"])[0]), 0), 5_000_000)
            except ValueError:
                return self._err("tail must be a whole number")
            g = ST.jobs.get(jid, log_tail=tail)
            return self._json({"ok": True, "job": g}) if g else self._err("no such job", 404)
        self._json({"ok": True, "jobs": ST.jobs.list()})

    def api_post_jobs_cancel(self, q):
        self._json({"ok": ST.jobs.cancel(json.loads(self._body())["id"])})

    def api_post_jobs_delete(self, q):
        self._json({"ok": ST.jobs.delete(json.loads(self._body())["id"])})

    def api_get_slurm(self, q):
        sq = J.squeue()
        ids = []
        for j in ST.jobs.list(100):
            ids += j.get("slurm_ids", [])
        ids = sorted(set(ids), key=int)
        sq["history"] = J.sacct(ids[-200:])
        sq["known_ids"] = ids
        self._json({"ok": True, **sq})

    def api_post_slurm_cancel(self, q):
        ids = [i for i in json.loads(self._body()).get("ids", []) if re.fullmatch(r"\d+(_\d+)?", str(i))]
        self._json({"ok": True, "output": J.scancel(ids) if ids else "no ids"})

    # ---------------------------------------------------------------- runs
    def api_post_runs_create(self, q):
        spec = json.loads(self._body())
        try:
            out = R.create_run(spec, ST.scripts)
        except ValueError as e:
            return self._err(str(e))
        rec = ST.prefs().get("recent_runs", [])
        ST.save_prefs({"recent_runs": [out["config"]] + [r for r in rec if r != out["config"]][:11]})
        self._json({"ok": True, **out})

    def api_get_runs_list(self, q):
        """run folders (campaign.toml) up to 2 levels below the start folder (or ?path=), newest first"""
        base = Path((q.get("path") or [str(ST.start_dir)])[0]).expanduser()
        out = []
        for r in R.find_runs(base, depth=2)[:200]:
            rec = r.get("record") or {}
            out.append({"folder": r["folder"], "config": r.get("config_path"), "workdir": r.get("workdir"),
                        "name": rec.get("name", ""), "starting_phase": rec.get("starting_phase", ""),
                        "created": rec.get("created", ""),
                        "mtime": Path(r["config_path"]).stat().st_mtime if r.get("config_path") else 0})
        out.sort(key=lambda x: -x["mtime"])
        self._json({"ok": True, "base": str(base), "runs": out})

    def api_get_runs_read(self, q):
        f = Path((q.get("folder") or [""])[0]).expanduser()
        if f.is_file():
            f = f.parent
        info = R.read_run(f)
        if info.get("config_path"):
            rec = ST.prefs().get("recent_runs", [])
            ST.save_prefs({"recent_runs": [info["config_path"]] + [r for r in rec if r != info["config_path"]][:11]})
        self._json({"ok": bool(info.get("config")), **info,
                    **({} if info.get("config") else {"error": info.get("error") or "no campaign config (*.toml) in this folder"})})

    def api_get_runs_check(self, q):
        f = Path((q.get("folder") or [""])[0]).expanduser()
        self._json(R.check_reproducible(f.parent if f.is_file() else f))

    def api_get_runs_defaults(self, q):
        self._json({"ok": True, "defaults": R.DEFAULTS})

    # ---------------------------------------------------------------- XRD page (scripts/xrd_average.py)
    def api_get_xrd_rows(self, q):
        src = (q.get("source") or [""])[0]
        if not src or not Path(src).expanduser().exists():
            return self._err("choose a results folder (run folder, workdir) or a result CSV")
        try:
            out = X.rows(str(Path(src).expanduser().resolve()), (q.get("T") or [None])[0] or None,
                         (q.get("comp") or [None])[0])
        except (ValueError, OSError) as e:
            return self._err(str(e))
        out["refs_dir"] = str(APP / "refs")
        self._json(out)

    def _xrd_job(self, b, mode):
        if not X.SCRIPT.is_file():
            raise BadRequest(f"{X.SCRIPT} not found")
        try:
            args, key = X.build_args(b, mode, ST.state_dir)
        except ValueError as e:
            raise BadRequest(str(e))
        (ST.state_dir / "xrd").mkdir(exist_ok=True)
        cmd = [sys.executable, str(X.SCRIPT)] + args
        title = (f"XRD {'sticks' if mode == 'sticks' else 'files'} - {b.get('comp')} {'/'.join(b.get('phases') or [])}"
                 + (f" {float(b['lo']):.2f}..{float(b['hi']):.2f} meV/atom" if b.get("lo") is not None and b.get("hi") is not None else ""))
        watch = [str(ST.state_dir / "xrd")] if mode == "sticks" else [str(Path(args[args.index("--out-dir") + 1]))]
        jid = ST.jobs.start(cmd, ST.state_dir / "xrd", title=title, module="xrd_" + mode, watch=watch,
                            meta={"kind": "xrd", "mode": mode, "input": b.get("source")})
        return jid, key, cmd

    def api_post_xrd_compute(self, q):
        b = json.loads(self._body() or b"{}")
        old = b.get("cancel")                    # the page's previous computation, if still running
        if old and re.fullmatch(r"[\w-]{1,64}", str(old)):
            ST.jobs.cancel(old)
        jid, key, cmd = self._xrd_job(b, "sticks")
        cli = None
        try:
            ea, _ = X.build_args({**b, "out_dir": b.get("out_dir") or "/OUT_DIR"}, "export", ST.state_dir)
            cli = shlex.join(["python3", "xrd_average.py"] + ea)
        except ValueError:
            pass
        self._json({"ok": True, "job": jid, "key": key, "cmd": shlex.join(cmd), "cli": cli})

    def api_post_xrd_export(self, q):
        b = json.loads(self._body() or b"{}")
        od = Path(str(b.get("out_dir") or "")).expanduser()
        if od.is_absolute():
            od.mkdir(parents=True, exist_ok=True)
        jid, _, cmd = self._xrd_job(b, "export")
        self._json({"ok": True, "job": jid, "cmd": shlex.join(cmd)})

    def api_get_xrd_sticks(self, q):
        p = X.sticks_path(ST.state_dir, (q.get("key") or [""])[0])
        if not p or not p.is_file():
            return self._err("no such result (computation not finished?)", 404)
        data = p.read_bytes()
        if "gzip" in (self.headers.get("Accept-Encoding") or ""):
            import gzip
            return self._send(200, gzip.compress(data, 5), "application/json", {"Content-Encoding": "gzip"})
        self._send(200, data)


class Server6(ThreadingHTTPServer):
    import socket as _s
    address_family = _s.AF_INET6


def serve(host="127.0.0.1", port=8765, scripts=None, state_dir=None, token=None, start_dir=None, open_browser=False):
    global ST
    scripts = Path(scripts or APP / "scripts").resolve()
    state_dir = Path(state_dir or Path.home() / ".fluorite_studio").expanduser().resolve()
    if not state_dir.exists():
        state_dir.mkdir(parents=True, mode=0o700)
    elif state_dir.stat().st_uid != os.getuid():
        sys.exit(f"state folder {state_dir} belongs to another user: choose your own with --state")
    elif state_dir.stat().st_mode & 0o077:
        os.chmod(state_dir, 0o700)              # the token and job logs live here
    tok_file = state_dir / "token"
    if tok_file.is_symlink():
        sys.exit(f"{tok_file} is a symbolic link: remove it")
    if token is None:
        token = tok_file.read_text().strip() if tok_file.is_file() else secrets.token_urlsafe(18)
    old = os.umask(0o077)
    try:
        tok_file.write_text(token)
    finally:
        os.umask(old)
    os.chmod(tok_file, 0o600)
    ST = State(scripts, state_dir, token, Path(start_dir or os.getcwd()).resolve())
    httpd = None
    for p in range(port, port + 20):
        try:
            httpd = (Server6 if ":" in host else ThreadingHTTPServer)((host, p), H)
            port = p
            break
        except OSError:
            continue
    if httpd is None:
        sys.exit(f"no free port in {port}-{port + 19}")
    ST.port = port
    shown = "127.0.0.1" if host in ("127.0.0.1", "0.0.0.0", "localhost") else f"[{host}]" if ":" in host else host
    url = f"http://{shown}:{port}/?token={token}"           # 127.0.0.1, not localhost: that may resolve to ::1
    node = os.uname().nodename
    print(f"Fluorite Ordering Studio {R.APP_VERSION}\n  scripts : {scripts}\n  state   : {state_dir}\n"
          f"  open    : {url}\n"
          f"  from your laptop (cluster): ssh -N -L {port}:127.0.0.1:{port} {os.environ.get('USER', 'you')}@{node}"
          f"   then open the link above\n  stop    : Ctrl+C (running jobs keep running)", flush=True)
    if open_browser:
        import webbrowser
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
