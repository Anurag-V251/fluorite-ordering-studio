"""
jobs.py - background jobs started from the web page (local processes) and a view of SLURM jobs.

Every job lives in STATE/jobs/<id>/: job.json (command, cwd, title, module, start time, watched folders),
log.txt, pid, exit, end.  Jobs survive a restart of the server (runner.py runs them in their own session).
"""
from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
SLURM_ID_RES = [r"submitted (?:merge )?job (\d+)", r"Submitted batch job (\d+)", r"submitted \w+ array (\d+)"]


def _wjson(path, obj):
    """atomic write (readers in other threads never see a half-written file)"""
    tmp = Path(str(path) + f".tmp{os.getpid()}-{uuid.uuid4().hex[:6]}")
    tmp.write_text(json.dumps(obj, indent=1), encoding="utf-8")
    os.replace(tmp, path)


def _rjson(path, default=None):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default
OUT_EXT = {".png", ".pdf", ".svg", ".html", ".csv", ".txt", ".eps", ".sh", ".json", ".toml"}


class Jobs:
    def __init__(self, state: Path):
        self.dir = Path(state) / "jobs"
        self.dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ start / stop
    def start(self, cmd, cwd, title="", module=None, watch=(), meta=None, env=None):
        jid = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
        jd = self.dir / jid
        jd.mkdir()
        cmd = [str(c) for c in cmd]
        info = {"id": jid, "cmd": cmd, "cwd": str(cwd), "title": title or Path(cmd[1] if len(cmd) > 1 else cmd[0]).name,
                "module": module, "start": time.time(), "watch": [str(w) for w in watch], "meta": meta or {}}
        _wjson(jd / "snapshot.json", _snapshot(info["watch"]))
        (jd / "cwd").write_text(str(cwd), encoding="utf-8")
        _wjson(jd / "job.json", info)
        env = dict(os.environ, PYTHONUNBUFFERED="1", **(env or {}))
        with open(jd / "runner_errors.txt", "wb") as err:
            subprocess.Popen([sys.executable, str(HERE / "runner.py"), str(jd.resolve()), "--", *cmd],
                             stdin=subprocess.DEVNULL, stdout=err, stderr=err,
                             start_new_session=True, env=env, cwd=str(cwd))
        return jid

    def cancel(self, jid):
        st = self.get(jid)
        if not st or st["state"] != "running" or not st.get("pid"):
            return False
        try:                                   # the runner leads its own process group (start_new_session)
            os.killpg(int(st["pid"]), signal.SIGTERM)
        except OSError:
            return False
        return True

    def delete(self, jid):
        jd = self._jd(jid)
        if jd and self.get(jid)["state"] != "running":
            shutil.rmtree(jd)
            return True
        return False

    # ------------------------------------------------------------------ read
    def _jd(self, jid):
        if not re.fullmatch(r"[\w\-]+", jid or ""):
            return None
        jd = self.dir / jid
        return jd if (jd / "job.json").is_file() else None

    def get(self, jid, log_tail=None, with_outputs=True):
        jd = self._jd(jid)
        if not jd:
            return None
        info = _rjson(jd / "job.json")
        if not isinstance(info, dict) or "start" not in info:
            return None
        rd = lambda n: (jd / n).read_text().strip() if (jd / n).is_file() else None
        pidtxt, ex, end = rd("pid"), rd("exit"), rd("end")
        if ex is not None and not re.fullmatch(r"-?\d+", ex):
            ex = None                          # caught mid-write (older runner): not ended yet
        pid, ticks = (pidtxt.split() + ["0"])[:2] if pidtxt else (None, None)
        if ex is not None:
            state = "cancelled" if (jd / "cancelled").is_file() else "done" if ex == "0" else "failed"
        elif pid and _alive(int(pid), ticks):
            state = "running"
        elif time.time() - info["start"] < 60:
            state = "starting"
        else:
            state = "lost"                      # process vanished without writing its exit code (machine reboot)
        info.update(state=state, pid=pid, exit=None if ex is None else int(ex),
                    end=float(end) if end else None)
        log = jd / "log.txt"
        if log_tail is not None and log.is_file():
            size = log.stat().st_size
            with log.open("rb") as f:
                f.seek(max(0, size - log_tail))
                info["log"] = f.read().decode(errors="replace")
            info["log_size"] = size
            info["log_path"] = str(log)
        ended = state in ("done", "failed", "cancelled")
        cache = jd / "slurm_ids.json"
        ids = _rjson(cache) if ended else None
        if not isinstance(ids, list):
            txt = ""
            if log.is_file():                  # a huge log: its first and last 2 MB hold the submit lines
                size = log.stat().st_size
                with log.open("rb") as f:
                    txt = f.read(2_000_000).decode(errors="replace")
                    if size > 4_000_000:
                        f.seek(size - 2_000_000)
                    txt += f.read().decode(errors="replace")
            ids = set()
            for rx in SLURM_ID_RES:
                ids |= set(re.findall(rx, txt))
            ids = sorted(ids, key=int)
            if ended:
                _wjson(cache, ids)
        info["slurm_ids"] = ids
        if ended and with_outputs:
            outs = _rjson(jd / "outputs.json")  # frozen at the first look after the end, so a later job that
            if not isinstance(outs, list):     # touches the same folder does not add its files here
                outs = self.outputs(info)
                _wjson(jd / "outputs.json", outs)
            info["outputs"] = outs
        return info

    def outputs(self, info):
        """files in the watched folders that are new or changed since the job started"""
        try:
            before = _rjson(self.dir / info["id"] / "snapshot.json", {})
        except Exception:
            before = {}
        t0 = info["start"] - 1
        out = []
        for w in info.get("watch", []):
            rec = w.startswith("R:")
            p = Path(w[2:] if rec else w)
            if not p.is_dir():
                continue
            try:
                for e in _files(p, 3 if rec else 1):
                    if not (Path(e.name).suffix.lower() in OUT_EXT or (rec and e.name.startswith("POSCAR"))):
                        continue
                    mt = e.stat().st_mtime
                    if (e.path in before and mt != before[e.path]) or (e.path not in before and mt >= t0):
                        out.append({"path": e.path, "name": os.path.relpath(e.path, p), "size": e.stat().st_size,
                                    "mtime": e.stat().st_mtime})
            except OSError:
                pass
        seen, uniq = set(), []
        for o in sorted(out, key=lambda x: x["name"]):
            if o["path"] not in seen:
                seen.add(o["path"])
                uniq.append(o)
        return uniq[:500]

    def running(self, module=None, config=None):
        out = []
        ids = sorted((p.name for p in self.dir.iterdir()
                      if (p / "job.json").is_file() and not (p / "exit").is_file()), reverse=True)   # all, not only recent
        for jid in ids:
            j = self.get(jid, with_outputs=False)
            if j and j["state"] in ("running", "starting") and (module is None or j.get("module") == module) and \
                    (config is None or (j.get("meta") or {}).get("config") == config):
                out.append(j)
        return out

    def list(self, limit=200):
        ids = sorted((p.name for p in self.dir.iterdir() if (p / "job.json").is_file()), reverse=True)[:limit]
        out = []
        for j in ids:
            g = self.get(j, with_outputs=False)
            if g:
                out.append(g)
        return out


def _files(d, depth):
    """files below d, at most `depth` levels deep (1 = only d itself), at most 20000"""
    out, stack = [], [(str(d), 1)]
    while stack and len(out) < 20000:
        cur, lv = stack.pop()
        try:
            for e in os.scandir(cur):
                if e.is_file():
                    out.append(e)
                elif e.is_dir() and lv < depth:
                    stack.append((e.path, lv + 1))
        except OSError:
            pass
    return out


def _snapshot(dirs):
    snap = {}
    for w in dirs:
        rec = w.startswith("R:")
        for e in _files(w[2:] if rec else w, 3 if rec else 1):
            try:
                snap[e.path] = e.stat().st_mtime
            except OSError:
                pass
    return snap


def _alive(pid, ticks=None):
    if ticks and ticks != "0":                 # same process as recorded (pid not reused)?
        try:
            with open(f"/proc/{pid}/stat") as f:
                st = f.read().rsplit(")", 1)[1].split()
            if st[19] != ticks:
                return False
            return st[0] != "Z"
        except (OSError, IndexError):
            pass
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # a zombie child of this server counts as finished
    try:
        with open(f"/proc/{pid}/stat") as f:
            return f.read().split(")")[-1].split()[0] != "Z"
    except OSError:
        return True


# ---------------------------------------------------------------------- SLURM
def username():
    import getpass
    try:
        return getpass.getuser()
    except Exception:                          # no passwd entry (containers)
        return str(os.getuid())


def slurm_available():
    return shutil.which("squeue") is not None


def squeue():
    if not slurm_available():
        return {"available": False, "jobs": []}
    user = username()
    fmt = "%i|%T|%M|%l|%D|%P|%R|%j"            # name last: it may contain "|"
    try:
        r = subprocess.run(["squeue", "-h", "-u", user, "-o", fmt], capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as e:
        return {"available": True, "error": str(e), "jobs": []}
    jobs = []
    for line in r.stdout.splitlines():
        p = line.split("|", 7)
        if len(p) == 8:
            jobs.append(dict(zip(["id", "state", "time", "limit", "nodes", "partition", "reason", "name"], p)))
    return {"available": True, "jobs": jobs, "error": r.stderr.strip() or None}


def sacct(ids):
    if not ids or shutil.which("sacct") is None:
        return []
    try:
        r = subprocess.run(["sacct", "-n", "-P", "-X", "-j", ",".join(ids), "-o", "JobID,JobName,State,Elapsed,ExitCode,Timelimit,NodeList"],
                           capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return []
    out = []
    for line in r.stdout.splitlines():
        p = line.split("|")
        if len(p) >= 5:
            out.append(dict(zip(["id", "name", "state", "elapsed", "exit", "limit", "nodes"], p)))
    return out


def scancel(ids):
    if shutil.which("scancel") is None:
        return "scancel not found"
    r = subprocess.run(["scancel", *ids], capture_output=True, text=True)
    return (r.stdout + r.stderr).strip() or "ok"
