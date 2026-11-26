"""
runner.py - run one command for the studio, detached from the web server, and record how it ended.

    python3 runner.py JOBDIR -- CMD ARG ...

Started in its own session (process group = this runner).  Writes JOBDIR/log.txt (stdout + stderr),
JOBDIR/pid ("<runner pid> <start ticks>" so a reused pid is not mistaken for the job), and when the command
ends JOBDIR/exit (exit code), JOBDIR/end (unix time) and JOBDIR/cancelled if it was stopped from the page.
Stop = SIGTERM to the whole process group: the command and its children get it, the runner records it.
If the command is still alive GRACE seconds later, the runner records the end and SIGKILLs the whole group.
"""
import os
import shlex
import signal
import subprocess
import sys
import time

GRACE = float(os.environ.get("FO_STUDIO_KILL_GRACE", "10"))


def start_ticks(pid):
    try:
        with open(f"/proc/{pid}/stat") as f:
            return f.read().rsplit(")", 1)[1].split()[19]
    except (OSError, IndexError):
        return "0"


def _put(jobdir, name, text):
    tmp = os.path.join(jobdir, f".{name}.tmp")
    with open(tmp, "w") as f:
        f.write(text)
    os.replace(tmp, os.path.join(jobdir, name))           # readers never see an empty file


def finish(jobdir, code, cancelled):
    if cancelled:
        _put(jobdir, "cancelled", "")
    _put(jobdir, "end", str(time.time()))
    _put(jobdir, "exit", str(code))                        # written last: its presence means "ended"


def main():
    cancelled = []
    signal.signal(signal.SIGTERM, lambda *a: cancelled.append(time.time()))   # first, before anything else
    signal.signal(signal.SIGHUP, signal.SIG_IGN)
    jobdir = sys.argv[1]
    cmd = sys.argv[sys.argv.index("--") + 1:]
    cwd = open(os.path.join(jobdir, "cwd")).read().strip() if os.path.isfile(os.path.join(jobdir, "cwd")) else None
    with open(os.path.join(jobdir, "pid"), "w") as f:
        f.write(f"{os.getpid()} {start_ticks(os.getpid())}")
    code = 127
    with open(os.path.join(jobdir, "log.txt"), "ab", buffering=0) as log:
        log.write(("$ " + shlex.join(cmd) + "\n").encode())
        try:
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
            while True:
                try:
                    code = p.wait(timeout=1)
                    break
                except (InterruptedError, subprocess.TimeoutExpired):
                    pass
                if cancelled and time.time() - cancelled[0] > GRACE:     # ignored SIGTERM: kill the group
                    log.write(f"\n[stopped from the page; still running after {GRACE:g} s, killed]\n".encode())
                    finish(jobdir, -9, True)
                    os.killpg(os.getpgrp(), signal.SIGKILL)
        except OSError as e:
            log.write(f"cannot start: {e}\n".encode())
        if cancelled:
            log.write(b"\n[stopped from the page]\n")
    # a job that finished by itself with exit 0 just as Stop was pressed is "done", not "cancelled"
    finish(jobdir, code, bool(cancelled) and code != 0)


if __name__ == "__main__":
    main()
