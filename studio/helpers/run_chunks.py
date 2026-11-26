#!/usr/bin/env python3
"""
run_chunks.py - run the unfinished chunks of a campaign one after the other on this machine (no SLURM).

    python3 run_chunks.py campaign.toml [--chunks 0 1 5]

Each chunk is `campaign.py CONFIG run --chunk K` (which uses all cores of this machine, cores_per_node /
threads structures in parallel; SLURM_CPUS_ON_NODE overrides).  Finished chunks are skipped, so the command can
be started again after an interruption.  Run `merge` afterwards.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_SCRIPTS = HERE.parent.parent / "scripts"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("config")
    ap.add_argument("--chunks", type=int, nargs="+", help="only these chunk ids")
    ap.add_argument("--scripts", default=None, help="scripts folder (default: [paths] scripts of the config)")
    a = ap.parse_args()
    scripts = a.scripts
    if scripts is None:                                          # the run's own (frozen) scripts, as campaign.py uses
        try:
            import tomllib
        except ModuleNotFoundError:
            import tomli as tomllib                              # type: ignore
        with open(a.config, "rb") as f:
            sp = tomllib.load(f).get("paths", {}).get("scripts")
        scripts = (Path(a.config).resolve().parent / sp) if sp else DEFAULT_SCRIPTS
    sys.path.insert(0, str(scripts))
    import campaign as C                                         # noqa: E402
    cfg = C.Cfg(a.config)
    plan = C.load_plan(cfg)
    prog = C.chunk_progress(cfg, plan)
    todo = [p[0] for p in prog if p[2] < p[1] and (a.chunks is None or p[0] in a.chunks)]
    print(f"{len(prog)} chunks, {len(todo)} to run: {C.compress(todo) if todo else '-'}", flush=True)
    fails = 0
    for i, k in enumerate(todo, 1):
        print(f"\n=== chunk {k} ({i}/{len(todo)})", flush=True)
        r = subprocess.call([sys.executable, str(cfg.scripts / "campaign.py"), str(cfg.path), "run", "--chunk", str(k)],
                            env=dict(os.environ, PYTHONUNBUFFERED="1"))
        fails += r != 0
    print(f"\ndone: {len(todo) - fails} chunks ok, {fails} with errors. Next: merge.")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
