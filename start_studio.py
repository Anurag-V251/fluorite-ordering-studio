#!/usr/bin/env python3
"""
Fluorite Ordering Studio - browser interface for the cation-ordering / NEP free-energy workflow.

    python3 start_studio.py                       # http://localhost:8765/?token=...  (printed)
    python3 start_studio.py --port 8800 --dir /path/to/results

On a cluster login node: start it there, then on your laptop  ssh -N -L 8765:localhost:8765 user@login-node
and open the printed link.  Only Python's standard library is needed for the app itself; the scripts it runs
need numpy / scipy / matplotlib (plots) and ase / calorine / phonopy (free energies), as before.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="127.0.0.1", help="interface (default 127.0.0.1 = this machine only)")
    ap.add_argument("--dir", help="folder the file pickers start in (default: current folder)")
    ap.add_argument("--scripts", help="scripts folder (default: scripts/ next to studio.py)")
    ap.add_argument("--state", help="where jobs/settings are kept (default ~/.fluorite_studio)")
    ap.add_argument("--new-token", action="store_true", help="make a new access token")
    ap.add_argument("--open", action="store_true", help="open a browser on this machine")
    a = ap.parse_args()
    from studio import server
    tok = None
    if a.new_token:
        import secrets
        tok = secrets.token_urlsafe(18)                 # written (mode 600) by the server
    if a.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"WARNING: listening on {a.host}: anyone who can reach this port and knows the token can run "
              "commands as you. Prefer the default 127.0.0.1 with ssh -L.", file=sys.stderr)
    server.serve(a.host, a.port, a.scripts, a.state, tok, a.dir, a.open)


if __name__ == "__main__":
    main()
