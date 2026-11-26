#!/usr/bin/env python3
"""
map_check.py - map or transform ONE structure onto other fluorite-derived phases and verify the cation ordering.

    python3 map_check.py POSCAR --method map --templates t=refs/POSCAR_t oIII=refs/POSCAR_oIII m=refs/POSCAR_m
    python3 map_check.py POSCAR --method transform --refs refs --targets t oIII m --polar-axis a --unique-axis a

  --method map        map_fluorite.py: direct smallest-distortion correspondence onto each template
  --method transform  transform_set.py: physical path through the tetragonal hub
                      (t main axis -> longest oIII / m axis; polar / unique axis in-plane)

Every output is checked with compare_ordering.py against the input (must be IDENTICAL).  Prints one line per
output (rms shift, strain S, equally good variants, verdict) and writes <out-dir>/map_check.csv.
Only calls the scripts in the scripts folder; changes none of them.
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SCRIPTS = os.path.normpath(os.path.join(HERE, "..", "..", "scripts"))
DEFAULT_REFS = os.path.normpath(os.path.join(HERE, "..", "..", "refs"))


def verdict(scripts, a, b):
    r = subprocess.run([sys.executable, os.path.join(scripts, "compare_ordering.py"), a, b],
                       capture_output=True, text=True)
    m = re.search(r"VERDICT: (\w+)", r.stdout)
    return m.group(1) if m else "ERROR " + (r.stderr.strip().split("\n") or [""])[-1][-200:]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("poscar", help="input structure (VASP5, diagonal repeat of the 12-atom cell)")
    ap.add_argument("--method", choices=["map", "transform"], default="map")
    ap.add_argument("--templates", nargs="+", metavar="LABEL=POSCAR", default=[],
                    help="map: target templates, e.g. t=refs/POSCAR_t (default: every POSCAR_<label> in --refs)")
    ap.add_argument("--refs", default=DEFAULT_REFS,
                    help="folder with POSCAR_c, POSCAR_t, POSCAR_oIII, POSCAR_m (default: the app's refs/); "
                         "map uses all of them when --templates is not given")
    ap.add_argument("--targets", nargs="+", default=["t", "oIII", "m"], help="transform: target phases")
    ap.add_argument("--main-axis", choices=list("abc"), help="transform: t main axis when several tie")
    ap.add_argument("--polar-axis", choices=list("abc"), help="transform: oIII polar axis")
    ap.add_argument("--unique-axis", choices=list("abc"), help="transform: m unique axis")
    ap.add_argument("--all-variants", action="store_true", help="write every equally good variant")
    ap.add_argument("--out-dir", help="output folder (default: $FO_STUDIO_OUT/map_check/<input name>, else "
                                      "./map_check_<input name>)")
    ap.add_argument("--scripts", default=DEFAULT_SCRIPTS, help="scripts folder (default: the app's scripts/)")
    a = ap.parse_args()
    if a.method == "map" and not a.templates:
        a.templates = sorted(glob.glob(os.path.join(a.refs, "POSCAR_*")))
        if not a.templates:
            sys.exit(f"--method map: no --templates given and no POSCAR_<label> files in {a.refs}")
    if not a.out_dir:
        base = os.environ.get("FO_STUDIO_OUT")
        name = os.path.basename(a.poscar) + "_" + a.method
        a.out_dir = os.path.join(base, "map_check", name) if base else "map_check_" + name
    os.makedirs(a.out_dir, exist_ok=True)
    print(f"output folder: {os.path.abspath(a.out_dir)}")
    stem = os.path.basename(a.poscar)
    rows = []
    if a.method == "map":
        for spec in a.templates:
            label, tpl = spec.split("=", 1) if re.match(r"[\w.+\-]+=", spec) else \
                (re.sub(r"^POSCAR[_.\-]?", "", os.path.basename(spec)), spec)
            out = os.path.join(a.out_dir, f"POSCAR_{label}_{stem}")
            cmd = [sys.executable, os.path.join(a.scripts, "map_fluorite.py"), a.poscar, "--template", tpl,
                   "--out", out, "--label", label] + (["--all-variants"] if a.all_variants else [])
            r = subprocess.run(cmd, capture_output=True, text=True)
            print(f"$ {' '.join(cmd)}\n{r.stdout}{r.stderr}")
            row = {"target": label, "method": "map_fluorite", "output": out, "rms_shift_A": "", "strain_S": "",
                   "equally_good_variants": "1", "axes": "", "verdict": ""}
            m = re.search(r"rms atomic shift ([\d.]+) A.*?S = ([\d.]+)", r.stdout)
            if m:
                row["rms_shift_A"], row["strain_S"] = m.groups()
            m = re.search(r"variants: (\d+) structures need exactly the same distortion", r.stdout)
            if m:
                row["equally_good_variants"] = m.group(1)
            m = re.search(r"axes\s*: (.*)", r.stdout)
            if m:
                row["axes"] = m.group(1).strip()
            if r.returncode != 0:
                row["verdict"] = "MAP FAILED"
                rows.append(row)
                continue
            outs = sorted(glob.glob(out + "_v[0-9]*"), key=lambda f: int(f.rsplit("_v", 1)[1])) or [out]
            for o in outs:
                rows.append(dict(row, output=o, verdict=verdict(a.scripts, a.poscar, o)))
    else:
        cmd = [sys.executable, os.path.join(a.scripts, "transform_set.py"), a.poscar, "--targets", *a.targets,
               "--out-dir", a.out_dir]
        if a.refs:
            cmd += ["--refs", a.refs]
        for k in ("main_axis", "polar_axis", "unique_axis"):
            if getattr(a, k):
                cmd += ["--" + k.replace("_", "-"), getattr(a, k)]
        if a.all_variants:
            cmd.append("--all-variants")
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(f"$ {' '.join(cmd)}\n{r.stdout}{r.stderr}")
        if not os.path.isfile(os.path.join(a.out_dir, "summary.csv")):
            sys.exit(f"transform_set.py wrote no summary.csv (exit code {r.returncode}); see its output above")
        with open(os.path.join(a.out_dir, "summary.csv")) as f:
            for s in csv.DictReader(f):
                o = s.get("output", "")
                if o and not os.path.isabs(o) and not os.path.isfile(o):
                    o = os.path.join(a.out_dir, o)
                rows.append({"target": s["target"], "method": "transform_set", "output": o,
                             "rms_shift_A": s.get("step_rms_shift_A", ""), "strain_S": s.get("step_strain_S", ""),
                             "equally_good_variants": "", "axes": f"{s.get('main_axis', '')} {s.get('detail', '')}".strip(),
                             "verdict": verdict(a.scripts, a.poscar, o) if s.get("status") == "ok" and os.path.isfile(o)
                             else "FAILED " + s.get("note", "")})
    print(f"\n{'target':<6} {'rms shift A':>11} {'strain S':>9} {'variants':>8}  {'verdict':<10} output")
    for r in rows:
        print(f"{r['target']:<6} {r['rms_shift_A']:>11} {r['strain_S']:>9} {r['equally_good_variants']:>8}  "
              f"{r['verdict']:<10} {os.path.basename(r['output'])}")
    if any(r["equally_good_variants"] not in ("", "1") for r in rows):
        print("\nEqually good variants cannot be ranked by distortion, only by energy (use --all-variants and relax them all).")
    with open(os.path.join(a.out_dir, "map_check.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["target"])
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if r["verdict"] != "IDENTICAL"]
    print(f"\n{len(rows) - len(bad)} of {len(rows)} outputs IDENTICAL to the input ordering")
    sys.exit(1 if bad or not rows else 0)


if __name__ == "__main__":
    main()
