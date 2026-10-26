#!/usr/bin/env python3
"""Convergence tests for the numerical settings used in the H2 scan.

Re-runs the five bond lengths used in the quadratic fit (0.68-0.84 A) for both
functionals with (i) a higher plane-wave cutoff and (ii) a larger box, and
refits b_eq, so the influence of the fixed numerical settings on b_eq can be
quantified instead of guessed.  Runs are independent and executed in parallel in
separate directories under conv/.  Writes conv/summary.dat and conv/conv.json.
"""
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_scan import XC, make_input, parse_output  # noqa: E402

os.environ.setdefault("HWLOC_COMPONENTS", "-x86")
BONDS = [0.68, 0.72, 0.76, 0.80, 0.84]
SETTINGS = {              # label: (ecutwfc / Ry, A / Angstrom)
    "ref_60Ry_12A": (60.0, 12.0),
    "ecut_100Ry_12A": (100.0, 12.0),
    "box_60Ry_16A": (60.0, 16.0),
    "ecut_100Ry_16A": (100.0, 16.0),
}


def modify(inp, ecut, box):
    inp = re.sub(r"ecutwfc\s*=\s*[\d.]+", f"ecutwfc = {ecut:.1f}", inp)
    inp = re.sub(r"\bA\s*=\s*[\d.]+", f"A = {box:.1f}", inp)
    return inp


def run_one(job):
    xc, label, b = job
    ecut, box = SETTINGS[label]
    d = os.path.join("conv", f"{xc}_{label}_b{b:.2f}")
    os.makedirs(d, exist_ok=True)
    for f in XC.values():
        if not os.path.exists(os.path.join(d, f)):
            shutil.copy(f, d)
    inp = modify(make_input(open("H2.in").read(), b, XC[xc]), ecut, box)
    open(os.path.join(d, "pw.in"), "w").write(inp)
    with open(os.path.join(d, "pw.out"), "w") as fo:
        subprocess.run(["pw.x", "-in", "pw.in"], cwd=d, stdout=fo, stderr=subprocess.STDOUT, check=True)
    conv, e, xcline = parse_output(os.path.join(d, "pw.out"))
    if not conv or e is None:
        sys.exit(f"ERROR: {d} did not converge")
    shutil.rmtree(os.path.join(d, "pwscf.save"), ignore_errors=True)
    print(f"{xc:3s} {label:15s} b={b:.2f}  E={e:.8f} Ry  [{xcline}]", flush=True)
    return (xc, label, b, e)


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    jobs = [(xc, lab, b) for xc in XC for lab in SETTINGS for b in BONDS]
    with ThreadPoolExecutor(nproc) as ex:
        res = list(ex.map(run_one, jobs))
    table = {}
    for xc, lab, b, e in res:
        table.setdefault(xc, {}).setdefault(lab, {})[b] = e
    out = {}
    lines = ["# xc  setting  ecut(Ry)  A(A)  E(0.76)(Ry)  b_eq(A)  gamma(Ry/A^2)  E_min(Ry)"]
    for xc in table:
        for lab in SETTINGS:
            bs = np.array(BONDS)
            es = np.array([table[xc][lab][b] for b in BONDS])
            g, be, al = np.polyfit(bs, es, 2)
            beq = -be / (2 * g)
            emin = al + be * beq + g * beq**2
            out[f"{xc}/{lab}"] = dict(ecut=SETTINGS[lab][0], A=SETTINGS[lab][1], E=dict(zip(map(str, BONDS), es.tolist())),
                                     beq=float(beq), gamma=float(g), emin=float(emin))
            lines.append(f"{xc:3s} {lab:15s} {SETTINGS[lab][0]:6.1f} {SETTINGS[lab][1]:5.1f} "
                         f"{table[xc][lab][0.76]:.8f} {beq:.4f} {g:.4f} {emin:.6f}")
    open("conv/summary.dat", "w").write("\n".join(lines) + "\n")
    json.dump(out, open("conv/conv.json", "w"), indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
