#!/usr/bin/env python3
"""Bond-length scan of H2 with Quantum ESPRESSO pw.x.

Usage: python3 run_scan.py lda|pbe
Generates one pw.x input per bond length from the template H2.in (only the
z-coordinate of the second H atom and, for PBE, the pseudopotential filename
are changed), runs pw.x, checks SCF convergence and writes <xc>.dat (b, E).
"""
import os
import re
import subprocess
import sys

XC = {"lda": "H.pz-vbc.UPF", "pbe": "H.pbe-vbc.UPF"}
BONDS = [round(0.60 + 0.04 * i, 2) for i in range(16)]   # 0.60 ... 1.20 Å

os.environ.setdefault("HWLOC_COMPONENTS", "-x86")        # OpenMPI/hwloc crash workaround on this VM


def make_input(template, b, upf):
    lines = template.splitlines()
    out = []
    in_pos = False
    n_atom = 0
    for line in lines:
        if line.strip().startswith("ATOMIC_POSITIONS"):
            in_pos = True
            out.append(line)
            continue
        if in_pos and line.strip():
            n_atom += 1
            if n_atom == 2:
                line = f" H 0.00 0.00 {b:.2f}"
            elif n_atom > 2:
                in_pos = False
        elif in_pos:
            in_pos = False
        if "H.pz-vbc.UPF" in line:
            line = line.replace("H.pz-vbc.UPF", upf)
        out.append(line)
    return "\n".join(out) + "\n"


def parse_output(path):
    txt = open(path).read()
    conv = "convergence has been achieved" in txt
    m = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", txt, re.M)
    xc = re.search(r"Exchange-correlation=.*", txt)
    return conv, (float(m[-1]) if m else None), (xc.group(0).strip() if xc else None)


def main():
    xc = sys.argv[1].lower()
    upf = XC[xc]
    template = open("H2.in").read()
    os.makedirs("runs", exist_ok=True)
    results = []
    for b in BONDS:
        tag = f"{xc}_b{b:.2f}"
        inp = os.path.join("runs", tag + ".in")
        outp = os.path.join("runs", tag + ".out")
        open(inp, "w").write(make_input(template, b, upf))
        with open(outp, "w") as fo:
            subprocess.run(["pw.x", "-in", inp], stdout=fo, stderr=subprocess.STDOUT, check=True)
        conv, e, xcline = parse_output(outp)
        if not conv or e is None:
            sys.exit(f"ERROR: run {tag} did not converge / no total energy found")
        print(f"{tag}: E = {e:.8f} Ry   converged   [{xcline}]")
        results.append((b, e))
    with open(f"{xc}.dat", "w") as f:
        f.write(f"# H2 bond-length scan, pseudopotential {upf}\n# b (Angstrom)   E_total (Ry)\n")
        for b, e in results:
            f.write(f"{b:.2f}  {e:.8f}\n")
    print(f"wrote {xc}.dat")


if __name__ == "__main__":
    main()
