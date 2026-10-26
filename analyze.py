#!/usr/bin/env python3
"""Analysis of the H2 bond-length scans (lda.dat, pbe.dat).

- plots E(b) for each functional (lda_curve.png, pbe_curve.png, both_curves.png)
- fits E(b) = alpha + beta*b + gamma*b^2 to the points within +/- FIT_HALF_WIDTH
  of the grid minimum and reports alpha, beta, gamma, b_eq = -beta/(2 gamma)
- compares b_eq with experiment (0.74 Angstrom) and writes results.md
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B_EXP = 0.74            # Angstrom, experimental H2 bond length (value given on the problem sheet)
FIT_HALF_WIDTH = 0.10   # Angstrom: fit window around the grid minimum (+/-)


def load(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


def fit(b, e, half_width=FIT_HALF_WIDTH):
    imin = int(np.argmin(e))
    bmin = b[imin]
    mask = np.abs(b - bmin) <= half_width + 1e-9
    # numpy.polyfit returns highest power first: [gamma, beta, alpha]
    gamma, beta, alpha = np.polyfit(b[mask], e[mask], 2)
    beq = -beta / (2.0 * gamma)
    emin_fit = alpha + beta * beq + gamma * beq**2
    resid = e[mask] - np.polyval([gamma, beta, alpha], b[mask])
    return dict(alpha=alpha, beta=beta, gamma=gamma, beq=beq, emin_fit=emin_fit,
                bmin_grid=bmin, emin_grid=e[imin], mask=mask,
                rms_resid=float(np.sqrt(np.mean(resid**2))))


def plot_single(b, e, r, label, fname):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(b, e, "o-", label=f"{label} pw.x data")
    bb = np.linspace(b[r["mask"]].min(), b[r["mask"]].max(), 200)
    ax.plot(bb, np.polyval([r["gamma"], r["beta"], r["alpha"]], bb), "--",
            label=f"quadratic fit (|b - {r['bmin_grid']:.2f}| <= {FIT_HALF_WIDTH:.2f} Å)")
    ax.axvline(r["beq"], color="k", ls=":", lw=1, label=f"$b_{{eq}}$ = {r['beq']:.4f} Å")
    ax.axvline(B_EXP, color="gray", ls="-.", lw=1, label=f"experiment {B_EXP:.2f} Å")
    ax.set_xlabel("H–H bond length $b$ (Å)")
    ax.set_ylabel("Total energy $E$ (Ry)")
    ax.set_title(f"H$_2$ total energy vs bond length ({label})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


def plot_both(data):
    """Single figure for the report: left, both E(b) curves; right, the 5-point
    fit window for each functional with the fitted parabola."""
    fig = plt.figure(figsize=(8.4, 4.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], hspace=0.45, wspace=0.3)
    ax = fig.add_subplot(gs[:, 0])
    colors = {"LDA": "C0", "PBE": "C1"}
    for label, (b, e, r) in data.items():
        c = colors[label[:3]]
        ax.plot(b, e, "o-", color=c, ms=4, label=f"{label}")
    ax.axvline(B_EXP, color="gray", ls="-.", lw=1, label=f"experiment {B_EXP:.2f} Å")
    ax.set_xlabel("H–H bond length $b$ (Å)")
    ax.set_ylabel("Total energy $E$ (Ry)")
    ax.set_title("(a) full scan, 0.60–1.20 Å", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    for i, (label, (b, e, r)) in enumerate(data.items()):
        c = colors[label[:3]]
        pax = fig.add_subplot(gs[i, 1])
        m = r["mask"]
        bb = np.linspace(b[m].min(), b[m].max(), 200)
        pax.plot(b[m], 1e3 * (e[m] - r["emin_fit"]), "o", color=c, ms=4, label="pw.x")
        pax.plot(bb, 1e3 * (np.polyval([r["gamma"], r["beta"], r["alpha"]], bb) - r["emin_fit"]), "--", color=c, lw=1,
                 label=r"$\alpha+\beta b+\gamma b^2$")
        pax.axvline(r["beq"], color="k", ls=":", lw=0.8)
        pax.axvline(B_EXP, color="gray", ls="-.", lw=0.8)
        pax.set_title(f"({'bc'[i]}) {label[:3]} fit: $b_{{eq}}$ = {r['beq']:.3f} Å", fontsize=9)
        pax.set_ylabel("$E-E_{min}$ (mRy)", fontsize=8)
        pax.tick_params(labelsize=7)
        pax.grid(alpha=0.3)
        pax.legend(fontsize=7, loc="upper center")
        if i == 1:
            pax.set_xlabel("$b$ (Å)", fontsize=8)
    fig.savefig("both_curves.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    data = {}
    for key, label in (("lda", "LDA (H.pz-vbc)"), ("pbe", "PBE (H.pbe-vbc)")):
        b, e = load(f"{key}.dat")
        r = fit(b, e)
        data[label] = (b, e, r)
        plot_single(b, e, r, label, f"{key}_curve.png")
        print(f"=== {label} ===")
        print(f"grid minimum: b = {r['bmin_grid']:.2f} Å, E = {r['emin_grid']:.8f} Ry")
        print(f"fit points  : b = {', '.join(f'{x:.2f}' for x in b[r['mask']])} Å ({r['mask'].sum()} points)")
        print(f"alpha = {r['alpha']:.6f} Ry   beta = {r['beta']:.6f} Ry/Å   gamma = {r['gamma']:.6f} Ry/Å^2")
        print(f"b_eq  = {r['beq']:.4f} Å   E(b_eq) = {r['emin_fit']:.6f} Ry   rms fit residual = {r['rms_resid']:.2e} Ry")
        print(f"error vs experiment ({B_EXP} Å): {100*(r['beq']-B_EXP)/B_EXP:+.2f} %")
        print()
    plot_both(data)

    print("| Method | b_eq (Å) | Error vs exp. (%) |")
    print("|---|---|---|")
    for label, (b, e, r) in data.items():
        print(f"| {label} | {r['beq']:.4f} | {100*(r['beq']-B_EXP)/B_EXP:+.2f} |")
    print(f"| Experiment | {B_EXP:.4f} | — |")
    # grid-induced uncertainty: spread of b_eq between 3- and 7-point windows
    for label, (b, e, r) in data.items():
        spread = [fit(b, e, hw)['beq'] for hw in (0.06, 0.10, 0.14)]
        print(f"{label}: b_eq window spread {min(spread):.4f}-{max(spread):.4f} Å -> +/- {(max(spread)-min(spread))/2:.4f} Å")


if __name__ == "__main__":
    main()
