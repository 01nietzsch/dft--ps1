# PSet 1 — H₂ bond length with Quantum ESPRESSO (LDA vs PBE)

Setup: `pw.x` 6.7 (Ubuntu `quantum-espresso` package), `H2.in` unchanged except for the
z-coordinate of the second H atom (part a/c) and the pseudopotential filename (part c).
`ecutwfc = 60 Ry`, cubic box `A = 12 Å`, Γ-point only (`1 1 1 0 0 0`) — not modified.
All 32 SCF runs converged (`convergence has been achieved` found in every output; the script
aborts otherwise). Raw inputs/outputs are in `runs/`.

Functional actually used, quoted from the `pw.x` output:

| Pseudopotential | `pw.x` output line |
|---|---|
| `H.pz-vbc.UPF` | `Exchange-correlation=  SLA  PZ   NOGX NOGC` (Slater exchange + Perdew–Zunger correlation, no gradient corrections → LDA) |
| `H.pbe-vbc.UPF` | `Exchange-correlation= SLA  PW   PBE  PBE` (Slater + Perdew–Wang LDA part with PBE gradient corrections for exchange and correlation → GGA-PBE) |

Input-file sanity check (requested before running the loop), e.g. `runs/lda_b0.60.in` lines 21–22:
```
 H 0.00 0.00 0.00
 H 0.00 0.00 0.60
```

## (a) LDA scan — `lda.dat`

| b (Å) | E (Ry) | b (Å) | E (Ry) |
|---|---|---|---|
| 0.60 | −2.22352148 | 0.92 | −2.25120042 |
| 0.64 | −2.24671146 | 0.96 | −2.24128285 |
| 0.68 | −2.26131788 | 1.00 | −2.23030085 |
| 0.72 | −2.26924915 | 1.04 | −2.21851510 |
| **0.76** | **−2.27200515** | 1.08 | −2.20613839 |
| 0.80 | −2.27075872 | 1.12 | −2.19334916 |
| 0.84 | −2.26642664 | 1.16 | −2.18030016 |
| 0.88 | −2.25972301 | 1.20 | −2.16712226 |

Lowest grid point: **b = 0.76 Å, E = −2.27200515 Ry**. Plot: `lda_curve.png`.

## (b) Quadratic fit near the minimum (LDA)

Fit E(b) = α + βb + γb² to the five points with |b − 0.76| ≤ 0.10 Å (b = 0.68, 0.72, 0.76, 0.80, 0.84 Å),
using `numpy.polyfit(b, E, 2)` (returns `[γ, β, α]`, highest power first).

| α (Ry) | β (Ry/Å) | γ (Ry/Å²) | b_eq = −β/2γ (Å) | E(b_eq) (Ry) |
|---|---|---|---|---|
| −1.514102 | −1.965223 | 1.273622 | **0.7715** | −2.272196 |

Choice of fit range: the curve is strongly anharmonic (steep repulsive wall on the left, soft
dissociation tail on the right), so a parabola is only a good model very close to the minimum.
±0.10 Å gives 5 points — enough to over-determine the 3 parameters and average out noise, while the
rms residual stays at 3×10⁻⁴ Ry. Widening the window biases b_eq towards larger values because the
asymmetric tail pulls the parabola to the right; narrowing it to 3 points makes the fit exact (no
redundancy, no error estimate). Sensitivity check:

| half-width (Å) | # points | LDA b_eq (Å) | rms residual (Ry) | PBE b_eq (Å) |
|---|---|---|---|---|
| 0.06 | 3 | 0.7675 | 0 (exact) | 0.7545 |
| **0.10** | **5** | **0.7715** | 3.0×10⁻⁴ | **0.7587** |
| 0.14 | 7 | 0.7773 | 9.9×10⁻⁴ | 0.7650 |
| 0.30 | 12 | 0.7990 | 4.2×10⁻³ | 0.7830 |

The 0.04 Å grid therefore limits b_eq to roughly ±0.005 Å; a finer grid or a Morse/Birch-type fit
would be needed for better precision.

## (c) PBE (GGA) scan — `pbe.dat`

| b (Å) | E (Ry) | b (Å) | E (Ry) |
|---|---|---|---|
| 0.60 | −2.28370102 | 0.92 | −2.30067764 |
| 0.64 | −2.30503490 | 0.96 | −2.28988421 |
| 0.68 | −2.31801191 | 1.00 | −2.27810958 |
| 0.72 | −2.32447971 | 1.04 | −2.26560794 |
| **0.76** | **−2.32589317** | 1.08 | −2.25258417 |
| 0.80 | −2.32340442 | 1.12 | −2.23920953 |
| 0.84 | −2.31792429 | 1.16 | −2.22563023 |
| 0.88 | −2.31016518 | 1.20 | −2.21197218 |

Lowest grid point: **b = 0.76 Å, E = −2.32589317 Ry**. Plot: `pbe_curve.png`.
Same fit window (|b − 0.76| ≤ 0.10 Å, 5 points):

| α (Ry) | β (Ry/Å) | γ (Ry/Å²) | b_eq (Å) | E(b_eq) (Ry) |
|---|---|---|---|---|
| −1.611497 | −1.883171 | 1.240985 | **0.7587** | −2.325916 |

## Comparison with experiment

| Method | b_eq (Å) | Error vs 0.7414 Å |
|---|---|---|
| LDA (`H.pz-vbc.UPF`) | 0.772 ± 0.005 | +4.1 % (+3.5 to +4.8 %) |
| GGA-PBE (`H.pbe-vbc.UPF`) | 0.759 ± 0.005 | +2.3 % (+1.7 to +3.0 %) |
| Experiment (Huber & Herzberg 1979) | 0.7414 | — |

Both functionals overestimate the H–H bond length; PBE is closer to experiment. Physics: LDA uses only
the local density (every volume element treated as uniform electron gas) and is known to overbind;
PBE adds a dependence on |∇n|, which matters for a small molecule whose density is strongly peaked at
the nuclei and decays exponentially outside, and corrects the LDA balance between bonding and atomic
regions. The PBE curvature γ is slightly smaller (softer bond). Absolute total energies of LDA and
PBE are *not* comparable (different functional and pseudopotential) — only differences within one
set-up are meaningful.

### Convergence test (`run_convergence.py`, `conv/summary.dat`)

The five fit points were recomputed with ecutwfc = 100 Ry and/or A = 16 Å and refitted:

| ecut (Ry) | A (Å) | LDA ΔE(0.76) (mRy) | LDA b_eq (Å) | PBE ΔE(0.76) (mRy) | PBE b_eq (Å) |
|---|---|---|---|---|---|
| 60 | 12 | 0 | 0.7715 | 0 | 0.7587 |
| 100 | 12 | −3.5 | 0.7707 | −2.3 | 0.7580 |
| 60 | 16 | −0.007 | 0.7715 | −0.006 | 0.7588 |
| 100 | 16 | −3.5 | 0.7707 | −2.3 | 0.7580 |

Box size and Γ-only sampling are fully converged (< 1e-5 Ry, < 1e-4 Å). The cutoff still lowers the
absolute energy by 2–3 mRy but shifts b_eq by only −0.0008 Å, far below the ±0.005 Å grid/fit
uncertainty. The remaining offset of the PBE result from the all-electron PBE value (~0.750 Å,
Ernzerhof & Scuseria JCP 110, 5029 (1999); NIST CCCBDB) is thus due to the von Barth–Car
pseudopotential and the coarse 0.04 Å grid + parabolic fit, not to the numerical settings.

Combined plot: `both_curves.png`. Full write-up: `report.pdf` (`report.tex`).

## Files

| File | Content |
|---|---|
| `run_scan.py` | generates inputs from `H2.in`, runs `pw.x`, checks convergence, writes `lda.dat` / `pbe.dat` |
| `analyze.py` | plots, quadratic fits, b_eq, comparison table (prints the numbers above) |
| `run_convergence.py`, `conv/` | cutoff / box-size convergence test at the five fit points |
| `report.tex`, `report.pdf` | LaTeX submission |
| `lda.dat`, `pbe.dat` | b (Å), E_total (Ry) |
| `lda_curve.png`, `pbe_curve.png`, `both_curves.png` | plots |
| `runs/*.in`, `runs/*.out` | all pw.x inputs/outputs |

Reproduce: `python3 run_scan.py lda && python3 run_scan.py pbe && python3 analyze.py`
(on this VM `HWLOC_COMPONENTS=-x86` is set by the script to work around an OpenMPI/hwloc crash).
