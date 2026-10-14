# Measurements of work package L2 on the beam (09/10/2026)

Beam r4 (111,843 unknowns), von Mises plasticity with linear hardening,
PCG + BoomerAMG (tolerance 1e-4), 10 steps. mfem-mgis 4a68166ca9,
`generic-parallel` behaviour compiled with the fixed options (patch 0093 of
TFEL), quiet machine. Supersedes the provisional measurements of
`results/L2-2026-10-09/` (loaded machine, badly placed threads).

| Case | Total | Integration | Residual | Tangent | Solve | Prediction |
|---|---|---|---|---|---|---|
| 8 processes, point by point | 15.62 | 0.49 | 0.43 | 2.22 | 9.55 | 2.87 |
| 8 processes, batches, `generic` | 15.56 | 0.50 | 0.43 | 2.21 | 9.52 | 2.84 |
| 8 processes, batches, `openmp` 1 thread | 15.53 | 0.49 | 0.43 | 2.21 | 9.49 | 2.83 |
| 1 process, point by point | 65.36 | 3.39 | 2.92 | 11.63 | 33.73 | 13.50 |
| 1 process, batches, `generic` | 65.06 | 3.19 | 2.91 | 11.55 | 33.63 | 13.56 |
| 1 process, batches, `openmp` 1 thread | 65.10 | 3.22 | 2.95 | 11.61 | 33.65 | 13.46 |
| 1 process, batches, `openmp` 8 threads | 64.82 | 2.63 | 2.94 | 11.55 | 33.96 | 13.49 |
| 1 process, batches, `stlpar` 8 cores | 65.11 | 2.68 | 2.95 | 11.57 | 34.15 | 13.54 |
| 2 processes, point by point | 36.62 | 1.67 | 1.53 | 6.53 | 19.17 | 7.57 |
| 2 processes, batches, `openmp` 4 threads | 36.34 | 1.35 | 1.56 | 6.54 | 19.26 | 7.49 |

Same solution everywhere (norm 48.83704586).

Finding: the behaviour weighs little (less than 5 % of the time). Several
threads only reduce its share, not the computation of the gradients element
by element, which stays in the integration phase. The gains on the beam will
come from the solve (more than 50 %), the tangent and the prediction
(work packages L3 to L5).
