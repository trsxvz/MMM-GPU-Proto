# Tangent by kernels end to end, 8 processes (11/10/2026)

Script: `cases/measure-positions.sh`. Stationary thermal r5 (279,873
unknowns, 5 steps) and beam r4 (10 steps), 8 single-threaded processes,
GMRES or conjugate gradient of hypre with BoomerAMG. Three computations per
case: point by point, by kernels with mfem-mgis installed in
`install/mfem-mgis` (8a6325335e, before), by kernels with mfem-mgis
installed in `install/mfem-mgis-registers` (4a007e85f2: positions of the
entries in the sparse matrix and matrices without forming B, after). The
line "mfem-mgis" of `summary.txt` gives the version of the sources, not
that of the installations. Quiet machine, GPU unused.

Times in seconds (`cases/summary-profiling.py`, number of calls in
parentheses):

| Case | Total | Integration | Residual | Tangent | Solution | Prediction |
|---|---|---|---|---|---|---|
| beam, point by point | 15.70 | 0.51 | 0.45 | 2.25 (31) | 9.56 | 2.87 |
| beam, kernels before | 13.05 | 0.18 | 0.08 | 1.40 (29) | 8.76 | 2.57 |
| beam, kernels after | 12.25 | 0.18 | 0.08 | 0.81 (29) | 8.71 | 2.42 |
| thermal, point by point | 9.44 | 2.15 | 2.30 | 3.12 (18) | 1.59 | |
| thermal, kernels before | 4.21 | 0.51 | 0.23 | 1.64 (18) | 1.56 | |
| thermal, kernels after | 3.65 | 0.51 | 0.22 | 1.06 (18) | 1.59 | |

- The tangent by kernels gains ×1.73 (beam) and ×1.55 (thermal); it is now
  2.8 to 2.9 times faster than the point by point tangent. The rest of the
  phase is mostly the PtAP product of hypre.
- Totals: ×1.07 (beam, the solution dominates) and ×1.15 (thermal).
- Same Newton iterations and same norms of the solution as the computation
  before these work packages.
