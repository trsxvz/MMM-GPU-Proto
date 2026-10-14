# Global matrix with pattern kept: end to end on 8 processes (11/10/2026)

Script: `cases/measure-global-matrix.sh`. Thermal r5 (5 steps) and beam r4
(10 steps), 8 single-threaded processes, tangent assembled by kernels.
Before: mfem-mgis 7af829b3a5 (`install/mfem-mgis-jacobians`); after:
e111df8454 (`install/mfem-mgis-global`, `batch::GlobalStiffnessMatrix`).
The line "mfem-mgis" of `summary.txt` gives the version of the sources, not
that of the installations. Quiet machine.

Times in seconds (`cases/summary-profiling.py`):

| Case | Loop over the steps | Tangent | Solution | Prediction |
|---|---|---|---|---|
| thermal, before | 3.38 | 1.04 (18) | 1.55 | |
| thermal, after | 2.74 | 0.43 (18) | 1.53 | |
| beam, before | 12.30 | 0.80 (29) | 8.77 | 2.44 |
| beam, after | 11.84 | 0.44 (29) | 8.74 | 2.38 |

- Tangent ×2.4 (thermal) and ×1.8 (beam); loop over the steps ×1.23 and
  ×1.04 (the solution dominates the beam).
- Same Newton iterations and same norms of the solution.
- Since the beginning of the day (point by point tangent, 3.12 s for the
  thermal, 2.25 s for the beam, `results/positions-2026-10-11-1008`), the
  tangent phase on 8 processes is 7 and 5 times faster.
