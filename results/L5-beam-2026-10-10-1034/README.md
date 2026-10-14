# Options of BoomerAMG for NewtonK0 (10/10/2026)

Beam r4 of the prototype (111,843 unknowns, 10 steps), all the kernels,
mfem-mgis a066c82baa. Best method of the previous measurements: NewtonK0
(conjugate gradient preconditioned by BoomerAMG built once on K0,
prediction by this conjugate gradient) with Eisenstat-Walker 1, initial
tolerance 0.1. Options of hypre passed by `-amgo`. Script:
`cases/measure-L5-amg.sh`. One RTX 5080 (`cuda` behaviour), 16 cores
(`generic` behaviour, 8 processes). Load 0.2 at the start.

| Options of BoomerAMG | GPU | Iterations | CPU, 8 processes | Iterations |
|---|---|---|---|---|
| default | 1.85 s | 33 | 8.78 s | 38 |
| `StrengthThreshold=0.5` | 1.87 s | 35 | 8.80 s | 38 |
| `StrengthThreshold=0.7` | 2.55 s | 38 | 10.05 s | 38 |
| `AggressiveCoarseningLevels=1` | **1.71 s** | 34 | 7.41 s | 34 |
| `StrengthThreshold=0.5`, `AggressiveCoarseningLevels=1` | 1.72 s | 37 | **7.22 s** | 34 |
| `CoarseningType=10` (HMIS) | 1.87 s | 38 | 9.74 s | 38 |
| `InterpolationType=14` | 1.79 s | 36 | 9.27 s | 35 |
| `RelaxationType=16` (Chebyshev) | 2.31 s | 37 | 9.15 s | 35 |
| `RelaxationSweeps=2` | 2.25 s | 33 | 10.45 s | 38 |
| `RelaxationType=7` (Jacobi) | failure | | failure | |

Time of the loop over the steps, one run per case; total Newton
iterations.

Three more runs on the GPU (less quiet machine, load 4): default 1.91,
2.05, 1.91 s; aggressive coarsening 1.87, 1.86, 1.73 s; threshold 0.5 and
aggressive coarsening 1.73, 1.81, 1.72 s; interpolation 14 1.93, 1.88,
1.82 s.

## Findings

- **Aggressive coarsening on one level is the best setting**: -16 to
  -18 % with 8 CPU processes (7.22 s with the threshold 0.5, against
  8.78 s), -5 to -8 % on the GPU. It reduces the cost of each cycle of
  BoomerAMG (smaller hierarchy) without increasing the number of solves.
- The more expensive smoothers (Chebyshev, two sweeps) and a strength
  threshold of 0.7 slow down. The plain Jacobi of hypre (7) makes the
  conjugate gradient fail from the first step (the preconditioner is no
  longer symmetric positive definite).
- Compared with the 05:00 measurement
  (`results/L5-beam-2026-10-10-0429`, NewtonK0 + EW1 with initial tolerance
  0.1: 1.96 s on GPU, 9.46 s on 8 cores), the kernels of the geometry and
  of the matrices on the host give 1.85 s and 8.78 s with the default
  options.
- On the GPU, the linear solve still dominates (0.94 to 1.07 s out of 1.7
  to 2.1 s, `cases/summary-profiling.py`), then the prediction (0.42 s).
