# Inexact Newton on the polycrystal (10/10/2026)

Polycrystal of mm-opera-hpc: 5 grains, single crystal behaviour of UO2,
periodic problem, uniaxial tension, 10 steps. mfem-mgis 4b72095ee0.
Script: `cases/measure-L5-polycrystal.sh` (options `--non-linear-solver`
and `--adaptive-linear-solver-tolerance` of the polycrystal). One RTX
5080, 16 cores, quiet machine (load 1.9).

All the kernels. Newton: PCG of hypre and BoomerAMG rebuilt at each
iteration. NewtonK0: conjugate gradient preconditioned by BoomerAMG built
once on the elastic operator. Tolerance of the linear solver of the
polycrystal (1e-12), or Eisenstat-Walker 1 with an initial tolerance of
0.1 (EW1).

| Case | CPU, 8 processes | GPU |
|---|---|---|
| Newton | 29.60 s (98 it.) | 11.88 s (98 it.) |
| Newton, EW1 | 33.96 s (156 it.) | 11.57 s (162 it.) |
| NewtonK0 | 25.02 s (98 it.) | 11.22 s (98 it.) |
| NewtonK0, EW1 | **23.46 s** (156 it.) | **9.58 s** (164 it.) |

Total times (`cases/summary-profiling.py`), total Newton iterations (one
linear solve each).

Detail on the GPU:

| Case | Integration | Solve |
|---|---|---|
| Newton | 5.09 s (158) | 5.11 s (98) |
| NewtonK0 | 5.16 s (159) | 4.45 s (98) |
| NewtonK0, EW1 | 7.12 s (225) | 0.70 s (164) |

## Findings

- **NewtonK0 with Eisenstat-Walker 1 is the best on both machines**: ×1.24
  against Newton on the GPU, ×1.26 on 8 cores. The linear solve becomes
  negligible on the GPU (0.70 s against 5.11 s).
- **The behaviour now dominates**: Eisenstat-Walker needs 66 more Newton
  iterations, hence as many integrations of the behaviour (7.1 s out of
  9.6 s on the GPU). The behaviour in single precision or a better
  behaviour kernel are the next levers on this case.
- With the classical Newton, Eisenstat-Walker does not pay off on CPU (each
  iteration rebuilds BoomerAMG), and barely on the GPU.
- Macroscopic stresses: differences from Newton of 9.6e-7 to 3.5e-6 with
  Eisenstat-Walker (2.4e-13 without). The microscopic Newton of the
  polycrystal converges up to 1e-4 in relative terms: with a linear
  tolerance of 1e-12, the quadratic convergence goes far beyond this
  criterion; with Eisenstat-Walker, the solution stops closer to the
  criterion. The differences are within the tolerance of the computation.
