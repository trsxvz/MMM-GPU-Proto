# Steady-state thermal by kernels (10/10/2026)

Nonlinear steady-state thermal problem in the beam (`cases/thermal`, r5,
279,873 unknowns, 1.57 million tetrahedra, 5 steps), conductivity
depending on the temperature, temperature imposed at both ends. Script:
`cases/measure-L7-thermal.sh`. One RTX 5080, 16 cores, load 2.4 at the
start. mfem-mgis 6cae30af99 plus the computation of the jacobians by a
kernel (8b0b47acfb, `summary.txt` gives the base commit).

The thermal tangent is negative definite and the diagonal of the imposed
temperatures is one: GMRES replaces the conjugate gradient (GMRES of hypre
and BoomerAMG, tolerance 1e-8, or Eisenstat-Walker 1; GMRES of MFEM and
Jacobi matrix free).

| Case | Before (0542) | After | Newton iterations |
|---|---|---|---|
| CPU 8 processes, point by point | 9.17 s | 9.03 s | 18 |
| CPU 8 processes, kernels | 4.89 s | 4.76 s | 18 |
| GPU, kernels | 3.01 s | **1.25 s** | 18 |
| GPU, kernels, Eisenstat-Walker 1 | 4.03 s | 1.36 s | 25 |
| GPU, kernels, matrix free | 183.40 s | 179.09 s | 18 |

Time of the loop over the steps. Same norm of the solution in all the
cases (306,232.626).

Detail (`cases/summary-profiling.py`, seconds and number of calls):

| Case | Integration | Residual | Tangent | Solve |
|---|---|---|---|---|
| CPU 8 processes, point by point | 2.11 (23) | 2.22 (23) | 3.05 (18) | 1.60 (18) |
| CPU 8 processes, kernels | 0.50 (23) | 0.23 (23) | 2.45 (18) | 1.52 (18) |
| GPU, kernels, before | 2.03 (23) | 0.03 (23) | 0.34 (18) | 0.54 (18) |
| GPU, kernels | 0.27 (23) | 0.03 (23) | 0.35 (18) | 0.53 (18) |

## Findings

- **The GPU does ×3.8 against 8 processes with the kernels, ×7.2 against
  point by point.** Before, the geometry of the elements cost 1.76 s on
  the host at the first integration (finite element space of order one and
  projection of the coordinates): 60 % of the time on the GPU. The kernel
  of the jacobians brings it down to a few milliseconds.
- **On CPU, the tangent dominates** (2.45 s out of 4.76 s with the
  kernels, 3.05 s point by point): the profiling does not yet separate the
  kernel of the element matrices, the sparse assembly and the parallel
  matrix of hypre. This is the next measurement on CPU.
- **Eisenstat-Walker does not pay off here**: the Newton converges in 3 or
  4 iterations, quadratically. The inexact Newton needs 25 iterations
  instead of 18, and each solve is already short.
- **Matrix free, Jacobi is not enough**: 178.7 s of solve, 10 s per solve,
  against 0.03 s with BoomerAMG. It would need a preconditioner built once
  (AMG on a reference tangent, like NewtonK0 in mechanics).
