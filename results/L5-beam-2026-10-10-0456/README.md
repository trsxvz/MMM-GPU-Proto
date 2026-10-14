# Tangent without assembly on the GPU (10/10/2026)

Beam r4 of the prototype (111,843 unknowns, 10 steps), `cuda` behaviour,
one RTX 5080, kernels for the gradients and the internal forces. mfem-mgis
4b72095ee0. Script: `cases/measure-L5-matrix-free.sh`. Quiet machine (load
1.7 at the start, decreasing).

Tangent applied without assembly (`TangentOperatorComputation` set to
`MatrixFree`, action by the fused kernel) with three possible
preconditioners:

- the Jacobi smoother;
- the Chebyshev smoother of MFEM (`OperatorChebyshevSmoother`, orders 2 to 5);
- NewtonK0: conjugate gradient preconditioned by BoomerAMG built once on
  K0, with or without Eisenstat-Walker 1 and an initial tolerance of 0.1.

These cases are compared with NewtonK0 with the tangent assembled by
kernels. Tolerance of the linear solver: 1e-4.

| Case | Time | Newton iterations | Solve |
|---|---|---|---|
| matrix free, CG + Jacobi | 21.04 s | 30 | 19.50 s |
| matrix free, CG + Chebyshev 2 | 23.13 s | 32 | 21.41 s |
| matrix free, CG + Chebyshev 3 | 25.50 s | 30 | 23.72 s |
| matrix free, CG + Chebyshev 4 | 28.22 s | 30 | 26.33 s |
| matrix free, CG + Chebyshev 5 | 31.43 s | 32 | 29.45 s |
| matrix free, NewtonK0 | 3.05 s | 33 | 2.24 s |
| matrix free, NewtonK0 + Eisenstat-Walker 1 | 2.05 s | 37 | 1.24 s |
| assembled, NewtonK0 + Eisenstat-Walker 1 | **1.99 s** | 37 | 0.97 s |

Time of the loop over the steps; solve according to
`cases/summary-profiling.py`.

## Findings

- **Chebyshev does not pay off here**: each additional order adds one
  action of the matrix free tangent per application of the preconditioner
  (tetrahedra of order 1, action on 3 x 4 = 12 unknowns per element), and
  the time grows with the order. The iterations of the conjugate gradient
  are not recorded by these measurements: the reduction brought by each
  order is not quantified here (on a Laplacian, 46 iterations with Jacobi,
  27 and 17 with Chebyshev of order 2 and 4, test
  `OperatorSmoothersTest`). The smoother stays useful as a building block
  (smoothing of a geometric or p multigrid, high orders): it is in the
  factory and tested.
- **Matrix free, NewtonK0 with Eisenstat-Walker 1 is within 3 % of the
  assembled tangent** (2.05 s against 1.99 s, one run per case): the
  tangent is no longer assembled during the iterations. The memory is not
  reduced for all that: K0 (copy of an assembled tangent), the matrix of
  the problem and the matrix of the prediction stay assembled.
- Eisenstat-Walker brings a bit more matrix free (3.05 s to 2.05 s,
  -33 %) than with the assembled tangent (2.57 s to 1.96 s, -24 %,
  `results/L5-beam-2026-10-10-0429/`).
