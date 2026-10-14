# Fixed point on the beam, CPU and GPU (10/10/2026)

Beam r4 (111,843 unknowns, tetrahedra of order 1), von Mises plasticity
with linear hardening, 10 steps, elastic prediction, relative tolerance of
the nonlinear solver 1e-8. mfem-mgis d05475d265. Script:
`cases/measure-L5.sh`. One RTX 5080, 16 cores. A review was running at the
same time (average load of 4.5 at the start): the times on CPU are
overestimated by about 20 % (Newton: 16.6 s on a quiet machine, 20.2 s
here).

All the computations use the kernels (gradients, residual, assembled
tangent); `generic` behaviour on CPU (8 processes), `cuda` behaviour on the
GPU.

- **Newton**: best Newton measured, PCG + BoomerAMG of hypre, linear
  tolerance 1e-4;
- **fixed point** on the elastic operator K0 (`FixedPoint`, method
  `ElasticOperator`, K0 assembled once), accelerated by Anderson (depth 8)
  or by the algorithm of Cast3M, K0 solved by one cycle of BoomerAMG per
  iteration (default settings of hypre, or Chebyshev smoother,
  `RelaxationType=16`), or by PCG + BoomerAMG; prediction solved up to
  1e-8;
- **modified Newton** accelerated by Anderson (tangent of the first
  iteration of the step).

| Case | CPU, 8 processes | GPU |
|---|---|---|
| Newton | 20.27 s (29 it.) | **3.43 s** (32 it.) |
| fixed point, Anderson 8, one AMG cycle | **13.33 s** (341 it.) | 4.80 s (1,231 it.) |
| same, Chebyshev smoother | 13.25 s (330 it.) | **3.44 s** (438 it.) |
| fixed point, Cast3M, one AMG cycle | 14.36 s (477 it.) | 4.32 s (1,162 it.) |
| fixed point, Anderson 8, PCG + AMG | 31.52 s (94 it.) | 8.49 s (96 it.) |
| modified Newton, Anderson 8 | 29.67 s (74 it.) | 7.82 s (76 it.) |

Total profiled time (maximum over the processes).

Breakdown on the GPU (seconds, number of calls):

| Case | Integration | Residual | Solve | Acceleration | Prediction |
|---|---|---|---|---|---|
| Newton | 0.05 (42) | 0.02 (42) | 2.05 (32) | | 0.78 (10) |
| fixed point, Chebyshev | 0.39 (449) | 0.23 (448) | 0.83 (438) | 0.12 (438) | 1.49 (10) |

## Findings

- **On CPU, the fixed point wins**: 13.3 s against 20.3 s for Newton
  (×1.5, to compare with the ×1.8 of the prototype of 04/10). The fixed
  point avoids the assembly of the tangent and the setup of BoomerAMG at
  each iteration. One iteration (behaviour, residual, one AMG cycle,
  Anderson) costs about 21 ms against 560 ms for one Newton iteration
  (times without the prediction divided by the number of iterations).
- **On the GPU, a tie** (3.44 s against 3.43 s), and only with the
  Chebyshev smoother: with the default settings of hypre on GPU (l1-Jacobi,
  PMIS coarsening), one cycle of BoomerAMG is much weaker (1,231 iterations
  instead of 438). The prototype (×4 on GPU) used its own smoothed
  aggregation AMG with the rigid body modes. Chebyshev does not benefit
  Newton (3.95 s).
- **The prediction weighs**: 1.49 s out of 3.44 s on the GPU (6.2 s out of
  13.3 s on CPU), because it is solved up to 1e-8 by PCG + BoomerAMG, with
  a new setup of BoomerAMG at each step. Without prediction, the fixed
  point is slower (1,394 iterations on the GPU, 23 s on CPU): Anderson and
  one AMG cycle solve the large elastic increment badly. Lead: do the
  prediction with K0 and its AMG already built (planned).
- One iteration of the fixed point on the GPU costs about 4 ms (behaviour
  0.9 ms, residual 0.5 ms, AMG cycle 1.9 ms, Anderson 0.3 ms): the latency
  of the launches still dominates (see `results/L6-beam-2026-10-09/`).
- The numbers of iterations of the fixed point on the GPU vary from one
  run to the other (from 1,128 to 2,133 for the same configuration with the
  default smoother): hypre on GPU and the atomic additions are not bitwise
  reproducible, and Anderson amplifies these differences.

## Next: prediction by the AMG of K0 (mfem-mgis 846551e0cb, 1ef05104d1)

Option `PredictionLinearSolver=ElasticOperator` of the fixed point
(`-pk0`): the elastic prediction is solved by the conjugate gradient
preconditioned by the BoomerAMG cycle of K0, already built, instead of a
new setup of BoomerAMG at each step. Fixed point, Anderson 8, one Chebyshev
AMG cycle, K0 kept, on the GPU:

| Tolerance of the prediction | Iterations | Total | Prediction | Solves |
|---|---|---|---|---|
| 1e-4 | 534 | **2.80 s** | 0.77 s | 1.01 s |
| 1e-5 | 535 | 3.13 s | 1.10 s | 1.02 s |
| 1e-6 | 458 | 2.94 s | 1.18 s | 0.88 s |
| 1e-7 | 605 | 3.60 s | 1.31 s | 1.15 s |
| 1e-8 | 490 | 3.31 s | 1.41 s | 0.96 s |

With a prediction at 1e-4, three runs of each: Newton 3.15 s, 3.16 s,
3.18 s; fixed point 3.32 s, 2.88 s, 2.97 s (566 to 674 iterations). **On
the GPU, the fixed point is on par with Newton or slightly faster** (median
2.97 s against 3.16 s). On 8 cores it stays clearly faster.

- The prediction costs about 40 ms per step in steady state (assembly by
  kernels 3.3 ms, global assembly 2.6 ms, conjugate gradient about 35 ms).
  The first step also carries the one-off setups (geometry of the elements
  214 ms, pattern of the matrix), which the averages per call mix in.
- One iteration of the fixed point costs about 3.5 ms on the GPU:
  integration of the behaviour 0.96 ms (of which 0.65 ms of the behaviour
  kernel), residual 0.55 ms (the kernel of the internal forces, one thread
  per node, inverts each jacobian 4 times: to be written like the fused
  kernel of the action, one thread per element), AMG cycle 1.9 ms, Anderson
  0.3 ms.
- Settings of BoomerAMG tried in addition to Chebyshev (aggressive
  coarsening, threshold 0.5, interpolation 18, limited levels): no gain
  beyond the scatter of the numbers of iterations.

## Next: NewtonK0 method (mfem-mgis 30beae4764)

Newton with the consistent tangent operator, solved by the conjugate
gradient of MFEM preconditioned by BoomerAMG built once on K0
(`Method=NewtonK0`, `-fpm newton-k0 -els amg -k0 -pk0`, tolerance of the
conjugate gradient 1e-4), no setup of BoomerAMG during the computation
after the first step. Quiet machine (load 3).

| Case | GPU | CPU, 8 processes |
|---|---|---|
| Newton (PCG + BoomerAMG rebuilt at each iteration) | 3.16 s (32 it.) | 13.49 s (29 it.) |
| NewtonK0 | **2.58 to 2.66 s** (31 to 32 it.) | 11.56 s (29 it.) |
| fixed point, Anderson 8, one AMG cycle, prediction by K0 | 2.88 to 3.32 s | **9.51 s** (518 it.) |

Three runs of NewtonK0 on the GPU: 2.604 s, 2.579 s, 2.655 s, same numbers
of iterations as Newton: unlike the fixed point, the result is stable.

- **On the GPU, NewtonK0 is the best solver measured** (×1.22 against
  Newton): the setup of BoomerAMG at each Newton iteration cost more than
  the additional iterations of the conjugate gradient with the AMG of K0.
- **On CPU, the accelerated fixed point stays the best** (×1.42 against
  Newton), NewtonK0 in between (×1.17).
- With the matrix free tangent (`TangentOperatorComputation=MatrixFree`),
  NewtonK0 no longer assembles anything during the iterations (3.10 s on
  the GPU): it is the first good preconditioner of the matrix free path,
  still slower than the assembled tangent on these tetrahedra of order 1
  (see `results/L4-beam-2026-10-10-0123/`).
