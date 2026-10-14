# Inexact Newton on the beam (10/10/2026)

Beam r4 of the prototype (111,843 unknowns, von Mises with linear
hardening, 10 steps), all the kernels, elastic prediction. mfem-mgis
517c0c0b08. Script: `cases/measure-L5-ew.sh`. One RTX 5080, 16 cores.
Load 3.7 at the start (end of compilations), results of Newton and of
NewtonK0 equal to those of `results/L5-beam-2026-10-10-0250/` within 1 %.

Tolerance of the linear solver fixed (1e-4, also for the predictions) or
adapted to the decrease of the residual (option `-alst`, option
`AdaptiveLinearSolverTolerance` of the Newton and of the `NewtonK0` method
of the fixed point: choices 1 and 2 of Eisenstat and Walker, default
settings of MFEM, initial tolerance 0.5). Newton: PCG of hypre and
BoomerAMG rebuilt at each iteration. NewtonK0: conjugate gradient of MFEM
preconditioned by BoomerAMG built once on K0, prediction by this conjugate
gradient (`-fpm newton-k0 -els amg -k0 -pk0`).

| Case | GPU | CPU, 8 processes |
|---|---|---|
| Newton, fixed tolerance | 3.12 s (32 it.) | 13.59 s (29 it.) |
| Newton, Eisenstat-Walker 1 | 3.01 s (53 it.) | 15.73 s (51 it.) |
| Newton, Eisenstat-Walker 2 | 3.32 s (63 it.) | 17.70 s (61 it.) |
| NewtonK0, fixed tolerance | 2.57 s (30 it.) | 11.56 s (29 it.) |
| NewtonK0, Eisenstat-Walker 1 | **2.14 s** (52 it.) | **10.93 s** (51 it.) |
| NewtonK0, Eisenstat-Walker 2 | 2.19 s (58 it.) | 12.30 s (61 it.) |

Time of the loop over the steps, total Newton iterations. For comparison,
the best fixed point measured on 8 cores (Anderson 8, one AMG cycle) takes
9.5 s (`results/L5-beam-2026-10-10-0250/`).

Solve (time and calls, `cases/summary-profiling.py`):

| Case | GPU | CPU, 8 processes |
|---|---|---|
| Newton, fixed tolerance | 2.01 s (32) | 8.73 s (29) |
| Newton, Eisenstat-Walker 1 | 1.78 s (53) | 9.37 s (51) |
| NewtonK0, fixed tolerance | 1.60 s (30) | 7.40 s (29) |
| NewtonK0, Eisenstat-Walker 1 | 1.05 s (52) | 5.21 s (51) |

## Findings

- **On the GPU, NewtonK0 with Eisenstat-Walker 1 is the best solver
  measured**: 2.14 s, ×1.46 against Newton and ×1.20 against NewtonK0 with
  a fixed tolerance. The solve goes from 1.60 s to 1.05 s despite 22 more
  Newton iterations: the solves are on average much cheaper (the
  iterations of the conjugate gradient are not recorded by these
  measurements), and one more iteration only costs one integration and one
  assembly by kernels (BoomerAMG is not rebuilt).
- **With the classical Newton, Eisenstat-Walker does not pay off**: each
  additional iteration rebuilds BoomerAMG and assembles the tangent. On 8
  cores, the computation is 16 % (choice 1) to 30 % (choice 2) slower.
- **Choice 1 against choice 2**: choice 1 (residual of the linear system,
  one more product by the tangent per iteration) does fewer Newton
  iterations (51 to 53 against 58 to 63) and wins everywhere.
- **Elastic steps**: with the fixed tolerance, the prediction and one
  iteration (8 processes) or two (GPU) are enough. With Eisenstat-Walker,
  the initial tolerance (0.5) gives 4 or 5 iterations to the elastic steps.
  A lower initial tolerance (option `InitialTolerance`) is a lead to gain
  more.
- **On CPU**, the accelerated fixed point stays the best (9.5 s), NewtonK0
  with Eisenstat-Walker 1 follows (10.9 s).

## Initial tolerance (NewtonK0, Eisenstat-Walker 1)

Option `-alsti` of the beam (`InitialTolerance`), same settings, outputs in
`initial-tolerance/` (run by hand, load 3.7 to 4.6 during the trials on 8
cores).

| Initial tolerance | GPU | CPU, 8 processes |
|---|---|---|
| 0.5 (default) | 2.19 s (52 it.) | 10.91 s (51 it.) |
| 0.1 | **1.96 to 2.02 s** (37 and 33 it.) | 9.66 s (38 it.) |
| 0.01 | 1.96 s (37 it.) | **9.46 s** (38 it.) |
| 0.001 | 2.07 s (36 it.) | |

- An initial tolerance of 0.1 or 0.01 gives the elastic steps back their
  two iterations: on the GPU, 1.96 s, that is ×1.59 against Newton
  (3.12 s); on 8 cores, 9.46 s, at the level of the best fixed point
  (9.5 s). On the GPU, two runs with 0.1 give 37 and 33 iterations (1.96 s
  and 2.02 s): hypre on the GPU is not bitwise reproducible, and the
  adapted tolerance depends on the residuals.
- With a fixed tolerance, the elastic steps do one iteration on 8
  processes and two on the GPU.
- The default value (0.5, the one of MFEM) stays conservative. The choice
  0.1 is a good setting for this case, to be confirmed on the polycrystal.
