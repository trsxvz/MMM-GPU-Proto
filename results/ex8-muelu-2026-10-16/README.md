# ex8 with MueLu (16/10/2026)

Example ex8 of mfem-mgis-examples (thermomechanical plate: Robin boundary
condition, swelling, iterative coupling), 5 steps, everything by kernels,
heat transfer solved by GMRES and BoomerAMG of hypre. The mechanics
(28,548 unknowns on the original mesh, 214,812 on the mesh refined once)
is solved by MUMPS or by a Krylov solver of MFEM preconditioned by MueLu
of Trilinos with the rigid body modes (`-pcMc MueLu`, options added to
ex8 on 16/10/2026: `-mls` chooses the smoother of MueLu, `-mld` the
damping factor of its Jacobi smoother). Build `build/examples-trilinos-ex8`
(ex8 linked to `install/mfem-mgis-trilinos`), GPU reserved, one run each.
Times: whole coupling over the steps
(`IterativeCouplingScheme::computeNextState`) and mechanics (`Mechanics`),
in seconds, from the profiler of each log.

The logs were produced by scripts of the working directory, partly with
intermediate versions of the options of ex8 (`-mlj`, replaced by `-mls`
and `-mld`); the options of each run are printed at the top of its log.
`cases/measure-ex8-muelu.sh` runs the same computations with the final
options.

## Original mesh

| target | solver of the mechanics | total | mechanics |
|---|---|---|---|
| CPU, 2 processes | MUMPS | 22.8 | 20.1 |
| CPU, 2 processes | CG, MueLu (symmetric Gauss-Seidel) | 57.2 | 54.3 |
| CPU, 2 processes | CG, MueLu (Chebyshev) | 65.1 | 62.2 |
| CPU, 2 processes | PCG of hypre, BoomerAMG, elasticity options (`../ex8-amg-elasticity-2026-10-11`) | 201.2 | |
| GPU, 1 process | MUMPS on the host (`-msMc SymmetricPositiveDefinite`) | 11.2 | 8.7 |
| GPU, 1 process | CG, MueLu (symmetric Gauss-Seidel) | 94.5 | 90.7 |
| GPU, 1 process | CG, MueLu (Chebyshev) | 16.8 | 13.1 |

- All the runs give the same statistics of the fields (minimum, maximum
  and mean of the temperature and of the norm of the displacement, to
  1e-14) and the same Newton iterations (114 lines `Newton iteration` in
  each log, heat transfer and mechanics).
- MueLu does not reduce the iterations on this plate: 480 iterations of
  the conjugate gradient per resolution on average (456 to 522, 29
  resolutions, `cpu-np2-cg-muelu-v1.txt`), against 441 for BoomerAMG with
  its elasticity options, and about 240 with its systems options (93
  seconds, measurement of 11/10/2026 noted in `JOURNAL.md`, raw output
  not kept). Its cycle is cheaper (57 seconds against 201 with the
  elasticity options), but MUMPS stays the fastest at this size, on the
  CPU and on the GPU.
- The Jacobi smoother damped by 0.6 (the choice of the beam) does not give
  a positive definite preconditioner here (`*-jacobi.txt`, option `-mlj`
  of an intermediate version: "PCG: The preconditioner is not positive
  definite", then failure of the first step). The runs `*-l1jacobi` used
  `relaxation: use l1`: in Ifpack2, it only adds half the absolute values
  of the couplings with the other processes to the small diagonal entries,
  hence an undamped Jacobi relaxation on one process: same failure. With
  a damping of 0.3, the GPU run fails, and the run on the coarse mesh
  (`check-coarse-jacobi-0.3.txt`) finishes after cut steps but differs
  from the reference statistics (1.1e-3 in relative value). The Chebyshev
  smoother of degree 2 (30 iterations of the power method, boost factor
  1.2) works, and is the one to use on a GPU: the symmetric Gauss-Seidel
  of MueLu is sequential and takes 90.7 seconds on the GPU.
- On the coarse mesh of the example (`check-coarse-*`, compared to its
  reference statistics), the conjugate gradient meets directions of
  negative curvature of the tangent operator in the first mechanical
  Newton iteration ("PCG: The operator is not positive definite", 3
  times, (Ad, d) down to -4e5), and continues: the tangent operator is
  not positive definite there, unlike on the original mesh (no such
  message). The test of the example (`rjh_plate-muelu`) therefore uses
  MINRES: on the coarse mesh, MINRES with MueLu gives the reference
  statistics without any such message (`check-coarse-minres-chebyshev`),
  and the test passes (`ctest-rjh_plate-muelu.log`, original mesh, the
  build being in Release mode).

## Mesh refined once (214,812 unknowns)

| target | solver of the mechanics | steps | total | mechanics |
|---|---|---|---|---|
| GPU, 1 process | CG, MueLu (Chebyshev) | 5 | 264.0 | 247.8 |
| GPU, 1 process | MINRES, MueLu (Chebyshev) | 5 | 205.1 | 192.5 |
| GPU, 1 process | GMRES, MueLu (Chebyshev) | 5 | failed | |
| GPU, 1 process | CG, MueLu (Chebyshev) | first | 39.7 | 37.4 |
| CPU, 2 processes | CG, MueLu (symmetric Gauss-Seidel) | first | 734.7 | 730.4 |

- The header of `cases/measure-ex8.sh` (11/10/2026) noted this mesh as out
  of reach: more than 20 minutes per step with MUMPS on 8 processes, 1,100
  seconds of PCG with BoomerAMG for 5 steps on 2 processes (about 220
  seconds per step for the linear solver alone). MINRES and MueLu on the
  GPU run the 5 steps in 205 seconds.
- The first step takes 39.7 seconds on the GPU with the Chebyshev
  smoother, 734.7 seconds on 2 CPU processes with the default smoother of
  MueLu. A run of the 5 steps on 2 CPU processes was stopped during its
  first step (not kept).
- The conjugate gradient meets directions of negative curvature of the
  tangent operator on this mesh: (Ad, d) from -4e-2 to -2e-12, 7 to 9
  times in the first step on the GPU (two runs), 8 on the CPU, the Newton
  iterations converging quadratically. Over the 5 steps, the third step
  (40,000 to 60,000 seconds) failed: its first mechanical resolution
  stopped after 10,000 iterations, with 1,377 of the 1,387 messages, and
  the step was cut in two (the log reports 40,000 to 50,000: the message
  of mfem-mgis described the sub-step instead of the failed step, fixed
  since). The conjugate gradient does not suit this tangent operator.
- MINRES (`-lsMc MINRESSolver`), which only needs a positive definite
  preconditioner, runs the 5 steps in 205 seconds without any failed
  resolution (110 lines `Newton iteration` against 136 for the conjugate
  gradient, whose cut step adds iterations). Its statistics of the fields
  differ from those of the conjugate gradient by 1.5e-6 at most in
  relative value (the cut step changes the time discretization of the
  conjugate gradient run). GMRES of MFEM, restarted with its default
  Krylov dimension, does not converge in 10,000 iterations on the first
  step. MINRES with MueLu is the choice for this mesh on the GPU.
