# ex8 with MueLu (16/10/2026)

Example ex8 of mfem-mgis-examples (thermomechanical plate: Robin boundary
condition, swelling, iterative coupling), 5 steps, everything by kernels,
heat transfer solved by GMRES and BoomerAMG of hypre. The mechanics
(28,548 unknowns on the original mesh, 214,812 on the mesh refined once)
is solved by MUMPS or by the conjugate gradient of MFEM preconditioned by
MueLu of Trilinos with the rigid body modes (`-lsMc CGSolver -pcMc
MueLu`, options added to ex8 on 16/10/2026: `-mls` chooses the smoother
of MueLu, `-mld` the damping factor of its Jacobi smoother). Build
`build/examples-trilinos-ex8` (ex8 linked to `install/mfem-mgis-trilinos`),
`cases/measure-ex8-muelu.sh`, GPU reserved, one run each. Times: whole
coupling over the steps (`IterativeCouplingScheme::computeNextState`) and
mechanics (`Mechanics`), in seconds, from the profiler of each log.

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
  its elasticity options. It is faster because its cycle is cheaper (57
  seconds against 201), but MUMPS stays the fastest at this size, on the
  CPU and on the GPU.
- The Jacobi smoother damped by 0.6 (the choice of the beam) or by 0.3
  does not give a positive definite preconditioner here (`*-jacobi*`:
  "PCG: The preconditioner is not positive definite", then failure of the
  first step). The runs `*-l1jacobi` used an intermediate version of the
  option (`relaxation: use l1`, which only adds the couplings with the
  other processes, hence an undamped Jacobi relaxation on one process):
  same failure. The Chebyshev smoother of degree 2 (30 iterations of the
  power method, boost factor 1.2) works, and is the one to use on a GPU:
  the symmetric Gauss-Seidel of MueLu is sequential and takes 90.7
  seconds on the GPU.
- The short test on the coarse mesh compared to the reference statistics
  of the example passes with MueLu (`check-coarse-*`, except the l1 and
  Jacobi 0.3 variants).

## Mesh refined once (214,812 unknowns)

| target | solver of the mechanics | steps | total | mechanics |
|---|---|---|---|---|
| GPU, 1 process | CG, MueLu (Chebyshev) | 5 | 264.0 | 247.8 |
| GPU, 1 process | CG, MueLu (Chebyshev) | first | 39.7 | 37.4 |
| CPU, 2 processes | CG, MueLu (symmetric Gauss-Seidel) | first | 734.7 | 730.4 |

- The header of `cases/measure-ex8.sh` (11/10/2026) noted this mesh as out
  of reach: more than 20 minutes per step with MUMPS on 8 processes, 1,100
  seconds of PCG with BoomerAMG for 5 steps on 2 processes. MueLu on the
  GPU runs the 5 steps in 264 seconds, and the first step 18.5 times
  faster than MueLu with its default smoother on 2 CPU processes.
- On this mesh, the conjugate gradient meets directions of negative
  curvature of the tangent operator ("PCG: The operator is not positive
  definite", with (Ad, d) from -4e-2 to -2e-12 in the first step on the
  GPU): 9 times in the first step on the GPU, 8 on the CPU (the Newton
  iterations converging quadratically), 1,387 times over the 5 steps, where
  one
  resolution of the mechanics failed (step from 40,000 to 50,000 seconds),
  the step being then cut and the computation finishing. The conjugate
  gradient does not suit a tangent operator which is not positive
  definite: GMRES or MINRES with MueLu are to be measured on this mesh.
