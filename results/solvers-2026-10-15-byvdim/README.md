# Linear solvers, unknowns ordered by vector dimension (15/10/2026)

Second pass of `cases/measure-solvers.sh` for the configurations that the
first pass (`results/solvers-cpu-2026-10-15`, `results/solvers-gpu-2026-10-15`)
did not measure correctly:

- `petsc-cg-gamg`: conjugate gradient of PETSc preconditioned by GAMG, now
  given the block size and the rigid body modes (mfem-mgis patches 0232
  and 0233), which requires the unknowns ordered by vector dimension (`-bv`);
- `cg-boomeramg-elasticity`: BoomerAMG with its elasticity options, which
  require the same ordering;
- `cg-boomeramg`: BoomerAMG with its systems options (default), as a
  reference measured in the same conditions.

Same problems and sizes as the first pass, 8 processes on the CPU, one
process on the GPU (RTX 5080). Build `solvers` of mfem-mgis 0399b20670
(patches up to 0237), quiet machine. Each time is the minimum of three
repetitions (`runs.txt`, `repeat-2`, `repeat-3`). `petsc-true-residual.txt`:
GAMG with the PETSc option `ksp_norm_type=unpreconditioned`, so that the
tolerance applies to the true residual (one repetition). Table:
`summary.md`.

## Main results (CPU, 8 processes)

| problem | configuration | iterations | true residual | first (s) | second operator (s) |
|---|---|---|---|---|---|
| elasticity | BoomerAMG, systems | 30 | 1.2e-7 | 0.203 | 0.203 |
| elasticity | BoomerAMG, elasticity | 31 | 1.1e-7 | 0.647 | 0.661 |
| elasticity | GAMG, rigid body modes | 15 | 5.4e-6 | 0.116 | 0.116 |
| elasticity | GAMG, true residual | 21 | 2.8e-9 | 0.149 | 0.119 |
| heterogeneous | BoomerAMG, systems | 47 | 2.5e-7 | 0.263 | 0.276 |
| heterogeneous | GAMG, true residual | 39 | 5.5e-9 | 0.144 | 0.134 |

- With the rigid body modes, GAMG needs 15 iterations instead of 177 on
  the elasticity, and 25 instead of 284 with the stiff inclusions.
- With the same relative tolerance, the default conjugate gradient of
  PETSc stops on the preconditioned residual: its true residual is 5e-6.
  Stopped on the true residual (21 iterations), GAMG is the fastest
  multigrid measured on the CPU for the elasticity: 0.12 seconds per
  operator against 0.20 for BoomerAMG, for a residual 40 times smaller.
  Its cycles are cheaper (about 1 ms per iteration against 4 ms). At this
  size, the conjugate gradient with the Jacobi preconditioner is still
  faster on the homogeneous beam (0.089 seconds per operator, 319
  iterations, `../solvers-cpu-2026-10-15`); with the stiff inclusions,
  GAMG is the fastest configuration (0.134 seconds against 0.25).
- The elasticity options of BoomerAMG do not reduce the iterations on
  this beam (31 against 30) and cost three times more. On the GPU, hypre
  does not support them: mfem-mgis falls back to the systems options,
  with a warning (53 iterations against 54).
- Laplace: BoomerAMG stays faster than GAMG (0.13 against 0.29 seconds per
  operator), whose setup dominates.

## Defects found by this pass (fixed in mfem-mgis)

- The first setup of the solvers of PETSc took 0.4 second more than the
  following ones, and each process opened a context on the GPU (262 MiB,
  2.1 GiB for 8 processes) although PETSc computes on the host in
  mfem-mgis: MFEM converts the matrices of hypre through the `MATHYPRE`
  format of PETSc, whose creation initializes the GPU when PETSc and hypre
  are built for it. mfem-mgis now copies the matrices itself (patches
  0235, 0236 and 0238): first setup 0.011 second instead of 0.41,
  following ones 0.010 instead of 0.015.
- Requesting the elasticity options of BoomerAMG for a problem with
  scalar unknowns made MFEM abort: they are refused (patches 0234, 0236, 0238).
- The generic test of the linear solvers failed for AmgX in parallel on
  the GPU: `mfem::ParNormlp` of MFEM 4.10 recomputes the norm from the
  values of the host without synchronizing them (patch 0237).
