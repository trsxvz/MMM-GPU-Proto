# Nonlinear beam with MueLu (16/10/2026)

`cases/measure-beam-muelu.sh`: the beam of the prototype (`cases/beam`,
tetrahedra, isotropic plasticity with linear hardening, 10 time steps,
Newton solver, linear tolerance 1e-4, behaviours integrated by kernels),
solved with the conjugate gradient preconditioned by BoomerAMG of hypre
(`pcg-amg`: conjugate gradient of hypre; `cg-amg`: conjugate gradient of
MFEM) or by MueLu (`cg-muelu`: two sweeps of a Jacobi smoother;
`cg-muelu-sgs`: symmetric Gauss-Seidel), the prolongators and the
restrictions of MueLu being kept from one tangent operator to the next
(`reuse: type` set to `RP`). Refinement 4 (111,843 unknowns) on the GPU
(one process) and on 8 CPU processes, refinement 5 (839,619 unknowns) on
the GPU. Cases built against mfem-mgis installed apart
(`install/mfem-mgis-trilinos`, branch `trilinos-host`), whole machine
reserved, one repetition. All the runs give the same norm of the
displacement (to 1e-8) and the same Newton iterations (28 to 35).

| refinement | target | linear solver | time (s) |
|---|---|---|---|
| 4 | GPU | `pcg-amg` | 2.99 |
| 4 | GPU | `cg-amg` | 3.08 |
| 4 | GPU | `cg-muelu`, factories of MueLu written with Kokkos | 5.91 |
| 4 | GPU | `cg-muelu`, factories not written with Kokkos (`-classic`) | 3.21 |
| 4 | CPU, 8 processes | `pcg-amg` | 12.09 |
| 4 | CPU, 8 processes | `cg-muelu-sgs` | 6.06 to 6.12 |
| 4 | CPU, 8 processes | `cg-muelu` | 6.53 to 6.76 |
| 5 | GPU | `pcg-amg` | 15.81 |
| 5 | GPU | `cg-amg` | 16.57 |
| 5 | GPU | `cg-muelu`, factories written with Kokkos | 41.27 |
| 5 | GPU | `cg-muelu`, factories not written with Kokkos (`-classic`) | 15.61 |

- On 8 CPU processes, MueLu halves the time of the computation (6.1
  seconds against 12.1 with BoomerAMG): the conjugate gradient needs 9.6
  iterations per resolution with MueLu, 25.2 with BoomerAMG (runs with the
  iterations printed, not stored).
- On the GPU, MueLu first needed 78 iterations per resolution (48 for
  BoomerAMG), against 11 for the same MueLu on one CPU process: the
  factories of MueLu written with Kokkos, which MueLu uses on a GPU, give a
  worse coarse space with the rigid body modes (cause: the local QR of
  their tentative prolongator, see `../muelu-factories-2026-10-16`). With
  the other factories
  (`use kokkos refactor` set to false, the default of mfem-mgis since
  d08224ba8f), MueLu is as fast as BoomerAMG on the GPU (3.21 against
  2.99 seconds at refinement 4, 15.61 against 15.81 at refinement 5).
- The runs `-pu4` (preconditioner built every 4 Newton iterations, option
  `-pu 4`) gave the same times as the runs without it, for BoomerAMG and
  for MueLu: the preconditioner `AssembledOperator` of mfem-mgis is set up
  again whenever the operator has been assembled again, which happens at
  each Newton iteration with the tangent operator assembled by kernels
  (`MFEMLinearSolvers.cxx`); the period only applies to the tangent
  operator applied without assembly.
- The first runs (before d08224ba8f) used mfem-mgis 44158ee740.

## Tangent operator applied without assembly (16/10/2026, 10 h 53)

GPU only, refinements 4 and 5, mfem-mgis 7ab30d43f4 (MueLu declared as an
approximate linear solver), tangent operator applied without assembly
(`TANGENT=MatrixFree`), two repetitions (the runs `-rep2`), GPU reserved:

- `-mf-k0-amg` and `-mf-k0-muelu`: Newton-K0 method of the fixed point
  solver (`-nls fixed-point -fpm newton-k0 -k0 -pk0 -alst
  EisenstatWalker1 -alsti 0.1`): conjugate gradient on the tangent
  operator, preconditioned by one cycle of a multigrid built once on the
  elastic operator, BoomerAMG (`-els amg`) or MueLu with two sweeps of a
  Jacobi smoother (`-els muelu`). The linear solver of the problem
  (`cg-amg` in the names of the runs) is not used.
- `-mf-pu4`: Newton, the preconditioner (`cg-amg` or `cg-muelu`, built on
  the operator assembled for it) set up again every 4 Newton iterations
  (`-pu 4`).

| refinement | method | multigrid | time (s) | Newton iterations |
|---|---|---|---|---|
| 4 | Newton-K0 | BoomerAMG | 1.87 to 2.06 | 33 |
| 4 | Newton-K0 | MueLu | 0.85 to 0.87 | 35 |
| 4 | Newton, `-pu 4` | BoomerAMG | 2.91 to 2.98 | 31 to 32 |
| 4 | Newton, `-pu 4` | MueLu | 1.46 to 1.48 | 31 |
| 5 | Newton-K0 | BoomerAMG | 11.24 to 11.58 | 38 to 39 |
| 5 | Newton-K0 | MueLu | 6.53 to 6.60 | 34 |
| 5 | Newton, `-pu 4` | BoomerAMG | 18.63 to 18.77 | 32 |
| 5 | Newton, `-pu 4` | MueLu | 10.80 to 10.84 | 32 |

- Newton-K0 with MueLu is the fastest configuration measured on this beam:
  6.5 seconds at refinement 5, against 11.2 with BoomerAMG and 15.6 to
  15.8 for Newton with the tangent operator assembled (first table), 2.4
  times faster. At refinement 4, 0.85 seconds against 1.87 and 3.0.
- With the tangent operator applied without assembly and the
  preconditioner set up every 4 Newton iterations, MueLu is 1.7 times
  faster than BoomerAMG (10.8 against 18.7 seconds at refinement 5).
- The norm of the displacement of Newton-K0 differs from that of Newton by
  1.4e-6 in relative value at refinement 5 with MueLu (133.69566 against
  133.69585), by less than 1e-8 with BoomerAMG: the linear resolutions of
  Newton-K0 are inexact (Eisenstat-Walker tolerances).
