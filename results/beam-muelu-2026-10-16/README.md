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
  worse coarse space with the rigid body modes. With the other factories
  (`use kokkos refactor` set to false, the default of mfem-mgis since
  d08224ba8f), MueLu is as fast as BoomerAMG on the GPU (3.21 against
  2.99 seconds at refinement 4, 15.61 against 15.81 at refinement 5).
- The runs `-pu4` (preconditioner built every 4 Newton iterations, option
  `-pu 4`) gave the same times as the runs without it, for BoomerAMG and
  for MueLu (not investigated).
- The first runs (before d08224ba8f) used mfem-mgis 44158ee740.
