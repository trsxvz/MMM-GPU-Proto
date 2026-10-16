# Solvers of Trilinos on the CPU, large systems (16/10/2026)

8 processes on the CPU, serial node of Tpetra (mfem-mgis 3a3909cec8, build
`trilinos-host`), `SIZE=32`: the beam of 421,443 unknowns, homogeneous and
with stiff inclusions, and the Laplace equation of 2,146,689 unknowns.
Whole machine reserved, two repetitions (the logs of the second one end
with `-rep2`); `summary.md` gives the minimum times. References from
`../elasticity` and `../laplace`.

| problem | configuration | second operator (s) | iterations | residual |
|---|---|---|---|---|
| elasticity | CG of MFEM + MueLu (default smoother) | 1.29 | 11 | 1.9e-6 |
| elasticity | CG of MFEM + MueLu, Jacobi smoother, hierarchy reused | 1.17 | 15 | 2.8e-6 |
| elasticity | CG + GAMG (PETSc), true residual | 1.76 | 26 | 5.8e-9 |
| elasticity | CG + BoomerAMG | 3.82 | 37 | 2.6e-7 |
| heterogeneous | CG of MFEM + MueLu | 2.12 | 26 | 3.8e-6 |
| heterogeneous | CG + GAMG, true residual | 3.03 | 59 | 7.4e-9 |
| heterogeneous | CG + BoomerAMG | 5.09 | 56 | 4.1e-7 |
| laplace | CG of MFEM + MueLu, Jacobi smoother, hierarchy reused | 1.39 | 12 | 9.7e-8 |
| laplace | CG of MFEM + MueLu | 1.55 | 9 | 9.7e-8 |
| laplace | CG + BoomerAMG | 1.93 | 14 | 4.3e-8 |

- MueLu is the fastest configuration measured on 8 cores at these sizes:
  1.29 seconds per operator on the elasticity (1.76 for GAMG, 3.82 for
  BoomerAMG), 2.12 with the inclusions (3.03 for GAMG), 1.55 on the
  Laplace equation (1.93 for BoomerAMG).
- The conjugate gradient of MFEM stops on the norm of the preconditioned
  residual: its true residual (1.9e-6 for MueLu on the elasticity, 2.6e-7
  for BoomerAMG) is larger than the one of GAMG stopped on the true
  residual (5.8e-9). The comparison with GAMG is therefore slightly
  favourable to MueLu; with the conjugate gradient of Belos, stopped on
  the true residual, MueLu needed 14 instead of 11 iterations at the
  smaller size (`../../solvers-2026-10-16-trilinos-cpu`).
- On the CPU, the default smoother of MueLu (symmetric Gauss-Seidel) is
  faster than the Jacobi smoother; reusing the hierarchy saves 0.25
  seconds on the elasticity and 0.36 on the Laplace equation.
