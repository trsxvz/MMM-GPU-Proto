# Solvers of Trilinos on the CPU, size of the first passes (16/10/2026)

Since the commit 89463ba63b of mfem-mgis, the solvers of Trilinos compute
where MFEM computes: these runs labelled `cpu` are 8 processes on the CPU
(serial node of Tpetra), unlike those of `../solvers-2026-10-16-trilinos`,
which shared the GPU. FROSch and Tacho still compute on the GPU (default
node of Tpetra, see the user guide of mfem-mgis): their runs are 8
processes sharing the GPU. Beam of 56,355 unknowns (homogeneous and with
stiff inclusions), Laplace equation of 274,625 unknowns. Build
`trilinos-host` of mfem-mgis 3a3909cec8, whole machine reserved (an
indexer of files used about a third of one core), one repetition. Table:
`summary.md`. The references on the CPU come from
`../solvers-cpu-2026-10-15` and `../solvers-2026-10-15-byvdim`.

| problem | configuration | second operator (s) | iterations | residual |
|---|---|---|---|---|
| elasticity | CG of Belos + MueLu | 0.124 | 14 | 8.9e-9 |
| elasticity | CG of MFEM + MueLu, Jacobi smoother, hierarchy reused | 0.106 | 14 | 1.3e-6 |
| elasticity | CG + GAMG (PETSc), true residual | 0.119 | 21 | 2.8e-9 |
| elasticity | CG + BoomerAMG | 0.199 | 30 | 1.2e-7 |
| elasticity | CG of Belos + Jacobi relaxation (Ifpack2) | 0.082 | 319 | 9.0e-9 |
| heterogeneous | CG of Belos + MueLu | 0.155 | 27 | 5.4e-9 |
| heterogeneous | CG + GAMG, true residual | 0.134 | 39 | 5.5e-9 |
| heterogeneous | CG + BoomerAMG | 0.272 | 47 | 2.5e-7 |
| laplace | CG of Belos + MueLu | 0.108 | 9 | 4.3e-9 |
| laplace | CG + BoomerAMG | 0.121 | 13 | 3.7e-8 |
| laplace | CG of Belos + Jacobi relaxation | 0.057 | 91 | 7.8e-9 |

- On 8 cores, MueLu is as fast as GAMG of PETSc on the elasticity (0.124
  against 0.119 seconds per operator for the same true residual) and
  faster than BoomerAMG; with the inclusions, GAMG stays slightly faster.
- At this size, the Jacobi relaxation remains the fastest on the
  homogeneous beam and on the Laplace equation, as in the first pass.
- FROSch: 6.9 seconds per operator on the elasticity with its default
  subdomain solver (KLU2), 4.8 with Tacho for the subdomains, the
  extensions and the coarse problem (on the GPU). On the Laplace equation,
  KLU2 did not finish in 150 seconds and Tacho ran out of memory on the
  GPU with 8 processes (`errors.txt`).
- Tacho alone (on the GPU shared by the 8 processes) factorizes the beam
  in 0.96 seconds and the Laplace equation in 5.8 seconds.
