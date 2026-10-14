# cuDSS, direct solver on the GPU (12/10/2026)

Solver `CuDSSSolver` of the factory of linear solvers of mfem-mgis
(branch cudss, 840683d96f, and 45afaa0653 after review), which wraps
`mfem::CuDSSSolver` (MFEM 4.10, cuDSS 0.7.1). Parameters `Symmetric`,
`PositiveDefinite` (the symmetric modes read the lower triangle, like
MUMPS) and `ReuseAnalysis` (analysis reused as long as the pattern does
not change on any process, as for MUMPS,
`results/mumps-analysis-2026-10-12`). The solver is refused without a
device on a CUDA GPU. ex8 selects it by `-lsMc CuDSSSolver`
(mfem-mgis-examples efeb8b8).

## Elasticity matrix of ex8 alone

Program `program/cudss-bench.cpp` (MFEM alone): linear elasticity on the
mesh of ex8 (28,548 unknowns, linear elements, ordering by nodes, all the
boundaries clamped), one process, six matrices with the same pattern and
different values. Average time of the factorisations after the first one,
in seconds (`elasticity-matrix-np1-hypre-gpu.txt`):

| Solver | Without reuse | With reuse |
|---|---|---|
| MUMPS on the host, symmetric positive definite | 0.253 | 0.203 |
| cuDSS, general matrix | 0.126 | 0.031 |
| cuDSS, symmetric positive definite (lower triangle) | 0.111 | 0.018 |

- With the reuse of the analysis, cuDSS factorises 11 times faster than
  MUMPS on the host.
- Relative residuals of the order of 6e-15 in all cases.
- With hypre on the host (`elasticity-matrix-np1-hypre-host.txt`), same
  orders of magnitude (0.022 s for symmetric cuDSS with reuse): MFEM
  copies the matrix to cuDSS from the memory of the host.
- This program hangs on two processes even before the construction of the
  finite element space, with the default configuration of MFEM (memory of
  the host through Umpire): not analysed, unrelated to cuDSS. The tests of
  mfem-mgis on two processes and the GPU pass.

## Full ex8 on the GPU

5 steps, original mesh, one process (cores 4 to 7), everything by
kernels, cuda behaviours, hypre on the host (thermal problem by GMRES and
BoomerAMG), symmetric positive definite mechanics. Build
`build/examples-mumps`, linked to `install/mfem-mgis-mumps` (branch
cudss). Two alternating passes (`summary.txt`, commands in
`measure-ex8-gpu.sh`). Times in seconds:

| Mechanics | Total | Coupling | Factorisations (29) |
|---|---|---|---|
| MUMPS on the host, analysis reused | 9.80 and 9.73 | 9.60 and 9.55 | 5.87 and 5.83 |
| cuDSS, without reuse | 7.12 and 7.14 | 6.93 and 6.95 | 3.47 and 3.47 |
| cuDSS, analysis reused | 4.37 and 4.38 | 4.19 and 4.19 | 0.69 and 0.70 |

- ex8 on the GPU goes from 11.1 s (MUMPS without reuse, before this work
  package and the previous one) to 4.4 s, that is ×2.5.
- After the review (45afaa0653: solutions checked by their residual,
  pattern compared without the values; `summary-review.txt`, commands in
  `measure-ex8-gpu-review.sh`, cores 0 to 3): 4.52 and 4.41 s with hypre
  on the host (the check costs at most 0.1 s); 5.63 and 5.58 s with hypre
  on the GPU (slower thermal problem, factorisations 0.63 s, the
  comparison of the patterns now only reads the indices); 8.38 and 8.45 s
  with hypre on the GPU without reuse. Same fields.
- Same statistics of the fields (temperature, norm of the displacement)
  and same numbers of iterations (114) in the three cases.
- The factorisations now weigh only 16 % of the time.

## ex8, thermal problem also by cuDSS

With the two direct solvers on the GPU (`-lsTh CuDSSSolver -lsMc
CuDSSSolver -mra`, thermal matrices given whole, mechanical ones by their
lower triangle; `summary-thermal.txt`, commands in
`measure-ex8-gpu-thermal.sh`, cores 0 to 3), in seconds:

| hypre | Total | Coupling | Thermal factorisations (35) | Thermal solves (35) | Mechanical factorisations (29) |
|---|---|---|---|---|---|
| on the host | 2.78 and 2.77 | 2.60 and 2.59 | 0.22 and 0.21 | 0.036 and 0.034 | 0.68 and 0.68 |
| on the GPU | 2.58 and 2.73 | 2.40 and 2.53 | 0.20 and 0.21 | 0.036 and 0.036 | 0.62 and 0.62 |

- The thermal solves by GMRES and BoomerAMG took 1.83 s (hypre on the
  host): cuDSS factorises the thermal matrix (9,516 unknowns) in 6 ms and
  solves it in 1 ms.
- ex8 on the GPU goes from 11.1 s this morning (MUMPS without reuse,
  GMRES and BoomerAMG) to 2.6 s, that is ×4.3, same results and same
  iterations (114).
- hypre now only carries the parallel matrices and vectors: keeping it on
  the GPU avoids the copies between the host and the GPU of the matrices
  given to cuDSS.

## ex7 with the direct solvers

Full ex7 (40 steps, 2,634 unknowns, periodic cell), everything by
kernels, symmetric positive definite matrices, analysis reused (options
`-ls`, `-ms` and `-mra` added to ex7, mfem-mgis-examples e98e259). Build
`build/examples-mumps/ex7`, commands in `measure-ex7.sh`, total durations
in `durations-ex7.txt`, in seconds:

| Case | GMRES and BoomerAMG | Direct solver |
|---|---|---|
| CPU, 2 processes (cores 0 to 3), generic behaviour | 6.17 and 6.08 | MUMPS: 5.41 and 5.36 |
| GPU, 1 process, cuda behaviour, hypre on the host | 4.61 and 4.32 | cuDSS: 2.08 and 2.02 |
| GPU, 1 process, cuda behaviour, hypre on the GPU | | cuDSS: 2.35 and 2.37 |

- On the GPU, cuDSS divides the time of ex7 by 2.1. The GPU is 2.6 times
  faster than the CPU with two cores.
- Same iterations (149) and same comparison with the reference of the
  average stresses in all the cases.
- Here, keeping hypre on the host stays faster (2.0 s against 2.4 s),
  unlike ex8 (cause not analysed).

## Limits of distributed cuDSS

With several processes on one GPU, the distributed analysis of cuDSS
0.7.1 fails (`CUDSS_STATUS_EXECUTION_FAILED`, code 5, or segmentation
fault) for some matrices that MUMPS solves:

- when rows only have their diagonal coefficient (isolated vertices of
  the graph of the matrix): from 2 such rows when the processes are
  coupled; without coupling, the threshold depends on the pattern;
- for a tridiagonal matrix coupling the processes, from 4 rows per
  process;
- for some uneven distributions of the rows of a matrix where each row is
  coupled to a row of another process (4, 2 and 6 rows on three
  processes; 4, 3 and 5 pass).

One process is not affected. Standalone reproducer:
`reproducers/cudss-distributed`. The diffusion, matching and coupling
matrices of the test pass. The tangent matrices of mfem-mgis keep the
coefficients of the rows of the imposed degrees of freedom (set to zero,
not removed). The test of mfem-mgis skips the matrices with diagonal rows
for cuDSS on several processes, and moves a single row between processes.

## Tests

`DirectSolverReuseAnalysisTest`, 10 variants: MUMPS on 1, 2 and 3
processes and on the GPU, cuDSS refused without GPU, cuDSS on the GPU on
1, 2 and 3 processes and with hypre on the host on 1 and 2 processes;
10/10. A negative definite matrix declared positive definite gives with
cuDSS a finite and wrong solution, rejected by the check of the residual.
Mutations: six detected (collective decision short-circuited, pattern not
compared, residual not checked, pattern forgotten after a reuse, buffers
not copied, global indices of the off-diagonal columns not compared). The
local decision without reuse is not detected, and it changes nothing here:
the compared layout includes the global number of rows and the first row,
and a process that keeps its own keeps valid descriptors. Independent
review Z9 addressed.
