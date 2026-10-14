# Reuse of the analysis of MUMPS (12/10/2026)

Parameter `ReuseAnalysis` of the `MUMPSSolver` solver of mfem-mgis (branch
mumps-reuse, c020804d68 and a114b6c870 after review): the analysis of MUMPS
(ordering and symbolic factorisation) of the previous matrix is reused
when the new matrix has the same pattern on all the processes. The
matrices of the successive Newton iterations have the same pattern: only
the first one is analysed. `mfem::MUMPSSolver::SetReorderingReuse` (MFEM
4.10) does not check the pattern: mfem-mgis compares it (rows, columns of
the diagonal and off-diagonal blocks, global indices of the off-diagonal
columns), and redoes the analysis on all the processes if the pattern has
changed on one of them. The function `getNumberOfMUMPSAnalyses` gives the
number of analyses done. Since the matrix is distributed, the analysis of
MUMPS only depends on the pattern: reusing it does not change the
factorisation. Option `-mra` of ex8 (mfem-mgis-examples 2720017).

## Elasticity matrix of ex8 alone

Program `program/mumps-reuse.cpp` (MFEM alone): linear elasticity on the
mesh of ex8 (28,548 unknowns, linear elements, ordering by nodes, all the
boundaries clamped), symmetric positive definite MUMPS, six factorisations
of matrices with the same pattern and different values. Mean time of the
factorisations after the first one, in seconds
(`variants-elasticity-matrix.txt`):

| Variant | Without reuse | With reuse |
|---|---|---|
| 1 process, automatic ordering | 0.258 | 0.207 |
| 1 process, METIS, Scotch or PORD | 0.255 to 0.256 | 0.205 to 0.206 |
| 1 process, AMF | 0.380 | 0.359 |
| 1 process, AMD | 0.496 | 0.476 |
| 1 process, BLR at 1e-10, 1e-8 or 1e-6 | 0.277 to 0.291 | 0.226 to 0.238 |
| 2 processes, automatic ordering | 0.189 | 0.136 |

- The reuse saves 20% on one process and 28% on two.
- The automatic ordering of MUMPS is the best; the approximate BLR
  factorisation is slower at this size.
- The parallel orderings (ParMETIS, PT-Scotch) are not in the MUMPS of the
  environment: MUMPS returns the error -38 and MFEM stops.

Only the reuse of the analysis is therefore offered in mfem-mgis.

## Complete ex8

5 steps, original mesh, everything by kernels, symmetric positive definite
MUMPS for the mechanics (`-msMc SymmetricPositiveDefinite`), GMRES and
BoomerAMG for the thermal. Build `build/examples-mumps`, linked to
`install/mfem-mgis-mumps`. Two alternated runs per variant
(`summary-np2.txt`, `summary-gpu.txt`). Times in seconds:

| Case | Variant | Total | Coupling | Factorisations (29) |
|---|---|---|---|---|
| CPU, 2 processes (cores 0 to 3), generic behaviours | without | 20.42 and 20.63 | 20.30 and 20.52 | 5.47 and 5.57 |
| | with `-mra` | 19.15 and 19.21 | 19.03 and 19.08 | 4.08 and 4.07 |
| GPU, 1 process (cores 4 to 7), cuda behaviours, hypre on the host | without | 11.15 and 11.08 | 10.96 and 10.90 | 7.23 and 7.22 |
| | with `-mra` | 9.73 and 9.73 | 9.55 and 9.55 | 5.87 and 5.86 |

- The factorisations take 26% less on CPU (2 processes) and 19% less on
  the GPU (one process, MUMPS on the host). ex8 on the GPU goes from 11.1
  to 9.7 s.
- Same statistics of the fields (temperature, norm of the displacement)
  and same numbers of Newton iterations (35 in thermal, 29 in mechanics)
  with and without reuse. The residuals differ by roundoff, as much as
  between two identical runs without reuse.
- On the GPU, MUMPS remains the first item (60% of the time).

## Tests

`MUMPSReuseAnalysisTest` (1 and 2 processes, and on the GPU): sequence of
matrices whose pattern is kept, changes on all the processes, changes on
the last process only, changes only through the columns of the
off-diagonal blocks, or only through the global indices of these columns,
with and without reuse, symmetric or not. After each matrix, the test
checks the residual and the number of analyses: a new analysis if and only
if the pattern has changed. Seven mutations detected (no comparison, no
reduction over the processes, global indices or off-diagonal columns not
compared, never any reuse, pattern forgotten after a reuse, analyses not
counted). Independent review Z8 handled.

## Files

- `ex8-np2-kernels-{without,with}-{1,2}.txt`, `ex8-gpu-kernels-{without,with}-{1,2}.txt`:
  raw outputs with the profiling tables.
- `measure-ex8-cpu.sh`, `measure-ex8-gpu.sh`: commands of the measurements.
- `program/`: measurement program on the elasticity matrix (`make` with
  `scripts/env.sh`, then
  `./mumps-reuse assemblage_hexa.msh 6 [ordering] [BLR tolerance]`).
