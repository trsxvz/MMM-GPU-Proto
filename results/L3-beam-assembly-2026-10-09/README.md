# Direct assembly of the tangent by kernels, beam (09/10/2026)

Intermediate measurements of work package L3: the tangent is no longer made of
element matrices stored and then copied by the loop of MFEM. A kernel
computes each `N x N` block of a pair of nodes of an element and adds it
directly into the local sparse matrix (position searched in the row, atomic
addition). The pattern of the matrix is the one of the first assembly of
MFEM. The end of `GetGradient` is the one of MFEM (diagonal block and
`PtAP` in parallel, `RAP` in sequential).

Beam r4 (111,843 unknowns, tetrahedra of order 1), PCG + BoomerAMG
(tolerance 1e-4). Script: `cases/measure-L3.sh beam` (10 steps, on CPU).
The directory `trials-3-step` contains computations of 3 steps with 1
process on CPU and on GPU (`--device cuda`), with provisional timers (lines
`TMP`).

Caveat: the measurements over 10 steps were made with a web browser open
(load of 1 to 2). The total times vary by 10 to 15 % between identical
computations. The summary gives the base commit (cf7f764628), since the
changes are not committed yet.

## On CPU (10 steps)

| Case | Total | Integration | Residual | Tangent | Solve |
|---|---|---|---|---|---|
| 8 processes, batches | 17.65 | 0.55 (41) | 0.48 (41) | 2.36 (31) | 10.99 |
| 8 processes, tangent by kernels | 20.13 | 0.71 (39) | 0.64 (39) | 3.25 (29) | 11.70 |
| 8 processes, everything by kernels | 18.46 | 0.22 (39) | 0.13 (39) | 3.04 (29) | 11.28 |
| 1 process, batches | 66.61 | 3.22 (38) | 3.00 (38) | 11.72 (28) | 34.50 |
| 1 process, tangent by kernels | 70.29 | 3.26 (38) | 3.04 (38) | 14.27 (28) | 35.35 |
| 1 process, everything by kernels | 67.11 | 0.96 (38) | 0.65 (38) | 14.31 (28) | 36.23 |

The direct assembly gains nothing on CPU: the tangent by kernels stays
22 % slower than element by element, as with the stored element matrices.
The timers show why: the kernel alone takes 0.48 s per call, against
0.34 s for the whole loop of MFEM (element matrices and `AddSubMatrix`).
One thread per pair of nodes recomputes the inverse of the jacobian and
reads the blocks of the tangent again for each pair: a good granularity
for a GPU, a bad one for a CPU core.

The residual by kernels is divided by 4.6 with 1 process (loop of MFEM
skipped), the gradients by 3.4.

## On GPU (3 steps, 1 process)

| Case | Total | Integration | Residual | Tangent | Solve | Prediction |
|---|---|---|---|---|---|---|
| CPU, batches | 24.53 | 1.14 (13) | 1.05 (13) | 4.47 (10) | 13.36 | 4.32 |
| GPU, batches | 29.05 | 6.70 (13) | 6.91 (13) | 9.32 (10) | 0.89 | 4.96 |
| GPU, everything by kernels, `EliminateRowsCols` | 9.05 | 0.47 (13) | 0.02 (13) | 2.38 (10) | 0.83 | 5.08 |
| GPU, everything by kernels, `EliminateBC` | 7.84 | 0.49 (14) | 0.02 (14) | 1.26 (11) | 0.90 | 4.92 |

- hypre is compiled for the GPU: the solve (PCG + BoomerAMG) goes from
  13.4 s to 0.9 s.
- Without kernels, the element loops of MFEM with vectors on the GPU are
  slower than on CPU (copies).
- The assembly kernel takes 7 ms per call (0.48 s on a CPU core).
- Timers of the end of `GetGradient` (`cuda-np1-batches-kernels-timers.txt`):
  diagonal block 1 ms, `PtAP` 1.3 ms, elimination of the imposed dofs
  122 ms. `HypreParMatrix::EliminateRowsCols` copies the matrix back to the
  host, eliminates while computing the eliminated part (unused), then sends
  it back to the GPU. `HypreParMatrix::EliminateBC` (`DIAG_ONE`) gives the
  same matrix, since the pattern is symmetric, in place and on the GPU: the
  tangent goes from 2.38 s to 1.26 s, of which about 0.9 s for the first
  call, which still goes through the loop of MFEM to build the pattern.
- The number of iterations varies from one computation to the other on GPU
  (10 or 11): atomic additions and reductions of the GPU, non deterministic
  order of the sums. Same norm of the solution up to 1e-10.

## Pattern built from the elements (mfem-mgis 799941c8fd)

The pattern of the local matrix is built from the dofs of the elements at
the first call: no more loop of MFEM. In each row, the diagonal term comes
first, then the other columns sorted: this is the order that
`HypreParMatrix` imposes in place on the rows of the local matrix it shares
(`hypre_CSRMatrixReorder`). The kernel searches the columns by bisection
instead of a scan of the row.

| Case (3 steps, 1 process) | Total | Integration | Residual | Tangent | Solve | Prediction |
|---|---|---|---|---|---|---|
| CPU, batches | 24.53 | 1.14 (13) | 1.05 (13) | 4.47 (10) | 13.36 | 4.32 |
| CPU, everything by kernels | 22.62 | 0.37 (13) | 0.22 (13) | 4.00 (10) | 13.32 | 4.53 |
| GPU, everything by kernels | 6.84 | 0.50 (13) | 0.02 (13) | 0.16 (10) | 0.89 | 4.99 |

- On CPU, the tangent by kernels becomes faster than the loop of MFEM
  (0.40 s per call against 0.45 s): the cost came mostly from the linear
  scan of the rows, not from the recomputation per pair of nodes.
- On GPU, 16 ms per tangent (against 126 ms with the first pattern of
  MFEM). The computation is dominated by the prediction (73 %), which
  assembles bilinear forms element by element on the host.

## Next steps

- prediction by kernels (work package L4);
- behaviour on the GPU (`generic-parallel` CUDA library);
- MPI test of the kernels in mfem-mgis;
- variant of the kernel with one thread per element, to compare on CPU.

## Behaviour integrated on the GPU (mfem-mgis ee1f31fa1c, MGIS a0b5140)

`generic-parallel` behaviour compiled for `cuda`, integrated by MGIS on the
copies of the arrays of the material on the GPU (view of MGIS on the copies
of MFEM), properties and external state variables copied to the GPU at each
call; all the kernels. Beam r4, 3 steps, 1 process:

| Case | Total | Integration | Residual | Tangent | Solve | Prediction |
|---|---|---|---|---|---|---|
| GPU, `generic` behaviour (host) | 6.51 | 0.47 (13) | 0.02 (13) | 0.15 (10) | 0.76 | 4.85 |
| GPU, `cuda` behaviour | 6.25 | 0.24 (13) | 0.01 (13) | 0.12 (10) | 0.77 | 4.80 |

Same solution. The prediction remains (77 % of the time), which assembles
its operators element by element on the host.
