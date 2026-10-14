# Notched bar ssna303 (11/10/2026)

Measurement campaign, case 2. Example ex4 of mfem-mgis-examples, variant
solved by hypre (Newton, FGMRES of hypre preconditioned by BoomerAMG),
finite strain plasticity (Hencky logarithmic strain), complete loading in
10 steps, mesh refined zero, one and two times. Script:
`cases/measure-ssna303.sh`. mfem-mgis 8a6325335e installed (FBar by
kernels, before the fixes of reviews F and G; `summary.txt` gives the
version of the sources, f0324f1d9c). One RTX 5080 (16 GB), 8 cores (16
threads), 30 GB of memory. Quiet machine during these computations.

INCOMPLETE series: stopped at 7:18 during the computation
`ssna303-r2-fbar-np8-kernels`, the memory of the machine being almost
exhausted (6 GiB of swap used afterwards). This computation and
`ssna303-r2-fbar-gpu-kernels` are missing.

| Refinement | Unknowns | 8 cores, point by point | 8 cores, kernels | GPU, kernels | Gain of the GPU |
|---|---|---|---|---|---|
| 0 | 27,945 | 5.97 s | 6.08 s | 6.10 s | ×1.0 |
| 1 | 196,047 | 66.8 s | 65.5 s | 27.1 s | ×2.4 |
| 2 | 1,461,915 | 683 s | 670 s | failure | |

FBar (`--use-fbar`, by kernels since the work package L7.9):

| Refinement | 8 cores, point by point | 8 cores, kernels | GPU, kernels | Gain of the GPU |
|---|---|---|---|---|
| 0 | 7.40 s | 6.97 s | 6.34 s | ×1.1 |
| 1 | 85.5 s | 78.7 s | 29.5 s | ×2.7 |
| 2 | 1420 s | stopped | not launched | |

Time of the loop over the steps; the gain compares the GPU to the kernels
on 8 cores.

- On the GPU, the solution (FGMRES and BoomerAMG of hypre) dominates:
  4.5 s out of 6.1 s at zero refinement, 18 s out of 27 s at one
  refinement. The gain comes mostly from the integration (1.2 s against
  11.2 s), from the tangent (3.0 s against 10.6 s) and from the
  prediction.
- On 8 cores, the kernels gain on the integration and the residual, but
  the tangent by kernels is slower than the element by element tangent
  (10.6 s against 8.4 s at one refinement), as already noted for the
  thermal (search of the columns in the sparse matrix).
- FBar costs 20% more than the standard formulation on 8 cores, 9% on the
  GPU, at one refinement.
- Two refinements, FBar, point by point: 1420 s, with a split step (48
  solutions against 25 without FBar) and perhaps memory swapping: a time
  not to be compared.
- Two refinements on the GPU, without FBar: failure from the first step on
  (10 splits, then stop), with no message in the output. ex4 has since
  printed the cause of the failure of each solution (mfem-mgis-examples
  db0ce2a). Relaunched alone at 7:35: the memory of the GPU is exhausted
  (at most 15.2 GiB out of 16); the allocation of 880 MB of the Jacobians
  of the integration rule (first computation of the geometric data, kept
  afterwards) fails at the first computation of the tangent. With 27
  points per hexahedron (rule of mfem-mgis for these elements), the mesh
  has about 12 million points: the blocks of the tangent (81 reals per
  point) alone take nearly 7.9 GiB (estimate from the size of the
  allocation, not measured). This refinement does not fit on this GPU
  without reducing the memory at the points (8-point rule, blocks of the
  tangent in single precision or computed on the fly).

Breakdown at one refinement (`cases/summary-profiling.py`, seconds, total
with the setup):

| Computation | Total | Integration | Tangent | Solution | Prediction |
|---|---|---|---|---|---|
| 8 cores, point by point | 66.8 | 12.6 | 8.4 | 29.4 | 14.5 |
| 8 cores, kernels | 65.5 | 11.2 | 10.6 | 28.6 | 14.3 |
| GPU, kernels | 27.2 | 1.2 | 3.0 | 18.1 | 3.6 |
| FBar, GPU, kernels | 29.6 | 1.2 | 3.2 | 20.0 | 3.8 |
