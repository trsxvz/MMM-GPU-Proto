# Matrix free tangent on the beam (10/10/2026)

Beam r4 (111,843 unknowns, tetrahedra of order 1), von Mises plasticity,
10 steps, elastic prediction. mfem-mgis 1bbcb877d8. Script:
`cases/measure-L4.sh`. One RTX 5080, 16 cores, quiet machine.

Two ways to apply the tangent, same kernels for the gradients, the
residual and the behaviour (`generic` behaviour on CPU, `cuda` on the
GPU):

- **assembled**: sparse matrix assembled by kernels
  (`TangentOperatorComputation=Kernels`), PCG + BoomerAMG of hypre
  (`-ls pcg-amg`);
- **matrix free**: action `Σ Bᵀ K B` by kernels
  (`TangentOperatorComputation=MatrixFree`), conjugate gradient of MFEM,
  BoomerAMG preconditioner built on the assembled tangent, rebuilt every 1,
  2 or 4 Newton iterations (`-ls cg-amg -pu n`, preconditioner
  `AssembledOperator`).

Same relative tolerance of the linear solver (1e-4) everywhere.

| Case | Total | Tangent | Preconditioner | Solve | Prediction |
|---|---|---|---|---|---|
| CPU 8 processes, assembled | 19.15 | 2.60 (29) | 0.01 | 12.35 | 3.75 |
| CPU 8 processes, matrix free, pu 1 | 27.52 | 0.00 | 2.63 | 20.79 | 3.65 |
| CPU 8 processes, matrix free, pu 2 | 20.84 | 0.00 | 1.26 | 16.19 | 3.00 |
| CPU 8 processes, matrix free, pu 4 | 25.28 | 0.00 | 0.75 | 21.00 | 3.13 |
| GPU, assembled | **4.08** | 0.24 (32) | 0.01 | 2.40 | 0.94 |
| GPU, matrix free, pu 1 | 7.38 | 0.00 | 0.23 | 5.69 | 0.98 |
| GPU, matrix free, pu 2 | 7.27 | 0.00 | 0.18 | 5.54 | 1.05 |
| GPU, matrix free, pu 4 | 6.83 | 0.00 | 0.11 | 5.39 | 0.85 |

Times in seconds (maximum over the processes), number of calls in
parentheses. "Preconditioner": assembly of the tangent for BoomerAMG (the
setup of BoomerAMG is counted in the solve). Same solution everywhere
(norm 48.83704586), 29 to 32 Newton iterations.

## Findings

- **The matrix free tangent is correct, but slower on this case**: 6.8 s
  against 4.1 s on the GPU, 20.8 s against 19.2 s on 8 cores at best. The
  gain on the assembly (0.24 s) does not make up for the cost of the action
  in the conjugate gradient.
- **Keeping BoomerAMG for several iterations** hardly helps: the solve
  goes from 5.7 s to 5.4 s on the GPU (pu 1 to 4). On CPU, pu 4 is slower
  than pu 2 (the conjugate gradient does more iterations with an old
  preconditioner).

## Nsight Systems profile (GPU, pu 1, 10 steps)

`nsys profile -t cuda,nvtx`, compute time of the GPU:

| Variant | Launches | GPU time | Total time |
|---|---|---|---|
| assembled | 318,768 | 2.62 s | 4.0 s |
| matrix free | 312,315 | 5.15 s | 6.4 s |

- The assembled variant is limited by the latency of the launches (2.6 s
  of compute out of 4 s), the matrix free variant by the compute.
- One action of the matrix free tangent costs about 1.3 ms on the GPU
  (1,788 actions), against 0.1 to 0.2 ms for a product by the assembled
  matrix. Detail per action: derivatives of the unknowns at the points by
  `QuadratureInterpolator::Derivatives` of MFEM, 0.89 ms (generic version
  for simplices, without tensor product: this is the main item); internal
  forces, 0.22 ms; gradients from the derivatives, 0.07 ms; product by the
  blocks of the tangent, 0.04 ms; restriction to the elements, 0.05 ms.
- Leads: a kernel fused per element (reading of the unknowns of the
  element, gradient at the points with the gradients of the shape
  functions and the inverse of the jacobian, product by K, `Bᵀ`,
  accumulation), which removes the intermediate arrays and the pass through
  the interpolator of MFEM. The same kernel serves for the residual. For
  tetrahedra of order 1, the data volume of a fused action (blocks of K at
  the points) stays of the same order as that of the assembled matrix: the
  benefit of the matrix free tangent will come mostly from order 2 and from
  large meshes (memory of the matrix).
