# Beam on GPU, minimal slice (09/10/2026)

Beam r4 (111,843 unknowns, tetrahedra of order 1), von Mises plasticity,
PCG + BoomerAMG (tolerance 1e-4), 10 steps. mfem-mgis 7159b6de4c, MGIS
a0b5140. Script: `cases/measure-L6.sh beam`. One RTX 5080, 16 cores.

Computations on GPU: `cuda` device of MFEM, hypre on the GPU, all the
kernels (gradients, internal forces, tangent assembled directly into the
sparse matrix, prediction), `generic` behaviour (integrated on the host)
or `generic-parallel` compiled for `cuda` (integrated on the GPU, on the
copies of the arrays of the material).

Caveat: the machine was coming out of a compilation (average load of 10 at
the start, decreasing). The point by point reference computation on 8
processes takes 17.5 s against 15.7 s on a quiet machine.

| Case | Total | Integration | Residual | Tangent | Solve | Prediction |
|---|---|---|---|---|---|---|
| CPU, 8 processes, point by point | 17.54 | 0.58 (41) | 0.48 (41) | 2.38 (31) | 10.89 | 3.12 (10) |
| CPU, 8 processes, everything by kernels | 14.68 | 0.20 (39) | 0.11 (39) | 1.91 (29) | 9.50 | 2.89 (10) |
| GPU, 1 process, `generic` behaviour | 5.01 | 0.79 (42) | 0.05 (42) | 0.33 (32) | 2.55 | 0.99 (10) |
| GPU, 1 process, `cuda` behaviour | **4.06** | 0.10 (42) | 0.03 (42) | 0.24 (32) | 2.47 | 0.87 (10) |
| GPU, 2 processes, `cuda` behaviour | 28.64 | 0.11 (41) | 0.04 (41) | 0.44 (31) | 22.90 | 4.91 (10) |

Times in seconds (maximum over the processes); in parentheses, number of
calls. Same solution everywhere (norm 48.83704586). The numbers of Newton
iterations vary from 29 to 32 with the rounding (order of the sums, atomic
additions on GPU).

## Findings

- **One GPU against 8 cores**: 4.06 s against 17.54 s (×4.3), and against
  14.68 s for the same kernels on CPU (×3.6).
- **What remains**: the linear solve (PCG + BoomerAMG of hypre on the GPU)
  takes 61 % of the time, the prediction 21 % (including its own linear
  solve). The tangent costs 7.5 ms per call, the integration of the
  behaviour 2.4 ms per call.
- **Behaviour on the GPU**: the integration goes from 0.79 s to 0.10 s (no
  more copies of the arrays between the host and the GPU).
- **Two processes on the same GPU**: the solve becomes 9 times slower
  (22.9 s). Not analysed. Most probable hypothesis: without the MPS service
  of CUDA, the contexts of the two processes share the GPU by time slices,
  and the many small kernels of BoomerAMG run one after the other with
  context switches. To be checked with MPS (`nvidia-cuda-mps-control -d`).
  One process per GPU stays the natural configuration.
- **Memory**: each process on the GPU reserves a pool of 4 GiB (Umpire).
  Four processes on the 16 GB GPU fail (seen when running the tests in
  parallel).

## Nsight Systems profile (3 steps, 1 process, `cuda` behaviour, mfem-mgis f0239ab807)

`nsys profile -t cuda,nvtx` (Nsight Systems 2025.5 of the HPC SDK 25.11;
the `nsys` of CUDA 13.0 does not read this report) on the beam r4, 3 steps
(10 Newton iterations, 13 residuals), about 1.8 s:

| Copies | Count | Volume | Time |
|---|---|---|---|
| host to GPU | 835 | 146 MB | 7.3 ms |
| GPU to host | 5,814 | 9 MB | 2.2 ms |
| within the GPU | 14,704 | 3.0 GB | 13.7 ms |

- **No notable hidden copy**: the large copies to the GPU are made only
  once or a few times (geometry, pattern, mesh). Per iteration, one or two
  vectors (0.9 MB) remain in each direction (norm and imposed dofs of the
  prediction, Dirichlet conditions) and scalars (norms, statuses of hypre:
  3,215 copies of 8 bytes).
- **Compute time of the GPU: 0.71 s out of 1.8 s**, in 95,339 kernel
  launches. The sparse matrix-vector products of hypre (`csrmv`, PCG and
  BoomerAMG) take 65 % of it (460 ms), the matrix product of the setup of
  the AMG 78 ms, the kernels of MFEM and of mfem-mgis (`CuKernel1D`,
  `CuKernel2D`) 69 ms, the behaviour (`apply`) 16 ms.
- The rest of the time (about 1.1 s) is spent on the host launching and
  waiting for these small kernels: the solve is limited by the latency of
  the launches of hypre, not by the compute. Leads: fewer PCG iterations
  (preconditioner reused from one iteration to the other, options of
  BoomerAMG for the GPU), larger problems, matrix free operator.

## Matrix free tangent (mfem-mgis 4d465460c3)

Option `TangentOperatorComputation=MatrixFree`: the tangent is applied by
kernels without assembly (`Σ Bᵀ K B` at the points), and so is its
diagonal; conjugate gradient of MFEM with a Jacobi preconditioner (`-ls
cg-jacobi`), same relative tolerance (1e-4) as PCG + BoomerAMG. `cuda`
behaviour, 10 steps, 1 process:

| Case | Total | Tangent | Solve | Prediction |
|---|---|---|---|---|
| assembled tangent, PCG + BoomerAMG | 3.61 | 0.22 (32) | 2.23 | 0.77 |
| matrix free tangent, CG + Jacobi | 95.31 | 0.00 (30) | 93.30 | 1.58 |

Same solution (48.83704586), same numbers of Newton iterations up to the
rounding. The matrix free path is correct, but Jacobi is not enough: 3.1 s
of conjugate gradient per Newton iteration (thousands of iterations on
this slender beam). The matrix free operator needs a real preconditioner:
AMG built on the assembled tangent (or on the elastic operator) and kept
for several iterations, or polynomial multigrid with Chebyshev smoothers
(work packages L4 and L5 of the plan).
