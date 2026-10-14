# Local assembly and global matrix (11/10/2026)

Bench `cases/assembly` in parallel (`--parallel`): assembly by kernels of
the stiffness matrices into the local matrix of each process, then
construction of the global matrix as mfem-mgis does
(`assembleTangentOperatorByKernels`): block diagonal matrix, product
`Pᵀ A P` of hypre, elimination of the degrees of freedom of a clamped end.
mfem-mgis 4a007e85f2 (separate installation `install/mfem-mgis-registers`).
Median times of 10 assemblies, maxima over the processes. Outputs:
`cpu-np8.txt`, `gpu.txt`.

| Case | Local assembly | Pᵀ A P | Elimination |
|---|---|---|---|
| 8 processes, beam r4 | 9.5 ms | 15.4 ms | 0.3 ms |
| 8 processes, thermal r5 | 14.8 ms | 36.5 ms | 0.2 ms |
| 8 processes, finite strain, Q1 hexahedra r4 | 6.5 ms | 3.1 ms | 0.2 ms |
| GPU, 1 process, beam r4 | 2.3 ms | 1.1 ms | 0.2 ms |
| GPU, 1 process, thermal r5 | 13.8 ms | 1.5 ms | 0.04 ms |
| GPU, 1 process, finite strain r4 | 4.5 ms | 0.4 ms | 0.5 ms |

- On 8 processes, the product `Pᵀ A P` now costs more than the local
  assembly for the linear elements: hypre recomputes the pattern and the
  communications at each tangent, while they do not change. Lead: pattern
  of the global matrix and exchange scheme kept, only the values being
  sent (item of plan L4).
- On the GPU with one process, this product is small.
- Nsight Compute (one launch of the kernel of the matrices): thermal and
  mechanics bound by the computation (85 to 89% of the throughput of the
  multiprocessors, memory almost idle), hence the next work package
  (inverses of the Jacobians computed once,
  `results/jacobians-inverses-2026-10-11`).
