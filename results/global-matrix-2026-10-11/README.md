# Global matrix with pattern and exchanges kept (11/10/2026)

Bench `cases/assembly --parallel` on 8 processes: local assembly by kernels
(with the positions), then global matrix by the product `Pᵀ A P` of hypre
(as before) and by `batch::GlobalStiffnessMatrix` (mfem-mgis branch
`global`, 28f2367b39, separate installation `install/mfem-mgis-global`):
pattern of the global matrix and exchanges between processes computed at
the first assembly, then only the values copied, exchanged
(`MPI_Alltoallv`) and added. Median times of 10 assemblies, maxima over
the processes. Quiet machine. Output: `cpu-np8.txt`.

| Case | Local assembly | Pᵀ A P of hypre | Pattern kept (first assembly) | Gain |
|---|---|---|---|---|
| beam, P1 tetrahedra (r4) | 9.0 ms | 15.4 ms | 2.5 ms (32 ms) | ×6.3 |
| thermal, P1 tetrahedra (r5) | 14.1 ms | 36.8 ms | 1.7 ms (46 ms) | ×21.6 |
| finite strain, Q1 hexahedra (r4) | 6.4 ms | 3.0 ms | 0.15 ms (7 ms) | ×20 |
| P2 tetrahedra, mechanics (r3) | 15.2 ms | 22.0 ms | 5.8 ms (65 ms) | ×3.8 |

Same matrices (relative gap of the products by a vector at most 1.4e-16).
The first assembly costs two to three products of hypre: it pays off from
the second tangent on. The global matrix now costs only a fraction of the
local assembly.
