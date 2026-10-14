# Positions of the entries in the sparse matrix (11/10/2026)

Benchmark `cases/assembly`: assembly by kernels of the stiffness matrices
of all the elements into the local sparse matrix of the tangent, with the
positions of the entries computed once
(`batch::computeStiffnessMatricesPositions`) and without them (search of
each entry in its row, by bisection), in the same program, alternating the
two, 10 times each. Arbitrary blocks of the tangent, one process.
mfem-mgis branch `positions` (d04acab66c), separate build installed in
`install/mfem-mgis-positions`. Quiet machine, one core (`taskset -c 4`) or
the RTX 5080.

Outputs: `host-2.txt` (one core) and `gpu-2.txt` (`--device cuda`). The
files `host.txt` and `gpu.txt` are the first series (b64d5dda96): on the
GPU, each timed assembly counted there the copy of the whole matrix to the
host (times inflated by 10 to 40 %), and the time of the positions the
copy of the positions; the ratios are the same.

## One core

| Case | Elements | Entries of the matrix | Positions | Without positions | With positions | Gain | Paid back after |
|---|---|---|---|---|---|---|---|
| beam, P1 tetrahedra, mechanics (r4) | 196,608 | 4.7 M | 57 MB, 419 ms | 235 ms | 117 ms | ×2.0 | 4 assemblies |
| thermal, P1 tetrahedra (r5) | 1,572,864 | 4.1 M | 50 MB, 272 ms | 382 ms | 133 ms | ×2.9 | 2 |
| Q1 hexahedra, 27 points, small strain (r4) | 4,096 | 1.1 M | 4.7 MB, 31 ms | 130 ms | 116 ms | ×1.1 | 3 |
| Q1 hexahedra, 27 points, finite strain (r4) | 4,096 | 1.1 M | 4.7 MB, 31 ms | 178 ms | 165 ms | ×1.1 | 3 |
| P2 tetrahedra, mechanics (r3) | 24,576 | 8.9 M | 44 MB, 357 ms | 337 ms | 216 ms | ×1.6 | 3 |
| Q2 hexahedra, mechanics (r2) | 64 | 0.3 M | 0.8 MB, 5 ms | 44 ms | 41 ms | ×1.1 | 2 |

Median times; the time of the positions is the one of a second
computation. Bitwise identical matrices, no entry searched with the
positions. The computation of the positions costs from 0.1 to 1.8
assemblies. It is paid back between the second and the fourth tangent.
The gain is large when the computation of the blocks is light (linear
elements, few points, thermal). With 27 points per hexahedron, the
computation of the blocks dominates (see
`results/registers-kernels-2026-10-11`, which reduces it). Memory: 2 n²
bytes per element, of the same order as the local matrix for linear
elements.

## GPU

| Case | Elements | Positions | Without positions | With positions | Gain |
|---|---|---|---|---|---|
| beam, P1 tetrahedra (r5) | 1,572,864 | 453 MB | 24.6 ms | 24.7 ms | ×0.99 |
| thermal, P1 tetrahedra (r5) | 1,572,864 | 50 MB | 13.0 ms | 13.0 ms | ×1.00 |
| Q1 hexahedra, 27 points (r5) | 32,768 | 38 MB | 42.9 ms | 42.5 ms | ×1.01 |
| P2 tetrahedra (r4) | 196,608 | 354 MB | 67.5 ms | 67.1 ms | ×1.01 |

No gain on the GPU: the search is hidden there by the computation of the
blocks and the atomic additions, and the positions would take up to
450 MB. Decision: `MultiMaterialNonLinearIntegrator` does not use them
with a GPU (nor with the `debug` device of MFEM). It keeps them for each
layout of the rows (columns sorted by a direct solver, or diagonal first),
recomputes them after an assembly where entries were searched, at most
three times, then gives them up.

The benchmark only measures the standard kernel. The kernels of FBar (rank
one) and of Faltus read the same positions.
