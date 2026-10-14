# Registers of the GPU kernels, stiffness matrices without forming B (11/10/2026)

## Finding

`cuobjdump -res-usage` on the objects of the kernels of mfem-mgis (CUDA
13.0, sm_120). Columns: file, registers per thread, stack (bytes), local
memory (bytes), kernel. `kind=N`: `GradientsKind` in the order of its
declaration (0 3D small strain, 1 2D small strain, 2 3D deformation
gradient, 3 2D deformation gradient, 4 and 5 temperature gradient, 6 and 7
temperature, 8 and 9 micromorphic).

- `registers.txt`: branch `positions` (b64d5dda96). The stiffness matrices
  in 3D finite strain (kind 2, ssna303 and polycrystal cases) took 254
  registers per thread, the maximum, and 576 bytes of stack: at most 8
  warps per multiprocessor. In 3D small strain: 122 to 126 registers, 576
  bytes of stack. The stack came from the arrays `Ba[9][3]`, `Bb[9][3]`,
  while `B` is sparse (in finite strain, `F_ij` only depends on `u_i`,
  through `dN/dX_j`).
- `registers-after.txt`: branch `registers` (4a007e85f2), products `K B_b`
  and `B_aᵀ (K B_b)` of the stiffness matrices, and products of the action
  of the tangent without assembly, written with a constexpr table of the
  nonzero entries of `B` (`getBOperatorEntries`, checked at compile time
  against `computeB` for each derivative taken alone), in the same order
  of operations; the kernels of the diagonals now only compute the
  diagonal of the blocks. 3D finite strain: 64 registers and 72 bytes of
  stack (the remaining stack is the block `m[3][3]`); 3D small strain: 80
  and 72; diagonals in 3D mechanics: 58 to 60; the 2D cases no longer have
  any stack. The 3D temperature gradient and the 3D micromorphic go from
  64 and 72 to 84 registers (94 and 92 for their diagonals): two blocks of
  256 threads per multiprocessor instead of four.

## Measurements

Bench `cases/assembly` (assembly by kernels of all the stiffness matrices
into the local sparse matrix, arbitrary blocks of the tangent, one
process), built against the separate installations of the branches
`positions` (d04acab66c, before) and `registers` (28842e74dd, after; the
fixes of review N, 4a007e85f2, came after these measurements: they reduce
the work memory of the host path and the computation of the diagonals,
without changing the values), same commands. Quiet machine, one core
(`taskset -c 4`) or the RTX 5080. Median times of 10 assemblies. Outputs:
`host-2.txt`, `gpu-2.txt` (the files `host.txt` and `gpu.txt` are a first
series where the time measured on the GPU included the copy of the whole
matrix to the host).

One core (with the positions of the entries, the default path without a
device):

| Case | Before | After | Gain | Total gain |
|---|---|---|---|---|
| finite strain, Q1 hexahedra 27 points (r4) | 164 ms | 52 ms | ×3.1 | ×3.4 |
| small strain, Q1 hexahedra 27 points (r4) | 116 ms | 48 ms | ×2.4 | ×2.7 |
| beam, P1 tetrahedra (r4) | 116 ms | 75 ms | ×1.5 | ×3.1 |
| finite strain, P2 tetrahedra (r3) | 277 ms | 133 ms | ×2.1 | ×3.0 |
| thermal, P1 tetrahedra (r5) | 137 ms | 121 ms | ×1.1 | ×3.2 |

"Total gain": with respect to the assembly without positions of the
branch `positions`, that is to the state of gpu-proto before these two
work packages. Bitwise identical values of the matrices between the two
branches (files written by `--output-values`).

GPU (without positions, the default path with a GPU):

| Case | Before | After | Gain |
|---|---|---|---|
| finite strain, Q1 hexahedra 27 points (r5) | 59.8 ms | 33.0 ms | ×1.81 |
| small strain, Q1 hexahedra 27 points (r5) | 42.5 ms | 30.5 ms | ×1.40 |
| beam, P1 tetrahedra (r5) | 25.0 ms | 18.1 ms | ×1.38 |
| finite strain, P2 tetrahedra (r4) | 102.2 ms | 52.9 ms | ×1.93 |
| thermal, P1 tetrahedra (r5) | 13.0 ms | 13.8 ms | ×0.94 |

On the GPU, the values differ from one computation to the next, with the
same program, because of the order of the atomic additions: maximum gap
1.4e-17 for values of 5.5e-2, the same between the two branches. The
thermal loses 6% on the GPU (84 registers instead of 64 for the 3D
temperature gradient); the kernels of the diagonals and of the action were
not timed.

Identical results: for finite values; with the multiplications and
additions contracted (default option of nvcc), only the sign of the zeros
coming from subnormal products may change. The components of the tangent
multiplied by null entries of `B` are no longer read: a non-finite value at
these places is no longer propagated.

The kernel of the action of the tangent without assembly (option
`MatrixFree`) uses the same table (second commit): 80 then 72 registers in
3D finite strain, identical results; not measured by this bench.
