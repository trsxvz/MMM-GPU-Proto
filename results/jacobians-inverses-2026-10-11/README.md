# Inverses of the Jacobians computed once (11/10/2026)

The kernels only used the Jacobians through their inverses, recomputed by
each thread at each point (for each pair of nodes in the kernel of the
stiffness matrices): 9 divisions and about thirty multiplications in
double precision, slow on this GPU. mfem-mgis (branch `jacobians`,
bf7673f4f8): the kernel of the geometry stores the inverses in place of
the Jacobians (same memory), and the kernels read them.

Bench `cases/assembly`, separate installations
`install/mfem-mgis-registers` (4a007e85f2, before) and
`install/mfem-mgis-jacobians` (after), same commands, one core
(`taskset -c 4`) or the RTX 5080, quiet machine. Median times of 10
assemblies. Outputs: `gpu.txt`, `host.txt`.

GPU (without positions, the default path with a GPU):

| Case | Before | After | Gain |
|---|---|---|---|
| finite strain, Q1 hexahedra 27 points (r5) | 31.8 ms | 19.0 ms | ×1.67 |
| small strain, Q1 hexahedra 27 points (r5) | 29.2 ms | 16.3 ms | ×1.79 |
| beam, P1 tetrahedra (r5) | 18.1 ms | 8.1 ms | ×2.2 |
| finite strain, P2 tetrahedra (r4) | 52.8 ms | 26.4 ms | ×2.0 |
| thermal, P1 tetrahedra (r5) | 13.8 ms | 2.3 ms | ×5.9 |

Gaps between the two versions of 1e-17 to 6e-17 for values of 4e-3 to
0.27: order of the atomic additions.

One core (with the positions): 52.3 against 51.5 ms (finite strain), 75.6
against 72.3 ms (beam), 116.7 against 109.7 ms (thermal); bitwise
identical values. The host path only inverted the Jacobians once per
point.

Registers (cuobjdump): stiffness matrices of the 3D thermal 56 instead of
84 (64 before the work package of the matrices without forming B);
gradients 40 to 46; the rest unchanged.
