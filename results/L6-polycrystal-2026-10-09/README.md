# Polycrystal on GPU, minimal slice (09/10/2026)

Polycrystal of mm-opera-hpc: 5 grains, single crystal behaviour of UO2
(orthotropic, finite strain, implicit), periodic problem, uniaxial
tension, 10 steps. PCG + BoomerAMG of hypre. mfem-mgis decb6fd2c1, MGIS
a0b5140. Script: `cases/measure-L6.sh polycrystal`. One RTX 5080, 16
cores. Quiet machine (load decreasing at the start, after a review).

Computations on GPU: `cuda` device of MFEM, hypre on the GPU, all the
kernels (the polycrystal uses no prediction), `generic` behaviour
(integrated on the host) or `generic-parallel` compiled for `cuda`
(integrated on the GPU, rotations of the grains included).

| Case | Total | Integration | Residual | Tangent | Solve |
|---|---|---|---|---|---|
| CPU, 8 processes, point by point | 31.65 | 14.86 (158) | 0.65 (158) | 1.78 (98) | 14.11 |
| CPU, 8 processes, everything by kernels | 30.55 | 14.03 (158) | 0.10 (158) | 2.09 (98) | 14.08 |
| CPU, 1 process, point by point | 151.34 | 105.62 (158) | 4.20 (158) | 8.70 (98) | 32.49 |
| GPU, 1 process, `generic` behaviour | 106.70 | 99.16 (158) | 0.15 (158) | 0.35 (98) | 5.56 |
| GPU, 1 process, `cuda` behaviour | **13.18** | 5.71 (158) | 0.07 (158) | 0.22 (98) | 5.60 |
| GPU, 2 processes, `cuda` behaviour | 91.57 | 12.18 (158) | 0.17 (158) | 0.95 (98) | 76.99 |

Times in seconds (maximum over the processes), number of calls in
parentheses. Macroscopic stresses of each computation equal to those of
the point by point computation on the same number of processes, up to
2.4e-13 to 4.4e-13 of the stress scale (`compare-outputs.py
--groups 0 1-3 4-6`).

## Findings

- **One GPU against 8 cores**: 13.18 s against 31.65 s (×2.4); against one
  core, ×11.5.
- **The behaviour dominates on CPU** (105 s out of 151 s with one
  process): on the GPU, the integration drops to 5.71 s (36 ms per call,
  against 14.9 s on 8 cores), on par with the linear solve (5.60 s).
- **Behaviour on the host with the GPU**: 106.7 s. The single crystal
  behaviour on a single core cancels the whole gain of the kernels.
- **Two processes on the same GPU**: solve 14 times slower (77 s), as for
  the beam; see `results/L6-beam-2026-10-09/`.
