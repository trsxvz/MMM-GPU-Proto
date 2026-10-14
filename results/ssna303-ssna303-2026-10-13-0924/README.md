# Notched bar ssna303 (13/10/2026)

Measurement campaign, case 2, resumption of the series of 11/10
(`results/ssna303-ssna303-2026-10-11-0625`, stopped for lack of memory).
Example ex4 of mfem-mgis-examples, variant solved by hypre (Newton, FGMRES
of hypre preconditioned by BoomerAMG), finite strain plasticity (Hencky
logarithmic strain), complete loading in 10 steps, mesh refined zero, one
and two times. Script: `cases/measure-ssna303.sh`, with
`examples_build=build/examples-f57072220f` (ex4 built against a separate
installation of mfem-mgis f57072220f, `install/mfem-mgis-f57072220f`). One
RTX 5080 (16 GB), 8 cores (16 threads), 30 GB of memory; another work
session on the machine, with no computation during the series.

| Refinement | Unknowns | 8 cores, point by point | 8 cores, kernels | GPU, kernels | Gain of the GPU |
|---|---|---|---|---|---|
| 0 | 27,945 | 7.13 s | 5.18 s | 6.89 s | ×0.75 |
| 1 | 196,047 | 71.0 s | 59.3 s | 26.0 s | ×2.3 |
| 2 | 1,461,915 | 693 s | 601 s | does not fit | |

FBar (`--use-fbar`):

| Refinement | 8 cores, point by point | 8 cores, kernels | GPU, kernels | Gain of the GPU |
|---|---|---|---|---|
| 0 | 7.64 s | 5.86 s | 7.28 s | ×0.80 |
| 1 | 90.9 s | 74.0 s | 27.9 s | ×2.7 |
| 2 | not launched | not launched | does not fit | |

Time of the loop over the steps; the gain compares the GPU to the kernels
on 8 cores.

- Compared to 11/10 (mfem-mgis 8a6325335e), the kernels on 8 cores gain
  10% at refinement 1 (65.5 → 59.3 s; FBar 78.7 → 74.0 s) and at
  refinement 2 (670 → 601 s); the GPU gains 4 to 5% (27.1 → 26.0 s; FBar
  29.5 → 27.9 s). These gains come from the work packages of 11/10
  (positions of the entries, matrices without forming B, inverses of the
  Jacobians).
- Refinement 2 on 8 processes: the memory of the machine is almost full
  during the computation by kernels (27.5 GB out of 30, 4.5 GB of swap);
  FBar is not launched there (it had exhausted the memory on 11/10).
- Refinement 2 on the GPU: with the 27-point rule of the hexahedra, the
  memory of the GPU is not enough (11/10). The order 3 rule (8 points,
  option `-iro 3`) was tried: it does not suit this case. On 8 cores,
  Newton fails 8 times and the steps concerned are split (810 s against
  601 s); on the GPU, Newton stagnates at a relative residual of 1.1e-4
  at the first step, then the splitting fails for lack of memory
  (allocation of 451 MB by hypre in the triple product). This refinement
  does not fit on this GPU with the complete rule; the memory at the
  points would have to be reduced (blocks of the tangent in single
  precision or computed on the fly), or several GPUs used.
- At refinement 0 (28,000 unknowns), the GPU remains slower than 8 cores:
  the system is too small for it.
