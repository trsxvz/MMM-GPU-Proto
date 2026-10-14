# Post-processings by kernels (11/10/2026)

Bubble of mm-opera-hpc, default mesh (`single_sphere.msh`, 111,580
quadratic tetrahedra, 467,109 unknowns), one step, PCG and BoomerAMG,
everything by kernels for the mechanics. Post-processings of the bubble:
export of the stresses at the nodes
(`ParaviewExportIntegrationPointResultsAtNodes`), principal stress and
hydrostatic pressure at the points, mean of the hydrostatic pressure
(`computeMeasure`, `computeIntegral`), L2 projection of the principal
stress, ParaView export.

Comparison of `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=ElementByElement`
(loops over the elements, the default) and `Kernels` (kernels of
mfem-mgis, branch post-processings, 0434dfd1de, standard host memory).
Executable: `build/mm-opera-hpc-post` (mm-opera-hpc 79c49a8, mean by
`computeMeasure` and `computeIntegral`). Outputs of the bubble identical to
the digits written in all cases.

## GPU (one process, cuda behaviour, two runs)

| Phase (s) | ElementByElement | Kernels |
|---|---|---|
| total | 14.77 and 14.86 | 5.32 and 5.33 |
| solution | 1.48 and 1.48 | 1.47 and 1.48 |
| post-processings | 11.32 and 11.40 | 1.86 and 1.87 |
| export at the nodes | 4.46 and 4.46 | 0.85 and 0.85 |
| hydrostatic mean | 1.90 and 1.91 | 0.008 and 0.008 |
| L2 projection | 4.35 and 4.42 | 0.39 and 0.40 |
| of which operator | 2.42 and 2.46 | 0.037 and 0.037 |
| of which right-hand side | 1.42 and 1.47 | 0.002 and 0.002 |
| of which solution | 0.49 and 0.46 | 0.32 and 0.33 |
| export of the principal stress | 0.49 and 0.49 | 0.49 and 0.49 |

## CPU (one core, one run)

| Phase (s) | ElementByElement | Kernels |
|---|---|---|
| total | 33.25 | 25.19 |
| post-processings | 12.16 | 5.35 |
| export at the nodes | 4.47 | 1.47 |
| hydrostatic mean | 1.89 | 0.33 |
| L2 projection | 5.20 | 2.95 |
| of which operator | 2.51 | 1.67 |
| of which right-hand side | 1.41 | 0.01 |
| of which solution | 1.25 | 1.07 |

- On the GPU, the post-processings run 6 times faster and the whole
  bubble 2.8 times faster. With the fix of the host memory
  (`results/host-memory-2026-10-11`), the bubble on the GPU went from 69 s
  to 5.3 s, against 25 s on one CPU core.
- The mass matrix of the tetrahedra is computed by a kernel and assembled
  by MFEM (`AssemblyLevel::FULL`): MFEM 4.10 does not assemble the mass of
  simplices on the device.
- What remains on the GPU: the solution of the L2 projection (conjugate
  gradient preconditioned by the ILU of hypre), the ParaView exports (on
  the host) and the gathering of the values at the points (on the host,
  where the functions of the user are stored).
- Final version after review (mfem-mgis affeddab21: mass computed by the
  kernel for all the elements, mesh not modified, collective decisions),
  same quiet machine: GPU 14.86 s element by element, 5.32 s by kernels
  (post-processings 11.41 and 1.87 s; files `bubble-gpu-*-final.txt`),
  identical outputs.
- The element by element path of the mean now goes through
  `computeMeasure` and `computeIntegral` (two loops): 1.9 s against 1.0 s
  for the former loop of the bubble.
