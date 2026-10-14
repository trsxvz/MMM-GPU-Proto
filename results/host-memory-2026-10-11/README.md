# Host memory of MFEM on the GPU (11/10/2026)

Script: `cases/measure-host-memory.sh <case>`, one process and the GPU,
quiet machine, one run.

## The defect

MFEM compiled with Umpire (this is the case here, as for hypre) takes
Umpire for the host memory as soon as a device is configured
(`HOST_UMPIRE`). Umpire records each allocation in its table. The loops
over the elements done on the host allocate and free temporaries at each
point (shape functions, element matrices): each of them goes through this
table, under a lock. Diagnosis by call stacks taken with gdb during the
post-processings of the bubble: all in `umpire::...::AllocationMap`, called
from `CalcShape` and `CalcDShape`. The CUDA calls only accounted for 1 s
out of 61.

Fix (mfem-mgis 0c61fbfcb4, branch host-memory): before creating the
device, `mfem_mgis::initialize` selects the standard memory for the host
and keeps Umpire for the GPU, unless the variable `MFEM_MEMORY` of MFEM is
given.

## Variants

- `umpire`: `MFEM_MEMORY=umpire`, Umpire for the host and the GPU (choice
  of MFEM, the one of all the previous GPU measurements of the prototype);
- `host`: `MFEM_MEMORY=host`, standard memory for the host and the GPU;
- `default`: choice of mfem-mgis after the fix, standard memory for the
  host and Umpire for the GPU.

Executables: mm-opera-hpc (`build/mm-opera-hpc-memory`) and beam
(`build/cases-memory`) linked to `install/mfem-mgis-memory` (fix); ex8 of
`build/examples-pointwise-ex8`, linked to an installation without the fix
(its `default` variant would be equal to `umpire`, it is not measured).
Cuda behaviours of `build/cases-memory`, unless point by point is mentioned
(`generic` behaviours, integrated on the host).

Identical results between variants: bubble to the digits written, cermet
within 5.4e-15, polycrystal within 2.5e-13 (70 values).

## Times (seconds, total of the profile)

| Case | umpire | host | default |
|---|---|---|---|
| bubble (467,109 unknowns, everything by kernels, post-processings included) | 69.17 | 13.83 | 13.97 |
| cermet, 1 step, everything by kernels | 3.52 | 1.53 | 1.53 |
| cermet, 1 step, point by point | 28.10 | 14.81 | 14.80 |
| ex8, 5 steps, everything by kernels, cuda behaviours | 18.64 | 16.89 | |
| ex8, 5 steps, point by point | 164.57 | 57.30 | |
| polycrystal, 10 steps, everything by kernels | 13.26 | 11.76 | 12.32 |
| beam r4, 10 steps, everything by kernels | 3.55 | 3.22 | 3.31 |

- The gain follows the share of the work done on the host by loops over
  the elements of MFEM. Bubble: post-processings (export of the stresses at
  the nodes, L2 projection) from 60.8 s to 10.4 s (`host`) and 10.5 s
  (`default`), against 11.3 s on one CPU core. ex8 point by point:
  integration, residual and tangent (assembled on the host) from 149 s to
  41.7 s. ex8 everything by kernels: the Robin condition, evaluated on the
  host, brings the residual from 1.47 to 0.31 s and the tangent from 1.10
  to 0.42 s.
- When everything is computed by kernels on the GPU (polycrystal, beam),
  the gap stays of the order of the noise of one run (from 3 to 11%).
- `host` and `default` cannot be told apart here: the GPU memory of Umpire
  brings and costs nothing visible. The fix keeps it, as MFEM does.
- Bubble on the GPU after the fix: 14 s, against 31 s on one CPU core,
  everything by kernels. Before, the GPU was twice as slow as the CPU.
- The earlier GPU measurements of the prototype all ran with `umpire`.
  Those where the GPU computes everything by kernels change little. Those
  that keep loops on the host were penalised: point by point integration
  on the GPU, post-processings, other integrators evaluated on the host.
