# Placements of the computations by the --execution option, beam (13/10/2026)

Beam r4 of the prototype (111,843 unknowns), 10 steps, Newton, PCG of hypre
preconditioned by BoomerAMG (tolerance 1e-4), von Mises behaviour with
linear hardening, built for the CPU (`generic` interface) and for the GPU
(`generic-parallel` interface, `cuda` backend). mfem-mgis 8728b7979f
(separate installation `install/mfem-mgis-9a`, build `build/cases-9a`).
Script: `cases/measure-execution-beam.sh`. Quiet machine (the neighbouring
session paused, GPU free).

Everything on the CPU: `--execution cpu`, 8 processes. The other
computations: one process and the GPU, `--execution gpu`, or
`--execution mixed` with the given placement of the finite elements (FE),
of the behaviours and of the linear solver (the beam has no
post-processings).

| Placement | FE | Behaviours | Solver | Time (s) | Copies host to GPU | Copies GPU to host |
|---|---|---|---|---|---|---|
| `cpu`, 8 processes | CPU | CPU | CPU | 12.17 | | |
| `gpu` | GPU | GPU | GPU | 2.96 | 120, 598 MB, 0.035 s | 120, 598 MB, 0.036 s |
| `mixed` | GPU | GPU | GPU | 2.87 | 120, 598 MB, 0.034 s | 120, 598 MB, 0.036 s |
| `mixed` | GPU | CPU | GPU | 3.73 | 184, 2,803 MB, 0.139 s | 82, 387 MB, 0.027 s |
| `mixed` | GPU | GPU | CPU | 42.74 | 120, 598 MB, 0.033 s | 120, 598 MB, 0.036 s |
| `mixed` | GPU | CPU | CPU | 44.66 | 176, 2,671 MB, 0.146 s | 78, 368 MB, 0.024 s |
| `mixed` | CPU | GPU | GPU | 30.62 | 624, 3,108 MB, 0.173 s | 372, 4,307 MB, 0.242 s |
| `mixed` | CPU | CPU | GPU | 31.17 | | |
| `mixed` | CPU | GPU | CPU | 66.00 | 576, 2,869 MB, 0.156 s | 348, 3,998 MB, 0.226 s |
| `mixed` | CPU | CPU | CPU | 66.90 | | |

- The ten placements give the same solution (norm of the displacements
  48.83704586). The number of Newton iterations varies from 28 to 33: the
  BoomerAMG of hypre is not the same on the CPU and on the GPU, and the
  tolerance of 1e-4 of PCG makes Newton inexact.
- The copies of the material arrays cost at most 0.42 s (finite elements
  on the CPU and behaviours on the GPU): the extra cost of the mixed
  placements comes from the processor that computes, not from the
  transfers. With a single process, a phase placed on the CPU has only one
  core: the linear solver on the CPU adds about 40 s, the element by
  element finite elements about 28 s, the integration of the behaviours on
  the CPU 0.86 s (including 0.17 s of copies).
- With everything on the GPU, 120 copies in each direction (0.07 s, 2.5%
  of the computation) remain, 12 per step. Their origin is the one found by
  the analysis of the transfers of 11/10 (`results/L7-transfers-2026-10-11`,
  same 80 and 40 copies per direction on the beam): the end of step update
  is done on the GPU, then the state arrays are brought back to the host,
  as required by the contract of `MaterialArrays` (the state may be
  modified on the host after an update), and the next integration sends
  them back to the GPU. Lead: an option that keeps the states on the GPU,
  the writes on the host being then synchronised explicitly.
- The counters only see the material arrays: the material properties and
  the external variables copied to the GPU at each integration, and the
  copies made by MFEM (vectors, matrices given to hypre) are not counted.
  The total extra cost of a placement is therefore read in the times, the
  counters giving the share of the arrays.
- With everything on the CPU, the finite elements are now computed by
  kernels on the host by default: 12.17 s on 8 processes, against 16.47 s
  element by element (first series, `results/execution-beam-2026-10-13-1418`,
  done before this change; 12.43 s by kernels with the environment
  variables). On 2 processes, element by element: 37.38 s.
