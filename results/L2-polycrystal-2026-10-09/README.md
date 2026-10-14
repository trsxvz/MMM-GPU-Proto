# Measurements of work package L2 on the polycrystal (09/10/2026)

Polycrystal of mm-opera-hpc: 5 grains, single crystal behaviour of UO2
(finite strain, orthotropic, implicit scheme), 20,553 unknowns, 10 steps,
PCG + BoomerAMG. mfem-mgis 4a68166ca9 (work package L2 reviewed),
`generic-parallel` behaviour compiled with the fixed options (patch 0093 of
TFEL). Quiet machine (load 0.9 at the start). Script:
`cases/measure-L2.sh polycrystal`; summary by phase:
`cases/summary-profiling.py *.txt`.

Integration: `NS::processNewUnknownsEstimate` (gradients element by
element, then call of the behaviour). Times in seconds, maximum over the
processes.

| Case | Total | Integration | Residual | Tangent | Solve |
|---|---|---|---|---|---|
| 1 process, point by point | 150.35 | 104.90 | 4.21 | 8.55 | 32.37 |
| 1 process, batches, `generic` | 148.31 | 103.07 | 4.21 | 8.48 | 32.23 |
| 1 process, batches, `openmp` 1 thread | 150.72 | 105.56 | 3.97 | 8.47 | 32.39 |
| 1 process, batches, `openmp` 8 threads | 63.03 | 17.84 | 4.03 | 8.39 | 32.45 |
| 1 process, batches, `stlpar` 8 cores | 59.13 | 13.94 | 4.08 | 8.33 | 32.47 |
| 2 processes, point by point | 85.64 | 53.35 | 2.27 | 5.65 | 24.09 |
| 2 processes, batches, `openmp` 4 threads | 49.44 | 17.23 | 2.19 | 5.53 | 24.22 |
| 8 processes, point by point | 30.52 | 14.51 | 0.60 | 1.72 | 13.45 |
| 8 processes, batches, `generic` | 30.38 | 14.33 | 0.62 | 1.69 | 13.50 |
| 8 processes, batches, `openmp` 1 thread | 30.73 | 14.72 | 0.58 | 1.70 | 13.49 |

Macroscopic stresses compared with the point by point computation on the
same number of processes: identical with the `generic` behaviour in
batches. Gap of at most 4e-9 (relative to the maximum of each column) with
the `generic-parallel` behaviour, compiled differently, at the level of the
tolerance of the iterative solver.

## Findings

- The behaviour weighs 70 % of the computation on one process. In batches,
  several threads divide it by 5.9 (`openmp`, 8 threads) and by 7.5
  (`stlpar`, TBB on the 16 hardware threads of the 8 cores). The full
  computation goes from 150 s to 59 s.
- One thread per process: the `generic-parallel` behaviour is as fast as
  the `generic` behaviour, once the compilation options are fixed (it was
  17 % slower before patch 0093 of TFEL).
- With 8 processes of one thread, the behaviour takes 14.5 s, like 8
  threads in one process. On CPU, the gain of the threads comes when the
  number of processes is limited (memory, GPU per process). The rest of the
  computation (mostly the solve) then comes first.
- The integration phase also includes the computation of the gradients
  element by element, which is sequential: this is the part that work package L3
  replaces with kernels.

## First trial (abandoned, kept for the record)

The first measurements (directory `scratch/`) gave a sequential `stlpar`
and an `openmp` with 8 threads only 4.3 times faster: `scripts/env.sh`
exports `OMP_PROC_BIND=spread OMP_PLACES=threads`, and libgomp bound the
main thread to the first processor (TBB saw only one processor, and
`spread` only used 4 cores). The measurement script now removes these
variables.
