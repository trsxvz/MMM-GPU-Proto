# Measurements of work package L2 on the beam (09/10/2026, provisional)

**Superseded** by `results/L2-beam-2026-10-09/` (quiet machine, threads
placed correctly, compilation options of the behaviour fixed).

Batch integration of the behaviours (option `BehaviourIntegration`,
mfem-mgis ef52adb445), beam r4 (111,843 unknowns), PCG + BoomerAMG,
tolerance 1e-4, 10 steps. Script: `cases/measure-L2.sh`; summary by phase:
`cases/summary-profiling.py *.txt`.

**Provisional measurements**: another compilation was running on the
machine (load of 5 to 9). The reference time of work package L0 for the same case,
on a quiet machine, is 15.45 s with 8 processes. Here we measure 20.02 s
for the same path (point by point). Differences below 10 to 20 % are not
significant.

## Results

All variants give the same solution (norm 48.83704586) and the same Newton
iterations (31 with 8 processes, 28 with 1 and 2 processes).

Time per phase (seconds, maximum over the processes; in parentheses, the
number of calls). Integration: `NS::processNewUnknownsEstimate`, which
includes the computation of the gradients element by element and the call
of the behaviour.

| Case | Total | Integration | Residual | Tangent | Solve | Prediction |
|---|---|---|---|---|---|---|
| 8 processes, point by point | 20.02 | 0.86 | 0.74 | 2.89 | 11.84 | 3.61 |
| 8 processes, batches, `generic` | 19.59 | 0.73 | 0.76 | 2.97 | 11.55 | 3.49 |
| 8 processes, batches, `openmp` 1 thread | 17.19 | 0.59 | 0.53 | 2.43 | 10.40 | 3.16 |
| 1 process, point by point | 67.06 | 3.36 | 2.97 | 11.72 | 34.92 | 13.87 |
| 1 process, batches, `generic` | 67.62 | 3.25 | 2.96 | 11.73 | 35.25 | 14.20 |
| 1 process, batches, `openmp` 1 thread | 67.55 | 3.24 | 3.01 | 11.76 | 35.16 | 14.16 |
| 1 process, batches, `openmp` 8 threads | 70.35 | 2.93 | 3.20 | 12.19 | 36.95 | 14.85 |
| 1 process, batches, `stlpar` 8 cores | 69.02 | 3.37 | 3.03 | 11.96 | 36.07 | 14.38 |
| 2 processes, batches, `openmp` 4 threads | 37.48 | 1.40 | 1.63 | 6.71 | 19.96 | 7.66 |

## Findings

- **Correctness**: the batch mode changes no result, whatever the library
  (`generic`, `generic-parallel` openmp or stlpar) and the number of
  threads.
- **The beam is not the right case for L2**: the integration of the
  behaviour (von Mises plasticity with linear hardening) weighs 5 % of the
  time. The linear solve (52 %), the prediction (21 %) and the assembly of
  the tangent (17 %) dominate. Several threads only reduce the call of the
  behaviour, not the computation of the gradients, which stays sequential
  element by element in the same phase.
- **Several threads per process**: no visible effect here (1 process with
  8 threads: 2.93 s of integration against 3.24 s with 1 thread). With 8
  processes of one thread, the `generic-parallel` library saves 0.27 s out
  of 0.86 s of integration, within the noise of the loaded machine.
- **Placement**: with `--bind-to none`, OpenMPI 5 binds all the processes
  to the first core (first trial abandoned: 8 processes on one core). The
  script now reserves cores per process
  (`--map-by slot:PE=<threads> --bind-to core`).
- **Case to measure next**: the polycrystal (55 % of the time in the
  integration in work package L0, orthotropic crystal behaviour). It needs the
  support of orthotropic `generic-parallel` behaviours (rotations given by
  fields), in progress.
