# Scaling of the beam (11/10/2026)

Measurement campaign, case 5. Beam of the prototype refined three, four
and five times, 10 steps, best method of work package L5: NewtonK0 (conjugate
gradient preconditioned by the multigrid of K0, built once), prediction by
K0, Eisenstat-Walker 1 tolerance, all the kernels. Script:
`cases/measure-beam-scaling.sh`. mfem-mgis e03fd927f7 (option of the pool
of hypre in the beam, since taken over by fd93974cdc). One RTX 5080
(16 GB, 15.9 GiB usable), 8 cores (16 threads). Quiet machine during these
computations.

On the GPU, the memory pool of hypre is reduced to 0.5 GiB (`-hdps 0.5`,
see below). Memory: maximum read by `nvidia-smi` every 200 ms, including
686 MiB used by the desktop; in GiB (1 GiB = 1024 MiB).

| Refinement | Unknowns | GPU, BoomerAMG | GPU, AmgX single | GPU, matrix free | 8 cores | Gain of the GPU |
|---|---|---|---|---|---|---|
| 3 | 15,795 | 0.67 s | **0.48 s** | 0.70 s | 0.85 s | ×1.8 |
| 4 | 111,843 | 1.81 s | **1.54 s** | 2.01 s | 7.02 s | ×4.5 |
| 5 | 839,619 | **9.86 s** | 11.98 s | 12.12 s | 107.9 s | ×10.9 |

Time of the loop over the steps; the gain compares the best time of the
GPU with the computation on 8 cores. Iterations of the fixed point: 33 to
38, almost independent of the size. Memory of the GPU at most: 2.29 GiB
with three and four refinements (2.73 and 2.87 GiB with AmgX), 6.43 GiB
with five (7.76 GiB with AmgX, 6.72 GiB matrix free).

- The gain of the GPU grows with the size: ×1.8 at 16,000 unknowns, ×10.9
  at 840,000. The GPU goes from 4 to 5 refinements in 5.5 times more time
  for 7.5 times more unknowns (better occupied), the CPU in 15 times more
  time (the solve goes from 3.4 s to 58 s: more data than the caches,
  hypothesis not measured).
- Tuned AmgX in single precision is the best up to 110,000 unknowns,
  aggressive BoomerAMG beyond (solve 4.4 s against 6.3 s with five
  refinements).
- The tangent applied without assembly does not pay off: the avoided
  assembly (2.0 s with five refinements) costs more in solve (8.0 s
  against 4.4 s), and it saves no memory (K0 stays assembled).

Breakdown (`cases/summary-profiling.py`, seconds, total with the setup):

| Computation | Total | Integration | Tangent | Solve | Prediction |
|---|---|---|---|---|---|
| 5 refinements, GPU, BoomerAMG | 12.0 | 1.1 | 2.0 | 4.4 | 2.6 |
| 5 refinements, 8 cores | 108.1 | 3.1 | 22.0 | 58.4 | 23.6 |

## Six refinements

6,502,275 unknowns: does not fit in the memory of the GPU. With the
default pool of hypre, the allocation of the matrix of K0 (2.14 GiB) is
refused at 13.46 GiB used (`beam-r6-gpu-amg-aggressive-assembled`, same
refusal matrix free for the tangent). With a pool of 0.5 GiB, an
allocation of 1.07 GiB is refused at 15.34 GiB, whether the tangent is
assembled or not (`beam-r6-gpu-amg-aggressive-reserved-hypre-0.5-*`). The
outputs start with the command line and give the error message, except
`beam-r6-gpu-amg-aggressive-matrix-free` (first trial, message only).
Memory of the host: 18.3 GiB.

## Memory pool of hypre

hypre 3.2 allocates on the GPU, through Umpire, a pool of 4 GiB from its
first allocation, whatever the size of the problem. With a pool of
0.5 GiB (which grows if needed), the memory used with three and four
refinements goes from 5.79 GiB to 2.29 GiB (3.5 GiB less). The time
changes by a few %: with BoomerAMG, 0.671 s against 0.644 s with three
refinements, 1.807 s against 1.776 s with four, 9.860 s against 9.469 s
with five (one more iteration), that is 2 to 4 % slower; with AmgX, less
than 1 %. First pass with the default pool:
`../scaling-beam-2026-10-11-0233`. Its computations on 8 cores, and their
rerun, were 13 to 43 % slower, for a cause not established (compilations
of another session were seen during these computations, but the rerun, at
a load of 1.1, gave the same times): only the series of this directory
serves as the reference for the CPU. The option has become an option of
mfem-mgis (`--hypre-device-pool-size` or variable
`MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`).
