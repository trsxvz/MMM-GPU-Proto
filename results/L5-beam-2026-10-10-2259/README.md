# AmgX against BoomerAMG for NewtonK0 (10/10/2026)

Beam r4 of the prototype (111,843 unknowns, 10 steps), all the kernels,
mfem-mgis de2e76a913 (wrapper of AmgX specific to mfem-mgis, AmgX patched
for CUDA 13). Script: `cases/measure-L5-amgx.sh`. One RTX 5080, 16 cores.
The sweep and the repetitions of 11/10 (`sweep-k0`, `repetitions-k0.txt`)
were made with mfem-mgis 309849dff3, which uploads the matrix of a single
process as a local matrix (and not a distributed one): same numbers of
iterations (110 and 111).

Parameters of AmgX passed by `-amgxo` with their names in AmgX. "Tuned":
`strength_threshold=0.8, interpolator=MULTIPASS` (best setting of a sweep
on K0, see below). "Single": computations of AmgX in single precision
(`-amgxs`, mode `dFFI`).

## Operator K0 alone (`beam-k0`, GPU)

Preconditioned conjugate gradient, relative tolerance 1e-8. Same
discretisation and same elastic behaviour as the beam. Time of one solve,
setup apart.

| Solver | Setup | Solve | Iterations |
|---|---|---|---|
| CG of MFEM + BoomerAMG | in the 1st solve (+25 ms) | 70 ms (69 to 89 ms) | 70 |
| CG of MFEM + aggressive BoomerAMG | in the 1st solve (+30 ms) | 70 ms | 80 |
| PCG of hypre + BoomerAMG | in the 1st solve (+22 ms) | 65 ms | 70 |
| CG of MFEM + default AmgX | 190 ms | 1.95 s | 1078 |
| CG of MFEM + default AmgX, single | 176 ms | 1.10 s | 1110 |
| CG of MFEM + tuned AmgX | 71 ms | 92 ms (75 to 92 ms) | 111 |
| CG of MFEM + tuned AmgX, 2 sweeps | 66 ms | 108 to 123 ms | 96 |
| **CG of MFEM + tuned AmgX, single** | 67 ms | **53 to 61 ms (54 to 72 ms)** | 110 |
| tuned PCG of AmgX | 63 ms | 78 ms | 116 |
| CG of MFEM + Jacobi | 1 ms | 227 ms | 1862 |

In parentheses: three new runs of 11/10 (`repetitions-k0.txt`). From one
solve to the other, the time varies by about 15 % on this GPU, for the
same number of iterations: only larger differences are significant. The
order does not change: tuned AmgX in single precision ahead of BoomerAMG
and of tuned AmgX in double.

The default setting (the one of MFEM: threshold 0.25, interpolation D2) is
a bad preconditioner for this vector elasticity, handled as a scalar
matrix. The threshold 0.8 and the MULTIPASS interpolation make it
comparable to BoomerAMG, whose "System" strategy takes the three components
into account. One iteration of tuned AmgX costs 0.83 ms in double and 0.48
to 0.55 ms in single, against 1.0 ms for BoomerAMG: the FP64 of the RTX
5080 is slow, and single precision speeds up the multigrid by 1.5 to 1.7
times without changing the number of iterations.

Sweep on K0 (`cases/sweep-amgx-k0.sh`, directory `sweep-k0`: JSON
configurations, raw outputs and `summary.txt`, redone on 11/10 with the
current version of `beam-k0`). Iterations and time of one solve:

| Configuration | Iterations | Solve |
|---|---|---|
| classical default (threshold 0.25, D2) | 1078 | 1.8 to 1.9 s |
| classical default, 2 sweeps | 779 | 1.9 s |
| classical default, aggressive (1 level) | 1283 | 0.75 s |
| threshold 0.5 | 1071 | 1.0 to 1.2 s |
| threshold 0.5, 2 sweeps | 769 | 1.15 to 1.18 s |
| threshold 0.6 | 128 | 165 ms |
| threshold 0.6, 2 sweeps | 102 | 197 to 203 ms |
| threshold 0.7 | 116 | 130 ms |
| threshold 0.8 | 118 | 107 ms |
| threshold 0.85 | 124 | 99 ms |
| threshold 0.9 | 133 | 102 ms |
| **threshold 0.8, MULTIPASS** | **111** | **75 to 77 ms** |
| MULTIPASS, threshold 0.6 or 0.7 | 129 to 130 | 90 to 93 ms |
| MULTIPASS, threshold 0.7, 2 sweeps | 107 | 118 to 119 ms |
| MULTIPASS, threshold 0.8, 2 sweeps | 96 | 101 to 104 ms |
| MULTIPASS, threshold 0.9 | 125 | 82 to 85 ms |
| threshold 0.8, 2 sweeps | 98 | 140 ms |
| threshold 0.8, HMIS | 117 | 178 ms |
| threshold 0.8, aggressive (1 level) | 252 | 170 ms |
| threshold 0.5, F cycle, W, dense LU or max_row_sum 1 | 1055 to 1588 | 1.0 to 4.7 s |
| classical, interpolation D1 | 1602 | 4.5 s |
| aggregation, multiple pairs | 925 | 1.15 s |
| aggregation SIZE_2, multicolour Gauss-Seidel | 458 to 483 | 5.7 to 6.0 s |
| aggregation, block Jacobi (SIZE_2, SIZE_8) | not converged | |
| classical, multicolour Gauss-Seidel | 806 to 6920 | 20 to 166 s |

The strength threshold decides everything: up to 0.5 the classical
multigrid is almost useless, from 0.6 on the number of iterations is
divided by eight. Two sweeps of the smoother remove 13 to 28 % of the
iterations but cost more than they bring. MULTIPASS, threshold 0.8, 2
sweeps gives 101 to 104 ms here, against 108 to 123 ms in the previous
table (same configuration passed by `-amgxo`): the gap is within the
variability of the times. With the block Jacobi, the aggregation gives a
preconditioner that the conjugate gradient does not tolerate (stop after 3
or 4 iterations, relative residual 0.58 to 0.69).

## Full computation, NewtonK0 + Eisenstat-Walker 1

| Multigrid of K0 | GPU | Iterations | CPU, 8 processes | Iterations |
|---|---|---|---|---|
| default BoomerAMG | 1.99 s | 33 | 9.23 s | 38 |
| aggressive BoomerAMG | 1.77 to 1.87 s | 34 to 36 | **7.40 s** | 34 |
| default AmgX | 51.5 s | 90 | 849 s | 92 |
| tuned AmgX | 1.95 to 2.04 s | 36 to 38 | 133 s | 40 |
| tuned AmgX, 2 sweeps | 2.39 s | 38 | 167 s | 44 |
| **tuned AmgX, single** | **1.58 to 1.60 s** | 37 | 128 s | 39 |

Time of the loop over the steps. On the GPU, four runs (the one of the
script and three repetitions, `repetitions-gpu.txt`) for the rows given as
a range and for tuned AmgX in single precision. Same solution in all the
cases (final norms within 1e-7 in relative terms).

- On the GPU, tuned AmgX in single precision is the best result on the
  beam: 1.58 to 1.60 s against 1.77 to 1.87 s for the best BoomerAMG (-10
  to -15 %).
- With 8 CPU processes, AmgX is unusable: each process has its AmgX
  resources on the same GPU. Probable cause, not measured: eight CUDA
  contexts share the GPU by time slices (no MPS), and the halo exchanges go
  through the host. AmgX is only of interest here with one GPU per process.

## Possible next steps

- Block matrices 3×3 (ordering by components): multigrid of AmgX adapted
  to systems, to compare with the "System" strategy of BoomerAMG.
- Single precision for BoomerAMG (hypre in mixed precision) and for the
  multigrid of Ginkgo.
- Setup of AmgX without going through the host, reuse
  (`AMGX_solver_resetup`) for an updated K0.
