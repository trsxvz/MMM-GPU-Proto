# Benchmark of the examples (CPU against GPU)

## Goal

For each example of mfem-mgis-examples (ex1 to ex8) and each case of
mm-opera-hpc (polycrystal, bubble, cermet), two wall times on the same
machine, for the same case and the same tolerances:

- T_cpu: the best time with everything on the CPU, in parallel (MPI), over
  the solvers and the numbers of processes tried;
- T_gpu: the best time with at least one phase on the GPU (behaviour
  integration, finite element kernels, linear solver), over the same kind
  of choices;
- the speedup T_cpu / T_gpu, reported as measured, also when it is below
  one.

## Machine and rules

- 8 cores with SMT (cores N and N+8 are twins), 30 GB, one RTX 5080
  (16 GB). CPU runs on 1, 2, 4 and 8 processes, one per physical core
  (cores 0 to 7), `OMP_NUM_THREADS=1`. They need the whole machine: no
  build or other measurement at the same time.
- GPU runs on one process: several processes on one GPU were much slower
  (polycrystal: 91.57 s on 2 processes, 13.18 s on 1,
  `results/L6-polycrystal-2026-10-09`). Sides: `--execution gpu`, and
  `--execution mixed --linear-solvers-execution cpu` when a direct solver
  of the CPU is the best choice.
- One build for both sides (-O3, double): mfem-mgis installed with all the
  solver libraries (hypre, PETSc, MUMPS, cuDSS, AmgX, Ginkgo, Trilinos with
  MueLu, FROSch, Teko and FastILU), and the examples and mm-opera-hpc built
  against it.
- CPU baseline with `--execution cpu`: behaviours integrated by batches and
  finite element kernels on the host. Without `--execution`, the behaviours
  are integrated point by point, a slower baseline which would flatter the
  GPU.
- Wall time of the whole program, and the phases of the profiler of
  mfem-mgis. Best of 3 runs for the retained configurations.
- Each retained run is checked: reference files of the examples when they
  exist (ex2, ex4, ex7, ex8), otherwise comparison with the CPU run
  (`cases/compare-outputs.py`).

## Driver

`cases/benchmark_slides.py` runs one case on one side with one linear
solver and, optionally, one nonlinear solver, and appends one line to
`runs.csv` in the results directory (output in `logs/`):

```
cases/benchmark_slides.py list
cases/benchmark_slides.py run RESULTS CASE SIDE SOLVER [--non-linear NAME]
                          [--np N] [--size SIZE] [--repeat R] [--dry-run]
cases/benchmark_slides.py best RESULTS
```

The solvers are given by the overrides of mfem-mgis, the same for all the
programs, without changing them: `--linear-solver-override` (linear solver
with its preconditioner, table `SOLVERS`: hypre with BoomerAMG, MUMPS,
cuDSS, Belos with MueLu, AmgX, PETSc with GAMG) and
`--non-linear-solver-override` (table `NON_LINEAR_SOLVERS`: Newton,
modified Newton, NewtonK0 whose elastic operator is solved by the linear
solver of the run). The program keeps its convergence criteria. `best`
prints, for each case and size, the best time of each side and the
speedup.

The runs need the libraries of the installation in `LD_LIBRARY_PATH` (TFEL
built with CUDA, MGIS, CUDA). The solvers of Trilinos initialize the GPU
through Kokkos: they only run with the GPU visible.

## Sizes

The default meshes of ex1, ex3, ex6 (about 200 unknowns), ex2 (3,726),
ex5 (7,442) and ex7 (2,634) are too small for a GPU to matter. All the
examples can be refined (`-r`), and ex5 runs in parallel. Proposed sizes
of the measurements (option `--size bench` of the driver; `--size small`
checks the scripts):

| Case | Size | Unknowns | Steps |
|---|---|---|---|
| ex1 | -r 6 | about 824k (estimate) | 10 |
| ex2 | -r 2 | about 60k (estimate) | 50 |
| ex3 | -r 4 | about 786k (estimate) | 1 solve |
| ex4 | -r 1 | 196,047 | 10 |
| ex5 | -r 2, parallel | about 119k (estimate) | 1 solve |
| ex6 | -r 3 | to measure | 1 solve |
| ex7 | -r 1 | to measure | 40 |
| ex8 | -r 1 | 214,812 (mechanics) | 5 |
| polycrystal | -r 1 | 168,885 | 10 |
| bubble | default | 467,109 | 1 solve |
| cermet | default | 33,015 | 10 |

Screening first (one or two steps per pair, the pairs which fail or are
more than twice slower than the best one being dropped), then full runs of
the best two or three pairs of each side.

## Trials at the small size (16/10/2026)

On the CPU (two processes), the driver and the overrides work on all the
cases. On the GPU (one process), seconds, X for a failure with its exit
status:

| Case | gpu, default | gpu, hypre-amg | gpu, cudss | gpu, amgx | gpu, belos-muelu | mixed, mumps |
|---|---|---|---|---|---|---|
| ex1 | X1 | X14 | 0.8 | 1.1 | 0.7 | X14 |
| ex2 | X-6 | 2.2 | 1.2 | 13.5 | 1.5 | 0.7 |
| ex3 | 0.7 | 0.8 | 0.9 | 0.8 | 0.8 | 1.2 |
| ex4 | 2.4 | 2.4 | 1.9 | 8.9 | X1 | 3.7 |
| ex5 | X-6 | 0.6 | 0.6 | 1.0 | 0.6 | 0.6 |
| ex6 | 0.7 | 0.6 | 0.7 | 0.7 | 0.8 | 0.6 |
| ex7 | 1.4 | 1.4 | 1.1 | 4.3 | 1.8 | 1.1 |
| ex8 | X-6 | 47.9 | 3.1 | X1 | X1 | 4.9 |
| polycrystal | 6.5 | 7.8 | 13.2 | 44.0 | 7.4 | 50.3 |
| bubble | 5.1 | 4.7 | X1 | 9.5 | X-6 | 276.8 |
| cermet | 2.7 | 2.2 | 5.3 | 3.7 | 11.0 | 32.8 |

These times only check the scripts: at this size, the setup dominates.
The failures:

- ex2, ex5 and ex8 use MUMPS by default, which is refused on the GPU with
  a message giving the options to use;
- ex1 with Newton: the Newton solver of ex1 is fragile at the onset of
  plasticity, on the CPU as well;
- ex4 with GMRES of Belos and MueLu: the linear solver does not converge;
- ex8 with AmgX and with MueLu: the linear solver override replaces the
  solvers of both problems, and the thermal operator is indefinite. A
  single global override does not suit the programs whose problems need
  different solvers;
- bubble with cuDSS: the factorization fails, probably for lack of memory
  on the GPU (467k unknowns of quadratic elements);
- bubble with MueLu: the conjugate gradient meets a non positive
  curvature, with the Jacobi smoother of MueLu.

## Outputs

`results/benchmark-slides-<date>/`: `runs.csv` and the logs of the
driver, the table of T_cpu, T_gpu and the speedups (`best`), and a
README with the configurations retained and the checks.

## Open questions

- The sizes above.
- Exclusive use of the machine for the CPU runs, and GPU slots for the GPU
  runs.
- Whole program time or time of the resolution only (both are recorded).
