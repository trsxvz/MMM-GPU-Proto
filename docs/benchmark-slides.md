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
  (16 GB). CPU runs on 8 processes (one per physical core) and on 16
  (all the logical cores, `--use-hwthread-cpus`), `OMP_NUM_THREADS=1`.
  They need the whole machine: they are made in exclusive windows, the
  builds of the other users being suspended.
- GPU runs on one process: several processes on one GPU were much slower
  (polycrystal: 91.57 s on 2 processes, 13.18 s on 1,
  `results/L6-polycrystal-2026-10-09`). Sides: `--execution gpu`, and
  `--execution mixed --linear-solvers-execution cpu` (behaviours and finite
  elements on the GPU, linear solver on the CPU). The screening of the GPU
  configurations may be made outside the exclusive windows, the GPU being
  shared by a lock taken around each run (option `--gpu-lock`); the final
  runs are made in exclusive windows.
- Two builds of the same sources: double and float (behaviours generated
  in the precision of MFEM), each with all the solver libraries which
  MFEM accepts in this precision (spack environments `mmm-trilinos-teko`
  and `mmm-trilinos-float`), against which the examples and mm-opera-hpc
  are built (`build/bench` and `build/bench-float`).
- CPU baseline with `--execution cpu`: behaviours integrated by batches and
  finite element kernels on the host. Without `--execution`, the behaviours
  are integrated point by point, a slower baseline which would flatter the
  GPU.
- Wall time of the whole program. The reference is the CPU with the
  solvers of each example as shipped; the best CPU configuration is always
  shown next to it. The CPU and the GPU runs of a case have the same mesh,
  order, time steps and tolerances.

## Driver

`cases/benchmark_slides.py` runs the cases and appends one line per run
to `runs.csv` in the results directory (output in `logs/`); its docstring
gives the usage of its subcommands:

- `list [CASE]`: the sizes and the configurations of the cases;
- `run RESULTS CASE CONFIGURATION`: one configuration;
- `calibrate RESULTS`: each case at the measurement size with the solvers
  of the program, on 8 processes and with the GPU;
- `screen RESULTS CASE`: the configurations of a case (sides, tiers) at the
  screening size (the mesh of the measurements, one or two steps of the
  same time step, or all the steps when the program has no option of the
  end time); the runs already made are skipped, so that an interrupted
  screening resumes;
- `finals RESULTS CASE`: at the measurement size, the configurations
  selected by the screening (the solvers of the program on the CPU with
  their best number of processes, the best other CPU configuration, the
  two best configurations with the GPU), the first ones being repeated
  (three runs) unless their first run took more than five minutes;
- `best RESULTS`: the best times and the speedups.

Common options: `--precision double|float`, `--size small|screen|bench`,
`--timeout S` (a run is stopped with all its processes, `mpirun` and its
ranks, after S seconds; its status is then `timeout`, and its time a lower
bound), `--gpu-lock FILE`, `--dry-run`.

`cases/benchmark_plots.py OUTPUT RESULTS...` draws, for each case, size and
precision, the best time with the solvers of the program on the CPU, the
best CPU configuration and the best configurations with the GPU, with
their speedups, and two overviews of the speedups (against the solvers of
the examples and against the best CPU configuration), with an index
`README.md`. The graphs of the measurements are in `slides/`.

The runs need the libraries of the installation in `LD_LIBRARY_PATH` (TFEL
built with CUDA, MGIS, CUDA). The solvers of Trilinos initialize the GPU
through Kokkos: they only run with the GPU visible.

## Configurations

The solvers are given by the overrides of mfem-mgis, the same for all the
programs, without changing them (`--linear-solver-override`,
`--non-linear-solver-override`), by the tangent operator applied without
assembly (`MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree`, kept with
`--execution`), and for ex8 by its own options (one linear solver per
problem: a global override would give the same one to the thermal and to
the mechanical problems). The programs keep their convergence criteria;
the relative tolerance of the linear solvers is the one of the program,
at least 1e-6 in float. The inexact Newton method and NewtonK0 may make up
to 50 iterations: they need more iterations than the exact Newton method
for which the programs set their maximum (10 for ex2, reached by the
inexact Newton method while it still converges).

- CPU: the solvers of the program on 8 and 16 processes; on 8 processes,
  the Krylov solver of hypre with BoomerAMG, Belos with MueLu (Gauss-Seidel
  smoother, hierarchy reused), PETSc with GAMG, MUMPS, and for the non
  linear computations the inexact Newton method and NewtonK0 (elastic
  operator by one cycle of BoomerAMG or of MueLu).
- GPU, tier 1 (screened first): the solvers of the program when they
  compute on the GPU, hypre with BoomerAMG, Belos with MueLu (Jacobi or
  Chebyshev smoother), AmgX (classical or aggregation multigrid), cuDSS, MUMPS on the CPU with the elements and the
  behaviours on the GPU, the solvers of the program in this mixed
  placement, the inexact Newton method, NewtonK0 with its elastic operator
  by BoomerAMG or MueLu, with or without the tangent operator applied
  without assembly, NewtonK0 with the elastic operator factorized by cuDSS,
  and the conjugate gradient with a Jacobi preconditioner on the tangent
  operator applied without assembly.
- GPU, tier 2 (variants dominated in the measurements of the solvers,
  `docs/solver-recommendations.md`): aggressive coarsening of BoomerAMG,
  Krylov solvers of MFEM with BoomerAMG, AmgX in single precision within the double build, the Jacobi
  preconditioner on the assembled operator, Ginkgo, Tacho, NewtonK0 with
  AmgX, the modified Newton method with cuDSS.

Not used: the sparse solvers of MAGMA (sequential only, the programs are
parallel), AmgX in a float build (refused by mfem-mgis), SuperLU_DIST,
STRUMPACK and SuiteSparse in float (refused by MFEM).

## Sizes

The measurement sizes fill the 8 cores and the GPU. The numbers of
unknowns are printed by the programs, or given by
`cases/tools/count_unknowns` (vector H1 space on the refined mesh).

| Case | Measurement size | Order | Unknowns | Steps |
|---|---|---|---|---|
| ex1 | -r 6, --prediction | 1 | 823,875 | 25 |
| ex2 | -r 4 | 1 | 908,226 | 50 |
| ex3 | -r 4 | 1 | 786,432 | 1 solve |
| ex4 | -r 1 | 1 | 196,047 | 10 |
| ex5 | -r 4, parallel | 2 | 1,847,042 | 1 solve |
| ex6 | -r 4 | 1 | 786,432 | 1 solve |
| ex7 | -r 1 | 2 | printed by the program | 40 |
| ex8 | -r 1 | 1 | 214,812 (mechanics) | 5 |
| polycrystal | -r 1 | 1 | 168,885 | 4 steps of 0.2 s |
| bubble | default mesh | 2 | 467,109 | 1 solve |
| cermet | -r 1 | 1 | 269,253 | 2 steps of 0.2 s |

- ex1: by default, the increment of the imposed displacement is
  concentrated in the elements next to the loaded face. With 4 refinements
  or more, their strain increment is far beyond the elastic limit and the
  Newton solver does not converge when the loading is reversed (with 6
  refinements, not even at the first step). The option `--prediction` of
  the example (elastic prediction, as in ex2 and ex4) removes the failure:
  41 Newton iterations for the 25 steps. ex1 has no option of the end
  time: it is screened with 5 refinements (107,811 unknowns).
- The polycrystal and the cermet run for short durations (the same in
  double and in float): their grain law of UO2 fails in float beyond about
  0.8 s (rounding of its elastic strain when the deformation gradient is
  close to the identity).

## Plan

| Stage | Where | Content |
|---|---|---|
| Calibration | GPU shared, CPU exclusive | each case at the measurement size, solvers of the program |
| Screening, GPU | GPU shared | tier 1 of each case at the screening size |
| Screening, CPU | exclusive | the CPU configurations at the screening size |
| Final runs | exclusive | `finals` of each case, graphs after each case |
| Float | same stages | once the float build is installed |
| Refinement | as available | tier 2, variants of the winners |

The full grid (all the solvers, processes, orders and precisions, three
repetitions) would take about 9,900 runs, three to five weeks of the
machine; the stages above select the configurations which may win on each
side.

## Single precision

The float build (`deps/spack/mmm-trilinos-float`) has hypre, PETSc and
MFEM built in single precision, Trilinos with the float
instantiations of Tpetra, MUMPS (its single precision library), cuDSS,
Ginkgo and oneMKL. Not available in float: SuperLU_DIST, STRUMPACK and
SuiteSparse (the wrappers of MFEM are written for double), SUNDIALS
(refused by the makefile of MFEM in single precision; not used by the
benchmark), GSLIB, AmgX
(refused by mfem-mgis in a float build; it is measured in double, also
computing in single precision within the double build), Teko (double
only). The relative tolerances of the linear solvers are at least 1e-6.

## Trials at the small size (16/10/2026, previous driver)

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
- ex1 with Newton: the time step of ex1 is the inverse of the number of
  steps, and the driver gave 5 steps instead of 100. With 20 steps or
  less, the Newton solver does not converge at the onset of plasticity.
  The driver now gives 25 steps, and ex1 passes on the CPU with all the
  solvers tried (to be checked again on the GPU);
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

- Overrides per problem: a single linear solver override replaces the
  solvers of all the problems of a program (ex8 uses its own options
  instead).
- Whole program time or time of the resolution only (the whole program is
  measured).
