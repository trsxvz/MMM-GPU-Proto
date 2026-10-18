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

## Second pass

The first pass (17/10/2026, `results/benchmark-slides-2026-10-17`,
`slides/v1`) showed several weaknesses, addressed by the second one
(`results/benchmark-slides-v2`, `slides/order-1`, `slides/order-2`):

- Orders 1 and 2 for all the cases, with sizes of about one million
  unknowns at both orders (one refinement less at the order 2) where the
  memory of the GPU allows it (see the sizes below). The first pass filled
  the GPU poorly at the order 1 for ex7, ex8 and the polycrystal.
- The ParaView exports of ex2, ex3, ex4, ex8 and of the polycrystal were
  timed with the computations: they are disabled (options added to ex2,
  ex3 and the bubble). The quantity of interest of each case is written in
  any case.
- Results checked: the quantity of interest of each run (resultant force,
  mean stresses, reaction, error to the analytical solution, statistics
  of the fields of ex8) is compared to the one of the solvers of the
  example on the CPU, or of the fastest configuration on the CPU when they
  did not finish (relative tolerance 1e-5 in double, 1e-3 in float). The
  runs in float are compared to the reference in double when it is
  measured: they differ beyond 1e-3 and beyond three times the median of
  the best differences of the configurations in float of the case, and
  always beyond 5 %, the accuracy reached in float depending on the case
  (at the measurement size, the crystal plasticity of the polycrystal and
  of the cermet stays 0.8 to 4 % off double at the order 1 and 0.5 to
  1.6 % at the order 2 with all the solvers, the solvers of the programs
  included; 0.15 to 0.8 % at the screening size). The solvers of a program may give
  wrong results in float: on ex8, MUMPS in single precision gives a norm of
  the displacement 9 % off the one of all the solvers in double (which
  agree to 1e-12; its maximum is 8.7 % below, its mean 11 %), while
  NewtonK0 and the conjugate gradient of hypre in float agree with double
  (at the order 2, only the solvers using MUMPS ran in float, 13 % off).
  No speedup is then given against the solvers of the program, and the
  best configurations are taken among the runs which agree. The values are compared by groups (the
  components of a tensor, the minimum, maximum and mean values of a field),
  relatively to the norm of their group in the reference, and the
  deformation gradients of the polycrystal and of the cermet by their
  strains (F - 1). The runs whose
  results differ keep their bar in the graphs, hatched, with the
  difference, and are listed in the indexes; the runs without a quantity
  of interest are dotted and listed as not checked. The output files of
  each run are kept in `outputs/<run>/`, a lock of the directory of the
  program preventing two drivers from exchanging them. The end of a
  program is waited for without polling: the times are not rounded.
- Analysis of each run from its log, without any additional run: the
  number of non linear iterations, the number of resolutions (time steps,
  sub-steps, iterations of the coupling scheme of ex8) and the times of the
  phases given by the time table of the profiler. Each section of the
  profiler gives its own time (without the ones of its sub-sections) to
  its phase (`PHASE_OF_SECTION` in the driver), so that the whole time of
  the profiler is counted once:
  - initialization: mesh, finite element spaces, boundary conditions,
    behaviours and solvers of the problems, and the time before the start
    of the profiler (initialization of MPI, of the device and of the
    libraries; up to 13 s of the 27 s of the polycrystal in float on the
    GPU);
  - integration of the behaviours (in the estimation of the new unknowns,
    in the fixed point methods and in the prediction), inner forces and
    update of the state at the end of each time step;
  - tangent operator and its assembly, with the elimination of the
    boundary conditions (also the assembly of the prediction, not timed by
    a sub-section in the coupled organization before patch `0273`);
  - setup of the linear solvers by `SetOperator` (factorizations,
    hierarchies of AMG outside hypre: MFEM sets the solvers of hypre up
    at their first resolution, counted in the linear solves);
  - linear solves, with the ones of the prediction (one per time step,
    each costing about one Newton iteration);
  - post-processings of mfem-mgis and of the programs;
  - the rest: the own time of the non linear solvers (with the Anderson
    acceleration), of the time loops, of the coupling scheme of ex8 and of
    the code of the programs outside the sections (up to 2 % at the
    measurement size, 4 % at the screening size).

  The first pass of the analysis left the prediction, the fixed point
  methods and the post-processings of the programs out of the phases
  (most of the time of ex1 in gray): `benchmark_slides.py reanalyze
  RESULTS...` computes the iterations, resolutions and phases of the runs
  again from their logs (`logs/`, or `logs.tar.xz` when `logs/` is
  absent), without any run, once no driver writes in the directory (a
  `runs.csv` modified during the analysis is not replaced). The profiler
  is enabled in all the programs (option `--profiling` added to ex1, ex2,
  ex3 and ex5); on the GPU, its sections synchronize the device.
- Organization of the computations on the CPU: the solvers of the program
  are also run without the option `--execution` (behaviours integrated
  point by point and computations element by element, the integration and
  the assembly being coupled in the loops over the elements, as in the
  master branch of mfem-mgis), against the organization of the prototype
  (`--execution cpu`: integration by batches and kernels).
- Mixed placement on 8 processes (the linear solvers compute on the 8
  cores): on one process, its solvers could not compete. Each process
  sharing the GPU reserves the memory of the integration of its
  behaviours: 6 or 8 processes exhaust the 16 GiB of the GPU for the
  polycrystal and the cermet (grain law of UO2) at both orders, run on 4
  processes, and for ex6 at the order 1 (screened on 8 and 4 processes,
  8 fitting at the order 2). The other cases fit on 8 processes at both
  orders (checked at the screening size). The best GPU
  configuration is also tried on 2 and 4 processes sharing the GPU
  (stopped beyond 1.5 times the best time).
- More configurations, all screened, a configuration being stopped beyond
  three times the best time of its side: the AMG built once per time step
  (modified Newton method, also accelerated by Anderson) or every 5
  iterations (assembled operator of the tangent operator applied without
  assembly), the fixed point on the elastic operator accelerated by
  Anderson, the aggregation multigrid of AmgX with aggregates of 8
  unknowns or with the multicolour Gauss-Seidel smoother, the multigrid
  (PGM) of Ginkgo, GAMG of PETSc with Chebyshev and Jacobi smoothers (CPU
  and mixed placement), Tacho (its parameters were given outside of its
  sublist in the first pass, where it always failed). The options of the
  solvers suited to the GPU were chosen by the trials below.
- PETSc computes on the GPU (`patches/mfem-mgis/0272-*`; the defect of
  MFEM it works around is described in `reproducers/petsc-gpu-vectors`).
  Three configurations use it: GAMG with the Krylov solver of PETSc
  (`gpu-petsc-gamg`), and NewtonK0 with one cycle of GAMG on the elastic
  operator, with or without the tangent operator assembled
  (`gpu-k0-petsc`, `gpu-k0-petsc-mf`). With the options of the CPU, the
  product PtAP of cuSPARSE exceeds the memory of the GPU at 786,432
  unknowns. Built for each resolution, the hierarchy takes most of the
  time (its setup is partly on the host). On the GPU, GAMG thus uses the
  plain aggregation, the hierarchy being built for the first resolution
  and reused by the next ones (trials below). Some failures are expected
  in the screening: out of memory on ex6 and the polycrystal (except
  `gpu-k0-petsc-mf` on the polycrystal), no convergence of GMRES on the
  cermet (also on the CPU) and on ex4 in float. The checks of these
  configurations through the driver are in
  `results/benchmark-slides-v2-petsc` (screening size, order 1).
- The final runs try the next configuration of the screening when one
  fails at the measurement size (in float, the selected configurations of
  ex2 and ex7 did not converge over the whole loading in the first pass),
  and measure first the configurations whose results agree. In float, the
  solvers of the program on the CPU (in both organizations) and the 5 best
  configurations of each side in double, at the same order, are
  screened.
- In float, ex3 and ex6 compare their solution to the analytical one with
  a tolerance calibrated on the coarse meshes of their tests (errors of
  7.9e-4 and 9.9e-4 with 786,432 unknowns): the comparison is disabled
  (`--no-check`), the error being printed and compared to the one of the
  solvers of the example on the CPU. The polycrystal runs 3 steps up to
  0.6 s (its grain law of UO2 fails in float beyond about 0.8 s).

Sizes of the second pass (the screening size has the mesh of the
measurement size, with one or two time steps):

| Case | Order 1 | Order 2 | Unknowns (order 1, order 2) |
|---|---|---|---|
| ex1 | -r 6 (screened with -r 5) | -r 5 (-r 4) | 823,875 at both orders |
| ex2 | -r 4 | -r 3 | 908,226 |
| ex3, ex6 | -r 4 | -r 3 | 786,432 |
| ex4 | -r 1 | default mesh | 196,047; printed by the program |
| ex5 | -r 5 | -r 4 | 1,847,042 |
| ex7 | -r 3 | -r 3 | printed by the program |
| ex8 | -r 1 | default mesh | 214,812 (mechanics) |
| polycrystal | -r 2 | -r 1 | 1,360,002 |
| bubble | -r 1 | default mesh | printed by the program |
| cermet | -r 1 | default mesh | 269,253 |

One refinement more for ex4, ex7 (order 1) and the cermet does not fit in
the memory of the GPU, even with the solvers of the programs: the
integration by batches keeps the tangent operator of each integration
point on the GPU (out of memory in `cudaMalloc` at the screening size).
The polycrystal at the order 1 fits with BoomerAMG, not with MueLu.

Trials of the options of the GPU solvers (17/10/2026, screening size
unless noted, `results/benchmark-slides-v2-tune`, with the scripts of the
trials; time of the resolutions, from the profiler, and total time). The
functional tests of all the configurations at the small size are in
`results/benchmark-slides-v2-smoke`, the checks of the memory of the GPU
in `results/benchmark-slides-v2-fit`:

| Solver | Trial | Result | Decision |
|---|---|---|---|
| BoomerAMG | aggressive coarsening (1 level) | ex2: resolutions 9.5 s to 4.8 s, total 11.7 s to 7.4 s; polycrystal: 1.85 s to 1.46 s; ex3, ex5: as fast; cermet at -r 2: the only one fitting in memory | kept (`hypre-amg-agg`) |
| BoomerAMG | Chebyshev, l1-Jacobi, threshold 0.5, interpolation 14, 6 elements per row, aggressive interpolations 5 and 7, transpose kept, and Chebyshev with aggressive coarsening and transpose kept | as fast as the default options or slower (up to +60 % for the resolutions) | not kept |
| Ginkgo, Jacobi | scalar (blocks of 1) against the default blocks (up to 32) | ex3: 396 iterations and 0.52 s against 652 and 0.99 s | blocks of 1 |
| Ginkgo, multigrid | smoother of two Jacobi iterations relaxed by 0.5 (default: one relaxed by 0.9) | ex3 at the small size: 61 iterations, 3.3 s, against 113 s; converges on the polycrystal; a cycle stays costly (ex3: 233 iterations, 34 s) | kept |
| Ginkgo, multigrid | W-cycle; K-cycle (suited to the pairwise aggregation) | 3.5 times slower; refused by the configuration of Ginkgo 1.11 | not kept |
| Ginkgo | incomplete factorizations (IC, ILU, with ISAI) | not a number from the first iteration on ex3 (breakdown of the factorization); IC with exact triangular solves: no convergence in 120 s | not screened |
| Ginkgo | CB-GMRES (Krylov basis compressed) | refused by mfem-mgis in parallel | not available |
| AmgX | aggregates of 8 | ex3: 0.65 s | kept (`amgx-agg8`) |
| AmgX | multicolour DILU smoother | fails on 5 cases of 11 at the small size | replaced by the aggregates of 8 |
| AmgX | l1-Jacobi, Chebyshev polynomial, multicolour ILU smoothers | timeout on the polycrystal; failure; "Configuration feature is not implemented" | not kept |
| MueLu | Chebyshev of degree 3, dropping tolerance 0.02 | no gain on ex3 and ex5 | not kept |
| hypre | FSAI, ILU(0) by blocks | ex3: 6 and 12 times slower than BoomerAMG | not kept |
| GAMG of PETSc | options of the CPU (smoothed aggregation, aggressive coarsening by the square of the graph) | ex3, polycrystal: out of memory in the product PtAP of cuSPARSE; ex1: 64.9 s (solvers of the program: 7.4 s) | replaced |
| GAMG of PETSc | prolongators kept from one resolution to the next (`pc_gamg_reuse_interpolation`, already the default of PETSc 3.26); no aggressive coarsening | ex1: 65.2 s, no effect (the square of the graph is computed again at each of the 41 setups: mfem-mgis gives PETSc a new matrix at each resolution); ex3: out of memory | not kept |
| GAMG of PETSc | plain aggregation (`pc_gamg_agg_nsmooths: 0`) | ex3: 17.2 s twice, the only options fitting in memory (BoomerAMG: 5.5 s); ex1: 63.9 s, 135 to 204 iterations per resolution, the 41 setups taking 72 % of the time and the square of the graph (on the host) 60 % (profile of PETSc) | kept |
| GAMG of PETSc | aggressive coarsening by MIS-2 instead of the square of the graph | ex1: 28.4 s, 29.9 s (plain aggregation), 19.7 s and 22.8 s (with the preconditioner reused), but four runs stopped at 300 s, the numbers of iterations of the linear solves being unchanged | not kept (not repeatable) |
| GAMG of PETSc | preconditioner built for the first resolution and reused by the next ones (`ksp_reuse_preconditioner`) | ex1: 20.1 and 20.2 s (smoothed aggregation), 20.3 s (plain aggregation) | kept with the plain aggregation (`petsc-gamg-gpu`) |
| GAMG of PETSc | ex6 (786,432 unknowns, non linear): plain aggregation with the preconditioner reused, NewtonK0 (smoothed or plain aggregation); polycrystal (1,360,002 unknowns): the options of the CPU, plain aggregation (with or without the preconditioner reused), NewtonK0 (smoothed aggregation) | out of memory in all the runs | screened, failures expected |
| GAMG of PETSc | non symmetric cases (GMRES): cermet, ex4 in float | cermet: first resolution stopped at 5,000 iterations without converging, 2,825 iterations after the time step was halved, run stopped at 300 s, also on the CPU; ex4 in float: no convergence (ex4 in double: 22.4 s) | screened, failures expected |
| NewtonK0 with GAMG | tangent operator applied without assembly | polycrystal: 47.0 s (fits in memory; NewtonK0 with BoomerAMG: 22.5 s); ex6: out of memory; ex1: 17.1 s | kept (`gpu-k0-petsc-mf`) |
| NewtonK0 with GAMG | one cycle of GAMG (Krylov solver of PETSc without iterations, `preonly`: the preconditioners of PETSc are not linear solvers of mfem-mgis), smoothed or plain aggregation | ex1: 12.0 s in both cases (one cycle of BoomerAMG: 9.4 s); ex4: 9.0 s, cermet: 8.7 s; the results differ from the ones of the CPU on ex4 (2.5e-4 in double) and on the cermet (1.2e-5), as with NewtonK0 and the other libraries on ex4 (AmgX: 1.8e-4, BoomerAMG: 2.7e-5 at the small size; 2e-6 to 5e-6 on the cermet) | kept with the plain aggregation (`k0-petsc`) |

Not available, with the reason:

| Solver | Reason | Kind |
|---|---|---|
| Elasticity options of BoomerAMG on the GPU | not supported by hypre on the GPU; MFEM aborts (`linalg/hypre.cpp`) | limit of hypre |
| Elasticity options of BoomerAMG and rigid body modes of GAMG on the CPU | require the ordering `byVDIM` of the unknowns, the programs use `byNODES` | not done: option of the programs |
| GAMG of PETSc on the GPU for the large non linear cases | out of memory on ex6, and on the polycrystal except with NewtonK0 and the tangent operator applied without assembly (products of cuSPARSE, operator copied to the format of PETSc in addition to the one of hypre) | limit of the backend CUDA of PETSc here; its backend Kokkos is not supported by mfem-mgis yet |
| Anderson acceleration of NewtonK0 | refused by the fixed point solver of mfem-mgis (an acceleration of a Newton method) | choice of mfem-mgis |
| Ginkgo, SUNDIALS, SuiteSparse, SuperLU_DIST, STRUMPACK, GSLIB, libCEED in float | refused by the makefile of MFEM in single precision | limit of MFEM |
| AmgX in float | refused by mfem-mgis in a float build (measured computing in single precision in the double build) | choice of mfem-mgis |
| Tacho in float | Tacho of Trilinos 17.2.1 is only instantiated in double | limit of Trilinos |
| MUMPS on the GPU | MUMPS computes on the CPU only (cuDSS is the direct solver of the GPU) | limit of MUMPS |

Measured differences between double and float to keep in mind: the
tolerances of mfem-mgis have floors in float (about 1e-5 relative), so
that a float run may make fewer iterations than a double one (ex8, first
pass, solvers of the example on the CPU: 20 Newton iterations and 28
resolutions in float against 62 and 48 in double, the factorizations of
MUMPS taking 13.6 s against 93 s), besides the faster arithmetic; the
iterations are given next to the times.

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
  screening size (the mesh of the measurements and one or two steps of
  the same time step; ex1 and ex7, without an option of the end time, are
  screened on a coarser mesh or with fewer steps, see Sizes); the runs
  already made are skipped, so that an interrupted screening resumes;
- `finals RESULTS CASE`: at the measurement size, the configurations
  selected by the screening (the solvers of the program on the CPU with
  their best number of processes, the best other CPU configuration, the
  two best configurations with the GPU), all but the second GPU one being
  repeated (three runs) unless their first run took more than two
  minutes; the configurations already measured are skipped, so that
  `finals` can be run again after a new screening;
- `best RESULTS`: the best times and the speedups;
- `reanalyze RESULTS...`: the numbers of iterations and of resolutions
  and the times of the phases of the runs, computed again from their logs
  (`logs/`, or `logs.tar.xz` when `logs/` is absent), without any run and
  once no driver writes in the directory.

Common options: `--precision double|float`, `--size small|screen|bench`,
`--timeout S` (a run is stopped with all its processes, `mpirun` and its
ranks, after S seconds; its status is then `timeout`, and its time a lower
bound), `--gpu-lock FILE`, `--dry-run`.

`cases/benchmark_plots.py OUTPUT RESULTS...` draws, for each case, size and
precision, the best time with the solvers of the program on the CPU, the
best CPU configuration and the best configurations with the GPU, with
their speedups, and two overviews of the speedups (against the solvers of
the examples and against the best CPU configuration), with an index
`README.md`. The graphs are generated in `slides/` on the machine of the
measurements, and not versioned: they are drawn again from the `runs.csv`
of the results. With `--agreeing-only`, the graphs only show the runs whose
results agree with the reference, without any hatched bar (the lower
bounds of the runs stopped by the time limit are given by the texts); with
`--by-precision`, the graphs of the cases of each order are in `double/`
and `float/`, the overviews and the index in the directory of the order. The final graphs
of the second pass are drawn with both options (`slides/final-fixed`).

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
inexact Newton method while it still converges). With the GPU, hypre
allocates at once a pool of 4 GiB, which left too little memory to
Trilinos, AmgX and cuDSS on ex3, ex4 and ex6: with the solvers of the
other libraries (Trilinos, AmgX, cuDSS, Ginkgo, PETSc), the initial pool
of hypre is 1 GiB (`MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`; it
grows if needed).

- CPU: the solvers of the program on 8 and 16 processes; on 8 processes,
  the Krylov solver of hypre with BoomerAMG, Belos with MueLu (Gauss-Seidel
  smoother, hierarchy reused), PETSc with GAMG, MUMPS, and for the non
  linear computations the inexact Newton method and NewtonK0 (elastic
  operator by one cycle of BoomerAMG or of MueLu).
- GPU, tier 1 (screened first): the solvers of the program when they
  compute on the GPU, hypre with BoomerAMG, Belos with MueLu (Jacobi or
  Chebyshev smoother), AmgX (classical or aggregation multigrid), cuDSS,
  MUMPS on the CPU with the elements and the behaviours on the GPU, the solvers of the program in this mixed
  placement, the inexact Newton method, NewtonK0 with its elastic operator
  by BoomerAMG or MueLu, with or without the tangent operator applied
  without assembly, NewtonK0 with the elastic operator factorized by cuDSS,
  and the conjugate gradient with a Jacobi preconditioner on the tangent
  operator applied without assembly; GAMG of PETSc (plain aggregation,
  hierarchy built once) and NewtonK0 with one cycle of GAMG, with or
  without the tangent operator assembled.
- GPU, tier 2 (variants dominated in the measurements of the solvers,
  `docs/solver-recommendations.md`): aggressive coarsening of BoomerAMG,
  Krylov solvers of MFEM with BoomerAMG, AmgX in single precision within
  the double build, the Jacobi preconditioner on the assembled operator, Ginkgo, Tacho, NewtonK0 with
  AmgX, the modified Newton method with cuDSS.

Left out for some cases: the direct solvers (MUMPS on the CPU, cuDSS)
for the large 3D cases (ex3, ex6, the bubble; cuDSS for ex1, MUMPS for the
cermet), MUMPS having saturated the memory of the machine on ex3; the solvers of the
program on the GPU when they are refused there (MUMPS of ex2 and ex5) and
in the mixed placement when they are refused on the CPU with the rest on
the GPU (conjugate gradient of MFEM of ex1 and ex3).

Not used: the sparse solvers of MAGMA (sequential only, the programs are
parallel), AmgX in a float build (refused by mfem-mgis), Ginkgo,
SuperLU_DIST, STRUMPACK and SuiteSparse in float (refused by MFEM).

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
| ex7 | -r 2 | 2 | printed by the program | 40 |
| ex8 | -r 1 | 1 | 214,812 (mechanics) | 5 |
| polycrystal | -r 1 | 1 | 168,885 | 4 steps of 0.2 s |
| ex3-o2 | -r 3 | 2 | 786,432 | 1 solve |
| ex6-o2 | -r 3 | 2 | 786,432 | 1 solve |
| bubble | default mesh | 2 | 467,109 | 1 solve |
| cermet | -r 1 | 1 | 269,253 | 2 steps of 0.2 s |
| polycrystal-r2 | -r 2 | 1 | 1,360,002 | 4 steps of 0.2 s |
| polycrystal-o2 | -r 1 | 2 | 1,360,002 | 4 steps of 0.2 s |
| cermet-r2 | -r 2 | 1 | 2,164,290 | 2 steps of 0.2 s |

- ex1: by default, the increment of the imposed displacement is
  concentrated in the elements next to the loaded face. With 4 refinements
  or more, their strain increment is far beyond the elastic limit and the
  Newton solver does not converge when the loading is reversed (with 6
  refinements, not even at the first step). The option `--prediction` of
  the example (elastic prediction, as in ex2 and ex4) removes the failure:
  25 predictions and 16 Newton iterations for the 25 steps. ex1 has no
  option of the end time: it is screened with 5 refinements (107,811
  unknowns).
- ex7: with 1 refinement, it has only 21,012 unknowns, which do not fill
  the GPU: it is measured with 2 refinements, and screened with 10 steps
  instead of 40 (no option of the end time).
- ex3-o2 and ex6-o2: ex3 and ex6 with quadratic elements and one
  refinement less, for the same number of unknowns: more work per unknown
  in the kernels of the elements and denser matrices.
- polycrystal-r2, polycrystal-o2 and cermet-r2: the cases of mm-opera-hpc
  on larger meshes, and the polycrystal with quadratic elements for the
  same number of unknowns, without the direct solvers (memory).
- The polycrystal and the cermet run for short durations (the same in
  double and in float): their grain law of UO2 fails in float beyond about
  0.8 s (rounding of its elastic strain when the deformation gradient is
  close to the identity).

## Plan

Second pass: four exclusive windows, in this order, so that the graphs
of the order 1 are complete in both precisions first: double at the order
1, float at the order 1 (`screen --reference
results/benchmark-slides-v2/double --keep 5`), double at the order 2,
float at the order 2. In each window, for each case: `screen` of the GPU
and mixed configurations, `screen` of the CPU configurations, `finals`,
then the graphs (`slides/order-1`, `slides/order-2`). Estimated
duration: 6 to 8 hours for a window in double (PETSc on the GPU adds
about 30 minutes), 2 to 3 hours in float.

First pass:

| Stage | Where | Content |
|---|---|---|
| Calibration | GPU shared, CPU exclusive | each case at the measurement size, solvers of the program |
| Screening, GPU | GPU shared | tier 1 of each case at the screening size |
| Screening, CPU | exclusive | the CPU configurations at the screening size |
| Final runs | exclusive | `finals` of each case, graphs after each case |
| Float | same stages | once the float build is installed; screening of the solvers of the program and of the 3 best configurations of each side in double (`screen --reference`) |
| Refinement | as available | tier 2, variants of the winners |

The full grid (all the solvers, processes, orders and precisions, three
repetitions) would take about 9,900 runs, three to five weeks of the
machine; the stages above select the configurations which may win on each
side.

## Validation of the second pass

All the runs of `results/benchmark-slides-v2` were checked after the
campaign (18/10/2026), without any additional run: 1,088 runs in double
and 561 in float, at the orders 1 and 2.

- Quantities of interest in double: all the final runs agree with the
  solvers of the programs (1e-5), at both orders. At the order 2 of ex1,
  the solvers of the program (conjugate gradient without preconditioner)
  did not finish within 1,800 s: the runs are checked against the fastest
  CPU configuration (BoomerAMG), and agree within 7e-13.
- Quantities of interest in float, against double:
  - ex8: the solvers of the program (MUMPS in single precision) are 9 %
    (order 1) and 13 % (order 2) off; at the order 1, NewtonK0 and the
    conjugate gradient of hypre agree, and give the speedup (x1.6).
  - ex7 at the order 1: the GPU runs on one process are 1.8e-3 off (the
    solvers of the program: 1.2e-4); the run on 2 processes sharing the
    GPU agrees and gives the speedup.
  - polycrystal and cermet: all the runs are 0.5 to 4 % off (crystal
    plasticity in float), the solvers of the programs included, within
    the accuracy reached in float on these cases.
  - ex7 at the order 2: no configuration converged at the screening size;
    at the measurement size, only the solvers of the program finished, on
    the CPU (1,315 s, coupled organization 1,726 s): no graph.
- Repeatability: the times of the repeated runs are within 10 % for most
  configurations. Some runs are not deterministic: one more non linear
  iteration or sub-step in some repetitions (ex5 at the order 2 with
  MUMPS: 4 or 5 Newton iterations, 19 to 24 s; NewtonK0 on the GPU, one
  iteration), more often in float (ex8 at the order 1: 39 to 41
  iterations, 40 to 54 s). The first run of a configuration may be slower
  (polycrystal at the order 2: 40 s, then 29 and 30 s). The final runs
  keep the best time of each configuration.
- Final runs which failed: the runs with several processes sharing the
  GPU, stopped at 1.5 times the best time or out of memory, the solvers
  of the program at the order 2 of ex1 (time limit), the solvers of
  NewtonK0 in float on ex2 (no convergence at the measurement size), and
  the solvers of the program of ex4 at the order 2 in float.
- Phases: the time left out of the phases (rest) is below 2 % of every
  final run.

## Single precision

The float build (`deps/spack/mmm-trilinos-float`) has hypre, PETSc and
MFEM built in single precision, Trilinos with the float
instantiations of Tpetra, MUMPS (its single precision library), cuDSS
and oneMKL. Not available in float: the packages refused by the makefile
of MFEM in single precision (SuperLU_DIST, STRUMPACK, SuiteSparse,
SUNDIALS, Ginkgo, GSLIB, libCEED, among others), AmgX (refused by
mfem-mgis in a float build; it is measured in double, also computing in
single precision within the double build), Teko (double only), Tacho (only
instantiated in double: Amesos2 of the float Trilinos refers to its float
instantiations, provided to the linker by a library of stubs which abort,
`deps/spack/mmm-trilinos-float/tacho-float-stubs.sh`). MGIS is built with
`MGIS_REAL_TYPE=float` (`install/h-interface/mgis-float`), and the driver
puts the library of MGIS of the precision of each run first in
`LD_LIBRARY_PATH`. In float, ptxas exceeds 15 GB on `fem/lor/lor_batched.cpp`
of MFEM, compiled as PTX only (`reproducers/ptxas-float-lor-batched`). The
relative tolerances of the linear solvers are at least 1e-6.
The tests of the solvers of PETSc of mfem-mgis (`PETScSolversTest`, linked
with the library of stubs of Tacho) fail in float on the CPU and on the GPU
alike, their tolerances being calibrated in double: the numbers of
iterations are the same on the CPU and on the GPU, up to two iterations of
GAMG.

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
