# Polycrystal on Adastra: convergence, solvers and series in number of grains

These are the Slurm jobs of the measurements of 13/11/2026 for the "fait
marquant" (see `JOURNAL.md`). They ran on the MI300A partition, from
`/lus/scratch/RES1/genden15/tchenaille`, the directory where they were
submitted (`SLURM_SUBMIT_DIR`).

The quota of files of the project on `/lus/work` was reached that day. The
results, the meshes and the versions of the scripts are therefore on
`/lus/scratch`.

## Common environment

- `common.sh`: environment of the prototype, results directory
  `fm/$FM_RUN`, journal `chain.log`, and the function `energy` (`sacct`).
  For the jobs which build or run tests: `load_build_environment`
  (compilers and `mpirun` of the machine profile), `behaviour_targets` and
  `test_settings` (variables of the runs of the tests by `ctest`: host or
  GPU, kernels forced or not).
- `deploy-overlay.sh <commit> <archive>`: a version of the scripts of
  `cases/`, with links to the rest of the deployed prototype.
  - The jobs take it with `FM_PROTO`.
  - This avoids overwriting the scripts of a running job: bash reads a script
    while it runs it.

## The jobs

All of them run `cases/measure-polycrystal-grains.sh`.

- **`conv.slurm`: convergence in mesh.**
  - 100 grains, P2, refinement `FM_REFINEMENT`, GPUs only.
  - `FM_MESHES=.../meshes/clean` for the cleaned meshes.
  - `FM_TAG`, `FM_NSTEPS` and `FM_DURATION` set the time steps, for instance
    20 steps of 0.25 s.
- **`solvers.slurm`: choice of the solvers.**
  - 100 grains cleaned, r1, 10 steps of 1 s.
  - Configurations `FM_CONFIGS` on the GPUs or on the cores
    (`FM_MODE=gpu` or `cpu`).
- **`memtest.slurm`: reading of a large mesh** by 96 processes on one node,
  giving the memory and the time of `Mesh::Constructor`.
- **`series.slurm`: series in number of grains.**
  - `FM_GRAINS` grains on the nodes of the allocation, refinement 1,
    `FM_NSTEPS` steps of 1 s.
  - The best solvers of each side: `FM_SOLVER_GPU`, `FM_SOLVER_CPU`.
  - `FM_GPU_REPS` computations on all the GPUs, then `FM_CPU_REPS` on all
    the cores.
- **`partition.slurm`: partitioned meshes** (`meshes/partition/`, see
  `meshes/README.md`).
  - For each `grains:processes` of `FM_PARTITIONS`, the cleaned mesh
    refined once and split in parts, in
    `meshes/part/v<grains>-np<processes>/mesh.<rank>`. The computations
    read them with `--mesh_mode Restart`, `FM_SERIES_REFINEMENT=0` and
    `FM_MESH_GPU` or `FM_MESH_CPU` of `series.slurm`.
  - By default `partition-mesh`, one process per part
    (`FM_PARTITION_TASKS_PER_NODE` processes per node: each one reads the
    whole mesh).
  - `FM_PARTITION_WORKERS=<w>`: `split-mesh` with `<w>` processes for any
    number of parts. One node is enough for the 6144 parts of 38,400
    grains.
  - Each partition is then read back by `split-mesh --check` on all the
    cores, which writes the file `checked` in its directory when its parts
    are consistent. Only a partition with this file is kept by the next
    jobs; the job fails if a partitioning or a check fails.
  - `FM_SPLIT_TEST=1`: tests of `split-mesh` first (400 grains, 96 parts,
    on the local disk of the node); the job stops if they fail.
  - `FM_TEST=1`: computations reading the parts against the ones building
    the mesh themselves.
  - `FM_REMOVE_PARTITIONS`: partitions removed first, the parts of
    computations already done. `FM_FILES_LIMIT`: a partition is not
    written, and the job fails, if its parts would bring the files of the
    group on the scratch above this number (quota read by `lfs`).
- **`geo.slurm`, `mesh.slurm`: meshes made on a node** (see
  `meshes/README.md`).
  - `geo.slurm`: geometries of Merope of `FM_GRAINS` grains for the seeds
    `FM_SEEDS`, in parallel.
  - `mesh.slurm`: correction of the plane surfaces, gmsh and cleaning. With
    `FM_SEEDS`, the geometries of `geo.slurm` are tried in this order until
    gmsh meshes one without error (the seed kept is in the file `seed`);
    the job fails without a cleaned mesh.
- **`ls-tests.slurm`: tests of the line search** (option `LineSearch` of
  mfem-mgis, 14/11/2026).
  - Separate build in `FM_FIX` (`build-fix.sh`), then the unit tests of
    mfem-mgis on the line search (`FM_TESTS`) and other tests which use
    the non linear solvers (`FM_OTHER_TESTS`).
  - Then the commands of the tests of the polycrystal
    (`polycrystal/CMakeLists.txt` of mm-opera-hpc), run by hand with the
    behaviour library of the deployed prototype.
- **`ls-identity.slurm`: options disabled, same outputs.**
  - CPU, 100 grains, 96 processes, one material for all the grains:
    executable `FM_EXECUTABLE_B` against `FM_EXECUTABLE_A`, cases
    `FM_CASES` (NewtonK0 with the two choices of Eisenstat and Walker,
    Newton solver, raw mesh with integration failures).
  - Compares the `.res` files, the averages per grain and the lines of the
    iterations of the logs.
  - `FM_VARIANTS=B` with `FM_REFERENCE_DIR` runs the second executable
    alone and compares it to the outputs of a previous job (the raw mesh
    takes 30 minutes per computation); `FM_TIMEOUT_CPU` is the time limit of
    a computation.
  - `FM_EXTRA_OPTIONS` replaces `--single-material`: empty for one material
    per grain, or with the options of the measurements (line search,
    prediction) to compare two executables with these options enabled.
- **`ls-stats.py <log>...`: counters of a computation**: attempts of time
  steps by time increment, Newton corrections, iterations whose step was
  reduced by the line search, refused trials, failed searches.
  - For the first resolution of each attempt, which starts from the
    prediction: corrections, steps reduced at the first correction, and
    ratio `||r||/||r_0||` at the iteration 0 (below 1 when the prediction
    by linear extrapolation is better than the default one).
  - Number of extrapolated estimates discarded because the behaviours
    could not be integrated at them (the resolution then starts again from
    the default prediction).
- **`ls-gpu-tests.slurm`: tests of the non linear solvers on a GPU.**
  - Tests of mfem-mgis which use the non linear solvers (`FM_TESTS`), with
    the mfem-mgis of the separate build of `FM_FIX`, on one MI300A.
  - Settings of `scripts/mfem-mgis-gpu-suite.sh`, one test at a time:
    kernels forced, computations selected by each test, then the tests
    whose name selects a device (`FM_MODES`).
- **`pred-tests.slurm`: tests of the prediction by linear extrapolation**
  (strategy `LINEAR_EXTRAPOLATION_PREDICTION` of mfem-mgis, option
  `--prediction` of the polycrystal, 14/11/2026).
  - Separate build in `FM_FIX`, unit tests of mfem-mgis.
  - Then the polycrystal of the tests of mm-opera-hpc in four time steps,
    with and without the extrapolation: axial stresses, lines of
    iterations, ratios at the iteration 0.
- **`pred-mutations.slurm`: mutations of the prediction by linear
  extrapolation.**
  - Separate build in `FM_FIX`, which must not be the one of an executable
    used by measurements: the library built last is a mutated one.
  - Eight mutations of the extrapolation and ten of its fallback on the
    default prediction, applied one at a time to the sources of mfem-mgis:
    each one must make `LinearExtrapolationPredictionTest` fail in at
    least one of its variants (sequential, 2 processes, batches and
    kernels on the host).
- **`pred-raw.slurm`: the extrapolation on the raw mesh of 100 grains**
  (5 steps of 1 s on all the cores, the case `i3-raw` of
  `ls-identity.slurm`, whose behaviour fails to integrate at some
  estimates).
  - Options `FM_OPTIONS` (default: the line search and the extrapolation),
    executable `FM_EXECUTABLE`, name of the computation `FM_LABEL`.
  - Logged: the line of the computation, the extrapolated estimates
    discarded, the estimates moved after an integration failure, the failed
    resolutions.
- **`review-tests.slurm`: tests added after the review of the patches**
  (14/11/2026).
  - Separate build in `FM_FIX`, then the unit tests of mfem-mgis on the
    weights of the integration points kept by the problems, on the kept
    rotation matrices and on the failure of the line search of NewtonK0
    (`FM_TESTS`), and other tests using the materials and the non linear
    solvers (`FM_OTHER_TESTS`).
  - Then the tests of the polycrystal of mm-opera-hpc by `ctest`, with the
    behaviours built there, and the comparison of `--single-material` to
    one material per grain run by hand (`cmake/run-and-compare.py` of
    mm-opera-hpc), on one process and on two.
- **`review-mutations.slurm`: mutations of the code covered by these
  tests.**
  - Separate build in `FM_FIX`, which must not be the one of an executable
    used by measurements.
  - Eleven mutations, one at a time, each of which must make its test
    fail: four of the kept weights, three of the kept rotation matrices,
    two of the line search of NewtonK0, two of the single material of the
    polycrystal. `FM_MUTATIONS` selects them.
- **`pred-ab.slurm`: the extrapolation measured on an RVE.**
  - `FM_GRAINS` grains, one material for all the grains, 20 s simulated,
    options `FM_COMMON` (the line search), without and with
    `--prediction LinearExtrapolation`, one after the other at each
    repetition, on the same nodes.
  - `FM_GPU_REPS` repetitions on all the GPUs, `FM_CPU_REPS` on all the
    cores; `FM_VARIANTS` selects the variants (`ls`, `ls-pred`).
- **`profile-trace.slurm`: trace of the final computation**, in the style
  of Nsight Systems (14/11/2026).
  - `FM_GRAINS` grains on all the GPUs of the allocation, final
    configuration (one material, line search, extrapolation, mesh
    partitioned beforehand), under `rocprofv3 --runtime-trace
    --marker-trace --output-format pftrace csv`.
  - The executable is built with the option `enable-roctx` of mfem-mgis
    (`build-roctx.sh`): the timed sections, the time steps and the
    resolutions are ranges of roctx, shown in the trace above the calls of
    the HIP runtime, the kernels and the memory copies.
  - `trace-ranks.sh` is the launcher placed before the executable: only
    the processes `FM_TRACE_RANKS` run under the profiler (default
    `node:0`, the processes of the first node; `all`; ranks `0-3,128`).
    The other processes run directly.
  - `FM_REFERENCE=1` runs the same computation without profiler first, for
    the overhead of the trace.
  - Results: one `trace_results.pftrace` per traced process, to open in
    <https://ui.perfetto.dev>, the same events as `csv` tables, and
    `trace-zones.txt`. Beyond 4 nodes, the analysis is left to the login
    node (the command is written in `chain.log`).
  - `FM_TRACE_TOOL=rocprof-sys` runs the traced processes under
    `rocprof-sys-sample` instead (formerly Omnitrace): one
    `perfetto-trace-<rank>.proto` per process, with the ranges of roctx,
    the calls of HIP and of MPI, the kernels, the sampled call stacks of
    the host and counters of the node, without `csv` tables.
  - The trace of Perfetto is limited by its buffer (1 GB per process by
    default, 0.3 GB for 20 steps); `--perfetto-buffer-size 4194304` is
    refused by Perfetto and the program hangs at its end.
- **`trace-zones.py <directory of the traces>`: what the GPU does in each
  zone.**
  - Tree of the zones of the loop over the time steps (ranges of roctx):
    calls, time, share of the loop, time during which a kernel runs and
    occupancy of the GPU.
  - Innermost zones (each instant counted once), calls of the HIP runtime
    and memory copies in the loop; `--steps` adds the time steps.
- **`trace-allocations.py <directory of the traces>`: allocations on the
  GPU in each zone.** Calls to `hipMalloc` of at least 2 MB (arrays of the
  size of the unknowns or of the points), per zone, over the loop and over
  the second half of the time steps: a loop which keeps its arrays, or
  takes them in a pool, shows none in its second half.
- **`pftrace-tracks.py <file.pftrace>`: tracks of a Perfetto trace**, with
  their numbers of slices, decoded without library: what Perfetto will
  show, and whether events were discarded (the number of slices must be
  the number of rows of the tables).
- **`profile-host.slurm`: profile of the host of the final computation**
  (15/11/2026): which functions of the host run while the GPU waits, which
  the traces of the GPU do not show.
  - `FM_GRAINS` grains on all the GPUs of the allocation, final
    configuration, one step by default (the setup and the first assembly
    of the tangent operator), the executable being `FM_EXECUTABLE`.
  - The processes `FM_TRACE_RANKS` (default `0`) run under Linux perf
    (`trace-ranks.sh` with `FM_TRACE_TOOL=perf`: `perf record -F 199
    --call-graph dwarf,16384`). The call stacks are unwound without frame
    pointers, and no debug information is needed for the names of the
    functions: the executables of the measurements are profiled as they
    are.
  - Results: `trace/rank<rank>-<host>/perf.data` (250 MB for one minute).
    The command of `perf report` written in `chain.log`, to run on the
    login node (a few seconds; it did not end in 20 minutes on a node),
    gives `report-rank<rank>.txt`: share of the samples of each function,
    with and without its callees. Only the user code is sampled (event
    `cycles:Pu`): the allocations and the waits are not in the profile.
- **`roctx-tests.slurm`: tests of the ranges of roctx.**
  - Separate build by `build-roctx.sh` of MGIS, mfem-mgis and the
    polycrystal from the archive `FM_FIX/src.tgz`, on the local disk of the
    node (the quota of files of the project on `/lus/scratch` was nearly
    reached): variants without the option `enable-roctx` (`off`), with it
    (`on`) and with the roctx of roctracer (`legacy`). Only the executables
    and the library of MGIS are kept, in `FM_FIX/bin` and `FM_FIX/lib`.
  - `ProfilingTest` of MGIS, `ProfilingMarkersTest` of mfem-mgis on the
    host, on a GPU and under `rocprofv3` (the ranges must be in the trace),
    the polycrystal of the tests with two traced processes out of four
    (and one under `rocprof-sys-sample` with `FM_ROCPROF_SYS=1`), the
    tests of the polycrystal by `ctest`, other tests of mfem-mgis
    (`FM_OTHER_TESTS`).
- **`roctx-ab.slurm`: cost of the ranges without profiler.**
  - The final computation with the executables `FM_EXECUTABLES`
    (`<name>=<path>`), one after the other at each repetition, on the same
    nodes.
- **`perf-tests.slurm`: build and tests of the performance changes**
  (MFEM with Umpire, work vectors, global matrix assembled on the GPU, mean
  stresses computed by kernels, options of BoomerAMG; `JOURNAL.md`,
  2026-11-15).
  - Separate build by `build-perf.sh` from the archive `FM_FIX/src.tgz`
    (sources of MGIS, mfem-mgis and mm-opera-hpc, and the machine profile
    `machine.sh` when its options of MFEM are not deployed yet), on the
    local disk of the node: MFEM with Umpire (kept in
    `FM_FIX/install/mfem`), then the variants `umpire` and `deployed` (the
    MFEM of the deployed prototype) of mfem-mgis, with `enable-roctx`, and
    of the polycrystal (`FM_FIX/bin/uniaxial-polycrystal-<variant>`).
    Nothing is written on `/lus/work`.
  - `FM_FIX` is the directory of one version of the sources (checksum of
    the archive in `sources.sha256`, another archive is refused): the
    executables load the library of MGIS of `FM_FIX/lib`, and their run
    path does not name the build directory. `FM_FIX/install` may be a link
    to the installation of MFEM of another directory.
  - Tests of the changed components and of the non linear solvers on the
    host and on a MI300A (kernels forced, default, tests selecting a
    device), tests of the polycrystal, then the polycrystal of the tests on
    4 GPUs: each process seeing its GPU, all the processes seeing the 4
    GPUs (allocators of Umpire `DEVICE::<id>`), mean stresses by kernels.
- **`perf-ab.slurm`: the changes measured on an RVE.**
  - The final computation in the variants `ref` (executable of the
    measurements of the FM), `deployed`, `umpire`, `umpire-kernels` (the
    final configuration) and `amg-<name>` (options of BoomerAMG of the
    elastic operator `FM_BOOMERAMG_<name>`), one after the other at each
    repetition, on the GPUs, and on the cores with `FM_CPU_REPS`.
  - The computations are compared to the first one: axial stress within
    `FM_TOLERANCE` on the GPUs, output file bit for bit on the cores.
  - `sections.txt`: times of the main sections and levels of the multigrid
    of each computation.
- **`perf-load.slurm`: more unknowns per GPU than the series** (1600 grains
  on 2 nodes: 6.6 M per GPU), with the memory of the first node read every
  2 s.
- **`perf-mutations.slurm`: mutations of the code of these changes**, each
  one having to make its tests fail in at least one mode, in a build of
  the local disk. The job fails if a mutation is not detected. It also
  holds the mutations of the iterations of the linear solver written in
  the messages and of the pattern of the local matrix (`FM_MUTATIONS="i1
  i2 s1 s2 s3"`), and the ones of the setup shared between the threads of
  OpenMP (16/11/2026): pattern of the local matrix (`s4 s5 s6`), ranges of
  the threads (`t1 t2 t3`, tested by `HostThreadsTest`) and setup of the
  global matrix (`g5 g6`).
- **`host-threads.slurm`: setup of the host shared between threads**
  (option `enable-openmp` of mfem-mgis, 16/11/2026).
  - `build-perf.sh` builds the variant `umpire` with this option, and the
    variant `umpire-serial` without it.
  - The job first writes `placement.txt`: the cores that the launcher gives
    to each process of a GPU computation (4 processes per node, 24 cores
    each) and of a CPU computation (96 processes per node, one core each),
    and the settings of OpenMP they see.
  - Then the final computation of `FM_GRAINS` grains, one step by default,
    with the executable `FM_EXECUTABLE` for each number of threads per
    process of `FM_THREADS` (`OMP_NUM_THREADS`, given to the computation
    by the variable `gpu_environment` of
    `cases/measure-polycrystal-grains.sh`), and with `FM_REFERENCE`
    before them, `FM_REPS` times.
  - `threads.txt`: time of the first assembly of the tangent operator, of
    the mesh and of the computation, Newton iterations and axial stress of
    each computation.
- **`host-threads-identity.slurm`: same outputs whatever the number of
  threads.** A computation on the cores (100 grains without refinement, one
  step by default) with 4 processes per node of 24 cores each (variable
  `cpu_threads` of `cases/measure-polycrystal-grains.sh`), for each number
  of threads per process of `FM_THREADS` (default `1 7 24`): the `.res`
  file, the averages per grain and the lines of the iterations must be the
  ones of the first number, bit for bit.
- **`perf-walltime.py`: wall clock time of the functions of one thread**,
  from the samples of `profile-host.slurm`:
  `perf script -i perf.data -F tid,time,ip,sym,dso | perf-walltime.py
  --root <function>`.
  - `perf report` counts cycles of the code of the user: the waits (GPU,
    MPI) and the first accesses to the pages of a large array are not in
    its shares. Here each sample is given the time elapsed since the
    previous sample of its thread: the times of the call tree add up to
    the wall clock time of its root.
  - `--timeline` prints the successive steps of a call instead (the
    function called by the root, as long as it stays the same), which
    separates the phases of a setup.
  - The perf of the nodes (5.14) does not unwind the call stacks through
    the library of MFEM (built by the clang of ROCm) nor through the HIP
    runtime: these stacks are cut, and gathered under `(cut call stacks)`.
- **Chains of jobs.** `perf-tests.slurm`, `ls-identity.slurm`,
  `perf-ab.slurm`, `perf-load.slurm`, `perf-mutations.slurm` and
  `host-threads.slurm` fail when a check fails and write
  the file `passed` in their results directory otherwise: the next job is
  submitted with `--dependency=afterok:<job>` and
  `FM_REQUIRE=<results>/passed` (see `common.sh`), so that a large job does
  not start after a failure. `FM_CHAIN_LOG` names a second journal
  receiving the lines of `chain.log`.

## Scalability studies

`scaling.sh` submits chains of jobs in the final configuration of the FM,
from the login node (15/11/2026). Each job depends on the previous one, so
that the cluster goes on without the login session.

- `scaling.sh mesh <grains>...`: `geo.slurm` then `mesh.slurm` for each
  number of grains. The orientations must already be in
  `meshes/series/v<grains>`.
- `scaling.sh gpu <grains> <nodes>...`: the parts of one RVE for all its
  numbers of nodes (4 processes per node) in one job, then one job of
  `FM_GPU_REPS` computations per number of nodes.
- `scaling.sh cpu <grains>:<nodes>[:keep]...`: the points one after the
  other (96 processes per node): parts, computation, and the parts removed
  by the partitioning job of the next point unless the point ends with
  `:keep`. The 6144 parts of a computation on 64 nodes are as many files:
  only one set at a time weighs on the quota of the group.
- `scaling.sh nodes <nodes>`: the points `FM_GPU_POINTS`
  (`grains[:computations]`) then `FM_CPU_POINTS` (`grains[:keep]`) of one
  number of nodes in one job (`scaling-nodes.slurm`): one wait in the queue
  for all of them. The parts of each CPU point are written in the job just
  before its computation unless they exist, and removed after it unless the
  point ends with `:keep`; the nodes wait meanwhile (2 to 8 minutes).
- `scaling.sh parts <grains>:<processes>...`: partitions written and
  checked, for instance the kept ones after a study which removed them.
- `scaling.sh check <grains>:<processes>...`: partitions made before the
  file `checked` existed, read back.
- A computation starts only after its parts were written and checked
  (`afterok` and `FM_REQUIRE`). In a chain of the CPU, a point which fails
  does not stop the next ones (`afterany`, jobs with a dependency which can
  no longer be satisfied cancelled).
- `FM_AFTER` is the dependency of the first job of a chain, for instance
  the mesh job.

`scaling-fm.sh` is the study of the FM: `meshes`, then `computations`,
with the executable writing the iterations of the linear solver.

- Weak series: 600 N grains on N nodes, from 1 to 64 nodes (4.94 M
  unknowns per GPU, 206,000 per core).
- Strong series: 2400 grains from 4 to 64 nodes, 4800 from 8, 9600 from
  16 and 19,200 from 32.
- GPU three times and CPU once per point, one job per number of nodes (two
  on 64 nodes).

`scaling-table.py <series directories>` prints one line per RVE and number
of nodes (login node). `scaling-figures.py <table> <directory> [grains...]`
writes the tables (efficiencies, speedups, where the time goes, the same
without the launch of the processes) and three figures (matplotlib, on a
workstation): the largest point, the weak series and the strong series of
the RVEs given. `scaling-linear.py <series directories>` prints the
iterations of the linear solver of each computation: per solve, per solve
at the first tolerance of the Newton method (the same on both sides), slow
solves, first time steps. They tell the difficulty of the problem from the
parallel cost.

## The report

`report.py <results directory>` prints the report from the login node:
- the times of the loop over the steps (`NLEP::solve`);
- the speedup, also per Newton iteration;
- the energies and the powers per node;
- the counters of the steps and of the iterations.
