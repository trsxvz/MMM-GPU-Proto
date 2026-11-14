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
- **`pftrace-tracks.py <file.pftrace>`: tracks of a Perfetto trace**, with
  their numbers of slices, decoded without library: what Perfetto will
  show, and whether events were discarded (the number of slices must be
  the number of rows of the tables).
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

## The report

`report.py <results directory>` prints the report from the login node:
- the times of the loop over the steps (`NLEP::solve`);
- the speedup, also per Newton iteration;
- the energies and the powers per node;
- the counters of the steps and of the iterations.
