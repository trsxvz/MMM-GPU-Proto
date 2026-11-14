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
- **`ls-stats.py <log>...`: counters of a computation**: attempts of time
  steps by time increment, Newton corrections, iterations whose step was
  reduced by the line search, refused trials, failed searches.
  - For the first resolution of each attempt, which starts from the
    prediction: corrections, steps reduced at the first correction, and
    ratio `||r||/||r_0||` at the iteration 0 (below 1 when the prediction
    by linear extrapolation is better than the default one).
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
- **`pred-ab.slurm`: the extrapolation measured on an RVE.**
  - `FM_GRAINS` grains, one material for all the grains, 20 s simulated,
    options `FM_COMMON` (the line search), without and with
    `--prediction LinearExtrapolation`, one after the other at each
    repetition, on the same nodes.
  - `FM_GPU_REPS` repetitions on all the GPUs, `FM_CPU_REPS` on all the
    cores; `FM_VARIANTS` selects the variants (`ls`, `ls-pred`).

## The report

`report.py <results directory>` prints the report from the login node:
- the times of the loop over the steps (`NLEP::solve`);
- the speedup, also per Newton iteration;
- the energies and the powers per node;
- the counters of the steps and of the iterations.
