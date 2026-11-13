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

## The report

`report.py <results directory>` prints the report from the login node:
- the times of the loop over the steps (`NLEP::solve`);
- the speedup, also per Newton iteration;
- the energies and the powers per node;
- the counters of the steps and of the iterations.
