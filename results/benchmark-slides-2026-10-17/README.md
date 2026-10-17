# Benchmark of the examples, 17/10/2026

Runs of `cases/benchmark_slides.py` (see `docs/benchmark-slides.md`) on
an AMD Ryzen 7 9800X3D (8 cores, 16 threads) and an NVIDIA GeForce RTX
5080, the machine being reserved for the screening of the CPU and for the
final runs.

- `calibration/`: each case at the measurement size with the solvers of
  the program, with the GPU (shared machine).
- `double/`: screening and final runs in double precision.
- `float/`: the same in single precision (float build), on the
  configurations which may win (the solvers of the program and the three
  best configurations of each side in double).

Each directory holds `runs.csv` (one line per run: case, size, order,
precision, configuration, side, processes, status, wall time of the whole
program, unknowns, command), `best.csv` (best times and speedups, from
`benchmark_slides.py best`) and `logs.tar.xz` (outputs of the runs). The
graphs are in `slides/`.
