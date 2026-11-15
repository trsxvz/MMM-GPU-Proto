#!/usr/bin/env python3
"""Iterations of the linear solver of the computations of a scalability
study, to tell the difficulty of the problem from the parallel cost.

usage: scaling-linear.py <series directory>...

Reads the logs of the computations of the given directories (series.slurm,
executable writing `<n> iterations of the linear solver` after each linear
solve) and prints one line per computation: grains, nodes, mode, label, then
- the linear solves and their iterations, over the computation;
- the iterations per solve;
- the iterations per solve at the first tolerance of the Newton method
  (iteration 0 of each resolution: the same tolerance, 0.01, for every size
  and on both sides, unlike the next ones, set by Eisenstat and Walker);
- the slow solves, which need more than SLOW_ITERATIONS iterations for a
  tolerance above SLOW_TOLERANCE (a loose tolerance is normally reached in
  a few iterations), and their iterations;
- the iterations of each of the first FIRST_STEPS time steps.
"""
import glob
import os
import re
import sys

FIRST_STEPS = 3
SLOW_ITERATIONS, SLOW_TOLERANCE = 30, 0.1
TOLERANCE = re.compile(r"relative tolerance of the linear solver (\S+)$")
SOLVE = re.compile(r"iteration +(\d+) : (\d+) iterations of the linear solver$")
NAME = re.compile(r"g(\d+)c-(.+)-r\d+-(gpu|cpu)-n(\d+)-")


def counts(path):
    """iterations of the solves, of the ones at the iteration 0 of a
    resolution, of the slow ones, and of each time step"""
    solves, first, slow, steps = [], [], [], []
    tolerance = 0
    with open(path, errors="replace") as log:
        for line in log:
            if line.startswith("Solving time step"):
                steps.append(0)
            m = TOLERANCE.search(line)
            if m:
                tolerance = float(m.group(1))
            m = SOLVE.search(line)
            if m and steps:
                n = int(m.group(2))
                solves.append(n)
                steps[-1] += n
                if m.group(1) == "0":
                    first.append(n)
                if n > SLOW_ITERATIONS and tolerance > SLOW_TOLERANCE:
                    slow.append(n)
    return solves, first, slow, steps


def main():
    print("grains nodes mode label solves iterations per_solve"
          " per_first_solve slow_solves slow_iterations",
          *(f"step{i + 1}" for i in range(FIRST_STEPS)))
    rows = []
    for directory in sys.argv[1:]:
        for path in glob.glob(os.path.join(directory, "*.txt")):
            m = NAME.match(os.path.basename(path))
            if not m:
                continue
            solves, first, slow, steps = counts(path)
            if not solves:
                continue
            rows.append((int(m.group(1)), int(m.group(4)), m.group(3),
                         m.group(2), len(solves), sum(solves),
                         "%.2f" % (sum(solves) / len(solves)),
                         "%.2f" % (sum(first) / len(first)),
                         len(slow), sum(slow),
                         *steps[:FIRST_STEPS]))
    for row in sorted(rows):
        print(*row)


if __name__ == "__main__":
    main()
