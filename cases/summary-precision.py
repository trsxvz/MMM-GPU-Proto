#!/usr/bin/env python3
"""Summary of the measurements float against double of cases/measure-precision-beam.sh.

For each computation (precision and Newton criterion): minimum and median
time of the repetitions, Newton iterations, norm of the solution and its
relative gap to the one in double with the same criterion, then the ratio
of the times double / float.

usage: cases/summary-precision.py <results directory>
"""
import re
import statistics
import sys
from pathlib import Path

LINE = re.compile(r"^dofs (\d+) \| nls (\S+) \| ls (\S+) \| \|\|u\|\| (\S+) \| "
                  r"iterations per step:([\d ]+)\| total (\d+)")
TIME = re.compile(r"^time ([0-9.]+) s")
NAME = re.compile(r"^(double|float)-(.+)-(\d+)\.txt$")


def read(path):
    """reads an output of the beam: unknowns, norm, iterations, time"""
    dofs = norm = total = time = None
    for line in path.read_text(errors="replace").splitlines():
        m = LINE.match(line)
        if m:
            dofs, norm, total = int(m.group(1)), float(m.group(4)), int(m.group(6))
        m = TIME.match(line)
        if m:
            time = float(m.group(1))
    return dofs, norm, total, time


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    directory = Path(sys.argv[1])
    runs = {}
    for path in sorted(directory.glob("*.txt")):
        m = NAME.match(path.name)
        if not m:
            continue
        dofs, norm, total, time = read(path)
        if (time is None) or (norm is None):
            print(f"{path.name}: no time or no norm "
                  "(failed computation?)")
            continue
        runs.setdefault((m.group(2), m.group(1)), []).append(
            (time, norm, total, dofs))
    criteria = sorted({c for c, _ in runs})
    print(f"{'computation':38s} {'n':>2s} {'min (s)':>9s} {'median':>9s} "
          f"{'iter.':>6s} {'||u||':>14s} {'gap':>9s}")
    for criterion in criteria:
        reference = runs.get((criterion, "double"))
        for precision in ("double", "float"):
            values = runs.get((criterion, precision))
            if not values:
                continue
            times = [v[0] for v in values]
            norm, total = values[0][1], values[0][2]
            gap = ""
            if reference and precision == "float":
                gap = f"{abs(norm - reference[0][1]) / reference[0][1]:.2e}"
            print(f"{precision + '-' + criterion:38s} {len(values):2d} "
                  f"{min(times):9.3f} {statistics.median(times):9.3f} "
                  f"{total:6d} {norm:14.10g} {gap:>9s}")
        if (criterion, "double") in runs and (criterion, "float") in runs:
            d = min(v[0] for v in runs[(criterion, "double")])
            f = min(v[0] for v in runs[(criterion, "float")])
            print(f"{'  double / float (minimum times)':38s} {d / f:9.2f}")


if __name__ == "__main__":
    main()
