#!/usr/bin/env python3
"""Points of a scalability study: one line per RVE and number of nodes.

usage: scaling-table.py <series directory>...

Reads the series.txt of the given directories (series.slurm) and prints, for
each number of grains (label g<N>c...) and of nodes, the unknowns and, for
the GPUs then for the cores, the mean over the completed computations of:
the duration of the Slurm step (time to solution, s), the wall clock time
from the end of the first step to the end of the last one (s), the energy
of the step (J), the Newton iterations and the failed resolutions, after
the number of computations. `-` without a completed computation. The table
is read by scaling-figures.py.
"""
import os
import re
import statistics
import sys

FIELDS = ["label", "mode", "nodes", "np", "status", "unknowns", "total", "szz",
          "steps", "resolutions", "newton", "failures", "solver", "mesh", "loop",
          "energy", "elapsed", "cauchy", "pp", "wall"]
MEANS = ["elapsed", "wall", "energy", "newton", "failures"]


def completed(directories):
    """(grains, nodes) -> mode -> completed computations"""
    points = {}
    for d in directories:
        try:
            lines = open(os.path.join(d, "series.txt")).read().splitlines()
        except OSError:
            continue
        for line in lines:
            row = dict(zip(FIELDS, line.split()))
            m = re.match(r"g(\d+)c", row.get("label", ""))
            if m and row.get("status") == "0":
                key = (int(m.group(1)), int(row["nodes"]))
                points.setdefault(key, {}).setdefault(row["mode"], []).append(row)
    return points


def main():
    points = completed(sys.argv[1:])
    print("grains nodes unknowns",
          *(f"{mode}_{name}" for mode in ("gpu", "cpu")
            for name in ["runs"] + MEANS))
    for (grains, nodes), modes in sorted(points.items()):
        unknowns = next(iter(modes.values()))[0]["unknowns"]
        values = []
        for mode in ("gpu", "cpu"):
            rows = modes.get(mode, [])
            values.append(len(rows))
            for name in MEANS:
                values.append("%.6g" % statistics.mean(
                    float(r[name]) for r in rows) if rows else "-")
        print(grains, nodes, unknowns, *values)


if __name__ == "__main__":
    main()
