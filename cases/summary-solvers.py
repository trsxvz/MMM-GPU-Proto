#!/usr/bin/env python3
"""Tables of the campaign of measurements of the linear solvers
(cases/measure-solvers.sh): one table per problem and per target (CPU or
GPU), sorted by the time of the first resolution, setup included.

usage: summary-solvers.py runs.txt > summary.md
"""

import sys


def read_runs(path):
    """Return the runs of a file of lines key=value."""
    runs = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            fields = dict(
                field.split("=", 1) for field in line.split() if "=" in field
            )
            if "label" in fields:
                runs.append(fields)
    return runs


def number(run, key):
    """Return a value of a run as a real, or None."""
    try:
        return float(run[key])
    except (KeyError, ValueError):
        return None


def main():
    runs = read_runs(sys.argv[1])
    groups = {}
    for run in runs:
        problem, _, target = run["label"].split("/")
        groups.setdefault((problem, target), []).append(run)
    print("# Linear solvers\n")
    print(
        "Times in seconds: setup (`SetOperator`) plus first resolution, "
        "following resolutions (mean), setup plus resolution for a second "
        "operator of the same sparsity. Relative residual of the first "
        "operator.\n"
    )
    for (problem, target), group in sorted(groups.items()):
        unknowns = group[0].get("unknowns", "?")
        processes = group[0].get("processes", "?")
        print(f"## {problem}, {target} ({unknowns} unknowns, "
              f"{processes} processes)\n")
        print("| configuration | setup + first | solve | second operator "
              "| iterations | residual | converged |")
        print("|---|---|---|---|---|---|---|")

        def first(run):
            s = number(run, "setup")
            f = number(run, "first")
            return (s or 0) + (f or 0)

        for run in sorted(group, key=first):
            label = run["label"].split("/")[1]
            second = (number(run, "resetup") or 0) + (number(run, "refirst") or 0)
            iterations = run.get("iterations", "")
            if iterations == "-1":
                iterations = "direct"
            print(f"| {label} | {first(run):.3g} | "
                  f"{number(run, 'solve') or 0:.3g} | {second:.3g} | "
                  f"{iterations} | {number(run, 'residual') or 0:.1e} | "
                  f"{'yes' if run.get('converged') == '1' else 'no'} |")
        print()


if __name__ == "__main__":
    main()
