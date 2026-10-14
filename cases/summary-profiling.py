#!/usr/bin/env python3
"""Summary per phase of the profiling tables of mfem-mgis.

usage: cases/summary-profiling.py file...

For each file (output of a case run with the profiling), sums the time of
each section (the maximum time over the processes in parallel), whatever
its position in the tree (the sections of the same name are not nested in
one another), and prints it per phase.
"""
import re
import sys

# phases: printed name, sections of mfem-mgis
PHASES = [
    ("total", ["root"]),
    ("integration", ["NS::processNewUnknownsEstimate", "FPS::integrate"]),
    ("residual", ["NS::computeResidual"]),
    ("tangent", ["NS::getJacobian"]),
    ("preconditioner", ["NS::setLinearSolverOperator"]),
    ("resolution", ["MFEM::Mult(r,c)"]),
    ("acceleration", ["FPS::accelerate"]),
    ("prediction", ["NLEPIB::computePrediction"]),
    ("end of step", ["NLEPIB::update"]),
]
# sequential table: calls, time, share; parallel: calls, min, average,
# max, share, imbalance
LINE = re.compile(r"^\s*\|[\s|]*(?:-->|>)\s*(\S.*?)\s*\|\s*(\d+)\s*\|"
                  r"((?:\s*[\d.]+%?\s*\|)+)\s*$")


def summarize(path):
    times = {}
    calls = {}
    with open(path, errors="replace") as f:
        for line in f:
            m = LINE.match(line)
            if m is None:
                continue
            name = m.group(1)
            values = [v.strip() for v in m.group(3).split("|") if v.strip()]
            # maximum time in parallel (4th column), time in sequential
            t = float(values[2] if len(values) > 3 else values[0])
            times[name] = times.get(name, 0.0) + t
            calls[name] = calls.get(name, 0) + int(m.group(2))
    return times, calls


def main():
    header = ["case"] + [p for p, _ in PHASES]
    print(" | ".join(header))
    print(" | ".join(["---"] * len(header)))
    for path in sys.argv[1:]:
        times, calls = summarize(path)
        if not times:
            continue
        row = [re.sub(r"^.*/|\.txt$", "", path)]
        for _, sections in PHASES:
            t = sum(times.get(s, 0.0) for s in sections)
            n = sum(calls.get(s, 0) for s in sections)
            row.append(f"{t:.2f} ({n})" if n else "-")
        print(" | ".join(row))


if __name__ == "__main__":
    main()
