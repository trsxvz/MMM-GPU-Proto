#!/usr/bin/env python3
"""Convergence in mesh of the polycrystal of mm-opera-hpc: compares the
computations of successive refinements of the same RVE.

usage: cases/polycrystal-convergence.py [--tolerance t] dir0 dir1 [dir2...]

Each directory is the one of a computation of
cases/measure-polycrystal-grains.sh (<results>/<computation>, refinements in
increasing order): its macroscopic stresses are read in <directory>.res and
the averages of the stresses per grain at each step in
<directory>/avgStressPolycrystal. For two successive refinements k and k + 1:

- macroscopic: largest relative difference of the axial stress Szz over
  the steps;
- per grain: root mean square over the grains of the difference of the
  average axial stress of the grain, relative to the macroscopic Szz, at
  the last step (and the largest one).

The retained refinement is the first k whose two differences with k + 1 are
below the tolerance (0.01 by default); printed on the last line as
`retained refinement: k` (or `retained refinement: none`), exit status 0 in
both cases.
"""
import math
import re
import sys


def read_res(path):
    """Axial stress (7th column) per step."""
    szz = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 7 and not line.startswith("#"):
                szz.append(float(fields[6]))
    return szz


def read_averages(path):
    """Average axial stress (component 2) per grain at the last step."""
    columns = {}
    last = None
    pattern = re.compile(r"# (\d+) column: (\d+) component .* material (\d+)")
    with open(path) as f:
        for line in f:
            m = pattern.match(line)
            if m:
                c, k, mi = (int(g) for g in m.groups())
                if k == 2:
                    columns[mi] = c - 1
            elif line.strip() and not line.startswith("#"):
                last = line.split()
    if last is None:
        return {}
    return {mi: float(last[c]) for mi, c in columns.items()}


def main():
    args = sys.argv[1:]
    tolerance = 0.01
    if len(args) >= 2 and args[0] == "--tolerance":
        tolerance = float(args[1])
        args = args[2:]
    if len(args) < 2:
        sys.exit(__doc__)
    levels = []
    for d in args:
        d = d.rstrip("/")
        try:
            szz = read_res(d + ".res")
        except OSError:
            szz = []
        try:
            grains = read_averages(d + "/avgStressPolycrystal")
        except OSError:
            grains = {}
        levels.append((d, szz, grains))
        print(f"{d}: {len(szz)} steps, {len(grains)} grains, "
              f"Szz at the end {szz[-1] if szz else '-'}")
    retained = None
    for k in range(len(levels) - 1):
        (_, s0, g0), (_, s1, g1) = levels[k], levels[k + 1]
        n = min(len(s0), len(s1))
        if n == 0 or not g0 or set(g0) != set(g1) or len(s0) != len(s1):
            print(f"refinements {k} and {k + 1}: incomplete computations")
            continue
        macro = max(abs(a - b) / abs(b) for a, b in zip(s0, s1))
        reference = abs(s1[-1])
        differences = [(g0[mi] - g1[mi]) / reference for mi in g1]
        rms = math.sqrt(sum(x * x for x in differences) / len(differences))
        largest = max(abs(x) for x in differences)
        ok = macro < tolerance and rms < tolerance
        print(f"refinements {k} and {k + 1}: Szz {macro:.2e}, "
              f"grains rms {rms:.2e} (largest {largest:.2e}) relative to Szz"
              f" -> {'converged' if ok else 'not converged'}"
              f" (tolerance {tolerance:g})")
        if ok and retained is None:
            retained = k
    print(f"retained refinement: {retained if retained is not None else 'none'}")


if __name__ == "__main__":
    main()
