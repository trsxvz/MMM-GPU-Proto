#!/usr/bin/env python3
"""Comparison of two files of numerical results.

usage: cases/compare-outputs.py reference file [--rtol r] [--atol a]
       [--groups g[:r]...] [--csv]

Compares one by one the numbers of the two files, arranged in lines of
columns (for instance: time, components of the gradient and of the
macroscopic stresses). The difference of each value must be lower than
atol + rtol * m, m being the largest absolute value of its group of
columns in the reference: a component close to zero (a null imposed
stress) is not compared to itself. By default, each column is a group;
with --groups, columns of the same nature share a scale, for instance
`--groups 0 1-3 4-6` for the time, the gradient and the stresses: a
transverse stress null in theory (numerical noise of 1e2 Pa under a
loading of 1e7 Pa) is then compared to the scale of the stresses. A group
may have its own relative tolerance (`1-3:1e-6`), instead of --rtol. A non
finite value (NaN, infinity) in one of the files is a difference. With
--csv, the values are separated by commas and the first line (the names
of the columns) is ignored. Prints the largest difference so scaled; non
zero exit code if the files differ.
"""
import argparse
import math
import sys


def read(path, csv=False):
    with open(path) as f:
        lines = f.readlines()
    if csv:
        lines = [line.replace(",", " ") for line in lines[1:]]
    return [[float(v) for v in line.split()] for line in lines if line.split()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("file")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--groups", nargs="*", default=[],
                        help="groups of columns sharing a scale (indices "
                        "from 0, intervals i-j), possibly followed by their "
                        "relative tolerance (i-j:r)")
    parser.add_argument("--csv", action="store_true",
                        help="values separated by commas, first line "
                        "ignored")
    args = parser.parse_args()
    a = read(args.reference, args.csv)
    b = read(args.file, args.csv)
    if len(a) != len(b) or not a or any(len(x) != len(y) for x, y in zip(a, b)):
        print("different shapes or empty file")
        return 1
    ncolumns = max(len(x) for x in a)
    scale = [0.0] * ncolumns
    rtol = [args.rtol] * ncolumns
    for x in a:
        for j, v in enumerate(x):
            scale[j] = max(scale[j], abs(v))
    for g in args.groups:
        interval, _, tolerance = g.partition(":")
        first, _, last = interval.partition("-")
        columns = range(int(first), int(last or first) + 1)
        if any(j >= ncolumns for j in columns):
            print(f"group {g}: missing column")
            return 1
        m = max(scale[j] for j in columns)
        for j in columns:
            scale[j] = m
            if tolerance:
                rtol[j] = float(tolerance)
    worst = 0.0
    ok = True
    n = 0
    for x, y in zip(a, b):
        for j, (u, v) in enumerate(zip(x, y)):
            n += 1
            d = abs(u - v)
            if (not (math.isfinite(u) and math.isfinite(v))) or \
                    d > args.atol + rtol[j] * scale[j]:
                ok = False
            if scale[j] != 0:
                worst = max(worst, d / scale[j])
    print(f"{n} values, largest difference relative to the maximum of its "
          f"group: {worst:.3e}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
