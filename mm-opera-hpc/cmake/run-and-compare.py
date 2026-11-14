#!/usr/bin/env python3
"""Runs a test case twice, the second time with additional options, then
compares the values of the output files of the two runs, with a relative
tolerance.

usage: run-and-compare.py --directory <directory>
           --compare <file> <columns> <tolerance> [--compare ...]
           --option=<option> [--option=<option> ...] -- <command...>

- the command is run in <directory>/reference then, followed by the
  additional options, in <directory>/compared: an output file whose name is
  fixed by the case is overwritten neither by the other run nor by another
  test. The files given to the command must thus be absolute paths, and its
  output files relative ones;
- <file>: output file, in the directory of each run;
- <columns>: columns compared (starting at 1), as ranges separated by
  commas, for instance `1,2-4,5-7` (`2-`: from the second column to the last
  one). The values of a range are of the same nature (the components of a
  stress, for instance): the difference between two values is relative to
  the largest absolute value of the range in the row, so that a component
  which is null up to the rounding errors can be compared;
- <tolerance>: relative tolerance on the values of this file;
- the directories of the two runs are emptied before running the command,
  so that a result of a previous run can't be compared.

The values are separated by spaces or by commas; the lines whose first
field is not a number are headers. Only the standard library is used.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_values(path):
    """returns the rows of values"""
    rows = []
    with open(path) as f:
        for line in f:
            sep = "," if "," in line else None
            fields = [v.strip() for v in line.split(sep)]
            fields = [v for v in fields if v]
            if fields and is_number(fields[0]):
                rows.append([float(v) for v in fields])
    return rows


def read_ranges(columns):
    """returns the ranges of columns as pairs (first, last), starting at 1,
    the last one being None for the last column of the file"""
    ranges = []
    for r in columns.split(","):
        first, sep, last = r.partition("-")
        if not sep:
            last = first
        ranges.append((int(first), int(last) if last else None))
    return ranges


def compare(name, reference, compared, columns, tolerance):
    """compares two files, returns true if their values are the same up to
    the tolerance"""
    a = read_values(reference)
    b = read_values(compared)
    if not a:
        print(f"no row in '{reference}': MISMATCH")
        return False
    if [len(r) for r in a] != [len(r) for r in b]:
        print(f"'{name}': {len(a)} and {len(b)} rows, or rows of different "
              "sizes: MISMATCH")
        return False
    success = True
    for first, last in read_ranges(columns):
        error = 0.
        for ra, rb in zip(a, b):
            va = ra[first - 1:last]
            vb = rb[first - 1:last]
            # a missing column or a value which is not finite is a mismatch
            if len(ra) < (last or first) or \
                    not all(math.isfinite(v) for v in va + vb):
                error = math.inf
                break
            scale = max(abs(v) for v in va + vb)
            if scale != 0:
                difference = max(abs(x - y) for x, y in zip(va, vb))
                error = max(error, difference / scale)
        ok = error <= tolerance
        print(f"'{name}', columns {first}-{last if last else ''}, "
              f"{len(a)} rows: maximum relative difference {error:.3g} "
              f"(tolerance {tolerance:g}): {'OK' if ok else 'MISMATCH'}")
        success = success and ok
    return success


def main():
    if "--" not in sys.argv:
        sys.exit("run-and-compare.py: no command given (after --)")
    i = sys.argv.index("--")
    command = sys.argv[i + 1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", required=True)
    parser.add_argument("--compare", nargs=3, action="append", required=True,
                        metavar=("FILE", "COLUMNS", "TOLERANCE"))
    parser.add_argument("--option", action="append", required=True)
    args = parser.parse_args(sys.argv[1:i])
    directories = {}
    for run, options in (("reference", []), ("compared", args.option)):
        d = os.path.join(args.directory, run)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        directories[run] = d
        # the outputs of the command precede the ones of this script
        sys.stdout.flush()
        status = subprocess.run(command + options, cwd=d).returncode
        if status != 0:
            sys.exit(f"run-and-compare.py: the {run} command failed "
                     f"(status {status})")
    success = True
    for name, columns, tolerance in args.compare:
        files = [os.path.join(directories[r], name)
                 for r in ("reference", "compared")]
        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            sys.exit(f"run-and-compare.py: no output file '{missing[0]}'")
        success = compare(name, files[0], files[1], columns,
                          float(tolerance)) and success
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # reported as the other failures (see the tests expecting a pattern
        # in the output, whose exit status is ignored by ctest)
        sys.exit(f"run-and-compare.py: {e!r}")
