#!/usr/bin/env python3
"""Runs a test case, then compares values of its output file to expected
values (reference solutions), with a relative tolerance.

usage: run-and-check.py --file <output file> --tolerance <relative>
           --check <row> <column> <expected> [--check ...] -- <command...>

- <row>: value of the first column of the row (for instance a time), or
  `first` (first row of values);
- <column>: index of the column (starting at 1) or name given by the header
  of the file;
- the output file is removed before running the command, so that a result
  of a previous run can't be checked.

The values are separated by spaces or by commas; the lines whose first
field is not a number are headers. Only the standard library is used.
"""

import argparse
import os
import subprocess
import sys


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_values(path):
    """returns the header (names of the columns) and the rows of values"""
    header = []
    rows = []
    with open(path) as f:
        for line in f:
            sep = "," if "," in line else None
            fields = [v.strip() for v in line.split(sep)]
            fields = [v for v in fields if v]
            if not fields:
                continue
            if is_number(fields[0]):
                rows.append(fields)
            elif not rows:
                header = fields
    return header, rows


def main():
    if "--" not in sys.argv:
        sys.exit("run-and-check.py: no command given (after --)")
    i = sys.argv.index("--")
    command = sys.argv[i + 1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--tolerance", type=float, required=True)
    parser.add_argument("--check", nargs=3, action="append", required=True,
                        metavar=("ROW", "COLUMN", "EXPECTED"))
    args = parser.parse_args(sys.argv[1:i])
    if os.path.exists(args.file):
        os.remove(args.file)
    status = subprocess.run(command).returncode
    if status != 0:
        sys.exit(f"run-and-check.py: the command failed (status {status})")
    if not os.path.exists(args.file):
        sys.exit(f"run-and-check.py: no output file '{args.file}'")
    header, rows = read_values(args.file)
    success = True
    for row, column, expected in args.check:
        if row == "first":
            selected = rows[0] if rows else None
        else:
            # the closest row, its time being equal to the given one up to
            # the rounding of the float type (the time 0.2 is written
            # 0.20000000298023 with float)
            closest = min(rows, key=lambda r: abs(float(r[0]) - float(row)),
                          default=None)
            selected = closest if (closest is not None and
                                   abs(float(closest[0]) - float(row)) <=
                                   1e-6 * max(1, abs(float(row)))) else None
        if selected is None:
            print(f"no row '{row}' in '{args.file}'")
            success = False
            continue
        c = int(column) - 1 if column.isdigit() else header.index(column)
        value = float(selected[c])
        reference = float(expected)
        error = abs(value - reference) / abs(reference)
        ok = error <= args.tolerance
        print(f"row {row}, column {column}: {value:.8g} (expected "
              f"{reference:.8g}, relative error {error:.3g}, tolerance "
              f"{args.tolerance:g}): {'OK' if ok else 'MISMATCH'}")
        success = success and ok
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
