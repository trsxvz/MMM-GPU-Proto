#!/usr/bin/env python3
"""Triage of the differences found by scripts/runtime-sweep.py.

For each law whose status is ROUNDOFF or MISMATCH in the summary of a
configuration (a GPU configuration, or the host configuration, whose two
libraries may be contracted differently by the compiler), the comparison
giving the largest difference is run again on the host by
scripts/roundoff-check.sh, with the OpenMP libraries of the
openmp-host-gcc sweep: the library of the batch is compiled with and
without FMA contraction, the reference without. The difference is
explained by round-off errors if:

- the host gives no difference without FMA contraction;
- the host gives a difference of the same order with FMA contraction: at
  least a hundredth of the difference of the configuration, or of one if
  the statuses differ (the law then turns round-off errors into
  differences of the order of one).

usage: scripts/roundoff-triage.py configuration [--summary name]
                                  [--pattern regex]

The results are written in
`scratch/runtime-<configuration>/roundoff-triage.tsv`.
"""

import argparse
import os
import re
import subprocess
import sys

PROTO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def difference(line):
    """largest relative difference of a line of results of `compare`,
    infinite if the statuses differ"""
    if "same statuses=0" in line:
        return float("inf")
    values = re.findall(r"(?:forces|isvs|K|energies|sound) ([0-9.e+-]+|nan)", line)
    return max((float(v) if v != "nan" else float("inf") for v in values), default=0.0)


def worst_comparison(log):
    """arguments and largest difference of the comparisons of a log"""
    worst, arguments, current = -1.0, None, None
    with open(log) as f:
        for line in f:
            if line.startswith("# "):
                current = line.split()[4:]
            elif line.startswith("K0=") and (not line.rstrip().endswith("OK")):
                d = difference(line)
                if (d > worst) and (current is not None):
                    worst, arguments = d, current
    return worst, arguments


def main():
    p = argparse.ArgumentParser()
    p.add_argument("configuration")
    p.add_argument("--summary", default="summary.tsv")
    p.add_argument("--pattern", default=".")
    args = p.parse_args()
    results = os.path.join(PROTO, "scratch", f"runtime-{args.configuration}")
    host = os.path.join(PROTO, "scratch", "runtime-openmp-host-gcc")
    output = os.path.join(results, "roundoff-triage.tsv")
    with open(os.path.join(results, args.summary)) as s, open(output, "w") as out:
        for line in s:
            rel, status = line.split("\t")[:2]
            if (status not in ("ROUNDOFF", "MISMATCH")) or (not re.search(args.pattern, rel)):
                continue
            name = rel[:-len(".mfront")].replace("/", "_")
            gpu, arguments = worst_comparison(os.path.join(results, name, "compare.log"))
            if arguments is None:
                continue
            r = subprocess.run([os.path.join(PROTO, "scripts", "roundoff-check.sh"),
                                os.path.join(host, name)] + arguments,
                               capture_output=True, text=True)
            off, fast, block = 0.0, 0.0, None
            for l in r.stdout.splitlines():
                if l.startswith("batch compiled with"):
                    block = "off" if "=off" in l else "fast"
                elif l.startswith("K0=") and (block is not None):
                    d = difference(l)
                    if block == "off":
                        off = max(off, d)
                    else:
                        fast = max(fast, d)
            explained = (off == 0) and (fast > 0) and (fast >= min(gpu, 1.0) / 100)
            verdict = "round-off" if explained else "UNEXPLAINED"
            result = (f"{rel}\t{status}\t{args.configuration} {gpu:.1e}\t"
                      f"host without FMA {off:.1e}\thost with FMA {fast:.1e}\t{verdict}")
            out.write(result + "\n")
            out.flush()
            print(result, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
