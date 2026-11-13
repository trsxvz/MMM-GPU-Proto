#!/usr/bin/env python3
"""Cost of the failed resolutions of the polycrystal: usage: failures.py out.txt...

For each output, Newton iterations (lines `||r||/||r_0||`) of the resolutions
which converged and of those which failed (`error: the resolution failed.`),
backtrackings after an integration failure, and the share of the iterations
spent in failed resolutions."""
import re
import sys

for path in sys.argv[1:]:
    its = 0
    ok = failed = 0
    nok = nfailed = 0
    backtracks = 0
    lengths = []
    with open(path, errors="replace") as f:
        for line in f:
            if "||r||/||r_0||" in line and line.startswith("Fixed point iteration"):
                its += 1
            elif line.startswith("Fixed point iteration") and "integration failure" in line:
                backtracks += 1
            elif line.startswith("error: the resolution failed"):
                failed += its
                nfailed += 1
                lengths.append(its)
                its = 0
            elif line.startswith("Fixed Point iteration"):
                # end of a converged resolution (macroscopic fixed point line)
                ok += its
                nok += 1
                its = 0
    total = ok + failed
    share = failed / total if total else 0
    print(f"{path}: {total} Newton iterations, {failed} in {nfailed} failed "
          f"resolutions ({100 * share:.0f} %), {ok} in {nok} converged ones "
          f"({ok / max(nok, 1):.1f} per resolution), {backtracks} backtrackings "
          f"after integration failures")
