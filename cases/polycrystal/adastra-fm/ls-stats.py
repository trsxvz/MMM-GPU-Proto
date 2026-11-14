#!/usr/bin/env python3
"""Statistics of the line search in the logs of the polycrystal.
usage: ls-stats.py <log>...
For each log: attempts of time steps (by time increment: succeeded/failed),
Newton corrections (lines of iteration k >= 1), iterations whose step was
reduced, refused trials, failed searches, failed resolutions by cause.
Then, for the first resolution of each attempt (the one which starts from
the prediction): corrections, ratio ||r||/||r_0|| at the iteration 0 (1
without prediction; with the prediction by linear extrapolation, below 1
when the extrapolation is better than the default prediction), steps
reduced at the first correction, and extrapolated estimates discarded
because the behaviours could not be integrated at them (the resolution then
starts again from the default prediction).
"""
import re
import statistics
import sys
from collections import Counter, defaultdict

step = re.compile(r"Solving time step from\s+(\S+)\s+to\s+(\S+)")
it = re.compile(r"(?:Fixed point|Newton) iteration\s+(\d+)\s*: \|\|r\|\| = (\S+), \|\|r\|\|/\|\|r_0\|\| = (\S+)")
ls = re.compile(r"iteration\s+(\d+)\s*: line search, step length (\S+) after (\d+) reductions")
lsf = re.compile(r"iteration\s+(\d+)\s*: line search failed after (\d+) reductions")
moved = re.compile(r"integration failure, the estimate is moved")
fallback = "the resolution starts again from the default prediction"


def analyse(path):
    attempts = []  # dt, ok, corrections, reduced, refused, lsfailed, maxit
    cur = None
    corrections = reduced = refused = failed_searches = moves = fallbacks = 0
    thetas = Counter()
    resolutions = 0
    last_k = 0
    for line in open(path, errors="replace"):
        if fallback in line:
            fallbacks += 1
            continue
        m = step.search(line)
        if m:
            if cur is not None:
                attempts.append(cur)
            dt = float(m.group(2)) - float(m.group(1))
            cur = dict(dt=round(dt, 6), ok=True, corrections=0, reduced=0,
                       refused=0, lsfailed=False, t=float(m.group(1)),
                       resolutions=0, first_corrections=0, first_ratio=None,
                       first_reduced=False)
            continue
        if cur is None:
            continue
        m = it.search(line)
        if m:
            k = int(m.group(1))
            if k == 0:
                resolutions += 1
                cur["resolutions"] += 1
                if cur["resolutions"] == 1:
                    cur["first_ratio"] = float(m.group(3))
            else:
                corrections += 1
                cur["corrections"] += 1
                if cur["resolutions"] == 1:
                    cur["first_corrections"] += 1
            last_k = k
            continue
        m = ls.search(line)
        if m:
            reduced += 1
            n = int(m.group(3))
            refused += n
            cur["reduced"] += 1
            cur["refused"] += n
            thetas[m.group(2)] += 1
            if cur["resolutions"] == 1 and int(m.group(1)) == 0:
                cur["first_reduced"] = True
            continue
        m = lsf.search(line)
        if m:
            failed_searches += 1
            refused += int(m.group(2)) + 1
            cur["refused"] += int(m.group(2)) + 1
            cur["lsfailed"] = True
            continue
        if moved.search(line):
            moves += 1
            continue
        if "the resolution failed" in line or "maximum number of iterations for the" in line:
            cur["ok"] = False
            cur["last_k"] = last_k
    if cur is not None:
        attempts.append(cur)
    by_dt = defaultdict(lambda: [0, 0])
    for a in attempts:
        by_dt[a["dt"]][0 if a["ok"] else 1] += 1
    nfail = sum(1 for a in attempts if not a["ok"])
    print(f"{path}")
    print(f"  attempts {len(attempts)}: {len(attempts) - nfail} accepted, {nfail} failed"
          f" ({sum(1 for a in attempts if not a['ok'] and a['lsfailed'])} by a failed line search)")
    print("  by time increment (accepted/failed): " +
          ", ".join(f"{dt:g} s: {v[0]}/{v[1]}" for dt, v in sorted(by_dt.items(), reverse=True)))
    print(f"  resolutions {resolutions}, Newton corrections {corrections},"
          f" of which {sum(a['corrections'] for a in attempts if not a['ok'])} in failed attempts")
    pct = 100 * reduced / corrections if corrections else 0
    print(f"  iterations with a reduced step {reduced} ({pct:.1f} %), refused trials {refused},"
          f" failed searches {failed_searches}, estimates moved (no line search) {moves}")
    if thetas:
        print("  step lengths: " + ", ".join(f"{t}: {n}" for t, n in sorted(thetas.items(), key=lambda x: -float(x[0]))))
    first = [a for a in attempts if a["first_ratio"] is not None]
    if first:
        c = [a["first_corrections"] for a in first]
        r = [a["first_ratio"] for a in first]
        print(f"  first resolutions {len(first)}: corrections {sum(c)}"
              f" (min {min(c)}, median {statistics.median(c):g}, max {max(c)}),"
              f" step reduced at the first correction {sum(1 for a in first if a['first_reduced'])}")
        print(f"  ratio ||r||/||r_0|| at their iteration 0: min {min(r):.3g},"
              f" median {statistics.median(r):.3g}, max {max(r):.3g},"
              f" below 1: {sum(1 for v in r if v < 1)}, above 1: {sum(1 for v in r if v > 1)}")
        print(f"  extrapolated estimates discarded after an integration failure: {fallbacks}")
    for a in attempts:
        if not a["ok"]:
            print(f"    failed: t = {a['t']:g}, dt = {a['dt']:g}, {a['corrections']} corrections,"
                  f" {a['reduced']} reduced steps, {'line search' if a['lsfailed'] else 'iterations'}")


for p in sys.argv[1:]:
    analyse(p)
