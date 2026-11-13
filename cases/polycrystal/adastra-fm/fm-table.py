#!/usr/bin/env python3
"""Table of the GPUs against the cores of the same nodes, by number of grains.

usage: fm-table.py <results directory> <series directory>...

The computations of the given series directories (series.txt and the
<computation>.steps of measure-polycrystal-grains.sh) are grouped by number
of grains (label g<N>c...), on the same number of nodes. For each mode, the
median (GPU, non deterministic) or the mean (CPU) of:
- the duration of the Slurm step (time to solution: start, mesh, setup,
  steps), only for the computations which completed;
- the wall clock time from the end of the first step to the time T, the
  last simulated time reached by all the computations of the group (the
  CPU computation of 6,400 grains was stopped at 18 s of 20 by its time
  limit), read in the timestamps of the steps; the computations cancelled
  (status other than 0 and 124) are left out;
- the energy of the Slurm step (completed computations).
The speedups are the CPU values over the GPU ones.
"""
import glob
import os
import re
import statistics
import sys

FIELDS = ["label", "mode", "nodes", "np", "status", "unknowns", "total", "szz",
          "steps", "resolutions", "newton", "failures", "solver", "mesh", "loop",
          "energy", "elapsed", "cauchy", "pp", "wall"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def steps(path):
    """simulated time -> timestamp of the end of the step"""
    t = {}
    try:
        with open(path) as f:
            for line in f:
                w = line.split()
                if len(w) >= 2:
                    t[round(float(w[1]), 6)] = float(w[0])
    except OSError:
        pass
    return t


def runs(directories):
    out = []
    for d in directories:
        try:
            lines = open(os.path.join(d, "series.txt")).read().splitlines()
        except OSError:
            continue
        for line in lines:
            row = dict(zip(FIELDS, line.split()))
            m = re.match(r"g(\d+)c", row.get("label", ""))
            if not m:
                continue
            name = glob.glob(os.path.join(
                d, f"{row['label']}-r*-{row['mode']}-n{row['nodes']}-"
                   f"{row['solver']}.steps"))
            row["grains"] = int(m.group(1))
            row["dir"] = os.path.basename(d)
            row["times"] = steps(name[0]) if name else {}
            out.append(row)
    return out


def summary(values, gpu):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return statistics.median(values) if gpu else statistics.mean(values)


def fmt(x, f="%.0f"):
    return (f % x) if x is not None else "-"


def main():
    rows = runs(sys.argv[2:])
    groups = {}
    for r in rows:
        groups.setdefault((r["grains"], r["nodes"]), []).append(r)
    print("grains | nodes | unknowns | T (s) | GPU runs | CPU runs |"
          " elapsed GPU / CPU (s) | speedup elapsed | wall to T GPU / CPU (s) |"
          " speedup wall | energy GPU / CPU (MJ) | energy ratio |"
          " Newton GPU / CPU | failures GPU / CPU")
    for (grains, nodes), group in sorted(groups.items(),
                                         key=lambda x: (x[0][0], int(x[0][1]))):
        # completed computations, or stopped by their time limit (124): the
        # computations stopped otherwise (cancelled) are left out
        kept = lambda m: [r for r in group if r["mode"] == m and r["times"]
                          and r.get("status") in ("0", "124")]
        g, c = kept("gpu"), kept("cpu")
        if not g or not c:
            continue
        T = min(max(r["times"]) for r in g + c)
        def wall(r):
            t = r["times"]
            if T not in t or 1.0 not in t:
                return None
            return t[T] - t[1.0]
        done = lambda rs: [r for r in rs if r.get("status") == "0"]
        eg = summary([num(r["elapsed"]) for r in done(g)], True)
        ec = summary([num(r["elapsed"]) for r in done(c)], False)
        wg = summary([wall(r) for r in g], True)
        wc = summary([wall(r) for r in c], False)
        jg = summary([num(r["energy"]) for r in done(g)], True)
        jc = summary([num(r["energy"]) for r in done(c)], False)
        ng = summary([num(r["newton"]) for r in done(g)], True)
        nc = summary([num(r["newton"]) for r in done(c)], False)
        fg = summary([num(r["failures"]) for r in done(g)], True)
        fc = summary([num(r["failures"]) for r in done(c)], False)
        ratio = lambda a, b: (b / a) if a and b else None
        print(f"{grains} | {nodes} | {g[0]['unknowns']} | {T:g} | {len(g)} |"
              f" {len(c)} | {fmt(eg)} / {fmt(ec)} |"
              f" {fmt(ratio(eg, ec), '%.1f')} | {fmt(wg)} / {fmt(wc)} |"
              f" {fmt(ratio(wg, wc), '%.1f')} |"
              f" {fmt(jg and jg / 1e6, '%.1f')} / {fmt(jc and jc / 1e6, '%.1f')} |"
              f" {fmt(ratio(jg, jc), '%.1f')} | {fmt(ng)} / {fmt(nc)} |"
              f" {fmt(fg)} / {fmt(fc)}")
    print("\nruns: directory label mode nodes status elapsed energy newton"
          " failures last time")
    for r in rows:
        print(f"  {r['dir']} {r['label']} {r['mode']} {r['nodes']}"
              f" {r.get('status')} {r.get('elapsed')} {r.get('energy')}"
              f" {r.get('newton')} {r.get('failures')}"
              f" {max(r['times']) if r['times'] else '-'}")


if __name__ == "__main__":
    main()
