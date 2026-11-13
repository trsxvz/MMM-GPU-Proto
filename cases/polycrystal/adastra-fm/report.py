#!/usr/bin/env python3
"""Report of the chain fm-chain: usage: report.py <results directory>.

Reads the lines of series.txt written by cases/measure-polycrystal-grains.sh:
label mode nodes processes status unknowns total Szz steps resolutions
newton failures solver mesh loop energy elapsed (the last five may be
missing in the older lines).
"""
import glob
import os
import sys

FIELDS = ["label", "mode", "nodes", "np", "status", "unknowns", "total", "szz",
          "steps", "resolutions", "newton", "failures", "solver", "mesh", "loop",
          "energy", "elapsed", "cauchy", "pp", "wall"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read_series(path):
    rows = []
    try:
        with open(path) as f:
            for line in f:
                fields = line.split()
                if len(fields) >= 8:
                    rows.append(dict(zip(FIELDS, fields)))
    except OSError:
        pass
    return rows


def fmt(x, f):
    return (f % x) if x is not None else "-"


def main():
    R = sys.argv[1]
    print(f"Report of {R}\n")
    print("== chain")
    try:
        print(open(os.path.join(R, "chain.log")).read())
    except OSError:
        print("no chain.log\n")
    directories = sorted(glob.glob(os.path.join(R, "*/series.txt")))
    print("== all the computations (label mode nodes processes status unknowns"
          " total Szz steps resolutions newton failures solver mesh loop energy"
          " elapsed cauchy pp wall)")
    for path in directories:
        print(f"-- {os.path.dirname(path)}")
        for row in read_series(path):
            print("   " + " ".join(row.get(k, "-") for k in FIELDS))
    print("\n== GPUs against cores of the same nodes (pairs gpu/cpu of the same"
          " label, nodes and solver label in a directory; times of the loop over"
          " the steps when known, total otherwise; elapsed: durations of the"
          " Slurm steps, start, mesh and setup included, the time to solution;"
          " wall: from the end of the first step to the end of the last one)")
    print("directory | label | nodes | solver | unknowns | GPU (s) | CPU (s) |"
          " speedup | speedup per Newton it. | speedup elapsed |"
          " speedup wall | energy CPU/GPU |"
          " power GPU/CPU (kW/node) | Szz GPU vs CPU |"
          " GPU steps/resol./Newton/failed | CPU steps/resol./Newton/failed")
    for path in directories:
        rows = read_series(path)
        pairs = {}
        for row in rows:
            k = (row["label"], row["nodes"], row.get("solver", "default"))
            pairs.setdefault(k, {})[row["mode"]] = row
        for (label, nodes, solver), modes in pairs.items():
            if "gpu" not in modes or "cpu" not in modes:
                continue
            g, c = modes["gpu"], modes["cpu"]
            tg = num(g.get("loop")) or num(g["total"])
            tc = num(c.get("loop")) or num(c["total"])
            speedup = tc / tg if tg and tc else None
            per_it = None
            ng, nc = num(g.get("newton")), num(c.get("newton"))
            if speedup and ng and nc:
                per_it = (tc / nc) / (tg / ng)
            # times to solution: both computations completed
            done = g.get("status") == "0" and c.get("status") == "0"
            de, ce = num(g.get("elapsed")), num(c.get("elapsed"))
            speedup_elapsed = ce / de if done and de and ce else None
            dw, cw = num(g.get("wall")), num(c.get("wall"))
            speedup_wall = cw / dw if done and dw and cw else None
            eg, ec = num(g.get("energy")), num(c.get("energy"))
            ratio = ec / eg if eg and ec else None
            n = num(nodes) or 1
            pg = (eg / num(g.get("elapsed")) / n / 1e3
                  if eg and num(g.get("elapsed")) else None)
            pc = (ec / num(c.get("elapsed")) / n / 1e3
                  if ec and num(c.get("elapsed")) else None)
            dszz = None
            if num(g["szz"]) and num(c["szz"]):
                dszz = abs(num(g["szz"]) - num(c["szz"])) / abs(num(c["szz"]))
            counts = lambda r: "/".join(r.get(k, "-") for k in
                                        ("steps", "resolutions", "newton",
                                         "failures"))
            print(f"{os.path.basename(os.path.dirname(path))} | {label} | {nodes}"
                  f" | {solver} | {g['unknowns']} | {fmt(tg, '%.1f')} | "
                  f"{fmt(tc, '%.1f')} | {fmt(speedup, '%.1f')} | "
                  f"{fmt(per_it, '%.1f')} | {fmt(speedup_elapsed, '%.1f')} | "
                  f"{fmt(speedup_wall, '%.1f')} | {fmt(ratio, '%.1f')} | "
                  f"{fmt(pg, '%.2f')}/{fmt(pc, '%.2f')} | {fmt(dszz, '%.1e')} | "
                  f"{counts(g)} | {counts(c)}")


if __name__ == "__main__":
    main()
