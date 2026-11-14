#!/usr/bin/env python3
"""What the GPU does in each zone of the computation.

usage: trace-zones.py <directory of the traces> [--steps] [--depth N]

Reads, for each traced process, the tables written by rocprofv3 with the
options --runtime-trace --marker-trace --output-format csv (directory of
fm-chain/trace-ranks.sh, one subdirectory per process):
- *marker_api_trace.csv: the ranges of roctx, that is the timed sections of
  mfem-mgis (`NLEP::solve`, `FPS::integrate`, ...), the time steps and the
  resolutions of the polycrystal (`time step from 3 to 4`, `resolution 0`);
- *kernel_trace.csv: the kernels run by the GPU;
- *hip_api_trace.csv, *memory_copy_trace.csv: the calls of the HIP runtime
  and the memory copies.

The loop over the time steps is the union of the ranges `time step ...`.
The GPU is busy when at least one kernel runs. Prints, as means over the
traced processes:
- the tree of the zones of the loop: calls, time, share of the loop, time
  during which the GPU is busy and occupancy of the GPU in the zone;
- the innermost zones: each instant of the loop is given to the deepest
  range containing it, so that the times add up to the loop;
- the calls of the HIP runtime and the memory copies in the loop;
- the steady regime, that is the second half of the time steps (the first
  ones hold the setup of the solvers and more resolutions);
- with --steps, the time steps of the first traced process.
--depth limits the depth of the tree (default 8)."""
import bisect
import csv
import glob
import os
import re
import sys
from collections import defaultdict


def read(path, name_column):
    """intervals (start, end, name) of a table of rocprofv3, sorted"""
    out = []
    try:
        with open(path, newline="") as f:
            rows = csv.reader(f)
            header = next(rows)
            s = header.index("Start_Timestamp")
            e = header.index("End_Timestamp")
            n = header.index(name_column) if name_column else None
            for r in rows:
                out.append((int(r[s]), int(r[e]),
                            r[n] if n is not None else ""))
    except (OSError, StopIteration, ValueError):
        pass
    out.sort()
    return out


class Busy:
    """union of intervals, and time covered by it in a given interval"""

    def __init__(self, intervals):
        self.starts, self.ends, self.before = [], [], [0]
        for s, e, _ in intervals:
            if self.ends and s <= self.ends[-1]:
                if e > self.ends[-1]:
                    self.before[-1] += e - self.ends[-1]
                    self.ends[-1] = e
            else:
                self.starts.append(s)
                self.ends.append(e)
                self.before.append(self.before[-1] + e - s)

    def until(self, t):
        i = bisect.bisect_right(self.starts, t)
        if i == 0:
            return 0
        return self.before[i] - max(self.ends[i - 1] - t, 0)

    def within(self, a, b):
        return self.until(b) - self.until(a)


def zone(name):
    """name of the zone of a range: the time steps and the resolutions are
    gathered, the first resolution of a time step apart from the next ones
    (it starts from the prediction)"""
    if name.startswith("time step from "):
        return "time step"
    m = re.fullmatch(r"resolution (\d+)", name)
    if m:
        return "resolution 0" if m.group(1) == "0" else "resolution 1+"
    return name


def analyse(prefix):
    """zones of the loop of one process"""
    ranges = [r for r in read(prefix + "marker_api_trace.csv", "Function")
              if r[1] > r[0]]
    kernels = read(prefix + "kernel_trace.csv", None)
    if not ranges:
        return None
    busy = Busy(kernels)
    # nesting of the ranges: sorted by start, the longest first (the time
    # steps and the resolutions before the sections, should two ranges have
    # the same times)
    ranges.sort(key=lambda r: (r[0], -r[1], not zone(r[2]).startswith(
        ("time step", "resolution")), not r[2].startswith("time step")))
    tree = defaultdict(lambda: [0, 0, 0])  # path -> calls, time, busy
    flat = defaultdict(lambda: [0, 0])  # zone -> self time, self busy
    steps = []  # start, end, busy time, resolutions
    resolutions = []  # starts of the resolutions
    # open ranges: end, start, path, in the loop, [time, busy] of the children
    stack = []

    def close():
        end, start, path, inside, children = stack.pop()
        if not inside:
            return
        t, b = end - start, busy.within(start, end)
        node = tree[path]
        node[0] += 1
        node[1] += t
        node[2] += b
        leaf = flat[path[-1]]
        leaf[0] += t - children[0]
        leaf[1] += b - children[1]
        if stack and stack[-1][3]:
            stack[-1][4][0] += t
            stack[-1][4][1] += b
        if len(path) == 1:
            steps.append([start, end, b, 0])
        elif len(path) == 2 and path[1].startswith("resolution"):
            resolutions.append(start)

    for s, e, n in ranges:
        while stack and stack[-1][0] <= s:
            close()
        z = zone(n)
        if stack and stack[-1][3]:
            path, inside = stack[-1][2] + (z,), True
        else:
            path, inside = (z,), z == "time step"
        stack.append((e, s, path, inside, [0, 0]))
    while stack:
        close()
    for st in steps:
        st[3] = sum(1 for r in resolutions if st[0] <= r < st[1])
    loop = Busy([(s, e, "") for s, e, _, _ in steps])
    total = sum(e - s for s, e, _, _ in steps)
    in_loop = [k for k in kernels if loop.within(k[0], k[1]) > 0]
    api = defaultdict(lambda: [0, 0])
    for s, e, n in read(prefix + "hip_api_trace.csv", "Function"):
        o = loop.within(s, e)
        if o > 0:
            api[n][0] += 1
            api[n][1] += o
    copies = defaultdict(lambda: [0, 0])
    for s, e, n in read(prefix + "memory_copy_trace.csv", "Direction"):
        o = loop.within(s, e)
        if o > 0:
            copies[n][0] += 1
            copies[n][1] += o
    durations = sorted(e - s for s, e, _ in in_loop)
    return {"tree": tree, "flat": flat, "steps": steps, "loop": total,
            "busy": sum(st[2] for st in steps), "kernels": len(in_loop),
            "median": durations[len(durations) // 2] if durations else 0,
            "ranges": len(ranges), "api": api, "copies": copies}


def main():
    args = sys.argv[1:]
    depth = 8
    if "--depth" in args:
        i = args.index("--depth")
        depth = int(args[i + 1])
        del args[i:i + 2]
    with_steps = "--steps" in args
    args = [a for a in args if a != "--steps"]
    if len(args) != 1:
        sys.exit(__doc__)
    directory = args[0]
    results = []
    for path in sorted(glob.glob(os.path.join(
            directory, "**", "*marker_api_trace.csv"), recursive=True)):
        prefix = path[: -len("marker_api_trace.csv")]
        r = analyse(prefix)
        name = os.path.basename(os.path.dirname(path))
        if r is None or not r["steps"]:
            print(f"{name}: no range of time step in {path}")
            continue
        results.append(r)
        print(f"{name}: loop {r['loop'] * 1e-9:.2f} s in {len(r['steps'])} "
              f"time steps and {sum(s[3] for s in r['steps'])} resolutions, "
              f"GPU busy {r['busy'] * 1e-9:.2f} s "
              f"({100 * r['busy'] / r['loop']:.1f} %), {r['ranges']} ranges, "
              f"{r['kernels']} kernels in the loop (median "
              f"{r['median'] * 1e-3:.0f} us)")
    if not results:
        sys.exit("no trace of markers with time steps")
    n = len(results)
    loop = sum(r["loop"] for r in results) / n
    busy = sum(r["busy"] for r in results) / n
    print(f"\nmean of {n} processes: loop {loop * 1e-9:.2f} s, GPU busy "
          f"{busy * 1e-9:.2f} s ({100 * busy / loop:.1f} %), idle "
          f"{(loop - busy) * 1e-9:.2f} s")
    # steady regime: second half of the time steps
    m = min(len(r["steps"]) for r in results)
    if m >= 4:
        t = sum(e - s for r in results for s, e, _, _ in r["steps"][m // 2:m])
        b = sum(st[2] for r in results for st in r["steps"][m // 2:m])
        print(f"steady regime (time steps {m // 2 + 1} to {m}): "
              f"{t * 1e-9 / n / (m - m // 2):.2f} s per time step, GPU busy "
              f"{100 * b / t:.1f} %")
    # tree of the zones
    tree = defaultdict(lambda: [0, 0, 0])
    for r in results:
        for path, v in r["tree"].items():
            for i in range(3):
                tree[path][i] += v[i] / n
    print("\nzones of the loop (mean per process): calls, time (s), share "
          "of the loop, GPU busy (s), occupancy of the GPU")
    children = defaultdict(list)
    for path in tree:
        children[path[:-1]].append(path)

    def show(path):
        c, t, b = tree[path]
        if len(path) > depth:
            return
        label = "  " * (len(path) - 1) + path[-1]
        print(f"  {label:<58} {c:8.0f} {t * 1e-9:9.2f} {100 * t / loop:6.1f} %"
              f" {b * 1e-9:9.2f} {100 * b / t if t else 0:6.1f} %")
        for child in sorted(children[path], key=lambda p: -tree[p][1]):
            show(child)

    for root in sorted(children[()], key=lambda p: -tree[p][1]):
        show(root)
    # innermost zones
    flat = defaultdict(lambda: [0, 0])
    for r in results:
        for z, v in r["flat"].items():
            flat[z][0] += v[0] / n
            flat[z][1] += v[1] / n
    print("\ninnermost zones (mean per process): time (s), share of the "
          "loop, GPU busy (s), occupancy, GPU idle (s), share of the idle "
          "time of the loop")
    idle = loop - busy
    for z, (t, b) in sorted(flat.items(), key=lambda x: -x[1][0])[:25]:
        if t < 5e6:
            continue
        print(f"  {z:<44} {t * 1e-9:9.2f} {100 * t / loop:6.1f} % "
              f"{b * 1e-9:9.2f} {100 * b / t:6.1f} % {(t - b) * 1e-9:9.2f} "
              f"{100 * (t - b) / idle if idle else 0:6.1f} %")
    print(f"  {'total':<44} {sum(v[0] for v in flat.values()) * 1e-9:9.2f}")
    # HIP runtime and memory copies
    api = defaultdict(lambda: [0, 0])
    copies = defaultdict(lambda: [0, 0])
    for r in results:
        for d, total in ((r["api"], api), (r["copies"], copies)):
            for k, v in d.items():
                total[k][0] += v[0] / n
                total[k][1] += v[1] / n
    print("\ncalls of the HIP runtime in the loop (mean per process): "
          "count, time (s), share of the loop")
    for k, (c, t) in sorted(api.items(), key=lambda x: -x[1][1])[:10]:
        print(f"  {k:<32} {c:10.0f} {t * 1e-9:9.2f} {100 * t / loop:6.1f} %")
    if copies:
        print("memory copies in the loop (mean per process): count, time (s)")
        for k, (c, t) in sorted(copies.items(), key=lambda x: -x[1][1]):
            print(f"  {k:<32} {c:10.0f} {t * 1e-9:9.2f}")
    if with_steps:
        print("\ntime steps of the first process: duration (s), "
              "resolutions, occupancy of the GPU")
        for i, (s, e, b, m) in enumerate(results[0]["steps"]):
            print(f"  {i + 1:4d} {(e - s) * 1e-9:8.2f} {m:3d} "
                  f"{100 * b / (e - s):6.1f} %")


main()
