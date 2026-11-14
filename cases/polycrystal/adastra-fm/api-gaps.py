#!/usr/bin/env python3
"""Where the GPUs wait: usage: api-gaps.py <directory of rocprofv3>

Reads, for each process, the kernel trace (*kernel_trace.csv), the trace of
the HIP runtime (*hip_api_trace.csv) and, if any, the trace of the memory
copies (*memory_copy_trace.csv) of rocprofv3 (options --kernel-trace
--hip-runtime-trace --memory-copy-trace). The idle intervals of the GPU of
a process are the intervals without kernel between its first and its last
kernel. Each idle interval is attributed to the HIP calls of the process
which overlap it (the part of the call inside the interval); the rest of
the idle time is spent outside the HIP runtime (computations on the host,
MPI, ...). Prints, per process and for all of them, the idle time and its
attribution, the calls of the HIP runtime taking the most time, and the
memory copies."""
import bisect
import csv
import glob
import os
import sys
from collections import defaultdict


def intervals(path, name_key=None):
    out = []
    try:
        with open(path) as f:
            for r in csv.DictReader(f):
                s, e = int(r["Start_Timestamp"]), int(r["End_Timestamp"])
                out.append((s, e, r.get(name_key, "") if name_key else ""))
    except OSError:
        pass
    out.sort()
    return out


def merge(ks):
    busy = []
    for s, e, _ in ks:
        if busy and s <= busy[-1][1]:
            busy[-1][1] = max(busy[-1][1], e)
        else:
            busy.append([s, e])
    return busy


d = sys.argv[1]
total_idle = defaultdict(float)
total_api = defaultdict(lambda: [0, 0.0])
copies = defaultdict(lambda: [0, 0.0, 0])
span_sum = busy_sum = idle_sum = 0.0
for kpath in sorted(glob.glob(os.path.join(d, "**", "*kernel_trace.csv"),
                              recursive=True)):
    prefix = kpath[: -len("kernel_trace.csv")]
    ks = intervals(kpath)
    if not ks:
        continue
    busy = merge(ks)
    span = (busy[-1][1] - busy[0][0]) * 1e-9
    busy_t = sum(e - s for s, e in busy) * 1e-9
    gaps = [(busy[i][1], busy[i + 1][0]) for i in range(len(busy) - 1)]
    api = intervals(prefix + "hip_api_trace.csv", "Function")
    for s, e, n in api:
        total_api[n][0] += 1
        total_api[n][1] += (e - s) * 1e-9
    starts = [a[0] for a in api]
    idle = defaultdict(float)
    # calls overlapping each gap: the calls are sorted by start; a call
    # starting before the gap may overlap it (calls are short, a window of
    # calls before the gap is enough)
    for gs, ge in gaps:
        g = ge - gs
        covered = 0
        i = max(bisect.bisect_left(starts, gs) - 64, 0)
        while i < len(api) and api[i][0] < ge:
            s, e, n = api[i]
            o = min(e, ge) - max(s, gs)
            if o > 0:
                idle[n] += o * 1e-9
                covered += o
            i += 1
        idle["(outside the HIP runtime)"] += max(g - covered, 0) * 1e-9
    idle_t = sum(ge - gs for gs, ge in gaps) * 1e-9
    span_sum += span
    busy_sum += busy_t
    idle_sum += idle_t
    for n, t in idle.items():
        total_idle[n] += t
    top = ", ".join(f"{n} {t:.1f} s" for n, t in
                    sorted(idle.items(), key=lambda x: -x[1])[:4])
    print(f"{os.path.basename(os.path.dirname(kpath))}/"
          f"{os.path.basename(prefix)}: span {span:.1f} s, kernels "
          f"{busy_t:.1f} s ({100 * busy_t / span:.0f} %), idle {idle_t:.1f} s: "
          f"{top}")
    for s, e, n in intervals(prefix + "memory_copy_trace.csv", "Direction"):
        copies[n][0] += 1
        copies[n][1] += (e - s) * 1e-9
if span_sum:
    print(f"\nall processes: span {span_sum:.1f} s, kernels {busy_sum:.1f} s "
          f"({100 * busy_sum / span_sum:.0f} %), idle {idle_sum:.1f} s")
    print("idle time of the GPUs attributed to (s, share of the idle time):")
    for n, t in sorted(total_idle.items(), key=lambda x: -x[1])[:15]:
        print(f"  {t:9.2f} {100 * t / idle_sum:5.1f} %  {n}")
    print("calls of the HIP runtime taking the most time (count, s):")
    for n, (c, t) in sorted(total_api.items(), key=lambda x: -x[1][1])[:15]:
        print(f"  {c:9d} {t:9.2f}  {n}")
    if copies:
        print("memory copies (count, s):")
        for n, (c, t, _) in sorted(copies.items(), key=lambda x: -x[1][1]):
            print(f"  {c:9d} {t:9.2f}  {n}")
