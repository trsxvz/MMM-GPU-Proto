#!/usr/bin/env python3
"""Fragmentation of the GPU work: usage: kernels.py <directory of rocprofv3>
[Newton iterations]

Reads the kernel traces (*kernel_trace.csv) of all the processes: per GPU
(agent), kernels launched, time spent in kernels (union of their intervals)
over the span from the first start to the last end (the rest is time
during which the GPU waits: launches, host work, MPI), size of the
kernels (work-items), and the kernels taking the most time."""
import csv
import glob
import os
import sys
from collections import defaultdict

d = sys.argv[1]
newton = float(sys.argv[2]) if len(sys.argv) > 2 else None
rows = []
for path in glob.glob(os.path.join(d, "**", "*kernel_trace.csv"), recursive=True):
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append((path, r))
if not rows:
    sys.exit("no kernel trace")
by_gpu = defaultdict(list)
names = defaultdict(lambda: [0, 0.0, 0])
for path, r in rows:
    s, e = int(r["Start_Timestamp"]), int(r["End_Timestamp"])
    key = (os.path.dirname(path), r.get("Agent_Id", "?"))
    grid = int(r.get("Grid_Size", r.get("Grid_Size_X", 0)) or 0)
    by_gpu[key].append((s, e, grid))
    n = r["Kernel_Name"][:90]
    names[n][0] += 1
    names[n][1] += (e - s) * 1e-9
    names[n][2] += grid
for key, ks in sorted(by_gpu.items()):
    ks.sort()
    busy, cur_s, cur_e = 0, None, None
    for s, e, _ in ks:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                busy += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    busy += cur_e - cur_s
    span = ks[-1][1] - ks[0][0]
    grids = sorted(g for _, _, g in ks)
    small = sum(1 for g in grids if g < 65536)
    durations = sorted(e - s for s, e, _ in ks)
    print(f"GPU {key[1]} ({os.path.basename(key[0])}): {len(ks)} kernels, "
          f"busy {busy * 1e-9:.1f} s of {span * 1e-9:.1f} s "
          f"({100 * busy / span:.0f} %), median duration "
          f"{durations[len(durations) // 2] * 1e-3:.0f} us, median size "
          f"{grids[len(grids) // 2]} work-items, {small} kernels below 65536 "
          f"work-items" + (f", {len(ks) / newton:.0f} kernels per Newton "
                           f"iteration" if newton else ""))
print("\nkernels taking the most time (all GPUs): count, time (s), mean "
      "work-items")
for n, (c, t, g) in sorted(names.items(), key=lambda x: -x[1][1])[:15]:
    print(f"{c:8d} {t:9.2f} {g / c:12.0f}  {n}")
