#!/usr/bin/env python3
# Allocations on the GPU in a trace of rocprofv3 (profile-trace.slurm): for
# each traced process, the calls to hipMalloc of at least a given size
# (default 2 MB: the arrays of the size of the unknowns or of the points),
# with their number, their time and their size per zone (the innermost range
# of roctx containing the call), over the whole loop over the time steps and
# over its second half (steady regime: the setups are done). A loop which
# allocates nothing shows no call in its second half: the arrays are kept,
# or taken in a pool (MFEM built with Umpire).
#
# usage: trace-allocations.py <directory of the traces> [minimum size in bytes]
import bisect
import csv
import os
import sys
from collections import defaultdict


def rows(path):
    with open(path, newline="") as f:
        yield from csv.DictReader(f)


def main():
    directory = sys.argv[1]
    minimum = int(sys.argv[2]) if len(sys.argv) > 2 else 2 << 20
    for rank in sorted(os.listdir(directory)):
        d = os.path.join(directory, rank)
        if not os.path.isfile(os.path.join(d, "trace_marker_api_trace.csv")):
            continue
        marks = sorted(
            (int(r["Start_Timestamp"]), int(r["End_Timestamp"]), r["Function"])
            for r in rows(os.path.join(d, "trace_marker_api_trace.csv")))
        steps = [m for m in marks if m[2].startswith("time step")]
        if not steps:
            print(f"{rank}: no time step in the trace of the markers")
            continue
        begin, end = steps[0][0], max(m[1] for m in steps)
        half = steps[len(steps) // 2][0]

        def zone(t):
            best = None
            for s, e, name in marks:
                if s > t:
                    break
                if t < e and (best is None or e - s < best[0]):
                    best = (e - s, name)
            return best[1] if best else "(none)"

        # sizes of the allocations, by start time
        allocations = sorted(
            (int(r["Start_Timestamp"]), int(r["Allocation_Size"]))
            for r in rows(os.path.join(d, "trace_memory_allocation_trace.csv"))
            if r["Operation"].endswith("ALLOCATE")
            and "DEALLOC" not in r["Operation"])
        starts = [a[0] for a in allocations]
        # [calls, time, bytes] per zone: whole loop, second half
        loop = defaultdict(lambda: [0, 0.0, 0])
        steady = defaultdict(lambda: [0, 0.0, 0])
        for r in rows(os.path.join(d, "trace_hip_api_trace.csv")):
            if r["Function"] != "hipMalloc":
                continue
            s, e = int(r["Start_Timestamp"]), int(r["End_Timestamp"])
            if not begin <= s < end:
                continue
            i = bisect.bisect_left(starts, s)
            size = 0
            while i < len(allocations) and allocations[i][0] < e:
                size += allocations[i][1]
                i += 1
            if size < minimum:
                continue
            z = zone(s)
            for table in (loop,) + ((steady,) if s >= half else ()):
                table[z][0] += 1
                table[z][1] += (e - s) * 1e-9
                table[z][2] += size
        print(f"{rank}: hipMalloc of at least {minimum} bytes, "
              f"{len(steps)} attempts of time steps")
        for title, table in (("loop over the time steps", loop),
                             ("second half of the time steps", steady)):
            total = [sum(v[k] for v in table.values()) for k in range(3)]
            print(f"  {title}: {total[0]} calls, {total[1]:.3f} s, "
                  f"{total[2] / 2**30:.2f} GiB")
            for z, v in sorted(table.items(), key=lambda kv: -kv[1][0]):
                print(f"    {z:42s} {v[0]:6d} calls {v[1]:9.3f} s "
                      f"{v[2] / 2**30:8.2f} GiB")


if __name__ == "__main__":
    main()
