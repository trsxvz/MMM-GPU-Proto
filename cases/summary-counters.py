#!/usr/bin/env python3
"""Summary of the hardware counters measured by cases/measure-layout.sh
(counters=1): for each computation (file <name>-counters.csv exported by
Nsight Compute) and each kernel, number of launches measured, mean duration,
sectors per request of the global loads and stores, bytes read and written
in the DRAM per launch, throughput of the DRAM and of the SMs relative to
their peak.

usage: cases/summary-counters.py <results directory>
"""

import csv
import glob
import os
import re
import sys

DURATION = "gpu__time_duration.sum"
LD_REQUESTS = "l1tex__t_requests_pipe_lsu_mem_global_op_ld.sum"
LD_SECTORS = "l1tex__t_sectors_pipe_lsu_mem_global_op_ld.sum"
ST_REQUESTS = "l1tex__t_requests_pipe_lsu_mem_global_op_st.sum"
ST_SECTORS = "l1tex__t_sectors_pipe_lsu_mem_global_op_st.sum"
DRAM_BYTES = "dram__bytes.sum"
DRAM_PEAK = "dram__throughput.avg.pct_of_peak_sustained_elapsed"
SM_PEAK = "sm__throughput.avg.pct_of_peak_sustained_elapsed"

# factors converting the units of Nsight Compute to seconds and bytes
FACTORS = {
    "ns": 1e-9, "us": 1e-6, "usecond": 1e-6, "ms": 1e-3, "msecond": 1e-3,
    "nsecond": 1e-9, "s": 1, "second": 1,
    "byte": 1, "Kbyte": 1e3, "Mbyte": 1e6, "Gbyte": 1e9,
    "KB": 1e3, "MB": 1e6, "GB": 1e9,
}


def short_name(name):
    """name of a kernel: function of mfem-mgis launching it, or integration
    of the behaviour"""
    if ("generic_parallel::gpu" in name) or ("mfront::apply" in name):
        return "behaviour integration"
    m = re.search(r"(launch\w+Kernel)", name)
    if m:
        return m.group(1)
    m = re.search(r"mfem_mgis::(?:batch::)?(\w+)", name)
    return m.group(1) if m else name[:60]


def read(path):
    """launches of a file exported by `ncu --csv --page raw`: list of
    (kernel, values), the values being converted to seconds and bytes"""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        return []
    header, units = rows[0], rows[1]
    launches = []
    for row in rows[2:]:
        if len(row) != len(header) or not row[0].strip().isdigit():
            continue
        values = {}
        for h, u, v in zip(header, units, row):
            # the separators of the thousands depend on the locale (comma,
            # space, no-break spaces)
            try:
                x = float(re.sub(r"[,\s  ]", "", v))
            except ValueError:
                continue
            values[h] = x * FACTORS.get(u, 1)
        launches.append((short_name(row[header.index("Kernel Name")]), values))
    return launches


def ratio(a, b):
    return a / b if b > 0 else float("nan")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    files = sorted(glob.glob(os.path.join(sys.argv[1], "*-counters.csv")))
    columns = ("launches", "time (us)", "ld sect/req", "st sect/req",
               "DRAM MB", "DRAM %peak", "SM %peak")
    for path in files:
        name = os.path.basename(path)[: -len("-counters.csv")]
        kernels = {}
        for kernel, values in read(path):
            kernels.setdefault(kernel, []).append(values)
        print(f"== {name}")
        print(f"   {'kernel':40s}" + "".join(f"{c:>13s}" for c in columns))
        for kernel, launches in sorted(kernels.items()):
            def total(m):
                return sum(v.get(m, 0) for v in launches)
            n = len(launches)
            line = (
                n,
                total(DURATION) / n * 1e6,
                ratio(total(LD_SECTORS), total(LD_REQUESTS)),
                ratio(total(ST_SECTORS), total(ST_REQUESTS)),
                total(DRAM_BYTES) / n / 1e6,
                total(DRAM_PEAK) / n,
                total(SM_PEAK) / n,
            )
            print(f"   {kernel:40s}{line[0]:13d}"
                  + "".join(f"{x:13.2f}" for x in line[1:]))


if __name__ == "__main__":
    main()
