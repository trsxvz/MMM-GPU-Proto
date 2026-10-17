#!/usr/bin/env python3
"""summary of the trials of tune.py: per case and size, the total time, the
Newton iterations, the time of the linear resolutions (with the number of
calls) and of the tangent operator from the time table of the profiler,
the linear iterations printed by Ginkgo, the agreement of the quantity"""
import csv
import pathlib
import re
import sys
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cases"))
import benchmark_slides as b  # noqa: E402

R = ROOT / "results"
TUNE = R / "benchmark-slides-v2-tune" / "double"
SMOKE = R / "benchmark-slides-v2-smoke" / "double" / "runs.csv"
SECTION = re.compile(r"^ \|[ |]*(?:-->|>) (.+?)\s*\|\s*(\d+)\s*\|\s*"
                     r"([0-9.]+)\s*\|", re.MULTILINE)
GINKGO = re.compile(r"^\|\s*(\d+)\|", re.MULTILINE)


def sections(text):
    d = {}
    for name, calls, t in SECTION.findall(text):
        c, s = d.get(name, (0, 0.))
        d[name] = (c + int(calls), s + float(t))
    return d


def ginkgo_iterations(text):
    """the number of iterations of each resolution logged by Ginkgo"""
    its = []
    for block in text.split("Iteration log with residual norms:")[1:]:
        found = GINKGO.findall(block.split("Newton iteration")[0])
        if found:
            its.append(int(found[-1]))
    return its


references = {}
for r in csv.DictReader(open(SMOKE)):
    if r["configuration"] == "cpu-native" and r["status"] == "0":
        references[(r["case"], r["size"])] = r["quantity"]
rows = list(csv.DictReader(open(TUNE / "runs.csv")))
for r in sorted(rows, key=lambda r: (r["case"], r["size"])):
    log = TUNE / "logs" / (f"{r['case']}-{r['size']}-o{r['order']}-double-"
                          f"{r['configuration']}-np{r['processes']}-"
                          f"{r['repeat']}.log")
    text = log.read_text(errors="replace") if log.exists() else ""
    s = sections(text)
    solve = s.get("MFEM::Mult(r,c)", (0, 0.))
    tangent = s.get("NS::getJacobian", (0, 0.))
    its = ginkgo_iterations(text)
    q = [float(x) for x in r["quantity"].split()] if r["quantity"] else None
    ref = references.get((r["case"], r["size"]))
    ref = [float(x) for x in ref.split()] if ref else None
    agree = b.check_quantity(r["case"], "double", q, ref)
    print(f"{r['case']:11s} {r['size']:6s} "
          f"{r['configuration'].replace('tune-', ''):22s} "
          f"{r['status']:>7s} {float(r['seconds']):7.2f} s  newton "
          f"{r['iterations']:>3s}  solves {solve[0]:3d} {solve[1]:7.2f} s"
          f"  tangent {tangent[1]:6.2f} s"
          f"  qoi {'-' if agree is None else ('ok' if agree else 'DIFF')}"
          + (f"  ginkgo its {its}" if its else ""))
