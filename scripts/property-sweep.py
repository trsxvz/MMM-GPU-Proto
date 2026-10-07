#!/usr/bin/env python3
"""Runtime comparisons of the test material properties of MFront.

For each file of `mfront/tests/properties`, the material property is built
with the `generic` interface (reference) and the `generic-parallel`
interface with a configuration file. The functions of the
`generic-parallel` interface are compared to one call of the reference per
point:

- the first prototype, all the arguments being given point by point;
- the first prototype, the first argument being uniform (null stride);
- the second prototype (contiguous values);
- the first prototype with values out of the bounds and the strict out of
  bounds policy: the batch must fail if the reference fails at one point.

The values of the arguments are drawn in their bounds (or physical bounds),
or in a range chosen from their names. With the CUDA and HIP backends, the
values are allocated in managed memory.

usage: scripts/property-sweep.py configuration [--jobs n] [--pattern regex]

The configuration directory `scratch/sweep-env/<configuration>` must
contain the configuration file `config-<configuration>.json` and the mfront
wrapper `mfront-dev.sh` (see scratch/dev-test.sh). The results are written
in `scratch/runtime-<configuration>/summary-properties.tsv`.

Statuses:
- OK: the batch gives the results of the reference;
- ROUNDOFF: same statuses, relative differences between 1e-10 and 1e-5
  (see scripts/runtime-sweep.py);
- MISMATCH: different statuses or larger differences;
- INCONCLUSIVE: the reference fails at all the points;
- CRASH: the comparison crashed or was stopped after a timeout;
- BUILD: the build of the libraries failed.
"""

import argparse
import ctypes
import os
import random
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

PROTO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

POLICY_WARNING = 1
POLICY_STRICT = 2


class OutputStatus(ctypes.Structure):
    """`mfront_gmp_OutputStatus`"""
    _fields_ = [("status", ctypes.c_int), ("c_error_number", ctypes.c_int),
                ("bounds_status", ctypes.c_int), ("msg", ctypes.c_char * 512)]


def default_range(name):
    """range of the values of an argument without bounds"""
    if name == "Temperature":
        return 293.15, 1273.15
    return 0.1, 1.0


def argument_range(lib, f, a):
    """range of the values of the argument `a` of the function `f`"""
    def bound(kind):
        for k in (kind + "Bound", kind + "PhysicalBound"):
            try:
                return float(ctypes.c_longdouble.in_dll(lib, f"{f}_{a}_{k}").value)
            except ValueError:
                pass
        return None
    lo, hi = bound("Lower"), bound("Upper")
    d0, d1 = default_range(a)
    if (lo is None) and (hi is None):
        return d0, d1
    if lo is None:
        return hi - (d1 - d0), hi
    if hi is None:
        return lo, max(d1, lo + (d1 - d0))
    return lo, hi


class Memory:
    """allocation of arrays of doubles, in managed memory if a GPU runtime
    is given"""

    def __init__(self, runtime):
        self.runtime = runtime
        if runtime is not None:
            self.cuda = ctypes.CDLL(runtime)

    def array(self, values):
        n = len(values)
        if self.runtime is None:
            return (ctypes.c_double * n)(*values)
        p = ctypes.c_void_p()
        if self.cuda.cudaMallocManaged(ctypes.byref(p), ctypes.c_size_t(8 * max(n, 1)),
                                       ctypes.c_uint(1)) != 0:
            raise RuntimeError("managed allocation failed")
        a = ctypes.cast(p, ctypes.POINTER(ctypes.c_double * n)).contents
        for i, v in enumerate(values):
            a[i] = v
        return a

    def synchronize(self):
        if self.runtime is not None:
            self.cuda.cudaDeviceSynchronize()


def compare_one(reference, batch, f, runtime, n=1009):
    """compare the functions `f` and `f2` of the batch library to the
    function `f` of the reference library. Returns a line of results."""
    ref = ctypes.CDLL(reference)
    lib = ctypes.CDLL(batch)
    nargs = ctypes.c_ushort.in_dll(ref, f + "_nargs").value
    names = ([s.decode() for s in (ctypes.c_char_p * nargs).in_dll(ref, f + "_args")]
             if nargs else [])
    rng = random.Random(42)
    values, wide = [], []
    for a in names:
        lo, hi = argument_range(ref, f, a)
        values.append([lo + (hi - lo) * rng.random() for _ in range(n)])
        # three times the range, centered on the range
        wide.append([lo - (hi - lo) + 3 * (hi - lo) * rng.random() for _ in range(n)])
    gf = ref[f]
    gf.restype = ctypes.c_double
    gf.argtypes = [ctypes.POINTER(OutputStatus), ctypes.POINTER(ctypes.c_double),
                   ctypes.c_size_t, ctypes.c_int]

    def reference_values(arguments, policy):
        r, statuses = [], []
        for i in range(n):
            s = OutputStatus()
            point = (ctypes.c_double * max(nargs, 1))(*[v[i] for v in arguments])
            r.append(gf(ctypes.byref(s), point, nargs, policy))
            statuses.append(s.status)
        return r, statuses

    memory = Memory(runtime)
    f1 = lib[f]
    f1.restype = None
    f2 = lib[f + "2"]
    f2.restype = None

    def batch_values(arguments, strides, contiguous, policy):
        arrays = [memory.array(v) for v in arguments]
        pointers = (ctypes.POINTER(ctypes.c_double) * max(nargs, 1))(
            *[ctypes.cast(a, ctypes.POINTER(ctypes.c_double)) for a in arrays])
        output = memory.array([0.0] * n)
        s = OutputStatus()
        if contiguous:
            f2(ctypes.byref(s), output, pointers, ctypes.c_size_t(nargs),
               ctypes.c_size_t(n), ctypes.c_int(policy))
        else:
            st = (ctypes.c_size_t * max(nargs, 1))(*strides)
            f1(ctypes.byref(s), output, ctypes.c_size_t(1), pointers, st,
               ctypes.c_size_t(nargs), ctypes.c_size_t(n), ctypes.c_int(policy))
        memory.synchronize()
        return list(output), s.status, s.msg.decode(errors="replace")

    results = []
    cases = [("by point", values, [1] * nargs, False, POLICY_WARNING),
             ("contiguous", values, [1] * nargs, True, POLICY_WARNING)]
    if nargs > 0:
        uniform = [[values[0][0]] * n] + values[1:]
        cases.append(("uniform", uniform, [0] + [1] * (nargs - 1), False,
                      POLICY_WARNING))
        cases.append(("out of bounds", wide, [1] * nargs, False, POLICY_STRICT))
    worst = "OK"
    for name, arguments, strides, contiguous, policy in cases:
        r, statuses = reference_values(arguments, policy)
        if name == "uniform":
            # only the first value of a uniform argument is given
            arguments = [[values[0][0]]] + values[1:]
        b, status, msg = batch_values(arguments, strides, contiguous, policy)
        valid = [i for i in range(n) if statuses[i] >= 0]
        if not valid:
            results.append(f"{name}: reference fails at all points, batch status {status}")
            if status >= 0:
                worst = "MISMATCH"
            elif name != "out of bounds":
                worst = "INCONCLUSIVE" if worst == "OK" else worst
            continue
        if (len(valid) != n) and (status < 0):
            # the batch failed, as the reference at some points: the values
            # computed by the batch are not meaningful
            results.append(f"{name}: batch status {status}, failures {n - len(valid)} ({msg})")
            continue
        same_status = (status < 0) == (len(valid) != n)
        m = max(abs(r[i]) for i in valid)
        d = max(abs(r[i] - b[i]) for i in valid)
        d = d / m if m != 0 else d
        line = f"{name}: batch status {status}, failures {n - len(valid)}, relative difference {d:.1e}"
        if msg:
            line += f" ({msg})"
        results.append(line)
        if (not same_status) or (d > 1e-5):
            worst = "MISMATCH"
        elif (d > 1e-10) and (worst in ("OK", "INCONCLUSIVE")):
            worst = "ROUNDOFF"
    return worst, "; ".join(results)


def check(args, law_path):
    rel = os.path.relpath(law_path, args.tests)
    directory = os.path.join(args.output, rel.replace("/", "_")[:-len(".mfront")])
    os.makedirs(directory, exist_ok=True)
    env_dir = os.path.join(PROTO, "scratch", "sweep-env", args.configuration)
    cmd = [os.path.join(env_dir, "mfront-dev.sh"), "--obuild",
           "--interface=generic", "--interface=generic-parallel",
           "--configuration-file=" + os.path.join(
               env_dir, f"config-{args.configuration}.json"),
           f"--search-path={os.path.dirname(law_path)}", law_path]
    with open(os.path.join(directory, "build.log"), "w") as log:
        try:
            r = subprocess.run(cmd, cwd=directory, stdout=log,
                               stderr=subprocess.STDOUT, timeout=3600)
        except subprocess.TimeoutExpired:
            return rel, "BUILD", "timeout"
    if r.returncode != 0:
        return rel, "BUILD", "see build.log"
    src = os.path.join(directory, "src")
    libs = sorted(l for l in os.listdir(src) if l.endswith(".so"))
    batch = [l for l in libs if l.startswith("libGenericParallel")]
    references = [l for l in libs if not l.startswith("libGenericParallel")]
    if (len(batch) != 1) or (len(references) != 1):
        return rel, "BUILD", "libraries not found"
    symbols = subprocess.run(["nm", "-D", "--defined-only",
                              os.path.join(src, references[0])],
                             capture_output=True, text=True).stdout.split()
    functions = sorted(s[:-len("_nargs")] for s in symbols if s.endswith("_nargs"))
    statuses, details = [], []
    for f in functions:
        cmd = [sys.executable, os.path.abspath(__file__), "--compare-one",
               os.path.join(src, references[0]), os.path.join(src, batch[0]), f]
        if args.runtime is not None:
            cmd += ["--runtime", args.runtime]
        try:
            with args.compare_slots:
                p = subprocess.run(cmd, cwd=directory, capture_output=True,
                                   text=True, timeout=600)
        except subprocess.TimeoutExpired:
            statuses.append("CRASH")
            details.append(f"{f}: timeout")
            continue
        out = p.stdout.strip().splitlines()
        if (p.returncode != 0) or (not out):
            statuses.append("CRASH")
            details.append(f"{f}: exit code {p.returncode}, {p.stderr.strip()[-300:]}")
            continue
        status, detail = out[-1].split("\t", 1)
        statuses.append(status)
        details.append(f"{f}: {detail}")
    with open(os.path.join(directory, "compare.log"), "w") as log:
        log.write("\n".join(details) + "\n")
    for s in ["CRASH", "MISMATCH", "ROUNDOFF", "INCONCLUSIVE", "OK"]:
        if s in statuses:
            return rel, s, f"{len(functions)} functions"
    return rel, "BUILD", "no function found"


def main():
    if (len(sys.argv) > 1) and (sys.argv[1] == "--compare-one"):
        p = argparse.ArgumentParser()
        p.add_argument("--compare-one", nargs=3)
        p.add_argument("--runtime", default=None)
        a = p.parse_args()
        status, detail = compare_one(*a.compare_one, a.runtime)
        print(f"{status}\t{detail}")
        return 0
    p = argparse.ArgumentParser()
    p.add_argument("configuration")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--compare-jobs", type=int, default=None)
    p.add_argument("--pattern", default=".")
    p.add_argument("--tests", default=os.path.join(PROTO, "tfel-dev", "mfront", "tests"))
    p.add_argument("--runtime", default=None,
                   help="CUDA runtime library used to allocate managed memory")
    args = p.parse_args()
    import threading
    args.compare_slots = threading.BoundedSemaphore(
        args.compare_jobs if args.compare_jobs is not None else args.jobs)
    args.output = os.path.join(PROTO, "scratch", f"runtime-{args.configuration}")
    os.makedirs(args.output, exist_ok=True)
    directory = os.path.join(args.tests, "properties")
    files = sorted(os.path.join(directory, f) for f in os.listdir(directory)
                   if f.endswith(".mfront") and re.search(args.pattern, f))
    summary = os.path.join(args.output, "summary-properties.tsv")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool, open(summary, "w") as s:
        for rel, status, detail in pool.map(lambda f: check(args, f), files):
            s.write(f"{rel}\t{status}\t{detail}\n")
            s.flush()
            print(f"{rel}\t{status}\t{detail}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
