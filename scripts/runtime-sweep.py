#!/usr/bin/env python3
"""Runtime comparisons of the test behaviours and models of MFront.

For each file of `mfront/tests` which compiles with the `generic-parallel`
interface, the file is built with the `generic` interface (reference) and
the `generic-parallel` interface with a configuration file, and the batch
integration is compared to the reference by the comparison program of the
tests of the `generic-parallel` interface (`compare.cxx`).

The material properties, the temperature and the initial values of the
internal state variables are read in the MTest files of the tests of the
`generic` interface (or of the `umat` interface of Cast3M). The time
increment is the one of the comparison program.

usage: scripts/runtime-sweep.py configuration [--jobs n] [--pattern regex]
                                [--tests directory] [--points n]

The configuration directory `scratch/sweep-env/<configuration>` must
contain the configuration file `config-<configuration>.json`, the mfront
wrapper `mfront-dev.sh` and the comparison program `compare` (see
scratch/dev-test.sh). The results are written in
`scratch/runtime-<configuration>/summary.tsv`.

Statuses:
- OK: the batch integration gives the results of the reference;
- MISMATCH: different statuses or results;
- INCONCLUSIVE: the integration fails at all the points of the first step,
  for all the tested types of tangent operator;
- CRASH: the comparison program crashed or was stopped after a timeout;
- NO_DATA: no MTest file gives the material properties;
- BUILD: the build of the libraries failed;
- UNSUPPORTED: modelling hypothesis or variable not handled by the
  comparison program.
"""

import argparse
import ctypes
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

PROTO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HYPOTHESES = ["Tridimensional", "PlaneStrain", "GeneralisedPlaneStrain",
              "Axisymmetrical", "PlaneStress",
              "AxisymmetricalGeneralisedPlaneStrain",
              "AxisymmetricalGeneralisedPlaneStress"]


def read_mtest_data(path):
    """return the data of an MTest file: real constants, material
    properties, external state variables, internal state variables and
    modelling hypothesis"""
    with open(path, errors="replace") as f:
        text = f.read()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    data = {"reals": {}, "mps": {}, "esvs": {}, "isvs": {},
            "hypothesis": None, "times": None, "loading": []}
    names = {"pi": math.pi}
    names.update({k: getattr(math, k) for k in
                  ["exp", "log", "sqrt", "sin", "cos", "tan", "pow"]})

    def value(v, t=0.0):
        v = v.strip()
        if v.startswith("{"):
            # evolution or array: first value
            items = [i for i in re.split(r"[,{}]", v) if i.strip()]
            first = items[0]
            if ":" in first:
                first = first.split(":", 1)[1]
            return value(first)
        v = v.strip("'\"")
        scope = dict(names)
        scope.update(data["reals"])
        scope["t"] = t
        return float(eval(v.replace("^", "**"), {"__builtins__": {}}, scope))

    for m in re.finditer(r"@(\w+)(<(\w+)>)?\s+([^;]*);", text):
        keyword, kind, rest = m.group(1), m.group(3), m.group(4).strip()
        try:
            if keyword == "Real":
                n, v = re.match(r"'([^']*)'\s+(.*)", rest, re.S).groups()
                data["reals"][n] = value(v)
            elif keyword == "ModellingHypothesis":
                data["hypothesis"] = rest.strip("'\" ")
            elif keyword == "Times":
                # first time and first time increment
                items = [i.strip() for i in rest.strip("{} ").split(",")]
                t0 = value(items[0])
                m1 = re.match(r"(.*)\bin\b(.*)", items[1])
                if m1 is None:
                    data["times"] = (t0, value(items[1]) - t0,
                                     value(items[-1].split("in")[0]))
                else:
                    t1, n = value(m1.group(1)), value(m1.group(2))
                    data["times"] = (t0, (t1 - t0) / n,
                                     value(items[-1].split("in")[0]))
            elif keyword in ("ImposedStrain", "ImposedDeformationGradient",
                             "ImposedGradient", "ImposedOpeningDisplacement"):
                n, v = re.match(r"'([^']*)'\s+(.*)", rest, re.S).groups()
                data["loading"].append((keyword, n, kind, v.strip()))
            elif keyword in ("MaterialProperty", "ExternalStateVariable",
                             "InternalStateVariable"):
                if kind not in (None, "constant", "function", "evolution"):
                    continue
                n, v = re.match(r"'([^']*)'\s+(.*)", rest, re.S).groups()
                key = {"MaterialProperty": "mps",
                       "ExternalStateVariable": "esvs",
                       "InternalStateVariable": "isvs"}[keyword]
                if keyword == "InternalStateVariable" and v.strip().startswith("{"):
                    items = [i for i in re.split(r"[,{}]", v) if i.strip()]
                    data[key][n] = [value(i) for i in items]
                else:
                    data[key][n] = value(v)
        except Exception:
            continue
    return data


def loading(data):
    """return the time increment of the first step of an MTest file and the
    amplitudes of the increments of the gradients tested: the maximal
    increment of the imposed gradients over this step, and ten times this
    value"""
    if data["times"] is None:
        return None, []
    t0, dt, t1 = data["times"]
    if (dt <= 0) or (t1 <= t0):
        return None, []
    increment = 0.0
    for keyword, name, kind, v in data["loading"]:
        try:
            if kind == "function":
                f = v.strip("'\"").replace("^", "**")
                scope = dict(data["reals"])
                scope.update({k: getattr(math, k) for k in ["exp", "log", "sqrt", "sin", "cos"]})
                g = lambda t: float(eval(f, {"__builtins__": {}}, dict(scope, t=t)))
                d = max(abs(g(t1) - g(t0)) / (t1 - t0) * dt,
                        abs(g(t0 + dt) - g(t0)))
            elif v.startswith("{"):
                pairs = [p.split(":") for p in v.strip("{} ").split(",") if ":" in p]
                ts = [float(eval(p[0], {"__builtins__": {}}, dict(data["reals"]))) for p in pairs]
                vs = [float(eval(p[1].strip().strip("'\""), {"__builtins__": {}},
                                 dict(data["reals"]))) for p in pairs]
                rate = max((abs(vs[i + 1] - vs[i]) / (ts[i + 1] - ts[i])
                            for i in range(len(ts) - 1) if ts[i + 1] > ts[i]),
                           default=0.0)
                d = rate * dt
            else:
                continue
            increment = max(increment, d)
        except Exception:
            continue
    if increment <= 0:
        return dt, []
    increment = min(max(increment, 1e-8), 1e-2)
    return dt, [f"{increment:.3g}", f"{min(10 * increment, 1e-2):.3g}"]


def find_mtest(tests, law):
    """return the MTest files testing the given behaviour, the tests of the
    `generic` interface first"""
    candidates = []
    patterns = [(os.path.join(tests, "behaviours", "generic"),
                 re.compile(r"@Behaviour<[^>]*>[^;]*'" + re.escape(law) + r"'")),
                (os.path.join(tests, "behaviours", "castem"),
                 re.compile(r"@Behaviour<[^>]*>[^;]*'umat" + re.escape(law.lower()) + r"'")),
                (os.path.join(tests, "behaviours", "aster"),
                 re.compile(r"@Behaviour<[^>]*>[^;]*'aster" + re.escape(law.lower()) + r"'")),
                (os.path.join(tests, "models"),
                 re.compile(r"@Model<[^>]*>[^;]*'" + re.escape(law) + r"'"))]
    for directory, pattern in patterns:
        if not os.path.isdir(directory):
            continue
        for f in sorted(os.listdir(directory)):
            if not f.endswith(".mtest"):
                continue
            p = os.path.join(directory, f)
            with open(p, errors="replace") as s:
                if pattern.search(s.read()):
                    candidates.append(p)
    return candidates


class Library:
    """symbols exported by a library generated by the `generic` interface"""

    def __init__(self, path):
        self.lib = ctypes.CDLL(path)

    def symbol(self, name, ctype):
        try:
            return ctype.in_dll(self.lib, name)
        except ValueError:
            return None

    def strings(self, prefixes, name):
        for p in prefixes:
            n = self.symbol(p + "n" + name, ctypes.c_ushort)
            if n is None:
                continue
            if n.value == 0:
                return []
            return [s.decode() for s in
                    (ctypes.c_char_p * n.value).in_dll(self.lib, p + name)]
        return None

    def ushort(self, prefixes, name):
        for p in prefixes:
            v = self.symbol(p + name, ctypes.c_ushort)
            if v is not None:
                return v.value
        return None


def elastic_names(orthotropic, hypothesis):
    if not orthotropic:
        return ["YoungModulus", "PoissonRatio"]
    names = ["YoungModulus1", "YoungModulus2", "YoungModulus3",
             "PoissonRatio12", "PoissonRatio23", "PoissonRatio13"]
    if hypothesis.startswith("AxisymmetricalGeneralised"):
        return names
    if hypothesis == "Tridimensional":
        return names + ["ShearModulus12", "ShearModulus23", "ShearModulus13"]
    return names + ["ShearModulus12"]


def get_behaviour_name(law_path):
    """return the name of the behaviour or of the model, prefixed by the
    name of the material if any"""
    with open(law_path, errors="replace") as f:
        text = re.sub(r"/\*.*?\*/", "", f.read(), flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    name = re.search(r"@(?:Behaviour|Model)\s+(\w+)\s*;", text)
    if name is None:
        return None
    material = re.search(r"@Material\s+(\w+)\s*;", text)
    if material is None:
        return name.group(1)
    return material.group(1) + "_" + name.group(1)


def build(args, law_path, directory):
    os.makedirs(directory, exist_ok=True)
    # data files read by some behaviours at runtime (orientations of the
    # grains of polycrystals), copied as by the tests of TFEL
    references = os.path.join(args.tests, "behaviours", "references")
    for f in os.listdir(references):
        if f.endswith(".txt"):
            shutil.copy(os.path.join(references, f), directory)
    tests = args.tests
    search = [f"--search-path={os.path.dirname(law_path)}"] + [
        f"--search-path={os.path.join(tests, d)}"
        for d in ["behaviours", "properties", "models"]]
    env_dir = os.path.join(PROTO, "scratch", "sweep-env", args.configuration)
    cmd = [os.path.join(env_dir, "mfront-dev.sh"), "--obuild",
           "--interface=generic", "--interface=generic-parallel",
           "--configuration-file=" + os.path.join(
               env_dir, f"config-{args.configuration}.json")] + search + [law_path]
    with open(os.path.join(directory, "build.log"), "w") as log:
        r = subprocess.run(cmd, cwd=directory, stdout=log,
                           stderr=subprocess.STDOUT, timeout=3600)
    return r.returncode == 0


def compare(args, directory, reference, batch, function, options):
    env = dict(os.environ)
    env["MFRONT_GPB_COMPARE_NUMBER_OF_POINTS"] = str(args.points)
    env["LD_LIBRARY_PATH"] = "src:" + env.get("LD_LIBRARY_PATH", "")
    compare_exe = os.path.join(PROTO, "scratch", "sweep-env",
                               args.configuration, "compare")
    cmd = [compare_exe, reference, batch, function] + options
    try:
        r = subprocess.run(cmd, cwd=directory, capture_output=True,
                           text=True, timeout=600, env=env)
    except subprocess.TimeoutExpired:
        return "CRASH", "timeout", ""
    out = r.stdout + r.stderr
    lines = [l for l in out.splitlines() if l.startswith("K0=")]
    if r.returncode not in (0, 1) or not lines:
        return "CRASH", f"exit code {r.returncode}", out[-500:]
    if all(l.endswith("OK") for l in lines):
        return "OK", "", out
    mismatch = False
    for l in lines:
        if l.endswith("OK"):
            continue
        if "same statuses=0" in l:
            mismatch = True
        for v in re.findall(r"(?:forces|isvs|K|energies|sound) ([0-9.e+-]+)", l):
            if float(v) > 1e-10:
                mismatch = True
    return ("MISMATCH" if mismatch else "INCONCLUSIVE"), "", out


def check(args, law_path):
    rel = os.path.relpath(law_path, args.tests)
    name = os.path.basename(law_path)[:-len(".mfront")]
    directory = os.path.join(args.output, rel.replace("/", "_")[:-len(".mfront")])
    try:
        if not build(args, law_path, directory):
            return rel, "BUILD", "see build.log"
    except subprocess.TimeoutExpired:
        return rel, "BUILD", "timeout"
    src = os.path.join(directory, "src")
    libs = os.listdir(src) if os.path.isdir(src) else []
    is_model = "libModel.so" in libs
    reference = os.path.join("src", "libModel.so" if is_model else "libBehaviour.so")
    batch = [l for l in libs if l.startswith("libGenericParallel") and l.endswith(".so")]
    if not os.path.exists(os.path.join(directory, reference)) or not batch:
        return rel, "BUILD", "libraries not found"
    batch = os.path.join("src", batch[0])
    # name of the behaviour (or of the model), which starts the names of
    # the generated functions
    behaviour = get_behaviour_name(law_path)
    lib = Library(os.path.join(directory, reference))
    if behaviour is None or lib.symbol(behaviour + "_mfront_interface",
                                       ctypes.c_char_p) is None:
        return rel, "BUILD", "no behaviour in the reference library"
    hypotheses = lib.strings([behaviour + "_"], "ModellingHypotheses") or []
    mtests = find_mtest(args.tests, behaviour) or find_mtest(args.tests, name)
    data = read_mtest_data(mtests[0]) if mtests else None
    hypothesis = None
    if data is not None and data["hypothesis"] in hypotheses:
        hypothesis = data["hypothesis"]
    for h in HYPOTHESES:
        if hypothesis is None and h in hypotheses:
            hypothesis = h
    if hypothesis is None:
        return rel, "UNSUPPORTED", "no supported modelling hypothesis"
    prefixes = [f"{behaviour}_{hypothesis}_", f"{behaviour}_"]
    mps = lib.strings(prefixes, "MaterialProperties") or []
    if lib.ushort(prefixes, "requiresStiffnessTensor"):
        orthotropic = lib.ushort(prefixes, "ElasticSymmetryType") == 1
        mps = elastic_names(orthotropic, hypothesis) + mps
    if lib.ushort(prefixes, "requiresThermalExpansionCoefficientTensor"):
        orthotropic = lib.ushort(prefixes, "ElasticSymmetryType") == 1
        mps = (["ThermalExpansion1", "ThermalExpansion2", "ThermalExpansion3"]
               if orthotropic else ["ThermalExpansion"]) + mps
    options = []
    if mps:
        if data is None:
            return rel, "NO_DATA", "no MTest file"
        values = []
        for n in mps:
            if n not in data["mps"]:
                return rel, "NO_DATA", f"material property '{n}' not found in {os.path.basename(mtests[0])}"
            values.append(repr(data["mps"][n]))
        options.append("--material-properties=" + ",".join(values))
    if data is not None and "Temperature" in data["esvs"]:
        options.append(f"--temperature={data['esvs']['Temperature']!r}")
    if lib.ushort(prefixes, "BehaviourType") == 2:
        options.append("--finite-strain=CAUCHY,DSIG_DF")
    if data is not None:
        esvs = lib.strings(prefixes, "ExternalStateVariables") or []
        if esvs:
            options.append("--external-state-variables=" + ",".join(
                repr(data["esvs"].get(n, data["esvs"].get("Temperature", 293.15)))
                for n in esvs))
        dt, amplitudes = loading(data)
        if dt is not None:
            options.append(f"--time-increment={dt!r}")
        if amplitudes:
            args_amplitudes = amplitudes
        else:
            args_amplitudes = args.amplitudes
    else:
        args_amplitudes = args.amplitudes
    function = f"{behaviour}_{hypothesis}"
    detail = os.path.basename(mtests[0]) if mtests else "no material property"
    status, out = "INCONCLUSIVE", ""
    for amplitude in args_amplitudes:
        for tangent in ["4", "1", "0"]:
            o = options + [f"--amplitude={amplitude}", f"--tangent-operator={tangent}"]
            status, why, out = compare(args, directory, reference, batch, function, o)
            if status != "INCONCLUSIVE":
                break
        if status != "INCONCLUSIVE":
            break
    with open(os.path.join(directory, "compare.log"), "w") as f:
        f.write(out)
    return rel, status, detail


def main():
    p = argparse.ArgumentParser()
    p.add_argument("configuration")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--pattern", default=".")
    p.add_argument("--tests", default=os.path.join(PROTO, "tfel-dev", "mfront", "tests"))
    p.add_argument("--points", type=int, default=1009)
    p.add_argument("--amplitudes", default="1e-4,1e-3")
    p.add_argument("--compile-sweep", default=None,
                   help="summary of a compilation sweep: only the files "
                   "compiled successfully are tested")
    args = p.parse_args()
    args.amplitudes = args.amplitudes.split(",")
    args.output = os.path.join(PROTO, "scratch", f"runtime-{args.configuration}")
    os.makedirs(args.output, exist_ok=True)
    files = []
    for d in ["behaviours", "models"]:
        for root, _, names in os.walk(os.path.join(args.tests, d)):
            files += [os.path.join(root, n) for n in names if n.endswith(".mfront")]
    files = sorted(f for f in files if re.search(args.pattern, f))
    if args.compile_sweep is not None:
        ok = set()
        with open(args.compile_sweep) as s:
            for l in s:
                c = l.rstrip("\n").split("\t")
                if len(c) >= 2 and c[1] == "OK":
                    ok.add(c[0])
        files = [f for f in files if os.path.relpath(f, args.tests) in ok]
    summary = os.path.join(args.output, "summary.tsv")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool, open(summary, "w") as s:
        for rel, status, detail in pool.map(lambda f: check(args, f), files):
            s.write(f"{rel}\t{status}\t{detail}\n")
            s.flush()
            print(f"{rel}\t{status}\t{detail}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
