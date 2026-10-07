#!/usr/bin/env python3
"""Runtime comparisons of the test behaviours and models of MFront.

For each file of `mfront/tests` which compiles with the `generic-parallel`
interface, the file is built with the `generic` interface (reference) and
the `generic-parallel` interface with a configuration file, and the batch
integration is compared to the reference by the comparison program of the
tests of the `generic-parallel` interface (`compare.cxx`).

The material properties, the temperature, the other external state
variables, the initial values of the internal state variables, the values
of the real parameters, the time increment and the amplitude of the
loading are read in the MTest files of the tests of the `generic`
interface (or of the `umat` interface of Cast3M or of the `aster`
interface). The material properties which are not given
by an MTest file take synthetic values chosen from their names (the summary
tells it): they are not physical, but allow to run and compare the code of
the law. All the amplitudes of the loading are compared, so that the
inelastic part of the behaviour is tested, and the worst conclusive result
is kept. If the integration fails at all the points for all of them,
smaller time increments and amplitudes are tried.

usage: scripts/runtime-sweep.py configuration [--jobs n] [--compare-jobs n]
                                [--pattern regex] [--tests directory]
                                [--points n] [--summary name]

The comparisons may be run one at a time (`--compare-jobs 1`) while the
libraries are built in parallel: the behaviours whose threads need a large
stack may use most of the memory of a GPU, so that concurrent comparisons
fail with "out of memory" errors.

The configuration directory `scratch/sweep-env/<configuration>` must
contain the configuration file `config-<configuration>.json`, the mfront
wrapper `mfront-dev.sh` and the comparison program `compare` (see
scratch/dev-test.sh). The results are written in
`scratch/runtime-<configuration>/summary.tsv` (or in the file given by the
`--summary` option, for instance to check again some laws). The command of
the comparisons of a law and their outputs are written in its `compare.log`
file.

Statuses:
- OK: the batch integration gives the results of the reference;
- ROUNDOFF: same statuses, relative differences between 1e-10 and 1e-5.
  Those differences are expected when the libraries are compiled by
  different compilers and the law amplifies round-off errors (adaptive
  sub-stepping, numerical jacobian). This shall be checked law by law with
  scripts/roundoff-check.sh, which compiles the library of the batch on the
  host with FMA contraction (as nvcc in device code) and the reference
  without it;
- MISMATCH: different statuses or larger differences;
- INCONCLUSIVE: the integration fails at all the points of the first step,
  for all the tested types of tangent operator and all the tested loadings;
- TIMEOUT: a comparison was stopped after 600 seconds (for instance, an
  adaptive integration with a large amplitude), the other comparisons
  being OK or inconclusive;
- CRASH: the comparison program crashed;
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
import threading
from concurrent.futures import ThreadPoolExecutor

PROTO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HYPOTHESES = ["Tridimensional", "PlaneStrain", "GeneralisedPlaneStrain",
              "Axisymmetrical", "PlaneStress",
              "AxisymmetricalGeneralisedPlaneStrain",
              "AxisymmetricalGeneralisedPlaneStress"]


def read_mtest_data(path):
    """return the data of an MTest file: real constants, material
    properties, external state variables, internal state variables,
    parameters, modelling hypothesis, times and imposed gradients. The
    values are the ones at the first time: the definitions are evaluated
    until all the names they use are known, whatever their order in the
    file (as in MTest, the expressions may use the constants, the constant
    material properties, the external state variables and the evolutions)."""
    with open(path, errors="replace") as f:
        text = f.read()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    data = {"reals": {}, "mps": {}, "esvs": {}, "isvs": {}, "parameters": {},
            "evolutions": {}, "hypothesis": None, "times": None, "loading": []}
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
        for k in ["evolutions", "esvs", "mps", "reals"]:
            scope.update(data[k])
        scope["t"] = t
        return float(eval(v.replace("^", "**"), {"__builtins__": {}}, scope))

    statements = [(m.group(1), m.group(3), m.group(4).strip()) for m in
                  re.finditer(r"@(\w+)(<(\w+)>)?\s+([^;]*);", text)]
    # definitions of named values: group, name, expression, list of values
    definitions = []
    groups = {"Real": "reals", "Parameter": "parameters",
              "MaterialProperty": "mps", "ExternalStateVariable": "esvs",
              "InternalStateVariable": "isvs", "Evolution": "evolutions"}
    for keyword, kind, rest in statements:
        if keyword not in groups:
            continue
        if kind not in (None, "constant", "function", "evolution"):
            continue
        m = re.match(r"'([^']*)'\s+(.*)", rest, re.S)
        if m is None:
            continue
        n, v = m.groups()
        is_list = (keyword == "InternalStateVariable") and v.strip().startswith("{")
        definitions.append((groups[keyword], n, v, is_list))
    while definitions:
        remaining = []
        for group, n, v, is_list in definitions:
            try:
                if is_list:
                    items = [i for i in re.split(r"[,{}]", v) if i.strip()]
                    data[group][n] = [value(i) for i in items]
                else:
                    data[group][n] = value(v)
            except Exception:
                remaining.append((group, n, v, is_list))
        if len(remaining) == len(definitions):
            break
        definitions = remaining
    for keyword, kind, rest in statements:
        try:
            if keyword == "ModellingHypothesis":
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
    return dt, [f"{increment:.3g}", f"{min(10 * increment, 1e-2):.3g}",
                f"{min(100 * increment, 1e-2):.3g}"]


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

    def ints(self, prefixes, name, n):
        for p in prefixes:
            v = self.symbol(p + name, ctypes.c_int * n)
            if v is not None:
                return list(v)
        return None

    def ushort(self, prefixes, name):
        for p in prefixes:
            v = self.symbol(p + name, ctypes.c_ushort)
            if v is not None:
                return v.value
        return None


def variable_size(t, hypothesis):
    """size of a variable whose type is encoded as described by the
    `SupportedTypes::getTypeIdentifier` method of MFront"""
    d = (t >> 3) & 3
    if d == 0:
        d = 3 if hypothesis == "Tridimensional" else (
            1 if hypothesis.startswith("AxisymmetricalGeneralised") else 2)
    return {0: 1, 1: {1: 3, 2: 4, 3: 6}[d], 2: d, 3: {1: 3, 2: 5, 3: 9}[d]}[t & 7]


def initial_internal_state_variables(lib, prefixes, hypothesis, data):
    """initial values of all the components of the internal state
    variables, read in an MTest file (null by default)"""
    names = lib.strings(prefixes, "InternalStateVariables") or []
    if not names:
        return []
    types = lib.ints(prefixes, "InternalStateVariablesTypes", len(names))
    values = []
    for n, t in zip(names, types):
        size = variable_size(t, hypothesis)
        v = data["isvs"].get(n, 0.0)
        # element of an array given by the name of the array in the MTest
        # file: one value per element, or one value for all of them
        element = re.fullmatch(r"(\w+)\[(\d+)\]", n)
        if (n not in data["isvs"]) and element and (element.group(1) in data["isvs"]):
            v = data["isvs"][element.group(1)]
            if isinstance(v, list):
                i = int(element.group(2))
                v = v[i] if i < len(v) else 0.0
        if isinstance(v, list):
            values += (v + [0.0] * size)[:size]
        else:
            values += [v] + [0.0] * (size - 1) if size > 1 else [v]
    return values


def synthetic_material_property(name, default):
    """value of a material property which is not given by an MTest file,
    chosen from its name (the given default value otherwise). Those values
    only allow to run and compare the code of the law: they are not
    physical."""
    n = re.sub(r"\[\d+\]$", "", name)
    if n.startswith("YoungModulus") or n == "E":
        return 150e9
    if n.startswith("PoissonRatio") or n == "nu":
        return 0.3
    if n.startswith("ShearModulus") or n in ("Mu", "G"):
        return 60e9
    if n.startswith("ThermalExpansion"):
        return 1e-5
    if n == "MassDensity":
        return 7800.0
    if n in ("s0", "R0", "sy", "YieldStress", "YieldStrength",
             "ReferenceEquivalentStress"):
        return 100e6
    if n == "Rinf":
        return 200e6
    return default


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
    header = "# " + " ".join(cmd) + "\n"
    try:
        with args.compare_slots:
            r = subprocess.run(cmd, cwd=directory, capture_output=True,
                               text=True, timeout=600, env=env)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", "", header + "timeout\n"
    out = header + r.stdout + r.stderr
    lines = [l for l in out.splitlines() if l.startswith("K0=")]
    if r.returncode not in (0, 1) or not lines:
        return "CRASH", f"exit code {r.returncode}", out[-500:]
    if all(l.endswith("OK") for l in lines):
        return "OK", "", out
    # largest relative difference, infinite if the statuses differ
    difference = 0.0
    for l in lines:
        if l.endswith("OK"):
            continue
        if "same statuses=0" in l:
            difference = math.inf
        for v in re.findall(r"(?:forces|isvs|K|energies|sound) ([0-9.e+-]+)", l):
            difference = max(difference, float(v))
    if difference > 1e-5:
        return "MISMATCH", "", out
    if difference > 1e-10:
        return "ROUNDOFF", "", out
    return "INCONCLUSIVE", "", out


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
    # one library per interface, named after the material if any (for
    # instance libB4C-generic.so and libGenericParallelB4C-openmp.so)
    libs = sorted(l for l in os.listdir(src) if l.endswith(".so")) if os.path.isdir(src) else []
    batch = [l for l in libs if l.startswith("libGenericParallel")]
    references = [l for l in libs if not l.startswith("libGenericParallel")]
    if (len(batch) != 1) or (len(references) != 1):
        return rel, "BUILD", "libraries not found"
    reference = os.path.join("src", references[0])
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
    synthetic = [n for n in mps if (data is None) or (n not in data["mps"])]

    def material_properties(default):
        """option giving the material properties, the synthetic values
        being computed with the given default value"""
        values = [data["mps"][n] if n not in synthetic
                  else synthetic_material_property(n, default) for n in mps]
        return ["--material-properties=" + ",".join(repr(v) for v in values)]
    if data is not None and "Temperature" in data["esvs"]:
        options.append(f"--temperature={data['esvs']['Temperature']!r}")
    if (data is not None) and data["parameters"]:
        # real parameters changed by the MTest file (the other parameters
        # can't be changed by the comparison program)
        names = lib.strings(prefixes, "Parameters") or []
        types = lib.ints(prefixes, "ParametersTypes", len(names)) or []
        reals = {n for n, t in zip(names, types) if t == 0}
        parameters = [f"{n}={v!r}" for n, v in data["parameters"].items() if n in reals]
        if parameters:
            options.append("--parameters=" + ",".join(parameters))
    if lib.ushort(prefixes, "BehaviourType") == 2:
        options.append("--finite-strain=CAUCHY,DSIG_DF")
    if (data is not None) and data["isvs"]:
        isvs = initial_internal_state_variables(lib, prefixes, hypothesis, data)
        if isvs:
            options.append("--initial-internal-state-variables=" +
                           ",".join(repr(v) for v in isvs))
    if data is not None:
        esvs = lib.strings(prefixes, "ExternalStateVariables") or []
        if esvs:
            options.append("--external-state-variables=" + ",".join(
                repr(data["esvs"].get(n, data["esvs"].get("Temperature", 293.15)))
                for n in esvs))
        dt, amplitudes = loading(data)
    else:
        dt, amplitudes = None, []
    if not amplitudes:
        amplitudes = args.amplitudes
    # time increment of the comparison program by default
    dt = 1.0 if dt is None else dt
    # if the integration fails at all the points (divergence of the local
    # solver, for instance), smaller time increments and amplitudes are tried
    # all the amplitudes of the first loading are compared, so that the
    # inelastic part of the behaviour is tested: the results of the largest
    # amplitudes are inconclusive if the local solver diverges
    first = [(dt, a) for a in amplitudes]
    smallest = f"{min(float(a) for a in amplitudes) / 10:.3g}"
    fallbacks = [(dt / 10, a) for a in amplitudes] + [(dt / 1000, a) for a in amplitudes]
    fallbacks += [(dt, smallest), (dt / 1000, smallest)]
    function = f"{behaviour}_{hypothesis}"
    detail = os.path.basename(mtests[0]) if mtests else "no MTest file"
    if synthetic:
        detail += f", synthetic values of {len(synthetic)} material properties"
    # the synthetic values of the material properties which are not given
    # a value from their names are first 1, then 1e-3
    defaults = [1.0, 1e-3] if synthetic else [None]

    def run(time_increment, amplitude, mps_options):
        """first conclusive comparison for the types of tangent operators"""
        result = ("INCONCLUSIVE", "")
        for tangent in ["4", "1", "0"]:
            o = options + mps_options + [
                f"--time-increment={time_increment!r}",
                f"--amplitude={amplitude}", f"--tangent-operator={tangent}"]
            status, why, out = compare(args, directory, reference, batch, function, o)
            result = (status, out)
            if status != "INCONCLUSIVE":
                break
        return result

    # statuses from the best to the worst
    ranks = ["INCONCLUSIVE", "OK", "TIMEOUT", "ROUNDOFF", "MISMATCH", "CRASH"]
    status, logs = "INCONCLUSIVE", []
    for default in defaults:
        mps_options = material_properties(default) if mps else []
        for time_increment, amplitude in first:
            s, out = run(time_increment, amplitude, mps_options)
            logs.append(out)
            if ranks.index(s) > ranks.index(status):
                status = s
        for time_increment, amplitude in fallbacks:
            if status != "INCONCLUSIVE":
                break
            status, out = run(time_increment, amplitude, mps_options)
            logs.append(out)
        if status != "INCONCLUSIVE":
            break
    with open(os.path.join(directory, "compare.log"), "w") as f:
        f.write("".join(logs))
    return rel, status, detail


def main():
    p = argparse.ArgumentParser()
    p.add_argument("configuration")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--compare-jobs", type=int, default=None,
                   help="maximal number of comparisons run at the same time "
                   "(number of jobs by default)")
    p.add_argument("--summary", default="summary.tsv")
    p.add_argument("--pattern", default=".")
    p.add_argument("--tests", default=os.path.join(PROTO, "tfel-dev", "mfront", "tests"))
    p.add_argument("--points", type=int, default=1009)
    p.add_argument("--amplitudes", default="1e-4,1e-3,1e-2")
    p.add_argument("--compile-sweep", default=None,
                   help="summary of a compilation sweep: only the files "
                   "compiled successfully are tested")
    args = p.parse_args()
    args.amplitudes = args.amplitudes.split(",")
    args.compare_slots = threading.BoundedSemaphore(
        args.compare_jobs if args.compare_jobs is not None else args.jobs)
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
    summary = os.path.join(args.output, args.summary)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool, open(summary, "w") as s:
        for rel, status, detail in pool.map(lambda f: check(args, f), files):
            s.write(f"{rel}\t{status}\t{detail}\n")
            s.flush()
            print(f"{rel}\t{status}\t{detail}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
