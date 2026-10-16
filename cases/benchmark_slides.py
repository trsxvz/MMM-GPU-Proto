#!/usr/bin/env python3
"""Benchmark of the examples of mfem-mgis-examples and of mm-opera-hpc for
the slides of the prototype.

For each case, the wall time of the whole program is measured on the CPU
(MPI, behaviours integrated by batches and finite element kernels on the
host: `--execution cpu`) and with at least one stage on the GPU (`--execution
gpu`, or `--execution mixed` with the linear solver on the CPU), for several
linear solvers given by the linear solver override of mfem-mgis
(`--linear-solver-override`, with their preconditioners) and several non
linear solvers given by the non linear solver override
(`--non-linear-solver-override`), the program keeping its convergence
criteria. Each run is written as one line of `runs.csv` in the results
directory, its output in `logs/`.

The programs are the ones of `build/bench/examples` and `build/bench/opera`,
built against `install/mfem-mgis-bench` (see `docs/benchmark-slides.md`).

usage:
  benchmark_slides.py list
  benchmark_slides.py run RESULTS CASE SIDE SOLVER [--non-linear NAME]
                      [--np N] [--size SIZE] [--repeat R] [--dry-run]
  benchmark_slides.py best RESULTS
"""

import argparse
import csv
import datetime
import os
import pathlib
import shlex
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "build" / "bench" / "examples"
OPERA = ROOT / "build" / "bench" / "opera"
OPENMPI = pathlib.Path(
    "/home/t/Bureau/spack-mmm/opt/spack/linux-x86_64_v3/"
    "openmpi-5.0.10-57cvs5dydumy2shauuhl64r4ybhyn2iw")


class Case:
    """a program, its directory and its options for each size"""

    def __init__(self, directory, program, options, sizes, tolerance,
                 symmetric, parallel_option=None, sequential_only=False):
        # directory of the program, where it runs (meshes and references)
        self.directory = directory
        # program, relative to the directory
        self.program = program
        # options of all the runs
        self.options = options
        # options of each size, the first one being the default size
        self.sizes = sizes
        # relative tolerance of the linear solver of the program, given to
        # all the linear solvers of the benchmark
        self.tolerance = tolerance
        # if true, the tangent operator is symmetric (conjugate gradient)
        self.symmetric = symmetric
        # option selecting the parallel computations, if the program
        # computes in sequential by default
        self.parallel_option = parallel_option
        # if true, the program only runs on one process
        self.sequential_only = sequential_only


# sizes: options of the size of the benchmark ("bench") and of a small
# size checking the scripts ("small"); steps shortened identically on the
# CPU and on the GPU
CASES = {
    # the time step of ex1 is 1 / nbsteps (100 steps by default): with 20
    # steps or less, the Newton solver does not converge at the onset of
    # plasticity
    "ex1": Case(EXAMPLES / "ex1", "./uniaxial_tensile_test",
                ["--no-post-processing"],
                {"small": ["-r", "2", "-ns", "25"],
                 "bench": ["-r", "6", "-ns", "25"]},
                1e-12, True, parallel_option="--parallel"),
    "ex2": Case(EXAMPLES / "ex2", "./ssna303", [],
                {"small": ["--nbsteps", "5"],
                 "bench": ["-r", "2", "--nbsteps", "50"]},
                1e-12, False, parallel_option="--parallel"),
    "ex3": Case(EXAMPLES / "ex3", "./two_layer_cube",
                ["--test-case", "1"],
                {"small": ["-r", "2"], "bench": ["-r", "4"]},
                1e-12, True, parallel_option="--parallel"),
    "ex4": Case(EXAMPLES / "ex4", "./ssna303_3d_hypre",
                ["--end-time", "1"],
                {"small": ["--nbsteps", "2"],
                 "bench": ["-r", "1", "--nbsteps", "10"]},
                1e-6, False),
    "ex5": Case(EXAMPLES / "ex5", "./satoh", ["--no-post-processing"],
                {"small": [], "bench": ["-r", "2"]},
                1e-12, True, parallel_option="--parallel"),
    "ex6": Case(EXAMPLES / "ex6", "./rve", ["--no-post-processing"],
                {"small": ["-r", "1"], "bench": ["-r", "3"]},
                1e-12, True),
    "ex7": Case(EXAMPLES / "ex7", "./mox2",
                ["--no-post-processing", "--order", "2", "--mesh",
                 "mesh/OneSphere.msh"],
                {"small": ["--nbsteps", "5"],
                 "bench": ["-r", "1", "--nbsteps", "40"]},
                1e-12, False),
    "ex8": Case(EXAMPLES / "ex8", "./rjh_plate",
                ["-m", "assemblage_hexa.msh"],
                {"small": ["-ns", "1"], "bench": ["-r", "1", "-ns", "5"]},
                1e-12, False),
    "polycrystal": Case(OPERA / "polycrystal", "./uniaxial-polycrystal", [],
                        {"small": ["--nstep", "2", "--duration", "2"],
                         "bench": ["-r", "1", "--nstep", "10",
                                   "--duration", "10"]},
                        1e-12, False),
    "bubble": Case(OPERA / "bubble", "./test-bubble",
                   ["--post-processing", "0"],
                   {"small": [], "bench": []}, 1e-12, True),
    "cermet": Case(OPERA / "cermet", "./cermet", ["--postprocessing", "0"],
                   {"small": ["--nstep", "1", "--duration", "0.4"],
                    "bench": ["--nstep", "10", "--duration", "4"]},
                   1e-12, True),
}

# linear solvers, given by the linear solver override of mfem-mgis (name
# and options, in the syntax of the data of TFEL), {tolerance} being
# replaced by the relative tolerance of the case, {hypre}, {belos} and
# {petsc} by the conjugate gradient or GMRES of each library, and {amgx}
# by the solver of AmgX
SOLVERS = {
    # Krylov solvers of hypre preconditioned by BoomerAMG
    "hypre-amg": '{{Name: "{hypre}", Options: {{Tolerance: {tolerance}, '
                 'MaximumNumberOfIterations: 5000, Preconditioner: {{Name: '
                 '"HypreBoomerAMG", Options: {{VerbosityLevel: 0}}}}}}}}',
    # direct solvers
    "mumps": '{{Name: "MUMPSSolver"}}',
    "cudss": '{{Name: "CuDSSSolver"}}',
    # Krylov solvers of Belos preconditioned by MueLu (Jacobi smoother,
    # faster than the default Gauss-Seidel smoother on a GPU)
    "belos-muelu": '{{Name: "BelosSolver", Options: {{Solver: "{belos}", '
                   'RelativeTolerance: {tolerance}, '
                   'MaximumNumberOfIterations: 5000, Preconditioner: {{Name: '
                   '"MueLu", Options: {{Parameters: {{"smoother: type": '
                   '"RELAXATION", "smoother: params": {{"relaxation: type": '
                   '"Jacobi", "relaxation: sweeps": 2, "relaxation: damping '
                   'factor": 0.6}}}}}}}}}}}}',
    # AmgX: its conjugate gradient preconditioned by its multigrid, or the
    # GMRES solver of MFEM preconditioned by one cycle of this multigrid
    "amgx": '{amgx}',
    # PETSc, preconditioned by its algebraic multigrid
    "petsc-gamg": '{{Name: "PETScLinearSolver", Options: {{Options: '
                  '{{ksp_type: "{petsc}", pc_type: "gamg", ksp_rtol: '
                  '{tolerance}, ksp_max_it: 5000}}}}}}',
}

# non linear solvers, given by the non linear solver override of mfem-mgis,
# {linear} being replaced by the linear solver of the run (the hypre-amg one
# for the default linear solver), which also solves the elastic operator of
# NewtonK0
NON_LINEAR_SOLVERS = {
    "newton": '{{Name: "Newton"}}',
    # modified Newton method: the elastic operator, assembled once, is the
    # operator of all the iterations
    "modified-newton": '{{Name: "FixedPoint", Options: {{Method: '
                       '"ModifiedNewton"}}}}',
    # Newton's method preconditioned by the elastic operator
    "newton-k0": '{{Name: "FixedPoint", Options: {{Method: "NewtonK0", '
                 'ElasticOperatorLinearSolver: {linear}}}}}',
}

# sides: execution options of mfem-mgis
SIDES = {
    "cpu": ["--execution", "cpu"],
    "gpu": ["--execution", "gpu"],
    # behaviours and finite elements on the GPU, linear solver on the CPU
    "mixed": ["--execution", "mixed", "--linear-solvers-execution", "cpu"],
}

FIELDS = ["date", "case", "size", "side", "solver", "non_linear",
          "processes", "repeat", "status", "seconds", "command"]


def amgx(case):
    """the solver of AmgX of a case (`AmgX` alone is one cycle of its
    multigrid, without convergence criteria)"""
    if case.symmetric:
        return ('{Name: "AmgXPCG", Options: {RelativeTolerance: '
                + repr(case.tolerance) + ', MaximumNumberOfIterations: 5000}}')
    return ('{Name: "GMRESSolver", Options: {RelativeTolerance: '
            + repr(case.tolerance) + ', MaximumNumberOfIterations: 5000, '
            'Preconditioner: {Name: "AmgX"}}}')


def override(case, solver):
    """the linear solver override of a case"""
    template = SOLVERS[solver]
    return template.format(
        tolerance=repr(case.tolerance),
        hypre="HyprePCG" if case.symmetric else "HypreGMRES",
        belos="CG" if case.symmetric else "GMRES",
        petsc="cg" if case.symmetric else "gmres",
        amgx=amgx(case))


def non_linear_override(case, solver, non_linear):
    """the non linear solver override of a case"""
    linear = override(case, "hypre-amg" if solver == "default" else solver)
    return NON_LINEAR_SOLVERS[non_linear].format(linear=linear)


def command(case, side, solver, non_linear, processes, size):
    """the command of a run"""
    c = [case.program] + case.options + case.sizes[size] + SIDES[side]
    if solver != "default":
        c += ["--linear-solver-override", override(case, solver)]
    if non_linear != "default":
        c += ["--non-linear-solver-override",
              non_linear_override(case, solver, non_linear)]
    # the parallel computations, also on one process (the path of the
    # applications, and of the linear solvers of hypre, Trilinos...)
    if case.parallel_option:
        c.append(case.parallel_option)
    if processes > 1:
        # one process per physical core, cores 0 to 7 (the twins 8 to 15
        # stay idle)
        c = [str(OPENMPI / "bin" / "mpirun"), "-n", str(processes),
             "--bind-to", "core", "--map-by", "core"] + c
    return c


def run(results, case_name, side, solver, non_linear, processes, size,
        repeat, dry_run=False):
    """run a case, append the result to results/runs.csv"""
    case = CASES[case_name]
    if case.sequential_only and processes > 1:
        raise RuntimeError(case_name + " only runs on one process")
    c = command(case, side, solver, non_linear, processes, size)
    if not dry_run:
        (results / "logs").mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    for r in range(repeat):
        name = (f"{case_name}-{size}-{side}-{solver}-{non_linear}"
                f"-np{processes}-{r}")
        print(shlex.join(c))
        if dry_run:
            continue
        start = time.perf_counter()
        with open(results / "logs" / (name + ".log"), "w") as log:
            p = subprocess.run(c, cwd=case.directory, stdout=log,
                               stderr=subprocess.STDOUT, env=env)
        seconds = time.perf_counter() - start
        row = {"date": datetime.datetime.now().isoformat(timespec="seconds"),
               "case": case_name, "size": size, "side": side,
               "solver": solver, "non_linear": non_linear,
               "processes": processes, "repeat": r,
               "status": p.returncode, "seconds": f"{seconds:.3f}",
               "command": shlex.join(c)}
        path = results / "runs.csv"
        new = not path.exists()
        with open(path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            if new:
                w.writeheader()
            w.writerow(row)
        print(f"{name}: status {p.returncode}, {seconds:.2f} s")


def best(results):
    """best time of each case on each side, among the successful runs"""
    rows = list(csv.DictReader(open(results / "runs.csv")))
    table = {}
    for row in rows:
        if row["status"] != "0":
            continue
        key = (row["case"], row["size"],
               "cpu" if row["side"] == "cpu" else "gpu")
        t = float(row["seconds"])
        # the files written before the non linear axis have no such column
        solver = row["solver"] + "/" + (row.get("non_linear") or "default")
        if (key not in table) or (t < table[key][0]):
            table[key] = (t, solver, row["side"], row["processes"])
    print("case,size,cpu_seconds,cpu_solver,cpu_processes,gpu_seconds,"
          "gpu_solver,gpu_side,speedup")
    for case, size in sorted({(k[0], k[1]) for k in table}):
        cpu = table.get((case, size, "cpu"))
        gpu = table.get((case, size, "gpu"))
        speedup = f"{cpu[0] / gpu[0]:.2f}" if (cpu and gpu) else ""
        print(",".join([
            case, size,
            f"{cpu[0]:.2f}" if cpu else "", cpu[1] if cpu else "",
            cpu[3] if cpu else "",
            f"{gpu[0]:.2f}" if gpu else "", gpu[1] if gpu else "",
            gpu[2] if gpu else "", speedup]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("results", type=pathlib.Path)
    r.add_argument("case", choices=sorted(CASES))
    r.add_argument("side", choices=sorted(SIDES))
    r.add_argument("solver", choices=["default"] + sorted(SOLVERS))
    r.add_argument("--non-linear", default="default",
                   choices=["default"] + sorted(NON_LINEAR_SOLVERS))
    r.add_argument("--np", type=int, default=1)
    r.add_argument("--size", default="small")
    r.add_argument("--repeat", type=int, default=1)
    r.add_argument("--dry-run", action="store_true")
    b = sub.add_parser("best")
    b.add_argument("results", type=pathlib.Path)
    args = parser.parse_args()
    if args.action == "list":
        for name, case in CASES.items():
            print(name, case.program, "sizes:", ", ".join(case.sizes))
        print("solvers:", ", ".join(["default"] + sorted(SOLVERS)))
        print("non linear solvers:",
              ", ".join(["default"] + sorted(NON_LINEAR_SOLVERS)))
        print("sides:", ", ".join(sorted(SIDES)))
    elif args.action == "run":
        run(args.results, args.case, args.side, args.solver, args.non_linear,
            args.np, args.size, args.repeat, args.dry_run)
    else:
        best(args.results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
