#!/usr/bin/env python3
"""Benchmark of the examples of mfem-mgis-examples and of mm-opera-hpc for
the slides of the prototype.

For each case, size and precision (double or float), the wall time of the
whole program is measured on the CPU (MPI on 8 or 16 processes, behaviours
integrated by batches and finite element kernels on the host: `--execution
cpu`) and with at least one major phase on the GPU (`--execution gpu`, or
`--execution mixed` with the linear solver on the CPU). Each case has a
list of configurations per side (table `CONFIGURATIONS`): the solvers of
the program (`native`), and the linear solvers, preconditioners and non
linear solvers of the libraries of mfem-mgis which may suit the case, given
by the solver overrides of mfem-mgis (`--linear-solver-override`,
`--non-linear-solver-override`), by options of the program, or by
environment variables (tangent operator applied without assembly). The
program keeps its convergence criteria.

Each run is written as one line of `runs.csv` in the results directory,
its output in `logs/`. The programs are the ones of `build/bench` (double)
and `build/bench-float` (float), built against `install/mfem-mgis-bench`
and `install/mfem-mgis-bench-float` (see `docs/benchmark-slides.md`).

Stages of the measurements: `calibrate` runs each case once with the
solvers of the program, on 8 processes and on the GPU; `screen` runs all
the configurations of a case at the screening size (the mesh of the
measurements and one or two steps of the same time step, or for ex1 and
ex7 a coarser mesh or fewer steps); `run` runs one
configuration, for the repetitions of the best ones; `best` and
`cases/benchmark_plots.py` summarize the results.

usage:
  benchmark_slides.py list [CASE]
  benchmark_slides.py run RESULTS CASE CONFIGURATION [--precision P]
                      [--np N] [--size SIZE] [--repeat R] [--timeout S]
                      [--dry-run]
  benchmark_slides.py calibrate RESULTS [--precision P] [--size SIZE]
                      [--cases LIST] [--timeout S] [--dry-run]
  benchmark_slides.py screen RESULTS CASE [--precision P] [--sides LIST]
                      [--size SIZE] [--timeout S] [--dry-run]
  benchmark_slides.py finals RESULTS CASE [--precision P] [--size SIZE]
                      [--screen-size SIZE] [--repeat R] [--timeout S]
                      [--dry-run]
  benchmark_slides.py best RESULTS
"""

import argparse
import contextlib
import csv
import datetime
import fcntl
import os
import pathlib
import re
import shlex
import signal
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
# builds of the examples and of mm-opera-hpc for each precision
BUILDS = {"double": ROOT / "build" / "bench",
          "float": ROOT / "build" / "bench-float"}
OPENMPI = pathlib.Path(
    "/home/t/Bureau/spack-mmm/opt/spack/linux-x86_64_v3/"
    "openmpi-5.0.10-57cvs5dydumy2shauuhl64r4ybhyn2iw")


class Size:
    """discretization and time steps of a run: order of the elements and
    options of the program"""

    def __init__(self, order, options):
        self.order = order
        self.options = options


class Case:
    """a program, its directory and its sizes"""

    def __init__(self, project, subdirectory, program, options, sizes,
                 tolerance, symmetric, parallel_option=None):
        # "examples" or "opera"; the directory of the program, where it
        # runs (meshes and references), is given by the precision
        self.project = project
        self.subdirectory = subdirectory
        # program, relative to the directory
        self.program = program
        # options of all the runs
        self.options = options
        # sizes (class Size), "small" checking the scripts, "screen" for
        # the screening and "bench" for the measurements
        self.sizes = sizes
        # relative tolerance of the linear solver of the program, given to
        # the linear solvers of the overrides
        self.tolerance = tolerance
        # if true, the tangent operator is symmetric (conjugate gradient)
        self.symmetric = symmetric
        # option selecting the parallel computations, if the program
        # computes in sequential by default
        self.parallel_option = parallel_option

    def directory(self, precision):
        """the directory of the program for a precision"""
        return BUILDS[precision] / self.project / self.subdirectory


def _sizes(small, screen, bench, order=1):
    """the sizes of a case, all at the same order"""
    return {"small": Size(order, small), "screen": Size(order, screen),
            "bench": Size(order, bench)}


# The sizes "bench" fill the 8 cores and the GPU (about one million
# unknowns, one or several waves of the multiprocessors of the GPU for the
# kernels of the elements), as given by the calibration; "screen" keeps
# the mesh and the time step of "bench", with one or two steps when the
# program has an option of the end time, except for ex1 and ex7 which
# have no such option (see below). The orders are the ones of the programs
# (2 for ex5, ex7 and the bubble), always given explicitly.
CASES = {
    # the time step of ex1 is 1 / nbsteps (100 steps by default); its
    # loading is cyclic, reversed at t = 0.3. Without --prediction, the
    # increment of the imposed displacement is concentrated in the elements
    # next to the loaded face: with 4 refinements or more, the Newton solver
    # does not converge at the reversal (with 6, not even at the first
    # step). ex1 has no option of the end time: "screen" uses the mesh of 5
    # refinements (107811 unknowns), "bench" the one of 6 (823875).
    "ex1": Case("examples", "ex1", "./uniaxial_tensile_test",
                ["--no-post-processing", "--prediction"],
                _sizes(["-r", "2", "-ns", "25"], ["-r", "5", "-ns", "25"],
                       ["-r", "6", "-ns", "25"]),
                1e-12, True, parallel_option="--parallel"),
    "ex2": Case("examples", "ex2", "./ssna303", [],
                _sizes(["--nbsteps", "5"],
                       ["-r", "4", "--nbsteps", "2", "--end-time", "0.04"],
                       ["-r", "4", "--nbsteps", "50"]),
                1e-12, False, parallel_option="--parallel"),
    "ex3": Case("examples", "ex3", "./two_layer_cube", ["--test-case", "1"],
                _sizes(["-r", "2"], ["-r", "4"], ["-r", "4"]),
                1e-12, True, parallel_option="--parallel"),
    "ex4": Case("examples", "ex4", "./ssna303_3d_hypre", [],
                _sizes(["--nbsteps", "2", "--end-time", "1"],
                       ["-r", "1", "--nbsteps", "2", "--end-time", "0.2"],
                       ["-r", "1", "--nbsteps", "10", "--end-time", "1"]),
                1e-6, False),
    "ex5": Case("examples", "ex5", "./satoh", ["--no-post-processing"],
                _sizes([], ["-r", "4"], ["-r", "4"], order=2),
                1e-12, True, parallel_option="--parallel"),
    "ex6": Case("examples", "ex6", "./rve", ["--no-post-processing"],
                _sizes(["-r", "1"], ["-r", "4"], ["-r", "4"]),
                1e-12, True),
    # ex7 has no option of the end time: with 1 refinement, it has only
    # 21012 unknowns; it is measured with 2 refinements, and screened with
    # 10 time steps instead of 40
    "ex7": Case("examples", "ex7", "./mox2",
                ["--no-post-processing", "--mesh", "mesh/OneSphere.msh"],
                _sizes(["--nbsteps", "5"], ["-r", "2", "--nbsteps", "10"],
                       ["-r", "2", "--nbsteps", "40"], order=2),
                1e-12, False),
    "ex8": Case("examples", "ex8", "./rjh_plate",
                ["-m", "assemblage_hexa.msh"],
                # "screen": the first of the 5 steps of the measurements
                _sizes(["-ns", "1"],
                       ["-r", "1", "-et", "2e4", "-tr", "1e5", "-ns", "1"],
                       ["-r", "1", "-ns", "5"]),
                1e-12, False),
    # the polycrystal and the cermet run for short durations, the same in
    # double and in float: the grain law of UO2 fails in float beyond about
    # 0.8 s (rounding of its elastic strain, results/precision-polycrystal)
    "polycrystal": Case("opera", "polycrystal", "./uniaxial-polycrystal",
                        [],
                        _sizes(["--nstep", "2", "--duration", "0.4"],
                               ["-r", "1", "--nstep", "1", "--duration",
                                "0.2"],
                               ["-r", "1", "--nstep", "4", "--duration",
                                "0.8"]),
                        1e-12, False),
    "bubble": Case("opera", "bubble", "./test-bubble",
                   ["--post-processing", "0"],
                   # the program solves with a relative tolerance of 1e-10
                   _sizes([], [], [], order=2), 1e-10, True),
    "cermet": Case("opera", "cermet", "./cermet", ["--postprocessing", "0"],
                   _sizes(["--nstep", "1", "--duration", "0.4"],
                          ["-r", "1", "--nstep", "1", "--duration", "0.2"],
                          ["-r", "1", "--nstep", "2", "--duration", "0.4"]),
                   # finite strain behaviours: tangent operator not
                   # symmetric in general (same grain law as the polycrystal)
                   1e-12, False),
}

# Linear solvers, given by the linear solver override of mfem-mgis (name
# and options, in the syntax of the data of TFEL), {tolerance} being
# replaced by the relative tolerance of the case and {hypre}, {mfem},
# {belos}, {petsc}, {ginkgo}, {magma} by the conjugate gradient or GMRES of
# each library (the case being symmetric or not). The value gives the
# description and the precisions in which the solver is available.
_JACOBI = ('{{"relaxation: type": "Jacobi", "relaxation: sweeps": 2, '
           '"relaxation: damping factor": 0.6}}')
_MUELU_JACOBI = ('{{Name: "MueLu", Options: {{Parameters: {{"smoother: '
                 'type": "RELAXATION", "smoother: params": ' + _JACOBI +
                 ', "reuse: type": "RP"}}}}}}')
_MUELU_CHEBYSHEV = ('{{Name: "MueLu", Options: {{Parameters: {{"smoother: '
                    'type": "CHEBYSHEV", "smoother: params": {{"chebyshev: '
                    'eigenvalue max iterations": 30, "chebyshev: boost '
                    'factor": 1.2}}, "reuse: type": "RP"}}}}}}')
_MUELU_DEFAULT = ('{{Name: "MueLu", Options: {{Parameters: {{"reuse: type": '
                  '"RP"}}}}}}')
_BOTH = ("double", "float")
_DOUBLE = ("double",)


def _krylov(name, preconditioner, tolerance_key="RelativeTolerance",
            iterations=5000):
    """a Krylov solver with a preconditioner"""
    return ('{{Name: "' + name + '", Options: {{' + tolerance_key +
            ': {tolerance}, MaximumNumberOfIterations: ' + str(iterations) +
            ((', Preconditioner: ' + preconditioner) if preconditioner
             else '') + '}}}}')


LINEAR_SOLVERS = {
    # Krylov solvers of hypre preconditioned by BoomerAMG (systems options
    # for vector unknowns, the default of mfem-mgis)
    "hypre-amg": (_krylov("{hypre}", '{{Name: "HypreBoomerAMG", Options: '
                          '{{VerbosityLevel: 0}}}}', "Tolerance"), _BOTH),
    # the same with one level of aggressive coarsening (cheaper setup)
    "hypre-amg-agg": (_krylov("{hypre}", '{{Name: "HypreBoomerAMG", Options: '
                              '{{VerbosityLevel: 0, '
                              'AggressiveCoarseningLevels: 1}}}}',
                              "Tolerance"), _BOTH),
    # Krylov solvers of MFEM preconditioned by BoomerAMG
    "mfem-amg": (_krylov("{mfem}", '{{Name: "HypreBoomerAMG", Options: '
                         '{{VerbosityLevel: 0}}}}'), _BOTH),
    # Krylov solvers of MFEM with the Jacobi preconditioner (also with the
    # tangent operator applied without assembly)
    "cg-jacobi": (_krylov("{mfem}", '{{Name: "OperatorJacobiSmoother"}}',
                          iterations=100000), _BOTH),
    # Krylov solvers of Belos preconditioned by MueLu with the rigid body
    # modes, its hierarchy reused between operators: Jacobi smoother (GPU),
    # Chebyshev smoother (GPU, thin plates), default smoother (symmetric
    # Gauss-Seidel, CPU)
    "muelu-jacobi": ('{{Name: "BelosSolver", Options: {{Solver: "{belos}", '
                     'RelativeTolerance: {tolerance}, '
                     'MaximumNumberOfIterations: 5000, Preconditioner: ' +
                     _MUELU_JACOBI + '}}}}', _BOTH),
    "muelu-chebyshev": ('{{Name: "BelosSolver", Options: {{Solver: '
                        '"{belos}", RelativeTolerance: {tolerance}, '
                        'MaximumNumberOfIterations: 5000, Preconditioner: ' +
                        _MUELU_CHEBYSHEV + '}}}}', _BOTH),
    "muelu-sgs": ('{{Name: "BelosSolver", Options: {{Solver: "{belos}", '
                  'RelativeTolerance: {tolerance}, '
                  'MaximumNumberOfIterations: 5000, Preconditioner: ' +
                  _MUELU_DEFAULT + '}}}}', _BOTH),
    # AmgX: its conjugate gradient preconditioned by its multigrid, or the
    # GMRES solver of MFEM preconditioned by one cycle of this multigrid;
    # in double precision, or computing in single precision within the
    # double build (mfem-mgis refuses AmgX in a float build)
    "amgx": ("{amgx}", _DOUBLE),
    "amgx-single": ("{amgx_single}", _DOUBLE),
    # aggregation multigrid of AmgX (block Jacobi smoother): the classical
    # one fails at its setup on ex3 (CUDA kernel launch error)
    "amgx-agg": ("{amgx_agg}", _DOUBLE),
    # direct solvers, the analysis being reused while the pattern of the
    # matrix does not change
    "cudss": ('{{Name: "CuDSSSolver", Options: {{ReuseAnalysis: true{spd}}}}}',
              _BOTH),
    "mumps": ('{{Name: "MUMPSSolver", Options: {{ReuseAnalysis: true{sym}}}}}',
              _BOTH),
    # Tacho: Cholesky factorization by default, LU for the non symmetric
    # cases
    "tacho": ('{{Name: "Amesos2Solver", Options: {{Solver: "Tacho", '
              'Parameters: {{method: "{tacho}"}}}}}}', _BOTH),
    "strumpack": ('{{Name: "STRUMPACKSolver", Options: {{ReuseAnalysis: '
                  'true}}}}', _DOUBLE),
    # PETSc with its algebraic multigrid, stopped on the true residual
    "petsc-gamg": ('{{Name: "PETScLinearSolver", Options: {{Options: '
                   '{{ksp_type: "{petsc}", pc_type: "gamg", ksp_rtol: '
                   '{tolerance}, ksp_max_it: 5000, ksp_norm_type: '
                   '"unpreconditioned"}}}}}}', _BOTH),
    # Ginkgo: Krylov solvers with a Jacobi preconditioner on the GPU (the
    # sparse solvers of MAGMA are only available in sequential, and all the
    # programs compute in parallel)
    "ginkgo-jacobi": (_krylov("{ginkgo}", '{{Name: "GinkgoJacobi"}}',
                              iterations=100000), _BOTH),
}

# Non linear solvers, given by the non linear solver override of
# mfem-mgis, {tangent} being replaced by the iterative solver of the tangent
# operator of NewtonK0 (conjugate gradient or GMRES); the linear solver of
# the elastic operator of NewtonK0 is set up once. The inexact methods need
# more iterations than the exact Newton method for which the programs set
# their maximum number of iterations (10 for ex2, which the inexact Newton
# method reaches while it still converges), hence the maximum of 50
_EW = ('AdaptiveLinearSolverTolerance: {{EisenstatWalker1: '
       '{{InitialTolerance: 0.1}}}}, MaximumNumberOfIterations: 50')


def _newton_k0(elastic_operator_solver):
    """Newton's method preconditioned by the elastic operator, computed once
    by integrating the behaviours (some behaviours, such as the ones of the
    cermet, do not provide the prediction operator)"""
    return ('{{Name: "FixedPoint", Options: {{Method: "NewtonK0", '
            'ElasticOperatorComputation: "Integration", '
            'ElasticOperatorUpdate: "Once", ElasticOperatorLinearSolver: ' +
            elastic_operator_solver + ', TangentLinearSolver: "{tangent}", ' +
            _EW + '}}}}')


NON_LINEAR_SOLVERS = {
    # inexact Newton method (Eisenstat and Walker)
    "newton-ew": ('{{Name: "Newton", Options: {{' + _EW + '}}}}', _BOTH),
    # NewtonK0, the elastic operator being inverted by one cycle of
    # BoomerAMG, of MueLu, of AmgX, or factorized by cuDSS
    "k0-amg": (_newton_k0('"HypreBoomerAMG"'), _BOTH),
    "k0-muelu": (_newton_k0('{{Name: "MueLu", Options: {{Parameters: '
                            '{{"smoother: type": "RELAXATION", "smoother: '
                            'params": ' + _JACOBI + '}}}}}}'), _BOTH),
    # the default smoother of MueLu (symmetric Gauss-Seidel), for the CPU
    "k0-muelu-sgs": (_newton_k0('{{Name: "MueLu"}}'), _BOTH),
    "k0-amgx": (_newton_k0('"AmgX"'), _DOUBLE),
    "k0-cudss": (_newton_k0('"CuDSSSolver"'), _BOTH),
    # modified Newton method: the tangent operator of the first iteration
    # of each resolution, given to the linear solver of the problem (a
    # factorization reused by the iterations with a direct solver)
    "modified-newton": ('{{Name: "FixedPoint", Options: {{Method: '
                        '"ModifiedNewton"}}}}', _BOTH),
}

# tangent operator applied without assembly (kept with --execution)
MATRIX_FREE = {"MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION": "MatrixFree"}

# With the GPU, hypre allocates at once a pool of 4 GiB (Umpire), which
# leaves too little memory to Trilinos, AmgX or cuDSS for the largest
# cases (failed allocations on ex3, ex4, ex6): with their solvers, the
# initial pool of hypre is 1 GiB (it grows if needed)
SMALL_HYPRE_POOL = {"MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE": "1"}
_NOT_HYPRE_LINEAR = {"muelu-jacobi", "muelu-chebyshev", "muelu-sgs", "amgx",
                     "amgx-single", "amgx-agg", "cudss", "tacho",
                     "ginkgo-jacobi"}
_NOT_HYPRE_NON_LINEAR = {"k0-muelu", "k0-muelu-sgs", "k0-amgx", "k0-cudss",
                         "modified-newton"}

# sides: execution options of mfem-mgis
SIDES = {
    "cpu": ["--execution", "cpu"],
    "gpu": ["--execution", "gpu"],
    # behaviours and finite elements on the GPU, linear solver on the CPU
    "mixed": ["--execution", "mixed", "--linear-solvers-execution", "cpu"],
}


class Configuration:
    """solvers of a run: linear solver override, non linear solver
    override, options of the program and environment variables"""

    def __init__(self, side, linear=None, non_linear=None, options=(),
                 environment=None, processes=None, tier=1):
        self.side = side
        # 1: screened first; 2: variants dominated in the measurements of
        # the solvers (docs/solver-recommendations.md), screened in a second
        # pass
        self.tier = tier
        self.linear = linear
        self.non_linear = non_linear
        self.options = list(options)
        self.environment = environment or {}
        # numbers of processes: 8 and 16 (SMT) on the CPU, 1 on the GPU
        self.processes = processes or ([8, 16] if side == "cpu" else [1])

    def precisions(self):
        """the precisions in which the solvers are available"""
        p = set(_BOTH)
        if self.linear is not None:
            p &= set(LINEAR_SOLVERS[self.linear][1])
        if self.non_linear is not None:
            p &= set(NON_LINEAR_SOLVERS[self.non_linear][1])
        return p

    def label(self):
        """a short description, for the tables and the graphs"""
        parts = [self.linear or "native"]
        if self.non_linear is not None:
            parts.append(self.non_linear)
        if self.environment == MATRIX_FREE:
            parts.append("matrix-free")
        if self.options:
            parts.append(" ".join(self.options))
        return ", ".join(parts)


def _cpu_configurations(nonlinear, direct_cpu=True):
    """the solvers of the program on 8 and 16 processes, strong alternatives
    on 8 processes, and the non linear solvers for the non linear
    computations; `direct_cpu`: if false, MUMPS is left out (its
    factorization of a large 3D case exceeds the memory of the machine)"""
    c = {"cpu-native": Configuration("cpu"),
         "cpu-hypre-amg": Configuration("cpu", "hypre-amg", processes=[8]),
         "cpu-muelu-sgs": Configuration("cpu", "muelu-sgs", processes=[8]),
         "cpu-petsc-gamg": Configuration("cpu", "petsc-gamg",
                                         processes=[8])}
    if direct_cpu:
        c["cpu-mumps"] = Configuration("cpu", "mumps", processes=[8])
    if nonlinear:
        # the non linear solvers tried with the GPU, on the CPU too
        c["cpu-hypre-amg-newton-ew"] = Configuration(
            "cpu", "hypre-amg", "newton-ew", processes=[8])
        c["cpu-hypre-amg-k0-amg"] = Configuration(
            "cpu", "hypre-amg", "k0-amg", processes=[8])
        c["cpu-k0-muelu"] = Configuration("cpu", "hypre-amg", "k0-muelu-sgs",
                                          processes=[8])
    return c


def _gpu_configurations(nonlinear, direct=True, native_on_cpu=True,
                        native_on_gpu=True, direct_cpu=True):
    """the configurations on the GPU: the solvers of the program, the
    linear solvers suited to the GPU, and for the non linear computations
    the inexact Newton method and NewtonK0 (also with the tangent operator
    applied without assembly); `direct`: if false, the size of the case is
    too large for the direct solvers on the GPU; `native_on_cpu`: if false,
    the linear solver of the program computes where MFEM computes (Krylov
    solver of MFEM), and can't be placed on the CPU with the elements on
    the GPU; `native_on_gpu`: if false, the linear solver of the program
    only computes on the CPU (MUMPS)"""
    c = {"gpu-native": Configuration("gpu")} if native_on_gpu else {}
    c.update({
         "gpu-hypre-amg": Configuration("gpu", "hypre-amg"),
         "gpu-hypre-amg-agg": Configuration("gpu", "hypre-amg-agg", tier=2),
         "gpu-mfem-amg": Configuration("gpu", "mfem-amg", tier=2),
         "gpu-muelu-jacobi": Configuration("gpu", "muelu-jacobi"),
         # the conjugate gradient of Belos stops on a non positive p.A.p
         # on ex1 with the Jacobi smoother (GPU) and the symmetric
         # Gauss-Seidel one (CPU), and on the bubble with the Jacobi one:
         # the Chebyshev smoother is screened too
         "gpu-muelu-chebyshev": Configuration("gpu", "muelu-chebyshev"),
         "gpu-amgx": Configuration("gpu", "amgx"),
         "gpu-amgx-agg": Configuration("gpu", "amgx-agg"),
         "gpu-amgx-single": Configuration("gpu", "amgx-single", tier=2),
         "gpu-cg-jacobi": Configuration("gpu", "cg-jacobi", tier=2),
         "gpu-ginkgo-jacobi": Configuration("gpu", "ginkgo-jacobi", tier=2)})
    if direct_cpu:
        c["mixed-mumps"] = Configuration("mixed", "mumps")
    if native_on_cpu:
        c["mixed-native"] = Configuration("mixed")
    if direct:
        c["gpu-cudss"] = Configuration("gpu", "cudss")
        c["gpu-tacho"] = Configuration("gpu", "tacho", tier=2)
    if nonlinear:
        c["gpu-hypre-amg-newton-ew"] = Configuration("gpu", "hypre-amg",
                                                     "newton-ew")
        c["gpu-cg-jacobi-mf"] = Configuration("gpu", "cg-jacobi",
                                              environment=MATRIX_FREE)
        for k0, tier in (("k0-amg", 1), ("k0-muelu", 1), ("k0-amgx", 2)):
            c["gpu-" + k0] = Configuration("gpu", "hypre-amg", k0,
                                           tier=tier)
            c["gpu-" + k0 + "-mf"] = Configuration("gpu", "hypre-amg", k0,
                                                   environment=MATRIX_FREE,
                                                   tier=tier)
        if direct:
            c["gpu-k0-cudss-mf"] = Configuration("gpu", "hypre-amg",
                                                 "k0-cudss",
                                                 environment=MATRIX_FREE)
            c["gpu-cudss-modified-newton"] = Configuration(
                "gpu", "cudss", "modified-newton", tier=2)
    return c


def _configurations(nonlinear, direct=True, native_on_cpu=True,
                    native_on_gpu=True, direct_cpu=True):
    """the configurations of a case, CPU then GPU"""
    c = _cpu_configurations(nonlinear, direct_cpu)
    c.update(_gpu_configurations(nonlinear, direct, native_on_cpu,
                                 native_on_gpu, direct_cpu))
    return c


def _ex8_configurations():
    """ex8 has a linear solver for each problem (options of the example):
    the global overrides would give the same solver to the thermal and to
    the mechanical problems"""
    muelu = ["-pcMc", "MueLu"]
    c = {"cpu-native": Configuration("cpu"),
         "cpu-minres-muelu": Configuration(
             "cpu", options=["-lsMc", "MINRESSolver"] + muelu,
             processes=[8]),
         "cpu-cg-amg": Configuration(
             "cpu", options=["-lsMc", "HyprePCG", "-pcMc", "HypreBoomerAMG"],
             processes=[8]),
         "cpu-newton-k0": Configuration("cpu", options=["-nlsMc",
                                                        "NewtonK0"],
                                        processes=[8]),
         "mixed-native": Configuration("mixed"),
         "gpu-minres-muelu-chebyshev": Configuration(
             "gpu", options=["-lsMc", "MINRESSolver"] + muelu +
             ["-mls", "chebyshev"]),
         "gpu-cg-muelu-chebyshev": Configuration(
             "gpu", options=["-lsMc", "CGSolver"] + muelu +
             ["-mls", "chebyshev"]),
         "gpu-minres-muelu-jacobi": Configuration(
             "gpu", options=["-lsMc", "MINRESSolver"] + muelu +
             ["-mls", "jacobi", "-mld", "0.3"]),
         "gpu-cudss": Configuration(
             "gpu", options=["-lsMc", "CuDSSSolver", "-mra"]),
         "gpu-hypre-pcg-amg": Configuration(
             "gpu", options=["-lsMc", "HyprePCG", "-pcMc",
                             "HypreBoomerAMG"])}
    return c


# configurations of each case: one solve of a linear problem (ex3, ex5,
# bubble) or non linear computations; the direct solvers are left out of
# the large 3D cases (on the GPU, its 16 GB; on the CPU, MUMPS on ex3 at
# 786,432 unknowns exceeded the 30 GB of the machine)
CONFIGURATIONS = {
    "ex1": _configurations(True, direct=False, native_on_cpu=False),
    "ex2": _configurations(True, native_on_gpu=False),
    "ex3": _configurations(False, direct=False, native_on_cpu=False,
                           direct_cpu=False),
    "ex4": _configurations(True),
    "ex5": _configurations(False, native_on_gpu=False),
    "ex6": _configurations(True, direct=False, direct_cpu=False),
    "ex7": _configurations(True),
    "ex8": _ex8_configurations(),
    "polycrystal": _configurations(True),
    "bubble": _configurations(False, direct=False, direct_cpu=False),
    "cermet": _configurations(True, direct_cpu=False),
}

# numbers of unknowns of the cases whose programs do not print them, given
# by cases/tools/count_unknowns (mesh refined uniformly, vector H1 space of
# the order of the size)
UNKNOWNS = {
    ("ex1", "small"): 375, ("ex1", "screen"): 107811,
    ("ex1", "bench"): 823875,
    ("ex2", "small"): 3726, ("ex2", "screen"): 908226,
    ("ex2", "bench"): 908226,
    ("ex3", "small"): 12288, ("ex3", "screen"): 786432,
    ("ex3", "bench"): 786432,
    ("ex5", "small"): 7442, ("ex5", "screen"): 1847042,
    ("ex5", "bench"): 1847042,
    ("ex6", "small"): 1536, ("ex6", "screen"): 786432,
    ("ex6", "bench"): 786432,
    ("polycrystal", "small"): 20553, ("polycrystal", "screen"): 168885,
    ("polycrystal", "bench"): 168885,
    ("cermet", "small"): 33015, ("cermet", "screen"): 269253,
    ("cermet", "bench"): 269253,
    # ex8 prints the unknowns of the thermal problem first: the ones of the
    # mechanical problem are given
    ("ex8", "small"): 28548, ("ex8", "screen"): 214812,
    ("ex8", "bench"): 214812,
}

FIELDS = ["date", "case", "size", "order", "precision", "configuration",
          "label", "side", "processes", "repeat", "status", "seconds",
          "unknowns", "command"]


def _amgx(case, single, precision="double", aggregation=False):
    """the solver of AmgX of a case (`AmgX` alone is one cycle of its
    multigrid, without convergence criteria), with the classical multigrid
    or the aggregation one; in single precision, the relative tolerance of
    AmgXPCG is at least 1e-6 (floor of AmgX)"""
    base = _tolerance(case, precision)
    tolerance = max(base, 1e-6) if single else base
    precision = ", Precision: \"Single\"" if single else ""
    if aggregation:
        precision += ", Multigrid: \"Aggregation\""
    if case.symmetric:
        return ('{Name: "AmgXPCG", Options: {RelativeTolerance: ' +
                repr(tolerance) + ', MaximumNumberOfIterations: 5000' +
                precision + '}}')
    return ('{Name: "GMRESSolver", Options: {RelativeTolerance: ' +
            repr(case.tolerance) + ', MaximumNumberOfIterations: 5000, '
            'Preconditioner: {Name: "AmgX", Options: {' +
            precision.lstrip(", ") + '}}}}').replace(
                "RelativeTolerance: " + repr(case.tolerance),
                "RelativeTolerance: " + repr(base))


# smallest relative tolerance given to the linear solvers in float
FLOAT_TOLERANCE = 1e-6


def _tolerance(case, precision):
    """the relative tolerance of the linear solvers of a case"""
    if precision == "float":
        return max(case.tolerance, FLOAT_TOLERANCE)
    return case.tolerance


def _format(template, case, precision="double"):
    """replaces the placeholders of a description"""
    s = case.symmetric
    return template.format(
        tolerance=repr(_tolerance(case, precision)),
        hypre="HyprePCG" if s else "HypreGMRES",
        mfem="CGSolver" if s else "GMRESSolver",
        belos="CG" if s else "GMRES",
        petsc="cg" if s else "gmres",
        ginkgo="GinkgoCG" if s else "GinkgoGMRES",
        tangent="CGSolver" if s else "GMRESSolver",
        tacho="chol" if s else "lu",
        spd=", Symmetric: true, PositiveDefinite: true" if s else "",
        sym=", Symmetric: true" if s else "",
        amgx=_amgx(case, False, precision),
        amgx_single=_amgx(case, True, precision),
        amgx_agg=_amgx(case, False, precision, aggregation=True))


def command(case, configuration, processes, size, precision="double"):
    """the command of a run, and its additional environment"""
    c = configuration
    sz = case.sizes[size]
    cmd = ([case.program] + case.options + sz.options +
           ["--order", str(sz.order)] + c.options + SIDES[c.side])
    if c.linear is not None:
        cmd += ["--linear-solver-override",
                _format(LINEAR_SOLVERS[c.linear][0], case, precision)]
    if c.non_linear is not None:
        cmd += ["--non-linear-solver-override",
                _format(NON_LINEAR_SOLVERS[c.non_linear][0], case,
                        precision)]
    # the parallel computations, also on one process (the path of the
    # applications, and of the linear solvers of hypre, Trilinos...)
    if case.parallel_option:
        cmd.append(case.parallel_option)
    if processes > 1:
        # one process per core: physical cores 0 to 7 for 8 processes,
        # all the logical cores (SMT) for 16
        unit = "core" if processes <= 8 else "hwthread"
        smt = ["--use-hwthread-cpus"] if processes > 8 else []
        cmd = ([str(OPENMPI / "bin" / "mpirun"), "-n", str(processes)] + smt +
               ["--bind-to", unit, "--map-by", unit] + cmd)
    environment = dict(c.environment)
    if c.side != "cpu" and (c.linear in _NOT_HYPRE_LINEAR or
                            c.non_linear in _NOT_HYPRE_NON_LINEAR):
        environment.update(SMALL_HYPRE_POOL)
    return cmd, environment


def _session_processes(sid):
    """the processes of a session (the ranks of mpirun are in their own
    process groups, but in the session of mpirun)"""
    pids = []
    for entry in pathlib.Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text()
        except OSError:
            continue
        # the fields after the name of the program, which may hold spaces
        fields = stat[stat.rindex(")") + 2:].split()
        if int(fields[3]) == sid:
            pids.append(int(entry.name))
    return pids


def _stop(p):
    """stops a program started in its own session: SIGTERM to mpirun (which
    stops its ranks), then SIGKILL to all the processes of the session"""
    try:
        os.killpg(p.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        p.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass
    for pid in _session_processes(p.pid):
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    if p.poll() is None:
        p.wait()


def execute(cmd, directory, log, env, timeout):
    """runs a command in its own session; returns the exit status, or the
    string timeout. On a timeout or on any exception (Ctrl-C, SIGTERM of
    the driver), the program is stopped with all its processes"""
    p = subprocess.Popen(cmd, cwd=directory, stdout=log,
                         stderr=subprocess.STDOUT, env=env,
                         start_new_session=True)
    try:
        return p.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        _stop(p)
        return "timeout"
    except BaseException:
        _stop(p)
        raise


# lines of the outputs giving the number of unknowns
_UNKNOWNS = re.compile(
    r"(?:number of (?:finite element )?unknowns|unknowns|dofs|"
    r"degrees of freedom)\s*[:=]?\s*(\d+)", re.IGNORECASE)


def _unknowns(path):
    """the number of unknowns printed by a program, if any"""
    try:
        with open(path, errors="replace") as f:
            for line in f:
                m = _UNKNOWNS.search(line)
                if m:
                    return m.group(1)
    except OSError:
        pass
    return ""


def _check_header(path):
    """refuses a runs.csv written by another version of the driver"""
    if not path.exists():
        return
    with open(path, newline="") as f:
        header = next(csv.reader(f), [])
    if header != FIELDS:
        raise RuntimeError(f"{path}: columns {header}, expected {FIELDS} "
                           "(older driver): use a new results directory")


# lock of the GPU shared with other users of the machine, taken around each
# run using the GPU when given (option --gpu-lock), the time of a run being
# measured once the lock is held
GPU_LOCK = None


@contextlib.contextmanager
def _gpu_lock(side):
    """holds the lock of the GPU during a run using it, if any"""
    if GPU_LOCK is None or side == "cpu":
        yield
        return
    with open(GPU_LOCK, "a") as f:
        # wait for the lock one hour at most
        for _ in range(3600):
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(1)
        else:
            raise RuntimeError(f"the lock {GPU_LOCK} was not obtained in "
                               "one hour")
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def run(results, case_name, configuration_name, precision, processes,
        size, repeat, dry_run=False, timeout=None, first_repeat=0):
    """runs a configuration of a case, appends the results to
    results/runs.csv; returns the statuses and the times of the runs"""
    case = CASES[case_name]
    c = CONFIGURATIONS[case_name][configuration_name]
    outcomes = []
    if precision not in c.precisions():
        print(f"{case_name}, {configuration_name}: not available in "
              f"{precision}")
        return outcomes
    cmd, extra = command(case, c, processes, size, precision)
    path = results / "runs.csv"
    if not dry_run:
        (results / "logs").mkdir(parents=True, exist_ok=True)
        _check_header(path)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env.update(extra)
    order = case.sizes[size].order
    prefix = " ".join(f"{k}={v}" for k, v in extra.items())
    for r in range(first_repeat, first_repeat + repeat):
        name = (f"{case_name}-{size}-o{order}-{precision}-"
                f"{configuration_name}-np{processes}-{r}")
        print((prefix + " " if prefix else "") + shlex.join(cmd))
        if dry_run:
            continue
        log_path = results / "logs" / (name + ".log")
        with _gpu_lock(c.side), open(log_path, "w") as log:
            start = time.perf_counter()
            status = execute(cmd, case.directory(precision), log, env,
                             timeout)
            seconds = time.perf_counter() - start
        row = {"date": datetime.datetime.now().isoformat(timespec="seconds"),
               "case": case_name, "size": size, "order": order,
               "precision": precision, "configuration": configuration_name,
               "label": c.label(), "side": c.side, "processes": processes,
               "repeat": r, "status": status, "seconds": f"{seconds:.3f}",
               "unknowns": (str(UNKNOWNS[(case_name, size)])
                            if (case_name, size) in UNKNOWNS
                            else _unknowns(log_path)),
               "command": (prefix + " " if prefix else "") +
               shlex.join(cmd)}
        new = not path.exists()
        with open(path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            if new:
                w.writeheader()
            w.writerow(row)
        print(f"{name}: status {status}, {seconds:.2f} s")
        outcomes.append((status, seconds))
    return outcomes


def calibrate(results, precision, size, cases, timeout, dry_run):
    """runs each case once with the solvers of the program on 8 processes,
    and on the GPU (hypre with BoomerAMG when the solver of the program is
    refused on the GPU)"""
    for name in cases:
        run(results, name, "cpu-native", precision, 8, size, 1, dry_run,
            timeout)
        gpu = ("gpu-minres-muelu-chebyshev" if name == "ex8" else
               "gpu-hypre-amg" if name in ("ex2", "ex5") else "gpu-native")
        run(results, name, gpu, precision, 1, size, 1, dry_run, timeout)


def _done(results, case_name, precision, size):
    """the configurations and numbers of processes already run (whatever
    their status), to resume an interrupted screening"""
    path = results / "runs.csv"
    if not path.exists():
        return set()
    with open(path, newline="") as f:
        return {(row["configuration"], int(row["processes"]))
                for row in csv.DictReader(f)
                if row["case"] == case_name and row["precision"] == precision
                and row["size"] == size}


def screen(results, case_name, precision, sides, size, timeout, dry_run,
           tiers=(1,)):
    """runs the configurations of a case of the given tiers on the given
    sides, except the ones already run"""
    done = _done(results, case_name, precision, size)
    for name, c in CONFIGURATIONS[case_name].items():
        if c.side not in sides or c.tier not in tiers:
            continue
        for processes in c.processes:
            if (name, processes) in done:
                print(f"{case_name}, {name}, {processes} processes: "
                      "already run")
                continue
            run(results, case_name, name, precision, processes, size, 1,
                dry_run, timeout)


def _screened(results, case_name, precision, size):
    """the best time of each configuration and number of processes in the
    screening of a case (successful runs only)"""
    path = results / "runs.csv"
    if not path.exists():
        return {}
    times = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if (row["case"] != case_name or row["precision"] != precision or
                    row["size"] != size or row["status"] != "0"):
                continue
            k = (row["configuration"], int(row["processes"]))
            times[k] = min(times.get(k, float("inf")), float(row["seconds"]))
    return times


def finals(results, case_name, precision, screen_size, size, repeat,
           timeout, dry_run, long_run=120):
    """runs at the measurement size the configurations selected by the
    screening: the solvers of the program on the CPU (best number of
    processes, or 8 if the screening did not finish), the best other CPU
    configuration, and the two best configurations with the GPU; the
    solvers of the program on the CPU, the best other CPU configuration and
    the best GPU configuration are repeated (`repeat` runs in all), unless
    their first run succeeds in more than `long_run` seconds or fails. The
    configurations of the results which are no longer defined are
    ignored, and the ones already run at this size are not run again"""
    times = _screened(results, case_name, precision, screen_size)
    configurations = CONFIGURATIONS[case_name]

    def ranked(predicate):
        return sorted((t, k) for k, t in times.items()
                      if k[0] in configurations and
                      predicate(configurations[k[0]], k[0]))

    native = ranked(lambda c, n: n == "cpu-native")
    cpu = ranked(lambda c, n: c.side == "cpu" and n != "cpu-native")
    gpu = ranked(lambda c, n: c.side != "cpu")
    selection = [(("cpu-native", native[0][1][1] if native else 8),
                  repeat)]
    if cpu:
        selection.append((cpu[0][1], repeat))
    for i, (_, k) in enumerate(gpu[:2]):
        selection.append((k, repeat if i == 0 else 1))
    # the configurations already measured are not run again, so that finals
    # can be run again after a new screening
    done = _done(results, case_name, precision, size)
    for (name, processes), r in selection:
        if (name, processes) in done:
            continue
        first = run(results, case_name, name, precision, processes, size, 1,
                    dry_run, timeout)
        if (r > 1 and first and first[0][0] == 0 and
                first[0][1] < long_run):
            run(results, case_name, name, precision, processes, size, r - 1,
                dry_run, timeout, first_repeat=1)


def best(results):
    """best time of each case on each side, among the successful runs: the
    CPU with the solvers of the program, the CPU with any configuration,
    and the GPU (at least the behaviours and the finite elements on the
    GPU)"""
    with open(results / "runs.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    table = {}
    for row in rows:
        if row["status"] != "0":
            continue
        key = (row["case"], row["size"], row["order"], row["precision"])
        t = float(row["seconds"])
        kinds = ["gpu"] if row["side"] != "cpu" else (
            ["cpu-native", "cpu"] if row["configuration"] == "cpu-native"
            else ["cpu"])
        for kind in kinds:
            k = key + (kind,)
            if (k not in table) or (t < table[k][0]):
                table[k] = (t, row["configuration"], row["processes"])
    print("case,size,order,precision,cpu_native_seconds,cpu_native_processes,"
          "cpu_best_seconds,cpu_best_configuration,gpu_best_seconds,"
          "gpu_best_configuration,speedup_native,speedup_best")

    def fmt(v):
        return f"{v[0]:.2f}" if v else ""

    for key in sorted({k[:4] for k in table}):
        nat = table.get(key + ("cpu-native",))
        cpu = table.get(key + ("cpu",))
        gpu = table.get(key + ("gpu",))
        print(",".join(list(key) + [
            fmt(nat), nat[2] if nat else "", fmt(cpu), cpu[1] if cpu else "",
            fmt(gpu), gpu[1] if gpu else "",
            f"{nat[0] / gpu[0]:.2f}" if (nat and gpu) else "",
            f"{cpu[0] / gpu[0]:.2f}" if (cpu and gpu) else ""]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="action", required=True)
    li = sub.add_parser("list")
    li.add_argument("case", nargs="?", choices=sorted(CASES))
    r = sub.add_parser("run")
    r.add_argument("results", type=pathlib.Path)
    r.add_argument("case", choices=sorted(CASES))
    r.add_argument("configuration")
    r.add_argument("--np", type=int, default=None)
    r.add_argument("--repeat", type=int, default=1)
    c = sub.add_parser("calibrate")
    c.add_argument("results", type=pathlib.Path)
    c.add_argument("--cases", default=",".join(CASES))
    sc = sub.add_parser("screen")
    sc.add_argument("results", type=pathlib.Path)
    sc.add_argument("case", choices=sorted(CASES))
    sc.add_argument("--sides", default="cpu,gpu,mixed")
    sc.add_argument("--tiers", default="1",
                    help="tiers of the configurations, for instance 1,2")
    fi = sub.add_parser("finals")
    fi.add_argument("results", type=pathlib.Path)
    fi.add_argument("case", choices=sorted(CASES))
    fi.add_argument("--screen-size", default="screen")
    fi.add_argument("--repeat", type=int, default=3)
    for a, size in ((r, "small"), (c, "bench"), (sc, "screen"),
                    (fi, "bench")):
        a.add_argument("--precision", default="double", choices=sorted(BUILDS))
        a.add_argument("--size", default=size)
        a.add_argument("--timeout", type=float, default=None,
                       help="seconds after which a run is stopped")
        a.add_argument("--dry-run", action="store_true")
        a.add_argument("--gpu-lock", default=None,
                       help="file locked around each run using the GPU")
    b = sub.add_parser("best")
    b.add_argument("results", type=pathlib.Path)
    args = parser.parse_args()
    # SIGTERM stops the run in progress, as Ctrl-C (SIGKILL can't be
    # handled: the program would then go on)
    signal.signal(signal.SIGTERM, lambda *a: sys.exit(143))
    global GPU_LOCK
    GPU_LOCK = getattr(args, "gpu_lock", None)
    if args.action == "list":
        for name, case in CASES.items():
            if args.case and name != args.case:
                continue
            sizes = ", ".join(f"{k} (order {v.order}: {' '.join(v.options)})"
                              for k, v in case.sizes.items())
            print(f"{name}: {case.program}, sizes {sizes}")
            for cn, cf in CONFIGURATIONS[name].items():
                print(f"  {cn}: {cf.side}, {cf.label()}, processes "
                      f"{cf.processes}, {'/'.join(sorted(cf.precisions()))}, "
                      f"tier {cf.tier}")
    elif args.action == "run":
        cf = CONFIGURATIONS[args.case].get(args.configuration)
        if cf is None:
            parser.error(f"unknown configuration '{args.configuration}' of "
                         f"{args.case}")
        if args.size not in CASES[args.case].sizes:
            parser.error(f"unknown size '{args.size}'")
        for n in ([args.np] if args.np else cf.processes):
            run(args.results, args.case, args.configuration, args.precision,
                n, args.size, args.repeat, args.dry_run, args.timeout)
    elif args.action == "calibrate":
        cases = args.cases.split(",")
        unknown = [x for x in cases if x not in CASES]
        if unknown:
            parser.error("unknown cases: " + ", ".join(unknown))
        calibrate(args.results, args.precision, args.size, cases,
                  args.timeout, args.dry_run)
    elif args.action == "screen":
        sides = args.sides.split(",")
        unknown = [x for x in sides if x not in SIDES]
        if unknown:
            parser.error("unknown sides: " + ", ".join(unknown))
        screen(args.results, args.case, args.precision, sides, args.size,
               args.timeout, args.dry_run,
               tuple(int(t) for t in args.tiers.split(",")))
    elif args.action == "finals":
        finals(args.results, args.case, args.precision, args.screen_size,
               args.size, args.repeat, args.timeout, args.dry_run)
    else:
        best(args.results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
