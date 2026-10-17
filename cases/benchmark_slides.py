#!/usr/bin/env python3
"""Benchmark of the examples of mfem-mgis-examples and of mm-opera-hpc for
the slides of the prototype.

For each case, order of the elements (1 or 2, with a size of each order),
size and precision (double or float), the wall time of the whole program
is measured on the CPU (MPI on 8 or 16 processes, behaviours
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
ex7 a coarser mesh or fewer steps), a configuration being stopped beyond
three times the best time of its side; `finals` runs at the measurement
size the configurations selected by the screening (the next one when one
fails, and the best GPU configuration with several processes sharing the
GPU); `run` runs one configuration; `best` and `cases/benchmark_plots.py`
summarize the results.

usage:
  benchmark_slides.py list [CASE]
  benchmark_slides.py run RESULTS CASE CONFIGURATION [--precision P]
                      [--order O] [--np N] [--size SIZE] [--repeat R]
                      [--timeout S] [--dry-run]
  benchmark_slides.py calibrate RESULTS [--precision P] [--order O]
                      [--size SIZE] [--cases LIST] [--timeout S] [--dry-run]
  benchmark_slides.py screen RESULTS CASE [--precision P] [--order O]
                      [--sides LIST] [--size SIZE] [--timeout S]
                      [--reference RESULTS --keep K] [--dry-run]
  benchmark_slides.py finals RESULTS CASE [--precision P] [--order O]
                      [--size SIZE] [--screen-size SIZE] [--repeat R]
                      [--timeout S] [--dry-run]

In dry-run mode, `finals` prints the runs of all the candidates (the runs
which would follow a failure included), but not the runs with several
processes sharing the GPU (they depend on the measured times)
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
import shutil
import signal
import subprocess
import sys
import threading
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
# builds of the examples and of mm-opera-hpc for each precision
BUILDS = {"double": ROOT / "build" / "bench",
          "float": ROOT / "build" / "bench-float"}
# libraries of MGIS (displacement gradient interface) of each build: the
# one of the precision of a run comes first in LD_LIBRARY_PATH, the other
# one is removed
MGIS_LIBRARIES = {"double": ROOT / "install" / "h-interface" / "mgis" / "lib",
                  "float": ROOT / "install" / "h-interface" / "mgis-float" /
                  "lib"}
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
                 tolerance, symmetric, parallel_option=None,
                 precision_options=None):
        # "examples" or "opera"; the directory of the program, where it
        # runs (meshes and references), is given by the precision
        self.project = project
        self.subdirectory = subdirectory
        # program, relative to the directory
        self.program = program
        # options of all the runs
        self.options = options
        # sizes (class Size) of each order (1 and 2): "small" checking the
        # scripts, "screen" for the screening and "bench" for the
        # measurements
        self.sizes = sizes
        # relative tolerance of the linear solver of the program, given to
        # the linear solvers of the overrides
        self.tolerance = tolerance
        # if true, the tangent operator is symmetric (conjugate gradient)
        self.symmetric = symmetric
        # option selecting the parallel computations, if the program
        # computes in sequential by default
        self.parallel_option = parallel_option
        # options of the runs of a precision
        self.precision_options = precision_options or {}

    def directory(self, precision):
        """the directory of the program for a precision"""
        return BUILDS[precision] / self.project / self.subdirectory


def _sizes(small, screen, bench, order=1):
    """the sizes of a case at an order"""
    return {"small": Size(order, small), "screen": Size(order, screen),
            "bench": Size(order, bench)}


def _orders(first, second):
    """the sizes of a case at the orders 1 and 2: `first` and `second` are
    the triples of options (small, screen, bench) of each order"""
    return {1: _sizes(*first, order=1), 2: _sizes(*second, order=2)}


# The ParaView exports are disabled (they would be timed with the
# computations); the quantity of interest of each case is written in any
# case (see QUANTITIES); the profiler of each program is enabled (its time
# table gives the phases of the computations, see PHASES; on the GPU, the
# sections synchronize the device). The sizes are given for the orders 1 and 2. The sizes "bench" fill the 8
# cores and the GPU (about one million unknowns, one or several waves of
# the multiprocessors of the GPU for the kernels of the elements); at the
# order 2, the mesh has one refinement less, for about the same number of
# unknowns as at the order 1, except for ex8 (direct solver of the program,
# MUMPS, whose memory limits the size). "screen" keeps the mesh and the
# time step of "bench", with one or two steps when the program has an
# option of the end time; ex1 and ex7 have no such option: "screen" uses a
# coarser mesh (ex1) or fewer steps (ex7).
CASES = {
    # the time step of ex1 is 1 / nbsteps (100 steps by default); its
    # loading is cyclic, reversed at t = 0.3. Without --prediction, the
    # increment of the imposed displacement is concentrated in the elements
    # next to the loaded face: with 4 refinements or more, the Newton solver
    # does not converge at the reversal
    "ex1": Case("examples", "ex1", "./uniaxial_tensile_test",
                ["--no-post-processing", "--prediction", "--profiling"],
                _orders((["-r", "2", "-ns", "25"], ["-r", "5", "-ns", "25"],
                         ["-r", "6", "-ns", "25"]),
                        (["-r", "1", "-ns", "25"], ["-r", "4", "-ns", "25"],
                         ["-r", "5", "-ns", "25"])),
                1e-12, True, parallel_option="--parallel"),
    "ex2": Case("examples", "ex2", "./ssna303",
                ["--no-paraview", "--profiling"],
                _orders((["--nbsteps", "5"],
                         ["-r", "4", "--nbsteps", "2", "--end-time", "0.04"],
                         ["-r", "4", "--nbsteps", "50"]),
                        (["--nbsteps", "5"],
                         ["-r", "3", "--nbsteps", "2", "--end-time", "0.04"],
                         ["-r", "3", "--nbsteps", "50"])),
                1e-12, False, parallel_option="--parallel"),
    # the threshold of the L2 error of ex3 to the analytical solution is
    # calibrated in float on the coarse mesh of the test (4e-5 with the
    # conjugate gradient): with 786432 unknowns, the error is 7.9e-4 in
    # float (2.6e-12 in double), the relative tolerance of the iterative
    # solvers being bounded to 1.2e-5 by mfem-mgis in float. The check is
    # disabled in float.
    "ex3": Case("examples", "ex3", "./two_layer_cube",
                ["--test-case", "1", "--no-post-processing", "--profiling"],
                _orders((["-r", "2"], ["-r", "4"], ["-r", "4"]),
                        (["-r", "1"], ["-r", "3"], ["-r", "3"])),
                1e-12, True, parallel_option="--parallel",
                precision_options={"float": ["--no-check"]}),
    "ex4": Case("examples", "ex4", "./ssna303_3d_hypre",
                ["--no-paraview", "--profiling"],
                _orders((["--nbsteps", "2", "--end-time", "1"],
                         ["-r", "1", "--nbsteps", "2", "--end-time", "0.2"],
                         ["-r", "1", "--nbsteps", "10", "--end-time", "1"]),
                        (["--nbsteps", "2", "--end-time", "1"],
                         ["--nbsteps", "2", "--end-time", "0.2"],
                         ["--nbsteps", "10", "--end-time", "1"])),
                1e-6, False),
    "ex5": Case("examples", "ex5", "./satoh",
                ["--no-post-processing", "--profiling"],
                _orders(([], ["-r", "5"], ["-r", "5"]),
                        ([], ["-r", "4"], ["-r", "4"])),
                1e-12, True, parallel_option="--parallel"),
    # the same for ex6: its check of the analytical solution tolerates
    # 4e-5 in float (error of 3.9e-6 on the coarse mesh), the error being
    # 9.9e-4 with 786432 unknowns
    "ex6": Case("examples", "ex6", "./rve", ["--no-post-processing"],
                _orders((["-r", "1"], ["-r", "4"], ["-r", "4"]),
                        (["-r", "1"], ["-r", "3"], ["-r", "3"])),
                1e-12, True, precision_options={"float": ["--no-check"]}),
    # ex7 has no option of the end time: 40 steps by default; "screen"
    # has 5 steps, "bench" 10
    "ex7": Case("examples", "ex7", "./mox2",
                ["--no-post-processing", "--mesh", "mesh/OneSphere.msh"],
                _orders((["--nbsteps", "5"], ["-r", "3", "--nbsteps", "5"],
                         ["-r", "3", "--nbsteps", "10"]),
                        (["--nbsteps", "5"], ["-r", "3", "--nbsteps", "5"],
                         ["-r", "3", "--nbsteps", "10"])),
                1e-12, False),
    # the mechanical problem of ex8 is solved by MUMPS: 214812 unknowns at
    # both orders; "screen" is the first of the 5 steps of the measurements
    "ex8": Case("examples", "ex8", "./rjh_plate",
                ["-m", "assemblage_hexa.msh", "--no-post-processing",
                 "--debug"],
                _orders((["-ns", "1"],
                         ["-r", "1", "-et", "2e4", "-tr", "1e5", "-ns", "1"],
                         ["-r", "1", "-ns", "5"]),
                        (["-ns", "1"],
                         ["-r", "0", "-et", "2e4", "-tr", "1e5", "-ns", "1"],
                         ["-r", "0", "-ns", "5"])),
                1e-12, False),
    # the polycrystal and the cermet run for short durations, the same in
    # double and in float: the grain law of UO2 fails in float beyond about
    # 0.8 s (rounding of its elastic strain, results/precision-polycrystal;
    # the solvers of the program failed in float with 4 steps up to 0.8 s
    # in the first pass): 3 steps up to 0.6 s
    "polycrystal": Case("opera", "polycrystal", "./uniaxial-polycrystal",
                        ["--disable-post-processing"],
                        _orders((["--nstep", "2", "--duration", "0.4"],
                                 ["-r", "2", "--nstep", "1", "--duration",
                                  "0.2"],
                                 ["-r", "2", "--nstep", "3", "--duration",
                                  "0.6"]),
                                (["--nstep", "2", "--duration", "0.4"],
                                 ["-r", "1", "--nstep", "1", "--duration",
                                  "0.2"],
                                 ["-r", "1", "--nstep", "3", "--duration",
                                  "0.6"])),
                        1e-12, False),
    # the program solves with a relative tolerance of 1e-10
    "bubble": Case("opera", "bubble", "./test-bubble",
                   ["--post-processing", "2"],
                   _orders(([], ["-r", "1"], ["-r", "1"]), ([], [], [])),
                   1e-10, True),
    # finite strain behaviours: tangent operator not symmetric in general
    # (same grain law as the polycrystal)
    "cermet": Case("opera", "cermet", "./cermet", ["--postprocessing", "0"],
                   _orders((["--nstep", "1", "--duration", "0.4"],
                            ["-r", "1", "--nstep", "1", "--duration", "0.2"],
                            ["-r", "1", "--nstep", "2", "--duration",
                             "0.4"]),
                           (["--nstep", "1", "--duration", "0.4"],
                            ["--nstep", "1", "--duration", "0.2"],
                            ["--nstep", "2", "--duration", "0.4"])),
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
    # (trials of the options of BoomerAMG on the GPU at the screening
    # sizes: the aggressive coarsening alone is the best one, 36 % faster
    # than the default options on ex2, as fast elsewhere, and uses less
    # memory; the Chebyshev smoother, the l1-Jacobi one, the strength
    # threshold 0.5, the interpolation 14, at most 6 elements per
    # row, the aggressive interpolations 5 or 7 and the transpose kept do
    # not improve the default options. The aggressive interpolation of type 6
    # crashes in hypre 3.2.0 on the GPU, and the options of elasticity, the
    # rigid body modes, are refused by hypre on the GPU)
    # Krylov solvers of MFEM preconditioned by BoomerAMG built every 5
    # iterations of the non linear solver on the assembled tangent operator
    # (only with the tangent operator applied without assembly, which
    # provides the assembled operator)
    "mfem-amg-lagged": (_krylov("{mfem}", '{{Name: "AssembledOperator", '
                                'Options: {{UpdatePeriod: 5, Preconditioner: '
                                '{{Name: "HypreBoomerAMG", Options: '
                                '{{VerbosityLevel: 0}}}}}}}}'), _BOTH),
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
    "amgx-agg8": ("{amgx_agg8}", _DOUBLE),
    "amgx-gs": ("{amgx_gs}", _DOUBLE),
    # direct solvers, the analysis being reused while the pattern of the
    # matrix does not change
    "cudss": ('{{Name: "CuDSSSolver", Options: {{ReuseAnalysis: true{spd}}}}}',
              _BOTH),
    "mumps": ('{{Name: "MUMPSSolver", Options: {{ReuseAnalysis: true{sym}}}}}',
              _BOTH),
    # Tacho: Cholesky factorization by default, LU for the non symmetric
    # cases (Tacho of Trilinos 17.2.1 is only instantiated in double); the
    # parameters of the solver are in the sublist named after it
    "tacho": ('{{Name: "Amesos2Solver", Options: {{Solver: "Tacho", '
              'Parameters: {{Tacho: {{method: "{tacho}"}}}}}}}}', _DOUBLE),
    "strumpack": ('{{Name: "STRUMPACKSolver", Options: {{ReuseAnalysis: '
                  'true}}}}', _DOUBLE),
    # PETSc with its algebraic multigrid, stopped on the true residual
    # PETSc computes on the host only (mfem-mgis copies the operator to an
    # AIJ matrix of the host: MFEM 4.10 with PETSc 3.26 and CUDA sees null
    # vectors); smoothed aggregation with Chebyshev and Jacobi smoothers.
    # The rigid body modes are only given with the ordering byVDIM (the
    # programs use byNODES)
    "petsc-gamg": ('{{Name: "PETScLinearSolver", Options: {{Options: '
                   '{{ksp_type: "{petsc}", pc_type: "gamg", '
                   'pc_gamg_threshold: 0.02, pc_gamg_agg_nsmooths: 1, '
                   'mg_levels_ksp_type: "chebyshev", mg_levels_pc_type: '
                   '"jacobi", ksp_rtol: {tolerance}, ksp_max_it: 5000, '
                   'ksp_norm_type: "unpreconditioned"}}}}}}', _BOTH),
    # Ginkgo: Krylov solvers with a Jacobi preconditioner on the GPU (the
    # sparse solvers of MAGMA are only available in sequential, and all the
    # programs compute in parallel)
    # (Ginkgo is refused by MFEM in single precision). The scalar Jacobi
    # preconditioner needs fewer iterations than the blocks of Ginkgo (up
    # to 32 unknowns, by default): 396 iterations instead of 652 on ex3
    "ginkgo-jacobi": (_krylov("{ginkgo}", '{{Name: "GinkgoJacobi", Options: '
                              '{{Configuration: {{max_block_size: 1}}}}}}',
                              iterations=100000), _DOUBLE),
    # Ginkgo with its algebraic multigrid (parallel graph match, the only
    # preconditioner of Ginkgo which is not local to each process): the
    # default smoother (one Jacobi iteration relaxed by 0.9) makes the
    # multigrid almost diverge on the elasticity (113 s instead of 3 s on
    # ex3 at the small size), hence two iterations relaxed by 0.5. The
    # incomplete factorizations (IC, ILU, with ISAI) break down on ex3 (not a
    # number from the first iteration): they are not screened
    "ginkgo-amg": (_krylov("{ginkgo}", '{{Name: "GinkgoMultigrid", Options: '
                           '{{Configuration: {{mg_level: {{type: '
                           '"multigrid::Pgm", deterministic: true}}, '
                           'smoother_relax: 0.5, smoother_iters: 2}}}}}}'),
                   _DOUBLE),
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
    # factorization reused by the iterations with a direct solver, the
    # hierarchy of an AMG built once per resolution), with or without the
    # Anderson acceleration (NewtonK0 can't be accelerated)
    "modified-newton": ('{{Name: "FixedPoint", Options: {{Method: '
                        '"ModifiedNewton", MaximumNumberOfIterations: 50}}}}',
                        _BOTH),
    "modified-newton-anderson": ('{{Name: "FixedPoint", Options: {{Method: '
                                 '"ModifiedNewton", '
                                 'MaximumNumberOfIterations: 50, '
                                 'Acceleration: {{Name: "Anderson", '
                                 'Options: {{Depth: 5, Regularization: '
                                 '1e-10}}}}}}}}', _BOTH),
    # fixed point on the elastic operator (computed once by integrating the
    # behaviours, inverted by one cycle of BoomerAMG) accelerated by
    # Anderson
    "k0-anderson": ('{{Name: "FixedPoint", Options: {{Method: '
                    '"ElasticOperator", ElasticOperatorComputation: '
                    '"Integration", ElasticOperatorUpdate: "Once", '
                    'ElasticOperatorLinearSolver: "HypreBoomerAMG", '
                    'MaximumNumberOfIterations: 200, Acceleration: {{Name: '
                    '"Anderson", Options: {{Depth: 8}}}}}}}}', _BOTH),
}

# tangent operator applied without assembly (kept with --execution)
MATRIX_FREE = {"MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION": "MatrixFree"}

# With the GPU, hypre allocates at once a pool of 4 GiB (Umpire), which
# leaves too little memory to Trilinos, AmgX or cuDSS for the largest
# cases (failed allocations on ex3, ex4, ex6): with their solvers, the
# initial pool of hypre is 1 GiB (it grows if needed)
SMALL_HYPRE_POOL = {"MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE": "1"}
# in the mixed placement, hypre computes on the host: the initial pool of
# hypre on the GPU, if any, is reduced since the GPU is shared by all the
# processes
MIXED_HYPRE_POOL = {"MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE": "0.25"}
_NOT_HYPRE_LINEAR = {"muelu-jacobi", "muelu-chebyshev", "muelu-sgs", "amgx",
                     "amgx-single", "amgx-agg", "amgx-agg8", "amgx-gs",
                     "cudss", "tacho", "ginkgo-jacobi", "ginkgo-amg"}
_NOT_HYPRE_NON_LINEAR = {"k0-muelu", "k0-muelu-sgs", "k0-amgx", "k0-cudss"}

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
                 environment=None, processes=None, tier=1, coupled=False):
        self.side = side
        # on the CPU, without the option --execution: the behaviours are
        # integrated point by point and the finite element computations are
        # made element by element, the integration and the assembly being
        # coupled in the loops over the elements (the organization of the
        # master branch of mfem-mgis), instead of the integration by
        # batches and the kernels of `--execution cpu`
        self.coupled = coupled
        # 1: screened first; 2: screened in a second pass (all the
        # configurations are of tier 1 since the screening stops the slow
        # ones, see `_adaptive_timeout`)
        self.tier = tier
        self.linear = linear
        self.non_linear = non_linear
        self.options = list(options)
        self.environment = environment or {}
        # numbers of processes: 8 and 16 (SMT) on the CPU, 1 on the GPU (the
        # best GPU configuration is also run with several processes sharing
        # the GPU by `finals`), 8 in the mixed placement (linear solver on
        # the CPU)
        self.processes = processes or {"cpu": [8, 16],
                                       "mixed": [8]}.get(side, [1])

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
        if self.coupled:
            parts.append("coupled")
        if self.options:
            parts.append(" ".join(self.options))
        return ", ".join(parts)


def _cpu_configurations(nonlinear, direct_cpu=True):
    """the solvers of the program on 8 and 16 processes, strong alternatives
    on 8 processes, and the non linear solvers for the non linear
    computations; `direct_cpu`: if false, MUMPS is left out (its
    factorization of a large 3D case exceeds the memory of the machine)"""
    c = {"cpu-native": Configuration("cpu"),
         "cpu-native-coupled": Configuration("cpu", coupled=True),
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
        c["cpu-hypre-amg-modified-newton"] = Configuration(
            "cpu", "hypre-amg", "modified-newton", processes=[8])
        c["cpu-k0-anderson"] = Configuration("cpu", "hypre-amg",
                                             "k0-anderson", processes=[8])
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
         "gpu-hypre-amg-agg": Configuration("gpu", "hypre-amg-agg"),
         "gpu-mfem-amg": Configuration("gpu", "mfem-amg"),
         "gpu-muelu-jacobi": Configuration("gpu", "muelu-jacobi"),
         # the conjugate gradient of Belos stops on a non positive p.A.p
         # on ex1 with the Jacobi smoother (GPU) and the symmetric
         # Gauss-Seidel one (CPU), and on the bubble with the Jacobi one:
         # the Chebyshev smoother is screened too
         "gpu-muelu-chebyshev": Configuration("gpu", "muelu-chebyshev"),
         "gpu-amgx": Configuration("gpu", "amgx"),
         "gpu-amgx-agg": Configuration("gpu", "amgx-agg"),
         "gpu-amgx-single": Configuration("gpu", "amgx-single"),
         "gpu-cg-jacobi": Configuration("gpu", "cg-jacobi"),
         "gpu-ginkgo-jacobi": Configuration("gpu", "ginkgo-jacobi"),
         "gpu-ginkgo-amg": Configuration("gpu", "ginkgo-amg"),
         "gpu-amgx-agg8": Configuration("gpu", "amgx-agg8"),
         "gpu-amgx-gs": Configuration("gpu", "amgx-gs")})
    # mixed placement: the elements and the behaviours on the GPU, shared by
    # the 8 processes, the linear solvers on the 8 cores
    if direct_cpu:
        c["mixed-mumps"] = Configuration("mixed", "mumps")
    if native_on_cpu:
        c["mixed-native"] = Configuration("mixed")
    c["mixed-hypre-amg"] = Configuration("mixed", "hypre-amg")
    c["mixed-petsc-gamg"] = Configuration("mixed", "petsc-gamg")
    if nonlinear:
        # NewtonK0 is refused in the mixed placement: its tangent solver (a
        # Krylov solver of MFEM) computes where MFEM computes, on the GPU
        c["mixed-hypre-amg-newton-ew"] = Configuration("mixed", "hypre-amg",
                                                       "newton-ew")
    if direct:
        c["gpu-cudss"] = Configuration("gpu", "cudss")
        c["gpu-tacho"] = Configuration("gpu", "tacho")
    if nonlinear:
        c["gpu-hypre-amg-newton-ew"] = Configuration("gpu", "hypre-amg",
                                                     "newton-ew")
        # the hierarchy of BoomerAMG built once per resolution (modified
        # Newton method), every 5 iterations (assembled operator of the
        # tangent operator applied without assembly), or once for the whole
        # computation (fixed point on the elastic operator, Anderson)
        c["gpu-hypre-amg-modified-newton"] = Configuration(
            "gpu", "hypre-amg", "modified-newton")
        c["gpu-hypre-amg-modified-newton-anderson"] = Configuration(
            "gpu", "hypre-amg", "modified-newton-anderson")
        c["gpu-amg-lagged-mf"] = Configuration("gpu", "mfem-amg-lagged",
                                               environment=MATRIX_FREE)
        c["gpu-k0-anderson"] = Configuration("gpu", "hypre-amg",
                                             "k0-anderson")
        c["gpu-cg-jacobi-mf"] = Configuration("gpu", "cg-jacobi",
                                              environment=MATRIX_FREE)
        for k0, tier in (("k0-amg", 1), ("k0-muelu", 1), ("k0-amgx", 1)):
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
                "gpu", "cudss", "modified-newton")
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
         "cpu-native-coupled": Configuration("cpu", coupled=True),
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
             ["-mls", "chebyshev"],
             environment=SMALL_HYPRE_POOL),
         "gpu-cg-muelu-chebyshev": Configuration(
             "gpu", options=["-lsMc", "CGSolver"] + muelu +
             ["-mls", "chebyshev"],
             environment=SMALL_HYPRE_POOL),
         "gpu-minres-muelu-jacobi": Configuration(
             "gpu", options=["-lsMc", "MINRESSolver"] + muelu +
             ["-mls", "jacobi", "-mld", "0.3"],
             environment=SMALL_HYPRE_POOL),
         "gpu-cudss": Configuration(
             "gpu", options=["-lsMc", "CuDSSSolver", "-mra"],
             environment=SMALL_HYPRE_POOL),
         "gpu-hypre-pcg-amg": Configuration(
             "gpu", options=["-lsMc", "HyprePCG", "-pcMc",
                             "HypreBoomerAMG"])}
    return c


# configurations of each case: one solve of a linear problem (ex3, ex5,
# bubble) or non linear computations; the direct solvers are left out of
# the large 3D cases (on the GPU, its 16 GB; on the CPU, MUMPS on ex3 at
# 786,432 unknowns exceeded the 30 GB of the machine), except when they are
# the solvers of the program (ex2, ex5, ex8)
CONFIGURATIONS = {
    "ex1": _configurations(True, direct=False, native_on_cpu=False,
                           direct_cpu=False),
    "ex2": _configurations(True, native_on_gpu=False),
    "ex3": _configurations(False, direct=False, native_on_cpu=False,
                           direct_cpu=False),
    "ex4": _configurations(True, direct=False, direct_cpu=False),
    "ex5": _configurations(False, native_on_gpu=False),
    "ex6": _configurations(True, direct=False, direct_cpu=False),
    "ex7": _configurations(True, direct=False, direct_cpu=False),
    "ex8": _ex8_configurations(),
    "polycrystal": _configurations(True, direct=False, direct_cpu=False),
    "bubble": _configurations(False, direct=False, direct_cpu=False),
    "cermet": _configurations(True, direct=False, direct_cpu=False),
}
# memory of the GPU in the mixed placement (16 GiB): each process sharing
# the GPU reserves the memory of the integration of its behaviours, and 6 or
# 8 processes exhaust it for the crystal plasticity behaviour (polycrystal,
# cermet) at both orders and for ex6 at the order 1 (the order 2 of ex6
# fits on 8 processes, the run on 8 processes fails at once at the order 1)
_MIXED_PROCESSES = {"polycrystal": [4], "cermet": [4], "ex6": [8, 4]}
for _case, _processes in _MIXED_PROCESSES.items():
    for _c in CONFIGURATIONS[_case].values():
        if _c.side == "mixed":
            _c.processes = list(_processes)

# numbers of unknowns of the cases whose programs do not print them, given
# by cases/tools/count_unknowns (mesh refined uniformly, vector H1 space of
# the order of the size)
UNKNOWNS = {
    ("ex1", 1, "small"): 375, ("ex1", 1, "screen"): 107811,
    ("ex1", 1, "bench"): 823875,
    ("ex1", 2, "small"): 375, ("ex1", 2, "screen"): 107811,
    ("ex1", 2, "bench"): 823875,
    ("ex2", 1, "small"): 3726, ("ex2", 1, "screen"): 908226,
    ("ex2", 1, "bench"): 908226,
    ("ex2", 2, "small"): 14522, ("ex2", 2, "screen"): 908226,
    ("ex2", 2, "bench"): 908226,
    ("ex3", 1, "small"): 12288, ("ex3", 1, "screen"): 786432,
    ("ex3", 1, "bench"): 786432,
    ("ex3", 2, "small"): 12288, ("ex3", 2, "screen"): 786432,
    ("ex3", 2, "bench"): 786432,
    ("ex5", 1, "screen"): 1847042, ("ex5", 1, "bench"): 1847042,
    ("ex5", 2, "small"): 7442, ("ex5", 2, "screen"): 1847042,
    ("ex5", 2, "bench"): 1847042,
    ("ex6", 1, "small"): 1536, ("ex6", 1, "screen"): 786432,
    ("ex6", 1, "bench"): 786432,
    ("ex6", 2, "small"): 12288, ("ex6", 2, "screen"): 786432,
    ("ex6", 2, "bench"): 786432,
    ("polycrystal", 1, "small"): 20553,
    ("polycrystal", 1, "screen"): 1360002,
    ("polycrystal", 1, "bench"): 1360002,
    ("polycrystal", 2, "small"): 168885,
    ("polycrystal", 2, "screen"): 1360002,
    ("polycrystal", 2, "bench"): 1360002,
    ("cermet", 1, "small"): 33015, ("cermet", 1, "screen"): 269253,
    ("cermet", 1, "bench"): 269253,
    ("cermet", 2, "small"): 269253, ("cermet", 2, "screen"): 269253,
    ("cermet", 2, "bench"): 269253,
    # ex8 prints the unknowns of the thermal problem first: the ones of the
    # mechanical problem are given
    ("ex8", 1, "small"): 28548, ("ex8", 1, "screen"): 214812,
    ("ex8", 1, "bench"): 214812,
    ("ex8", 2, "small"): 214812, ("ex8", 2, "screen"): 214812,
    ("ex8", 2, "bench"): 214812,
}

# Quantities of interest of the cases, compared to the ones of the solvers
# of the program on the CPU to check the results of each configuration:
# read in the output files of a run (moved to outputs/<run>/ after it), or
# in its log. A quantity is a list of values: the last row of a file
# (without its first column, the time, if `skip` is true), a column of a
# CSV file, or the values found by a regular expression in the log. Its
# kind is "value" (compared by the relative difference of the norms of the
# differences) or "error" (an error to an analytical solution, which must
# not exceed ten times the one of the reference, or an absolute floor).
class Quantity:
    """the quantity of interest of a case"""

    def __init__(self, kind, files=(), log=None, skip=True, column=None,
                 groups=None, offsets=None):
        self.kind = kind
        # sizes of the groups of values compared together (the components
        # of a tensor, the minimum, maximum and mean values of a field): the
        # difference of each value is relative to the norm of its group in
        # the reference. By default, each value is compared alone (values
        # of different scales, as a displacement and a force)
        self.groups = groups
        # values subtracted before the comparison (the identity from the
        # components of a deformation gradient, close to 1, whose
        # differences are the ones of the strains)
        self.offsets = offsets
        # output files of the program (names or patterns in its directory)
        self.files = list(files)
        # function giving the values from the text of the log
        self.log = log
        self.skip = skip
        self.column = column

    def read(self, log_path, outputs):
        """the values of the quantity of a run, or None"""
        try:
            if self.log is not None:
                return self.log(log_path.read_text(errors="replace"))
            for pattern in self.files:
                paths = sorted(outputs.glob(pattern))
                if not paths:
                    continue
                lines = [line for line in paths[-1].read_text().splitlines()
                         if line.strip() and not line.startswith("#")]
                if not lines:
                    return None
                if self.column is not None:
                    header = lines[0].rstrip(",").split(",")
                    row = lines[-1].rstrip(",").split(",")
                    return [float(row[header.index(self.column)])]
                values = [float(v) for v in lines[-1].split()]
                return values[1:] if self.skip else values
        except (OSError, ValueError, IndexError):
            return None
        return None


def _last_number(label):
    """a reader of the last number printed after `label` in a log"""
    expression = re.compile(re.escape(label) +
                            r"\s*([-+]?(?:nan|inf)|[-+0-9.eE]+)",
                            re.IGNORECASE)

    def read(text):
        found = expression.findall(text)
        return [float(found[-1])] if found else None
    return read


def _statistics(text):
    """the minimum, maximum and mean values of the fields printed by ex8
    (option --debug) after the line "# field, minimum, maximum and mean
    values" """
    header = "# field, minimum, maximum and mean values"
    if header not in text:
        return None
    values = []
    for line in text.split(header)[-1].splitlines()[1:]:
        parts = line.split()
        if len(parts) != 4:
            break
        values += [float(v) for v in parts[1:]]
    return values or None


QUANTITIES = {
    # uniaxial tension: gradients, stress and plastic strain at the first
    # integration point at the last step (no time column)
    "ex1": Quantity("value", ["UniaxialTensileTest-*.txt"], skip=False),
    # mean stress of the material at the last step
    "ex2": Quantity("value", ["force.txt"], groups=[2]),
    "ex3": Quantity("error", log=_last_number("L2Error:")),
    # resultant force on the upper boundary at the last step (the
    # horizontal components are the noise of the resolutions)
    "ex4": Quantity("value", ["force_*.txt"], groups=[3]),
    "ex5": Quantity("value",
                    log=_last_number("reaction of the left boundary:")),
    "ex6": Quantity("error",
                    log=_last_number("error to the analytical solution:")),
    # mean stresses of the materials at the last step
    "ex7": Quantity("value", ["avgStress"], groups=[6, 6]),
    # minimum, maximum and mean values of the fields (option --debug)
    "ex8": Quantity("value", log=_statistics, groups=[3] * 4),
    # diagonal components of the macroscopic gradient and of the
    # macroscopic stress at the last step
    "polycrystal": Quantity("value", ["uniaxial-polycrystal.res"],
                            groups=[3, 3], offsets=[1, 1, 1, 0, 0, 0]),
    # average hydrostatic stress of the cell
    "bubble": Quantity("value", ["sig_hydro.txt"], column="Average"),
    "cermet": Quantity("value", ["cermet.res"], groups=[3, 3],
                       offsets=[1, 1, 1, 0, 0, 0]),
}

# relative tolerances of the comparison of the quantities of interest to the
# ones of the solvers of the program on the CPU, and absolute floors of the
# errors to the analytical solutions
QUANTITY_TOLERANCES = {"double": 1e-5, "float": 1e-3}
ERROR_FLOORS = {"double": 1e-8, "float": 1e-2}


def _norm(v):
    return sum(x * x for x in v) ** 0.5


def _relative_difference(case_name, values, reference):
    """the largest difference of the values to the reference, relative to
    the norm of their group in the reference (see `Quantity.groups`)"""
    q = QUANTITIES[case_name]
    if q.offsets is not None and len(q.offsets) == len(reference):
        values = [v - o for v, o in zip(values, q.offsets)]
        reference = [v - o for v, o in zip(reference, q.offsets)]
    groups = q.groups or [1] * len(reference)
    if sum(groups) != len(reference):
        groups = [len(reference)]
    worst = 0.
    start = 0
    for n in groups:
        a = values[start:start + n]
        r = reference[start:start + n]
        start += n
        difference = _norm([x - y for x, y in zip(a, r)])
        scale = _norm(r)
        if difference == 0:
            continue
        if scale == 0 or difference != difference:
            return float("inf")
        worst = max(worst, difference / scale)
    return worst


def quantity_difference(case_name, values, reference):
    """the difference of the quantity of a run to the one of the reference:
    the relative difference of the values (the largest one over the groups
    of values), or for an error the ratio of the errors; None if one of
    them is missing"""
    if values is None or reference is None or len(values) != len(reference):
        return None
    if QUANTITIES[case_name].kind == "error":
        r = max(reference)
        return max(values) / r if r > 0 else None
    return _relative_difference(case_name, values, reference)


def check_quantity(case_name, precision, values, reference):
    """if the quantity of a run agrees with the one of the reference (the
    solvers of the program on the CPU); None if one of them is missing"""
    if values is None or reference is None or len(values) != len(reference):
        return None
    q = QUANTITIES[case_name]
    if q.kind == "error":
        return max(values) <= max(ERROR_FLOORS[precision],
                                  10 * max(reference))
    return (_relative_difference(case_name, values, reference) <=
            QUANTITY_TOLERANCES[precision])


def _collect_outputs(case_name, precision, results, name, start):
    """moves the output files written by a run (modified after `start`, a
    time of the epoch) from the directory of the program to
    results/outputs/<run>; returns this directory"""
    outputs = results / "outputs" / name
    outputs.mkdir(parents=True, exist_ok=True)
    directory = CASES[case_name].directory(precision)
    q = QUANTITIES.get(case_name)
    for pattern in (q.files if q else []):
        for f in directory.glob(pattern):
            if f.is_file() and f.stat().st_mtime >= start:
                shutil.move(str(f), str(outputs / f.name))
    return outputs


# analysis of a run from its log, without any additional run: the number of
# iterations of the non linear solvers (corrections of Newton's method or
# of the fixed point methods, the iteration 0 being the initial residual),
# and the times of the phases given by the time table of the profiler of
# mfem-mgis (mean over the processes, summed over the sections of the same
# name)
# (the lines of the residuals only: the inexact methods also print the
# tolerance of the linear solver of each iteration, and the outer loop of
# the polycrystal and of the cermet on the macroscopic stress prints
# "Fixed Point iteration" with |res|)
_ITERATION = re.compile(r"^(?:Newton|Fixed point) iteration\s+(\d+)\s*:\s*"
                        r"\|\|r\|\|", re.MULTILINE)
# a line of the time table: name, number of calls and the times (one
# process: time and part; several processes: minimum, mean, maximum, part
# and imbalance)
_SECTION = re.compile(r"^ \|[ |]*(?:-->|>) (.+?)\s*\|\s*(\d+)\s*\|(.*)$",
                      re.MULTILINE)
# phases: name in runs.csv and sections of the profiler (the behaviours are
# integrated in the estimation of the new unknowns, before the computation
# of the residual)
PHASES = {
    "total": ["root"],
    "nonlinear_solve": ["NLEP::solve"],
    "residual": ["NS::processNewUnknownsEstimate", "NS::computeResidual"],
    "tangent": ["NS::getJacobian"],
    "linear_setup": ["NS::setLinearSolverOperator"],
    "linear_solve": ["MFEM::Mult(r,c)"],
    "post_processings": ["NLEP::executePostProcessings"],
}


def _analysis(log_path):
    """the number of non linear iterations, the number of resolutions (calls
    to the solve method of the problems: time steps, sub-steps, iterations
    of a coupling scheme) and the times of the phases of a run"""
    try:
        text = log_path.read_text(errors="replace")
    except OSError:
        return "", "", ""
    iterations = sum(1 for i in _ITERATION.findall(text) if int(i) > 0)
    sections = {}
    calls = {}
    for name, count, rest in _SECTION.findall(text):
        try:
            values = [float(v.strip().rstrip("%"))
                      for v in rest.split("|") if v.strip()]
        except ValueError:
            continue
        if not values:
            continue
        # mean over the processes, or the time of a single process
        t = values[1] if len(values) >= 5 else values[0]
        sections[name] = sections.get(name, 0.0) + t
        calls[name] = calls.get(name, 0) + int(count)
    phases = {k: sum(sections.get(n, 0.0) for n in names)
              for k, names in PHASES.items()
              if any(n in sections for n in names)}
    resolutions = calls.get("NLEP::solve", "")
    return (str(iterations), str(resolutions),
            ";".join(f"{k}={v:.4g}" for k, v in phases.items()))


FIELDS = ["date", "case", "size", "order", "precision", "configuration",
          "label", "side", "processes", "repeat", "status", "seconds",
          "unknowns", "quantity", "iterations", "resolutions", "phases",
          "command"]


def _amgx(case, single, precision="double", aggregation=False,
          parameters=None):
    """the solver of AmgX of a case (`AmgX` alone is one cycle of its
    multigrid, without convergence criteria), with the classical multigrid
    or the aggregation one, and parameters of AmgX added to the default
    configuration of mfem-mgis (`AmgXParameters`, a parameter without scope
    being given to the multigrid); in single precision, the relative
    tolerance of AmgXPCG is at least 1e-6 (floor of AmgX)"""
    base = _tolerance(case, precision)
    tolerance = max(base, 1e-6) if single else base
    precision = ", Precision: \"Single\"" if single else ""
    if aggregation:
        precision += ", Multigrid: \"Aggregation\""
    if parameters:
        precision += ", AmgXParameters: {" + parameters + "}"
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


# options of the aggregation multigrid of AmgX: aggregates of 8 unknowns
# (fewer levels; the multicolour DILU smoother fails on 5 cases of 11),
# multicolour Gauss-Seidel smoother in a new scope (redefining the scope
# "smoother" of the default configuration is refused by AmgX)
_AMGX_AGG8 = 'selector: "SIZE_8"'
_AMGX_GS = ('selector: "SIZE_2", "smoother(gs)": "MULTICOLOR_GS", '
            '"gs:symmetric_GS": 1')


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
        amgx_agg=_amgx(case, False, precision, aggregation=True),
        amgx_agg8=_amgx(case, False, precision, aggregation=True,
                        parameters=_AMGX_AGG8),
        amgx_gs=_amgx(case, False, precision, aggregation=True,
                      parameters=_AMGX_GS))


def command(case, configuration, processes, size, precision="double",
            order=1):
    """the command of a run, and its additional environment"""
    c = configuration
    sz = case.sizes[order][size]
    cmd = ([case.program] + case.options +
           case.precision_options.get(precision, []) + sz.options +
           ["--order", str(sz.order)] + c.options +
           ([] if c.coupled else SIDES[c.side]))
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
    # the initial pool of hypre on the GPU is reduced when the GPU is shared
    # by several processes (one pool each) or used by other libraries
    if c.side == "gpu" and (processes > 1 or
                            c.linear in _NOT_HYPRE_LINEAR or
                            c.non_linear in _NOT_HYPRE_NON_LINEAR):
        environment.update(SMALL_HYPRE_POOL)
    if c.side == "mixed":
        environment.update(MIXED_HYPRE_POOL)
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
    # the end of the program is waited for without polling (`Popen.wait`
    # with a timeout polls every 50 ms, which would round the times), a
    # timer stopping it after `timeout` seconds
    expired = threading.Event()

    def stop():
        expired.set()
        _stop(p)
    timer = threading.Timer(timeout, stop) if timeout is not None else None
    try:
        if timer is not None:
            timer.start()
        status = p.wait()
    except BaseException:
        if timer is not None:
            timer.cancel()
        _stop(p)
        raise
    if timer is not None:
        timer.cancel()
        timer.join()
    return "timeout" if expired.is_set() else status


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
            # the other users of the GPU waiting for the lock get it before
            # the next run of the driver
            time.sleep(3)


@contextlib.contextmanager
def _case_lock(directory):
    """holds a lock of the directory of a program during a run and the
    collection of its output files: two drivers running the same program at
    the same time (screening with the GPU while the CPU is measured) would
    exchange their output files"""
    with open(directory / ".benchmark-slides.lock", "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def run(results, case_name, configuration_name, precision, processes,
        size, repeat, dry_run=False, timeout=None, first_repeat=0,
        order=1):
    """runs a configuration of a case, appends the results to
    results/runs.csv; returns the statuses and the times of the runs"""
    case = CASES[case_name]
    c = CONFIGURATIONS[case_name][configuration_name]
    outcomes = []
    if precision not in c.precisions():
        print(f"{case_name}, {configuration_name}: not available in "
              f"{precision}")
        return outcomes
    cmd, extra = command(case, c, processes, size, precision, order)
    path = results / "runs.csv"
    if not dry_run:
        (results / "logs").mkdir(parents=True, exist_ok=True)
        _check_header(path)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env.update(extra)
    others = {str(p) for p in MGIS_LIBRARIES.values()}
    env["LD_LIBRARY_PATH"] = ":".join(
        [str(MGIS_LIBRARIES[precision])] +
        [p for p in env.get("LD_LIBRARY_PATH", "").split(":")
         if p and p not in others])
    prefix = " ".join(f"{k}={v}" for k, v in extra.items())
    for r in range(first_repeat, first_repeat + repeat):
        name = (f"{case_name}-{size}-o{order}-{precision}-"
                f"{configuration_name}-np{processes}-{r}")
        print((prefix + " " if prefix else "") + shlex.join(cmd))
        if dry_run:
            continue
        log_path = results / "logs" / (name + ".log")
        with _case_lock(case.directory(precision)):
            with _gpu_lock(c.side), open(log_path, "w") as log:
                epoch = time.time()
                start = time.perf_counter()
                status = execute(cmd, case.directory(precision), log, env,
                                 timeout)
                seconds = time.perf_counter() - start
            outputs = _collect_outputs(case_name, precision, results, name,
                                       epoch)
        values = QUANTITIES[case_name].read(log_path, outputs) \
            if (case_name in QUANTITIES and status == 0) else None
        analysis = _analysis(log_path)
        row = {"date": datetime.datetime.now().isoformat(timespec="seconds"),
               "case": case_name, "size": size, "order": order,
               "precision": precision, "configuration": configuration_name,
               "label": c.label(), "side": c.side, "processes": processes,
               "repeat": r, "status": status, "seconds": f"{seconds:.3f}",
               "unknowns": (str(UNKNOWNS[(case_name, order, size)])
                            if (case_name, order, size) in UNKNOWNS
                            else _unknowns(log_path)),
               "quantity": (" ".join(repr(v) for v in values)
                            if values else ""),
               "iterations": analysis[0], "resolutions": analysis[1],
               "phases": analysis[2],
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


def calibrate(results, precision, size, cases, timeout, dry_run, order=1):
    """runs each case once with the solvers of the program on 8 processes,
    and on the GPU (hypre with BoomerAMG when the solver of the program is
    refused on the GPU)"""
    for name in cases:
        run(results, name, "cpu-native", precision, 8, size, 1, dry_run,
            timeout, order=order)
        gpu = ("gpu-minres-muelu-chebyshev" if name == "ex8" else
               "gpu-hypre-amg" if name in ("ex2", "ex5") else "gpu-native")
        run(results, name, gpu, precision, 1, size, 1, dry_run, timeout,
            order=order)


def _done(results, case_name, precision, size, order=1):
    """the configurations and numbers of processes already run (whatever
    their status), to resume an interrupted screening"""
    path = results / "runs.csv"
    if not path.exists():
        return set()
    with open(path, newline="") as f:
        return {(row["configuration"], int(row["processes"]))
                for row in csv.DictReader(f)
                if row["case"] == case_name and row["precision"] == precision
                and row["size"] == size and row["order"] == str(order)}


def _best_of_reference(reference, case_name, precision, size, keep,
                       reference_precision="double", order=1):
    """the configurations to screen in `precision` from the screening made
    in `reference_precision` in the results `reference`: the solvers of the
    program on the CPU and, among the configurations available in
    `precision`, the `keep` best ones on the CPU and the `keep` best ones
    with the GPU"""
    times = _screened(reference, case_name, reference_precision, size,
                      order, differing=False)
    if not times:
        print(f"{case_name}, order {order}: no successful run in the "
              f"reference {reference}: all the configurations are screened")
        return None
    configurations = CONFIGURATIONS[case_name]
    # the solvers of the program on the CPU, also in the organization of
    # the master branch (coupled loops), for the comparison of the
    # organizations
    selected = {"cpu-native", "cpu-native-coupled"} & set(configurations)
    for on_cpu in (True, False):
        ranked = []
        for (name, _), t in sorted(times.items(), key=lambda x: x[1]):
            c = configurations.get(name)
            if (c is None or name in ranked or name in selected or
                    (c.side == "cpu") != on_cpu or
                    precision not in c.precisions()):
                continue
            ranked.append(name)
        selected.update(ranked[:keep])
    return selected


def screen(results, case_name, precision, sides, size, timeout, dry_run,
           tiers=(1,), only=None, order=1):
    """runs the configurations of a case of the given tiers on the given
    sides, except the ones already run (and, if `only` is given, the ones
    which are not in it)"""
    configurations = CONFIGURATIONS[case_name]
    for name, c in configurations.items():
        if c.side not in sides or c.tier not in tiers:
            continue
        if only is not None and name not in only:
            continue
        for processes in c.processes:
            # read again before each run: another driver may have run it
            done = _done(results, case_name, precision, size, order)
            if (name, processes) in done:
                print(f"{case_name}, {name}, {processes} processes: "
                      "already run")
                continue
            run(results, case_name, name, precision, processes, size, 1,
                dry_run, _adaptive_timeout(results, case_name, precision,
                                           size, order, c, name, timeout),
                order=order)


# a configuration is stopped in the screening when it takes more than
# ADAPTIVE_FACTOR times the best time already measured on its side (the CPU,
# or the GPU), and ADAPTIVE_MINIMUM seconds at least: it can't be selected
ADAPTIVE_FACTOR = 3
ADAPTIVE_MINIMUM = 60


def _adaptive_timeout(results, case_name, precision, size, order,
                      configuration, name, timeout):
    """the time limit of a run of the screening: `timeout`, or a smaller one
    given by the best time already measured on the same side; the solvers
    of the program on the CPU always have the time limit `timeout`"""
    if name == "cpu-native":
        return timeout
    configurations = CONFIGURATIONS[case_name]
    on_cpu = configuration.side == "cpu"
    # (a fast run whose results differ does not shorten the limit)
    times = [t for (n, _), t in
             _screened(results, case_name, precision, size, order,
                       differing=False).items()
             if n in configurations and
             (configurations[n].side == "cpu") == on_cpu]
    if not times:
        return timeout
    limit = max(ADAPTIVE_MINIMUM, ADAPTIVE_FACTOR * min(times))
    return limit if timeout is None else min(timeout, limit)


def _quantity(row):
    """the values of the quantity of interest of a row of runs.csv"""
    q = row.get("quantity", "")
    return [float(v) for v in q.split()] if q else None


def _rows(results, case_name, precision, size, order):
    """the rows of runs.csv of a case, precision, size and order; the key
    "agrees" of a row is true when its quantity of interest agrees with the
    one of the solvers of the program on the CPU, false when it differs,
    None when one of them is missing"""
    path = results / "runs.csv"
    if not path.exists():
        return []
    with open(path, newline="") as f:
        rows = [row for row in csv.DictReader(f)
                if row["case"] == case_name and
                row["precision"] == precision and row["size"] == size and
                row["order"] == str(order)]
    reference = next((_quantity(r) for r in rows
                      if r["configuration"] == "cpu-native" and
                      r["status"] == "0" and _quantity(r) is not None), None)
    for row in rows:
        row["agrees"] = (True if row["configuration"] == "cpu-native" and
                         row["status"] == "0" else
                         check_quantity(case_name, precision, _quantity(row),
                                        reference))
    return rows


def _screened(results, case_name, precision, size, order=1,
              agreeing=False, differing=True):
    """the best time of each configuration and number of processes in the
    screening of a case (successful runs; only the ones whose quantity of
    interest agrees with the one of the solvers of the program on the CPU if
    `agreeing` is true; without the ones whose quantity differs if
    `differing` is false)"""
    times = {}
    for row in _rows(results, case_name, precision, size, order):
        if (row["status"] != "0" or
                (agreeing and row["agrees"] is not True) or
                (not differing and row["agrees"] is False)):
            continue
        k = (row["configuration"], int(row["processes"]))
        times[k] = min(times.get(k, float("inf")), float(row["seconds"]))
    return times


def finals(results, case_name, precision, screen_size, size, repeat,
           timeout, dry_run, long_run=120, order=1, tries=3,
           gpu_processes=(2, 4)):
    """runs at the measurement size the configurations selected by the
    screening: the solvers of the program on the CPU (best number of
    processes, or 8 if the screening did not finish), the best other CPU
    configuration and the two best configurations with the GPU. When a
    selected configuration fails at the measurement size, the next one of
    the screening is tried (`tries` configurations at most for each
    selection). The best configuration computing on the GPU with one
    process is also run with `gpu_processes` processes sharing the GPU (the
    parts computed on the host are then parallel). The solvers of the
    program on the CPU, the best other CPU configuration, the best GPU
    configuration and its best number of processes are repeated (`repeat`
    runs in all), unless their first run succeeds in more than `long_run`
    seconds or fails. The configurations of the results which are no
    longer defined are ignored, and the ones already run at this size are
    not run again"""
    times = _screened(results, case_name, precision, screen_size, order)
    agreeing = _screened(results, case_name, precision, screen_size, order,
                         agreeing=True)
    configurations = CONFIGURATIONS[case_name]

    def ranked(predicate):
        """the configurations by increasing time of the screening, the ones
        whose results agree with the solvers of the program first"""
        return [k for _, _, k in sorted((k not in agreeing, t, k)
                                        for k, t in times.items()
                                        if k[0] in configurations and
                                        predicate(configurations[k[0]],
                                                  k[0]))]

    def measure(key, r):
        """runs a configuration (resumed: the runs already made are
        counted), returns its best successful time at the measurement size,
        None if it failed"""
        name, processes = key
        made = _runs(results, case_name, precision, size, order).get(key, [])
        if not made:
            made = [o[1] if o[0] == 0 else None
                    for o in run(results, case_name, name, precision,
                                 processes, size, 1, dry_run, timeout,
                                 order=order)]
        # repeated while the first run succeeded in less than `long_run`
        if (r > len(made) and made and made[0] is not None and
                made[0] < long_run):
            run(results, case_name, name, precision, processes, size,
                r - len(made), dry_run, timeout, first_repeat=len(made),
                order=order)
        return _measured(results, case_name, precision, size,
                         order).get(key)

    def first_success(candidates, r, count):
        """measures the candidates in order until `count` of them succeed
        (`tries` attempts at most for each one); returns the successful
        ones"""
        succeeded = []
        for key in candidates[:tries * count]:
            if len(succeeded) == count:
                break
            if measure(key, r if not succeeded else 1) is not None:
                succeeded.append(key)
        return succeeded

    native = ranked(lambda c, n: n == "cpu-native")
    measure(native[0] if native else ("cpu-native", 8), repeat)
    # the same solvers with the coupled organization of the computations
    # (as in the master branch of mfem-mgis), with the same number of
    # processes
    if "cpu-native-coupled" in configurations:
        measure(("cpu-native-coupled", native[0][1] if native else 8),
                repeat)
    first_success(ranked(lambda c, n: c.side == "cpu" and
                         n not in ("cpu-native", "cpu-native-coupled")),
                  repeat, 1)
    gpu = first_success(ranked(lambda c, n: c.side != "cpu"), repeat, 2)
    # several processes sharing the GPU, for the best configuration
    # computing on the GPU with one process
    if gpu and configurations[gpu[0][0]].side == "gpu" and gpu[0][1] == 1:
        name = gpu[0][0]
        # several processes on one GPU were much slower in the first trials
        # (results/L6-polycrystal-2026-10-09): each number of processes is
        # stopped beyond 1.5 times the best time, and tried only if the
        # previous one was faster than one process
        for n in gpu_processes:
            measured = _measured(results, case_name, precision, size, order)
            times = {m: t for (c, m), t in measured.items()
                     if c == name and t is not None}
            if 1 not in times:
                break
            best_previous = min(t for m, t in times.items() if m < n)
            if n != gpu_processes[0] and best_previous >= times[1]:
                break
            if (name, n) not in measured:
                limit = 1.5 * best_previous
                run(results, case_name, name, precision, n, size, 1,
                    dry_run, limit if timeout is None
                    else min(timeout, limit), order=order)
        measured = _measured(results, case_name, precision, size, order)
        best_n = min((n for n in (1,) + tuple(gpu_processes)
                      if measured.get((name, n)) is not None),
                     key=lambda n: measured[(name, n)], default=1)
        # the best number of processes is repeated as the others
        if best_n != 1:
            measure((name, best_n), repeat)


def _runs(results, case_name, precision, size, order=1):
    """the times of the runs made of each configuration and number of
    processes at a size, in their order (None for a failed run; the runs
    whose quantity of interest differs are kept, and marked in the
    graphs)"""
    runs = {}
    for row in _rows(results, case_name, precision, size, order):
        k = (row["configuration"], int(row["processes"]))
        runs.setdefault(k, []).append(
            float(row["seconds"]) if row["status"] == "0" else None)
    return runs


def _measured(results, case_name, precision, size, order=1):
    """the best successful time of each configuration and number of
    processes at a size (the configurations which only failed are mapped
    to None)"""
    times = {}
    for k, ts in _runs(results, case_name, precision, size, order).items():
        ok = [t for t in ts if t is not None]
        times[k] = min(ok) if ok else None
    return times


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
        # (the coupled organization of the solvers of the program is not a
        # CPU configuration of the comparison)
        kinds = ["gpu"] if row["side"] != "cpu" else (
            ["cpu-native", "cpu"] if row["configuration"] == "cpu-native"
            else [] if row["configuration"] == "cpu-native-coupled"
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
    sc.add_argument("--reference", type=pathlib.Path, default=None,
                    help="results of the screening in double: only the "
                    "solvers of the program on the CPU and the best "
                    "configurations of this screening are run")
    sc.add_argument("--keep", type=int, default=5,
                    help="number of the best configurations of the "
                    "reference kept on each side (CPU, GPU)")
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
        a.add_argument("--order", type=int, default=1, choices=(1, 2),
                       help="order of the elements (sizes of each order)")
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
                              for o in sorted(case.sizes)
                              for k, v in case.sizes[o].items())
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
        if args.size not in CASES[args.case].sizes[args.order]:
            parser.error(f"unknown size '{args.size}'")
        for n in ([args.np] if args.np else cf.processes):
            run(args.results, args.case, args.configuration, args.precision,
                n, args.size, args.repeat, args.dry_run, args.timeout,
                order=args.order)
    elif args.action == "calibrate":
        cases = args.cases.split(",")
        unknown = [x for x in cases if x not in CASES]
        if unknown:
            parser.error("unknown cases: " + ", ".join(unknown))
        calibrate(args.results, args.precision, args.size, cases,
                  args.timeout, args.dry_run, args.order)
    elif args.action == "screen":
        sides = args.sides.split(",")
        unknown = [x for x in sides if x not in SIDES]
        if unknown:
            parser.error("unknown sides: " + ", ".join(unknown))
        only = (_best_of_reference(args.reference, args.case,
                                   args.precision, args.size, args.keep,
                                   order=args.order)
                if args.reference is not None else None)
        screen(args.results, args.case, args.precision, sides, args.size,
               args.timeout, args.dry_run,
               tuple(int(t) for t in args.tiers.split(",")), only,
               args.order)
    elif args.action == "finals":
        finals(args.results, args.case, args.precision, args.screen_size,
               args.size, args.repeat, args.timeout, args.dry_run,
               order=args.order)
    else:
        best(args.results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
