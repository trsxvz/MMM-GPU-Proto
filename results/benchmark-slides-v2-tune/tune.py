#!/usr/bin/env python3
"""Trials of the parameters of the GPU solvers before the second pass: each
variant is a linear solver (and optionally a non linear solver) description,
run on the GPU at the screening size through the driver (lock of the GPU,
quantity of interest, iterations); the results are appended to
results/benchmark-slides-v2-tune/<precision>/runs.csv.

usage: tune.py CASE VARIANT [VARIANT...] [--order O] [--precision P]
       tune.py --list
The lock of the GPU shared with other users is given by the environment
variable BENCHMARK_GPU_LOCK, if any; the libraries of the installation must
be in LD_LIBRARY_PATH (see docs/benchmark-slides.md).
"""

import argparse
import os
import pathlib
import signal
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cases"))
import benchmark_slides as b  # noqa: E402

LOCK = os.environ.get("BENCHMARK_GPU_LOCK") or None


def amg(options, krylov="{hypre}", key="Tolerance"):
    """a Krylov solver of hypre preconditioned by BoomerAMG with options"""
    return b._krylov(krylov, '{{Name: "HypreBoomerAMG", Options: '
                     '{{VerbosityLevel: 0' + (", " + options if options
                                              else "") + '}}}}', key)


def ginkgo(preconditioner, iterations=100000):
    return b._krylov("{ginkgo}", preconditioner, iterations=iterations)


# variants: name -> (linear solver, non linear solver or None, environment)
VARIANTS = {
    # BoomerAMG on the GPU
    "amg-default": (amg(""), None, {}),
    "amg-cheb": (amg("RelaxationType: 16, ChebyshevOrder: 2"), None, {}),
    "amg-cheb3": (amg("RelaxationType: 16, ChebyshevOrder: 3"), None, {}),
    "amg-l1j2": (amg("RelaxationType: 18, RelaxationSweeps: 2"), None, {}),
    "amg-agg": (amg("AggressiveCoarseningLevels: 1"), None, {}),
    "amg-theta07": (amg("StrengthThreshold: 0.7"), None, {}),
    "amg-theta025": (amg("StrengthThreshold: 0.25"), None, {}),
    "amg-interp14": (amg("InterpolationType: 14"), None, {}),
    "amg-interp18": (amg("InterpolationType: 18"), None, {}),
    "amg-pmax6": (amg("MaximumElementsPerRow: 6"), None, {}),
    "amg-theta05": (amg("StrengthThreshold: 0.5"), None, {}),
    "amg-cheb-theta05": (amg("RelaxationType: 16, ChebyshevOrder: 2, "
                             "StrengthThreshold: 0.5"), None, {}),
    "amg-cheb-keept": (amg("RelaxationType: 16, ChebyshevOrder: 2, "
                           "KeepTranspose: true"), None, {}),
    "amg-keept": (amg("KeepTranspose: true"), None, {}),
    "amg-cheb-agg-keept": (amg("RelaxationType: 16, ChebyshevOrder: 2, "
                               "AggressiveCoarseningLevels: 1, "
                               "KeepTranspose: true"), None, {}),
    # aggressive coarsening with the interpolations computed on the GPU
    # (two-stage extended+i, extended)
    "amg-agg5": (amg("AggressiveCoarseningLevels: 1, "
                     "AggressiveInterpolationType: 5"), None, {}),
    "amg-agg7": (amg("AggressiveCoarseningLevels: 1, "
                     "AggressiveInterpolationType: 7"), None, {}),
    # other preconditioners of hypre suited to the GPU: factorized sparse
    # approximate inverse (symmetric positive definite operators), ILU(0)
    # by blocks
    "hypre-fsai": (b._krylov("{hypre}", '{{Name: "HypreFSAI"}}',
                             "Tolerance", 100000), None, {}),
    "hypre-ilu0": (b._krylov("{hypre}", '{{Name: "HypreILU"}}',
                             "Tolerance", 100000), None, {}),
    # Ginkgo
    "gko-jacobi": (ginkgo('{{Name: "GinkgoJacobi"}}'), None, {}),
    # block Jacobi whose blocks are stored in the lowest precision keeping
    # their conditioning (adaptive precision)
    "gko-jacobi-adaptive": (ginkgo('{{Name: "GinkgoJacobi", Options: {{'
                                   'Configuration: {{storage_optimization: '
                                   '"autodetect"}}}}}}'), None, {}),
    "gko-jacobi1": (ginkgo('{{Name: "GinkgoJacobi", Options: {{'
                           'Configuration: {{max_block_size: 1}}}}}}'),
                    None, {}),
    "gko-jacobi3": (ginkgo('{{Name: "GinkgoJacobi", Options: {{'
                           'Configuration: {{max_block_size: 3}}}}}}'),
                    None, {}),
    "gko-mg": (ginkgo('{{Name: "GinkgoMultigrid"}}', 5000), None, {}),
    "gko-mg-relax06": (ginkgo('{{Name: "GinkgoMultigrid", Options: {{'
                              'Configuration: {{smoother_relax: 0.6, '
                              'smoother_iters: 2}}}}}}', 5000), None, {}),
    "gko-mg-levels": (ginkgo('{{Name: "GinkgoMultigrid", Options: {{'
                             'Configuration: {{smoother_relax: 0.6, '
                             'smoother_iters: 2, max_levels: 25, '
                             'min_coarse_rows: 500}}}}}}', 5000), None, {}),
    "gko-mg-w": (ginkgo('{{Name: "GinkgoMultigrid", Options: {{'
                        'Configuration: {{smoother_relax: 0.6, '
                        'smoother_iters: 2, max_levels: 25, '
                        'min_coarse_rows: 500, cycle: "w"}}}}}}', 5000),
                 None, {}),
    "gko-mg-relax05": (ginkgo('{{Name: "GinkgoMultigrid", Options: {{'
                              'Configuration: {{smoother_relax: 0.5, '
                              'smoother_iters: 2}}}}}}', 5000), None, {}),
    "gko-mg-relax06-i1": (ginkgo('{{Name: "GinkgoMultigrid", Options: {{'
                                 'Configuration: {{smoother_relax: 0.6}}}}}}',
                                 5000), None, {}),
    "gko-mg-relax06-i3": (ginkgo('{{Name: "GinkgoMultigrid", Options: {{'
                                 'Configuration: {{smoother_relax: 0.6, '
                                 'smoother_iters: 3}}}}}}', 5000), None, {}),
    # K-cycles (Krylov iterations on the coarse levels, the cycle suited
    # to the pairwise aggregation of Pgm): variable preconditioner, hence
    # the flexible solvers FCG or GCR ({flexible})
    "gko-kcycle": (b._krylov("{flexible}", '{{Name: "GinkgoMultigrid", '
                             'Options: {{Configuration: {{cycle: "{kcycle}", '
                             'smoother_relax: 0.6, smoother_iters: 2, '
                             'max_levels: 25, min_coarse_rows: 500}}}}}}'),
                   None, {}),
    "gko-kcycle-relax09": (b._krylov("{flexible}", '{{Name: '
                                     '"GinkgoMultigrid", Options: {{'
                                     'Configuration: {{cycle: "{kcycle}", '
                                     'max_levels: 25, min_coarse_rows: '
                                     '500}}}}}}'), None, {}),
    # GMRES with a compressed Krylov basis (vectors stored in a reduced
    # precision: less memory traffic in the orthogonalization)
    "gko-cbgmres-jacobi": (b._krylov("GinkgoCBGMRES",
                                     '{{Name: "GinkgoJacobi"}}',
                                     iterations=100000), None, {}),
    "gko-cbgmres-mg06": (b._krylov("GinkgoCBGMRES", '{{Name: '
                                   '"GinkgoMultigrid", Options: {{'
                                   'Configuration: {{smoother_relax: 0.6, '
                                   'smoother_iters: 2}}}}}}'), None, {}),
    "gko-ic": (ginkgo('{{Name: "GinkgoIC"}}'), None, {}),
    "gko-ilu-isai": (ginkgo('{{Name: "GinkgoILUISAI"}}'), None, {}),
    "gko-ic-isai": (ginkgo('{{Name: "GinkgoICISAI"}}'), None, {}),
    # MueLu
    "muelu-cheb": (b.LINEAR_SOLVERS["muelu-chebyshev"][0], None, {}),
    "muelu-cheb3": ('{{Name: "BelosSolver", Options: {{Solver: "{belos}", '
                    'RelativeTolerance: {tolerance}, '
                    'MaximumNumberOfIterations: 5000, Preconditioner: '
                    '{{Name: "MueLu", Options: {{Parameters: {{"smoother: '
                    'type": "CHEBYSHEV", "smoother: params": {{"chebyshev: '
                    'degree": 3, "chebyshev: eigenvalue max iterations": 30, '
                    '"chebyshev: boost factor": 1.2}}, "reuse: type": '
                    '"RP"}}}}}}}}}}', None, {}),
    "muelu-cheb-drop": ('{{Name: "BelosSolver", Options: {{Solver: '
                        '"{belos}", RelativeTolerance: {tolerance}, '
                        'MaximumNumberOfIterations: 5000, Preconditioner: '
                        '{{Name: "MueLu", Options: {{Parameters: {{"smoother: '
                        'type": "CHEBYSHEV", "smoother: params": '
                        '{{"chebyshev: eigenvalue max iterations": 30, '
                        '"chebyshev: boost factor": 1.2}}, "aggregation: '
                        'drop tol": 0.02, "reuse: type": "RP"}}}}}}}}}}',
                        None, {}),
}


# aggregation multigrid of AmgX with other smoothers suited to the GPU:
# l1-Jacobi, Chebyshev polynomial, multicolour ILU (built for each case,
# the solver depending on the symmetry)
AMGX_PARAMETERS = {
    "amgx-l1": 'selector: "SIZE_2", "smoother(l1)": "JACOBI_L1"',
    "amgx-cheb": ('selector: "SIZE_2", "smoother(cheb)": "CHEBYSHEV_POLY", '
                  '"cheb:chebyshev_polynomial_order": 2'),
    "amgx-ilu": 'selector: "SIZE_4", "smoother(ilu)": "MULTICOLOR_ILU"',
    "amgx-agg8": 'selector: "SIZE_8"',
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", nargs="?")
    parser.add_argument("variants", nargs="*")
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--precision", default="double")
    parser.add_argument("--size", default="screen")
    parser.add_argument("--non-linear", default=None,
                        help="non linear solver of the driver")
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    if args.list:
        print("\n".join(list(VARIANTS) + list(AMGX_PARAMETERS)))
        return 0
    # SIGTERM stops the run in progress (see the driver)
    signal.signal(signal.SIGTERM, lambda *a: sys.exit(143))
    b.GPU_LOCK = LOCK
    results = ROOT / "results" / "benchmark-slides-v2-tune" / args.precision
    for v in args.variants:
        if v in AMGX_PARAMETERS:
            # the description is formatted again by the driver
            d = b._amgx(b.CASES[args.case], False, args.precision,
                        aggregation=True, parameters=AMGX_PARAMETERS[v])
            VARIANTS[v] = (d.replace("{", "{{").replace("}", "}}"), None,
                           {})
        linear, non_linear, environment = VARIANTS[v]
        symmetric = b.CASES[args.case].symmetric
        linear = (linear.replace("{flexible}", "GinkgoFCG" if symmetric
                                 else "GinkgoGCR")
                  .replace("{kcycle}", "kfcg" if symmetric else "kgcr"))
        name = "tune-" + v + ("-" + args.non_linear if args.non_linear
                              else "")
        if v.startswith("gko-"):
            # number of iterations of each resolution printed by MFEM
            linear = linear.replace("Options: {{RelativeTolerance",
                                    "Options: {{VerbosityLevel: 1, "
                                    "RelativeTolerance", 1)
        b.LINEAR_SOLVERS[name] = (linear, b._BOTH)
        if not (v.startswith("amg-") or v.startswith("hypre-")):
            # the pool of hypre on the GPU is reduced for the other libraries
            b._NOT_HYPRE_LINEAR.add(name)
        b.CONFIGURATIONS[args.case][name] = b.Configuration(
            "gpu", name, args.non_linear, environment=environment)
        b.run(results, args.case, name, args.precision, 1, args.size, 1,
              False, args.timeout, order=args.order)
    return 0


if __name__ == "__main__":
    sys.exit(main())
