#!/usr/bin/env python3
"""Graphs of the benchmark of the examples (cases/benchmark_slides.py) for
the slides: for each case, size, order and precision, the best time on the
CPU with the solvers of the program, the best time on the CPU with any
configuration, and the best times with the GPU, with the speedups; and an
overview of the speedups of all the cases. The graphs are written in the
output directory with an index (README.md) giving their context.

usage: benchmark_plots.py OUTPUT RESULTS [RESULTS...] [--threshold S]
                          [--agreeing-only] [--by-precision]
"""

import argparse
import csv
import pathlib
import subprocess
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

# the sizes of the cases (options and order), from the driver
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import benchmark_slides  # noqa: E402

# reference palette (light surface): categorical slots 1 and 2, a neutral
# gray for the reference (solvers of the program on the CPU), ink and grid
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
GPU_COLOR = "#2a78d6"
CPU_COLOR = "#eb6834"
REFERENCE_COLOR = "#a8a69f"
# overview: the precisions (categorical slots 1 and 3, values labelled)
DOUBLE_COLOR = "#2a78d6"
FLOAT_COLOR = "#1baf7a"
# phases of the computations (categorical slots 1 to 6 in their fixed
# order, and a neutral gray for the rest), from the time table of the
# profiler (column phases of runs.csv, see `_phases` in the driver): the
# prediction of each time step is counted in the integration, the tangent
# operator and the linear solves
PHASE_COLORS = {
    "residual": ("behaviours, inner forces, update of the state", "#2a78d6"),
    "tangent": ("tangent operator and its assembly", "#eb6834"),
    "linear_setup": ("setup of the linear solvers", "#1baf7a"),
    "linear_solve": ("linear solves (setup of hypre included)", "#eda100"),
    "post_processings": ("post-processings", "#e87ba4"),
    "initialization": ("initialization (MPI, device, mesh, problems)",
                       "#008300"),
    "other": ("rest", "#a8a69f"),
}

# the bars of the runs stopped by the time limit (lower bounds), whose
# results differ or which are not checked are hatched (not with the option
# --agreeing-only, whose graphs have no mark: the texts give the lower
# bounds)
HATCH = True

# solvers of each program (the reference of the speedups)
NATIVE_NAMES = {
    "ex1": "CG of MFEM without preconditioner",
    "ex2": "MUMPS",
    "ex3": "CG of MFEM without preconditioner",
    "ex3-o2": "CG of MFEM without preconditioner",
    "ex4": "FGMRES of hypre + BoomerAMG",
    "ex5": "MUMPS",
    "ex6": "PCG of hypre + BoomerAMG",
    "ex6-o2": "PCG of hypre + BoomerAMG",
    "ex7": "GMRES of hypre + BoomerAMG",
    "ex8": "thermal: GMRES of hypre + BoomerAMG, mechanics: MUMPS",
    "polycrystal": "PCG of hypre + BoomerAMG",
    "polycrystal-r2": "PCG of hypre + BoomerAMG",
    "polycrystal-o2": "PCG of hypre + BoomerAMG",
    "bubble": "PCG of hypre + BoomerAMG",
    "cermet": "PCG of hypre + BoomerAMG",
    "cermet-r2": "PCG of hypre + BoomerAMG",
}

CASES = {
    "ex1": "ex1, uniaxial tension of a cube, J2 plasticity",
    "ex2": "ex2, notched bar (ssna303), plane strain, finite strain "
           "plasticity",
    "ex3": "ex3, two-layer cube, elasticity",
    "ex3-o2": "ex3, two-layer cube, elasticity, quadratic elements",
    "ex4": "ex4, notched bar (ssna303) in 3D, finite strain plasticity",
    "ex5": "ex5, plate of Satoh",
    "ex6": "ex6, RVE made of two layers, non linear elasticity",
    "ex6-o2": "ex6, RVE made of two layers, non linear elasticity, "
              "quadratic elements",
    "ex7": "ex7, RVE of a MOX fuel",
    "ex8": "ex8, thermomechanical plate (thermal and mechanical problems)",
    "polycrystal": "mm-opera-hpc, polycrystal of UO2, crystal plasticity",
    "polycrystal-r2": "mm-opera-hpc, polycrystal of UO2, crystal "
                      "plasticity, refined mesh",
    "polycrystal-o2": "mm-opera-hpc, polycrystal of UO2, crystal "
                      "plasticity, quadratic elements",
    "bubble": "mm-opera-hpc, pressurized bubbles in a periodic cell, "
              "elasticity",
    "cermet": "mm-opera-hpc, cermet (UO2 grains, chromium interface)",
    "cermet-r2": "mm-opera-hpc, cermet (UO2 grains, chromium interface), "
                 "refined mesh",
}


def machine():
    """the processor and the GPU of the machine"""
    cpu = gpu = ""
    try:
        for line in subprocess.run(["lscpu"], capture_output=True, text=True,
                                   env={"LC_ALL": "C"}).stdout.splitlines():
            if line.startswith("Model name:"):
                cpu = line.split(":", 1)[1].strip()
        gpu = subprocess.run(["nvidia-smi", "--query-gpu=name",
                              "--format=csv,noheader"], capture_output=True,
                             text=True).stdout.splitlines()[0].strip()
    except (OSError, IndexError):
        pass
    return cpu, gpu


# readable names of the solvers of the configurations
LINEAR_NAMES = {
    "native": "solvers of the example",
    "hypre-amg": "Krylov of hypre + BoomerAMG",
    "hypre-amg-agg": "Krylov of hypre + BoomerAMG (aggressive coarsening)",
    "mfem-amg-lagged": "Krylov of MFEM + BoomerAMG built every 5 "
                       "iterations",
    "mfem-amg": "Krylov of MFEM + BoomerAMG",
    "cg-jacobi": "Krylov of MFEM + Jacobi",
    "muelu-jacobi": "Krylov of Belos + MueLu (Jacobi smoother)",
    "muelu-chebyshev": "Krylov of Belos + MueLu (Chebyshev smoother)",
    "muelu-sgs": "Krylov of Belos + MueLu (Gauss-Seidel smoother)",
    "amgx": "AmgX",
    "amgx-single": "AmgX in single precision",
    "amgx-agg": "AmgX (aggregation)",
    "amgx-agg8": "AmgX (aggregation, aggregates of 8)",
    "amgx-gs": "AmgX (aggregation, multicolour Gauss-Seidel smoother)",
    "cudss": "cuDSS (direct)",
    "mumps": "MUMPS (direct)",
    "tacho": "Tacho of Trilinos (direct)",
    "strumpack": "STRUMPACK (direct)",
    "petsc-gamg": "Krylov of PETSc + GAMG",
    "petsc-gamg-gpu": "Krylov of PETSc + GAMG (plain aggregation, built "
                      "once)",
    "ginkgo-jacobi": "Krylov of Ginkgo + Jacobi",
    "ginkgo-amg": "Krylov of Ginkgo + its multigrid (PGM)",
    "magma-jacobi": "Krylov of MAGMA + Jacobi",
}
NON_LINEAR_NAMES = {
    "newton-ew": "inexact Newton (Eisenstat-Walker)",
    "k0-amg": "Newton-K0, K0 by one cycle of BoomerAMG",
    "k0-muelu": "Newton-K0, K0 by one cycle of MueLu",
    "k0-amgx": "Newton-K0, K0 by one cycle of AmgX",
    "k0-cudss": "Newton-K0, K0 factorized by cuDSS",
    "k0-petsc": "Newton-K0, K0 by one cycle of GAMG (PETSc)",
    "k0-muelu-sgs": "Newton-K0, K0 by one cycle of MueLu (Gauss-Seidel "
                    "smoother)",
    "modified-newton": "modified Newton",
    "modified-newton-anderson": "modified Newton, Anderson acceleration",
    "k0-anderson": "fixed point on K0 (one cycle of BoomerAMG), Anderson "
                   "acceleration",
}
# configurations of ex8, given by the options of the example (mechanics)
OPTION_NAMES = {
    "-lsMc MINRESSolver -pcMc MueLu": "mechanics: MINRES + MueLu",
    "-lsMc CGSolver -pcMc MueLu": "mechanics: CG + MueLu",
    "-mls chebyshev": "(Chebyshev smoother)",
    "-mls jacobi -mld 0.3": "(Jacobi smoother)",
    "-lsMc CuDSSSolver -mra": "mechanics: cuDSS (direct)",
    "-lsMc HyprePCG -pcMc HypreBoomerAMG": "mechanics: PCG of hypre + "
                                           "BoomerAMG",
    "-oMc byVDIM -amgMc Elasticity": "(elasticity options)",
    "-nlsMc NewtonK0": "mechanics: Newton-K0",
}


def solvers(label, case=None):
    """readable names of the solvers of a configuration (its label in
    runs.csv: linear solver, non linear solver, matrix-free, options)"""
    parts = [p.strip() for p in label.split(",")]
    if parts == ["native"] and case in NATIVE_NAMES:
        return f"solvers of the example ({NATIVE_NAMES[case]})"
    linear, rest = parts[0], parts[1:]
    non_linear = [p for p in rest if p in NON_LINEAR_NAMES]
    matrix_free = "matrix-free" in rest
    options = [p for p in rest if p not in NON_LINEAR_NAMES and
               p != "matrix-free"]
    if non_linear and non_linear[0].startswith("k0"):
        # the linear solver of the problem only solves the predictions
        text = NON_LINEAR_NAMES[non_linear[0]]
    else:
        text = LINEAR_NAMES.get(linear, linear)
        if non_linear:
            text += ", " + NON_LINEAR_NAMES[non_linear[0]]
    if matrix_free:
        text += ", tangent operator without assembly"
    for o in options:
        for k, v in OPTION_NAMES.items():
            o = o.replace(k, v)
        text = (o if linear == "native" else text + ", " + o)
    return text


def describe(row):
    """the side and the solvers of a run, for the labels"""
    n = int(row["processes"])
    if row["side"] == "cpu":
        where = f"CPU, {n} processes"
    elif row["side"] == "mixed":
        where = f"GPU shared by {n} processes, linear solver on the CPU"
    else:
        where = "GPU" if n == 1 else f"GPU shared by {n} processes"
    return f"{where}: {solvers(row['label'], row['case'])}"


def load(paths):
    """the successful runs, grouped by case, size, order and precision, and
    the runs of the solvers of the program on the CPU stopped by the time
    limit (lower bounds of their time)"""
    groups = {}
    stopped = {}
    # the directory of the results of each group
    directory = {}
    for p in paths:
        path = pathlib.Path(p) / "runs.csv"
        if not path.exists():
            continue
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            required = {"case", "size", "order", "precision",
                        "configuration", "label", "side", "processes",
                        "status", "seconds", "quantity"}
            if not required.issubset(reader.fieldnames or []):
                print(f"{path}: written by an older driver, skipped")
                continue
            for row in reader:
                key = (row["case"], row["size"], row["order"],
                       row["precision"])
                directory[key] = pathlib.Path(p)
                if row["status"] == "0":
                    groups.setdefault(key, []).append(row)
                elif (row["status"] == "timeout" and
                      row["configuration"] == "cpu-native"):
                    stopped.setdefault(key, []).append(row)
    # the runs whose quantity of interest differs from the one of the
    # solvers of the program on the CPU are kept, marked (key "differs":
    # the relative difference, or the ratio of the errors), and listed. The
    # runs in float are checked against the reference in double (same
    # case, size and order), when it is measured: the solvers of a program
    # may give wrong results in float (MUMPS in single precision on ex8),
    # then marked too
    differs = {}
    for key, rows in groups.items():
        case, size, order, precision = key
        if case not in benchmark_slides.QUANTITIES:
            continue
        reference, kind = None, None
        if precision == "float":
            double = groups.get((case, size, order, "double"))
            if double is None:
                # the results in double next to the ones in float, as the
                # driver
                double = [r for r in benchmark_slides._read_rows(
                    directory[key].resolve().parent / "double", case,
                    "double", size, int(order)) if r["status"] == "0"]
            reference, kind = benchmark_slides._reference(double)
        in_double = reference is not None
        if not in_double:
            reference, kind = benchmark_slides._reference(rows)
        benchmark_slides.annotate_quantities(case, precision, rows,
                                             reference, in_double)
        for r in rows:
            if r["configuration"] == "cpu-native" and not in_double:
                continue
            r["reference"] = ({"native": "the solvers of the example",
                               "fastest CPU": "the fastest CPU "
                               "configuration"}.get(kind, "the reference") +
                              (" in double" if in_double else ""))
            if r["agrees"] is False:
                r["differs"] = r["difference"]
                differs.setdefault(key, {})[r["configuration"]] = \
                    r["differs"]
            elif r["agrees"] is None:
                # no quantity for this run or for the reference: not
                # checked, also marked and listed
                r["unchecked"] = True
                differs.setdefault(key, {})[r["configuration"]] = \
                    "not checked"
    return groups, stopped, differs


def keep_agreeing(groups):
    """leaves out of the groups the runs whose results differ from the
    reference or are not checked: the graphs then only show runs whose
    results agree. When all the runs of the solvers of the program are left
    out, the other runs of the case are marked, so that no speedup is given
    against them; returns the runs left out of each group (configuration
    and number of processes)"""
    left_out = {}
    for key, rows in groups.items():
        wrong = [r for r in rows if "differs" in r or "unchecked" in r]
        rows[:] = [r for r in rows
                   if "differs" not in r and "unchecked" not in r]
        native_wrong = (any(r["configuration"] == "cpu-native"
                            for r in wrong) and
                        not any(r["configuration"] == "cpu-native"
                                for r in rows))
        for r in rows:
            if native_wrong:
                r["native_left_out"] = True
        if wrong:
            left_out[key] = sorted(
                {f"{r['configuration']} ({r['processes']} process"
                 f"{'es' if r['processes'] != '1' else ''}"
                 f"{', not checked' if 'unchecked' in r else ''})"
                 for r in wrong})
    return left_out


def best_by_configuration(rows):
    """the best run of each configuration, whatever its number of processes
    (given in the labels)"""
    def order(r):
        # the runs whose results agree first
        return ("differs" in r, float(r["seconds"]))

    best = {}
    for row in rows:
        k = row["configuration"]
        if (k not in best) or (order(row) < order(best[k])):
            best[k] = row
    return sorted(best.values(), key=order)


def summarize(rows):
    """the reference (best of the solvers of the program on the CPU), the
    best CPU configuration and the three best configurations with the GPU"""
    runs = best_by_configuration(rows)
    native = [r for r in runs if r["configuration"] == "cpu-native"]
    # the coupled organization of the solvers of the program is shown in
    # the graphs of the organization, not as a CPU configuration; the best
    # configurations are taken among the runs whose results agree first
    agree_first = (lambda r: ("differs" in r, float(r["seconds"])))
    cpu = sorted((r for r in runs if r["side"] == "cpu" and
                  r["configuration"] != "cpu-native-coupled"),
                 key=agree_first)
    gpu = sorted((r for r in runs if r["side"] != "cpu"), key=agree_first)
    return (native[0] if native else None, cpu[0] if cpu else None,
            gpu[:3])


def style(ax):
    """recessive axes and grid"""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, length=0)
    ax.xaxis.grid(True, color=GRID, linewidth=1)
    ax.set_axisbelow(True)


def size_description(case, size, order):
    """the size of a run: its options in the driver"""
    sz = benchmark_slides.CASES.get(case, None)
    order = int(order)
    if sz is None or order not in sz.sizes or size not in sz.sizes[order]:
        return f"size {size}"
    kind = {"screen": "screening size", "bench": "measurement size",
            "small": "small size"}.get(size, f"size {size}")
    options = " ".join(sz.sizes[order][size].options) or "default mesh"
    return f"{kind} ({options})"


def case_graph(key, rows, stopped, output, hardware):
    """the graph of a case; returns its speedups, or the reason why there is
    no graph. If the solvers of the program did not finish on the CPU
    within the time limit, their time is a lower bound (hatched bar,
    speedups given as lower bounds)"""
    native, cpu, gpu = summarize(rows)
    bound = False
    if native is None and stopped:
        native = max(stopped, key=lambda r: float(r["seconds"]))
        bound = True
    if not gpu:
        return "no successful run with the GPU"
    if native is None and cpu is None:
        return "no successful run on the CPU"
    # the solvers of the program give wrong results (in float, against the
    # reference in double), shown or left out (see `keep_agreeing`): no
    # speedup against them
    native_wrong = ((native is not None and "differs" in native) or
                    any(r.get("native_left_out") for r in rows))
    bars = []
    if native is not None:
        bars.append((native, REFERENCE_COLOR))
    if cpu is not None and cpu is not native:
        bars.append((cpu, CPU_COLOR))
    bars += [(g, GPU_COLOR) for g in gpu]
    reference = native or cpu
    t_ref = float(reference["seconds"])
    t_cpu = float(cpu["seconds"]) if cpu else t_ref
    # the best CPU time is a lower bound if it is the one of the solvers of
    # the program stopped by the time limit
    cpu_bound = bound and (cpu is None or t_ref < t_cpu)
    t_cpu = min(t_cpu, t_ref) if bound else t_cpu
    # long descriptions on several lines
    labels = [textwrap.fill(describe(r), 70) for r, _ in bars]
    lines = sum(label.count("\n") for label in labels)
    fig, ax = plt.subplots(figsize=(11, 1.2 + 0.62 * len(bars) +
                                    0.18 * lines), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    times = [float(r["seconds"]) for r, _ in bars]
    y = range(len(bars))[::-1]
    patches = ax.barh(list(y), times, height=0.5, color=[c for _, c in bars],
                      edgecolor=SURFACE, linewidth=2)
    if bound and HATCH:
        patches[0].set_hatch("//")
        patches[0].set_edgecolor(INK_SECONDARY)
    # results which differ from the ones of the solvers of the example
    for patch, (r, _) in zip(patches, bars):
        if HATCH and ("differs" in r or "unchecked" in r):
            patch.set_hatch("xx" if "differs" in r else "..")
            patch.set_edgecolor(INK_SECONDARY)
    ax.set_yticks(list(y), labels, color=INK, fontsize=9)
    ax.set_xlabel("wall time of the whole program (s), lower is better",
                  color=INK_SECONDARY, fontsize=9)
    right = max(times)
    ge = ">= " if bound else ""
    for yi, (r, color), t in zip(y, bars, times):
        text = f"{t:.1f} s"
        if r.get("iterations"):
            text += f" ({r['iterations']} iterations)"
        if bound and r is native:
            text = f"> {t:.0f} s (stopped)"
        if r["side"] != "cpu":
            if native is not None and not native_wrong:
                text += f"   {ge}x{t_ref / t:.1f} vs CPU with the " \
                        "solvers of the example"
                if cpu is not None and cpu is not native:
                    text += f", {'>= ' if cpu_bound else ''}" \
                            f"x{t_cpu / t:.1f} vs best CPU"
            elif native_wrong:
                text += f"   x{t_cpu / t:.1f} vs best CPU"
            else:
                # the solvers of the program failed on the CPU
                text += f"   x{t_cpu / t:.1f} vs best CPU (the solvers " \
                        "of the example failed)"
        if "differs" in r:
            d = r["differs"]
            text += (f"   results differ from {r['reference']}" +
                     (f" ({d:.1e})" if d is not None else ""))
        elif "unchecked" in r:
            text += "   results not checked (no quantity of interest)"
        ax.text(t + right * 0.01, yi, text, va="center", fontsize=8.5,
                color=INK)
    ax.set_xlim(0, right * 2.4)
    case, size, order, precision = key
    unknowns = next((r["unknowns"] for r in rows if r["unknowns"]), "")
    title = CASES.get(case, case)
    detail = [f"order {order}", precision,
              size_description(case, size, order)]
    if unknowns:
        detail.insert(0, f"{int(unknowns):,} unknowns".replace(",", " "))
    fig.suptitle(title, x=0.01, ha="left", color=INK, fontsize=11.5,
                 fontweight="bold")
    ax.set_title(", ".join(detail), loc="left", color=INK_SECONDARY,
                 fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(0.01, 0.01, hardware, color=INK_SECONDARY, fontsize=7.5)
    name = f"{case}-{size}-o{order}-{precision}.png"
    fig.savefig(output / name, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    t_gpu = float(gpu[0]["seconds"])
    return {"name": name, "case": case, "size": size, "order": order,
            "bound": bound, "cpu_bound": cpu_bound,
            "precision": precision, "unknowns": unknowns,
            "speedup_native": (t_ref / t_gpu if native is not None and
                               not native_wrong else None),
            "speedup_best": t_cpu / t_gpu, "gpu": describe(gpu[0]),
            "cpu": describe(cpu) if cpu else "",
            "t_native": t_ref if native is not None else None,
            "t_cpu": t_cpu, "t_gpu": t_gpu,
            "gpu_differs": "differs" in gpu[0],
            "native_wrong": native_wrong,
            "iterations": "/".join(
                r.get("iterations", "") or "-"
                for r in (cpu if native_wrong else native or cpu, gpu[0]))}


def overview(summaries, output, hardware, field, title, name):
    """the speedups of all the cases, double and float side by side"""
    rows = [s for s in summaries if s[field] is not None]
    if not rows:
        return None
    cases = sorted({(s["case"], s["size"], s["order"]) for s in rows},
                   key=lambda k: list(CASES).index(k[0]) if k[0] in CASES
                   else 99)
    fig, ax = plt.subplots(figsize=(11, 1.4 + 0.55 * len(cases)), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    height = 0.36
    labelled = set()
    for i, (precision, color) in enumerate((("double", DOUBLE_COLOR),
                                            ("float", FLOAT_COLOR))):
        for j, k in enumerate(cases):
            s = next((x for x in rows if (x["case"], x["size"], x["order"])
                      == k and x["precision"] == precision), None)
            if s is None:
                continue
            yj = len(cases) - 1 - j + (0.2 if i == 0 else -0.2)
            # bars from 1 (same time on both sides) on the log scale
            ax.barh(yj, s[field] - 1, left=1, height=height, color=color,
                    edgecolor=SURFACE, linewidth=2,
                    label=None if precision in labelled else precision)
            labelled.add(precision)
            ge = ">= " if ((s["bound"] and field == "speedup_native") or
                           (s["cpu_bound"] and field == "speedup_best")) \
                else ""
            ax.text(s[field] * 1.04, yj, f"{ge}x{s[field]:.1f}",
                    va="center", fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ticks = [t for t in (0.25, 0.5, 1, 2, 5, 10, 20, 50, 100)
             if t <= 1.5 * max(s[field] for s in rows)]
    ax.set_xticks(ticks, [f"x{t:g}" for t in ticks])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(min(0.5, 0.8 * min(s[field] for s in rows)),
                2.2 * max(s[field] for s in rows))
    ax.axvline(1, color=INK_SECONDARY, linewidth=1)
    ax.set_yticks(range(len(cases))[::-1],
                  [f"{c} (order {o}{', screening' if s == 'screen' else ''})"
                   for c, s, o in cases], color=INK, fontsize=9)
    ax.set_xlabel("speedup (time on the CPU / time with the GPU), log "
                  "scale; 1 = same time", color=INK_SECONDARY, fontsize=9)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.suptitle(title, x=0.01, ha="left", color=INK, fontsize=11.5,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(0.01, 0.01, hardware, color=INK_SECONDARY, fontsize=7.5)
    fig.savefig(output / name, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return name


def precision_overview(summaries, output, hardware,
                       name="overview-float-vs-double.png"):
    """the gain of the float build over the double one for each case
    measured in both precisions: best time on the CPU in double / best time
    on the CPU in float, and the same with the GPU"""
    keys = sorted({(s["case"], s["size"], s["order"]) for s in summaries},
                  key=lambda k: list(CASES).index(k[0]) if k[0] in CASES
                  else 99)
    pairs = []
    for k in keys:
        by = {s["precision"]: s for s in summaries
              if (s["case"], s["size"], s["order"]) == k}
        if "double" in by and "float" in by:
            pairs.append((k, by["double"], by["float"]))
    if not pairs:
        return None
    fig, ax = plt.subplots(figsize=(11, 1.4 + 0.55 * len(pairs)), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    gains = []
    for j, (k, d, f) in enumerate(pairs):
        for i, (side, color, field) in enumerate(
                (("best CPU", CPU_COLOR, "t_cpu"),
                 ("best with the GPU", GPU_COLOR, "t_gpu"))):
            g = d[field] / f[field]
            gains.append(g)
            yj = len(pairs) - 1 - j + (0.2 if i == 0 else -0.2)
            ax.barh(yj, g - 1, left=1, height=0.36, color=color,
                    edgecolor=SURFACE, linewidth=2,
                    label=side if j == 0 else None)
            ax.text(max(g, 1) * 1.04, yj, f"x{g:.2f}", va="center",
                    fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ticks = [t for t in (0.5, 0.75, 1, 1.5, 2, 3, 4)
             if 0.6 * min(gains) <= t <= 1.5 * max(gains)]
    ax.set_xticks(ticks, [f"x{t:g}" for t in ticks])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(min(0.8, 0.8 * min(gains)), 1.6 * max(gains))
    ax.axvline(1, color=INK_SECONDARY, linewidth=1)
    ax.set_yticks(range(len(pairs))[::-1],
                  [f"{c} (order {o}{', screening' if sz == 'screen' else ''})"
                   for (c, sz, o), _, _ in pairs], color=INK, fontsize=9)
    ax.set_xlabel("time in double / time in float (best configuration of "
                  "each precision), log scale; 1 = same time",
                  color=INK_SECONDARY, fontsize=9)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.suptitle("Gain of single precision (float) over double precision",
                 x=0.01, ha="left", color=INK, fontsize=11.5,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(0.01, 0.01, hardware, color=INK_SECONDARY, fontsize=7.5)
    fig.savefig(output / name, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return name


def phases(row):
    """the times of the phases of a run (column phases of runs.csv), the
    rest of the total time being the phase "other"; None without the time
    table"""
    text = row.get("phases", "")
    if not text:
        return None
    values = dict((k, float(v)) for k, v in
                  (item.split("=") for item in text.split(";")))
    # rows of an older analysis (without the initialization): see
    # `benchmark_slides.py reanalyze`
    if "total" not in values or "initialization" not in values:
        return None
    # the time table starts after the initialization of MPI, of the device
    # and of the libraries: the rest of the wall time of the program is
    # counted in the initialization
    wall = max(float(row["seconds"]), values["total"])
    parts = {k: values[k] for k in PHASE_COLORS
             if k in values and k != "other"}
    parts["initialization"] = (parts.get("initialization", 0.0) + wall -
                               values["total"])
    parts["other"] = max(wall - sum(parts.values()), 0.0)
    return parts


def coupling_graph(key, rows, output, hardware):
    """the graph of the organization of the computations on the CPU for a
    case: the solvers of the program with the coupled organization (option
    --execution not given: integration point by point and computations
    element by element, as in the master branch of mfem-mgis) and with the
    organization of the prototype (--execution cpu: integration by batches
    and kernels), with the times of the phases; returns a summary, or None
    if one of them is missing"""
    best = {c: min((r for r in rows if r["configuration"] == c),
                   key=lambda r: float(r["seconds"]), default=None)
            for c in ("cpu-native-coupled", "cpu-native")}
    if None in best.values():
        return None
    case, size, order, precision = key
    fig, ax = plt.subplots(figsize=(11, 3.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    labels = []
    used = set()
    for y, (c, text) in enumerate(
            (("cpu-native", "organization of the prototype (batches and "
              "kernels)"),
             ("cpu-native-coupled", "coupled organization (as the master "
              "branch)"))):
        r = best[c]
        t = float(r["seconds"])
        parts = phases(r) or {"other": t}
        left = 0.0
        for k, (name, color) in PHASE_COLORS.items():
            if k not in parts or parts[k] <= 0:
                continue
            ax.barh(y, parts[k], left=left, height=0.5, color=color,
                    edgecolor=SURFACE, linewidth=2,
                    label=name if k not in used else None)
            used.add(k)
            left += parts[k]
        iterations = r.get("iterations", "")
        ax.text(max(t, left) * 1.01, y,
                f"{t:.1f} s" + (f", {iterations} iterations"
                                if iterations else ""),
                va="center", fontsize=8.5, color=INK)
        labels.append(f"CPU, {r['processes']} processes: {text}")
    ax.set_yticks([0, 1], labels, color=INK, fontsize=9)
    ax.set_xlabel("wall time of the whole program (s), lower is better",
                  color=INK_SECONDARY, fontsize=9)
    right = max(float(best[c]["seconds"]) for c in best)
    ax.set_xlim(0, right * 1.5)
    gain = float(best["cpu-native-coupled"]["seconds"]) / \
        float(best["cpu-native"]["seconds"])
    fig.suptitle(CASES.get(case, case) + ": organization of the "
                 f"computations on the CPU (x{gain:.2f})", x=0.01,
                 ha="left", color=INK, fontsize=11.5, fontweight="bold")
    ax.set_title(f"solvers of the example, order {order}, {precision}, "
                 f"{size_description(case, size, order)}", loc="left",
                 color=INK_SECONDARY, fontsize=9)
    # the legend under the axes (it hides no bar), in a band of the figure
    # reserved above the line of the hardware
    handles, names = ax.get_legend_handles_labels()
    rows = -(-len(names) // 3)
    bottom = 0.08 + 0.055 * rows
    fig.tight_layout(rect=(0, bottom, 1, 1))
    fig.legend(handles, names, frameon=False, fontsize=7.5, ncol=3,
               loc="upper left",
               bbox_to_anchor=(ax.get_position().x0, bottom - 0.005))
    fig.text(0.01, 0.01, hardware, color=INK_SECONDARY, fontsize=7.5)
    name = f"{case}-{size}-o{order}-{precision}-cpu-organization.png"
    fig.savefig(output / name, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return {"name": name, "case": case, "size": size, "order": order,
            "precision": precision, "gain": gain,
            "t_coupled": float(best["cpu-native-coupled"]["seconds"]),
            "t_prototype": float(best["cpu-native"]["seconds"])}


def coupling_overview(summaries, output, hardware,
                      name="overview-cpu-organization.png"):
    """the gain of the organization of the prototype on the CPU (time of
    the coupled organization / time of the prototype) of all the cases,
    double and float side by side"""
    if not summaries:
        return None
    keys = sorted({(s["case"], s["size"]) for s in summaries},
                  key=lambda k: list(CASES).index(k[0]) if k[0] in CASES
                  else 99)
    fig, ax = plt.subplots(figsize=(11, 1.4 + 0.55 * len(keys)), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    gains = [s["gain"] for s in summaries]
    labelled = set()
    for i, (precision, color) in enumerate((("double", DOUBLE_COLOR),
                                            ("float", FLOAT_COLOR))):
        for j, k in enumerate(keys):
            s = next((x for x in summaries if (x["case"], x["size"]) == k
                      and x["precision"] == precision), None)
            if s is None:
                continue
            yj = len(keys) - 1 - j + (0.2 if i == 0 else -0.2)
            ax.barh(yj, s["gain"] - 1, left=1, height=0.36, color=color,
                    edgecolor=SURFACE, linewidth=2,
                    label=None if precision in labelled else precision)
            labelled.add(precision)
            ax.text(max(s["gain"], 1) * 1.03, yj, f"x{s['gain']:.2f}",
                    va="center", fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ticks = [t for t in (0.5, 0.75, 1, 1.5, 2, 3, 4)
             if 0.6 * min(gains) <= t <= 1.6 * max(gains)]
    ax.set_xticks(ticks, [f"x{t:g}" for t in ticks])
    ax.set_xlim(min(0.8, 0.8 * min(gains)), 1.6 * max(gains))
    ax.axvline(1, color=INK_SECONDARY, linewidth=1)
    ax.set_yticks(range(len(keys))[::-1], [c for c, _ in keys], color=INK,
                  fontsize=9)
    ax.set_xlabel("time of the coupled organization (as the master branch) / "
                  "time of the prototype, log scale; 1 = same time",
                  color=INK_SECONDARY, fontsize=9)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.suptitle("Gain of the organization of the prototype on the CPU "
                 "(solvers of the examples)", x=0.01, ha="left", color=INK,
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(0.01, 0.01, hardware, color=INK_SECONDARY, fontsize=7.5)
    fig.savefig(output / name, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return name


def order_graphs(order, keys, groups, stopped, output, hardware,
                 threshold, differs, by_precision=False, left_out=None):
    """the graphs, the overviews and the index of an order, in `output`
    (the graphs of each case in a directory of its precision if
    `by_precision` is true, the overviews and the index in `output`);
    `left_out`: the configurations left out of each group because their
    results differ (see `keep_agreeing`), None if they are shown; returns
    the number of graphs"""
    output.mkdir(parents=True, exist_ok=True)

    def place(precision):
        """the directory of the graphs of a precision"""
        d = output / precision if by_precision else output
        d.mkdir(parents=True, exist_ok=True)
        return d

    def link(precision, name):
        """the link to a graph from the index"""
        return f"{precision}/{name}" if by_precision else name

    summaries = []
    missing = []
    # a case measured at the measurement size is not shown at the screening
    # size any more (preliminary graph, removed)
    measured = {(c, p) for c, sz, o, p in keys if sz == "bench"}
    for key in keys:
        case, size, _, precision = key
        if size != "bench" and (case, precision) in measured:
            (place(precision) /
             f"{case}-{size}-o{order}-{precision}.png").unlink(missing_ok=True)
            continue
        s = case_graph(key, groups.get(key, []), stopped.get(key, []),
                       place(precision), hardware)
        if isinstance(s, str):
            if left_out and key in left_out:
                s += " whose results agree"
            missing.append((key, s))
        else:
            summaries.append(s)
    elements = "linear" if order == "1" else "quadratic"
    # the overviews only show the speedups of the runs with the GPU whose
    # results agree (the others are in the graphs of the cases, hatched,
    # and listed in the index)
    valid = [s for s in summaries if not s["gpu_differs"]]
    o1 = overview(valid, output, hardware, "speedup_native",
                  "Speedup with the GPU against the CPU with the solvers of "
                  f"each example, order {order}",
                  "overview-speedup-native.png")
    o2 = overview(valid, output, hardware, "speedup_best",
                  "Speedup with the GPU against the best configuration on "
                  f"the CPU, order {order}", "overview-speedup-best-cpu.png")
    o3 = precision_overview(valid, output, hardware)
    organization = [c for c in (coupling_graph(k, groups.get(k, []),
                                               place(k[3]), hardware)
                                for k in keys if k[1] == "bench")
                    if c is not None]
    o4 = coupling_overview(organization, output, hardware)
    lines = [f"# Graphs of the benchmark of the examples, order {order} "
             f"({elements} elements)", "",
             f"{hardware}. Generated by `cases/benchmark_plots.py` from "
             "the runs of `cases/benchmark_slides.py` (see "
             "`docs/benchmark-slides.md`). Wall time of the whole program. "
             "The CPU runs use 8 or 16 MPI processes; the GPU runs one "
             "process, or several sharing the GPU; in the mixed placement, "
             "8 processes share the GPU and the linear solvers compute on "
             "the 8 cores. The reference is the CPU with the solvers of "
             "each example (as shipped); the best CPU configuration is "
             "shown too.", ""]
    for o in (o1, o2, o3, o4):
        if o:
            lines += [f"![{o}]({o})", ""]
    lines += ["The graphs at the screening size are preliminary (one run of "
              "each configuration, fewer steps): they are replaced by the "
              "graphs at the measurement size once the final runs of the "
              "case are made.", "",
              "| graph | case | size | unknowns | precision | GPU | speedup "
              "vs solvers of the example | speedup vs best CPU | non linear "
              "iterations (CPU / GPU) |",
              "|---|---|---|---|---|---|---|---|---|"]
    # sorted and highlighted by the speedup against the best CPU
    # configuration, the fairest comparison
    for s in sorted(summaries, key=lambda s: -s["speedup_best"]):
        mark = " **(highlight)**" if s["speedup_best"] >= threshold else ""
        native = ((">= " if s["bound"] else "") +
                  f"x{s['speedup_native']:.1f}" if s["speedup_native"]
                  else "solvers of the example wrong in float"
                  if s["native_wrong"] else "solvers of the example failed")
        best_cpu = (">= " if s["cpu_bound"] else "") + \
            f"x{s['speedup_best']:.1f}"
        size = {"bench": "measurement",
                "screen": "screening (preliminary)"}.get(s["size"],
                                                          s["size"])
        gpu = s["gpu"] + (" (results differ)" if s["gpu_differs"] else "")
        lines.append(f"| [{s['name']}]({link(s['precision'], s['name'])})"
                     f"{mark} | {s['case']} | "
                     f"{size} | {s['unknowns']} | {s['precision']} | "
                     f"{gpu} | {native} | {best_cpu} | "
                     f"{s['iterations']} |")
    if organization:
        lines += ["", "Organization of the computations on the CPU (solvers "
                  "of the example): coupled organization (option --execution "
                  "not given: integration point by point and computations "
                  "element by element, as in the master branch of mfem-mgis) "
                  "against the organization of the prototype (--execution "
                  "cpu: integration by batches and kernels):", "",
                  "| graph | case | precision | coupled | prototype | gain |",
                  "|---|---|---|---|---|---|"]
        lines += [f"| [{c['name']}]({link(c['precision'], c['name'])}) | "
                  f"{c['case']} | "
                  f"{c['precision']} | {c['t_coupled']:.1f} s | "
                  f"{c['t_prototype']:.1f} s | x{c['gain']:.2f} |"
                  for c in organization]
    if missing:
        lines += ["", "Cases and sizes without a graph:", ""]
        lines += [f"- {c}, {sz}, order {o}, {pr}: {reason}"
                  for (c, sz, o, pr), reason in missing]
    lines += ["", "The results of each run are checked: its quantity of "
              "interest (resultant force, mean stresses, error to the "
              "analytical solution...) is compared to the one of the "
              "solvers of the example on the CPU (or, when they did not "
              "finish, of the fastest configuration on the CPU), with a "
              "relative tolerance of 1e-5 in double and 1e-3 in float (an "
              "error to the analytical solution may be ten times the one of "
              "the reference). The runs in float are compared to the "
              "reference in double when it is measured: they differ beyond "
              "1e-3 and three times the median of the best differences of "
              "the configurations in float of the case (the accuracy "
              "reached in float), and always beyond 5 %. The best "
              "configurations are taken among the runs whose results agree "
              "first; when the solvers of the example give wrong results in "
              "float, no speedup is given against them. "]
    if left_out is not None:
        lines[-1] += ("Only the runs whose results agree are shown (option "
                      "--agreeing-only); the runs left out because their "
                      "results differ (or are not checked):")
        lines.append("")
        # the runs of the graphs: at the measurement size, or at the
        # screening size for a case not measured
        shown = {(c, pr) for c, sz, _, pr in keys if sz == "bench"}
        lines += ([f"- {c}, {sz}, {pr}: " + ", ".join(names)
                   for (c, sz, _, pr), names in sorted(left_out.items())
                   if sz == "bench" or (c, pr) not in shown]
                  or ["- none"])
        (output / "README.md").write_text("\n".join(lines) + "\n")
        return len(summaries)
    lines[-1] += ("The runs whose results differ are kept in the graphs, "
                  "hatched, with their relative difference (or the ratio of "
                  "the errors), to be examined; the runs without a quantity "
                  "of interest (or whose reference has none) are dotted, not "
                  "checked; the overviews leave out the cases whose best run "
                  "with the GPU gives results which differ:")
    lines.append("")
    if differs:
        lines += [f"- {c}, {sz}, {pr}: " + ", ".join(
                  f"{n} ({d})" if isinstance(d, str) else
                  f"{n} ({d:.1e})" if d is not None else n
                  for n, d in sorted(names.items()))
                  for (c, sz, _, pr), names in sorted(differs.items())]
    else:
        lines.append("- none")
    (output / "README.md").write_text("\n".join(lines) + "\n")
    return len(summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("results", nargs="+")
    parser.add_argument("--threshold", type=float, default=1.5,
                        help="speedup above which a case is highlighted")
    parser.add_argument("--agreeing-only", action="store_true",
                        help="only show the runs whose results agree with "
                        "the reference (no hatched bar)")
    parser.add_argument("--by-precision", action="store_true",
                        help="the graphs of each case in a directory of its "
                        "precision (double, float), the overviews and the "
                        "index in the directory of the order")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cpu_name, gpu_name = machine()
    hardware = f"CPU: {cpu_name} (8 cores, 16 threads); GPU: {gpu_name}"
    groups, stopped, differs = load(args.results)
    left_out = keep_agreeing(groups) if args.agreeing_only else None
    global HATCH
    HATCH = not args.agreeing_only
    keys = sorted(set(groups) | set(stopped))
    # one directory for each order of the elements
    orders = sorted({k[2] for k in keys})
    lines = ["# Graphs of the benchmark of the examples", "",
             "One directory for each order of the elements, each with its "
             "overviews (speedups in double and in float, gain of the float "
             "build) and its index:", ""]
    for order in orders:
        n = order_graphs(order, [k for k in keys if k[2] == order], groups,
                         stopped, args.output / f"order-{order}", hardware,
                         args.threshold,
                         {k: v for k, v in differs.items() if k[2] == order},
                         args.by_precision,
                         None if left_out is None else
                         {k: v for k, v in left_out.items() if k[2] == order})
        lines.append(f"- [order {order}](order-{order}/README.md): {n} "
                     "graphs")
        print(f"{n} graphs in {args.output / f'order-{order}'}")
    (args.output / "README.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
