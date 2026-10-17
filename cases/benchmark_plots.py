#!/usr/bin/env python3
"""Graphs of the benchmark of the examples (cases/benchmark_slides.py) for
the slides: for each case, size, order and precision, the best time on the
CPU with the solvers of the program, the best time on the CPU with any
configuration, and the best times with the GPU, with the speedups; and an
overview of the speedups of all the cases. The graphs are written in the
output directory with an index (README.md) giving their context.

usage: benchmark_plots.py OUTPUT RESULTS [RESULTS...] [--threshold S]
"""

import argparse
import csv
import pathlib
import subprocess
import sys

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

# solvers of each program (the reference of the speedups)
NATIVE_NAMES = {
    "ex1": "CG of MFEM without preconditioner",
    "ex2": "MUMPS",
    "ex3": "CG of MFEM without preconditioner",
    "ex4": "FGMRES of hypre + BoomerAMG",
    "ex5": "MUMPS",
    "ex6": "PCG of hypre + BoomerAMG",
    "ex7": "GMRES of hypre + BoomerAMG",
    "ex8": "thermal: GMRES of hypre + BoomerAMG, mechanics: MUMPS",
    "polycrystal": "PCG of hypre + BoomerAMG",
    "bubble": "PCG of hypre + BoomerAMG",
    "cermet": "PCG of hypre + BoomerAMG",
}

CASES = {
    "ex1": "ex1, uniaxial tension of a cube, J2 plasticity",
    "ex2": "ex2, notched bar (ssna303), plane strain, finite strain "
           "plasticity",
    "ex3": "ex3, two-layer cube, elasticity",
    "ex4": "ex4, notched bar (ssna303) in 3D, finite strain plasticity",
    "ex5": "ex5, plate of Satoh",
    "ex6": "ex6, RVE made of two layers, non linear elasticity",
    "ex7": "ex7, RVE of a MOX fuel",
    "ex8": "ex8, thermomechanical plate (thermal and mechanical problems)",
    "polycrystal": "mm-opera-hpc, polycrystal of UO2, crystal plasticity",
    "bubble": "mm-opera-hpc, pressurized bubbles in a periodic cell, "
              "elasticity",
    "cermet": "mm-opera-hpc, cermet (UO2 grains, chromium interface)",
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
    "mfem-amg": "Krylov of MFEM + BoomerAMG",
    "cg-jacobi": "Krylov of MFEM + Jacobi",
    "muelu-jacobi": "Krylov of Belos + MueLu (Jacobi smoother)",
    "muelu-chebyshev": "Krylov of Belos + MueLu (Chebyshev smoother)",
    "muelu-sgs": "Krylov of Belos + MueLu (Gauss-Seidel smoother)",
    "amgx": "AmgX",
    "amgx-single": "AmgX in single precision",
    "cudss": "cuDSS (direct)",
    "mumps": "MUMPS (direct)",
    "tacho": "Tacho of Trilinos (direct)",
    "strumpack": "STRUMPACK (direct)",
    "petsc-gamg": "Krylov of PETSc + GAMG",
    "ginkgo-jacobi": "Krylov of Ginkgo + Jacobi",
    "magma-jacobi": "Krylov of MAGMA + Jacobi",
}
NON_LINEAR_NAMES = {
    "newton-ew": "inexact Newton (Eisenstat-Walker)",
    "k0-amg": "Newton-K0, K0 by one cycle of BoomerAMG",
    "k0-muelu": "Newton-K0, K0 by one cycle of MueLu",
    "k0-amgx": "Newton-K0, K0 by one cycle of AmgX",
    "k0-cudss": "Newton-K0, K0 factorized by cuDSS",
    "modified-newton": "modified Newton",
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
    if row["side"] == "cpu":
        where = f"CPU, {row['processes']} processes"
    elif row["side"] == "mixed":
        where = "GPU, linear solver on the CPU"
    else:
        where = "GPU"
    return f"{where}: {solvers(row['label'], row['case'])}"


def load(paths):
    """the successful runs, grouped by case, size, order and precision, and
    the runs of the solvers of the program on the CPU stopped by the time
    limit (lower bounds of their time)"""
    groups = {}
    stopped = {}
    for p in paths:
        path = pathlib.Path(p) / "runs.csv"
        if not path.exists():
            continue
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            required = {"case", "size", "precision", "configuration",
                        "label", "side", "processes", "status", "seconds"}
            if not required.issubset(reader.fieldnames or []):
                print(f"{path}: written by an older driver, skipped")
                continue
            for row in reader:
                # the files of the previous driver have no column order
                key = (row["case"], row["size"], row.get("order", ""),
                       row["precision"])
                if row["status"] == "0":
                    groups.setdefault(key, []).append(row)
                elif (row["status"] == "timeout" and
                      row["configuration"] == "cpu-native"):
                    stopped.setdefault(key, []).append(row)
    return groups, stopped


def best_by_configuration(rows):
    """the best run of each configuration and number of processes"""
    best = {}
    for row in rows:
        k = (row["configuration"], row["processes"])
        if (k not in best) or (float(row["seconds"]) <
                               float(best[k]["seconds"])):
            best[k] = row
    return sorted(best.values(), key=lambda r: float(r["seconds"]))


def summarize(rows):
    """the reference (best of the solvers of the program on the CPU), the
    best CPU configuration and the three best configurations with the GPU"""
    runs = best_by_configuration(rows)
    native = [r for r in runs if r["configuration"] == "cpu-native"]
    cpu = [r for r in runs if r["side"] == "cpu"]
    gpu = [r for r in runs if r["side"] != "cpu"]
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


def size_description(case, size):
    """the size of a run: its options in the driver"""
    sz = benchmark_slides.CASES.get(case, None)
    if sz is None or size not in sz.sizes:
        return f"size {size}"
    kind = {"screen": "screening size", "bench": "measurement size",
            "small": "small size"}.get(size, f"size {size}")
    return f"{kind} ({' '.join(sz.sizes[size].options) or 'default mesh'})"


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
    fig, ax = plt.subplots(figsize=(11, 1.2 + 0.62 * len(bars)), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    labels = [describe(r) for r, _ in bars]
    times = [float(r["seconds"]) for r, _ in bars]
    y = range(len(bars))[::-1]
    patches = ax.barh(list(y), times, height=0.5, color=[c for _, c in bars],
                      edgecolor=SURFACE, linewidth=2)
    if bound:
        patches[0].set_hatch("//")
        patches[0].set_edgecolor(INK_SECONDARY)
    ax.set_yticks(list(y), labels, color=INK, fontsize=9)
    ax.set_xlabel("wall time of the whole program (s), lower is better",
                  color=INK_SECONDARY, fontsize=9)
    right = max(times)
    ge = ">= " if bound else ""
    for yi, (r, color), t in zip(y, bars, times):
        text = f"{t:.1f} s"
        if bound and r is native:
            text = f"> {t:.0f} s (stopped)"
        if r["side"] != "cpu":
            if native is not None:
                text += f"   {ge}x{t_ref / t:.1f} vs CPU with the " \
                        "solvers of the example"
                if cpu is not None and cpu is not native:
                    text += f", {'>= ' if cpu_bound else ''}" \
                            f"x{t_cpu / t:.1f} vs best CPU"
            else:
                # the solvers of the program failed on the CPU
                text += f"   x{t_cpu / t:.1f} vs best CPU (the solvers " \
                        "of the example failed)"
        ax.text(t + right * 0.01, yi, text, va="center", fontsize=8.5,
                color=INK)
    ax.set_xlim(0, right * 2.4)
    case, size, order, precision = key
    unknowns = next((r["unknowns"] for r in rows if r["unknowns"]), "")
    title = CASES.get(case, case)
    detail = [f"order {order}", precision, size_description(case, size)]
    if unknowns:
        detail.insert(0, f"{int(unknowns):,} unknowns".replace(",", " "))
    fig.suptitle(title, x=0.01, ha="left", color=INK, fontsize=11.5,
                 fontweight="bold")
    ax.set_title(", ".join(detail), loc="left", color=INK_SECONDARY,
                 fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.text(0.01, 0.01, hardware, color=INK_SECONDARY, fontsize=7.5)
    name = f"{case}-{size}-o{order}-{precision}.png"
    fig.savefig(output / name, facecolor=SURFACE)
    plt.close(fig)
    t_gpu = float(gpu[0]["seconds"])
    return {"name": name, "case": case, "size": size, "order": order,
            "bound": bound, "cpu_bound": cpu_bound,
            "precision": precision, "unknowns": unknowns,
            "speedup_native": t_ref / t_gpu if native is not None else None,
            "speedup_best": t_cpu / t_gpu, "gpu": describe(gpu[0]),
            "cpu": describe(cpu) if cpu else "",
            "t_native": t_ref if native is not None else None,
            "t_cpu": t_cpu, "t_gpu": t_gpu}


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
    fig.savefig(output / name, facecolor=SURFACE)
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
    fig.savefig(output / name, facecolor=SURFACE)
    plt.close(fig)
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("results", nargs="+")
    parser.add_argument("--threshold", type=float, default=1.5,
                        help="speedup above which a case is highlighted")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cpu_name, gpu_name = machine()
    hardware = f"CPU: {cpu_name} (8 cores, 16 threads); GPU: {gpu_name}"
    summaries = []
    missing = []
    groups, stopped = load(args.results)
    keys = sorted(set(groups) | set(stopped))
    # a case measured at the measurement size is not shown at the screening
    # size any more (preliminary graph, removed)
    measured = {(c, o, p) for c, sz, o, p in keys if sz == "bench"}
    for key in keys:
        case, size, order, precision = key
        if size != "bench" and (case, order, precision) in measured:
            (args.output /
             f"{case}-{size}-o{order}-{precision}.png").unlink(missing_ok=True)
            continue
        s = case_graph(key, groups.get(key, []), stopped.get(key, []),
                       args.output, hardware)
        if isinstance(s, str):
            missing.append((key, s))
        else:
            summaries.append(s)
    o1 = overview(summaries, args.output, hardware, "speedup_native",
                  "Speedup with the GPU against the CPU with the solvers of "
                  "each example", "overview-speedup-native.png")
    o2 = overview(summaries, args.output, hardware, "speedup_best",
                  "Speedup with the GPU against the best configuration on "
                  "the CPU", "overview-speedup-best-cpu.png")
    o3 = precision_overview(summaries, args.output, hardware)
    lines = ["# Graphs of the benchmark of the examples", "",
             f"{hardware}. Generated by `cases/benchmark_plots.py` from "
             "the runs of `cases/benchmark_slides.py` (see "
             "`docs/benchmark-slides.md`). Wall time of the whole program; "
             "the CPU runs use 8 or 16 MPI processes, the GPU runs one "
             "process. The reference is the CPU with the solvers of each "
             "example (as shipped); the best CPU configuration is shown "
             "too.", ""]
    for o in (o1, o2, o3):
        if o:
            lines += [f"![{o}]({o})", ""]
    lines += ["The graphs at the screening size are preliminary (one run of "
              "each configuration, fewer steps): they are replaced by the "
              "graphs at the measurement size once the final runs of the "
              "case are made.", "",
              "| graph | case | size | unknowns | precision | GPU | speedup "
              "vs solvers of the example | speedup vs best CPU |",
              "|---|---|---|---|---|---|---|---|"]
    # sorted and highlighted by the speedup against the best CPU
    # configuration, the fairest comparison
    for s in sorted(summaries, key=lambda s: -s["speedup_best"]):
        mark = " **(highlight)**" if s["speedup_best"] >= \
            args.threshold else ""
        native = ((">= " if s["bound"] else "") +
                  f"x{s['speedup_native']:.1f}" if s["speedup_native"]
                  else "solvers of the example failed")
        best_cpu = (">= " if s["cpu_bound"] else "") + \
            f"x{s['speedup_best']:.1f}"
        size = {"bench": "measurement",
                "screen": "screening (preliminary)"}.get(s["size"],
                                                          s["size"])
        lines.append(f"| [{s['name']}]({s['name']}){mark} | {s['case']} | "
                     f"{size} | {s['unknowns']} | {s['precision']} | "
                     f"{s['gpu']} | "
                     f"{native} | {best_cpu} |")
    if missing:
        lines += ["", "Cases and sizes without a graph:", ""]
        lines += [f"- {c}, {sz}, order {o}, {pr}: {reason}"
                  for (c, sz, o, pr), reason in missing]
    (args.output / "README.md").write_text("\n".join(lines) + "\n")
    print(f"{len(summaries)} graphs in {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
