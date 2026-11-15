#!/usr/bin/env python3
"""Tables and figures of a scalability study from the points of
scaling-table.py.

usage: scaling-figures.py <table> <output directory> [grains...]

The weak series is made of the points at the largest load, in grains per
node; each RVE with points on more nodes than its base one gives a strong
series. Efficiencies: weak T(1) / T(N), strong T(N0) N0 / (T(N) N), on the
time to solution (duration of the Slurm step). Writes tables.md (Markdown,
in French as the reports), with where the time goes, and, in PNG
(matplotlib), one graph each:
- largest.*: time to solution and energy of the largest point of the weak
  series, GPUs against cores;
- weak.*: times of the weak series, with the ideal ones, the GPUs above
  the cores, each with its linear axis;
- strong.*: the same for the strong series of the RVEs given (default:
  all).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator

# colours by entity: the GPUs and the cores, or the RVEs of the strong series
GPU, CPU = "#2a78d6", "#eb6834"
SIZES = ["#1baf7a", "#e87ba4", "#6250d6", "#eda100"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
NAMES = {"gpu": "GPU (4 MI300A par nœud)", "cpu": "CPU (96 cœurs par nœud)"}
STYLES = {"gpu": dict(marker="o", linestyle="-"),
          "cpu": dict(marker="s", linestyle=":")}
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "lines.linewidth": 2, "lines.markersize": 7, "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.bbox": "tight", "savefig.dpi": 200})


def read(path):
    """(grains, nodes) -> values of the line, `-` as None"""
    with open(path) as f:
        names = f.readline().split()
        points = {}
        for line in f:
            row = {n: None if w == "-" else float(w)
                   for n, w in zip(names, line.split())}
            points[int(row["grains"]), int(row["nodes"])] = row
    return points


def french(x, digits=0):
    return "-" if x is None else f"{x:.{digits}f}".replace(".", ",")


def ratio(a, b):
    return a / b if a and b else None


def series(points):
    """weak series and strong series {grains: points}, as lists of
    (nodes, row) sorted by nodes"""
    base = max(g / n for g, n in points)
    weak = sorted((n, r) for (g, n), r in points.items() if g == base * n)
    strong = {}
    for (g, n), r in sorted(points.items()):
        strong.setdefault(g, []).append((n, r))
    return weak, {g: s for g, s in strong.items()
                  if len(s) > 1 and s[0][0] * base == g}


def efficiency(rows, mode, strong, time="elapsed"):
    """nodes -> efficiency against the first point, on the time to solution
    or on another time of the points"""
    n0, first = rows[0]
    t0 = first[f"{mode}_{time}"]
    return [(n, ratio(t0 * (n0 if strong else n), r[f"{mode}_{time}"] * n))
            for n, r in rows if r[f"{mode}_{time}"] and t0]


def percent(e):
    return "-" if e is None else f"{french(100 * e)} %"


def table(rows, strong, out):
    out.write(
        "| Grains | DDL | Nœuds | DDL par GPU | DDL par cœur |"
        " Temps jusqu'à la solution GPU / CPU (s) | Efficacité GPU / CPU |"
        " Speedup | Temps mural GPU / CPU (s) | Speedup |"
        " Énergie GPU / CPU (MJ) | Rapport | Newton GPU / CPU |\n"
        + "|---" * 13 + "|\n")
    eg = dict(efficiency(rows, "gpu", strong))
    ec = dict(efficiency(rows, "cpu", strong))
    for n, r in rows:
        mj = lambda m: r[f"{m}_energy"] and r[f"{m}_energy"] / 1e6
        pair = lambda a, b, d=0: f"{french(a, d)} / {french(b, d)}"
        out.write(
            f"| {r['grains']:.0f} | {french(r['unknowns'] / 1e6, 1)} M | {n} |"
            f" {french(r['unknowns'] / (4e6 * n), 2)} M |"
            f" {french(r['unknowns'] / (96 * n))} |"
            f" {pair(r['gpu_elapsed'], r['cpu_elapsed'])} |"
            f" {percent(eg.get(n))} / {percent(ec.get(n))} |"
            f" {french(ratio(r['cpu_elapsed'], r['gpu_elapsed']), 1)} |"
            f" {pair(r['gpu_wall'], r['cpu_wall'])} |"
            f" {french(ratio(r['cpu_wall'], r['gpu_wall']), 1)} |"
            f" {pair(mj('gpu'), mj('cpu'), 2)} |"
            f" {french(ratio(r['cpu_energy'], r['gpu_energy']), 1)} |"
            f" {pair(r['gpu_newton'], r['cpu_newton'])} |\n")


def breakdown(rows, mode, out):
    """where the time to solution goes: launch of the processes (Slurm step
    without the program), setup (program without its loop), first step
    (loop without the wall clock time of the next steps), next steps"""
    out.write("| Nœuds | Lancement (s) | Mise en place (s) | Premier pas (s)"
              " | 19 pas suivants (s) | Total (s) |\n" + "|---" * 6 + "|\n")
    for n, r in rows:
        v = {k: r[f"{mode}_{k}"] for k in ("elapsed", "total", "loop", "wall")}
        if None in v.values():
            continue
        parts = (v["elapsed"] - v["total"], v["total"] - v["loop"],
                 v["loop"] - v["wall"], v["wall"], v["elapsed"])
        out.write(f"| {n} | " + " | ".join(french(x, 1) for x in parts) + " |\n")


def program(rows, strong, out):
    """the same comparison without the launch of the processes, which
    depends on what the nodes ran before: time of the program"""
    out.write("| Nœuds | Temps du programme GPU / CPU (s) | Efficacité GPU / CPU"
              " | Speedup |\n" + "|---" * 4 + "|\n")
    eg = dict(efficiency(rows, "gpu", strong, "total"))
    ec = dict(efficiency(rows, "cpu", strong, "total"))
    for n, r in rows:
        out.write(
            f"| {n} | {french(r['gpu_total'])} / {french(r['cpu_total'])} |"
            f" {percent(eg.get(n))} / {percent(ec.get(n))} |"
            f" {french(ratio(r['cpu_total'], r['gpu_total']), 1)} |\n")


def node_axis(ax, nodes):
    ax.set_xscale("log", base=2)
    ax.set_xticks(nodes)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Nœuds")


def unknowns(row):
    """unknowns of a point, in millions (M) or in billions, in words"""
    u = row["unknowns"]
    if u >= 1e9:
        return f"{french(u / 1e9, 2)} milliard de DDL"
    return f"{french(u / 1e6)} M DDL"


def curve(ax, rows, mode, strong, **style):
    """time to solution of a series and its ideal time, a thin grey line
    from its first point"""
    xy = [(n, r[f"{mode}_elapsed"]) for n, r in rows if r[f"{mode}_elapsed"]]
    if xy:
        n0, t0 = xy[0]
        nodes = [n0 * (rows[-1][0] / n0) ** (i / 50) for i in range(51)]
        ax.plot(nodes, [t0 * n0 / n if strong else t0 for n in nodes],
                color=MUTED, linewidth=0.8)
        ax.plot(*zip(*xy), **STYLES[mode], **style)
    return xy


def save(fig, directory, name):
    fig.savefig(os.path.join(directory, f"{name}.png"))
    plt.close(fig)


def largest_figure(nodes, row, directory):
    """bars of the GPUs and of the cores for the time to solution and for
    the energy of one point, with their ratios"""
    fig, axes = plt.subplots(1, 2, figsize=(7, 4))
    panels = (("elapsed", 1, "Temps jusqu'à la solution (s)", 0),
              ("energy", 1e6, "Énergie (MJ)", 1))
    for ax, (key, unit, name, digits) in zip(axes, panels):
        gpu, cpu = row[f"gpu_{key}"] / unit, row[f"cpu_{key}"] / unit
        ax.bar(["CPU", "GPU"], [cpu, gpu], width=0.5, color=[CPU, GPU])
        for x, value in enumerate((cpu, gpu)):
            ax.annotate(french(value, digits), (x, value), (0, 4),
                        textcoords="offset points", ha="center")
        ax.annotate(f"÷ {french(cpu / gpu, 1)}", (1, gpu), (0, 24),
                    textcoords="offset points", ha="center", fontsize=13,
                    fontweight="bold")
        ax.set_title(name, loc="left", fontsize=10)
        ax.grid(axis="x", visible=False)
        ax.set_axisbelow(True)
        ax.margins(y=0.12)
    fig.suptitle(
        f"{row['grains']:,.0f}".replace(",", " ")
        + f" grains, {unknowns(row)}, {nodes} nœuds :"
        f" {4 * nodes} MI300A contre {96 * nodes} cœurs", x=0.02, ha="left")
    fig.tight_layout()
    save(fig, directory, "largest")


def panels(nodes, title, height):
    """two graphs, the GPUs above the cores, each with its linear axis of
    the times from zero, on the same axis of the nodes"""
    fig, axes = plt.subplots(2, 1, figsize=(6.5, height), sharex=True)
    for ax, mode in zip(axes, ("gpu", "cpu")):
        ax.set_title(NAMES[mode], loc="left", fontsize=10)
        ax.set_ylabel("Temps jusqu'à la solution (s)")
    node_axis(axes[1], nodes)
    fig.suptitle(title, x=0.02, ha="left")
    return fig, axes


def finish(fig, axes, directory, name):
    for ax in axes:
        ax.set_ylim(0, 1.15 * ax.get_ylim()[1])
        ax.plot([], [], color=MUTED, linewidth=0.8, label="idéal")
        ax.legend(loc="lower left" if name == "weak" else "upper right")
    fig.tight_layout()
    save(fig, directory, name)


def weak_figure(rows, directory):
    n, first = rows[0]
    fig, axes = panels(
        [n for n, _ in rows],
        f"Scalabilité faible : {first['grains'] / n:.0f} grains par nœud", 6.5)
    for ax, mode, colour in zip(axes, ("gpu", "cpu"), (GPU, CPU)):
        for n, t in curve(ax, rows, mode, False, color=colour):
            ax.annotate(f"{french(t)} s", (n, t), (0, 8),
                        textcoords="offset points", ha="center", fontsize=9)
    # sizes of the first and of the last problem, under the axis
    for (n, r), side in ((rows[0], "left"), (rows[-1], "right")):
        axes[1].annotate(unknowns(r), (n, 0), (0, -30), ha=side, va="top",
                         xycoords=("data", "axes fraction"),
                         textcoords="offset points", color=MUTED)
    finish(fig, axes, directory, "weak")


def strong_figure(strong, directory):
    """one colour per RVE"""
    fig, axes = panels(
        sorted({n for rows in strong.values() for n, _ in rows}),
        "Scalabilité forte : chaque problème part du plus petit nombre\n"
        "de nœuds que permet la mémoire des GPU", 7.5)
    for ax, mode in zip(axes, ("gpu", "cpu")):
        for (grains, rows), colour in zip(strong.items(), SIZES):
            name = f"{grains} grains, {unknowns(rows[0][1])}"
            for n, t in curve(ax, rows, mode, True, color=colour, label=name):
                ax.annotate(french(t), (n, t), (0, 7), ha="center",
                            textcoords="offset points", fontsize=9)
    finish(fig, axes, directory, "strong")


def main():
    points = read(sys.argv[1])
    directory = sys.argv[2]
    os.makedirs(directory, exist_ok=True)
    weak, strong = series(points)
    with open(os.path.join(directory, "tables.md"), "w") as out:
        out.write("## Scalabilité faible\n\n")
        table(weak, False, out)
        out.write("\nSans le lancement des processus :\n\n")
        program(weak, False, out)
        for mode in ("gpu", "cpu"):
            out.write(f"\nDécomposition du temps {mode.upper()} :\n\n")
            breakdown(weak, mode, out)
        for grains, rows in strong.items():
            out.write(f"\n## Scalabilité forte, {grains} grains\n\n")
            table(rows, True, out)
            out.write(f"\nSans le lancement des processus, {grains} grains :\n\n")
            program(rows, True, out)
            for mode in ("gpu", "cpu"):
                out.write(f"\nDécomposition du temps {mode.upper()},"
                          f" {grains} grains :\n\n")
                breakdown(rows, mode, out)
    largest_figure(*weak[-1], directory)
    weak_figure(weak, directory)
    shown = [int(g) for g in sys.argv[3:]] or list(strong)
    strong_figure({g: strong[g] for g in shown if g in strong}, directory)


if __name__ == "__main__":
    main()
