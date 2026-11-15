#!/usr/bin/env python3
"""Tables and figures of a scalability study from the points of
scaling-table.py.

usage: scaling-figures.py <table> <output directory> [base grains per node]

The weak series is made of the points at the base load (default 600 grains
per node); each RVE with points on more nodes than its base one gives a
strong series. Efficiencies: weak T(1) / T(N), strong T(N0) N0 / (T(N) N),
on the time to solution (duration of the Slurm step). Writes tables.md
(Markdown, in French as the reports) and, in PNG and PDF, weak.*, strong.*
and speedup.* (matplotlib).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator

# colours by entity: the GPUs, the cores, then the series of the speedups
GPU, CPU = "#2a78d6", "#eb6834"
SERIES = ["#52514e", "#1baf7a", "#e87ba4", "#6250d6"]
MARKERS = ["o", "s", "^", "D"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
NAMES = {"gpu": "GPU (4 MI300A par nœud)", "cpu": "CPU (96 cœurs par nœud)"}
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


def series(points, base):
    """weak series and strong series {grains: points}, as lists of
    (nodes, row) sorted by nodes"""
    weak = sorted((n, r) for (g, n), r in points.items() if g == base * n)
    strong = {}
    for (g, n), r in sorted(points.items()):
        strong.setdefault(g, []).append((n, r))
    return weak, {g: s for g, s in strong.items()
                  if len(s) > 1 and s[0][0] * base == g}


def efficiency(rows, mode, strong):
    """nodes -> efficiency of the time to solution against the first point"""
    n0, first = rows[0]
    t0 = first[f"{mode}_elapsed"]
    return [(n, ratio(t0 * (n0 if strong else n), r[f"{mode}_elapsed"] * n))
            for n, r in rows if r[f"{mode}_elapsed"] and t0]


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
        percent = lambda e: "-" if e is None else f"{french(100 * e)} %"
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


def node_axis(ax, nodes):
    ax.set_xscale("log", base=2)
    ax.set_xticks(nodes)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Nœuds")


def plain_log(axis):
    axis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    axis.set_minor_formatter(FuncFormatter(lambda y, _: ""))


def curve(ax, rows, key, **style):
    xy = [(n, r[key]) for n, r in rows if r[key]]
    if xy:
        ax.plot(*zip(*xy), **style)
    return xy


def save(fig, directory, name):
    for extension in ("png", "pdf"):
        fig.savefig(os.path.join(directory, f"{name}.{extension}"))
    plt.close(fig)


def time_and_efficiency(axes, rows, strong):
    """time to solution (with the ideal one for a strong series) above the
    efficiency, for the GPUs and the cores"""
    top, bottom = axes
    nodes = [n for n, _ in rows]
    for mode, colour, marker in (("gpu", GPU, "o"), ("cpu", CPU, "s")):
        xy = curve(top, rows, f"{mode}_elapsed", color=colour, marker=marker,
                   label=NAMES[mode])
        if strong and xy:
            n0, t0 = xy[0]
            top.plot(nodes, [t0 * n0 / n for n in nodes], color=colour,
                     linewidth=1, linestyle="--")
        e = efficiency(rows, mode, strong)
        if e:
            bottom.plot([n for n, _ in e], [100 * v for _, v in e],
                        color=colour, marker=marker)
    if strong:
        top.plot([], [], color=MUTED, linewidth=1, linestyle="--",
                 label="idéal")
    top.set_yscale("log")
    plain_log(top.yaxis)
    bottom.axhline(100, color=MUTED, linewidth=1, linestyle="--")
    bottom.set_ylim(0, 115)
    for ax in axes:
        node_axis(ax, nodes)
    top.set_xlabel("")


def weak_figure(rows, directory):
    fig, axes = plt.subplots(2, 1, figsize=(6, 6.5), sharex=True,
                             height_ratios=[3, 2])
    time_and_efficiency(axes, rows, False)
    axes[0].set_title("Scalabilité faible : 600 grains (19,8 M DDL) par nœud",
                      loc="left")
    axes[0].set_ylabel("Temps jusqu'à la solution (s)")
    axes[1].set_ylabel("Efficacité T(1) / T(N) (%)")
    axes[0].legend()
    save(fig, directory, "weak")


def strong_figure(strong, directory):
    fig, axes = plt.subplots(2, len(strong), figsize=(4 * len(strong), 6.5),
                             sharex="col", height_ratios=[3, 2], squeeze=False)
    for column, (grains, rows) in enumerate(strong.items()):
        time_and_efficiency(axes[:, column], rows, True)
        unknowns = rows[0][1]["unknowns"] / 1e6
        axes[0, column].set_title(
            f"{grains} grains, {french(unknowns, 1)} M DDL", loc="left")
    axes[0, 0].set_ylabel("Temps jusqu'à la solution (s)")
    axes[1, 0].set_ylabel("Efficacité (%)")
    axes[0, 0].legend()
    fig.suptitle("Scalabilité forte", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, directory, "strong")


def speedup_figure(weak, strong, directory):
    """speedup of the GPUs over the cores of the same nodes, against the
    nodes and against the unknowns per GPU"""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    named = [("faible (600 grains par nœud)", weak)] + [
        (f"forte, {g} grains", rows) for g, rows in strong.items()]
    for (name, rows), colour, marker in zip(named, SERIES, MARKERS):
        xy = [(n, r["unknowns"] / (4e6 * n),
               ratio(r["cpu_elapsed"], r["gpu_elapsed"])) for n, r in rows]
        xy = [p for p in xy if p[2]]
        if xy:
            nodes, load, s = zip(*xy)
            axes[0].plot(nodes, s, color=colour, marker=marker, label=name)
            axes[1].plot(load, s, color=colour, marker=marker)
    node_axis(axes[0], [n for n, _ in weak])
    axes[1].set_xscale("log")
    plain_log(axes[1].xaxis)
    axes[1].set_xlabel("DDL par GPU (millions)")
    axes[0].set_ylabel("Speedup GPU / CPU à nœuds égaux")
    axes[0].set_ylim(bottom=0)
    axes[0].set_title("Temps jusqu'à la solution CPU / GPU", loc="left")
    axes[1].set_title("Les mêmes points, selon la charge des GPU", loc="left")
    axes[0].legend(title="Scalabilité")
    fig.tight_layout()
    save(fig, directory, "speedup")


def main():
    points = read(sys.argv[1])
    directory = sys.argv[2]
    base = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    os.makedirs(directory, exist_ok=True)
    weak, strong = series(points, base)
    with open(os.path.join(directory, "tables.md"), "w") as out:
        out.write("## Scalabilité faible\n\n")
        table(weak, False, out)
        for grains, rows in strong.items():
            out.write(f"\n## Scalabilité forte, {grains} grains\n\n")
            table(rows, True, out)
    weak_figure(weak, directory)
    if strong:
        strong_figure(strong, directory)
    speedup_figure(weak, strong, directory)


if __name__ == "__main__":
    main()
