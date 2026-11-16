#!/usr/bin/env python3
"""Wall clock time of the functions of one thread, from the samples of Linux
perf (fm-chain/profile-host.slurm).

usage: perf script -i perf.data -F tid,time,ip,sym,dso |
         perf-walltime.py [--root NAME] [--skip N] [--depth N] [--min S]
                          [--thread TID] [--timeline]

`perf report` counts cycles of the processor spent in the code of the user:
the waits (GPU, MPI, files) and the time spent in the kernel (first accesses
to the pages of a large array) are not in its shares. Here each sample is
given the time elapsed since the previous sample of its thread, a wait or a
system call being given to the function in which the thread resumes. The
times of the tree therefore add up to the wall clock time of its root.

Prints the call tree of the thread --thread (default: the one of the first
sample, the main thread of the process) below the function whose name
contains --root (default `main`): time during which each function was on the
call stack, in s and in share of the root, and its own time. --skip N starts
at the (N+1)-th call of the root (calls separated by samples outside of it),
--depth limits the depth of the tree (default 6), --min the time of the
functions printed (default 0.1 s). The functions without symbol are named
after their library.

A call stack deeper than the copy made by `perf record --call-graph
dwarf,<size>`, or going through a library without unwinding tables (the HIP
runtime), is cut: its outermost functions are missing. The cut stacks sampled
during a call of the root are gathered below `(cut call stacks)`, from their
outermost known function.

--timeline prints instead the successive periods of the first call of the
root kept (see --skip), in the order of the time: a period lasts as long as
the function called by the root stays the same (`(own)` for the code of the
root itself), which separates the steps of a function that does one thing
after the other."""
import argparse
import re
import sys

FRAME = re.compile(r"^\t\s*[0-9a-f]+\s+(.*?)(?:\s+\(([^()]*)\))?$")
HEADER = re.compile(r"^(?:.*\s)?(\d+)\s+(\d+\.\d+):")
# outermost functions of a complete call stack
OUTERMOST = ("_start", "main", "__libc_start_call_main", "clone3",
             "__clone3", "start_thread")
CUT = "(cut call stacks)"


def samples(stream):
    """(thread, time, call stack from the outermost function) of each
    sample"""
    tid, time, stack = None, None, []
    for line in stream:
        if line.startswith("\t"):
            m = FRAME.match(line.rstrip("\n"))
            if m:
                name, dso = m.group(1), m.group(2) or ""
                if name == "[unknown]":
                    name = "[" + dso.rsplit("/", 1)[-1] + "]"
                stack.append(name)
            continue
        if tid is not None and line.strip() == "":
            yield tid, time, stack[::-1]
            tid, stack = None, []
            continue
        m = HEADER.match(line)
        if m:
            tid, time, stack = int(m.group(1)), float(m.group(2)), []
    if tid is not None:
        yield tid, time, stack[::-1]


class Node:
    def __init__(self):
        self.total, self.own, self.children = 0.0, 0.0, {}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--root", default="main")
    p.add_argument("--skip", type=int, default=0)
    p.add_argument("--depth", type=int, default=6)
    p.add_argument("--min", type=float, default=0.1)
    p.add_argument("--thread", type=int)
    p.add_argument("--timeline", action="store_true")
    a = p.parse_args()
    tree, thread, previous = Node(), a.thread, None
    call, inside, root, cut = -1, False, None, []
    periods = []

    def add(dt, stack):
        if call == a.skip:
            step = " < ".join(stack[1:3][::-1]) if len(stack) > 1 else "(own)"
            if periods and periods[-1][0] == step:
                periods[-1][1] += dt
            else:
                periods.append([step, dt])
        node = tree
        node.total += dt
        for f in stack:
            node = node.children.setdefault(f, Node())
            node.total += dt
        node.own += dt

    for tid, time, stack in samples(sys.stdin):
        if thread is None:
            thread = tid
        if tid != thread:
            continue
        dt = 0.0 if previous is None else time - previous
        previous = time
        at = next((i for i, f in enumerate(stack) if a.root in f), None)
        if at is None:
            # a cut stack belongs to the call if the root is sampled again
            # before a complete stack without it
            if inside and not (stack and stack[0] in OUTERMOST):
                cut.append((dt, stack))
            else:
                inside, cut = False, []
            continue
        if not inside:
            inside, call = True, call + 1
        if call >= a.skip:
            root = stack[at]
            for d, s in cut:
                add(d, [root, CUT] + s)
            add(dt, stack[at:])
        cut = []

    def show(name, node, depth, total):
        if node.total < a.min or depth > a.depth:
            return
        print(f"{node.total:9.2f} s {100 * node.total / total:5.1f} % "
              f"(own {node.own:6.2f} s) {'  ' * depth}{name[:150]}")
        for n, c in sorted(node.children.items(), key=lambda x: -x[1].total):
            show(n, c, depth + 1, total)

    if tree.total == 0.0:
        sys.exit(f"no sample of the thread {thread} below `{a.root}`")
    if a.timeline:
        t = 0.0
        for step, dt in periods:
            if dt >= a.min:
                print(f"{t:9.2f} s  + {dt:6.2f} s  {step[:170]}")
            t += dt
        return
    for name, node in tree.children.items():
        show(name, node, 0, tree.total)


if __name__ == "__main__":
    main()
