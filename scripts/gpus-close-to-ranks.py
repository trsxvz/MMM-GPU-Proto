#!/usr/bin/env python3
"""Order of the GPUs of a node for CUDA_VISIBLE_DEVICES such that the GPU of
the rank r (the r-th visible GPU, see `getDeviceIdentifier` in mfem-mgis) is
close to the cores of this rank: on its NUMA node if possible, otherwise on
its socket. All the GPUs stay visible (the GPU-aware MPI of Open MPI with
PSM2 requires the GPUs of the other processes).

usage: scripts/gpus-close-to-ranks.py processes cores-per-process

The rank r runs on the cores r * c to (r + 1) * c - 1, as srun places the
processes of `--cpus-per-task=c`. The NUMA nodes and sockets of the cores
are read with `lscpu -p=CPU,NODE,SOCKET`, the NUMA node of each GPU in
`nvidia-smi topo -m`. Prints nothing if the topology is unknown.
"""
import re
import subprocess
import sys


def run(command):
    try:
        return subprocess.run(command, capture_output=True, text=True,
                              check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return ""


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    processes, cores = int(sys.argv[1]), int(sys.argv[2])
    node_of, socket_of = {}, {}
    for line in run(["lscpu", "-p=CPU,NODE,SOCKET"]).splitlines():
        if line.startswith("#"):
            continue
        cpu, node, socket = (int(v) for v in line.split(","))
        node_of[cpu], socket_of[cpu] = node, socket
    # NUMA node of each GPU: the field after the CPU affinity
    gpu_node = []
    cpus = re.compile(r"^\d+(-\d+)?(,\d+(-\d+)?)*$")
    for line in run(["nvidia-smi", "topo", "-m"]).splitlines():
        if not re.match(r"^GPU\d+\t", line):
            continue
        fields = [f.strip() for f in line.split("\t")]
        for i in range(1, len(fields) - 1):
            if cpus.match(fields[i]) and fields[i + 1].isdigit():
                gpu_node.append(int(fields[i + 1]))
                break
        else:
            return
    if not gpu_node or not node_of:
        return
    # socket of a NUMA node
    socket_of_node = {node_of[c]: socket_of[c] for c in node_of}
    ranks = range(min(processes, len(gpu_node)))
    nodes = [{node_of[c] for c in range(r * cores, (r + 1) * cores)
              if c in node_of} for r in ranks]
    sockets = [{socket_of_node[n] for n in nodes[r]} for r in ranks]
    free = list(range(len(gpu_node)))
    chosen = {}
    # GPUs on the NUMA node of the cores of the rank, then on their socket,
    # then any free GPU
    for close in (lambda r, g: gpu_node[g] in nodes[r],
                  lambda r, g: socket_of_node.get(gpu_node[g]) in sockets[r],
                  lambda r, g: True):
        for r in ranks:
            if r not in chosen:
                g = next((g for g in free if close(r, g)), None)
                if g is not None:
                    free.remove(g)
                    chosen[r] = g
    order = [chosen[r] for r in ranks]
    print(",".join(str(g) for g in order + free))


if __name__ == "__main__":
    main()
