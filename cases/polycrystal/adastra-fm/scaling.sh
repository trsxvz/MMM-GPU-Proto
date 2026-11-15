#!/bin/bash
# Chains of jobs of a scalability study of the polycrystal in the final
# configuration of the FM (one material, NewtonK0 with Eisenstat-Walker 1
# from 0.01, line search, extrapolation, partitioned meshes, 20 steps of
# 1 s): 4 processes per node on the GPUs, 96 on the cores. Run on the login
# node, in the directory holding fm-chain/:
#   FM_RUN=<name> fm-chain/scaling.sh mesh <grains>...
#   FM_RUN=<name> fm-chain/scaling.sh check <grains>:<processes>...
#   FM_RUN=<name> fm-chain/scaling.sh gpu <grains> <nodes>...
#   FM_RUN=<name> fm-chain/scaling.sh cpu <grains>:<nodes>[:keep]...
# - mesh: geometries of the seeds FM_SEEDS (default 1 2) then mesh of the
#   first one gmsh accepts, for each number of grains. The orientations
#   must be in <scratch>/meshes/series/v<grains>.
# - check: partitions already made, read back and marked as checked.
# - gpu: parts of one RVE for all its numbers of nodes in one job, then one
#   job of FM_GPU_REPS computations (default 1) per number of nodes.
# - cpu: the points one after the other: parts, computation, and the parts
#   removed by the next job unless the point ends with `:keep` (the parts of
#   the CPU are many files: only one set at a time weighs on the quota).
# A computation only starts after its parts were written and checked; a
# point which fails does not stop the next ones of the chain.
# FM_AFTER: dependency of the first job of each chain (for instance
# afterok:<mesh job>). FM_FILES_LIMIT: files of the group on the scratch
# that the parts must not exceed (default 1490000). FM_TIME_CPU,
# FM_TIMEOUT_CPU: time limits of a CPU job and of its computation in s
# (default 02:00:00 and 6600). FM_TIMEOUT_GPU: time limit of a GPU
# computation in s (default 1500). FM_EXECUTABLE: executable of both sides.
# Prints the last job of each chain; the jobs are listed in
# <results>/jobs.txt.
set -eu
W=/lus/scratch/RES1/genden15/tchenaille
R=$W/fm/${FM_RUN:?FM_RUN required}
mkdir -p "$R"
solver="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01"
common="ALL,FM_RUN=$FM_RUN"
partitioning="$common,FM_PARTITION_WORKERS=16,FM_FILES_LIMIT=${FM_FILES_LIMIT:-1490000}"
computation="$common,FM_PROTO=$W/proto/7eb41368"
computation+=",executable=${FM_EXECUTABLE:-$W/perfimpl4/bin/uniaxial-polycrystal-umpire}"
computation+=",extra_options=--single-material --line-search Backtracking --prediction LinearExtrapolation --mesh_mode Restart"
computation+=",FM_SERIES_REFINEMENT=0,duration=20,FM_NSTEPS=20"
computation+=",FM_SOLVER_GPU=$solver,FM_SOLVER_CPU=$solver"
# usage: submit dependency options...
# job cancelled if its dependency can no longer be satisfied
submit() {
  local after=$1
  shift
  sbatch --parsable -A genden15 --constraint=MI300 --kill-on-invalid-dep=yes \
    ${after:+--dependency=$after} "$@"
}
note() { echo "$(date '+%F %T') $*" >> "$R/jobs.txt"; }
parts() { echo "$W/meshes/part/v$1-np$2"; }
command=${1:?mesh, check, gpu or cpu}
shift
case $command in
  mesh)
    seeds=${FM_SEEDS:-1 2}
    for n; do
      g=$(submit "${FM_AFTER:-}" --job-name="fm-scal-geo-g$n" \
        --export="$common,FM_GRAINS=$n,FM_SEEDS=$seeds" fm-chain/geo.slurm)
      m=$(submit "afterok:$g" --job-name="fm-scal-mesh-g$n" \
        --export="$common,FM_GRAINS=$n,FM_SEEDS=$seeds,FM_GMSH=$W/gmsh-pip/gmsh" \
        fm-chain/mesh.slurm)
      note "mesh of $n grains (seeds $seeds): geometries $g, mesh $m"
      echo "$m"
    done
    ;;
  check)
    c=$(submit "${FM_AFTER:-}" --nodes=1 --job-name=fm-scal-check \
      --export="$common,FM_CHECK=$*" fm-chain/partition.slurm)
    note "check of the partitions $*: $c"
    echo "$c"
    ;;
  gpu)
    n=$1
    shift
    partitions=
    for nodes; do
      partitions+="$n:$((4 * nodes)) "
    done
    p=$(submit "${FM_AFTER:-}" --nodes=1 --job-name="fm-scal-part-g$n-gpu" \
      --export="$partitioning,FM_PARTITIONS=$partitions" fm-chain/partition.slurm)
    note "GPU, $n grains: partitions $partitions$p"
    for nodes; do
      d=$(parts "$n" $((4 * nodes)))
      j=$(submit "afterok:$p" --nodes="$nodes" --time=01:30:00 \
        --job-name="fm-scal-g$n-n$nodes-gpu" \
        --export="$computation,MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels,FM_REQUIRE=$d/checked,FM_GRAINS=$n,FM_MESH_GPU=$d/mesh.,FM_GPU_REPS=${FM_GPU_REPS:-1},FM_CPU_REPS=0,FM_TIMEOUT_GPU=${FM_TIMEOUT_GPU:-1500}" \
        fm-chain/series.slurm)
      note "GPU, $n grains on $nodes nodes x${FM_GPU_REPS:-1}: $j"
    done
    echo "$p"
    ;;
  cpu)
    after=${FM_AFTER:-}
    remove=
    for point; do
      IFS=: read -r n nodes keep <<< "$point"
      np=$((96 * nodes))
      d=$(parts "$n" "$np")
      p=$(submit "$after" --nodes=1 --job-name="fm-scal-part-g$n-n$nodes-cpu" \
        --export="$partitioning,FM_PARTITIONS=$n:$np,FM_REMOVE_PARTITIONS=$remove" \
        fm-chain/partition.slurm)
      # setting of the Slingshot network for the large numbers of ranks
      environment=
      [ "$np" -le 1000 ] || environment=FI_CXI_RX_MATCH_MODE=hybrid
      c=$(submit "afterok:$p" --nodes="$nodes" --time="${FM_TIME_CPU:-02:00:00}" \
        --job-name="fm-scal-g$n-n$nodes-cpu" \
        --export="$computation,FM_REQUIRE=$d/checked,FM_GRAINS=$n,FM_MESH_CPU=$d/mesh.,FM_GPU_REPS=0,FM_CPU_REPS=1,FM_TIMEOUT_CPU=${FM_TIMEOUT_CPU:-6600},FM_CPU_ENVIRONMENT=$environment" \
        fm-chain/series.slurm)
      note "CPU, $n grains on $nodes nodes${keep:+ (parts kept)}: partition $p${remove:+ (removes $remove)}, computation $c"
      after=afterany:$c
      remove=
      [ "$keep" = keep ] || remove=$n:$np
    done
    if [ -n "$remove" ]; then
      c=$(submit "$after" --nodes=1 --job-name=fm-scal-remove \
        --export="$common,FM_REMOVE_PARTITIONS=$remove" fm-chain/partition.slurm)
      note "CPU: removal of $remove $c"
    fi
    echo "$c"
    ;;
  *)
    echo "unknown command $command" >&2
    exit 1
    ;;
esac
