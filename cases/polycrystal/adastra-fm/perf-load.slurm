#!/bin/bash
# Load per GPU: the final computation of the polycrystal with more unknowns
# per GPU than the series (JOURNAL.md, 2026-11-15). Polycrystal of
# FM_GRAINS grains (default 1600: 53.0 M unknowns after one uniform
# refinement, 6.6 M per GPU on 2 nodes), mesh read by each process, FM_NSTEPS
# steps of 1 s (default 8), with the executable FM_EXECUTABLE (default: the
# one of FM_FIX built with MFEM with Umpire) and, if FM_KERNELS=1 (default),
# the mean stresses computed by kernels. The memory of the first node is
# read every 2 s (<results>/memory-samples.txt: free memory of each NUMA
# node, memory used on each GPU, in GiB; the memory of the GPU of a MI300A
# is counted as GTT). Reference: 8 steps in 193.5 s, 4.36 s in
# NS::computeResidual and 75 to 92 GiB per GPU with the executable of
# 14/11/2026, whose MFEM has no pool of memory (job 5989504).
# The job fails, so that a job submitted with --dependency=afterok does not
# start, if the computation fails or does not end its steps; otherwise the
# file <results>/passed is written (FM_REQUIRE of the next jobs, see
# common.sh).
#SBATCH --job-name=fm-perf-load
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:-1600}
out=$R/perf-load-$SLURM_JOB_ID
mkdir -p "$out"
(
  while :; do
    printf '%s' "$(date +%s.%N)"
    for i in 0 1 2 3; do
      awk -v i=$i '$3 == "MemFree:" {printf " n%d.free:%.2f", i, $4 / 1048576}' \
        /sys/devices/system/node/node$i/meminfo
    done
    for c in /sys/class/drm/card*/device; do
      [ -r "$c/mem_info_gtt_used" ] && printf ' %s.gtt:%.2f' \
        "$(basename "$(dirname "$c")")" \
        "$(awk '{print $1 / 2^30}' "$c/mem_info_gtt_used")"
    done
    echo
    sleep 2
  done
) > "$out/memory-samples.txt" 2>&1 &
sampler=$!
export mesh=$W/meshes/clean/v$n/${n}crystals.msh
export vectors=$MESHES/v$n/vectors_${n}crystals.txt
export refinement=${FM_REFINEMENT:-1} nsteps=${FM_NSTEPS:-8}
export duration=$nsteps modes=gpu run_timeout_gpu=${FM_TIMEOUT_GPU:-1500}
export solver_gpu="--non-linear-solver NewtonK0 --adaptive-linear-solver-tolerance EisenstatWalker1 --adaptive-linear-solver-initial-tolerance 0.01"
export executable=${FM_EXECUTABLE:-${FM_FIX:?}/bin/uniaxial-polycrystal-umpire}
export extra_options="--single-material --line-search Backtracking --prediction LinearExtrapolation"
[ "${FM_KERNELS:-1}" = 1 ] && export MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels
log "start, $n grains, refinement $refinement, $SLURM_JOB_NUM_NODES nodes," \
  "$nsteps steps, executable $executable, kernels of the mean stresses ${FM_KERNELS:-1}"
label=g${n}c solver_label=load "$P/cases/measure-polycrystal-grains.sh" "$out" \
  > "$out/g${n}c.out" 2>&1
kill "$sampler"
line=$(tail -1 "$out/series.txt")
log "$line"
name=g${n}c-r$refinement-gpu-n$SLURM_JOB_NUM_NODES-load
log "sections: $(grep -E -- '--> (NLEP::solve|FPS::integrate|NS::computeResidual|PRED::residual|NS::getJacobian|MFEM::Mult\(r,c\)|computeMacroscopicCauchyStress) ' "$out/$name.txt" | awk -F'|' '{gsub(/ +/, " "); print $(NF - 7), $(NF - 6), "calls, min", $(NF - 5), "mean", $(NF - 4), "max", $(NF - 3)}' | tr '\n' ';')"
log "largest memory used on a GPU of the first node: $(tr ' ' '\n' < "$out/memory-samples.txt" | sed -n 's/.*\.gtt://p' | sort -g | tail -1) GiB;" \
  "smallest free memory of a NUMA node: $(tr ' ' '\n' < "$out/memory-samples.txt" | sed -n 's/.*\.free://p' | sort -g | head -1) GiB"
if [ "$(awk '{print $5, $9}' <<< "$line")" = "0 $nsteps" ]; then
  touch "$out/passed"
  log end
else
  log "end, FAILURE (the next jobs of the chain do not start)"
  exit 1
fi
