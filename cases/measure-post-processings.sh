#!/bin/bash
# Cost of the post-processings (work package L8): polycrystal of
# mm-opera-hpc, 10 steps, all the kernels, Newton, on the GPU and on 8 CPU
# processes, without then with the post-processings of the case (ParaView
# export of the results and averages of the stresses at each step),
# computed on the host then by kernels
# (MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels). The variable
# `executable` gives another executable than the one of the build
# directory (absolute path). The profiling sections of mfem-mgis separate
# the copy of the materials to the host (NLEP::PP::synchronizeOnHost) from
# the execution of the post-processings (NLEP::PP::execute). The ParaView
# files are written in the build directory of the case, then deleted.
#
# usage: cases/measure-post-processings.sh [results directory]
work_package=L8
source "$(dirname "$0")/measure-common.sh" polycrystal "$@"
batch=MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
kernels=(MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
  MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
  MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
cuda=$gp/generic-parallel-cuda/src/libGenericParallelBehaviour-cuda.so
options=""

# waits at most 10 minutes for the GPU to be free of other processes
wait_for_gpu() {
  local i
  for ((i = 0; i != 120; ++i)); do
    if [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
      return 0
    fi
    sleep 5
  done
  echo "   GPU shared with another process" | tee -a "$out/summary.txt"
}

for device in cuda ""; do
  if [ -z "$device" ]; then
    np=8
    library=$generic
    where=np8
  else
    np=1
    library=$cuda
    where=np1-gpu
  fi
  post_processing=--disable-post-processing
  reference=""
  [ -n "$device" ] && wait_for_gpu
  run "polycrystal-$where-without-post-processings" 1 $np "$library" "$batch" \
    "${kernels[@]}"
  post_processing=""
  reference="polycrystal-$where-without-post-processings"
  [ -n "$device" ] && wait_for_gpu
  run "polycrystal-$where-post-processings" 1 $np "$library" "$batch" \
    "${kernels[@]}"
  [ -n "$device" ] && wait_for_gpu
  run "polycrystal-$where-post-processings-kernels" 1 $np "$library" \
    "$batch" "${kernels[@]}" MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels
  rm -rf "$B"/OutputFile-Uniaxial-polycrystal* "$B"/PolycrystalParaviewOutput
done
python3 "$P/cases/summary-profiling.py" "$out"/polycrystal-*.txt > "$out/profiling.txt" 2>&1
echo "end: $(date '+%F %T')" | tee -a "$out/summary.txt"
