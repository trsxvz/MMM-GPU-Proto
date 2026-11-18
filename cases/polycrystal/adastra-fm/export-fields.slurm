#!/bin/bash
# Computation of the polycrystal on the GPUs of the allocation with the
# ParaView export of its fields at each time step, for a video: one GPU
# computation of series.slurm (FM_GRAINS, FM_MESH_GPU, ...), the export
# being asked by the variables of the measurement script (paraview_export=yes
# and the options of the fields in extra_options, for instance
# --enable-export-von-Mises-stress --enable-paraview-single-precision).
# The files of a time step of an export (one directory CycleNNNNNN per step)
# are gathered in one tar archive, CycleNNNNNN.tar, as soon as the step is
# complete (the next one exists), the last step at the end: the number of
# files of the group stays small (one or two steps at a time instead of all
# of them). Each archive is checked (same number of files) before its files
# are removed. The .pvd files are kept: after `tar xf CycleNNNNNN.tar` in
# the directory of the export, ParaView reads the export unchanged.
# usage, from the submission directory of fm-chain:
#   sbatch -A genden15 --nodes=<n> --time=<t> --job-name=fm-export-fields \
#     --export="ALL,FM_RUN=<run>,FM_PROTO=<scripts>,executable=<binary>,\
#   paraview_export=yes,extra_options=<options>,FM_SOLVER_GPU=<solver>,\
#   FM_SOLVER_CPU=<solver>,FM_SERIES_REFINEMENT=0,duration=20,FM_NSTEPS=20,\
#   MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels,FM_GRAINS=<n>,\
#   FM_MESH_GPU=<parts>/mesh.,FM_REQUIRE=<parts>/checked" \
#     fm-chain/export-fields.slurm
#SBATCH --constraint=MI300
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
computation=$R/series-g${FM_GRAINS:?}-$SLURM_JOB_ID/g${FM_GRAINS}c-rep1-r0-gpu-n$SLURM_JOB_NUM_NODES-best

# archive_cycles <directory of the computation> [final]: archives the
# complete steps of each export, and also the last ones with `final`
archive_cycles() {
  local export_dir cycles cycle name n i
  for export_dir in "$1"/*/; do
    mapfile -t cycles < <(find "$export_dir" -maxdepth 1 -type d \
                               -name 'Cycle*' | sort)
    n=${#cycles[@]}
    [ -n "${2:-}" ] || n=$((n - 1))
    for ((i = 0; i < n; i++)); do
      cycle=${cycles[i]}
      name=$(basename "$cycle")
      if tar cf "$cycle.tar.part" -C "$export_dir" "$name" &&
         [ "$(tar tf "$cycle.tar.part" | grep -vc '/$')" = \
           "$(find "$cycle" -type f | wc -l)" ]; then
        mv "$cycle.tar.part" "$cycle.tar" && rm -rf "$cycle"
      else
        rm -f "$cycle.tar.part"
        log "archive of $cycle failed, its files are kept"
      fi
    done
  done
}

done_file=$R/.export-fields-$SLURM_JOB_ID-done
(
  while [ ! -e "$done_file" ]; do
    [ -d "$computation" ] && archive_cycles "$computation"
    sleep 20
  done
) &
archiver=$!
FM_GPU_REPS=1 FM_CPU_REPS=0 bash "$SLURM_SUBMIT_DIR/fm-chain/series.slurm"
touch "$done_file"
wait "$archiver"
[ -d "$computation" ] && archive_cycles "$computation" final
rm -f "$done_file"
log "export of the fields: $(find "$computation" -name 'Cycle*.tar' | wc -l)" \
  "archives, $(find "$computation" -type d -name 'Cycle*' | wc -l) steps" \
  "left unarchived, $(du -sh "$computation" | cut -f1)"
