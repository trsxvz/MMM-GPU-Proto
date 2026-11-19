#!/bin/bash
# Computation of the polycrystal on the GPUs of the allocation with the
# ParaView export of its fields at each time step, for a video: one GPU
# computation of series.slurm (FM_GRAINS, FM_MESH_GPU, ...), the export
# being asked by the variables of the measurement script (paraview_export=yes
# and the options of the fields in extra_options, for instance
# --enable-export-von-Mises-stress --enable-paraview-single-precision).
# Once the computation is over, the files of each time step of each export
# (one directory CycleNNNNNN per step) are gathered in one tar archive,
# CycleNNNNNN.tar, by archive-cycles.slurm: one process per step, spread over
# the nodes of the allocation (archiving the steps one at a time during the
# computation, with the default records of tar, ran at about 7 MB/s on
# Lustre, far behind the computation). The number of files of the group thus
# grows by the files of all the steps during the computation (about 5 per
# step and per 4 GPUs), then falls to one per step. The .pvd files are kept:
# after `tar xf CycleNNNNNN.tar` in the directory of the export, ParaView
# reads the export unchanged.
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

FM_GPU_REPS=1 FM_CPU_REPS=0 bash "$SLURM_SUBMIT_DIR/fm-chain/series.slurm"
mapfile -t exports < <(find "$computation" -mindepth 1 -maxdepth 1 -type d)
if [ ${#exports[@]} -gt 0 ]; then
  ARCHIVE_DIRS="${exports[*]}" ARCHIVE_SCRIPTS=$SLURM_SUBMIT_DIR/fm-chain \
    bash "$SLURM_SUBMIT_DIR/fm-chain/archive-cycles.slurm"
fi
log "export of the fields: $(find "$computation" -name 'Cycle*.tar' | wc -l)" \
  "archives, $(find "$computation" -type d -name 'Cycle*' | wc -l) steps" \
  "left unarchived, $(du -sh "$computation" | cut -f1)"
