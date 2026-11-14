#!/bin/bash
# Meshes in binary Gmsh 2.2 (meshes/msh22_binary.py) against the same meshes
# in ASCII: on the CPU (FM_CPU_GRAINS grains, default 100), the outputs must
# be identical; on the GPU (FM_GPU_GRAINS, default 400), time of the reading
# of the mesh (Mesh::LoadMesh). FM_NSTEPS steps of 1 s (default 3), the
# options of the FM (variables executable and extra_options).
#SBATCH --job-name=fm-mesh-bin
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
out=$R/mesh-bin-$SLURM_JOB_ID
mkdir -p "$out"
log "start, scripts $P"
export refinement=1 nsteps=${FM_NSTEPS:-3} duration=${FM_NSTEPS:-3} solver_label=best
export solver_gpu=${FM_SOLVER_GPU:?} solver_cpu=${FM_SOLVER_CPU:?}
export run_timeout_gpu=2400 run_timeout_cpu=2400
load() {
  awk -F'|' '$3 ~ /--> Mesh::LoadMesh/ {gsub(/ /, "", $7); print $7; exit}' "$1"
}
for mode in cpu gpu; do
  if [ $mode = cpu ]; then n=${FM_CPU_GRAINS:-100}; else n=${FM_GPU_GRAINS:-400}; fi
  export vectors=$MESHES/v$n/vectors_${n}crystals.txt
  for f in ascii bin; do
    if [ $f = ascii ]; then m=${n}crystals.msh; else m=${n}crystals-bin.msh; fi
    mesh=$W/meshes/clean/v$n/$m label=g${n}c-$f modes=$mode \
      "$P/cases/measure-polycrystal-grains.sh" "$out"
    log "$mode $f: $(tail -1 "$out/series.txt"); LoadMesh $(load "$out/g${n}c-$f-r1-$mode-n1-best.txt") s"
  done
  a=$out/g${n}c-ascii-r1-$mode-n1-best b=$out/g${n}c-bin-r1-$mode-n1-best
  if cmp -s "$a.res" "$b.res"; then
    log "$mode: .res of the ASCII and binary meshes identical"
  else
    log "$mode: .res of the ASCII and binary meshes differ"
  fi
done
log end
