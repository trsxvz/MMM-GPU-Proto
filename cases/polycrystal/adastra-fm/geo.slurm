#!/bin/bash
# Geometries of the polycrystal (Laguerre tessellation of Merope, as
# meshes/make-mesh.sh) of FM_GRAINS grains for the seeds FM_SEEDS (default
# 1 to 8), computed in parallel on a node of the MI300A partition with the
# Merope built in <scratch>/merope (Merope is sequential and needs about
# 0.3 GB per 1000 grains: the geometries of more than 25,000 grains do not
# fit in the memory of a workstation). Each geometry goes to
# <scratch>/meshes/geo/v<grains>-seed<seed>; the number of warnings "An
# entire face should be merged" of each seed is logged (0: the geometry is
# meshed by gmsh as is). With FM_CHECK=1, the geometry of 6,400 grains of
# seed 2 is also made, to compare with the one made on the workstation.
#SBATCH --job-name=fm-geo
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
M=$W/merope
# library path set, not extended: the inherited one may give the libgomp of
# ROCm (LLVM libomp), a second OpenMP runtime next to the libiomp5 of MKL
export LD_LIBRARY_PATH=$M/extra-lib:/opt/software/gaia/prod/6.0.0/intel-oneapi-mkl-2025.2.0-gcc-14.2.1-njwr/mkl/2025.2/lib
unset LD_PRELOAD
export PYTHONPATH=$M/INSTALL-DIR/lib:$M/tools/python
template=$D/MMM-GPU-Proto/mm-opera-hpc/polycrystal/mesh/1000crystals.py
geometry() {
  local grains=$1 seed=$2 d=$W/meshes/geo/v$1-seed$2 L
  mkdir -p "$d"
  cd "$d" || return
  L=$(python3 -c "print(($grains*1000)**(1/3))")
  sed -e "s/L = \[100, 100, 100\]/L = [$L, $L, $L]/" \
    -e "s/nbSpheres = 1000 /nbSpheres = $grains/" -e 's/MeshSize = 1$/MeshSize = 2/' \
    -e "s/randomSeed = 1$/randomSeed = $seed/" \
    -e "s/1000crystals.geo/${grains}crystals.geo/" "$template" > "${grains}crystals.py"
  start=$(date +%s)
  /usr/bin/time -v python3 "${grains}crystals.py" > geo.log 2>&1
  echo "status $?, $(( $(date +%s) - start )) s, $(grep -c 'entire face should be merged' geo.log) warnings, $(grep 'Maximum resident' geo.log | sed 's/^ *//')" > geo.done
}
log "start, $n grains, seeds ${FM_SEEDS:-1 2 3 4 5 6 7 8}"
for s in ${FM_SEEDS:-1 2 3 4 5 6 7 8}; do
  geometry "$n" "$s" &
done
if [ "${FM_CHECK:-0}" = 1 ]; then
  geometry 6400 2 &
fi
wait
for s in ${FM_SEEDS:-1 2 3 4 5 6 7 8}; do
  log "$n grains, seed $s: $(cat "$W/meshes/geo/v$n-seed$s/geo.done")"
done
if [ "${FM_CHECK:-0}" = 1 ]; then
  log "check 6400 grains seed 2: $(cat "$W/meshes/geo/v6400-seed2/geo.done"); md5 $(md5sum < "$W/meshes/geo/v6400-seed2/6400crystals.geo" | cut -c1-12)"
fi
log end
