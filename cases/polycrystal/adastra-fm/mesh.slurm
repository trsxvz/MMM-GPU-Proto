#!/bin/bash
# Mesh of the polycrystal of FM_GRAINS grains on a node of the MI300A
# partition (the memory of a workstation is too small for gmsh and
# clean_mesh.py above 6,400 grains): the geometry (<grains>crystals.geo) and
# the orientations (vectors_<grains>crystals.txt), made by make-mesh.sh up to
# gmsh, are in <scratch>/meshes/series/v<grains>; gmsh (msh 2.2, mesh size of
# the .geo), then clean_mesh.py (eps 0.3) into
# <scratch>/meshes/clean/v<grains>/<grains>crystals.msh, as make-mesh.sh.
#SBATCH --job-name=fm-mesh
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
s=$W/meshes/series/v$n
c=$W/meshes/clean/v$n
gmsh=${FM_GMSH:-/opt/software/gaia/prod/6.0.0/gmsh-4.13.1-gcc-14.2.1-ftsu/bin/gmsh}
mkdir -p "$c"
log "start, $n grains, $gmsh"
cd "$s" || exit 1
start=$(date +%s)
/usr/bin/time -v "$gmsh" -3 -nt 96 -format msh22 "${n}crystals.geo" \
  -o "${n}crystals-raw.msh" > gmsh.log 2>&1
if grep -qE '^Error|No elements in volume' gmsh.log; then
  log "gmsh errors (see $s/gmsh.log)"
  exit 1
fi
log "gmsh: $(( $(date +%s) - start )) s, $(grep -E 'Maximum resident' gmsh.log | sed 's/^ *//'), $(ls -la "${n}crystals-raw.msh" | awk '{print $5}') bytes"
start=$(date +%s)
/usr/bin/time -v python3 -I "$SLURM_SUBMIT_DIR/fm-chain/clean_mesh.py" \
  "${n}crystals-raw.msh" 0.3 "$c/${n}crystals.msh" > clean.log 2>&1
log "clean: status $?, $(( $(date +%s) - start )) s, $(grep -E 'Maximum resident' clean.log | sed 's/^ *//'); $(grep -v '^\s' clean.log | tail -3 | tr '\n' ' ' | cut -c1-300)"
log end
