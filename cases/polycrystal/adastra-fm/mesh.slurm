#!/bin/bash
# Mesh of the polycrystal of FM_GRAINS grains on a node of the MI300A
# partition (the memory of a workstation is too small for gmsh and
# clean_mesh.py above 6,400 grains): the geometry (<grains>crystals.geo) and
# the orientations (vectors_<grains>crystals.txt), made by make-mesh.sh up to
# gmsh, are in <scratch>/meshes/series/v<grains>; gmsh (msh 2.2, mesh size of
# the .geo), then clean_mesh.py (eps 0.3) into
# <scratch>/meshes/clean/v<grains>/<grains>crystals.msh, as make-mesh.sh.
# The plane surfaces made non-planar by the merging of a tiny face by Merope
# (refused by OpenCASCADE) are first defined as surfaces filling their loop
# (fix_plane_surfaces.py, from the distances of their vertices to their
# plane).
# FM_SEEDS: the geometry is one of geo.slurm
# (<scratch>/meshes/geo/v<grains>-seed<seed>): the seeds are tried in this
# order until gmsh meshes one without error, and the seed kept is written in
# the file `seed` next to the geometry.
#SBATCH --job-name=fm-mesh
#SBATCH --constraint=MI300
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
source "$SLURM_SUBMIT_DIR/fm-chain/common.sh"
n=${FM_GRAINS:?}
s=$W/meshes/series/v$n
c=$W/meshes/clean/v$n
gmsh=${FM_GMSH:-/opt/software/gaia/prod/6.0.0/gmsh-4.13.1-gcc-14.2.1-ftsu/bin/gmsh}
mkdir -p "$c" "$s"
log "start, $n grains, $gmsh"
cd "$s" || exit 1
geo=${n}crystals.geo
# usage: mesh_geometry
# fix_plane_surfaces.py then gmsh on the geometry of the current directory;
# fails if gmsh reports an error or a volume without elements
mesh_geometry() {
  local start
  # plane surfaces made non-planar by Merope, defined as surfaces (distance
  # of the vertices to their plane above 1e-8; OpenCASCADE refuses 1e-7)
  start=$(date +%s)
  python3 -I "$SLURM_SUBMIT_DIR/fm-chain/fix_plane_surfaces.py" "$geo" \
    "${n}crystals-fixed.geo" > fix.log 2>&1
  log "fix: $(head -1 fix.log), $(( $(date +%s) - start )) s"
  start=$(date +%s)
  /usr/bin/time -v "$gmsh" -3 -nt 96 -format msh22 "${n}crystals-fixed.geo" \
    -o "${n}crystals-raw.msh" > gmsh.log 2>&1
  if grep -qE '^Error|No elements in volume' gmsh.log; then
    log "gmsh errors (see $s/gmsh.log)"
    return 1
  fi
  log "gmsh: $(( $(date +%s) - start )) s, $(grep -E 'Maximum resident' gmsh.log | sed 's/^ *//'), $(ls -la "${n}crystals-raw.msh" | awk '{print $5}') bytes"
}
if [ -z "${FM_SEEDS:-}" ]; then
  mesh_geometry || exit 1
else
  rm -f seed
  for seed in $FM_SEEDS; do
    log "geometry of seed $seed"
    cp "$W/meshes/geo/v$n-seed$seed/$geo" . || continue
    if mesh_geometry; then
      echo "$seed" > seed
      break
    fi
  done
  [ -e seed ] || exit 1
fi
start=$(date +%s)
/usr/bin/time -v python3 -I "$SLURM_SUBMIT_DIR/fm-chain/clean_mesh.py" \
  "${n}crystals-raw.msh" 0.3 "$c/${n}crystals.msh" > clean.log 2>&1
status=$?
log "clean: status $status, $(( $(date +%s) - start )) s, $(grep -E 'Maximum resident' clean.log | sed 's/^ *//'); $(grep -v '^\s' clean.log | tail -3 | tr '\n' ' ' | cut -c1-300)"
# the jobs of a chain only go on after a cleaned mesh
[ "$status" = 0 ] || exit 1
log end
