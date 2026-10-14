#!/bin/bash
# Reference times of the current path of mfem-mgis (work package L0): cases
# of the measurement campaign, on the CPU (no device), one OpenMP thread per
# process. To be run on a quiet machine, after build-mfem-mgis.sh.
#
# usage: cases/measure-reference.sh [results directory]
#
# Raw outputs (with the profiling tables) in the results directory (by
# default results/reference-L0-<date>), one per case, and a summary in
# summary.txt: command, status, total duration of the process, summary
# lines of the case.
set -u
P=$(cd "$(dirname "$0")/.." && pwd)
source "$P/scripts/env.sh" > /dev/null
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
  location -i openmpi)
export PATH=$MPI/bin:$PATH
export OMP_NUM_THREADS=1
out=${1:-$P/results/reference-L0-$(date +%Y-%m-%d-%H%M)}
mkdir -p "$out"
B=$P/build
# usage: run name directory command...
run() {
  local name=$1 dir=$2
  shift 2
  echo "== $name: $*" | tee -a "$out/summary.txt"
  local start end status
  start=$(date +%s.%N)
  (cd "$dir" && timeout 3600 "$@") < /dev/null > "$out/$name.txt" 2>&1
  status=$?
  end=$(date +%s.%N)
  echo "   status $status, $(echo "$end - $start" | bc) s" |
    tee -a "$out/summary.txt"
  /usr/bin/grep -E "^dofs|^time" "$out/$name.txt" | tee -a "$out/summary.txt"
}
{
  echo "date: $(date '+%F %T')"
  echo "machine: $(nproc --all) cores, $(free -g | awk '/^Mem/ {print $2}') GB"
  echo "mfem-mgis: $(git -C "$P/mfem-mgis" log -1 --format=%h)"
} | tee -a "$out/summary.txt"

# 1. beam in J2 plasticity (small strains, 3D), Newton of mfem-mgis:
# 15,795 unknowns (r3) on 4 processes, 111,843 (r4) on 8
beam=$B/cases/beam
for ls in "pcg-amg" "pcg-amg -lt 1e-4" "mumps"; do
  name=beam-r3-np4-$(echo "$ls" | tr ' ' '_')
  run "$name" "$beam" mpirun -n 4 ./beam -r 3 -n 10 -mi 3000 -ls $ls \
    --profiling
done
for ls in "pcg-amg -lt 1e-4" "mumps"; do
  name=beam-r4-np8-$(echo "$ls" | tr ' ' '_')
  run "$name" "$beam" mpirun -n 8 ./beam -r 4 -n 10 -mi 3000 -ls $ls \
    --profiling
done

# 2. ssna303 in 3D (finite strain), FGMRES and BoomerAMG of hypre, 10 steps
# until the end of the loading
ex4=$B/mfem-mgis-examples/ex4
run ssna303-3d-hypre-r0-np4 "$ex4" mpirun -n 4 ./ssna303_3d_hypre \
  -ns 10 -et 1 --profiling
run ssna303-3d-hypre-r1-np8 "$ex4" mpirun -n 8 ./ssna303_3d_hypre \
  -ns 10 -et 1 -r 1 --profiling

# 3. polycrystal (5 grains, single crystal behaviour), without
# post-processing
poly=$B/mm-opera-hpc/polycrystal
run polycrystal-r0-np4 "$poly" mpirun -n 4 ./uniaxial-polycrystal \
  --disable-post-processing --nstep 10 --duration 10
run polycrystal-r1-np8 "$poly" mpirun -n 8 ./uniaxial-polycrystal \
  --disable-post-processing --nstep 10 --duration 10 -r 1

# 4. thermomechanics (coupled non linear heat transfer), plate of the
# example ex8, complete mesh, MUMPS for the mechanics
ex8=$B/mfem-mgis-examples/ex8
run thermomechanical-np4 "$ex8" mpirun -n 4 ./rjh_plate \
  -m assemblage_hexa.msh -lsMc MUMPSSolver
echo "results in $out"
