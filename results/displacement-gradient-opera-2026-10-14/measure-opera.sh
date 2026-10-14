#!/bin/bash
# Cermet and polycrystal of mm-opera-hpc, with the deformation gradient (F)
# and with the displacement gradient (H), in double and in float, two
# processes on the cores 0 and 1. Builds in build/h-interface (TFEL, MGIS,
# mfem-mgis and mm-opera-hpc installed in install/h-interface, MFEM of the
# environments mmm-double-cpu and mmm-float-cpu). Results (.res: time, F11,
# F22, F33, then sigma11, sigma22, sigma33) and logs in the directory given
# as argument.
# usage: measure-opera.sh <output directory>
set -u
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
O=$(realpath -m "${1:?output directory}")
mkdir -p "$O"
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
T=$P/install/h-interface/tfel
run() {  # precision, case, gradient, name, options...
  local p=$1 c=$2 g=$3 n=$4
  shift 4
  local E=$P/deps/spack/mmm-$p-cpu
  local MPI
  MPI=$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$E" location -i mpi)
  local M=$P/install/h-interface/mgis
  [ "$p" = float ] && M=$P/install/h-interface/mgis-float
  local B=$P/build/h-interface/mm-opera-hpc-$p/$c
  local exe=cermet
  [ "$c" = polycrystal ] && exe=uniaxial-polycrystal
  local opt=""
  [ "$g" = H ] && opt=--displacement-gradient
  local f="$O/$p-$c-$n-$g"
  (cd "$B" && PATH=$MPI/bin:$PATH LD_LIBRARY_PATH=$T/lib:$M/lib \
    OMPI_MCA_btl=self,tcp \
    taskset -c 0,1 mpirun -n 2 --bind-to none ./$exe "$@" $opt \
      --macroscopic-stress-output-file "$f.res" > "$f.log" 2>&1)
  echo "$p $c $n $g exit $? $(tail -1 "$f.res" 2>/dev/null)"
}
for p in double float; do
  for g in F H; do
    run $p cermet $g 2x0.2 --duration 0.4 --nstep 2 --postprocessing 0
    run $p cermet $g 4x0.1 --duration 0.4 --nstep 4 --postprocessing 0
    for d in 0.1 0.4 1; do
      run $p polycrystal $g 1x$d --disable-post-processing --nstep 1 \
        --duration $d
    done
  done
done
