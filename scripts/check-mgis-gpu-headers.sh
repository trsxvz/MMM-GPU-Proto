#!/bin/bash
# Compile and run reproducers/mgis-gpu-headers/check.cxx with each GPU
# programming model: the headers of MGIS compiled by the calling code (field
# operations, memory resources) are checked on the device.
# The compilers and their flags are the ones of the configurations of the
# generic-parallel interface of TFEL (substitutions of build/tfel-<env>).
# usage: scripts/check-mgis-gpu-headers.sh [cuda] [hip] [sycl] [kokkos] [openmp]
set -u
PROTO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$PROTO/scripts/env.sh" > /dev/null
SRC="$PROTO/reproducers/mgis-gpu-headers/check.cxx"
OUT="$PROTO/scratch/mgis-gpu-headers"
MGIS_LIB="$PROTO/build/mgis/src"
mkdir -p "$OUT"
# substitutions of the TFEL build
subst() {
  grep '^tfel-substitutions:STRING=' "$PROTO_TFEL_BUILD/CMakeCache.txt" |
    tr ';' '\n' | sed -n "s|^$1%||p"
}
INC="-I$PROTO/mgis/include"
LIBS="-L$MGIS_LIB -lMFrontGenericInterface"
export LD_LIBRARY_PATH="$MGIS_LIB:$LD_LIBRARY_PATH"
BACKENDS=${*:-cuda hip sycl kokkos openmp}
status=0
for b in $BACKENDS; do
  exe="$OUT/check-$b"
  case $b in
    cuda)
      cmd="$(subst TFEL_NVCC_COMPILER) -std=c++20 -O2 -arch=native -x cu -ccbin g++-14 -DMGIS_CHECK_CUDA $INC $SRC -o $exe $LIBS" ;;
    hip)
      cmd="$(subst TFEL_HIPCC_COMPILER) $(subst TFEL_HIPCC_FLAGS) -std=c++20 -O2 -DMGIS_CHECK_HIP $INC $SRC -o $exe $LIBS" ;;
    sycl)
      cmd="$(subst TFEL_ACPP_COMPILER) -std=c++20 -O2 -DMGIS_CHECK_SYCL $INC $SRC -o $exe $LIBS" ;;
    kokkos)
      cmd="$(subst TFEL_KOKKOS_CXX_COMPILER) -std=c++20 -O2 $(subst TFEL_KOKKOS_CXX_FLAGS) -DMGIS_CHECK_KOKKOS $INC $SRC -o $exe $(subst TFEL_KOKKOS_LINKER_FLAGS) $LIBS" ;;
    openmp)
      cmd="$(subst TFEL_NVHPC_COMPILER) -std=c++20 -O2 -mp=gpu -DMGIS_CHECK_OPENMP $INC $SRC -o $exe $LIBS" ;;
    *)
      echo "unknown backend '$b'" >&2; status=1; continue ;;
  esac
  if ! $cmd > "$OUT/$b.log" 2>&1; then
    echo "$b: compilation failed (see $OUT/$b.log)"; status=1; continue
  fi
  if "$exe" >> "$OUT/$b.log" 2>&1; then
    echo "$b: success"
  else
    echo "$b: failure (see $OUT/$b.log)"; status=1
  fi
done
exit $status
