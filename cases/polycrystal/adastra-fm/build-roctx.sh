#!/bin/bash
# Separate build of MGIS, of mfem-mgis and of the polycrystal of
# mm-opera-hpc from the complete sources of the three repositories (archive
# FM_FIX/src.tgz holding mgis/, mfem-mgis/ and mm-opera-hpc/), for the
# profiling sections declared as ranges of roctx to the profilers of ROCm.
# MGIS is built once, with its test of the profiling; mfem-mgis and the
# polycrystal are built in the variants FM_BUILD_VARIANTS (default
# "off on"):
# - off: without the option enable-roctx of mfem-mgis (the default);
# - on: with enable-roctx (roctx of rocprofiler-sdk);
# - legacy: with enable-roctx and the roctx of roctracer (libroctx64).
# The sources, the builds and the installations are in FM_BUILD_DIR, by
# default a directory of the local disk of the node, which is lost at the
# end of the job: the quota of files of the project on /lus/scratch was
# nearly reached on 14/11/2026 (a build tree holds thousands of files). Only
# the executables (FM_FIX/bin/uniaxial-polycrystal-<variant>) and the
# library of MGIS which they load (FM_FIX/lib, first in their run path) are
# kept in FM_FIX. The sources and the builds of the prototype on /lus/work
# are not modified.
# To run on a node of the MI300A partition (fm-chain/roctx-tests.slurm), as
# the build of the prototype: CMake detects the architecture of the GPUs
# (GPU_TARGETS gfx942) and -march=native is the one of the nodes.
set -eu
D=/lus/work/RES1/genden15/tchenaille/mmm
F=${FM_FIX:?}
export PROTO_DEPLOY=$D PROTO_MACHINE=adastra-mi300
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
# as build-mgis.sh and build-mfem-mgis.sh, with the installations of the
# machine profile
MFEM=$PROTO_MFEM
MPI_CC=${PROTO_MPI_CC:-cc}
MPI_CXX=${PROTO_MPI_CXX:-CC}
J=${JOBS:-48}
B=${FM_BUILD_DIR:-${TMPDIR:-/tmp}/fix-roctx-${SLURM_JOB_ID:-$$}}
echo "MFEM $MFEM, MPI $MPI_CC $MPI_CXX, compilers $PROTO_CC $PROTO_CXX"
echo "start $(date +%T) on $(hostname), build in $B"
rm -rf "$B" "$F/bin" "$F/lib"
mkdir -p "$B/src" "$F/bin" "$F/lib"
tar xzf "$F/src.tgz" -C "$B/src"
cmake -S "$B/src/mgis" -B "$B/build/mgis" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$PROTO_CC" -DCMAKE_CXX_COMPILER="$PROTO_CXX" \
  -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL" \
  -Denable-testing=ON -Denable-mgis-function=ON "${PROTO_MGIS_OPTIONS[@]}" \
  -Denable-hdf5-automatic-support=OFF \
  -DCMAKE_INSTALL_PREFIX="$B/install/mgis" > "$F/cmake-mgis.log"
cmake --build "$B/build/mgis" -j "$J"
cmake --build "$B/build/mgis" -j "$J" --target ProfilingTest
cmake --install "$B/build/mgis" > /dev/null
cp "$B/install/mgis/lib/libMFrontGenericInterface.so" "$F/lib/"
echo "mgis built $(date +%T)"
common=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$PROTO_CC"
  -DCMAKE_CXX_COMPILER="$PROTO_CXX" -DMPI_C_COMPILER="$MPI_CC"
  -DMPI_CXX_COMPILER="$MPI_CXX" "${PROTO_CMAKE_OPTIONS[@]}" -DMFEM_DIR="$MFEM")
for v in ${FM_BUILD_VARIANTS:-off on}; do
  case $v in
    off) roctx=(-Denable-roctx=OFF) ;;
    on) roctx=(-Denable-roctx=ON) ;;
    legacy) roctx=(-Denable-roctx=ON
      -DROCTX_LIBRARY="$ROCM_PATH/lib/libroctx64.so") ;;
    *)
      echo "unknown variant $v" >&2
      exit 1
      ;;
  esac
  cmake -S "$B/src/mfem-mgis" -B "$B/build/mfem-mgis-$v" "${common[@]}" \
    "${roctx[@]}" \
    -DMFrontGenericInterface_DIR="$B/install/mgis/share/mgis/cmake" \
    -DTFEL_DIR="$PROTO_TFEL_INSTALL/share/tfel/cmake" \
    -DCMAKE_INSTALL_PREFIX="$B/install/mfem-mgis-$v" \
    > "$F/cmake-mfem-mgis-$v.log"
  if ! grep -q '^GPU_TARGETS:STRING=gfx942' \
    "$B/build/mfem-mgis-$v/CMakeCache.txt"; then
    echo "the architecture of the GPUs (gfx942) was not detected" >&2
    exit 1
  fi
  grep -E 'mfem-mgis roctx|Found ROCTX' "$F/cmake-mfem-mgis-$v.log" || true
  cmake --build "$B/build/mfem-mgis-$v" -j "$J"
  cmake --install "$B/build/mfem-mgis-$v" > /dev/null
  # the executable loads the library of MGIS kept in FM_FIX/lib
  cmake -S "$B/src/mm-opera-hpc" -B "$B/build/mm-opera-hpc-$v" \
    "${common[@]}" -DCMAKE_BUILD_RPATH="$F/lib" \
    -DMFEMMGIS_DIR="$B/install/mfem-mgis-$v/share/mfem-mgis/cmake" \
    > "$F/cmake-opera-$v.log"
  cmake --build "$B/build/mm-opera-hpc-$v" --target uniaxial-polycrystal -j 8
  cp "$B/build/mm-opera-hpc-$v/polycrystal/uniaxial-polycrystal" \
    "$F/bin/uniaxial-polycrystal-$v"
  echo "$v built $(date +%T)"
done
echo "end $(date +%T)"
