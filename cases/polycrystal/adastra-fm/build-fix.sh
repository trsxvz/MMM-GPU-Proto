#!/bin/bash
# Separate build, on the scratch, of mfem-mgis with the weights of the
# integration points cached (mean values of the thermodynamic forces) and of
# the polycrystal of mm-opera-hpc with the macroscopic stress timed:
# - A: timed polycrystal against the installed mfem-mgis of the prototype;
# - B: timed polycrystal against the modified mfem-mgis.
# FM_FIX: directory of the build (default <scratch>/fix), holding this
# script and the modified files in patched/; FM_BUILD_VARIANTS: variants
# built (default "A B").
# The sources and the builds of the prototype on /lus/work are not modified.
# To run on a node of the MI300A partition (fm-chain/fix-ab.slurm), as the
# build of the prototype: CMake detects the architecture of the GPUs
# (GPU_TARGETS gfx942) and -march=native is the one of the nodes.
set -eu
D=/lus/work/RES1/genden15/tchenaille/mmm
F=${FM_FIX:-/lus/scratch/RES1/genden15/tchenaille/fix}
export PROTO_DEPLOY=$D PROTO_MACHINE=adastra-mi300
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null 2>&1
set -u
P=$D/MMM-GPU-Proto
# as build-mfem-mgis.sh, with the installations of the machine profile
MFEM=$PROTO_MFEM
MPI_CC=${PROTO_MPI_CC:-cc}
MPI_CXX=${PROTO_MPI_CXX:-CC}
echo "MFEM $MFEM, MPI $MPI_CC $MPI_CXX, compilers $PROTO_CC $PROTO_CXX"
echo "start $(date +%T) on $(hostname)"
rm -rf "$F/src" "$F/build" "$F/install" "$F/bin"
for r in mfem-mgis mm-opera-hpc; do
  rm -rf "$F/src/$r"
  mkdir -p "$F/src"
  cp -r "$P/$r" "$F/src/$r"
  cp -r "$F/patched/$r/." "$F/src/$r/"
done
common=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$PROTO_CC"
  -DCMAKE_CXX_COMPILER="$PROTO_CXX" -DMPI_C_COMPILER="$MPI_CC"
  -DMPI_CXX_COMPILER="$MPI_CXX" "${PROTO_CMAKE_OPTIONS[@]}" -DMFEM_DIR="$MFEM")
cmake -S "$F/src/mfem-mgis" -B "$F/build/mfem-mgis" "${common[@]}" \
  -DMFrontGenericInterface_DIR="$P/install/mgis/share/mgis/cmake" \
  -DTFEL_DIR="$PROTO_TFEL_INSTALL/share/tfel/cmake" \
  -DCMAKE_INSTALL_PREFIX="$F/install/mfem-mgis" > "$F/cmake-mfem-mgis.log"
if ! grep -q '^GPU_TARGETS:STRING=gfx942' "$F/build/mfem-mgis/CMakeCache.txt"; then
  echo "the architecture of the GPUs (gfx942) was not detected" >&2
  exit 1
fi
for v in ${FM_BUILD_VARIANTS:-A B}; do
  if [ $v = A ]; then i=$P/install/mfem-mgis; else i=$F/install/mfem-mgis; fi
  if [ $v = B ]; then
    cmake --build "$F/build/mfem-mgis" -j "${JOBS:-24}"
    cmake --install "$F/build/mfem-mgis" > /dev/null
  fi
  cmake -S "$F/src/mm-opera-hpc" -B "$F/build/mm-opera-hpc-$v" "${common[@]}" \
    -DMFEMMGIS_DIR="$i/share/mfem-mgis/cmake" > "$F/cmake-opera-$v.log"
  cmake --build "$F/build/mm-opera-hpc-$v" --target uniaxial-polycrystal -j 8
  mkdir -p "$F/bin"
  cp "$F/build/mm-opera-hpc-$v/polycrystal/uniaxial-polycrystal" \
    "$F/bin/uniaxial-polycrystal-$v"
  echo "$v built $(date +%T)"
done
echo "end $(date +%T)"
