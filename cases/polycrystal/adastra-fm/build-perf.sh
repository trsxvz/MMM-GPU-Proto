#!/bin/bash
# Separate build of MFEM, of MGIS, of mfem-mgis and of the polycrystal of
# mm-opera-hpc for the performance of the polycrystal on the GPUs
# (JOURNAL.md, 2026-11-15), from the archive FM_FIX/src.tgz holding the
# complete sources of the three repositories (mgis/, mfem-mgis/ and
# mm-opera-hpc/) and, optionally, the machine profile machine.sh
# (deps/machines/adastra-mi300/machine.sh), which then replaces the one of
# the deployed prototype (options of MFEM not deployed yet). MGIS is built
# once; mfem-mgis (with the option enable-roctx, which costs nothing without
# profiler) and the polycrystal are built in the variants FM_BUILD_VARIANTS
# (default "umpire"), which differ by their MFEM:
# - umpire: the MFEM built by this script with the options of the machine
#   profile (PROTO_MFEM_OPTIONS: with Umpire), installed in
#   FM_FIX/install/mfem and built only if it is not there;
# - deployed: the MFEM of the deployed prototype (built without Umpire on
#   14/11/2026).
# The sources and the builds are in FM_BUILD_DIR, by default a directory of
# the local disk of the node, lost at the end of the job (quota of files of
# the project on /lus/scratch); a second call in the same directory only
# builds what changed. Kept in FM_FIX: the installation of MFEM, the
# executables (bin/uniaxial-polycrystal-<variant>) and the library of MGIS
# which they load (lib/). The run path of the executables only names
# FM_FIX/lib and the installations of TFEL and of MFEM, not the build
# directory. FM_FIX is the directory of one version of the sources: the
# checksum of the archive is kept in FM_FIX/sources.sha256, and another
# archive is refused, since the executables already built load the library
# of MGIS of FM_FIX. FM_FIX/install may be a link to the installation of
# MFEM of another directory. Nothing is written on /lus/work.
# To run on a node of the MI300A partition (fm-chain/perf-tests.slurm), as
# the build of the prototype: CMake detects the architecture of the GPUs
# (GPU_TARGETS gfx942) and -march=native is the one of the nodes.
set -eu
D=/lus/work/RES1/genden15/tchenaille/mmm
F=${FM_FIX:?}
export PROTO_DEPLOY=$D PROTO_MACHINE=adastra-mi300
# the errors of the machine profile are kept (libraries not found)
set +u
source "$D/MMM-GPU-Proto/scripts/env.sh" > /dev/null
set -u
J=${JOBS:-48}
B=${FM_BUILD_DIR:-${TMPDIR:-/tmp}/fix-perf-${SLURM_JOB_ID:-$$}}
echo "start $(date +%T) on $(hostname), build in $B"
mkdir -p "$B/src" "$F/bin" "$F/lib"
sources=$(sha256sum "$F/src.tgz" | cut -d' ' -f1)
if [ -f "$F/sources.sha256" ] && [ "$(cat "$F/sources.sha256")" != "$sources" ]; then
  echo "$F holds executables built from other sources" \
    "($(cut -c1-16 "$F/sources.sha256"), archive ${sources:0:16}):" \
    "use another directory" >&2
  exit 1
fi
echo "$sources" > "$F/sources.sha256"
tar xzf "$F/src.tgz" -C "$B/src"
deployed_mfem=$PROTO_MFEM
if [ -f "$B/src/machine.sh" ]; then
  # P is the deployed prototype (installation of TFEL, launchers)
  set +u
  source "$B/src/machine.sh" > /dev/null
  set -u
fi
MPI_CC=${PROTO_MPI_CC:-cc}
MPI_CXX=${PROTO_MPI_CXX:-CC}
echo "MPI $MPI_CC $MPI_CXX, compilers $PROTO_CC $PROTO_CXX"
cmake -S "$B/src/mgis" -B "$B/build/mgis" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$PROTO_CC" -DCMAKE_CXX_COMPILER="$PROTO_CXX" \
  -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL" \
  -Denable-testing=ON -Denable-mgis-function=ON "${PROTO_MGIS_OPTIONS[@]}" \
  -Denable-hdf5-automatic-support=OFF \
  -DCMAKE_INSTALL_PREFIX="$B/install/mgis" > "$F/cmake-mgis.log"
cmake --build "$B/build/mgis" -j "$J"
cmake --install "$B/build/mgis" > /dev/null
cp "$B/install/mgis/lib/libMFrontGenericInterface.so" "$F/lib/"
echo "mgis built $(date +%T)"
# MFEM with the options of the machine profile, as `deploy.sh mfem`, in
# another prefix
build_mfem() {
  local src=$B/mfem-src
  source "$D/MMM-GPU-Proto/scripts/offline/versions.sh"
  local archive=$D/mfem-$PROTO_MFEM_VERSION.tar.gz
  echo "$PROTO_MFEM_SHA256  $archive" | sha256sum -c --quiet -
  rm -rf "$src"
  mkdir -p "$src"
  tar -xzf "$archive" -C "$src" --strip-components=1
  (cd "$src" && make config PREFIX="$1" "${PROTO_MFEM_OPTIONS[@]}") \
    > "$F/config-mfem.log"
  make -C "$src" -j "$J" lib > "$F/build-mfem.log" 2>&1
  make -C "$src" install >> "$F/build-mfem.log" 2>&1
  echo "MFEM $PROTO_MFEM_VERSION built $(date +%T)"
}
for v in ${FM_BUILD_VARIANTS:-umpire}; do
  case $v in
    umpire)
      MFEM=$F/install/mfem
      [ -f "$MFEM/lib/libmfem.so" ] || build_mfem "$MFEM"
      ;;
    deployed) MFEM=$deployed_mfem ;;
    *)
      echo "unknown variant $v" >&2
      exit 1
      ;;
  esac
  echo "$v: MFEM $MFEM," \
    "$(grep -E '^MFEM_USE_UMPIRE|^MFEM_CXXFLAGS' "$MFEM/share/mfem/config.mk" | tr -s ' ' | tr '\n' ';')"
  cmake -S "$B/src/mfem-mgis" -B "$B/build/mfem-mgis-$v" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$PROTO_CC" \
    -DCMAKE_CXX_COMPILER="$PROTO_CXX" -DMPI_C_COMPILER="$MPI_CC" \
    -DMPI_CXX_COMPILER="$MPI_CXX" "${PROTO_CMAKE_OPTIONS[@]}" \
    -DMFEM_DIR="$MFEM" -Denable-roctx=ON \
    -DMFrontGenericInterface_DIR="$B/install/mgis/share/mgis/cmake" \
    -DTFEL_DIR="$PROTO_TFEL_INSTALL/share/tfel/cmake" \
    -DCMAKE_INSTALL_PREFIX="$B/install/mfem-mgis-$v" \
    > "$F/cmake-mfem-mgis-$v.log"
  if ! grep -q '^GPU_TARGETS:STRING=gfx942' \
    "$B/build/mfem-mgis-$v/CMakeCache.txt"; then
    echo "the architecture of the GPUs (gfx942) was not detected" >&2
    exit 1
  fi
  cmake --build "$B/build/mfem-mgis-$v" -j "$J"
  cmake --install "$B/build/mfem-mgis-$v" > /dev/null
  # the executable loads the library of MGIS kept in FM_FIX/lib: its run
  # path is given, instead of the one of the build tree, which names the
  # build directory
  cmake -S "$B/src/mm-opera-hpc" -B "$B/build/mm-opera-hpc-$v" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$PROTO_CC" \
    -DCMAKE_CXX_COMPILER="$PROTO_CXX" -DMPI_C_COMPILER="$MPI_CC" \
    -DMPI_CXX_COMPILER="$MPI_CXX" "${PROTO_CMAKE_OPTIONS[@]}" \
    -DMFEM_DIR="$MFEM" -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
    -DCMAKE_INSTALL_RPATH="$F/lib;$PROTO_TFEL_INSTALL/lib;$MFEM/lib" \
    -DMFEMMGIS_DIR="$B/install/mfem-mgis-$v/share/mfem-mgis/cmake" \
    > "$F/cmake-opera-$v.log"
  cmake --build "$B/build/mm-opera-hpc-$v" --target uniaxial-polycrystal -j 8
  e=$B/build/mm-opera-hpc-$v/polycrystal/uniaxial-polycrystal
  if readelf -d "$e" | grep -E 'RPATH|RUNPATH' | grep -q "$B"; then
    echo "the run path of the executable names the build directory $B" >&2
    exit 1
  fi
  if ldd "$e" | grep "not found" >&2; then
    echo "libraries of the executable not found" >&2
    exit 1
  fi
  cp "$e" "$F/bin/uniaxial-polycrystal-$v"
  echo "$v built $(date +%T)"
done
echo "end $(date +%T)"
