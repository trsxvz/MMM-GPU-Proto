#!/usr/bin/env bash
# Construction d'un MFEM compilé avec HIP pour un GPU AMD (gfx90a), pour
# vérifier que mfem-mgis compile et se lie avec lui. La machine n'a pas de
# GPU AMD : rien ne s'exécute sur le GPU, seules la compilation et l'édition
# de liens sont vérifiées.
#
# ROCm vient des paquets Python de ROCm (TheRock), installés sans droits
# d'administrateur dans un environnement virtuel :
#   python3 -m venv $ROCM_VENV
#   $ROCM_VENV/bin/pip install --index-url \
#     https://stable.repo.amd.com/rocm/whl-next/ \
#     "rocm[libraries,devel,device-gfx90a]"
#   $ROCM_VENV/bin/rocm-sdk init
#
# Pièges :
# - un paquet système (libamdhip64-dev d'Ubuntu) installe d'anciens en-têtes
#   de HIP dans /usr/include/hip, que le compilateur de ROCm lit avant les
#   siens : les en-têtes de ROCm sont donnés explicitement (-isystem) à
#   hipcc par HIPCC_COMPILE_FLAGS_APPEND et au compilateur HIP de CMake par
#   CMAKE_HIP_FLAGS ;
# - scripts/env.sh exporte HIP_PLATFORM=nvidia (HIP sur GPU NVIDIA pour
#   TFEL), que hipcc et hipconfig suivent : MFEM est compilé avec
#   HIP_PLATFORM=amd, et mfem-mgis choisit lui-même la plateforme AMD quand
#   MFEM utilise HIP.
#
# hypre est compilé pour le CPU (MFEM avec HIP accepte un hypre sur l'hôte).
# MPI, METIS, SuiteSparse et OpenBLAS sont ceux de l'environnement spack
# mmm-gpu. Résultat : deps/hip-amd/install/{hypre,mfem}.
#
# L'étape mfem-mgis construit mfem-mgis (dossier mfem-mgis du prototype)
# avec ce MFEM, MGIS et TFEL du prototype, dans build/mfem-mgis-hip, avec
# ses tests (cible build-tests), puis lance les tests qui ne demandent pas
# de GPU (ctest -E gpu ; il n'y a pas de GPU AMD).
#
# L'étape examples installe ce mfem-mgis dans install/mfem-mgis-hip et
# construit les exemples de mfem-mgis-examples avec lui (configuration
# installée, lois compilées pour le GPU AMD), dans build/examples-hip, puis
# lance les tests qui ne demandent pas de GPU.
#
# Utilisation : scripts/build-hip-amd.sh [hypre|mfem|mfem-mgis|examples|all]

set -euo pipefail

P=$(cd "$(dirname "$0")/.." && pwd)
ROCM_VENV=${ROCM_VENV:-$HOME/Bureau/rocm-sdk-therock/venv}
ROCM=$("$ROCM_VENV/bin/rocm-sdk" path --root)
HIP_ARCH=${HIP_ARCH:-gfx90a}
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
JOBS=${JOBS:-4}
# cœurs utilisés par la compilation
CPUS=${CPUS:-0,1,8,9}

location() {
  SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
    location -i "$1"
}
MPI=$(location mpi)
METIS=$(location metis)
SUITESPARSE=$(location suite-sparse)
OPENBLAS=$(location openblas)

D=$P/deps/hip-amd
mkdir -p "$D"

build_hypre() {
  local version=3.2.0
  # archive mise en cache par spack (somme de contrôle de la recette)
  local sha=5273205a310fb6aa3ae506ce216760fb67b30e02024874f3cdb8b811e4801de7
  local archive=$HOME/Bureau/spack-mmm/var/spack/cache/_source-cache/archive/${sha:0:2}/$sha.tar.gz
  echo "$sha  $archive" | sha256sum -c -
  rm -rf "$D/hypre-src" "$D/build-hypre"
  mkdir -p "$D/hypre-src"
  tar -xf "$archive" -C "$D/hypre-src" --strip-components=1
  cmake -S "$D/hypre-src/src" -B "$D/build-hypre" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$MPI/bin/mpicc" \
    -DCMAKE_CXX_COMPILER="$MPI/bin/mpic++" \
    -DHYPRE_ENABLE_MPI=ON -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX="$D/install/hypre"
  taskset -c "$CPUS" cmake --build "$D/build-hypre" -j "$JOBS"
  cmake --install "$D/build-hypre"
  echo "hypre $version installé dans $D/install/hypre"
}

build_mfem() {
  # MFEM 4.10, comme l'environnement mmm-gpu ; les sources déjà extraites
  # sont gardées (la compilation reprend où elle s'est arrêtée)
  if [ ! -f "$D/mfem-src/makefile" ]; then
    mkdir -p "$D/mfem-src"
    git -C "$HOME/Documents/GitHub/mfem-amgx-pr" archive v4.10 |
      tar -x -C "$D/mfem-src"
  fi
  local blas="-Wl,-rpath,$OPENBLAS/lib -L$OPENBLAS/lib -lopenblas"
  export HIPCC_COMPILE_FLAGS_APPEND="-isystem $ROCM/include"
  # plateforme AMD, même si scripts/env.sh choisit NVIDIA pour TFEL
  export HIP_PATH=$ROCM ROCM_PATH=$ROCM HIP_PLATFORM=amd
  cd "$D/mfem-src"
  # `make config` récrit config/config.mk, dont dépendent tous les objets :
  # il n'est relancé que pour une nouvelle extraction des sources
  [ -f config/config.mk ] || make config \
    PREFIX="$D/install/mfem" \
    MFEM_USE_MPI=YES MFEM_USE_HIP=YES MFEM_USE_METIS_5=YES \
    MFEM_USE_SUITESPARSE=YES MFEM_USE_ZLIB=YES \
    MFEM_SHARED=YES MFEM_STATIC=NO \
    HIP_CXX="$ROCM/bin/hipcc" HIP_DIR="$ROCM" HIP_ARCH="$HIP_ARCH" \
    MPICXX="$MPI/bin/mpic++" \
    HYPRE_DIR="$D/install/hypre" \
    HYPRE_LIB="-Wl,-rpath,$D/install/hypre/lib -L$D/install/hypre/lib -lHYPRE" \
    METIS_DIR="$METIS" \
    METIS_LIB="-Wl,-rpath,$METIS/lib -L$METIS/lib -lmetis" \
    SUITESPARSE_DIR="$SUITESPARSE" \
    SUITESPARSE_OPT="-I$SUITESPARSE/include -I$SUITESPARSE/include/suitesparse" \
    LAPACK_LIB="$blas"
  taskset -c "$CPUS" make -j "$JOBS" lib
  make install
  echo "MFEM installé dans $D/install/mfem"
}

# environnement propre, sans scripts/env.sh : la vue de l'environnement
# spack mmm-gpu contient le MFEM compilé avec CUDA, dont les en-têtes
# passeraient avant ceux du MFEM compilé avec HIP (le dossier include de la
# vue est celui de MPI). Seuls le TFEL du prototype et MPI sont dans les
# chemins, aucune autre installation de TFEL ou de MFEM.
build_environment() {
  tfel=$P/install/tfel-mmm-gpu
  unset CMAKE_PREFIX_PATH PKG_CONFIG_PATH HIP_PLATFORM HIP_COMPILER
  export PATH=$tfel/bin:$MPI/bin:/usr/local/bin:/usr/bin:/bin
  export LD_LIBRARY_PATH=$tfel/lib
  export HIPCC_COMPILE_FLAGS_APPEND="-isystem $ROCM/include"
}

build_mfem_mgis() {
  build_environment
  local b=$P/build/mfem-mgis-hip
  cmake -S "$P/mfem-mgis" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -DCMAKE_HIP_FLAGS="-isystem $ROCM/include" \
    -DMFEM_DIR="$D/install/mfem" \
    -DMFrontGenericInterface_DIR="$P/install/mgis/share/mgis/cmake" \
    -DTFELHOME="$tfel" -DTFEL_DIR="$tfel/share/tfel/cmake" \
    -DCMAKE_INSTALL_PREFIX="$P/install/mfem-mgis-hip"
  taskset -c "$CPUS" cmake --build "$b" -j "$JOBS"
  taskset -c "$CPUS" cmake --build "$b" -j "$JOBS" --target build-tests
  run_tests "$b"
}

# lance les tests de la construction donnée qui ne demandent pas de GPU.
# Un MFEM construit avant l'ajout des chemins de hypre et de METIS à
# HYPRE_LIB et METIS_LIB ne les trouve pas seul à l'exécution.
run_tests() {
  (cd "$1" &&
    LD_LIBRARY_PATH="$D/install/hypre/lib:$METIS/lib:$LD_LIBRARY_PATH" \
    taskset -c "$CPUS" ctest -E gpu --output-on-failure)
}

build_examples() {
  build_environment
  cmake --install "$P/build/mfem-mgis-hip"
  local b=$P/build/examples-hip
  cmake -S "$P/mfem-mgis-examples" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -DCMAKE_HIP_FLAGS="-isystem $ROCM/include" \
    -DMFEMMGIS_DIR="$P/install/mfem-mgis-hip/share/mfem-mgis/cmake" \
    -DMFEM_DIR="$D/install/mfem"
  taskset -c "$CPUS" cmake --build "$b" -j "$JOBS"
  run_tests "$b"
}

case "${1:-all}" in
  hypre) build_hypre ;;
  mfem) build_mfem ;;
  mfem-mgis) build_mfem_mgis ;;
  examples) build_examples ;;
  all) build_hypre; build_mfem; build_mfem_mgis; build_examples ;;
  *) echo "usage: $0 [hypre|mfem|mfem-mgis|examples|all]" >&2; exit 1 ;;
esac
