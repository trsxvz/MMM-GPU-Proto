#!/usr/bin/env bash
# Construction d'un MFEM compilé avec HIP pour un GPU AMD, puis de mfem-mgis
# et des exemples avec lui. L'architecture est donnée par HIP_ARCH :
# - gfx90a (par défaut, MI210 et MI250) : compilation et édition de liens
#   seulement, la machine de développement n'ayant pas ce GPU ;
# - gfx1036 : le GPU intégré de la machine de développement (Ryzen), sur
#   lequel les tests s'exécutent vraiment avec GPU_TESTS=ON. /dev/kfd est
#   réservé au groupe render : les tests passent par `sg render` si la
#   session a été ouverte avant l'ajout à ce groupe.
# Les dossiers d'une architecture autre que gfx90a portent son nom en
# suffixe (deps/hip-amd/mfem-src-gfx1036, build/mfem-mgis-hip-gfx1036...) ;
# hypre, compilé pour le CPU, est commun.
#
# ROCm vient des paquets Python de ROCm (TheRock), installés sans droits
# d'administrateur dans un environnement virtuel :
#   python3 -m venv $ROCM_VENV
#   $ROCM_VENV/bin/pip install --index-url \
#     https://stable.repo.amd.com/rocm/whl-next/ \
#     "rocm[libraries,devel,device-gfx90a]==10.1.0"
#   $ROCM_VENV/bin/rocm-sdk init
# et, pour une autre architecture, son paquet de noyaux des bibliothèques :
#   $ROCM_VENV/bin/pip install --index-url \
#     https://stable.repo.amd.com/rocm/whl-next/ \
#     "rocm-sdk-device-gfx1036==10.1.0"
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
# hypre est compilé pour le CPU (MFEM avec HIP accepte un hypre sur l'hôte),
# ou avec HIP pour le GPU si HYPRE_HIP=ON (nécessaire aux tests sur le GPU :
# les réglages d'exécution refusent un solveur linéaire de hypre qui calcule
# sur l'hôte quand il est placé sur le GPU), depuis l'archive mise en cache
# par spack. Avec HYPRE_HIP=ON, MFEM et la suite portent le suffixe
# supplémentaire -hypre-hip. MFEM 4.10 vient du dépôt officiel
# (étiquette v4.10). MPI, METIS, SuiteSparse et OpenBLAS sont ceux de
# l'environnement spack mmm-gpu. Résultat : deps/hip-amd/install/{hypre,mfem}.
#
# L'étape mfem-mgis construit mfem-mgis (dossier mfem-mgis du prototype)
# avec ce MFEM, MGIS et TFEL du prototype, dans build/mfem-mgis-hip, avec
# ses tests (cibles build-tests et BehaviourTest), ceux qui demandent un GPU
# seulement avec GPU_TESTS=ON, puis les lance.
#
# L'étape examples installe ce mfem-mgis dans install/mfem-mgis-hip et
# construit les exemples de mfem-mgis-examples avec lui (configuration
# installée, lois compilées pour le GPU AMD), dans build/examples-hip, puis
# lance ses tests (ceux sur GPU avec GPU_TESTS=ON, option
# MFEM_MGIS_EXAMPLES_GPU_TESTS).
#
# Utilisation : [HIP_ARCH=gfx1036] [GPU_TESTS=ON] [HYPRE_HIP=ON] \
#   scripts/build-hip-amd.sh [hypre|mfem|mfem-mgis|examples|all]

set -euo pipefail

P=$(cd "$(dirname "$0")/.." && pwd)
ROCM_VENV=${ROCM_VENV:-$HOME/Bureau/rocm-sdk-therock/venv}
ROCM=$("$ROCM_VENV/bin/rocm-sdk" path --root)
HIP_ARCH=${HIP_ARCH:-gfx90a}
GPU_TESTS=${GPU_TESTS:-OFF}
HYPRE_HIP=${HYPRE_HIP:-OFF}
# suffixe des dossiers propres à l'architecture
if [ "$HIP_ARCH" = gfx90a ]; then
  A=
else
  A=-$HIP_ARCH
fi
# suffixes de hypre et des dossiers qui en dépendent
if [ "$HYPRE_HIP" = ON ]; then
  H=-hip-$HIP_ARCH
  A=$A-hypre-hip
else
  H=
fi
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
JOBS=${JOBS:-4}
# cœurs utilisés par la compilation et par les tests
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
  local spack_root
  spack_root=$(dirname "$(dirname "$SPACK")")
  local archive=$spack_root/var/spack/cache/_source-cache/archive/${sha:0:2}/$sha.tar.gz
  echo "$sha  $archive" | sha256sum -c -
  local src=$D/hypre-src$H b=$D/build-hypre$H
  rm -rf "$src" "$b"
  mkdir -p "$src"
  tar -xf "$archive" -C "$src" --strip-components=1
  local gpu=()
  if [ "$HYPRE_HIP" = ON ]; then
    # hypre compilé par le compilateur HIP de ROCm (qui n'est pas l'enveloppe
    # de MPI : les en-têtes de MPI lui sont donnés), bibliothèques rocsparse,
    # rocrand et rocblas de ROCm ; sans Umpire, que hypre active d'office
    # pour un GPU (non installé ici : hypre alloue alors lui-même)
    export HIP_PATH=$ROCM ROCM_PATH=$ROCM HIP_PLATFORM=amd
    gpu=(-DHYPRE_ENABLE_HIP=ON -DCMAKE_HIP_ARCHITECTURES="$HIP_ARCH"
      -DHYPRE_ENABLE_UMPIRE=OFF
      -DCMAKE_HIP_COMPILER="$ROCM/lib/llvm/bin/clang++"
      -DCMAKE_HIP_FLAGS="-isystem $ROCM/include -isystem $MPI/include"
      -DCMAKE_PREFIX_PATH="$ROCM")
  fi
  cmake -S "$src/src" -B "$b" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$MPI/bin/mpicc" \
    -DCMAKE_CXX_COMPILER="$MPI/bin/mpic++" \
    -DHYPRE_ENABLE_MPI=ON -DBUILD_SHARED_LIBS=ON \
    "${gpu[@]}" \
    -DCMAKE_INSTALL_PREFIX="$D/install/hypre$H"
  taskset -c "$CPUS" cmake --build "$b" -j "$JOBS"
  cmake --install "$b"
  echo "hypre $version installé dans $D/install/hypre$H"
}

build_mfem() {
  # MFEM 4.10, comme l'environnement mmm-gpu ; les sources déjà extraites
  # sont gardées (la compilation reprend où elle s'est arrêtée)
  local src=$D/mfem-src$A
  if [ ! -f "$src/makefile" ]; then
    git clone --quiet --depth 1 --branch v4.10 \
      https://github.com/mfem/mfem.git "$src"
  fi
  local blas="-Wl,-rpath,$OPENBLAS/lib -L$OPENBLAS/lib -lopenblas"
  export HIPCC_COMPILE_FLAGS_APPEND="-isystem $ROCM/include"
  # plateforme AMD, même si scripts/env.sh choisit NVIDIA pour TFEL
  export HIP_PATH=$ROCM ROCM_PATH=$ROCM HIP_PLATFORM=amd
  cd "$src"
  # `make config` récrit config/config.mk, dont dépendent tous les objets :
  # il n'est relancé que pour une nouvelle extraction des sources
  [ -f config/config.mk ] || make config \
    PREFIX="$D/install/mfem$A" \
    MFEM_USE_MPI=YES MFEM_USE_HIP=YES MFEM_USE_METIS_5=YES \
    MFEM_USE_SUITESPARSE=YES MFEM_USE_ZLIB=YES \
    MFEM_SHARED=YES MFEM_STATIC=NO \
    HIP_CXX="$ROCM/bin/hipcc" HIP_DIR="$ROCM" HIP_ARCH="$HIP_ARCH" \
    MPICXX="$MPI/bin/mpic++" \
    HYPRE_DIR="$D/install/hypre$H" \
    HYPRE_LIB="-Wl,-rpath,$D/install/hypre$H/lib -L$D/install/hypre$H/lib -lHYPRE" \
    METIS_DIR="$METIS" \
    METIS_LIB="-Wl,-rpath,$METIS/lib -L$METIS/lib -lmetis" \
    SUITESPARSE_DIR="$SUITESPARSE" \
    SUITESPARSE_OPT="-I$SUITESPARSE/include -I$SUITESPARSE/include/suitesparse" \
    LAPACK_LIB="$blas"
  taskset -c "$CPUS" make -j "$JOBS" lib
  make install
  echo "MFEM installé dans $D/install/mfem$A"
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
  local b=$P/build/mfem-mgis-hip$A
  cmake -S "$P/mfem-mgis" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -DCMAKE_HIP_FLAGS="-isystem $ROCM/include" \
    -DMFEM_DIR="$D/install/mfem$A" \
    -DMFrontGenericInterface_DIR="$P/install/mgis/share/mgis/cmake" \
    -DTFELHOME="$tfel" -DTFEL_DIR="$tfel/share/tfel/cmake" \
    -Denable-gpu-tests="$GPU_TESTS" \
    -DCMAKE_INSTALL_PREFIX="$P/install/mfem-mgis-hip$A"
  taskset -c "$CPUS" cmake --build "$b" -j "$JOBS"
  # la bibliothèque de lois des tests (BehaviourTest) dépend de la cible
  # check, pas de build-tests
  taskset -c "$CPUS" cmake --build "$b" -j "$JOBS" \
    --target build-tests BehaviourTest
  run_tests "$b"
}

# lance les tests de la construction donnée ; sur le GPU, par `sg render`
# si la session n'a pas encore le groupe render (sg, programme setuid,
# retire LD_LIBRARY_PATH de l'environnement : il est redonné)
#
# Sur le GPU, les noyaux et les copies sont sérialisés : avec ROCm 10.1
# (HIP 7.16, paquets de TheRock) sur le GPU intégré gfx1036, un noyau encore
# en attente est perdu quand un noyau d'une bibliothèque partagée dont le
# code n'est pas encore chargé est lancé pour la première fois (reproducteur
# reproducers/hip-lost-kernel) : résidus faux (composantes imposées non
# annulées) et fautes mémoire de hypre. Le chargement immédiat du code
# (HIP_ENABLE_DEFERRED_LOADING=0) échoue, libhipblaslt n'ayant pas de code
# pour gfx1036. La sérialisation ne suffit pas toujours : rve et
# ssna303_3d_hypre restent alors bloqués dans hipLaunchKernel, au premier
# lancement d'un noyau de la bibliothèque des lois (attente d'un signal qui
# ne vient jamais).
#
# L'OpenMPI de l'environnement mmm-gpu est construit avec CUDA : sans
# OMPI_MCA_accelerator=null, chaque processus ouvre un contexte sur le GPU
# NVIDIA de la machine, qui peut servir à d'autres travaux.
run_tests() {
  local serialize="export OMPI_MCA_accelerator=null; "
  if [ "$GPU_TESTS" = ON ]; then
    serialize="${serialize}export AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3; "
  fi
  # délai de 600 s par test : un test bloqué ne garde pas les cœurs 1500 s
  local cmd="${serialize}cd '$1' && taskset -c '$CPUS' ctest --timeout 600 --output-on-failure"
  if [ "$GPU_TESTS" = ON ] && ! id -nG | grep -qw render; then
    sg render -c "export LD_LIBRARY_PATH='$LD_LIBRARY_PATH'; $cmd"
  else
    bash -c "$cmd"
  fi
}

build_examples() {
  build_environment
  cmake --install "$P/build/mfem-mgis-hip$A"
  local b=$P/build/examples-hip$A
  cmake -S "$P/mfem-mgis-examples" -B "$b" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DMPI_C_COMPILER="$MPI/bin/mpicc" -DMPI_CXX_COMPILER="$MPI/bin/mpicxx" \
    -DCMAKE_HIP_FLAGS="-isystem $ROCM/include" \
    -DMFEMMGIS_DIR="$P/install/mfem-mgis-hip$A/share/mfem-mgis/cmake" \
    -DMFEM_DIR="$D/install/mfem$A" \
    -DMFEM_MGIS_EXAMPLES_GPU_TESTS="$GPU_TESTS"
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
