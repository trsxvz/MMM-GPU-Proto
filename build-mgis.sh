#!/bin/bash
# Compile, teste et installe le fork de MGIS du prototype (dépôt `mgis/`,
# branche `gpu-proto`), contre le TFEL installé par build-tfel.sh, dans
# l'environnement de scripts/env.sh (mmm-gpu par défaut).
#
# usage : ./build-mgis.sh [full]
#
# Sans argument : construction dans build/mgis (Release, tests CUDA des
# bibliothèques de lots), cible `check` (tests compris), puis installation
# dans install/mgis, où mfem-mgis trouve MGIS.
#
# Avec `full` : construction dans build/mgis-full, avec HDF5 et les liaisons
# C, Fortran et Python (HDF5, pybind11, Python et numpy de l'installation
# spack de la machine), cible `check`, sans installation.
#
# Un dossier de construction existant est réutilisé, et reconfiguré avec
# les options de ce script à chaque appel. Pour repartir de zéro, le
# supprimer avant.
set -e
P=$(cd "$(dirname "$0")" && pwd)
source "$P/scripts/env.sh"
case "${1:-default}" in
  default)
    BUILD=$P/build/mgis
    OPTIONS=(-Denable-generic-parallel-cuda-tests=ON
      -Denable-hdf5-automatic-support=OFF
      -DCMAKE_INSTALL_PREFIX="$P/install/mgis")
    ;;
  full)
    SPACK_OPT=/home/t/Bureau/spack/opt/spack/linux-x86_64_v4
    HDF5=$SPACK_OPT/hdf5-1.14.5-fipxiy3en2unepc3jrn44pgxx2l7uiw2
    PYBIND11=$SPACK_OPT/py-pybind11-3.0.1-aktpzf47fftmlfpa7p5d4bwffzzvwjf5
    PYTHON=$SPACK_OPT/python-3.14.0-4hccxfwtgv4h6shmhneyyivlki4yk3mh
    NUMPY=$SPACK_OPT/py-numpy-2.3.4-5arojcw3zrhljmwvetbaly2buuxcuvlh
    # numpy, HDF5 and Python found by CMake and by the tests
    export PYTHONPATH=$NUMPY/lib/python3.14/site-packages
    export LD_LIBRARY_PATH=$HDF5/lib:$PYTHON/lib:$LD_LIBRARY_PATH
    BUILD=$P/build/mgis-full
    OPTIONS=(-DCMAKE_Fortran_COMPILER=gfortran-14
      -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL;$HDF5;$PYBIND11"
      -DPython_EXECUTABLE="$PYTHON/bin/python3"
      -Denable-hdf5-support=ON -Denable-c-bindings=ON
      -Denable-fortran-bindings=ON -Denable-python-bindings=ON)
    ;;
  *)
    echo "usage: $0 [full]" >&2
    exit 1
    ;;
esac
cmake -S "$P/mgis" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
  -DCMAKE_CUDA_COMPILER="$PROTO_CUDA/bin/nvcc" \
  -DCMAKE_CUDA_HOST_COMPILER=g++-14 \
  -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL" \
  -Denable-testing=ON -Denable-mgis-function=ON "${OPTIONS[@]}"
cmake --build "$BUILD" -j8
# the tests of the C and Fortran bindings are built by `check`, not by
# `build-tests`
cmake --build "$BUILD" -j8 --target check
if [ "${1:-default}" = default ]; then
  cmake --install "$BUILD"
fi
