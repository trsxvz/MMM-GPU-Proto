#!/bin/bash
# Validation of the TFEL branch with all the configurations of build-tfel.sh:
# reconfiguration and installation, regression tests of TFEL, matrix of
# the generic-parallel tests at -O3 (VALIDATE_LEVELS, for instance
# "O3 O2 O1", selects the levels), compilation sweeps of the
# behaviours, models and material properties of MFront for CUDA with nvcc
# and clang.
# MAKEFLAGS is lowered to -j4 to leave memory to other builds.
# usage: scratch/validate3.sh tag   (results/generic-parallel-tests-<tag>-*.txt)
tag=${1:-$(date +%F)}
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
R=$P/results
export MAKEFLAGS=-j4
cd "$P/build/tfel" || exit 1
CONFS=$(grep '^CONFIGURATIONS=' "$P/build-tfel.sh" | sed 's/^CONFIGURATIONS="\(.*\)"$/\1/' | sed 's/\$CONFIGURATIONS;//' | paste -sd';')
(export PATH=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view/bin:$PATH HIP_PLATFORM=nvidia CUDA_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/12.9
 cmake -Dgeneric-parallel-configurations="$CONFS" \
   -DTFEL_CLANGXX_COMPILER=/home/t/Bureau/spack/opt/spack/linux-x86_64_v4/llvm-20.1.8-v3gg3rqravuex5zixx6xglsqudvhcnjw/bin/clang++ \
   -DTFEL_ICPX_COMPILER=/opt/intel/oneapi/compiler/2026.1/bin/icpx \
   -DTFEL_CLANG_CUDA_PATH=/home/t/Bureau/spack-mmm-externals/cuda-12.9.1 \
   -U "tfel-substitutions*" . > "$P/scratch/validate3-cmake.log" 2>&1); echo "cmake=$? ($CONFS)"
cd "$P" || exit 1
./build-tfel.sh > scratch/validate3-install.log 2>&1; echo "install=$?"
cd "$P/build/tfel" || exit 1
CLEANPATH=$(echo "$PATH" | tr ':' '\n' | grep -v spack | paste -sd:)
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" make -k -j8 build-tests > "$P/scratch/validate3-build-tests.log" 2>&1; echo "build-tests=$?"
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" ctest -j8 > "$P/scratch/validate3-ctest.log" 2>&1; echo "ctest=$?"
grep -E "tests passed|tests failed" "$P/scratch/validate3-ctest.log"
cd "$P" || exit 1
for level in ${VALIDATE_LEVELS:-O3}; do
  if [ $level = O2 ]; then
    scripts/run-generic-parallel-tests.sh > "$R/generic-parallel-tests-$tag-$level.txt" 2>&1
  else
    GPB_OPTIMIZATION_LEVEL=-$level GPB_TESTS_DIRECTORY=$P/scratch/tc-$level \
      scripts/run-generic-parallel-tests.sh > "$R/generic-parallel-tests-$tag-$level.txt" 2>&1
  fi
  echo "matrix-$level=$?"
done
for compiler in nvcc clang; do
  scripts/compile-sweep.sh config-cuda-$compiler.json 6 > scratch/validate3-sweep-$compiler.log 2>&1
  echo "sweep-$compiler=$?"
  tail -12 scratch/validate3-sweep-$compiler.log
done
