#!/bin/bash
# Validation of the TFEL branch with all the configurations of build-tfel.sh:
# reconfiguration and installation, regression tests of TFEL, matrices of
# the generic-parallel tests at -O2, -O1 and -O3, compilation sweep of the
# behaviours of MFront for CUDA.
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
   -DTFEL_CLANGXX_COMPILER=/usr/bin/clang++-20 \
   -DTFEL_ICPX_COMPILER=/opt/intel/oneapi/compiler/2026.1/bin/icpx \
   -U "tfel-substitutions*" . > "$P/scratch/validate3-cmake.log" 2>&1); echo "cmake=$? ($CONFS)"
cd "$P" || exit 1
./build-tfel.sh > scratch/validate3-install.log 2>&1; echo "install=$?"
cd "$P/build/tfel" || exit 1
CLEANPATH=$(echo "$PATH" | tr ':' '\n' | grep -v spack | paste -sd:)
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" make -k -j8 build-tests > "$P/scratch/validate3-build-tests.log" 2>&1; echo "build-tests=$?"
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" ctest -j8 > "$P/scratch/validate3-ctest.log" 2>&1; echo "ctest=$?"
grep -E "tests passed|tests failed" "$P/scratch/validate3-ctest.log"
cd "$P" || exit 1
for level in O2 O1 O3; do
  if [ $level = O2 ]; then
    scripts/run-generic-parallel-tests.sh > "$R/generic-parallel-tests-$tag-$level.txt" 2>&1
  else
    GPB_OPTIMIZATION_LEVEL=-$level GPB_TESTS_DIRECTORY=$P/scratch/tc-$level \
      scripts/run-generic-parallel-tests.sh > "$R/generic-parallel-tests-$tag-$level.txt" 2>&1
  fi
  echo "matrix-$level=$?"
done
scripts/compile-sweep.sh config-cuda-nvcc.json 6 > scratch/validate3-sweep.log 2>&1; echo "sweep=$?"
tail -12 scratch/validate3-sweep.log
