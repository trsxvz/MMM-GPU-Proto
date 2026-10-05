#!/bin/bash
# reconfiguration (new configurations), installation, TFEL regression, then
# backend matrices at -O2, -O1 and -O3 and compile sweep of all the tests
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
cd "$P/build/tfel" || exit 1
CONFS="cuda-nvcc;hip-hipcc;stlpar-parunseq-gcc;stlpar-parunseq-nvhpc-gpu;kokkos;sycl-default-acpp;openmp-host-gcc"
(export PATH=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view/bin:$PATH HIP_PLATFORM=nvidia CUDA_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/12.9; cmake -Dgeneric-parallel-configurations="$CONFS" -U "tfel-substitutions*" . > ../../scratch/validate2-cmake.log 2>&1); echo "cmake=$?"
cd "$P" || exit 1
./build-tfel.sh > scratch/validate2-install.log 2>&1; echo "install=$?"
cd "$P/build/tfel" || exit 1
CLEANPATH=$(echo "$PATH" | tr ':' '\n' | grep -v spack | paste -sd:)
touch mfront/src/mfront
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" make -k -j12 build-tests > ../../scratch/validate2-build-tests.log 2>&1; echo "build-tests=$?"
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" ctest -j12 > ../../scratch/validate2-ctest.log 2>&1; echo "ctest=$?"
grep -E "tests passed|tests failed" ../../scratch/validate2-ctest.log
cd "$P" || exit 1
scripts/run-generic-parallel-tests.sh > results/generic-parallel-tests-2026-10-05-d-O2.txt 2>&1; echo "matrix-O2=$?"
GPB_OPTIMIZATION_LEVEL=-O1 GPB_TESTS_DIRECTORY=$P/scratch/tc-O1 scripts/run-generic-parallel-tests.sh > results/generic-parallel-tests-2026-10-05-d-O1.txt 2>&1; echo "matrix-O1=$?"
GPB_OPTIMIZATION_LEVEL=-O3 GPB_TESTS_DIRECTORY=$P/scratch/tc-O3 scripts/run-generic-parallel-tests.sh > results/generic-parallel-tests-2026-10-05-d-O3.txt 2>&1; echo "matrix-O3=$?"
scripts/compile-sweep.sh config-cuda-nvcc.json 10 > scratch/sweep-full.log 2>&1; echo "sweep=$?"
tail -12 scratch/sweep-full.log
