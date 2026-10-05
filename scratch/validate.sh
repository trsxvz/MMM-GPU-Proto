#!/bin/bash
# installation, régression TFEL (environnement propre) et matrice complète
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
cd "$P" || exit 1
./build-tfel.sh > scratch/validate-install.log 2>&1; echo "install=$?"
cd "$P/build/tfel" || exit 1
CLEANPATH=$(echo "$PATH" | tr ':' '\n' | grep -v spack | paste -sd:)
touch mfront/src/mfront
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" make -k -j12 build-tests > ../../scratch/validate-build-tests.log 2>&1; echo "build-tests=$?"
env -u LD_LIBRARY_PATH PATH="$CLEANPATH" ctest -j12 > ../../scratch/validate-ctest.log 2>&1; echo "ctest=$?"
grep -E "tests passed|tests failed" ../../scratch/validate-ctest.log
cd "$P" || exit 1
scripts/run-generic-parallel-tests.sh > results/generic-parallel-tests-2026-10-05-c.txt 2>&1; echo "matrix=$?"
cat results/generic-parallel-tests-2026-10-05-c.txt
