#!/bin/bash
# Build ilu.cxx against the MFEM and hypre of the prototype and run it on
# the CPU and on the GPU (two processes), with the preconditioner of MFEM
# and with a preconditioner kept from one operator to the other, for the
# levels of fill 0 and 1.
# usage: reproducers/hypre-ilu-gpu/run.sh
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
M=$(grep "^MFEM_DIR" "$P/build/mfem-mgis/CMakeCache.txt" | cut -d= -f2)
V=$PROTO_SPACK_VIEW
W=$P/scratch/reproducers/hypre-ilu-gpu
mkdir -p "$W" && cd "$W" || exit 1
# only the headers of MFEM needed (mfem.hpp includes AmgX)
mpicxx -O2 -std=c++20 -I"$M/include" -I"$V/include" -I"$PROTO_CUDA/include" \
  "$R/ilu.cxx" -o ilu -L"$M/lib" -lmfem -L"$V/lib" -lHYPRE \
  -Wl,-rpath,"$M/lib" -Wl,-rpath,"$V/lib" || exit 1
for args in "cpu 0" "cpu 1 1" "cuda 1 0" "cuda 1 1" "cuda 0"; do
  echo "== $args"
  timeout 120 mpirun -n 2 ./ilu $args 2>&1 |
    grep -E "resolution|Segmentation fault" | sort -u | head -4
done
