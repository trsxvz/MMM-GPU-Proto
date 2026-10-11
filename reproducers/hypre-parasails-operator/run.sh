#!/bin/bash
# Build parasails.cxx against the MFEM and hypre of the prototype and run
# it on two processes with the preconditioner of MFEM, with the first
# operator kept alive, and with the fixed preconditioner. The allocator
# overwrites the start of the memory freed, where the communicator of the
# destroyed operator was stored. On a GPU (cuda), all the modes fail at the
# first setup: ParaSails only computes on the host.
# usage: reproducers/hypre-parasails-operator/run.sh [cpu|cuda]
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
M=$(grep "^MFEM_DIR" "$P/build/mfem-mgis/CMakeCache.txt" | cut -d= -f2)
V=$PROTO_SPACK_VIEW
W=$P/scratch/reproducers/hypre-parasails-operator
D=${1:-cpu}
mkdir -p "$W" && cd "$W" || exit 1
# only the headers of MFEM needed (mfem.hpp includes AmgX)
mpicxx -O2 -std=c++20 -I"$M/include" -I"$V/include" -I"$PROTO_CUDA/include" \
  "$R/parasails.cxx" -o parasails -L"$M/lib" -lmfem -L"$V/lib" -lHYPRE \
  -Wl,-rpath,"$M/lib" -Wl,-rpath,"$V/lib" || exit 1
for mode in mfem kept fixed; do
  echo "== $D $mode"
  timeout 120 mpirun -n 2 ./parasails "$D" "$mode" 2>&1 |
    grep -aE "resolution|Segmentation fault|MPI_ERR" | sort -u | head -4
done
