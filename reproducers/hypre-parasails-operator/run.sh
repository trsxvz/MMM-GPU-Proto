#!/bin/bash
# Build parasails.cxx and newton.cxx against the MFEM and hypre of the
# prototype and run them on two processes.
# - parasails: preconditioner of MFEM, first operator kept alive, fixed
#   preconditioner. The allocator overwrites the start of the memory freed,
#   where the communicator of the destroyed operator was stored. On a GPU
#   (cuda), all the modes fail at the first setup: ParaSails only computes
#   on the host.
# - newton (CPU only): NewtonSolver and ParNonlinearForm of MFEM, with the
#   preconditioner of MFEM and with the fixed one.
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
for p in parasails newton; do
  mpicxx -O2 -std=c++20 -w -I"$M/include" -I"$V/include" -I"$PROTO_CUDA/include" \
    "$R/$p.cxx" -o $p -L"$M/lib" -lmfem -L"$V/lib" -lHYPRE \
    -Wl,-rpath,"$M/lib" -Wl,-rpath,"$V/lib" || exit 1
done
for mode in mfem kept fixed; do
  echo "== parasails $D $mode"
  timeout 120 mpirun -n 2 ./parasails "$D" "$mode" 2>&1 |
    grep -aE "resolution|Segmentation fault|MPI_ERR" | sort -u | head -4
done
for mode in mfem fixed; do
  echo "== newton cpu $mode"
  timeout 120 mpirun -n 2 ./newton "$mode" 2>&1 |
    grep -aE "newton:|Segmentation fault|MPI_ERR" | sort -u | head -4
done
