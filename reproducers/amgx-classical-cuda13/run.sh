#!/bin/bash
# Runs the amgx_capi example of AmgX (conjugate gradient preconditioned by
# its algebraic multigrid) on a Poisson matrix, with the classical and the
# aggregation multigrids, for each installation of AmgX given (by default,
# the one of the spack environment). Each installation is run with its own
# example (the examples have a RPATH to their library).
# usage: reproducers/amgx-classical-cuda13/run.sh [prefix of AmgX...]
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
SPACK=${SPACK:-$HOME/Bureau/spack-mmm/bin/spack}
W=$P/scratch/reproducers/amgx-classical-cuda13
mkdir -p "$W" && cd "$W" || exit 1
python3 "$R/poisson.py" 20 poisson.mtx
prefixes=("$@")
if [ ${#prefixes[@]} -eq 0 ]; then
  prefixes=("$(SPACK_DISABLE_LOCAL_CONFIG=1 "$SPACK" -e "$P/deps/spack/mmm-gpu" \
    location -i amgx)")
fi
for A in "${prefixes[@]}"; do
  echo "== $A"
  for c in PCG_CLASSICAL_V_JACOBI PCG_AGGREGATION_JACOBI; do
    r=$(timeout 120 "$A/lib/examples/amgx_capi" -m poisson.mtx \
      -c "$A/lib/configs/$c.json" 2>&1 |
      grep -E "Thrust failure|Solver Failed|Total Iterations|Total Time" |
      tr -s ' \t' ' ' | tr '\n' ' ')
    echo "$c: $r"
  done
done
