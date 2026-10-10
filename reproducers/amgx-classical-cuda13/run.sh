#!/bin/bash
# Runs the amgx_capi example of AmgX (conjugate gradient preconditioned by
# its algebraic multigrid) on a Poisson matrix, with the classical and the
# aggregation multigrids, for the AmgX of the spack environment and, if
# built, for the AmgX of build/amgx-stub (same commit, compiled with
# -static-global-template-stub=false).
# usage: reproducers/amgx-classical-cuda13/run.sh
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
W=$P/scratch/reproducers/amgx-classical-cuda13
mkdir -p "$W" && cd "$W" || exit 1
python3 "$R/poisson.py" 20 poisson.mtx
# AmgX linked by MFEM
A=$(grep -o "/[^ ]*/amgx-[^/ ]*/include" \
  "$(grep '^MFEM_DIR' "$P/build/mfem-mgis/CMakeCache.txt" | cut -d= -f2)/share/mfem/config.mk" |
  head -1 | sed 's|/include$||')
# the examples have a RPATH to their library: each library is run with its
# own example
for d in "$A/lib|$A/lib/examples" "$P/build/amgx-stub/build|$P/build/amgx-stub/build/examples"; do
  IFS='|' read -r lib examples <<< "$d"
  if [ ! -x "$examples/amgx_capi" ]; then
    continue
  fi
  echo "== $lib"
  for c in PCG_CLASSICAL_V_JACOBI PCG_AGGREGATION_JACOBI; do
    r=$(timeout 120 "$examples/amgx_capi" -m poisson.mtx \
      -c "$A/lib/configs/$c.json" 2>&1 |
      grep -E "Thrust failure|Solver Failed|Total Iterations|Total Time" |
      tr -s ' \t' ' ' | tr '\n' ' ')
    echo "$c: $r"
  done
done
