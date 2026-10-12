#!/bin/bash
# Compile static.cu with nvcc and the headers of TFEL installed by the
# prototype, for each mode (compilation only, no GPU needed).
# usage: reproducers/mfront-static-quantity-cuda/run.sh
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
P=$(cd "$R/../.." && pwd)
source "$P/scripts/env.sh" > /dev/null 2>&1
W=$P/scratch/reproducers/mfront-static-quantity-cuda
mkdir -p "$W" && cd "$W" || exit 1
for mode in 1 2 3 4; do
  if nvcc -std=c++20 -arch=sm_120 --expt-relaxed-constexpr -DMODE=$mode \
    -I"$PROTO_TFEL_INSTALL/include" "$R/static.cu" -o static-$mode \
    > log-$mode 2>&1; then
    echo "mode $mode: OK"
  else
    echo "mode $mode: $(grep -m1 'error' log-$mode)"
  fi
done
