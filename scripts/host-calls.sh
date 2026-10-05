#!/bin/bash
# Generate the given test behaviours with the generic-parallel interface
# (CUDA backend) and list the host functions called from device code.
# usage: scripts/host-calls.sh law1 law2 ...
P=$(cd "$(dirname "$0")/.." && pwd)
V=/home/t/Bureau/spack/var/spack/environments/tfelgpu_local/.spack-env/view
export PATH=$P/install/tfel/bin:$V/bin:$PATH
export LD_LIBRARY_PATH=$P/install/tfel/lib:${LD_LIBRARY_PATH:-}
SRC=$P/tfel-gpu-proto/mfront/tests/behaviours
OUT=$P/scratch/host-calls
mkdir -p "$OUT"
run() {
  n=$1
  d=$OUT/$n
  rm -rf "$d" && mkdir -p "$d" && cd "$d" || exit 1
  mfront --obuild --interface=generic-parallel --configuration-file=config-cuda-nvcc.json \
    --search-path="$SRC" --search-path="$SRC/../properties" --search-path="$SRC/../models" \
    "$SRC/$n.mfront" > build.log 2>&1
  r=$(cuobjdump --dump-resource-usage src/*-cuda.so 2>/dev/null | sed -n 's/.*REG:\([0-9]*\).*/\1/p' | sort -n | head -1)
  calls=$(grep -o 'calling a __host__ function("[^"]*")' build.log | sed -E 's/calling a __host__ function\("//; s/"\)$//' | sed -E 's/\(.*//' | sort -u | tr '\n' ';' | cut -c1-400)
  echo -e "$n\tminREG=${r:-none}\t$calls"
}
export -f run
export OUT SRC
printf "%s\n" "$@" | xargs -P 12 -I{} bash -c 'run {}'
