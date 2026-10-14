#!/bin/bash
# runs the mfront of the build tree (build/tfel) with the headers of the sources
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
source "$P/scripts/env.sh"
LD_LIBRARY_PATH=$(find "$P/build/tfel" -name "*.so" -printf "%h\n" | sort -u | paste -sd:):${LD_LIBRARY_PATH:-} \
  exec "$P/build/tfel/mfront/src/mfront" "$@"
