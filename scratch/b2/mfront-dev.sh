#!/bin/bash
# lance le mfront de l'arbre de build (build/tfel) avec les en-têtes des sources
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
source "$P/scripts/env.sh"
LD_LIBRARY_PATH=$(find "$P/build/tfel" -name "*.so" -printf "%h\n" | sort -u | paste -sd:):${LD_LIBRARY_PATH:-} \
  exec "$P/build/tfel/mfront/src/mfront" "$@"
