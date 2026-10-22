#!/bin/sh
# Build the checks and run them for each order of declaration of the base
# class and each default visibility of the symbols. The other arguments are
# given to cmake (compiler, generator, flags).
# usage: reproducers/execution-policies/run.sh [build directory] [cmake options...]
set -e
R=$(dirname "$0")
W=${1:-$R/../../scratch/reproducers/execution-policies}
[ $# -gt 0 ] && shift
for first in OFF ON; do
  for hidden in ON OFF; do
    B=$W/identifier-first-$first-hidden-$hidden
    cmake -S "$R" -B "$B" -DEXECUTION_POLICIES_IDENTIFIER_FIRST=$first \
      -DEXECUTION_POLICIES_HIDDEN_VISIBILITY=$hidden "$@"
    cmake --build "$B" --config Release
    ctest --test-dir "$B" -C Release -V
  done
done
