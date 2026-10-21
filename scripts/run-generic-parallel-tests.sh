#!/bin/bash
# Runs the installed tests of the generic-parallel interface (behaviours
# and material properties) for all the configurations, in a disposable
# copy (scratch/tc by default), with the environment of scripts/env.sh. The
# copy has the layout of the installed test-suite
# (<kind>/generic-parallel/<configuration>): the comparison of the material
# properties includes a header of the tests of the behaviours.
# Usage: scripts/run-generic-parallel-tests.sh [configuration pattern]
# Variables:
# - GPB_OPTIMIZATION_LEVEL (for example -O1 or -O3): replaces -O2 in the
#   configuration files (copied in each test directory, where MFront looks
#   for them first) and in the commands of the tests;
# - GPB_TESTS_DIRECTORY: directory of the copy of the tests.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
# MFront calls make, which reads MAKEFLAGS: the behaviours are compiled in
# parallel (6 simultaneous GPU compilations fit in 30 GB of memory)
export MAKEFLAGS=${MAKEFLAGS:--j6}
T=${GPB_TESTS_DIRECTORY:-$PROTO/scratch/tc}
level=${GPB_OPTIMIZATION_LEVEL:-}
configurations=$PROTO_TFEL_INSTALL/share/tfel/mfront/configuration-files
# an existing directory is only deleted if it is in scratch/ or if this
# script created it
case "$T" in
  "$PROTO"/scratch/?*) ;;
  *)
    if [ -e "$T" ] && [ ! -f "$T/.generic-parallel-tests-copy" ]; then
      echo "$T exists and is not a copy made by this script" >&2
      exit 1
    fi
    ;;
esac
rm -rf -- "$T"
mkdir -p "$T/behaviours" "$T/material-properties"
touch "$T/.generic-parallel-tests-copy"
cp -r "$PROTO_TFEL_INSTALL/share/tfel/test-suite/behaviours/generic-parallel" "$T/behaviours/"
cp -r "$PROTO_TFEL_INSTALL/share/tfel/test-suite/material-properties/generic-parallel" "$T/material-properties/"
pattern=${1:-.}
status=0
for kind in behaviours material-properties; do
  for d in "$T/$kind/generic-parallel"/*/; do
    [ -f "$d/test.check" ] || continue
    name=$(basename "$d")
    echo "$name" | grep -qE "$pattern" || continue
    if [ -n "$level" ]; then
      # the test of the configuration and its sub-tests (TDLS)
      for t in $(find "$d" -name test.check); do
        td=$(dirname "$t")
        for c in $(grep -o -- '--configuration-file=[^ "]*' "$t" | sed 's/--configuration-file=//' | sort -u); do
          sed "s/-O2/$level/g" "$configurations/$c" > "$td/$c"
        done
        sed -i "s/-O2/$level/g" "$t"
      done
    fi
    start=$(date +%s)
    out=$(cd "$d" && tfel-check 2>&1)
    end=$(date +%s)
    # a configuration may have sub-tests (TDLS): one failure is a failure.
    # tfel-check reports a test whose requirements are not met as a success
    # whose commands are all skipped: such a test is counted as skipped.
    counts=$(echo "$out" | awk '
      /beginning of test/ { ran = 0 }
      /\*\* Exec/ && /SUCCESS/ { ran = 1 }
      /end of test/ { if ($0 ~ /FAILED/) f++; else if (ran) p++; else k++ }
      END { print p + 0, f + 0, k + 0 }')
    read -r passed failed skipped <<< "$counts"
    if [ "$failed" -ne 0 ]; then
      r=FAILED; status=1
    elif [ "$passed" -ne 0 ]; then
      r=SUCCESS
    else
      r=SKIPPED
    fi
    printf "%-20s %-32s %-8s %4d s (%d passed, %d skipped)\n" "$kind" "$name" "$r" $((end - start)) "$passed" "$skipped"
  done
done
exit $status
