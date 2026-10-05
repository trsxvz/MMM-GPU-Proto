#!/bin/bash
# Lance les tests installés de l'interface generic-parallel (lois et
# propriétés matériau) pour toutes les configurations, dans une copie
# jetable (scratch/tc par défaut), avec l'environnement de scripts/env.sh.
# Usage : scripts/run-generic-parallel-tests.sh [motif de configuration]
# Variables :
# - GPB_OPTIMIZATION_LEVEL (par exemple -O1 ou -O3) : remplace -O2 dans les
#   fichiers de configuration (copiés dans chaque dossier de test, où MFront
#   les cherche en premier) et dans les commandes des tests ;
# - GPB_TESTS_DIRECTORY : dossier de la copie des tests.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
# MFront appelle make, qui lit MAKEFLAGS : compilation des lois en parallèle
# (6 compilations GPU simultanées tiennent dans 30 Go de mémoire)
export MAKEFLAGS=${MAKEFLAGS:--j6}
T=${GPB_TESTS_DIRECTORY:-$PROTO/scratch/tc}
level=${GPB_OPTIMIZATION_LEVEL:-}
configurations=$PROTO/install/tfel/share/tfel/mfront/configuration-files
rm -rf "$T" && mkdir -p "$T"
cp -r "$PROTO/install/tfel/share/tfel/test-suite/behaviours/generic-parallel" "$T/behaviours"
cp -r "$PROTO/install/tfel/share/tfel/test-suite/material-properties/generic-parallel" "$T/material-properties"
pattern=${1:-.}
status=0
for kind in behaviours material-properties; do
  for d in "$T/$kind"/*/; do
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
