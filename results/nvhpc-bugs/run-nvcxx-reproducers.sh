#!/bin/bash
# Compile and run the reproducers of the defects of nvc++ with a given
# compiler, at -O0 to -O3, and print one line per case (PASS, FAIL or
# COMPILE-ERROR).
# usage: run-nvcxx-reproducers.sh path/to/nvc++ [work directory]
set -u
NVCXX=$1
W=${2:-$(mktemp -d)}
R=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$W" && cd "$W" || exit 1
echo "# $("$NVCXX" --version | grep -m1 -o 'nvc++ [0-9.-]*')"
# check: name, source, flags, expected output (grep -E pattern)
check() {
  local name=$1 src=$2 flags=$3 expected=$4
  if ! "$NVCXX" $flags -std=c++20 "$R/$src" -o t > log 2>&1; then
    printf "%-28s %-40s COMPILE-ERROR (%s)\n" "$name" "$flags" \
      "$(grep -m1 -i -E 'error|fatal|-[FS]-[0-9]' log | cut -c1-100)"
    return
  fi
  local out
  out=$(timeout 60 ./t 2>&1)
  if echo "$out" | grep -qE "$expected"; then
    printf "%-28s %-40s PASS\n" "$name" "$flags"
  else
    printf "%-28s %-40s FAIL (%s)\n" "$name" "$flags" \
      "$(echo "$out" | grep -v '^$' | head -1 | cut -c1-80)"
  fi
}
for o in -O0 -O1 -O2 -O3; do
  # signed loop counters: all the lines must give 4, 3, 4 and 16
  if "$NVCXX" -cuda $o -std=c++20 "$R/signed-loop-trip-count.cu" -o t > log 2>&1; then
    bad=$(./t | grep -v -E 'sum 4: 4 \| arr 3: 3 \| post 4: 4 \| mul 16: 16$' | awk '{print $1}' | paste -sd, -)
    if [ -n "$bad" ]; then bad="FAIL ($bad)"; else bad=PASS; fi
    printf "%-28s %-40s %s\n" signed-loop-trip-count "-cuda $o" "$bad"
  else
    printf "%-28s %-40s COMPILE-ERROR\n" signed-loop-trip-count "-cuda $o"
  fi
  check empty-kernel-unified empty-kernel-stdpar.cxx "$o -stdpar=gpu -gpu=mem:unified" '^23\.6 '
  check empty-kernel-default empty-kernel-stdpar.cxx "$o -stdpar=gpu" '^23\.6 '
  check aggregate-zero-init aggregate-zero-init.cxx "$o -stdpar=gpu -gpu=mem:managed" '^2 '
  check member-function-pointer member-function-pointer.cxx "$o -stdpar=gpu -gpu=mem:managed" '^2 '
  check std-variant-visit std-variant-visit.cxx "$o -stdpar=gpu -gpu=mem:managed" '^5 '
done
