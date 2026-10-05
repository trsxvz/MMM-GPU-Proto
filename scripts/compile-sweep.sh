#!/bin/bash
# Compile all the test behaviours and models of MFront (sub-directories
# included) with the generic-parallel interface and the given configuration
# file, and summarize the results.
# usage: scripts/compile-sweep.sh [configuration file] [number of jobs]
#        [pattern selecting the files]
# The MFRONT variable selects the mfront executable (for instance
# scratch/b2/mfront-dev.sh to use the build tree).
# The results are written in scratch/sweep-<configuration>/summary.tsv: one
# line per file (path relative to mfront/tests, status, first error).
# Statuses:
# - OK: compiled, no call to a host function from device code;
# - HOST_CALLS: calls to host functions from device code (nvcc removes them
#   with a warning, so the results would be wrong);
# - SUSPICIOUS: a kernel uses very few registers (code removed silently);
# - USER_CODE: the code written by the user can't be used on GPUs
#   (exceptions, standard streams, etc.);
# - COMPILATION: compilation error in the generated code or in TFEL;
# - UNSUPPORTED: feature reported as not supported by MFront;
# - FRAGMENT: file meant to be imported by other files;
# - INVALID: the `generic` interface fails too (for instance the unit tests
#   of MFront which are expected to fail);
# - MFRONT: other errors reported by MFront.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
P=$PROTO
CONF=${1:-config-cuda-nvcc.json}
JOBS=${2:-12}
PATTERN=${3:-.}
TESTS=$P/tfel/mfront/tests
MFRONT=${MFRONT:-mfront}
OUT=$P/scratch/sweep-$(basename "${CONF%.json}")
rm -rf "$OUT" && mkdir -p "$OUT"
compile() {
  f=$1
  r=${f#"$TESTS"/}
  n=${r%.mfront}
  d=$OUT/${n//\//_}
  mkdir -p "$d" && cd "$d" || exit 1
  search="--search-path=$(dirname "$f") --search-path=$TESTS/behaviours"
  search="$search --search-path=$TESTS/properties --search-path=$TESTS/models"
  if $MFRONT --obuild --interface=generic-parallel --configuration-file="$CONF" \
      $search "$f" > build.log 2>&1; then
    # calls to host functions from device code are removed by nvcc (with a
    # warning), so the results would be wrong. A kernel with very few
    # registers is also suspicious: some constructs (for instance building a
    # std::string_view from a C string) make nvcc remove the code silently.
    s=OK
    e=""
    if grep -q "calling a __host__ function" build.log; then
      s=HOST_CALLS
      e=$(grep -o 'calling a __host__ function("[^"]*")' build.log | sort -u | head -3 | tr '\n' ' ' | cut -c1-300)
    elif [ -n "$(ls src/*-cuda.so 2>/dev/null)" ]; then
      g=$(cuobjdump --dump-resource-usage src/*-cuda.so 2>/dev/null |
            sed -n 's/.*REG:\([0-9]*\).*/\1/p' | sort -n | head -1)
      if [ -n "$g" ] && [ "$g" -lt 16 ]; then
        s=SUSPICIOUS
        e="a kernel uses only $g registers"
      fi
    fi
  else
    # calls to host functions may be reported as errors (`-diag-error`)
    if grep -q "calling a __host__ function" build.log; then
      s=HOST_CALLS
    elif grep -q "not supported yet" build.log; then
      s=UNSUPPORTED
    elif grep -q "no code block integrator defined" build.log; then
      s=FRAGMENT
    elif grep -q "make -C src" build.log; then
      if grep -qE "\.mfront\([0-9]+\): error" build.log; then
        s=USER_CODE
      else
        s=COMPILATION
      fi
    elif ! (mkdir -p generic && cd generic &&
              $MFRONT --interface=generic $search "$f" > build.log 2>&1); then
      s=INVALID
    else
      s=MFRONT
    fi
    if [ "$s" = HOST_CALLS ]; then
      e=$(grep -o 'calling a __host__ function("[^"]*")' build.log | sort -u | head -3 | tr '\n' ' ' | cut -c1-300)
    else
      e=$(grep -m1 -E "error|Error|ERROR|raise|unsupported" build.log | cut -c1-300)
    fi
  fi
  echo -e "$r\t$s\t$e" > status
}
export -f compile
export OUT CONF TESTS MFRONT
find "$TESTS/behaviours" "$TESTS/models" -name "*.mfront" | sort |
  grep -E "$PATTERN" | xargs -P "$JOBS" -I{} bash -c 'compile {}'
cat "$OUT"/*/status | sort > "$OUT/summary.tsv"
cut -f2 "$OUT/summary.tsv" | sort | uniq -c
