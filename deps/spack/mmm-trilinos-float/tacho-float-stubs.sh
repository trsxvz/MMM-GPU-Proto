#!/bin/bash
# Library of stubs of the single precision instantiations of
# Tacho::Driver: Amesos2 of the float build of Trilinos 17.2.1 refers to
# them, but Tacho is only instantiated in double, so that the programs
# linked with Amesos2 do not link. Each stub aborts with a message: Tacho
# is not used in single precision.
# usage: tacho-float-stubs.sh libamesos2.so output-directory
set -eu
amesos2=${1:?usage: tacho-float-stubs.sh libamesos2.so output-directory}
output=${2:?usage: tacho-float-stubs.sh libamesos2.so output-directory}
mkdir -p "$output"
{
  echo '#include <stdio.h>'
  echo '#include <stdlib.h>'
  echo 'static void tacho_float_unavailable(void) {'
  echo '  fputs("Tacho is not available in single precision\n", stderr);'
  echo '  abort();'
  echo '}'
  i=0
  for s in $(nm -D --undefined-only "$amesos2" | awk '{print $2}' |
             /usr/bin/grep '^_ZN5Tacho6DriverIf'); do
    echo "void tacho_float_stub_$i(void) __asm__(\"$s\");"
    echo "void tacho_float_stub_$i(void) { tacho_float_unavailable(); }"
    i=$((i + 1))
  done
} > "$output/tacho-float-stubs.c"
gcc -O2 -fPIC -shared -o "$output/libtacho-float-stubs.so" \
  "$output/tacho-float-stubs.c"
echo "$output/libtacho-float-stubs.so"
