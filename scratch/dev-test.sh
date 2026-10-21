#!/bin/bash
# Runs the test of a generic-parallel configuration with the
# development tree (tfel-dev: mfront of scratch/build-dev2 and headers of
# tfel-dev), without installing it. The substitutions of tfel-check are read
# in the cache of the build directory of TFEL (see scripts/env.sh).
# usage: scratch/dev-test.sh configuration[/sub-test] [pattern of compare-all.sh]
# (for instance cuda-clang or cuda-clang/tdls). DEVTEST_KIND selects the
# tests: behaviours (default) or material-properties. DEVTEST_OPTIMIZATION_LEVEL
# (for instance -O1 or -O3) replaces -O2 in the configuration file and in the
# commands of the test. DEVTEST_BEHAVIOURS (names of the test behaviours)
# restricts the built behaviours, for targeted tests. DEVTEST_TREE and
# DEVTEST_BUILD select another development tree and its build directory.
# The tests are copied with the layout of the installed test-suite
# (<kind>/generic-parallel/<configuration>): the comparison of the material
# properties includes a header of the tests of the behaviours.
set -u
c=$1
conf=${c%%/*}
pattern=${2:-.}
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
D=${DEVTEST_TREE:-$P/tfel-dev}
BUILD=${DEVTEST_BUILD:-$P/scratch/build-dev2}
kind=${DEVTEST_KIND:-behaviours}
case "$kind" in
  behaviours | material-properties) ;;
  *)
    echo "invalid kind of tests '$kind'" >&2
    exit 1
    ;;
esac
# the copy is deleted with a literal prefix, before sourcing scripts/env.sh
# (the slashes of the configuration are replaced)
rm -rf -- "$HOME/Documents/GitHub/mfem-mgis-gpu-proto/scratch/devtest/$kind-${c//\//-}"
T=$P/scratch/devtest/$kind-${c//\//-}
source "$P/scripts/env.sh" > /dev/null 2>&1
export MAKEFLAGS=${MAKEFLAGS:--j3}
mkdir -p "$T/behaviours" "$T/material-properties"
cp -r "$D/test-suite/behaviours/generic-parallel" "$T/behaviours/"
cp -r "$D/test-suite/material-properties/generic-parallel" "$T/material-properties/"
cd "$T/$kind/generic-parallel/$c" || exit 1
# substitutions of tfel-check, as KEY%value pairs
subst=$(grep '^tfel-substitutions:STRING=' "$PROTO_TFEL_BUILD/CMakeCache.txt" | cut -d= -f2-)
# additional substitutions, as KEY%value pairs separated by semicolons
subst="$subst${DEVTEST_SUBSTITUTIONS:+;$DEVTEST_SUBSTITUTIONS}"
apply() {
  local line=$1 pair
  IFS=';' read -ra pairs <<< "$subst"
  for pair in "${pairs[@]}"; do
    line=${line//@${pair%%%*}@/${pair#*%}}
  done
  echo "$line"
}
dev_includes="-I$D/mfront/include -I$D/include"
# configuration files with the headers of the development tree first: the
# one of the configuration and the ones used by the test (for instance,
# sycl-gpu-acpp uses config-sycl-default-acpp.json)
configurations=$( (echo "$conf"; grep -o 'config-[A-Za-z0-9_-]*\.json' test.check |
  sed 's/^config-//; s/\.json$//') | sort -u)
for cf in $configurations; do
  [ -f "$D/mfront/configuration-files/config-$cf.json" ] || continue
  apply "$(cat "$D/mfront/configuration-files/config-$cf.json")" |
    sed "s|compilation_flags: \({\?\)\"|compilation_flags: \1\"$dev_includes |g" > "config-$cf.json"
done
level=${DEVTEST_OPTIMIZATION_LEVEL:-}
if [ -n "$level" ]; then
  sed -i "s/-O2/$level/g" config-*.json test.check
fi
# the headers of the development tree are also used by the sources compiled
# by the C++ compiler (INCLUDES precedes the other include paths)
cat > mfront-dev.sh << EOF
#!/bin/bash
INCLUDES="$dev_includes \${INCLUDES:-}" \\
LD_LIBRARY_PATH=\$(find "$BUILD" -name "*.so" -printf "%h\n" | sort -u | paste -sd:):\${LD_LIBRARY_PATH:-} \\
  exec "$BUILD/mfront/src/mfront" "\$@"
EOF
chmod +x mfront-dev.sh
# reference library built with the headers of the development tree too (see
# build-libraries.sh)
cat > config-reference.json << EOF
compilation_options : {
  cxx : {
    compilation_flags: "$dev_includes -O2 -std=c++20"
  }
}
EOF
export MFRONT_GPB_REFERENCE_CONFIGURATION_FILE=$PWD/config-reference.json
export LD_LIBRARY_PATH=src:reference/src:$LD_LIBRARY_PATH
grep '^@Command' test.check | sed 's/^@Command "//; s/";$//' | while read -r cmd; do
  cmd=$(apply "$cmd")
  cmd=${cmd//@mfront@/./mfront-dev.sh}
  # headers of the development tree first (substituted by tfel-check)
  cmd=${cmd//-I @TFEL_INCLUDE_DIR@/$dev_includes -I$PROTO_TFEL_INSTALL/include}
  cmd=${cmd//-I@TFEL_INCLUDE_DIR@/$dev_includes -I$PROTO_TFEL_INSTALL/include}
  cmd=${cmd//@TFEL_LIBRARY_DIR@/$PROTO_TFEL_INSTALL/lib}
  case "$cmd" in
    *compare-all.sh*) cmd="$cmd '$pattern'" ;;
    # DEVTEST_BEHAVIOURS restricts the built behaviours (all by default)
    *build-libraries.sh*) cmd="$cmd ${DEVTEST_BEHAVIOURS:-}" ;;
    # the comparisons called directly by the test may use behaviours which
    # are not built in a targeted test
    ./compare*) [ -n "${DEVTEST_BEHAVIOURS:-}" ] && continue ;;
  esac
  if ! sh -c "$cmd" > cmd.log 2>&1; then
    echo "$c: FAILED: ${cmd:0:120}"
    grep -E "FAILED|error|Error" cmd.log | head -5
    cat cmd.log >> all.log
    exit 1
  fi
  cat cmd.log >> all.log
done || exit 1
if [ "$kind" = behaviours ]; then
  echo "$c: OK=$(grep -c 'OK$' all.log) FAILED=$(grep -c 'FAILED$' all.log) by-variable-OK=$(grep 'by variable' all.log | grep -c 'OK$')"
else
  echo "$c: $kind OK"
fi
