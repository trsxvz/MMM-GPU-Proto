#!/bin/bash
# Runs the behaviour test of a generic-parallel configuration with the
# development tree (tfel-dev: mfront of scratch/build-dev2 and headers of
# tfel-dev), without installing it. The substitutions of tfel-check are read
# in the cache of build/tfel.
# usage: scratch/dev-test.sh configuration[/sub-test] [pattern of compare-all.sh]
# (for instance cuda-clang or cuda-clang/tdls)
set -u
c=$1
conf=${c%%/*}
pattern=${2:-.}
P=$HOME/Documents/GitHub/mfem-mgis-gpu-proto
D=$P/tfel-dev
T=$P/scratch/devtest/$c
source "$P/scripts/env.sh" > /dev/null 2>&1
export MAKEFLAGS=${MAKEFLAGS:--j3}
rm -rf "$T" && mkdir -p "$(dirname "$T")"
cp -r "$D/test-suite/behaviours/generic-parallel" "$T"
cd "$T/$c" || exit 1
# substitutions of tfel-check, as KEY%value pairs
subst=$(grep '^tfel-substitutions:STRING=' "$P/build/tfel/CMakeCache.txt" | cut -d= -f2-)
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
# configuration file with the headers of the development tree first
apply "$(cat "$D/mfront/configuration-files/config-$conf.json")" |
  sed "s|compilation_flags: \({\?\)\"|compilation_flags: \1\"$dev_includes |g" > "config-$conf.json"
cat > mfront-dev.sh << EOF
#!/bin/bash
LD_LIBRARY_PATH=\$(find "$P/scratch/build-dev2" -name "*.so" -printf "%h\n" | sort -u | paste -sd:):\${LD_LIBRARY_PATH:-} \\
  exec "$P/scratch/build-dev2/mfront/src/mfront" "\$@"
EOF
chmod +x mfront-dev.sh
export LD_LIBRARY_PATH=src:$LD_LIBRARY_PATH
grep '^@Command' test.check | sed 's/^@Command "//; s/";$//' | while read -r cmd; do
  cmd=$(apply "$cmd")
  cmd=${cmd//@mfront@/./mfront-dev.sh}
  # headers of the development tree first (substituted by tfel-check)
  cmd=${cmd//-I @TFEL_INCLUDE_DIR@/$dev_includes -I $P/install/tfel/include}
  case "$cmd" in
    *compare-all.sh*) cmd="$cmd $pattern" ;;
  esac
  if ! sh -c "$cmd" > cmd.log 2>&1; then
    echo "$c: FAILED: ${cmd:0:120}"
    grep -E "FAILED|error|Error" cmd.log | head -5
    cat cmd.log >> all.log
    exit 1
  fi
  cat cmd.log >> all.log
done || exit 1
echo "$c: OK=$(grep -c 'OK$' all.log) FAILED=$(grep -c 'FAILED$' all.log) by-variable-OK=$(grep 'by variable' all.log | grep -c 'OK$')"
