# Check the libraries of the `generic-parallel` interface which also contain
# the functions of the `generic` interface (`generic_functions` option): the
# libraries of the test material properties are built in the
# `generic-functions` directory with the tested configuration and this
# option, together with the libraries of the `generic` interface, and
# compared by compare-all.sh with the `--generic-functions=true` option of
# compare.cxx.
#
# usage: sh ../check-generic-functions.sh mfront configuration_file
#        compare_executable backend [options]
set -e
mfront=$1
configuration=$2
compare=$3
backend=$4
shift 4
for v in mfront compare; do
  eval "p=\$$v"
  case "$p" in
    /*) ;;
    */*) eval "$v=\$PWD/\$p" ;;
  esac
done
# the configuration file is looked for in the current directory, then in
# the configuration files installed with TFEL, as MFront does
find_configuration() {
  if [ -f "$1" ]; then
    echo "$PWD/$1"
    return
  fi
  prefix=$(cd "$(dirname "$(command -v "$mfront")")/.." && pwd)
  for d in "$prefix"/share/tfel*/mfront/configuration-files; do
    if [ -f "$d/$1" ]; then
      echo "$d/$1"
      return
    fi
  done
  echo "configuration file '$1' not found" >&2
  exit 1
}
configuration=$(find_configuration "$configuration")
directory=$(cd "$(dirname "$0")" && pwd)
mkdir -p generic-functions
cd generic-functions
cat > config-generic-functions.json << 'END'
interfaces_options: {
  generic-parallel: {
    generic_functions: true
  }
}
END
files=""
for m in PoissonRatioTest Inconel600_YoungModulus Inconel600_YoungModulus-qt \
  UO2_YoungModulus_Martin1989 LinearDataInterpolationTest MaterialLawCallTest \
  BoundsOrderTest; do
  files="$files $directory/$m.mfront"
done
"$mfront" --obuild --search-path="$directory" --interface=generic \
  --interface=generic-parallel --configuration-file="$configuration" \
  --configuration-file="$PWD/config-generic-functions.json" $files
# a library of this interface has its own parameters: the parameters
# handler is not exported as a unique symbol of the process, which would be
# shared with the libraries of other backends (ELF libraries)
if [ "$(uname)" = Linux ] && command -v nm > /dev/null 2>&1; then
  for l in src/libGenericParallel*.so; do
    if nm -DC "$l" | awk '$2 == "u"' | grep -q 'ParametersHandler'; then
      echo "the parameters handler of '$l' is shared by the process" >&2
      exit 1
    fi
  done
fi
sh "$directory/compare-all.sh" "$compare" "$backend" --generic-functions=true "$@"
# runtime checks disabled, as with the `generic` interface (no check of the
# number of arguments, of the parameters, of the bounds and of the results).
# Some configuration files set this option to `false`, which conflicts with
# the option given on the command line: a copy of the configuration file
# setting it to `true` is then used (see build-no-runtime-checks.sh).
mkdir -p no-runtime-checks
cd no-runtime-checks
if grep -q 'disable_runtime_checks' "$configuration"; then
  sed 's/disable_runtime_checks: *false/disable_runtime_checks: true/' \
    "$configuration" > no-runtime-checks.json
  option=""
else
  cp "$configuration" no-runtime-checks.json
  option="--dsl-option=disable_runtime_checks:true"
fi
"$mfront" --obuild --search-path="$directory" --interface=generic \
  --interface=generic-parallel --configuration-file=no-runtime-checks.json \
  --configuration-file="$PWD/../config-generic-functions.json" $option \
  "$directory/Inconel600_YoungModulus.mfront"
"$compare" src/libGenericInconel600.so \
  "src/libGenericParallelInconel600-$backend.so" Inconel600_YoungModulus \
  --generic-functions=true "$@"
