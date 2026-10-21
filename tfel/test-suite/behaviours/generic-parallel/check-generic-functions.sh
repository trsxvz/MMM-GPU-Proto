# Check the libraries of the `generic-parallel` interface which also contain
# the functions of the `generic` interface (`generic_functions` option):
# libraries of a few test behaviours are built in the `generic-functions`
# directory with the tested configuration and this option, and compared to
# the reference built by build-libraries.sh, with the
# `--generic-functions=true` option of compare.cxx (interface `Generic`,
# functions of a point compared to the functions of batches of points,
# shared parameters and out of bounds policy).
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
for b in IsotropicLinearHardeningPlasticity ImplicitNortonWithEnergies \
  OrthotropicElastic InitializeFunctionTest PostProcessingTest \
  GenericParallelModelTest; do
  files="$files $directory/$b.mfront"
done
"$mfront" --obuild --search-path="$directory" --interface=generic-parallel \
  --configuration-file="$configuration" \
  --configuration-file="$PWD/config-generic-functions.json" $files
# a library has its own out of bounds policies: the functions giving them are
# not exported as unique symbols of the process, which would be shared with
# the libraries of other backends (ELF libraries)
if [ "$(uname)" = Linux ] && command -v nm > /dev/null 2>&1; then
  for l in src/libGenericParallel*.so; do
    if nm -DC "$l" | awk '$2 == "u"' | grep -q '_getOutOfBoundsPolicy'; then
      echo "an out of bounds policy of '$l' is shared by the process" >&2
      exit 1
    fi
  done
fi
gf="--generic-functions=true $*"
# compare-functions.sh reads the comparison program and the batch library in
# the positional parameters (`.` does not take arguments in POSIX shells)
set -- "$compare" "src/libGenericParallelBehaviour-$backend.so"
. "$directory/compare-functions.sh"
# reference built by build-libraries.sh
reference=../reference/src/libBehaviour.so
run IsotropicLinearHardeningPlasticity_Tridimensional --amplitude=3e-3 $gf
run IsotropicLinearHardeningPlasticity_Tridimensional --amplitude=3e-3 \
  --parameters=YoungModulus=150e9,PoissonRatio=0.3,HardeningSlope=5e9,YieldStrength=200e6 $gf
run ImplicitNortonWithEnergies_Tridimensional --material-properties=150e9,0.3 \
  --perturbation-of-internal-state-variables=1e-4 --out-of-bounds-policy=Strict $gf
run OrthotropicElastic_Tridimensional \
  --material-properties=200e9,150e9,100e9,0.3,0.25,0.2,70e9,60e9,50e9 $gf
run OrthotropicElastic_Tridimensional_rotateGradients $gf
run InitializeFunctionTest_Tridimensional_InitializeFunction_StressFromInitialPressure $gf
run PostProcessingTest_Tridimensional_PostProcessing_PrincipalStrain \
  --material-properties=150e9,0.3 $gf
# standalone model: failures at some points, parameter changed in the batch
# library
reference=../reference/src/libModel.so
batch=src/libGenericParallelModel-$backend.so
run GenericParallelModelTest_PlaneStress --material-properties=500 \
  --amplitude=0 --temperature=600 \
  --perturbation-of-internal-state-variables=1e-4 $gf
run GenericParallelModelTest_Tridimensional --material-properties=500 \
  --amplitude=0 --temperature=600 --parameters=SwellingCoefficient=1e-3 $gf
exit $status
