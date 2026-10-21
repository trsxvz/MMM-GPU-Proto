# Comparisons of the `generic-parallel` interface to the `generic` interface
# for the test material properties, shared by the tests of all the backends
# (see compare.cxx).
#
# usage: sh ../compare-all.sh compare_executable backend [options] [pattern]
#
# The libraries are the ones built by build-libraries.sh in the `src`
# directory. The options (arguments starting with `--`, for instance
# `--tolerance=1e-13`) are passed to all the comparisons. The pattern (an
# extended regular expression, `.` by default) selects the compared
# material properties.
set -e
compare=$1
backend=$2
shift 2
pattern=.
options=""
for a in "$@"; do
  case "$a" in
    --*) options="$options $a" ;;
    *) pattern=$a ;;
  esac
done
run() {
  if echo "$2" | grep -qE "$pattern"; then
    library=$1
    function=$2
    shift 2
    "$compare" "src/libGeneric$library.so" \
      "src/libGenericParallel$library-$backend.so" "$function" $options "$@"
  fi
}
run Inconel600 Inconel600_YoungModulus
run Inconel600 Inconel600_YoungModulus_qt
run UO2 UO2_YoungModulus_Martin1989
run MaterialLaw PoissonRatioTest
run MaterialLaw LinearDataInterpolationTest
run MaterialLaw BoundsOrderTest
# the `generic` interface does not check the physical bounds of a material
# property called by another one, whereas the implementation of the called
# material property used in device code returns NaN if they are not
# satisfied: the temperatures are kept in the physical bounds of the called
# material property
run MaterialLaw MaterialLawCallTest --domain=0:1e300
