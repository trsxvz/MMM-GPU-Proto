# Generation and compilation of the libraries compared by the tests of a
# configuration, shared by the tests of all the configurations.
#
# usage: sh ../build-libraries.sh mfront configuration_file [behaviours]
#
# - the reference library of the `generic` interface is built in the
#   `reference` directory with the default compiler of MFront, or with the
#   configuration file given by the MFRONT_GPB_REFERENCE_CONFIGURATION_FILE
#   environment variable. The reference then does not depend on the tested
#   compiler: the compilers translating host and device code in a single
#   pass (nvc++, AdaptiveCpp) define `TFEL_DEVICE_CODE` in host code too;
# - the libraries of the `generic` and `generic-parallel` interfaces are built
#   in the current directory with the configuration file of the tested
#   configuration. Using both interfaces in the same call checks that they
#   can be used together.
#
# The behaviours (names of the MFront files of the directory of this script,
# without extension) are all the test behaviours by default.
set -e
mfront=$1
configuration=$2
shift 2
case "$mfront" in
  /*) ;;
  */*) mfront=$PWD/$mfront ;;
esac
directory=$(cd "$(dirname "$0")" && pwd)
behaviours=${*:-"
  IsotropicLinearHardeningPlasticity ImplicitNorton
  HosfordPerfectPlasticity StrainHardeningCreep3 HenckyImplicitNorton
  GreenLagrangeImplicitNorton SaintVenantKirchhoffElasticity
  DefaultGenericSaintVenantKirchhoffElasticity IsotropicMisesCreepTest
  MultipleIsotropicMisesFlowsTest ImplicitGenericBehaviourTest
  ImplicitModelTest ParametersAsStaticVariablesTest
  InitializeFunctionTest PostProcessingTest OrthotropicElastic
  OrthotropicSaintVenantKirchhoffElasticity
  OrthotropicStationaryHeatTransfer SingleCrystal_DD_FCC
  FiniteStrainSingleCrystal ImplicitNortonWithEnergies Norton_rk42
  Tvergaard BehaviourVariableTest ExternalModelCallerTest
  StandardElastoViscoPlasticityTest
  OrthotropicStandardElastoViscoPlasticityTest
  PorousStandardElastoViscoPlasticityTest GenericParallelModelTest
  GenericParallelModelWithQuantitiesTest
"}
files=$(for b in $behaviours; do echo "$directory/$b.mfront"; done)
reference_configuration=""
if [ -n "${MFRONT_GPB_REFERENCE_CONFIGURATION_FILE:-}" ]; then
  reference_configuration="--configuration-file=$MFRONT_GPB_REFERENCE_CONFIGURATION_FILE"
fi
mkdir -p reference
(cd reference &&
  "$mfront" --obuild --search-path="$directory" --interface=generic \
    $reference_configuration $files)
"$mfront" --obuild --search-path="$directory" --interface=generic \
  --interface=generic-parallel --configuration-file="$configuration" $files
