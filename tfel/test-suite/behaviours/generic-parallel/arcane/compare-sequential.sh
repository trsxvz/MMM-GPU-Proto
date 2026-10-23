# Comparisons of the `generic-parallel` interface to the `generic` interface
# on a `Sequential` queue of Arccore, restricted to an integration, a
# rotation, an initialize function, a post-processing and shared data (see
# compare-all.sh and the `MFRONT_GPB_COMPARE_ARCANE_QUEUE` environment
# variable of runtime.hxx).
#
# usage: sh compare-sequential.sh compare_executable batch_library
MFRONT_GPB_COMPARE_ARCANE_QUEUE=sequential sh "$(dirname "$0")/../compare-all.sh" "$1" "$2" \
  '^(ImplicitNorton_Tridimensional|PolyCrystalSharedDataTest|OrthotropicElastic_Tridimensional_rotateGradients|InitializeFunctionTest_Tridimensional|ImplicitNorton_PlaneStress_PostProcessing)'
