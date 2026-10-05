# Comparisons of the `generic-parallel` interface to the `generic` interface
# for behaviours using the TDLS linear system solver, shared by the tests of
# all the backends.
#
# usage: sh ../../compare-tdls.sh compare_executable batch_library [pattern]
#
# See compare-functions.sh for the arguments.
. "$(dirname "$0")/compare-functions.sh"
# consistent tangent operator computed from the inverse of the jacobian
run ImplicitNortonTDLS_Tridimensional 6 6 2 7 1 36 1e-3 150e9,0.3 293.15
run ImplicitNortonTDLS_PlaneStress 4 4 2 6 1 16 1e-3 150e9,0.3 293.15
exit $status
