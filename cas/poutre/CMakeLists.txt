add_executable(poutre poutre.cxx)
target_link_libraries(poutre PRIVATE mfem-mgis::MFEMMGIS)
add_mfront_library(poutre-behaviours
  "${MFEM_MGIS_SOURCE_DIR}/tests/IsotropicLinearHardeningPlasticity.mfront")
configure_file("${MFEM_MGIS_SOURCE_DIR}/tests/beam-tet.mesh"
  "${CMAKE_CURRENT_BINARY_DIR}/beam-tet.mesh" COPYONLY)

# short tests: two steps on a coarse mesh; the three linear solvers give the
# same solution
set(norm 3.264041277)
add_test(NAME poutre-pcg-amg
  COMMAND mpirun -n 2 ./poutre -r 1 -n 2 -ls pcg-amg -en ${norm})
add_test(NAME poutre-mumps
  COMMAND mpirun -n 2 ./poutre -r 1 -n 2 -ls mumps -en ${norm})
add_test(NAME poutre-umfpack
  COMMAND ./poutre -r 1 -n 2 -ls umfpack -no-p -en ${norm})
# the comparison detects a wrong solution
add_test(NAME poutre-wrong-norm
  COMMAND ./poutre -r 1 -n 2 -ls umfpack -no-p -en 3.3)
set_tests_properties(poutre-wrong-norm PROPERTIES WILL_FAIL TRUE)
