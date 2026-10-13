# Copy the law. Unlike file(COPY), configure_file updates the copy when the law
# is modified.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/StationaryLinearHeatTransfer.mfront
               ${CMAKE_CURRENT_BINARY_DIR}/StationaryLinearHeatTransfer.mfront
               COPYONLY)

# Library of the law for the CPU and, if MFEM uses a GPU, for the GPU, chosen
# at run time by mfem_mgis::getBehavioursLibrary from the option --execution
mfem_mgis_add_behaviours_library(BehavioursRobin
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/StationaryLinearHeatTransfer.mfront)

add_executable(robin_test robin_test.cxx ../src/RobinBC.cxx)
target_link_libraries(robin_test PRIVATE MPI::MPI_CXX mfem-mgis::MFEMMGIS)
target_compile_definitions(robin_test PRIVATE
    MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursRobin_BEHAVIOURS_DIRECTORY}")

add_dependencies(check robin_test BehavioursRobin)

add_test(NAME robin_test
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        $<TARGET_FILE:robin_test>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(robin_test PROPERTIES PROCESSORS 2)
# the same computation on a GPU (see ../../config/TestMode.cmake)
if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS)
  add_test(NAME robin_test-gpu
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          $<TARGET_FILE:robin_test> --execution gpu
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(robin_test-gpu PROPERTIES
      PROCESSORS 2 RESOURCE_LOCK gpu)
endif()
