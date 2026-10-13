cmake_minimum_required(VERSION 3.10)
project(ProjetSimuThermoMeca VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MFEMMGIS REQUIRED)

# Test mode of the examples, see ../config/TestMode.cmake
if(NOT DEFINED MFEM_MGIS_EXAMPLES_FULL_TESTS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/../config/TestMode.cmake)
endif()

include(CTest)
enable_testing()
if(CMAKE_CONFIGURATION_TYPES)
if (NOT TARGET check)
add_custom_target(check COMMAND
${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
if (NOT TARGET check)
add_custom_target(check COMMAND
${CMAKE_CTEST_COMMAND} -T test )
endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)


get_filename_component(MFRONT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MFront" ABSOLUTE)

file(GLOB MFRONT_LOIS CONFIGURE_DEPENDS "${MFRONT_INCLUDE_DIR}/*.mfront")
# headers included by the MFront files
file(GLOB MFRONT_HEADERS CONFIGURE_DEPENDS "${MFRONT_INCLUDE_DIR}/*.hxx")

# Libraries of the U3Si2 and ALFENI materials for the CPU and, if MFEM uses
# a GPU, for the GPU, chosen at run time by mfem_mgis::getBehavioursLibrary
# from the option --execution
mfem_mgis_add_behaviours_library(BehavioursEx8
    SOURCES ${MFRONT_LOIS}
    LIBRARIES U3SI2 ALFENI
    INCLUDE_DIRECTORY ${MFRONT_INCLUDE_DIR}
    DEPENDS ${MFRONT_HEADERS})

add_executable(rjh_plate
    Thermomechanical.cxx
    src/BoundaryConditions.cxx
    src/RobinBC.cxx
)

# Copy the files needed at run time next to the executable. Unlike
# file(COPY), configure_file updates the copies when these files are modified.
foreach(input mesh/assemblage_hexa.msh mesh/assemblage_hexa_coarse.msh
              assemblage_hexa-statistics.ref
              assemblage_hexa_coarse-statistics.ref)
  get_filename_component(name ${input} NAME)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${input}
                 ${CMAKE_CURRENT_BINARY_DIR}/${name} COPYONLY)
endforeach()

target_include_directories(rjh_plate PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/headers")
target_compile_definitions(rjh_plate PRIVATE
    MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursEx8_BEHAVIOURS_DIRECTORY}")

target_link_libraries(rjh_plate PRIVATE
    MPI::MPI_CXX
    mfem-mgis::MFEMMGIS
)

add_dependencies(check rjh_plate BehavioursEx8)

# The full test mode runs the simulation with the default mesh. The restricted
# test mode uses the coarser mesh. The statistics of the temperature and of the
# displacement are compared to the reference values computed with the same
# mesh. The swelling is always compared to its exact value.
if(MFEM_MGIS_EXAMPLES_FULL_TESTS)
  set(rjh_plate_mesh assemblage_hexa)
else()
  set(rjh_plate_mesh assemblage_hexa_coarse)
endif()
# the mechanical problem is solved by MUMPS
if(MFEM_USE_MUMPS)
  add_test(NAME rjh_plate
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          $<TARGET_FILE:rjh_plate>
          -m ${rjh_plate_mesh}.msh
          -lsMc MUMPSSolver
          -rf ${rjh_plate_mesh}-statistics.ref
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  # one OpenMP thread per process, otherwise MUMPS oversubscribes the cores
  set_tests_properties(rjh_plate PROPERTIES
      PROCESSORS 2
      ENVIRONMENT OMP_NUM_THREADS=1)
  # the same computation by batches and by kernels (integration of the
  # behaviours by batches, gradients, inner forces and tangent operator
  # computed by kernels, the Robin boundary condition being added to them),
  # compared to the same reference
  add_test(NAME rjh_plate-kernels
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          $<TARGET_FILE:rjh_plate>
          -m ${rjh_plate_mesh}.msh
          -lsMc MUMPSSolver
          -rf ${rjh_plate_mesh}-statistics.ref
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set(rjh_plate_kernels_environment
      "OMP_NUM_THREADS=1;MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch;MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels;MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels;MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels")
  set_tests_properties(rjh_plate-kernels PROPERTIES
      PROCESSORS 2
      ENVIRONMENT "${rjh_plate_kernels_environment}")
  # the tangent operator of the mechanical problem given to MUMPS as a
  # symmetric positive definite matrix, compared to the same reference
  add_test(NAME rjh_plate-mumps-spd
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          $<TARGET_FILE:rjh_plate>
          -m ${rjh_plate_mesh}.msh
          -lsMc MUMPSSolver
          -msMc SymmetricPositiveDefinite
          -rf ${rjh_plate_mesh}-statistics.ref
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(rjh_plate-mumps-spd PROPERTIES
      PROCESSORS 2
      ENVIRONMENT "${rjh_plate_kernels_environment}")
  # the mechanical problem solved by the NewtonK0 method (tangent systems
  # solved by GMRES preconditioned by the elastic operator, factorized by
  # MUMPS at each resolution of the mechanical problem), all the computations
  # being made by kernels
  add_test(NAME rjh_plate-newton-k0
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          $<TARGET_FILE:rjh_plate>
          -m ${rjh_plate_mesh}.msh
          -lsMc MUMPSSolver
          -nlsMc NewtonK0
          -rf ${rjh_plate_mesh}-statistics.ref
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(rjh_plate-newton-k0 PROPERTIES
      PROCESSORS 2
      ENVIRONMENT "${rjh_plate_kernels_environment}")
  # the same computation on a GPU, the direct solver MUMPS computing on the
  # CPU (option --execution mixed --linear-solver-device cpu), and with the
  # finite elements computed on the CPU too (see ../config/TestMode.cmake)
  if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS AND BehavioursEx8_HAVE_GPU_BEHAVIOURS)
    add_test(NAME rjh_plate-gpu
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:rjh_plate>
            -m ${rjh_plate_mesh}.msh
            -lsMc MUMPSSolver
            -rf ${rjh_plate_mesh}-statistics.ref
            --execution mixed --linear-solver-device cpu
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME rjh_plate-mixed-gpu
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:rjh_plate>
            -m ${rjh_plate_mesh}.msh
            -lsMc MUMPSSolver
            -rf ${rjh_plate_mesh}-statistics.ref
            --execution mixed --linear-solver-device cpu
            --finite-elements cpu
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(rjh_plate-gpu rjh_plate-mixed-gpu PROPERTIES
        PROCESSORS 2
        RESOURCE_LOCK gpu
        ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
endif()

# comparison of the swelling model to its exact value. MTest computes in
# double: with another numerical type, the behaviours are also generated in
# double for this test (a mfem-mgis without numerical type uses double).
if((NOT BehavioursEx8_REAL_TYPE) OR (BehavioursEx8_REAL_TYPE STREQUAL "double"))
  set(U3SI2_MTEST_LIBRARY
      "BehavioursEx8/cpu/libU3SI2${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  set(_mtest_directory "${CMAKE_CURRENT_BINARY_DIR}/BehavioursEx8-mtest")
  set(U3SI2_MTEST_LIBRARY
      "${_mtest_directory}/src/libU3SI2-generic${CMAKE_SHARED_LIBRARY_SUFFIX}")
  file(MAKE_DIRECTORY "${_mtest_directory}")
  add_custom_command(
    OUTPUT "${U3SI2_MTEST_LIBRARY}"
    COMMAND "${MFRONT}" --obuild --interface=generic
            "--include=${MFRONT_INCLUDE_DIR}" ${MFRONT_LOIS}
    WORKING_DIRECTORY "${_mtest_directory}"
    DEPENDS ${MFRONT_LOIS} ${MFRONT_HEADERS}
    COMMENT "MFront: behaviours of ex8 in double for MTest")
  add_custom_target(BehavioursEx8-mtest ALL DEPENDS "${U3SI2_MTEST_LIBRARY}")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mtest/Swelling.mtest.in
               ${CMAKE_CURRENT_BINARY_DIR}/Swelling.mtest @ONLY)
add_test(NAME u3si2_swelling
    COMMAND mtest Swelling.mtest
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# comparison of the Robin boundary condition to an exact solution
add_subdirectory(Robin)
