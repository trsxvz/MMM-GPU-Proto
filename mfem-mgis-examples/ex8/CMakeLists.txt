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

set(MFRONT_OUTPUTS
    "${CMAKE_CURRENT_BINARY_DIR}/src/libU3SI2-generic.so"
    "${CMAKE_CURRENT_BINARY_DIR}/src/libALFENI-generic.so"
)

# Compiling MFront laws
add_custom_command(
    OUTPUT ${MFRONT_OUTPUTS}
    COMMAND mfront --obuild --interface=generic --include=${MFRONT_INCLUDE_DIR} ${MFRONT_LOIS}
    DEPENDS ${MFRONT_LOIS} ${MFRONT_HEADERS}
    COMMENT "Generating laws for thermomechanical simulation."
    VERBATIM
)

add_executable(rjh_plate
    Thermomechanical.cxx
    src/BoundaryConditions.cxx
    src/RobinBC.cxx
)

# Copy the files needed at run time next to the executable. Unlike
# file(COPY), configure_file updates the copies when these files are modified.
foreach(input mesh/assemblage_hexa.msh mesh/assemblage_hexa_coarse.msh
              assemblage_hexa-statistics.ref
              assemblage_hexa_coarse-statistics.ref mtest/Swelling.mtest)
  get_filename_component(name ${input} NAME)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${input}
                 ${CMAKE_CURRENT_BINARY_DIR}/${name} COPYONLY)
endforeach()

target_include_directories(rjh_plate PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/headers")

add_custom_target(
    generate-libBehaviour ALL
    DEPENDS ${MFRONT_OUTPUTS}
)

target_link_libraries(rjh_plate PRIVATE
    MPI::MPI_CXX
    mfem-mgis::MFEMMGIS
)

add_dependencies(check rjh_plate generate-libBehaviour)

# The full test mode runs the simulation with the default mesh. The restricted
# test mode uses the coarser mesh. The statistics of the temperature and of the
# displacement are compared to the reference values computed with the same
# mesh. The swelling is always compared to its exact value.
if(MFEM_MGIS_EXAMPLES_FULL_TESTS)
  set(rjh_plate_mesh assemblage_hexa)
else()
  set(rjh_plate_mesh assemblage_hexa_coarse)
endif()
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
# the mechanical problem solved by the NewtonK0 method (tangent systems
# solved by GMRES preconditioned by the elastic operator, factorized by
# MUMPS at each time step), all the computations being made by kernels
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
# the same computations on a GPU (see ../config/TestMode.cmake)
if(MFEM_MGIS_EXAMPLES_CUDA_TESTS)
  foreach(_variant "" "-kernels")
    add_test(NAME rjh_plate${_variant}-cuda
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:rjh_plate>
            -m ${rjh_plate_mesh}.msh
            -lsMc MUMPSSolver
            -rf ${rjh_plate_mesh}-statistics.ref
            --device cuda
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(rjh_plate${_variant}-cuda PROPERTIES
        PROCESSORS 2
        RESOURCE_LOCK gpu)
  endforeach()
  set_tests_properties(rjh_plate-cuda PROPERTIES
      ENVIRONMENT OMP_NUM_THREADS=1)
  set_tests_properties(rjh_plate-kernels-cuda PROPERTIES
      ENVIRONMENT "${rjh_plate_kernels_environment}")
endif()

# comparison of the swelling model to its exact value
add_test(NAME u3si2_swelling
    COMMAND mtest Swelling.mtest
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# comparison of the Robin boundary condition to an exact solution
add_subdirectory(Robin)
