cmake_minimum_required(VERSION 3.13)
project(SatohTestExample)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
# testing
# (must be placed *before* any add_subdirectory, cmake bug ?)
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

# Require MFEM-MGIS
if (NOT TARGET mfem-mgis::MFEMMGIS)
  find_package(MFEMMGIS REQUIRED)
endif(NOT TARGET mfem-mgis::MFEMMGIS)

# Define a test that will be called during the "make check"

if(MFEM_USE_SUITESPARSE)
  # Copy the files needed at run time. Unlike file(COPY), configure_file updates
  # the copies when these files are modified.
  foreach(input IsotropicLinearThermoElasticity.mfront square.msh)
    configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                   COPYONLY)
  endforeach()

  # Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU,
  # chosen at run time by mfem_mgis::getBehavioursLibrary from the
  # option --execution
  mfem_mgis_add_behaviours_library(BehavioursEx5
    SOURCES ${PROJECT_BINARY_DIR}/IsotropicLinearThermoElasticity.mfront)

  add_executable(satoh SatohTest.cxx)
  target_compile_definitions(satoh PRIVATE
    MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursEx5_BEHAVIOURS_DIRECTORY}")
  target_link_libraries(satoh mfem-mgis::MFEMMGIS)
  target_include_directories(satoh
    SYSTEM
    PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
    PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

  add_dependencies(check satoh BehavioursEx5)

  # The test is the same in the full and restricted test modes, since it is a
  # single computation on a coarse mesh, whose reference reaction depends on
  # the mesh (it is only checked with the original mesh and order: options
  # --refinement and --order).
  set(test "satoh")
  add_test(NAME ${test} COMMAND ${test})
  # the same computation in parallel, MUMPS replacing UMFPACK
  if(MFEM_USE_MUMPS AND MFEM_USE_MPI)
    add_test(NAME ${test}-parallel
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                     $<TARGET_FILE:satoh> --parallel)
    set_tests_properties(${test}-parallel PROPERTIES
                         PROCESSORS 2 ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
  # the same computation on a GPU, UMFPack computing on the CPU (see
  # ../config/TestMode.cmake)
  if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS AND BehavioursEx5_HAVE_GPU_BEHAVIOURS)
    add_test(NAME ${test}-gpu
             COMMAND ${test} --execution mixed --linear-solvers-execution cpu)
    set_tests_properties(${test}-gpu PROPERTIES RESOURCE_LOCK gpu)
  endif()
else()
  message(STATUS "Disabling satoh. Compile MFEM with SuiteSparse support to enable it")
endif()


