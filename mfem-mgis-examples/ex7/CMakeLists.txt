cmake_minimum_required(VERSION 3.13)
project(MixedOxideFuelsExample)

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

# Test mode of the examples, see ../config/TestMode.cmake
if(NOT DEFINED MFEM_MGIS_EXAMPLES_FULL_TESTS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/../config/TestMode.cmake)
endif()

# Define a test that will be called during the "make check"

# Copy the files needed at run time. Unlike file(COPY), configure_file updates
# the copies when these files are modified.
# The contents of the directories are listed again when files are added
file(GLOB_RECURSE directory_files CONFIGURE_DEPENDS
     RELATIVE ${PROJECT_SOURCE_DIR}
     ${PROJECT_SOURCE_DIR}/mesh/*
     ${PROJECT_SOURCE_DIR}/file_geo/*
     ${PROJECT_SOURCE_DIR}/script_merope/*)
foreach(input NortonViscoplasticityWithThreshold.mfront petscrc
              OneSphere-avgStress-5steps.ref OneSphere-avgStress-40steps.ref
              ${directory_files})
  configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                 COPYONLY)
endforeach()

# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU, chosen
# at run time by mfem_mgis::getBehavioursLibrary from the option --execution
mfem_mgis_add_behaviours_library(BehavioursEx7
  SOURCES ${PROJECT_BINARY_DIR}/NortonViscoplasticityWithThreshold.mfront)

add_executable(mox2 MixedOxideFuels.cxx)
target_compile_definitions(mox2 PRIVATE
  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursEx7_BEHAVIOURS_DIRECTORY}")
target_link_libraries(mox2 mfem-mgis::MFEMMGIS)
target_include_directories(mox2
  SYSTEM
  PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
  PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

add_dependencies(check mox2 BehavioursEx7)

# The full test mode computes the complete loading in 40 time steps. The
# restricted test mode computes it in 5 time steps only. The mean stresses are
# compared to the reference values computed with the same number of time
# steps.
if(MFEM_MGIS_EXAMPLES_FULL_TESTS)
  set(mox2_nbsteps 40)
else()
  set(mox2_nbsteps 5)
endif()
set(mox2_test_options --no-post-processing --order 2 --mesh mesh/OneSphere.msh
    --nbsteps ${mox2_nbsteps}
    --reference-file OneSphere-avgStress-${mox2_nbsteps}steps.ref)

set(test "mox2")
add_test(NAME ${test} COMMAND mox2 ${mox2_test_options})

# a direct solver (MUMPS, symmetric positive definite matrices) gives the
# same mean stresses
if(MFEM_USE_MUMPS)
  add_test(NAME ${test}-mumps
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:mox2> ${mox2_test_options}
                   -ls MUMPSSolver -ms SymmetricPositiveDefinite)
  set_tests_properties(${test}-mumps PROPERTIES
                       RESOURCE_LOCK avgStress PROCESSORS 2
                       ENVIRONMENT OMP_NUM_THREADS=1)
endif()

if(MFEM_USE_PETSC)
  add_test(NAME ${test}-petsc
           COMMAND mox2 ${mox2_test_options}
                   --use-petsc --petsc-configuration-file petscrc)
  # both tests write the mean stresses in the same file
  set_tests_properties(${test} ${test}-petsc PROPERTIES RESOURCE_LOCK avgStress)
endif()

# the same computation on a GPU, and with the behaviours integrated on the
# CPU (see ../config/TestMode.cmake)
if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS AND BehavioursEx7_HAVE_GPU_BEHAVIOURS)
  add_test(NAME ${test}-gpu
           COMMAND mox2 ${mox2_test_options} --execution gpu)
  add_test(NAME ${test}-mixed-gpu
           COMMAND mox2 ${mox2_test_options} --execution mixed
                   --behaviours-execution cpu)
  set_tests_properties(${test}-gpu ${test}-mixed-gpu PROPERTIES
                       RESOURCE_LOCK "avgStress;gpu")
endif()
