cmake_minimum_required(VERSION 3.13)
project(UniaxialTensileTestExample)

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

# The sources of this example are installed by mfem-mgis
set(MFEMMGIS_EX1_DIR "${MFEMMGIS_DIR}/../examples/ex1")

# Define the executable to build for this example
add_executable(uniaxial_tensile_test ${MFEMMGIS_EX1_DIR}/UniaxialTensileTest.cxx)
target_link_libraries(uniaxial_tensile_test mfem-mgis::MFEMMGIS)
target_include_directories(uniaxial_tensile_test
    SYSTEM
    PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
    PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

# Copy the files needed at run time. Unlike file(COPY), configure_file updates
# the copies when these files are modified.
foreach(input IsotropicLinearHardeningPlasticity.mfront cube.mesh
              IsotropicLinearHardeningPlasticity.ref)
  configure_file(${MFEMMGIS_EX1_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                 COPYONLY)
endforeach()

# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU,
# chosen at run time by mfem_mgis::getBehavioursLibrary from the
# option --execution
mfem_mgis_add_behaviours_library(BehavioursEx1
  SOURCES ${PROJECT_BINARY_DIR}/IsotropicLinearHardeningPlasticity.mfront)
target_compile_definitions(uniaxial_tensile_test PRIVATE
  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursEx1_BEHAVIOURS_DIRECTORY}")

add_dependencies(check uniaxial_tensile_test BehavioursEx1)

# Define the test that will be called during the "make check", which compares
# the results to the reference values. It is the same in the full and
# restricted test modes, since it is a single small computation.
add_test(NAME uniaxial_tensile_test
         COMMAND uniaxial_tensile_test
                 --reference-file IsotropicLinearHardeningPlasticity.ref)
# the same computation on a GPU (see ../config/TestMode.cmake)
if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS)
  add_test(NAME uniaxial_tensile_test-gpu
           COMMAND uniaxial_tensile_test
                   --reference-file IsotropicLinearHardeningPlasticity.ref
                   --execution gpu)
  set_tests_properties(uniaxial_tensile_test-gpu PROPERTIES
                       RESOURCE_LOCK gpu)
endif()
