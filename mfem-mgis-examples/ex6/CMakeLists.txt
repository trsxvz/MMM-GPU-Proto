cmake_minimum_required(VERSION 3.13)
project(RveNonLinearElasticExample)

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

# Copy the files needed at run time. Unlike file(COPY), configure_file updates
# the copies when these files are modified.
foreach(input SaintVenantKirchhoffElasticity.mfront cube_2mat_per.mesh
              inclusions_49.geo)
  configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                 COPYONLY)
endforeach()

# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU,
# chosen at run time by mfem_mgis::getBehavioursLibrary from the
# option --execution
mfem_mgis_add_behaviours_library(BehavioursEx6
  SOURCES ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront)

add_executable(rve RveNonLinearElastic.cxx)
target_compile_definitions(rve PRIVATE
  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursEx6_BEHAVIOURS_DIRECTORY}")
target_link_libraries(rve mfem-mgis::MFEMMGIS)
target_include_directories(rve
  SYSTEM
  PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
  PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")

add_dependencies(check rve BehavioursEx6)

# The test is the same in the full and restricted test modes, since it is a
# single computation on the two-layer cube, compared to its analytical
# solution.
set(test "rve")
add_test(NAME ${test} COMMAND rve --no-post-processing)
# the same computation on a GPU, and with the behaviours integrated on the
# CPU (see ../config/TestMode.cmake)
if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS AND BehavioursEx6_HAVE_GPU_BEHAVIOURS)
  add_test(NAME ${test}-gpu
           COMMAND rve --no-post-processing --execution gpu)
  add_test(NAME ${test}-mixed-gpu
           COMMAND rve --no-post-processing --execution mixed
                   --behaviours-execution cpu)
  set_tests_properties(${test}-gpu ${test}-mixed-gpu PROPERTIES
                       RESOURCE_LOCK gpu)
endif()
