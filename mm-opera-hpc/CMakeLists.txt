cmake_minimum_required(VERSION 3.22)
project(opera-hpc
        VERSION 0.1.0
        DESCRIPTION "Evolution of a population of bubble in a brittle material")

# ===========================
# === CMake customization ===
# ===========================
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Coverage RelWithDebInfo Debug)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
include(cmake/modules/tfel.cmake)

find_package(MFEMMGIS REQUIRED)

include(CTest)
# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
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

# Tests on a GPU, selected by the MM_OPERA_HPC_GPU_TESTS option (off by
# default, a GPU being required when the tests are run): the cases are also
# tested with the option --execution gpu, if MFEM is compiled with CUDA or
# HIP (as described by mfem-mgis). The MM_OPERA_HPC_RUN_GPU_TESTS variable is
# then true.
option(MM_OPERA_HPC_GPU_TESTS "Test the cases on a GPU too" OFF)
set(MM_OPERA_HPC_RUN_GPU_TESTS OFF)
if(MM_OPERA_HPC_GPU_TESTS)
  get_property(_mm_opera_hpc_mfem_use_cuda GLOBAL
               PROPERTY MFEM_MGIS_MFEM_USE_CUDA)
  get_property(_mm_opera_hpc_mfem_use_hip GLOBAL
               PROPERTY MFEM_MGIS_MFEM_USE_HIP)
  if((_mm_opera_hpc_mfem_use_cuda STREQUAL "YES") OR
     (_mm_opera_hpc_mfem_use_hip STREQUAL "YES"))
    set(MM_OPERA_HPC_RUN_GPU_TESTS ON)
  endif()
endif()

add_subdirectory(src)
add_subdirectory(polycrystal)
add_subdirectory(bubble)
add_subdirectory(cermet)

