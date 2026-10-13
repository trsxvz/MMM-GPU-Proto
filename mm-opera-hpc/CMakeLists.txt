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

# Tests comparing results to reference solutions: the case is run by
# cmake/run-and-check.py, which then compares values of its output file to
# the expected ones. Without Python 3, the cases are only run.
#
# usage: mm_opera_hpc_add_checked_test(NAME <name> FILE <output file>
#            TOLERANCE <relative tolerance>
#            CHECKS <row> <column> <expected> [<row> <column> <expected>...]
#            COMMAND <command...>)
set(MM_OPERA_HPC_RUN_AND_CHECK
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run-and-check.py")
find_package(Python3 COMPONENTS Interpreter QUIET)
if(NOT Python3_Interpreter_FOUND)
  message(STATUS "Python 3 not found: the results of the cases are not "
                 "compared to the reference solutions")
endif()
function(mm_opera_hpc_add_checked_test)
  cmake_parse_arguments(_t "" "NAME;FILE;TOLERANCE" "CHECKS;COMMAND" ${ARGN})
  if(NOT Python3_Interpreter_FOUND)
    add_test(NAME ${_t_NAME} COMMAND ${_t_COMMAND})
    return()
  endif()
  set(_checks)
  list(LENGTH _t_CHECKS _n)
  math(EXPR _last "${_n} - 1")
  foreach(_i RANGE 0 ${_last} 3)
    math(EXPR _j "${_i} + 1")
    math(EXPR _k "${_i} + 2")
    list(GET _t_CHECKS ${_i} _row)
    list(GET _t_CHECKS ${_j} _column)
    list(GET _t_CHECKS ${_k} _expected)
    list(APPEND _checks --check ${_row} ${_column} ${_expected})
  endforeach()
  add_test(NAME ${_t_NAME}
           COMMAND ${Python3_EXECUTABLE}
                   ${MM_OPERA_HPC_RUN_AND_CHECK}
                   --file ${_t_FILE} --tolerance ${_t_TOLERANCE} ${_checks}
                   -- ${_t_COMMAND})
endfunction()

add_subdirectory(src)
add_subdirectory(polycrystal)
add_subdirectory(bubble)
add_subdirectory(cermet)

