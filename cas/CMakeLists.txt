# Cases of the measurement campaign which are not in mfem-mgis-examples or
# mm-opera-hpc. Built against the mfem-mgis installed by build-mfem-mgis.sh
# (MFEMMGIS_DIR), with mfront in the PATH.
cmake_minimum_required(VERSION 3.18)
project(mfem-mgis-proto-cases LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(MFEMMGIS REQUIRED)
find_program(MFRONT mfront REQUIRED)

include(CTest)
enable_testing()
if(NOT TARGET check)
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -T test)
endif()

# sources of mfem-mgis, whose tests provide meshes and behaviours
set(MFEM_MGIS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../mfem-mgis"
  CACHE PATH "sources of mfem-mgis")

# usage: add_mfront_library(target sources...)
# Behaviours compiled by MFront (generic interface) in
# <current binary dir>/src/libBehaviour.so
function(add_mfront_library target)
  set(wd "${CMAKE_CURRENT_BINARY_DIR}")
  add_custom_command(
    OUTPUT "${wd}/src/libBehaviour.so"
    COMMAND ${MFRONT} --obuild --interface=generic ${ARGN}
    WORKING_DIRECTORY "${wd}"
    DEPENDS ${ARGN}
    COMMENT "MFront: behaviours of ${target}")
  add_custom_target(${target} ALL DEPENDS "${wd}/src/libBehaviour.so")
endfunction()

add_subdirectory(poutre)
