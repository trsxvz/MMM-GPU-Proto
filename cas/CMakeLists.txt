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

# usage: add_generic_parallel_mfront_libraries(target sources...)
# Behaviours compiled by MFront with the generic-parallel interface on the
# host, for each backend whose configuration file is provided by TFEL
# (openmp, stlpar), in
# <current binary dir>/generic-parallel-<backend>/src/libGenericParallelBehaviour-<backend>.so
# (targets <target>-<backend>). The backends are returned in the variable
# <target>_BACKENDS.
function(add_generic_parallel_mfront_libraries target)
  get_filename_component(_mfront_bin "${MFRONT}" DIRECTORY)
  set(_backends)
  foreach(backend_configuration openmp:config-openmp-host-gcc.json
                                stlpar:config-stlpar-parunseq-gcc.json)
    string(REPLACE ":" ";" _bc "${backend_configuration}")
    list(GET _bc 0 _backend)
    list(GET _bc 1 _configuration)
    if(EXISTS
       "${_mfront_bin}/../share/tfel/mfront/configuration-files/${_configuration}")
      set(_directory "${CMAKE_CURRENT_BINARY_DIR}/generic-parallel-${_backend}")
      set(_library "${_directory}/src/libGenericParallelBehaviour-${_backend}.so")
      file(MAKE_DIRECTORY "${_directory}")
      add_custom_command(
        OUTPUT "${_library}"
        COMMAND ${MFRONT} --obuild --interface=generic-parallel
                "--configuration-file=${_configuration}" ${ARGN}
        WORKING_DIRECTORY "${_directory}"
        DEPENDS ${ARGN}
        COMMENT "MFront: behaviours of ${target}, generic-parallel (${_backend})")
      add_custom_target(${target}-${_backend} ALL DEPENDS "${_library}")
      list(APPEND _backends ${_backend})
    endif()
  endforeach()
  set(${target}_BACKENDS ${_backends} PARENT_SCOPE)
endfunction()

find_package(Python3 COMPONENTS Interpreter)

add_subdirectory(poutre)
add_subdirectory(polycristal)
