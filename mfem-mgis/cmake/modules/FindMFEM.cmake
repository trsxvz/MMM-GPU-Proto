# Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-745557. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the ParElag library. For more information and source code
# availability see http://github.com/LLNL/parelag.
#
# ParElag is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
# Modified version 2023, CEA, France

# Defines the following variables
#   - MFEM_FOUND
#   - MFEM_LIBRARIES
#   - MFEM_INCLUDE_DIRS
#   - MFEM_LIBRARY_DIRS
# and the variables of the `config.mk` file of MFEM (`MFEM_USE_CUDA`, ...).
#
# When MFEM is compiled with CUDA or HIP, the sources containing kernels of
# MFEM must be compiled as CUDA or HIP code: see `MFEMMGISDeviceCode.cmake`,
# loaded by this file (`mfem_mgis_enable_device_language`,
# `mfem_mgis_compile_as_device_code`). CMake 3.18 (CUDA) or 3.21 (HIP) is
# then required.

# FIXME: Perhaps I should try to assert that we're building against
# the same hypre library MFEM built on?

## We require parallel MFEM, which has explicit dependency on METIS and HYPRE

# Start by finding the MFEM config.mk file: in the directory given by
# MFEM_DIR (variable, then environment variable) first, then in the
# prefixes of the environment variable CMAKE_PREFIX_PATH, which must not
# hide the MFEM explicitly given

set(INPUT_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}")
STRING(REPLACE ":" ";" MYSEARCH_PATH "${INPUT_PREFIX_PATH}")
find_file(MFEM_CONFIG_FILE config.mk
  PATHS ${MFEM_DIR}/build $ENV{MFEM_DIR}/build
  ${MFEM_DIR}
  $ENV{MFEM_DIR}
  # MFEM_DIR generally points to MFEM_INSTALL_DIR/lib/share/mfem
  ${MFEM_DIR}/../../../
  $ENV{MFEM_DIR}/../../../
  ${MYSEARCH_PATH}
  PATH_SUFFIXES share/mfem
  NO_DEFAULT_PATH
  DOC "The MFEM configuration file")

if (MFEM_CONFIG_FILE)
  # Extract the directory name
  message(STATUS "Found ${MFEM_CONFIG_FILE}")
  get_filename_component(MFEM_CONFIG_PATH ${MFEM_CONFIG_FILE} DIRECTORY)

  # Extract the relevant list of lines
  file(STRINGS ${MFEM_CONFIG_FILE} MFEM_CONFIG_STRINGS REGEX " = ")
  # Now I need to parse this file; it's not pretty. Mayber there's a better way?
  foreach(MY_LINE IN LISTS MFEM_CONFIG_STRINGS)
    # The name of the variable is before the first equal sign, its value
    # after it. The value may contain equal signs (for instance the flags
    # `-Xcompiler=-fopenmp` of a build with CUDA).
    if(MY_LINE MATCHES "^([^=]*)=(.*)$")
      string(STRIP "${CMAKE_MATCH_1}" MY_VAR_NAME)
      string(STRIP "${CMAKE_MATCH_2}" MY_VAR_VALUE)
      set(${MY_VAR_NAME} "${MY_VAR_VALUE}")
    endif()
  endforeach(MY_LINE)
  set(MFEM_DIR ${MFEM_INSTALL_DIR})
  
else(MFEM_CONFIG_FILE)
  message(FATAL_ERROR "MFEM configuration file not found!")
endif(MFEM_CONFIG_FILE)
  
# Find the header
find_path(MFEM_INCLUDE_DIRS mfem.hpp
  HINTS ${MFEM_DIR}/build $ENV{MFEM_DIR}/build ${MFEM_CONFIG_DIR} ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "MFEM header location.")
find_path(MFEM_INCLUDE_DIRS mfem.hpp)

# When MFEM is compiled with CUDA, `MFEM_CXX` is nvcc and the flags of
# config.mk are given in the syntax of nvcc (`-Xcompiler=...`,
# `-Xlinker=...`). The sources of mfem-mgis compiled by the C++ compiler get
# their equivalent for this compiler. The other flags are unchanged.
function(mfem_mgis_host_flags out)
  set(_flags)
  foreach(_flag IN LISTS ARGN)
    if(_flag MATCHES "^-Xcompiler=(.*)$")
      string(REPLACE "," ";" _options "${CMAKE_MATCH_1}")
      list(APPEND _flags ${_options})
    elseif(_flag MATCHES "^-Xlinker=(.*)$")
      list(APPEND _flags "-Wl,${CMAKE_MATCH_1}")
    else()
      list(APPEND _flags "${_flag}")
    endif()
  endforeach(_flag)
  set(${out} ${_flags} PARENT_SCOPE)
endfunction(mfem_mgis_host_flags)

separate_arguments(MFEM_EXTRA_INC_DIRS UNIX_COMMAND "${MFEM_TPLFLAGS}")
# `-isystem` is followed by its directory: both are given together, since
# `add_compile_options` removes duplicated options. The flags written in
# the syntax of nvcc are given as such to the sources compiled as CUDA code
# by nvcc (see `MFEMMGISDeviceCode.cmake`) and translated for the others.
if(MFEM_USE_CUDA STREQUAL "YES")
  if(CMAKE_VERSION VERSION_LESS 3.18)
    message(FATAL_ERROR "mfem-mgis: CMake 3.18 or later is required with "
            "MFEM compiled with CUDA")
  endif()
endif()
set(_isystem OFF)
set(_nvcc "$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:NVIDIA>>")
foreach(FLAG IN LISTS MFEM_EXTRA_INC_DIRS)
  if(_isystem)
    add_compile_options("SHELL:-isystem \"${FLAG}\"")
    set(_isystem OFF)
  elseif(FLAG STREQUAL "-isystem")
    set(_isystem ON)
  else()
    mfem_mgis_host_flags(_host_flags ${FLAG})
    if("${_host_flags}" STREQUAL "${FLAG}")
      add_compile_options(${FLAG})
    else()
      foreach(_host_flag IN LISTS _host_flags)
        add_compile_options("$<$<NOT:${_nvcc}>:${_host_flag}>")
      endforeach(_host_flag)
      add_compile_options("$<${_nvcc}:${FLAG}>")
    endif()
  endif()
endforeach(FLAG)
if(_isystem)
  message(WARNING "mfem-mgis: the last flag of MFEM_TPLFLAGS is "
          "'-isystem', without directory: ignored")
endif()

# Description of MFEM used to compile its kernels (see
# `MFEMMGISDeviceCode.cmake`)
include("${CMAKE_CURRENT_LIST_DIR}/MFEMMGISDeviceCode.cmake")
mfem_mgis_set_mfem_device_properties("${MFEM_USE_CUDA}" "${MFEM_USE_HIP}"
  "${MFEM_CXX}" "${MFEM_HOST_CXX}" "${MFEM_CXXFLAGS}")


# Find the library
find_library(MFEM_LIBRARY mfem
  HINTS ${MFEM_DIR}/build $ENV{MFEM_DIR}/build ${MFEM_CONFIG_DIR} ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The MFEM library.")

# Check if we have shared or static libraries
include(ParELAGCMakeUtilities)
parelag_determine_library_type(${MFEM_LIBRARY} MFEM_LIB_TYPE)
add_library(MFEM::mfem ${MFEM_LIB_TYPE} IMPORTED)
 
# Set the library location
set_property(TARGET MFEM::mfem
  PROPERTY IMPORTED_LOCATION ${MFEM_LIBRARY})

set_property(TARGET MFEM::mfem APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MFEM_INCLUDE_DIRS})

# Set external libraries (lapack, hypre, metis ...)
separate_arguments(MFEM_EXT_LIBS_LIST UNIX_COMMAND "${MFEM_EXT_LIBS}")
mfem_mgis_host_flags(MFEM_EXT_LIBS_LIST ${MFEM_EXT_LIBS_LIST})
# MFEM compiled with CUDA: its headers include the headers of the CUDA
# toolkit (`cusparse.h`, ...), and its external libraries name the libraries
# of the toolkit (`-lcusparse`, ...) without their directory, which nvcc
# knows but not the C++ compiler
if(MFEM_USE_CUDA STREQUAL "YES")
  # the toolkit of the nvcc used to compile MFEM
  if((NOT DEFINED CUDAToolkit_ROOT) AND (MFEM_CXX MATCHES "nvcc$"))
    get_filename_component(_mfem_cuda_bin "${MFEM_CXX}" DIRECTORY)
    get_filename_component(CUDAToolkit_ROOT "${_mfem_cuda_bin}" DIRECTORY)
  endif()
  # `CUDAToolkit_ROOT` is used (policy CMP0074, not set by the minimum
  # version of CMake required by mfem-mgis)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0074 NEW)
  find_package(CUDAToolkit REQUIRED)
  cmake_policy(POP)
  list(PREPEND MFEM_EXT_LIBS_LIST "-L${CUDAToolkit_LIBRARY_DIR}")
  set_property(TARGET MFEM::mfem APPEND
    PROPERTY INTERFACE_LINK_LIBRARIES CUDA::cudart)
endif()
# MFEM compiled with HIP: its headers include the headers of HIP
# (`hip/hip_runtime.h`, ...), found by the C++ compiler through the target
# `hip::host` of ROCm, which also gives the HIP runtime used by MFEM (linked
# by hipcc, so not named in its external libraries)
if(MFEM_USE_HIP STREQUAL "YES")
  # the ROCm of the hipcc used to compile MFEM
  if((NOT DEFINED hip_ROOT) AND (MFEM_CXX MATCHES "hipcc$"))
    get_filename_component(_mfem_hip_bin "${MFEM_CXX}" DIRECTORY)
    get_filename_component(hip_ROOT "${_mfem_hip_bin}" DIRECTORY)
  endif()
  # MFEM compiled with HIP uses the libraries of ROCm (rocsparse, ...): the
  # platform is AMD, whatever the environment variable HIP_PLATFORM (which
  # may select the NVIDIA platform for other codes)
  if(NOT DEFINED HIP_PLATFORM)
    set(HIP_PLATFORM amd)
  endif()
  # `hip_ROOT` is used (policy CMP0074, not set by the minimum version of
  # CMake required by mfem-mgis)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0074 NEW)
  find_package(hip CONFIG REQUIRED)
  cmake_policy(POP)
  set_property(TARGET MFEM::mfem APPEND
    PROPERTY INTERFACE_LINK_LIBRARIES hip::host)
endif()
set_property(TARGET MFEM::mfem APPEND
    PROPERTY INTERFACE_LINK_LIBRARIES ${MFEM_EXT_LIBS_LIST})

# Set MPI library
if(MFEM_USE_MPI)
  # The following include directories are automatically filled within FindMFEM.cmake
  find_package(MPI REQUIRED)
  set_property(TARGET MFEM::mfem APPEND
    PROPERTY INTERFACE_LINK_LIBRARIES
    ${MPI_CXX_LIBRARIES})
    include_directories(${MPI_INCLUDE_PATH})
endif(MFEM_USE_MPI)

# Optional MUMPS
if(MFEM_USE_MUMPS)
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(MFEM_USE_MUMPS)

if (NOT MFEM_USE_SUITESPARSE)
  message(FATAL_ERROR "You shall have SUITESPARSE support activated in MFEM.CMake will exit.\n" )
endif()

# Set the include directories
set(MFEM_INCLUDE_DIRS ${MFEM_INCLUDE_DIRS}
  CACHE PATH
  "Directories in which to find headers for MFEM.")
mark_as_advanced(FORCE MFEM_INCLUDE_DIRS)

# Set the libraries
set(MFEM_LIBRARIES MFEM::mfem)
set(MFEM_LIBRARY ${MFEM_LIBRARY} CACHE PATH
  "The primary library file for MFEM")

mark_as_advanced(FORCE MFEM_LIBRARY)
mark_as_advanced(FORCE MFEM_CONFIG_FILE)


# This handles "REQUIRED" etc keywords
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MFEM
  DEFAULT_MSG
  MFEM_CONFIG_FILE MFEM_INCLUDE_DIRS)
