# Find roctx, the annotation library of ROCm: the ranges declared by a
# program are shown by the profilers of ROCm (rocprofv3...) with the calls
# of the HIP runtime and the kernels
#
# The library of rocprofiler-sdk (`librocprofiler-sdk-roctx`, header
# `rocprofiler-sdk-roctx/roctx.h`, ROCm 6.2 or later) is searched first,
# then the one of roctracer which it replaces (`libroctx64`, header
# `roctracer/roctx.h`), both having the same functions. They are searched in
# the installation of ROCm given by the variable `ROCTX_DIR` or by the
# environment variable of the same name, then in the one of HIP (`hip_DIR`)
# and in the one given by the environment variable `ROCM_PATH`. This module
# defines:
#   - ROCTX_FOUND
#   - ROCTX_INCLUDE_DIR, directory of `rocprofiler-sdk-roctx/roctx.h` or of
#     `roctracer/roctx.h`
#   - ROCTX_LIBRARY, shared library of roctx

set(_roctx_hints ${ROCTX_DIR} $ENV{ROCTX_DIR})
if(hip_DIR)
  # <ROCm>/lib/cmake/hip
  get_filename_component(_roctx_rocm "${hip_DIR}/../../.." ABSOLUTE)
  list(APPEND _roctx_hints "${_roctx_rocm}")
endif()
list(APPEND _roctx_hints $ENV{ROCM_PATH})
# these installations are searched first, before the prefixes of CMake
# (`CMAKE_PREFIX_PATH`), which may contain another ROCm
foreach(_roctx_paths "NO_DEFAULT_PATH" "")
  find_library(ROCTX_LIBRARY NAMES rocprofiler-sdk-roctx roctx64
    HINTS ${_roctx_hints}
    PATH_SUFFIXES lib lib64
    DOC "Shared library of roctx"
    ${_roctx_paths})
  find_path(ROCTX_INCLUDE_DIR
    NAMES rocprofiler-sdk-roctx/roctx.h roctracer/roctx.h
    HINTS ${_roctx_hints}
    PATH_SUFFIXES include
    DOC "Directory of the headers of roctx"
    ${_roctx_paths})
endforeach(_roctx_paths)
unset(_roctx_hints)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ROCTX
  REQUIRED_VARS ROCTX_LIBRARY ROCTX_INCLUDE_DIR)
mark_as_advanced(ROCTX_INCLUDE_DIR ROCTX_LIBRARY)
