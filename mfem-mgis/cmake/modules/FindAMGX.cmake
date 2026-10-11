# Find the AmgX library of NVIDIA (https://github.com/NVIDIA/AMGX)
#
# The installation directory of AmgX may be given by the variable
# `AMGX_DIR` or by the environment variable of the same name. This module
# defines:
#   - AMGX_FOUND
#   - AMGX_INCLUDE_DIR, directory of `amgx_c.h`
#   - AMGX_LIBRARY, shared library of AmgX (`libamgxsh`)
#   - the imported target `AMGX::amgxsh`, linked with the runtime of CUDA
#     (`CUDA::cudart`) if this target is defined

# the installation directory given is searched first, before the prefixes
# of CMake (`CMAKE_PREFIX_PATH`), which may contain another AmgX
foreach(_amgx_paths "NO_DEFAULT_PATH" "")
  find_path(AMGX_INCLUDE_DIR amgx_c.h
    HINTS ${AMGX_DIR} $ENV{AMGX_DIR}
    PATH_SUFFIXES include
    DOC "Directory of the headers of AmgX"
    ${_amgx_paths})
  find_library(AMGX_LIBRARY amgxsh
    HINTS ${AMGX_DIR} $ENV{AMGX_DIR}
    PATH_SUFFIXES lib lib64
    DOC "Shared library of AmgX"
    ${_amgx_paths})
endforeach(_amgx_paths)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMGX
  REQUIRED_VARS AMGX_LIBRARY AMGX_INCLUDE_DIR)
if(AMGX_FOUND AND (NOT TARGET AMGX::amgxsh))
  add_library(AMGX::amgxsh SHARED IMPORTED)
  set_target_properties(AMGX::amgxsh PROPERTIES
    IMPORTED_LOCATION "${AMGX_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${AMGX_INCLUDE_DIR}")
  if(TARGET CUDA::cudart)
    set_property(TARGET AMGX::amgxsh APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES CUDA::cudart)
  endif()
endif()
mark_as_advanced(AMGX_INCLUDE_DIR AMGX_LIBRARY)
