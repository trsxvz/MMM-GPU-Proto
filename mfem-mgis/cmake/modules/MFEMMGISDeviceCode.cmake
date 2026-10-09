# Compilation of the sources containing kernels of MFEM (`mfem::forall`).
#
# When MFEM is compiled with CUDA, its kernels must be compiled as CUDA code
# (MFEM checks it at compile time). This module provides:
#
# - `mfem_mgis_enable_device_language()`: enables the CUDA language if MFEM
#   is compiled with CUDA, and does nothing otherwise. The CUDA compiler and
#   its host compiler default to the ones used to compile MFEM
#   (`MFEM_CXX` and `MFEM_HOST_CXX` in `config.mk`), unless
#   `CMAKE_CUDA_COMPILER` and `CMAKE_CUDA_HOST_COMPILER` are given. This
#   macro enables a language: it must be called at file scope (not in a
#   function), before `mfem_mgis_compile_as_device_code`.
# - `mfem_mgis_compile_as_device_code(target sources...)`: compiles the
#   given sources of the target as CUDA code if MFEM is compiled with CUDA,
#   and does nothing otherwise.
#
# Both rely on global properties describing MFEM, set by `FindMFEM.cmake`
# or by the configuration file of mfem-mgis:
#
# - `MFEM_MGIS_MFEM_USE_CUDA`: `YES` if MFEM is compiled with CUDA;
# - `MFEM_MGIS_MFEM_USE_HIP`: `YES` if MFEM is compiled with HIP (not
#   supported yet);
# - `MFEM_MGIS_MFEM_CUDA_COMPILER`, `MFEM_MGIS_MFEM_CUDA_HOST_COMPILER`:
#   compilers used to compile MFEM;
# - `MFEM_MGIS_MFEM_CUDA_ARCHITECTURES`: GPU architectures of MFEM, in the
#   syntax of the `CUDA_ARCHITECTURES` property;
# - `MFEM_MGIS_MFEM_CXX_STANDARD`: C++ standard used to compile MFEM.

# usage: mfem_mgis_set_mfem_device_properties(use_cuda use_hip cxx host_cxx
#                                             cxxflags)
# Sets the global properties describing MFEM from the values of its
# `config.mk` file (`MFEM_USE_CUDA`, `MFEM_USE_HIP`, `MFEM_CXX`,
# `MFEM_HOST_CXX`, `MFEM_CXXFLAGS`).
function(mfem_mgis_set_mfem_device_properties use_cuda use_hip cxx host_cxx
         cxxflags)
  set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_CUDA "${use_cuda}")
  set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_HIP "${use_hip}")
  foreach(_property CUDA_COMPILER CUDA_HOST_COMPILER CUDA_ARCHITECTURES
                    CXX_STANDARD)
    set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_${_property} "")
  endforeach(_property)
  if(NOT use_cuda STREQUAL "YES")
    return()
  endif()
  set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_CUDA_COMPILER "${cxx}")
  set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_CUDA_HOST_COMPILER
               "${host_cxx}")
  separate_arguments(_flags UNIX_COMMAND "${cxxflags}")
  # GPU architectures: `-arch=X`, `--gpu-architecture=X` or
  # `-gencode arch=compute_X,code=[sm_X,compute_X]` (or `code=sm_X`), the
  # value of each option being possibly given in the next argument
  set(_archs)
  set(_codes)
  set(_option "")
  foreach(_flag IN LISTS _flags)
    if(_option)
      set(_flag "${_option}=${_flag}")
      set(_option "")
    endif()
    if(_flag MATCHES "^(-arch|--gpu-architecture|-gencode|--generate-code)$")
      set(_option "${_flag}")
    elseif(_flag MATCHES "^(-arch|--gpu-architecture)=?(.+)$")
      list(APPEND _archs "${CMAKE_MATCH_2}")
    elseif(_flag MATCHES "^(-gencode|--generate-code)=?(.*)code=\\[?([^]]*)\\]?$")
      string(REPLACE "," ";" _c "${CMAKE_MATCH_3}")
      list(APPEND _codes ${_c})
    elseif(_flag MATCHES "^-std=(c|gnu)\\+\\+([0-9]+)$")
      set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_CXX_STANDARD
                   "${CMAKE_MATCH_2}")
    endif()
  endforeach(_flag)
  # translation in the syntax of the `CUDA_ARCHITECTURES` property: `sm_X`
  # gives real code (`X-real`) and `compute_X` virtual code (`X-virtual`),
  # `-arch=sm_X` both (`X`)
  set(_cuda_archs)
  foreach(_arch IN LISTS _archs)
    if(_arch MATCHES "^(native|all|all-major)$")
      list(APPEND _cuda_archs "${_arch}")
    elseif(_arch MATCHES "^(sm|compute)_([0-9]+[af]?)$")
      list(APPEND _cuda_archs "${CMAKE_MATCH_2}")
    endif()
  endforeach(_arch)
  foreach(_code IN LISTS _codes)
    if(_code MATCHES "^sm_([0-9]+[af]?)$")
      list(APPEND _cuda_archs "${CMAKE_MATCH_1}-real")
    elseif(_code MATCHES "^compute_([0-9]+[af]?)$")
      list(APPEND _cuda_archs "${CMAKE_MATCH_1}-virtual")
    endif()
  endforeach(_code)
  list(REMOVE_DUPLICATES _cuda_archs)
  if(NOT _cuda_archs)
    message(WARNING "mfem-mgis: the GPU architectures of MFEM are not "
            "found in its compilation flags (${cxxflags}): the sources "
            "containing kernels of MFEM are compiled for the default "
            "architectures of the CUDA compiler")
  endif()
  set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_CUDA_ARCHITECTURES
               "${_cuda_archs}")
endfunction(mfem_mgis_set_mfem_device_properties)

# usage: mfem_mgis_enable_device_language()
macro(mfem_mgis_enable_device_language)
  get_property(_mfem_mgis_use_cuda GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_CUDA)
  if(_mfem_mgis_use_cuda STREQUAL "YES")
    if(CMAKE_VERSION VERSION_LESS 3.18)
      message(FATAL_ERROR "mfem-mgis: CMake 3.18 or later is required to "
              "compile the kernels of MFEM compiled with CUDA")
    endif()
    if(NOT CMAKE_CUDA_COMPILER AND NOT DEFINED ENV{CUDACXX})
      get_property(CMAKE_CUDA_COMPILER GLOBAL
                   PROPERTY MFEM_MGIS_MFEM_CUDA_COMPILER)
    endif()
    if(NOT CMAKE_CUDA_HOST_COMPILER AND NOT DEFINED ENV{CUDAHOSTCXX})
      get_property(CMAKE_CUDA_HOST_COMPILER GLOBAL
                   PROPERTY MFEM_MGIS_MFEM_CUDA_HOST_COMPILER)
    endif()
    enable_language(CUDA)
  endif()
endmacro(mfem_mgis_enable_device_language)

# usage: mfem_mgis_compile_as_device_code(target sources...)
function(mfem_mgis_compile_as_device_code target)
  get_property(_use_hip GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_HIP)
  if(_use_hip STREQUAL "YES")
    message(FATAL_ERROR "mfem-mgis: the compilation of the kernels of MFEM "
            "compiled with HIP is not supported yet")
  endif()
  get_property(_use_cuda GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_CUDA)
  if(NOT _use_cuda STREQUAL "YES")
    return()
  endif()
  get_property(_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  if(NOT "CUDA" IN_LIST _languages)
    message(FATAL_ERROR "mfem-mgis: MFEM is compiled with CUDA and the "
            "sources of target '${target}' contain kernels of MFEM: the "
            "CUDA language must be enabled (see "
            "`mfem_mgis_enable_device_language`)")
  endif()
  set_source_files_properties(${ARGN}
    TARGET_DIRECTORY ${target} PROPERTIES LANGUAGE CUDA)
  # same C++ standard as the other sources of the target, or as MFEM
  get_target_property(_standard ${target} CXX_STANDARD)
  if(NOT _standard)
    set(_standard "${CMAKE_CXX_STANDARD}")
  endif()
  if(NOT _standard)
    get_property(_standard GLOBAL PROPERTY MFEM_MGIS_MFEM_CXX_STANDARD)
  endif()
  if(_standard)
    set_target_properties(${target} PROPERTIES
      CUDA_STANDARD "${_standard}" CUDA_STANDARD_REQUIRED ON)
  endif()
  get_property(_archs GLOBAL PROPERTY MFEM_MGIS_MFEM_CUDA_ARCHITECTURES)
  if(_archs)
    set_target_properties(${target} PROPERTIES CUDA_ARCHITECTURES "${_archs}")
  endif()
  # MFEM is linked to the shared CUDA runtime (`CUDA::cudart`): a second,
  # static, runtime would not share its state (device, streams, memory)
  set_target_properties(${target} PROPERTIES
    CUDA_RUNTIME_LIBRARY Shared LINKER_LANGUAGE CXX)
  # options required by the kernels of MFEM (nvcc only)
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    target_compile_options(${target} PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
  endif()
  # the definitions and visibility options of `CMAKE_CXX_FLAGS` (and of its
  # variant for the build type) do not apply to the CUDA sources: they are
  # given to them, since the C++ and CUDA sources of the target share
  # classes (for instance `_GLIBCXX_DEBUG` changes the standard containers)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type)
  separate_arguments(_cxx_flags UNIX_COMMAND
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}}")
  foreach(_flag IN LISTS _cxx_flags)
    if(_flag MATCHES "^-[DU].")
      target_compile_options(${target} PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:${_flag}>")
    elseif((_flag MATCHES "^-fvisibility") AND
           (CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA"))
      target_compile_options(${target} PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_flag}>")
    endif()
  endforeach(_flag)
endfunction(mfem_mgis_compile_as_device_code)
