# Compilation of the sources containing kernels of MFEM (`mfem::forall`).
#
# When MFEM is compiled with CUDA or HIP, its kernels must be compiled as
# CUDA or HIP code (MFEM checks it at compile time). This module provides:
#
# - `mfem_mgis_enable_device_language()`: enables the CUDA or HIP language
#   if MFEM is compiled with CUDA or HIP, and does nothing otherwise. With
#   CUDA, the CUDA compiler and its host compiler default to the ones used
#   to compile MFEM (`MFEM_CXX` and `MFEM_HOST_CXX` in `config.mk`), unless
#   `CMAKE_CUDA_COMPILER` and `CMAKE_CUDA_HOST_COMPILER` are given. With
#   HIP, the HIP compiler defaults to the clang of ROCm used by the hipcc
#   which compiled MFEM (CMake does not accept hipcc as HIP compiler),
#   unless `CMAKE_HIP_COMPILER` is given, the platform is AMD (the one of
#   the libraries of ROCm used by MFEM) unless `CMAKE_HIP_PLATFORM` is given,
#   and `CMAKE_HIP_ARCHITECTURES`, unless given, is set to the GPU
#   architectures of MFEM (CMake would otherwise detect the GPU of the
#   machine, or use the default architecture of clang), the sources being
#   compiled for the architectures of MFEM anyway (as with CUDA). This
#   macro enables a language: it must be called at file scope (not in a
#   function), before `mfem_mgis_compile_as_device_code`. A warning is
#   printed if the compiler is not the one of MFEM (for instance a build
#   directory configured before a change of MFEM).
# - `mfem_mgis_compile_as_device_code(target sources...)`: compiles the
#   given sources of the target as CUDA or HIP code if MFEM is compiled with
#   CUDA or HIP, and does nothing otherwise.
# - `mfem_mgis_disable_device_floating_point_contraction(target
#   sources...)`: compiles the device code of the given sources of the
#   target without contraction of the floating point operations into fused
#   multiply-add operations.
#
# The first two rely on global properties describing MFEM, set by
# `FindMFEM.cmake` or by the configuration file of mfem-mgis:
#
# - `MFEM_MGIS_MFEM_USE_CUDA`: `YES` if MFEM is compiled with CUDA;
# - `MFEM_MGIS_MFEM_USE_HIP`: `YES` if MFEM is compiled with HIP;
# - `MFEM_MGIS_MFEM_CUDA_COMPILER`, `MFEM_MGIS_MFEM_CUDA_HOST_COMPILER`:
#   compilers used to compile MFEM with CUDA;
# - `MFEM_MGIS_MFEM_CUDA_ARCHITECTURES`: GPU architectures of MFEM, in the
#   syntax of the `CUDA_ARCHITECTURES` property;
# - `MFEM_MGIS_MFEM_HIP_COMPILER`: compiler used to compile MFEM with HIP
#   (generally hipcc);
# - `MFEM_MGIS_MFEM_HIP_ARCHITECTURES`: GPU architectures of MFEM, in the
#   syntax of the `HIP_ARCHITECTURES` property (`gfx90a`, ...).

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
                    HIP_COMPILER HIP_ARCHITECTURES)
    set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_${_property} "")
  endforeach(_property)
  if(use_hip STREQUAL "YES")
    set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_HIP_COMPILER "${cxx}")
    # GPU architectures: `--offload-arch=X` (X possibly being a list
    # separated by commas) or the former `--amdgpu-target=X`, the value of
    # each option being possibly given in the next argument
    separate_arguments(_flags UNIX_COMMAND "${cxxflags}")
    set(_hip_archs)
    set(_option "")
    foreach(_flag IN LISTS _flags)
      if(_option)
        set(_flag "${_option}=${_flag}")
        set(_option "")
      endif()
      if(_flag MATCHES "^(--offload-arch|--amdgpu-target)$")
        set(_option "${_flag}")
      elseif(_flag MATCHES "^(--offload-arch|--amdgpu-target)=(.+)$")
        string(REPLACE "," ";" _a "${CMAKE_MATCH_2}")
        list(APPEND _hip_archs ${_a})
      endif()
    endforeach(_flag)
    list(REMOVE_DUPLICATES _hip_archs)
    if(NOT _hip_archs)
      message(WARNING "mfem-mgis: the GPU architectures of MFEM are not "
              "found in its compilation flags (${cxxflags}): the sources "
              "containing kernels of MFEM are compiled for the "
              "architectures given by CMAKE_HIP_ARCHITECTURES")
    endif()
    set_property(GLOBAL PROPERTY MFEM_MGIS_MFEM_HIP_ARCHITECTURES
                 "${_hip_archs}")
  endif()
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

# usage: mfem_mgis_hip_compiler(out hipcc)
# Sets `out` to the clang of ROCm used by the given HIP compiler of MFEM:
# the compiler itself if it is a clang, the clang reported by the
# `hipconfig` installed next to it if it is hipcc, or an empty value if it
# is not found (CMake then searches the HIP compiler itself).
function(mfem_mgis_hip_compiler out hipcc)
  set(${out} "" PARENT_SCOPE)
  get_filename_component(_name "${hipcc}" NAME)
  if(_name MATCHES "clang")
    set(${out} "${hipcc}" PARENT_SCOPE)
  elseif(_name MATCHES "^hipcc")
    get_filename_component(_bin "${hipcc}" DIRECTORY)
    if(EXISTS "${_bin}/hipconfig")
      # AMD platform (the answer of hipconfig depends on HIP_PLATFORM)
      execute_process(COMMAND ${CMAKE_COMMAND} -E env HIP_PLATFORM=amd
                              "${_bin}/hipconfig" --hipclangpath
                      OUTPUT_VARIABLE _path RESULT_VARIABLE _status
                      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if((_status EQUAL 0) AND (EXISTS "${_path}/clang++"))
        set(${out} "${_path}/clang++" PARENT_SCOPE)
      endif()
    endif()
  endif()
endfunction(mfem_mgis_hip_compiler)

# usage: mfem_mgis_enable_device_language()
macro(mfem_mgis_enable_device_language)
  get_property(_mfem_mgis_use_hip GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_HIP)
  if(_mfem_mgis_use_hip STREQUAL "YES")
    if(CMAKE_VERSION VERSION_LESS 3.21)
      message(FATAL_ERROR "mfem-mgis: CMake 3.21 or later is required to "
              "compile the kernels of MFEM compiled with HIP")
    endif()
    get_property(_mfem_mgis_hipcc GLOBAL PROPERTY MFEM_MGIS_MFEM_HIP_COMPILER)
    mfem_mgis_hip_compiler(_mfem_mgis_hip_clang "${_mfem_mgis_hipcc}")
    if(NOT CMAKE_HIP_COMPILER AND ("$ENV{HIPCXX}" STREQUAL "")
       AND _mfem_mgis_hip_clang)
      set(CMAKE_HIP_COMPILER "${_mfem_mgis_hip_clang}")
    endif()
    # MFEM compiled with HIP uses the libraries of ROCm: AMD platform
    if(NOT DEFINED CMAKE_HIP_PLATFORM)
      set(CMAKE_HIP_PLATFORM amd)
    endif()
    if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
      get_property(_mfem_mgis_hip_archs GLOBAL
                   PROPERTY MFEM_MGIS_MFEM_HIP_ARCHITECTURES)
      if(_mfem_mgis_hip_archs)
        set(CMAKE_HIP_ARCHITECTURES "${_mfem_mgis_hip_archs}")
      endif()
    endif()
    enable_language(HIP)
    # the programs and libraries linked to mfem-mgis (a static library
    # containing HIP code) are linked by the C++ compiler, as with CUDA, and
    # not by the HIP compiler (whose linker preference is greater than the
    # one of C++); the HIP runtime is given by the target `hip::host`
    set(CMAKE_HIP_LINKER_PREFERENCE_PROPAGATES 0)
    if(_mfem_mgis_hip_clang)
      get_filename_component(_mfem_mgis_c1 "${CMAKE_HIP_COMPILER}" REALPATH)
      get_filename_component(_mfem_mgis_c2 "${_mfem_mgis_hip_clang}" REALPATH)
      if(NOT _mfem_mgis_c1 STREQUAL _mfem_mgis_c2)
        message(WARNING "mfem-mgis: the HIP compiler "
                "(${CMAKE_HIP_COMPILER}) is not the one of MFEM "
                "(${_mfem_mgis_hip_clang}, used by ${_mfem_mgis_hipcc})")
      endif()
    endif()
  endif()
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
    get_property(_mfem_mgis_cuda_compiler GLOBAL
                 PROPERTY MFEM_MGIS_MFEM_CUDA_COMPILER)
    if(_mfem_mgis_cuda_compiler)
      get_filename_component(_mfem_mgis_c1 "${CMAKE_CUDA_COMPILER}" REALPATH)
      get_filename_component(_mfem_mgis_c2 "${_mfem_mgis_cuda_compiler}"
                             REALPATH)
      if(NOT _mfem_mgis_c1 STREQUAL _mfem_mgis_c2)
        message(WARNING "mfem-mgis: the CUDA compiler "
                "(${CMAKE_CUDA_COMPILER}) is not the one of MFEM "
                "(${_mfem_mgis_cuda_compiler})")
      endif()
    endif()
  endif()
endmacro(mfem_mgis_enable_device_language)

# usage: mfem_mgis_compile_as_device_code(target sources...)
function(mfem_mgis_compile_as_device_code target)
  get_property(_use_cuda GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_CUDA)
  get_property(_use_hip GLOBAL PROPERTY MFEM_MGIS_MFEM_USE_HIP)
  if(_use_cuda STREQUAL "YES")
    set(_language CUDA)
  elseif(_use_hip STREQUAL "YES")
    set(_language HIP)
  else()
    return()
  endif()
  # the language must be enabled in the directory of the call (or in one of
  # its parents)
  if(NOT CMAKE_${_language}_COMPILER_LOADED)
    message(FATAL_ERROR "mfem-mgis: MFEM is compiled with ${_language} and "
            "the sources of target '${target}' contain kernels of MFEM: the "
            "${_language} language must be enabled (see "
            "`mfem_mgis_enable_device_language`)")
  endif()
  set_source_files_properties(${ARGN}
    TARGET_DIRECTORY ${target} PROPERTIES LANGUAGE ${_language})
  # same C++ standard as the other sources of the target, at least C++20
  # (required by the headers of mfem-mgis)
  get_target_property(_standard ${target} CXX_STANDARD)
  if(NOT _standard)
    set(_standard "${CMAKE_CXX_STANDARD}")
  endif()
  if((NOT _standard) OR (_standard LESS 20) OR (_standard EQUAL 98))
    set(_standard 20)
  endif()
  set_target_properties(${target} PROPERTIES
    ${_language}_STANDARD "${_standard}" ${_language}_STANDARD_REQUIRED ON)
  get_property(_archs GLOBAL PROPERTY MFEM_MGIS_MFEM_${_language}_ARCHITECTURES)
  if(_archs)
    set_target_properties(${target} PROPERTIES
      ${_language}_ARCHITECTURES "${_archs}")
  endif()
  # the target is linked by the C++ compiler. The CUDA sources use the
  # shared CUDA runtime (`CUDA::cudart`), linked to all the code of
  # mfem-mgis by the target `MFEM::mfem`, rather than a static copy per
  # target (MFEM itself embeds its own static runtime). The HIP sources use
  # the HIP runtime of `hip::host`, linked in the same way
  set_target_properties(${target} PROPERTIES LINKER_LANGUAGE CXX)
  if(_language STREQUAL "CUDA")
    set_target_properties(${target} PROPERTIES CUDA_RUNTIME_LIBRARY Shared)
  endif()
  # options required by the kernels of MFEM (nvcc only)
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    target_compile_options(${target} PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
  endif()
  # the definitions and visibility options of `CMAKE_CXX_FLAGS` (and of its
  # variant for the build type) do not apply to the CUDA or HIP sources:
  # they are given to them, since the C++ and device sources of the target
  # share classes (for instance `_GLIBCXX_DEBUG` changes the standard
  # containers)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type)
  separate_arguments(_cxx_flags UNIX_COMMAND
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}}")
  foreach(_flag IN LISTS _cxx_flags)
    if(_flag MATCHES "^-[DU].")
      target_compile_options(${target} PRIVATE
        "$<$<COMPILE_LANGUAGE:${_language}>:${_flag}>")
    elseif(_flag MATCHES "^-fvisibility")
      if(CMAKE_${_language}_COMPILER_ID STREQUAL "NVIDIA")
        target_compile_options(${target} PRIVATE
          "$<$<COMPILE_LANGUAGE:${_language}>:-Xcompiler=${_flag}>")
      elseif(_language STREQUAL "HIP")
        target_compile_options(${target} PRIVATE
          "$<$<COMPILE_LANGUAGE:HIP>:${_flag}>")
      endif()
    endif()
  endforeach(_flag)
endfunction(mfem_mgis_compile_as_device_code)

# usage: mfem_mgis_disable_device_floating_point_contraction(target
#                                                           sources...)
# Compiles the device code of the given sources of the target without
# contraction of the floating point operations into fused multiply-add
# operations: `-fmad=false` for nvcc, `-ffp-contract=off` for clang (HIP, or
# CUDA compiled by clang, the option then also applying to the host code of
# the sources). Does nothing if the sources are not compiled as CUDA or HIP
# code (see `mfem_mgis_compile_as_device_code`).
function(mfem_mgis_disable_device_floating_point_contraction target)
  set(_options)
  if(CMAKE_CUDA_COMPILER_LOADED)
    if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
      list(APPEND _options "$<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>")
    else()
      list(APPEND _options "$<$<COMPILE_LANGUAGE:CUDA>:-ffp-contract=off>")
    endif()
  endif()
  if(CMAKE_HIP_COMPILER_LOADED)
    list(APPEND _options "$<$<COMPILE_LANGUAGE:HIP>:-ffp-contract=off>")
  endif()
  if(NOT _options)
    return()
  endif()
  # the CUDA and HIP languages require CMake 3.18 or later (see
  # `mfem_mgis_enable_device_language`)
  set_property(SOURCE ${ARGN} TARGET_DIRECTORY ${target}
    APPEND PROPERTY COMPILE_OPTIONS ${_options})
endfunction(mfem_mgis_disable_device_floating_point_contraction)
