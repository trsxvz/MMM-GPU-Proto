# Optional libraries of MFEM used by mfem-mgis (solvers and preconditioners).
#
# For each library, the cache variable `enable-<name>` selects its use:
#
# - `AUTO` (default): the library is used if MFEM is built with it;
# - `ON`: the library is required, an error is reported if MFEM is not built
#   with it;
# - `OFF`: the library is not used, even if MFEM is built with it.
#
# The other boolean constants of CMake (`TRUE`, `YES`, `1`, `FALSE`, `NO`,
# `0`...) are accepted for `ON` and `OFF`.
#
# The function `mfem_mgis_select_optional_library` gives the same choice for
# the other optional libraries (packages of Trilinos, see `CMakeLists.txt`).
#
# The variable `MFEM_MGIS_HAVE_<NAME>` (`ON` or `OFF`) gives the result, and
# the target `MFEMMGIS` defines the macro of the same name (see
# `src/CMakeLists.txt`), which the sources use instead of the macros
# `MFEM_USE_<NAME>` of MFEM. A package manager maps each library to one
# switch (for instance the variants `+petsc` and `~petsc` of spack to
# `-Denable-petsc=ON` and `-Denable-petsc=OFF`).

# names of the options, variables of MFEM and descriptions
set(MFEM_MGIS_OPTIONAL_MFEM_LIBRARIES
  "petsc:PETSC:Krylov solvers, preconditioners and nonlinear solvers of PETSc"
  "ginkgo:GINKGO:solvers and preconditioners of Ginkgo"
  "mumps:MUMPS:direct solver MUMPS"
  "superlu:SUPERLU:direct solver SuperLU_DIST"
  "strumpack:STRUMPACK:direct solver and preconditioner STRUMPACK"
  "suitesparse:SUITESPARSE:direct solvers UMFPACK and KLU of SuiteSparse"
  "sundials:SUNDIALS:nonlinear solver KINSOL of SUNDIALS"
  "mkl-pardiso:MKL_PARDISO:direct solver PARDISO of oneMKL"
  "mkl-cpardiso:MKL_CPARDISO:direct solver Cluster Sparse Solver of oneMKL"
  "cudss:CUDSS:direct solver cuDSS"
  "magma:MAGMA:sparse solvers of MAGMA")

#! \brief selects the use of an optional library (option `enable-<name>`)
#! \param name: name of the option
#! \param available: true if the library can be used
#! \param reason: why the library can not be used (error message of `ON`)
#! \param description: description of the option
#! \param result: variable set to `ON` or `OFF` in the scope of the caller
function(mfem_mgis_select_optional_library name available reason description
         result)
  set(enable-${name} "AUTO" CACHE STRING
      "${description}: AUTO (if available), ON or OFF")
  set_property(CACHE enable-${name} PROPERTY STRINGS AUTO ON OFF)
  # the boolean constants of CMake are accepted for ON and OFF
  string(TOUPPER "${enable-${name}}" _choice)
  if(_choice MATCHES "^(ON|TRUE|YES|Y|1)$")
    set(_choice "ON")
  elseif(_choice MATCHES "^(OFF|FALSE|NO|N|0)$")
    set(_choice "OFF")
  endif()
  if(available)
    set(_available ON)
  else()
    set(_available OFF)
  endif()
  if(_choice STREQUAL "AUTO")
    set(_have ${_available})
  elseif(_choice STREQUAL "ON")
    if(NOT _available)
      message(FATAL_ERROR "mfem-mgis: enable-${name} is ON, but ${reason}")
    endif()
    set(_have ON)
  elseif(_choice STREQUAL "OFF")
    set(_have OFF)
  else()
    message(FATAL_ERROR "mfem-mgis: invalid value '${enable-${name}}' of "
            "enable-${name} (expected AUTO, ON or OFF)")
  endif()
  set(${result} ${_have} PARENT_SCOPE)
endfunction(mfem_mgis_select_optional_library)

set(_mfem_mgis_optional_libraries_summary)
foreach(_library ${MFEM_MGIS_OPTIONAL_MFEM_LIBRARIES})
  string(REPLACE ":" ";" _fields "${_library}")
  list(GET _fields 0 _name)
  list(GET _fields 1 _variable)
  list(GET _fields 2 _description)
  # the variables of MFEM are YES or NO (config.mk), or ON or OFF
  mfem_mgis_select_optional_library(${_name} "${MFEM_USE_${_variable}}"
    "MFEM is not built with it (MFEM_USE_${_variable})" "${_description}"
    MFEM_MGIS_HAVE_${_variable})
  list(APPEND _mfem_mgis_optional_libraries_summary
       "${_name}=${MFEM_MGIS_HAVE_${_variable}}")
endforeach(_library)
string(REPLACE ";" " " _mfem_mgis_optional_libraries_summary
       "${_mfem_mgis_optional_libraries_summary}")
message(STATUS "mfem-mgis optional libraries of MFEM: "
        "${_mfem_mgis_optional_libraries_summary}")
