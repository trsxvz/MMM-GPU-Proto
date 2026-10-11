# Cermet of mm-opera-hpc: polycrystal of UO2 (5 grains, single crystal
# behaviour, orthotropic, finite strain, implicit) with a metallic interface
# (Norton behaviour), periodic problem, one time step on two processes. The
# executable is the one of mm-opera-hpc; the behaviours are also generated
# by the generic-parallel interface (one library for both). The macroscopic
# stresses are compared to the ones of the integration point by point.
set(MM_OPERA_HPC_CERMET_DIR "${MM_OPERA_HPC_BINARY_DIR}/cermet"
  CACHE PATH "build directory of the cermet of mm-opera-hpc")
add_generic_parallel_mfront_libraries(cermet-behaviours CUDA
  "${MM_OPERA_HPC_SOURCE_DIR}/cermet/law/NortonCr.mfront"
  "${MM_OPERA_HPC_SOURCE_DIR}/cermet/law/Mono_UO2_CosH_Jaco3.mfront")

if(NOT (Python3_Interpreter_FOUND AND
        EXISTS "${MM_OPERA_HPC_CERMET_DIR}/cermet"))
  return()
endif()
set(_run mpirun -n 2 "${MM_OPERA_HPC_CERMET_DIR}/cermet"
  --mesh "${MM_OPERA_HPC_CERMET_DIR}/mesh/5grains.msh"
  --file "${MM_OPERA_HPC_CERMET_DIR}/vectors_5grains.txt"
  --duration 0.4 --nstep 1 --postprocessing 0)
set(_kernels "MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch;MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels;MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels;MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels")
set(_reference "${CMAKE_CURRENT_BINARY_DIR}/cermet-point-par-point.res")

# usage: add_cermet_test(name environment [options...])
# Test cermet-<name>, writing <current binary dir>/cermet-<name>.res, and
# the comparison of its macroscopic stresses to the ones of the integration
# point by point
function(add_cermet_test name environment)
  set(_output "${CMAKE_CURRENT_BINARY_DIR}/cermet-${name}.res")
  add_test(NAME cermet-${name} COMMAND ${_run} ${ARGN}
    --macroscopic-stress-output-file "${_output}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties(cermet-${name} PROPERTIES
    ENVIRONMENT "${environment}"
    FIXTURES_SETUP cermet-${name})
  if(name STREQUAL "point-par-point")
    return()
  endif()
  add_test(NAME cermet-comparison-${name}
    COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/../comparer-sorties.py"
            "${_reference}" "${_output}" --rtol 1e-8 --groupes 0 1-3 4-6)
  set_tests_properties(cermet-comparison-${name} PROPERTIES
    FIXTURES_REQUIRED "cermet-point-par-point;cermet-${name}")
endfunction()

set(_generic "${MM_OPERA_HPC_CERMET_DIR}/src/libBehaviour.so")
add_cermet_test(point-par-point "MFEM_MGIS_BEHAVIOUR_INTEGRATION=PointByPoint"
  -l "${_generic}")
add_cermet_test(noyaux "${_kernels}" -l "${_generic}")
foreach(_backend ${cermet-behaviours_BACKENDS})
  set(_library "${CMAKE_CURRENT_BINARY_DIR}/generic-parallel-${_backend}/src/libGenericParallelBehaviour-${_backend}.so")
  if(_backend STREQUAL "cuda")
    # behaviours integrated on the GPU, all the computations by kernels
    add_cermet_test(batch-${_backend} "${_kernels}" --device cuda
      -l "${_library}")
    set_tests_properties(cermet-batch-${_backend} PROPERTIES
      RESOURCE_LOCK gpu)
  else()
    add_cermet_test(batch-${_backend}
      "MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch;OMP_NUM_THREADS=2"
      -l "${_library}")
  endif()
endforeach()
