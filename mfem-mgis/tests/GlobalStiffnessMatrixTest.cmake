# Runs GlobalStiffnessMatrixTest on the host with one thread, then in other
# settings, and checks that the checksums of the global matrices are the
# same, that is that their values are the same bit for bit:
# - with DEVICE, on a GPU: the values assembled on the GPU are the ones of
#   the host. With MPI=gpu-aware-mpi, the values exchanged between the
#   processes are given to MPI on the GPU;
# - with THREADS (mfem-mgis built with the option enable-openmp), with
#   this number of threads of OpenMP, on the host and on the GPU if any:
#   the setup shared between the threads gives the same matrices.
#
# usage: cmake -DNP=<processes> [-DDEVICE=<cuda|hip>]
#          [-DMPI=<host|gpu-aware-mpi>] [-DTHREADS=<threads>]
#          -DTEST=<executable> -DMESH=<mesh>
#          -P GlobalStiffnessMatrixTest.cmake
if(MPI STREQUAL "gpu-aware-mpi")
  set(ENV{MFEM_MGIS_GPU_AWARE_MPI} 1)
endif()
# settings: device and number of threads, the first one being the reference
set(_settings cpu:1)
if(DEVICE)
  list(APPEND _settings ${DEVICE}:1)
endif()
if(THREADS)
  list(APPEND _settings cpu:${THREADS})
  if(DEVICE)
    list(APPEND _settings ${DEVICE}:${THREADS})
  endif()
endif()
foreach(_setting ${_settings})
  string(REPLACE ":" ";" _fields ${_setting})
  list(GET _fields 0 _device)
  list(GET _fields 1 _threads)
  set(ENV{MFEM_MGIS_DEVICE} ${_device})
  set(ENV{OMP_NUM_THREADS} ${_threads})
  execute_process(
    COMMAND mpirun -n ${NP} "${TEST}" --mesh "${MESH}"
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE _output)
  message("${_output}")
  if(NOT _status EQUAL 0)
    message(FATAL_ERROR "GlobalStiffnessMatrixTest failed on the device "
            "${_device} with ${_threads} thread(s) (${_status})")
  endif()
  string(REGEX MATCHALL "[^\n]*: checksum [0-9a-f]+" _checksums "${_output}")
  if(NOT DEFINED _reference)
    set(_reference "${_checksums}")
    list(LENGTH _reference _n)
    if(_n EQUAL 0)
      message(FATAL_ERROR
              "no checksum in the output of GlobalStiffnessMatrixTest")
    endif()
  elseif(NOT "${_checksums}" STREQUAL "${_reference}")
    message(FATAL_ERROR "the global matrices assembled on the host with one "
            "thread and on the device ${_device} with ${_threads} thread(s) "
            "differ:\n${_reference}\n${_checksums}")
  endif()
endforeach()
message(STATUS "${_n} global matrices: same values in the settings "
        "${_settings} (device:threads)")
