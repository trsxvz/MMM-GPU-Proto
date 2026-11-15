# Runs GlobalStiffnessMatrixTest on the host, then on a GPU, and checks that
# the checksums of the global matrices are the same: the values assembled on
# the GPU are the ones of the host, bit for bit. With MPI=gpu-aware-mpi, the
# values exchanged between the processes are given to MPI on the GPU.
#
# usage: cmake -DNP=<processes> -DDEVICE=<cuda|hip> -DMPI=<host|gpu-aware-mpi>
#          -DTEST=<executable> -DMESH=<mesh>
#          -P GlobalStiffnessMatrixTest.cmake
set(ENV{OMP_NUM_THREADS} 1)
if(MPI STREQUAL "gpu-aware-mpi")
  set(ENV{MFEM_MGIS_GPU_AWARE_MPI} 1)
endif()
foreach(_device cpu ${DEVICE})
  set(ENV{MFEM_MGIS_DEVICE} ${_device})
  execute_process(
    COMMAND mpirun -n ${NP} "${TEST}" --mesh "${MESH}"
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE _output)
  message("${_output}")
  if(NOT _status EQUAL 0)
    message(FATAL_ERROR "GlobalStiffnessMatrixTest failed on the device "
            "${_device} (${_status})")
  endif()
  string(REGEX MATCHALL "[^\n]*: checksum [0-9a-f]+" _checksums_${_device}
         "${_output}")
endforeach()
list(LENGTH _checksums_cpu _n)
if(_n EQUAL 0)
  message(FATAL_ERROR "no checksum in the output of GlobalStiffnessMatrixTest")
endif()
if(NOT "${_checksums_cpu}" STREQUAL "${_checksums_${DEVICE}}")
  message(FATAL_ERROR "the global matrices assembled on the host and on the "
          "device ${DEVICE} differ:\n${_checksums_cpu}\n"
          "${_checksums_${DEVICE}}")
endif()
message(STATUS "${_n} global matrices: same values on the host and on the "
        "device ${DEVICE}")
