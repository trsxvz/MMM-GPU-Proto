# Runs ProfilingMarkersTest on a GPU under rocprofv3, the profiler of ROCm,
# and checks that the ranges of roctx declared by the test (three steps,
# each one with three sections) are in the trace of the markers.
#
# usage: cmake -DROCPROFV3=<rocprofv3> -DTEST=<executable>
#          -DDIRECTORY=<directory of the trace> -P ProfilingMarkersTest.cmake
file(REMOVE_RECURSE "${DIRECTORY}")
execute_process(
  COMMAND "${ROCPROFV3}" --marker-trace --output-format csv
          -d "${DIRECTORY}" -o trace --
          "${TEST}" --device hip --expected yes
  RESULT_VARIABLE _status)
if(NOT _status EQUAL 0)
  message(FATAL_ERROR "ProfilingMarkersTest failed under rocprofv3 "
          "(${_status})")
endif()
set(_trace "${DIRECTORY}/trace_marker_api_trace.csv")
if(NOT EXISTS "${_trace}")
  message(FATAL_ERROR "no trace of the markers (${_trace})")
endif()
file(STRINGS "${_trace}" _lines)
foreach(_range "ProfilingMarkersTest: step 0=1"
               "ProfilingMarkersTest: step 1=1"
               "ProfilingMarkersTest: step 2=1"
               "ProfilingMarkersTest::Step=3"
               "ProfilingMarkersTest::Fill=3"
               "ProfilingMarkersTest::Norm=3")
  string(REGEX REPLACE "=[0-9]+$" "" _name "${_range}")
  string(REGEX REPLACE "^.*=" "" _expected "${_range}")
  set(_count 0)
  foreach(_line IN LISTS _lines)
    string(FIND "${_line}" "\"${_name}\"" _position)
    if(NOT _position EQUAL -1)
      math(EXPR _count "${_count} + 1")
    endif()
  endforeach()
  if(NOT _count EQUAL _expected)
    message(FATAL_ERROR "range '${_name}': ${_count} in the trace of the "
            "markers, ${_expected} expected")
  endif()
  message(STATUS "range '${_name}': ${_count} in the trace of the markers")
endforeach()
