# Test mode of the examples, selected by the MFEM_MGIS_EXAMPLES_TEST_MODE
# option:
#
# - full: the tests compute the complete simulations;
# - restricted: the tests only compute the beginning of the simulations, or
#   their smallest meaningful case, which is much faster;
# - auto, the default: restricted in the Debug and Coverage builds, which are
#   much slower, and full otherwise.
#
# The MFEM_MGIS_EXAMPLES_FULL_TESTS variable is true in the full test mode.
# Tests which are already small are the same in both modes.

set(MFEM_MGIS_EXAMPLES_TEST_MODE "auto" CACHE STRING
    "Test mode of the examples: auto, full or restricted")
set_property(CACHE MFEM_MGIS_EXAMPLES_TEST_MODE
             PROPERTY STRINGS auto full restricted)

if(MFEM_MGIS_EXAMPLES_TEST_MODE STREQUAL "auto")
  if(CMAKE_BUILD_TYPE MATCHES "^(Debug|Coverage)$")
    set(MFEM_MGIS_EXAMPLES_FULL_TESTS OFF)
  else()
    set(MFEM_MGIS_EXAMPLES_FULL_TESTS ON)
  endif()
elseif(MFEM_MGIS_EXAMPLES_TEST_MODE STREQUAL "full")
  set(MFEM_MGIS_EXAMPLES_FULL_TESTS ON)
elseif(MFEM_MGIS_EXAMPLES_TEST_MODE STREQUAL "restricted")
  set(MFEM_MGIS_EXAMPLES_FULL_TESTS OFF)
else()
  message(FATAL_ERROR "invalid test mode '${MFEM_MGIS_EXAMPLES_TEST_MODE}', "
                      "MFEM_MGIS_EXAMPLES_TEST_MODE must be auto, full or "
                      "restricted")
endif()

if(MFEM_MGIS_EXAMPLES_FULL_TESTS)
  message(STATUS "Test mode of the examples: full")
else()
  message(STATUS "Test mode of the examples: restricted")
endif()

# Tests on a GPU, selected by the MFEM_MGIS_EXAMPLES_GPU_TESTS option (off by
# default, a GPU being required when the tests are run): the examples which
# support it are also tested with the option --execution gpu (or mixed), if
# MFEM is compiled with CUDA or HIP (as described by mfem-mgis, which must be
# found before this file is included). The MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS
# variable is then true.
option(MFEM_MGIS_EXAMPLES_GPU_TESTS "Test the examples on a GPU too" OFF)
set(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS OFF)
if(MFEM_MGIS_EXAMPLES_GPU_TESTS)
  get_property(_mfem_mgis_examples_mfem_use_cuda GLOBAL
               PROPERTY MFEM_MGIS_MFEM_USE_CUDA)
  get_property(_mfem_mgis_examples_mfem_use_hip GLOBAL
               PROPERTY MFEM_MGIS_MFEM_USE_HIP)
  if((_mfem_mgis_examples_mfem_use_cuda STREQUAL "YES") OR
     (_mfem_mgis_examples_mfem_use_hip STREQUAL "YES"))
    set(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS ON)
    message(STATUS "Tests of the examples on a GPU: yes")
  else()
    message(STATUS "Tests of the examples on a GPU: none, MFEM is not "
                   "compiled with CUDA nor HIP")
  endif()
endif()
