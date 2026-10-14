# Configures the Doxyfile of mfem-mgis without configuring the project (the
# documentation does not depend on its dependencies). Variables given by -D:
# PROJECT_SOURCE_DIR (sources of mfem-mgis), PROJECT_BINARY_DIR (output).
# usage: cmake -DPROJECT_SOURCE_DIR=... -DPROJECT_BINARY_DIR=... -P
#        scripts/site-doxyfile.cmake
file(READ "${PROJECT_SOURCE_DIR}/CMakeLists.txt" _cmakelists)
string(REGEX MATCH
  "project\\(([a-z-]+)[ \n]+VERSION ([0-9.]+)[ \n]+DESCRIPTION \"([^\"]*)\""
  _project "${_cmakelists}")
if(NOT _project)
  message(FATAL_ERROR "project() of mfem-mgis not found")
endif()
set(CMAKE_PROJECT_NAME "${CMAKE_MATCH_1}")
set(CMAKE_PROJECT_VERSION "${CMAKE_MATCH_2}")
set(CMAKE_PROJECT_DESCRIPTION "${CMAKE_MATCH_3}")
configure_file("${PROJECT_SOURCE_DIR}/Doxyfile.in"
  "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
# output directory of the Doxyfile, which doxygen does not create
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/docs/doxygen")
