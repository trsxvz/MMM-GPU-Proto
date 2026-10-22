# Checks of the execution policies of the generic-parallel interfaces across
# libraries (see README.md). No dependency: CUDA, SYCL and Kokkos are
# replaced by stand-ins.
cmake_minimum_required(VERSION 3.20)
project(execution-policies LANGUAGES C CXX)

option(EXECUTION_POLICIES_IDENTIFIER_FIRST
  "declare getIdentifier before the destructor of the base class" OFF)
option(EXECUTION_POLICIES_HIDDEN_VISIBILITY
  "hide the symbols by default, as the flags given by tfel-config" ON)
set(EXECUTION_POLICIES_OTHER_GENERATED_LIBRARY "" CACHE FILEPATH
  "generated library built by another toolchain (test cross)")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
if(EXECUTION_POLICIES_HIDDEN_VISIBILITY)
  set(CMAKE_C_VISIBILITY_PRESET hidden)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
# all the libraries and programs in one directory (one per configuration
# with Visual Studio), where Windows looks for the libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra)
endif()
if(EXECUTION_POLICIES_IDENTIFIER_FIRST)
  add_compile_definitions(MFRONT_GENERICPARALLEL_IDENTIFIER_FIRST)
endif()
include_directories(include stand-ins src)

# runtime library: the base class of the policies
add_library(MFrontGenericParallelRuntime SHARED
  src/AbstractBackendExecutionPolicy.cxx)

# generated library, not linked to the runtime library, as the libraries
# generated today: the link fails if a symbol of this library is needed
add_library(generated MODULE src/generated.cxx)
if(CMAKE_SYSTEM_NAME MATCHES "Linux|FreeBSD")
  target_link_options(generated PRIVATE "LINKER:--no-undefined")
endif()
# the same, linked to the runtime library
add_library(generated-linked MODULE src/generated.cxx)
target_link_libraries(generated-linked PRIVATE MFrontGenericParallelRuntime)
# the same checks by dynamic_cast, for information
add_library(dynamic-cast MODULE src/dynamic-cast.cxx)
target_link_libraries(dynamic-cast PRIVATE MFrontGenericParallelRuntime)

# caller (MGIS, a user code), linked to the runtime library
add_library(caller SHARED src/caller.cxx)
target_link_libraries(caller PRIVATE MFrontGenericParallelRuntime
  ${CMAKE_DL_LIBS})
add_executable(execution-policies-linked src/main-linked.cxx)
target_link_libraries(execution-policies-linked PRIVATE caller)
add_executable(execution-policies-loaded src/main-loaded.cxx)
target_link_libraries(execution-policies-loaded PRIVATE ${CMAKE_DL_LIBS})
add_dependencies(execution-policies-loaded caller generated dynamic-cast)
add_executable(execution-policies-c src/main-c.c)
target_link_libraries(execution-policies-c PRIVATE ${CMAKE_DL_LIBS})
add_dependencies(execution-policies-c generated)

enable_testing()
add_test(NAME linked
  COMMAND execution-policies-linked $<TARGET_FILE:generated>
          $<TARGET_FILE:dynamic-cast>)
add_test(NAME linked-generated-linked
  COMMAND execution-policies-linked $<TARGET_FILE:generated-linked>)
add_test(NAME loaded
  COMMAND execution-policies-loaded $<TARGET_FILE:caller>
          $<TARGET_FILE:generated> $<TARGET_FILE:dynamic-cast>)
add_test(NAME c
  COMMAND execution-policies-c $<TARGET_FILE:generated>)
if(EXECUTION_POLICIES_OTHER_GENERATED_LIBRARY)
  add_test(NAME cross
    COMMAND execution-policies-linked
            "${EXECUTION_POLICIES_OTHER_GENERATED_LIBRARY}" - c-types)
endif()
