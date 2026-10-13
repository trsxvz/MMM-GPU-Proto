cmake_minimum_required(VERSION 3.13)
project(ssna303_3dExample)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
# testing
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()
if(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND
      ${CMAKE_CTEST_COMMAND} -T test )
   endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)

# Require MFEM-MGIS
if (NOT TARGET mfem-mgis::MFEMMGIS)
  find_package(MFEMMGIS REQUIRED)
endif(NOT TARGET mfem-mgis::MFEMMGIS)

# Test mode of the examples, see ../config/TestMode.cmake
if(NOT DEFINED MFEM_MGIS_EXAMPLES_FULL_TESTS)
  include(${CMAKE_CURRENT_SOURCE_DIR}/../config/TestMode.cmake)
endif()


# Copy the files needed at run time. Unlike file(COPY), configure_file updates
# the copies when these files are modified.
foreach(input IsotropicLinearHardeningPlasticity.mfront ssna303_3d.msh rc_ex10p
              ssna303_3d-force-1step.ref ssna303_3d-force-10steps.ref
              ssna303_3d-force-10steps-fbar.ref)
  configure_file(${PROJECT_SOURCE_DIR}/${input} ${PROJECT_BINARY_DIR}/${input}
                 COPYONLY)
endforeach()

# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU, chosen
# at run time by mfem_mgis::getBehavioursLibrary from the option --execution
mfem_mgis_add_behaviours_library(BehavioursEx4
  SOURCES ${PROJECT_BINARY_DIR}/IsotropicLinearHardeningPlasticity.mfront)

# Define a test that will be called during the "make check"
function(add_ssna303_3d_test_ex solver nsteps end_time reference_file use_fbar)
 set(test "ssna303_3d_${solver}")
 add_test(NAME ${test} COMMAND ${test} --nbsteps ${nsteps} --end-time ${end_time}
                                  --reference-file ${reference_file})
 # one OpenMP thread per process, otherwise MUMPS oversubscribes the cores
 set_tests_properties(${test} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
 if(${use_fbar})
   if(MGIS_HAVE_TFEL)
     string(REPLACE ".ref" "-fbar.ref" fbar_reference_file ${reference_file})
     add_test(NAME ${test}-fbar COMMAND ${test} --use-fbar
                                                 --nbsteps ${nsteps} --end-time ${end_time}
                                                 --reference-file ${fbar_reference_file}
                                                 --standard-reference-file ${reference_file})
     set_tests_properties(${test}-fbar PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
   endif(MGIS_HAVE_TFEL)
 endif()
endfunction(add_ssna303_3d_test_ex)

# The full test mode computes the complete loading in 10 time steps. The
# restricted test mode only computes its first time step, during which the
# plastic flow starts. In both modes, the results are compared to the first
# values of the reference files. The PETSc case always computes a single small
# time step.
if(MFEM_MGIS_EXAMPLES_FULL_TESTS)
  set(ssna303_3d_nsteps 10)
  set(ssna303_3d_end_time 1)
else()
  set(ssna303_3d_nsteps 1)
  set(ssna303_3d_end_time 0.1)
endif()

foreach(tg IN ITEMS mumps hypre petsc)
  # Define the executable to build for this example
  if(((NOT MFEM_USE_PETSC) AND (${tg} STREQUAL "petsc")) OR
      ((NOT MFEM_USE_MUMPS) AND (${tg} STREQUAL "mumps")))
  else()
    add_executable("ssna303_3d_${tg}" ssna303_${tg}.cxx)
    target_compile_definitions(ssna303_3d_${tg} PRIVATE
      MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursEx4_BEHAVIOURS_DIRECTORY}")
    target_link_libraries(ssna303_3d_${tg} mfem-mgis::MFEMMGIS)
    target_include_directories(ssna303_3d_${tg}
      SYSTEM
      PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
      PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
    add_dependencies(check ssna303_3d_${tg} BehavioursEx4)
    if(${tg} STREQUAL "petsc")
      add_ssna303_3d_test_ex(${tg} 1 0.02 ssna303_3d-force-1step.ref OFF)
    else()
      if(${tg} STREQUAL "mumps")
        add_ssna303_3d_test_ex(${tg} ${ssna303_3d_nsteps} ${ssna303_3d_end_time}
                               ssna303_3d-force-10steps.ref ON)
      else()
        add_ssna303_3d_test_ex(${tg} ${ssna303_3d_nsteps} ${ssna303_3d_end_time}
                               ssna303_3d-force-10steps.ref OFF)
      endif()
    endif()
  endif()
endforeach()

# the same computations on a GPU (see ../config/TestMode.cmake): hypre on the
# GPU, and MUMPS on the CPU (option --linear-solver-device cpu)
if(MFEM_MGIS_EXAMPLES_RUN_GPU_TESTS AND BehavioursEx4_HAVE_GPU_BEHAVIOURS)
  add_test(NAME ssna303_3d_hypre-gpu
           COMMAND ssna303_3d_hypre --nbsteps ${ssna303_3d_nsteps}
                   --end-time ${ssna303_3d_end_time}
                   --reference-file ssna303_3d-force-10steps.ref
                   --execution gpu)
  set_tests_properties(ssna303_3d_hypre-gpu PROPERTIES RESOURCE_LOCK gpu)
  if(MFEM_USE_MUMPS)
    add_test(NAME ssna303_3d_mumps-gpu
             COMMAND ssna303_3d_mumps --nbsteps ${ssna303_3d_nsteps}
                     --end-time ${ssna303_3d_end_time}
                     --reference-file ssna303_3d-force-10steps.ref
                     --execution mixed --linear-solver-device cpu)
    set_tests_properties(ssna303_3d_mumps-gpu PROPERTIES
                         RESOURCE_LOCK gpu ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
endif()
