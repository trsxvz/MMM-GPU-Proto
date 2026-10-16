project(UniaxialPolycrystalExample)

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco.mfront
  ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco3.mfront
  DESTINATION ${PROJECT_BINARY_DIR})

file(COPY law mesh results plot_polycrystal_results.py check_polycrystal_restults.py DESTINATION ${PROJECT_BINARY_DIR})

# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU, chosen
# at run time by mfem_mgis::getBehavioursLibrary from the option --execution
mfem_mgis_add_behaviours_library(BehavioursPolycrystal
  SOURCES ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco.mfront
          ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco3.mfront)

add_executable(uniaxial-polycrystal UniaxialPolycrystal.cxx)
target_link_libraries(uniaxial-polycrystal MMOPERAHPC)
target_include_directories(uniaxial-polycrystal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(uniaxial-polycrystal PRIVATE
  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursPolycrystal_BEHAVIOURS_DIRECTORY}")
    
# testing
add_dependencies(check uniaxial-polycrystal BehavioursPolycrystal)
# one time step of 0.2 s, the axial stress Szz (column 7 of the output file)
# being compared to the one computed by Cast3M (results/polycrystal_cast3m.csv)
set(_polycrystal_checks CHECKS 0.2 7 1.7994e7 TOLERANCE 0.01)
set(_polycrystal_options --disable-post-processing --nstep 1 --duration 0.2)
mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest
  FILE polycrystal-test.res ${_polycrystal_checks}
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test.res)
# the same with the displacement gradient given to the behaviour
mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest-DisplacementGradient
  FILE polycrystal-test-dg.res ${_polycrystal_checks}
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-dg.res
          --displacement-gradient)
if(MM_OPERA_HPC_RUN_GPU_TESTS AND BehavioursPolycrystal_HAVE_GPU_BEHAVIOURS)
  mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest-gpu
    FILE polycrystal-test-gpu.res ${_polycrystal_checks}
    COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
            --macroscopic-stress-output-file polycrystal-test-gpu.res
            --execution gpu)
  mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest-mixed-gpu
    FILE polycrystal-test-mixed-gpu.res ${_polycrystal_checks}
    COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
            --macroscopic-stress-output-file polycrystal-test-mixed-gpu.res
            --execution mixed --behaviours-execution cpu)
  set_tests_properties(UniaxialPolycrystalTest-gpu UniaxialPolycrystalTest-mixed-gpu PROPERTIES RESOURCE_LOCK gpu)
endif()
mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTestPar
  FILE polycrystal-test-par.res ${_polycrystal_checks}
  COMMAND mpirun -n 2 $<TARGET_FILE:uniaxial-polycrystal>
          ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-par.res)
add_test(NAME PolycrystalTestWithPostProPar COMMAND mpirun -n 2 ./uniaxial-polycrystal --enable-export-von-Mises-stress --enable-export-first_eigen_stress --nstep 1 --duration 0.01)
