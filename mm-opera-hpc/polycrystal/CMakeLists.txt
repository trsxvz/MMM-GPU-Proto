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
# the maximum number of iterations of the non linear solver is taken into
# account. With NewtonK0 and the tolerance of Eisenstat and Walker, whose
# first linear solve is loose (the residual is about halved), a resolution
# of this time step takes four iterations: with one iteration at most, the
# first resolution fails and the time step is divided, which is expected
# (in double, no resolution converges, and the temporal sequence fails
# after its ten divisions; in float, the absolute tolerance on the residual
# accepts the smaller time steps). A second iteration would reveal an
# ignored option.
add_test(NAME UniaxialPolycrystalTest-MaximumIterations
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-iterations.res
          --non-linear-solver NewtonK0
          --adaptive-linear-solver-tolerance EisenstatWalker1
          --non-linear-solver-maximum-iterations 1)
set_tests_properties(UniaxialPolycrystalTest-MaximumIterations PROPERTIES
  PASS_REGULAR_EXPRESSION "dividing time step by 2"
  FAIL_REGULAR_EXPRESSION "Fixed point iteration  2")
add_test(NAME PolycrystalTestWithPostProPar COMMAND mpirun -n 2 ./uniaxial-polycrystal --enable-export-von-Mises-stress --enable-export-first_eigen_stress --nstep 1 --duration 0.01)
