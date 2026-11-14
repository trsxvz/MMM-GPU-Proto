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
# the same with one material for all the grains (option --single-material),
# whose number is written when the option is taken into account (expected in
# the output, which reveals an ignored option; ctest then ignores the exit
# status: a mismatch, a failure of the command or an error of the scripts
# of cmake/ is recognized by its message)
mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest-SingleMaterial
  FILE polycrystal-test-single-material.res ${_polycrystal_checks}
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-single-material.res
          --single-material)
# one material for all the grains gives the results of one material per
# grain: two time steps of 0.2 s with the averages of the stresses per grain
# (post-processings without the ParaView export), with and without the
# option, on one process and on two. Compared: the main output file
# (time, diagonal components of the deformation gradient, diagonal
# components of the stress, relatively to the axial one) and the averages
# per grain (file avgStressPolycrystal, relatively to the largest component
# at each time).
# In double, the two computations only differ by the order of some sums:
# relative differences of 1e-16 were measured on the main output file of
# this case (4e-17 on one process, 1e-16 on two), and 3e-13 for 100 grains
# on 96 processes. Its 14 digits are compared with a tolerance of 1e-9,
# which still accepts linear solves (relative tolerance of 1e-12) stopping
# one iteration apart. The averages per grain are written with 6 digits,
# hence the tolerance of 2e-5 (twice a unit of the last digit); they were
# identical. In float, the stresses are only known within 3e4 Pa or so (see
# the tolerances of the simulation), for an axial stress of 1.8e7 Pa at the
# first time step: tolerance of 1 %, the one of the comparison to Cast3M.
if("${MFEMMGIS_REAL_TYPE}" STREQUAL "float")
  set(_polycrystal_comparisons polycrystal.res 1,2-4,5-7 1e-2
                               avgStressPolycrystal 1,2- 1e-2)
else()
  set(_polycrystal_comparisons polycrystal.res 1,2-4,5-7 1e-9
                               avgStressPolycrystal 1,2- 2e-5)
endif()
set(_polycrystal_comparison_options
    --disable-paraview-export --nstep 2 --duration 0.4
    --mesh ${CMAKE_CURRENT_SOURCE_DIR}/mesh/5crystals.msh
    --vect ${CMAKE_CURRENT_SOURCE_DIR}/mesh/vectors_5crystals.txt
    --macroscopic-stress-output-file polycrystal.res)
mm_opera_hpc_add_compared_test(
  NAME UniaxialPolycrystalTest-SingleMaterial-Comparison
  COMPARISONS ${_polycrystal_comparisons}
  OPTIONS --single-material
  COMMAND $<TARGET_FILE:uniaxial-polycrystal>
          ${_polycrystal_comparison_options})
mm_opera_hpc_add_compared_test(
  NAME UniaxialPolycrystalTestPar-SingleMaterial-Comparison
  COMPARISONS ${_polycrystal_comparisons}
  OPTIONS --single-material
  COMMAND mpirun -n 2 $<TARGET_FILE:uniaxial-polycrystal>
          ${_polycrystal_comparison_options})
set_tests_properties(UniaxialPolycrystalTest-SingleMaterial
  UniaxialPolycrystalTest-SingleMaterial-Comparison
  UniaxialPolycrystalTestPar-SingleMaterial-Comparison PROPERTIES
  PASS_REGULAR_EXPRESSION "Nombre de grains : +5"
  FAIL_REGULAR_EXPRESSION "MISMATCH;run-and-check.py: ;run-and-compare.py: ;no row ;Traceback")
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
# line search of the non linear solver (NewtonK0, tolerance of Eisenstat and
# Walker): same axial stress
mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest-LineSearch
  FILE polycrystal-test-line-search.res ${_polycrystal_checks}
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-line-search.res
          --non-linear-solver NewtonK0
          --adaptive-linear-solver-tolerance EisenstatWalker1
          --line-search Backtracking)
# an unknown line search is refused
add_test(NAME UniaxialPolycrystalTest-LineSearch-Invalid
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-line-search-invalid.res
          --line-search Unknown)
set_tests_properties(UniaxialPolycrystalTest-LineSearch-Invalid PROPERTIES
  PASS_REGULAR_EXPRESSION "error: invalid line search")
# prediction by linear extrapolation of the previous time step (NewtonK0,
# tolerance of Eisenstat and Walker, line search), 4 time steps of 0.05 s:
# same axial stress. From the second time step, the first resolution starts
# from the extrapolation, whose residual is lower than the one of the default
# prediction, the reference of the relative criterion: the ratio written at
# the iteration 0 is lower than one, whereas it is always one without
# prediction (expected in the output, which reveals an ignored option; a
# mismatch of the axial stress or a failure of the command is reported by
# cmake/run-and-check.py)
mm_opera_hpc_add_checked_test(NAME UniaxialPolycrystalTest-LinearExtrapolation
  FILE polycrystal-test-linear-extrapolation.res ${_polycrystal_checks}
  COMMAND $<TARGET_FILE:uniaxial-polycrystal>
          --disable-post-processing --nstep 4 --duration 0.2
          --macroscopic-stress-output-file polycrystal-test-linear-extrapolation.res
          --non-linear-solver NewtonK0
          --adaptive-linear-solver-tolerance EisenstatWalker1
          --line-search Backtracking
          --prediction LinearExtrapolation)
set_tests_properties(UniaxialPolycrystalTest-LinearExtrapolation PROPERTIES
  PASS_REGULAR_EXPRESSION "iteration  0 : \\|\\|r\\|\\| = [^,]+, \\|\\|r\\|\\|/\\|\\|r_0\\|\\| = (0\\.[0-9]|[0-9.]+e-)"
  FAIL_REGULAR_EXPRESSION "MISMATCH;run-and-check.py: ;no row ")
# an unknown prediction is refused
add_test(NAME UniaxialPolycrystalTest-Prediction-Invalid
  COMMAND $<TARGET_FILE:uniaxial-polycrystal> ${_polycrystal_options}
          --macroscopic-stress-output-file polycrystal-test-prediction-invalid.res
          --prediction Unknown)
set_tests_properties(UniaxialPolycrystalTest-Prediction-Invalid PROPERTIES
  PASS_REGULAR_EXPRESSION "error: invalid prediction")
add_test(NAME PolycrystalTestWithPostProPar COMMAND mpirun -n 2 ./uniaxial-polycrystal --enable-export-von-Mises-stress --enable-export-first_eigen_stress --nstep 1 --duration 0.01)
