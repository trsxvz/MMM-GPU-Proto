project(OperaHPCCermetExample)

# Copy
file(COPY law mesh results DESTINATION ${PROJECT_BINARY_DIR})
file(COPY mesh/vectors_5grains.txt check_cermet_restults.py plot_cermet_results.py DESTINATION ${PROJECT_BINARY_DIR})

# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU, chosen
# at run time by mfem_mgis::getBehavioursLibrary from the option --execution
mfem_mgis_add_behaviours_library(BehavioursCermet
  SOURCES ${PROJECT_SOURCE_DIR}/law/NortonCr.mfront
          ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco.mfront
          ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco3.mfront)

## Build executable
add_executable(cermet cermet.cxx)
target_link_libraries(cermet
 PUBLIC MMOPERAHPC)
add_dependencies(cermet BehavioursCermet)
target_compile_definitions(cermet PRIVATE
  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursCermet_BEHAVIOURS_DIRECTORY}")

## Test
add_dependencies(check cermet BehavioursCermet)
# two time steps of 0.2 s, the axial stress Szz (column 7 of the output
# file) being compared to the one computed by Cast3M at each step
# (results/cermet_cast3m.csv). In single precision, the stress at 0.2 s
# (axial strain of 1e-4) depends on the compiler of the behaviour (gap to
# the solution in double: 0.15 % with gcc 13, 2.7 % with gcc 14), the
# deformation gradient, stored in float close to the identity, being only
# known within one ulp (1.2e-7, 1.2e-3 of the strain): only the stress at
# 0.4 s is checked (gap below 0.25 % with both compilers).
get_property(_real_type GLOBAL PROPERTY MFEM_MGIS_REAL_TYPE)
if(_real_type STREQUAL "float")
  set(_cermet_checks CHECKS 0.4 7 2.9462e7 TOLERANCE 0.01)
else()
  set(_cermet_checks CHECKS 0.2 7 1.7246e7 0.4 7 2.9462e7 TOLERANCE 0.01)
endif()
set(_cermet_options --duration 0.4 --nstep 2 --postprocessing 0)
mm_opera_hpc_add_checked_test(NAME CermetNullStrainXXYY
  FILE cermet-test.res ${_cermet_checks}
  COMMAND $<TARGET_FILE:cermet> ${_cermet_options}
          --macroscopic-stress-output-file cermet-test.res)
add_test(NAME CermetImposeGradXXYY COMMAND cermet --duration 0.4 --nstep 1 --postprocessing 0)
if(MM_OPERA_HPC_RUN_GPU_TESTS AND BehavioursCermet_HAVE_GPU_BEHAVIOURS)
  mm_opera_hpc_add_checked_test(NAME CermetNullStrainXXYY-gpu
    FILE cermet-test-gpu.res ${_cermet_checks}
    COMMAND $<TARGET_FILE:cermet> ${_cermet_options}
            --macroscopic-stress-output-file cermet-test-gpu.res
            --execution gpu)
  set_tests_properties(CermetNullStrainXXYY-gpu PROPERTIES RESOURCE_LOCK gpu)
endif()

## Install
install(TARGETS cermet DESTINATION bin)
