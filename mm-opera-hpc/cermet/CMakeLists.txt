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
add_test(NAME CermetNullStrainXXYY COMMAND cermet --duration 0.4 --nstep 1 --postprocessing 0)
add_test(NAME CermetImposeGradXXYY COMMAND cermet --duration 0.4 --nstep 1 --postprocessing 0)
if(MM_OPERA_HPC_RUN_GPU_TESTS AND BehavioursCermet_HAVE_GPU_BEHAVIOURS)
  add_test(NAME CermetNullStrainXXYY-gpu COMMAND cermet --duration 0.4 --nstep 1 --postprocessing 0 --execution gpu)
  set_tests_properties(CermetNullStrainXXYY-gpu PROPERTIES RESOURCE_LOCK gpu)
endif()

## Install
install(TARGETS cermet DESTINATION bin)
