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
add_test(NAME UniaxialPolycrystalTest COMMAND uniaxial-polycrystal --disable-post-processing --nstep 1 --duration 0.1)
if(MM_OPERA_HPC_RUN_GPU_TESTS)
  add_test(NAME UniaxialPolycrystalTest-gpu COMMAND uniaxial-polycrystal --disable-post-processing --nstep 1 --duration 0.1 --execution gpu)
  add_test(NAME UniaxialPolycrystalTest-mixed-gpu COMMAND uniaxial-polycrystal --disable-post-processing --nstep 1 --duration 0.1 --execution mixed --behaviours cpu)
  set_tests_properties(UniaxialPolycrystalTest-gpu UniaxialPolycrystalTest-mixed-gpu PROPERTIES RESOURCE_LOCK gpu)
endif()
add_test(NAME UniaxialPolycrystalTestPar COMMAND mpirun -n 2 ./uniaxial-polycrystal --disable-post-processing --nstep 1 --duration 0.1)
add_test(NAME PolycrystalTestWithPostProPar COMMAND mpirun -n 2 ./uniaxial-polycrystal --enable-export-von-Mises-stress --enable-export-first_eigen_stress --nstep 1 --duration 0.01)
