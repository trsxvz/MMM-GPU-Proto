project(OperaHPCBubble)

# behavior law
file(COPY verification/mmm_vs_analytical.py DESTINATION ${PROJECT_BINARY_DIR})


# Behaviour library for the CPU and, if MFEM uses a GPU, for the GPU, chosen
# at run time by mfem_mgis::getBehavioursLibrary from the option --execution
mfem_mgis_add_behaviours_library(BehavioursBubble
  SOURCES ${PROJECT_SOURCE_DIR}/Elasticity.mfront)

# build executable

file( GLOB_RECURSE SRCS *.cxx )

add_executable(test-bubble
  ${SRCS}
)

add_dependencies(test-bubble BehavioursBubble)
target_compile_definitions(test-bubble PRIVATE
  MFEM_MGIS_EXAMPLE_BEHAVIOURS_DIRECTORY="${BehavioursBubble_BEHAVIOURS_DIRECTORY}")

target_link_libraries(test-bubble
 PRIVATE tfel::TFELMath tfel::TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(test-bubble
   PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# copy and install

file(COPY Elasticity.mfront mesh DESTINATION ${PROJECT_BINARY_DIR})

# the stress at the bubble (column Stress of the output file) is compared to
# the analytical hoop stress at the surface of a spherical cavity under the
# internal pressure p (1 by default), p / 2, with a tolerance of 20 % for
# the coarse mesh of the test (a null stress field fails). Both tests write
# the same output file.
set(_bubble_options --mesh mesh/mesh_sphere.msh
    --bubble-file mesh/single_bubble_ci.txt)
mm_opera_hpc_add_checked_test(NAME BubbleTest
  FILE bubbles_and_stresses_selected.txt CHECKS first Stress 0.5
  TOLERANCE 0.2
  COMMAND $<TARGET_FILE:test-bubble> ${_bubble_options})
set_tests_properties(BubbleTest PROPERTIES RESOURCE_LOCK bubble_output)
# the mean hydrostatic stress of the solid (column Average of sig_hydro.txt)
# is compared to the one of the composite sphere assemblage with a null
# macroscopic strain, -3 K f p / (3 K f + 4 mu), f being the volume
# fraction of the bubble (2.68e-4): -4.354e-4 p (written -0.0004354, which
# argparse reads as a number). This checks the sign and
# the magnitude of the pressure, the Poisson ratio and the volume of the
# cavity (-0.5 % on the coarse mesh of the test)
mm_opera_hpc_add_checked_test(NAME BubbleTest-hydrostatic-stress
  FILE sig_hydro.txt CHECKS first Average -0.0004354
  TOLERANCE 0.02
  COMMAND $<TARGET_FILE:test-bubble> ${_bubble_options})
set_tests_properties(BubbleTest-hydrostatic-stress PROPERTIES
                     RESOURCE_LOCK bubble_output)
if(MM_OPERA_HPC_RUN_GPU_TESTS AND BehavioursBubble_HAVE_GPU_BEHAVIOURS)
  mm_opera_hpc_add_checked_test(NAME BubbleTest-gpu
    FILE bubbles_and_stresses_selected.txt CHECKS first Stress 0.5
    TOLERANCE 0.2
    COMMAND $<TARGET_FILE:test-bubble> ${_bubble_options} --execution gpu)
  set_tests_properties(BubbleTest-gpu PROPERTIES RESOURCE_LOCK
                       "gpu;bubble_output")
endif()


install(TARGETS test-bubble DESTINATION bin)
