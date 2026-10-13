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

add_test(NAME BubbleTest COMMAND test-bubble --mesh mesh/mesh_sphere.msh --bubble-file mesh/single_bubble_ci.txt)
if(MM_OPERA_HPC_RUN_GPU_TESTS)
  add_test(NAME BubbleTest-gpu COMMAND test-bubble --mesh mesh/mesh_sphere.msh --bubble-file mesh/single_bubble_ci.txt --execution gpu)
  set_tests_properties(BubbleTest-gpu PROPERTIES RESOURCE_LOCK gpu)
endif()


install(TARGETS test-bubble DESTINATION bin)
