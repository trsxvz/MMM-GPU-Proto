if(MFEM_MGIS_HAVE_TFEL)
mfront_behaviours_check_library(BehaviourBenchmarks  IsotropicLinearElasticity)

  file(COPY
      ${PROJECT_SOURCE_DIR}/benchmarks/IsotropicLinearElasticity.mfront
      ${PROJECT_SOURCE_DIR}/benchmarks/beam-tet.mesh
      DESTINATION ${PROJECT_BINARY_DIR})

# MFEM-MGIS example
  add_executable(MFEMMGISLinearElasticityBenchmark
      MFEMMGISLinearElasticityBenchmark.cxx)
  add_custom_command(
      OUTPUT "src/libBehaviour.so"
      COMMAND "${MFRONT}" --obuild --interface=generic
              ${MFEM_MGIS_MFRONT_REAL_TYPE_OPTION}
              ${PROJECT_BINARY_DIR}/IsotropicLinearElasticity.mfront
      DEPENDS ${PROJECT_BINARY_DIR}/IsotropicLinearElasticity.mfront
              ${MFEM_MGIS_MFRONT_DEPENDENCIES}
      COMMENT "Performing mfront code generation")

  target_link_libraries(MFEMMGISLinearElasticityBenchmark
      PUBLIC MFEMMGIS)


  add_custom_target(generate-libBehaviour-benchmark ALL
      DEPENDS "src/libBehaviour.so")

add_dependencies(MFEMMGISLinearElasticityBenchmark generate-libBehaviour-benchmark)

# MFEM example
  add_executable(MFEMLinearElasticityBenchmark
      MFEMLinearElasticityBenchmark.cxx)

  target_link_libraries(MFEMLinearElasticityBenchmark
      PUBLIC MFEMMGIS)

endif()

# benchmark of the linear solvers of the factory on assembled systems
# (Laplace equation, linear elasticity of a beam, possibly heterogeneous)
if(MFEM_USE_MPI)
  add_executable(LinearSolversBenchmark LinearSolversBenchmark.cxx)
  target_link_libraries(LinearSolversBenchmark PUBLIC MFEMMGIS)
endif()
