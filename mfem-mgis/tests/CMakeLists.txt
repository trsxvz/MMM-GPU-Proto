# Helper function to create a test executable and add it to build-tests
function(add_test_executable_base name libs)
  add_executable(${name} EXCLUDE_FROM_ALL ${name}.cxx)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${name} PRIVATE "${libs}")
  add_dependencies(build-tests ${name})
endfunction()

function(add_test_executable name)
  add_test_executable_base("${name}" "mfem-mgis::MFEMMGIS")
endfunction()

function(add_tfel_test_executable name)
  if(MFEM_MGIS_ENABLE_TFEL_TESTS)
    add_test_executable_base("${name}" "mfem-mgis::MFEMMGIS;tfel::TFELTests")
  endif(MFEM_MGIS_ENABLE_TFEL_TESTS)
endfunction()

# Helper function to add a simple test case (serial only)
function(add_simple_test name)
  add_test_executable(${name})
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# Helper function to add a simple test case based on TFEL::tests (serial only)
function(add_simple_tfel_test name)
  if(MFEM_MGIS_ENABLE_TFEL_TESTS)
    add_tfel_test_executable(${name})
    add_test(NAME ${name} COMMAND ${name})
  endif(MFEM_MGIS_ENABLE_TFEL_TESTS)
endfunction()

# Helper function to add a test case with MPI support Usage: add_mpi_test(name
# parameters [NO_PARALLEL])
function(add_mpi_test name parameters)
  if(MFEM_USE_MPI)
    add_test_executable(${name})
    add_test(NAME ${name} COMMAND ${name} ${${parameters}} "--parallel" "0")
    add_test(NAME ${name}P COMMAND mpirun -n 2 ./${name} ${${parameters}}
                                   "--parallel" "1")
  endif(MFEM_USE_MPI)
endfunction()

# Simple test case (no TFEL dependency)
add_simple_test(UnitTestParameters)
add_simple_tfel_test(ParametersValidatorTest)
add_simple_tfel_test(TimeDiscretizationTest)
add_tfel_test_executable(MaximumNumberOfTimeStepsTest)
add_test(NAME MaximumNumberOfTimeStepsTest COMMAND ./MaximumNumberOfTimeStepsTest
         "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh")
add_tfel_test_executable(UnreliableResultsTest)
add_test(NAME UnreliableResultsTest COMMAND ./UnreliableResultsTest
         "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh")
add_simple_tfel_test(QPEvaluatorTest)


if(MFEM_MGIS_HAVE_TFEL)
  # Compile MFront behaviours library
  mfront_behaviours_check_library(
    BehaviourTest
    IsotropicLinearElasticity
    OrthotropicElasticity
    IsotropicLinearHardeningPlasticity
    IsotropicLinearHardeningThermoPlasticity
    Mazars
    SaintVenantKirchhoffElasticity
    SaintVenantKirchhoffElasticityV2
    MiehePhaseFieldDamage
    MicromorphicDamageI_SpectralSplit
    MicromorphicDamageII
    StationaryNonLinearHeatTransfer
    OrthotropicStationaryNonLinearHeatTransfer
    HeatStorage
    OrthotropicMicromorphicDamage
    UO2Shrinkage_RAPHAEL2008
    FailingIsotropicLinearElasticity
    OrthotropicSaintVenantKirchhoffElasticity)

  # TwoLayerCube executable
  add_test_executable(TwoLayerCube)

  # Helper function for TwoLayerCube cases
  function(add_two_layer_cube_test ncase solver)
    set(test "TwoLayerCube-${ncase}-${solver}")
    set(TwoLayerCubeParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--test-case"
        "${ncase}"
        "--linearsolver"
        "${solver}")
    if(solver STREQUAL "CGSolver" OR solver STREQUAL "GMRESSolver")
      add_test(NAME ${test} COMMAND ./TwoLayerCube ${TwoLayerCubeParameters}
                                    "--no-parallel")
      set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    endif()
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./TwoLayerCube
                                     ${TwoLayerCubeParameters} "--parallel"
                                     "--refinement" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # Generate all TwoLayerCube cases. The hypre solvers and MUMPS are only
  # tested in parallel.
  set(TwoLayerCubeSolvers CGSolver GMRESSolver HypreFGMRES HyprePCG HypreGMRES)
  if(MFEM_USE_MUMPS)
    list(APPEND TwoLayerCubeSolvers MUMPSSolver)
  endif(MFEM_USE_MUMPS)
  foreach(ncase RANGE 0 5)
    foreach(solver IN LISTS TwoLayerCubeSolvers)
      add_two_layer_cube_test(${ncase} ${solver})
    endforeach()
  endforeach()

  # UniaxialTensileTest executable
  add_test_executable(UniaxialTensileTest)

  # Helper function for UniaxialTensileTest cases
  function(add_uniaxial_tensile_test behaviour internal_state_variable)
    set(test "UniaxialTensileTest-${behaviour}")
    set(UniaxialTensileTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        "--internal-state-variable"
        "${internal_state_variable}")
    add_test(NAME ${test}
             COMMAND ./UniaxialTensileTest ${UniaxialTensileTestParameters}
                     "--no-parallel")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./UniaxialTensileTest
                       ${UniaxialTensileTestParameters} "--parallel")
    endif(MFEM_USE_MPI)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
  endfunction()

  # UniaxialTensileTest cases
  add_uniaxial_tensile_test(OrthotropicElasticity EquivalentStrain)
  add_uniaxial_tensile_test(IsotropicLinearHardeningPlasticity
                            EquivalentPlasticStrain)
  add_uniaxial_tensile_test(Mazars Damage)
  add_uniaxial_tensile_test(SaintVenantKirchhoffElasticity EquivalentStrain)

  # ImposedPressureTest executable
  add_test_executable(ImposedPressureTest)

  # Helper function for ImposedPressureTest cases
  function(add_imposed_pressure_test behaviour)
    set(test "ImposedPressureTest-${behaviour}")
    set(ImposedPressureTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref")
    add_test(NAME ${test}
             COMMAND ./ImposedPressureTest ${ImposedPressureTestParameters}
                     "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./ImposedPressureTest
                       ${ImposedPressureTestParameters} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_imposed_pressure_test(IsotropicLinearElasticity)

  # ImposedPressureTest2 executable
  add_test_executable(ImposedPressureTest2)

  # Helper function for ImposedPressureTest2 cases
  function(add_imposed_pressure_test2 behaviour)
    set(test "ImposedPressureTest2-${behaviour}")
    set(ImposedPressureTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref")
    add_test(NAME ${test}
             COMMAND ./ImposedPressureTest2 ${ImposedPressureTestParameters}
                     "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./ImposedPressureTest2
                       ${ImposedPressureTestParameters} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_imposed_pressure_test2(IsotropicLinearElasticity)

  # FiniteElementSpaceOrderingTest executable
  add_test_executable(FiniteElementSpaceOrderingTest)
  set(FiniteElementSpaceOrderingTestParameters
      "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
      "--library" "$<TARGET_FILE:BehaviourTest>")
  add_test(NAME FiniteElementSpaceOrderingTest
           COMMAND ./FiniteElementSpaceOrderingTest
                   ${FiniteElementSpaceOrderingTestParameters} "--parallel" "0")
  set_property(TEST FiniteElementSpaceOrderingTest PROPERTY DEPENDS BehaviourTest)
  if(MFEM_USE_MPI)
    add_test(NAME FiniteElementSpaceOrderingTestP
             COMMAND mpirun -n 2 ./FiniteElementSpaceOrderingTest
                     ${FiniteElementSpaceOrderingTestParameters} "--parallel" "1")
    set_property(TEST FiniteElementSpaceOrderingTestP
                 PROPERTY DEPENDS BehaviourTest)
  endif(MFEM_USE_MPI)

  # PartialQuadratureSpaceTest executable
  add_tfel_test_executable(PartialQuadratureSpaceTest)

  # Helper function for PartialQuadratureSpaceTest cases
  function(add_partial_quadrature_space_test behaviour internal_state_variable)
    set(test "PartialQuadratureSpaceTest-${behaviour}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        "--internal-state-variable"
        "${internal_state_variable}")
    add_test(NAME ${test} COMMAND ./PartialQuadratureSpaceTest ${params}
                                  "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PartialQuadratureSpaceTest
                                     ${params} "--parallel" 1)
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_partial_quadrature_space_test(OrthotropicElasticity EquivalentStrain)

  # PartialQuadratureSpaceTest2 executable
  add_tfel_test_executable(PartialQuadratureSpaceTest2)

  # Helper function for PartialQuadratureSpaceTest2 cases
  function(add_partial_quadrature_space_test2)
    set(test "PartialQuadratureSpaceTest2")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh")
    add_test(NAME ${test}
             COMMAND ./PartialQuadratureSpaceTest2 ${params} "--parallel" 0)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./PartialQuadratureSpaceTest2 ${params} "--parallel" 1)
    endif(MFEM_USE_MPI)
  endfunction()

  add_partial_quadrature_space_test2()

  # NonLinearEvolutionProblemPostProcessingsTest executable
  add_tfel_test_executable(NonLinearEvolutionProblemPostProcessingsTest)

  # Helper function for NonLinearEvolutionProblemPostProcessingsTest cases
  function(add_nonlinear_evolution_problem_post_processings_test)
    set(test "NonLinearEvolutionProblemPostProcessingsTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh")
    add_test(NAME ${test}
             COMMAND ./NonLinearEvolutionProblemPostProcessingsTest ${params} "--parallel" 0)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./NonLinearEvolutionProblemPostProcessingsTest ${params} "--parallel" 1)
    endif(MFEM_USE_MPI)
  endfunction()

  add_nonlinear_evolution_problem_post_processings_test()

  # UseMultiMaterialNonLinearIntegratorTest executable
  add_tfel_test_executable(UseMultiMaterialNonLinearIntegratorTest)

  # Helper function for UseMultiMaterialNonLinearIntegratorTest cases
  function(add_use_multi_material_non_linear_integrator_test)
    set(test "UseMultiMaterialNonLinearIntegratorTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh")
    add_test(NAME ${test}
             COMMAND ./UseMultiMaterialNonLinearIntegratorTest ${params} "--parallel" 0)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./UseMultiMaterialNonLinearIntegratorTest ${params} "--parallel" 1)
    endif(MFEM_USE_MPI)
  endfunction()

  add_use_multi_material_non_linear_integrator_test()

  # AddBehaviourIntegratorTest executable
  add_tfel_test_executable(AddBehaviourIntegratorTest)

  # Helper function for AddBehaviourIntegratorTest cases
  function(add_add_behaviour_integrator_test)
    set(test "AddBehaviourIntegratorTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh")
    add_test(NAME ${test}
             COMMAND ./AddBehaviourIntegratorTest ${params} "--parallel" 0)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./AddBehaviourIntegratorTest ${params} "--parallel" 1)
    endif(MFEM_USE_MPI)
  endfunction()

  add_add_behaviour_integrator_test()

  # MaterialQuantityProviderSearchTest executable
  add_tfel_test_executable(MaterialQuantityProviderSearchTest)

  # Helper function for MaterialQuantityProviderSearchTest cases
  function(add_material_quantity_provider_search_test)
    set(test "MaterialQuantityProviderSearchTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
               "--library" "$<TARGET_FILE:BehaviourTest>")
    add_test(NAME ${test}
             COMMAND ./MaterialQuantityProviderSearchTest ${params} "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./MaterialQuantityProviderSearchTest ${params} "--parallel" 1)
      set_property(TEST Parallel${test} PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_material_quantity_provider_search_test()

  # BehaviourIntegratorBaseTest executable
  add_tfel_test_executable(BehaviourIntegratorBaseTest)

  # Helper function for BehaviourIntegratorBaseTest cases
  function(add_behaviour_integrator_base_test)
    set(test "BehaviourIntegratorBaseTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
               "--library" "$<TARGET_FILE:BehaviourTest>")
    add_test(NAME ${test}
             COMMAND ./BehaviourIntegratorBaseTest ${params} "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./BehaviourIntegratorBaseTest ${params} "--parallel" 1)
      set_property(TEST Parallel${test} PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_behaviour_integrator_base_test()

  # BehaviourIntegratorsSelectionTest executable
  add_test_executable(BehaviourIntegratorsSelectionTest)

  # Helper function for BehaviourIntegratorsSelectionTest
  function(add_behaviour_integrators_selection_test)
    set(test "BehaviourIntegratorsSelectionTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh"
               "--library" "$<TARGET_FILE:BehaviourTest>"
               "--behaviour" "SaintVenantKirchhoffElasticityV2")
    add_test(NAME ${test} COMMAND ./BehaviourIntegratorsSelectionTest ${params}
                                  "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./BehaviourIntegratorsSelectionTest ${params}
                       "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_behaviour_integrators_selection_test()

  # StateManagerTest executable
  add_tfel_test_executable(StateManagerTest)

  # Helper function for StateManager cases
  function(add_state_manager_test)
    set(test "StateManagerTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh")
    add_test(NAME ${test}
             COMMAND ./StateManagerTest ${params} "--parallel" 0)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./StateManagerTest ${params} "--parallel" 1)
    endif(MFEM_USE_MPI)
  endfunction()

  add_state_manager_test()

  # MeshDiscretizationTest executable
  add_tfel_test_executable(MeshDiscretizationTest)

  # Helper function for MeshDiscretizationTest cases
  function(add_mesh_discretization_test)
    set(test "MeshDiscretizationTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh")
    add_test(NAME ${test} COMMAND ./MeshDiscretizationTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./MeshDiscretizationTest ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_mesh_discretization_test()

  # FiniteElementDiscretizationTest executable
  add_tfel_test_executable(FiniteElementDiscretizationTest)

  # Helper function for FiniteElementDiscretizationTest cases
  function(add_finite_element_discretization_test)
    set(test "FiniteElementDiscretizationTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh")
    add_test(NAME ${test} COMMAND ./FiniteElementDiscretizationTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./FiniteElementDiscretizationTest ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_finite_element_discretization_test()

  # PartialQuadratureFunctionTest executable
  add_tfel_test_executable(PartialQuadratureFunctionTest)

  # Helper function for PartialQuadratureFunctionTest cases
  function(add_partial_quadrature_function_test behaviour
           internal_state_variable)
    set(test "PartialQuadratureFunctionTest-${behaviour}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        "--internal-state-variable"
        "${internal_state_variable}")
    add_test(NAME ${test} COMMAND ./PartialQuadratureFunctionTest ${params}
                                  "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./PartialQuadratureFunctionTest ${params}
                       "--parallel" 1)
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_partial_quadrature_function_test(OrthotropicElasticity EquivalentStrain)

  # StationaryNonLinearHeatTransferTest executable
  add_test_executable(StationaryNonLinearHeatTransferTest)

  # Helper function for StationaryNonLinearHeatTransferTest cases
  function(add_stationary_non_linear_heat_transfer_test nsolver)
    set(test "StationaryNonLinearHeatTransferTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library" "$<TARGET_FILE:BehaviourTest>"
        "--behaviour" "StationaryNonLinearHeatTransfer"
        "--linearsolver" "${nsolver}")
    add_test(NAME ${test} COMMAND ./StationaryNonLinearHeatTransferTest
                                  ${params} "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./StationaryNonLinearHeatTransferTest
                       ${params} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # StationaryNonLinearHeatTransferTest cases
  add_stationary_non_linear_heat_transfer_test(1)
  if(MFEM_USE_SUITESPARSE)
    add_stationary_non_linear_heat_transfer_test(2)
  endif()

  # PointWiseModelTest executable
  add_test_executable(PointWiseModelTest)

  # Helper function for PointWiseModelTest
  function(add_point_wise_model_test)
    set(test "PointWiseModelTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
    add_test(NAME ${test} COMMAND ./PointWiseModelTest ${params} "--parallel"
                                  "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PointWiseModelTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_point_wise_model_test()

  # PointWiseModelTest-2 executable
  add_test_executable(PointWiseModelTest-2)

  # Helper function for PointWiseModelTest-2
  function(add_point_wise_model_test_2)
    set(test "PointWiseModelTest-2")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
    add_test(NAME ${test} COMMAND ./PointWiseModelTest-2 ${params} "--parallel"
                                  "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PointWiseModelTest-2
                                     ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_point_wise_model_test_2()

  # L2ProjectionTest executable
  add_test_executable(L2ProjectionTest)

  # Helper function for L2ProjectionTest cases
  function(add_l2_projection_test nsolver)
    set(test "L2ProjectionTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver"
               "${nsolver}")
    add_test(NAME ${test} COMMAND ./L2ProjectionTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test} COMMAND mpirun -n 2 ./L2ProjectionTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_l2_projection_test(0)

  # PointsSetCurveTest executable
  add_tfel_test_executable(PointsSetCurveTest)

  # Helper function for PointsSetCurveTest cases
  function(add_points_set_curve_test nsolver)
    set(test "PointsSetCurveTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver"
               "${nsolver}")
    add_test(NAME ${test} COMMAND ./PointsSetCurveTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test} COMMAND mpirun -n 2 ./PointsSetCurveTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_points_set_curve_test(0)

  # GridFunctionTest executable
  add_test_executable(GridFunctionTest)

  # Helper function for GridFunctionTest
  function(add_grid_function_test)
    set(test "GridFunctionTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
    add_test(NAME ${test} COMMAND ./GridFunctionTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./GridFunctionTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_grid_function_test()

  # GridFunctionOnSubMeshTest executable
  add_tfel_test_executable(GridFunctionOnSubMeshTest)

  # Helper function for GridFunctionOnSubMeshTest
  function(add_grid_function_on_submesh_test)
    set(test "GridFunctionOnSubMeshTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh")
    add_test(NAME ${test} COMMAND ./GridFunctionOnSubMeshTest ${params}
                                  "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./GridFunctionOnSubMeshTest
                                     ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_grid_function_on_submesh_test()

  # ImplicitGradientRegularizationTest executable
  add_test_executable(ImplicitGradientRegularizationTest)

  # Helper function for ImplicitGradientRegularizationTest cases
  function(add_implicit_gradient_regularization_test nsolver)
    set(test "ImplicitGradientRegularizationTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver"
               "${nsolver}")
    add_test(NAME ${test} COMMAND ./ImplicitGradientRegularizationTest
                                  ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./ImplicitGradientRegularizationTest
                       ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_implicit_gradient_regularization_test(0)

  # MicromorphicDamage2DTest executable
  add_test_executable(MicromorphicDamage2DTest)

  # Helper function for MicromorphicDamage2DTest cases
  function(add_micromorphic_damage_2d_test nsolver)
    set(test "MicromorphicDamage2DTest-${nsolver}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "MiehePhaseFieldDamage"
        "--linearsolver"
        "${nsolver}")
    add_test(NAME ${test} COMMAND ./MicromorphicDamage2DTest ${params}
                                  "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage2DTest
                                     ${params} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # MicromorphicDamage2DTest2 executable
  add_test_executable(MicromorphicDamage2DTest2)

  # Helper function for MicromorphicDamage2DTest2 cases
  function(add_micromorphic_damage_2d_test2 nsolver)
    set(test "MicromorphicDamage2DTest2-${nsolver}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "MicromorphicDamageII"
        "--linearsolver"
        "${nsolver}")
    add_test(NAME ${test} COMMAND ./MicromorphicDamage2DTest2 ${params}
                                  "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage2DTest2
                                     ${params} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # MicromorphicDamage3DTest executable
  add_test_executable(MicromorphicDamage3DTest)

  # Helper function for MicromorphicDamage3DTest cases
  function(add_micromorphic_damage_3d_test nsolver)
    set(test "MicromorphicDamage3DTest-${nsolver}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/bar-with-defect-3D.msh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "MicromorphicDamageII"
        "--linearsolver"
        "${nsolver}")
    add_test(NAME ${test} COMMAND ./MicromorphicDamage3DTest ${params}
                                  "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage3DTest
                                     ${params} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # Micromorphic damage tests
  add_micromorphic_damage_2d_test2(0)
  add_micromorphic_damage_2d_test(0)
  add_micromorphic_damage_3d_test(0)

  # Helper function for PredictionTest cases
  function(add_prediction_test test)
    add_test_executable(${test})
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh" "--library"
               "$<TARGET_FILE:BehaviourTest>")
    add_test(NAME ${test} COMMAND "./${test}" ${params} "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 "./${test}" ${params}
                                     "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # PredictionTest cases
  add_prediction_test(PredictionTest)
  add_prediction_test(PredictionTest2)

  if(MGIS_HAVE_TFEL)
    add_tfel_test_executable(GeometryTest)
    add_test(NAME GeometryTest COMMAND ./GeometryTest)
  endif(MGIS_HAVE_TFEL)

  if(MFEM_USE_GSLIB AND MGIS_HAVE_TFEL)
    # GridFunctionInterpolatorTest executable
    add_tfel_test_executable(GridFunctionInterpolatorTest)
    
    # Helper function for GridFunctionInterpolatorTest
    function(add_grid_function_interpolator_test)
      set(test "GridFunctionInterpolatorTest")
      set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
      add_test(NAME ${test} COMMAND ./GridFunctionInterpolatorTest ${params} "--serial")
      if(MFEM_USE_MPI)
        add_test(NAME ${test}P
                 COMMAND mpirun -n 2 ./GridFunctionInterpolatorTest ${params} "--parallel")
      endif(MFEM_USE_MPI)
    endfunction()
    
    add_grid_function_interpolator_test()
  endif(MFEM_USE_GSLIB AND MGIS_HAVE_TFEL)

  # ParallelReadMode executable (MPI only)
  if(MFEM_USE_MPI)
    add_test_executable(ParallelReadMode)
    set_target_properties(ParallelReadMode PROPERTIES COMPILE_FLAGS
                                                      "-DDO_USE_MPI")
    # Commented out test - can be enabled if needed add_test(NAME
    # Restart_Parallel_Read_Mode_MPI_2 COMMAND mpirun -n 2 ./ParallelReadMode
    # --mesh ${CMAKE_CURRENT_SOURCE_DIR}/mesh-test. --library
    # "$<TARGET_FILE:BehaviourTest>" --read-mesh-mode Restart)
  endif()

  # PostProSubMesh executable (MPI only)
  if(MFEM_USE_MPI)
    add_test_executable(PostProSubMesh)
    add_test(
      NAME PostProSubMeshPar
      COMMAND
        mpirun -n 2 ./PostProSubMesh --mesh
        ${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh --library
        "$<TARGET_FILE:BehaviourTest>")
  endif()

  # BatchIntegrationTest executable: integration of the behaviours by
  # batches compared to the integration point by point, with the behaviours
  # generated by the generic interface and, if the configuration files of
  # TFEL are available, by the generic-parallel interface on the host
  # (backends openmp and stlpar)
  add_test_executable(BatchIntegrationTest)
  set(batch_integration_test_libraries)
  set(batch_integration_test_behaviours)
  foreach(_behaviour IsotropicLinearHardeningPlasticity
                     IsotropicLinearHardeningThermoPlasticity
                     IsotropicLinearElasticity
                     FailingIsotropicLinearElasticity
                     OrthotropicElasticity
                     OrthotropicSaintVenantKirchhoffElasticity
                     StationaryNonLinearHeatTransfer
                     OrthotropicStationaryNonLinearHeatTransfer
                     HeatStorage
                     MicromorphicDamageII
                     OrthotropicMicromorphicDamage)
    list(APPEND batch_integration_test_behaviours
         "${CMAKE_CURRENT_SOURCE_DIR}/${_behaviour}.mfront")
  endforeach()
  get_filename_component(_mfront_bin "${MFRONT}" DIRECTORY)
  # the behaviours of the cuda backend are only given to the tests on the GPU
  set(batch_integration_test_backends openmp:config-openmp-host-gcc.json
                                      stlpar:config-stlpar-parunseq-gcc.json)
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    list(APPEND batch_integration_test_backends cuda:config-cuda-nvcc.json)
  endif()
  set(batch_integration_test_cuda_library)
  foreach(backend_configuration ${batch_integration_test_backends})
    string(REPLACE ":" ";" _bc "${backend_configuration}")
    list(GET _bc 0 _backend)
    list(GET _bc 1 _configuration)
    set(_configuration_file
        "${_mfront_bin}/../share/tfel/mfront/configuration-files/${_configuration}")
    if(EXISTS "${_configuration_file}")
      set(_directory "${CMAKE_CURRENT_BINARY_DIR}/generic-parallel-${_backend}")
      set(_library
          "${_directory}/src/libGenericParallelBehaviour-${_backend}.so")
      file(MAKE_DIRECTORY "${_directory}")
      add_custom_command(
        OUTPUT "${_library}"
        COMMAND "${MFRONT}" --obuild --interface=generic-parallel
                "--configuration-file=${_configuration}"
                ${batch_integration_test_behaviours}
        WORKING_DIRECTORY "${_directory}"
        DEPENDS ${batch_integration_test_behaviours}
        COMMENT "MFront: generic-parallel behaviour (${_backend})")
      add_custom_target(GenericParallelBehaviour-${_backend}
                        DEPENDS "${_library}")
      add_dependencies(build-tests GenericParallelBehaviour-${_backend})
      if(_backend STREQUAL "cuda")
        set(batch_integration_test_cuda_library
            "--cuda-library" "${_library}")
      else()
        list(APPEND batch_integration_test_libraries
             "--${_backend}-library" "${_library}")
      endif()
    endif()
  endforeach()
  set(batch_integration_test_arguments
      --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
      --mesh-2d ${CMAKE_CURRENT_SOURCE_DIR}/bar.msh
      --square-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh
      --two-materials-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube-two-materials.mesh
      --interface-mesh
      ${CMAKE_CURRENT_SOURCE_DIR}/cube-two-materials-interface.mesh
      --library $<TARGET_FILE:BehaviourTest>)
  add_test(NAME BatchIntegrationTest
           COMMAND BatchIntegrationTest ${batch_integration_test_arguments}
                   ${batch_integration_test_libraries})
  # the variants of the longest tests on the debug device of MFEM, which
  # take up to one minute and a half in an optimised build, are only
  # registered in the optimised builds: they would take much longer in the
  # Debug and Coverage builds. Without device, the stiffness matrices and
  # their diagonals are computed by the version of the kernel for the host
  # (one iteration per element): in these builds, the kernel itself (one
  # thread per pair of nodes) is run by the short variants on the debug
  # device (GradientsKernelsTest-debug, FixedPointSolverTest-kernels-debug)
  # and by the variants on the GPU.
  set(long_debug_device_tests ON)
  if((CMAKE_BUILD_TYPE STREQUAL "Debug") OR
     (CMAKE_BUILD_TYPE STREQUAL "Coverage"))
    set(long_debug_device_tests OFF)
  endif()
  # computations by kernels on the debug device of MFEM and on the GPU; on
  # the debug device, the behaviours generated by the generic-parallel
  # interface on the host are integrated on the copies of the arrays on the
  # device, which are in the memory of the host
  if(long_debug_device_tests)
    add_test(NAME BatchIntegrationTest-debug
             COMMAND BatchIntegrationTest --device debug
                     ${batch_integration_test_arguments}
                     ${batch_integration_test_libraries})
  endif()
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME BatchIntegrationTest-cuda
             COMMAND BatchIntegrationTest --device cuda
                     ${batch_integration_test_arguments}
                     ${batch_integration_test_cuda_library})
  endif()
  # computations on two processes (the comparisons only), with MUMPS, one
  # thread per process
  if(MFEM_USE_MPI AND MFEM_USE_MUMPS)
    add_test(NAME BatchIntegrationTestP
             COMMAND mpirun -n 2 ./BatchIntegrationTest --parallel
                     ${batch_integration_test_arguments}
                     ${batch_integration_test_libraries})
    set_tests_properties(BatchIntegrationTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME BatchIntegrationTestP-cuda
               COMMAND mpirun -n 2 ./BatchIntegrationTest --parallel
                       --device cuda ${batch_integration_test_arguments}
                       ${batch_integration_test_cuda_library})
      set_tests_properties(BatchIntegrationTestP-cuda
                           PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
    endif()
  endif()

  # HeatTransferKernelsTest executable: stationary non linear heat transfer
  # by kernels compared to the integration point by point, on the host, on
  # the debug device of MFEM and on the GPU, sequential and on two processes
  add_test_executable(HeatTransferKernelsTest)
  set(heat_transfer_kernels_test_arguments
      --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
      --two-materials-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube-two-materials.mesh
      --square-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh
      --library $<TARGET_FILE:BehaviourTest>)
  # the sequential computations use UMFPack
  if(MFEM_USE_SUITESPARSE)
    add_test(NAME HeatTransferKernelsTest
             COMMAND HeatTransferKernelsTest
                     ${heat_transfer_kernels_test_arguments}
                     ${batch_integration_test_libraries})
    # the debug device checks the memory accesses of the kernels: the
    # behaviours generated by the generic interface are enough
    if(long_debug_device_tests)
      add_test(NAME HeatTransferKernelsTest-debug
               COMMAND HeatTransferKernelsTest --device debug
                       ${heat_transfer_kernels_test_arguments})
    endif()
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME HeatTransferKernelsTest-cuda
               COMMAND HeatTransferKernelsTest --device cuda
                       ${heat_transfer_kernels_test_arguments}
                       ${batch_integration_test_cuda_library})
      set_tests_properties(HeatTransferKernelsTest-cuda
                           PROPERTIES RESOURCE_LOCK gpu)
    endif()
  endif()
  if(MFEM_USE_MPI AND MFEM_USE_MUMPS)
    add_test(NAME HeatTransferKernelsTestP
             COMMAND mpirun -n 2 ./HeatTransferKernelsTest --parallel
                     ${heat_transfer_kernels_test_arguments})
    set_tests_properties(HeatTransferKernelsTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
  endif()

  # MicromorphicKernelsTest executable: micromorphic damage by kernels
  # compared to the integration point by point, on the host, on the debug
  # device of MFEM and on the GPU, sequential and on two processes
  add_test_executable(MicromorphicKernelsTest)
  set(micromorphic_kernels_test_arguments
      --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
      --square-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh
      --library $<TARGET_FILE:BehaviourTest>)
  # the sequential computations use UMFPack
  if(MFEM_USE_SUITESPARSE)
    add_test(NAME MicromorphicKernelsTest
             COMMAND MicromorphicKernelsTest
                     ${micromorphic_kernels_test_arguments}
                     ${batch_integration_test_libraries})
    # the debug device checks the memory accesses of the kernels: the
    # behaviours generated by the generic interface are enough
    if(long_debug_device_tests)
      add_test(NAME MicromorphicKernelsTest-debug
               COMMAND MicromorphicKernelsTest --device debug
                       ${micromorphic_kernels_test_arguments})
    endif()
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME MicromorphicKernelsTest-cuda
               COMMAND MicromorphicKernelsTest --device cuda
                       ${micromorphic_kernels_test_arguments}
                       ${batch_integration_test_cuda_library})
      set_tests_properties(MicromorphicKernelsTest-cuda
                           PROPERTIES RESOURCE_LOCK gpu)
    endif()
  endif()
  if(MFEM_USE_MPI AND MFEM_USE_MUMPS)
    add_test(NAME MicromorphicKernelsTestP
             COMMAND mpirun -n 2 ./MicromorphicKernelsTest --parallel
                     ${micromorphic_kernels_test_arguments})
    set_tests_properties(MicromorphicKernelsTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
  endif()

  # FixedPointSolverTest executable: fixed point solver (methods,
  # accelerations, linear solvers of the elastic operator) compared to the
  # Newton solver
  add_test_executable(FixedPointSolverTest)
  add_test(NAME FixedPointSolverTest
           COMMAND FixedPointSolverTest
                   --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                   --library $<TARGET_FILE:BehaviourTest>)
  if(MFEM_USE_MPI AND MFEM_USE_MUMPS)
    add_test(NAME FixedPointSolverTestP
             COMMAND mpirun -n 2 ./FixedPointSolverTest --parallel
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>)
    set_tests_properties(FixedPointSolverTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
  # integration by batches and computations by kernels, on the host and on
  # the debug device of MFEM
  set(fixed_point_solver_test_kernels
      MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch
      MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels
      MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels
      MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels)
  foreach(_device cpu debug)
    add_test(NAME FixedPointSolverTest-kernels-${_device}
             COMMAND FixedPointSolverTest --device ${_device}
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>)
    set_tests_properties(FixedPointSolverTest-kernels-${_device}
                         PROPERTIES ENVIRONMENT
                         "${fixed_point_solver_test_kernels}")
  endforeach()

  # MaterialArraysTest executable: arrays of the materials in the memory of
  # MFEM, on the host, on the debug device of MFEM and on the GPU
  if(NOT MFEM_USE_HIP STREQUAL "YES")
    add_test_executable(MaterialArraysTest)
    mfem_mgis_compile_as_device_code(MaterialArraysTest MaterialArraysTest.cxx)
    add_test(NAME MaterialArraysTest
             COMMAND MaterialArraysTest
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --multi-material-mesh
                     ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                     --mesh-2d ${CMAKE_CURRENT_SOURCE_DIR}/bar.msh
                     --library $<TARGET_FILE:BehaviourTest>)
    add_test(NAME MaterialArraysTest-debug
             COMMAND MaterialArraysTest --device debug
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --multi-material-mesh
                     ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                     --mesh-2d ${CMAKE_CURRENT_SOURCE_DIR}/bar.msh
                     --library $<TARGET_FILE:BehaviourTest>)
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME MaterialArraysTest-cuda
               COMMAND MaterialArraysTest --device cuda
                       --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                       --multi-material-mesh
                       ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                       --mesh-2d ${CMAKE_CURRENT_SOURCE_DIR}/bar.msh
                       --library $<TARGET_FILE:BehaviourTest>)
    endif()
  endif()

  # GradientsKernelsTest executable: gradients of the behaviours computed
  # by kernels, on the host, on the debug device of MFEM and on the GPU
  if(NOT MFEM_USE_HIP STREQUAL "YES")
    add_test_executable(GradientsKernelsTest)
    set(gradients_kernels_test_arguments
        --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
        --tetrahedral-mesh ${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh
        --two-materials-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
        --square-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh
        --library $<TARGET_FILE:BehaviourTest>)
    add_test(NAME GradientsKernelsTest
             COMMAND GradientsKernelsTest ${gradients_kernels_test_arguments})
    add_test(NAME GradientsKernelsTest-debug
             COMMAND GradientsKernelsTest --device debug
                     ${gradients_kernels_test_arguments})
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME GradientsKernelsTest-cuda
               COMMAND GradientsKernelsTest --device cuda
                       ${gradients_kernels_test_arguments})
    endif()
    # GlobalStiffnessMatrixTest executable: global matrix assembled from the
    # local matrices with a pattern and exchanges computed once, compared to
    # the product of hypre, on one, two and four processes
    if(MFEM_USE_MPI)
      add_test_executable(GlobalStiffnessMatrixTest)
      foreach(_np 1 2 4)
        add_test(NAME GlobalStiffnessMatrixTest-np${_np}
                 COMMAND mpirun -n ${_np} ./GlobalStiffnessMatrixTest
                         --mesh
                         ${CMAKE_CURRENT_SOURCE_DIR}/cube-two-materials.mesh)
        set_tests_properties(GlobalStiffnessMatrixTest-np${_np}
                             PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
      endforeach()
    endif()
    # IntegrationRuleOrderTest executable: order of the integration rules
    # of the behaviour integrators given by the user, point by point and by
    # kernels, on the host and on the debug device of MFEM (UMFPack)
    if(MFEM_USE_SUITESPARSE)
      add_test_executable(IntegrationRuleOrderTest)
      set(integration_rule_order_test_arguments
          --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
          --library $<TARGET_FILE:BehaviourTest>)
      add_test(NAME IntegrationRuleOrderTest
               COMMAND IntegrationRuleOrderTest
                       ${integration_rule_order_test_arguments})
      add_test(NAME IntegrationRuleOrderTest-debug
               COMMAND IntegrationRuleOrderTest --device debug
                       ${integration_rule_order_test_arguments})
    endif()
  endif()

  # HypreBoomerAMGTest executable
  if(MFEM_USE_MPI)
    add_test_executable(HypreBoomerAMGTest)
    add_test(NAME HypreBoomerAMGTestP COMMAND mpirun -n 2 ./HypreBoomerAMGTest)
    # HypreILUTest executable: options of the incomplete factorisation of
    # hypre kept after a new operator, ILU(0) when hypre computes on a GPU.
    # The levels of fill change the numbers of iterations more on one
    # process than on two (block Jacobi factorisations)
    add_test_executable(HypreILUTest)
    add_test(NAME HypreILUTest COMMAND HypreILUTest)
    add_test(NAME HypreILUTestP COMMAND mpirun -n 2 ./HypreILUTest)
  endif(MFEM_USE_MPI)

endif(MFEM_MGIS_HAVE_TFEL)

# tests of the linear solvers, without behaviours
# OperatorSmoothersTest executable: Jacobi and Chebyshev smoothers, which
# only use products by the operator and its diagonal, as preconditioners
# of the conjugate gradient
add_test_executable(OperatorSmoothersTest)
add_test(NAME OperatorSmoothersTest COMMAND OperatorSmoothersTest)
if(MFEM_USE_MPI)
  add_test(NAME OperatorSmoothersTestP
           COMMAND mpirun -n 2 ./OperatorSmoothersTest --parallel)
  set_tests_properties(OperatorSmoothersTestP
                       PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
endif()

# AdaptiveLinearSolverToleranceTest executable: relative tolerance of the
# iterative linear solvers adapted to the decrease of the residual
add_test_executable(AdaptiveLinearSolverToleranceTest)
add_test(NAME AdaptiveLinearSolverToleranceTest
         COMMAND AdaptiveLinearSolverToleranceTest)
if(MFEM_USE_MPI)
  add_test(NAME AdaptiveLinearSolverToleranceTestP
           COMMAND mpirun -n 2 ./AdaptiveLinearSolverToleranceTest)
  set_tests_properties(AdaptiveLinearSolverToleranceTestP
                       PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
endif()

# HypreKrylovAbsoluteToleranceTest executable
if(MFEM_USE_MPI)
  add_test_executable(HypreKrylovAbsoluteToleranceTest)
  add_test(NAME HypreKrylovAbsoluteToleranceTestP
           COMMAND mpirun -n 2 ./HypreKrylovAbsoluteToleranceTest)
endif(MFEM_USE_MPI)

# AmgXSolverTest executable: solvers of AmgX built by the factory of linear
# solvers. AmgX computes on the GPU whatever the device of MFEM: the tests
# require one, and are not run at the same time as the other tests on the
# GPU
if(MFEM_MGIS_HAVE_AMGX AND enable-gpu-tests)
  add_test_executable(AmgXSolverTest)
  add_test(NAME AmgXSolverTest-gpu COMMAND AmgXSolverTest)
  add_test(NAME AmgXSolverTest-cuda COMMAND AmgXSolverTest --device cuda)
  # no GPU visible: the solvers of AmgX are refused
  add_test(NAME AmgXSolverTest-no-gpu COMMAND AmgXSolverTest --no-gpu)
  set_tests_properties(AmgXSolverTest-no-gpu
                       PROPERTIES ENVIRONMENT "CUDA_VISIBLE_DEVICES=")
  set(amgx_tests AmgXSolverTest-gpu AmgXSolverTest-cuda)
  if(MFEM_USE_MPI)
    add_test(NAME AmgXSolverTestP-gpu
             COMMAND mpirun -n 2 ./AmgXSolverTest --parallel)
    add_test(NAME AmgXSolverTestP-cuda
             COMMAND mpirun -n 2 ./AmgXSolverTest --parallel --device cuda)
    list(APPEND amgx_tests AmgXSolverTestP-gpu AmgXSolverTestP-cuda)
  endif(MFEM_USE_MPI)
  # a deadlock fails quickly
  set_tests_properties(${amgx_tests} PROPERTIES RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# HypreDevicePoolSizeTest executable: size of the pool of memory of hypre on
# the GPU, given by the option of mfem-mgis or by its environment variable
# (no memory allocated on the GPU)
if(MFEM_USE_MPI)
  add_test_executable(HypreDevicePoolSizeTest)
  add_test(NAME HypreDevicePoolSizeTest
           COMMAND HypreDevicePoolSizeTest --hypre-device-pool-size 0.25
                   --expected 0.25)
  add_test(NAME HypreDevicePoolSizeTest-short
           COMMAND HypreDevicePoolSizeTest -hdps 0.125 --expected 0.125)
  add_test(NAME HypreDevicePoolSizeTest-environment
           COMMAND HypreDevicePoolSizeTest --expected 0.5)
  add_test(NAME HypreDevicePoolSizeTest-option-and-environment
           COMMAND HypreDevicePoolSizeTest -hdps 0.25 --expected 0.25)
  set_tests_properties(HypreDevicePoolSizeTest-environment
                       HypreDevicePoolSizeTest-option-and-environment
                       PROPERTIES ENVIRONMENT
                       MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=0.5)
  add_test(NAME HypreDevicePoolSizeTest-negative
           COMMAND HypreDevicePoolSizeTest -hdps -1 --expected error)
  add_test(NAME HypreDevicePoolSizeTest-invalid
           COMMAND HypreDevicePoolSizeTest -hdps 1x --expected error)
  add_test(NAME HypreDevicePoolSizeTest-missing
           COMMAND HypreDevicePoolSizeTest --expected error -hdps)
  # less than a byte, more bytes than the range of std::size_t
  add_test(NAME HypreDevicePoolSizeTest-too-small
           COMMAND HypreDevicePoolSizeTest -hdps 1e-10 --expected error)
  add_test(NAME HypreDevicePoolSizeTest-too-large
           COMMAND HypreDevicePoolSizeTest -hdps 1e30 --expected error)
  # actual size of the pool created on the GPU by a first allocation
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME HypreDevicePoolSizeTest-gpu
             COMMAND HypreDevicePoolSizeTest --device cuda -hdps 0.25
                     --expected 0.25 --allocate)
    set_tests_properties(HypreDevicePoolSizeTest-gpu
                         PROPERTIES RESOURCE_LOCK gpu)
  endif()
endif(MFEM_USE_MPI)

# DeviceForallTest executable: a kernel of MFEM, compiled as CUDA code
# when MFEM is compiled with CUDA (not supported yet with HIP). The test on
# the GPU requires one.
if(NOT MFEM_USE_HIP STREQUAL "YES")
  add_test_executable(DeviceForallTest)
  mfem_mgis_compile_as_device_code(DeviceForallTest DeviceForallTest.cxx)
  add_test(NAME DeviceForallTest-cpu
           COMMAND DeviceForallTest --expected cpu)
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME DeviceForallTest-cuda
             COMMAND DeviceForallTest --device cuda --expected cuda)
    # the processes of a node use its GPUs in turn
    if(MFEM_USE_MPI)
      add_test(NAME DeviceForallTestP-cuda
               COMMAND mpirun -n 2 ./DeviceForallTest --device cuda
                       --expected cuda)
    endif(MFEM_USE_MPI)
  endif()
endif()

# ElementsGeometryTest executable: jacobians of the element transformations
# used by the kernels of the batch path, computed by a kernel for the meshes
# without nodes (not compiled with HIP), on the host, on the debug device,
# on the GPU and on two processes
if(NOT MFEM_USE_HIP STREQUAL "YES")
  add_test_executable(ElementsGeometryTest)
  add_test(NAME ElementsGeometryTest COMMAND ElementsGeometryTest)
  add_test(NAME ElementsGeometryTest-debug
           COMMAND ElementsGeometryTest --device debug)
  if(MFEM_USE_MPI)
    add_test(NAME ElementsGeometryTestP
             COMMAND mpirun -n 2 ./ElementsGeometryTest)
    set_tests_properties(ElementsGeometryTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME ElementsGeometryTest-cuda
             COMMAND ElementsGeometryTest --device cuda)
    set_tests_properties(ElementsGeometryTest-cuda
                         PROPERTIES RESOURCE_LOCK gpu)
  endif()
endif()

# MFEMFieldOperationsTest executable: copies and fills of views by kernels
# of MFEM, on the host, on the debug device and on the GPU
if(NOT MFEM_USE_HIP STREQUAL "YES")
  add_test_executable(MFEMFieldOperationsTest)
  mfem_mgis_compile_as_device_code(MFEMFieldOperationsTest
                                   MFEMFieldOperationsTest.cxx)
  add_test(NAME MFEMFieldOperationsTest COMMAND MFEMFieldOperationsTest)
  add_test(NAME MFEMFieldOperationsTest-debug
           COMMAND MFEMFieldOperationsTest --device debug)
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME MFEMFieldOperationsTest-cuda
             COMMAND MFEMFieldOperationsTest --device cuda)
  endif()
endif()

# DeviceOptionTest executable
add_test_executable(DeviceOptionTest)
add_test(NAME DeviceOptionTest-none
         COMMAND DeviceOptionTest --expected none)
add_test(NAME DeviceOptionTest-debug
         COMMAND DeviceOptionTest --device debug --expected debug)
add_test(NAME DeviceOptionTest-debug-short
         COMMAND DeviceOptionTest -dev debug --expected debug)

# Install the example 1, based on UniaxialTensileTest

set(EX1_DIR "share/mfem-mgis/examples/ex1")
install(FILES UniaxialTensileTest.cxx IsotropicLinearHardeningPlasticity.mfront
              cube.mesh references/IsotropicLinearHardeningPlasticity.ref
        DESTINATION ${EX1_DIR})

install(
  FILES ex1CMakeLists.txt
  DESTINATION ${EX1_DIR}
  RENAME CMakeLists.txt)

set(EX3_DIR "share/mfem-mgis/examples/ex3")
install(FILES TwoLayerCube.cxx IsotropicLinearElasticity.mfront
              cube_2mat_per.mesh Box.med Box.per
        DESTINATION ${EX3_DIR})

install(
  FILES ex3CMakeLists.txt
  DESTINATION ${EX3_DIR}
  RENAME CMakeLists.txt)

# each process on the GPU may reserve a large memory pool (4 GiB for hypre
# with Umpire): the tests on the GPU are not run at the same time
get_property(_tests DIRECTORY PROPERTY TESTS)
foreach(_test ${_tests})
  if(_test MATCHES "-cuda$")
    set_tests_properties(${_test} PROPERTIES RESOURCE_LOCK gpu)
  endif()
endforeach()
