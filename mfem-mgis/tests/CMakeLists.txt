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
add_tfel_test_executable(SubSteppingMessageTest)
add_test(NAME SubSteppingMessageTest COMMAND ./SubSteppingMessageTest
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
    PlaneStressSaintVenantKirchhoffElasticity
    MiehePhaseFieldDamage
    MicromorphicDamageI_SpectralSplit
    MicromorphicDamageII
    StationaryNonLinearHeatTransfer
    OrthotropicStationaryNonLinearHeatTransfer
    HeatStorage
    OrthotropicMicromorphicDamage
    UO2Shrinkage_RAPHAEL2008
    FailingIsotropicLinearElasticity
    OrthotropicSaintVenantKirchhoffElasticity
    LogarithmicStrainElasticity)

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
  if(MFEM_MGIS_HAVE_MUMPS)
    list(APPEND TwoLayerCubeSolvers MUMPSSolver)
  endif(MFEM_MGIS_HAVE_MUMPS)
  foreach(ncase RANGE 0 5)
    foreach(solver IN LISTS TwoLayerCubeSolvers)
      add_two_layer_cube_test(${ncase} ${solver})
    endforeach()
  endforeach()

  # linear solver override, given by the variable
  # MFEM_MGIS_LINEAR_SOLVER_OVERRIDE (sequential test) or by the option
  # --linear-solver-override (parallel test, with a preconditioner whose
  # options are empty braces), read by the utilities of TFEL: the GMRES solver
  # requested by the program is replaced by a conjugate gradient, the
  # replacement being reported (see ExecutionSettingsTest for the refusals)
  if(MFEM_MGIS_HAVE_TFEL_UTILITIES)
    set(_override_options "AbsoluteTolerance: 1e-12, RelativeTolerance: 1e-12, MaximumNumberOfIterations: 5000")
    set(_override_test_options
        --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
        --library $<TARGET_FILE:BehaviourTest> --test-case 1
        --linearsolver GMRESSolver)
    set(_override_regex
        "linear solver override 'CGSolver' replaces the linear solver 'GMRESSolver'.*Error is lower than threshold")
    set(_override_failures
        "Error is greater than threshold;Process received signal;terminate called")
    add_test(NAME TwoLayerCube-1-LinearSolverOverride
             COMMAND ./TwoLayerCube ${_override_test_options} --no-parallel)
    set_tests_properties(TwoLayerCube-1-LinearSolverOverride PROPERTIES
      DEPENDS BehaviourTest
      ENVIRONMENT "MFEM_MGIS_LINEAR_SOLVER_OVERRIDE={Name: \"CGSolver\", Options: {${_override_options}}}"
      PASS_REGULAR_EXPRESSION "${_override_regex}"
      FAIL_REGULAR_EXPRESSION "${_override_failures}")
    if(MFEM_USE_MPI)
      add_test(NAME TwoLayerCube-1-LinearSolverOverrideP
               COMMAND mpirun -n 2 ./TwoLayerCube ${_override_test_options}
                       --parallel --linear-solver-override
                       "{Name: \"CGSolver\", Options: {${_override_options}, Preconditioner: {Name: \"HypreDiagScale\", Options: {}}}}")
      set_tests_properties(TwoLayerCube-1-LinearSolverOverrideP PROPERTIES
        DEPENDS BehaviourTest
        PASS_REGULAR_EXPRESSION "${_override_regex}"
        FAIL_REGULAR_EXPRESSION "${_override_failures}")
      # preconditioner override alone: the GMRES solver of the program keeps
      # its options and receives the preconditioner
      add_test(NAME TwoLayerCube-1-PreconditionerOverrideP
               COMMAND mpirun -n 2 ./TwoLayerCube ${_override_test_options}
                       --parallel --preconditioner-override HypreDiagScale)
      set_tests_properties(TwoLayerCube-1-PreconditionerOverrideP PROPERTIES
        DEPENDS BehaviourTest
        PASS_REGULAR_EXPRESSION
        "preconditioner override 'HypreDiagScale' is given to the linear solver 'GMRESSolver'.*Error is lower than threshold"
        FAIL_REGULAR_EXPRESSION "${_override_failures}")
      # the preconditioner HypreDiagScale of the conjugate gradient of hypre
      # requested by the program is replaced by BoomerAMG
      add_test(NAME TwoLayerCube-1-PreconditionerOverrideReplacementP
               COMMAND mpirun -n 2 ./TwoLayerCube
                       --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                       --library $<TARGET_FILE:BehaviourTest> --test-case 1
                       --linearsolver HyprePCG --parallel
                       --preconditioner-override
                       "{Name: \"HypreBoomerAMG\", Options: {VerbosityLevel: 0}}")
      set_tests_properties(TwoLayerCube-1-PreconditionerOverrideReplacementP
        PROPERTIES
        DEPENDS BehaviourTest
        PASS_REGULAR_EXPRESSION
        "preconditioner override 'HypreBoomerAMG' is given to the linear solver 'HyprePCG'.*Error is lower than threshold"
        FAIL_REGULAR_EXPRESSION "${_override_failures}")
      # None, given to a direct solver which accepts no preconditioner
      if(MFEM_MGIS_HAVE_MUMPS)
        add_test(NAME TwoLayerCube-1-PreconditionerOverrideNoneP
                 COMMAND mpirun -n 2 ./TwoLayerCube
                         --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                         --library $<TARGET_FILE:BehaviourTest> --test-case 1
                         --linearsolver MUMPSSolver --parallel
                         --preconditioner-override None)
        set_tests_properties(TwoLayerCube-1-PreconditionerOverrideNoneP
          PROPERTIES
          DEPENDS BehaviourTest
          PASS_REGULAR_EXPRESSION
          "preconditioner override 'None' is given to the linear solver 'MUMPSSolver'.*Error is lower than threshold"
          FAIL_REGULAR_EXPRESSION "${_override_failures}")
      endif(MFEM_MGIS_HAVE_MUMPS)
    endif(MFEM_USE_MPI)
    # both overrides, given by the variables: the preconditioner override is
    # given to the linear solver override
    add_test(NAME TwoLayerCube-1-LinearSolverAndPreconditionerOverrides
             COMMAND ./TwoLayerCube ${_override_test_options} --no-parallel)
    set_tests_properties(TwoLayerCube-1-LinearSolverAndPreconditionerOverrides
      PROPERTIES
      DEPENDS BehaviourTest
      ENVIRONMENT "MFEM_MGIS_LINEAR_SOLVER_OVERRIDE={Name: \"CGSolver\", Options: {${_override_options}}};MFEM_MGIS_PRECONDITIONER_OVERRIDE=DSmoother"
      PASS_REGULAR_EXPRESSION
      "linear solver override 'CGSolver' replaces the linear solver 'GMRESSolver'.*preconditioner override 'DSmoother' is given to the linear solver 'CGSolver'.*Error is lower than threshold"
      FAIL_REGULAR_EXPRESSION "${_override_failures}")
    # non linear solver override, given by the option
    # --non-linear-solver-override (sequential test) or by the variable
    # MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE (parallel test): the Newton solver
    # of the program is replaced by the modified Newton method of the fixed
    # point solver, whose verbosity level, given by the override (-1, no
    # output), takes precedence over the one set by the program (0, which
    # prints the iterations)
    set(_nl_override
        "{Name: \"FixedPoint\", Options: {Method: \"ModifiedNewton\", VerbosityLevel: -1}}")
    set(_nl_override_regex
        "non linear solver override 'FixedPoint' replaces the non linear solver 'Newton'.*Error is lower than threshold")
    set(_nl_override_failures "${_override_failures};Fixed point iteration")
    add_test(NAME TwoLayerCube-1-NonLinearSolverOverride
             COMMAND ./TwoLayerCube ${_override_test_options} --no-parallel
                     --non-linear-solver-override "${_nl_override}")
    set_tests_properties(TwoLayerCube-1-NonLinearSolverOverride PROPERTIES
      DEPENDS BehaviourTest
      PASS_REGULAR_EXPRESSION "${_nl_override_regex}"
      FAIL_REGULAR_EXPRESSION "${_nl_override_failures}")
    if(MFEM_USE_MPI)
      add_test(NAME TwoLayerCube-1-NonLinearSolverOverrideP
               COMMAND mpirun -n 2 ./TwoLayerCube ${_override_test_options}
                       --parallel)
      set_tests_properties(TwoLayerCube-1-NonLinearSolverOverrideP PROPERTIES
        DEPENDS BehaviourTest
        ENVIRONMENT "MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE=${_nl_override}"
        PASS_REGULAR_EXPRESSION "${_nl_override_regex}"
        FAIL_REGULAR_EXPRESSION "${_nl_override_failures}")
    endif(MFEM_USE_MPI)
  endif(MFEM_MGIS_HAVE_TFEL_UTILITIES)

  # the value MatrixFree of MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION chooses
  # an algorithm, not a placement: it is kept with --execution, without the
  # warning of the ignored placement variables. The conjugate gradient
  # without preconditioner converges with the tangent operator applied
  # without assembly, while UMFPack, which requires an assembled matrix,
  # fails, which shows that the operator is not assembled
  add_test(NAME TwoLayerCube-1-MatrixFreeExecution
           COMMAND ./TwoLayerCube
                   --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                   --library $<TARGET_FILE:BehaviourTest> --test-case 1
                   --no-parallel --execution cpu)
  set_tests_properties(TwoLayerCube-1-MatrixFreeExecution PROPERTIES
    DEPENDS BehaviourTest
    ENVIRONMENT "MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree"
    PASS_REGULAR_EXPRESSION "Error is lower than threshold"
    FAIL_REGULAR_EXPRESSION
    "is ignored;Error is greater than threshold;Process received signal;terminate called")
  if(MFEM_MGIS_HAVE_SUITESPARSE)
    add_test(NAME TwoLayerCube-1-MatrixFreeExecution-UMFPack
             COMMAND ./TwoLayerCube
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
                     --library $<TARGET_FILE:BehaviourTest> --test-case 1
                     --no-parallel --execution cpu
                     --linearsolver UMFPackSolver)
    set_tests_properties(TwoLayerCube-1-MatrixFreeExecution-UMFPack
      PROPERTIES
      DEPENDS BehaviourTest
      ENVIRONMENT "MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=MatrixFree"
      PASS_REGULAR_EXPRESSION "not a SparseMatrix")
  endif(MFEM_MGIS_HAVE_SUITESPARSE)

  # UniaxialTensileTest executable
  add_test_executable(UniaxialTensileTest)

  # Helper function for UniaxialTensileTest cases. The optional arguments are
  # a suffix of the name of the test and additional options of the test, the
  # reference file being the one of the behaviour. The internal state
  # variable may be empty (none extracted): it is given apart from the list
  # of the parameters, which would drop an empty value.
  function(add_uniaxial_tensile_test behaviour internal_state_variable)
    set(test "UniaxialTensileTest-${behaviour}")
    set(options)
    if(ARGC GREATER 2)
      set(test "${test}-${ARGV2}")
      list(SUBLIST ARGN 1 -1 options)
    endif()
    set(UniaxialTensileTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        ${options})
    add_test(NAME ${test}
             COMMAND ./UniaxialTensileTest ${UniaxialTensileTestParameters}
                     "--internal-state-variable" "${internal_state_variable}"
                     "--no-parallel")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./UniaxialTensileTest
                       ${UniaxialTensileTestParameters}
                       "--internal-state-variable" "${internal_state_variable}"
                       "--parallel")
    endif(MFEM_USE_MPI)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
  endfunction()

  # UniaxialTensileTest cases
  add_uniaxial_tensile_test(OrthotropicElasticity EquivalentStrain)
  add_uniaxial_tensile_test(IsotropicLinearHardeningPlasticity
                            EquivalentPlasticStrain)
  add_uniaxial_tensile_test(Mazars Damage)
  add_uniaxial_tensile_test(SaintVenantKirchhoffElasticity EquivalentStrain)
  add_uniaxial_tensile_test(LogarithmicStrainElasticity EquivalentStrain)
  # finite strain behaviours given the displacement gradient H = F - I
  # (parameter DisplacementGradient of the behaviour integrators), compared
  # to the same reference values as with the deformation gradient: a
  # behaviour written with the deformation gradient (which receives
  # F = I + H) and a behaviour in the Hencky strain framework (whose strain
  # is computed from H)
  add_uniaxial_tensile_test(SaintVenantKirchhoffElasticity EquivalentStrain
                            DisplacementGradient --displacement-gradient)
  add_uniaxial_tensile_test(LogarithmicStrainElasticity EquivalentStrain
                            DisplacementGradient --displacement-gradient)
  # orthotropic finite strain behaviour (rotated by the test): the gradients
  # are rotated in the material frame, H as F. In double only: the stress
  # extracted is a lateral one in the material frame, zero within the
  # convergence criterion of the Newton solver, about 5e4 Pa in float
  if(NOT MFEM_MGIS_REAL_TYPE STREQUAL "float")
    add_uniaxial_tensile_test(OrthotropicSaintVenantKirchhoffElasticity "")
    add_uniaxial_tensile_test(OrthotropicSaintVenantKirchhoffElasticity ""
                              DisplacementGradient --displacement-gradient)
  endif()
  # convergence of the Newton solver at the floor of the residual, in float
  # only: with an unreachable relative tolerance, a stagnation of the
  # residual is accepted as convergence (option StagnationRelativeTolerance
  # of the Newton solver, whose message is required), and the resolution
  # fails when this option is disabled. The relative tolerance is bounded
  # below by 100 times the machine epsilon: the floor of the residual is
  # reached with the finite strain behaviour (about 3e-4 to 8e-4), not with
  # the small strain one. In double, the residual of these tests keeps
  # decreasing slowly (factor 0.5 per iteration at 5e-19 of the initial
  # residual): the floor is not reached in ten iterations
  if(MFEM_MGIS_REAL_TYPE STREQUAL "float")
    set(_svk_options
        --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
        --library $<TARGET_FILE:BehaviourTest>
        --behaviour SaintVenantKirchhoffElasticity
        --internal-state-variable EquivalentStrain
        --relative-tolerance 1e-30 --no-parallel --no-post-processing)
    add_test(NAME UniaxialTensileTest-Stagnation-SaintVenantKirchhoffElasticity
             COMMAND ./UniaxialTensileTest ${_svk_options}
                     --reference-file
                     ${CMAKE_CURRENT_SOURCE_DIR}/references/SaintVenantKirchhoffElasticity.ref)
    set_tests_properties(
      UniaxialTensileTest-Stagnation-SaintVenantKirchhoffElasticity
      PROPERTIES PASS_REGULAR_EXPRESSION "accepted as convergence"
                 FAIL_REGULAR_EXPRESSION "the resolution failed;test failed")
    add_test(NAME UniaxialTensileTest-NoStagnation-SaintVenantKirchhoffElasticity
             COMMAND ./UniaxialTensileTest ${_svk_options}
                     --stagnation-relative-tolerance 0)
    set_tests_properties(
      UniaxialTensileTest-NoStagnation-SaintVenantKirchhoffElasticity
      PROPERTIES PASS_REGULAR_EXPRESSION "the resolution failed")
    # the displacement gradient keeps the precision of small strains: in the
    # Hencky strain framework, the Newton solver converges at every step
    # with a relative tolerance of 1e-5, below the floor of the residual
    # reached with the deformation gradient (about 2e-5)
    add_test(NAME UniaxialTensileTest-Precision-LogarithmicStrainElasticity-DisplacementGradient
             COMMAND ./UniaxialTensileTest
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>
                     --behaviour LogarithmicStrainElasticity
                     --internal-state-variable EquivalentStrain
                     --reference-file
                     ${CMAKE_CURRENT_SOURCE_DIR}/references/LogarithmicStrainElasticity.ref
                     --displacement-gradient --relative-tolerance 1e-5
                     --stagnation-relative-tolerance 0
                     --no-parallel --no-post-processing)
    set_tests_properties(
      UniaxialTensileTest-Precision-LogarithmicStrainElasticity-DisplacementGradient
      PROPERTIES FAIL_REGULAR_EXPRESSION "the resolution failed;test failed"
                 DEPENDS BehaviourTest)
  endif()
  # a parameter specific to the Newton solver of the program (option
  # StagnationRelativeTolerance) is refused by the non linear solver override,
  # which only receives the parameters accepted by all the non linear
  # solvers: the refused parameter is reported, and the verbosity level of
  # the program (0) is kept (the iterations are printed)
  if(MFEM_MGIS_HAVE_TFEL_UTILITIES)
    add_test(NAME UniaxialTensileTest-NonLinearSolverOverride-OrthotropicElasticity
             COMMAND ./UniaxialTensileTest
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>
                     --behaviour OrthotropicElasticity
                     --internal-state-variable EquivalentStrain
                     --reference-file
                     ${CMAKE_CURRENT_SOURCE_DIR}/references/OrthotropicElasticity.ref
                     --stagnation-relative-tolerance 0.1
                     --no-parallel --no-post-processing
                     --non-linear-solver-override
                     "{Name: \"FixedPoint\", Options: {Method: \"ModifiedNewton\"}}")
    set_tests_properties(
      UniaxialTensileTest-NonLinearSolverOverride-OrthotropicElasticity
      PROPERTIES
      PASS_REGULAR_EXPRESSION
      "non linear solver override 'FixedPoint' refuses the parameters given by the program to its non linear solver: the parameters 'StagnationRelativeTolerance' are ignored.*Fixed point iteration"
      FAIL_REGULAR_EXPRESSION
      "the resolution failed;test failed;Process received signal;terminate called"
      DEPENDS BehaviourTest)
  endif(MFEM_MGIS_HAVE_TFEL_UTILITIES)

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
  # the linear solver built by the program (setLinearSolver taking a solver
  # handler) is replaced by the linear solver override given by the variable
  # MFEM_MGIS_LINEAR_SOLVER_OVERRIDE, the replacement being reported once
  if(MFEM_MGIS_HAVE_TFEL_UTILITIES)
    add_test(NAME ImposedPressureTest-IsotropicLinearElasticity-LinearSolverOverride
             COMMAND ./ImposedPressureTest
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>
                     --behaviour IsotropicLinearElasticity
                     --reference-file
                     ${CMAKE_CURRENT_SOURCE_DIR}/references/IsotropicLinearElasticity.ref
                     --parallel 0)
    set_tests_properties(
      ImposedPressureTest-IsotropicLinearElasticity-LinearSolverOverride
      PROPERTIES
      DEPENDS BehaviourTest
      ENVIRONMENT "MFEM_MGIS_LINEAR_SOLVER_OVERRIDE={Name: \"CGSolver\", Options: {RelativeTolerance: 1e-12, MaximumNumberOfIterations: 5000}}"
      PASS_REGULAR_EXPRESSION
      "linear solver override 'CGSolver' replaces the linear solver built by the program"
      FAIL_REGULAR_EXPRESSION
      "replaces the linear solver '';did not converge;valid parameters refused;misspelled parameter;Process received signal;terminate called")
  endif(MFEM_MGIS_HAVE_TFEL_UTILITIES)

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
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
               "--library" "$<TARGET_FILE:BehaviourTest>")
    add_test(NAME ${test}
             COMMAND ./AddBehaviourIntegratorTest ${params} "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test}
               COMMAND mpirun -n 2 ./AddBehaviourIntegratorTest ${params} "--parallel" 1)
      set_property(TEST Parallel${test} PROPERTY DEPENDS BehaviourTest)
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
  if(MFEM_MGIS_HAVE_SUITESPARSE)
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

  # Helper function declaring the variant of a test computing the
  # post-processings by kernels when possible (the other cases falling back
  # to the computations element by element)
  function(add_post_processings_kernels_variant test)
    add_test(NAME ${test}-kernels COMMAND ${ARGN})
    set_tests_properties(${test}-kernels PROPERTIES ENVIRONMENT
                         MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels)
  endfunction()

  # L2ProjectionTest executable
  add_test_executable(L2ProjectionTest)

  # Helper function for L2ProjectionTest cases
  function(add_l2_projection_test nsolver)
    set(test "L2ProjectionTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver"
               "${nsolver}")
    add_test(NAME ${test} COMMAND ./L2ProjectionTest ${params} "--parallel" "0")
    add_post_processings_kernels_variant(${test} ./L2ProjectionTest ${params}
                                         "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME Parallel${test} COMMAND mpirun -n 2 ./L2ProjectionTest ${params}
                                     "--parallel" "1")
      add_post_processings_kernels_variant(Parallel${test} mpirun -n 2
                                           ./L2ProjectionTest ${params}
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
    add_post_processings_kernels_variant(${test} ./GridFunctionTest ${params}
                                         "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./GridFunctionTest ${params}
                                     "--parallel" "1")
      add_post_processings_kernels_variant(${test}P mpirun -n 2
                                           ./GridFunctionTest ${params}
                                           "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_grid_function_test()

  # GridFunctionOnSubMeshTest executable
  add_tfel_test_executable(GridFunctionOnSubMeshTest)

  # Helper function for GridFunctionOnSubMeshTest
  function(add_grid_function_on_submesh_test test mesh)
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/${mesh}")
    add_test(NAME ${test} COMMAND ./GridFunctionOnSubMeshTest ${params}
                                  "--parallel" "0")
    add_post_processings_kernels_variant(${test} ./GridFunctionOnSubMeshTest
                                         ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./GridFunctionOnSubMeshTest
                                     ${params} "--parallel" "1")
      add_post_processings_kernels_variant(${test}P mpirun -n 2
                                           ./GridFunctionOnSubMeshTest
                                           ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_grid_function_on_submesh_test(GridFunctionOnSubMeshTest beam-tet.mesh)
  add_grid_function_on_submesh_test(GridFunctionOnSubMeshTest2D
                                    quad-2-materials.mesh)

  # AssignValuesTest executable
  add_tfel_test_executable(AssignValuesTest)

  # Helper function for AssignValuesTest
  function(add_assign_values_test test mesh)
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/${mesh}")
    add_test(NAME ${test} COMMAND ./AssignValuesTest ${params} "--parallel"
                                  "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./AssignValuesTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_assign_values_test(AssignValuesTest beam-tet.mesh)
  add_assign_values_test(AssignValuesTest2D quad-2-materials.mesh)

  # PostProcessingsKernelsTest executable: measures, integrals, values at the
  # nodes and L2 projections computed by kernels and element by element, on
  # hexahedra (two materials), tetrahedra and quadrangles, on the host, on the
  # debug device of MFEM, on two processes and on the GPU
  add_test_executable(PostProcessingsKernelsTest)
  function(add_post_processings_kernels_test name mesh hypothesis)
    set(test "PostProcessingsKernelsTest-${name}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/${mesh}"
               "--library" "$<TARGET_FILE:BehaviourTest>"
               "--hypothesis" "${hypothesis}" ${ARGN})
    add_test(NAME ${test} COMMAND ./PostProcessingsKernelsTest ${params})
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    add_test(NAME ${test}-debug
             COMMAND ./PostProcessingsKernelsTest ${params} --device debug)
    set_property(TEST ${test}-debug PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PostProcessingsKernelsTest
                                     ${params} --parallel)
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME ${test}-cuda
               COMMAND ./PostProcessingsKernelsTest ${params} --device cuda)
      set_tests_properties(${test}-cuda PROPERTIES
                           DEPENDS BehaviourTest RESOURCE_LOCK gpu)
      if(MFEM_USE_MPI)
        add_test(NAME ${test}P-cuda
                 COMMAND mpirun -n 2 ./PostProcessingsKernelsTest ${params}
                         --parallel --device cuda)
        set_tests_properties(${test}P-cuda PROPERTIES
                             DEPENDS BehaviourTest RESOURCE_LOCK gpu)
      endif(MFEM_USE_MPI)
    endif()
  endfunction()

  add_post_processings_kernels_test(hexahedra cube_2mat_per.mesh Tridimensional)
  add_post_processings_kernels_test(tetrahedra beam-tet.mesh Tridimensional)
  add_post_processings_kernels_test(quadrangles bar.msh PlaneStrain)
  # degrees of freedom ordered by vector dimension; elements of order 3,
  # whose values at the nodes are computed element by element (64 degrees of
  # freedom per hexahedron)
  add_post_processings_kernels_test(hexahedra-byVDIM cube_2mat_per.mesh
                                    Tridimensional --ordering byVDIM)
  add_post_processings_kernels_test(hexahedra-order3 cube_2mat_per.mesh
                                    Tridimensional --order 3)

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
    # MeshReadModeTest executable: errors of the mesh read modes reported to
    # the calling code instead of aborting the program
    add_test_executable(MeshReadModeTest)
    add_test(NAME MeshReadModeTestP
             COMMAND mpirun -n 2 ./MeshReadModeTest
                     ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh)
    set_tests_properties(MeshReadModeTestP PROPERTIES
                         PROCESSORS 2 ENVIRONMENT OMP_NUM_THREADS=1)
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
                     OrthotropicMicromorphicDamage
                     UO2Shrinkage_RAPHAEL2008)
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
      # the model UO2Shrinkage_RAPHAEL2008 (material UO2) is in its own
      # library
      set(_uo2_library
          "${_directory}/src/libGenericParallelUO2-${_backend}.so")
      file(MAKE_DIRECTORY "${_directory}")
      add_custom_command(
        OUTPUT "${_library}" "${_uo2_library}"
        COMMAND "${MFRONT}" --obuild --interface=generic-parallel
                ${MFEM_MGIS_MFRONT_REAL_TYPE_OPTION}
                "--configuration-file=${_configuration}"
                ${batch_integration_test_behaviours}
        WORKING_DIRECTORY "${_directory}"
        DEPENDS ${batch_integration_test_behaviours}
                ${MFEM_MGIS_MFRONT_DEPENDENCIES}
        COMMENT "MFront: generic-parallel behaviour (${_backend})")
      add_custom_target(GenericParallelBehaviour-${_backend}
                        DEPENDS "${_library}" "${_uo2_library}")
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
  if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS)
    add_test(NAME BatchIntegrationTestP
             COMMAND mpirun -n 2 ./BatchIntegrationTest --parallel
                     ${batch_integration_test_arguments}
                     ${batch_integration_test_libraries})
    set_tests_properties(BatchIntegrationTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
    # the debug device of MFEM in parallel, without GPU: the memory of the
    # device is checked (the parallel spaces have a prolongation, whose
    # products are made on the device)
    if(long_debug_device_tests)
      add_test(NAME BatchIntegrationTestP-debug
               COMMAND mpirun -n 2 ./BatchIntegrationTest --parallel
                       --device debug ${batch_integration_test_arguments}
                       ${batch_integration_test_libraries})
      set_tests_properties(BatchIntegrationTestP-debug
                           PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
    endif()
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME BatchIntegrationTestP-cuda
               COMMAND mpirun -n 2 ./BatchIntegrationTest --parallel
                       --device cuda ${batch_integration_test_arguments}
                       ${batch_integration_test_cuda_library})
      # behaviours, kernels and assembly on the GPU, hypre on the host
      add_test(NAME BatchIntegrationTestP-cuda-hypre-on-host
               COMMAND mpirun -n 2 ./BatchIntegrationTest --parallel
                       --device cuda --hypre-on-host
                       ${batch_integration_test_arguments}
                       ${batch_integration_test_cuda_library})
      set_tests_properties(BatchIntegrationTestP-cuda
                           BatchIntegrationTestP-cuda-hypre-on-host
                           PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
      set_tests_properties(BatchIntegrationTestP-cuda-hypre-on-host
                           PROPERTIES RESOURCE_LOCK gpu)
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
  # the sequential computations use UMFPack, or the linear solver of
  # SequentialLinearSolver.hxx without SuiteSparse
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
  if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS)
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
  # the sequential computations use UMFPack, or the linear solver of
  # SequentialLinearSolver.hxx without SuiteSparse
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
  if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS)
    add_test(NAME MicromorphicKernelsTestP
             COMMAND mpirun -n 2 ./MicromorphicKernelsTest --parallel
                     ${micromorphic_kernels_test_arguments})
    set_tests_properties(MicromorphicKernelsTestP
                         PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
  endif()

  # KINSOLSolverTest executable: non linear solver based on KINSOL
  # (SUNDIALS), compared to the Newton solver
  if(MFEM_MGIS_HAVE_SUNDIALS)
    add_test_executable(KINSOLSolverTest)
    add_test(NAME KINSOLSolverTest
             COMMAND KINSOLSolverTest
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>)
    if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS)
      add_test(NAME KINSOLSolverTestP
               COMMAND mpirun -n 2 ./KINSOLSolverTest --parallel
                       --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                       --library $<TARGET_FILE:BehaviourTest>)
      set_tests_properties(KINSOLSolverTestP
                           PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
    endif()
    # with a GPU, KINSOL computes on the device (vectors of CUDA or HIP)
    if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES"))
       AND enable-gpu-tests)
      if(MFEM_USE_CUDA STREQUAL "YES")
        set(kinsol_solver_test_device cuda)
      else()
        set(kinsol_solver_test_device hip)
      endif()
      add_test(NAME KINSOLSolverTest-gpu
               COMMAND KINSOLSolverTest --device ${kinsol_solver_test_device}
                       --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                       --library $<TARGET_FILE:BehaviourTest>)
      set_tests_properties(KINSOLSolverTest-gpu
                           PROPERTIES RESOURCE_LOCK gpu TIMEOUT 300)
    endif()
  endif()

  # SNESSolverTest executable: non linear solver based on the SNES of
  # PETSc, compared to the Newton solver, in parallel
  if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS AND MFEM_MGIS_HAVE_PETSC)
    add_test_executable(SNESSolverTest)
    add_test(NAME SNESSolverTestP
             COMMAND mpirun -n 2 ./SNESSolverTest
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>)
    set_tests_properties(SNESSolverTestP PROPERTIES
      ENVIRONMENT OMP_NUM_THREADS=1 FAIL_REGULAR_EXPRESSION "check failed")
    if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES")) AND
       enable-gpu-tests)
      if(MFEM_USE_CUDA STREQUAL "YES")
        set(snes_solver_test_device cuda)
      else()
        set(snes_solver_test_device hip)
      endif()
      add_test(NAME SNESSolverTestP-gpu
               COMMAND mpirun -n 2 ./SNESSolverTest
                       --device ${snes_solver_test_device}
                       --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                       --library $<TARGET_FILE:BehaviourTest>)
      set_tests_properties(SNESSolverTestP-gpu PROPERTIES
        ENVIRONMENT OMP_NUM_THREADS=1 FAIL_REGULAR_EXPRESSION "check failed"
        RESOURCE_LOCK gpu TIMEOUT 600)
    endif()
  endif()

  # FixedPointSolverTest executable: fixed point solver (methods,
  # accelerations, linear solvers of the elastic operator) compared to the
  # Newton solver
  add_test_executable(FixedPointSolverTest)
  add_test(NAME FixedPointSolverTest
           COMMAND FixedPointSolverTest
                   --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                   --library $<TARGET_FILE:BehaviourTest>)
  if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS)
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
  # stiffness matrices added in a fixed order (`TangentOperatorAssembly`):
  # the same results at each run, also on the GPU, where the atomic additions
  # make the accelerations of the fixed point solver stagnate in rare runs
  set(fixed_point_solver_test_ordered ${fixed_point_solver_test_kernels}
      MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered)
  set(fixed_point_solver_test_ordered_devices cpu)
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    list(APPEND fixed_point_solver_test_ordered_devices cuda)
  endif()
  foreach(_device ${fixed_point_solver_test_ordered_devices})
    add_test(NAME FixedPointSolverTest-kernels-ordered-${_device}
             COMMAND FixedPointSolverTest --device ${_device}
                     --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
                     --library $<TARGET_FILE:BehaviourTest>)
    set_tests_properties(FixedPointSolverTest-kernels-ordered-${_device}
                         PROPERTIES ENVIRONMENT
                         "${fixed_point_solver_test_ordered}")
    if(_device STREQUAL "cuda")
      set_tests_properties(FixedPointSolverTest-kernels-ordered-${_device}
                           PROPERTIES RESOURCE_LOCK gpu)
    endif()
  endforeach()

  # MaterialArraysTest executable: arrays of the materials in the memory of
  # MFEM, on the host, on the debug device of MFEM and on the GPU. The
  # point-wise model is also integrated by batches with the library
  # generated by the generic-parallel interface for the device of the test
  # (openmp on the host and on the debug device, cuda on the GPU), if any
  add_test_executable(MaterialArraysTest)
  mfem_mgis_compile_as_device_code(MaterialArraysTest MaterialArraysTest.cxx)
  set(material_arrays_test_arguments
      --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
      --multi-material-mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh
      --mesh-2d ${CMAKE_CURRENT_SOURCE_DIR}/bar.msh
      --library $<TARGET_FILE:BehaviourTest>)
  set(material_arrays_test_host_library)
  if(TARGET GenericParallelBehaviour-openmp)
    set(material_arrays_test_host_library
        --generic-parallel-library
        "${CMAKE_CURRENT_BINARY_DIR}/generic-parallel-openmp/src/libGenericParallelUO2-openmp.so")
  endif()
  add_test(NAME MaterialArraysTest
           COMMAND MaterialArraysTest ${material_arrays_test_arguments}
                   ${material_arrays_test_host_library})
  add_test(NAME MaterialArraysTest-debug
           COMMAND MaterialArraysTest --device debug
                   ${material_arrays_test_arguments}
                   ${material_arrays_test_host_library})
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    set(material_arrays_test_cuda_library)
    if(TARGET GenericParallelBehaviour-cuda)
      set(material_arrays_test_cuda_library
          --generic-parallel-library
          "${CMAKE_CURRENT_BINARY_DIR}/generic-parallel-cuda/src/libGenericParallelUO2-cuda.so")
    endif()
    add_test(NAME MaterialArraysTest-cuda
             COMMAND MaterialArraysTest --device cuda
                     ${material_arrays_test_arguments}
                     ${material_arrays_test_cuda_library})
    # the point-wise model generated by the kokkos backend, integrated on
    # the GPU, Kokkos (built for CUDA) being initialized by
    # mfem_mgis::initialize
    get_filename_component(_mfront_bin "${MFRONT}" DIRECTORY)
    if(MFEM_MGIS_HAVE_KOKKOS AND ("CUDA" IN_LIST Kokkos_DEVICES) AND
       EXISTS "${_mfront_bin}/../share/tfel/mfront/configuration-files/config-kokkos.json")
      set(_directory "${CMAKE_CURRENT_BINARY_DIR}/generic-parallel-kokkos")
      set(_library "${_directory}/src/libGenericParallelUO2-kokkos.so")
      file(MAKE_DIRECTORY "${_directory}")
      add_custom_command(
        OUTPUT "${_library}"
        COMMAND "${MFRONT}" --obuild --interface=generic-parallel
                --configuration-file=config-kokkos.json
                ${MFEM_MGIS_MFRONT_REAL_TYPE_OPTION}
                "${CMAKE_CURRENT_SOURCE_DIR}/UO2Shrinkage_RAPHAEL2008.mfront"
        WORKING_DIRECTORY "${_directory}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/UO2Shrinkage_RAPHAEL2008.mfront"
                ${MFEM_MGIS_MFRONT_DEPENDENCIES}
        COMMENT "MFront: generic-parallel model (kokkos)")
      add_custom_target(GenericParallelBehaviour-kokkos DEPENDS "${_library}")
      add_dependencies(build-tests GenericParallelBehaviour-kokkos)
      add_test(NAME MaterialArraysTest-kokkos
               COMMAND MaterialArraysTest --device cuda --kokkos
                       ${material_arrays_test_arguments}
                       --generic-parallel-library "${_library}")
      set_tests_properties(MaterialArraysTest-kokkos PROPERTIES
                           RESOURCE_LOCK gpu)
    endif()
  endif()

  # GradientsKernelsTest executable: gradients of the behaviours computed
  # by kernels, on the host, on the debug device of MFEM and on the GPU
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
  # kernels, on the host and on the debug device of MFEM (UMFPack, or the
  # linear solver of SequentialLinearSolver.hxx without SuiteSparse)
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
  # ExecutionSettingsTest executable: execution settings (option
  # --execution and locations of the phases in the mixed mode), uniaxial
  # tension computed where the settings place each phase, refusal of the
  # behaviour libraries and of the linear solvers computing elsewhere, and
  # refusal of the contradictory options by initialize
  if(MFEM_USE_MPI)
    add_test_executable(ExecutionSettingsTest)
    # behaviour libraries for the CPU and, if MFEM uses a GPU, for the
    # GPU, chosen by getBehavioursLibrary
    include(${PROJECT_SOURCE_DIR}/cmake/modules/MFEMMGISBehaviours.cmake)
    # the library for the GPU is only built for the tests on the GPU
    set(_execution_settings_no_gpu NO_GPU)
    if(enable-gpu-tests)
      set(_execution_settings_no_gpu)
    endif()
    mfem_mgis_add_behaviours_library(ExecutionSettingsBehaviours
      SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/IsotropicLinearHardeningPlasticity.mfront
      EXCLUDE_FROM_ALL ${_execution_settings_no_gpu})
    add_dependencies(build-tests ExecutionSettingsBehaviours)
    set(execution_settings_test_arguments
        --mesh ${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh
        --mixed-mesh
        ${CMAKE_CURRENT_SOURCE_DIR}/square-quadrangle-triangles.mesh
        --behaviours-directory
        ${ExecutionSettingsBehaviours_BEHAVIOURS_DIRECTORY})
    # library generated by the generic-parallel interface for the host
    # (openmp backend), used too when the behaviours are on the CPU
    list(FIND batch_integration_test_libraries --openmp-library _openmp)
    if(NOT _openmp EQUAL -1)
      math(EXPR _openmp "${_openmp} + 1")
      list(GET batch_integration_test_libraries ${_openmp} _openmp_library)
      list(APPEND execution_settings_test_arguments
           --host-library ${_openmp_library})
    endif()
    add_test(NAME ExecutionSettingsTest
             COMMAND ExecutionSettingsTest
                     ${execution_settings_test_arguments})
    add_test(NAME ExecutionSettingsTest-cpu
             COMMAND ExecutionSettingsTest --execution cpu
                     --expected-settings cpu,cpu,cpu,cpu,cpu
                     ${execution_settings_test_arguments})
    # the environment variables placing the computations are ignored
    # when the execution mode is given
    add_test(NAME ExecutionSettingsTest-cpu-environment
             COMMAND ExecutionSettingsTest --execution cpu
                     --expected-settings cpu,cpu,cpu,cpu,cpu
                     ${execution_settings_test_arguments})
    set_tests_properties(ExecutionSettingsTest-cpu-environment PROPERTIES
      ENVIRONMENT "MFEM_MGIS_DEVICE=cuda;MFEM_MGIS_HYPRE_ON_HOST=1;MFEM_MGIS_BEHAVIOUR_INTEGRATION=PointByPoint;MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels;MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels;MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels;MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels")
    # the debug device of MFEM computes on the host
    add_test(NAME ExecutionSettingsTest-cpu-debug
             COMMAND ExecutionSettingsTest --execution cpu --device debug
                     --expected-settings cpu,cpu,cpu,cpu,cpu
                     ${execution_settings_test_arguments})
    function(add_execution_settings_refusal_test name message)
      add_test(NAME ExecutionSettingsTest-refused-${name}
               COMMAND ExecutionSettingsTest ${ARGN})
      set_tests_properties(ExecutionSettingsTest-refused-${name}
                           PROPERTIES PASS_REGULAR_EXPRESSION "${message}")
    endfunction()
    add_execution_settings_refusal_test(
      location "option '--finite-elements-execution' requires '--execution mixed'"
      --execution cpu --finite-elements-execution gpu)
    # an invalid linear solver override (or, without the utilities of TFEL,
    # any linear solver override) is refused by mfem_mgis::initialize
    add_execution_settings_refusal_test(
      linear-solver-override "invalid linear solver override"
      --linear-solver-override "{Name: 1}")
    add_execution_settings_refusal_test(
      linear-solver-override-twice "--linear-solver-override is given twice"
      --linear-solver-override CGSolver --linear-solver-override CGSolver)
    add_execution_settings_refusal_test(
      preconditioner-override "invalid preconditioner override"
      --preconditioner-override "{Name: \"HypreBoomerAMG\", Unknown: 1}")
    # the names are only checked once the descriptions are read by the
    # utilities of TFEL
    if(MFEM_MGIS_HAVE_TFEL_UTILITIES)
      add_execution_settings_refusal_test(
        preconditioner-override-unknown "unknown name 'UnknownPreconditioner'"
        --preconditioner-override UnknownPreconditioner)
      add_execution_settings_refusal_test(
        preconditioner-override-none "no options expected for 'None'"
        --preconditioner-override "{Name: None, Options: {VerbosityLevel: 0}}")
      add_execution_settings_refusal_test(
        non-linear-solver-override-unknown "unknown name 'UnknownSolver'"
        --non-linear-solver-override UnknownSolver)
    endif(MFEM_MGIS_HAVE_TFEL_UTILITIES)
    add_execution_settings_refusal_test(
      non-linear-solver-override "invalid non linear solver override"
      --non-linear-solver-override "{Name: 1}")
    add_execution_settings_refusal_test(
      mode "invalid value 'xpu' for option '--execution'" --execution xpu)
    add_execution_settings_refusal_test(
      twice "option '--execution' given twice" --execution cpu -x cpu)
    add_execution_settings_refusal_test(
      cpu-device "the device 'cuda' uses a GPU, which contradicts"
      --execution cpu --device cuda)
    add_execution_settings_refusal_test(
      missing-value "no value given to '--execution'"
      --execution --behaviours-execution cpu)
    add_execution_settings_refusal_test(
      cpu-pool-size "can't be given with '--execution cpu'"
      --execution cpu --hypre-device-pool-size 1)
    if((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES"))
      # refused before any use of the GPU
      add_execution_settings_refusal_test(
        location-value "invalid value 'tpu' for option '--behaviours-execution'"
        --execution mixed --behaviours-execution tpu)
      add_execution_settings_refusal_test(
        host-device "the device 'debug' does not use a GPU"
        --execution gpu --device debug)
      add_execution_settings_refusal_test(
        hypre-on-host
        "contradicts the location of the linear solvers"
        --execution gpu --hypre-on-host)
      add_execution_settings_refusal_test(
        no-hypre-on-host
        "contradicts the location of the linear solvers"
        --execution mixed --linear-solvers-execution cpu --no-hypre-on-host)
    endif()
    if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES"))
       AND enable-gpu-tests)
      set(execution_settings_test_gpu_arguments
          ${execution_settings_test_arguments})
      # the variants placing the behaviours on the GPU require their
      # library for the GPU
      set(_bh_locations cpu)
      if(ExecutionSettingsBehaviours_HAVE_GPU_BEHAVIOURS)
        list(APPEND _bh_locations gpu)
        add_test(NAME ExecutionSettingsTest-gpu
                 COMMAND ExecutionSettingsTest --execution gpu
                         --expected-settings gpu,gpu,gpu,gpu,gpu
                         ${execution_settings_test_gpu_arguments})
        add_test(NAME ExecutionSettingsTest-mixed-gpu
                 COMMAND ExecutionSettingsTest --execution mixed
                         --expected-settings mixed,gpu,gpu,gpu,gpu
                         ${execution_settings_test_gpu_arguments})
        add_test(NAME ExecutionSettingsTest-mixed-environment-gpu
                 COMMAND ExecutionSettingsTest --execution mixed
                         --expected-settings mixed,gpu,gpu,gpu,gpu
                         ${execution_settings_test_gpu_arguments})
        set_tests_properties(ExecutionSettingsTest-mixed-environment-gpu
          PROPERTIES ENVIRONMENT
          "MFEM_MGIS_DEVICE=debug;MFEM_MGIS_HYPRE_ON_HOST=1;MFEM_MGIS_BEHAVIOUR_INTEGRATION=PointByPoint;MFEM_MGIS_GRADIENTS_COMPUTATION=ElementByElement;MFEM_MGIS_INNER_FORCES_COMPUTATION=ElementByElement;MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=ElementByElement;MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=ElementByElement")
      endif()
      # the 16 locations of the four phases in the mixed mode
      foreach(_fe cpu gpu)
        foreach(_bh ${_bh_locations})
          foreach(_lsd cpu gpu)
            foreach(_ppd cpu gpu)
              set(_name "${_fe}-${_bh}-${_lsd}-${_ppd}")
              add_test(NAME ExecutionSettingsTest-mixed-${_name}-gpu
                       COMMAND ExecutionSettingsTest --execution mixed
                               --finite-elements-execution ${_fe}
                               --behaviours-execution ${_bh}
                               --linear-solvers-execution ${_lsd}
                               --post-processings-execution ${_ppd}
                               --expected-settings
                               mixed,${_fe},${_bh},${_lsd},${_ppd}
                               ${execution_settings_test_gpu_arguments})
            endforeach()
          endforeach()
        endforeach()
      endforeach()
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
    # hypre on the host with a device on the GPU: the levels of fill of the
    # host are available
    if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
      add_test(NAME HypreILUTestP-cuda-hypre-on-host
               COMMAND mpirun -n 2 ./HypreILUTest --device cuda
                       --hypre-on-host)
      set_tests_properties(HypreILUTestP-cuda-hypre-on-host
                           PROPERTIES RESOURCE_LOCK gpu)
    endif()
    # HypreParaSailsTest executable: ParaSails preconditioner given a new
    # operator after the destruction of the previous one, refused when hypre
    # computes on a GPU. In float, hypre must be built with its internal
    # LAPACK: with an external LAPACK, hypre 3.2.0 in single precision calls
    # the LAPACK routines of the double type (dpotrf, dpotrs) on float arrays
    # and ParaSails does not converge
    add_test_executable(HypreParaSailsTest)
    add_test(NAME HypreParaSailsTest COMMAND HypreParaSailsTest)
    add_test(NAME HypreParaSailsTestP
             COMMAND mpirun -n 2 ./HypreParaSailsTest)
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

# DirectSolverReuseAnalysisTest executable: analysis of the direct solvers
# (MUMPS, cuDSS) reused for a new operator only if its pattern is the same
# on all processes. A wrong reuse may block the solver: the tests have a
# time limit
if(MFEM_USE_MPI AND (MFEM_MGIS_HAVE_MUMPS OR MFEM_MGIS_HAVE_CUDSS))
  add_test_executable(DirectSolverReuseAnalysisTest)
endif()
if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MUMPS)
  add_test(NAME DirectSolverReuseAnalysisTest
           COMMAND DirectSolverReuseAnalysisTest)
  add_test(NAME DirectSolverReuseAnalysisTestP
           COMMAND mpirun -n 2 ./DirectSolverReuseAnalysisTest)
  # three processes: layout of the first process kept when the one of the
  # others changes
  add_test(NAME DirectSolverReuseAnalysisTestP3
           COMMAND mpirun -n 3 ./DirectSolverReuseAnalysisTest)
  set_tests_properties(DirectSolverReuseAnalysisTest
                       DirectSolverReuseAnalysisTestP
                       DirectSolverReuseAnalysisTestP3
                       PROPERTIES TIMEOUT 120 ENVIRONMENT OMP_NUM_THREADS=1)
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME DirectSolverReuseAnalysisTestP-cuda
             COMMAND mpirun -n 2 ./DirectSolverReuseAnalysisTest
                     --device cuda)
    set_tests_properties(DirectSolverReuseAnalysisTestP-cuda
                         PROPERTIES RESOURCE_LOCK gpu TIMEOUT 120
                                    ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
endif()
if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_CUDSS)
  # cuDSS refused without a device on a CUDA GPU, even when the test suite
  # is run on a GPU (variable MFEM_MGIS_DEVICE)
  add_test(NAME DirectSolverReuseAnalysisTest-cudss-refused
           COMMAND DirectSolverReuseAnalysisTest --solver CuDSSSolver)
  set_tests_properties(DirectSolverReuseAnalysisTest-cudss-refused
                       PROPERTIES TIMEOUT 120 ENVIRONMENT MFEM_MGIS_DEVICE=
                                  ENVIRONMENT_MODIFICATION MFEM_DEVICE=unset:)
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME DirectSolverReuseAnalysisTest-cudss-cuda
             COMMAND DirectSolverReuseAnalysisTest --solver CuDSSSolver
                     --device cuda)
    add_test(NAME DirectSolverReuseAnalysisTestP-cudss-cuda
             COMMAND mpirun -n 2 ./DirectSolverReuseAnalysisTest
                     --solver CuDSSSolver --device cuda)
    add_test(NAME DirectSolverReuseAnalysisTestP3-cudss-cuda
             COMMAND mpirun -n 3 ./DirectSolverReuseAnalysisTest
                     --solver CuDSSSolver --device cuda)
    # hypre on the host: the matrices are stored on the host (accepted on a
    # GPU reading the pageable memory of the host, refused otherwise)
    add_test(NAME DirectSolverReuseAnalysisTest-cudss-cuda-hypre-on-host
             COMMAND DirectSolverReuseAnalysisTest --solver CuDSSSolver
                     --device cuda --hypre-on-host)
    add_test(NAME DirectSolverReuseAnalysisTestP-cudss-cuda-hypre-on-host
             COMMAND mpirun -n 2 ./DirectSolverReuseAnalysisTest
                     --solver CuDSSSolver --device cuda --hypre-on-host)
    set_tests_properties(DirectSolverReuseAnalysisTest-cudss-cuda
                         DirectSolverReuseAnalysisTestP-cudss-cuda
                         DirectSolverReuseAnalysisTestP3-cudss-cuda
                         DirectSolverReuseAnalysisTest-cudss-cuda-hypre-on-host
                         DirectSolverReuseAnalysisTestP-cudss-cuda-hypre-on-host
                         PROPERTIES RESOURCE_LOCK gpu TIMEOUT 120
                                    ENVIRONMENT OMP_NUM_THREADS=1)
    # three pools of hypre of the default size (4 GiB each), with the
    # memory of cuDSS, nearly fill a GPU of 16 GiB: the test failed for
    # lack of memory when other processes used the GPU
    set_tests_properties(DirectSolverReuseAnalysisTestP3-cudss-cuda
                         PROPERTIES ENVIRONMENT
                         "OMP_NUM_THREADS=1;MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1")
  endif()
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

# LinearSolverFactoriesTest executable: each linear solver and each
# preconditioner registered in the factories (generic test, see the test)
add_test_executable(LinearSolverFactoriesTest)
add_test(NAME LinearSolverFactoriesTest COMMAND LinearSolverFactoriesTest)
# the same on a beam in linear elasticity (in single precision, with looser
# bounds and the entries which do not converge listed by the test)
add_test(NAME LinearSolverFactoriesTest-elasticity
         COMMAND LinearSolverFactoriesTest --elasticity)
set_tests_properties(LinearSolverFactoriesTest
                     LinearSolverFactoriesTest-elasticity PROPERTIES
  FAIL_REGULAR_EXPRESSION "check failed")
if(MFEM_USE_MPI)
  add_test(NAME LinearSolverFactoriesTestP
           COMMAND mpirun -n 2 ./LinearSolverFactoriesTest --parallel)
  add_test(NAME LinearSolverFactoriesTestP-elasticity
           COMMAND mpirun -n 2 ./LinearSolverFactoriesTest --parallel
                   --elasticity)
  set_tests_properties(LinearSolverFactoriesTestP
                       LinearSolverFactoriesTestP-elasticity PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
endif(MFEM_USE_MPI)
# on a CUDA GPU, the entries which only compute on the GPU (AmgX, cuDSS) are
# also checked
if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
  add_test(NAME LinearSolverFactoriesTest-cuda
           COMMAND LinearSolverFactoriesTest --device cuda)
  add_test(NAME LinearSolverFactoriesTest-cuda-elasticity
           COMMAND LinearSolverFactoriesTest --device cuda --elasticity)
  set(linear_solver_factories_gpu_tests LinearSolverFactoriesTest-cuda
      LinearSolverFactoriesTest-cuda-elasticity)
  if(MFEM_USE_MPI)
    add_test(NAME LinearSolverFactoriesTestP-cuda
             COMMAND mpirun -n 2 ./LinearSolverFactoriesTest --parallel
                     --device cuda)
    list(APPEND linear_solver_factories_gpu_tests
         LinearSolverFactoriesTestP-cuda)
  endif(MFEM_USE_MPI)
  set_tests_properties(${linear_solver_factories_gpu_tests} PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# GinkgoSolversTest executable: solvers and preconditioners of Ginkgo, on
# the device of MFEM
if(MFEM_MGIS_HAVE_GINKGO)
  add_test_executable(GinkgoSolversTest)
  add_test(NAME GinkgoSolversTest COMMAND GinkgoSolversTest)
  set(ginkgo_tests GinkgoSolversTest)
  if(MFEM_USE_MPI)
    add_test(NAME GinkgoSolversTestP
             COMMAND mpirun -n 2 ./GinkgoSolversTest --parallel)
    list(APPEND ginkgo_tests GinkgoSolversTestP)
  endif(MFEM_USE_MPI)
  set_tests_properties(${ginkgo_tests} PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
  if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
    add_test(NAME GinkgoSolversTest-cuda
             COMMAND GinkgoSolversTest --device cuda)
    set(ginkgo_gpu_tests GinkgoSolversTest-cuda)
    if(MFEM_USE_MPI)
      add_test(NAME GinkgoSolversTestP-cuda
               COMMAND mpirun -n 2 ./GinkgoSolversTest --parallel
                       --device cuda)
      list(APPEND ginkgo_gpu_tests GinkgoSolversTestP-cuda)
    endif(MFEM_USE_MPI)
    set_tests_properties(${ginkgo_gpu_tests} PROPERTIES
      FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
  endif()
endif()

# MFEMSolversTest executable: smoothers and incomplete factorization of the
# sparse matrices of MFEM, direct solvers of SuiteSparse used as
# preconditioners or solvers, orthogonalizing wrapper of a preconditioner
add_test_executable(MFEMSolversTest)
add_test(NAME MFEMSolversTest COMMAND MFEMSolversTest)
set_tests_properties(MFEMSolversTest PROPERTIES
  FAIL_REGULAR_EXPRESSION "check failed")
if(MFEM_USE_MPI)
  add_test(NAME MFEMSolversTestP COMMAND mpirun -n 2 ./MFEMSolversTest -p)
  set_tests_properties(MFEMSolversTestP PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
endif(MFEM_USE_MPI)
# on a GPU, the matrices and the vectors are copied on the host
if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES"))
   AND enable-gpu-tests)
  if(MFEM_USE_CUDA STREQUAL "YES")
    set(mfem_solvers_test_device cuda)
  else()
    set(mfem_solvers_test_device hip)
  endif()
  add_test(NAME MFEMSolversTest-gpu
           COMMAND MFEMSolversTest --device ${mfem_solvers_test_device})
  set(mfem_solvers_gpu_tests MFEMSolversTest-gpu)
  if(MFEM_USE_MPI)
    add_test(NAME MFEMSolversTestP-gpu
             COMMAND mpirun -n 2 ./MFEMSolversTest -p
                     --device ${mfem_solvers_test_device})
    list(APPEND mfem_solvers_gpu_tests MFEMSolversTestP-gpu)
  endif(MFEM_USE_MPI)
  set_tests_properties(${mfem_solvers_gpu_tests} PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# MKLSolversTest executable: sparse direct solvers of oneMKL, PARDISO in
# sequential and the Cluster Sparse Solver in parallel
if(MFEM_MGIS_HAVE_MKL_PARDISO)
  add_test_executable(MKLSolversTest)
  add_test(NAME MKLSolversTest COMMAND MKLSolversTest)
  set_tests_properties(MKLSolversTest PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
endif()
if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MKL_CPARDISO)
  if(NOT MFEM_MGIS_HAVE_MKL_PARDISO)
    add_test_executable(MKLSolversTest)
  endif()
  add_test(NAME MKLSolversTestP COMMAND mpirun -n 2 ./MKLSolversTest -p)
  set_tests_properties(MKLSolversTestP PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
endif()
# on a GPU, the matrices and the vectors are copied on the host
if((MFEM_MGIS_HAVE_MKL_PARDISO OR
    MFEM_MGIS_HAVE_MKL_CPARDISO) AND
   ((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES")) AND
   enable-gpu-tests)
  if(MFEM_USE_CUDA STREQUAL "YES")
    set(mkl_solvers_test_device cuda)
  else()
    set(mkl_solvers_test_device hip)
  endif()
  set(mkl_solvers_gpu_tests)
  if(MFEM_MGIS_HAVE_MKL_PARDISO)
    add_test(NAME MKLSolversTest-gpu
             COMMAND MKLSolversTest --device ${mkl_solvers_test_device})
    list(APPEND mkl_solvers_gpu_tests MKLSolversTest-gpu)
  endif()
  if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_MKL_CPARDISO)
    add_test(NAME MKLSolversTestP-gpu
             COMMAND mpirun -n 2 ./MKLSolversTest -p
                     --device ${mkl_solvers_test_device})
    list(APPEND mkl_solvers_gpu_tests MKLSolversTestP-gpu)
  endif()
  set_tests_properties(${mkl_solvers_gpu_tests} PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# NodalBlockJacobiTest executable: Jacobi preconditioner by nodal blocks,
# with the backends of the batched linear algebra of MFEM available
add_test_executable(NodalBlockJacobiTest)
add_test(NAME NodalBlockJacobiTest COMMAND NodalBlockJacobiTest)
set(nodal_block_jacobi_tests NodalBlockJacobiTest)
if(MFEM_USE_MPI)
  add_test(NAME NodalBlockJacobiTestP
           COMMAND mpirun -n 2 ./NodalBlockJacobiTest -p)
  list(APPEND nodal_block_jacobi_tests NodalBlockJacobiTestP)
endif()
set_tests_properties(${nodal_block_jacobi_tests} PROPERTIES
  FAIL_REGULAR_EXPRESSION "check failed")
if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES")) AND
   enable-gpu-tests)
  if(MFEM_USE_CUDA STREQUAL "YES")
    set(nodal_block_jacobi_test_device cuda)
  else()
    set(nodal_block_jacobi_test_device hip)
  endif()
  add_test(NAME NodalBlockJacobiTest-gpu
           COMMAND NodalBlockJacobiTest
                   --device ${nodal_block_jacobi_test_device})
  set(nodal_block_jacobi_gpu_tests NodalBlockJacobiTest-gpu)
  if(MFEM_USE_MPI)
    add_test(NAME NodalBlockJacobiTestP-gpu
             COMMAND mpirun -n 2 ./NodalBlockJacobiTest -p
                     --device ${nodal_block_jacobi_test_device})
    list(APPEND nodal_block_jacobi_gpu_tests NodalBlockJacobiTestP-gpu)
  endif()
  set_tests_properties(${nodal_block_jacobi_gpu_tests} PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# HypreMGRTest executable: multigrid reduction of hypre on a beam in linear
# elasticity, with one and two levels of reduction, the components of the
# unknowns being given to hypre by markers, whatever their ordering
if(MFEM_USE_MPI)
  add_test_executable(HypreMGRTest)
  add_test(NAME HypreMGRTest COMMAND HypreMGRTest)
  add_test(NAME HypreMGRTestP COMMAND mpirun -n 2 ./HypreMGRTest)
  set_tests_properties(HypreMGRTest HypreMGRTestP PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
  if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES")) AND
     enable-gpu-tests)
    if(MFEM_USE_CUDA STREQUAL "YES")
      set(hypre_mgr_test_device cuda)
    else()
      set(hypre_mgr_test_device hip)
    endif()
    add_test(NAME HypreMGRTestP-gpu
             COMMAND mpirun -n 2 ./HypreMGRTest
                     --device ${hypre_mgr_test_device})
    set_tests_properties(HypreMGRTestP-gpu PROPERTIES
      FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
  endif()
endif()

# MAGMASolversTest executable: sparse solvers of MAGMA, which always
# compute on a GPU: Krylov solvers in sequential, preconditioners applied by
# the conjugate gradient of MFEM in sequential and in parallel; the vectors
# are copied on the GPU with the device cpu of MFEM, and given to MAGMA
# with the device of the GPU
if(MFEM_MGIS_HAVE_MAGMA_SPARSE AND enable-gpu-tests)
  add_test_executable(MAGMASolversTest)
  if(MFEM_USE_CUDA STREQUAL "YES")
    set(magma_solvers_test_device cuda)
  else()
    set(magma_solvers_test_device hip)
  endif()
  add_test(NAME MAGMASolversTest COMMAND MAGMASolversTest)
  add_test(NAME MAGMASolversTest-gpu
           COMMAND MAGMASolversTest --device ${magma_solvers_test_device})
  set(magma_solvers_tests MAGMASolversTest MAGMASolversTest-gpu)
  if(MFEM_USE_MPI)
    add_test(NAME MAGMASolversTestP COMMAND mpirun -n 2 ./MAGMASolversTest -p)
    add_test(NAME MAGMASolversTestP-gpu
             COMMAND mpirun -n 2 ./MAGMASolversTest -p
                     --device ${magma_solvers_test_device})
    list(APPEND magma_solvers_tests MAGMASolversTestP MAGMASolversTestP-gpu)
  endif()
  set_tests_properties(${magma_solvers_tests} PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# HypreSolversTest executable: solvers and preconditioners of hypre which
# are not wrapped by MFEM, smoothers of hypre and options of BoomerAMG
if(MFEM_USE_MPI)
  add_test_executable(HypreSolversTest)
  add_test(NAME HypreSolversTest COMMAND HypreSolversTest)
  add_test(NAME HypreSolversTestP COMMAND mpirun -n 2 ./HypreSolversTest)
  set_tests_properties(HypreSolversTest HypreSolversTestP PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
endif(MFEM_USE_MPI)

# SparseDirectSolversTest executable: direct solvers of MUMPS, SuperLU_DIST
# and STRUMPACK, preconditioners of MUMPS and STRUMPACK
if(MFEM_USE_MPI AND (MFEM_MGIS_HAVE_MUMPS OR
                     MFEM_MGIS_HAVE_SUPERLU OR
                     MFEM_MGIS_HAVE_STRUMPACK))
  add_test_executable(SparseDirectSolversTest)
  add_test(NAME SparseDirectSolversTest COMMAND SparseDirectSolversTest)
  add_test(NAME SparseDirectSolversTestP
           COMMAND mpirun -n 2 ./SparseDirectSolversTest)
  set_tests_properties(SparseDirectSolversTest SparseDirectSolversTestP
    PROPERTIES FAIL_REGULAR_EXPRESSION "check failed")
endif()

# PETScSolversTest executable: Krylov solver and preconditioner of PETSc
# configured by the options of PETSc
if(MFEM_USE_MPI AND MFEM_MGIS_HAVE_PETSC)
  add_test_executable(PETScSolversTest)
  add_test(NAME PETScSolversTest COMMAND PETScSolversTest)
  add_test(NAME PETScSolversTestP COMMAND mpirun -n 2 ./PETScSolversTest)
  set_tests_properties(PETScSolversTest PETScSolversTestP PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed")
  # with MFEM on a GPU, the solvers of PETSc compute on the GPU if PETSc is
  # built for it, on the host otherwise
  if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES")) AND
     enable-gpu-tests)
    if(MFEM_USE_CUDA STREQUAL "YES")
      set(petsc_solvers_test_device cuda)
    else()
      set(petsc_solvers_test_device hip)
    endif()
    add_test(NAME PETScSolversTest-gpu
             COMMAND PETScSolversTest --device ${petsc_solvers_test_device})
    add_test(NAME PETScSolversTestP-gpu
             COMMAND mpirun -n 2 ./PETScSolversTest
                     --device ${petsc_solvers_test_device})
    set_tests_properties(PETScSolversTest-gpu PETScSolversTestP-gpu
      PROPERTIES FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu
                 TIMEOUT 600)
  endif()
endif()

# TrilinosSolversTest executable: solvers and preconditioners of Trilinos.
# With MFEM on a GPU, Tpetra computes on the GPU if Trilinos is built with
# CUDA or HIP; Trilinos built with CUDA or HIP uses the GPU in all the tests
# (lock gpu). With three processes, a process in the middle owns no
# unknown
if(MFEM_MGIS_HAVE_TRILINOS)
  add_test_executable(TrilinosSolversTest)
  if(MFEM_MGIS_HAVE_TEKO)
    target_compile_definitions(TrilinosSolversTest PRIVATE MFEM_MGIS_HAVE_TEKO)
  endif()
  add_test(NAME TrilinosSolversTest COMMAND TrilinosSolversTest)
  add_test(NAME TrilinosSolversTestP COMMAND mpirun -n 2 ./TrilinosSolversTest)
  add_test(NAME TrilinosSolversTestP3
           COMMAND mpirun -n 3 ./TrilinosSolversTest)
  set_tests_properties(TrilinosSolversTest TrilinosSolversTestP
                       TrilinosSolversTestP3 PROPERTIES
    FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu TIMEOUT 600)
  if(((MFEM_USE_CUDA STREQUAL "YES") OR (MFEM_USE_HIP STREQUAL "YES")) AND
     enable-gpu-tests)
    if(MFEM_USE_CUDA STREQUAL "YES")
      set(trilinos_solvers_test_device cuda)
    else()
      set(trilinos_solvers_test_device hip)
    endif()
    add_test(NAME TrilinosSolversTest-gpu
             COMMAND TrilinosSolversTest
                     --device ${trilinos_solvers_test_device})
    add_test(NAME TrilinosSolversTestP-gpu
             COMMAND mpirun -n 2 ./TrilinosSolversTest
                     --device ${trilinos_solvers_test_device})
    set_tests_properties(TrilinosSolversTest-gpu TrilinosSolversTestP-gpu
      PROPERTIES FAIL_REGULAR_EXPRESSION "check failed" RESOURCE_LOCK gpu
                 TIMEOUT 600)
  endif()
endif()

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
    # parallel matrices on one process (uploaded as sequential matrices)
    add_test(NAME AmgXSolverTestP1-gpu
             COMMAND mpirun -n 1 ./AmgXSolverTest --parallel)
    list(APPEND amgx_tests AmgXSolverTestP-gpu AmgXSolverTestP-cuda
         AmgXSolverTestP1-gpu)
  endif(MFEM_USE_MPI)
  # a deadlock fails quickly
  set_tests_properties(${amgx_tests} PROPERTIES RESOURCE_LOCK gpu TIMEOUT 300)
endif()

# HypreDevicePoolSizeTest executable: size of the pool of memory of hypre on
# the GPU, given by the option of mfem-mgis or by its environment variable
# (no memory allocated on the GPU)
if(MFEM_USE_MPI)
  add_test_executable(HypreDevicePoolSizeTest)
  # with a hypre built with Umpire, the test reads the size of the pool
  # created by hypre with a function of Umpire, which is linked explicitly
  # (MFEM may not link it, for instance when it is built without Umpire)
  find_library(MFEM_MGIS_UMPIRE_LIBRARY umpire)
  if(MFEM_MGIS_UMPIRE_LIBRARY)
    target_link_libraries(HypreDevicePoolSizeTest
      PRIVATE ${MFEM_MGIS_UMPIRE_LIBRARY})
  endif()
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

# DeviceForallTest executable: a kernel of MFEM, compiled as CUDA or HIP
# code when MFEM is compiled with CUDA or HIP. The test on the GPU requires
# one.
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

# ElementsGeometryTest executable: jacobians of the element transformations
# used by the kernels of the batch path, computed by a kernel for the meshes
# without nodes, on the host, on the debug device, on the GPU and on two
# processes
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

# MFEMFieldOperationsTest executable: copies and fills of views by kernels
# of MFEM, on the host, on the debug device and on the GPU
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

# DeviceOptionTest executable
add_test_executable(DeviceOptionTest)
add_test(NAME DeviceOptionTest-none
         COMMAND DeviceOptionTest --expected none)
add_test(NAME DeviceOptionTest-debug
         COMMAND DeviceOptionTest --device debug --expected debug)
add_test(NAME DeviceOptionTest-debug-short
         COMMAND DeviceOptionTest -dev debug --expected debug)
# a device on the host or on a GPU, with the standard memory of the host
# (MFEM compiled with Umpire allocates it with Umpire by default), unless
# the variable MFEM_MEMORY of MFEM gives it, which prevails (these tests
# fail if MFEM_MEMORY is set for the whole test suite)
add_test(NAME DeviceOptionTest-cpu
         COMMAND DeviceOptionTest --device cpu --expected cpu)
if(MFEM_USE_UMPIRE STREQUAL "YES")
  add_test(NAME DeviceOptionTest-cpu-umpire
           COMMAND DeviceOptionTest --device cpu --expected cpu
                   --expected-host-memory umpire)
  set_tests_properties(DeviceOptionTest-cpu-umpire PROPERTIES
                       ENVIRONMENT MFEM_MEMORY=umpire)
endif()
if((MFEM_USE_CUDA STREQUAL "YES") AND enable-gpu-tests)
  add_test(NAME DeviceOptionTest-cuda
           COMMAND DeviceOptionTest --device cuda --expected gpu
                   --expected-hypre device)
  # hypre on the host with a device on the GPU, requested by the option or
  # by the environment; the size of the pool of memory of hypre on the GPU
  # is then refused
  add_test(NAME DeviceOptionTest-cuda-hypre-on-host
           COMMAND DeviceOptionTest --device cuda --hypre-on-host
                   --expected gpu --expected-hypre host)
  add_test(NAME DeviceOptionTest-cuda-hypre-on-host-environment
           COMMAND DeviceOptionTest --device cuda --expected gpu
                   --expected-hypre host)
  set_tests_properties(DeviceOptionTest-cuda-hypre-on-host-environment
                       PROPERTIES ENVIRONMENT MFEM_MGIS_HYPRE_ON_HOST=1)
  add_test(NAME DeviceOptionTest-cuda-hypre-on-host-cancelled
           COMMAND DeviceOptionTest --device cuda --no-hypre-on-host
                   --expected gpu --expected-hypre device)
  set_tests_properties(DeviceOptionTest-cuda-hypre-on-host-cancelled
                       PROPERTIES ENVIRONMENT MFEM_MGIS_HYPRE_ON_HOST=1)
  add_test(NAME DeviceOptionTest-cuda-hypre-on-host-pool-size
           COMMAND DeviceOptionTest --device cuda --hypre-on-host
                   --hypre-device-pool-size 1 --expected gpu
                   --expected-error "when hypre computes on the host")
  # device configured before the initialization of hypre, by the variable
  # MFEM_DEVICE of MFEM (with the memory of the host chosen by MFEM), or
  # by the program after the initialization of mfem-mgis
  add_test(NAME DeviceOptionTest-cuda-hypre-on-host-mfem-device
           COMMAND DeviceOptionTest --hypre-on-host --expected gpu
                   --expected-host-memory any --expected-hypre host)
  set_tests_properties(DeviceOptionTest-cuda-hypre-on-host-mfem-device
                       PROPERTIES ENVIRONMENT MFEM_DEVICE=cuda)
  add_test(NAME DeviceOptionTest-cuda-hypre-on-host-program-device
           COMMAND DeviceOptionTest --hypre-on-host --create-device cuda
                   --expected gpu --expected-host-memory any
                   --expected-hypre host)
  add_test(NAME DeviceOptionTest-cuda-program-device
           COMMAND DeviceOptionTest --create-device cuda --expected gpu
                   --expected-host-memory any --expected-hypre device)
  set_tests_properties(DeviceOptionTest-cuda
                       DeviceOptionTest-cuda-hypre-on-host
                       DeviceOptionTest-cuda-hypre-on-host-environment
                       DeviceOptionTest-cuda-hypre-on-host-cancelled
                       DeviceOptionTest-cuda-hypre-on-host-pool-size
                       DeviceOptionTest-cuda-hypre-on-host-mfem-device
                       DeviceOptionTest-cuda-hypre-on-host-program-device
                       DeviceOptionTest-cuda-program-device
                       PROPERTIES RESOURCE_LOCK gpu)
endif()
# KokkosOptionTest executable: the option --kokkos and the variable
# MFEM_MGIS_KOKKOS initialize Kokkos if mfem-mgis is built with Kokkos;
# otherwise, the option is refused and the variable is ignored with a
# warning. Kokkos is finalized by mfem_mgis::finalize
add_test_executable(KokkosOptionTest)
add_test(NAME KokkosOptionTest-none
         COMMAND KokkosOptionTest --expected none)
add_test(NAME KokkosOptionTest-environment-zero
         COMMAND KokkosOptionTest --expected none)
set_tests_properties(KokkosOptionTest-environment-zero PROPERTIES
                     ENVIRONMENT MFEM_MGIS_KOKKOS=0)
add_test(NAME KokkosOptionTest-environment-overridden
         COMMAND KokkosOptionTest --no-kokkos --expected none)
set_tests_properties(KokkosOptionTest-environment-overridden PROPERTIES
                     ENVIRONMENT MFEM_MGIS_KOKKOS=1)
if(MFEM_MGIS_HAVE_KOKKOS)
  # Kokkos built for a GPU requires a GPU at its initialization
  set(_kokkos_gpu OFF)
  foreach(_kokkos_device CUDA HIP SYCL)
    if(_kokkos_device IN_LIST Kokkos_DEVICES)
      set(_kokkos_gpu ON)
    endif()
  endforeach()
  if((NOT _kokkos_gpu) OR enable-gpu-tests)
    add_test(NAME KokkosOptionTest-option
             COMMAND KokkosOptionTest --kokkos --expected initialized)
    add_test(NAME KokkosOptionTest-short
             COMMAND KokkosOptionTest -kk --expected initialized)
    add_test(NAME KokkosOptionTest-environment
             COMMAND KokkosOptionTest --expected initialized)
    set_tests_properties(KokkosOptionTest-environment PROPERTIES
                         ENVIRONMENT MFEM_MGIS_KOKKOS=1)
    if(_kokkos_gpu)
      set_tests_properties(KokkosOptionTest-option KokkosOptionTest-short
                           KokkosOptionTest-environment PROPERTIES
                           RESOURCE_LOCK gpu)
    endif()
  endif()
else()
  add_test(NAME KokkosOptionTest-option
           COMMAND KokkosOptionTest --kokkos --expected error)
  add_test(NAME KokkosOptionTest-short
           COMMAND KokkosOptionTest -kk --expected error)
  add_test(NAME KokkosOptionTest-environment
           COMMAND KokkosOptionTest --expected none)
  set_tests_properties(KokkosOptionTest-environment PROPERTIES
                       ENVIRONMENT MFEM_MGIS_KOKKOS=1)
endif()

# Install the example 1, based on UniaxialTensileTest

set(EX1_DIR "share/mfem-mgis/examples/ex1")
install(FILES UniaxialTensileTest.cxx Tolerance.hxx
              IsotropicLinearHardeningPlasticity.mfront
              cube.mesh references/IsotropicLinearHardeningPlasticity.ref
        DESTINATION ${EX1_DIR})

install(
  FILES ex1CMakeLists.txt
  DESTINATION ${EX1_DIR}
  RENAME CMakeLists.txt)

set(EX3_DIR "share/mfem-mgis/examples/ex3")
install(FILES TwoLayerCube.cxx Tolerance.hxx IsotropicLinearElasticity.mfront
              cube_2mat_per.mesh Box.med Box.per
        DESTINATION ${EX3_DIR})

install(
  FILES ex3CMakeLists.txt
  DESTINATION ${EX3_DIR}
  RENAME CMakeLists.txt)

# each process on the GPU may reserve a large memory pool (4 GiB for hypre
# with Umpire): the tests on the GPU, whose names end with -cuda, -hip or
# -gpu, are not run at the same time (the other tests on the GPU are given
# the lock where they are declared)
get_property(_tests DIRECTORY PROPERTY TESTS)
foreach(_test ${_tests})
  if(_test MATCHES "-(cuda|hip|gpu)$")
    set_tests_properties(${_test} PROPERTIES RESOURCE_LOCK gpu)
  endif()
endforeach()
