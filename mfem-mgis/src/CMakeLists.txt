set(MFEMMGIS_SOURCES)
list(APPEND MFEMMGIS_SOURCES
  Config.cxx
  ExitStatus.cxx
  MPI.cxx
  Profiler.cxx
  Parameter.cxx
  Parameters.cxx
  ParametersValidator.cxx
  MeshDiscretization.cxx
  FiniteElementDiscretization.cxx
  FiniteElementSpacesManager.cxx
  GridFunctionUtilities.cxx
  PartialQuadratureSpace.cxx
  PartialQuadratureSpaceIdentifiersManager.cxx
  PartialQuadratureFunction.cxx
  PartialQuadratureFunctionUtilities.cxx
  AbstractQPEvaluator.cxx
  AbstractQPEvaluatorGenerator.cxx
  QPEvaluatorBase.cxx
  QPEvaluators.cxx
  QPEvaluatorDescription.cxx
  QPEvaluatorsFactory.cxx
  Behaviour.cxx
  Material.cxx
  MaterialArrays.cxx
  MaterialQuantityProviderSearch.cxx
  AbstractBehaviourIntegrator.cxx
  PostProcessingsComputation.cxx
  ExecutionSettings.cxx
  BehaviourIntegratorBase.cxx
  BehaviourIntegratorFactory.cxx
  MatrixFreeTangentOperator.cxx
  MultiMaterialNonLinearIntegrator.cxx
  AbstractDirichletBoundaryCondition.cxx
  DirichletBoundaryConditionBase.cxx
  UniformDirichletBoundaryCondition.cxx
  AbstractBoundaryCondition.cxx
  UniformHeatSourceBoundaryCondition.cxx
  UniformImposedPressureBoundaryCondition.cxx
  ImposedDirichletBoundaryConditionAtClosestNode.cxx
  AbstractNonLinearEvolutionProblemPostProcessing.cxx
  NonLinearEvolutionProblemPostProcessingBase.cxx
  PostProcessingFactory.cxx
  ParaviewExportIntegrationPointResultsAtNodes.cxx
  ComputeResultantForceOnBoundary.cxx
  AbstractNonLinearEvolutionProblem.cxx
  NonLinearResolutionOutput.cxx
  NonLinearEvolutionProblemImplementationBase.cxx
  NonLinearEvolutionProblemImplementation.cxx
  NonLinearEvolutionProblem.cxx
  PeriodicNonLinearEvolutionProblem.cxx
  SolverUtilities.cxx
  LinearSolverHandler.cxx
  LinearSolverFactory.cxx
  PreconditionerFactory.cxx
  LinearSolverModules.cxx
  MFEMLinearSolvers.cxx
  HypreLinearSolvers.cxx
  PETScLinearSolvers.cxx
  AmgXSolvers.cxx
  HypreMatrixPattern.cxx
  MUMPSSolvers.cxx
  CuDSSSolvers.cxx
  AbstractNonLinearSolver.cxx
  NonLinearSolverBase.cxx
  NonLinearSolverFactory.cxx
  NewtonSolver.cxx
  AccelerationAlgorithm.cxx
  AndersonAccelerationAlgorithm.cxx
  Cast3MAccelerationAlgorithm.cxx
  SecantAccelerationAlgorithms.cxx
  FixedPointSolver.cxx
  AdaptiveLinearSolverTolerance.cxx
  AnalyticalTests.cxx
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  OrthotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  OrthotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator.cxx
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator.cxx
  IsotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator.cxx
  OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator.cxx
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator.cxx
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator.cxx
  BidimensionalMicromorphicDamageBehaviourIntegrator.cxx
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator.cxx
  TridimensionalMicromorphicDamageBehaviourIntegrator.cxx
  TransientHeatTransferBehaviourIntegrator.cxx
  Faltus2026RegularizedBehaviourIntegrators.cxx
  L2Projection.cxx
  AbstractTimeStepValidator.cxx
  TimeStepValidatorBase.cxx
  DefaultTimeStepValidator.cxx
  AbstractConvergenceFailureHandler.cxx
  ConvergenceFailureHandlerBase.cxx
  DefaultConvergenceFailureHandler.cxx
  AbstractTimeIncrementComputer.cxx
  TimeIncrementComputerBase.cxx
  DefaultTimeIncrementComputer.cxx
  AbstractSimulationMonitor.cxx
  ComputeNextStateOutput.cxx
  Provider.cxx
  Dependency.cxx
  DependenciesManager.cxx
  StateManager.cxx
  AbstractCouplingItem.cxx
  AbstractCouplingScheme.cxx
  AbstractCouplingSchemeConvergenceCriterion.cxx
  FirstIterationConvergenceCriterion.cxx
  CouplingSchemeBase.cxx
  CouplingSchemeConvergenceCriterionBase.cxx
  LoopCouplingScheme.cxx
  IterativeCouplingScheme.cxx
  CouplingSchemeConvergenceCriterionBase.cxx
  AbstractModel.cxx
  ModelBase.cxx
  NonLinearModel.cxx
  PointWiseModel.cxx
  PhysicalSystem.cxx
  Simulation.cxx
  DataFileUtilities.cxx
  AbstractPostProcessing.cxx
  PostProcessingBase.cxx
  AbstractCurve.cxx
  MultipleCurves.cxx
  GridFunctionValuesCurve.cxx
  CurvesWriter.cxx
  CurvesPostProcessing.cxx
  PointsSetCurves.cxx
  PointsSetCurvesWriter.cxx
  PointsSetCurvesPostProcessing.cxx
)

if(MGIS_HAVE_TFEL)
  list(APPEND MFEMMGIS_SOURCES
       Geometry.cxx
       FBarBehaviourIntegrators.cxx
       FBarIsotropicPlaneStrainBehaviourIntegrator.cxx
       FBarIsotropicTridimensionalBehaviourIntegrator.cxx
       FBarOrthotropicPlaneStrainBehaviourIntegrator.cxx
       FBarOrthotropicTridimensionalBehaviourIntegrator.cxx)
endif(MGIS_HAVE_TFEL)

if(MFEM_USE_GSLIB AND MGIS_HAVE_TFEL)
list(APPEND MFEMMGIS_SOURCES
     GridFunctionInterpolator.cxx)
endif(MFEM_USE_GSLIB AND MGIS_HAVE_TFEL)

if(MGIS_FUNCTION_SUPPORT AND MGIS_HAVE_TFEL)
  list(APPEND MFEMMGIS_SOURCES
       QPEvaluator.cxx
       MechanicalPostProcessings.cxx
       PartialQuadratureFunctionsSet.cxx)
endif(MGIS_FUNCTION_SUPPORT AND MGIS_HAVE_TFEL)

# batch processing (sources containing kernels of MFEM, compiled as CUDA or
# HIP code when MFEM is compiled with CUDA or HIP)
set(MFEMMGIS_DEVICE_SOURCES
  Batch/MFEMFieldOperations.cxx
  Batch/GradientsKernels.cxx
  Batch/InnerForcesKernels.cxx
  Batch/IntegrationPointsFunctionsKernels.cxx
  Batch/StiffnessMatricesKernels.cxx
  Batch/FBarKernels.cxx
  Batch/FaltusKernels.cxx
  Batch/ElementsGeometryKernels.cxx)
# scaling and conversions of the vectors exchanged with AmgX (single
# precision)
if(MFEM_MGIS_HAVE_AMGX)
  list(APPEND MFEMMGIS_DEVICE_SOURCES AmgXKernels.cxx)
endif()
list(APPEND MFEMMGIS_SOURCES ${MFEMMGIS_DEVICE_SOURCES})
# host part of the computations by kernels (geometry of the elements,
# gradients, inner forces, values at the integration points, stiffness
# matrices, FBar, Faltus et al., global stiffness matrix)
list(APPEND MFEMMGIS_SOURCES
  Batch/ElementsGeometry.cxx
  Batch/Gradients.cxx
  Batch/InnerForces.cxx
  Batch/IntegrationPointsFunctions.cxx
  Batch/StiffnessMatrices.cxx
  Batch/FBar.cxx
  Batch/Faltus.cxx
  Batch/GlobalStiffnessMatrix.cxx)

mfem_mgis_library(MFEMMGIS ${MFEMMGIS_SOURCES})
mfem_mgis_compile_as_device_code(MFEMMGIS ${MFEMMGIS_DEVICE_SOURCES})
# in single precision, the stiffness matrices are computed on the device
# without contraction into fused multiply-add operations. With the
# contractions of nvcc, the Newton algorithm diverged, or converged to less
# accurate solutions, on a refined beam in plasticity, probably by drifting
# along its soft modes (bending). The kernels of the FBar formulation and of
# the regularization of Faltus et al., which compute stiffness matrices in
# the same way, are compiled likewise. In double precision, the contractions
# are kept: they did not change the convergence, and the assembly of the
# beam was 13 to 16 % slower without them.
if(MFEM_MGIS_REAL_TYPE STREQUAL "float")
  mfem_mgis_disable_device_floating_point_contraction(MFEMMGIS
    Batch/StiffnessMatricesKernels.cxx
    Batch/FBarKernels.cxx
    Batch/FaltusKernels.cxx)
endif()
# the kernels of the library use the runtime of CUDA or HIP: the dependency
# is carried by the library itself, not only by the target of MFEM, which a
# project may define on its own
if(MFEM_USE_CUDA STREQUAL "YES")
  target_link_libraries(MFEMMGIS PUBLIC CUDA::cudart)
endif()
if(MFEM_USE_HIP STREQUAL "YES")
  target_link_libraries(MFEMMGIS PUBLIC hip::host)
endif()
# AmgX precedes MFEM in the link line: a MFEM compiled with AmgX contains
# its own copy of AmgX (static library), which the library of AmgX given to
# mfem-mgis replaces
if(MFEM_MGIS_HAVE_AMGX)
  target_compile_definitions(MFEMMGIS PRIVATE MFEM_MGIS_HAVE_AMGX)
  foreach(_property LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
    get_target_property(_libraries MFEMMGIS ${_property})
    set_property(TARGET MFEMMGIS PROPERTY ${_property}
                 AMGX::amgxsh ${_libraries})
  endforeach(_property)
endif()
# Kokkos: the only source including Kokkos is compiled by the compiler of
# Kokkos in a shared library privately linked to Kokkos, so that neither the
# other sources nor the projects linked to mfem-mgis inherit the flags of
# Kokkos
if(MFEM_MGIS_HAVE_KOKKOS)
  add_library(MFEMMGISKokkos SHARED KokkosRuntime.cxx)
  target_link_libraries(MFEMMGISKokkos PRIVATE Kokkos::kokkos)
  # compilation and link redirected to the compiler of Kokkos by its
  # launcher, which does not expect another launcher (such as ccache)
  kokkos_compilation(TARGET MFEMMGISKokkos)
  # the installed Kokkos may not support the extensions of the compilers;
  # the symbols are exported despite the flag -fvisibility=hidden of the
  # compilation; the installed library finds the one of Kokkos
  set_target_properties(MFEMMGISKokkos PROPERTIES
    CXX_COMPILER_LAUNCHER ""
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET default
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    INSTALL_RPATH_USE_LINK_PATH ON)
  install(TARGETS MFEMMGISKokkos EXPORT MFEMMGIS
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib${LIB_SUFFIX}
          ARCHIVE DESTINATION lib${LIB_SUFFIX})
  target_link_libraries(MFEMMGIS PRIVATE MFEMMGISKokkos)
  target_compile_definitions(MFEMMGIS PRIVATE MFEM_MGIS_HAVE_KOKKOS)
endif()
if(MFEM_USE_GSLIB AND MGIS_HAVE_TFEL)
 target_compile_definitions(MFEMMGIS PUBLIC 
   MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR)
endif(MFEM_USE_GSLIB AND MGIS_HAVE_TFEL)

target_compile_definitions(MFEMMGIS PUBLIC 
  $<$<CONFIG:Debug>:MFEM_MGIS_DEBUG>)

if(OpenMP_FOUND)
  target_link_libraries(MFEMMGIS
    PUBLIC OpenMP::OpenMP_CXX)
endif(OpenMP_FOUND)

