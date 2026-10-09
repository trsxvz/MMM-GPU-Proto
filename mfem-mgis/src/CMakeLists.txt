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
  BehaviourIntegratorBase.cxx
  BehaviourIntegratorFactory.cxx
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
  AbstractNonLinearSolver.cxx
  NonLinearSolverBase.cxx
  NonLinearSolverFactory.cxx
  NewtonSolver.cxx
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

# batch processing (sources containing kernels of MFEM, compiled as CUDA
# code when MFEM is compiled with CUDA). The compilation of kernels of MFEM
# compiled with HIP is not supported yet: those sources are then left out.
set(MFEMMGIS_DEVICE_SOURCES
  Batch/MFEMFieldOperations.cxx
  Batch/GradientsKernels.cxx
  Batch/InnerForcesKernels.cxx
  Batch/StiffnessMatricesKernels.cxx)
if(MFEM_USE_HIP STREQUAL "YES")
  message(STATUS "mfem-mgis: MFEM compiled with HIP, the sources containing "
          "kernels of MFEM are not compiled (${MFEMMGIS_DEVICE_SOURCES})")
  set(MFEMMGIS_DEVICE_SOURCES)
endif()
list(APPEND MFEMMGIS_SOURCES ${MFEMMGIS_DEVICE_SOURCES})
# host part of the computations by kernels (gradients, inner forces), which
# reports an error if the kernels are not compiled
list(APPEND MFEMMGIS_SOURCES
  Batch/ElementsGeometry.cxx
  Batch/Gradients.cxx
  Batch/InnerForces.cxx
  Batch/StiffnessMatrices.cxx)

mfem_mgis_library(MFEMMGIS ${MFEMMGIS_SOURCES})
if(MFEMMGIS_DEVICE_SOURCES)
  mfem_mgis_compile_as_device_code(MFEMMGIS ${MFEMMGIS_DEVICE_SOURCES})
  target_compile_definitions(MFEMMGIS PRIVATE MFEM_MGIS_HAVE_DEVICE_SOURCES)
endif()
# the kernels of the library use the runtime of CUDA: the dependency is
# carried by the library itself, not only by the target of MFEM, which a
# project may define on its own
if(MFEM_USE_CUDA STREQUAL "YES")
  target_link_libraries(MFEMMGIS PUBLIC CUDA::cudart)
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

