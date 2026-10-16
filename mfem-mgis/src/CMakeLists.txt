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
  SolverOverrides.cxx
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
  HypreSolvers.cxx
  NodesCoordinates.cxx
  PETScSupport.cxx
  PETScLinearSolvers.cxx
  SuperLUSolvers.cxx
  STRUMPACKSolvers.cxx
  MKLSparseSolvers.cxx
  MAGMASparseSolvers.cxx
  NodalBlockJacobiPreconditioner.cxx
  GinkgoLinearSolvers.cxx
  AmgXSolvers.cxx
  TrilinosLinearSolvers.cxx
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
  KINSOLSolver.cxx
  PETScSNESSolver.cxx
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
  Batch/ElementsGeometryKernels.cxx
  NodalBlockJacobiKernels.cxx)
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

# the solver overrides are read by the utilities of TFEL (found
# before the library, whose installed configuration finds them for the
# projects using the static library)
set(MFEM_MGIS_HAVE_TFEL_UTILITIES OFF)
if(MGIS_HAVE_TFEL)
  if(NOT TARGET tfel::TFELUtilities)
    find_package(TFELUtilities QUIET
                 HINTS "${TFEL_DIR}" "${_TFEL_MODULE_PATH}")
  endif()
  if(TARGET tfel::TFELUtilities)
    set(MFEM_MGIS_HAVE_TFEL_UTILITIES ON)
  endif()
endif()
message(STATUS "mfem-mgis solver overrides (utilities of TFEL): "
        "${MFEM_MGIS_HAVE_TFEL_UTILITIES}")
# for the tests
set(MFEM_MGIS_HAVE_TFEL_UTILITIES ${MFEM_MGIS_HAVE_TFEL_UTILITIES}
    PARENT_SCOPE)
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
# optional libraries of MFEM (see cmake/modules/MFEMMGISOptionalLibraries.cmake):
# public macros, also read by the headers and by the programs using
# mfem-mgis
foreach(_library ${MFEM_MGIS_OPTIONAL_MFEM_LIBRARIES})
  string(REPLACE ":" ";" _fields "${_library}")
  list(GET _fields 1 _variable)
  if(MFEM_MGIS_HAVE_${_variable})
    target_compile_definitions(MFEMMGIS PUBLIC MFEM_MGIS_HAVE_${_variable})
  endif()
endforeach(_library)
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
# sparse part of MAGMA, whose interface is only used by the library
if(MFEM_MGIS_HAVE_MAGMA_SPARSE)
  target_compile_definitions(MFEMMGIS PRIVATE MFEM_MGIS_HAVE_MAGMA_SPARSE)
  target_include_directories(MFEMMGIS PRIVATE
                             ${MFEM_MGIS_MAGMA_SPARSE_INCLUDE_DIR})
  target_link_libraries(MFEMMGIS PRIVATE ${MFEM_MGIS_MAGMA_SPARSE_LIBRARY})
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
# Trilinos: the sources including Trilinos are compiled by the compiler of
# Kokkos in a shared library privately linked to Trilinos (as the library
# MFEMMGISKokkos), whose interface only uses standard types (see
# TrilinosBackend.hxx)
if(MFEM_MGIS_HAVE_TRILINOS)
  add_library(MFEMMGISTrilinos SHARED TrilinosBackend.cxx
              TrilinosTpetraBackendDefaultNode.cxx
              TrilinosTpetraBackendHostNode.cxx)
  if(MFEM_MGIS_HAVE_TEKO)
    # block preconditioners of Teko (option enable-teko)
    target_sources(MFEMMGISTrilinos PRIVATE TrilinosTekoBackend.cxx)
  endif()
  target_link_libraries(MFEMMGISTrilinos PRIVATE
                        Trilinos::all_selected_libs MPI::MPI_CXX)
  kokkos_compilation(TARGET MFEMMGISTrilinos)
  if(Kokkos_ENABLE_CUDA AND (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
    # the host code generated by nvcc writes the functional casts of the
    # sources and of the headers of Trilinos as old-style casts
    target_compile_options(MFEMMGISTrilinos PRIVATE -Wno-old-style-cast)
  endif()
  set_target_properties(MFEMMGISTrilinos PROPERTIES
    CXX_COMPILER_LAUNCHER ""
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET default
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    INSTALL_RPATH_USE_LINK_PATH ON)
  install(TARGETS MFEMMGISTrilinos EXPORT MFEMMGIS
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib${LIB_SUFFIX}
          ARCHIVE DESTINATION lib${LIB_SUFFIX})
  target_link_libraries(MFEMMGIS PRIVATE MFEMMGISTrilinos)
  target_compile_definitions(MFEMMGIS PRIVATE MFEM_MGIS_HAVE_TRILINOS)
  # numerical type of the interface between mfem-mgis and the library
  # MFEMMGISTrilinos (see TrilinosBackend.hxx)
  foreach(_target MFEMMGIS MFEMMGISTrilinos)
    target_compile_definitions(${_target} PRIVATE
      MFEM_MGIS_TRILINOS_REAL_TYPE=${MFEM_MGIS_REAL_TYPE})
    if(MFEM_MGIS_HAVE_TEKO)
      target_compile_definitions(${_target} PRIVATE MFEM_MGIS_HAVE_TEKO)
    endif()
  endforeach()
endif()
if(MFEM_MGIS_HAVE_TFEL_UTILITIES)
  target_link_libraries(MFEMMGIS PRIVATE tfel::TFELUtilities)
  target_compile_definitions(MFEMMGIS PRIVATE MFEM_MGIS_HAVE_TFEL_UTILITIES)
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

