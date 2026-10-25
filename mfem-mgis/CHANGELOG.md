# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- update of partial quadrature functions and grid functions defined on
  different meshes (whole mesh, submeshes of materials or of boundaries),
  see `MeshDiscretization::getElementsMapping` and `getInverseMapping`
- copy by `assign_values` between partial quadrature functions defined on
  the same location of different meshes
- `PartialQuadratureSpace::getLocation`
- `FiniteElementSpaceOrdering` parameter, to order the degrees of freedom by
  vector dimension (`byVDIM`) instead of by nodes (`byNODES`, the default). The
  `Elasticity` strategy of the `HypreBoomerAMG` preconditioner requires
  `byVDIM`.
- `isMPIMainProcess`, `getMPIRank` and `getMPISize` for a communicator or a
  mesh discretization, and the same functions for the default communicator
- `synchronize` gathers the error messages of a context on all the
  processes of a communicator or of a mesh discretization
- `IntegrationRuleOrder` parameter of the finite element discretization, to
  choose the order of the integration rules of the behaviour integrators (8
  points instead of 27 on the hexahedra of order 1 with the order 3).
- execution settings: option `--execution cpu|gpu|mixed` and, in the mixed
  mode, `--finite-elements-execution`, `--behaviours-execution`,
  `--linear-solvers-execution` and `--post-processings-execution` placing
  each phase on the CPU or on the GPU
  (`getExecutionSettings`); behaviour libraries and linear solvers checked
  against them; parameter `BehaviourExecutionLocation` of
  `addBehaviourIntegrator` giving the location expected for the behaviour
  of a material.
- CMake function `mfem_mgis_add_behaviours_library` (module
  `MFEMMGISBehaviours.cmake`): behaviour libraries for the CPU (generic
  interface) and for the GPU of MFEM (generic-parallel interface), chosen at
  run time by `getBehavioursLibrary` from the execution settings
- `getMaterialArraysTransfers`: copies of the arrays of the materials
  between the host and the device, printed after the profiling table;
  `MaterialArrays::read` and `MaterialArrays::readTangentOperatorBlocks`
  (the latter also available in `BehaviourIntegratorBase`).
- integration of the behaviours by batches and computations by kernels of
  the gradients, of the inner forces and of the tangent operator, on the
  host or on a GPU: parameters `BehaviourIntegration`,
  `GradientsComputation`, `InnerForcesComputation`,
  `TangentOperatorComputation` and `TangentOperatorAssembly` (`Atomic` or
  `Ordered`) of the non linear evolution problems, and environment
  variables `MFEM_MGIS_BEHAVIOUR_INTEGRATION`, `MFEM_MGIS_*_COMPUTATION`
  and `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`; FBar and Faltus
  regularizations computed by kernels
- tangent operator applied without assembly (`TangentOperatorComputation`
  set to `MatrixFree`) and preconditioner `AssembledOperator`
- post-processings computed by kernels (environment variable
  `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION`)
- `FixedPoint` non linear solver (methods `ElasticOperator`,
  `ModifiedNewton` and `NewtonK0`) and its acceleration algorithms
  (`Anderson`, `UAnderson`, `FAnderson`, secant methods, `IronsTuck`,
  `Steffensen`, `Cast3M`...)
- `AdaptiveLinearSolverTolerance` option (Eisenstat-Walker) of the Newton
  solver and of the `NewtonK0` method
- solvers of AmgX (`AmgX`, `AmgXPCG`, CMake option `enable-amgx`) and of
  cuDSS (`CuDSSSolver`); parameter `ReuseAnalysis` of `MUMPSSolver`;
  `AbsoluteTolerance` of the Krylov solvers of hypre; options of
  `HypreBoomerAMG` (`StrongThreshold`, `CoarseningType`,
  `AggressiveCoarseningLevels`, `InterpolationType`, `RelaxationType`...);
  preconditioners `OperatorJacobiSmoother` and `OperatorChebyshevSmoother`
- options `--device`, `--hypre-device-pool-size`, `--hypre-on-host` and
  `--kokkos` read by `mfem_mgis::initialize` (environment variables
  `MFEM_MGIS_DEVICE`, `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`,
  `MFEM_MGIS_HYPRE_ON_HOST` and `MFEM_MGIS_KOKKOS`);
  CMake options `enable-gpu-tests` and `enable-kokkos`
- option `--gpu-aware-mpi` read by `mfem_mgis::initialize` (environment
  variable `MFEM_MGIS_GPU_AWARE_MPI`): MFEM, hypre and AmgX give the arrays
  on the GPU to MPI, without copies on the host
- CMake module `MFEMMGISDeviceCode.cmake`, compiling the sources containing
  kernels of MFEM as CUDA or HIP code; MFEM compiled with HIP is supported
  (compiled and linked, not run on an AMD GPU yet)

### Changed

- use std::span instead of mgis::span
- the arrays of the materials are stored in the memory of MFEM: after a
  computation on a GPU, a direct access to `s0`, `s1` or `K` on the host
  requires `synchronizeOnHost`; `Material` privately inherits from
  `MaterialArrays` and its constructor takes a layout (`mgis::FieldLayout`)
- partial quadrature spaces are associated with a location (material or
  boundary). Partial quadrature spaces defined on a boundary using a finite
  element discretization are built on the boundary elements of the mesh
- `QPDependency`, `PartialQuadratureFunctionsSet` and
  `PartialQuadratureSpace::treatInvalidElementIndex` use locations instead
  of material identifiers
- `LocationIdentifier` is no longer an aggregate and can not be assigned
- `PartialQuadratureSpace::getLocationName` returns a description of the
  location, e.g. `material 'BEAM' (12)`
- the static `buildQuadratureSpace` methods of the behaviour integrators take
  `attributes::Throwing` as first argument
- `setMaterialProperty` and `setExternalStateVariable` called with an
  evaluator require an equivalent partial quadrature space and are collective
- `updateL2Projection` reports an error for functions defined outside the
  mesh of the result
- `updateGridFunction` requires a grid function whose finite element space is
  managed by the finite element spaces manager of the partial quadrature
  functions
- `updateGridFunction` takes a `Context` as first argument and returns a
  boolean
- `assign_values` requires both views to be defined on the same mesh
  discretization
- partial quadrature spaces are handled by value: copies share the same
  data and compare equal. Constructors and functions taking a
  `std::shared_ptr<const PartialQuadratureSpace>` take a `const
  PartialQuadratureSpace&`, and the `buildQuadratureSpace` methods of
  the behaviour integrators return a `PartialQuadratureSpace`
- partial quadrature function views can no longer be default constructed nor
  assigned
- `getPartialQuadratureSpacePointer`, replaced by `getPartialQuadratureSpace`

### Fixed

- the message of a failed time step described its first sub-step instead
  of the failed time step
- with Umpire, the processes of a node allocated the device memory of MFEM
  on the GPU 0, whatever their GPU (MFEM 4.10, mfem/mfem#5696)

### Removed

- deprecated methods of `AbstractNonLinearEvolutionProblem` and deprecated
  functions `getGradient`, `getThermodynamicForce` and
  `getInternalStateVariable` which did not take a `Context`
- `PartialQuadratureSpace::getId` and the constructors of
  `PartialQuadratureSpace` taking a material identifier
- the overloads of `updateGridFunction` taking a mesh or a submesh
- `QPDependency::getMaterialIdentifier`,
  `PartialQuadratureFunctionsSet::getMaterialIdentifiers` and the overload of
  `PartialQuadratureFunctionsSet::get` taking a material identifier
- deprecated functions `getMPIrank` and `getMPIsize`

## [1.0.0] - 2024-08-02

### Added

- Support for multiple materials with distinct constitutive equations:
  - Micromorphic laws
  - Orthotropic laws
- Support for internal state variables (defined at integration points):
  - Use of MGIS data structures
  - Use of mechanical behavior laws generated by MFront
- Post-processing features for nonlinear mechanics:
  - `ComputeResultantForceOnBoundary`: Compute the resultant of the inner forces on a boundary
  - `ComputeMeanThermodynamicForcesValues`: Compute macroscopic stress and strain for each material
  - `ParaviewExportIntegrationPointResultsAtNodes`: Generate Paraview post-processing files with internal state variables using MGIS
- Flexible C++ API for defining MFEM simulation features
- CI workflow with examples in `test` and 58% coverage
- Profiler for performance analysis
- Support for MFEM solvers:
  - Direct solvers: 
    - MUMPS
    - UMFPackSolver
  - Krylov solvers:
    - CGSolver
    - GMRESSolver
    - BiCGSTABSolver
    - MINRESSolver
    - SLISolver
		- HyprePCG
    - HypreGMRES
		- HypreFGMRES
    - Solver from Petsc
  - Preconditioner:
    - HypreILU
    - HypreEuclid
    - HypreBoomerAMG
    - HypreParaSails
    - HypreDiagScale
    - Preconditioner from Petsc
