# Detailed plan for mfem-mgis

Status: plan drawn up on 09/10/2026, reviewed independently the same day,
then checked against the whole code of mfem-mgis (82 features
inventoried; omissions and partial coverages integrated). Progress on
the evening of 09/10/2026: L0, L1 and L2 done and reviewed (see the
[log](../JOURNAL.md)). L3 advanced: gradients, residual and assembled
tangent by kernels, through the entry points of
`NonLinearEvolutionProblemImplementation` (`Mult`, `GetGradient`), action
and diagonal of the tangent without assembly (`MatrixFree` option, with a
Jacobi preconditioner only); remaining: the nodal blocks, the complete
dedicated operator and the partial assembly of MFEM. Minimal slice of L6
done (`cuda` behaviour integrated on the GPU arrays, end of step on the
GPU, validation by the `debug` device, profile without hidden copy),
prediction of L4 by kernels done. On a GPU: beam r4 in 3.3 s after the
end of step on the GPU (log of 09/10; 4.06 s before,
`results/L6-beam-2026-10-09`) against 17.5 s on 8 cores, polycrystal in
13.2 s against 31.7 s (`results/L6-polycrystal-2026-10-09`). On
10/10/2026: action of the tangent without assembly by a fused kernel and
derivatives of the unknowns by a dedicated kernel (`results/L4-*`),
preconditioner built on the assembled tangent (`AssembledOperator`); L5
wave 1 largely done: `FixedPoint` solver (methods `ElasticOperator` and
`ModifiedNewton`, all the accelerations of MTest, prediction by the AMG
of K0), options of BoomerAMG; on 8 cores the fixed point beats Newton
(×1.5 to ×1.7), on the GPU a tie with the Chebyshev smoother
(`results/L5-*`). `NewtonK0` done (method of the `FixedPoint` solver, any
iterative solver of MFEM from the factory preconditioned by the K0
solver): on the GPU, ×1.2 on the beam and ×1.06 on the polycrystal, on 8
cores ×1.19 on the polycrystal, same iterations as Newton. Inexact Newton
done (`AdaptiveLinearSolverTolerance` option, Eisenstat-Walker): with
NewtonK0, ×1.59 against Newton on the GPU (beam). Remaining in L5: mixed
precision, line search, matrix free smoothers, waves 2 and 3. On
13/10/2026: work packages L0 to L10 done (status in
[overview.md](overview.md#progress)); matrix free smoothers done
(`OperatorJacobiSmoother`, `OperatorChebyshevSmoother`), as well as
cuDSS and AmgX (own wrapper) of wave 2. Remaining: mixed precision (task
9b, [plan-precision.md](plan-precision.md)), line search, Jacobi by nodal
blocks, Ginkgo and the rest of waves 2 and 3. Tasks done next on 13/10:
9a, execution settings ([plan-execution.md](plan-execution.md)), and 9c,
kernels of mfem-mgis compiled for HIP. The fundamental choices are in
[plan-mfem-mgis.md](plan-mfem-mgis.md); the list of solvers in
[gpu-solvers.md](gpu-solvers.md). This document splits the work into
work packages, says which files change, what is tested and when a work
package is finished.

Base: mfem-mgis master `a011d97d6f` (06/10/2026), repository
`mfem-mgis/`. The `file:line` references refer to this base. On
13/10/2026, the branch of the prototype was rebased onto master
`8fac67336e`, which contains the rework of `PartialQuadratureSpace`
(point 7 below). TFEL and MGIS are those of the prototype (MGIS installed
in `install/mgis` by `build-mgis.sh`), MFEM 4.10 with CUDA (and the
OpenMP and libCEED backends) comes from the `mmm-gpu` environment. The
plan applies to mfem-mgis, to mfem-mgis-examples and to mm-opera-hpc.

## Contents

1. [The code today](#the-code-today)
2. [Principles](#principles)
3. [Proposed parameters](#proposed-parameters)
4. [Overview of the work packages](#overview-of-the-work-packages)
5. [L0. Working base](#l0-working-base)
6. [L1. Data at the points in the memory of MFEM](#l1-data-at-the-points-in-the-memory-of-mfem)
7. [L2. Behaviours integrated by batches on CPU](#l2-behaviours-integrated-by-batches-on-cpu)
8. [L3. Batch kernels and dedicated operator](#l3-batch-kernels-and-dedicated-operator)
9. [L4. Tangent operator and prediction](#l4-tangent-operator-and-prediction)
10. [L5. Solvers](#l5-solvers)
11. [L6. GPU](#l6-gpu)
12. [L7. All of mfem-mgis](#l7-all-of-mfem-mgis)
13. [L8. Outputs and post-processings](#l8-outputs-and-post-processings)
14. [L9. Examples and measurements](#l9-examples-and-measurements)
15. [L10. Documentation and upstream proposals](#l10-documentation-and-upstream-proposals)
16. [Test strategy](#test-strategy)
17. [Measurement campaign](#measurement-campaign)
18. [Questions to settle](#questions-to-settle)
19. [Risks](#risks)

## The code today

What the batch path must replace or work around, as found in the code
(files of `src/` and `include/MFEMMGIS/`).

### Solution loop

- `NonLinearEvolutionProblemImplementation<parallel>` inherits from
  `(Par)NonlinearForm` and from `NonLinearEvolutionProblemImplementationBase`.
  It holds a `MultiMaterialNonLinearIntegrator* const` (`final` class),
  added to the form by `AddDomainIntegrator`
  (`NonLinearEvolutionProblemImplementation.cxx:335` and `:557`).
- A Newton iteration makes **three element by element passes on the
  host**:
  1. `integrate(u, it, dt)` (`...Implementation.cxx:485-514` and
     `:715-745`): for each element, element values, then integration of
     the behaviour point by point;
  2. `Mult` of the form: `AssembleElementVector`, which reads the
     thermodynamic forces;
  3. `GetGradient` of the form: `AssembleElementGrad`, which reads `K`.
- **`Mult` does not integrate the behaviour.** The integration is
  triggered by an action of the nonlinear solver of mfem-mgis, called at
  each new estimate (`NonLinearSolverBase.cxx:39-58`). An external
  nonlinear solver (MFEM, SUNDIALS, PETSc), which only calls `Mult`,
  would compute the residual with outdated stresses. The only exception
  is the PETSc mode: `AssembleElementVector` integrates the behaviour,
  and a failure of the behaviour stops the program (`MFEM_VERIFY`,
  `MultiMaterialNonLinearIntegrator.cxx:290-297`).
- The residual goes through `computeResidual`
  (`NonLinearSolverBase.cxx:135-140`), the tangent through `getJacobian`
  (`:142-152`), that is through `oper->Mult` and `oper->GetGradient`: a
  dedicated operator can be plugged in there.
- The Newton (`mfem_mgis::NewtonSolver`) is specific to mfem-mgis: it
  derives from `mfem::IterativeSolver` through `NonLinearSolverBase`, not
  from `mfem::NewtonSolver`. It divides the correction by 2 when the
  behaviour fails (`NewtonSolver.cxx:166-195`), without any other line
  search.
- The linear solver receives the tangent at each iteration
  (`NewtonSolver.cxx:245`, `prec` is the linear solver there); its
  `SetOperator` rebuilds the preconditioner (Krylov of MFEM and of
  hypre).
- The only nonlinear solver of the factory is `Newton`
  (`NonLinearSolverFactory.cxx`), plus `PetscNonlinearSolver` when PETSc
  is requested.
- The prediction (`...Implementation.cxx:98-233`) goes through a
  `BilinearForm` and a `LinearForm` (`StiffnessMatrixIntegrator`,
  `InternalForcesIntegrator`), element by element. It is set in C++
  only (`setPredictionPolicy`), along four axes: strategy (including
  `CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`), operator, integration
  operator, zero time increment
  (`AbstractNonLinearEvolutionProblem.hxx:39-205`).
- The Dirichlet conditions write the imposed values into `u` on the host
  (`UniformDirichletBoundaryCondition.cxx:64-79`).
- End of step: `update` and `revert` of MGIS on the host, material by
  material (`BehaviourIntegratorBase.cxx:779-785`).
- Profiling: measurement sections already exist (`CatchTimeSection`,
  `ProfilingSection`, "NS::..." sections).

### Integrators, materials, data at the points

- A behaviour integrator **is** the material: `BehaviourIntegratorBase`
  derives from `Material`, which derives from
  `mgis::behaviour::MaterialDataManager`.
- An element attribute can carry **several** integrators, whose
  contributions are summed (`MultiMaterialNonLinearIntegrator.cxx:272-276`
  and `:301-318`).
- `Material` builds the MGIS manager without an initializer
  (`Material.cxx:48-61`): arrays allocated by MGIS on the host, laid out
  point by point.
- The integration is done point by point
  (`BehaviourIntegratorBase::performsLocalBehaviourIntegration`,
  `BehaviourIntegratorBase.cxx:672-743`), with a `BehaviourDataView` and
  shared buffers. It requires a `generic` library.
- The real code of the integrators is in
  `StandardBehaviourIntegratorCRTPBase.ixx`: gradients (`:41-142`),
  internal forces (`:154-219`), tangent (`:222-290`), with
  `CalcPhysDShape` and `DenseMatrix` members. The `B` matrices of the
  concrete classes are **generated** (`code-generation/behaviour-integrator.cxx`,
  GiNaC). Element indexing `u[ni + c * nnodes]`, native order of the
  degrees of freedom of MFEM.
- **Integration rule**: `IntRules.Get(geom, 2 * t.OrderGrad(&e))`
  (generated classes, for example
  `IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.cxx:16-20`),
  that is 27 points on the hexahedron of order 1. It is chosen element by
  element: a mixed mesh has several rules.
- **Rotations of orthotropic behaviours**: the gradients are rotated
  before the integration and the blocks of the tangent after; the
  thermodynamic forces stay in the material frame (they serve as the
  state at the next step) and are only rotated at assembly, in a copy
  (`StandardBehaviourIntegratorCRTPBase.ixx:134-139` and `:208-210`).
  Matrices per point given by a `variant` (`Material.hxx:99-106`),
  rotations by function pointers of MGIS.
- `PartialQuadratureSpace`: contiguous points per element, elements of a
  material in increasing order, but offsets stored in a
  `std::unordered_map` (`PartialQuadratureSpace.hxx:230-236`).
- `PartialQuadratureFunction`: `std::vector` or borrowed memory,
  `std::span` views, access `values[o * stride + begin]`.
  `getIntegrationPointValues` returns a contiguous `std::span`
  (`PartialQuadratureFunction.hxx:198-228` and `:323-352`), called about
  forty times; the integrators read the data through `o * stride`. **The
  whole current path assumes the point by point layout.**
- Material properties and external variables: evaluators
  (`AbstractQPEvaluator`), evaluated on the host in `setup`
  (`BehaviourIntegratorBase.cxx:360-670`).

### Linear solvers

- Factory `LinearSolverFactory.cxx`: Krylov of MFEM (sequential and
  parallel), `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `MUMPSSolver`,
  `UMFPackSolver`; preconditioners `HypreBoomerAMG` (strategies
  `Elasticity`, `System`, `None`; `System` by default when the space has
  several components), `HypreEuclid`, `HypreILU`, `HypreParaSails`,
  `HypreDiagScale`. All the preconditioners and direct solvers require an
  assembled matrix.
- `SystemHypreBoomerAMG` calls `SetSystemsOptions` again at each
  `SetOperator` (`LinearSolverFactory.cxx:46-50`), which resets two
  settings of the AMG to their default value.
- Two defects noticed along the way: `HypreGMRES` tests an
  `AbsoluteTolerance` parameter that is not in the allowed list (dead
  code), `HypreFGMRES` too, and calls `SetTol` instead of an absolute
  tolerance (`mfem::HypreFGMRES` has no `SetAbsTol` in MFEM 4.9 or 4.10:
  go through `HYPRE_ParCSRFlexGMRESSetAbsoluteTol`).

### GPU

- No use of `mfem::Device`, `AssemblyLevel`, `forall`,
  `ElementRestriction`, `QuadratureInterpolator`, `GeometricFactors` or
  `DofToQuad`. `mfem_mgis::initialize` does not create a device.
- `mfem::forall` only launches a CUDA kernel if the file is compiled by
  nvcc (`__CUDACC__`, `general/forall.hpp:789-794`); otherwise the body
  runs on the host. The files of mfem-mgis that contain kernels must
  therefore be compiled as CUDA when MFEM is.

### Extension points and features to preserve

Found by checking the plan against the whole code:

- **Integrators of the form other than the one of mfem-mgis.** The
  problem inherits publicly from `(Par)NonlinearForm`: the user can add
  their own integrators to it. mfem-mgis itself puts there the imposed
  pressure (boundary integrator) and the heat source (domain integrator
  with markers, `UniformHeatSourceBoundaryCondition.cxx:208` and `:223`);
  the ex8 example of mfem-mgis-examples puts there a Robin condition
  (boundary face integrator whose tangent depends on `u`,
  `ex8/src/RobinBC.cxx`). The
  `UseMultiMaterialNonLinearIntegrator=false` mode
  (`NonLinearEvolutionProblemImplementationBase.cxx:18-42`) does without
  the integrator of mfem-mgis.
- **Imposed degrees of freedom**: those of the form
  (`SetEssentialTrueDofs`), not only those of the Dirichlet conditions.
  The periodic problem sets them directly
  (`PeriodicNonLinearEvolutionProblem.cxx:100`, `:192`, `:235`, `:291`):
  `getEssentialDegreesOfFreedom()` is empty there. The interface of the
  Dirichlet conditions is public and writes `u` on the host.
- **Data of the materials read or written on the host by the public
  API**: `getGradient`, `getThermodynamicForce`,
  `getInternalStateVariable` (`Material.cxx:269-370`),
  `MechanicalPostProcessings`, `makeGridFunction`, `computeL2Projection`,
  `StateManager`, `makeGradientEvaluator`, post-processings given by
  `std::function`, initialization of the internal variables by the user.
  All these views rely on a stride between points (`data_stride`).
- **Material properties and external variables set by the user in their
  own host memory**, through the API of MGIS (tests
  `StationaryNonLinearHeatTransferTest`, `TwoLayerCube`; mm-opera-hpc
  `cermet/cermet.cxx:118-135`).
- **Convergence of the linear solvers**: `hasConverged` and
  `getNumberOfIterationsAtConvergence` (`SolverUtilities.cxx:103-163`)
  declare converged any solver they do not know.
- **Failures**: an attribute without an integrator stops the computation
  today (`MultiMaterialNonLinearIntegrator.cxx:30-53`); the integration
  stops at the first failing element; status 0 ("unreliable results")
  counts as a success.
- **Prediction**: it adds the contribution of the boundary conditions
  (`addLinearFormIntegrators`, `...Implementation.cxx:160`) and respects
  `DiscardLinearSolverFailure`. `getLinearizedOperators` is a public API
  (used by `PredictionTest`).
- **Sequential**: the assembled tangent is a `SparseMatrix` (UMFPACK,
  Krylov of MFEM); the hypre preconditioners are refused there
  (`LinearSolverFactory.cxx:245-300`).
- **Dependencies between problems**:
  `resolveBehaviourIntegratorsDependencies`
  (`NonLinearEvolutionProblem.cxx:795-854`) links the integrators of two
  problems, which must keep the same integration rule.
- **Implicit gradient regularization** (`L2Projection.cxx:624-702`) and
  integrals (`computeMeasure`, `computeIntegral`, energies, resultant):
  host loops over the points.
- **What needs no work**: the evaluator factory (almost entirely
  commented out), the dependencies between models (not implemented), the
  restart of the state of the materials and the Python bindings (absent),
  the reading of the meshes, the time management.

### What the partial assembly of MFEM 4.10 imposes

Checked in the sources of MFEM 4.9 and 4.10 (`fem/nonlinearform.cpp`,
`fem/nonlinearform_ext.cpp`, `fem/pnonlinearform.cpp`):

- `NonlinearForm::Mult` returns right after the extension: the boundary
  integrators (`AddBoundaryIntegrator`) are **ignored without a message**
  (imposed pressure, `UniformImposedPressureBoundaryCondition`). A domain
  integrator with markers without partial assembly stops the computation
  (`AssemblePA` not implemented): this is the case of the heat source
  (`UniformHeatSourceBoundaryCondition`).
- The element vectors of the extension are in lexicographic order, not
  in the native order of the current kernels.
- `GetGradient` returns a matrix free operator, Dirichlet conditions
  included; `ParNonlinearForm::GetLocalGradient` refuses the extension:
  no `HypreParMatrix` through this path.

Consequence: the dedicated operator (entry point B of the plan) is the
main path. The partial assembly of MFEM (entry point A) remains an
option, refused with a clear message if the problem has boundary
integrators or integrators with markers.

## Principles

- **The current path does not change** and remains the default path. The
  existing tests do not change.
- **Work in `mfem-mgis/`, `mfem-mgis-examples/` and `mm-opera-hpc/`,
  branch `gpu-proto` of each** (bases `a011d97d6f`, `8397eec`,
  `8bd56eb`, then `8fac67336e`, `600da70` and `49147f2` since the rebase
  of 13/10/2026), copied into this repository by
  `scripts/sync-repository.sh`. Commit messages in English. A finished
  work package is a set of reviewed commits.
- **New code grouped together**: `include/MFEMMGIS/Batch/` and
  `src/Batch/` (name to be discussed). The existing classes gain
  extension points, no GPU code.
- **Kernels written once**, with `mfem::forall`, generic over the
  hypothesis and the gradient type. They serve on CPU and on GPU, and for
  both entry points. The files of `src/Batch/` are compiled as CUDA
  (CUDA language of CMake) when MFEM is compiled with CUDA, as MFEM does
  for its own sources.
- **All the variants, chosen by option**, compared on the same cases
  (layout, entry point, tangent operator, integration rule, solver).
- **Each option has a test that fails if the option is ignored.**
- **No silent failure**: an unsupported configuration (boundary
  integrator with the partial assembly, `Elasticity` strategy on GPU,
  `generic` behaviour on a GPU device, component by component layout on
  the current path...) is refused with a message.

## Proposed parameters

Names to be discussed, in the style of the existing parameters.

Names chosen by the code (mfem-mgis b5e7d12fc4, 13/10/2026), which
replace those of the table:

- problem: `BehaviourIntegration` (`PointByPoint`, `Batch`),
  `GradientsComputation`, `InnerForcesComputation` and
  `TangentOperatorComputation` (`ElementByElement`, `Kernels`, and
  `MatrixFree` for the tangent), `TangentOperatorAssembly` (`Atomic`,
  `Ordered`), each with its `MFEM_MGIS_*` variable; post-processings
  through the `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION` variable;
- discretization: `IntegrationRuleOrder`;
- `FixedPoint` solver: `Method` (`ElasticOperator`, `ModifiedNewton`,
  `NewtonK0`), `ElasticOperatorUpdate` (`EachResolution`, `Once`),
  `ElasticOperatorComputation`, `ElasticOperatorLinearSolver`,
  `PredictionLinearSolver`, `TangentLinearSolver`, `Acceleration` (with
  the options of each acceleration, for example `Depth` and `Period`) and
  `AdaptiveLinearSolverTolerance`, also an option of `Newton`;
- preconditioners: `OperatorJacobiSmoother`,
  `OperatorChebyshevSmoother`, `AssembledOperator`;
- application: `--device`, `--execution` and the placement options
  ([plan-execution.md](plan-execution.md)).

Not parameters anywhere in the code: `AssemblyStrategy`,
`QuadratureDataLayout`, `TangentOperator`, `Prediction`,
`PreconditionerUpdate`, `PreconditionerOperator`, `LineSearch`,
`NodalBlockJacobi`. The preconditioner kept on `K0` goes through the
`NewtonK0` method.

| Where | Parameter | Values |
|---|---|---|
| problem | `AssemblyStrategy` | `Legacy` (default, current path), `BatchedOperator` (dedicated operator), `PartialAssembly` (partial assembly of MFEM) |
| problem | `QuadratureDataLayout` | `PointByPoint` (default), `ComponentByComponent` (refused with `Legacy`) |
| problem | `TangentOperator` | `Assembled` (option c, default), `MatrixFree` (option a), `ElementMatrices` (element assembly), `ElasticOperator` (option b) |
| problem | `IntegrationRuleOrder` | order of the integration rule; by default the current one (`2 * OrderGrad`) |
| problem | `Prediction` | dictionary taking up the four current axes (`Strategy`, `Operator`, `IntegrationOperator`, `NullTimeIncrement`), plus the new `Extrapolation` strategy |
| nonlinear solver | `Newton` | new options: `AdaptiveLinearSolverTolerance` (Eisenstat-Walker), `LineSearch`, `PreconditionerUpdate` |
| nonlinear solver | `FixedPoint` | `Method` (`ElasticOperator`, `ModifiedNewton`), `ElasticOperatorUpdate`, `LinearSolver` (exact or number of AMG cycles), `Acceleration` (`None`, `Anderson`, `Cast3M`, secants, `IronsTuck`, `Steffensen`, delta squared), `AccelerationDepth`, `AccelerationDamping`, `AccelerationTrigger`, `AccelerationPeriod`, `OperatorPrecision` (see [L5](#fixed-point-and-related-methods)) |
| nonlinear solver | `Newton` | `PreconditionerOperator`: `Tangent` (default) or `Elastic` (frozen AMG of `K0`, `NewtonK0` method) |
| nonlinear solver | `LBFGS` (preconditioned, to be written), `KINSOL` | their own options |
| preconditioner | `OperatorJacobi`, `OperatorChebyshev`, `NodalBlockJacobi`, `LOR` | matrix free |
| preconditioner | `HypreBoomerAMG` | new options: coarsening, interpolation, smoother, sweeps, threshold, levels |
| preconditioner | `ElasticOperator` | wrapper: preconditioner built once on the assembled elastic operator, whatever the operator of the Krylov |
| solver, preconditioner | `Ginkgo*`, `AmgX`, `CuDSS` | according to [gpu-solvers.md](gpu-solvers.md) |
| application | `--device` | configuration of `mfem::Device` (`cpu`, `cuda`, `debug`...), added to `declareDefaultOptions` |

`PreconditionerUpdate`: `EachIteration` (current behaviour),
`EachTimeStep`, `Never`, `EveryNIterations`. With `Never`, the
preconditioner is wrapped in an object whose `SetOperator` rebuilds
nothing.

**The options of the tangent operator, precisely:**

- **option c**: `TangentOperator=Assembled`, Newton, any solver and
  preconditioner on the assembled tangent;
- **option a**: `TangentOperator=MatrixFree`, Newton, Krylov of MFEM,
  `ElasticOperator` preconditioner (AMG on the elastic operator) or
  matrix free preconditioner (`OperatorJacobi`, `OperatorChebyshev`,
  `NodalBlockJacobi`);
- **option b**: `TangentOperator=ElasticOperator` and `FixedPoint`
  solver: the linear operator is the elastic operator, assembled once;
  each iteration applies an approximate solve (one AMG cycle, or a short
  Krylov) to `F(u)`, with Anderson acceleration. The Newton with
  `TangentOperator=ElasticOperator` is the same method without Anderson
  and with a complete Krylov: kept as a variant;
- `TangentOperator=ElementMatrices`: variant of option a, with the
  element matrices kept.

## Overview of the work packages

| Package | Content | Depends on | Size | Completion |
|---|---|---|---|---|
| L0 | build of the three repositories, compatibility with MFEM 4.10, reference, device, small fixes | `mmm-gpu` environment | S | existing tests green with MFEM 4.10 and CUDA, reference times |
| L1 | data at the points in the memory of MFEM, two layouts | L0 | M | green suite, component by component layout on the data |
| L2 | behaviours integrated by batches on CPU | L1 | M | same results, gain measured with `stlpar` and `openmp` |
| L3 | batch kernels, dedicated operator, partial assembly | L2 | L | residual and tangent equal to the current path, green suite |
| L4 | tangent operator (a, b, c, elements), prediction, integration rule | L3 | M | the options on the reference cases |
| L5 | solvers in three waves, fixed point and its variants | L3, L4 | L | wave 1 and fixed point complete |
| L6 | GPU | L3 (minimal slice), then L4, L5 | L | same results as on CPU, no hidden copy with one process |
| L7 | all the hypotheses, kinematics and models | L3 to L6 | L | complete suite of mfem-mgis, completed with the missing tests, on the batch path |
| L8 | outputs and post-processings | L1, L6 | S | identical outputs, cost measured |
| L9 | examples and measurements | L6 | M | measurement campaign published |
| L10 | documentation, upstream proposals | all | M | documents and branches ready |

Sizes: S (a few days), M (one to two weeks), L (several weeks). Critical
path: L0, L1, L2, L3, GPU slice of L6, L4, L5. The minimal GPU slice (3D
small strains, one behaviour, option c then a) comes early to remove the
risks of the GPU before generalizing.

## L0. Working base

Objective: build mfem-mgis in the environment of the prototype, with the
time reference of the current path.

Work:

- Branch `gpu-proto` in `mfem-mgis/`. Script `build-mfem-mgis.sh` at the
  root (like `build-mgis.sh`): MFEM 4.10 with CUDA from `mmm-gpu`, MGIS
  from `install/mgis`, TFEL from `install/tfel-mmm-gpu`, build,
  `make check`.
- CMake of mfem-mgis: CUDA language enabled when MFEM is compiled with
  CUDA, for the files that will contain kernels.
- **Compatibility with MFEM 4.10**: mfem-mgis has only been used with 4.8
  and 4.9; fix what 4.10 breaks (removed or deprecated APIs), in separate
  commits, that can be proposed upstream.
- mfem-mgis-examples and mm-opera-hpc built in the same way
  (`build-mfem-mgis.sh` builds the three), their tests green.
- Check that the current suite passes with MFEM compiled for CUDA but
  with the `cpu` device. The MFEM of the environment is compiled with
  Umpire: its default memory types are then those of Umpire.
- Integers of hypre: the spack recipe of mfem-mgis requests hypre
  `+int64` by default, which hypre refuses with CUDA ("use +mixedint");
  the environment uses hypre with 32-bit integers. To be documented.
- `--device` option: `mfem::Device` created by `mfem_mgis::initialize`
  (or by the application) from a string (`Config.cxx:174-222`). Test the
  order with `mfem::Hypre::Init`. Any hypre preconditioner must be
  created after the device (the default values of BoomerAMG depend on
  the GPU, `linalg/hypre.cpp:5029`).
- Measurements per phase: extend the existing profiling sections
  (`CatchTimeSection`, `ProfilingSection`) to all the phases
  (integration, residual, tangent, preconditioner, solve, end of step),
  with a synchronization of the GPU before each reading of a time.
- Fix the two defects of `HypreGMRES` and `HypreFGMRES` (small,
  independent, can be proposed upstream right away).

Tests: existing suite; a test of the `--device` parameter (the device is
the one requested).

Completion: green suite, reference times noted for the cases of the
[measurement campaign](#measurement-campaign).

## L1. Data at the points in the memory of MFEM

Objective: the arrays of the materials live in a memory chosen by
mfem-mgis (`mfem::Memory`, host or GPU), in one of the two layouts. MGIS
uses them without copy.

Work:

- `Material`: build the `MaterialDataManager` with a
  `MaterialDataManagerInitializer` (`Material.cxx:48-61`). `Material`
  derives directly from `MaterialDataManager`, built before the members:
  the memory must be held by a base declared before
  `MaterialDataManager` (base-from-member idiom), or created by a factory
  before the material. The arrays (gradients, thermodynamic forces,
  internal variables, energies of the two states, blocks of the tangent)
  are `mfem::Vector`; their pointers are given to MGIS with the location
  (`HOST` or `DEVICE`, from `mfem::Device`) and the layout
  (`QuadratureDataLayout`). Variant to compare: arrays allocated by MGIS
  with a memory resource (`MGIS/MemoryResources`), seen by MFEM.
- **Validity of the copies**: an `mfem::Memory` has a host copy and a
  device copy, with validity flags that MFEM keeps up to date. MGIS
  writes through raw pointers: the pointers given to MGIS must be
  obtained by `Write` or `ReadWrite` (on the device) just before each
  call, and any host access must go through `HostRead` or
  `HostReadWrite`. This policy is written and tested in this work
  package.
- `PartialQuadratureSpace`: add, next to the hash table,
  `mfem::Array<int>` arrays of the elements of the material and of the
  offsets, readable on the device. The table stays for the current code.
  A rework of this class is in progress upstream (see
  [Risks](#risks)).
- **Component by component layout**: the current path assumes the point
  by point layout everywhere (contiguous accesses, `o * stride`). This
  layout is therefore refused with `AssemblyStrategy=Legacy`, and
  reserved for the batch path. Done in L1: the layout is chosen when a
  material is built; the integrators of the current path always build
  their materials point by point, and the stride-based accessors and
  evaluators refuse the other layout with a message. The parameter of
  the problem (`QuadratureDataLayout`) and its refusal with
  `AssemblyStrategy=Legacy` come with the batch path (L2, L3). In this
  work package, only the data of the materials and their description
  change; the accesses of the current path remain point by point.
- **Host access through the public API** (accessors of `Material`,
  post-processings, projections, `StateManager`, gradient evaluators,
  initialization of the internal variables): each accessor synchronizes
  to the host (`HostRead`); each host write is declared
  (`HostReadWrite`) so that the device copy is invalidated. With the
  component by component layout, these views (based on a stride between
  points) are refused, or served by a transposed copy.
- **Material properties and external variables in the host memory of
  the user** (API of MGIS): copied to the memory of the behaviour at each
  `setup` when the behaviour computes on GPU, or refused with a message
  that says how to proceed.
- End of step: `update` and `revert` go through the overloads of MGIS
  that take `FieldOperations`. MGIS provides `HostFieldOperations` and
  CUDA, HIP, SYCL, Kokkos and OpenMP operations in headers, to be
  compiled by mfem-mgis with the corresponding toolchain. Variant to
  compare: `MFEMFieldOperations`, written with `mfem::forall` and the
  copies of `mfem::Memory` (compiled as CUDA when MFEM is).

Tests: complete suite unchanged (point by point layout, data in
`mfem::Memory`); refusal of the component by component layout by the
accessors and the evaluators based on a stride between points (the
refusal with `Legacy` comes with the parameter of the problem, in L2 or
L3); host accessors after a computation on the device (up to date
values); unit tests of the end of step operations in both layouts; test
of the validity policy under the `debug` device of MFEM.

Completion: green suite, data in the memory of MFEM.

## L2. Behaviours integrated by batches on CPU

Objective: replace the point by point integration by one call per
material, without touching the assembly. First measurable gain: the
`stlpar` and `openmp` backends of TFEL on CPU.

Work:

- Split `integrate` (`...Implementation.cxx:485-514`):
  1. element by element pass, which only computes the gradients (current
     code of `implementIntegrate` without the call to the behaviour),
     including the macroscopic gradients and the update of the external
     variables from the unknowns (trait
     `updateExternalStateVariablesFromUnknownsValues`), then rotation of
     the gradients;
  2. for each material, one call to `mgis::behaviour::integrate` on a
     view of the material (`make_view`), with the chosen execution
     policy;
  3. rotation of the blocks of the tangent only, on views. **The
     thermodynamic forces are not rotated**: they stay in the material
     frame, the current assembly rotates them in a copy.
- Several integrators on the same attribute: one call per integrator
  (each one is a material).
- The `generic` libraries remain accepted: with a workspace, MGIS
  integrates them point by point (reference). The `generic-parallel`
  libraries are chosen by the path of the library.
- Material properties and external variables: the values evaluated in
  `setup` are given to the view (uniform or arrays).
- Failures: MGIS returns the first failing point and its status; the
  current behaviour (division of the Newton correction, substeps) must
  be kept. Change of behaviour to be documented: by batches, all the
  points are integrated (today the integration stops at the first
  failing element), so the state changes beyond the failing point; the
  iteration is rejected anyway. Status 0 ("unreliable results") and the
  reduction of the status between MPI processes are kept as they are.
- **Material properties at the beginning of the step**: the `generic`
  behaviours only read those at the end of the step, but the batch
  processing of MGIS checks both states (`mgis/src/Integrate.cxx`,
  `checkMaterialStateView`). A problem that only sets those at the end
  of the step would be refused: mfem-mgis also sets those at the
  beginning of the step (same values at `setup`), or MGIS accepts their
  absence when the behaviour does not read them (to be decided, see
  [Questions](#questions-to-settle)).

Tests: complete suite with a `generic` library (identical to the current
path), then with `stlpar` and `openmp` (identical up to rounding);
orthotropic behaviour (rotations); several integrators on one attribute;
test of a failing point.

Completion: same results, gain measured on the beam and the polycrystal.

Implemented (09/10/2026), differences from the plan:

- provisional parameter `BehaviourIntegration` (`PointByPoint`,
  `Batch`), default value given by the environment variable
  `MFEM_MGIS_BEHAVIOUR_INTEGRATION`; it will join `AssemblyStrategy`;
- material properties at the beginning of the step (question 6):
  mfem-mgis gives the behaviour, for the call, views on those at the end
  of the step, like the point by point path; MGIS is not modified;
- thermodynamic forces: still rotated point by point by the residual
  when the behaviour provides the rotation of one point; for an
  orthotropic `generic-parallel` behaviour (batch rotations only),
  rotated for all the points after the integration, in a separate array
  (the forces of the material stay in its frame). A cache invalidated by
  versions of the arrays was tried then removed after the review: the
  writes of MGIS through its views escape it;
- failures: MGIS stops at the first failing point for a `generic`
  behaviour and integrates all the points of a `generic-parallel`
  behaviour;
- the gain depends on the case: negligible on the beam (cheap
  behaviour), significant on the polycrystal (see the measurements).

## L3. Batch kernels and dedicated operator

Objective: gradients, residual, action of the tangent and element
matrices computed by batch kernels, without an element by element pass
on the host.

Kernels (`src/Batch/`, `mfem::forall`, 3D small strains first):

- **geometry**: geometric factors (`Mesh::GetGeometricFactors`,
  Jacobians and determinants) and basis functions (`DofToQuad`) for the
  integration rule of each material. By default, the same rule as the
  current path (`2 * OrderGrad`). The kernels of MFEM require one rule
  per geometry type: a mixed mesh is handled geometry by geometry. Done:
  for a mesh without nodes, the Jacobians are computed by a kernel from
  the vertices (the geometric factors of MFEM required a finite element
  space of order one and a projection of the coordinates on the host,
  1.8 s for 1.6 million tetrahedra);
- **gradients**: element values (`ElementRestriction`), derivatives at
  the points, conversion to the gradients of the behaviour (strain in
  Voigt notation with √2, like the current `B` matrices), macroscopic
  gradients added, writing into the end of step gradients of the
  material, in its layout, for the elements of the material only;
- **residual**: `Bᵀ σ w det J` per point (with the rotation of the forces
  in a copy for the orthotropic behaviours), sum per element, then
  transpose of the restriction;
- **action of the tangent**: `Σ Bᵀ K B v`;
- **diagonal** and **3×3 nodal blocks** of the tangent (smoothers);
- **element matrices** `Ke = Σ Bᵀ K B w det J`.

The `B` matrices are generated by GiNaC today. Two variants to compare:
kernels written by hand, generic over the hypothesis, or generated by an
extension of `code-generation/behaviour-integrator.cxx`.

Entry points:

- **B, dedicated operator** (`AssemblyStrategy=BatchedOperator`, first).
  `Mult` and `GetGradient` of `NonLinearEvolutionProblemImplementation<p>`
  overridden. The operator redoes what `(Par)NonlinearForm` does today:
  - prolongation `P` and `Pᵀ` (on GPU,
    `DeviceConformingProlongationOperator`);
  - imposed degrees of freedom of the form (`SetEssentialTrueDofs`,
    periodic included): zero residual, tangent eliminated or
    `ConstrainedOperator`, diagonal set to 1 for the smoothers;
  - **all the other lists of integrators of the form** (domain with
    markers, boundary, interior faces, boundary faces: pressure, heat
    source, Robin of ex8, integrators of the user), except the integrator
    of mfem-mgis: called on the current path for the residual and the
    assembled tangent; for option a, matrix free action of these terms,
    or separate assembly; refusal with a message otherwise, never a
    silent omission;
  - the current domain integrator is no longer added to the form
    (otherwise the domain would be counted twice);
  - if `Mult` integrates the behaviour, the "new estimate" action of the
    solver of mfem-mgis no longer integrates it (otherwise two
    integrations per iteration).

  `Mult` integrates the behaviour then assembles the residual (option),
  which makes external nonlinear solvers possible; a failure of the
  behaviour is reported without stopping the program. `mgis_integrator`
  (`final` type) becomes an interface common to both paths, or the batch
  path has its own object. Explicit refusals on the batch path: element
  attribute without an integrator (otherwise zero contribution and
  singular matrix), `UseMultiMaterialNonLinearIntegrator=false` mode,
  NURBS meshes and collections other than H1, integrators other than the
  18 standard integrators of the generator as long as they have no
  kernel.
- **A, partial assembly of MFEM** (`AssemblyStrategy=PartialAssembly`):
  a `NonlinearFormIntegrator` that implements `AssemblePA`, `AddMultPA`,
  `AssembleGradPA`, `AddMultGradPA`, `AssembleGradDiagonalPA`, with the
  element vectors in lexicographic order. Explicit refusal if the problem
  has boundary integrators or integrators with markers.

Tests:

- per element: element residual and matrix of the kernels against
  `AssembleElementVector` and `AssembleElementGrad`, random
  displacements, hexahedral and tetrahedral meshes, orders 1 and 2,
  several materials, byNODES and byVDIM numbering, orthotropic
  behaviour;
- global: residual, action of the tangent, diagonal against the current
  path, sequential and on 2 processes;
- problem: `UniaxialTensileTest`, `TwoLayerCube`, `PredictionTest`,
  `ImposedPressureTest`, `ImposedDirichletBoundaryConditionAtClosestNode`,
  a periodic problem, with both entry points (the second one must refuse
  the pressure): same iterations, same results; a test of each explicit
  refusal.

Completion: suite of the 3D small strains green on both entry points and
both layouts.

## L4. Tangent operator and prediction

Objective: the options of the tangent operator (see
[Proposed parameters](#proposed-parameters)), the choice of the
integration rule, and a prediction that no longer uses element by
element bilinear forms.

Work:

- `TangentOperator=Assembled` (c): element matrices of the kernels, then
  `HypreParMatrix`. On CPU first. Sparsity pattern, hypre matrix and
  communication scheme kept from one iteration to the next (only the
  values change). Done (11/10): positions of the entries of the element
  matrices in the rows of the sparse matrix computed once, without GPU
  (`computeStiffnessMatricesPositions`, assembly ×1.1 to ×2.9 on one
  core, no gain on the GPU), and products `K B` and `Bᵀ K B` written
  without forming `B` (constexpr table of the nonzero entries: 254 then
  64 registers in 3D finite strain, assembly ×1.1 to ×3.1 on one core and
  ×1.4 to ×1.9 on the GPU, identical results). Inverses of the Jacobians
  computed once by the geometry kernel: assembly on the GPU ×1.7 to
  ×5.9. Global matrix with pattern and exchanges kept
  (`batch::GlobalStiffnessMatrix`, conforming spaces, hypre on the host):
  product `Pᵀ A P` of 15 to 37 ms replaced by 2 ms on 8 processes.
  Remaining: the prediction (same class) and the case of hypre on the
  GPU.
- `TangentOperator=MatrixFree` (a): action of the kernels, diagonal and
  nodal blocks for the smoothers, AMG built once on the assembled
  elastic operator.
- `TangentOperator=ElementMatrices`: element matrices kept, element by
  element action.
- `TangentOperator=ElasticOperator` (b): with `FixedPoint` (L5).
- `IntegrationRuleOrder`: rule chosen by the user, for example 8 points
  on the hexahedron of order 1 instead of 27; check of the absence of
  spurious modes on the reference cases. Done (11/10): parameter of the
  discretization, same order for all the geometries, rules with
  non-positive weights refused, quadrature spaces and integrations by
  the same rule (`applyIntegrationRuleOrder`).
- Prediction (`...Implementation.cxx:98-233`) rewritten with the
  kernels: elastic, secant or tangent operator by the prediction
  integration of MGIS, right-hand side by the residual kernel plus the
  contribution of the boundary conditions (`addLinearFormIntegrators`),
  respect of `DiscardLinearSolverFailure`, and the four current axes
  kept. New strategy: extrapolation of the increment of the previous
  step. Parameter `Prediction`. `getLinearizedOperators` (public API):
  returned by the kernels, or refused on the batch path with a message.
- Sequential: option c gives a `SparseMatrix` (UMFPACK, Krylov of MFEM),
  the parallel case a `HypreParMatrix`.
- Solver or preconditioner provided by the user that requires a matrix
  with `MatrixFree` or `ElasticOperator`: explicit refusal.
- `IntegrationRuleOrder` changes the size of the quadrature spaces, the
  functions at the points of the user, the links between problems
  (`resolveBehaviourIntegratorsDependencies`) and the references of the
  tests: the integrals (`computeMeasure`, `computeIntegral`, energies,
  resultant) follow the chosen rule; two linked problems keep the same
  rule.

Tests: each option on the reference cases (same solutions up to the
Newton tolerance); a test per option that fails if it is ignored (number
of iterations, type of the returned operator, number of integration
points); predictions against `PredictionTest`, `PredictionTest2` and
`ImposedPressureTest2` (the only test that combines prediction and
pressure); option c in sequential with UMFPACK.

Completion: the options and the predictions on the reference cases.

## L5. Solvers

Objective: connect the solvers of [gpu-solvers.md](gpu-solvers.md), in
waves, in the factories of mfem-mgis.

**Wave 1** (required for options a and b and for the GPU):

- `PreconditionerUpdate` in the Newton of mfem-mgis
  (`NewtonSolver.cxx:245`): the linear solver still receives the new
  tangent, but its preconditioner is wrapped in an object whose
  `SetOperator` rebuilds nothing when it should not (with the Krylov
  solvers of hypre, the setup phase of the preconditioner must also be
  neutralized); `ElasticOperator` wrapper;
- `FixedPoint` solver with Anderson acceleration;
- inexact Newton (Eisenstat-Walker, done: option
  `AdaptiveLinearSolverTolerance` of the Newton and of `NewtonK0`) and
  line search, written in the Newton of mfem-mgis (the hooks of
  `mfem::NewtonSolver` do not apply there);
- matrix free smoothers: `OperatorJacobi`, `OperatorChebyshev`,
  `NodalBlockJacobi`;
- options of BoomerAMG exposed, without being overwritten by
  `SystemHypreBoomerAMG::SetOperator`; `Elasticity` strategy refused on
  GPU (MFEM stops the program in byVDIM), `System` proposed instead;
- `hasConverged` and `getNumberOfIterationsAtConvergence` extended to
  all the added solvers (today any unknown solver is declared converged:
  a failure would go unnoticed), with a test that makes each solver fail;
  `DiscardLinearSolverFailure` and the iteration information kept in the
  new solvers;
- fixed point and related methods: see the next section.

**Wave 2** (comparisons on GPU):

- preconditioned L-BFGS, to be written (`LBFGSSolver` of MFEM ignores
  the preconditioner; `SNESQN` of PETSc as the reference);
- Ginkgo in the factory: distributed solvers and preconditioners of
  MFEM 4.10 (Schwarz), local preconditioners (block Jacobi with adaptive
  precision, ILU, ISAI); cuDSS through `CuDSSSolver` (reuse of the
  symbolic analysis; done on 12/10, `CuDSSSolver` solver of the factory,
  ex8 on the GPU 9.8 -> 4.4 s);
- AmgX through `AmgXSolver`, with `UpdateOperator` for the complete
  reuse;
- KINSOL (`KINSolver`), with the operator that integrates the behaviour
  (L3);
- repair of the PETSc path (same operator, no more stop on a failure of
  the behaviour), then GAMG with near null space;
- mixed precision (Ginkgo, then AmgX `dDFI` through its own wrapper).

**Wave 3** (depending on the measurements): p-multigrid and LOR (order 2
and higher), SuperLU_DIST and STRUMPACK on GPU, MGR of hypre for the
coupled problems, FSAI and ILU of hypre, pipelined Krylov, multigrid of
Ginkgo (through its API), partial reuse of the AMG hierarchy (GAMG, AmgX
through its own wrapper).

Tests: each solver on a case where it converges, against a reference
solution; each option with a test that fails if it is ignored (number of
builds of the preconditioner for `PreconditionerUpdate`, number of
iterations for Anderson, settings of BoomerAMG read back after several
iterations...).

With an iterative coupling scheme (several solves in the same step), the
reuse policy of the preconditioner also applies between coupling
iterations.

Completion: wave 1 complete and measured; next waves along with L9.

### Fixed point and related methods

Preliminary study (04/10/2026, outside this repository, prototype
connected to the solver factory of mfem-mgis): on a beam in J2
plasticity (112,000 unknowns), the fixed point on the elastic operator
`K0`, accelerated by Anderson, with a single frozen AMG cycle per
iteration, is 2.15 times faster than the best Newton on CPU, and about 4
times on GPU (2 times faster than a Newton with a frozen AMG of `K0`).
On a polycrystal (expensive behaviour, weak hardening), it is 3.3 times
slower on CPU: 8 times more iterations, hence 8 times more integrations
of the behaviour. The GPU shifts exactly this ratio. Hence all the variants
below, compared on the same cases.

**Methods** (`FixedPoint` nonlinear solver, `Method` parameter):

- `ElasticOperator`: fixed point on `K0` (initial stiffness, as in
  Cast3M). `K0` is assembled once, for the whole computation or once per
  step (`ElasticOperatorUpdate`). Each iteration: integration of the
  behaviour, residual, one approximate solve with `K0`, acceleration;
- `ModifiedNewton`: tangent assembled and preconditioned (or factorized)
  at the first iteration of the step only;
- `NewtonK0` (done as a method of the `FixedPoint` solver): consistent
  tangent operator, Krylov preconditioned by the AMG of `K0` built once.
  Same iterations as Newton, without building an AMG at each iteration;
  it is the real competitor of the fixed point on GPU. The Krylov is an
  iterative solver of MFEM from the factory (`CGSolver` by default,
  `GMRESSolver`, `FGMRESSolver`, `BiCGSTABSolver`, `MINRESSolver`,
  `SLISolver`): the conjugate gradient assumes a symmetric positive
  definite tangent.

**Solve with `K0`** (`LinearSolver` parameter):

- exact: a single factorization (MUMPS or UMFPACK on CPU, cuDSS on GPU);
- inexact: `n` cycles of the frozen AMG of `K0` per iteration (a single
  one is enough with Anderson), or a short Krylov preconditioned by this
  AMG.

**Accelerations** (`Acceleration` parameter), all ported from the
accelerators of MTest (TFEL, `mtest/src/*AccelerationAlgorithm.cxx`,
unified framework published in Computers & Mathematics with Applications
70(9), 2015) to the vectors of MFEM: scalar products and combinations of
vectors, hence on GPU without dedicated code:

- `None` (reference: in J2 with linear hardening `H`, the rate is
  `3μ/(3μ+H)` in the plastic zones, predictable);
- `Anderson` (on the iterates `UAnderson` and on the residuals
  `FAnderson`; depth, damping, trigger, period);
- `Cast3M`;
- secants: `Secant`, `CrossedSecant`, `AlternateSecant`;
- `IronsTuck`, `Steffensen`;
- delta squared: `CrossedDelta2`, `Crossed2Delta`, `Crossed2Deltabis`,
  `AlternateDelta2`, `Alternate2Delta`.

The prototype had only ported Cast3M and Anderson. The literature review
of the study shows that Anderson is not always the best: hence the
porting of the whole family.

**Precision** (`OperatorPrecision` parameter): `K0` or its AMG in single
precision, residual and acceleration in double. On the RTX 5080, the
study measures a gain of about 30 % for the same number of iterations.

**Prediction**: the elastic prediction (essential in plasticity) reuses
`K0` already factorized or its AMG: no more factorization per step (it
cost 35 s out of 44 with MUMPS in the study).

**Measurements specific to these methods**: time per phase
(integration, residual, assembly, build of the preconditioner, solves,
acceleration), as in the prototype; iterations; and the "best against
best" comparison (never against a Newton with a direct solver only,
which distorts the gains).

**Tests**: convergence rate of the non accelerated fixed point equal to
`3μ/(3μ+H)` on a J2 case with linear hardening (exact test); each
accelerator on a linear problem (Anderson is equivalent there to GMRES
preconditioned by `K0`); same solutions as Newton up to the tolerance on
the reference cases; a test per parameter that fails if it is ignored
(number of builds of the AMG, number of factorizations).

## L6. GPU

Objective: the whole computation of an iteration on the GPU, data
included.

Minimal slice (as soon as L3 is finished): 3D small strains, one J2
behaviour compiled for `cuda`, option c then a.

Work:

- `cuda` device of MFEM, arrays of the materials on the GPU (L1),
  kernels of L3 on the GPU (files compiled by nvcc), integration by the
  `cuda` (or `kokkos`) library of the behaviour;
- refusal of a `generic` library on a GPU device, except for an explicit
  mixed configuration (material in shared memory);
- Dirichlet conditions applied on the GPU (`Vector::SetSubVector` with a
  value runs on the device since MFEM 4.9); no more host write into `u`
  for the conditions of mfem-mgis. Conditions written by the user
  (public interface that writes `u` on the host): `u` is synchronized to
  the host before the call and to the device after (correct, but
  costly), and an interface for the device is proposed;
- solvers that only work on the host (UMFPACK, MUMPS) refused with a GPU
  device, except for an explicit mixed configuration (copy of the
  matrix);
- end of step on the GPU by the CUDA operations of MGIS (or the
  `MFEMFieldOperations` variant);
- option c on GPU: global matrix assembled on the GPU from the element
  matrices. The sparse matrix filling of MFEM on the device
  (`ElementRestriction::FillJAndData`) works on scalar degrees of
  freedom: for 3 components, to be checked, or to be done by blocks of
  components or through the IJ interface of hypre on GPU. Elimination of
  the imposed degrees of freedom on the device;
- several materials: asynchronous calls to MGIS then a single
  synchronization;
- check that no hidden copy takes place: GPU without HMM and Nsight
  Systems profile on the real GPU; validity policy of the copies checked
  by the `debug` device of MFEM, which runs the kernels on the host (with
  an `openmp` or `stlpar` behaviour). With several processes, without
  GPU-aware MPI, the exchanges of the prolongation and of hypre go
  through the host: the criterion only holds with one process.

Tests: the tests of L3 and L4 on the `cuda` device, compared with the
CPU up to rounding; the same under the `debug` device with an `openmp`
behaviour.

Completion: same results as on CPU on the reference cases, no hidden
copy with one process, first measurements.

## L7. All of mfem-mgis

Objective: the batch path covers what the current path covers.

Proposed order:

1. plane strain and plane stress;
2. finite strains (deformation gradient, Piola-Kirchhoff stress,
   conversions of MGIS on views);
3. stationary and transient nonlinear thermal (the transient one
   requires a kernel on the values of the shape functions, not only on
   their derivatives; the temperature read from the unknowns as an
   external variable must live on the GPU). Stationary done (temperature
   gradient, flux, block of the derivative of the flux with respect to
   the temperature, values of the shape functions, temperature at the
   points copied into the external variable on the host). Transient done
   (gradient = temperature at the point, operator B = value of the shape
   function, `TransientHeatTransfer` integrator also declared in 2D);
   remaining: the external variable kept on the GPU;
4. orthotropic behaviours: rotation matrices per point in an array,
   rotations of MGIS on views (no more function pointers);
5. material properties and external variables in all their forms,
   including the `generic-parallel` properties evaluated by MGIS on the
   GPU; non uniform evaluators copied to the GPU at each `setup` after
   synchronization of their source;
6. boundary conditions on boundaries (pressure, heat source) by
   boundary kernels, for the GPU. Done differently (11/10): these
   conditions do not depend on the unknowns, their contribution is
   assembled once for a unit value, then multiplied by the value at the
   end of the step and added on the GPU to the internal forces computed
   by kernels (`AbstractBoundaryCondition::isIndependentOfTheUnknowns`).
   The nonlinear form only contains the integrator of mfem-mgis: no more
   residual assembled on the host, tangent assembled by kernels without
   the assembly of MFEM, tangent without assembly possible with a
   pressure. Remaining: the conditions that depend on the unknowns (Robin
   of ex8), which would require real boundary kernels;
7. models and coupling schemes (`Simulation`, `IterativeCouplingScheme`,
   `LoopCouplingScheme`, `NonLinearModel`, `PointWiseModel`,
   `QPEvaluators`), dependencies between integrators of two problems
   (`resolveBehaviourIntegratorsDependencies`): the thermomechanical
   example ex8 depends on them, and the exchanges of fields at the points
   between models are done on the host. `PointWiseModel` (behaviour at
   the points without finite elements) is a first simple case of batch
   behaviour on GPU;
8. micromorphic;
9. regularizations (FBar, Faltus): they require the current solution at
   assembly, to be handled last. FBar done (11/10): gradient at the
   center of the elements by a rule of order 0, `F̄ = α F`, standard
   tangent modified plus a rank one term per element (see
   `batch::FBarData`), in 3D and in plane strain, assembled matrix,
   diagonal and matrix free action. Faltus done (11/10): gradient at the
   center with the same rule, residual `Σ w Dᵀ α (F - F0)`, constant
   tangent `α Σ w Dᵀ D` computed on the fly by pairs of nodes (block
   `d_a · d_b` times the identity, `d` difference of the derivatives of
   the shape functions at the point and at the center), in 3D and in
   plane strain (plane stress refused);
10. periodic problems.

Tests: the suite of mfem-mgis run on the batch path, case by case. The
current suite tests neither FBar, nor Faltus, nor plane stress, nor
transient thermal, nor the heat source, nor the iterative coupling, nor
PETSc: these tests are to be written in this work package (on the
current path first, as the reference), or the corresponding examples of
mfem-mgis-examples enter the validation.

Completion: complete suite, completed with these tests, on the batch
path, on CPU and on GPU.

## L8. Outputs and post-processings

- First: copy to the host before each post-processing (`HostRead`),
  current behaviour unchanged. With the component by component layout,
  the accesses of the post-processings are to be adapted or the data to
  be transposed into a point by point copy.
- Then, depending on the measurements: average thermodynamic forces,
  resultant on a boundary, energies and integrals (`computeMeasure`,
  `computeIntegral`) by reductions on the GPU, with the chosen
  integration rule; L2 projections and implicit gradient regularization
  (`L2Projection.cxx:624-702`) by the partial assembly of MFEM.

Tests: identical outputs on CPU and on GPU, in both layouts.

Status on 13/10/2026: done. Post-processings by kernels (variable
`MFEM_MGIS_POST_PROCESSINGS_COMPUTATION`), measured in
`results/post-processings-kernels-2026-10-11` and
`results/L8-polycrystal-2026-10-13-1027`.

## L9. Examples and measurements

- mfem-mgis-examples and mm-opera-hpc: options `--device`,
  `AssemblyStrategy`, `TangentOperator`, `QuadratureDataLayout`,
  `IntegrationRuleOrder`, solver, passed on the command line.
- Special cases of these repositories to handle: the Robin condition of
  ex8 (boundary face integrator whose tangent depends on `u`) and its
  thermomechanical coupling; the material properties set in the host
  memory of the user (`cermet`); the periodic polycrystals (degrees of
  freedom imposed by the periodic problem); the bubble (`bubble`).
- Measurement campaign (next section), results in `results/`, scripts
  in `scripts/`.

Status on 13/10/2026: done, with other names. The examples and the cases
take the options of mfem-mgis (`--device`, `--execution` and the
placement), the library of the behaviour, the solvers and the profiling
(`patches/mfem-mgis-examples/`, `patches/mm-opera-hpc/`);
`AssemblyStrategy`, `TangentOperator` and `QuadratureDataLayout` do not
exist (see [Proposed parameters](#proposed-parameters)). The measurement
scripts are in `cases/`.

## L10. Documentation and upstream proposals

- Documentation of mfem-mgis (in English, in its tree): batch path,
  parameters, tested configurations.
- Split into independent upstream proposals: fixes (L0), data at the
  points (L1), batch behaviours (L2), kernels and operator (L3), solvers
  (L5), GPU (L6). See [upstream-integration.md](upstream-integration.md).

Status on 13/10/2026: documentation of mfem-mgis in its tree (for
example `docs/source/user_guide/batch_computations.rst`) and split in
[upstream-integration.md](upstream-integration.md) done; nothing is
proposed yet to the official repositories.

## Test strategy

- A single reference: the current path, with a `generic` library.
- Three levels of comparison: element (element residual and matrix),
  global (residual, action of the tangent, diagonal), problem (Newton
  iterations, displacements, stresses, internal variables). Relative
  tolerance 1e-12 at the first two levels, up to rounding at the third.
- Each new option has a test that fails if the option is ignored.
- GPU: same tests on the `cuda` device, compared with the CPU; same tests
  under the `debug` device of MFEM with an `openmp` behaviour, for the
  validity of the copies.
- Targeted pass (about twenty minutes) at each commit, complete suite at
  the end of a work package, independent review before declaring a work
  package finished.

## Measurement campaign

Cases:

1. beam in J2 plasticity, 3D, small strains (case of the preliminary
   trials);
2. ssna303 (mfem-mgis-examples ex2, finite strains);
3. polycrystal (mm-opera-hpc, several materials, single crystal
   behaviour);
4. nonlinear thermal;
5. cube refined until the memory of the GPU is full (scaling).

Measurements: time per step and per phase (gradients, behaviour,
residual, tangent, build of the preconditioner, solve, end of step,
outputs), Newton and Krylov iterations, memory of the GPU, transfers
between host and GPU (including those of the MPI exchanges without
GPU-aware MPI).

Comparisons:

- CPU (1 to 16 cores, MPI) against GPU;
- current path against batch path on CPU;
- point by point layout against component by component;
- options a, b, c and elements of the tangent operator;
- integration rule (27 points against 8 on the hexahedron of order 1);
- solvers, wave by wave;
- `cuda` library against `kokkos` of the same behaviour, and `generic`
  on the host.

## Questions to settle

1. **Names of the parameters** proposed above. Settled by the code (see
   [Proposed parameters](#proposed-parameters)).
2. **Order of the entry points**: dedicated operator first (proposed),
   because of the boundary integrators ignored or refused by the partial
   assembly of MFEM.
3. **MFEM 4.10**: settled on 09/10/2026, the whole prototype moves to
   4.10 (see [gpu-solvers.md](gpu-solvers.md#mfem-410)).
4. **Kernels written by hand or generated** (extension of the GiNaC
   generator of mfem-mgis).
5. **Owner of the memory**: mfem-mgis (proposed) or MGIS with a memory
   resource; both are planned in L1, to be compared.
6. **Material properties at the beginning of the step** (see L2):
   settled in L2, mfem-mgis gives views on those at the end of the step,
   like the point by point path.
7. **Working base**: `a011d97d6f`, or the rework of
   `PartialQuadratureSpace` in progress upstream once merged. Done on
   13/10/2026: branch rebased onto master `8fac67336e`.
8. **Upstream proposals**: when and in which order.

## Risks

- **Size**: many integrators and special cases. The batch path starts
  with the most used cases and remains optional.
- **Upstream branches in progress**: the rework of
  `PartialQuadratureSpace` (branch
  `466-refactor-partial-quadrature-space-to-better-support-spaces-on-boundaries`,
  89 files) touches the classes that L1 modifies; the branch
  `452-post-processings-only-consider-the-first-behaviour-integrator-of-a-material`
  touches the post-processings of L1 and L8. Conflicts to be expected.
- **Partial assembly of MFEM**: limits described above; entry point B
  does not depend on the extension.
- **CUDA compilation in mfem-mgis**: the kernels require nvcc;
  compilation time and memory to be monitored, as for MFEM.
- **Cost of the AMG on GPU**: setup more expensive than the solve; hence
  the reuse of the preconditioner in wave 1.
- **Double precision on a consumer GPU**: the RTX 5080 computes in double
  precision 64 times slower than in single; the measurements must be
  confirmed on a compute GPU (A100 or H100).
- **MFEM 4.10**: mfem-mgis has only been used with 4.8 and 4.9;
  compatibility checked in L0. Compilation with CUDA 13: `ptxas` takes
  22 GB for one file of MFEM (builds kept in check by
  `scripts/memory-governor.sh`).
- **Local failures of the integration**: on GPU, a failing point must be
  reported without copying the whole state; MGIS returns the first
  failing point and a global status.
- **Several processes without GPU-aware MPI**: copies to the host at each
  exchange, unavoidable with the local environment.

- **Several GPUs per node and two CUDA runtimes**: the kernels of
  mfem-mgis go through the shared CUDA runtime, those of MFEM through the
  static copy of the runtime contained in `libmfem.so`. Both follow the
  current context of the driver, hence a priori the same GPU, but this
  cannot be checked on a machine with a single GPU. To be tested on a
  node with several GPUs (one GPU per process, operations of MFEM and
  mfem-mgis on the same arrays).
- **External storage between materials**: a material can read its
  external variables in the arrays of another one (`EXTERNAL_STORAGE`
  mode, example ex8: swelling read in the internal variables of the
  swelling model). The synchronization of one material does not cover
  the arrays of the other: to be handled when these arrays are written
  on the GPU.
