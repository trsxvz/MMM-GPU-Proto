# What was added

This page lists what the prototype adds to each library. The sources
are `docs/tfel.md`, `docs/mgis.md`,
`docs/plan-mfem-mgis-detailed.md`, `docs/plan-execution.md`,
`docs/plan-precision.md` and `docs/gpu-solvers.md`. The documentation
of the interfaces, in English, is in the tree of each library:
`tfel/docs/web/generic-parallel-behaviour-interface.md`,
`tfel/docs/web/generic-parallel-material-property-interface.md`,
`mgis/docs/web/generic-parallel.md` and `mgis/docs/web/views.md`.

## TFEL and MFront

### The `generic-parallel` interface for behaviours and models

MFront had a `generic` interface that integrates a behaviour at one
point. The new `generic-parallel` interface generates, for the same
behaviour, functions that process a whole batch of points.

- **Generated functions.** For each behaviour and each modelling
  hypothesis: the integration of a batch, one function per initialize
  function, one function per post-processing, and, for an orthotropic
  behaviour, the rotations of the gradients, of the thermodynamic
  forces and of all the blocks of the tangent operator.
- **Description symbols**, read at load time: interface
  (`GenericParallel`), backend, API version and memory space. The
  caller can check everything when loading, without compiling against
  TFEL.
- **Versioned C data contract** (`MFront/GenericParallelBehaviour/BatchData.h`).
  A batch is described by a structure. Each group of variables is a
  field: a pointer, a stride between points and a stride between
  components. Two strides cover both useful layouts: point by point
  (the layout of MGIS) and component by component (coalesced accesses
  on the GPU). Material properties and external variables can also be
  given variable by variable, each one uniform or given per point.
- **Work at one point written once**, as templates usable on the host
  and on the GPU (`MFront/GenericParallelBehaviour/Integrate.hxx`). The
  generated behaviour class is the same for both interfaces: the
  `generic-parallel` interface only adds a constructor that receives
  the parameters. The integration drivers of the `generic` interface
  were split into a host and device core and the historical host
  functions, so both interfaces share the same cores.
- **Six backends**, each with one primitive: apply a function object to
  the points of the batch and reduce the statistics.

| Backend | Primitive | Execution queue | Asynchronous |
|---|---|---|---|
| `cuda`, `hip` | one thread per point, warp vote, atomics only on an event | stream | yes |
| `stlpar` | `std::transform_reduce` over a counting iterator | none | no |
| `kokkos` | `Kokkos::parallel_for`, atomics only on an event | execution space instance | yes |
| `raja` | `RAJA::forall`, atomics only on an event | resource of camp | yes |
| `sycl` | `parallel_for` and atomic updates | `sycl::queue`, required | yes, if no value is copied and the queue is in order |
| `openmp` | `target teams distribute parallel for` or `parallel for` | none | no |
| `openacc` | `parallel loop` with reductions | none | no |

- **Execution options** with the same meaning for all backends: queue,
  statistics of an asynchronous computation, and a flag that forbids
  changing the stack size of GPU threads. An empty batch is accepted
  without any check.
- **Memory space symbol** (`host`, `device` or `any`), fixed at compile
  time from the backend, the compiler and the memory mode.
- **Statuses and errors**: success, success with a time step factor
  below 0.99 at some points, failure at some points (with the number of
  failing points and the index of the first one, independent of the
  order of computation), invalid arguments, backend error. An optional
  array gives the status of each point.
- **Error reporting at no cost on success**: kernels write into the
  statistics only when an event occurs at a point. In synchronous mode,
  the `cuda`, `hip`, `sycl` and `kokkos` backends keep the statistics in
  reused blocks accessible by host and GPU. Measured on an elastic
  behaviour and 1000 points: the synchronous CUDA call goes from 135 to
  29 µs.
- **Synchronous or asynchronous computation**: with statistics given in
  the options, the `cuda`, `hip`, `kokkos` and `sycl` backends return as
  soon as the computation is submitted; after synchronising the queue,
  `mfront_gpb_report` fills the output status. Several batches and the
  assembly can be submitted before a single synchronisation.
- **Stack size of GPU threads**, computed from the attributes of the
  kernel (CUDA, HIP), from a probe kernel (Kokkos), or set through the
  CUDA runtime (`nvc++`, AdaptiveCpp).
- **Coverage**: everything the `generic` interface can do. Small strain
  and general behaviours (thermal, cohesive zones), finite strain with
  the Green-Lagrange and Hencky strategies, initialize functions and
  post-processings, rotations of orthotropic behaviours, single crystals,
  behaviour variables and their parameters, the
  `StandardElastoViscoPlasticity` brick, the TDLS linear solver in
  implicit behaviours, shared data (`@SharedData`), and standalone
  models, treated as behaviours without gradients.

### The `generic-parallel` interface for material properties

This interface already existed in TFEL. The prototype extends it to
offer the same possibilities as the interface for behaviours:

- `kokkos`, `raja`, `openmp` and `openacc` backends, added to `cuda`,
  `hip`, `stlpar` and `sycl`;
- properties defined by data or calling other laws, usable on CUDA and
  HIP;
- one function per property, `<name>_generic_parallel`, taking the data
  of a batch (version 3 of the functions): number of points, output and
  arguments with a stride, out-of-bounds policy, optional statuses of the
  points, and execution options with the same members and meaning as for
  behaviours. The signature is the same for all backends;
- the same evaluation at each point as the `generic` interface: same
  checks, same status per point, `NaN` where the evaluation fails;
- output status and statistics modelled on behaviours (numbers of
  failing points and of points out of bounds with the `Warning` policy,
  first point of each kind), merged only on an event;
- errors of the batch shared by the three kinds of objects: `-5`
  (invalid arguments), `-6` (invalid parameters), `-7` (runtime error of
  the backend);
- asynchronous evaluation, with an exported report function
  `<name>_generic_parallel_report`;
- tests written once for all the backends, and a comparison with the
  `generic` interface at each point;
- symbols for backend, function version and memory space;
- statistics of synchronous evaluations taken from a pool of reused
  blocks, and waiting on the stream or Kokkos instance used instead of
  the whole device;
- number of points of type `size_t` in CUDA, HIP and Kokkos kernels and
  with AdaptiveCpp, at most 2^48 - 1 points (the limit of 2^31 points is
  lifted).

### GPU code in TFEL and MFront

Code running on a GPU has no exceptions, no virtual calls on objects
created by the host, and no access to the global variables of the host.
The following changes make TFEL and the generated code compatible
without changing CPU results:

- host and device annotations of the mathematical and material headers
  (`TFEL_HOST_DEVICE`), checked by compiling each header with clang
  CUDA (`scripts/clang-cuda-header-scan.sh`);
- `TFEL_DEVICE_CODE`, defined when device code is compiled, including
  by single pass compilers (`nvc++`, AdaptiveCpp);
- behaviours without virtual functions in `MechanicalBehaviour`;
- `TrivialGenTypeBase` instead of `std::variant` for the tangent operator
  of finite strain behaviours;
- constant data built at compile time (slip systems) through
  `tfel::utilities::getConstantObject`;
- eigenvalues without recursion nor assertion on the GPU (Jacobi
  algorithm for multiple eigenvalues);
- portable errors: `tfel::fail<Exception>(message)`, `abortResolution`
  and `abortIntegration`. The 36 test behaviours whose user code was not
  usable on the GPU were made portable with these tools, without
  changing their results;
- shared data (`@SharedData`), computed on the host and copied once per
  GPU and per context;
- parameters of behaviour variables passed by value to the kernels;
- material laws called by behaviours, generated in the namespace
  `mfront::device`;
- static variables written `TFEL_STATIC_VARIABLE(v)`: the variable
  itself in host code, a copy initialised by a constant expression in
  device code.

### Single precision

- DSL option `real_type` (`double` by default, or `float`), and the
  macros `MFRONT_GB_REAL_TYPE` and `MFRONT_GMP_REAL_TYPE` written at the
  top of the sources of the `generic` and `generic-parallel`
  interfaces.
- Symbol `<f>_real_size` exported by each function, and API versions
  shifted by 1000 for a type other than `double`, so that an older
  caller refuses the library instead of reading it wrongly.
- In `float`, the criterion given by `@Epsilon` is bounded to ε/10 by
  MFront at generation time, with a warning, and
  `numerical_jacobian_epsilon` to ε.
- MTest refuses a behaviour, a material property or a `generic`
  evolution that is not in `double`.
- The other interfaces and the models of the DSL `Model` refuse
  `float`; models written with a behaviour DSL (`RungeKuttaModelDSL`)
  are generated in `float`.
- `IsAssignableTo` accepts a floating point type into a narrower one
  (decision D3 (a)).
- In `float`, MFront warns when the criterion of a behaviour
  (`@Epsilon`) is looser than 1e-5.
- **Displacement gradient** (`docs/displacement-gradient.md`): the caller
  may give H = F - I instead of F to a finite strain behaviour (stress
  measure offset by 10 in the `K` array, symbol
  `<f>_AcceptsDisplacementGradient`). The Green-Lagrange and Hencky
  strains are then computed from H without cancellation, which keeps the
  relative precision of small strains in `float`.
- The `LogarithmicStrainHandler` computes the divided differences of the
  logarithm with `log1p`, without cancellation when two eigenvalues are
  close: in `float`, the error on the derivative of the second
  Piola-Kirchhoff stress drops from 4.6e2 to 4e-6 in the worst case
  measured.

### Configurations

TFEL is built and tested with 17 configurations of the behaviour
interface (`build-tfel.sh`): CUDA with nvcc and clang, HIP with hipcc
and clang (NVIDIA platform of HIP), `stlpar` with g++, clang, icpx and
`nvc++` on the GPU, Kokkos, RAJA (nvcc, GPU and OpenMP on the host),
SYCL with AdaptiveCpp (CPU and GPU) and icpx, OpenMP on the host with g++,
OpenMP offloading with `nvc++` and OpenACC with `nvc++` (GPU and cores of
the host).
The `hip` backend also compiles behaviours for an AMD GPU (gfx90a,
compilation and linking only).

## MGIS

- **Loading.** MGIS recognises a `generic-parallel` library, checks the
  version of the data contract, and reads the backend and the memory
  space. It keeps its own copy of the C structures of TFEL, with a test
  that checks their memory layout is identical.
- **Automatic switch.** The existing functions (`integrate`,
  `executeInitializeFunction`, `executePostProcessing`, rotations)
  switch to batch computation by themselves. The calling code does not
  change. Data are never copied: the batch is described directly on the
  arrays of the `MaterialDataManager`.
- **Models** of the `generic-parallel` interface are used exactly like
  behaviours.
- **Material properties**: loading of `generic-parallel` properties
  (version 3 of the functions, two independent tests: functions of
  batches if the library exports their version, function of one point
  for the `generic` interface) and an `evaluate` function on a batch of
  points, each input uniform or given per point.
- **Memory.** The caller declares where its arrays are
  (`HOST`, `DEVICE` or `SHARED`). MGIS compares this location with the
  memory space of the library, so it needs no list of backends. Functions
  that read the arrays on the host refuse arrays placed on a GPU instead
  of crashing. A `MaterialDataManager` can allocate its arrays with a
  standard memory resource (`std::pmr::memory_resource*`). MGIS provides
  four resources as headers compiled by the caller: CUDA managed memory,
  HIP managed memory, SYCL shared memory and `Kokkos::SharedSpace`.
- **Execution policy** (`mgis::BatchExecutionPolicy`), common to
  behaviours, models and material properties: queue, and protection of
  the GPU stack size. Synchronous by default.
- **Asynchronous computation** (`mgis::AsynchronousBatchExecutionPolicy`)
  for every batch function. It returns a pending computation
  (`PendingBatchComputation`, `PendingBatchEvaluation`), whose
  `getResult` gives the result once the caller has synchronised the
  queue.
- **Views and layouts.** `FieldView` (values with a stride between
  points and a stride between components), a layout per material data
  manager (`POINT_BY_POINT` or `COMPONENT_BY_COMPONENT`), and
  `MaterialDataView`, the data of a material in the arrays of the
  caller. All functions of the manager exist on a view, with a single
  code path.
- **End of step on the GPU.** `update`, `revert` and `updateValues` take
  copy operations (`FieldOperations`), with headers for CUDA, HIP, SYCL,
  Kokkos, RAJA, OpenMP, OpenACC and the parallel algorithms of C++ (on the
  host, or on the GPU with `nvc++ -stdpar=gpu`) compiled by the caller.
  `make_unchecked_view` and `updateInputValues` let a caller such as
  mfem-mgis keep the data of its materials on the GPU.
- **No heap allocation** in the usual synchronous calls once the arrays
  are allocated, checked by a test that counts the calls to
  `operator new`.
- **Bindings**: the layout is exposed in C, Fortran and Python.
- **Profiling**: `setProfilingSynchronization` gives a function called
  at the start and end of each timed section, so that asynchronous GPU
  work is counted in the section that launches it.
- **Single precision**: CMake option `MGIS_REAL_TYPE`, generated header
  `MGIS/RealType.h`, and a check of the real type of each entry point
  at load time.
- **Displacement gradient**: option `gradient` of
  `FiniteStrainBehaviourOptions` and its bindings; the options of the
  behaviour are now given to the initialize functions and
  post-processings.
- **MGIS has no dependency on any GPU library.** The GPU tests need an
  NVIDIA GPU and are disabled by default
  (`enable-generic-parallel-cuda-tests` and
  `enable-generic-parallel-stlpar-nvhpc-gpu-tests`).

## mfem-mgis

The work is organised in work packages L0 to L10
(`docs/plan-mfem-mgis-detailed.md`), all done. The default path of
mfem-mgis does not change: each new path is chosen by a parameter or
an environment variable.

- **Material data in the memory of MFEM** (L1). The arrays of the
  materials are `mfem::Vector` objects, given to MGIS with their
  location. Host accesses through the public API synchronise with the
  host first. Elements and offsets of the quadrature spaces are readable
  on the GPU.
- **Batch integration** of behaviours (L2), on the host and on the GPU:
  one MGIS call per material, rotations of orthotropic behaviours,
  `update` and `revert` on the GPU, point-wise models (`PointWiseModel`).
  Parameter `BehaviourIntegration` (`PointByPoint`, `Batch`) or variable
  `MFEM_MGIS_BEHAVIOUR_INTEGRATION`.
- **Mechanics by kernels** (L3, L4), written once with `mfem::forall`
  and used on CPU and GPU: element geometry (Jacobians and their
  inverses computed once), gradients, internal forces, stiffness
  matrices. Parameters `GradientsComputation`, `InnerForcesComputation`
  and `TangentOperatorComputation` (`ElementByElement`, `Kernels`, and
  `MatrixFree` for the tangent).
- **Direct assembly into the sparse matrix**: pattern built from the
  elements, positions of the entries computed once, exchanges kept, by
  atomic additions or in a fixed order that gives identical results from
  one run to the next (parameter `TangentOperatorAssembly`, `Atomic` or
  `Ordered`).
- **Matrix free tangent operator**: action and diagonal without
  assembly, with a preconditioner built on the assembled operator
  (`AssembledOperator`) or matrix free smoothers
  (`OperatorJacobiSmoother`, `OperatorChebyshevSmoother`).
- **Prediction by kernels** (L4), and the integration rule order chosen
  by the user (`IntegrationRuleOrder`: with order 3, 8 points instead of
  27 on first order hexahedra).
- **Other integrators and boundary conditions**: boundary conditions
  that do not depend on the unknowns (imposed pressure, heat source) are
  assembled once for a unit value, then scaled and added on the GPU to
  the internal forces computed by kernels.
- **Other formulations by kernels** (L7): stationary and transient
  thermal problems, micromorphic damage, FBar, Faltus et al.
  regularisation.
- `TransientHeatTransfer` also in plane strain and plane stress (a
  feature without GPU code).
- **Post-processings by kernels** (L8): integrals, L2 projections,
  nodal values (variable `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION`).
- **Solvers** (L5):
  - `FixedPoint` nonlinear solver with the methods `ElasticOperator`,
    `ModifiedNewton` and `NewtonK0`, and all accelerations ported from
    MTest (Anderson, Cast3M, secant methods, Irons-Tuck, Steffensen,
    delta squared methods);
  - `NewtonK0`: consistent tangent, Krylov solver preconditioned by the
    solver of the elastic operator K0 built once;
  - inexact Newton (`AdaptiveLinearSolverTolerance`, Eisenstat and
    Walker), for Newton and `NewtonK0`;
  - BoomerAMG options (strength threshold, coarsening, aggressive
    coarsening, interpolation, relaxation, cycles, levels), reapplied
    after each new operator (`ConfiguredHypreBoomerAMG`);
  - reuse of the MUMPS analysis while the pattern of the matrix does not
    change on any process (`ReuseAnalysis`);
  - cuDSS direct solver on the GPU (`CuDSSSolver`), with the same reuse
    of the analysis;
  - hypre on the host while the device is on the GPU
    (`--hypre-on-host`), for small linear systems;
  - optional AmgX solvers (`enable-amgx`), using the C interface of
    AmgX directly: `AmgX` preconditioner and solver, `AmgXPCG`;
  - optional Kokkos initialised by mfem-mgis (`enable-kokkos`), for
    behaviours compiled with the Kokkos backend.
- **Solver libraries** (15 and 16/10/2026, `docs/plan-solver-libraries.md`
  and `docs/solver-recommendations.md`): the linear solvers and
  preconditioners are registered by one module per library, with a
  description (where they compute, kind, accepted preconditioners,
  options) and a generic test:
  - PETSc (Krylov solvers and preconditioners by the options of PETSc,
    GAMG with the rigid body modes, SNES), Ginkgo (solvers and
    preconditioners described by configurations of Ginkgo), the solvers
    of hypre not wrapped by MFEM (BiCGSTAB, COGMRES, LGMRES, hybrid
    solver, FSAI, MGR, smoothers), SuperLU_DIST, STRUMPACK, KLU, the
    preconditioners of the sparse matrices of MFEM, the options of MUMPS,
    KINSOL, oneMKL PARDISO and Cluster Sparse Solver, the sparse solvers
    of MAGMA;
  - Trilinos (optional, `enable-trilinos`): Krylov solvers of Belos,
    direct solvers of Amesos2 (Tacho on the GPU), Ifpack2, MueLu with the
    rigid body modes, also as an approximate solver (one V-cycle), with
    the reuse of its hierarchy between operators, FROSch. They compute
    where MFEM computes; on the GPU, the matrices of hypre are converted
    by kernels of Kokkos.
- **Solver overrides** (16/10/2026), to compare solvers on existing
  programs without changing them: `--linear-solver-override` (`-lso`),
  `--preconditioner-override` (`-pco`) and `--non-linear-solver-override`
  (`-nlso`), or the variables `MFEM_MGIS_LINEAR_SOLVER_OVERRIDE`,
  `MFEM_MGIS_PRECONDITIONER_OVERRIDE` and
  `MFEM_MGIS_NON_LINEAR_SOLVER_OVERRIDE`, replace the linear solvers, the
  preconditioners and the nonlinear solvers set by the program. Each one
  is described as a factory argument (a name, or `Name` and `Options`) in
  the syntax of the data of TFEL. The options of the nonlinear solver
  override take precedence over the parameters of the program; the
  parameters it refuses are ignored, except the convergence criteria.
- **GPU build and configuration**: CUDA enabled when MFEM is compiled
  with CUDA, and the sources containing kernels compiled as CUDA, or as
  HIP when MFEM uses HIP; `--device` option of `initialize`; hypre
  memory pool on the GPU; standard host memory for the device configured
  by `initialize`; profiling sections that wait for the launched kernels.
- **Execution settings** (`docs/plan-execution.md`):
  - `--execution cpu|gpu|mixed`, and, in `mixed` mode, the placement of
    each phase: `--finite-elements-execution`, `--behaviours-execution`,
    `--linear-solvers-execution`, `--post-processings-execution`;
  - the placement of the behaviour of each material (parameter
    `BehaviourExecutionLocation` of `addBehaviourIntegrator`);
  - counting of the copies of the material arrays
    (`getMaterialArraysTransfers`), shown after the profiling table;
  - the CMake function `mfem_mgis_add_behaviours_library`, which builds
    the `generic` library and, when MFEM uses a GPU, the library of the
    GPU backend, and `getBehavioursLibrary`, which chooses the library
    from the placement of the behaviours.
- **Single precision** (`docs/plan-precision.md`): the real type is the
  one of the installed MFEM. Configuration checks the type of MFEM and
  of MGIS and refuses AmgX in float; parameters are stored as `double`;
  MPI types follow `real`; relative tolerances of the solvers are
  bounded; in float, the GPU stiffness kernels are compiled without
  floating point contraction (CMake function
  `mfem_mgis_disable_device_floating_point_contraction`). The
  regularisation of the Anderson acceleration is raised to 10 ε in
  float. The tests are adapted to float (tolerances of
  `tests/Tolerance.hxx`, from the measured errors): 223/223 in float on
  the CPU and 280/280 on the GPU. Ten tests stay in double only: UMFPACK,
  GSLIB, Umpire, cuDSS and CUDA, absent from MFEM in float. ParaSails
  runs in float since hypre is built with its internal LAPACK (see the
  defect of hypre in single precision in "Defects found and corrected").
  The Newton solver accepts a stagnation of the residual at its floor
  (option `StagnationRelativeTolerance`, 10⁴ ε by default): the floor of
  the relative residual grows with the size of the problem (1.3e-5 for
  the beam r4 on the CPU).
- **Whole test suite on the GPU** (variable `MFEM_MGIS_DEVICE`), except
  the two tests that require the absence of a GPU.
- **Displacement gradient**: parameter `DisplacementGradient` of the
  behaviour integrators of finite strain behaviours; the macroscopic
  gradients are then displacement gradients. A test checks that the
  Newton solver converges in `float` below the floor reached with the
  deformation gradient.
- **Documentation** in the tree of mfem-mgis, for example
  `docs/source/user_guide/batch_computations.rst`.

## mfem-mgis-examples and mm-opera-hpc

- **mfem-mgis-examples** (38 patches): options of the examples
  (behaviour library, device, solvers, profiling), synchronisations on
  the host for the GPU, choice of the direct solvers, behaviour
  libraries for the CPU and the GPU chosen by `--execution`, tests on
  the GPU. Tests on the GPU are registered only when the GPU behaviour
  library is built.
- **mm-opera-hpc** (40 patches): for the polycrystal, choice of the
  nonlinear solver (Newton, NewtonK0), of the adaptive tolerance and of
  AmgX, behaviour integrated on the GPU; for the bubble and the cermet,
  options of mfem-mgis and of the behaviour, synchronisations on the
  host, profiling of the post-processings, direct solvers, behaviour
  libraries chosen by `--execution`, tests on the GPU, and cases in
  single precision. Tests are compared with reference solutions (Cast3M,
  analytical solution of the bubble).
- All examples and all cases (polycrystal, cermet, bubble) run by
  batches and on the GPU.
- The cermet and the polycrystal accept `--displacement-gradient`; their
  macroscopic unknowns are displacement gradients, which removes the
  quantization of the macroscopic loading in `float`.

## Prototype repository

- Build scripts: `build-tfel.sh`, `build-mgis.sh`, `build-mfem-mgis.sh`,
  and `scripts/build-hip-amd.sh` for an AMD GPU (compilation and
  linking only).
- Measurement cases in `cases/` (beam, thermal, polycrystal, bubble,
  cermet, ssna303, examples) and their results in `results/`. Their
  test suite also runs in float: an executable of the examples or of
  mm-opera-hpc is only used if it is built against the same installation
  of mfem-mgis as the cases, the tolerances of the comparisons depend on
  the real type, and the UMFPACK, AmgX and ex5 tests are absent in float
  (134 tests in float, 143 in double).
- The website (`site/`, `scripts/build-site.sh`,
  `.github/workflows/pages.yml`) and the continuous integration
  (`.github/workflows/ci.yml`: tests on the CPU in double and float,
  compilation with CUDA and with HIP), both started at each push on main
  since the repository became public (15/10/2026).
- Tools in `scripts/`: patch synchronisation and checking, test matrix,
  test suite of mfem-mgis on the GPU (`scripts/mfem-mgis-gpu-suite.sh`),
  host and device annotation checks, compile and run time sweeps,
  rounding checks, and a memory guard for long compilations.
- Small standalone reproducers of specific points in `reproducers/`,
  and compiler defects with their reproducers in `results/nvhpc-bugs/`.
- Spack environments in `deps/spack/`: `mmm-gpu`; `mmm-float-cpu` and
  its double precision twin `mmm-double-cpu`; `mmm-float-gpu`; `ci-cpu`
  and `ci-cpu-float` for the continuous integration.
