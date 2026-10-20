# Upstream integration

This document proposes a split of the changes into contributions for the
official repositories of TFEL, MGIS and mfem-mgis, then of the examples and
of mm-opera-hpc. Each contribution must be reviewable and testable on its
own. The fixes already proposed to MFEM and to AmgX are listed at the end.

Some commits of the prototype mix several topics. The contributions will
therefore be rebuilt from the final state of the branches, file by file,
rather than by selecting commits.

## TFEL

The order follows the dependencies: each contribution builds on the
previous ones. The technical details are in [tfel.md](tfel.md).

1. **Fixes of the existing code**, without GPU code. These are the fixes
   listed in [tfel.md](tfel.md#fixes-of-the-existing-code). They can
   change the results on CPU, since they are fixes. Each one is covered by
   an existing or new test.
2. **GPU code in TFEL** (mathematical and material headers): host and GPU
   annotations, `TFEL_DEVICE_CODE` for the single pass compilers,
   behaviours without virtual functions, `TrivialGenTypeBase`,
   `getConstantObject`, eigenvalues without recursion, compiler workarounds
   (for example the iterators of the cubic splines). The results on CPU do
   not change. The incompatibilities are listed in the release notes.
3. **GPU compatible generated code in MFront**: annotations of the
   generated classes, GPU compatible behaviour class, drivers of the
   `generic` interface split into a host and GPU core and the historical
   functions, rotation functions and method calls usable on GPU, type of
   the tangent operator checked before its conversion, parameters passed to
   the behaviour class and to its behaviour variables, portable errors
   (`tfel::fail`, `abortResolution`, `abortIntegration`), shared data
   (`@SharedData`), static variables written `TFEL_STATIC_VARIABLE(v)`
   (copy in the GPU code, the variable itself on the host). The results of
   the existing interfaces do not change.
4. **`generic-parallel` interface of the behaviours and of the models**,
   with a first set of backends: CUDA, HIP, parallel algorithms of C++ and
   OpenMP on the host. Versioned data contract, execution options,
   synchronous or asynchronous computation, errors without cost in case of
   success, memory space symbol, shared data copied on each GPU, tests and
   documentation.
5. **Additional backends and compilers**: Kokkos, RAJA, SYCL, offload of
   `nvc++` (`stdpar` and OpenMP), clang CUDA and HIP. Most of the compiler
   workarounds go into this contribution.
6. **Extensions of the interface of the material properties**: Kokkos, RAJA
   and OpenMP backends, properties defined by data or calling other laws on
   CUDA and HIP, execution options as last argument with the same meaning
   for all the backends (version 2 of the functions), asynchronous
   evaluation and report function, symbols of backend, version and memory
   space, code of the kernels shared by the backends, statuses in reused
   blocks, status 1 with the `Warning` policy.
7. **Float precision** (DSL option `real_type`, patches 0112 to 0120):
   numerical type of the `generic` and `generic-parallel` interfaces,
   symbol `<f>_real_size` and shifted API versions, refusal of float by the
   other interfaces and by the models of the `Model` DSL, `IsAssignableTo`
   towards a narrower floating point type, MTest that checks the type,
   convergence criteria bounded at generation with warnings, documentation.
   In double, the generated code only changes by the size symbol.

The data contract of the behaviours
(`MFront/GenericParallelBehaviour/BatchData.h`) is versioned from the
fourth contribution on, so that MGIS can check it.

### Checks for each contribution

- the complete test suite of TFEL (`ctest`);
- the comparisons of the `generic-parallel` tests for the configurations
  of the contribution, at `-O3`;
- no change of the results of the existing interfaces, except the fixes of
  the first contribution;
- formatting by `clang-format` limited to the modified lines.

## MGIS

The contributions of MGIS follow the fourth contribution of TFEL, whose
contract they read. The details are in [mgis.md](mgis.md).

1. **Fixes**: axes of the rotation matrices, material property given by a
   `std::span` of one value, failure status of the initialisation
   functions and of the point by point post-processings.
2. **Batches of behaviours and models**: loading and checks, copies of the
   structures of TFEL and test of their layout, integration without copy,
   initialize functions, post-processings, rotations, automatic switch
   of the existing functions, execution policy, refusal of the one point
   functions.
3. **Memory**: location of the arrays checked against the memory space of
   the libraries, refusal of host accesses to the arrays placed on a GPU,
   memory resources and their headers.
4. **Material properties**: loading of the `generic-parallel` properties
   and function `evaluate`.
5. **Layouts and copy operations**: layout of the arrays of a manager
   (`MaterialDataManagerInitializer::layout`), copy operations
   (`FieldOperations`), view without check of the values kept on the host
   (`make_unchecked_view`) and update of the inputs only of a state
   (`updateInputValues`), on which the data of the materials of mfem-mgis
   on the GPU depend, HDF5 files and C, Fortran and Python bindings that
   follow the layout.
6. **Profiling**: synchronisation function called around the measured
   sections, to wait for the kernels launched on a GPU.
7. **Documentation**: web page of the batches and release notes.
8. **Float precision** (option `MGIS_REAL_TYPE`, patches 0064 to 0070):
   type checked when loading each entry point, tests in float on the CPU
   and with CUDA, documentation. It follows contribution 7 of TFEL.

The GPU tests of MGIS need an NVIDIA GPU. They are disabled by default
(options `enable-generic-parallel-cuda-tests` and
`enable-generic-parallel-stlpar-nvhpc-gpu-tests`), so that the continuous
integration of MGIS stays on CPU.

## mfem-mgis

The branch `gpu-proto` of mfem-mgis starts from the official master
`20dcdf0a9b` (16/10/2026, merge of PR #487). The order follows the
dependencies:

- contributions 1 and 2 need neither the TFEL nor the MGIS of the
  prototype;
- contribution 3 needs the profiling synchronisation of MGIS (sixth
  contribution of MGIS);
- contribution 5 needs the location of the arrays and the layouts of MGIS
  (third and fifth contributions of MGIS);
- from contribution 6 on, the `generic-parallel` interface of TFEL (fourth
  contribution of TFEL) and the batches of MGIS (second contribution of
  MGIS).

The overall plan is in [plan-mfem-mgis.md](plan-mfem-mgis.md).

1. **Fixes of the existing code**, without GPU code:
   - ParaSails preconditioner given again to a new operator after the
     destruction of the previous one (defect of `mfem::HypreParaSails`);
   - settings of the ILU of hypre applied again after each new operator
     (MFEM resets them to the default values);
   - absolute tolerance of the Krylov solvers of hypre (refused, or given
     as relative tolerance to FGMRES);
   - short names of the options declared by `declareDefaultOptions`
     (`-up`, `-no-up`, `petscrc_file`) ignored by `initialize`;
   - exception thrown by the tangent operator during the Newton
     correction: failure of the resolution in sequential, message then
     `MPI_Abort` in parallel, instead of `std::terminate`;
   - exceptions of the L2 projections, declared `noexcept`, reported in the
     execution context;
   - sequential prediction: diagonal set to one at the imposed degrees of
     freedom, as in parallel;
   - CMake: several directories returned by `tfel-config`, one line of
     `env.sh` per directory;
   - documentation of the constructors that take an execution context.
2. **Features without GPU**, useful on CPU:
   - order of the integration rules chosen by the user
     (`IntegrationRuleOrder`: with order 3, 8 points instead of 27 on the
     hexahedra of order 1);
   - BoomerAMG options (strength threshold, coarsening, interpolation,
     relaxation, cycles), applied again after each new operator;
   - inexact Newton (tolerance of the linear solver adapted to the
     residual, Eisenstat and Walker);
   - Jacobi and Chebyshev smoothers of MFEM as preconditioners;
   - analysis of MUMPS reused as long as the pattern of the matrix does not
     change on any process (`ReuseAnalysis`, number of analyses given by
     `getNumberOfAnalyses`, for MUMPS and cuDSS; test
     `DirectSolverReuseAnalysisTest`);
   - parameters of the problem passed by the constructors of
     `PeriodicNonLinearEvolutionProblem`;
   - `TransientHeatTransfer` in plane strain and in plane stress;
   - profiling sections (setup of the linear solver, prediction, end of the
     time step, post-processings).
3. **Build and configuration for the GPU**:
   - MFEM compiled with CUDA, sources containing kernels of MFEM compiled
     as CUDA code, or as HIP code when MFEM uses HIP (AMD platform,
     architectures of MFEM; compilation checked, not run on an AMD GPU);
   - `FindMFEM`: `MFEM_DIR` takes precedence over `CMAKE_PREFIX_PATH`;
   - option `--device` of `initialize` (one GPU per process, device already
     configured by the program, device created before PETSc);
   - memory pool of hypre on the GPU;
   - standard host memory for the device configured by `initialize`: with
     MFEM compiled with Umpire, the temporary vectors of the host loops on
     the elements went through Umpire (2 to 5 times slower);
   - hypre on the host with the device on the GPU (`--hypre-on-host`), for
     the small linear systems;
   - cuDSS direct solver on the GPU (`CuDSSSolver`, wrapper of
     `mfem::CuDSSSolver` that recreates it when the size or the pattern
     changes), with the reuse of the analysis of contribution 2;
   - profiling sections that wait for the launched kernels.
4. **GPU path of the existing code**: the whole test suite on the GPU
   (variable `MFEM_MGIS_DEVICE`), and the fixes found this way:
   - host accesses declared to MFEM (`HostRead`, `HostReadWrite`):
     components of the L2 projections, imposed values of the Dirichlet
     boundary conditions, ParaView export, interpolation of the fields,
     analytical L2 error;
   - ILU of hypre limited to level 0 and ParaSails refused when hypre
     computes on the GPU (failures of hypre).
5. **Data of the materials in the memory of MFEM**: arrays of the
   integration points stored component by component, synchronised on the
   host before each host access, elements and offsets of the quadrature
   spaces readable on the GPU, copies of the states by kernels.
6. **Batch integration of the behaviours**, on the host and on the GPU:
   batches of MGIS, rotations of the orthotropic behaviours, `update` and
   `revert` on the GPU, point by point models (`PointWiseModel`). Without
   the option `--execution` (contribution 11), point by point integration
   remains the default; the parameter `BehaviourIntegration` (or the
   variable `MFEM_MGIS_BEHAVIOUR_INTEGRATION`) selects the batches.
7. **Mechanics by kernels**:
   - geometry of the elements (jacobians and their inverses computed
     once), gradients, internal forces, stiffness matrices;
   - direct assembly in the sparse matrix (pattern built from the
     elements, positions of the entries computed once, exchanges kept), by
     atomic additions or in a fixed order, which makes the results
     identical from one pass to the next (parameter
     `TangentOperatorAssembly`, variable
     `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`);
   - prediction, action and diagonal of the tangent operator without
     assembly (`MatrixFree`), preconditioner on the assembled operator
     (`AssembledOperator`);
   - other integrators of the form (Robin, springs), boundary conditions
     independent of the unknowns (imposed pressure, heat source) assembled
     once and added to the internal forces computed by kernels.

   Without the option `--execution`, element by element computation
   remains the default; the parameters `GradientsComputation`,
   `InnerForcesComputation` and `TangentOperatorComputation` (or the
   corresponding variables) select the kernels.
8. **Other formulations by kernels**: stationary and transient thermal,
   micromorphic damage, FBar, regularisation of Faltus et al.
9. **Fixed point solver** (`FixedPoint`): accelerations of MTest, NewtonK0
   method (factorised elastic operator, adapted tolerance), elastic
   prediction, iterative solvers of the factory.
10. **Post-processings by kernels** (variable
    `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION`): integrals, L2 projections,
    values at the nodes.
11. **Execution settings** (see [plan-execution.md](plan-execution.md)):
    option `--execution cpu|gpu|mixed` and placement of each phase in
    `mixed` mode (`--finite-elements-execution`, `--behaviours-execution`,
    `--linear-solvers-execution`, `--post-processings-execution`), default
    values deduced from the placement (batches and kernels, except with PETSc),
    placement variables ignored with a warning, placement of the behaviour
    of a material (parameter `BehaviourExecutionLocation`), copies of the
    arrays of the materials counted, libraries of behaviours for the CPU
    and the GPU (`mfem_mgis_add_behaviours_library`,
    `getBehavioursLibrary`). It depends on the batches, the kernels and the
    post-processings.
12. **Float precision** (see [plan-precision.md](plan-precision.md)): real
    type given by MFEM, configuration (refusal of an MGIS of another type
    and of AmgX in float), parameters stored as `double`, MPI types,
    tolerances of the tests adapted to float, Anderson regularisation
    raised to 10 ε, and, in float, stiffness kernels of the GPU compiled
    without contraction into fused multiply-adds. It follows the precision
    contributions of TFEL and MGIS.
13. **Optional, to be proposed later**:
    - Kokkos initialised by mfem-mgis (option `enable-kokkos`), for the
      behaviours compiled with the Kokkos backend;
    - AmgX solvers, which use the C interface of AmgX directly: after the
      merge of the fix proposed to AmgX and of its spack recipe.

Many commits of `gpu-proto` touch several contributions. For example, the
move of the whole suite to the GPU carries fixes (contribution 1),
configuration (3) and the GPU path (4), and `LinearSolverFactory.cxx`
carries pieces of six contributions. As for TFEL, the contributions will
be rebuilt from the final state, by pieces of files.

### Checks for each contribution

- the complete test suite of mfem-mgis on CPU, in sequential and on two
  processes;
- from contribution 4 on, the suite on the GPU (`MFEM_MGIS_DEVICE=cuda`);
- from contribution 6 on, the suite with the variables that force the
  batches and the kernels (`MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch`,
  `MFEM_MGIS_*_COMPUTATION=Kernels`), on CPU and on GPU;
- the examples of mfem-mgis-examples and the cases of mm-opera-hpc
  followed by the prototype, with the same results on CPU and on GPU up to
  roundoff;
- a test for each fix and each feature (some do not have one yet, for
  example the parameters of `PeriodicNonLinearEvolutionProblem`);
- formatting by `clang-format` limited to the modified lines (the
  repository is not formatted as a whole).

## mfem-mgis-examples and mm-opera-hpc

The changes of the examples (31 patches since `600da70`) and of
mm-opera-hpc (34 patches since `49147f2`) follow the contributions of
mfem-mgis that they use. They will be proposed after them, by topic:

- options of mfem-mgis and of the behaviour in each example and each case
  (library of behaviours, device, solvers, profiling), then `--execution`
  and the libraries of behaviours for the CPU and the GPU (contribution
  11);
- synchronisations on the host before the host reads (contribution 5);
- choice of the non linear solver (Newton, NewtonK0) and of the adapted
  tolerance (contributions 2 and 9), direct solvers and reuse of the
  analysis (contributions 2 and 3), AmgX (contribution 13);
- examples and cases in float: MPI types, typed options, adapted criteria
  and tolerances (contribution 12);
- tests on the GPU, registered only when the behaviour for the GPU is
  built;
- in mm-opera-hpc, tests compared to the reference solutions (Cast3M,
  analytical solution of the bubble) and fixes found by the review: help
  `-h` without failure, unit of the pressure of the bubble, output files
  opened after the options and written by the main process only.

## Already proposed upstream

Fixes of MFEM and AmgX proposed during the prototype to their official
repositories (open on 13/10/2026):

- MFEM #5604: shared resources of `AmgXSolver` (lifetime, sequential
  initialisation);
- MFEM #5605: order of the vertices of the boundary elements of a
  `SubMesh`, taken from the parent boundary element;
- MFEM #5615: `HypreParaSails::SetOperator` takes the communicator of the
  new operator. mfem-mgis works around this defect in contribution 1;
- AmgX #386: classical multigrid with CUDA 13 ("invalid device
  function", reproducer `reproducers/amgx-classical-cuda13`).

The rebuild of MFEM and AmgX in the environment waits for their merge (see
`deps/spack/README.md`).
