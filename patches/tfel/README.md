# Patches

Series of patches giving this tree from TFEL master, commit `4d50821c9`
(23/10/2026):

```bash
git clone <official TFEL repository> tfel && cd tfel
git checkout 4d50821c9
git am --keep-cr /path/to/patches/*.patch
```

`--keep-cr` keeps the CRLF line endings of some files of TFEL (test
behaviours), without which some patches do not apply.

Patches, by subject (independently of their numbers; the details and
their justifications are in `docs/tfel.md` of the prototype):

- device code in the headers of TFEL: mathematical and material headers,
  finite strain tools, rotations, constant objects (slip systems),
  trivially copyable tensors, eigenvalues and eigenvectors without
  recursion nor assertion (Jacobi algorithm for multiple eigenvalues),
  tensor accessors, cubic splines, consistent host and device annotations
  (checked by compiling each header with clang CUDA), device code
  selection for single pass compilers, portable error reporting,
  workarounds of compiler defects;
- generated code compatible with device code: behaviours without virtual
  functions, integration drivers of the `generic` interface split into
  device capable cores, rotation functions, behaviour variables, models
  and porous flows of the StandardElastoViscoPlasticity brick, parameters
  of the behaviour variables, data shared by all the integration points
  (`@SharedData`), static variables written `TFEL_STATIC_VARIABLE(v)`
  (Runge-Kutta DSLs and material properties included), type of the
  tangent operator checked before its conversion;
- `generic-parallel` interface for behaviours and models: CUDA, HIP, C++
  standard parallelism, Kokkos, RAJA, SYCL, OpenMP (host and offloading)
  and OpenACC backends, standalone models, finite strain behaviours, initialize
  functions, post-processings, rotation of all the tangent operator
  blocks, versioned C data structures, symbols of the interface, of the
  backend and of the memory space, execution options shared by the
  backends, asynchronous execution, error reporting without cost on
  success, shared data copied on each device, statistics;
- `generic-parallel` interface for material properties: fixes, Kokkos, RAJA,
  OpenMP and OpenACC backends, properties using data or other laws,
  asynchronous evaluation, statistics, symbols, statuses in reused blocks,
  status of the `Warning` policy, physical bounds. Patch 0124 gives
  version 3 of the functions: one function per property taking the data
  of a batch, evaluation at each point identical to the `generic`
  interface (statuses of the points, `NaN` where it fails), output status
  and statistics modelled on the behaviours, errors of the batch `-5`,
  `-6` and `-7`. Patch 0125 writes the tests once for all the backends
  and compares them with the `generic` interface at each point. Patch
  0123: the device implementation of a material law returns `NaN` for a
  result which is not finite;
- `generic-parallel` interface for behaviours (patches 0121 and 0122):
  models executed when a prediction is computed, one definition of the
  time step scaling factor of the statistics;
- tests of both interfaces against the `generic` interface for all the
  configurations (TDLS, isotropic DSLs, parameters changed at run time,
  static variables), and documentation (web pages, keyword pages,
  incompatibilities, effect of the single pass compilers);
- build and configuration files: linker flags after the object files,
  headers of HIP given to hipcc, configuration files of the backends and
  of the compilers (`nvc++ -stdpar=gpu` without relocatable device code,
  optimisation flags);
- fixes of the existing code, which may change the results on CPU: bounds
  checks of the generic models and messages of the bounds of the material
  properties, Poisson ratio of the stress update algorithm of the
  isotropic DSLs (and the test Plasticity8, which depends on it), DDIF2
  brick (non constant fracture stresses and energies, checks of the
  options), error normalization factor of `rkCastem`. Patch 0132: line of
  the errors of the files of parameters, `errno` restored after reading
  the file of parameters of a material property, option `shall_fail` of
  `tfel-check` refusing a command which succeeds;
- single precision (patches 0105 to 0113): DSL option `real_type`
  (numerical type of the `generic` and `generic-parallel` interfaces,
  symbols `<f>_real_size`, API versions shifted for `float`), other
  interfaces refusing `float`, `IsAssignableTo` for floating point types,
  MTest checking the numerical type, convergence criteria bounded at the
  generation with warnings, documentation;
- displacement gradient (patches 0114 to 0118): Green-Lagrange tensor
  computed from the displacement gradient, `LogarithmicStrainHandler`
  built from the Green-Lagrange tensor and divided differences without
  cancellation, displacement gradient given to the finite strain
  behaviours by the `generic` interfaces, axial deformation gradient of
  plane stress computed by the behaviour;
- RAJA backend (patch 0119): execution policy given as an option and
  checked before being copied in the generated code, resource of camp as
  queue, `raja` and `raja-host` configurations, tests;
- OpenACC backend (patch 0120): `parallel loop` with reductions, memory
  space deduced from the compilation macros, stack size set through the
  CUDA driver, `openacc-nvhpc-gpu` and `openacc-nvhpc-multicore`
  configurations, tests;
- common contract of the behaviours and of the models (patches 0126 and
  0128): functions named after the ones of the `generic` interface
  followed by `_generic_parallel`, version 3 of the data structures,
  report function exported by each library, errors of the batch `-5`,
  `-6` and `-7`;
- interfaces built on the `generic` interfaces and `generic_functions`
  option (patches 0127 and 0129 to 0131): code evaluating a material
  property at a point written once, libraries also containing the
  functions of the `generic` interface, with one set of parameters and,
  for a behaviour, one out of bounds policy (release notes: patch 0133);
- workaround of a defect of nvc++ `-stdpar=gpu` (patch 0134): the
  execution of an initialize function at a point is not inlined in the
  reduction kernel;
- execution policies (patches 0135 to 0137): one C++ class per backend
  deriving from `AbstractBackendExecutionPolicy` (identifier checked by
  the generated functions), common options (preserve the stack size,
  asynchronous computation), `MFrontGenericParallelRuntime` library,
  statistics and execution policy in the data of the batches (version 4
  of the data structures), identifier of the expected policy exported by
  each library (initialized statically, except with Kokkos), tests of the
  rules for all the functions;
- `arcane` backend (patch 0138): accelerator API of Arccore, the core
  component of the Arcane framework, with a queue selecting the GPU or the
  cores of the host at runtime, native atomic operations, unified memory
  of Arccore, `arcane` configuration (Arccore built with CUDA, clang as a
  CUDA compiler), tests on CUDA, `Thread` and `Sequential` queues;
- documentation: the pages and comments of the changed code, checked
  against it.
