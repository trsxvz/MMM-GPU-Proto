# Design review of the branch generic-parallel-behaviours (TFEL)

> Historical document. Independent review made on 06/10, on the state of
> the branch at that date. All its findings have been handled (phase R1 of
> the [log](../../JOURNAL.md)). The current state is described in
> [tfel.md](../tfel.md).

Date: 2026-10-06. Diff reviewed: fec0ade99..generic-parallel-behaviours (40 commits, 364 files, +18421/-4875).

## FILES READ / REMAINING

READ:
- DESIGN.md, docs/web/generic-parallel-behaviour-interface.md
- GenericParallelBehaviour/BatchData.h, Integrate.hxx, GPU.hxx, ParallelSTL.hxx, Kokkos.hxx, OpenMP.hxx, SYCL.hxx, FiniteStrainIntegrators.hxx, ModelIntegrator.hxx
- diff of GenericBehaviour/Integrate.hxx
- GenericParallelBehaviourInterface.{hxx,cxx}, GenericParallelModelInterface.cxx
- GenericParallel/Behaviour/* (generator backends), DeviceStackSize.hxx
- diff MaterialProperty/* (ParallelSTL, BackendBase, CxxProgrammingModelBackendBase)
- diff AbstractBehaviourInterface, BehaviourCodeGeneratorBase (.hxx, .cxx for the most part), GenericModelInterface.hxx
- MechanicalBehaviour.hxx, TFELConfig.hxx
- config-*.json, test.check, compare-all.sh, compare-functions.sh, compare-tdls.sh, compare.cxx (structure)

- BehaviourRotationFunctions.hxx, CodeGeneratorUtilities.cxx (diff), MFrontMaterialPropertyInterface.cxx (diff, mfront::device)
- GenericBehaviour/StressMeasures.hxx, outline of GreenLagrangeStrainIntegrate.hxx
- include/TFEL: TFELConfig, ConstantObject, TrivialGenTypeBase, GenTypeBase, BoundsCheck, MechanicalBehaviour, GenericFixedSizeArray, TensorViewFromStensor, StensorComputeEigenVectors (commit 9cd913eb2), FiniteStrainBehaviourTangentOperator, tmatrix (change_basis), LinearHomogenizationSchemes, DDIF2Base (samples)
- cmake/modules/tfel.cmake, ImplicitCodeGeneratorBase (diff), PerformanceProfiling, RungeKuttaDSLBase, BehaviourDSLCommon (diffs)
- messages of the workaround commits (9cd913eb2, a22a3cff6, 58f4ced56, f008859df, b0d52da8f, 494929641, 1b621828d)

REMAINING: none (ranked summary written at the end of the file)

## FINDINGS (in reading order)


### Group 1: execution headers (GenericParallelBehaviour/*)

- [MODERATE, duplication] Each backend rewrites the same four entry
  functions (integrate, executeInitializeFunction, executePostProcessing,
  rotate), which only bind `makeApply`: GPU.hxx:280-364,
  ParallelSTL.hxx:103-203, OpenMP.hxx (same structure), Kokkos.hxx (same).
  About 100 lines per backend, documentation included. The generated code
  could call `mfront::gpb::integrateBatch<...>(s, d, p,
  policy, gpu::makeApply(...))` directly: each backend would only provide
  its `apply`/`makeApply`. The announced design (DESIGN.md, "each backend
  provides a single primitive") is right, but the implementation does not
  go all the way.
- [MODERATE, duplication] SYCL.hxx:369-581 keeps its own entry points
  (4 functions of 40 to 50 lines), each of which repeats: check, test of the
  queue, creation of DeviceBatchData, transfer, dispatch, copyOutputs. The
  only reason is the transfer of the fields. One could transfer first, then
  call `integrateBatch` with the transferred BatchData. Minor
  inconsistency: `std::unique_ptr<DeviceBatchData>` for three functions,
  `std::optional<DeviceCopies>` for `rotate` (SYCL.hxx:563).
- [MODERATE, duplication] Integrate.hxx: `dispatch`, `dispatchInitializeFunction`,
  `dispatchPostProcessing` (Integrate.hxx:1081-1092, 1241-1258, 1286-1303)
  repeat the same direct/copy switch. `PointIntegrator`,
  `InitializeFunctionExecutor`, `PostProcessingExecutor` repeat the same
  members (BatchData, InputFields, parameters, policy). A generic switch
  (`std::bool_constant<direct>`) would be enough.
- [MINOR, duplication] Integrate.hxx:110-210: `gather`/`gatherValues` and
  `scatter`/`scatterValues` do the same thing, one version on
  `LocalValues<N>`, the other on a raw pointer.
- [MODERATE, over-engineering / boilerplate] FiniteStrainIntegrators.hxx:
  three structures (FiniteStrainIntegrator, GreenLagrangeStrainIntegrator,
  LogarithmicStrainIntegrator), plus StandardIntegrator in Integrate.hxx,
  which only differ by the namespace called. 300 lines, 80 % of which is
  repeated documentation. The style "strategy = struct of static template
  functions" is reasonable, but these structures would belong in
  `mfront::gb::*` (DESIGN.md actually mentions
  `mfront::gb::StandardIntegrator`, which does not exist: the type is in
  `mfront::gpb`).
- [IMPORTANT, MGIS/mfem-mgis integration] GPU.hxx:178-245: each call does
  cudaGetDevice, cudaDeviceGetAttribute, cudaMallocAsync,
  cudaFuncGetAttributes, cudaDeviceGetLimit (under a mutex), launch, copy,
  cudaFreeAsync then `cudaStreamSynchronize`. The call is therefore always
  blocking and the report (OutputStatus) is always brought back to the
  host. For mfem-mgis with assembly on GPU, this imposes one
  synchronisation per call and per material part, with no way to avoid it.
  Same thing for SYCL (`wait_and_throw` 3 times, SYCL.hxx:338-343) and
  Kokkos (reduction into a scalar). No asynchronous variant is planned in
  the C contract.
- [IMPORTANT, integration] GPU.hxx:49-83 and DeviceStackSize.hxx: the
  generated library modifies a global limit of the CUDA context
  (`cudaDeviceSetLimit(cudaLimitStackSize)`) without the caller knowing.
  This changes the local memory reserved for all the kernels of the
  process (MFEM, hypre). It is technically justified but must be
  documented as a side effect and the caller must be able to disable it
  (only the stlpar/openmp/sycl backend has an option, CUDA/HIP always
  applies it).
- [MODERATE, C contract] BatchData.h: no version number nor structure
  size. DESIGN.md:55-57 says that the structure "can evolve without
  breaking the prototype", but an MGIS compiled against an older version
  will silently pass a structure of the wrong size. A `version` field (or
  an exported symbol giving the version of the contract) costs little and
  protects MGIS.
- [MINOR, C contract] `statuses` is an `int*` while everything else is a
  field with strides (BatchData.h:152). `rdts` is a field, `statuses` is
  not. Status conventions inverted with respect to `generic` for the
  initialisation and post-processing functions (1 = success against 0),
  documented but a source of errors for MGIS.
- [MINOR] The functors passed to the kernels copy the whole `BatchData`
  (about 25 fields) plus `InputFields` (4 x N fields of 24 bytes plus
  sizes) plus the parameters. For a behaviour with many material
  properties, the size of the kernel arguments grows fast (limit 4 KB
  before CUDA 12.1). Nothing blocking, but to be watched.
- [GOOD] Integrate.hxx: clean separation between "view of a point"
  (`PointView`), "function applied at a point" (trivially copyable
  functors), "primitive of the backend" (`apply`). No extended lambda
  required. The direct/copy choice is made once on the host. The comments
  explain the choices and the compiler workarounds.
- [GOOD] Reduction per warp then atomics only if needed
  (GPU.hxx:121-156), ordered key for the min of rdt (Integrate.hxx:765-778).
  Simple and correct.

### Group 2: drivers of `generic` (GenericBehaviour/Integrate.hxx)

- [GOOD] The split "TFEL_HOST_DEVICE core + historical function that decodes
  K" is clean and avoids any duplication of the physics between `generic`
  and `generic-parallel`. `ComputationOptions` replaces scattered ad hoc
  decodings of `K[0]`. `MethodCaller<m>` keeps the API by member pointer on
  the host.
- [MODERATE, process] Commit 494929641 mixes the refactoring and six fixes
  of the behaviour of `generic` on CPU (status `if (r)` of the finite
  strain drivers, DTAU_DDF in Hencky, internal variables of the end of the
  time step in the post-processings, etc.). They are real fixes, well
  described in the message, but they contradict principle 4 of DESIGN.md
  ("No change of the CPU behaviour of existing interfaces") and should be
  separate commits (or even PRs) for upstream: they concern all the users
  of `generic`.
- [MINOR] `callAndReportExceptions(d, f, on_exception)`: the third
  argument is `[] {}` in 3 calls out of 4. An overload with two arguments
  would be more readable.

### Group 3: generators (interfaces and backends)

- [MODERATE, architectural consistency] Two philosophies coexist in the
  same `generic-parallel` family. The material properties (upstream) write
  the loops and kernels as strings in the generator
  (MaterialProperty/ParallelSTLBackend.cxx:68-120, OpenMPBackend.cxx and
  KokkosBackend.cxx added by the branch in this style). The behaviours put
  all the logic in template headers and the generator only writes a call
  (Behaviour/*Backend.cxx, 40 to 70 lines each). The second approach is
  clearly better, but the branch extended the first one (Kokkos, OpenMP for
  the material properties) instead of converging. Two hierarchies of
  backends, two factories
  (`generic_parallel::material_property::BackendFactory` and
  `generic_parallel::behaviour::BackendFactory`), same option key
  `backend` in the configuration files.
- [MODERATE, duplication] Parsing of the options duplicated between the
  two hierarchies: `execution_policy` (MaterialProperty/ParallelSTLBackend.cxx:
  29-46 and Behaviour/ParallelSTLBackend.cxx:27-49, same strings, slightly
  different messages), `device_stack_size` (5 copies), `data_transfer_policy`
  (SYCL). A small common function `parseExecutionPolicy` would be enough.
- [MINOR, duplication] GenericParallelBehaviourInterface.cxx and
  GenericParallelModelInterface.cxx repeat `setOptions` (identical except
  for the class name), `getLibraryName`, `getHeaderFileName`,
  `getSourceFileName`, `getTargetsDescription`, and the writing of the
  entry points (check of the null pointer, `handleExceptions`, parameters,
  `writeCall`).
- [MODERATE, fragile design] BehaviourCodeGeneratorBase.cxx,
  `writeInterfacesCode`: the compatibility `generic` + `generic-parallel`
  relies on the deduplication of code blocks by text equality. One space
  or one comment of difference is enough to produce duplicate definitions
  (compilation error hard to read for the user). An explicit relation
  ("this interface reuses the members of `generic`", for example a method
  `getBaseInterfaceForGeneratedMembers`) would be more robust. It is
  pragmatic and short, so acceptable for a prototype, but it is the kind
  of implicit mechanism that upstream will probably refuse.
- [IMPORTANT, design] `areDynamicallyAllocatedVectorsAllowed()` and
  `getParametersSourceInBehaviourConstructor()` are properties of an
  interface (AbstractBehaviourInterface.hxx:93-111) while the class of the
  behaviour is common to all the interfaces. Adding `generic-parallel`
  therefore changes the class generated for `generic`, `castem`, etc. in
  the same call of mfront (no dynamic vectors, additional constructor).
  It is intended, but it should be documented as such and named as a
  constraint on the class (for example "requiresDeviceCompatibleClass").
- [IMPORTANT, semantics] Parameters treated as static variables
  (BehaviourCodeGeneratorBase.cxx, `writeBehaviourParameters`): in device
  code, `static constexpr` becomes a non static `const` member under
  `#ifdef TFEL_DEVICE_CODE`. The class therefore has neither the same
  memory layout nor the same definition in the host and device passes of
  nvcc, clang and hipcc (ODR violation). It works because the behaviour
  object is never passed from the host to the device, but it is fragile:
  the day MGIS or mfem-mgis wants to build or copy a behaviour on the host
  side for the device, it will break silently. Simpler and safe: the same
  declaration in both passes (`const` member everywhere, at the cost of a
  few bytes), or a device constant of namespace scope on the model of
  `getConstantObject`. Lines: BehaviourCodeGeneratorBase.cxx:3809-3814.
- [IMPORTANT, semantics and tests] TFELConfig.hxx:383-390:
  `TFEL_DEVICE_CODE` is defined in the whole translation unit with
  `nvc++ -stdpar=gpu`, `nvc++ -mp=gpu` and AdaptiveCpp SSCP. With the
  provided configuration files, the `generic` library compiled in the same
  call (`--interface=generic --interface=generic-parallel`) is also
  compiled with these options: it loses the exceptions, the bounds
  warnings, the parameters read from the initialisers in some paths, etc.
  In the tests of these configurations
  (stlpar-parunseq-nvhpc-gpu/test.check, openmp-offload-nvhpc-gpu,
  sycl-default-acpp), the "reference" is therefore not independent: the
  device paths are compared to the device paths, with the same compiler.
  The name `TFEL_DEVICE_CODE` is also misleading in this case ("may be
  device code"). The reference should be compiled with the usual host
  compiler, or at least the tests should say so.
- [MINOR, contract bug] Generated entry points
  (GenericParallelBehaviourInterface.cxx, `write_function`, and
  GenericParallelModelInterface.cxx:318-323, GenericParallelBehaviourInterface.cxx:725-730): if `mfront_batch_data` is
  null, `report(-2)` is called without `initialize`: `number_of_failures`,
  `number_of_unreliable_results` and `rdt` of the status remain
  uninitialised. The rotation functions, on the other hand, call
  `initialize` before.
- [MINOR, disorder] GenericParallelBehaviourInterface.cxx:102-111: orphan
  documentation block (the one of `getPointDataSizes`) followed by another
  block; the documented function is at line 149. Inconsistent prefixes of
  the generated names: `mgb_d`, `mgpb_parameters`, `mfront_parameters`,
  `mfront_data`.
- [MINOR] GenericParallelModelInterface.cxx:344-346: the out of bounds
  policy passed to the models is `None` and ignored (ModelIntegrator.hxx:59,
  anonymous parameter). The comment says so, but the API passes a useless
  value.
- [MINOR, API] KokkosBackend: the string `execution_space` is copied as is
  as a template argument in the generated code
  (Behaviour/KokkosBackend.cxx:39-43), without validation.
- [GOOD] `GPUBackend` factors CUDA and HIP (a single class, two names, two
  extensions). The factory and the validation of the options follow
  exactly the model of the material properties. The generator of the
  interface is 840 readable lines, without physics.
- [GOOD] The models reuse `writeGenericModelClass` with an option
  `device` rather than a copy of the `generic` generator.

### Group 4: tests

- [MODERATE, maintainability] The list of the 29 .mfront files is copied
  identically in the 15 test.check (checked: same md5 hash for the 15
  configurations). Adding a test case requires modifying 16 files.
  compare-all.sh centralises the comparisons but not the generation.
- [MODERATE, maintainability] compare-all.sh passes by hand the sizes of
  each behaviour (`6 6 2 7 1 36 ...`) while the `generic` library exports
  this information (symbols read by MGIS). Any modification of a test
  behaviour requires recomputing these numbers. For the rotations, 8 dummy
  arguments (`0 0 0 0 0 0 0 - 0 --rotation=6`) show that the command line
  interface of compare.cxx has reached its limits.
- [MINOR] compare.cxx (1235 lines) mixes five execution environments by
  macros (Kokkos, SYCL, SYCL host memory, CUDA/HIP, CPU). It is reasonable
  for a test bench, but a future test through MGIS would advantageously
  replace a large part of this code.
- [GOOD] The reference is the `generic` interface point by point, on two
  time steps, three memory layouts, three types of computation, with
  partial failures, energies, speed of sound, initialize functions,
  post-processings and rotations. The functional coverage is very good.

### Group 5: TFEL headers (include/TFEL)

Measurement: 163 files, +2362/-1836 lines; 757 added lines out of 2362
contain a TFEL_HOST_DEVICE annotation. The rest: induced reformatting,
replacement of exceptions, and a handful of real changes.

- [GOOD] The compiler workarounds in the TFEL headers are few and all
  commented with version and symptom: ConstantObject.hxx:29,44 (nvcc,
  nvc++), BoundsCheck.hxx:49 (nvcc), TensorViewFromStensor.hxx:57
  (nvc++ 25.11, switch), GenericFixedSizeArray.ixx:36 (nvc++ 25.11, value
  initialisation). The commit messages (58f4ced56, a22a3cff6, f008859df)
  give the reproducers. No uncontrolled accumulation.
- [GOOD] BoundsCheck.hxx: change from `std::string` to `const char*` and
  from the specialisations per dimension to a loop: -365/+199 lines,
  simpler than before. Boolean return on device, exception on the host:
  clear and documented compromise.
- [GOOD] TrivialGenTypeBase (tagged union, 249 lines, unit test
  TrivialGenTypeBaseTest.cxx): justified (std::variant not supported by
  nvc++ on device), isolated, same interface as GenTypeBase.
- [MODERATE, invasive] API/ABI breaks in the core of TFEL, all justified
  but to be negotiated with upstream:
  * MechanicalBehaviour.hxx: removal of the virtual functions (commit
    1b621828d), replaced by a concept. Any third party code that writes
    `override` on these methods no longer compiles;
  * BoundsCheck: signature `const char*` instead of `std::string`;
  * src/Material/LogarithmicStrainHandler.cxx removed: two exported
    symbols disappear from libTFELMaterial;
  * GenericFixedSizeArray: copy constructors and assignment operators
    removed (to be trivially copyable).
- [MODERATE, change of CPU results] StensorComputeEigenVectors.hxx
  (commit 9cd913eb2): for a multiple eigenvalue in 3D, the eigenvectors
  are now computed by Jacobi (`fses::syevj3`) instead of the analytical
  algorithm, on the host too. It is a numerical change (degenerate case)
  for all the users of TFEL, motivated by the device. It could have been
  limited to the device code, or at least flagged as a change of results.
- [MINOR, churn] GenTypeBase.hxx: all the methods receive
  TFEL_HOST_DEVICE and `std::get` becomes `*std::get_if`, while
  GenTypeBase (std::variant) is no longer used in device code since the
  introduction of TrivialGenTypeBase in the same commit (329a740e9).
  Unnecessary modification of a central header.
- [MINOR] ConstantObject.hxx:57-70: the return type of
  `getConstantObject` depends on the compiler (`T` with nvc++, `const T&`
  elsewhere). The rule "bind to a constant reference" is documented, but
  it is an API trap (an `auto x = ...` copies everywhere, an
  `auto& x = ...` only compiles elsewhere than with nvc++).
- [MINOR, cosmetic] Part of the diff is clang-format reformatting
  (`requires` clauses, order `[[nodiscard]] TFEL_HOST_DEVICE`). The order
  of the attributes is a real fix (C++ grammar), but the reindentation of
  the `requires` (LinearHomogenizationSchemes.ixx, etc.) makes the
  upstream review heavier.

### Group 6: points that will hinder MGIS and mfem-mgis

- [IMPORTANT] Calls always synchronous and report always brought back to
  the host (see group 1). mfem-mgis on GPU will want to chain integration,
  assembly and solver on the same stream without synchronisation.
- [IMPORTANT, silent semantics] On device, several values silently differ
  from the host: parameters of the called material laws (`mfront::device`,
  default values, MFrontMaterialPropertyInterface.cxx,
  `writeDeviceMaterialPropertyBody`), parameters of the behaviour
  variables (default values), standard bounds of the material laws not
  checked, bounds warnings removed. Documented in
  generic-parallel-behaviour-interface.md:453-474, but no detection at
  runtime: a user of MGIS who modifies a parameter of a material law will
  get GPU results different from the CPU results without any message.
- [MODERATE] Two libraries (`libBehaviour.so` and
  `libGenericParallelBehaviour-<backend>.so`) export the same C names
  (integration function, `setParameter`, `setOutOfBoundsPolicy`) and each
  contains its own copy of the compiled object of the behaviour (generated
  Makefile: `X.o` linked in both). The state (parameters, out of bounds
  policy) is therefore duplicated: if MGIS falls back on `generic` for a
  part of the points, it will have to forward each `setParameter` to both
  libraries. With an RTLD_GLOBAL loading (or a direct link), the identical
  names can even interpose.
- [MODERATE] No version of the C contract (see group 1).
- [MINOR] The error messages per point are lost
  (`bd.error_message = nullptr`, Integrate.hxx:531, `PointView`): in case
  of failure on CPU (stlpar, host openmp), the message of the exception of
  the behaviour is never reported, only "integration failed at some
  integration points". Keeping the first message (at least on the host)
  would help debugging in MGIS a lot.
- [MINOR] GPU.hxx:54: the mutex of `ensureStackSize` is static in a
  function template, so one mutex per kernel. The read-compare-write
  sequence of the stack limit is not protected between different kernels
  nor between libraries: with concurrent calls from several host threads
  (different materials in mfem-mgis), a lower limit can overwrite a higher
  limit just before a launch.
- [GOOD] The C contract by fields (pointer + 2 strides), with material
  properties and external variables "per variable", matches exactly the
  `MaterialDataManager` of MGIS (uniform values or per point): MGIS will be
  able to pass its data without copy. The symbol
  `<behaviour>_generic_parallel_backend` and the reuse of the symbols of
  `generic` simplify the loading. The role of each library
  (DESIGN.md:234-265) is clear and sound: TFEL does not allocate, MGIS
  depends on no GPU tool.

### Miscellaneous checked (cosmetic)

- cmake/modules/tfel.cmake:554: `_tfel-generic-parallel_behaviours_only_configurations`
  is an empty list, used in a test: dead scaffolding.
- RungeKuttaDSLBase.cxx:753-760 and RungeKuttaCodeGeneratorBase.cxx:2254-2265
  still write their `#ifdef TFEL_DEVICE_CODE` and their `useQt` branches
  by hand instead of using `makeHostAndDeviceCode` and
  `getInvalidStandardTangentOperatorFlagCondition` introduced by the
  branch: incomplete harmonisation.
- DESIGN.md:94-97 quotes `mfront::gb::StandardIntegrator` and
  `finite_strain::Integrator`, which do not exist under these names
  (`mfront::gpb::StandardIntegrator`, `mfront::gpb::FiniteStrainIntegrator`).

### Complements (checked)

- [MODERATE, consistency] Different status conventions in the same
  `generic-parallel` family: material properties (`mfront_gmp_OutputStatus`,
  GenericMaterialProperty/OutputStatus.h:26-75): 0 success, 1 caution,
  -2 exception, -4 invalid result, -6 stack (added by the branch);
  behaviours (`mfront_gpb_OutputStatus`, BatchData.h:187-209): 1 success,
  0 unreliable, -1 failure, -2 invalid arguments, -3 runtime error or
  exception. The choice "as `generic`" is defensible for the behaviours,
  but -2 does not have the same meaning from one structure to the other
  and MGIS will handle both.
- [MODERATE, scope] Six backends and fifteen tested configurations. Most
  of the workarounds and of the particular paths come from nvc++ (21
  mentions in the added lines outside tests: switch, value init, reference
  to a local flag, member pointers, assertions, imposed `managed` memory,
  stack via dlopen) and from AdaptiveCpp (no reductions, no recursion, no
  virtual functions). For upstream and for MGIS, a first batch CUDA/HIP +
  host stlpar/OpenMP would be much easier to get accepted and to maintain;
  the others could follow.
- [MODERATE, compilation cost] Each entry point is instantiated twice
  (`direct` true or false, Integrate.hxx:1081-1092), so two complete
  kernels per function, initialisations and post-processings included. The
  gain of the direct mode is measured nowhere (JOURNAL.md gives no
  comparison). For heavy behaviours (single crystals), this doubles an
  already long compilation time. To be measured, or to be made optional.

## RANKED SUMMARY

Overall verdict: the core design is good (a single implementation per
point shared with `generic`, logic in headers, simple C contract suited to
MGIS). The serious defects are side effects on the existing code and
silent host/device divergences, not the structure of the new code. The
duplication is real but local and easy to remove.

### A. Serious problems

A1. `TFEL_DEVICE_CODE` defined in the whole translation unit for nvc++
-stdpar/-mp and AdaptiveCpp SSCP (TFELConfig.hxx:383-390). With the
provided configuration files, the `generic` library of the same call is
compiled with these options (generated Makefile: `X.o` linked in
libBehaviour.so and in libGenericParallelBehaviour-stlpar.so). It then
loses exceptions and warnings, and the default constructor of the
behaviour variables takes the default values of the parameters even on
the host (BehaviourCodeGeneratorBase.cxx:2446-2455). The tests of these
configurations (stlpar-parunseq-nvhpc-gpu, openmp-offload-nvhpc-gpu,
sycl-default-acpp) therefore compare device path against device path, with
the same compiler: the reference is not independent.

A2. Silent host/device divergences, documented but never detected:
parameters of the called material laws (`mfront::device`, default
values, MFrontMaterialPropertyInterface.cxx,
`writeDeviceMaterialPropertyBody`), parameters of the behaviour variables,
standard bounds not checked, warnings removed. A `setParameter` on the
MGIS side will give GPU results different from the CPU results without
any message. At a minimum, there should be a refusal (or a warning) at
runtime when a modified parameter is not passed to the device.

A3. Calls always blocking, global side effects. GPU.hxx:178-245:
cudaMallocAsync, copy, launch, copy, cudaFreeAsync then
cudaStreamSynchronize at each call; SYCL.hxx:338-343: three
`wait_and_throw`; Kokkos: reduction into a scalar. The report is always
brought back to the host. The stack limit of the CUDA context is modified
without any possibility of refusal for CUDA/HIP (GPU.hxx:49-83), and via
dlopen of libcudart for the others (DeviceStackSize.hxx). For mfem-mgis
with assembly on GPU, it is the first thing to revisit: an asynchronous
variant (report left in device memory, synchronisation left to the caller)
and an explicit stack setting.

A4. The behaviour class does not have the same definition in the host and
device passes: parameters treated as static variables declared `static
constexpr` on the host and non static `const` members on the device
(BehaviourCodeGeneratorBase.cxx:3809-3814). ODR violation for nvcc, clang
and hipcc. It works as long as the behaviour object is never built on the
host then copied to the device; it will be wrong the day MGIS or
mfem-mgis do so.

A5. Changes of behaviour of the existing code mixed with the prototype,
contrary to principle 4 of DESIGN.md ("No change of the CPU behaviour of
existing interfaces"): six fixes of the `generic` drivers in commit
494929641, fixes of the rotations (79eacfee1), eigenvectors of a multiple
eigenvalue computed by Jacobi on the host too (9cd913eb2), API breaks
(MechanicalBehaviour without virtual functions, 1b621828d; BoundsCheck as
`const char*`), ABI break (LogarithmicStrainHandler.cxx removed from
libTFELMaterial). Each change is justified and well explained in its
commit message, but the whole must be split into distinct PRs for
upstream and the principle of DESIGN.md corrected.

### B. Moderate problems

B1. Two architectures in the same `generic-parallel` family: material
properties generated as strings (the branch adds Kokkos and OpenMP to them
in this style), behaviours by template headers. Two factories, options
parsed twice (`execution_policy`:
MaterialProperty/ParallelSTLBackend.cxx:29-46 and
Behaviour/ParallelSTLBackend.cxx:22-49; `device_stack_size`: six copies),
two status conventions (`mfront_gmp_OutputStatus`: 0 success,
-2 exception; `mfront_gpb_OutputStatus`: 1 success, -2 invalid
arguments).

B2. Duplication in the execution headers: four entry points copied per
backend while they only bind `makeApply` (GPU.hxx:280-364,
ParallelSTL.hxx:103-203, OpenMP.hxx, Kokkos.hxx); four SYCL entry points
of 40 to 50 lines repeating the same sequence (SYCL.hxx:369-581); three
identical `dispatch*` (Integrate.hxx:1081-1092, 1241-1258, 1286-1303);
four integration structures that only differ by the namespace called
(FiniteStrainIntegrators.hxx, 300 lines). The generated code could call
`mfront::gpb::integrateBatch` with the `makeApply` of the backend: 300 to
400 fewer lines.

B3. Implicit mechanisms in the generator: compatibility
`generic` + `generic-parallel` obtained by deduplicating code blocks by
text equality (`writeInterfacesCode`, BehaviourCodeGeneratorBase.cxx);
interface properties (`areDynamicallyAllocatedVectorsAllowed`,
AbstractBehaviourInterface.hxx:93-111) that modify the class shared by
all the interfaces of the same call. Pragmatic, but fragile and hard to
read.

B4. C contract and MGIS integration: no version of the contract (while
DESIGN.md:53-57 counts on the evolution of the structure); state
duplicated between libBehaviour.so and libGenericParallelBehaviour-*.so
(parameters, out of bounds policy, same C names), to be forwarded by MGIS
in case of fallback on `generic`; error message per point lost
(Integrate.hxx:531).

B5. Tests hard to maintain: the list of the 29 .mfront files copied
identically in 15 test.check; sizes of each behaviour passed by hand in
compare-all.sh while the `generic` library exports them; dummy positional
arguments for the rotations (`0 0 0 0 0 0 0 - 0 --rotation=6`).

B6. Scope: six backends and fifteen configurations. Most of the
workarounds come from nvc++ (21 mentions in the added lines) and from
AdaptiveCpp. A first batch CUDA/HIP + host stlpar/OpenMP would be easier
to get accepted.

B7. Double direct/copy instantiation of each entry point, without
measurement of the gain: compilation time doubled for the heavy
behaviours.

### C. Minor and cosmetic

- Integrate.hxx:110-213: `gather`/`gatherValues`, `scatter`/`scatterValues`
  duplicated.
- Null batch pointer: `report(-2)` without `initialize`, fields of the
  status uninitialised (GenericParallelBehaviourInterface.cxx:725-730,
  GenericParallelModelInterface.cxx:318-323).
- Orphan documentation block (GenericParallelBehaviourInterface.cxx:102-111);
  inconsistent generated prefixes (`mgb_d`, `mgpb_parameters`, `mfront_parameters`).
- GenTypeBase.hxx annotated TFEL_HOST_DEVICE and modified (`*std::get_if`)
  while it is no longer used on device.
- `getConstantObject`: return type that depends on the compiler
  (ConstantObject.hxx:57-69).
- GPU.hxx:54: static mutex per template instantiation, the
  read-compare-write sequence of the stack is not protected between
  kernels.
- Models: out of bounds policy `None` passed then ignored
  (GenericParallelModelInterface.cxx:344-346, ModelIntegrator.hxx:59).
- KokkosBackend: `execution_space` copied without validation in the
  generated code.
- tfel.cmake:554: empty list `_tfel-generic-parallel_behaviours_only_configurations`.
- RungeKutta: `#ifdef TFEL_DEVICE_CODE` and `useQt` branches written by
  hand instead of the new utilities (RungeKuttaDSLBase.cxx:753-760,
  RungeKuttaCodeGeneratorBase.cxx:2254-2265).
- DESIGN.md:94-97: names of types that do not exist
  (`mfront::gb::StandardIntegrator`).
- Size of the kernel arguments (BatchData + InputFields + parameters
  copied by value) to be watched for the behaviours with many properties.
- clang-format reformatting of `requires` clauses that makes the diff
  heavier.

### D. What is done well

- A single implementation per point: the `generic` drivers are split into
  a TFEL_HOST_DEVICE core (explicit options, arguments of the constructor
  passed) and the historical functions. No duplicated physics.
- Clean separation view of a point / trivially copyable functor /
  primitive `apply` of the backend; no extended lambdas; direct/copy
  choice made once on the host.
- C contract by fields (pointer + two strides) with material properties
  and external variables "per variable": matches exactly
  `MaterialDataManager`, MGIS will be able to pass its data without copy.
  Roles of TFEL, MGIS and mfem-mgis clear (DESIGN.md:234-265).
- Reduction per warp and atomics only if needed; ordered key for the
  minimum of rdt.
- Few compiler workarounds in the TFEL headers (four), all commented with
  version and symptom; reproducers in results/nvhpc-bugs.
- BoundsCheck simplified (-365/+199 lines); TrivialGenTypeBase justified
  and tested; `change_basis` added with a test; rotations shared between
  `generic` and `generic-parallel` (BehaviourRotationFunctions) with fixes
  of `generic`.
- Models made by reusing `writeGenericModelClass` (option `device`); CUDA
  and HIP backends factored in `GPUBackend`.
- Exemplary commit messages (motivation, compiler defect, fixes);
  complete user documentation; very good functional coverage of the tests
  (two steps, three layouts, three types of computation, partial
  failures, energies, speed of sound, initialize functions,
  post-processings, rotations, models).

### E. Suggested work order before MGIS

1. Asynchronous variant and explicit stack setting (A3).
2. Detection of the parameters not passed to the device (A2).
3. Same class definition on host and device (A4).
4. `generic` reference compiled by the host compiler in the tests (A1).
5. Factor the entry points and the `dispatch` (B2), version of the C
   contract (B4).
6. Split the branch into upstream PRs: fixes of `generic`, TFEL device,
   interface (A5).
