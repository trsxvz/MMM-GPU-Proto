# D1, MGIS: integration of a batch with a `generic-parallel` library

> Historical document. Design proposal written before the work on MGIS,
> then completed during the first steps. The current state is described
> in [mgis.md](../mgis.md), the path that led to it in the
> [log](../../JOURNAL.md).

Status at the time of writing: steps 1 to 3 done on 2026-10-08 (branch
`gpu-proto` of `mgis/`, `patches/mgis`). Plan revised on 2026-10-08.
Bases: MGIS master `818f0f9`, TFEL of the proto (`tfel/`, R1 to R8:
structures in version 2, index of the first failed point).

Differences with the plan, found while coding:

- the uniform values kept by MGIS (`real` of the `FieldHolder`, rotation
  matrices built from axes) are in host memory: refused with a message for
  `cuda` and `hip`. The caller gives them by a `std::span` in accessible
  memory (material properties now accept a `std::span` of one value, like
  external state variables);
- the first failed point (`n`) is known without an array of statuses
  (index provided by TFEL);
- the rotations do not yet have an overload taking an execution policy
  (default queue of the backend);
- the libraries are no longer unloaded at the end of the process (TBB
  threads).
The version of 2026-10-06 (detailed analysis of MGIS, line numbers) stays
in the git history of this file.

## 1. What TFEL now provides

- Versioned C contract: `MFront/GenericParallelBehaviour/BatchData.h`
  (`MFRONT_GPB_API_VERSION`) for behaviours and models,
  `MFront/GenericParallelMaterialProperty/MaterialProperty.h`
  (`MFRONT_GPMP_API_VERSION`) for material properties.
- Symbols of each library: `<b>_mfront_interface` =
  `GenericParallel`, `<b>_generic_parallel_backend`,
  `<b>_generic_parallel_api_version` (behaviours, models and material
  properties).
- Material properties and external state variables given group by group
  or variable by variable: each `FieldHolder` of MGIS becomes a field,
  without copy (uniform: zero stride; array: `{p + b * s, s, 1}`).
- Empty batch accepted even with null pointers; invalid computation
  options refused for the batch (`-2`); output status always complete.
- Execution: queue or stream (`queue`), asynchronous statistics
  (`mfront_gpb_BatchStatistics`, C decoding functions), in-order SYCL queue
  required in asynchronous mode, stack of the GPU threads set by the
  library.
- Shared data (`@SharedData`) copied by the library, once per GPU: nothing
  to do in MGIS, several GPUs per process possible.
- Parameters: set in the `generic-parallel` library itself
  (`<b>_setParameter`), read on the host at each call.

## 2. What MGIS already provides and simplifies everything

- **Storage provided by the caller.** `MaterialDataManagerInitializer` and
  `MaterialStateManagerInitializer` accept external arrays for the
  gradients, the forces, the internal state variables, the energies, `K`
  and the speed of sound; the `FieldHolder` of material properties and
  external state variables accept `std::span`. A code that allocates GPU
  memory itself (mfem-mgis with the memory of MFEM) can therefore already
  place everything: **no new memory API is needed for mfem-mgis**.
- Private internal storage (`std::vector`): it can change allocator
  without breaking the API.
- `LibrariesManager::getInterface` exists (never called).
- `MGIS_HOST_DEVICE` and a `mgis::function` module with CUDA algorithms
  already exist: D1 does not depend on them, but does not prevent
  connecting to them later.

What mfem-mgis (D2) expects: today the integration is done point by point
in the middle of the assembly loop
(`BehaviourIntegratorBase::performsLocalBehaviourIntegration`). On GPU,
the assembly becomes three phases: computation of the gradients at all
points (on the GPU), integration of the batch, assembly. MGIS must
therefore:

1. integrate a batch (or a range `[b, e)`) on arrays that stay on the GPU;
2. never read or write these arrays on the host during the integration;
3. allow several materials to be submitted before a single
   synchronisation;
4. update the state (`update`, `revert`) without bringing the data back.

## 3. Principles

1. **MGIS stays a host library without GPU dependency** (compiled by
   g++). It loads, checks, describes the data and calls; the computation
   is in the generated library.
2. **Zero copy** between `MaterialDataManager` and `BatchData`, material
   properties and external state variables included.
3. **MGIS never allocates in the memory of a GPU.** The memory comes from
   the caller (external arrays), or from a memory resource that it gives
   (section 5, convenience for simple codes and tests).
4. **The batch path does not touch the content of the arrays on the host**
   (no copy, no reading of the statuses per point): the memory can be
   purely GPU memory.
5. **Explicit errors instead of undefined behaviours.** MGIS knows the
   backend (symbol) and the location of the arrays, declared at
   construction (section 4.2): a GPU backend with ordinary host memory is
   refused with a clear message. The per point entry points refuse a
   `generic-parallel` library, and conversely.
6. **Automatic switch.** The existing functions (`integrate`,
   `executeInitializeFunction`, `executePostProcessing`, rotations)
   recognise a `generic-parallel` library: an existing code works without
   change as soon as its memory suits the backend.
7. **Statuses.** `1`, `0`, `-1`: usual integration status; `-2`
   (arguments) and `-3` (backend error): exceptions, or errors of the
   `Context` in the `noexcept` variants. A bug is never taken for a
   convergence failure (which would make the time step decrease endlessly).

## 4. API

### 4.1 Loading

- `BehaviourDescription`: `interface` (`""`/`Generic` or
  `GenericParallel`), `backend`; API version checked at loading time.
- `Behaviour`: batch pointers (`batch_b`, batch initialize functions and
  post-processings, three rotations). For a `generic-parallel` library,
  the per point pointers stay null.
- `LibrariesManager`: `getBatchBehaviour`, `getBatchInitializeFunction`,
  `getBatchPostProcessing`, `getBatchRotation`; the per point functions
  check the interface.
- Models (`libGenericParallel<...>-<backend>.so`): same path, without
  additional work.
- MGIS keeps its own copy of the C structures (as for the `generic`
  interface, without compile-time dependency on TFEL); a test checks the
  memory layout against `BatchData.h`.

### 4.2 Location of the data and execution policy

The location of the data is a property of the data, declared at
construction; the execution policy only describes the execution.

```cpp
namespace mgis::behaviour {
  //! \brief location of the arrays of a material data manager
  enum struct MemoryLocation {
    HOST,    //!< ordinary memory of the host (default)
    DEVICE,  //!< memory of the device, not accessible by the host
    SHARED   //!< accessible by the host and the device (managed, unified)
  };
  // MaterialDataManagerInitializer: MemoryLocation location = HOST;
  //! \brief execution of a batch (generic-parallel libraries)
  struct BatchExecutionPolicy {
    //! \brief queue, stream or execution space (backend dependent)
    void* queue = nullptr;
    //! \brief see mfront_gpb_ExecutionOptions
    bool preserve_stack_size = false;
  };
}
```

MGIS crosses the location and the backend: `HOST` memory refused by a GPU
backend (except SYCL with automatic transfers); host functions
(extraction, saving, local values) refused for `DEVICE`.

`MaterialDataManager` carries a default policy
(`m.batch_execution_policy`), used by the automatic switch; the overloads
that take a policy replace it for one call.

### 4.3 Functions

```cpp
// synchronous integration of all the points, or of [b, e)
BehaviourIntegrationResult integrate(const BatchExecutionPolicy&,
                                     MaterialDataManager&,
                                     const BehaviourIntegrationOptions&,
                                     real dt);
BehaviourIntegrationResult integrate(const BatchExecutionPolicy&,
                                     MaterialDataManager&,
                                     const BehaviourIntegrationOptions&,
                                     real dt, size_type b, size_type e);
// Context variants (noexcept), initialize functions, post-processings
// and rotations with the same pattern
```

- Result: existing `BehaviourIntegrationResult`, extended with the number
  of failures and of unreliable results. The first failed point (`n`) is
  known only if the caller provides an array of statuses.
- Statuses and time step factors per point: optional arrays of the
  `MaterialDataManager`, provided like `K`
  (`useExternalArrayOfIntegrationStatuses`, `...OfTimeStepScalingFactors`).
- Initial time step factor: the same as the per point path (identical
  results).
- `ThreadPool&`: ignored for a `generic-parallel` library (the parallelism
  is the one of the backend).
- `integrate_debug`: refused for a `generic-parallel` library (use the
  `generic` library of the same behaviour).

### 4.4 Asynchronous (for several materials)

```cpp
// submits the integration and returns at once (CUDA, HIP, Kokkos, SYCL
// with an in-order queue); the statistics must be accessible by the device
void submit(const BatchExecutionPolicy&, MaterialDataManager&,
            const BehaviourIntegrationOptions&, real dt,
            mfront_gpb_BatchStatistics&);
// once the queue is synchronized by the caller
BehaviourIntegrationResult getIntegrationResult(
    const mfront_gpb_BatchStatistics&);
```

### 4.5 Update of the state

`update` and `revert` copy `s1` into `s0` (or the reverse) on the host.
The owner of the memory provides the copy primitive, MGIS keeps the
algorithm (which blocks to copy):

```cpp
// copy: callable copying n reals (the copy of MFEM for mfem-mgis)
template <typename CopyFunction>
void update(MaterialDataManager&, const CopyFunction&);
template <typename CopyFunction>
void revert(MaterialDataManager&, const CopyFunction&);
```

Without a primitive, `update`/`revert` copy on the host and refuse
`DEVICE` memory.

## 5. Memory resource (convenience)

For the codes that let MGIS allocate: optional `std::pmr::memory_resource*`
in the initializers, used for all the internal storage (states, `K`, speed
of sound, statuses). Resources provided as header-only files, compiled by
the calling code (MGIS does not depend on them): CUDA and HIP managed
memory, SYCL shared USM, Kokkos shared space. The `FieldHolder` keep
`std::vector` (public API used by mfem-mgis): with a GPU backend, the
material properties and external state variables must be uniform or
`std::span` in accessible memory, otherwise an error is raised.

## 6. Tests

Test behaviours of MGIS compiled with `generic` (reference) and
`generic-parallel` for the configurations chosen when configuring MGIS
(`stlpar-parunseq-gcc` by default; `cuda-nvcc`, `kokkos`,
`sycl-default-acpp` as options). The GPU tests allocate managed memory
through the C API of the runtime (`cudaMallocManaged` called from g++ code,
without nvcc) and run **without HMM**, to catch any access to host memory.

1. Loading: interface, backend, version, pointers, cross refusals of the
   per point and batch APIs.
2. Memory layout of the copied structures.
3. Integration: batch against loop per point on identical states (rounding
   tolerance, loadings far from the convergence limits); three types of
   integration; behaviours of MGIS (Plasticity, Norton, Gurson, tensorial
   external state variables, speed of sound, energies, finite strain,
   Runge-Kutta, single crystal).
4. Fields: uniform, in-place and mixed material properties and external
   state variables; ranges `[b, e)`; empty batch.
5. Failures and errors: failure at a few points, missing SYCL queue, host
   memory refused on GPU, `-2`/`-3` never taken for `-1`.
6. Initialize functions, post-processings, rotations.
7. Switch: the existing integration tests replayed on `generic-parallel`.
8. Asynchronous mode and update on the GPU.
9. Memory resource: complete suite with the managed resource.

## 7. Steps

Each step: code, tests, complete MGIS test suite, commit.

1. Loading, detection, guards; C structures and layout test.
2. Synchronous integration (policy, mapping of the data, memory checks,
   statuses, ranges) and automatic switch; CPU tests (`stlpar`, `openmp`)
   then GPU tests (`cuda`, managed memory).
3. Initialize functions, post-processings, rotations.
4. Memory resource and header-only resources.
5. Asynchronous mode, `update`/`revert` on the GPU.
6. `generic-parallel` material properties in `mgis::material_property`
   (evaluation of an array of points), if mfem-mgis needs them.
7. Documentation. C, Python and Fortran bindings: later.

Steps 1 to 3 give an MGIS usable by mfem-mgis on CPU and on GPU (memory
provided by the caller); 4 and 5 serve performance and convenience; 6
depends on D2.

## 8. Risks

- Identical entry point names, different signatures: guards and test 1.
- Drift of the copied structures: version checked at loading time, test 2.
- Managed memory touched by the host between two time steps: page
  migrations (performance); hence the copy on the GPU of section 4.5 and
  the rule "the batch path does not touch the arrays".
- Different semantics from the per point path: all points are integrated
  (the per point path stops at the first failure), messages per point not
  transmitted.
- GPU stack of polycrystals (nearly 15 GB reserved): a single process per
  GPU for these behaviours.

## 9. Decisions to take

1. Automatic switch of the existing functions (recommended: yes).
2. Default policy carried by `MaterialDataManager` (recommended: yes).
3. Asynchronous mode and copy on the GPU from D1 on (step 5) or after the
   first trials of mfem-mgis (recommended: after step 3, depending on the
   measurements).
4. `generic-parallel` material properties in MGIS (step 6): depending on
   D2.
