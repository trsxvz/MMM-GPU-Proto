# D2, MGIS: API on views (proposal, version 4)

> Historical document. Design proposal written before the work on the
> views of MGIS, after three reviews (appendix). The current state is
> described in [mgis.md](../mgis.md#views-and-layouts), the path that led
> to it in the [log](../../JOURNAL.md).

Differences with the proposal, found while coding:

- the CUDA, HIP and Kokkos memory resources could not serve as upstream
  resource of a pool (alignment refused beyond the one of their memory):
  fixed by over-allocation, without which the remedy of 3.5 did not work;
- the extraction of an internal state variable on a state view goes
  through operations (`extractInternalStateVariable(output, state, name, ops)`)
  and therefore also works on GPU. Per variable accessors on a state view
  (`getGradientView`, etc.) are added, as on the manager;
- the overloads of the rotations and of `evaluate` on arrays build views
  and call the functions on views: a single code path;
- a copy of a view into itself does nothing (the copies of the standard
  library and of the GPU runtimes forbid overlaps).

Status: proposal, no code written. First step of the
[plan for mfem-mgis](../plan-mfem-mgis.md). This version integrates three
independent reviews (version 1: findings A1 to C5; version 2: findings N1
to N10; version 3: findings F1 to F6), all checked in the code. Their
handling is in the appendix. The last review concludes that the design
holds and that no remaining point requires changing the API.

## 1. Need

mfem-mgis wants to lay out the data at the Gauss points in two ways, as
chosen, on CPU as on GPU:

- **point by point**: the components of a point are contiguous
  (`value(i, j) = v[i * nc + j]`);
- **component by component**: the values of a component are contiguous
  over all the points of the material (`value(i, j) = v[j * n + i]`),
  coalesced accesses on GPU. This is the *strided coalesced* layout of
  `mgis::function` and of TFEL (`tfel::math::StridedCoalescedView`).

MGIS must compute on these data without copy and without overhead, with
all its functions, and finish a time step on GPU.

**TFEL has nothing to change** (checked by the two reviews). Its contract
describes each group by a pointer and two strides. Its kernels thus read
and write all the groups, the inputs of the initialize functions, the
outputs of the post-processings and the rotations. The SYCL backend copies
the strided fields. The arrays of per variable fields are read on the host
during the call, in asynchronous mode too, for all backends (functors
passed by value; OpenMP and `stlpar` synchronous).

## 2. What MFEM 4.9 and mfem-mgis show

- The ordering of the degrees of freedom (`FiniteElementSpaceOrdering`,
  `byNODES` or `byVDIM`) is independent of the layout at the Gauss points.
  `ElementRestriction` always returns element vectors `(ND, VDIM, NE)`.
- The `byNODES` layout of `QuadratureInterpolator` is by element blocks
  (`q + NQ * (c + K * e)`). It cannot be used as is: the gradients of the
  behaviour require a conversion kernel, and the points of a material are
  not contiguous in it. This kernel writes in the chosen layout.
- mfem-mgis: `Material` derives from `MaterialDataManager`. The current
  path assumes the point by point layout in about forty places and keeps
  this layout. The component by component layout is used by the new batch
  paths. mfem-mgis gives its material properties and external state
  variables with the `NOUPDATE` policy.

## 3. Proposal

### 3.1 `FieldView`: a strided field, with its extents

The two-dimensional counterpart of `std::span`, with the meaning of
`std::mdspan` and `layout_stride` (absent from C++20 and from gcc 14):

```cpp
namespace mgis {
  //! layout of the arrays of a manager; COMPONENT_BY_COMPONENT is
  //! the "strided coalesced" layout of mgis::function
  enum struct FieldLayout { POINT_BY_POINT, COMPONENT_BY_COMPONENT };

  //! value(i, j) = data()[i * getPointStride() + j * getComponentStride()]
  template <typename ValueType>  // real or const real
  struct FieldView {
    //! empty view: never a uniform field by default
    constexpr FieldView() noexcept = default;
    //! preconditions: sufficient extent, consistent strides (see below)
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             std::span<ValueType>, size_type n,
                                             size_type nc, size_type ps,
                                             size_type cs);
    constexpr FieldView(std::span<ValueType>, size_type n, size_type nc,
                        size_type ps, size_type cs);
    template <bool doPreconditionsCheck>
    constexpr FieldView(const PreconditionsCheck<doPreconditionsCheck>&,
                        std::span<ValueType>, size_type n, size_type nc,
                        size_type ps, size_type cs);
    //! free conversion of a mutable view into a constant view
    template <typename OtherValueType>
    requires(std::is_same_v<const OtherValueType, ValueType>)
    constexpr FieldView(const FieldView<OtherValueType>&) noexcept;
    // data(), getNumberOfPoints(), getNumberOfComponents(),
    // getPointStride(), getComponentStride(), empty(), operator()(i, j)
   private:
    std::span<ValueType> values;
    size_type n = 0, nc = 0, point_stride = 0, component_stride = 0;
  };

  FieldView<T> makeFieldView(std::span<T>, FieldLayout, size_type n,
                             size_type nc);                        // dense
  FieldView<T> makePointByPointFieldView(std::span<T>, size_type n,
                                         size_type nc, size_type ps); // ps >= nc
  FieldView<T> makeComponentByComponentFieldView(std::span<T>, size_type n,
                                                 size_type nc, size_type ld); // ld >= n
  FieldView<T> makeUniformFieldView(std::span<T>, size_type n); // zero point stride
  FieldView<T> getComponents(FieldView<T>, size_type offset, size_type count);
}
```

- Error handling as in `mgis::function`: static `checkPreconditions`, and
  a constructor taking `PreconditionsCheck` that calls
  `check_preconditions<FieldView>` (check that can be bypassed when the
  caller has already done it). Without deriving from
  `PreconditionsChecker`, which has neither a default constructor nor
  assignment: `FieldView` stays a regular type, like `std::span`.
- Preconditions: an empty set (`n == 0` or `nc == 0`) is accepted first,
  without computing the extent; otherwise the extent covers
  `(n - 1) * ps + (nc - 1) * cs + 1` values.
- **Outputs** (checked at the call): no value may serve two pairs
  `(i, j)`. Exact, constant time check of a sufficient condition: when the
  two dimensions with an extent greater than one are sorted by increasing
  stride, the small stride is non-zero and the large stride is at least
  the small stride times its extent. Exotic interleavings are therefore
  refused for outputs. Inputs may be uniform.
- **Empty groups** (models without gradients, behaviours without internal
  state variables, energies not computed): the view must be empty, without
  any other check.
- Interoperability with `mgis::function`, in one direction only: a
  `FunctionView` and a *strided coalesced* view (to which a public access
  to its data is added) convert into a `FieldView`. The
  `CoalescedMemoryAccess` views (one pointer per component) have no
  equivalent: a single pointer per field, like the contract of TFEL.

### 3.2 `MaterialDataView`: the description of a material on arrays owned by others

`MaterialDataManager` owns its arrays. `MaterialDataView` describes arrays
of the caller, with any strides. These are regular types (copyable and
assignable): the behaviour is kept by pointer, as `BehaviourDataView` keeps
its arrays.

```cpp
namespace mgis::behaviour {
  struct MaterialStateView {
    MaterialStateView(const Behaviour&, size_type n,
                      std::pmr::memory_resource* = nullptr);
    FieldView<real> gradients, thermodynamic_forces,
        internal_state_variables, stored_energies, dissipated_energies;
    // one view per variable, in the order of the behaviour, empty at first
    std::pmr::vector<FieldView<const real>> material_properties;
    FieldView<const real> mass_density;
    std::pmr::vector<FieldView<const real>> external_state_variables;
    const Behaviour* b;  // never null
    size_type n;
  };
  // by name, with Context variants, as for MaterialStateManager
  void setMaterialProperty(MaterialStateView&, std::string_view,
                           FieldView<const real>);
  void setExternalStateVariable(MaterialStateView&, std::string_view,
                                FieldView<const real>);
  void setMassDensity(MaterialStateView&, FieldView<const real>);

  struct MaterialDataView {
    MaterialDataView(const Behaviour&, size_type n, MemoryLocation,
                     std::pmr::memory_resource* = nullptr);
    MaterialStateView s0, s1;
    FieldView<real> K;               // if an operator is requested
    FieldView<real> speed_of_sound;  // if requested
    //! policy of the functions that take none, like the manager
    BatchExecutionPolicy batch_execution_policy;
    const Behaviour* b;  // never null
    size_type n;
    MemoryLocation memory_location;
  };
  //! view on the arrays of a manager, in its layout, with its
  //! policy; the per variable vectors are allocated by the resource
  MaterialDataView make_view(MaterialDataManager&,
                             std::pmr::memory_resource* = nullptr);
}
```

- `std::pmr` vectors: copying a view allocates with the default resource,
  moving keeps the resource. A view built on a stack buffer must not leave
  its scope, and two views with different resources are not swapped
  (`swap`). Only the functions of the manager use a stack buffer, for a
  view that does not leave the call. `make_view` without a resource uses
  the default resource.
- `make_view` checks that the values kept by the manager on the host
  (scalars, local storage) are readable by the computations of the
  behaviour, as today: a view obtained by `make_view` bypasses no check.
- Public members, like the manager. Each function checks the consistency
  of the view (same behaviour and same number of points in the view and
  its states, number of per variable views, extents of each field,
  location). One check per variable, on the host, per call.
- Three levels: `MaterialDataManager` (owns), `MaterialDataView`
  (describes, lasts for the whole computation), and at each call the C
  structure `BatchDataView` (contract of TFEL). Copying a view copies
  descriptors, never values.

**Functions on a view.** The same as on the manager. The "without policy"
form of the manager, which switches between batch and point by point
computation, takes on a view an **explicit workspace**:

| Function | Forms on a view |
|---|---|
| `integrate` | `(policy, v, opts, dt[, b, e])` and `(asynchronous policy, ...)`: batch library; `(workspace, v, opts, dt[, b, e])`: any library (batch with the policy of the view, or point by point on the host); `integrate_debug(workspace, v, ...)` |
| `executeInitializeFunction` | same forms, inputs as `FieldView` or absent |
| `executePostProcessing` | same forms, outputs as `FieldView` |
| rotations (gradients, forces, blocks of K) | matrices as `FieldView` or from the axes (`RotationMatrix2D/3D`), in place or not; synchronous, asynchronous, without policy |
| `evaluate` (material properties) | output and inputs as `FieldView`, with or without `OutOfBoundsPolicy`; synchronous, asynchronous, without policy |

- The number of points comes from the extents: no separate `n` any more.
- **A single code path**: the functions of the manager call those of the
  view. Without a policy: `make_view(m)` and the workspace of the calling
  thread. With `ThreadPool`: one workspace per thread, on ranges. What is
  specific to the manager stays there (allocation of K on demand, check of
  the values kept on the host).
- **Workspace** (`BehaviourIntegrationWorkSpace`): one per thread,
  provided by the caller. It keeps a pointer to the behaviour for which it
  was built (checked at each call). It receives, allocated once at its
  construction, the buffers of one point for `s0` and `s1`: gradients,
  forces, internal state variables, energies, K, largest input of an
  initialize function, largest output of a post-processing.
- The functions of the manager without a policy route a `generic-parallel`
  behaviour to the batch **before** taking a workspace, as today: neither
  lock nor allocation on this path.
- **`generic` library**: point by point computation on the host. Two
  versions of the loop, chosen once per call: direct pointers if all the
  components are contiguous (the current case, unchanged), otherwise
  values gathered in the workspace then copied back. The gathering
  reproduces the current computation (integration options written in
  `K[0]`, properties and external state variables read with their
  strides). The rotations of a `generic` library on non-contiguous fields
  do the same with a buffer of one point (on the stack, on the heap
  beyond).
- No `ThreadPool` on a view: CPU parallelism goes through the `openmp` and
  `stlpar` backends, or through the caller with its workspaces.
- The old forms with `IntegrationType`, kept for compatibility (and used
  by the C bindings), stay on the manager only: views are a new API.
- Host functions on a view, as on the manager:
  `extractInternalStateVariable(output, state, name)`,
  `convertFiniteStrainStress` and `convertFiniteStrainTangentOperator`
  (output as `FieldView`), readable by the host only.

### 3.3 Layout of a `MaterialDataManager`

```cpp
auto i = MaterialDataManagerInitializer{};
i.layout = FieldLayout::COMPONENT_BY_COMPONENT;  // POINT_BY_POINT by default
```

The layout is stored in `MaterialDataManager::layout` and
`MaterialStateManager::layout` (new parameter of its constructor), like
the memory location. **Everything that is given per point to the manager
follows this layout**:

- state arrays, K, speed of sound;
- external state variables with several components given by a
  `std::span` (material properties are scalars);
- inputs of the initialize functions and outputs of the post-processings
  given by a `std::span`;
- outputs of `extractInternalStateVariable` and of the finite strain
  conversions, which also receive overloads taking `FieldView`.

What changes in MGIS:

- **F set to the identity** (`F-Cauchy` behaviours) with the strides of
  the layout;
- `updateValues(o, i)` checks that the two states have the same layout;
- `save` and `restore` always write and read point by point: a restart
  does not depend on the layout;
- access in the layout of the manager: `getGradientsView(s)`,
  `getGradientView(s, name)`, and likewise for the forces, the internal
  state variables and the blocks of K;
- the rotation overloads taking `std::span` stay point by point (they do
  not see the manager): the documentation says so and refers to the
  overloads taking `FieldView`;
- `MaterialFunctionManager`: constructor with an initializer; laid out
  component by component, it returns *strided coalesced* views through new
  fixed-size accessors (`getStridedCoalescedGradient<N>`, etc.). Its
  current accessors stay point by point (clear error otherwise):
  `mgis::function` has no *strided coalesced* view of dynamic size;
- bindings: the layout is exposed in the initializers and by an accessor
  in C, Fortran and Python. Python returns numpy arrays with the strides of
  the layout (pybind11). FEniCS always builds its manager point by point:
  nothing to do. Julia exposes no array. The `*_stride` accessors keep
  their meaning (number of values per point).

`update` and `revert` copy whole arrays: they do not depend on the layout.

### 3.4 End of time step on GPU: the same procedure, copies done by the caller

`update` and `revert` do more than copying (`updateValues` creates or
erases entries, changes a scalar into an array, resets the mass density to
zero, respects `shall_be_updated`). To keep exactly this meaning, MGIS
keeps its procedure and only delegates the copies of values to operations
provided by the caller, on the model of `std::pmr::memory_resource`:

```cpp
namespace mgis {
  struct FieldOperations {
    //! copy (same extents), possibly asynchronous on the queue of the
    //! caller; both operands are in the location of the
    //! manager or of the view
    void copy(FieldView<real>, FieldView<const real>);
    void fill(FieldView<real>, const real);
    virtual ~FieldOperations();
   private:
    // non-virtual public interface, like std::pmr::memory_resource
    virtual void do_copy(FieldView<real>, FieldView<const real>) = 0;
    virtual void do_fill(FieldView<real>, const real) = 0;
  };
  struct HostFieldOperations final : FieldOperations { /* std::copy */ };
}
namespace mgis::behaviour {
  void update(MaterialDataManager&, FieldOperations&);
  void revert(MaterialDataManager&, FieldOperations&);
  void updateValues(MaterialStateManager&, const MaterialStateManager&,
                    FieldOperations&);  // with Context variants
  void update(const MaterialDataView&, FieldOperations&);  // states and K
  void revert(const MaterialDataView&, FieldOperations&);
  void updateValues(const MaterialStateView&, const MaterialStateView&,
                    FieldOperations&);  // state arrays
}
```

- **A single procedure**: `update(m)` keeps its host access check, then
  calls `update(m, HostFieldOperations)`. Likewise for `revert(m)` and
  `updateValues(o, i)`. The current meaning is preserved case by case
  (resizing, mass density, insertions, erasures, `shall_be_updated`;
  checked by the review).
- Rule for the values: an array of the caller to an array of the caller
  goes through the operations. Everything else is done on the host by
  MGIS: values it keeps (scalars, local storage), and any copy between a
  value it keeps and an array of the caller, **in both directions**. These
  mixed copies are refused, with a message, if the location of the manager
  is not readable by the host. The operations therefore only receive
  arrays in the location of the manager: `FieldView` does not need to
  carry a location.
- Order and asynchronous mode: MGIS does its part on the host during the
  call, the copies delegated to the operations may start later on the
  queue. Preconditions, written in the documentation: no pending
  computation on the manager at the time of the call; synchronise the
  queue before reading on the host a copied array, or before a new
  `update` or `revert` if mixed copies are possible (shared location).
- View: state arrays (`updateValues` between two states of a view) and
  reset of K to zero. The material properties and the external state
  variables of a view belong to the caller.
- Operations provided as headers, compiled by the caller with its GPU
  tool, next to the memory resources (`MGIS/FieldOperations/`): CUDA and
  HIP (a stream), SYCL (a queue), Kokkos (an instance), OpenMP offload
  (synchronous, device number). Direct contiguous copy when both layouts
  are dense, a kernel otherwise.
- No swap of the states: it is not equivalent to `update` (a prediction
  does not write `s1`, nor do the failed points, mfem-mgis reads `s1` for
  the prediction residual), for a gain of one copy per time step.

### 3.5 Fixes of the existing code (step D1)

- Each asynchronous call allocates its status through the resource of the
  policy: with a CUDA managed resource, one `cudaMallocManaged` and one
  `cudaFree` per call. Standard remedy, without a new API: a
  `std::pmr::synchronized_pool_resource` on top of the managed resource,
  recommended by the documentation, used by the tests and the examples.
  The documentation recalls the precondition already written (the pending
  object lives until the end of the computation): with a pool, violating
  it corrupts the next call instead of being slow. The large blocks
  (matrices built from the axes) exceed the block size of the pool: they
  must be built once with `buildRotationMatrices`, as the documentation
  already advises.
- `evaluate` allocates two `std::vector` per call: replaced by a buffer on
  the stack.

### 3.6 Cost

- **TFEL kernels unchanged**: only the computation of the addresses
  changes with the layout.
- **Calls of the manager: no allocation in the common case.** `make_view`
  and the arrays of fields of the C structure are built in a buffer on the
  stack (`std::pmr::monotonic_buffer_resource`). Today, the manager
  allocates four `std::vector` per call.
- **View built by the caller**: allocation once, at its construction.
  Nothing per call.
- **Point by point computation**: nothing more for the current case;
  buffers of the workspace allocated once.
- **Checks**: one per variable and per field, on the host, per call.
- **End of time step**: one copy per state array and per time step (not
  per iteration), and one virtual call per copy.
- **To be measured, this is not an overhead**: point by point against
  component by component, depending on the behaviour and the machine.

## 4. Examples

mfem-mgis, all-GPU path, manager laid out component by component:

```cpp
auto i = MaterialDataManagerInitializer{};
i.layout = FieldLayout::COMPONENT_BY_COMPONENT;
i.location = MemoryLocation::DEVICE;
i.s0.gradients = std::span(g0, n * gsize);  // memory of MFEM
// ... all the arrays given
auto m = MaterialDataManager(b, n, i);
setExternalStateVariable(m.s0, "Temperature", std::span(T0, n),
                         MaterialStateManager::EXTERNAL_STORAGE);
setExternalStateVariable(m.s1, "Temperature", std::span(T1, n),
                         MaterialStateManager::EXTERNAL_STORAGE);
auto p = BatchExecutionPolicy{};
p.queue = stream;
const auto r = integrate(p, m, opts, dt);  // unchanged
auto ops = mgis::field_operations::CUDAFieldOperations{stream};
update(m, ops);  // end of time step on the GPU, in the order of the stream
```

Code that manages its data itself:

```cpp
auto v = MaterialDataView(b, n, MemoryLocation::DEVICE);
v.s0.gradients = makeComponentByComponentFieldView(std::span(g0, gsize * ld), n, gsize, ld);
v.s1.gradients = makeComponentByComponentFieldView(std::span(g1, gsize * ld), n, gsize, ld);
// ... forces, internal state variables, K
const auto E = makeUniformFieldView(std::span<const real>(E_d, 1), n);
setMaterialProperty(v.s0, "YoungModulus", E);
setMaterialProperty(v.s1, "YoungModulus", E);
auto p = BatchExecutionPolicy{};
p.queue = stream;
const auto r = integrate(p, v, opts, dt);
```

## 5. Symmetry to check before saying "done"

| Axis | Values |
|---|---|
| Object | behaviour, model, material property |
| Function | integration, initialisation, post-processing, rotations (3 kinds, 3 ways of giving the matrices), evaluation, `update`, `revert`, `updateValues`, extraction, finite strain conversions |
| Support | view, manager |
| Layout | point by point, component by component, arbitrary strides and `ld > n` (views) |
| Execution | without policy (workspace), synchronous, asynchronous, `ThreadPool` (manager) |
| Library | `generic-parallel` (host and GPU), `generic` |
| Points | all, range `[b, e)`, empty set, empty groups |
| Location | host, GPU, shared, expected refusals |
| Bindings and restarts | C, Fortran, Python; `save` and `restore` |

## 6. Set aside, with the reason

- **Layout by element blocks** (`byNODES` of `QuadratureInterpolator`):
  third stride and TFEL contract v3, without use (section 2).
- **Statuses per point**: the contract of TFEL has them for behaviours,
  not for material properties, with different conventions, and the point
  by point computation stops at the first failure. The current limitation
  stays for both. Possible addition later, together.
- **Swap of the states**: section 3.4.
- **`ThreadPool` on a view**: section 3.2.
- **`Context` variants of the computation functions**: the computation
  functions of the manager do not have them. Those that have them
  (assignments, `updateValues`) also have them on views.
- **Overlap between two different outputs**: not checked (costly in
  general), documented contract. The overlap of an output with itself is
  checked (3.1).
- **OpenMP memory resource**: OpenMP has no portable shared memory; the
  OpenMP end of time step operations exist, the resource does not
  (current situation unchanged).

## 7. What falls to mfem-mgis

- a component stride in `PartialQuadratureFunction`, or the *strided
  coalesced* views of `mgis::function`, for the outputs of a material laid
  out component by component;
- arrays readable by the GPU for `PartialQuadratureSpace`;
- in partial assembly, the handling of the elements of each material by
  the integrator (MFEM ignores the markers).

## 8. Steps

1. `FieldView`, `MaterialDataView`, functions on views, functions of the
   manager rewritten on top of them, explicit workspace, buffers on the
   stack. The 45 current tests unchanged and green.
2. Layout of the manager (3.3), `MaterialFunctionManager`, bindings.
3. `FieldOperations`, `update`, `revert` and `updateValues` with
   operations, CUDA, HIP, SYCL, Kokkos and OpenMP offload headers. Fixes
   of 3.5.
4. Tests of each cell of the symmetry matrix with a few behaviours
   (elasticity, Norton, orthotropic behaviour, `F-Cauchy` finite strain,
   model, behaviour with an initialize function and a post-processing,
   material property), on CPU (`stlpar`, `openmp`, `generic`) and CUDA,
   including `ld > n`, prediction, end of time step on GPU, and GPU without
   HMM. Bit for bit identical results between view and manager, and
   between the two layouts. A second build of MGIS with the bindings and
   HDF5 (disabled today) tests what is touched there. Web and Doxygen
   documentation of MGIS, documentation of the prototype, log.
5. Independent review, fixes, second review of the diff.

## Appendix: handling of the reviews

### Version 1

| Finding | Handling |
|---|---|
| A1 F initialised point by point | 3.3 |
| A2 workspace of a view | explicit argument (3.2) |
| A3 swap not equivalent to `update` | swap removed (3.4) |
| A4 `update` of the properties and external state variables | same procedure, copies through operations, written rule (3.4) |
| B1 TFEL and statuses per point | statuses removed for both (section 6) |
| B2 missing overloads | rotations by axes and in place, `evaluate` without bounds policy, policy of the view, setters and `updateValues` with `Context` |
| B3 statuses per point | removed |
| B4 layout of the outputs | everything given per point follows the layout (3.3) |
| B5 pitfalls of `FieldView` | empty view by default, extents, conversion to `const` |
| B6 allocations | buffers on the stack, pool for the asynchronous mode (3.5, 3.6) |
| B7 duplicates and analogies | three levels (3.2), regular types |
| B8 interoperability with `mgis::function` | one direction only, public access to the data, named accessors |
| B9 code that assumes point by point | layout exposed or refused in each binding (3.3) |
| B10 superfluous overloads | no `ThreadPool` on a view |
| C1 to C5 | empty set first, overlap of an output, layout in the manager, `make…` factories, `getComponents(view, start, count)`, `n` taken from the views |

### Version 2

| Finding | Handling |
|---|---|
| N1 workspace pointer in the view | removed: explicit argument, one per thread, behaviour checked, buffers sized at construction |
| N2 allocation of `make_view` | `std::pmr` vectors, buffer on the stack in the functions of the manager |
| N3 end of time step against `updateValues` | same procedure with copy operations, rule for the values kept by MGIS, `updateValues` with operations, check of the layout between states, K for a view |
| N4 overlap | exact check of a sufficient condition |
| N5 rotations | `std::span` overloads documented as point by point; buffer of one point for `generic` |
| N6 empty groups | written rule (3.1) |
| N7 pool | `synchronized_pool_resource`, precondition recalled, large blocks |
| N8 statements | confirmed true |
| N9 non-regular types | pointer to the behaviour, non-constant members, consistency checked at each call |
| N10 minor | `memory_location`, scalar material properties, `PreconditionsChecker`, parameter of the constructor of `MaterialStateManager`, numpy with strides, public host operations (`HostFieldOperations`), no more `get…` function that modifies, examples fixed |

### Version 3

| Finding | Handling |
|---|---|
| F1 end of time step: copy of a value kept by MGIS to an array of the caller, order in asynchronous mode, access check of `update(m)` | mixed copies on the host in both directions, refused outside memory readable by the host; order preconditions written; access check kept (3.4) |
| F2 `FieldView` derived from `PreconditionsChecker` not regular | no more inheritance, same mechanism through `check_preconditions` (3.1) |
| F3 public `make_view` bypasses `checkHostValues` | check in `make_view` (3.2) |
| F4 semantics of the `std::pmr` vectors | written (3.2) |
| F5 asymmetries: `updateValues` on a view, extraction and conversions on a view, `IntegrationType` forms | added, or written as reserved to the manager (3.2, 3.4) |
| F6 implementation details: routing before the workspace, `K[0]`, buffers of `s0` and `s1`, strides of the evaluators, buffer with a size that depends on the behaviour, pointer to the behaviour in the workspace, non-virtual interface of `FieldOperations` | written (3.2, 3.4) |
