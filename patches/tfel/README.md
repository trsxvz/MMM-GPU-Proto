# Patches

Series of patches giving this tree from TFEL master, commit `ff0959ece`
(12/10/2026):

```bash
git clone <official TFEL repository> tfel && cd tfel
git checkout ff0959ece
git am --keep-cr /path/to/patches/*.patch
```

`--keep-cr` keeps the CRLF line endings of some files of TFEL (test
behaviours), without which some patches do not apply.

Patches, by subject (independently of their numbers; the details and
their justifications are in `docs/tfel.md` of the prototype, in French):

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
  standard parallelism, Kokkos, SYCL and OpenMP (host and offloading)
  backends, standalone models, finite strain behaviours, initialize
  functions, post-processings, rotation of all the tangent operator
  blocks, versioned C data structures, symbols of the interface, of the
  backend and of the memory space, execution options shared by the
  backends, asynchronous execution, error reporting without cost on
  success, shared data copied on each device, statistics;
- `generic-parallel` interface for material properties: fixes, Kokkos and
  OpenMP backends, properties using data or other laws, execution options
  as last argument (version 2 of the functions), asynchronous evaluation,
  statistics, symbols, statuses in reused blocks, status of the `Warning`
  policy, physical bounds;
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
  options), error normalization factor of `rkCastem`.
