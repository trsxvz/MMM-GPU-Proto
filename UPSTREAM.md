# Upstream integration plan

The branch of TFEL (`tfel/`) mixes fixes of existing code, device code
support and the new interface. It will be proposed upstream as a series of
pull requests, each one reviewable and testable on its own. Some commits mix
several concerns: the pull requests will be rebuilt from the final state of
the branch, file by file, rather than by cherry-picking.

## Order of the pull requests

1. **Fixes of existing code**, without device code:
   - linker flags after the object files;
   - fixes of the `generic-parallel` interface for material properties:
     code generated with the runtime checks disabled, device memory of
     the `SYCL` backend freed after each call, pointers to the arguments
     copied in the kernels (the array given by the calling code may not be
     accessible by the device), extents of the strided fields and values
     between the strided outputs preserved (`SYCL`), data transfers inside
     the `try` block, backend in the header guard, bounds and default
     values of the parameters written with 14 digits, grid size of the
     `CUDA` kernels;
   - fixes of the drivers of the `generic` interface: status tested after
     a finite strain integration, moduli of the `DTAU_DDF` operator with the
     Hencky strain, state at the end of the time step in the
     post-processings, initialize functions of the Green-Lagrange and
     Hencky strategies, speed of sound with the prediction operator,
     pointers restored after an exception, stress of the initialize
     functions for all stress measures, exceptions of the constructor of
     the behaviour reported as failures;
   - fixes of the rotation functions of the `generic` interface (offsets,
     scalars, vectors);
   - bounds checks of the generic models in all compilation modes;
   - generated classes of the generic models given internal linkage: the
     parameters of a model were shared by all the libraries of this model
     loaded in a process (static variable of an inline function, emitted
     as a unique symbol);
   - flow rule of the `DigitSeparatorTest` test behaviour, which read an
     uninitialized derivative;
   - physical bounds of the arguments of the material properties exported
     even if those arguments have no standard bounds.

   Those fixes may change the results on CPUs. Each one is covered by an
   existing test or by a new one.

2. **Device code support in TFEL** (mathematical and material headers):
   host and device annotations, `TFEL_DEVICE_CODE` for single pass
   compilers, behaviours without virtual functions, `TrivialGenTypeBase`,
   `getConstantObject`, eigenvectors without recursion in device code,
   compiler workarounds (for instance, iterators of the cubic splines
   advanced without `std::advance` and `std::prev`, which nvcc compiles as
   functions doing nothing in device code). The results on CPUs are unchanged. The API
   changes are listed in the release notes.

3. **Device compatible generated code in MFront**: annotations of the
   generated classes, device compatible behaviour class (interface
   property and DSL option), drivers of the `generic` interface split in a
   host and device core and the historical host functions, rotation
   functions and method callers usable on GPUs, type of the tangent
   operator computed by the user checked before its conversion (the
   `get` method of `TrivialGenTypeBase` only checks it on the host). The
   results of the existing interfaces are unchanged.

4. **The `generic-parallel` interface for behaviours** with a first set of
   backends: CUDA, HIP, C++ parallel algorithms and OpenMP on the host.
   Tests and documentation.

5. **Additional backends and compilers**: Kokkos, SYCL, nvc++ offloading
   (`stdpar` and OpenMP), clang CUDA and HIP. Most compiler workarounds
   belong to this pull request.

6. **Extensions of the interface for material properties**: Kokkos and
   OpenMP backends, material laws using data or other material laws on
   CUDA and HIP, options shared with the behaviours.

The data contract (`MFront/GenericParallelBehaviour/BatchData.h`) is
versioned from the fourth pull request on, so that MGIS can check it.

## Checks for each pull request

- the full test suite of TFEL (`ctest`);
- the comparisons of the `generic-parallel` tests for the configurations
  of the pull request, at `-O3`;
- no change of the results of the existing interfaces, except for the
  fixes of the first pull request;
- formatting by `clang-format` limited to the modified lines.
