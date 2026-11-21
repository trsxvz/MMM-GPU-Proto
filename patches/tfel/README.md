# Patches

Series of patches giving this tree from TFEL master, commit `753cd76da`
(20/11/2026):

```bash
git clone <official TFEL repository> tfel && cd tfel
git checkout 753cd76da
git am --keep-cr /path/to/patches/*.patch
```

`--keep-cr` keeps the CRLF line endings of some files of TFEL (test
behaviours), without which some patches do not apply.

The series is written by subject: each patch gathers the changes of one
subject, whatever the order in which they were developed (the history of
the developments is in the journal of the prototype; the details and their
justifications are in `docs/tfel.md`). Each changed file belongs to one
patch, except the lists of `mfront/include/CMakeLists.txt` and
`mfront/src/CMakeLists.txt`, whose lines are given by the patch of the
files they name.

The `generic-parallel` interface for material properties of the master of
TFEL is the reviewed version of the one of the prototype (pull requests
#1088 and #1415: version 1 of the functions, runtime of the backends in
`MFront/GenericParallel`). The series is based on it and only adds the
behaviours and the models on the same runtime, the additional backends,
the libraries also containing the functions of the `generic` interface,
the numerical type of the generic interfaces and the material laws called
in device code.

1. TFEL, device code support in the mathematical and material libraries:
   `TFEL_DEVICE_CODE` macro and device code selection for single pass
   compilers, host and device annotations (checked by compiling each
   header with clang CUDA), portable error reporting (`tfel::fail`,
   contract violations), eigenvalues and eigenvectors without recursion
   nor assertion (Jacobi algorithm for multiple eigenvalues), finite strain
   tools, rotations, cubic splines, behaviours without virtual functions,
   slip systems built at compile time (`getConstantObject`), trivially
   copyable tensors (`TrivialGenTypeBase`, values built in place),
   `IsAssignableTo` for floating
   point types, size of the numerical type of a library read by
   `ExternalLibraryManager`, workarounds of compiler defects.
2. MFront, generated code compatible with device code and numerical type
   of the generic interfaces: integration drivers of the `generic`
   interface split into device capable cores, rotation functions,
   behaviour variables, models and porous flows of the
   StandardElastoViscoPlasticity brick, data shared by all the integration
   points (`@SharedData`), static variables written
   `TFEL_STATIC_VARIABLE(v)`, comparison to a numerical jacobian and
   profiling on the host only; DSL option `real_type` (symbols
   `<f>_real_size`, API versions shifted for `float`, other interfaces and
   MTest refusing `float`, convergence criteria bounded at the generation
   with warnings); displacement gradient given to the finite strain
   behaviours by the `generic` interfaces and axial deformation gradient
   of plane stress computed by the behaviour; fixes of the existing code,
   which may change the results on CPU (bounds checks of the generic
   models, Poisson ratio of the stress update algorithm of the isotropic
   DSLs, DDIF2 brick, error normalization factor of `rkCastem`).
3. `generic-parallel`, runtime shared by the behaviours, the models and
   the material properties: the backends of the interface for material
   properties also write the code of the behaviours and of the models
   (runtime of the interface, shared data copied on each device, names of
   the generated files, function object of the `CUDA` kernels declared
   `__grid_constant__` instead of being copied in the stack of each thread),
   libraries of material properties also containing the function of the
   `generic` interface (`generic_functions` option).
4. MFront, `generic-parallel` interface for behaviours and models:
   integration of a batch, initialize functions, post-processings,
   rotation of all the tangent operator blocks, standalone models, finite
   strain behaviours, values of a point copied with all the reads before the
   writes, versioned C data structures (`MFRONT_GPB_API_VERSION`
   4), output status, statistics and execution policies of the interface
   for material properties, report function exported by each library,
   errors of the batch `-5`, `-6` and `-7`.
5. Kokkos backend: `Kokkos::parallel_for` in the execution space given as
   an option, statistics in the memory space of the execution space,
   `kokkos` and `kokkos-host` configurations.
6. RAJA backend: execution policy given as an option and checked before
   being copied in the generated code, resource of camp as queue,
   statistics in device memory, `raja` and `raja-host` configurations.
7. OpenMP backend: offloading (`target teams distribute parallel for`) or
   threads of the host, `openmp-offload-nvhpc-gpu` and `openmp-host-gcc`
   configurations.
8. OpenACC backend: `parallel loop` with reductions, memory space deduced
   from the compilation macros, stack size set through the CUDA driver,
   `openacc-nvhpc-gpu` and `openacc-nvhpc-multicore` configurations.
9. Arcane backend: accelerator API of Arccore, the core component of the
   Arcane framework, with a queue selecting the GPU or the cores of the
   host at runtime, atomic operations of Arccore, statistics in the memory
   of the device of the queue, made the current device during the call,
   `arcane` configuration (Arccore of the main branch of Arcane from its
   commit `51b90b8`, built with CUDA, clang as a CUDA compiler).
10. CMake, configurations of the `generic-parallel` interfaces: detection
    of the compilers and of the libraries (Kokkos, RAJA, Arccore, nvc++,
    path of CUDA for clang), flags of the configuration files (`nvc++
    -stdpar=gpu` without relocatable device code, optimisation flags,
    precise floating point model of `icpx`), components of `tfel-check`.
11. Tests of the `generic-parallel` interfaces against the `generic`
    interface for all the configurations (behaviours, models, TDLS,
    isotropic DSLs, parameters changed at run time, static variables,
    libraries also containing the functions of the `generic` interface,
    material property calling another material property, additional
    backends).
12. Documentation: page of the `generic-parallel` interface for behaviours,
    additional backends and options of the interface for material
    properties, pages of the keywords and of the DSL options, release
    notes.
