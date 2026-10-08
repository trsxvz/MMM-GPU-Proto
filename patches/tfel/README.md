# Patches

Series of patches giving this tree from TFEL master, commit `fec0ade99`
(merge of pull request #1040, 01/10/2026):

```bash
git clone <official TFEL repository> tfel && cd tfel
git checkout fec0ade99
git am /path/to/patches/*.patch
```

Patches:

1. device code support in the mathematical and material headers of TFEL;
2. linker flags after the object files in the generated Makefiles;
3. generated code of the behaviours compatible with device code;
4. `generic-parallel` interface for behaviours (CUDA, HIP and C++ standard
   parallelism backends);
5. tests and documentation;
6. fixes of the `generic-parallel` interface for material properties;
7. HIP backend tested on NVIDIA GPUs, `hip-hipcc` configuration and test;
8. behaviours without virtual functions (`MechanicalBehaviourConcept`);
9. device code selection for single pass compilers (`nvc++ -stdpar=gpu`,
   AdaptiveCpp);
10. no recursion in the computation of the eigenvectors of stensors;
11. Kokkos and SYCL backends, GPU offloading with `nvc++`, stack size of the
    threads;
12. device capable integration drivers of the `generic` interface (cores
    taking the computation options, used by both interfaces);
13. device code support of the finite strain tools, rotations of vectors and
    matrices, constant objects (slip systems) and trivially copyable
    tensors (`TrivialGenTypeBase` instead of `std::variant`);
14. rotation functions shared by the `generic` and `generic-parallel`
    interfaces, initialize functions and post-processings called through
    function objects, slip systems built at compile time;
15. `generic-parallel`: finite strain behaviours, initialize functions,
    post-processings, rotations, extended tests.
16. consistent host and device annotations in TFEL (checked by compiling
    each header with clang CUDA), device code of the porous criteria and
    of the non linear solvers;
17. behaviour variables, models and porous flows of the
    StandardElastoViscoPlasticity brick usable in device code, fixes of the
    generic model interface;
18. `generic-parallel`: standalone models, tests of TDLS, behaviour
    variables and of the StandardElastoViscoPlasticity brick.
19. device code selected in the whole translation unit with `nvc++ -mp=gpu`;
20. bounds checks of the generic models reported in all compilation modes;
21. `generic-parallel`: OpenMP backend, Kokkos backend for material
    properties, entry points shared by the backends;
22. no assertion in the computation of the eigenvectors of 3D stensors;
23. function objects of the tangent operator of the implicit DSLs owning
    their status flag;
24. `generic-parallel`: OpenMP offloading at `-O2`, tests of the
    behaviours with clang and icpx;
25. tag names of the C structures of the batch data;
26. material properties and external state variables given variable by
    variable;
27. symbols identifying the interface and the backend;
28. integration drivers of the generic interface working on copies of the
    behaviour data;
29. stress computed by the initialize functions exported for all stress
    measures;
30. formatting.

Latest patches (0070 to 0078):

- shared data copied on each device, backend and version symbols and
  shared kernel code of the material properties;
- rotation of all the tangent operator blocks, parameters of the behaviour
  variables in device code, `nvc++ -stdpar=gpu` without relocatable device
  code;
- error reporting without cost on success (index of the first point where
  the computations failed, reused status blocks), status of the warning
  policy of the material properties.
