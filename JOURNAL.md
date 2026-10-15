# Journal

This log tells the work day by day: what was done, the defects found, the
decisions and their reasons. Each entry describes the state at the time it
was written. Some choices may have changed afterwards: the current state is
described in [docs/](docs/overview.md).

## Phases

| Marker | Topic |
|---|---|
| Milestones 1 to 3 | generated code that compiles for the GPU, `generic-parallel` interface of the behaviours, first tests |
| B1 to B7 | coverage of MFront: finite strains, initialize functions, post-processings, rotations, single crystals, models, TDLS |
| C1 to C5 | backends: HIP, `stlpar` with nvc++, Kokkos, SYCL, OpenMP |
| R1 | fixes from the review of 06/10 |
| R2 | test behaviours portable between CPU and GPU |
| R3 | execution trials of all the test behaviours |
| R4 | validations at `-O3`, spack environment |
| R5 | clang-tidy and independent review |
| R6, R7 | data shared per GPU, parameters of the behaviour variables, rotation of all the blocks |
| R8 | errors at no cost on success |
| R9 | memory space of the libraries |
| R10 | execution options and asynchronous evaluation of the material properties |
| R11 | documentation of TFEL and MGIS up to date |
| R12 | statistics of the material properties modelled on the behaviours |
| R13 | final audit of TFEL and MGIS, fixes, second review |
| D1 | MGIS: batches of behaviours, models and material properties |
| D2 | mfem-mgis: choices made, work to come ([docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md)) |

## Conventions

- The paths quoted are those of the time, with the names given on
  14/10/2026, when the repository was translated into English: files,
  directories, scripts, executables, tests, environment variables, local
  build directories and local branches were renamed then (see
  `cases/README.md`). Several documents have also been moved:
  - `DESIGN.md` became [docs/tfel.md](docs/tfel.md);
  - `UPSTREAM.md` became [docs/upstream-integration.md](docs/upstream-integration.md);
  - `scratch/review-tfel-2026-10-06.md` is in
    [docs/history/](docs/history/review-tfel-2026-10-06.md);
  - `scratch/d1-mgis-design.md` is in
    [docs/history/](docs/history/design-mgis-2026-10-08.md);
  - `scratch/d2-mfem-mgis-design.md` became
    [docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md).
- `tfel/` is the TFEL repository of the prototype, `tfel-dev/` its
  development copy (same repository, another working tree).
- The time measurements are those of the development machine (RTX 5080
  GPU).

## Index of entries

| Entry |
|---|
| [2026-10-05, state of TFEL master](#2026-10-05-state-of-tfel-master) |
| [2026-10-05, J2 probe (default DSL) in a CUDA kernel](#2026-10-05-j2-probe-default-dsl-in-a-cuda-kernel) |
| [2026-10-05, milestone 1 done: generated code that compiles for the GPU](#2026-10-05-milestone-1-done-generated-code-that-compiles-for-the-gpu) |
| [2026-10-05, milestone 2: `generic-parallel` interface of the behaviours](#2026-10-05-milestone-2-generic-parallel-interface-of-the-behaviours) |
| [2026-10-05, milestone 3 and sweep of the test behaviours](#2026-10-05-milestone-3-and-sweep-of-the-test-behaviours) |
| [2026-10-05, material laws, eigenvalues and stack on GPU](#2026-10-05-material-laws-eigenvalues-and-stack-on-gpu) |
| [2026-10-05, reorganisation, local commits, environment pitfall](#2026-10-05-reorganisation-local-commits-environment-pitfall) |
| [2026-10-05, HIP on NVIDIA GPU (C1, commit 0007)](#2026-10-05-hip-on-nvidia-gpu-c1-commit-0007) |
| [2026-10-05, `stlpar` on GPU with nvc++ 25.11 (C2, in progress)](#2026-10-05-stlpar-on-gpu-with-nvc-2511-c2-in-progress) |
| [2026-10-05, Kokkos and SYCL backends, factorisation, CMake fixes](#2026-10-05-kokkos-and-sycl-backends-factorisation-cmake-fixes) |
| [2026-10-05, full matrix of the backends, fixes](#2026-10-05-full-matrix-of-the-backends-fixes) |
| [2026-10-05 (evening), B5 finite strains, B2 initialisations and post-processings](#2026-10-05-evening-b5-finite-strains-b2-initialisations-and-post-processings) |
| [2026-10-05 (evening, continued), B3 rotations of the orthotropic behaviours](#2026-10-05-evening-continued-b3-rotations-of-the-orthotropic-behaviours) |
| [2026-10-05 (night), B4 single crystals, defects found by nvcc](#2026-10-05-night-b4-single-crystals-defects-found-by-nvcc) |
| [2026-10-05 (night, continued), B1 numerical coverage, new nvcc and clang defects](#2026-10-05-night-continued-b1-numerical-coverage-new-nvcc-and-clang-defects) |
| [2026-10-05 (night, continued), B6 standalone models](#2026-10-05-night-continued-b6-standalone-models) |
| [2026-10-05 (night, continued), B7 TDLS solver on GPU](#2026-10-05-night-continued-b7-tdls-solver-on-gpu) |
| [2026-10-05 (night, continued), behaviour variables and full sweep](#2026-10-05-night-continued-behaviour-variables-and-full-sweep) |
| [2026-10-05 (night, continued), StandardElastoViscoPlasticity brick and headers of TFEL](#2026-10-05-night-continued-standardelastoviscoplasticity-brick-and-headers-of-tfel) |
| [2026-10-05 (night, continued), commits 0016 to 0018](#2026-10-05-night-continued-commits-0016-to-0018) |
| [2026-10-05 (evening), C5 OpenMP backend, factored entry points](#2026-10-05-evening-c5-openmp-backend-factored-entry-points) |
| [2026-10-05 (evening, continued), material properties: Kokkos and OpenMP backends](#2026-10-05-evening-continued-material-properties-kokkos-and-openmp-backends) |
| [2026-10-05 (evening, continued), matrix c, OpenMP, models compiled by nvc++](#2026-10-05-evening-continued-matrix-c-openmp-models-compiled-by-nvc) |
| [2026-10-06 (night), versions of NVHPC and CUDA, defects of OpenMP offload](#2026-10-06-night-versions-of-nvhpc-and-cuda-defects-of-openmp-offload) |
| [2026-10-06 (night, continued), fields per variable, symbols, roles](#2026-10-06-night-continued-fields-per-variable-symbols-roles) |
| [2026-10-06 (night, continued), drivers of the generic interface, OpenMP in direct mode](#2026-10-06-night-continued-drivers-of-the-generic-interface-openmp-in-direct-mode) |
| [2026-10-06 (morning), clang CUDA and HIP configurations](#2026-10-06-morning-clang-cuda-and-hip-configurations) |
| [2026-10-06 (morning, continued), clang CUDA sweep of all the tests of MFront](#2026-10-06-morning-continued-clang-cuda-sweep-of-all-the-tests-of-mfront) |
| [2026-10-06 (morning, continued), DSLs not run on GPU and two defects of nvc++ OpenMP](#2026-10-06-morning-continued-dsls-not-run-on-gpu-and-two-defects-of-nvc-openmp) |
| [2026-10-06 (morning, end), validation of the merged branch](#2026-10-06-morning-end-validation-of-the-merged-branch) |
| [2026-10-06 (afternoon), fixes from the review (R1)](#2026-10-06-afternoon-fixes-from-the-review-r1) |
| [2026-10-06 (evening), portable test behaviours (R2)](#2026-10-06-evening-portable-test-behaviours-r2) |
| [2026-10-07 (night), execution trials of all the test behaviours (R3)](#2026-10-07-night-execution-trials-of-all-the-test-behaviours-r3) |
| [2026-10-07 (morning), large validations at -O3 (R4, first part)](#2026-10-07-morning-large-validations-at--o3-r4-first-part) |
| [2026-10-07 (afternoon), switch to the mmm-gpu spack environment (R4, in progress)](#2026-10-07-afternoon-switch-to-the-mmm-gpu-spack-environment-r4-in-progress) |
| [2026-10-07 (evening), mmm-gpu environment validated (end of R4)](#2026-10-07-evening-mmm-gpu-environment-validated-end-of-r4) |
| [2026-10-08 (night), clang-tidy and independent review (R5)](#2026-10-08-night-clang-tidy-and-independent-review-r5) |
| [2026-10-08 (night, continued), shared data, symbols, duplication (R6)](#2026-10-08-night-continued-shared-data-symbols-duplication-r6) |
| [2026-10-08 (night, end), parameters of the behaviour variables, rotations (R7)](#2026-10-08-night-end-parameters-of-the-behaviour-variables-rotations-r7) |
| [2026-10-08 (morning), uniform error handling (R8)](#2026-10-08-morning-uniform-error-handling-r8) |
| [2026-10-08 (morning, continued), MGIS: loading and computation of batches (D1, steps 1 to 3)](#2026-10-08-morning-continued-mgis-loading-and-computation-of-batches-d1-steps-1-to-3) |
| [2026-10-08 (afternoon), memory space (R9) and review of MGIS (D1, step 4b)](#2026-10-08-afternoon-memory-space-r9-and-review-of-mgis-d1-step-4b) |
| [2026-10-08 (afternoon, continued), MGIS: material properties and batch models (D1, step 6)](#2026-10-08-afternoon-continued-mgis-material-properties-and-batch-models-d1-step-6) |
| [2026-10-08 (evening), execution options of the material properties (R10)](#2026-10-08-evening-execution-options-of-the-material-properties-r10) |
| [2026-10-08 (evening, continued), asynchronous evaluation of the material properties (R10, continued)](#2026-10-08-evening-continued-asynchronous-evaluation-of-the-material-properties-r10-continued) |
| [2026-10-08 (evening, end), SYCL queues and thread blocks of the material properties (R10, end)](#2026-10-08-evening-end-sycl-queues-and-thread-blocks-of-the-material-properties-r10-end) |
| [2026-10-08 (night), documentation of TFEL and MGIS (R11)](#2026-10-08-night-documentation-of-tfel-and-mgis-r11) |
| [2026-10-08 (night, continued), asynchronous computation in MGIS (D1, step 5)](#2026-10-08-night-continued-asynchronous-computation-in-mgis-d1-step-5) |
| [2026-10-08 (night, end), first failing point of the material properties (R12)](#2026-10-08-night-end-first-failing-point-of-the-material-properties-r12) |
| [2026-10-08 (night, closing), final audit of TFEL and MGIS (R13)](#2026-10-08-night-closing-final-audit-of-tfel-and-mgis-r13) |
| [2026-10-08 (end), choices for mfem-mgis (D2)](#2026-10-08-end-choices-for-mfem-mgis-d2) |

## Entries

### 2026-10-05, state of TFEL master

- `generic-parallel` interface of the material properties (#966):
  CUDA, HIP, standard C++, SYCL backends; two prototypes (with and without
  memory stride); parameters in a POD structure passed by value to the kernel;
  errors in an array; exported metadata; tests in the installation
  directory (`test-suite/material-properties/generic-parallel`).
- Behaviours: the generated code is already almost GPU compatible
  (`TFEL_HOST_DEVICE` methods, no exception). Branch
  `411-introduce-a-derived-generic-interface-for-gpus`: draft from 2023,
  abandoned.

### 2026-10-05, J2 probe (default DSL) in a CUDA kernel

Errors or warnings found, in order:
1. `integrate` generated by the default DSL not annotated `TFEL_HOST_DEVICE`.
2. Constructors from `mfront::gb::BehaviourData` (behaviour, data, integration
   data) not annotated.
3. Virtual destructor of `MechanicalBehaviour` (TFEL) not annotated, so the
   generated destructor (`= default`) is host only.
4. Parameters read through the singleton `...ParametersInitializer::get()` (host).
5. Checks `raise_if(smflag != STANDARDTANGENTOPERATOR, ...)` generated
   in `integrate` and `computePredictionOperator`.

`tfel::reportContractViolation` does not fit for 5: it calls `abort`
on the host, which would change the CPU behaviour. Decision: macro
`TFEL_DEVICE_CODE` (compilation of the GPU code only) and a branch in the
generated code, the CPU staying identical.

### 2026-10-05, milestone 1 done: generated code that compiles for the GPU

- `integrate`, `computePredictionOperator`, constructors and destructor
  annotated `TFEL_HOST_DEVICE`; "impossible" errors (`raise`) replaced
  by `writeDeviceCompatibleRaise`: exception on the host (unchanged), failure
  returned on the GPU.
- Parameters: the constructor of an interface can read the parameters from
  an argument (`getParametersSourceInBehaviourConstructor`) instead of the
  host singletons.
- Bounds (`BoundsCheck.hxx` rewritten, same API): `TFEL_HOST_DEVICE`
  functions that return `bool`. Host: exception or warning
  as before. GPU: `false` if the policy is strict. `checkBounds()`
  returns `bool`, the physical bounds checked after integration do
  `return FAILURE;` on GPU. The three specialisations `BoundsCheck<1,2,3>`
  replaced by a single template.
- The generator writes only once the identical blocks produced by the
  interfaces (constructors of the data, `exportStateData`...): `generic` and
  `generic-parallel` can be used together.

### 2026-10-05, milestone 2: `generic-parallel` interface of the behaviours

- Design (see DESIGN.md): same skeleton as the interface of the
  material properties (factory of backends, options in the
  configuration files, library `GenericParallel<lib>-<backend>`,
  symbols exported with the `GenericParallel` interface). Difference: the
  per-point logic is written once in `TFEL_HOST_DEVICE` headers
  (`MFront/GenericParallelBehaviour/Integrate.hxx`, `GPU.hxx`,
  `ParallelSTL.hxx`); the code generated per modelling hypothesis is about
  fifteen lines (structure of the parameters, sizes, call).
- `Integrate.hxx`: equivalent without exception of `mfront::gb::integrate`.
  Two paths: direct access if all the fields are stored point by point,
  otherwise copy into local arrays then copy back (storage by
  component, coalesced accesses on GPU).
- `GPU.hxx`: one thread per point; summary (failures, unreliable results,
  minimal `rdt`) reduced per warp, atomics only if a point of the warp has
  something to report; `cudaMallocAsync` on the stream of the caller.
- Return values: 1 success, 0 `rdt` < 0.99 somewhere, -1 failure at
  least at one point, -2 invalid arguments, -3 execution error.

Pitfalls found:
1. nvcc cannot generate the host stub of a kernel whose non-type template
   parameter is a class (C++20): sizes passed through a type.
2. **`std::string_view` built from a C string in GPU code
   (libstdc++ 14 + nvcc 13) empties the whole kernel, without error or
   warning** (`__builtin_strlen` treated as unreachable). The
   Norton kernels used 4 registers and computed nothing.
   `BoundsCheck` now takes `const char*`. Only the numerical
   comparison detected it: a check to keep (number of registers, tests).
3. libstdc++: `std::for_each(par_unseq, ...)` on `std::views::iota`
   runs on **a single thread** (the C++17 category of the iterator is not
   "random access"). The `stlpar` backend of the material properties of master
   is therefore sequential with gcc. Backend of the behaviours: counting iterator.
4. MFront placed `$(LDFLAGS)` (hence `-ltbb`) before the objects: TBB ignored
   with `--as-needed` (Ubuntu), undefined TBB symbols. Fixed in
   `MakefileGenerator.cxx` (master has the same defect).
5. CUDA test of the material properties of master: one UO2 case fails here
   (`errno` non zero after the call), unrelated to our changes.
6. nvcc warnings 2810 (`nodiscard` on `operator=` of the types of
   TFEL/Math, false positives), 20208, 3189: suppressed in
   `config-cuda-nvcc.json`, with `-arch=native` added.

Indicative measurements (RTX 5080, double precision capped at 1/64), 100,003
points, full call (allocation, copy, synchronisation): J2 0.3-0.4 ms,
Norton 6-7 ms; `stlpar` 16 threads: J2 1.6-2.6 ms, Norton 17-21 ms.

### 2026-10-05, milestone 3 and sweep of the test behaviours

- `tfel-check` tests: `compare.cxx` compilable by nvcc (managed memory)
  or g++; both interfaces generated by the same call of `mfront`
  (also checks that they coexist). `tfel-check` components added in
  `cmake/modules/tfel.cmake` (`mfront::behaviour::generic-parallel_...`).
- Sweep: each behaviour compiled with the CUDA backend, classified OK if no
  call of a host function from the GPU code (nvcc warning 20011) and if
  no kernel has fewer than 16 registers. **A host call is not a
  mere warning: nvcc removes the corresponding code path**
  (empty kernel or wrong computation without error). 99 behaviours clean on
  the first pass, 120 after fixes.
- Fixes of the generator: Runge-Kutta DSLs (`integrate`, derivatives,
  divergence exceptions → `return FAILURE` on GPU, minimal time
  step), `MultipleIsotropicMisesFlows`, isotropic DSLs (`computeFlow`,
  `NewtonIntegration`, tangent operator), free expansion
  (`computeStressFreeStrain/Expansion`), Hooke potential (elastic
  prediction), accessors of blocks of the jacobian and of the tangent operator,
  wrappers of the internal behaviours (`@BehaviourVariable`), qualified
  `STANDARDTANGENTOPERATOR` flag (a regression I had introduced on the
  general behaviours), profiler disabled on GPU, dynamic arrays forbidden if
  an interface does not support them (new method
  `areDynamicallyAllocatedVectorsAllowed` of the interfaces), "static
  variables" parameters: non static members in GPU code (nvcc refuses
  the ODR-use of a `static constexpr` member, and the expressions of TFEL/Math
  take the scalars by reference). `FSAlgorithm` annotated.
- Non portable user code (19 behaviours): `throw`, `std::cout`,
  `raise_if`, `static const` initialised at run time. On GPU, use
  `return false;` or `return FAILURE;`. To be documented.

### 2026-10-05, material laws, eigenvalues and stack on GPU

- Material laws called by the behaviours (`@MaterialLaw`,
  `@ElasticMaterialProperties`, `@ComputeThermalExpansion`...):
  the `mfront` interface of the material properties also writes, in its
  header, a `TFEL_HOST_DEVICE inline` implementation in the namespace
  `mfront::device` (same body, constant parameters, no exception
  and no `errno`, `NaN` if a physical bound is not respected).
  `writeMaterialLaws` writes `using mfront::device::Law;` in GPU code; the
  blocks that check the bounds of the law (exceptions, `std::cerr`) are
  compiled for the host only. Models (`@AxialGrowth`...): class already
  in the header, annotated. Check: StrainHardeningCreep3 (Lamé
  coefficients computed by `Inconel600_YoungModulus(T)`, physical quantities)
  bit for bit identical to the `generic` interface.
- Eigenvalues: the 8 solvers of `stensor` (TFEL, FSES ×5, GTE,
  Harari) annotated and checked on GPU (gap to the CPU ≤ 3e-15,
  reconstruction error ≤ 1e-13). **New silent pitfall: `std::fpclassify`
  in GPU code empties the kernel** (like `std::string_view` from a C
  string). FSES used it: replaced by `tfel::math::ieee754::fpclassify`,
  which works on GPU.
- Stack of the CUDA threads: 1024 bytes by default, not increased
  automatically; Hosford (real power, `pow` of libdevice) needs
  2160 bytes (static estimate, the real usage is a bit larger:
  threshold measured between 2160 and 2220). The GPU backend reads the need of the kernel
  (`cudaFuncGetAttributes`, `localSizeBytes`) and increases the limit if
  needed (+512 bytes, rounded to the kB, under a lock because the limit is
  global to the context). Hosford identical to `generic` (≤ 6e-13).
- `tfel-check` tests extended: + Hosford, + StrainHardeningCreep3 (material
  law on GPU, quantities). Both configurations pass (CUDA
  1 min 31, CPU 32 s). The bench now also compares the statuses and the
  summary of the first step (instead of requiring its success).
- Documentation: `docs/web/generic-parallel-behaviour-interface.md`.

### 2026-10-05, reorganisation, local commits, environment pitfall

- Directories `tfel/` (formerly `tfel-gpu-proto`), `mgis/`, `mfem-mgis/`, each
  with `patches/`. 5 local commits in `tfel/`, series
  `tfel/patches/0001..0005`. Full recompilation of TFEL after the
  renaming (absolute paths in the CMake cache).
- **Environment pitfall**: the `LD_LIBRARY_PATH` and the `PATH` of the shell
  contain the spack view `tfelgpu_local` (older TFEL master,
  `5696f406f`). `LD_LIBRARY_PATH` comes before the RUNPATH of the binaries:
  the `mfront` of the build tree loaded the libTFELMFront of
  spack. The first pass of the test suite of TFEL therefore tested
  the old MFront (errors BarBehaviour2/4 and TDLS unrelated to our
  changes). Rule: for the build tree,
  `env -u LD_LIBRARY_PATH PATH=<without spack>`; for the installation,
  `LD_LIBRARY_PATH=$P/install/tfel/lib:...` and `PATH=$P/install/tfel/bin:...`
  at the front (which is what the scripts of the proto do).
- Seen along the way: `~/Documents/GitHub/mgis-gpu` (prototype
  from 01-02/2026): CUDA kernels written by hand on the
  `StridedCoalescedMemoryAccess` views of MGIS. The format of the fields with two
  strides of the interface matches these views (useful for D1).
- Test suite of TFEL (validated command, ~2 min 15 after
  compilation): in `build/tfel`,
  `CLEANPATH=$(echo "$PATH" | tr ':' '\n' | grep -v spack | paste -sd:)`,
  `touch mfront/src/mfront` (if MFront changed),
  `env -u LD_LIBRARY_PATH PATH="$CLEANPATH" make -k -j16 build-tests`, then
  `env -u LD_LIBRARY_PATH PATH="$CLEANPATH" ctest -j16`. Result:
  6963/6963.
- Open question (`generic` interface, not modified): the post-processing
  only replaces the thermodynamic forces of the initial state by those
  of the end of the step, whereas the comment also mentions the internal
  variables; and it passes `*(nullptr)` as the initial state when it is
  not used (undefined behaviour). The GPU version will pass a
  valid reference and will keep the semantics of `generic`.
- Error fixed: TFEL has a tracked `patches/` directory (homebrew patches);
  the series of the prototype had replaced it. Restored; the series of the proto is in
  `proto-patches/` (in the three libraries).

### 2026-10-05, HIP on NVIDIA GPU (C1, commit 0007)

- `hipcc` on the NVIDIA platform calls `nvcc` and defines neither `__HIP__`
  nor `__HIPCC__`: `GPU.hxx` silently took the CUDA path. The
  `hip` backend now includes `hip/hip_runtime.h` before `GPU.hxx`
  (like the HIP backend of the material properties), which defines
  `__HIPCC__`; selection in `GPU.hxx`: `__HIPCC__ || __HIP__`.
- CUDA/HIP differences fixed: warp size attribute
  (`cudaDevAttrWarpSize` / `hipDeviceAttributeWarpSize`), only the untyped
  overload for `hipFuncGetAttributes`.
- `compare.cxx`: HIP runtime if `__HIP__` or `MFRONT_GPB_COMPARE_USE_HIP`.
- Existing `hip-hipcc` configuration: `nvcc` flags moved into
  `TFEL_HIPCC_FLAGS` (nvidia platform only; `-fPIC` on amd);
  `-isystem @TFEL_HIP_INCLUDE_DIR@` removed (HIP installed in `/usr`:
  `-isystem /usr/include` breaks the `#include_next` of libstdc++, and hipcc
  already knows its directory). The same defect remains in the
  `hip-clang-nvidia-platform` configuration (not enabled here).
- Running the installed tests: copy
  `install/tfel/share/tfel/test-suite/{behaviours,material-properties}/generic-parallel`
  into `scratch/tc`, then `tfel-check` in the subdirectory, with
  `PATH=install/tfel/bin:...` and, for HIP,
  `HIP_PLATFORM=nvidia HIP_COMPILER=nvcc CUDA_PATH=<SDK>/cuda/12.9`.

### 2026-10-05, `stlpar` on GPU with nvc++ 25.11 (C2, in progress)

- `nvc++ -stdpar=gpu` compiles host and GPU in one pass, without `__CUDA_ARCH__`:
  `TFEL_DEVICE_CODE` is now also defined if `__NVCOMPILER_STDPAR_GPU`
  (any function may be offloaded, so the "device" code is chosen
  in the whole translation unit; it is valid on the host).
- **nvc++ 25.11 bug no. 1 (proven, reproducers in
  `results/nvhpc-bugs/`)**: in GPU code at `-O2`/`-O3`, a loop
  `for (int k = 0; k != N; ++k)` (**signed** counter, `!=`, constant
  bound) does N+1 iterations (`int`, `long`, `short`; `unsigned`,
  `size_t`, `unsigned short` and `<` correct; `-O1` correct; also in
  CUDA through `nvc++ -cuda`; nvcc 13 correct). The PTX already contains the
  wrong constant (22.2 instead of 21.6): bug of the optimiser of nvc++,
  not of ptxas. Variant: entirely empty kernel (`empty-kernel-stdpar.cxx`).
  TFEL: only the loops of `BoundsCheck` added by the prototype (stensor overloads) and
  `fses::fill` had this pattern → `unsigned short` counters (consistent with
  the rest of TFEL). The generators emit unsigned counters. The
  user code of the `.mfront` files remains exposed (to be documented).
- **nvc++ 25.11 bug no. 2**: in `unified` memory mode (default on an
  HMM machine like this one), a `STG` (global) write to a stack
  address (`R212 = cvta.local(SP)` in the SASS) in the constructor of the
  behaviour → `cudaErrorInvalidAddressSpace`, whatever `-O`. Not reduced to a
  minimal reproducer (self-referencing objects, non inlined constructor:
  correct in isolation). **Clean workaround: `-gpu=mem:managed`**, which
  is anyway the portable setting of stdpar (does not depend on the
  presence of HMM; on an HMM machine the GPU still accesses all the
  memory).
- Stack: Hosford needs 6 to 8 KB with nvc++ (2.2 KB with nvcc);
  `NV_ACC_CUDA_STACKSIZE` (or `NVCOMPILER_ACC_CUDA_STACKSIZE`) is taken into
  account by stdpar.
- Managed mode, `-O2`: 5/6 behaviours OK straight away, Hosford OK with
  `NV_ACC_CUDA_STACKSIZE=8192`. `-O1`: Hosford still failing at 8 KB (to be
  looked at).

### 2026-10-05, Kokkos and SYCL backends, factorisation, CMake fixes

- Factorisation: `PointIntegrator` (trivially copyable function object
  that integrates the behaviour at one point and returns the statistics),
  `CombineStatistics` and `dispatch` (choice of the direct path or
  gather/scatter depending on the storage of the fields) in
  `Integrate.hxx`. `stlpar`, Kokkos and SYCL use them; no more duplicated
  lambdas, no extended lambdas required (nvcc, Kokkos). CUDA/HIP
  keep their kernel (reductions per warp).
- Kokkos (`Kokkos.hxx`): reduction by functor (`value_type`, `init`,
  `join`), no `KOKKOS_LAMBDA`; on CUDA/HIP, probe kernel
  (`stackSizeProbe`, never launched) to size the stack, since Kokkos does not
  do it. The calling code initialises Kokkos; `queue` can point
  to an instance of the execution space.
- `kokkos` configuration of TFEL: compiler and options extracted from the
  imported targets of Kokkos (`find_package(Kokkos COMPONENTS
  separable_compilation)`, otherwise Kokkos redirects the whole compilation of
  TFEL to its compiler); `-ccbin` = compiler of TFEL;
  `-cudart shared` (Kokkos uses the shared CUDA runtime: a single
  runtime in the process). Pitfalls of `nvcc_wrapper`: `-I dir` with a
  space taken for a source file (use `-Idir`); double `-std`
  → warning on the standard output that corrupts the `.d` files
  (tfel-config already adds `-std=c++20`).
- CMake fix of TFEL: the substitution of a compiler
  (`@TFEL_NVHPC_COMPILER@`, `@TFEL_HIPCC_COMPILER@`, ...) was only
  recorded if the compiler was searched by CMake (missing if it
  is given by the user or already in the cache); the substitutions are in the
  cache: changing a value requires `cmake -U "tfel-substitutions*" .`.
- `compare.cxx`: Kokkos execution layers (`ScopeGuard`,
  `kokkos_malloc<SharedSpace>`) and SYCL (global queue, shared USM, queue
  passed through `BatchData::queue`).
- `scripts/env.sh`: full environment (TFEL, CUDA 12.9, LLVM 20,
  AdaptiveCpp, nvc++, Kokkos) to reproduce all the tests.
- SYCL / AdaptiveCpp: the generic compiler (SSCP) is single pass →
  `TFEL_DEVICE_CODE` defined if `__ACPP_ENABLE_LLVM_SSCP_TARGET__`, as
  for nvc++. On GPU: `ptxas fatal: Unresolved extern variable
  _ZTVN10__cxxabiv117__class_type_infoE`: the behaviour is polymorphic (pure
  virtual interface of `MechanicalBehaviour`), its vtable references the RTTI;
  AdaptiveCpp requires the RTTI, SYCL does not allow virtual functions in
  device code. **Decision**: `MechanicalBehaviour` keeps the types and
  constants but loses its virtual functions (nobody uses the
  behaviours through this base in TFEL; the generated behaviours are `final`, even the
  "behaviour variables" whose wrappers only add friendship);
  the interface is described by the C++20 concept
  `tfel::material::MechanicalBehaviourConcept`, checked by `static_assert`
  in the `generic` and `generic-parallel` drivers; the generators
  no longer emit `final`/`override` on these methods nor on the
  destructor. On CPU (OpenMP): crash in the JIT of AdaptiveCpp
  (`SubCfgFormationPass`), to be looked at again after these changes.

### 2026-10-05, full matrix of the backends, fixes

- A single test matrix: `scripts/run-generic-parallel-tests.sh`
  (copy of the installed tests into `scratch/tc`, `tfel-check` per
  configuration, environment `scripts/env.sh`).
- Pitfall: MFront reads the environment variable `HIP_COMPILER` to
  choose the HIP compiler; `HIP_COMPILER=nvcc` (meant for hipcc) made it
  call the nvcc 13 of the PATH. With `HIP_PLATFORM=nvidia`, hipcc
  already uses nvcc: `HIP_COMPILER` removed from `env.sh` and `build-tfel.sh`.
- CMake of TFEL: `HIP_PLATFORM` put in the cache (an automatic
  reconfiguration without the environment variable failed).
- Stack of the GPU threads: the CUDA driver enlarges the stack to the static need
  of the kernel, but this need ignores the internal routines of the compiler
  (`__cuda_sm20_div_rn_f64_full`, double division) → stack overflow
  (seen with compute-sanitizer for Hosford in SYCL/AdaptiveCpp). New
  shared header `MFront/GenericParallel/DeviceStackSize.hxx`:
  `setMinimalCUDAStackSize` acts through the CUDA runtime already loaded in the
  process (`dlopen(..., RTLD_NOLOAD)`, no CUDA dependency at
  compile time, no effect if CUDA is not used). Option
  `device_stack_size` of the `stlpar` and `sycl` backends, for the behaviours AND the
  material properties (the configuration files are shared).
- Contract violation reports: disabled in all device code
  (`TFEL_DEVICE_CODE`), not only CUDA/HIP or after the inclusion of
  SYCL (AdaptiveCpp on CPU: `reportContractViolation` not defined).
- **Recursion in TFEL**: `StensorComputeEigenVectors<3>::computeEigenVector`
  (eigenvalue assumed simple) fell back, in the degenerate case, on the
  full analytical decomposition which calls `computeEigenVector` again:
  cycle in the call graph (detected by an analysis of the LLVM IR of the
  kernel), which makes a pass of AdaptiveCpp on CPU loop forever
  (`fillTransitiveSplitterCallers`). Fallback replaced by the Jacobi
  algorithm (`fses::syevj3`), without recursion. Only the degenerate case changes
  (non unique eigenvectors).
- Generators: the last two `override` (`IsotropicMisesPlasticFlow`,
  `IsotropicStrainHardeningMisesCreep`) removed; the compilation of the tests
  of TFEL found them.
- Local commits: `TFEL: behaviours without virtual functions`,
  `TFEL: device code selection for single pass compilers`; the backends
  will be committed after the test suite of TFEL.
- **Defects found in the `generic` code (to be fixed in the
  B5 redesign)**:
  1. `logarithmic_strain::integrate`, operator `DTAU_DDF` at the end of the step:
     uses `lgh0` (beginning of the step) with `T1`, whereas the
     `DSIG_DF` branch uses `lgh1`;
  2. `executePostProcessing`: saves/restores
     `d.s0.internal_state_variables` without ever putting the end of step
     values in it (the comment says the opposite): the post-processings see
     the initial internal variables with the final stresses;
  3. `executePostProcessing` without initial state: `*(nullptr)` passed by
     reference (undefined behaviour).
  4. (fixed in commit 494929641, redesign of the `generic` drivers into
     cores callable on GPU): status `if (r)` of the finite strain
     drivers, GL/Hencky initialize functions with swapped
     arguments (compilation error), decoding of the prediction that
     ignored the offset of the speed of sound, pointers not restored if
     an initialisation function throws an exception.

### 2026-10-05 (evening), B5 finite strains, B2 initialisations and post-processings

- **B5** (not committed at the time of writing): integration strategies
  `FiniteStrainIntegrator`, `GreenLagrangeStrainIntegrator`,
  `LogarithmicStrainIntegrator` (`FiniteStrainIntegrators.hxx`), sizes
  per point depending on the options (forces: T for PK1, S otherwise;
  operator: S·S for DS_DEGL, T·T for DPK1_DF, S·T otherwise), generated
  generic traits (`writeGenericBehaviourTraits`), finite strain tools
  annotated for the device (`FiniteStrainBehaviourTangentOperator`,
  `LogarithmicStrainHandler` header only, LU without exception on device,
  `TinyMatrixInvert` returns a boolean).
- **B2**: `_InitializeFunction_<name>` and `_PostProcessing_<name>` per batch,
  same statuses as the integration (1 on success, whereas
  the `generic` interface returns 0). Common primitive of the backends
  `apply(status, n, queue, rdt0, f, ...)`: integration, initialisation
  functions and post-processings = check + choice of the kernel +
  `apply`. `PointView` (`generic` data of one point, local copies if the
  fields are not stored point by point). Null inputs accepted for
  an initialisation function without input (like MGIS).
- Comparison harness (`compare.cxx`): named options
  `--finite-strain=SM,TO`, `--initialize-function=n`, `--post-processing=n`,
  random states for initialisations and post-processings, check that the
  reference does change something (otherwise an empty test would pass),
  documented relative tolerance 1e-10 (different compilers for the
  two libraries: different FMA contraction, amplified by the local
  Newton solvers; SYCL/OpenMP gave 1.1e-12 on Hencky plane stress). The
  comparisons are listed only once in `compare-all.sh` (optional selection
  pattern), called by the 8 `test.check` files.
- Defects found by these tests and fixed:
  1. storage by component: the local view did not copy the end of step
     state, read by the post-processings and only partially
     rewritten by some initialize functions → `gatherState()`
     (mutation test: removing it makes the new test fail);
  2. **nvc++ 25.11** (`-stdpar=gpu`): three code generation defects,
     reproducers in `results/nvhpc-bugs/`:
     - visit of `std::variant` (tangent operator in finite
       strains): "Unexpected branch type" then invalid LLVM aliases
       for the internal union of libstdc++ → `std::variant` excluded from
       device code: new class `tfel::utilities::TrivialGenTypeBase`
       (recursive union, same interface as `GenTypeBase`, trivially
       copyable, unit test `TrivialGenTypeBaseTest`), used by
       `FiniteStrainBehaviourTangentOperator`; default constructor of
       the empty union (otherwise nvc++ still emits an alias);
     - default initialiser `= {}` of an array member:
       "parse integer constant must have integer type"
       (`aggregate-zero-init.cxx`) → member initialised in the constructor;
     - **call through a pointer to member function silently ignored on
       GPU** (`member-function-pointer.cxx`, wrong result at -O0..-O3, without
       diagnostic) → the drivers receive a function object; MFront generates
       a named functor per initialisation function and post-processing;
       the `generic` API with a member pointer stays (host wrapper, `MethodCaller`).
  3. `GenericFixedSizeArray` (base of all the fixed size tensors):
     copy and move written by hand → types not trivially
     copyable; switched to `= default` (same semantics on the logical
     elements), which makes `stensor`, `st2tost2`, etc. trivially
     copyable (memcpy, `variant` without visit, SYCL arguments).
  4. nvcc only warns (#20011) when device code calls a
     host function, and produces a wrong library → `-diag-error 20011`
     in the configurations using nvcc (`cuda-nvcc`, `hip-hipcc` on
     NVIDIA).
- Upstream finding not fixed (semantics): in finite strains, the
  stress modified by an initialisation function is written into a
  temporary buffer and lost (except the Cauchy measure in standard finite
  strains); the FS test behaviours therefore initialise internal
  variables.
- Upstream defects in the generator of the `generic` rotations (B3, to be fixed
  in the batch version): "array" variant of an order 2 tensor
  that is not the first variable written without the offset of the point;
  scalars not copied into `dest`; vectors not supported;
  `_dtau_ddF` uses the arguments of `_dsig_dF` (no effect).
- CUDA sweep to redo after installation (the 33 "UNSUPPORTED"
  were finite strain behaviours).

### 2026-10-05 (evening, continued), B3 rotations of the orthotropic behaviours

- Shared module `MFront/BehaviourRotationFunctions.{hxx,cxx}`:
  description of the rotation functions of an orthotropic behaviour (name, rotated
  objects with type and offset, size per point, direction of rotation,
  possible reason for not being supported) and writing of the rotation code
  of one point, usable on GPU:
  `map<T>(dest + o) = change_basis(T{src + o}, r);`, scalars copied.
- `generic` interface: its ~900 lines of generation of the rotations
  replaced by ~60 lines based on this module (functions and names
  unchanged). Upstream defects fixed along the way: offset of the point forgotten
  for the order 2 tensors that are not the first variable ("array"
  variant), scalars not copied into `dest`, vectors and
  vector/vector blocks not supported (`std::exit(-1)` at run time!):
  the orthotropic thermal integrators of mfem-mgis depend on them.
  The blocks still not supported (arrays, vector/tensor) throw
  an exception (`generic` interface) or return the status -2 (batches).
- TFEL/Math: `change_basis` for `tvector<N>` and `tmatrix<N, N>`, same
  convention as the tensors (test `ChangeBasisTest`: consistency with
  `change_basis` of the symmetric tensors in 2D and 3D).
- `generic-parallel` interface: `mfront_gpb_RotationData` (destination,
  source, matrices per point or uniform if the point stride is zero,
  queue), functions `<behaviour>_<H>_rotate...` per batch, rotation functor
  of one point generated per hypothesis, `RotationExecutor`, `rotate` in the 4
  backends; SYCL: handling of the copies extracted into `DeviceCopies`,
  reused by `DeviceBatchData` and by the rotations.
- Tests: behaviour `OrthotropicStationaryHeatTransfer` (vectors, conductivity
  matrix), copies of `OrthotropicElastic` and
  `OrthotropicSaintVenantKirchhoffElasticity`; harness `--rotation=n`
  (random matrices per point or uniform, two storages, reference =
  `generic` function in place); 15 rotation functions and 5 integrations
  of orthotropic behaviours added to `compare-all.sh`.
- **Defect found by these tests**: with `@RequireStiffnessTensor` (or the
  thermal expansion), the `generic` code reads the elastic properties at the
  head of the material properties (MGIS adds them), but the per-point sizes
  of the batch excluded them → reading of uninitialised values in storage
  by component. Fixed (total size read by the behaviour).
- stlpar/g++: 208/208 comparisons.
- Matrix: `MAKEFLAGS=-j6` in the script (MFront calls `make`), option
  `GPB_OPTIMIZATION_LEVEL` (-O1, -O3: configurations copied into each
  test directory, where MFront looks for them first) and
  `GPB_TESTS_DIRECTORY`.

### 2026-10-05 (night), B4 single crystals, defects found by nvcc

- **Generated methods not annotated**: `execute<Name>InitializeFunction` and
  `execute<Name>PostProcessing` were not `TFEL_HOST_DEVICE` (invisible
  in stlpar/g++, refused by nvcc thanks to `-diag-error 20011,20014`) →
  annotated in `GenericBehaviourInterface`. The matrix launched before this
  fix (CUDA, HIP, Kokkos failing in 24 s) was stopped.
- **Generated slip systems** (`@CrystalStructure`,
  `@SlidingSystem`): singleton with a dynamically initialised local `static`,
  forbidden on GPU. `constexpr` constructor (all the values are
  literals), `TFEL_HOST_DEVICE` methods, `areCoplanar` without
  `std::min/max`, `offset` error without exception on device, access through
  `tfel::utilities::getConstantObject<T>()` (new header
  `TFEL/Utilities/ConstantObject.hxx`):
  - nvcc fills with **zeros, without diagnostic**, the device copy of a local
    `static constexpr` object whose constructor is called (reproducer
    `results/nvhpc-bugs/nvcc-static-constexpr-object.cu`) → namespace
    `__device__ constexpr` variable in CUDA/HIP;
  - SYCL (AdaptiveCpp, GPU and CPU): `inline constexpr` variable read
    correctly;
  - nvc++ `-stdpar=gpu` explicitly refuses global variables and
    local statics in GPU code → copy by value (to be bound to a
    constant reference, as the generated code already does).
- Historical classes of TFEL: `MetallicCFCSlidingSystems` (same
  treatment) and `MetallicCFCGenericSlidingSystemsInteractionMatrix`
  (matrix returned by value, computed at each call: removes the cache
  of the first call, a documented pitfall that ignored the coefficients of the
  following calls; impossible on GPU).
- TFEL/Math: expressions `StensorProductExpr{1,2,3}D`,
  `BuildFromRotationMatrix`, matrix access of the tensors, scalar product
  of vectors, derivatives of fixed size arrays: annotated
  `TFEL_HOST_DEVICE`.
- nvcc: warning #20014 (same defect as #20011, without name) also
  turned into an error.
- Harness: `MFRONT_GPB_COMPARE_NUMBER_OF_POINTS` (debugging),
  `--initial-internal-state-variables=v` or full list, failure if all
  the points fail at the first step (empty test).
- Tests added: `SingleCrystal_DD_FCC` (generated systems, initial dislocation
  density) and `FiniteStrainSingleCrystal` (TFEL classes,
  initial `Fe` = identity). CUDA: gaps 1e-15 and 5e-12; the single crystal
  integrates 100,003 points in 65 ms on GPU against 640 ms for the sequential
  reference.
- stlpar/g++: 220/220 comparisons. Reinstallation, regression and matrix
  launched again.

### 2026-10-05 (night, continued), B1 numerical coverage, new nvcc and clang defects

- Harness: stored and dissipated energies compared at each integration
  (`--energies=true` requires them to be non zero), speed of sound
  (`--speed-of-sound=<density>`: offset +100 of `K[0]`, uniform
  density), type of tangent operator (`--tangent-operator`),
  random perturbation of the initial internal variables and bounds
  policy in both libraries (`--out-of-bounds-policy=<behaviour>,Strict`)
  for batches with mixed statuses.
- Behaviours added: `ImplicitNortonWithEnergies` (`@InternalEnergy`,
  `@DissipatedEnergy`, `@SpeedOfSound`, standard bounds), `Norton_rk42`
  (Runge-Kutta), `Tvergaard` (cohesive zone). stlpar: everything matches; with
  `Strict`, 97 % failures against 50 % (physical bounds only), identical
  statuses.
- **clang CUDA** (`cuda-clang` configuration) compiled no behaviour:
  `TFEL_HOST_DEVICE [[nodiscard]]` is refused (a C++11 attribute must
  come before `__host__ __device__`) → 285 fixes in 41 headers and
  in the generators of MFront; definition of `fses::is_negligible` without
  the annotation of its declaration (refused by clang).
- **nvcc 13.0, NVVM optimiser**: the stress written by the initialisation
  function of SVK (finite strains, Cauchy measure) was
  lost in storage by component (CUDA and HIP, which uses nvcc).
  Correct with clang CUDA (same code, same GPU), with `-G` and with
  `-Xcicc -O0`; the result changes if `printf` calls are added.
  Reduced (not minimal) reproducer in `results/nvhpc-bugs/`. Workaround
  chosen: the end of step state is copied in the constructor of
  `PointView` (like the other inputs) instead of a method
  `gatherState()` called after the construction; the outputs of the view are
  read back through the pointers of the data of the behaviour (`PointView` now owns
  `rdt`, accessor `getTimeStepScalingFactor()`).

### 2026-10-05 (night, continued), B6 standalone models

- New `generic-parallel` interface for models
  (`GenericParallelModelInterface`, registered next to `generic`):
  same backends and same options as the behaviours, library
  `GenericParallel<library|material|Model>-<backend>`, header
  `MFront/GenericParallelModel/<name>-generic-parallel-<backend>.hxx`, one
  function per hypothesis (same name as the `generic` interface), symbols of
  the `generic` interface, `<name>_setParameter`.
- A model is treated as a behaviour without gradients (outputs = internal
  variables, inputs = external variables, constant material properties =
  material properties): `mfront::gpb::ModelIntegrator` (integration
  strategy, like `StandardIntegrator`) builds the model from
  the data of the point and the parameters, then calls `execute`. No
  change in the kernels of the backends.
- Generator of the model class factored into `GenericModelInterface`
  (`writeGenericModelClass`, options `device` and `parameters_type`):
  host output identical to the old one on 7 models of the tests of TFEL; in
  device mode, `TFEL_HOST_DEVICE` constructor that reads the parameters from
  a structure passed by value, functions `bool execute_<f>` (bounds:
  `return false;`), `execute()` that chains them, material laws of
  the `mfront::device` namespace (`writeMaterialLaws`).
- Upstream defect fixed: a model with material properties and
  parameters treated as static variables
  (`parameters_as_static_variables`) did not compile with the `generic`
  interface (the static parameters were initialised in the
  constructor); checked with the installed mfront.
- Semantics: failure → `rdt = 0.1` like `generic`; host exceptions
  (physical bounds, material laws) caught point by point
  (`callAndReportExceptions`, otherwise `std::terminate` in the parallel
  algorithm); no prediction operator: in prediction, the model
  is not executed and the state is not updated (like a behaviour). The
  harness skips the prediction when the tangent operator is empty and
  no longer writes into an empty `K`.
- Tests: `GenericParallelModelTest` (parameter, material property, material
  law, bounds, physical bound for the failures, two chained
  functions) and `GenericParallelModelWithQuantitiesTest` (`@UseQt`),
  added to the 8 configurations and to `compare-all.sh`. stlpar and CUDA:
  zero gaps, identical statuses (50,172 failures out of 100,000 with
  perturbation). Models of TFEL: 6/7 compile for CUDA; the 7th
  (`SiC_IrradiationSwellingModel_GoFaster`) throws an exception in its
  code (not supported on GPU, documented).
- Second upstream defect fixed: a static variable of a type with
  units (`@StaticVariable temperature Tref = 300;`) of a `@UseQt` model
  did not compile (`static constexpr temperature Tref = 300;`); it is
  now written `static constexpr auto Tref = temperature{300};`
  (tested by `GenericParallelModelWithQuantitiesTest`, correct value on
  GPU).

### 2026-10-05 (night, continued), B7 TDLS solver on GPU

- `@LinearSystemSolver TDLS` (TDLS library, built in
  `~/Documents/GitHub/tdls/build-install`, TFEL configured with
  `enable-tdls`) works as is with `generic-parallel`: the
  library is callable on GPU and the generated code only passes
  TFEL types. `ImplicitNortonTDLS` (implicit Norton, consistent tangent
  operator through the partial inverse of the jacobian, hence `factorize`,
  `substitute_inplace` and `solve_inplace` of TDLS): CUDA, gaps 1e-14 to
  2e-13 with the `generic` interface; stlpar g++ OK.
- Tests: `tfel-check` component `mfront::linear_system_solver::tdls`
  declared by CMake if `enable-tdls`; subdirectory `tdls/` in each
  of the 8 configurations (`@Requires` configuration + TDLS), script
  `compare-tdls.sh`. The `run` function and the common variables are
  moved into `compare-functions.sh`, read by `compare-all.sh` and
  `compare-tdls.sh`.
- `scripts/run-generic-parallel-tests.sh`: optimisation level applied
  to the subtests too; result FAILED if one of the tests fails (before, a
  single SUCCESS was enough); number of passed tests printed. New
  version kept in `scratch/run-generic-parallel-tests.sh.new` while
  the running matrix finishes (bash reads the script as it goes).
- Open point (TDLS, outside GPU): g++ 13 and 14 emit a `-Wattributes`
  warning ("declared with greater visibility than the type of its
  field") on lambdas of `tfel::math` when compiling
  `ImplicitNortonTDLS` with the `generic` interface, even without the options
  of `tfel-config`; absent without TDLS and in the TDLS tests of TFEL
  (`StandardElasticity` brick). Harmless; cause not isolated (lead:
  class-type template parameters of TDLS and lambdas in
  initialisers of static members).

### 2026-10-05 (night, continued), behaviour variables and full sweep

- Finding: the compilation sweep (`scripts/compile-sweep.sh`) only
  covered the 226 files at the root of `mfront/tests/behaviours`.
  Left out: `StandardElastoViscoPlasticity` (140 files, the
  brick of the real plasticity behaviours), `StandardElasticity` (21),
  `IsotropicDSLs` (10), `AuxiliaryModels` (7), `DDIF2` (6),
  `BehaviourVariables` (4), `FiniteStrainSingleCrystal` (2), `tdls` (12),
  `unit-tests` (117, of which 96 must fail) and the models (7). The script
  now goes through all of `behaviours` and `models`; a failure of MFront is
  classified `INVALID` if the `generic` interface fails too (unit tests
  that must fail).
- Behaviour variables (`@BehaviourVariable`,
  `@BehaviourVariableFactory`, `@Model` in a behaviour): host calls from
  the GPU (`ExternalModelCaller2/4` in the old sweep). Fixed in
  `BehaviourCodeGeneratorBase`: constructor and `make()` of the factory,
  default constructors of the internal behaviour and of its classes
  `BehaviourData`/`IntegrationData` made `TFEL_HOST_DEVICE`. The parameters
  of an internal behaviour (non static) were read from its host singleton
  (`ParametersInitializer::get()`, not exposed by the library of the main
  behaviour): in device code, they take their default values
  (method `writeBehaviourParameterInitialisationToDefaultValues`, extracted
  from the existing code of the parameters treated as static variables).
  CUDA and stlpar: `BehaviourVariable1`, `BehaviourVariableFactory1/2`,
  `ExternalModelCaller2/4` identical to `generic`, and a variant with
  non static parameters.
- Tests added to the 8 configurations: `BehaviourVariableTest` (factory,
  internal behaviour `IsotropicPlasticMisesFlow` with units and non static
  parameters) and `ExternalModelCallerTest` (Runge-Kutta rk54 + external
  model with a parameter).

### 2026-10-05 (night, continued), StandardElastoViscoPlasticity brick and headers of TFEL

- CUDA sweep (development mfront) of the 140 behaviours of
  `StandardElastoViscoPlasticity`: 118 OK on the first pass, 140/140
  after fixes:
  - GTN and Rousselier criteria, elastic prediction of the axial
    strain (logarithmic finite strains): `throw` in case of
    non-convergence → in device code, NaN (the integration fails), like the
    material properties;
  - staggered scheme of the porosity (`tfel::raise<DivergenceException>` in
    `@AdditionalConvergenceChecks`, which cannot return a failure): in
    device code, the iteration counter of the solver is set to `iterMax`
    to stop the resolution on a failure. `TinyNonLinearSolverBase`
    now compares with `>=`/`<` instead of `==`/`!=` (identical in
    normal operation, well defined if the counter is forced);
  - generated functions `mfront_computeYieldRadius*` (isotropic hardening
    defined by data) without `TFEL_HOST_DEVICE`; **upstream defect**: the
    name of `mfront_computeYieldRadiusAndDerivative` was doubled
    (`m2 + local_id`) in the spline variant;
  - `AitkenAccelerationAlgorithm`, linear and cubic spline
    interpolations: annotations;
  - **upstream defect**: 2-argument declarations of the Chu-Needleman
    nucleation functions (strain) never defined, the 4-argument definitions
    used by the generated code not being declared (hence
    host only);
  - **nvcc 13.0 defect** (cudafe++): `auto i = decltype(s){};` with
    `const auto s = c.size();` is rewritten for the host compiler as
    `const unsigned long{}`, invalid → `size_type{}` in
    `LinearInterpolation.ixx` (more readable anyway).
- **clang CUDA sweep of the headers** (`scratch/hdscan`): each header
  of `TFEL/Math` and `TFEL/Material` (335) compiled alone by clang in CUDA
  mode. clang refuses a definition without `__host__ __device__` whose
  declaration has them (nvcc accepts it): 60 + 25 definitions fixed
  (criteria of Barlat, Cazacu, Drucker, Mohr-Coulomb, Michel-Suquet,
  Rousselier, GTN, pi plane, nucleations, `DDIF2Base`, solvers of Broyden,
  Levenberg-Marquardt, Powell, bisection, scalar Newton, PCW
  homogenisation schemes) by a script (`fix_hd.py`) then checked;
  attribute `[[gnu::visibility]]` placed after `__host__ __device__` in
  `DDIF2Base`; 19 headers using `TFEL_HOST_DEVICE` without including
  `TFEL/Config/TFELConfig.hxx`. No defect of this kind left.
- **False GPU defect explained**: `OrthotropicNortonTest` (Hill) gave
  a 66 % gap on GPU (nvcc and clang) at 3e-3 strain, but exact
  on CPU, including with `-DTFEL_DEVICE_CODE` and
  `-ftrivial-auto-var-init=pattern`. The same gap appears on CPU with
  `-march=native -ffp-contract=fast` (FMA): ill-conditioned problem
  (`@Epsilon 1e-16`, a few points failing), Newton converges elsewhere
  depending on the roundings. It is not a defect of the GPU. Criterion chosen for
  the test behaviours: the `generic` interface and `generic-parallel` compiled
  with FMA on CPU must agree before any GPU comparison.
- Tests added: `StandardElastoViscoPlasticityTest` (plastic, Mises,
  isotropic hardening by data, Armstrong-Frederick + Chaboche 2012),
  `OrthotropicStandardElastoViscoPlasticityTest` (Norton, Hill, Voce +
  linear), `PorousStandardElastoViscoPlasticityTest` (GTN, Chu-Needleman
  nucleation, staggered scheme, initial porosity 0.005): robust (FMA)
  and identical within 1e-15 on CUDA.

### 2026-10-05 (night, continued), commits 0016 to 0018

- TFEL regression: 6965/6965 with all of the above (fresh
  installation). Local commits: 7b1ed663c (TFEL: consistent annotations,
  criteria and solvers), 4aa867593 (MFront: behaviour variables,
  models, porous flows), bd5528869 (generic-parallel: models and
  tests). Bundle `backups/tfel-generic-parallel-behaviours-2026-10-05.bundle`
  checked, `proto-patches/` regenerated (0001 to 0018).
- O2 matrix on the fresh installation in progress
  (`results/generic-parallel-tests-2026-10-05-c.txt`).
- Script `scripts/clang-cuda-header-scan.sh`: check of the host/device
  annotations of all the headers by clang CUDA (0 inconsistency).

### 2026-10-05 (evening), C5 OpenMP backend, factored entry points

- Common entry points in `Integrate.hxx`: `integrateBatch`,
  `executeInitializeFunctionOnBatch`, `executePostProcessingOnBatch`,
  `rotateBatch` (check of the data, choice of the function object, call of
  the `apply` primitive of the backend passed as a callable object). The
  stlpar, CUDA/HIP and Kokkos backends now only bind their primitive
  to their options (`internals::makeApply`); SYCL keeps its entry points
  (copies to the device before the choice). stlpar and CUDA: 46/46
  identical comparisons after the factorisation.
- `openmp` backend (`OpenMP.hxx`, `OpenMPBackend`): `target teams
  distribute parallel for map(to: f)` with member by member reductions of the
  statistics (`+` for the numbers of failures and of unreliable results,
  `min` for `rdt`) or `parallel for` on the host (option `offload`).
  `TFEL_DEVICE_CODE` defined for `nvc++ -mp=gpu` (`__NVCOMPILER_OPENMP_GPU`,
  single pass compiler like `-stdpar=gpu`), `getConstantObject` returns
  a copy as with `-stdpar=gpu`.
- **nvc++ 25.11 defect (OpenMP, -O2 and -O3)**: the behaviours of the implicit DSLs
  fail as soon as a tangent operator is requested (even elastic); the
  Newton solver converges on GPU, but `computeConsistentTangentOperator` returns
  false (LU decomposition of the jacobian). Correct at -O0 and -O1; a `printf`
  in the failure branch is enough to fix the result (defect of
  the optimiser). Not reproduced out of context (7×7 LU alone in a
  `target` region: correct), nor fixed by `-Mnounroll`, `-Mnovect`, `-Mnoinline`,
  `-Mnofma`, nor by `map(to:)` instead of `firstprivate`. The same code is
  correct with `nvc++ -stdpar=gpu` at -O2. Documented workaround: the
  `openmp-offload-nvhpc-gpu` configuration is compiled at -O1. clang cannot
  serve as a reference compiler here: the LLVM 20 of spack
  is built without OpenMP offload support (`offload=False`).

### 2026-10-05 (evening, continued), material properties: Kokkos and OpenMP backends

- The `generic-parallel` interface of the material properties
  only had CUDA, HIP, stlpar and SYCL. Added `kokkos` (option
  `execution_space`, `Kokkos::parallel_for`, statuses of the bounds in a
  `Kokkos::View` in `SharedSpace`) and `openmp` (options `offload`,
  `device_stack_size`, same constructs as for the behaviours).
- Factorisation in `CxxProgrammingModelBackendBase`: choice of the kernel
  depending on the strides (`writeKernelSelection`, previously copied in stlpar and
  SYCL), writing of the kernels (`writeKernels`), specifiers of the lambdas
  (`getKernelSpecifiers`, `KOKKOS_FUNCTION` for Kokkos). nvcc forbids
  initialisers in the captures of extended lambdas (`__host__
  __device__`): the pointer to the statuses of the bounds and the parameters
  are now declared before the kernels and captured by name (container
  `mfront_bounds_statuses_storage` in each backend).
- Common test `test-cxx.cxx` (formerly `test-stlpar.cxx`) for stlpar,
  openmp and kokkos: backend chosen by `MFRONT_GPMP_TEST_BACKEND`, arrays
  in Kokkos shared memory if `MFRONT_GPMP_TEST_USE_KOKKOS`. Results:
  stlpar g++, OpenMP host, OpenMP offload (nvc++ -O2, 19 GPU kernels),
  Kokkos CUDA and Kokkos host: all the tests pass. Test directories
  `kokkos`, `kokkos-host`, `openmp-host-gcc`, `openmp-offload-nvhpc-gpu`
  added; no configuration reserved to the behaviours any more.

### 2026-10-05 (evening, continued), matrix c, OpenMP, models compiled by nvc++

- O2 matrix (`results/generic-parallel-tests-2026-10-05-c.txt`, state of
  commits 0016-0018): CUDA, HIP, Kokkos, Kokkos host, stlpar g++ pass
  with the TDLS subtests; `stlpar-parunseq-nvhpc-gpu`: 280 correct
  comparisons and 4 wrong ones, all on the failure case of the models
  (`GenericParallelModelTest_PlaneStress` with perturbation). Cause: the
  reference library (`generic`) is also compiled by
  `nvc++ -stdpar=gpu`, where `TFEL_DEVICE_CODE` is defined in the whole
  compilation unit: the bounds checks then return a boolean instead
  of throwing an exception, and the host functions of the models (`void`)
  ignored it (physical bounds silently ignored). Fixed:
  `writeGenericModelClass` always writes `bool` functions that
  return false if a check fails and an `execute` method, which
  the host implementation uses (failure reported, `rdt = 0.1`). On a normal
  host, nothing changes (the exceptions are thrown). 12/12 with nvc++.
- OpenMP host (g++ `-fopenmp`): `offload` becomes a template parameter
  (`if constexpr`): the host code no longer depends on the offload
  entry points of the OpenMP runtime (`GOMP_teams4` missing from the
  `libgomp.so.1` libraries of NVHPC and LLVM present in
  `LD_LIBRARY_PATH`). Full suite: 284/284.
- OpenMP offload (nvc++ -O1): each behaviour is correct on its own, but the
  library of the whole suite does not load on GPU ("Could not find
  symbol _st_... on the device"). Bisection by relinking the
  existing objects: `HosfordPerfectPlasticity` alone at fault; its device
  code (`.gpu` IR kept by `-gpu=keep`) defines an internal static block
  of 317 bytes never used, removed by the optimiser, that the
  runtime looks for at load time. Defect of nvc++ 25.11; with the one at
  -O2, the `openmp-offload-nvhpc-gpu` configuration remains provided but is
  not enabled (documented in the page of the interface and in
  `results/nvhpc-bugs/README.md`).
- End of matrix c: SYCL GPU and SYCL CPU (AdaptiveCpp, single pass
  generic compiler) fail on the same case of the models as nvc++
  (same cause, fixed by the `bool` functions); all the material
  properties pass; regression 6965/6965.
- `scripts/run-generic-parallel-tests.sh`: tfel-check reports as passed
  a test whose requirements (`@Requires`) are not met (all its
  commands skipped); the script now counts separately the passed
  and skipped tests (configurations not enabled, TDLS missing).
- Validation chain launched (`scratch/validate2.sh`): reconfiguration
  with `openmp-host-gcc`, installation, regression, matrices -O2, -O1, -O3
  (`results/generic-parallel-tests-2026-10-05-d-O*.txt`), CUDA sweep of
  all the tests (`scratch/sweep-config-cuda-nvcc/summary.tsv`).
- Regression after installation (validate2): 6965/6965. Local commits
  4da85b2c5 (TFEL: `TFEL_DEVICE_CODE` with nvc++ -mp=gpu), b0290dce3
  (MFront: bounds checks of the models), b401968ce (generic-parallel:
  OpenMP, Kokkos for the material properties, common entry points);
  backup bundle and `proto-patches/` (0001 to 0021) up to date.
- MGIS study (D1, read only): `scratch/d1-mgis-design.md`. Important
  finding: MGIS today loads without error an isotropic `generic-parallel`
  library (same function names, different
  signatures) and the call is then undefined behaviour; the interface must be
  detected (`_mfront_interface` == "GenericParallel").
  Recommendation: option B (MGIS assembles one field per group, without copy
  in the common cases), memory through `std::pmr::memory_resource`, new
  overload `integrate(BatchExecutionPolicy, MaterialDataManager&, ...)`.

### 2026-10-06 (night), versions of NVHPC and CUDA, defects of OpenMP offload

- Release notes of NVHPC 26.1 to 26.9: no list of fixed defects.
  Reproducers replayed with NVHPC 26.3 and 26.9 installed in
  `deps/nvhpc/` (archives checked against the checksums of the spack package):
  only `aggregate-zero-init` is fixed (26.9); `run-nvcxx-reproducers.sh`
  and the new reproducer `std-variant-visit.cxx` in
  `results/nvhpc-bugs/`. Defects of nvcc present from 13.0.88 to 13.4.92
  (compilers of the pip packages `nvidia-cuda-nvcc`).
- **Hosford under OpenMP offload** (symbol `_st_…__T1_0` not found): the
  317 byte block is the `__PRETTY_FUNCTION__` text of the `assert` calls of
  `StensorComputeEigenVectors<3>::computeEigenVectors`; nvc++ removes
  the `assert` from the device code but keeps the registration of the text. The three
  `assert` calls were redundant (implied by the conditions of the branches),
  one was wrong at equality (`<` instead of `<=`): replaced by
  comments stating the invariants. These are the only `assert` calls of the
  headers of TFEL used in device code.
- **Wrong tangent operator at -O2/-O3 under OpenMP offload**: neither the LU, nor
  the stack (identical frames at -O1 and -O2, ptxas), nor an uninitialised
  variable (`-ftrivial-auto-var-init`). Distinct markers of `rdt` on
  each failure path: `computeConsistentTangentOperator` returns false
  without going through any of its `return false`. Trigger: the local object
  `GetPartialJacobianInvert` generated by MFront, which keeps a `bool&`
  reference to the local flag `mfront_success` (its construction alone
  is enough). Generator fixed (`ImplicitCodeGeneratorBase`): each object
  (`GetPartialJacobianInvert`, `GetIntegrationVariablesDerivatives_*`)
  owns its flag, combined into `mfront_success` after the code of
  the user; the substitution goes through a local variable (a failure
  is no longer erased by a following successful call).
- icpx (oneAPI 2026.1.1, apt) installed: SYCL on CPU only (NVIDIA
  plugin removed since 2025.3); `libintelocl` requires the TBB of oneAPI
  (2023.1): paths `compiler/2026.1/lib`, `umf/1.1/lib`,
  `tbb/2023.1/lib/intel64/gcc4.8` to be added at the end of `LD_LIBRARY_PATH`.
- clang CUDA: clang 20 fails with CUDA 13.0 (`texture_fetch_functions.h`
  removed) and with the CUDA 12.9 of NVHPC (curand headers in
  `math_libs`). clang configurations on GPU postponed to the new
  environment (LLVM 23, CUDA 13); a CUDA prefix grouping `math_libs`
  is needed if the CUDA comes from NVHPC.
- Tests of the behaviours added for `stlpar-parunseq-clang`,
  `stlpar-parunseq-icpx` and `sycl-default-icpx` (to be validated).

### 2026-10-06 (night, continued), fields per variable, symbols, roles

- `BatchData`: the material properties and the external variables can
  be given per group (one field) or per variable (one field per
  variable, optional members `material_properties_by_variable` and
  `external_state_variables_by_variable`). A structure initialised
  as before keeps its meaning. Both forms are converted on the host side
  into `InputFields` (fields and numbers of components), passed by value
  to the kernels: a single code path per point. The material properties
  and the external variables of a point are always copied, the direct
  mode only depends on the other fields now (a uniform property no
  longer disables it). SYCL copies each variable if needed.
  `PointTraits` gives the number of components of each variable (behaviours
  and models), read on the host side only.
- `compare.cxx`: third storage "per variable" (uniform
  values, zero stride) for all the comparisons: 396 instead of 284
  with stdpar g++.
- Symbols: `<behaviour>_generic_parallel_backend` (behaviours and models).
  Defect fixed: the `generic-parallel` models exported
  `mfront_interface = "Generic"`.
- Roles of TFEL/MFront, MGIS and mfem-mgis: section of `DESIGN.md`.
  MGIS design (`scratch/d1-mgis-design.md`, 2.3): MGIS passes one field
  per variable, without workspace.
- Development in the `tfel-dev` working tree (branch
  `batchdata-by-variable`) during the validation of `tfel/`; development
  MFront in `scratch/build-dev2`, tests through
  `scratch/dev-test.sh <configuration>`.

### 2026-10-06 (night, continued), drivers of the generic interface, OpenMP in direct mode

- Development test of OpenMP on GPU (-O2) of the fields per variable: a
  single failure, the initialisation function of SaintVenantKirchhoffElasticity
  in direct mode (fields point by point). The written stress is lost
  (the batch keeps the initial stress). Correct at -O1. Moving a copy
  into `PointView` makes the three storages fail: a defect of
  the optimiser sensitive to the shape of the code, like the one of nvcc 13.0 on the
  same function.
- Trigger: the drivers of the `generic` interface (integration,
  initialize functions, post-processings of the behaviours in small
  strains, in finite strains, in Green-Lagrange strain
  and logarithmic strain) temporarily modified the pointers of the data
  of the behaviour given by the caller, then restored them. They now work
  on a local copy of the data whose pointers designate
  the local values; for the Cauchy stress, the finite strain drivers
  directly call the small strain
  ones. Results unchanged; OpenMP on GPU -O2: 396/396.
- Latent defect fixed along the way: the stress computed by an
  initialisation function was only exported for the Cauchy
  stress of the standard finite strain behaviours. It is now
  converted into the measure of the caller (deformation gradient
  at the beginning of the step) for all the measures and the behaviours written in
  Green-Lagrange or logarithmic strain. Common conversion in
  `MFront/GenericBehaviour/StressMeasures.hxx`. Tests added (PK2, PK1).
- C structures of `BatchData.h` named (`typedef struct name {...} name;`):
  cudafe (nvcc 13.0) crashed on the destructor of a `std::vector` of
  anonymous fields.
- -O2 matrix of the night (commits of the evening of 05/10 and of the OpenMP fixes):
  OpenMP on GPU and stdpar icpx pass; stdpar clang fails because
  of the environment: the system clang loads the library of the LLVM
  of spack (`LD_LIBRARY_PATH`), built without the AArch64 target. The clang
  of this LLVM is now used.

### 2026-10-06 (morning), clang CUDA and HIP configurations

- The behaviours are now tested with `cuda-clang` and `hip-clang`
  (NVIDIA platform). Only the material properties were before.
- clang 20 only knows the CUDA versions released before it (up to 12.8). The
  CMake variable `TFEL_CLANG_CUDA_PATH` chooses the CUDA of clang (by
  default the one of nvcc). Here: CUDA 12.9 of the HPC SDK (merged prefix from
  `deps/spack/make-nvhpc-cuda-prefix.sh`).
- clang refuses any call of a host function from GPU code. nvcc
  accepts `constexpr` host functions with `--expt-relaxed-constexpr`.
  Two holes found this way: `pull_back` and `invert(st2tost2)`, called by
  the conversions of the finite strain behaviours. They are now
  `TFEL_HOST_DEVICE` (`invert` relies on `TinyMatrixInvert`, already
  usable on GPU). Only one `constexpr` host function remains in
  TFEL: `tensor::operator()(i, j)`, which throws an exception.
- System HIP in `/usr`: `-isystem /usr/include` placed this directory
  before the headers of libstdc++ and broke their `#include_next`. The
  clang configurations pass the directory of HIP with `-I`, which the
  compilers ignore for a system directory.
- TDLS subtests added to the five configurations that had none
  (cuda-clang, hip-clang, stlpar clang, stlpar icpx, SYCL icpx).
- Results (development tree, -O2): cuda-clang and hip-clang
  402/402 each, of which 114 with the fields per variable; TDLS 18/18 for
  the five configurations.
- `scripts/compile-sweep.sh` recognises the errors of clang. Sweep of
  all the behaviours, models and properties of MFront with cuda-clang launched.
- Branch `clang-cuda` of `tfel-dev` (3 commits), to be merged after the
  validation in progress.

### 2026-10-06 (morning, continued), clang CUDA sweep of all the tests of MFront

- Sweep of the 608 behaviours, models and properties of `mfront/tests` with
  cuda-clang (development MFront and headers, 21 min): 458 OK on the
  first pass. clang refuses any host call from GPU code, so it
  complements the nvcc sweep.
- Host calls (38): 37 in the code written by the user of the test
  behaviours (`throw`, `std::cout`, `tfel::raise`). A single hole of TFEL: the
  `tensor` view of a `stensor`, used by any product or sum mixing
  `tensor` and `stensor` (for example `F * S`), and the `(i, j)` accessors
  of the tensors were reserved to the host, because they throw an exception
  on an invalid index. New function
  `internals::reportInvalidTensorIndex`: exception on the host, NaN on
  GPU, like `TinyMatrixInvert`. No `constexpr` function
  reserved to the host left in `TFEL/Math`.
- Test behaviour added to all the configurations:
  `DefaultGenericSaintVenantKirchhoffElasticity` (gradient F, force PK1,
  product `tensor` × `stensor`, conversion of the tangent operator).
  nvcc and clang: 18/18.
- Compilation errors (32):
  - 26 unit tests `WarningTest*`, which must not compile. The
    sweep now also compiles the `generic` interface before
    concluding: they are classified `INVALID`;
  - material properties defined by interpolation of data (`@Data`,
    5 files): the CUDA/HIP kernel named its input fields
    `mfront_<input>_values`, the name of the abscissae of the data. Renamed
    `mfront_<input>_field`. Also failed with nvcc;
  - material property calling another law (`@MaterialLaw`): clang
    also analyses the body of a `__global__` kernel during the host pass,
    where only the host versions of the called laws were visible. The
    kernels now always import the `mfront::device`
    implementations.
- Tests added to the 16 configurations of the material properties:
  `LinearDataInterpolationTest` and `MaterialLawCallTest`, checked by a
  common function in each test program. Latent defect fixed:
  `hipMallocManaged` receives a size in bytes (4 bytes allocated for
  4 reals). Six configurations tried (nvcc, clang, hipcc, hip-clang,
  OpenMP, stdpar): all pass.
- "Suspects" (39): all material properties, whose kernels are
  small with clang. The registers criterion now only applies to the
  behaviours and the models.
- Branch `clang-cuda`: 8 commits above the validated branch.

### 2026-10-06 (morning, continued), DSLs not run on GPU and two defects of nvc++ OpenMP

- Run-time coverage: four DSLs of MFront were only compiled
  by the sweep, never run on GPU. Test behaviours added to all
  the configurations: `IsotropicMisesCreepTest`,
  `MultipleIsotropicMisesFlowsTest` (Voce hardening and creep),
  `ImplicitGenericBehaviourTest` and `ImplicitModelTest`. CPU, nvcc, clang,
  hipcc: 60/60; SYCL AdaptiveCpp: 84/84.
- Upstream limitation, unrelated to the GPU: the
  `MultipleIsotropicMisesFlows` DSL does not provide a consistent tangent
  operator, and the driver of the `generic` interface then refuses any
  request for an operator, even elastic. No upstream test runs a behaviour
  of this DSL (they are only compiled). The comparison covers the
  integrations without operator and the prediction.
- nvc++ 25.11 OpenMP on GPU, -O2 and -O3, two defects found by the
  behaviour `DefaultGenericSaintVenantKirchhoffElasticity` (product `F * S`):
  - internal error "Computed goto control expression must be a 32-bit
    integer": the `switch` of the `tensor` view of a `stensor`, inlined into
    the one of the product. Workaround: view written with an index table
    and a comparison, without `switch`;
  - kernel code removed, no result written: zeroing an
    array member (`v()`, `v{}`, `= {}`) is done by copying an
    uninitialised `.common` global array. Standalone reproducer of 30
    lines. Workaround: zeroing by a loop in the
    default constructor of `GenericFixedSizeArray`.
  Reproducers and workarounds in `results/nvhpc-bugs/README.md`.
  Method: the sentinel matters. With `map(from:)`, the entries not
  written return the GPU memory of a previous trial; only
  `map(tofrom:)` with a sentinel value distinguishes an empty kernel.
- After the workarounds, the five new behaviours on nvc++: OpenMP on GPU
  78/78 at -O1, -O2 and -O3, stdpar 78/78. `Math` unit tests of TFEL
  (host): 80/80.
- `scratch/dev-test.sh`: variable `DEVTEST_OPTIMIZATION_LEVEL`, like the
  matrix.
- Branch `clang-cuda`: 10 commits above the validated branch.

### 2026-10-06 (morning, end), validation of the merged branch

- The matrices now run at -O3 only: when -O3 passes, -O1 and -O2
  pass in the vast majority of cases. `VALIDATE_LEVELS` chooses
  other levels if needed.
- Branch `clang-cuda` merged into `generic-parallel-behaviours`
  (9cae5c235). `scratch/validate3.sh 2026-10-06c`:
  - regression of TFEL: 6965/6965;
  - -O3 matrix: the 31 tests pass (15 configurations of behaviours,
    all with their TDLS subtest, 16 of material properties), none
    skipped (`results/generic-parallel-tests-2026-10-06c-O3.txt`);
  - compilation sweeps of the 608 behaviours, models and properties,
    identical with nvcc and clang: 504 OK, 36 user code, 60
    invalid (of which 26 unit tests that must not compile), 7
    fragments, 1 host call (`results/compile-sweep-*-2026-10-06.tsv`).
- The last host call came from MFront: the comparison to the numerical
  jacobian (`@CompareToNumericalJacobian`) prints its results with
  `std::cout`. It is now generated for the host only (53d30ff6c).
- The 36 cases of user code are all in the test behaviours:
  `throw`, `tfel::raise`, `std::cout`, statics initialised at
  run time, slip systems of the polycrystals read from a
  file (documented limitation).

### 2026-10-06 (afternoon), fixes from the review (R1)

Branch `review-fixes` of TFEL, on top of 53d30ff6c:

- A4: the behaviour class has the same definition on the host and the GPU
  (interface property `requiresDeviceCompatibleBehaviourClass`, DSL
  option `device_compatible_class` propagated to the behaviour variables).
- A1: the reference of the tests is compiled by the default compiler
  of MFront (`build-libraries.sh`), independently of the tested
  configuration.
- A2: the parameters of the behaviour variables are compared on
  the host to their default value, used on the GPU (status -2 otherwise).
  The parameters of the called material laws are constants on the host
  and the GPU: this point of the review was wrong.
- B2: a single implementation of the entry points. The generated code
  calls `integrateBatch` (and the three others) with the `makeApply`
  primitive of the backend; SYCL adds its transfers. 677 fewer lines.
- A3: execution options (`mfront_gpb_ExecutionOptions`): queue,
  asynchronous statistics (CUDA, HIP, Kokkos, SYCL without copy),
  `preserve_stack_size`. A single stack mutex per library. The second
  step of the comparisons is asynchronous.
- B4: `MFRONT_GPB_API_VERSION`, exported by the libraries. The
  per-point messages are not transmitted (no exceptions and no buffer
  per thread on GPU): the `generic` behaviour gives the message of a point.
- B1: options common to both interfaces (`BackendOptions`), Kokkos execution
  space validated, status conventions documented.
- B5: `compare` reads the sizes from the reference library and
  deduces the tested function from its name. The sizes read are identical
  to the sizes given by hand for the 75 comparisons.
- A5: eigenvectors of a multiple eigenvalue by Jacobi on GPU
  only (CPU results unchanged), incompatibilities in the release
  notes, effect of the single pass compilers documented, principle 4 of
  DESIGN.md fixed, plan to split into PRs (`UPSTREAM.md`).
- B3 and C: deduplication of the generator documented, status initialised
  for an empty batch, out of bounds policy of the models, `GenTypeBase`
  restored, overload of `callAndReportExceptions`, host/GPU utility
  in the Runge-Kutta DSL, empty CMake list removed.
  `getConstantObject` keeps its return type: copy imposed by nvc++,
  reference elsewhere, API unchanged for the other compilers.

- B7: direct accesses to the fields stored point by point against local
  copies only (`scratch/b7-measure.sh`, point by point layer):

  | | nvcc compilation | gcc compilation | execution |
  |---|---|---|---|
  | SingleCrystal_DD_FCC | 7.9 s against 5.8 s | 2.6 s against 2.6 s | identical |
  | ImplicitNorton | 16.8 s against 10.6 s | 5.3 s against 4.4 s | identical |

  The execution time is the same within 1 % on GPU (RTX 5080); on
  CPU, the noise of OpenMP (16 threads) hides any difference. The direct mode
  is removed: a single variant of the kernels.

Tests at -O3 (`scratch/dev-test.sh`): 489/489 for cuda-nvcc, cuda-clang,
hip-hipcc, kokkos, kokkos-host, openmp-host-gcc,
openmp-offload-nvhpc-gpu, stlpar-parunseq-gcc,
stlpar-parunseq-nvhpc-gpu; 513/513 for sycl-default-acpp (with the
copies from host memory); TDLS 18/18; material properties OK on
five configurations.

### 2026-10-06 (evening), portable test behaviours (R2)

Branch `portable-test-laws` of TFEL, on top of 7c1427f57.

- `tfel::fail<Exception>(...)`: throws the exception on the host as before,
  returns `false` on GPU. The code block returns this value.
- Returning `false` does not make the integration fail everywhere: in
  `@Integrator` of the implicit DSLs, the Newton correction is divided
  by two; in `@Derivative` of the Runge-Kutta DSLs, the time step is
  reduced. An exception, on the other hand, stops everything. On GPU, the same effect
  is obtained by `this->abortResolution()` (non linear solver of TFEL,
  iteration counter at the maximum) and `this->abortIntegration()` (Runge-Kutta
  behaviours, generated by MFront).
- `@SharedData Type name = expression;`: data built once on
  the host (orientations of the grains of a polycrystal read from a file).
  The behaviour accesses it through a reference member. The generic-parallel interface
  copies it at the first call into memory accessible to the GPU (memory of the
  GPU for CUDA, HIP, SYCL and Kokkos, heap otherwise) and passes its address to the
  kernels through the structure of the parameters. `PolyCrystalsSlidingSystems` in
  fixed size arrays. Forbidden in the behaviour variables.
- 36 behaviours of `mfront/tests` converted without changing their results on
  the host (same exceptions, same messages): 30 with exceptions or `cout`,
  6 polycrystals. The messages built at run time and the `cout`
  stay on the host (`#ifndef TFEL_DEVICE_CODE`).
- New comparison test: Berveiller-Zaoui polycrystal with four
  grains (`PolyCrystalSharedDataTest`).
- Environment defects of NVHPC 25.11 found along the way (in
  `results/nvhpc-bugs/README.md`): `-mp=gpu -gpu=mem:managed` does not link
  `libnvhpcmanaux` (unlike `-stdpar=gpu`), hence a
  `free(): invalid pointer` when g++ code frees managed memory;
  `libnvcpumath.so` does not depend on `libm`: a program that does not load
  `libm` cannot load an nvc++ library (no effect on
  MGIS).
- Compilation sweeps of the 608 files with nvcc and clang: 541 OK,
  0 user code, 0 host call (against 504 OK, 36 user code and
  1 host call on the morning of 06/10), 60 invalid and 7 fragments unchanged.
- The driver of the `generic` interface built the behaviour outside the zone
  that converts the exceptions into failures. With `@SharedData`, the reading
  of the orientations file happens in the constructor: a missing
  file would have made the exception leave the generated function instead
  of returning a failure as before. The behaviour is now built in
  the protected zone (checked: status -1 and message of the exception).
- `compare` accepts the time step and the values of the external variables
  (for R3).
- Validation at -O3: 498/498 on cuda-nvcc, cuda-clang, hip-hipcc, kokkos,
  kokkos-host, openmp-host-gcc, openmp-offload-nvhpc-gpu,
  stlpar-parunseq-gcc, stlpar-parunseq-nvhpc-gpu; 522/522 on
  sycl-default-acpp; TDLS 18/18; material properties OK on five
  configurations. After the fix of the `generic` driver: openmp-host-gcc,
  cuda-nvcc and openmp-offload-nvhpc-gpu again 498/498,
  sycl-default-acpp 522/522.

### 2026-10-07 (night), execution trials of all the test behaviours (R3)

The behaviours of `mfront/tests` that compile with the `generic-parallel` interface
are now compared at run time, and no longer only compiled:

- `scripts/runtime-sweep.py` (behaviours and models): for
  each behaviour, `generic` (reference) and `generic-parallel` libraries,
  data read from the MTest file of the behaviour (material properties,
  external variables, initial values of the internal variables,
  parameters, time step, amplitude of the loading; definitions evaluated
  in any order, evolutions included), then `compare`. Without
  MTest, the material properties take synthetic values
  chosen from their name. If the reference fails at all the points,
  smaller time steps and amplitudes are tried;
- `scripts/property-sweep.py` (material properties): `generic` function
  called point by point against the two `generic-parallel` prototypes,
  with a uniform argument, and with out of bounds values and the
  strict policy;
- `scripts/roundoff-check.sh`: confirms that a small gap comes from
  roundings (library of the batch compiled on the host with FMA contraction,
  like nvcc on the GPU, against a reference without FMA).

Configurations: openmp-host-gcc and cuda-nvcc (comparisons one by one
on the GPU).

Defects of TFEL found and fixed:

- arrays of gradients and of thermodynamic forces (`@Gradient Stensor
  a[2]`): the size of a point ignored the size of the arrays, only the
  first component was copied (`getPointDataSizes`);
- conversion of the tangent operator of the finite strain behaviours:
  the generated code reads the operator computed by the user with
  `Dt.template get<T>()`, which checks the type on the host only. If the
  block of the user does not compute the right operator, the host fails
  (GenTypeCastError exception), but the GPU returned a success with a
  wrong matrix. The type is now tested before the reading: same
  exception on the host, failure on the GPU;
- parameters of the models shared between libraries (upstream code of
  the `generic` interface): the singleton of the parameters is the static
  variable of an inline function, which GCC emits as a unique symbol
  (STB_GNU_UNIQUE) and which the loader unifies between libraries, even
  loaded with RTLD_LOCAL. The `generic` and
  `generic-parallel` libraries of the same model shared their parameters. The
  generated code of the models is now in an anonymous namespace;
- physical bounds of the arguments of the material properties (upstream code,
  all the interfaces): not exported when the argument had no
  standard bounds;
- cubic spline interpolation on GPU (nvcc): NaN between and after the
  collocation points, for the material properties as well as for the
  hardenings defined by data (as soon as there is plasticity).
  Cause: nvcc 13.0 compiles without diagnostic, into an empty function, any
  `constexpr` function using `__builtin_constant_p`, which libstdc++ 14
  calls in `std::advance`; `std::advance`, `std::next` and `std::prev`
  therefore do nothing on GPU (reproducer
  `results/nvhpc-bugs/nvcc-builtin-constant-p.cu`). Internal functions
  in `CubicSpline.ixx`, the only device code of TFEL concerned;
- DigitSeparatorTest: the flow rule did not compute `df_dseq`,
  which the Newton solver read without initialisation. The result depended on the
  leftovers of the stack (gap in OpenMP with several threads, with compare compiled
  at -O3). Shown with `-ftrivial-auto-var-init=zero` and `pattern`.

Defects of the comparison program, fixed: vector or tensor external
variables given variable by variable (one field per
component instead of one field per variable). New: option
`--parameters` (parameters changed in both libraries).

Coverage added to the `generic-parallel` suite:

- `TensorialVariablesTest`: arrays of gradients and of forces, external
  and internal variables of all the types (scalars, vectors, symmetric
  tensors, tensors, arrays);
- parameters changed at run time for a behaviour and a model;
- `CubicSplineHardeningTest`: hardening interpolated by a cubic
  spline, plastic loading.

Each new test was checked both ways: it passes with the
fixes and fails without them (old comparison program, old
sizes, parameters set in the reference only).

The remaining gaps come from roundings: `scripts/roundoff-triage.py`
replays on the host the worst comparison of each gap, the
library of the batch compiled with and without FMA contraction. For each one,
no gap without FMA, a gap of the same order with it. At large amplitudes,
some behaviours are at the convergence limit of their local solver:
a rounding difference there makes a point fail on one side and not on
the other (for example Chaboche2 at amplitude 1e-2: 101 failures out of 1009
in the reference, 2 more on GPU, 4 different statuses on the host with
FMA). Even on the host, g++ contracts the two libraries of a behaviour
differently: a few gaps also appear in OpenMP, and disappear
without FMA.

Remaining inconclusive cases: intended failures (MaterialPropertiesBoundsCheck,
ExternalStateVariablesBoundsCheck, StateVariablesBoundsCheck,
InitializationFailureTest: same status in the reference and the batch), and
behaviours whose synthetic values make the reference diverge (some
are covered by the suite with physical values). The tests of the
warnings of MFront WarningTest80 to 82 name `f` and `df_dseq` without
computing them: results undefined by construction.

Sharing a GPU between processes: PolyCrystals needs 113 to 121 KB of
stack per thread, that is 14.8 to 15.9 GB reserved on the RTX 5080 whatever the
number of points. Another process (even its CUDA context alone) is enough
to run out of memory: the interface cleanly returns -3 "out of
memory". Six simultaneous processes on an ordinary behaviour (stack of 3.4 KB,
100,003 points each): no error. Open question for MGIS: 7
"illegal memory access" in about 1000 launches, only when 6
processes each reserved about 15 GB of stack out of 16 GB; never
alone, never under compute-sanitizer. A behaviour whose stack takes almost
the whole GPU must have it to itself.

Results:

- final execution pass, amplitudes 1e-4, 1e-3 and 1e-2 (105 min):
  openmp-host-gcc 467 OK out of 485, cuda-nvcc 401 OK out of 485; all the
  other results explained (76 rounding gaps confirmed by the triage,
  72 on GPU and 4 on the host,
  WarningTest80 to 82, 10 inconclusive cases identical in both
  configurations, one timeout at amplitude 1e-2, one transient lack of
  GPU memory);
- material properties: 56/56 on openmp-host-gcc and cuda-nvcc;
- `generic-parallel` suite at -O3: openmp-host-gcc, stlpar-parunseq-gcc
  and kokkos-host 531/531 (before the addition of CubicSplineHardeningTest),
  cuda-nvcc 540/540 (with it);
- compilation sweeps of the 608 files with AdaptiveCpp (28 min),
  nvc++ -stdpar=gpu (329 min) and nvc++ -mp=gpu (335 min): 541 OK, 60
  invalid files and 7 fragments, like nvcc and clang, file by
  file. With nvc++ -mp=gpu, the polycrystals of 30 grains
  (BerveillerZaouiPolyCrystal, PolyCrystals) need about 4 h and
  5 GB in the front-end: `compile-sweep.sh` now has a maximal time
  per behaviour (`COMPILE_TIMEOUT`, 6 h by default, status TIMEOUT).

### 2026-10-07 (morning), large validations at -O3 (R4, first part)

TFEL of the branch (all the fixes of R3) reconfigured, installed and
validated by `scratch/validate3.sh`:

- full ctest of TFEL: 6965 tests out of 6965;
- matrix of the `generic-parallel` tests at -O3
  (`results/generic-parallel-tests-2026-10-07-r4-O3.txt`): 15
  configurations of behaviours (cuda-clang, cuda-nvcc, hip-clang-nvidia-platform,
  hip-hipcc, kokkos, kokkos-host, openmp-host-gcc,
  openmp-offload-nvhpc-gpu, stlpar-parunseq-clang, stlpar-parunseq-gcc,
  stlpar-parunseq-icpx, stlpar-parunseq-nvhpc-gpu, sycl-default-acpp,
  sycl-default-acpp-omp, sycl-default-icpx) and 16 of material properties,
  all passed, TDLS included;
- nvcc and clang compilation sweeps: 541 OK, 60 invalid files
  and 7 fragments, identical file by file to the reference sweep.

Rest of R4: switch to the mmm-gpu spack environment (everything is
installed in it except MFEM, postponed), then a new validation.

### 2026-10-07 (afternoon), switch to the mmm-gpu spack environment (R4, in progress)

`scripts/env.sh` chooses the environment according to `PROTO_ENVIRONMENT`:
`tfelgpu_local` (default, validated) or `mmm-gpu` (deps/spack/mmm-gpu:
Kokkos 5.2.2, AdaptiveCpp 25.10, HIP 7.2.3 and LLVM 20.1.8 from spack, CUDA
13.0.2 from the HPC SDK for nvcc and HIP). TFEL is built in
`build/tfel-mmm-gpu` and installed in `install/tfel-mmm-gpu`;
`build-tfel.sh`, `scripts/run-generic-parallel-tests.sh`,
`scratch/dev-test.sh` and `scratch/validate3.sh` follow this variable. View
of mmm-gpu regenerated (`spack env view regenerate`: `acpp` was missing). MFEM
remains postponed.

Validation in mmm-gpu (interrupted at 14:40, machine shut down):

- full ctest of TFEL: 6965 tests out of 6965;
- matrix at -O3 (`results/generic-parallel-tests-2026-10-07-mmm-gpu-O3.txt`,
  partial): cuda-nvcc, kokkos, kokkos-host, openmp-host-gcc,
  openmp-offload-nvhpc-gpu, stlpar-parunseq-clang, stlpar-parunseq-gcc and
  stlpar-parunseq-icpx passed;
- cuda-clang and hip-clang: `libcudart.so.12` not found (clang 20 compiles
  with CUDA 12.9). Fixed in `scripts/env.sh` (libraries of CUDA
  12.9 after those of CUDA 13);
- hip-hipcc: headers of the system HIP 5.7 taken instead of those of
  HIP 7.2 (hipcc does not give the headers of HIP to nvcc on the NVIDIA
  platform). Defect of TFEL fixed (`cmake/modules/tfel.cmake`, patch 0066):
  directory of the headers of HIP added to the options of hipcc if it is not an
  implicit directory of the C++ compiler (unchanged for HIP in /usr).

Resuming:

```bash
PROTO_ENVIRONMENT=mmm-gpu ./build-tfel.sh install -j12   # reconfiguration (patch 0066)
PROTO_ENVIRONMENT=mmm-gpu scratch/validate3.sh 2026-10-07-mmm-gpu
```

Then, if everything passes: mmm-gpu by default in `scripts/env.sh`, R4 ticked.
Next: performance campaign (large batches, multicore CPU against GPU,
several families of behaviours), clang-tidy and independent review, then D1
(MGIS).

### 2026-10-07 (evening), mmm-gpu environment validated (end of R4)

Full validation in the mmm-gpu spack environment, after the
fixes of the afternoon (libraries of CUDA 12.9 in
`LD_LIBRARY_PATH`, headers of HIP given to hipcc):

- full ctest of TFEL: 6965 tests out of 6965;
- matrix of the `generic-parallel` tests at -O3
  (`results/generic-parallel-tests-2026-10-07-mmm-gpu-O3.txt`): 15
  configurations of behaviours and 16 of material properties, all passed,
  including hip-hipcc and hip-clang with HIP 7.2.3, kokkos with Kokkos 5.2.2 and
  sycl-default-acpp with the AdaptiveCpp of spack;
- nvcc and clang compilation sweeps: 541 OK, identical file by
  file to the reference sweep.

Additional upstream fix of TFEL (patch 0067): argument of the `endif`
of the hip-clang configuration in `tfel.cmake` (CMake warning at
each configuration). Note: the substitutions of TFEL are in the
CMake cache; a reconfiguration that changes one of their values requires
`cmake -U "tfel-substitutions*" .` (done by `scratch/validate3.sh`).

mmm-gpu is now the default environment of `scripts/env.sh` and
of `build-tfel.sh`; `PROTO_ENVIRONMENT=tfelgpu_local` keeps the old one
(also validated on the morning of 07/10). TFEL is therefore validated in both
environments.

Proposed next steps: performance campaign (large batches, multicore CPU
against GPU, several families of behaviours), clang-tidy and independent
review, then D1 (MGIS).

### 2026-10-08 (night), clang-tidy and independent review (R5)

clang-tidy 20 (bugprone, clang-analyzer, performance checks and
a few modernize, readability and cppcoreguidelines ones) on the sources of
the interface in MFront and on the generated code of five behaviours (behaviours
in small and finite strains, model, material properties), with the
headers of the interface. No defect in the headers executed during the
computation. In the generator: dead variables (upstream code of the material
properties) and two useless copies, fixed. The rest is cosmetic.

Two independent reviews, one of the headers executed during the computation,
the other of the generator. Defects confirmed and fixed:

- behaviours and models:
  - an empty batch whose fields are null (data of an empty
    `std::vector`, process without any point of a material) was refused (`-2`);
  - invalid computation options of a finite strain behaviour
    gave a failure at each point (`-1`), which pushes a finite
    element code to reduce the step endlessly. They are now
    refused for the batch (`-2`);
  - the initialize functions and the post-processings overwrote the
    step factors `rdts` proposed by the integration;
  - the output status was not complete after an exception thrown
    before the check of the data;
  - SYCL, `automatic` policy: the copied extent of a field used
    the maximal sizes and not those of the options. For a finite
    strain behaviour with the Cauchy stress, the copy read
    then rewrote values after the end of the array of the caller;
  - SYCL, `automatic` policy: two outputs stored in the same area
    (stresses and internal variables of a point stored one after the other)
    overwrote each other on return. Only the values of one field are copied back;
  - SYCL: the field of a group given variable by variable was copied
    whereas the contract says it is ignored;
  - SYCL: no ordering between the kernel and the work already submitted to an
    out of order queue. The asynchronous computations now require an in-order
    queue (`-2` otherwise); in synchronous mode, an out of order queue is
    synchronised before the kernel;
- material properties (upstream code):
  - with `disable_runtime_checks`, the variant with strides did not compile for
    any backend (`mfront_errno_old` not declared);
  - SYCL: the device memory was never freed (leak at
    each call); extent of the fields with strides too large; values
    located between the values of an output with strides overwritten on return;
  - the array of the pointers of the arguments, often on the stack of
    the caller, was read in the kernels (Kokkos, OpenMP, stdpar, SYCL).
    The pointers are now copied into the kernel. This defect was
    hidden by the development card: the driver enables HMM on it, which
    makes all the host memory readable by the GPU;
  - SYCL transfers done outside the `try` block; include guard without the
    name of the backend; bounds and default values written with 6 digits;
    one CUDA block too many; wrong error message; stray file.

New defect of nvc++ 25.11 found while fixing: a lambda that captures
a `std::array` of pointers, passed through `map(to : ...)` to an OpenMP
region, reads wrong pointers on GPU (-O2, -O3), or does not compile (-O0,
-O1). Reproducer `results/nvhpc-bugs/omp-lambda-pointer-array.cxx`. The
pointers are captured one by one.

Tests added: empty batch with null fields and invalid options (comparison
program, all the backends), compilation of a material property with
the checks disabled. Each one fails with the old
libraries and passes with the new ones.

Documented without fix: the `@SharedData` are copied once, on the
device of the first call (a process with several GPUs cannot
share the library); outside CUDA, HIP, SYCL and Kokkos, the heap must be
accessible by the device.

Still open, less serious: `generic_parallel_backend` symbols
and API version missing for the material properties, duplication between
the behaviour and model interfaces and between the two kernels of the
material properties, encoding of `rdt` written twice (C and C++),
`device_stack_size` ignored without message if the CUDA runtime is missing.
The reviews are in `scratch/review-r5` (not versioned).

Validation at -O3 with the development tree:

- behaviours: openmp-host-gcc, cuda-nvcc and kokkos 601/601 (540
  comparisons, 50 empty batches, 11 invalid options), sycl-default-acpp
  627/627 (RTX 5080, host memory: `automatic` copies exercised);
- material properties: stlpar-parunseq-gcc, openmp-host-gcc, kokkos-host,
  kokkos, cuda-nvcc, hip-hipcc, sycl-default-acpp, sycl-gpu-acpp,
  stlpar-parunseq-nvhpc-gpu and openmp-offload-nvhpc-gpu passed.

`scratch/dev-test.sh` now produces all the configuration files
named by a test (sycl-gpu-acpp uses the one of
sycl-default-acpp).

Limit of the validation on this machine: the NVIDIA driver enables HMM on it,
which makes all the host memory accessible to the GPU. An access of the GPU to
the stack or to the heap of the host works there, whereas it crashes on a
machine without HMM. To make sure, the GPU tests would have to be run again
with the `nvidia_uvm` module loaded with `uvm_disable_hmm=1` (administrator
rights).

### 2026-10-08 (night, continued), shared data, symbols, duplication (R6)

Three points left open by R5.

- `@SharedData`: one copy per GPU (current GPU in CUDA and HIP), per
  device and per context in SYCL. Kokkos keeps one copy (one GPU per
  process). The copy is stored under the address of the function that
  provides the data (`getSharedDataCopy<&Behaviour::mfront_get_shared_data_x>`,
  cache `SharedDataCopies<T, Device>`). The integration, the initialization
  functions and the post-processings therefore share the same copy.
  Cost: one query of the current GPU per call. New SYCL test: the
  polycrystalline behaviour is integrated on the GPU, then on the CPU.
- material properties: symbols `<name>_generic_parallel_backend` and
  `<name>_generic_parallel_api_version`, as for the behaviours. The version
  is defined in the new C header
  `MFront/GenericParallelMaterialProperty/MaterialProperty.h`. The tests of
  all the backends check these symbols.
- duplication: a single code for the material property kernels of all
  the backends (1385 lines before, 1013 after). The two copies had
  diverged: order of the checks, messages, restoration of `errno`.
  The behaviours and the models share the beginning and the end of the
  generated functions, the file and library names, and the choice of the
  backend. An unknown backend gives the list of the available backends. On
  the host side, the time step factor is encoded by the C functions of the
  contract.

The generated code of all the test behaviours was compared before and
after, on 6 backends: only the intended changes appear.

Full validation at -O3, with HMM disabled (`uvm_disable_hmm=Y`):

- installation, TFEL tests: 6965/6965;
- behaviours: 15 configurations out of 15;
- material properties: 14 out of 16. cuda-nvcc and cuda-clang failed
  because of the upstream test `test-cuda.cu`: the second call of `test9`
  wrote the result into host memory from a kernel (compute-sanitizer:
  invalid write). HMM hid this defect. The test now uses device memory:
  cuda-nvcc and cuda-clang pass, compute-sanitizer reports no error. The
  HIP test already uses managed memory;
- compilation sweeps: nvcc 541 OK, 60 INVALID, 7 FRAGMENT, as before;
  clang identical (541 OK).

Without HMM, no hidden access of the GPU to host memory appears in the
code of the interface.

### 2026-10-08 (night, end), parameters of the behaviour variables, rotations (R7)

Three remaining limits.

- String parameters of the models: nothing to fix for the GPU. A model
  with a string parameter already fails at parsing on CPU, with the
  `generic` interface. The reproducer `reproducers/device-strings` shows
  in passing that `std::string_view` does not work in GPU code: wrong
  results without any message with nvcc and hipcc, errors with nvc++
  (calls to `std::char_traits`). A small type without `char_traits` works
  with the six compilers.
- Parameters of the behaviour variables: each GPU compatible class has a
  structure `mfront_parameters_values`, with its parameters and those of
  its behaviour variables (recursively). It is filled on CPU from the
  initializers (`mfront_get_parameters_values`), then applied on the
  device (`mfront_set_parameters_values`). The factories of behaviour
  variables receive it too. The parameter structure of the
  `generic-parallel` interface derives from it: no more duplicated types,
  no more refusal of changed values. New test: the parameter
  `hardening_factor` of `BehaviourVariableTestPlasticity` is changed by
  file, then the CPU and GPU results are compared. A version that ignores
  the value fails on GPU (gap of 1.5e-3).
- Rotation of the tangent operator: all the blocks are rotated, on CPU
  (`generic`) as on GPU. The array blocks are rotated sub-block by
  sub-block, because MFront stores the sub-blocks one after the other. The
  pairs without a TFEL type (vector and tensor, tensor and symmetric
  tensor: TFEL has no `change_basis` for `st2tot2`) go through a generic
  product `Q1 K Q2ᵀ`, usable on GPU
  (`GenericBehaviour/RotateTangentOperatorBlock.hxx`). The generated code
  of the existing cases does not change. New test behaviour
  `OrthotropicTangentOperatorBlocksTest`: 9 blocks of all types.

Checks limited to the touched code:

- reproducer `reproducers/rotations`, independent of the code that rotates
  the blocks: `rotateTangentOperatorBlocks(K) x =
  rotateThermodynamicForces(K rotateGradients(x))` on the 5 hypotheses
  (gap of the order of 1e-15), in-place rotation identical. A deliberately
  wrong version is detected;
- targeted comparisons of the orthotropic behaviours
  (`scratch/targeted-test.sh`), CPU and GPU: 137 OK;
- full suites of the behaviours, which also contain the comparisons of
  the rotations of the new behaviour: openmp-host-gcc and
  cuda-nvcc 637/637, sycl-default-acpp 664/664, stlpar-parunseq-nvhpc-gpu
  637/637 with the workaround below.

New defect of nvc++ 25.11 (`-stdpar=gpu`), revealed by the
stlpar-parunseq-nvhpc-gpu suite: linking the device code of the behaviour
library fails ("nvlink error : Multiple definition of
'__staticinit_110045_…'", between `HosfordPerfectPlasticity` and
`InitializeFunctionTest`). nvc++ stores some local constants of the device
code (here the empty lambda of `assign` and `addAndAssign` of the TFEL
arrays) in visible global data, named by a counter specific to each file.
Two files that give the same counter to the same inline function
conflict. The parameter code only shifted the counters: the R6 validation
passed by chance. Without the rotation header, the conflict remains
(counters shifted by 8). The reproducer
`results/nvhpc-bugs/stdpar-staticinit-link.cxx` (two copies of a four-line
file) fails from -O0 to -O3. Workaround: the device code of each file is
linked alone (`-gpu=nordc` in `config-stlpar-parunseq-nvhpc-gpu.json`).
This is possible because the generated sources do not call any device
function defined in another file. Offloaded OpenMP does not emit these
symbols.

### 2026-10-08 (morning), uniform error handling (R8)

Principle, the same everywhere: the path without error costs nothing. The
kernels write into the statistics or the statuses only if an event
occurs at a point (failure, unreliable result, decrease of the time step
factor, violated bound, invalid result). No character string on the GPU:
the detailed message of a point is obtained by replaying it with the
`generic` interface.

State before (listed in `scratch/r8-map.md`, not versioned):

- behaviours, CUDA/HIP synchronous path: allocation, two copies and
  deallocation at each call; SYCL: allocation, three waits and
  deallocation at each call;
- material properties: status array allocated at each call by all the
  backends (`cudaMallocManaged` then `cudaFree`, which synchronizes the
  GPU; `std::vector`, hence `cudaMallocManaged` with nvc++; shared
  `Kokkos::View`; shared USM vector), then `cudaDeviceSynchronize`;
- no index of the failing point: the status array had to be scanned point
  by point.

After:

- C contract of the behaviours in version 2: `first_failure_index` in
  `mfront_gpb_BatchStatistics` and `mfront_gpb_OutputStatus`. It is the
  smallest failing index, whatever the order of the computations (atomic
  minimum or reduction). The message gives the number of failures and this
  point;
- CUDA/HIP: one vote per warp; the warps without an event no longer do any
  reduction or atomic operation;
- new shared header `MFront/GenericParallel/StatusBlockPool.hxx`: pool of
  blocks accessible by the host and the GPU, allocated at the first call
  on a GPU then reused, one block per ongoing call (concurrent calls are
  possible). It serves the synchronous path of the behaviours (CUDA/HIP:
  managed memory, SYCL: shared USM, a single wait) and the statuses of the
  material properties of all the backends (managed, heap, shared USM,
  `Kokkos::SharedSpace` freed at `finalize`). The host initializes the
  block and reads it after the synchronization: without an event, no data
  passes between the host and the GPU;
- CUDA/HIP material properties: synchronization of the stream instead of
  the whole GPU;
- synchronous Kokkos, OpenMP and parallel algorithms: reductions
  unchanged, with the index of the first failure in addition.

Test defect revealed: the HIP test of the material properties assumed
that the memory from `hipMallocManaged` was zero. The allocations freed at
each call hid this. The memory is now set to zero, like the
`thrust::device_vector` of the CUDA test.

The contract of the material properties (`mfront_gmp_OutputStatus`, shared
with the `generic` interface) has no field for a point index: they say
which bound is violated, not where.

Measurement (`reproducers/sync-call-overhead`, RTX 5080, elastic
behaviour, data in managed memory): a synchronous CUDA call goes from 135
to 29 µs for 1000 points, and from 546 to 442 µs for 100000 points. The
fixed cost per call decreases by about 105 µs.

Material properties, at the same time: with the `Warning` policy, a
violated bound now gives the status 1 (computed result, to be used with
caution), like the `generic` interface. The status used to stay zero. The
codes specific to the interface (-2, -5, -6, -7) are documented. New test
`test14` (C++ and CUDA tests).

Checks:

- targeted tests of the behaviours with failures
  (`ImplicitNortonWithEnergies` with the strict policy, model
  `GenericParallelModelTest`, initialize functions, post-processings):
  50/50 on cuda-nvcc, hip-hipcc, sycl-default-acpp, kokkos,
  openmp-host-gcc, stlpar-parunseq-gcc. A version without the atomic
  minimum of the index is detected (12 failures);
- material properties: cuda-nvcc, hip-hipcc, stlpar-parunseq-gcc,
  openmp-host-gcc, kokkos, sycl-default-acpp;
- full validation without HMM (`validate3.sh 2026-10-08-r8`): TFEL tests
  6965/6965, -O3 matrix 31/31 (stlpar-parunseq-nvhpc-gpu included, with
  `-gpu=nordc`), nvcc and clang sweeps 541 OK, as before.

In targeted mode, `scratch/dev-test.sh` skips the comparisons called
directly by a test (behaviours not built).

### 2026-10-08 (morning, continued), MGIS: loading and computation of batches (D1, steps 1 to 3)

Branch `gpu-proto` of `mgis/` (base `818f0f9`). The patches are in
`patches/mgis`, regenerated by `scripts/sync-mgis.sh`. MGIS is built in
`build/mgis` against the TFEL of the prototype. The spack TFEL must be
kept out: `CMAKE_PREFIX_PATH` found it first.

- Step 1, loading. A `generic-parallel` library was loaded without error
  and its batch functions were called with the arguments of one point
  (undefined behaviour). `BehaviourDescription` now reads the interface,
  the backend and the version of the structures (2). `Behaviour` holds the
  batch pointers (`batch_b`, `batch_f`, three rotations).
  `LibrariesManager` refuses the one-point functions for a batch library,
  and conversely. MGIS keeps its own copy of the C structures
  (`BatchDataView.hxx`), checked against `BatchData.h` by a memory layout
  test.
- Step 2, integration. `integrate(BatchExecutionPolicy, MaterialDataManager,
  options, dt[, b, e])` builds the view of the batch on the arrays of the
  manager, without copy. The material properties and the external
  variables are given variable by variable: uniform value or value per
  point. The existing functions switch by themselves, with the policy of
  the manager (`batch_execution_policy`). The statuses -2 and -3 are
  raised as exceptions, never taken for a convergence failure. A failure
  gives the index of the first point (`n`). `integrate_debug` is refused:
  the `generic` library is the one used for diagnosis.
- Step 3: initialize functions, post-processings and rotations of
  batches, with the same switch (thread pool included).

GPU: with the `cuda` backend, the caller provides all the memory as
managed memory (MGIS initializers, `std::span`). The uniform values kept
by MGIS are in host memory: they are refused with a clear message for
`cuda` and `hip`. Test `GenericParallelCUDATest` (option
`enable-generic-parallel-cuda-tests`): Norton over 10 steps and 1000
points, elasticity with a Young's modulus per point, HMM disabled. Results
equal to those of the point by point computation up to rounding,
compute-sanitizer without error, kernels seen by Nsight (0.5 ms for 1000
points).

MGIS defects found along the way, fixed in separate commits:

- the axes of the rotation matrices were refused when their number of
  values was a multiple of the dimension (inverted test): no
  `RotationMatrix2D` nor `RotationMatrix3D` could be built. The number of
  points was also compared to a remainder instead of a quotient;
- intermittent crash at the end of the process: the libraries were
  unloaded while the TBB threads (parallel algorithms of libstdc++) were
  still alive. They are no longer unloaded (`RTLD_NODELETE`);
- a material property could not be given by a `std::span` of a single
  value, while the evaluation already treated it as uniform. It is the
  only way to place a uniform value in memory accessible by the GPU.

Noted without a fix: `executeInitializeFunction` point by point does not
set `exit_status` to -1 on failure (only `n` is filled). The batch path
reports the failure.

Checks: 44 MGIS tests (41 original ones, loading, batch/point comparison
on CPU, comparison on GPU). Each comparison case was checked against a
deliberately wrong version (range shift, forced uniform rotations, wrong
stress measure suffix). An orthotropic finite-strain behaviour was added
to the MGIS tests for the rotations with a stress measure and a tangent
operator.

Step 4a: the location of the arrays (`MemoryLocation`: host, GPU, shared)
is declared by the initializer of the manager and checked. A `cuda` or
`hip` backend refuses data in host memory. Outside host memory, all the
arrays of the states must be given, as well as those of the tangent
operators and of the speed of sound if they are used: otherwise MGIS would
allocate them in host memory. `update` and `revert`, which copy on the
host, refuse data on the GPU. Each refusal is tested by its message. The
header `BatchExecutionPolicy.hxx` was not installed: fixed.

Rest of the plan: memory resource (step 4b, which will make automatic what
the GPU test does by hand), asynchronous execution and `update`/`revert`
on the GPU (step 5), `generic-parallel` material properties and
documentation.

### 2026-10-08 (afternoon), memory space (R9) and review of MGIS (D1, step 4b)

TFEL (R9), two commits plus one for dates:

- each library exports the memory space of its fields,
  `<f>_generic_parallel_memory_space`: `host`, `device` or `any`. The
  value is fixed at compile time. `cuda` and `hip`: `device`. `stlpar`
  and `openmp` on the host: `host`. nvc++ (`-stdpar=gpu`, `-mp=gpu`):
  `device`, or `any` with unified memory. `kokkos`: depends on whether the
  execution space can access host memory. `sycl`: `any` with automatic
  transfers, `device` otherwise. Header
  `MFront/GenericParallel/MemorySpace.hxx`. The material properties
  export the same symbol. `compare` and the property tests check its
  value;
- the `stlpar` configurations linked to TBB (gcc, clang, icpx) add
  `-Wl,-z,nodelete`: the library stays loaded until the end of the
  process. The TBB threads were still running when it was unloaded at
  exit. MGIS no longer needs `RTLD_NODELETE` for all the libraries.

MGIS: independent review of the branch (19 points) and three points noted
before. Fixed:

- the list of backends written in MGIS (`cuda`, `hip`) is replaced by the
  TFEL symbol: `isMemoryLocationAccessible(behaviour, location)`.
  The cases missed before (Kokkos on GPU, offloaded OpenMP, SYCL without
  transfers, nvc++ without unified memory, data on GPU with a host
  backend) are covered;
- the one-point functions (integration, `integrate_debug`, initialization,
  post-processing, rotations, C bindings) called a null pointer for a
  batch library. They raise an explicit error;
- the values kept by MGIS (uniform value, `LOCAL_STORAGE`) are refused
  when host memory is not accessible to the computations. Before, only the
  scalar was checked: a `std::vector` passed and the CUDA kernel crashed,
  with the CUDA context lost;
- the mass density is passed only for the speed of sound. It can be given
  by a `std::span` of one value. A uniform property in `LOCAL_STORAGE` is
  stored as a scalar;
- arrays of variables: the review believed they were refused. MFront
  exports each element as a variable (`v2_esv[0]`, `v2_esv[1]`) and
  refuses array types for the material properties and the external
  variables. One field per variable, as in TFEL: the per-element split,
  added then tested, was wrong for an array type, so it is removed;
- initialization, post-processings and rotations of batches take an
  execution policy (queue), like the integration. The `sycl` rotations
  always failed: null queue. The batch rotation functions are loaded by
  typed accessors, and the finite-strain suffixes are written once (an
  unknown value raises an error instead of falling back to PK1);
- the state managers know the location of their arrays. `update`,
  `revert`, `updateValues`, the extraction, `save`, `restore` and the
  point by point computation refuse arrays on the GPU. A manager without
  any point, or a model without a tangent operator, is accepted outside
  host memory;
- a failure gives the absolute index of the first point in the message,
  and an unexpected status raises an exception. The unused asynchronous
  statistics functions are removed: they will come back with step 5;
- CMake: CUDA architecture `native` only with CMake 3.24; batch sources
  regenerated when `mfront` changes (the old sources, without the new
  symbol, remained).

Step 4b: the manager initializer takes a memory resource
(`std::pmr::memory_resource*`). MGIS uses it to allocate the arrays that
are not given. Their values are initialized by the host: the resource must
allocate memory accessible by the host (managed or shared). For memory on
the GPU only, all the arrays must still be given. Resources provided as
headers, compiled by the calling code: CUDA and HIP managed memory, SYCL
and Kokkos shared memory. They check the requested alignment. Cost: one
pointer per array, one virtual call per allocation, nothing per
computation.

Tests added: range with modified gradients, post-processing on a range,
partial failures caused by an external variable outside its physical
bounds (absolute index), all the variables of
`TensorialExternalStateVariableTest` checked value by value, finite strain
with four pairs of stress measure and tangent operator, speed of sound, 2D
and 3D rotations per point, functions with a policy, memory resource and
refused locations. On GPU: managed memory resource with no array given,
then integration, initialization, post-processing and rotation on a CUDA
stream that is not the default stream. Two deliberately wrong versions
(stride of the points of the variables, relative index) are detected. The
44 MGIS tests pass, HMM disabled. The HIP, SYCL and Kokkos headers
compile.

An older backup in `backups/` was overwritten by a rename. Its commits
remain in `tfel-dev` (branches and reflog).

Remaining for MGIS: asynchronous mode and `update`/`revert` on the GPU
(step 5, after mfem-mgis), `generic-parallel` material properties (step
6), documentation (step 7). Per-point statuses are not exposed. The
rotations given by `RotationMatrix2D`/`3D` build the matrices on the host
at each call: for repeated rotations, give them as an array.

### 2026-10-08 (afternoon, continued), MGIS: material properties and batch models (D1, step 6)

Everything that the `generic-parallel` interface of TFEL produces is now
usable by MGIS, as for the behaviours.

- Material properties. At loading, MGIS reads the interface, the version
  of the functions (1), the backend and the memory space. It loads the
  function for a set of points (`batch_fct`) instead of the one for a
  single point. Before, a `generic-parallel` library gave a pointer with
  the wrong signature. Each function refuses the other type of library.
- `evaluate(values, property, inputs[, location, policy])` evaluates
  a property on a set of points. The inputs are uniform or given per
  point. A single call for a `generic-parallel` property, point by point
  computation on the host otherwise: same switch as `integrate`. The
  location is checked against the memory space. The statuses keep the
  meaning of the `generic` interface (0, 1 out of bounds with the Warning
  policy, negative on failure).
- Models: the `generic-parallel` interface of the models exports the same
  symbols as the behaviours. MGIS already loaded them through `load`: a
  test checks it, with a parameter modified by `setParameter`.
- `MemoryLocation` moves into the `mgis` namespace
  (`MGIS/MemoryLocation.hxx`): the access rule is written once for the
  behaviours and the properties. The CMake functions of the batch
  libraries (behaviours, models, properties) share their code
  (`cmake/modules/generic-parallel.cmake`).

Asymmetry that comes from TFEL: the material property functions take no
execution options (queue, asynchronous statistics). They compute on the
default queue of their backend, synchronously. `evaluate` therefore has no
overload with a policy. Adding one would require a change of the contract
of the material properties in TFEL (version 2).

Tests (minimal): properties on the host (`stlpar`) and on the GPU
(`cuda`, managed memory), per-point and uniform values, physical bounds,
statuses of the three bound policies, locations, invalid arguments;
version and layout of the status checked against TFEL; model `ode_rk54`
as a batch on the host and on the GPU. Two deliberately wrong versions of
`evaluate` and an ignored parameter are detected. 45 MGIS tests, all
green, HMM disabled.

The message `Error at line ... of file 'src/targets.lst'` during the build
comes from parallel calls of `mfront` in the same directory. It also
exists in the original MGIS build. It has no effect on the sources.

Remaining for MGIS: asynchronous mode and `update`/`revert` on the GPU
(step 5, after mfem-mgis), documentation (step 7).

### 2026-10-08 (evening), execution options of the material properties (R10)

The functions of the `generic-parallel` material properties had no
execution options. The SYCL backend added the queue as a specific
argument, in the middle of the others: the signature depended on the
backend. The MGIS `evaluate`, written for the common signature, would have
called a SYCL library with shifted arguments. CUDA launched on the default
stream, Kokkos on the default instance with a global `Kokkos::fence()`.

TFEL, version 2 of the functions (`MFRONT_GPMP_API_VERSION`):

- a last argument `const mfront_gpmp_ExecutionOptions*` (null for the
  default options): queue and preservation of the stack size, same
  meaning as for the behaviours. CUDA and HIP stream, `sycl::queue*`
  mandatory (status -5 otherwise), instance of the Kokkos execution space
  (wait on this instance only). The parallel algorithms and OpenMP have no
  queue;
- the SYCL-specific argument and the mechanism of specific arguments
  disappear: a single signature for all the backends;
- the stack size code (`stlpar`, offloaded OpenMP, SYCL) is written once
  and respects `preserve_stack_size`;
- the evaluation stays synchronous. Asynchronous evaluation would require
  another way of reporting the bounds (the message names the variable, it
  is built after the computation): not done.

Tests: the six test files pass the options. A non-blocking CUDA and HIP
stream is checked by a kernel that writes the arguments on this stream
after a delay, without the host waiting: on the default stream, the test
fails (checked). Kokkos instance (`kokkos-host` receives the host space of
its configuration), missing SYCL queue refused. Ten configurations, one
per backend: `stlpar-parunseq-gcc`, `openmp-host-gcc`, `kokkos-host`,
`kokkos`, `cuda-nvcc`, `hip-hipcc`, `sycl-default-acpp`, `sycl-gpu-acpp`,
`stlpar-parunseq-nvhpc-gpu`, `openmp-offload-nvhpc-gpu`. All pass.
The matrix and the full TFEL suite are not rerun.

MGIS: `BatchExecutionPolicy` moves into the `mgis` namespace, shared by
the behaviours, the models and the properties. `evaluate(policy, ...)`
exists like `integrate(policy, ...)`. The statuses for invalid arguments
and backend errors (-5, -6, -7) are raised as exceptions, as for the
behaviours. 45 MGIS tests, all green, including the evaluation on a
non-default CUDA stream.

Behaviours, models and material properties now have the same execution
interface in TFEL and in MGIS. Only difference: in TFEL, the behaviours
can compute asynchronously (statistics on the GPU), the material
properties cannot. MGIS uses asynchronous mode for neither of them (step
5).

### 2026-10-08 (evening, continued), asynchronous evaluation of the material properties (R10, continued)

In TFEL, the behaviours could compute asynchronously (statistics given by
the execution options), the material properties could not. This is
solved: synchronous or asynchronous at the choice of the caller, on all
the backends, for both.

- The execution options of the properties take the statuses of an
  asynchronous evaluation (`statuses`, array of
  `<name>_generic_parallel_number_of_statuses` integers set to zero by the
  caller). The function then returns as soon as the computation is
  submitted to the queue (CUDA, HIP, Kokkos, SYCL if no value is copied).
  The output status only reports invalid arguments and backend errors.
- The analysis of the statuses is generated once, in an exported function
  (`<name>_generic_parallel_report`). The synchronous evaluation calls it
  itself; after an asynchronous evaluation, the caller calls it once the
  queue is synchronized. The report is the same in both cases, including
  the message that names the variable.
- The kernels do not change: they write the statuses only for a violated
  bound or an invalid result. Zero cost on success, as before.
- The structure is only published in the prototype: the version stays 2.

Tests limited to what changes, one configuration per backend (ten, all
green): asynchronous evaluation with a temperature outside its physical
bounds, report after synchronization, values of the other points. With
CUDA and HIP, a slow kernel writes the arguments on the stream before the
evaluation: the stream is still busy when the function returns
(asynchrony checked) and the values are right (order checked).

MGIS: the copy of the options structure follows TFEL (checked by the
layout test). MGIS uses asynchronous mode neither for the behaviours nor
for the properties (step 5): the two are symmetric.

### 2026-10-08 (evening, end), SYCL queues and thread blocks of the material properties (R10, end)

Two gaps remained between the behaviours and the material properties.
They are aligned on the behaviours.

- **Out-of-order SYCL queue.** The kernel must run after the computations
  already submitted to the queue, which may produce its arguments. An
  in-order queue guarantees it. The behaviours already handled the
  out-of-order queue: asynchronous computation refused, queue drained
  before the kernel in synchronous mode. The material properties did not:
  their kernel could run before the computations of its arguments. Same
  rule now (status `-5` for an asynchronous evaluation).
- **Threads per block with CUDA and HIP.** 64, fixed, for the material
  properties; 128 rounded to a multiple of the warp size for the
  behaviours. Same rule now: 128 by default, warp size read on the GPU at
  run time.

Tests: `test-sycl.cxx` checks the refusal of asynchronous evaluation on an
out-of-order queue (test 16). The test of the asynchronous evaluation uses
an in-order queue. Installed tests of the material properties, after
reinstalling TFEL: `cuda-nvcc`, `hip-hipcc`, `sycl-default-acpp`,
`sycl-gpu-acpp`, all green. MGIS: 45 tests out of 45 against this TFEL.

### 2026-10-08 (night), documentation of TFEL and MGIS (R11)

Part of the documentation dated from the official versions. It was
reviewed in full and updated, in three parts.

**Prototype documentation** (`docs/`). Restructured: overview, TFEL
details, MGIS details, validation, split for the official repositories,
plan for mfem-mgis. The old working documents are in `docs/history/`.
This log received an index and the table of phases.

**TFEL** (`tfel-dev`). Three independent reviews, then fixes:

- web pages of the two `generic-parallel` interfaces, release notes,
  pages of the `generic` interface and of `tfel-config`, menu, list of the
  documentations;
- code comments: memory space of the batches, SYCL queues, lists of
  backends, conventions of the material property backends, documentation
  written in the generated code (strides of the arrays, constructor
  parameters, order of the arguments);
- old errors in the touched files: misplaced comments, wrong parameter
  names, template parameters documented as parameters;
- section `@SharedData` added to the keyword pages of the behaviour DSLs.
  These pages are not regenerated in full: TFEL fixed typos directly in
  them, not in the help files they are drawn from. The page of the DSL
  `IsotropicStrainHardeningMisesCreep`, which contained an error message,
  is replaced by the generated page.

Each statement of the documentation was checked in the code. Nine
findings of the reviews were wrong or imprecise: they were corrected
before being applied. A CMake defect is noted, not fixed: the
`sycl-gpu-acpp` tests use the configuration files of `sycl-default-acpp`,
which are installed only if this configuration is also chosen.

**MGIS** (`mgis`). New page `generic-parallel.md`, linked from the home
page, the menu and the pages that mention it. Release notes, installation
and README updated. Each public function says what it does for a
`generic-parallel` library: switch, refusal, location of the arrays. Old
errors fixed: `update` described backwards, names of the rotations of the
finite-strain behaviours, badly annotated parameters.

Only comments changed in the code: checked by comparing the sources
without their comments (MGIS), by reading the diff and by compiling (TFEL
and MGIS).

### 2026-10-08 (night, continued), asynchronous computation in MGIS (D1, step 5)

TFEL offered asynchronous computation for everything. MGIS did not expose
it. This is done, and synchronous mode stays the default.

- An asynchronous policy (`mgis::AsynchronousBatchExecutionPolicy`) takes
  the queue and the stack of the synchronous policy. It adds the memory
  resource that allocates the status of the computation, and its location.
- Each batch function accepts it: integration, initialize functions,
  post-processings, rotations, evaluation of the material properties. It
  returns a pending computation (`PendingBatchComputation`,
  `PendingBatchEvaluation`).
- After the caller has synchronized the queue, `getResult` returns the
  result of the synchronous function. The code is common to both modes,
  except for obtaining the status. MGIS synchronizes nothing: it does not
  depend on any GPU library.
- The pending computation keeps the memory that the GPU writes or reads:
  the status, and the rotation matrices built from the material axes.
  These matrices are allocated by the resource of the policy: in
  asynchronous mode, they also work for a behaviour that computes on the
  GPU.
- MGIS checks that the status is accessible by the host (initialization,
  reading) and by the computations.
- The material properties load the report function and the number of
  statuses of their library.

Tests: on CPU (`stlpar`), same results as in synchronous mode, including
the index of the first failing point, refusal of a `generic` behaviour and
of a misplaced status. On GPU (CUDA, non-default stream, managed memory,
without HMM): integration, rotation from the material axes, bounds of a
material property. MGIS: 45 tests out of 45.

TFEL: the `sycl-gpu-acpp` tests use the configuration files of
`sycl-default-acpp`. They are now installed if either of the two
configurations is chosen (checked by a configuration with `sycl-gpu-acpp`
alone).

Still to discuss: the index of the first failing point for the material
properties (three options in [docs/tfel.md](docs/tfel.md#comparison-between-behaviours-and-material-properties)).

### 2026-10-08 (night, end), first failing point of the material properties (R12)

The `generic-parallel` material properties reused the status of the
`generic` interface, designed for a single point. They said neither how
many points had failed, nor which one failed first. It was the last
difference with the behaviours. Option chosen: a dedicated status,
modelled on the one of the behaviours.

TFEL:

- statistics (`mfront_gpmp_Statistics`) and output status
  (`mfront_gpmp_OutputStatus`) specific to the interface: number of
  failing points, the first of them, number of points out of bounds with
  the `Warning` policy, first failing check (to name the variable in the
  message);
- execution options with the members of those of the behaviours: queue,
  statistics (asynchronous), stack. The status array and the symbol
  `<name>_generic_parallel_number_of_statuses` disappear;
- each kernel computes the status of its point in a lambda, and merges it
  only if an event occurs: atomic operations for CUDA, HIP, SYCL and
  Kokkos, reduction for the parallel algorithms and OpenMP, as for the
  behaviours. Zero cost on success;
- new defect of `nvc++` 25.11 (`-mp=gpu`, -O0 to -O3): a nested lambda
  that captures by reference a pointer captured by the kernel does not
  write through this pointer. Workaround: capture by copy. Reproducer:
  `results/nvhpc-bugs/omp-lambda-reference-capture.cxx`.

The sign of the rank of the out-of-bounds variable was examined. The doc
of the `generic` interface says it is negative for a physical bound. The
code makes it negative as soon as the evaluation fails, including for a
bound with the `Strict` policy, and a test of the MGIS master expects
this. This behaviour of the code is kept, the `generic-parallel`
interface follows it, and the doc is fixed. Two defects of the `generic`
interface are fixed: message of the upper physical bound, `errno` not
restored.

MGIS: copies of the new structures (layouts checked at compile time),
`evaluate` returns a `MaterialPropertyEvaluationResult` (status, rank,
index of the first failing point, message), the same for both interfaces,
like `BehaviourIntegrationResult`. The index is the same with `generic`
and with `generic-parallel` (tested).

Tests: TFEL material properties on the 16 configurations (all the
backends and compilers), with the number of failures, the first failing
point and the unreliable results; TFEL `generic` tests (6412); MGIS 45
out of 45, CUDA included.

One difference remains, set aside: the optional per-point status array
of the behaviours does not exist for the material properties. MGIS
exposes it for neither of them.

### 2026-10-08 (night, closing), final audit of TFEL and MGIS (R13)

Two independent read-only reviews went through all the `generic-parallel`
code of TFEL (behaviours, models, material properties, six backends) and
of MGIS. Each finding was checked in the code before fixing. The findings
and how they were handled:

TFEL:

- Kokkos, synchronous behaviours: reduction, whereas the material
  properties used atomics. The behaviours move to atomics on a block of
  `Kokkos::SharedSpace`, like the properties. Both interfaces require this
  space (static assertion).
- SYCL, material properties: class of copies generated separately,
  different from the one of the behaviours. Statistics never copied (the
  GPU wrote into host memory at the first failure), allocation not
  checked, undefined output at the failing points, `wait()` instead of
  `wait_and_throw()`, different notion of accessibility. Both interfaces
  now share `SYCLDeviceCopies.hxx`.
- Thread stack: CUDA, HIP and Kokkos did not check it for the material
  properties, and `preserve_stack_size` had no effect there. Both
  interfaces share `GPUStackSize.hxx` and `KokkosStackSize.hxx`.
- The kernel becomes a template parameter of the function that reads its
  attributes. With an argument, the function was instantiated identically
  in several libraries. The dynamic linker resolved it to a single one.
  The CUDA runtime of that library, linked statically, did not know the
  kernel of the others: `invalid resource handle`. Found by the CUDA test
  of the material properties, which loads three libraries.
- Runtime exceptions (SYCL, Kokkos, parallel algorithms, OpenMP): status
  `-2` for the material properties, `-6` for the CUDA and HIP errors. All
  `-6` now.
- Asynchronous report: C function `mfront_gpb_report` for the behaviours,
  counterpart of `<name>_generic_parallel_report` of the properties. The
  synchronous report of the behaviours uses it: a single message text.
- Empty batch: accepted without any check in both interfaces (before:
  `-2` with SYCL without a queue, `-7` for a uniform output). `errno`
  restored at each return (RAII guard).
- SYCL behaviours with `@SharedData` and without a queue: `-3` instead of
  `-2`.
- Number of points as `size_t` in the CUDA, HIP and Kokkos kernels and
  with AdaptiveCpp (`int` before, limited to 2^31 points). AdaptiveCpp
  does offload the loop with a `size_t` index (checked).
- Doc: cost of the reductions (a few integers copied per call),
  statistics of `stlpar` and `openmp` accessible by the host.
- Tests added: `preserve_stack_size` (CUDA, HIP and Kokkos properties;
  CUDA and HIP behaviours), statistics in host memory with SYCL, SYCL
  behaviours with a null queue and an out-of-order queue, SYCL behaviours
  with the transfer policy `none` (configuration `sycl-default-acpp-2`),
  execution of a property compiled without checks, empty set.

MGIS:

- Move assignment of the pending objects removed: with different
  resources, it copied the statistics while the GPU was writing them
  (write after free). The move constructor remains.
- Material properties: a `generic-parallel` evaluation only returns `1`,
  `0`, `-1`, `-4`. Any other status raises an exception, already at
  submission in asynchronous mode (the `-2` was lost).
- Finite-strain conversions on a manager: check of the host access.
- Batch rotations: location of the arrays as an argument (`HOST` by
  default), like `evaluate`. `buildRotationMatrices` builds the matrices
  once in a given resource. In-place rotation documented.
- Location of the values given by `std::span` to the manager and of the
  arrays of the batch functions: the same as the arrays of the manager
  (documented contract). `allocatePostProcessingVariables` with a
  resource.
- Tests added: partial failures on GPU (first failure by atomics,
  synchronous, asynchronous, range), non-finite result (`-4`, host and
  CUDA), refused locations, rotation matrices built once, rotation without
  a policy on GPU, move of the pending objects.

Limitations kept, with their reason (see [docs/mgis.md](docs/mgis.md#current-limitations)):
point by point layout of the `MaterialDataManager` (an API on views is to
be defined with mfem-mgis), no variant with `Context`,
`MaterialFunctionManager` on the host, `update`/`revert` on the host, no
per-point status array, `generic-parallel` tests absent from the MGIS CI
as long as TFEL master does not have the interface.

Second review, on the diff of the fixes, with the report of the first
audit as input (status of each finding, defects introduced, remaining
asymmetries):

- TFEL: all the findings fixed, no class A or B defect introduced. Minor
  points handled: comments on the cost of the reductions, success status
  of an empty set (`1` or `0` depending on the convention), propagation of
  the SYCL errors by the asynchronous handler of the queue, possible leak
  if the description of a SYCL copy cannot be stored, C functions of the
  CUDA and HIP properties without `try`.
- MGIS: two holes found and fixed. The rotations without a policy passed
  the `HOST` location to the batch functions: they refused any `device`
  library, a regression from the previous step. They now take the
  location, like the others. `allocatePostProcessingVariables` without a
  resource gave host memory to a behaviour that cannot read it: now
  refused. A null resource has the same meaning everywhere (default
  resource).

Tests:

- material properties: the 16 configurations;
- behaviours: 3 configurations in full (`openmp-host-gcc`, `stlpar-parunseq-gcc`,
  `sycl-default-acpp-2`), the 13 others in a targeted pass (integration,
  failures, initialization function, post-processing, rotations, shared
  data for SYCL, SYCL queues, preserved stack on CUDA and HIP);
- MGIS: 45 tests out of 45 against the installed TFEL, CUDA included.

`scratch/dev-test.sh` did not copy the configuration file when
`test.check` passes it to `build-libraries.sh` (`sycl-default-acpp-omp`):
fixed.

When the MGIS test libraries were rebuilt with `-Wshadow`, a warning
appeared in the batch post-processing code (a local variable shadowed a
member): fixed.

### 2026-10-08 (end), choices for mfem-mgis (D2)

The plan for mfem-mgis is decided. It is described in
[docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md). In short:

- two goals: explore (measurements, variants) and propose clean code to
  the mfem-mgis developers. Everything can be selected by options, and the
  current path stays the default path;
- two entry points in MFEM: partial assembly of `NonlinearForm` and an
  operator specific to mfem-mgis, on the same batch kernels;
- two layouts of the data at the Gauss points: point by point and
  component by component. Today MGIS accepts only the first one: an API on
  views is the first step;
- three options for the tangent operator (matrix free, elastic fixed
  point, assembled tangent) and a nonlinear fixed point solver;
- Ginkgo in the solver factory, direct solvers kept on CPU;
- MFEM 4.9, the version used so far by mfem-mgis, to be compiled with
  CUDA. The specification of `deps/spack/mmm-gpu/spack.yaml` (4.10) will be
  changed before the compilation.

The sources of mfem-mgis (master `a011d97d6f`), mfem-mgis-examples
(master `a242567`) and mm-opera-hpc (main `c31c4bc`) are tracked in this
repository from their base, without modification
(`scripts/sync-repository.sh`).

### 2026-10-09, views API of MGIS (D2, step 1)

First step of the mfem-mgis plan: MGIS accepts arrays stored component by
component, or owned by the caller, described by views with arbitrary
strides. The proposal (version 4, after three reviews) is versioned in
[docs/history/design-views-mgis-2026-10-08.md](docs/history/design-views-mgis-2026-10-08.md).
The current state is described in [docs/mgis.md](docs/mgis.md#views-and-layouts).

What is done, in `mgis/` (patches `0024` to `0033`):

- `FieldView` (two strides, regular type, contract checked at
  construction), `MaterialStateView` and `MaterialDataView`, `make_view`;
- all the functions of the manager on a view, without a policy (explicit
  workspace), with a synchronous or asynchronous policy. The functions of
  the manager call those of the view: a single path;
- layout of a manager (`layout` of the initializer), followed by
  everything that is given per point. `save` and `restore` stay point by
  point;
- `MaterialFunctionManager` with an initializer and *strided coalesced*
  views; layout exposed in C, Fortran and Python;
- `FieldOperations`: `update`, `revert` and `updateValues` on the arrays
  of the caller, including on GPU. CUDA, HIP, SYCL, Kokkos and OpenMP
  headers, compiled by the caller;
- rotations and `evaluate` on views; the overloads on arrays go through
  them;
- no more heap allocation per call in the common case (buffers on the
  stack), pool recommended for the status of the asynchronous
  computations.

Defects found along the way and fixed: HDF5 restoration of an external
variable with several components (write outside the array), CUDA, HIP and
Kokkos memory resources unusable under a pool (alignment refused),
initializer without a constructor in Python, error message of the
rotations of the tangent. Found and not fixed: `MaterialAxisStorage` not
exported (copying a `RotationMatrix3D` outside MGIS does not link).

Checks: 48 tests of the main build (CUDA included), 94 of the build with
HDF5 and the C, Fortran and Python bindings, and the GPU headers compiled
and run with the five toolchains (`scripts/check-mgis-gpu-headers.sh`).
The results on views are bitwise identical to those of the manager, in
all layouts. Deliberate mutations of the code (copy omitted, stride
ignored, matrices read wrongly, incomplete update) are detected by the
tests.

MGIS documentation: page `docs/web/views.md`, batch page and release
notes updated.

### 2026-10-09 (continued), reviews of the views API and fixes

Before declaring the step finished, a coverage audit and three
independent reviews (correctness, design, then review of the diff of the
fixes). Patches `0035` to `0053` (commit `a2f8b21` wrongly announced
`0035` to `0057`).

- **Coverage**: the audit listed the untested cases (predictions, speed of
  sound, model, energies, thread pool, ranges, component by component
  manager beyond the simple case, views on GPU, expected refusals). All
  have a test, including a new test behaviour that computes the energies.
  Each computation on a view now starts from a copy of the state before
  the computation: a computation that does nothing is detected.
- **Cost**: the checks built the error messages at each call, and these
  were never displayed (about thirty allocations per call). A test counts
  the calls to `operator new`: no allocation left in the calling thread
  for the integration, the initialization, the post-processing, the end of
  step, the rotations and `evaluate`.
- **Defects fixed**: thread pool (a failing task let the others read
  destroyed objects), outputs of a failing point copied from the previous
  point, error message read from a destroyed exception, finite-strain
  conversions without a check of the output, copies of empty views, mass
  density (rule written), GPU kernels without a division per value, SYCL
  queue kept by reference.
- **API completed**: `getPoints`, tangent blocks of a view, `updateValues`
  with a context, extraction by operations on the manager, conversions on
  views of the manager, dedicated header for the workspace, refusal of a
  temporary behaviour.
- **Old MGIS defects fixed**: copies of the rotation matrices, shape of
  the tangent operators in Fortran, strided numpy arrays read as
  contiguous, context allocated by `update` and `revert`.
- **Left aside**: views in the C, Fortran and Python bindings (the
  prototype only works in C++), the check of the overlap between two
  different outputs (documented).

Lesson: the tests of the bindings are built by the `check` target, not by
`build-tests`. A mutation of the Fortran code went unnoticed with a stale
executable; it is detected once the `check` target is used.

Final checks: 50 tests (main build, CUDA included), 96 tests (HDF5 and C,
Fortran, Python bindings, through `check`), the five GPU toolchains
(`scripts/check-mgis-gpu-headers.sh`), and deliberate mutations detected
by the new tests.

### 2026-10-09 (night), MGIS frozen, list of solvers, spack environment

- **MGIS frozen.** New script `build-mgis.sh`: fresh build of
  `build/mgis` (50 tests, CUDA included) and installation in
  `install/mgis`, where mfem-mgis will find it; `build-mgis.sh full` for
  `build/mgis-full` (HDF5, C, Fortran, Python bindings: 96 tests). Both
  rebuilt from empty directories, green.
- **Environment leak fixed.** The first fresh build failed: the spack
  environment `tfelgpu_local`, activated in the shell, provided HDF5, and
  its view contains an old TFEL whose headers came before those of the
  prototype TFEL. `scripts/env.sh` and `build-tfel.sh` now neutralize this
  environment for `mmm-gpu`. A fresh configuration of TFEL gives the same
  cache as `build/tfel-mmm-gpu`, except for unused CUDA libraries (cuFile,
  nvfatbin) and the git hash displayed by TFEL, frozen at configuration
  time (`9f10d65a4` instead of `631019468`, cosmetic).
- **MGIS, patch 0054**: the CUDA behaviours of the tests are compiled with
  the options of the `cuda-nvcc` configuration of TFEL (nvcc warnings on
  the TFEL headers hidden, calls of host functions from the GPU turned
  into errors). No nvcc warning left.
- **List of solvers**: [docs/gpu-solvers.md](docs/gpu-solvers.md).
  Notable points, checked in the sources: the "Elasticity" strategy of
  BoomerAMG in mfem-mgis stops the program with GPU hypre
  (`SetElasticityOptions`); MFEM 4.10 wraps cuDSS and the distributed
  solvers of Ginkgo; AmgX 2.5.0 predates the fixes for CUDA 13; the AmgX
  wrapper of MFEM is double precision only.
- **Spack environment**: MFEM 4.9 with CUDA, Ginkgo, PETSc (Kokkos),
  SUNDIALS, AmgX, SuperLU_DIST, STRUMPACK, MUMPS, UMFPACK; cuDSS on its
  own. Problems met and solved: spack cache shared between two
  installations (index of the patches); BLAS of the HPC SDK chosen by
  default (OpenBLAS imposed); Kokkos 5 without the `cuda_lambda` variant
  required by the PETSc recipe (recipe modified in `mmm_gpu`); AmgX 2.5.0
  too old for CUDA 13 (commit of the main branch); OpenMPI with CUDA
  impossible without rebuilding LLVM and TBB (OpenMPI without CUDA, PETSc
  run with `-use_gpu_aware_mpi 0`); CUDA prefix whose `bin` was a link, so
  that nvcc took the headers of the SDK without cuBLAS
  (`make-nvhpc-cuda-prefix.sh` fixed, prefixes 13.0 and 12.9 redone).

### 2026-10-09 (night, continued), detailed plan of mfem-mgis and review

- **Map of mfem-mgis** (read-only): solve loop, integrators, data at the
  points, solvers, parameters, tests, with `file:line` references. It is
  the basis of the plan.
- **Detailed plan**: [docs/plan-mfem-mgis-detailed.md](docs/plan-mfem-mgis-detailed.md),
  eleven work packages (L0 to L10) with touched files, tests and completion
  criteria, proposed parameters, test strategy, measurement campaign.
- **Independent review** of the list of solvers and of the plan: 50
  findings, all integrated. The most important ones:
  - the `Mult` of mfem-mgis does not integrate the behaviour (it is an
    action of the mfem-mgis Newton): the nonlinear solvers of MFEM,
    SUNDIALS and PETSc require an operator that integrates it; in PETSc
    mode, a failure of the behaviour stops the program;
  - the Newton of mfem-mgis is not the one of MFEM (inexact Newton and
    line search would have to be written in it); the L-BFGS of MFEM
    ignores the preconditioner;
  - `mfem::forall` launches a CUDA kernel only in a file compiled by
    nvcc;
  - 27-point integration rule on the order 1 hexahedron (rule option
    added: about 3.4 times less computation with 8 points);
  - the thermodynamic forces are not rotated at the integration;
  - the current path assumes the point by point layout everywhere;
  - validity of the host and GPU copies of `mfem::Memory` around the
    calls to MGIS;
  - `SystemHypreBoomerAMG` overwrites two BoomerAMG settings at each
    `SetOperator`; `AmgXSolver::UpdateOperator` is a complete reuse, not a
    partial one.

### 2026-10-09 (night, end), compilation of the environment

- Second pass: everything installs (hypre, Ginkgo, AmgX, SuperLU_DIST,
  STRUMPACK, SUNDIALS, MUMPS, cuDSS, kokkos-kernels) except PETSc: its
  Kokkos files were compiled without OpenMP, which the headers of a Kokkos
  compiled with OpenMP require. Fixed by `petsc +openmp`.
- Third pass: PETSc installs; MFEM already fails at configuration,
  compiled in C++17 while Umpire 2026 (through camp) requires C++20. Fixed
  by `mfem cxxstd=20`.
- Kokkos was rebuilt (`+pic`, required by the shared kokkos-kernels): the
  Kokkos tests of TFEL pass (behaviours and material properties, `kokkos`
  and `kokkos-host`, 4/4), and so do the MGIS GPU headers with Kokkos.
- Tooling lesson: in the log monitors, the shell `grep` is a function
  (ugrep) that reported nothing; wait instead for the end of a process by
  its PID.

### 2026-10-09 (morning), MFEM 4.10, GPU-aware OpenMPI, monitored builds

- **The whole prototype moves to MFEM 4.10** (cuDSS, distributed Ginkgo,
  `GroupCommunicator` on the device), with the OpenMP and libCEED
  backends. The facts checked on 4.9 were checked again on 4.10.
- **OpenMPI with CUDA.** Its internal hwloc enabled the OpenCL plugin with
  an invalid path; the `openmpi` recipe of `mmm_gpu` builds it with
  `--disable-opencl`. The spack hwloc with CUDA would have rebuilt LLVM and
  TBB.
- **Robust builds.** The memory regulator suspends all the processes of
  the build except the largest one (which holds the memory and must
  finish); each phase of `install-mmm-gpu.sh` is restarted up to three
  times; an `mfem` recipe keeps the date of `config.mk` so that an
  interrupted build resumes without recompiling everything. A supervisor
  outside the session restarts the installation until it succeeds.
- **Two oversized MFEM 4.10 files with CUDA 13** (sm_120):
  `fem/lor/lor_batched.cpp` (`ptxas` at 23 GB) and
  `fem/integ/bilininteg_convection_ea.cpp` (`cicc` for more than two
  hours: 3D kernel specialized up to orders 8×9, unrolled loops). With two
  jobs, the Linux kernel killed `ptxas` for lack of memory despite the
  pause of the regulator; the phase restarts with one job, keeping the
  compiled objects.
- **Plan**: checked against all the code of mfem-mgis (82 features) and
  completed; fixed point techniques from the preliminary study added
  (fixed point on `K0`, ported MTest accelerators, modified Newton, frozen
  AMG of `K0`, single precision, prediction).

### 2026-10-09 (noon), MFEM 4.10 installed, work package L0 started

- **MFEM 4.10 installed** (12:32), with CUDA, OpenMP, libCEED, Ginkgo,
  AmgX, cuDSS, PETSc, SUNDIALS, SuperLU_DIST, STRUMPACK, MUMPS, SuiteSparse,
  Umpire and gslib. The example `ex1p` runs on CPU, on GPU with partial
  assembly, with BoomerAMG on GPU with two processes (CUDA-aware OpenMPI)
  and with the `ceed-cuda` backend.
- **Last build difficulties**:
  - `fem/integ/bilininteg_convection_ea.cpp`: cicc grew by 8 GB every 20
    minutes (17 GB after 2 h 30); patch that keeps the 3D specializations
    up to order 3 (the generic kernel takes over), 1 min 42;
  - link: `-lnvToolsExt` (absent from CUDA 13), requested by the AmgX
    options of MFEM, removed;
  - `libgs.a` of gslib not relocatable: recompiled with `-fPIC`.
  The current installation was finished by hand (same content as the
  recipes); a rebuild entirely through the recipes remains to be done
  (`deps/spack/README.md`).
- **Work package L0**: the three repositories build and pass their tests against
  MFEM 4.10 with CUDA (CPU device): mfem-mgis 128/128, mfem-mgis-examples
  20/20, mm-opera-hpc 6/6. Fixes:
  - mfem-mgis: `FindMFEM.cmake` truncated the values of `config.mk` that
    contain an equal sign; nvcc options translated for the C++ compiler;
    CUDA toolkit added; `tfel.cmake` reads several directories returned by
    `tfel-config` (TDLS); absolute tolerance of the hypre Krylov solvers
    (dead code and wrong call for FGMRES) with a test that detects the
    regression; option `--device` with a test;
  - bases of the examples and of mm-opera-hpc: reference repositories
    (`latug0/mfem-mgis-examples` master `8397eec`,
    `rprat-pro/mm-opera-hpc` main `8bd56eb`); the old bases, taken from
    older mirrors, no longer compiled with mfem-mgis;
  - mm-opera-hpc: the bubble loop, on host data, used `mfem::forall`,
    which since MFEM 4.10 requires a CUDA compilation when MFEM has CUDA
    (and would launch the loop on the GPU): `mfem::hypre_forall_cpu`.

### 2026-10-09 (afternoon), review of work package L0 and fixes

An independent review of the commits of work package L0 gave 25 findings. All
are handled:

- **Option `--device`**: the short name `-dev` was not recognized by
  `initialize` (the computation ran on the host without a message), like
  the short names of the PETSc options. The device is now created before
  the initialization of PETSc. Each process of a node takes a GPU of the
  node in turn. The number of GPUs is queried only for a `cuda` or `hip`
  configuration: otherwise MFEM stops when it has no GPU, or on a machine
  without a GPU.
- **Sources containing MFEM kernels**: installed CMake module
  `MFEMMGISDeviceCode.cmake`, also loaded when a project finds MFEM before
  mfem-mgis. nvcc and its host compiler are by default those of MFEM
  (`mpic++` of OpenMPI, which finds `mpi.h`). Architectures read in all
  the forms of the nvcc options. Explicit errors for HIP and when the CUDA
  language is not enabled. Shared CUDA runtime (`CUDA::cudart`, like the
  rest of the mfem-mgis code; MFEM embeds its own static copy of the
  runtime). C++ definitions and visibility options given to the CUDA
  sources. GPU tests outside the TFEL block, option `enable-gpu-tests`,
  GPU test with two processes.
- **TFEL**: only the `-I` and `-L` options of `tfel-config` are read;
  `env.sh` adds each library directory (mfem-mgis and mm-opera-hpc).
- **Documentation**: `AbsoluteTolerance` of the hypre Krylov solvers;
  32-bit integers of hypre and absence of `+gpu-aware-mpi`
  (`deps/spack/README.md`).
- **Prototype scripts**: LLVM and the other roots independent of MFEM
  are built with 4 jobs (and not 1 or 2); the memory regulator resumes the
  processes it suspended, even if they left the tree and even when it
  stops; the build directories are reconfigured at each call; convection
  patch limited to CUDA 13; PETSc recipe described as a full copy.
- **Bubble loop of mm-opera-hpc**: `hypre_forall_cpu` adds an OpenMP loop
  when hypre is compiled with OpenMP, which the old `mfem::forall` without
  a device did not do. No effect here (hypre without OpenMP), to be
  watched if hypre moves to OpenMP (oversubscription with several
  processes per node).

Full rebuild in a fresh directory: mfem-mgis 132/132,
mfem-mgis-examples 20/20, mm-opera-hpc 6/6, including the MFEM kernel on
the GPU with one and two processes.

### 2026-10-09 (afternoon), profiling and GPU

- MGIS: `mgis::setProfilingSynchronization` sets a function called at the
  beginning and at the end of each active profiling section, before the
  time is read (none by default). Without it, the duration of the kernels
  launched asynchronously on the GPU fell into the first following section
  that waits for them. Exact test of the number of calls.
- mfem-mgis: with a CUDA or HIP device, `initialize` installs a GPU
  synchronization there.
- ssna303 examples in 3D (ex4): option `--profiling`. Over one step, with
  two processes, the integration of the behaviours (in
  `NS::processNewUnknownsEstimate`) already takes 25 % of the time, the
  tangent and the solve 24 % each.
- Suites: MGIS 50/50, mfem-mgis 132/132, examples 20/20, mm-opera-hpc
  6/6.

### 2026-10-09 (afternoon), second review of the fixes of work package L0

A second independent review of the fixes confirmed 18 findings as fixed,
found 6 partially fixed and 1 not fixed, plus 9 new minor defects
(severity 1 or 2). All are handled:

- choice of the GPU of each process also for the `gpu` alias of MFEM (and
  `raja-gpu`, `ceed-gpu`...), which put all the processes on the first
  GPU; error if the MFEM device is already configured; the test checks the
  GPU of each process;
- CUDA toolkit taken next to the nvcc of MFEM; CMake 3.18 checked as soon
  as MFEM uses CUDA; warning if the CUDA compiler is not the one of MFEM;
  CUDA language checked in the directory of the target; at least C++20 for
  the CUDA sources; module documented in the installation guide; with a
  HIP MFEM, the configuration no longer fails;
- `env.sh`: one line per TFEL library directory, written at configuration
  time (zsh compatible);
- memory regulator: a suspended process that receives a signal (for
  example SIGTERM from spack) is resumed immediately so that it receives
  it (tried);
- `build-mfem-mgis.sh` removes from the cache the CUDA compilers of an old
  directory.

New directory `cases/`: campaign cases specific to the prototype. `beam`
(J2 plasticity, Newton of mfem-mgis); its tests compare the norm of the
final solution, identical for PCG and BoomerAMG, MUMPS and UMFPACK, and a
test checks that the comparison detects a wrong solution.
`cases/measure-reference.sh` measures the reference times of work package L0.

Suites: mfem-mgis 132/132, examples 20/20, mm-opera-hpc 6/6, cases 4/4.

### 2026-10-09 (afternoon), reference times, end of work package L0

Reference times of the current path on CPU for the cases of the campaign:
`results/reference-L0-2026-10-09/README.md` (beam, ssna303 in 3D,
polycrystal, thermomechanics). The linear solve (up to 61 %) and the
integration of the behaviours (up to 55 %) dominate; with MUMPS, the
prediction refactorizes the elastic operator at each step (40 % of the
time of the refined beam).

Comparison with the old stack (MFEM 4.9), same case and same Newton: PCG
and BoomerAMG a bit faster (15.4 s versus 17.5 s), MUMPS 1.55 times slower
(65.6 s versus 42.2 s), of which about 10 % is due to OpenMPI compiled
with CUDA (`--mca accelerator null`). The rest is still to be clarified
(MUMPS 5.9.1, OpenBLAS with OpenMP, MFEM 4.10); MUMPS with METIS is to be
tried when MFEM is rebuilt.

Work package L0 is finished: three repositories built and tested against MFEM
4.10 with CUDA, option `--device`, compilation of the MFEM kernels,
profiling sections synchronized with the GPU, reference times.

### 2026-10-09 (afternoon), work package L1: data at the points in the memory of MFEM

Steps done (mfem-mgis, patches 0016 to 0020):

- **Material arrays in `mfem::Vector` objects** (`MaterialArrays`,
  private base of `Material` declared before the MGIS manager): gradients,
  forces, internal variables and energies of both states, tangent blocks.
  MGIS uses them without copy (host location). MGIS no longer sets F to
  the identity when it is given the gradients: mfem-mgis does it, in both
  layouts.
- **Validity of the copies**: MGIS writes through raw pointers, which
  MFEM does not see. Each host access of mfem-mgis first goes through
  `synchronizeOnHost` (integration, end of step, revert, post-processings,
  accessors; read only access through `synchronizeOnHostForReading`).
  Without a GPU, nothing changes.
- **Elements of a material readable on the GPU**: arrays of the elements
  and of the offsets of their points in `PartialQuadratureSpace`.
- **Copies of the states by MFEM kernels** (`batch::MFEMFieldOperations`,
  first source of the library compiled as CUDA). The end of step and the
  revert go through `Material`, including for the point-wise models, which
  did not synchronize.
- **Component by component layout**: chosen when the material is built;
  the accessors based on a stride between points refuse it.

Postponed: material properties and external variables in user memory
(moot as long as MGIS receives host pointers), parameter
`QuadratureDataLayout` of the problem (no batch path yet), variant with
the arrays allocated by MGIS.

Tests: `MaterialArraysTest` and `MFEMFieldOperationsTest` on the `cpu`,
`debug` and `cuda` devices. Suites: mfem-mgis 138/138, examples 20/20,
mm-opera-hpc 6/6, cases 4/4. (Correction after review: this entry said
that each synchronization was checked by mutation; this was false for the
integration, the post-processings and the revert. See the next entry.)

### 2026-10-09 (afternoon), review of work package L1 and fixes

An independent review of work package L1 gave 27 findings, none of severity 3 or
more on the current path (nothing writes the arrays on the GPU yet), but
four that would become wrong results as soon as a kernel wrote them. All
are handled:

- **Missing host access points**: integration by PETSc (in the assembly
  of the residual, without going through `integrate`): synchronization at
  the beginning of `solve`; point-wise models (`computeNextState`); public
  functions that compute the mean forces and the resultants; evaluators of
  rotated stresses and gradients; loops of mm-opera-hpc (bubbles) and of
  the example ex8 that read the arrays directly.
- **Read-only post-processings**: they kept the GPU copies valid, instead
  of invalidating them at each step.
- **Complete validity rule**: a host access phase starts with the
  synchronization, and no access on the GPU takes place during this phase;
  the views returned by the accessors are only valid until then.
- **Private arrays**: no reallocation is possible any more behind the back
  of MGIS; access through pointers (`Read`, `Write`, `ReadWrite`).
- **End of step on the host** through `mgis::HostFieldOperations`, which
  keeps the dense copies (the generic kernels were 3.5 to 12 times
  slower); the operations by MFEM kernels remain for arrays on the GPU.
  With an MFEM compiled with HIP, their source is no longer compiled (the
  configuration failed).
- **CMake**: CUDA runtime carried by the library itself, check of the
  MGIS API at configuration time.
- **Types**: sizes as `std::size_t` (the "too large for MFEM" checks could
  not trigger), exceptions of the synchronizations reported in the context
  instead of stopping the program.
- **Component by component layout** also refused by the rotated
  evaluators.
- **Tests**: each entry point (accessors, integration, post-processings,
  revert, end of step) is called after all the arrays have been written on
  the GPU; end of step and revert in both layouts; six accessors refused;
  behaviour that computes an energy; material without any point; copies
  and fills in both directions. Checked by mutation: the synchronizations
  of the integration, the revert, the end of step, the tangent blocks and
  the accessors make the test fail when they are removed (`debug` device,
  and `cuda` when values are compared). Exception: the global
  synchronization before the post-processings duplicates the one that the
  functions they call already do; removing it makes nothing fail.

### 2026-10-09 (late afternoon), second review of work package L1

A second review of the fixes confirmed 22 findings as fixed, 5 as
partially fixed, and found 10 new minor defects (severity 1 or 2) and one
older omission. Handled:

- synchronization targeted on one array: the accessors, the evaluators
  and the computation of the mean forces now copy back only the array that
  is read (before: the eleven arrays, about 15 times too much to read the
  stresses);
- public linearized operators (`getLinearizedOperators`): synchronized
  when they are created;
- exceptions: `integrate` and the evaluators report them instead of
  stopping the program; the point-wise model synchronizes its status
  between processes even after an exception;
- the rotated evaluators use the common layout check;
- CMake: the MGIS API is checked by compiling a test (and not by the
  presence of a header); the installed configuration looks for the CUDA
  toolkit used at build time; with HIP, the header of the kernel
  operations is no longer installed;
- method names following the rules of mfem-mgis (`readWrite`, `write`);
- tests: point-wise model, rotation of the forces of an orthotropic
  behaviour (and its refusal with the component by component layout),
  mean forces called directly, four more accessors, energies, tangent
  blocks after the end of step, material without any point, distinct
  output files per device.

Mutations (each synchronization removed in turn): all make the test fail
under `debug`, and almost all under `cuda`. Two known exceptions: the
global synchronization before the post-processings (redundant) and the
one at the beginning of `solve`, which only serves the PETSc path (not
tested on the `debug` device).

Known limitation, as in the rest of mfem-mgis: an error local to one
process (here a synchronization exception, possible only with MFEM
compiled with exceptions) is reported by this process alone; the others
may wait in a collective operation.

### 2026-10-09 (evening), work package L2: behaviours integrated by batches on CPU

Parameter `BehaviourIntegration` of the nonlinear problems:
`PointByPoint` (default, current path) or `Batch`. The default value can
be changed by the environment variable
`MFEM_MGIS_BEHAVIOUR_INTEGRATION`, to run a test suite or an application
in this mode without touching the code.

By batches:

1. the gradients are computed element by element (code of
   `implementIntegrate` without the call of the behaviour), except for a
   prediction;
2. each integrator integrates its behaviour in a single MGIS call
   (`mgis::behaviour::integrate` on the data manager of the material):
   point by point inside MGIS for a `generic` library, as one batch for a
   `generic-parallel` library (backends `openmp`, `stlpar`...);
3. for an orthotropic behaviour, the gradients are rotated into the
   material frame before the call, and the tangent blocks into the global
   frame after it, at once for all the points (MGIS rotation functions on
   views; batch rotation function for a `generic-parallel` behaviour,
   which has no one-point rotation function).

Stresses in the global frame for the residual: rotated point by point, as
before, when the behaviour provides the one-point rotation; for an
orthotropic `generic-parallel` behaviour, rotated for all the points by
`integrateBatch`, right after the integration. The post-processings (mean
forces, rotated evaluators) compute them for all the points at each call.

Material properties at the beginning of the step: for the call, the
end-of-step ones are substituted for them as views (as in the point by
point path), then their own are put back in place, without copy.

Explicit refusals, with a message: `Batch` with PETSc, integrators
without batch support (micromorphic damage), `generic-parallel`
behaviour point by point.

**Independent review**: one serious defect found and reproduced. The
first version kept the rotated stresses in a cache of the material,
invalidated through versions of the arrays; MGIS writes through its views
without changing the versions, and the assembly of the residual with
PETSc (integration then residual, element by element) read stale stresses
for any orthotropic behaviour. Cache and versions removed (design
above); the sequence found by the review became a test. Also fixed: copy
of the material properties at each call (views), documentation, message
of the environment variable, 2D tests, orthotropic FBar, prediction.
Documented without change: MGIS gives the behaviour an initial time step
factor of 10 (1 point by point), without effect for mfem-mgis.

Tests:

- same results, to the bit, point by point and by batches: plasticity,
  finite strain, FBar, end-of-step material properties varying in space,
  orthotropic behaviours in small and finite strain (new test behaviour
  `OrthotropicSaintVenantKirchhoffElasticity`), uniform rotation or axes
  varying from one point to another, plane strain, prediction;
- `generic-parallel` libraries (`openmp`, `stlpar`) compared with the
  point by point computation with the one-point rotations: to the bit
  (1e-13 with prediction in finite strain), which checks the rotations of
  all the points independently;
- failure of the behaviour at some points;
- residual element by element right after the integration of each
  element;
- mfem-mgis suite: 139/139; in `Batch` mode (environment variable):
  133/139, the 6 failures are the expected refusals (micromorphic
  integrators); examples: 18/20 in `Batch` (refusal of PETSc);
  mm-opera-hpc: 6/6 in `Batch`;
- prototype cases: beam and polycrystal of mm-opera-hpc by batches with
  their `generic-parallel` behaviours, compared with the point by point
  computation.

**Two defects found by the measurements, outside mfem-mgis:**

- `scripts/env.sh` exports `OMP_PROC_BIND=spread OMP_PLACES=threads`
  (for Kokkos and AdaptiveCpp). libgomp, loaded by the libraries of the
  cases, then binds the main thread to the first processor: the `stlpar`
  backend (TBB) saw only one processor and stayed sequential, and
  `openmp` with 8 threads, spread over the hardware threads, used only 4
  cores. The measurement script removes these variables.
- The `g++` configuration files of the TFEL fork compiled the
  `generic-parallel` behaviours with `-O2` alone: at `-O2`, g++ 14
  vectorizes only with a very restrictive cost model. The polycrystal
  behaviour was 17 % slower than with the `generic` interface (78.4 s
  versus 67.1 s of integration over 2 steps). Fixed (TFEL patch 0093: the
  options of `tfel-config --oflags`, including `-ftree-vectorize`). The
  other configuration files (GPU) do not have `-DNDEBUG` either: to be
  checked in L6.

Measurements on a quiet machine (`results/L2-polycrystal-2026-10-09/`,
`results/L2-beam-2026-10-09/`): on the polycrystal, the behaviour takes
70 % of the computation on one process; by batches, it goes from 104.9 s
to 17.8 s with `openmp` with 8 threads and to 13.9 s with `stlpar` (full
computation: 150 s, then 59 s). With 8 single-thread processes, 14.5 s in
all modes. On the beam, the behaviour takes less than 5 %: the gain is
invisible.

### 2026-10-09 (night), work package L3 started: gradients and internal forces by kernels

First building blocks of the batch path without an element by element
pass on the host (`src/Batch/`, sources compiled as CUDA when MFEM is):

- **gradients** (`batch::computeUnknownsDerivatives`,
  `batch::computeGradients`): restriction to the elements (native
  ordering), physical derivatives at the points by the MFEM interpolator
  (`QuadratureInterpolator::PhysDerivatives`, layout `byVDIM`), once per
  integration rule for all the materials, then one kernel per material
  that writes the gradients of the behaviour (strain with √2, or
  deformation gradient) into the array of the material, in its layout,
  macroscopic gradients included;
- **internal forces** (`batch::computeElementsGeometry`,
  `batch::addInnerForces`): gradients of the shape functions on the
  reference element (`DofToQuad`), jacobians and determinants
  (`GeometricFactors`), one thread per element node, `Bᵀ σ w det J` with
  the same operator `B` as the generated integrators, written into a
  per-element vector with the layout of the MFEM restriction;
- interface of the integrators: `supportsKernels`,
  `computeGradientsByKernels`, `addInnerForcesByKernels` (standard CRTP
  base: 3D mechanics, plane strain, plane stress, small and finite strain,
  without FBar, integration rule common to all the elements of the
  material, mesh with a single geometry);
- provisional parameter `GradientsComputation` (`ElementByElement` by
  default, `Kernels`), with `BehaviourIntegration=Batch`, environment
  variable `MFEM_MGIS_GRADIENTS_COMPUTATION`.

`QuadratureInterpolator::MultTranspose` is not implemented in MFEM
4.10: the kernel of the internal forces is written here.

Tests: gradients and internal forces of the kernels equal to those of the
element by element computation up to 1e-16 (hexahedra and tetrahedra,
orders 1 and 2, small and finite strain, orthotropic, plane strain, two
materials with macroscopic gradients), on the host, on the `debug` device
of MFEM (which found a host read without synchronization in the test
itself) and on the GPU; full computations with kernels equal to the point
by point computation up to 1e-12. mfem-mgis suite 142/142; with all the
kernels enabled by the environment variables, only the expected refusals
fail (thermal, micromorphic damage).

Next step: the dedicated operator (residual without the MFEM loop over
the elements, tangent matrix free or by element matrices computed by
kernels), then the thermal problem and FBar.

### 2026-10-09 (night), work package L3: tangent by kernels, geometry, review

- **Geometry without nodes**: `Mesh::GetGeometricFactors` (and
  `QuadratureInterpolator::PhysDerivatives`, which calls it) gives nodes
  to the mesh (`EnsureNodes`). All the element transformations then go
  through them, and the loops still done element by element were 50 %
  slower (beam). The geometric factors of the kernels are computed on a
  copy of the vertices (module `Batch/ElementsGeometry`, kept per
  integration rule), and the gradients start from the reference
  derivatives, converted by the inverse of the jacobians in the kernel.
- **Element matrices of the tangent by kernels**
  (`TangentOperatorComputation`): one thread per pair of nodes of an
  element, `Σ Bᵀ K B w det J`, copied into the MFEM assembly. Equal to the
  element by element computation up to 1e-16 (CPU, `debug`, GPU). Slower
  on CPU (+24 % per call): they will only pay off with an assembly on the
  device, still to be written (`ElementRestriction::
  FillJAndData` of MFEM only handles scalar spaces), or matrix free.
- **Residual by kernels without the MFEM loop** when the form has no
  other integrator.
- **Independent review** of the kernels: MFEM conventions, operator B,
  layouts checked and correct; two serious defects reproduced and fixed:
  on the GPU, the Dirichlet conditions wrote the unknowns on the host
  without declaring it to MFEM (stale copy read by the gradient kernel,
  Newton failing); the regularized integrators (`Faltus2026`) were
  accepted and lost their regularization term without a message. Also
  fixed: refusals at `setup` (mixed mesh, unsupported space, attribute
  without an integrator), Newton exceptions reported instead of stopping
  the program, undeclared host accesses in the prediction (crash on the
  `debug` device), geometry cache without partial entries, options of the
  shared MFEM interpolator restored. New tests: refusal of `Faltus2026`,
  `BatchIntegrationTest` on `debug` and `cuda`, beam and polycrystal
  (periodic) in kernel mode on 2 processes.

Measurements (`results/L3-beam-2026-10-09/`, `results/L3-polycrystal-2026-10-09/`):
gradients by kernels, integration phase of the beam divided by 2.5 to
3.4; residual by kernels divided by 2.1 (beam) and 2.8 (polycrystal),
before the removal of the MFEM loop; element matrices by kernels slower
on CPU (+24 to +38 %). Full computations with all the kernels on the GPU
equal to the point by point computation (`BatchIntegrationTest-cuda`).

### 2026-10-09 (night), work package L3: direct assembly of the tangent, first computations on GPU

- **Direct assembly** (mfem-mgis b47f65c3cf): the `N x N` blocks of the
  element matrices computed by kernels are added directly into the local
  sparse matrix (position searched in the row, atomic addition, entries
  outside the sparsity pattern reported), without storing element
  matrices. The pattern comes from the first MFEM assembly; the end of
  `GetGradient` is the one of MFEM (diagonal block and `PtAP` in
  parallel, `RAP` in sequential).
- **Elimination of the imposed dofs**: `HypreParMatrix::EliminateRowsCols`
  copies the matrix back to the host and computes the eliminated part,
  useless here (122 ms per call on GPU for 111,843 unknowns). Replaced by
  `EliminateBC` (`DIAG_ONE`), in place and on the GPU, giving the same
  matrix since the pattern is symmetric.
- **First full computations on GPU** (beam r4, 3 steps, 1 process,
  `--device cuda`): 7.8 s with all the kernels, versus 24.5 s on one CPU
  core and 29.1 s on GPU without kernels. The assembly kernel takes 7 ms
  per call (0.48 s on one core). hypre is compiled for the GPU (solve 0.9
  s instead of 13.4 s). Still on the host: the prediction (linearized
  operators assembled by MFEM, 4.9 s), the behaviour (`generic`) and the
  first call of the tangent (MFEM loop for the pattern, 0.9 s).
- **On CPU**, the direct assembly gains nothing: one thread per pair of
  nodes recomputes the inverse of the jacobian and rereads the tangent for
  each pair (0.48 s per call, versus 0.34 s for the MFEM loop).
- **Second review** of the kernels: severity 5 defects checked as fixed.
  Handled: failure of a kernel on a single process made collective before
  the following collective operations (otherwise deadlock), options of the
  interpolator restored even on exception, force vector kept between two
  residuals, Dirichlet conditions without dofs on a process without a
  useless transfer. One finding of the review was wrong: with the `debug`
  device of MFEM, any access on the device, even read-only, protects the
  host copy (`MemoryManager::GetDevicePtr`); the host reread of the
  imposed dofs remains necessary (crash reproduced without it). The host
  rereads of the unknowns, on the other hand, were useless (checked on
  `debug` and `cuda`, sequential and 2 processes) and are removed.
- **Test pitfall**: the test executables of mfem-mgis are outside `all`
  (target `build-tests`), and the library is static; without rebuilding
  `build-tests`, ctest reruns the old binaries. `scripts/env.sh` must also
  be sourced (the `LD_LIBRARY_PATH` of the profile loads another TFEL).

Measurements: `results/L3-beam-assembly-2026-10-09/`. Next: pattern
built from the dofs of the elements (first call without the MFEM loop,
sorted columns), variant of the kernel adapted to the CPU (one thread per
element), MPI test of the kernels in mfem-mgis.

### 2026-10-09 (night), work package L3: tangent pattern built from the elements

- mfem-mgis 799941c8fd: the local matrix is created from the dofs of the
  elements at the first call of `GetGradient` (no more MFEM loop). First
  attempt with sorted columns: wrong in parallel. `HypreParMatrix`
  (diagonal block constructor) shares the arrays of the local matrix and
  puts the diagonal term at the head of each row, in place, without
  touching the `isSorted` flag. The pattern therefore follows the hypre
  order (diagonal, then sorted columns), which this reordering leaves
  intact, and the search in the kernel accepts any order (first column,
  bisection, then a scan).
- The defect was only seen by the prototype cases and by the suite in
  kernel mode (tests `...P`): no mfem-mgis test on the kernels runs in
  parallel. MPI test to be added.
- Tangent of the beam (3 steps): 16 ms per call on GPU, 0.40 s on one CPU
  core versus 0.45 s for the MFEM loop. On GPU, the prediction now takes
  73 % of the time.

### 2026-10-09 (evening), minimal GPU slice: behaviour integrated on the GPU

- MGIS a0b5140: `make_unchecked_view`, view on a manager that does not
  check that the values it keeps on the host (uniform values, local
  storage) are accessible to the computations of the behaviour; the
  caller replaces the views of the arrays and copies these values.
  `make_view` unchanged (refusal). CUDA test extended.
- mfem-mgis ee1f31fa1c: a `generic-parallel` behaviour is integrated on
  the device copies of the material arrays if its computations can access
  them (`cuda` behaviour on GPU; any behaviour on the `debug` device,
  whose copies are in host memory). MGIS view moved onto the copies made
  valid by MFEM (`Read`, `ReadWrite`, `Write`), properties, external
  variables and mass density copied to the device at each call, rotations
  of the orthotropic behaviours on the device, as well as the change to
  the global frame when the behaviour cannot access the host. No more
  synchronization on the host before the integration when everything is on
  the device; read synchronization after the integration if the forces or
  the tangent stay element by element.
- Validation: `openmp` and `stlpar` behaviours on the `debug` device
  (accesses to invalid copies detected), identical to point by point to
  the last bit without kernels; `cuda` behaviours on GPU, sequential and
  with 2 processes, with and without kernels, orthotropic ones included.
  Case `beam-gpu`.
- Beam r4, 3 steps, GPU: 6.25 s. Integration 0.24 s (0.47 s with the
  `generic` behaviour on the host). The prediction takes 77 % of the time:
  next work package (work package L4).

### 2026-10-09 (evening), prediction by kernels, beam in 4 s on GPU

- mfem-mgis 7159b6de4c: with the forces and the tangent by kernels, the
  operator of the prediction is assembled by kernels in a local matrix
  kept from one prediction to the next, and so are the internal forces at
  the beginning of the step: no more element by element bilinear forms.
  System formed like `FormLinearSystem` of MFEM, elimination in place
  (`B = Pᵀ b - K X`, then `B[e] = X[e]`). The boundary conditions that
  contribute to the right-hand side (pressure) are still assembled by
  MFEM; fallback to the old path if they contribute to the operator.
- Defect found along the way: the list of the imposed dofs of the
  prediction wrapped the memory of a `std::vector`. Read on the `debug`
  device, it made the test crash once in a few runs (depending on the
  memory layout). Copy allocated by MFEM.
- Beam r4, 10 steps (`results/L6-beam-2026-10-09/`): **4.06 s on one
  GPU** against 17.54 s on 8 cores (point by point) and 14.68 s on 8
  cores with the same kernels. What remains is the hypre solve (61 %) and
  the prediction (21 %, including its solve). Two processes on the same
  GPU: 28.6 s (solve 9 times slower, not analysed).

### 2026-10-09 (night), independent review of work packages L3, L4 (prediction) and L6 (GPU slice)

Adversarial review of mfem-mgis b47f65c3cf to 7159b6de4c and of MGIS
a0b5140 (complete reading, compared with MFEM 4.10, comparison programs on
cpu, `debug` and `cuda`, 1 to 3 processes including one without any
element, order 2 elements, `byVDIM` ordering). **No defect of severity 3
to 5.** Checked: prediction by kernels, elimination in place, pattern in
the order of hypre, and synchronisations after the integration on the
device.

Handled (mfem-mgis 9398a9e33c, decb6fd2c1):
- tests (severity 2): the prediction of each step is recorded and
  compared; case with imposed pressure (the internal forces at the
  beginning of the step count; forgetting them in the prediction by
  kernels, which the tests did not see, is now detected); two materials,
  one integrated on the device, the other on the host (new mesh
  `cube-two-materials.mesh`: `cube_2mat_per.mesh` is periodic, its faces
  x = 0 and x = 1 are identified and a Dirichlet condition loads nothing
  there); behaviour with properties varying from one point to another in
  `generic-parallel`; kernels with prediction for these behaviours;
- documentation that had become wrong (computation in the global frame on
  the device, comments of `solve` and of `MaterialArrays`);
- useless copy to the host of the rotated forces at each residual;
  possible leak of the prediction matrix; loading function of the
  Dirichlet conditions evaluated on all processes; integer overflow
  refused; GPU tests never run simultaneously.

Remaining (severity 1, noted): `update` and `revert` on the host (round
trip of the arrays at each step; MGIS has the CUDA operations to do them
on the GPU); synchronisation of all the arrays when only the gradients are
computed on the host; the prediction matrix duplicates the tangent; each
process reserves about 5 GB on the GPU (Umpire pool of hypre), and hypre
would be compiled without GPU-aware MPI, which could explain the slowness
of two processes on one GPU.

### 2026-10-09 (night), polycrystal on GPU

Polycrystal of mm-opera-hpc (single crystal behaviour of UO2 compiled for
`cuda`, 10 steps, `results/L6-polycrystal-2026-10-09/`): **13.18 s on one
GPU against 31.65 s on 8 cores** (151 s on one core), same macroscopic
stresses to within 4e-13 of the stress scale. The behaviour, which
dominates on CPU, takes 5.7 s on the GPU, as much as the solve. With the
behaviour on the host, the GPU gains nothing (107 s). Two processes on the
same GPU: solve 14 times slower, as for the beam.

### 2026-10-09 (night), end of step on the GPU

mfem-mgis (commit "Update and revert on the device"): the states of a
material whose behaviour is integrated on the device are updated and
restored on the device (`update` and `revert` of MGIS on a view of the
copies, through `MFEMFieldOperations`), and the synchronisation on the
host at the beginning of each solve is no longer done when all the
behaviours are integrated on the device (except with PETSc). The arrays
stay on the GPU from one step to the next. Beam r4, 10 steps, GPU: 3.3 s
instead of 4.1 s (8 cores: 17.5 s).

### 2026-10-09 (night), Nsight Systems profile of the beam on GPU

No notable hidden copy between the host and the GPU (146 MB to the GPU
over 3 steps, mostly setup data; one or two vectors per iteration). The
GPU computes 0.71 s out of 1.8 s, in 95,000 kernel launches, 65 % of which
are matrix-vector products of hypre: the solve is limited by the latency
of the launches. Details in `results/L6-beam-2026-10-09/README.md`.

### 2026-10-10 (night), action of the tangent without assembly

- mfem-mgis cbe693a9a8: action of the tangent `Σ Bᵀ K B x w det J` and its
  diagonal by kernels, without assembly (increments of the gradients by
  the gradients kernel, product by the blocks of the tangent at each
  point, internal forces kernel; diagonal by the elementary matrices
  kernel restricted to the pairs of a node with itself). Equal to the
  product by the assembled matrix and to its diagonal to within 1e-16
  (all the cases of `GradientsKernelsTest`, cpu, `debug`, `cuda`). Common
  checks of the computations by kernels factored out (`getKernelsData`).
- mfem-mgis 4d465460c3: `TangentOperatorComputation=MatrixFree`, the
  tangent returned by `GetGradient` is a matrix free operator
  (`MatrixFreeTangentOperator`, imposed dofs by `ConstrainedOperator`),
  with the preconditioner `OperatorJacobiSmoother` (diagonal by
  `AssembleDiagonal`, or `GetDiag` for a `SparseMatrix` of MFEM, which
  does not have `AssembleDiagonal`). Complete computations equal to point
  by point (small and finite strains, orthotropic, two materials,
  `generic-parallel` behaviours, `debug`, GPU, 2 processes). Case
  `beam-matrix-free`.
- Beam r4 on GPU: 95 s with CG + Jacobi against 3.4 s with the assembled
  tangent and BoomerAMG. Correct, but a suitable preconditioner is needed
  (AMG kept over several iterations, polynomial multigrid).

### 2026-10-10 (night), preconditioner on the assembled tangent

- mfem-mgis b129009b0c: preconditioner `AssembledOperator` for the matrix
  free tangent. It asks the operator for its assembled version (tangent
  assembled by kernels, imposed dofs eliminated) and builds on it the
  preconditioner given as option (BoomerAMG, Jacobi...), rebuilt every
  `UpdatePeriod` Newton iterations and kept in between. Solver `cg-amg`
  of the beam case.

### 2026-10-10 (night), fixes from the second review of work packages L3 and L6

Independent review of cbe693a9a8 and 4d465460c3
(`scratch/night-2026-10-09/review-L3-L6-suite.md`). Fixed in
mfem-mgis 1bbcb877d8 (and MGIS 860b693):

- F1: `update` and `revert` of a behaviour integrated on the device no
  longer carried the external variables of s1 over to s0 (temperature of
  s0 frozen: thermal expansion counted again at each step, without any
  message). New MGIS function `updateInputValues` (material properties,
  external variables and mass density, same rules as `updateValues`,
  without the arrays), called on the host after the copy of the arrays to
  the device.
- F2: after `update` on the device, the host copies of the arrays were
  stale. They are now synchronised for reading at the end of `update` and
  of `revert` (one copy from the GPU to the host per step).
- F3: an exception of the tangent, of its diagonal or of its action in
  `computeNewtonCorrection` (`noexcept`) ended in `std::terminate`. It is
  recorded in the context of the solve, which fails cleanly.
- F4: diagonal on a non-conforming mesh: `|Pᵀ| diag(A)`, like MFEM.
- F5, F6, A2, Q2: copy of the list of imposed dofs of the constrained
  operator (rebuilt when it changes), list of the Jacobi smoother kept as
  a member (smoother reused), guard on the number of threads per node for
  the diagonal, shared functions instead of copied ones, parameter
  `TangentOperatorComputation` read once (enumeration), documentation of
  the synchronisations.
- Tests (T1) in `BatchIntegrationTest`: behaviour with thermal expansion
  and temperature varying in time, states read back on the host without
  synchronisation after each `update` and `revert`, diagonal of the
  matrix free operator compared to that of the assembled tangent, clean
  refusal of MatrixFree with a boundary integrator. Each test was checked
  by mutation: without the call to `updateInputValues`, the
  "temperature" cases fail on `debug`; without the synchronisation,
  `debug` stops at the read; with a diagonal falsified by 1e-6, all the
  MatrixFree cases fail; without the `try/catch`, `std::terminate`.
- Suites: 146/146 (including 7 on GPU), forced kernels 129/139 (10
  expected refusals: heat, micromorphic), cases 19/19.

### 2026-10-10 (night), kernels of the derivatives and fused action

- mfem-mgis a2602d117a: the derivatives of the unknowns at the points are
  computed by a dedicated kernel (gradients of the shape functions of
  `mfem::DofToQuad`) instead of `QuadratureInterpolator::Derivatives` of
  MFEM, whose generic version for simplices cost 890 µs per call on the
  beam (39 µs now). The action of the matrix free tangent is computed by
  a single kernel per material, one thread per element (gradient of the
  increment, product by K, `Bᵀ`, no array at the points): 128 µs per
  material instead of 350 µs. GPU compute time of the matrix free beam
  (10 steps): 2.88 s instead of 5.15 s.
- The fused kernel reads the blocks of K point by point (non coalesced
  accesses, about 45 % of the bandwidth): the component by component
  layout (parameter `QuadratureDataLayout`, planned in L1) remains to be
  plugged into the integrators.

### 2026-10-10 (night), fixed point solver (work package L5)

mfem-mgis f3d8ed1e55: non-linear solver `FixedPoint`, taken from the
prototype of 04/10 (`mfem-mgis-gpu`) and written in mfem-mgis:

- methods `ElasticOperator` (K0 assembled once per solve or once for the
  whole computation, solved by the solver of the problem or by a
  dedicated solver whose setup is kept) and `ModifiedNewton`;
- accelerations in an extensible factory (`AccelerationAlgorithmFactory`):
  `Anderson` (type II, matrix of the scalar products updated along the
  iterations, two scalar products per iteration and per kept difference)
  and `Cast3M` (like MTest);
- `HypreBoomerAMG` also declared as a linear solver (one AMG cycle per
  iteration).

Cube in tension (plasticity with linear hardening), 2 steps: Newton 8
iterations, fixed point 103 (observed rate 0.642, bound 3μ/(3μ+H) = 0.887),
Anderson 37, Cast3M 56, modified Newton 59, one AMG cycle + Anderson 46.
Same solutions as Newton to within 1e-12, on CPU, `debug` and `cuda`, with
1 and 2 processes. On the GPU, a BoomerAMG cycle is much weaker (101
iterations instead of 46): the settings of hypre on GPU (l1-Jacobi
smoother, PMIS coarsening) are to be exposed.

Independent review of b129009b0c and 1bbcb877d8 finished
(`scratch/night-2026-10-10/review-1bbcb877d8.md`): F1 to F6 confirmed
fixed; new findings C3 (direct write after `update` ignored on the
device), C7 (exception on a single process: deadlock), C4 and C5 (matrix
of the preconditioner `AssembledOperator`), C9 (tests), C8 (copy of K).
Fixes in progress.

### 2026-10-10 (night), fixes from the review of 1bbcb877d8

mfem-mgis bae47cb272:

- C3: `update` and `revert` of a behaviour integrated on the device end
  with the synchronisation of the state arrays on the host (device copies
  invalidated), as on the host: a direct write into the state after
  `update` is taken into account. C8: the blocks of K, which are zero, are
  set to zero on the host instead of being copied.
- C7: a failure of the assembly of the tangent or of the setup of the
  linear solver on one process is made collective before the solve;
  during the solve, only the failures reported on all processes
  (`CollectiveFailure`, raised by `raiseOnFailureOnAnyProcess`) are
  caught in parallel.
- C4 and C5: the preconditioner `AssembledOperator` knows the number of
  assemblies of the tangent (counter shared by the problem and its
  providers): it sets itself up again if its matrix was reassembled by
  another call, and keeps a valid setup for the solve of the prediction
  instead of resetting the period to zero. Beam r3, 4 processes, 10
  steps, 25 Newton iterations: 25, 13 and 7 setups of BoomerAMG for a
  period of 1, 2 and 4 (before: one more per step, whatever the period).
- Tests: direct write after `update` (detects C3 by mutation on `debug`);
  diagonal of the matrix free operator also checked with the assembled
  preconditioner, and case without prediction (detects C4 by mutation:
  stop in `MPI_Comm_size` with 2 processes).
- Suites: 141/141, GPU 7/7, forced kernels 131/141 (10 expected
  refusals), cases 19/19.
- Remaining (severity 1): C11 (in sequential, AssembledOperator only
  wraps Jacobi), C12 (preconditioner kept after a Newton failure), C14
  (one MPI reduction per application of the matrix free operator), tests
  of F4, F5, F6.

### 2026-10-10 (night), accelerations of MTest

mfem-mgis: port of all the accelerations of MTest into the factory of the
fixed point: `UAnderson` and `FAnderson` (Anderson method on the residuals
of the fixed point or on those of the problem, acceleration every two
iteration steps by default), `Secant`, `CrossedSecant`, `AlternateSecant`,
`IronsTuck`, `Steffensen`, `CrossedDelta2`, `Crossed2Delta`,
`Crossed2Deltabis`, `AlternateDelta2`, `Alternate2Delta`. Cube in tension
(2 steps, Newton 8 iterations, plain fixed point 103): Anderson 37,
UAnderson and FAnderson 38, CrossedSecant 41, Alternate2Delta 46, IronsTuck
48, Cast3M 56, AlternateDelta2 60, AlternateSecant 61, Secant 65,
Crossed2Delta 70, CrossedDelta2 140, Crossed2Deltabis 374; Steffensen,
applied component by component, diverges on this coupled problem (same
formula as MTest). All the others give the Newton solution to within
1e-12.

### 2026-10-10 (night), BoomerAMG options, prediction by K0, review of L5

- mfem-mgis d05475d265: options of BoomerAMG (`StrengthThreshold`,
  `CoarseningType`, `AggressiveCoarseningLevels`, `InterpolationType`,
  `RelaxationType`, `RelaxationSweeps`, `CycleType`, `NumberOfCycles`,
  `MaximumNumberOfLevels`), applied again after each `SetOperator` (which
  keeps only part of the settings of hypre). On the GPU, the Chebyshev
  smoother (`RelaxationType=16`) makes a BoomerAMG cycle much more
  efficient for the fixed point (438 iterations instead of 1,231).
- L5 measurements (`results/L5-beam-2026-10-10-0250/`): on 8 cores, the
  fixed point (Anderson 8, one AMG cycle, K0 kept) beats Newton (13.3 s
  against 20.3 s on a loaded machine; 9.4 s against 16.6 s on a quiet
  machine, ×1.7); on the GPU, a tie (3.4 s) with Chebyshev.
- mfem-mgis 846551e0cb: `getElasticPredictionLinearSolver`: the fixed
  point can provide the solver of the elastic prediction (conjugate
  gradient preconditioned by the AMG of K0 already built); modest gain
  (prediction 1.34 s instead of 1.52 s on the GPU), the prediction remains
  the main cost item of the fixed point.
- Independent review of a2602d117a, f3d8ed1e55, bae47cb272
  (`scratch/night-2026-10-10/review-L5.md`): kernels and fixed point
  correct (42 configurations per device, Anderson identical to a direct
  implementation, Cast3M faithful to MTest); C7 not fixed (the collective
  agreement could pair with another reduction). Fixed in 1d5adf1a32: in
  parallel, only the `CollectiveFailure` exceptions are handled, any other
  exception stops the program; options of the fixed point validated
  before being applied, actions called once, NaN guard, information of
  the linear solver, bounds on the numbers of iterations in the test (a
  wrong Cast3M is detected), kernels and `debug` variants.

### 2026-10-10 (night), NewtonK0, internal forces, short review

- mfem-mgis 3e414bc352: internal forces with one thread per element (81 µs
  instead of 217 µs per call on the GPU).
- mfem-mgis 1ef05104d1: profiling sections of the prediction
  (`PRED::integrate`, `forms`, `assemble`, `solve`); the one-off setups
  (geometry 214 ms, pattern) skewed the averages per call.
- mfem-mgis 30beae4764: method `NewtonK0` of the fixed point (Newton,
  conjugate gradient preconditioned by BoomerAMG of K0 built once); fixes
  from the short review of 846551e0cb and 1d5adf1a32
  (`scratch/night-2026-10-10/review-L5-suite.md`): prediction solver of
  the fixed point limited to the elastic prediction at the beginning of
  the step (N1), message then `MPI_Abort` for a non-collective exception
  (N4), K0 and its solver replaced without intermediate state (N2), test
  of the prediction solver with a linear solver that counts its solves
  (mutation detected).
- Measurements (`results/L5-beam-2026-10-10-0250/README.md`): on the GPU,
  NewtonK0 2.6 s against 3.16 s for Newton (stable); on 8 cores,
  accelerated fixed point 9.5 s, NewtonK0 11.6 s, Newton 13.5 s.

### 2026-10-10 (night), NewtonK0 on the polycrystal, review K0

- mfem-mgis 9249013150: `PeriodicNonLinearEvolutionProblem` takes the
  parameters of the problem (for example `NonLinearSolver`). mm-opera-hpc
  6fdd03c then 155202f: option `--non-linear-solver` of the polycrystal
  (`Newton` or `NewtonK0`, K0 obtained by an integration of the
  behaviour, which does not provide a prediction operator).
- Measurements (`results/L5-polycrystal-2026-10-10-0340/README.md`, script
  `cases/measure-L5-polycrystal.sh`): on 8 cores, NewtonK0 24.9 s against
  29.6 s (×1.19, solve -36 %); on the GPU, 11.3 s against 11.9 s (×1.06),
  the behaviour dominates. Same Newton iterations, stresses equal to
  within 2e-13.
- Short review of 30beae4764 (`scratch/night-2026-10-10/review-K0.md`,
  36 NewtonK0 configurations in sequential and on 2 processes, within
  5e-11 of Newton): nothing above severity 2. Fixed in 9355256a9f:
  `TangentLinearSolver` and `PredictionLinearSolver` name an iterative
  solver of MFEM from the factory (conjugate gradient by default; GMRES,
  FGMRES, BiCGSTAB, MINRES, SLI), the conjugate gradient requires a
  symmetric positive definite tangent (documented); `None` for
  `Acceleration`; `FGMRESSolver` in the factory; Anderson regularisation
  bounded at 1e-14; tests: each iterative solver, tangent without
  assembly, options changed between steps, prediction with constant
  gradients; six mutations detected
  (`scratch/night-2026-10-10/mutations-K0/`). Option `-tls` of the beam
  (iterative solver of NewtonK0).

### 2026-10-10 (night), inexact Newton, review K0b

- Review of 9355256a9f (`scratch/night-2026-10-10/review-K0b.md`):
  nothing above severity 2. Fixed in 82087f9f83: sequential prediction
  without kernels with the diagonal set to one at the imposed dofs (an
  iterative solver preconditioned by K0 solves it in one iteration: GMRES
  took 2 21 fixed point iterations, 2 2 now), base options applied only if
  all are valid, messages, tests (replacement of the solver of K0 between
  two steps with a solver that counts the calls after its destruction,
  default solver and tolerance of NewtonK0, prediction by GMRES); all the
  mutations detected.
- mfem-mgis 72f2e6e29b and 517c0c0b08: inexact Newton, option
  `AdaptiveLinearSolverTolerance` of the Newton and of the method
  `NewtonK0` (choices 1 and 2 of Eisenstat and Walker, settings of MFEM,
  Kelley safeguard against oversolving, tolerance of the linear solver
  restored after each solve, norm of the linear residual in the iteration
  information); unit test `AdaptiveLinearSolverToleranceTest`.
- Measurements (`results/L5-beam-2026-10-10-0429/README.md`, script
  `cases/measure-L5-ew.sh`, options `-alst` and `-alsti` of the beam): on
  the GPU, NewtonK0 with Eisenstat-Walker 1 and an initial tolerance of
  0.1 takes 1.96 s against 3.12 s for Newton (×1.59) and 2.57 s for
  NewtonK0 with a fixed tolerance; on 8 cores, 9.46 s, at the level of the
  best fixed point. With the classical Newton, Eisenstat-Walker does not
  pay off (each additional iteration rebuilds BoomerAMG).
- Review of the inexact Newton (`scratch/night-2026-10-10/review-EW.md`):
  nothing above severity 3 (coverage). Fixed in 4b72095ee0: safeguard of
  choice 1 of Eisenstat and Walker as in the paper and in PETSc, initial
  tolerance bounded by the maximal one, tests that recompute the sequence
  of tolerances from the iteration information (ten mutations detected),
  glossary section.
- mfem-mgis b4127a6c12: Chebyshev smoother of MFEM
  (`OperatorChebyshevSmoother`, orders 1 to 5) in the factory, test
  `OperatorSmoothersTest` (Laplacian: 46 CG iterations with Jacobi, 27
  and 17 with Chebyshev of order 2 and 4). Matrix free measurements on the
  GPU (`results/L5-beam-2026-10-10-0456/README.md`, option
  `-ls cg-chebyshev` of the beam): Chebyshev slower than Jacobi (each
  order adds a matrix free action); NewtonK0 + Eisenstat-Walker 1 matrix
  free within 3 % of the assembled tangent (2.05 s against 1.99 s).
- Polycrystal (mm-opera-hpc, options `--adaptive-linear-solver-tolerance`
  and `--adaptive-linear-solver-initial-tolerance`), measurements
  `results/L5-polycrystal-2026-10-10-0502/README.md`: NewtonK0 with
  Eisenstat-Walker 1 is the best on the GPU (9.58 s against 11.88 s for
  Newton, ×1.24) and on 8 cores (23.46 s against 29.60 s); the solve
  becomes negligible on the GPU and the behaviour dominates (7.1 s out of
  9.6 s).

### 2026-10-10 (morning), heat transfer by kernels, review C

- Review of b4127a6c12, 4b72095ee0 and of the polycrystal
  (`scratch/night-2026-10-10/review-C.md`): nothing above severity 2.
  Fixed: mfem-mgis 6cae30af99 (invariance of the smoothers under a
  diagonal scaling, which detects an ignored diagonal or a frozen
  smoother; choice 1 without linear residual; glossary), mm-opera-hpc
  976a60b (errors of the options), README of the measurements (memory not
  reduced by the matrix free mode, rounded figures).
- mfem-mgis 0a26a38be0: non-linear stationary heat transfer by kernels
  (work package L7.3): 2D and 3D temperature gradient kinds, number of
  components of the unknowns separated from the dimension in all the
  kernels, block of the derivative of the flux with respect to the
  temperature (matrix, diagonal, matrix free action), values of the shape
  functions, temperature at the points copied into the external variable
  (new method of the generated thermal integrators and of the generator).
  Test `HeatTransferKernelsTest`: 3D and plane strain, assembled and
  matrix free, with and without prediction, generic and generic-parallel
  behaviours, sequential, 2 processes, `debug` and GPU, differences to
  point by point from 1e-16 to 3e-14; three mutations detected. Suite with
  forced kernels: only the 6 micromorphic tests are still refused (10
  before). The thermal tangent is negative definite: the conjugate
  gradient does not apply matrix free (GMRES in the test). Independent
  review in progress.

### 2026-10-10 (morning), thermal measurements, Jacobians by kernel

- Measurements of the stationary heat transfer by kernels (`cases/thermal`,
  r5, 279,873 unknowns, script `cases/measure-L7-thermal.sh`,
  `results/L7-thermal-2026-10-10-0542/README.md`): on the GPU, the first
  integration cost 1.76 s, one-time computation of the geometry on the
  host (first order finite element space: 1.0 s; projection of the
  coordinates: 0.83 s; geometric factors of MFEM on the GPU: 0.02 s).
- mfem-mgis 8b0b47acfb: Jacobians of the meshes without nodes computed by
  a kernel from the vertices of the elements and the finite element of the
  transformations of MFEM (`GetTransformationFEforElementType`). Test
  `ElementsGeometryTest`: segments, triangles, quadrilaterals, tetrahedra,
  hexahedra, prisms, pyramids, regular and distorted, unknowns of order
  one and two, curved mesh, host, `debug` and GPU, differences to the
  transformations of MFEM at most 2e-15; six mutations detected. Heat
  transfer on the GPU: 3.01 s down to 1.25 s, ×3.8 against 8 processes
  with the kernels (`results/L7-thermal-2026-10-10-0606/README.md`). Beam
  r4: NewtonK0 + Eisenstat-Walker 1 from 2.14 s to 1.90 s on the GPU
  (`results/L5-beam-2026-10-10-0610/README.md`); polycrystal unchanged
  (`results/L5-polycrystal-2026-10-10-0612/README.md`). Independent
  review done (`scratch/night-2026-10-10/review-G.md`, nothing above
  severity 3).
- mfem-mgis fbac8c024a: fixes from review H of the heat transfer by
  kernels. The test only saw the x component (one-dimensional solution):
  third imposed temperature on y = 1. New cases: plane stress, rotated
  anisotropic behaviour (3D, plane strain and plane stress), two
  materials. The anisotropic test behaviour has a non-symmetric
  conductivity: the transposition of the block of the derivative of the
  flux with respect to the gradient, which neither this test nor
  BatchIntegrationTest saw (all tangents symmetric), is now detected.
  Thermal behaviours recognised by their structure and not by the names of
  the glossary; checks of the copy of the temperature. Ten mutations
  detected.
- mfem-mgis e56f87ddd9: fixes from review G. G1 (severity 3): the finite
  element of the transformations of MFEM is a global object whose
  `GetDofToQuad` cache is indexed by the address of the rule; two rules
  allocated one after the other at the same address gave wrong Jacobians.
  The gradients of the shape functions are now computed for each rule.
  Jacobians always owned by the geometry (including for a mesh with
  nodes). Test: inverted elements (negative determinants), parallel
  meshes refined in a non-conforming way (variant on two processes),
  rules at the same address, refusals. Eight mutations detected.
- Review HG2 of fixes H and G (`scratch/night-2026-10-10/review-HG2.md`):
  nothing above severity 2, to be handled (comparison of the temperatures
  at the points that compares nothing, cache by rule address still there
  for the meshes with nodes and the DofToQuad of the unknowns).
- mfem-mgis 62908a890d: stiffness matrices by kernels, one iteration per
  element on the host (inverse of the Jacobian, operators B and products
  K B once per point); heat transfer r5 on 8 processes: kernel from 100 to
  66 ms per tangent, the search of the columns of the sparse matrix now
  dominates (40 ms against 12 ms of computation), PtAP of hypre 37 ms.
- mfem-mgis a394212d92: transient term of the heat transfer by kernels
  (kinds `TEMPERATURE_3D/2D`, gradient = temperature at the point,
  operator B = value of the shape function). `TransientHeatTransfer` was
  neither tested nor declared in 2D: test behaviour HeatStorage,
  integrator declared in plane strain and plane stress, transient cases in
  HeatTransferKernelsTest (differences to point by point at most 2e-15,
  six mutations detected). Open question: the rule of
  `TransientHeatTransfer` is of order 2 OrderGrad, that is a single point
  for first order tetrahedra, which under-integrates the mass term.
- mfem-mgis 3e3cc5b0f2: fixes from review HG2. Jacobians of all the meshes
  by the kernel, from the nodes (any order, ordering by nodes or by
  components, discontinuous nodes of periodic meshes) or from the
  vertices; the geometric factors of MFEM are no longer used (cache by
  rule address, rules assumed to be tensor products: 57 % error
  otherwise). Check of the shape functions cached by rule address. The
  thermal test compared the initial temperatures (storage mode by copy):
  external mode and non-uniform reference temperature. Seven mutations
  detected. Polycrystal unchanged
  (`results/L5-polycrystal-2026-10-10-0752/README.md`, CPU tangent from
  2.01 to 1.60 s).
- mfem-mgis 9b93bc8c36: micromorphic by kernels (work package L7.8). The three
  micromorphic integrators move to the CRTP base (point by point results
  identical to the last digit), kinds `MICROMORPHIC_DAMAGE_3D/2D` (value
  and gradient of the unknown, block diagonal tangent).
  MicromorphicKernelsTest (isotropic, and rotated anisotropic with a
  non-symmetric tangent, 3D, plane strain and plane stress, non-linear
  problem): differences to point by point at most 5e-14, eight mutations
  detected. **The suite passes entirely with all the computations forced
  to the kernels (156/156).**
- Review T (`scratch/night-2026-10-10/review-T.md`): T1 (severity 4) the
  test behaviour HeatStorage had the wrong sign for the thermal convention
  of mfem-mgis (residual ∫ grad N . j with j = -k grad T): the test
  validated the backward heat equation; to be fixed, with a physical check
  of the lag of the transient.
- mfem-mgis 4c61469c33: fixes from review T. HeatStorage returns
  r = -C dT/dt (sign convention of the heat transfer of mfem-mgis,
  documented in TransientHeatTransfer with the initialisation of the
  gradients of s0); physical check independent of the path: the transient
  temperature lags behind the stationary one (16 to 37 K on average).
  Diagonal of the matrix free tangent checked in the thermal and
  micromorphic tests (before the update of the materials, otherwise the
  check was empty). NURBS meshes refused. Long debug variants only in an
  optimised build. Three mutations detected.
- Review M of 9b93bc8c36 and 4c61469c33 (`scratch/night-2026-10-10/review-M.md`):
  no wrong result, identical point by point results confirmed; fixes in
  a066c82baa (dead check removed, documentation, saturated damage in
  MicromorphicKernelsTest).
- mfem-mgis 3410a92b3b: **the whole suite on the GPU**. Variable
  `MFEM_MGIS_DEVICE`. The first pass revealed defects of the GPU path: the
  ParaView export overwrote the unknowns computed on the GPU with the
  stale values of the host (reference `MakeTRef`), L2 error and GSLIB
  interpolation read without synchronisation, `Elasticity` strategy of
  BoomerAMG refused by hypre on GPU (fallback to `System`), fill level of
  `HypreILU` lost at the second operator (also on CPU) and ILU(1) crashing
  in hypre on GPU (reproducer `reproducers/hypre-ilu-gpu`, ILU(0) by
  default on GPU). Tests made portable. Result: 145/147 on the GPU (the two
  tests that require the absence of a device fail), with and without the
  forced kernels; 157/157 on CPU. Procedure in `docs/validation.md`.
- Sweep of the BoomerAMG options for NewtonK0 + Eisenstat-Walker 1
  (`cases/measure-L5-amg.sh`, `results/L5-beam-2026-10-10-1034/README.md`):
  aggressive coarsening on one level is the best setting, 7.22 s on 8
  cores (-18 %, with a threshold of 0.5) and 1.71 s on the GPU (-5 to
  -8 %); the plain Jacobi of hypre makes the conjugate gradient fail.

- Review D of 3410a92b3b and a066c82baa
  (`scratch/night-2026-10-10/review-D.md`): no wrong result
  introduced. An earlier defect of severity 5 (D6): the L2 projection and
  the implicit gradient regularisation of a function with several
  components gave zeros on the GPU (copy of the components through the
  operator `[]`, which reads the host without making it valid), invisible
  because the tests exported the result without checking it.
- mfem-mgis 5246498277: fix of D6 (`HostRead`, `HostReadWrite`).
  L2ProjectionTest compares the projection of linear functions with one
  and three components, which is exact, to their interpolant: the test
  failed on the GPU before the fix (difference equal to the largest
  value) and now passes on CPU and GPU, in sequential and in parallel.
- mfem-mgis 080c20d1e3: other findings of review D.
  - Option `HypreILUType` of `HypreILU`: ILUT works on the GPU (7
    iterations against 10 for ILU(0) on one process), the refusal of the
    fill levels greater than 0 mentions it; a level is refused for ILUT.
    The version of hypre is 3.2.0 (the message of 3410a92b3b wrongly says
    3.0.2).
  - HypreILUTest: variant with one process (10, 7 and 5 iterations for
    levels 0, 1 and 2, against 18, 17 and 17 on two), default level
    checked on CPU, ILUT checked; two mutations detected (level 0 by
    default on CPU, type not passed to hypre).
  - FixedPointSolverTest: Crossed2Deltabis is no longer bounded (240 to
    562 iterations depending on the rounding, even on CPU); the margin of
    the other cases on GPU goes from 50 % to two iterations (they differ
    from the CPU by one iteration at most); three green GPU passes.
  - Check `ndof` of the cache of the shape functions removed (it could
    not fail). Documentation: `isHypreUsingGPU`, `MFEM_DEVICE` (takes
    precedence over `MFEM_MGIS_DEVICE`), mechanism of `MakeTRef` (alias
    with frozen validity flags, with or without prolongation: the message
    of 3410a92b3b wrongly restricts it to the case without prolongation),
    pairwise kernel also run on the debug device in Debug and Coverage.
- Reproducer `reproducers/hypre-ilu-gpu`: mode 2, level 0 given once to
  `mfem::HypreILU`: CPU 18 then 17, 17 iterations (level lost at the new
  operator), GPU 18 then segmentation fault.
- AmgX in the linear solver factory of mfem-mgis (preconditioner and
  solver `AmgX`, cycles of the AmgX multigrid, configuration by string or
  by file; solver `AmgXPCG`, conjugate gradient of AmgX with relative and
  absolute tolerances, convergence judged on the number of iterations as
  for hypre): local branch `amgx-wip` (2755f98996), not pushed. The AmgX
  of the environment (commit 91a8413e, CUDA 13.0.2, sm_120) crashes in the
  setup of the classical multigrid, even with its own example
  `amgx_capi`: "invalid device function" in `cub::DeviceFor::Bulk` called
  by `thrust::sequence` (Galerkin product, csr_multiply.cu). The
  aggregation multigrid works (Poisson 20³: 22 iterations, 18 ms). Cause
  located: two objects of libamgx.a define `DeviceFor::Bulk` for this
  kernel; the one of amg.cu.o (call from `Matrix::colorMatrix`, inline
  virtual function of a class template) has the host stub of the kernel,
  internal with the CUDA 13 default `-static-global-template-stub=true`,
  but its cubin does not contain the kernel. Trial: AmgX rebuilt with
  `-static-global-template-stub=false` (build/amgx-stub).
- Spack reinstallation of the patched amgx, of gslib +shared and of MFEM
  (rerun of the README): amgx and gslib installed, build of MFEM stopped
  around 14:30 (system memory too low), to be restarted later. The
  previous spack.lock is put back (MFEM zxpdstv), the new one is kept in
  `backups/mmm-gpu-spack.lock-amgx-patch-a-installer`.
- Review X of the AmgX work package (`scratch/night-2026-10-10/review-X.md`):
  X1 (severity 4), destroying the first AmgX solver makes the others
  fail; X2 to X5 (severity 3): verbose path that aborts, sequential then
  parallel solver, GPU 0 imposed, tolerances poorly tested. X1, X3, X4
  and the convergence status come from the `AmgXSolver` wrapper of MFEM.
- mfem-mgis de2e76a913 (branch `amgx-wip`): own wrapper over the C
  interface of AmgX, without `mfem::AmgXSolver`.
  - The AmgX resources carry global objects (memory pool, cuSPARSE and
    cuBLAS handles) released by the destruction of any resources, and
    keep the address of their communicator: a single resource per
    process, created for the first solver on `MPI_COMM_WORLD`, destroyed
    with the last one.
  - `exception_handling=1` makes AmgX abort (AmgX itself handles the
    error): without it, the errors come back as return codes; an invalid
    configuration is refused by the factory.
  - Convergence status of AmgX; single precision (`Precision`: `Single`,
    mode `dFFI`, conversions on the GPU); the mixed precision `dDFI` is no
    longer supported by AmgX 2.5.
  - `AMGX_vector_set_zero` resizes the vector without the halo of a
    distributed matrix: the zero initial solution is uploaded instead.
  - Option `enable-amgx` and `FindAMGX.cmake` (`AMGX_DIR`, searched before
    `CMAKE_PREFIX_PATH`), AmgX placed before MFEM at link time (the MFEM
    of the environment contains its own copy of AmgX, not patched),
    dependency passed on by `MFEMMGISConfig.cmake`.
  - AmgXSolverTest green in the 4 variants (sequential, 2 processes,
    without and with `--device cuda`), 10 mutations detected in total.
- Beam case: options `-amgxo` (AmgX parameters), `-amgxs` (single
  precision), tool `beam-k0` (operator K0 alone, same discretisation and
  same behaviour as the beam), 4 AmgX tests checked by the norm (15/15).
  Reproducer `amgx-classical-cuda13`: compares the given AmgX
  installations (without the patch, classical fails; with the patch
  installed by spack, 14 iterations).
- AmgX measurements (`cases/measure-L5-amgx.sh`, results/L5-beam-2026-10-10-2259):
  on K0 alone, tuned AmgX (threshold 0.8, MULTIPASS interpolation) in
  single precision solves in 53 to 61 ms (110 iterations) against 70 ms
  for BoomerAMG (70 iterations); default setting of MFEM 1078 iterations.
  Complete beam on GPU: tuned AmgX single 1.58 to 1.60 s, aggressive
  BoomerAMG 1.77 to 1.87 s, tuned AmgX double 1.95 to 2.04 s (four runs).
  On 8 CPU processes with AmgX on the shared GPU: 128 to 849 s, unusable.
- Upstream PRs opened: MFEM #5604 (shared resource keeping its own
  configuration and a duplicated communicator, `InitSerial` only once;
  checked before/after) and AmgX #386 (patch of `amg.cu`, closes #380).
  The MFEM rule on AI is respected by a line of the description (the label
  is reserved to the maintainers). The rebuild of MFEM and AmgX will wait
  for the upstream merge, then for the merge in spack-packages; until
  then, `scripts/env.sh` sets `AMGX_DIR` to the patched AmgX.
- Polycrystal with AmgX (mm-opera-hpc 838144c: options
  `--elastic-operator-linear-solver`, `--amgx-parameters`,
  `--amgx-single-precision`; `cases/measure-L5-polycrystal-amgx.sh`,
  results/L5-polycrystal-2026-10-10-2341): BoomerAMG remains better
  (9.04 s against 10.70 s for tuned AmgX in single precision). The AmgX
  multigrid adds Newton iterations (245 integrations against 225), and
  the integration of the behaviour dominates (75 %). The gain of AmgX
  depends on the case: it pays off when the linear solve dominates
  (beam).

### 2026-10-11 (night), review Y of the AmgX wrapper

- Review Y (independent) of the AmgX wrapper of mfem-mgis: Y1 and Y2 of
  severity 3, the rest of 1 or 2. Correction of a cause given above: the
  AmgX defect with CUDA 13 does not come from the internal host stub
  (`-static-global-template-stub=false` does not fix it); the weak copy of
  `cub::DeviceFor::Bulk` (thrust::sequence<int>, from the inlined
  `Matrix::colorMatrix`) in amg.cu.o launches a null kernel, fixed by the
  patch of src/amg.cu.
- Fixes in the wrapper:
  - Y1: the resources kept a copy of `MPI_COMM_WORLD` made at the first
    solver, which is a collective operation: a sequential solver created
    on a single process deadlocked. AmgX only keeps the address of the
    communicator and duplicates it itself for its exchanges: the handle
    of `MPI_COMM_WORLD` is kept as is, creation and destruction are no
    longer collective.
  - Y2: a hypre matrix on a communicator of a single process, in a
    computation with several processes, was uploaded as a distributed
    matrix on the resources of `MPI_COMM_WORLD` (read outside the array
    of offsets): it is uploaded as a local matrix.
  - Y11: collective checks in `SetOperator` (each process owns rows, at
    most 2^31 - 1 rows).
  - Y7: after a failure of `SetOperator`, the solver refuses to solve (it
    kept the old matrix); sizes of the vectors checked.
  - Y5: convergence floors of AmgX (1e-12 of the initial residual in
    double, 1e-6 in single): `AmgXPCG` refuses a smaller relative
    tolerance and requires a relative tolerance in single precision.
  - Y6: in verbose mode, the residual computed for the statistics could
    stop the cycle before it starts (small right-hand side, zero
    solution): negative tolerance for the multigrid.
  - Y13: copy of the solver forbidden.
- Tests (Y8): verbose AmgXPCG and default options, right-hand side
  multiplied by 1e-20 (norm 5e-23) in verbose mode, failure of `SetOperator` then solve,
  vector of wrong size, floors, sequential solver on the first process
  only, communicator of one process, processes in reverse order, process
  without rows, complete reinitialisation of AmgX in parallel, test
  without a visible GPU (`AmgXSolverTest-no-gpu`, `CUDA_VISIBLE_DEVICES=`
  empty). The failures of the factory raise an exception instead of
  calling `std::exit`, to destroy the solvers before exiting (double copy
  of AmgX, Y4).
- Documentation (Y3, Y4, Y9, Y10, Y12): header of the wrapper and
  docs/gpu-solvers.md (link order, destruction before `MPI_Finalize`,
  exit with a live solver, floors), figures of the README of the beam
  measurements fixed, warning of build-mfem-mgis.sh when the AmgX used
  does not have the patch for CUDA 13 (fingerprint of the patch in the
  spack spec of the installation).
- Review Z (independent) of these fixes: Y1, Y2, Y6, Y7, Y9 and Y11
  correct. Z1 (severity 2): the relative criteria of AmgX also have an
  absolute floor, convergence is declared as soon as the norm of the
  residual goes below 1e-20 (and immediately if the initial norm is
  already below); AmgXPCG stopped too early with a success for a
  right-hand side of norm lower than 1e-20 divided by the tolerance
  (probe: actual relative residual 1.5e-8 for a tolerance of 1e-10 and a
  norm of 5e-13, zero solution for 5e-22). Fixed: without an absolute
  tolerance, the right-hand side is multiplied by a power of 2 that brings
  its norm between 1 and 2, and the solution by the (exact) inverse; same
  for the preconditioner in single precision (range of single precision);
  absolute tolerance lower than 1e-20 refused. Z2: the collective checks
  raise a `CollectiveFailure`. Z5: failure status as long as AmgX has not
  given its own (an exception left the status of the previous solve).
  Tests: right-hand sides multiplied by 1e-15, 1e-25, 1e-40 and 1e45
  (double and single precision), zero right-hand side, right-hand side of
  wrong size, status after a failure.
- Short review of dbc963185e (three findings of severity 1): the factor
  taken from the norm of the right-hand side overflows at the extremes
  (zero norm below about 1e-162, infinite beyond 1e154: zero solution
  reported as converged); it is now taken from the maximum of the
  absolute values of the components (reduction on the GPU, maximum over
  the processes). In single precision, an absolute tolerance prevented the
  scaling and the scalar products of AmgX overflowed beyond 1e19: it is
  refused. Documentation of the statistics of the scaled system and of
  the tolerances passed through `AmgXParameters`. Tests: factors 1e-158
  and 1e158, right-hand side with different scales depending on the
  process (common factor), residuals of the tests computed on scaled
  vectors; three mutations detected.

### 2026-10-11 (night), transfers between the host and the GPU (L7.3)

- Profile of the copies (`cases/transfers-profile.sh`,
  results/L7-transfers-2026-10-11): 8 % of the time of the heat transfer
  r5, 3 % of the beam r4. Origin given by gdb: mostly the end of step,
  which brings the state arrays back to the host and invalidates their
  copy on the GPU (sent again at the next integration), then the
  temperature at the points (values of the unknowns copied back to the
  host, external variables copied to the GPU at each integration).
- Trial of a read-only synchronisation after the update (GPU copy kept
  valid): it breaks the contract "direct modification of the state on the
  host after an update" of `MaterialArrays`, checked by
  `BatchIntegrationTest` on `debug` and `cuda`. Cancelled.
- Temperature kept on the GPU: gain of at most 2.3 % on the heat
  transfer, not done. Work package L7.3 closed on these measurements.

### 2026-10-11 (night), scaling of the beam

- `cases/measure-beam-scaling.sh`, results/scaling-beam-2026-10-11-0247:
  beam refined 3, 4 and 5 times (16,000 to 840,000 unknowns), NewtonK0 and
  Eisenstat-Walker 1. Gain of the GPU over 8 cores: ×1.8, ×4.5 then ×10.9
  (9.86 s against 107.9 s at 840,000 unknowns). Tuned AmgX in single
  precision is the best up to 110,000 unknowns, aggressive BoomerAMG
  beyond. The tangent without assembly does not pay off with NewtonK0.
- Six refinements (6.5 million unknowns): GPU memory exhausted (Umpire
  allocations refused around 13.5 GiB, then around 15.3 GiB with a
  reduced hypre pool). The driver of the beam now prints the error
  message of a non-converged step (before: "not converged" without a
  cause) and has a verbosity option (`-vl`).
- hypre 3.2 allocates a pool of 4 GiB on the GPU at once: mfem-mgis
  fd93974cdc, option `--hypre-device-pool-size` (or `-hdps`, or variable
  `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`) applied by `initialize`; at
  0.5 GiB, the beam with four refinements uses 2.29 GiB instead of
  5.79 GiB, for a time 2 to 4 % longer with BoomerAMG (and not "at the
  same time", in translation, as the commit message says).
  HypreDevicePoolSizeTest (7 cases, 5 mutations detected); suites 165/165
  CPU, 15/15 GPU, cases 33/33.
- The measurements on 8 cores of the first pass
  (results/scaling-beam-2026-10-11-0233) are 13 to 43 % slower than
  those of the second, for a cause not established (compilations of
  another session were running, but the rerun at load 1.1 gave the same
  times).
- Review E (scratch, review-E.md): option of mfem-mgis checked in the
  sources of hypre and by a probe on the GPU (exact size of the pool,
  growth, one pool per process). Fixed: times and causes in the README,
  units (GiB), number of unknowns of six refinements (6,502,275); in
  mfem-mgis, option declared as a string (like `--device`), bounded size
  (at least one byte, in the range of `std::size_t`), global error flag of
  hypre ignored, environment variable ignored with a warning if hypre does
  not use Umpire, documentation (one pool per process), test on the GPU
  that allocates and checks the actual size of the pool (the only one to
  detect a setting applied too late).

### 2026-10-11 (night), FBar by kernels (L7.9)

- FBar formulation computed by kernels (`batch::FBarData`,
  `src/Batch/FBar*.cxx`): the gradient at the centre of each element comes
  from an integration rule of order 0, computed with the other rules (new
  method `getAdditionalIntegrationRulesForKernels` of the integrators);
  the standard gradients kernel gives F, then a kernel replaces F by
  `F̄ = α F`, `α = (J0/J)^(1/3)` (1/2 in plane strain, `F̄33 = 1`), keeping
  F and α. Internal forces: those of the standard finite strain
  formulation. Tangent: the standard tangent with
  `K' = α K + (K F) ⊗ dα/dF`, plus a rank one matrix per element `a ⊗ b`
  (`a`: internal forces of a pseudo-stress `α/(ie J0) K F`,
  `b = dJ0/du`), for the assembled matrix, the diagonal and the matrix
  free action. The common helpers of the kernels (data, blocks of the
  tangent, internal forces) move from the standard base to
  `BehaviourIntegratorBase`, shared with the FBar base.
- Tests: `BatchIntegrationTest` compares with the FBar point by point
  computation, isotropic and orthotropic, in 3D and in plane strain, with
  and without prediction, assembled and matrix free (diagonal included),
  generic, openmp, stlpar and cuda behaviours, CPU, `debug`, GPU, 2
  processes: differences from 1e-15 to 1e-12. Seven mutations detected
  (rank one term of the matrix, of the action and of the diagonal, term
  `v ⊗ c`, forgotten α, pseudo-stress, F̄33 in plane strain). Faltus is
  still refused (GradientsKernelsTest).
- Review F (scratch, review-F.md), fixed in mfem-mgis (639db23788 to
  ad1d81db9e): FBar data marked invalid when the gradients are computed
  element by element, the tangent by kernels being then refused with a
  message; data of the beginning of the step kept by `update` and
  restored by `revert` (a prediction after a revert read data of another
  state); exact number of threads per element in the overflow check;
  vectors `a` and `b` of the rank one matrices stored per element of the
  material (new function `addMaterialInnerForces`), `centers` reduced to
  J0. Tests: face x = 0 clamped (non-homogeneous solution), FBar case with
  pressure (elementary matrices), two materials where the studied
  material is the second one of the mesh (different indices in the mesh
  and in the material), quadratic elements, refusal of the tangent by
  kernels without gradients by kernels, two predictions separated by a
  revert. In parallel, the residual of the FBar cases with pressure
  stagnates around 1e-12 times its reference (MUMPS): relative tolerance
  raised to 1e-11 for these cases (option `relative_tolerance`). CPU,
  `debug`, GPU, 2 processes (CPU and GPU), GradientsKernelsTest: 8/8.
- Review G (scratch, review-G.md) of the fixes: nothing of severity 3
  or more. Handled: the FBar tangent by kernels without gradients by
  kernels is refused at the setup of the solve (before, it was refused at
  the first assembly, seen as a non-convergence and thus cut again in a
  `Simulation`), tangent already computed forgotten together with the
  data, documentation of `TangentOperatorComputation`, refusal also
  tested without assembly, prediction in the order 2 case. The finding on
  the prediction without revert (iterative coupling) is not specific to
  FBar: the prediction integrates the gradients of the last iterate, and
  the FBar data are those of these gradients.

### 2026-10-11 (night), boundary conditions independent of the unknowns (L7.6)

- Imposed pressure and heat source: their residual is `s(t + dt) v`, `v`
  independent of the unknowns, zero tangent. New methods of
  `AbstractBoundaryCondition`: `isIndependentOfTheUnknowns` (false by
  default) and `addContributionToTheResidual`. `v` is assembled only
  once, for a unit value, by the linear form integrator already used by
  the prediction, then kept on the GPU. With the internal forces computed
  by kernels, these conditions are no longer integrators of the
  non-linear form: their contribution is added to the local vector of the
  internal forces before the passage to the true degrees of freedom. The
  form then only contains the integrator of mfem-mgis: no more loop over
  the elements on the host for the residual (copies of `u` and of the
  residual), tangent assembled by kernels without the assembly of MFEM,
  tangent without assembly possible with a pressure (refused before).
- Tests: BatchIntegrationTest (pressure with the tangent assembled by
  kernels or without assembly, and as an integrator of the form when the
  internal forces are computed element by element; refusal of the
  no-assembly mode with a test condition that depends on the unknowns),
  HeatTransferKernelsTest (heat source, stationary and transient,
  differences of 1e-16). With orthotropic FBar, a pressure and the tangent
  without assembly, the conjugate gradient of the test does not converge
  in 10,000 iterations (non-symmetric tangent): this case stays isotropic.
  Four mutations detected (doubled pressure, forgotten source, pressure
  declared dependent, boundary integrators ignored).

### 2026-10-11 (morning), ssna303 and stop for lack of memory

- Measurements of ssna303 (case 2 of the campaign,
  results/ssna303-ssna303-2026-10-11-0625): GPU ×2.4 at 200,000 unknowns
  (×2.7 with FBar), the hypre solve dominating on the GPU. Series stopped
  at 7:18 during `ssna303-r2-fbar-np8-kernels` (1.46 million unknowns,
  FBar, 8 processes): the memory of the machine (30 GB) was almost
  exhausted, and the background tasks of the session (measurements, and
  the chain leading to the installation and the measurements of the
  pressure and of the post-processings) were stopped. Nothing was
  restarted. The computation with two refinements on the GPU fails from
  the first step; ex4 now prints the cause of the solve failures
  (mfem-mgis-examples db0ce2a).

### 2026-10-11 (morning), Faltus regularisation by kernels (end of L7.9)

- Module `batch/Faltus` of mfem-mgis (ed19e652f2, tests 8e4a341cfd):
  gradient at the centre F0 computed with the gradients by kernels (rule
  of order 0 shared with FBar, computation of F0 factored out), kept by
  `update`, restored by `revert`; residual `Σ w Dᵀ α (F - F0)`,
  `D = B - B0`; tangent `α Σ w Dᵀ D`, independent of the unknowns,
  computed on the fly by pairs of nodes: for a deformation gradient, the
  block of nodes a and b is `d_a · d_b` times the identity (`d`:
  difference of the derivatives of the shape function at the point and at
  the centre), that is about 20 times less computation than term by term.
  Internal forces by kernels without gradients by kernels refused at
  setup.
- Tests: distorted meshes (non-affine elements, all different, option
  `distorted_mesh`), also for FBar; action of the no-assembly mode
  compared to the assembled operator; GradientsKernelsTest compares FBar
  and Faltus term by term (distorted hexahedra of orders 1 and 2,
  tetrahedra, distorted quadrilaterals), with a penalty coefficient of
  10: the differences between the element by element computation and the
  kernels are proportional to the coefficient (cancellation of the terms
  in B and B0 of the reference, rounding errors of F multiplied by α).
  Suite 183/183, BatchIntegrationTest CPU, debug, GPU, 2 processes CPU
  and GPU green, 15 mutations detected (9 then 6 on the pairwise kernels,
  including a confusion of the Jacobians of the point and of the centre
  that only the distorted meshes see). Review I handled (no wrong
  result), review J in progress.

### 2026-10-11 (morning), order of the integration rules (L4)

- Parameter `IntegrationRuleOrder` of the discretisation (mfem-mgis
  c941836702): order of the rules of all the behaviour integrators; by
  default, 27 points on the order 1 hexahedron, 8 with order 3 (exact
  integration of the elastic stiffness of an affine hexahedron), that is
  3.4 times less memory and work at the points. The 26 integrators build
  their quadrature space with this rule and their `getIntegrationRule`
  returns the one of the space. Negative order and rules with
  non-positive weights refused. IntegrationRuleOrderTest (host and
  `debug` device), suite 184/184, 4 mutations detected, review K handled
  (documentation of the risky cases: under-integration of order 2
  hexahedra, one point at orders 0 and 1).
- Motivation: ssna303 with two refinements does not fit on the GPU with
  27 points per element; the measurement with 8 points remains to be done
  (series of measurements stopped for lack of memory, to be restarted
  after approval).

### 2026-10-11 (morning), FBar by kernels without storing K'

- mfem-mgis b1a8c682ae: the modified tangent of FBar,
  `K' = α K + v ⊗ c`, is no longer stored at the points (81 reals per
  point in 3D); the kernels of the stiffness matrices apply it on the
  fly, `K' x = α K x + v (c · x)`, from `v` and `c` (18 reals per point):
  about 6 GB less for ssna303 with two refinements (12 million points),
  probable cause of the memory stop. New optional argument of the
  stiffness matrix functions (`TangentOperatorModification`, compiled
  for the deformation gradients only); the FBar functions check that they
  receive the blocks of the tangent that were used to compute `v` and
  `c`. Suite 185/185, 3 mutations detected, review L handled (shadowing
  of the member `K`, documentation).

### 2026-10-11 (morning), positions of the entries and matrices without forming B

- mfem-mgis d04acab66c: `batch::computeStiffnessMatricesPositions`
  computes once the position of each entry of the element stiffness
  matrices in its row of the sparse matrix (`std::uint16_t`, relative to
  the beginning of the row). The three assembly kernels (standard, rank
  one of FBar, Faltus) read the position and only use it if it designates
  the column of the entry; otherwise the entry is searched and counted
  (`AssembledStiffnessMatrix::searched_entries`). The multi-material
  integrator keeps positions for each layout of the rows (columns sorted
  by a direct solver, or diagonal first), without GPU, recomputes them
  after an assembly where entries were searched, at most three times, and
  switches to the search if their computation fails. On one core:
  assembly ×1.1 to ×2.9; on the GPU, no gain (positions not used).
  Reviews M2 and M3 handled (nothing above 3/5; the detection of the
  benchmark by CMake always failed), 14 mutations detected out of 15 (the
  15th only concerns rows of more than 65,535 entries).
- mfem-mgis 4a007e85f2: products `K B_b` and `B_aᵀ (K B_b)` of the
  stiffness matrices, and action of the tangent without assembly, written
  with a constexpr table of the non-zero entries of `B`
  (`getBOperatorEntries`, compared to `computeB` at compile time), in the
  same order of operations: values identical to the bit on the host.
  Kernel of the matrices in 3D finite strain: 254 registers and 576 bytes
  of stack before, 64 and 72 after. Assembly ×1.1 to ×3.1 on one core,
  ×1.4 to ×1.9 on the GPU (heat transfer ×0.94). Review N handled
  (nothing above 2/5; exact compile-time check, diagonals only), 12
  mutations detected.
- Benchmark `cases/assembly` (assembly of all the matrices, with and
  without positions, mechanics in small and finite strain, heat transfer,
  host or GPU, writing of the values to compare two versions). Results:
  `results/assembly-positions-2026-10-11-0901`,
  `results/registers-kernels-2026-10-11`. Suites 185/185.

### 2026-10-11 (morning), global matrix and inverses of the Jacobians

- Benchmark `cases/assembly` in parallel (`--parallel`): local assembly by
  kernels, then global matrix as in mfem-mgis (diagonal block, `Pᵀ A P`
  of hypre, elimination). On 8 processes, the product of hypre costs more
  than the local assembly (heat transfer r5: 36.5 against 14.8 ms; beam
  r4: 15.4 against 9.5 ms); on the GPU with one process, it is small
  (`results/global-assembly-2026-10-11`). Lead: pattern and exchanges of
  the global matrix kept.
- Nsight Compute: the kernels of the matrices on the GPU are compute
  bound (85 to 89 % of the throughput of the multiprocessors). Each thread
  inverted the Jacobian at each point. mfem-mgis 7af829b3a5: the geometry
  kernel stores the inverses of the Jacobians in place of the Jacobians
  (same memory), all the kernels read them. Assembly on the GPU ×1.7 to
  ×2.2 in mechanics, ×5.9 in heat transfer
  (`results/jacobians-inverses-2026-10-11`); on the host, same values to
  the bit. End to end on the GPU: heat transfer ×1.22, beam ×1.05, the
  solve dominating (`results/jacobians-gpu-2026-10-11`). Suite 185/185, 9
  mutations detected, review O handled (nothing above 1/5).

### 2026-10-11 (late morning), global matrix with kept pattern and exchanges

- mfem-mgis e111df8454: `batch::GlobalStiffnessMatrix`. For a conforming
  space, each entry of the local matrix goes to a fixed entry of the
  global matrix: the first assembly computes the pattern (product of
  hypre on a matrix of ones), the position of each local entry or its
  place in the message to the process that owns its row; the following
  ones copy, exchange (`MPI_Alltoallv`) and add the values. If the local
  pattern (arrays) or the space changes on one process, all of them redo
  the setup. Used by the tangent assembled by kernels when hypre does not
  use the GPU.
- Benchmark on 8 processes (`results/global-matrix-2026-10-11`): beam
  15.4 -> 2.5 ms, heat transfer 36.8 -> 1.7 ms. End to end
  (`results/global-matrix-end-to-end-2026-10-11`): tangent ×2.4 (heat
  transfer) and ×1.8 (beam), same solutions. Suite 188/188,
  GlobalStiffnessMatrixTest on 1, 2 and 4 processes, 9 mutations detected
  (the mutation that does not use the class only changes the time),
  review P handled (nothing above 2/5).
- Pitfall: the behaviour library of the tests (`BehaviourTest`) is
  neither in `all` nor in `build-tests`; a separate build must build it
  explicitly, otherwise a dozen parallel tests fail.

### 2026-10-11 (noon), documentation of the batch path (L10)

- mfem-mgis: guide page `batch_computations.rst` (selection of the batch
  and kernel computations, parameters and environment variables, checks,
  integration rules, what the kernels cover, assembly and application of
  the tangent, prediction, GPU, AmgX, example, limits), section on the
  FixedPoint solver in `nonlinear_algorithms.rst`, preconditioners
  `AssembledOperator` and `AmgX` in the glossary. Each statement checked
  in the code; review Q (independent fact checking) handled, including the
  name of the behaviour library produced by the generic-parallel
  interface in the example. Reduced Sphinx build without warnings on these
  pages.

### 2026-10-11 (early afternoon), other integrators of the form by kernels

- mfem-mgis 517e33cde4: when the non-linear form contains other
  integrators than the one of mfem-mgis (Robin condition on boundary
  faces, springs depending on the unknowns), the residual and the tangent
  by kernels no longer go through the loops of MFEM over all the
  elements. The domain, boundary and boundary face integrators are
  evaluated by loops copied from MFEM and restricted to them (markers
  included), on the host, then added to the internal forces and to the
  local matrix assembled by kernels; the positions of the entries and the
  global matrix with kept pattern are also used in this case. Fallback to
  the old path with interior face integrators or an extension of MFEM.
  The prolongation of the unknowns is done on all processes before the
  computations that can fail on a single one. Guide of mfem-mgis up to
  date.
- Tests: boundary, boundary face and domain springs with markers,
  interior boundary (faces skipped as by MFEM), fallback through a null
  interior face integrator, counters of the evaluations of the residual
  and of the assemblies of the tangent. Suite 188/188, 10 mutations
  detected, reviews R and S handled (one 3/5 defect: possible deadlock of
  the residual in parallel if a single process failed before the
  prolongation).
- Thermal case: option `-hc` (convective exchange at the second end,
  boundary face integrator). On 8 processes (heat transfer r5,
  `results/other-integrators-2026-10-11`): loop over the steps 4.17 ->
  2.70 s, tangent ×3.3, residual ×2, same solutions and iterations.

### 2026-10-11 (afternoon), ex8 entirely by kernels, on CPU and on GPU

- ex8 of mfem-mgis-examples (heat transfer with Robin condition,
  mechanics with swelling, iterative coupling) runs entirely by kernels
  with mfem-mgis 517e33cde4: same statistics as point by point, on CPU and
  on GPU. On the GPU, the example read its results in stale host copies
  (comparison to the reference failing, solves correct): fixed in the
  example (mfem-mgis-examples 3a907f1), with the tests rjh_plate-kernels
  and, under the option `MFEM_MGIS_EXAMPLES_GPU_TESTS`, rjh_plate-cuda and
  rjh_plate-kernels-cuda.
- Measurement (`results/ex8-2026-10-11`, original mesh): CPU 2 processes
  31.0 -> 22.3 s (×1.39), GPU 1 process 167.7 -> 38.3 s (×4.4). What
  remains is the integration of the behaviours on the host and MUMPS; the
  mesh refined once is out of reach of the solvers tried (MUMPS, PCG and
  BoomerAMG).

### 2026-10-11 (mid-afternoon), ex8 entirely on the GPU

- The behaviours of ex8 compile with the generic-parallel interface
  (openmp, stlpar, cuda), once `computeHardeningModulus` is made
  `constexpr` (mfem-mgis-examples 89718dd). New case `cases/ex8`
  (libraries GenericParallelU3SI2 and GenericParallelALFENI of each
  backend, tests ex8-batch-*), `add_generic_parallel_mfront_libraries`
  accepts several libraries and an include directory.
- mfem-mgis 8011fe9eaa: a model at the points (`PointWiseModel`, swelling
  of ex8) generated for the GPU is integrated by batches on the GPU, then
  its host copies are made readable (other models read them through
  external storage). The MGIS view on the device copies of the arrays
  (`Material::makeDeviceView`) and the criterion
  (`isIntegratedByBatchesOnDevice`) are shared with the integrators.
  MaterialArraysTest: two steps with update, on the host, on `debug`
  (device path exercised, host copies protected) and on the GPU. Suite
  188/188, mutations 3/4 (the copy of the inputs to the GPU cannot be
  detected on this machine, whose GPU reads the host memory), review U
  handled.
- ex8: the heat transfer reads the mechanical gradients on the host; when
  the mechanics is integrated on the GPU, they must be synchronised
  (FieldUpdaterModel, mfem-mgis-examples 6505b78; without this,
  difference of 3e-5 to the reference).
- Measurement (`results/ex8-behaviours-cuda-2026-10-11`): ex8 on the GPU
  in 18.5 s with the cuda behaviours, against 37.8 s with the behaviours
  on the host, 166 s point by point on the GPU and 22.1 s entirely by
  kernels on 2 CPU processes. The integration goes from 19.6 to 0.66 s;
  MUMPS now takes 15 s out of 18.5.

### 2026-10-11 (late afternoon), validations and solvers of the mechanics of ex8

- Validations of `docs/validation.md` on gpu-proto after the work
  packages of the day: forced kernels 163/163; the whole suite on the GPU 159/161 (the
  two tests that require the absence of a device fail).
  GlobalStiffnessMatrixTest assumed that the global matrix is accepted
  when hypre computes on the GPU, whereas the class refuses it in that
  case (mfem-mgis 76c2b50b65); the tests with 4 processes on the GPU
  require `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1` (hypre reserves 4 GiB per
  process). FixedPointSolverTest failed once in 26 runs, not reproduced.
- Mechanics of ex8, to reduce MUMPS (15 s out of 18.5 s on the GPU):
  PCG + BoomerAMG requires about 240 iterations per solve (93 s against
  22 s with MUMPS on 2 processes); NewtonK0 (option `-nlsMc NewtonK0`,
  mfem-mgis-examples d8287e4) converges like Newton but GMRES
  preconditioned by the elastic operator costs too much (28.8 s against
  13.1 s); ModifiedNewton diverges from the first step. The lever remains
  a preconditioner suited to this thin plate with contrasted materials.

### 2026-10-11 (late afternoon), GPU backends of the generic-parallel interface

- In mfem-mgis, these work: cuda (nvcc, clang), hip (clang, hipcc on the
  NVIDIA platform) and kokkos if the program initialises Kokkos
  (`cases/common/KokkosRuntime`, option `--kokkos` of the beam, Kokkos
  optional). These fail: sycl (the program must provide a SYCL queue),
  nvhpc stdpar and openmp offload (OpenMP runtime of nvhpc in conflict
  with the libgomp of MFEM and mfem-mgis compiled by GCC).
- `cases/CMakeLists.txt`: `add_generic_parallel_mfront_library` (one
  configuration, explicit name) under
  `add_generic_parallel_mfront_libraries`; the beam and the polycrystal
  have the libraries and the short tests of these backends (beam-gpu-*,
  polycrystal-batch-* and comparisons).
- Beam r4 (`results/backends-gpu-2026-10-11`): behaviour at 1 ms per
  integration, no difference. Polycrystal
  (`results/backends-polycrystal-2026-10-11`): integration 5.6 s with
  nvcc or hipcc, 7.7 s with clang (cuda or hip), 102 s on the host;
  stresses identical to within 2.4e-13.

### 2026-10-11 (evening), Kokkos initialised by mfem-mgis (option enable-kokkos)

- mfem-mgis (branch kokkos): CMake option `enable-kokkos` (disabled by
  default). `mfem_mgis::initialize` initialises Kokkos on request (option
  `--kokkos` of `declareDefaultOptions`, or variable `MFEM_MGIS_KOKKOS`),
  after the device of MFEM and on its GPU; `finalize` finalises it before
  MPI. Without the CMake option, `--kokkos` is an error and the variable
  is ignored with a warning. Tests: `KokkosOptionTest` (6 cases),
  `MaterialArraysTest-kokkos` on the GPU. Kokkos takes the GPU of the
  MFEM device, including when it is configured by the program before
  `initialize`.
- Design: the only source that includes Kokkos is compiled by the Kokkos
  compiler in a small shared library `MFEMMGISKokkos`, privately linked to
  Kokkos. A first version linked `MFEMMGIS`, static, to Kokkos with
  `$<LINK_ONLY:...>`: the link options of Kokkos
  (`-DKOKKOS_DEPENDENCE -arch=sm_120`) were then passed to all the client
  programs, and the `check_cxx_source_compiles` of `cases/` failed
  outside the Kokkos launcher. The clients (mm-opera-hpc, `cases/`) no
  longer look for Kokkos. The Kokkos launcher only applies to this target
  (component `separable_compilation`), compiled without
  `CMAKE_CXX_COMPILER_LAUNCHER`: ccache remains usable for the rest.
  Kokkos must be a shared library (checked at configuration).
- `cases/` detects the option of mfem-mgis
  (`MFEM_MGIS_HAS_KOKKOS_OPTION`): the beam then uses the one of
  mfem-mgis instead of its own. The polycrystal has a kokkos test
  (`MM_OPERA_HPC_POLYCRYSTAL_WITH_KOKKOS`).
- Polycrystal (`results/kokkos-polycrystal-2026-10-11`): the kokkos
  behaviour integrates in 5.62 and 5.69 s, the cuda behaviour in 5.61 s
  (both go through nvcc); stresses identical to within 2.8e-13.
  mm-opera-hpc was not modified.

### 2026-10-11 (evening), cermet and bubble of mm-opera-hpc, host memory on the GPU

- Cermet (UO2 polycrystal and metallic interface) and bubble (cavity
  under pressure, quadratic elements, post-processings) of mm-opera-hpc
  on the batch path: point by point, by batches, entirely by kernels, on
  the CPU (1 and 2 processes) and on the GPU (`generic` behaviours on the
  host and cuda behaviours): same results (cermet within 8.6e-15 at
  most). Tests in `cases/bubble` and `cases/cermet`, like the
  polycrystal; `compare-outputs.py --csv` for the outputs of the bubble.
  mm-opera-hpc: options of mfem-mgis (`--device`...) in both programs,
  option `-l` of the cermet, option `--linear-solver-preconditioner` of
  the bubble: its default case (467,109 unknowns) failed even point by
  point (PCG stopped at 1000 iterations with the diagonal, Newton limited
  to one iteration) and converges with BoomerAMG.
- On the GPU, the bubble took 68 s against 31 s on one CPU core, of which
  60 s of post-processings. Cause: MFEM compiled with Umpire also
  allocates the host memory with Umpire as soon as a device is
  configured, and each temporary of the element loops done on the host
  goes through the allocation table of Umpire (stacks taken with gdb; the
  CUDA calls only accounted for 1 s). Fix in mfem-mgis (branch
  host-memory): standard memory for the host, Umpire kept for the GPU,
  except with `MFEM_MEMORY`.
- Measurements (`results/host-memory-2026-10-11`): bubble 69.2 -> 14.0 s,
  ex8 point by point 164.6 -> 57.3 s, cermet point by point 28.1 ->
  14.8 s; little effect when everything is computed by kernels (beam,
  polycrystal, ex8 with the cuda behaviours: 18.6 -> 16.9 s, the Robin
  condition being evaluated on the host).

### 2026-10-11 (evening), post-processings by kernels (work package L8)

- Starting measurement (bubble of mm-opera-hpc, 467,109 unknowns, after
  the fix of the host memory): on the GPU, 10.5 s of post-processings out
  of 14 s; export of the stresses at the nodes 4.5 s, L2 projection 4.4 s,
  hydrostatic average 1.0 s.
- mfem-mgis (branch post-processings): with
  `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels`, integrals and measures
  by reduction, right-hand side of the L2 projection by kernel and
  operator assembled by MFEM on the device (mass of the simplices computed
  by kernel, since MFEM 4.10 does not do it), values at the nodes by local
  projections and averaging by the restriction to the elements and the
  prolongation (decision taken by all the processes before any exchange).
  The other cases remain element by element; a counter makes the choice
  observable. Defect fixed along the way: the L2 projection kept a
  reference to a destroyed matrix in full assembly.
- Comparison test of the two paths (hexahedra with two materials,
  tetrahedra, quadrilaterals; host, debug, two processes, GPU),
  differences of the order of 1e-14; complete suites 181/181 (forced
  kernels) and 174/174; six mutations detected.
- Bubble (`results/post-processings-kernels-2026-10-11`): GPU 14.8 ->
  5.3 s (post-processings 11.3 -> 1.9 s), one CPU core 33.2 -> 25.2 s.
  With the fix of the host memory, the bubble on the GPU went from 69 s to
  5.3 s. mm-opera-hpc: profiling sections of the post-processings of the
  bubble, average by `computeMeasure` and `computeIntegral`.
- Independent review of the work package (one 4/5 defect, five 3/5) handled:
  fallback to element by element beyond 32 degrees of freedom per element
  (the export stopped), mass computed by the kernel for all the elements
  (the assembly of MFEM gave nodes to the meshes that had none),
  collective decision of the L2 projections (possible deadlock on a mixed
  mesh), geometry of the elements of the material only for the integrals
  (1000 materials: 49.8 s against 0.04 s). Suites 180/180 (forced kernels
  and default), mutations 6/6, GPU variants 15/15; bubble unchanged
  (5.32 s on the GPU).

### 2026-10-11 (evening), ex8: BoomerAMG with the elasticity options

- Negative result (`results/ex8-amg-elasticity-2026-10-11`): for the
  mechanics of ex8, PCG and BoomerAMG with the elasticity options (byVDIM
  ordering) require 441 iterations on average, 201 s against 23.5 s with
  MUMPS (2 processes). Options `-oMc` and `-amgMc` added to ex8.

### 2026-10-11 (evening), ex8: symmetric matrices for MUMPS

- The tangent of the mechanics of ex8 is symmetric: given to MUMPS as
  symmetric positive definite (option `-msMc` of ex8, parameters
  `Symmetric` and `PositiveDefinite` of the MUMPS solver of mfem-mgis),
  same iterations and results, factorisations 29 % faster on the CPU (2
  processes) and 36 % on the GPU (on the host). ex8 entirely on the GPU:
  16.5 -> 12.4 s (`results/ex8-mumps-symmetric-2026-10-11`).
- Short review handled: option `-msMc` also for the elastic operator of
  NewtonK0, warnings for the options without effect, README of ex8, test
  `rjh_plate-mumps-spd` (mfem-mgis-examples 03c347e). On the GPU, the
  second cost item of ex8 is now the thermal solve (GMRES and BoomerAMG,
  3.15 s for 35 solves).

### 2026-10-11 (evening), ex8: solver of the heat transfer on the GPU, ParaSails

- Heat transfer of ex8 on the GPU (`results/ex8-thermal-gpu-2026-10-11`):
  MUMPS on the host (`-lsTh MUMPSSolver`) is the fastest, ex8 12.6 ->
  11.5 s (heat transfer 4.26 -> 3.03 s). GMRES of hypre on the GPU takes
  100 ms per solve for 9,516 unknowns, with BoomerAMG as with DiagScale.
  ILU(0), imposed on the GPU, does not converge (stopped at 600 s); block
  ILU(1) converges on the CPU.
- ParaSails crashed everywhere. Two defects:
  `mfem::HypreParaSails::SetOperator` (MFEM 4.10 and main branch) reads
  the communicator of the previous operator, which the non-linear solvers
  have already destroyed (`reproducers/hypre-parasails-operator`);
  ParaSails only computes on the host and crashes when hypre computes on
  the GPU. mfem-mgis 607f895a0d (branch parasails): preconditioner of the
  factory that points to the new operator before the call, refusal on the
  GPU, test `HypreParaSailsTest` (mutation detected). After the fix,
  ParaSails remains too weak for the heat transfer of ex8 (2,000 to 5,600
  GMRES iterations, then stagnation).

### 2026-10-11 (evening), upstream integration of mfem-mgis

- `docs/upstream-integration.md`: mfem-mgis section, split of the 115
  commits of `gpu-proto` into eleven contributions ordered by
  dependencies (fixes, features without GPU, GPU configuration, GPU path
  of the existing code, material data, batches, kernels, fixed point,
  post-processings, optional Kokkos and AmgX), dependencies on TFEL and
  MGIS, checks. MGIS section completed (layouts, copy operations,
  profiling). Independent review handled (classification of the boundary
  conditions, dependencies on MGIS, forgotten fixes).

### 2026-10-12 (night), examples ex1, ex2, ex3, ex5, ex6 and ex7 by batches and on the GPU

- Six examples of mfem-mgis-examples followed by the prototype
  (`cases/ex1`, `cases/ex2`, `cases/ex3`, `cases/ex5`, `cases/ex6`,
  `cases/ex7`): point by point, by kernels, by openmp and stlpar batches,
  on the GPU; CPU 48/48, GPU 12/12 (`results/examples-gpu-2026-10-12`).
  Options added to the examples: `-l` (ex2, ex5), `--device` (ex5, ex6).
- Complete ex7: CPU 2 processes 10.4 s point by point, 6.2 s by kernels;
  GPU 6.5 to 6.8 s, limited by GMRES and BoomerAMG on a system of 2,634
  unknowns.
- Behaviour of ex5 refused by nvcc: static variable of a quantity
  referenced in the GPU code. TFEL 837338523, then f1ff18259, d22416824
  and bf143cf46 after two reviews: MFront writes the static variables
  `TFEL_STATIC_VARIABLE(v)` (the variable itself on the host, a copy in
  the GPU code), in the code blocks, the inputs of the material
  properties, the behaviour variables and the Runge-Kutta algorithms
  (including a static lambda that kept the name of the first behaviour).
  Test behaviour `StrainHardeningCreep3` (comparisons openmp 10/10,
  cuda-nvcc 11/11), reproducer `reproducers/mfront-static-quantity-cuda`.
- Measurement pitfall: with `scripts/env.sh` (`OMP_PROC_BIND=spread`),
  `mpirun --bind-to none taskset -c ...` puts all the processes on the
  same core. The CPU measurements of ex8 of 11/10 were not affected; the
  difference with 19.9 s came from the element by element computation. An
  abandoned measurement script, left running in the background, also
  skewed a series: measurements redone.

### 2026-10-12 (night), fixes of MFront found by the reviews

- TFEL c9dc10510: Poisson ratio of the isotropic DSLs with the stress
  update algorithm converted to `real` (did not compile with
  `@UseQt true`), test behaviour `Plasticity8` (same results as
  `Plasticity6`).
- TFEL 4180d14b2: DDIF2 brick, non-constant fracture stresses and
  energies (inverted test, first energy used everywhere); checked by
  pairs of equal constant and non-constant behaviours, under MTest in
  each direction.
- TFEL fb5c52bc5: rkCastem, normalisation by a static variable written
  with the macro, debug message.
- Review Z6 handled: 58890e6f5 (test of `Plasticity8` that really depends
  on the Poisson ratio, mutation detected), 6d11373d0 (DDIF2: a single
  option of slope or energy, sign of the constant components), 31c5008d2
  (rkCastem: `young` used as normalisation factor, as announced by the
  release notes of 4.1; only `OrthotropicCreep_rkCastem` affected, same
  value; error messages put back in order).

### 2026-10-12 (night), hypre on the host with the device on the GPU

- mfem-mgis 279b5583d4 and 3ad387b002 (review Z7): option
  `--hypre-on-host` of `initialize` (and variable
  `MFEM_MGIS_HYPRE_ON_HOST`): behaviours, kernels and MFEM on the GPU,
  hypre on the host. The request is read before the initialisation of
  hypre; the flag of MFEM and the host policy are both needed
  (mutations). GPU tests 18/18.
- Measurements (`results/hypre-host-2026-10-12`): complete ex7 6.7 ->
  4.4 s on the GPU (CPU 2 cores 6.3 s), ex8 12.5 -> 11.6 s (heat transfer
  4.0 -> 2.8 s, as with MUMPS). For small systems, the GPU gets ahead of
  the CPU.

### 2026-10-12 (night), validations

- `gpu-proto` of mfem-mgis (3ad387b002): CPU suite 182/182, with forced
  batches and kernels 182/182, on the GPU 180/182 (only the two expected
  failures).
- TFEL: the 5905 tests of the directory `mfront/tests/behaviours/generic`
  pass with the behaviours regenerated by the MFront of the night (static
  variables, DDIF2, Poisson ratio of the isotropic DSLs, Runge-Kutta), and
  the whole TFEL suite (6970 tests, CPU) passes.
- mfem-mgis suite on the GPU with the forced kernels: 179/182, the two
  expected failures and `FixedPointSolverTest`, whose acceleration
  `Crossed2Deltabis` sometimes stagnates around 1e-10 on one step (1 time
  in 20 in isolated runs). The sums of the GPU are not reproducible; the
  algorithm is the one of MTest, faithfully ported; the test always
  passes on CPU.

### 2026-10-12 (morning), reuse of the MUMPS analysis

- On the GPU, ex8 spends 60 % of its time in the MUMPS factorisations of
  the mechanics (on the host), the Robin condition less than 5 %.
- mfem-mgis (branch mumps-reuse): parameter `ReuseAnalysis` of
  `MUMPSSolver`. The analysis of the previous matrix is reused if the
  pattern of the new one is the same on all processes; mfem-mgis compares
  the pattern, which `SetReorderingReuse` of MFEM does not check
  (c020804d68). Review Z8 handled (a114b6c870): solver in an internal
  header, pattern compared without allocation before the reduction and
  forgotten during the factorisation, number of analyses given by
  `getNumberOfMUMPSAnalyses`. Test `MUMPSReuseAnalysisTest`: pattern
  kept, changed everywhere, changed on the last process, or changed only
  by the off-diagonal columns or their global indices, number of analyses
  checked after each matrix; seven mutations detected.
  mfem-mgis-examples: option `-mra` of ex8.
- Measurements (`results/mumps-analysis-2026-10-12`): on the elasticity
  matrix of ex8, factorisations −20 % on one process, −28 % on two;
  orderings and approximate BLR factorisation without gain. ex8 on the GPU
  11.1 -> 9.7 s, on CPU (2 processes) 20.5 -> 19.2 s, same results and
  same iterations.

### 2026-10-12 (morning), cuDSS on the GPU

- mfem-mgis (branch cudss, 840683d96f): solver `CuDSSSolver` of the
  factory (parameters `Symmetric`, `PositiveDefinite`, `ReuseAnalysis`,
  refused without a CUDA GPU). The wrapper `CuDSSDirectSolver` recreates
  the MFEM solver when the size or the pattern changes (MFEM only chooses
  the reuse before the first matrix and stops on a different size).
  Comparison of the patterns shared with MUMPS (`HypreMatrixPattern`);
  `getNumberOfAnalyses` covers both solvers; test
  `DirectSolverReuseAnalysisTest` (option `--solver`), 6/6 on CPU and
  GPU. mfem-mgis-examples efeb8b8: `-lsMc CuDSSSolver` in ex8.
- Measurements (`results/cudss-2026-10-12`): on the elasticity matrix of
  ex8, factorisation 0.018 s with cuDSS (analysis reused) against 0.203 s
  with MUMPS on the host; ex8 on the GPU 9.8 -> 4.4 s (11.1 s at the
  beginning of the morning), same results and same iterations.
- The distributed analysis of cuDSS 0.7.1 fails for some matrices (rows
  reduced to their diagonal, tridiagonal coupling the processes, uneven
  distributions; several processes on one GPU):
  `reproducers/cudss-distributed`. One process is not affected.
- Review Z9 handled (45afaa0653): collective decision on the layout and
  the pattern (deadlock and stale descriptors when a process kept its
  local size), `HypreParMatrix` required, refusal with hypre on the host
  except for a GPU that reads the pageable memory of the host, solutions
  checked by their residual (cuDSS returns a finite and wrong solution
  for a negative definite matrix declared positive definite), pattern of
  a matrix on the GPU compared without its values. Test with 10 variants
  (1, 2 and 3 processes, hypre on the host), six mutations detected. ex8
  on the GPU: 4.4 to 4.5 s with hypre on the host, 5.6 s with hypre on
  the GPU.
- The heat transfer of ex8 also by cuDSS (`-lsTh CuDSSSolver`, option
  already there): its 35 solves by GMRES and BoomerAMG (1.83 s) become 35
  factorisations of 6 ms and solves of 1 ms. ex8 on the GPU: 2.6 s
  (11.1 s at the beginning of the morning, ×4.3), same results. Example of
  a run on the GPU in the README of ex8 (mfem-mgis-examples 7578d86).
- ex7: choice of the linear solver (`-ls`, `-ms`, `-mra`;
  mfem-mgis-examples e98e259). On the GPU, cuDSS 2.0 s against 4.3 to
  4.6 s with GMRES and BoomerAMG; on CPU (2 processes), MUMPS 5.4 s
  against 6.1 s; same iterations and same comparison to the reference.
- Beam r4 (111,843 unknowns) with cuDSS (options `-ls cudss`, `-mra`,
  `-els cudss` of the case, `results/cudss-beam-2026-10-12`): the direct
  solver does not pay off at this size. Newton with cuDSS 9.6 s against
  3.1 s with PCG and BoomerAMG (factorisation 0.235 s against 68 ms per
  solve); NewtonK0 with K0 factorised once by cuDSS 2.7 to 3.0 s against
  1.85 s with one BoomerAMG cycle (67 ms against 22 ms per tangent
  solve). cuDSS suits small systems (ex7, ex8).
- Review Z10 of the examples and of the beam case handled (no serious
  defect): mfem-mgis-examples 4ed849e (long names of the options of ex8
  without "mumps", options validated right after reading with messages
  that name them, warnings for the options without effect, README of ex7
  and ex8, complete GPU command of ex8 and comparable reference, test
  `mox2-mumps`); beam case: `-mra` passed to the elastic operator
  (`-els mumps|cudss`), symmetries described in the help, warning without
  a direct solver, complete list of the solvers in `cases/README.md`.
  Tests of ex7 and ex8 on CPU 6/6.
- Polycrystal of mm-opera-hpc (20,553 unknowns, periodic cell,
  non-symmetric) with cuDSS (`--linear-solver CuDSSSolver --reuse-analysis`,
  mm-opera-hpc 556b3e5, `results/cudss-polycrystal-2026-10-12`): 19.1 s
  against 11.4 s with PCG and BoomerAMG and 8.8 s with NewtonK0 and
  Eisenstat-Walker. Factorisation of 132 ms (reuse effective: 222 ms
  without), against 24 ms for the matrix of ex8 of similar size: the
  periodicity and the non-symmetry cost the direct solver a lot.

### 2026-10-12 and 13, ordered assembly of the tangent

- On the GPU, the computations were not reproducible from one run to the
  next (`FixedPointSolverTest` failed about once in 20). The source is the
  assembly by kernels: the element matrices are added to the sparse
  matrix by atomic additions, whose order changes at each run. The
  subtest without assembly diverged for the same reason: the elastic
  operator of its preconditioner is assembled by kernels; the matrix free
  action and the reduction of MFEM, for their part, are deterministic.
- mfem-mgis 0dbae0fcc7, 456c298dd6 and f57072220f: parameter
  `TangentOperatorAssembly` (`Atomic` by default, `Ordered`), variable
  `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`. With a device, the element
  matrices are stored, then added by one thread per value of the matrix,
  in a fixed order, from the list of the contributions of each value
  (computed once per row layout, checked column by column, recomputed if
  it no longer matches the matrix). Without a device, the host path,
  already sequential, is kept. Applies to the tangent, to the elastic
  operator of the fixed point, to FBar and to the Faltus regularisation.
- Result: four runs of `FixedPointSolverTest` on the GPU give histories
  identical bit for bit; three runs in atomic mode give three different
  ones. New tests: `GradientsKernelsTest` (assembly value by value
  compared to the atomic one, identical at each assembly, other row
  layout detected, contributions computed once per layout),
  `BatchIntegrationTest` (ordered cases, two computations identical bit
  for bit, number of ordered assemblies checked),
  `FixedPointSolverTest-kernels-ordered-cpu` and `-cuda`.
- First version with one thread per row: 7.1 ms per assembly for this
  kernel against 0.48 ms for the atomic one (Nsight Systems profile).
  Version per value: tangent of the beam r4 in 6 ms against 4.3 ms in
  atomic mode, plus 0.13 s at the first assembly; NewtonK0 2.2 to 2.35 s
  against 1.9 s. On CPU, same cost as the atomic one
  (`results/ordered-assembly-2026-10-12`).

### 2026-10-13, measurements of the pressure and of the post-processings, ordered GPU suite

- Imposed pressure on the beam r4 (work package L7.6,
  `results/pressure-beam-2026-10-13-1030`): on the GPU, with pressure,
  17.6 s before this work package (boundary integrator on the host, tangent
  assembled by MFEM on the host) against 1.70 s after; on 8 cores, the
  extra cost of the pressure disappears (9.88 against 5.68 s).
- Post-processings of the polycrystal (work package L8,
  `results/L8-polycrystal-2026-10-13-1027`): 1.5 s out of 11 s on the
  GPU, without gain by kernels (the ParaView export dominates), against
  ×2.8 on the total time of the bubble (projections).
- mfem-mgis suite on the GPU with the ordered assembly
  (`scripts/mfem-mgis-gpu-suite.sh`): 184/186 with the forced kernels,
  only the two tests that require the absence of a device failing.
- ssna303 series redone (`results/ssna303-ssna303-2026-10-13-0924`):
  refinement 1, GPU ×2.3 (FBar ×2.7); refinement 2 out of the GPU memory.

### 2026-10-14, float on the GPU: diagnosis of the divergence of the beam

- On the GPU in float, at the criterion 1e-4, the beam r4 diverged in one
  computation out of three to one out of two with atomic assembly
  (`results/precision-beam-gpu-2026-10-14-0604`). With ordered assembly at
  n = 10, all the computations by kernels with an iterative solver
  diverged (`results/diagnostic-float-gpu-2026-10-14`). With PCG and
  BoomerAMG, only the element by element tangent (assembled on the host)
  converged.
- Cause: the contractions into fused multiply-adds by nvcc in the
  stiffness kernels. The tangents by kernels and element by element only
  differ by 2e-7. But float works here close to its limit: at step 8, the
  Newton direction is only known to within 16 to 32 %. With the
  contractions, the Newton probably drifted along the soft modes of the
  beam, which the criterion on the residual does not see.
- Without contraction, the rigid translations are better respected
  (4.2e-8 against 5.1e-8 of the maximal term). The true linear residual is
  between 5e-3 and 9e-3 on the contracted trajectory (K), and 5e-4 on the
  trajectory without contraction.
- Fix in mfem-mgis: in float, stiffness kernels compiled for the GPU
  without contraction (`mfem_mgis_disable_device_floating_point_contraction`),
  same assembly cost; contractions kept in double (assembly 13 to 16 %
  slower without them). 31 computations out of 31 converge, within 4.2e-5
  at most of double; test `beam-gpu-r4` in cases/beam. Measurement redone
  (`results/precision-beam-gpu-2026-10-14-0843`): float 1.37 times faster
  than double, within 3.4e-5 at most of double.

### 2026-10-14, cases in float, polycrystal, site and CI on the CPU

- Suite of the cases of the prototype in float (root repository
  5d91ce18): 134 tests green, double 143/143. An executable of the
  examples or of mm-opera-hpc is only taken if it is built against the
  same installation of mfem-mgis as the cases.
- Polycrystal in float (`results/precision-polycrystal-2026-10-14`): axial
  stress at +6.1e-3 from double at 5e-5 of strain, +1.0e-3 at 2e-4; on a
  step of 1 s, the behaviour fails. Probable cause: F, stored in float
  close to the identity, only moves by one ulp at a time (argument for an
  interface in H = F - I).
- Draft of the site (`site/`, Sphinx and MyST, pages in English).
- `cpu` job of the CI simulated locally in an empty environment: the
  behaviour library of the tests of mfem-mgis was not built (126 failures
  out of 208), fixed. Then a matrix on the precision
  (`deps/spack/ci-cpu-float`), simulated in double and in float.
- The simulation in float found a tolerance too tight in mm-opera-hpc: at
  1e4 Pa, the fixed point on the lateral stresses of the first step of the
  cermet stagnated around 1.1e4 Pa or oscillated between 1.5e4 and
  1.7e4 Pa (five computations out of six, depending on the rounding of
  the whole computation, behaviour included). When the two lateral
  stresses are close to S, the correction of the lateral components of F,
  about (1 + ν) (1 - 2 ν) S / E, rounds to zero below half the gap between
  two neighbouring floats (6e-8 just below 1), that is 1.54e4 Pa of
  residual. The cut step gave a 2.7 % difference to double at 0.2 s, an
  effect of the time step (0.24 to 0.27 % at the same steps).
- Three successive fixes, each one caught by a review: fdfdafd reduced
  the check in float to 0.4 s, on a wrong explanation (the rounding of
  F); 2e4881c took 5e4 Pa from an underestimated floor and measurements of
  mixed provenance; 19d64ba took 3e4 Pa, a constant tied to UO2 and to
  tension, with an explanation of the difference wrong by a factor of ten
  (the rounding of 1 + 1e-4 in float only makes +0.017 % on the strain).
  Sweep redone with a single executable, with one and two processes
  (`results/cermet-float-2026-10-14`). Final rule (05358bd, then ae5768b
  for the gap between floats taken on the smaller side): the tolerance
  stays 1e4 Pa, and the fixed point takes at least twice its resolution,
  computed from F and the moduli (30,804 Pa here, twice as much in
  compression, negligible in double). Difference to double from 0.10 to
  0.16 % at 0.2 s, of which about 0.07 % is due to the absolute tolerance
  of the microscopic equilibrium in float (measured in double with this
  tolerance), and from 0.06 to 0.07 % at 0.4 s; same checks in both
  precisions; double unchanged to the bit; polycrystal unchanged at 0.1
  and 0.4 s, still with a failure of the behaviour with a step of 1 s
  (results/precision-polycrystal-2026-10-14).

### 2026-10-14, ParaSails in float and Newton criterion at the floor of the residual

- ParaSails in float: a defect of hypre 3.2.0 (also on master, not yet
  reported). With an external LAPACK in single precision, hypre calls the
  LAPACK routines of the double type on float arrays: the LAPACK header
  has no single precision names, unlike the BLAS header. Two fixes
  compared: a fixed header, and the internal LAPACK of hypre, both
  correct, same times. Chosen: `hypre ~lapack` in the float environments
  (`mmm-float-cpu` rebuilt). mfem-mgis 3ba3b6e729: warning at the
  configuration, ParaSails tests in float. Draft of the upstream report
  ready, not published.
- Newton criterion adapted to the precision: mfem-mgis 11fda6426e, option
  `StagnationRelativeTolerance` (10⁴ ε by default). In float, the beam r4
  converges with the criterion of the case (1e-8) in 27 iterations, 7.0e-5
  away from double; in double, nothing changes (mm-opera-hpc bit for bit,
  no stagnation in the test suites).
- CI simulated in double: mfem-mgis 208/208, examples 22/22, mm-opera-hpc
  6/6.

### 2026-10-14, review of the Newton criterion, translation into English

- Independent review of the stagnation criterion (no grave defect, three
  medium ones) handled in mfem-mgis 9c2dc8c4e5: the residual must not
  decrease by more than 10 % during two iterations (the first version
  compared each iteration with the smallest residual, and accepted two
  iterations with a factor 0.6); the best iterate of the resolution is
  kept (the residual at the floor is noisy, from 615 to 7594 on a step of
  the finite strain test in float); the test with the small strain
  behaviour is removed (it converged without the criterion), and the
  finite strain tests check the messages; `cmake_push_check_state` for the
  check of hypre, whose list of affected solvers is completed (Gaussian
  elimination 98 and 198, GSMG, AME). Float suite of mfem-mgis on the CPU:
  223/223 (91 stagnations accepted, 53 best iterates restored). Double
  unchanged: CI simulation 208/208, examples 22/22, mm-opera-hpc 6/6 with
  the same results bit for bit, no stagnation in the logs. Beam r4 in float
  with the criterion of the case (1e-8): same iterations and solution as
  with the first version (27 iterations, 7.0e-5 away from double, 30.3 s;
  `results/newton-stagnation-beam-2026-10-14`).
- Translation of the repository into English, before it becomes public:
  files and directories renamed (`cas/` became `cases/`, 1,150 files),
  contents translated (documents, log, results, scripts,
  CMake, test names, environment variables of the scripts, names of the
  local build directories and of the local branches of mfem-mgis). The
  commit messages of the history and the French text of upstream files
  (licence of MGIS, pages of TFEL) are kept. An inverted statement was
  fixed on the way (`docs/gpu-solvers.md`: 8 bytes per coefficient in single
  precision instead of 12). Checks: cases in float 101/101, in double
  105/105 on the CPU and 41/41 on the GPU, examples 21/21, mm-opera-hpc 6/6,
  site built without warning, no broken link; independent review (no
  grave defect) handled.
- Float chain on the GPU after the rebuild of `mmm-float-gpu` (hypre with
  its internal LAPACK): mfem-mgis 280/280 on the GPU, with the ParaSails
  and stagnation tests, after a fix of a test (mfem-mgis 73055c8ca6: the
  numbers of iterations of BatchIntegrationTest may differ by two in
  float; on the GPU, a time step converged once in 1 iteration instead of
  3, one run out of ten); examples 21/21, mm-opera-hpc 6/6, cases 134/134.

### 2026-10-14 (evening), defects of hypre and of ROCm reported

- hypre: issue [hypre-space/hypre#1634](https://github.com/hypre-space/hypre/issues/1634)
  and fix [#1635](https://github.com/hypre-space/hypre/pull/1635) (single
  precision names of the LAPACK routines). Checked on the master branch of
  hypre with OpenBLAS: the `ij` driver (ParaSails and PCG, 100³ unknowns, two
  processes) crashes in `ParaSailsSetupValues` without the fix and converges
  in 211 iterations with it. On master, single precision with an external
  BLAS does not even compile with GCC 14: `ParaSails/ConjGrad.c` calls
  `dcopy_`, `dscal_` and `daxpy_` without prototype (spack patches it); the
  same pull request fixes it.
- ROCm on the integrated GPU gfx1036: the diagnosis of 13/10 (a kernel lost
  at the first launch of a kernel of a shared library) was wrong. The
  kernel runs (7 to 15 µs after its launch), but its writes to device
  memory are lost when the GPU stays idle about 300 µs before the next
  operation: 0 trials out of 20 with 100 µs, 20 out of 20 from 300 µs on;
  none with a synchronisation right after the launch, with
  `AMD_SERIALIZE_KERNEL=3` or with pinned host buffers; rare when a host
  thread keeps the GPU busy. Loading the code of libmfem (2 ms) and large
  allocations give this idle time. Plain HIP reproducer in
  `reproducers/hip-idle-write-loss` (renamed from `hip-lost-kernel`),
  reported as [ROCm/rocm-systems#13469](https://github.com/ROCm/rocm-systems/issues/13469).
  It looks like the modified lines of the L2 cache are dropped when the GPU
  powers down on idle; not checked with GFXOFF disabled (reboot needed).
  No fix on our side beyond the serialisation of the tests.

### 2026-10-14 (night), displacement gradient for the finite strain behaviours

- Decided on 14/10: the solver may give H = F - I instead of F to the
  finite strain behaviours (`docs/displacement-gradient.md`). Implemented
  in TFEL (generic interfaces, Green-Lagrange strain from H, logarithmic
  strain handler built from F and E), MGIS (option `gradient`, bindings),
  mfem-mgis (parameter `DisplacementGradient`) and mm-opera-hpc (option
  `--displacement-gradient`, macroscopic unknowns stored as H).
- Three independent reviews. Found and fixed: Hencky strain in float with
  two close eigenvalues (57 % of error on the first Piola-Kirchhoff
  stress, tolerance of 1e-14 applied to the differences of eigenvalues),
  plane stress in the Hencky framework (NaN), options never given to the
  initialize functions and post-processings by MGIS (and a null `K` array
  in the material data managers), macroscopic loading of mm-opera-hpc
  still rounded as F.
- Found on the way, in TFEL as it was: in float, the logarithmic strain
  handler built from F loses all precision when two eigenvalues are close
  (rotated uniaxial strain of 10 %, stress not coaxial with the strain:
  1e-1 of error on the second Piola-Kirchhoff stress, 4.6e2 on its
  derivative). Fixed by the same divided differences with `log1p`
  (below 4e-6), the special cases for equal eigenvalues being removed.
- Results: `results/displacement-gradient-cube-2026-10-14` and
  `results/displacement-gradient-opera-2026-10-14`. Test suites: MGIS
  50/50 in double and in float, with the bindings 99/99, mfem-mgis 227/227
  in double and 230/230 in float, mm-opera-hpc 6/6 in double and in
  float.

### 2026-10-14 (night), CUDA and HIP levels of the continuous integration

- The jobs `cuda` and `hip` of `.github/workflows/ci.yml` (renamed from
  `cpu.yml`) compile the whole chain with CUDA 13.0 and with ROCm 10.1 at
  each run, with the build cache of spack on GHCR and ccache, as the job
  `cpu`. Simulated locally (`scratch/ci-sim/run-gpu.sh`): both green.
- First simulation of `cuda`: two defects. The environment `ci-cuda` did
  not use the repository `mmm_gpu`, hence not its patch of MFEM for CUDA 13
  (cicc never ends on `bilininteg_convection_ea.cpp`); and ptxas needs more
  than 12 GB for `fem/lor/lor_batched.cpp` with sm_120, even with `-O0`,
  more than the memory of the machines of GitHub. Fixed by the repository
  `mmm_gpu` in `ci-cuda`, a variant `cuda_ptx_only` of its package `mfem`
  (PTX code only: ptxas is not run, enough for a CI which only compiles),
  and one job for hypre and MFEM (cicc needs up to 8.8 GB for
  `lor_batched.cpp`, 12.5 GB with two jobs).
- Found by this job: the link of `HypreDevicePoolSizeTest` of mfem-mgis
  failed with a hypre built with Umpire and an MFEM built without it
  (mfem-mgis aa045a7853: Umpire linked explicitly).

### 2026-10-14 (late night), plane stress of the finite strain behaviours

- Found while testing the displacement gradient, and older than it: in
  plane stress, TFEL expected an axial component of zero in the
  deformation gradient given by the caller, whereas MGIS, mfem-mgis and
  MTest give one. The behaviours written with F and in the Green-Lagrange
  strain framework added the axial component they compute to the given
  one: with one, the Piola-Kirchhoff stresses of the former and the
  Cauchy stress of the latter were wrong by a factor of two. The Hencky
  strain framework used the given axial component in its conversions:
  the first and second Piola-Kirchhoff stresses were wrong by the axial
  strain, in relative value (with F as with H; mfem-mgis uses the first
  Piola-Kirchhoff stress).
- Fixed in TFEL: the axial component given by the caller is not used,
  the computed one is set in the deformation gradient and used by all
  the conversions, and the axial strain given in the Green-Lagrange and
  Hencky frameworks is zero (an independent review found that the
  Green-Lagrange framework still computed it from the given component).
  Test in MGIS against the analytical solution (Cauchy and first
  Piola-Kirchhoff stresses, F and H, axial component shifted by the
  caller), which detects these defects before the fix. Not fixed, older:
  the tangent operator DSIG_DF in plane stress misses the derivative of
  the axial component (1e-2 at 3 % of strain; DPK1_DF is exact).
- Test of an orthotropic finite strain behaviour given the displacement
  gradient in mfem-mgis (the gradients are rotated in the material frame,
  H as F).
- Audit of the whole history of the repository before it becomes public:
  no secret (tokens, keys) and no internal link.
- Initial value of the axial deformation gradient in plane stress: the
  behaviours written with F store it in the `AxialDeformationGradient`
  internal state variable, which TFEL now uses in the conversions. MGIS
  and mfem-mgis initialized it to zero, as the other internal state
  variables: the deformation gradient was singular at the first step
  (Piola-Kirchhoff stresses converted to NaN, Newton on this component
  divided by zero). MGIS (`getAxialDeformationGradientOffset`, also
  bound in C, Fortran and Python) and mfem-mgis (`MaterialArrays`) set
  it to one. Test behaviours in plane stress written with F, whose
  Newton starts from this value; the tests fail without the fix. An
  independent review found three gaps in the tests (offset zero, branch
  of the initializer, layout component by component), filled. In plane
  stress, MFront converts no tangent operator of the behaviours written
  with F: a behaviour used by mfem-mgis must compute `DPK1_DF`.
- Axisymmetrical generalised plane stress (axial component second, not
  third) added to the test behaviours of MGIS: analytical solution with
  F and with H in the three frameworks, integration with F and with H for
  all the stress measures and tangent operators, initialize functions and
  post-processings. TFEL was right in this hypothesis (never tested
  before).

### 2026-10-15, linear solvers of the major libraries: phase 0

- Plan [docs/plan-solver-libraries.md](docs/plan-solver-libraries.md):
  the whole coverage of what MFEM 4.10 wraps (MFEM, hypre, PETSc, Ginkgo,
  MUMPS, SuperLU_DIST, STRUMPACK, cuDSS, SuiteSparse, oneMKL, AmgX,
  KINSOL), plus MAGMA and Trilinos, PaStiX through PETSc.
- Phase 0 in mfem-mgis: each linear solver and preconditioner is
  registered with a description (library, kind, location, accepted
  preconditioners, options); a factory of preconditioners replaces the
  chain of tests on their names; the location is checked from the
  description instead of the name; one module per library (MFEM, hypre,
  MUMPS, cuDSS, AmgX); the other libraries report their convergence
  through an interface. Moved without change of behaviour (independent
  review); placements locked by the generic test, which solves the
  Laplace equation with every entry, in sequential, parallel and on the
  GPU. Suites: 228 (GPU MFEM, on the host), 233 (CPU double), 232 (CPU
  float); Doxygen and Sphinx without warnings.

### 2026-10-15, linear solvers of the major libraries: PETSc

- Phase 1 of the plan, first library: `PETScLinearSolver` (a Krylov
  solver of PETSc, `KSP`) and `PETScPreconditioner` (`PC`), configured by
  the options of PETSc given by their `Options` parameter. The options are
  written in the options database with a prefix unique to each solver and
  removed when it is destroyed; those which PETSc never read are reported.
  A preconditioner of another library is a shell preconditioner of PETSc,
  which receives the operators; the options of the preconditioner of PETSc
  are then refused in the options of the solver, since PETSc would replace
  it. PETSc is initialized when the first solver is built.
- The interface of the solvers of the other libraries now gives access to
  their relative tolerance, used by the adaptive tolerance of the Newton
  solver.
- Dedicated test (sequential and parallel): on a 16³ mesh in parallel,
  CG with Jacobi 25 iterations, with GAMG 9, with BoomerAMG as shell
  preconditioner 11, CG of MFEM with the GAMG of PETSc 10. Independent
  review: no defect of severity 4 or 5, its points handled. Suites: 230
  (GPU MFEM, on the host), targeted CPU double and float; Doxygen and
  Sphinx without warnings.

### 2026-10-15, linear solvers of the major libraries: Ginkgo

- Each solver and preconditioner of Ginkgo is a configuration of Ginkgo
  (`gko::config::parse`), completed by the `Configuration` option of the
  entry, on the classes of MFEM wrapping Ginkgo. Ginkgo computes where the
  device of MFEM computes. Solvers: CG, FCG, BiCGSTAB, CGS, GMRES, GCR,
  MINRES, IR in sequential and parallel; IDR, PipeCG, BiCG, CB-GMRES and
  the sparse direct solver in sequential (the distributed versions of IDR
  and PipeCG fail with the vectors wrapped by MFEM); `GinkgoSolver` for
  any iterative solver. Preconditioners: Jacobi, ILU, IC, their ISAI
  variants, SOR, Gauss-Seidel, direct solver (local, in the Schwarz
  preconditioner of Ginkgo in parallel, optional L1 smoother), algebraic
  multigrid (Pgm), `GinkgoPreconditioner` for any preconditioner.
- A preconditioner of another library is wrapped by MFEM and receives the
  operators. The stopping criteria use the norm of the residual (MFEM's
  wrappers used preconditioned or shadow norms), the relative tolerance
  can be changed. Each key of a configuration is checked alone: Ginkgo
  1.11 terminates the program on some errors in a configuration whose
  other keys are not read yet.
- All the entries pass the generic test and the dedicated test on the CPU
  and on the GPU (RTX 5080, CUDA), sequential and with 2 processes; on a
  12³ mesh, CG with the multigrid of Ginkgo 16 iterations against 31 with
  its block Jacobi. Independent review: its points handled (one 4/5 in
  the test in single precision; 3/5: preconditioner given twice, implicit
  criteria, termination on configuration errors). Suites: 232 (GPU MFEM,
  on the host), targeted CPU double and float; Doxygen and Sphinx without
  warnings.

### 2026-10-15, linear solvers of the major libraries: hypre completions, SuperLU_DIST and STRUMPACK

- hypre: the BiCGSTAB, COGMRES (host only) and LGMRES Krylov solvers,
  the hybrid solver (diagonal scaling, then BoomerAMG when the
  convergence is slow) and the factorized sparse approximate inverse,
  which MFEM does not wrap, are implemented after `mfem::HypreGMRES`;
  the smoothers of hypre applied by MFEM (`HypreSmoother`); new options
  of BoomerAMG with the ranges of hypre (truncation, aggressive
  interpolation, size of the coarsest level, Chebyshev order, ...).
  The convergence of the Krylov solvers is given by the indicator of
  hypre and by the code of the resolution: a breakdown of BiCGSTAB is
  not a convergence. The hybrid solver only computes its final residual
  when its logging is enabled: without it, every resolution was reported
  as converged.
- SuperLU_DIST: the matrix is given to the solver of MFEM just before
  the resolution (a solver of MFEM destroyed with an operator never
  factorized is avoided); the parallel symbolic factorization is refused
  on one process and with the reuse of the analysis, where SuperLU_DIST
  9.2.1 crashes; ParMETIS with the sequential symbolic factorization
  works in both cases.
- STRUMPACK: options given by its command line, double and mixed
  precision, compression of the frontal matrices with a Krylov solver of
  STRUMPACK, approximate factorization as a preconditioner. STRUMPACK
  reports a success when its Krylov solver stops at the maximum number of
  iterations (relative residual of 2e-2 measured): mfem-mgis reports it
  as not converged.
- Each new entry solves a second operator of the same pattern (the
  Newton path). Found on the way: on a GPU, `mfem::Add` sorts the columns
  of each row, while the diagonal scaling of hypre reads the diagonal as
  the first entry: NaN, then a crash in MPI with 2 processes (upstream
  note written, nothing posted).
- Independent review: one 5/5 (hybrid convergence) and the 3/5 handled
  during the integration. Suites: 236 (GPU MFEM, on the host), targeted
  CPU double (72) and float (71); Doxygen and Sphinx without warnings.

### 2026-10-15, linear solvers of the major libraries: MFEM completions and MUMPS

- The Gauss-Seidel and Jacobi smoothers of the sparse matrices of MFEM,
  its block ILU(0) factorization and the direct solvers UMFPACK and KLU
  are preconditioners, applied to the diagonal block of each process in
  parallel; KLU is also a direct solver in sequential; `OrthoSolver`
  orthogonalizes the vectors to the constant vector (scalar problems with
  null flux conditions). MUMPS accepts a reordering strategy, the
  tolerance of its block low-rank factorization and a print level.
- With `--device cuda`, `mfem::KLUSolver` and `mfem::BlockILU` gave wrong
  results (residual 0.9): they read the memory of the host without
  synchronizing it, while the vectors live on the GPU, and KLU reads the
  arrays of the host after sorting the columns on the GPU. They are given
  a copy of the matrix on the host only and copies of the vectors. MFEM
  does not check the status of KLU: a singular matrix now fails on all the
  processes.
- Independent review: two 3/5 (the two points above), all the points
  handled (`review-mfem-mumps.md`). The generic test now runs these host
  entries when hypre computes on the GPU. Suites: 238 (GPU MFEM, on the
  host), targeted CPU double (74) and float (73).

### 2026-10-15, non linear solver of KINSOL (end of phase 1)

- The `KINSOL` non linear solver uses KINSOL (SUNDIALS 7.9) through
  `mfem::KINSolver`, whose functions given to KINSOL are replaced by the
  ones of mfem-mgis: they report the failures of the integration of the
  behaviours and of the linear solver to KINSOL (MFEM's functions always
  report a success), and give the product by the tangent operator
  (otherwise approximated by differences of residuals, each one
  integrating the behaviours).
- On the uniaxial tension of a plastic cube: Newton's method of KINSOL,
  with or without line search, makes the 8 iterations of the Newton
  solver and gives the same solution; modified Newton 14, Picard 59, with
  Anderson's acceleration 32. Sequential, parallel, CPU and GPU.
- Found on the way: with `--device cuda`, MFEM wraps CUDA vectors of
  SUNDIALS, and a direct solver which writes on the host left the array of
  the GPU stale (null steps); KINSOL 7.9 asks for a setup at each
  iteration of Picard's method whatever its maximum number of setup calls;
  a null tolerance becomes KINSOL's default 6e-6, and the default maximum
  Newton step, 1000 times the norm of the initial guess, applies to
  unknowns which are not scaled. All handled; upstream notes written,
  nothing posted.
- Independent review: three 3/5, all handled (`review-kinsol.md`). The
  configuration of Doxygen now documents the modules of SuperLU_DIST,
  STRUMPACK, Ginkgo and SUNDIALS. Suites: 240 (GPU MFEM, on the host),
  targeted CPU double (74) and float (73).


### 2026-10-15, separate environment of the solver libraries (phase 2)

- The libraries which change the build of MFEM or of PETSc go into a new
  spack environment, `mmm-solvers`, so that the MFEM of `mmm-gpu` is not
  rebuilt before the merge of the pending upstream pull requests: a spack
  repository `mmm_solvers` derives the packages `mfem` (variants
  `mkl-pardiso`, `mkl-cpardiso`, `magma`) and `petsc` (variant `pastix`)
  from the ones of `mmm_gpu`, and adds `intel-oneapi-mkl-solvers`, oneMKL
  without the BLAS and LAPACK providers: in a unified environment, oneMKL
  would otherwise have replaced OpenBLAS for every package. PETSc with
  MUMPS requires a Fortran compiler and is not built. Only MFEM, PETSc,
  MUMPS, MAGMA, PaStiX, Scotch and oneMKL are built; the other packages
  are the ones of `mmm-gpu`. The build took 2 h 40 (MFEM alone more than 2 h, the
  CUDA kernels of `lor_batched` needing 23 GB with ptxas).
- `MKLPardisoSolver` (sequential) and `MKLClusterSparseSolver`
  (parallel) call oneMKL directly. The classes of MFEM check the errors
  of oneMKL with `MFEM_ASSERT` only, read the vectors on the host without
  synchronizing them, and give the whole matrix to the symmetric types,
  which only read the upper triangle. Found on the way: PARDISO perturbs
  the null pivots and reports a success on a singular matrix; a
  resolution is now reported as converged only if its residual is lower
  than 1e-3 times the right hand side, as for cuDSS.
- The test of PETSc solves with the external packages PETSc is built
  with: LU factorization of PaStiX (relative residual 1e-14), GCRODR of
  HPDDM, which needs 13 iterations instead of 25 for a second resolution
  thanks to the recycled subspace.
- Independent review of the oneMKL module: two 3/5, handled
  (`review-mkl.md`). Suite of the `mmm-solvers` build: 257 tests, all
  green; oneMKL also on the CUDA device, sequential and parallel.

### 2026-10-15, MAGMA: sparse solvers and nodal blocks (phase 3)

- The sparse part of MAGMA 2.10 (environment `mmm-solvers`) is wrapped by
  a module of mfem-mgis: seven Krylov solvers (sequential, one GPU) and
  the Jacobi, incomplete LU and incomplete Cholesky preconditioners of
  MAGMA, the preconditioners being applied to the diagonal block of each
  process in parallel and usable by any solver. MAGMA is called directly:
  its option parser prints its usage at each call.
- Defects of MAGMA found on the way: the preconditioned TFQMR reports its
  breakdowns as successes (convergence now checked on the residual
  recomputed by MAGMA); QMR applies the transposed incomplete Cholesky
  factors in the wrong order (refused, the incomplete LU factorization
  being the same for a symmetric matrix); QMR without look-ahead breaks
  down on the Laplace equation; a null diagonal entry goes undetected by
  the incomplete factorizations (now refused, on all the processes).
- Independent review: no finding rated 4 or 5, two rated 3 (the TFQMR
  breakdown, and a test which did not run the path giving the vectors of
  MFEM to MAGMA without copy), all the findings handled.
- The batched dense algebra of MFEM (native kernels, cuBLAS or hipBLAS,
  MAGMA) inverts the blocks of the new `NodalBlockJacobi` preconditioner:
  the 3 x 3 blocks which couple the components of the displacement at each
  node, extracted from the assembled matrix, or computed by kernels for the
  tangent operator applied without assembly (standard integrators, FBar
  formulation, regularization of Faltus et al., essential degrees of
  freedom). The nodal blocks computed without assembly match the ones of
  the assembled operator (relative difference 6e-17 at most over 120
  comparisons), and Newton's method with this preconditioner reproduces
  the reference solutions.
- Found on the way: the MAGMA backend of the batched algebra of MFEM
  computes on its own queue without waiting for the kernels of MFEM: the
  GMRES solver of MFEM diverged in parallel on the GPU; the device is now
  synchronized around the external backends. Also, on a GPU, reading on
  the device the diagonal block given by `HypreParMatrix::GetDiag` ended
  with an invalid `cudaFree`: the arrays of hypre are now read directly.
  Upstream notes written, nothing posted.

### 2026-10-15, end of phase 1: MUMPS, AmgX, MGR of hypre, SNES of PETSc

- The remaining items of phase 1 are done: the factorization of MUMPS as a
  preconditioner (`MUMPSPreconditioner`, the block low-rank factorization
  being checked on a cube large enough for MUMPS to compress fronts: 10
  iterations of the conjugate gradient instead of one), the aggregation
  multigrid of AmgX in the default configurations (`Multigrid` option),
  the generic test of all the linear solvers on a beam in linear
  elasticity (double and single precision, with the entries which do not
  converge on this problem listed with their reason), the multigrid
  reduction of hypre (`HypreMGR`) and the non linear solvers of PETSc in
  the factory of non linear solvers (`SNES`). `AMGFSolver`, `LORSolver`
  and the p-multigrid of MFEM are not applicable as is to the tangent
  operators of mfem-mgis: AMGF needs a filtered subspace (contact), LOR a
  high order bilinear form to refine, the p-multigrid the operators of each
  level.
- MGR: with the Jacobi F-relaxation of hypre (default), MGR is slower than
  the Jacobi preconditioner on the beam; with an F-relaxation by algebraic
  multigrid, GMRES needs about 35 iterations on the host and 50 on a GPU,
  against 150 to 300 with Jacobi. The points which are not kept at a level
  must be eliminated explicitly (`HYPRE_MGRSetNonCpointsToFpoints`): by
  default hypre added coarse points, and two levels needed 87 iterations
  on a GPU instead of 41 and did not converge in single precision. Several
  options of hypre 3.2.0 crash or give wrong results (relaxations of the
  host applied to the matrices on the GPU, block interpolations and
  restrictions): they are refused with a reason, and a workaround avoids
  the crash of the Jacobi relaxation with the block interpolation on a GPU.
- SNES: Newton's method with and without line search, the trust region
  method, the linear solver of the problem or the solvers of PETSc, the
  tangent operator assembled or applied without assembly, on the host and
  on a GPU; the failures of the integration are domain errors of SNES
  (the line search reduces its step). On the uniaxial tension of a plastic
  cube, SNES gives the solution of the Newton solver of mfem-mgis in the
  same number of iterations.
- Found on the way: the solvers of PETSc of the factory of linear solvers
  gave wrong results as soon as MFEM computes on a GPU (MFEM 4.10 sets an
  internal flag of the vectors of PETSc 3.26, which then saw null vectors),
  and the vectors of PETSc of the GPU made `PetscFinalize` abort the
  program: with a GPU, the solvers of PETSc now compute on the host. The
  test of the solvers of PETSc now also runs on a GPU.
- Independent reviews: the lot of the morning (no finding rated 4 or 5),
  MGR (one 4/5 reproduced, one 3/5, handled), SNES (one 4/5 already fixed
  by the GPU commit, the 3/5 and 2/5 handled). Upstream notes written for
  hypre, MFEM and PETSc, nothing posted.

### 2026-10-15, Trilinos (phase 4, in progress) and benchmark of the linear solvers

- Environment `mmm-trilinos` (deps/spack): Trilinos 17.2.1 with CUDA,
  Kokkos and Kokkos Kernels 5.2.1 (the versions Trilinos 17.2 requires),
  hence PETSc and MFEM built again. The `shylu` variant of the package of
  spack only enables FROSch: the local repository `mmm_trilinos` adds the
  variant `tacho` (Tacho in Amesos2; GoogleTest disabled, otherwise TriBITS
  stops). Trilinos took 1 h 29 to build with nvcc.
- mfem-mgis: module of Trilinos written and reviewed, not yet built in
  mfem-mgis. The calls to Trilinos are compiled by the compiler of Kokkos
  in a shared library, `MFEMMGISTrilinos`, whose interface only uses
  standard types; it compiles with nvcc without warning. Entries:
  `BelosSolver`, `Amesos2Solver`, `Ifpack2`, `MueLu` (rigid body modes),
  `FROSch` (two level Schwarz), `Amesos2Preconditioner`. The unknowns
  ordered by nodes are seen by Trilinos ordered by vector dimension, as
  MueLu expects. The static review found one defect rated 5 (the matrix of
  hypre left on the host after the conversion), one rated 4 (a test of
  MueLu on a system too small for a multigrid) and several rated 3 and 2
  (options of Amesos2 ignored, failures of Belos aborting the run, Dirichlet
  rows not seen by FROSch): all handled before the first build.
- Blocked: the MFEM of `mmm-trilinos` is not built. One of its files,
  `fem/lor/lor_batched.cpp`, needs about 22 GB for ptxas with sm_120; with
  another large compilation on the machine, the memory ran out twice, and
  the compilation was stopped. `install-mmm-trilinos.sh` builds MFEM alone
  with two jobs.
- mfem-mgis 8d829b4c26 (patch 0231): `LinearSolversBenchmark` measures any
  entry of the factory of linear solvers, its parameters being given on the
  command line in a small text syntax, on assembled systems (Laplace
  equation, beam in linear elasticity, stiff inclusions): setup, first and
  following resolutions, second operator of the same sparsity.
  `cases/measure-solvers.sh` runs the campaign of phase 5 with it, on 8
  processes and on the GPU, and `cases/summary-solvers.py` makes the tables.
