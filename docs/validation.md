# Validation

This document describes how the changes of TFEL, of MGIS, of mfem-mgis
and of the cases are checked, with which tools, and where to find the results.

## Principle

The `generic` interface serves as the reference. For each tested behaviour, the same
MFront source is compiled twice: with the `generic` interface by the
host compiler, and with the `generic-parallel` interface for the
tested configuration. The two libraries compute the same points, and
the results are compared. The accepted differences are rounding differences.

A test that cannot fail is worthless. The important tests are
therefore checked with a deliberately wrong version of the code (relative index
instead of absolute, wrong stride of a field, default stream instead of the stream of
the caller, etc.): the test must detect it.

## Environment

- Machine: NVIDIA RTX 5080 GPU (architecture `sm_120`); the AMD GPU
  integrated in the processor (gfx1036, access through the render group)
  only checks that HIP runs (see the end of "AMD GPU" below).
- Spack environment `mmm-gpu` (`deps/spack/`): NVIDIA HPC SDK 25.11
  (nvc++ and CUDA 13.0.2), LLVM 20.1.8, AdaptiveCpp 25.10, Kokkos 5.2.2,
  HIP 7.2.3 (NVIDIA platform), MFEM 4.10 with hypre 3.2.0, and the
  studied solvers (PETSc, SUNDIALS, SuperLU_DIST, STRUMPACK, Ginkgo,
  AmgX, cuDSS, MUMPS), OpenMPI 5. See `deps/spack/README.md`.
- Spack environment `mmm-raja`: RAJA 2026.07.0 and camp 2026.07.1, built
  with CUDA 13.0.2 (`sm_120`) and OpenMP, for the `raja` configuration of
  TFEL. It is separate from `mmm-gpu`, whose camp has no OpenMP.
- `scripts/env.sh` gives the complete environment of all the scripts.
- The heterogeneous memory management (HMM) of the NVIDIA drivers lets the
  GPU read the ordinary memory of the host. It therefore hides the location
  errors. The final validations are done without HMM
  (`uvm_disable_hmm=1`).

## TFEL

### Tests of the interface of the behaviours and of the models

`tfel/test-suite/behaviours/generic-parallel/<configuration>/`: one directory
per configuration (see [tfel.md](tfel.md#configurations-and-compilers)).
Each test compiles the test behaviours with the two interfaces, then
the tool `compare.cxx` compares them on varied batches:

- three types of computation (integration only, tangent operator, prediction)
  and two memory layouts;
- uniform values and values per point, given group by group or
  variable by variable;
- initialize functions, post-processings, rotations, models,
  parameters modified at run time, deliberate failures (index of the first
  failing point, statuses per point), asynchronous computation;
- report function of the library after each asynchronous computation,
  first unreliable point, initial time step factor of the batch (option
  `--time-step-scaling-factor`, case `0.5` in `compare-all.sh`);
- invalid file of parameters of a behaviour and of a behaviour variable
  (`-6`, option `--invalid-parameters`);
- refusal of the names containing `_generic_parallel`, for a behaviour and
  for a model (`check-refused-name.sh`);
- libraries built with the `generic_functions` option
  (`check-generic-functions.sh`, option `--generic-functions` of
  `compare.cxx`): interface `Generic`, functions of a point compared with
  the functions of batches, shared parameters and out of bounds policy.
  All the configurations run it, except `raja-options` and
  `sycl-default-acpp-2`;
- TDLS solver in the implicit behaviours (sub-test `tdls`).

### Tests of the interface of the material properties

`tfel/test-suite/material-properties/generic-parallel/<configuration>/`:
values compared with the expected values, uniform and strided arguments,
bounds and policies, properties calling other properties or defined
by data, execution options (non-blocking CUDA and HIP stream,
Kokkos instance, resource of camp on a stream created by the test for
RAJA, missing SYCL queue), asynchronous evaluation.

The tests of the semantics are written once (`tests.hxx`, data of a batch
built by `batch.hxx`) and shared by the backends whose arrays are
accessible by the host: statuses of the points, failures and unreliable
results in the same batch, invalid results, errors of the batch (null
pointers, wrong number of arguments, uniform output with several points,
more than 2^48 - 1 points, null statistics given to the report function),
empty batch, uniform and strided values. Each backend adds its own tests
(queues, asynchronous evaluation, stack, memory of the device). A
comparator (`compare.cxx`, run by `compare-all.sh`) evaluates each test
property with the `generic` interface at each point and with the batch,
for the three policies, contiguous, strided and uniform values,
synchronous and asynchronous evaluations: values, statuses of the points
and summary. The values are compared exactly with `g++`, `clang++` and
`nvc++` on the host, up to 1e-13 with `icpx` and on GPUs. The
property with bounds and physical bounds on two arguments and on its
output (`BoundsOrderTest`) checks the retained check, in `tests.hxx` and in
the comparisons. The property generated with the runtime checks disabled
is run by `test-no-runtime-checks.cxx` (memory, queue and synchronization
of `runtime.hxx`) on `stlpar-parunseq-gcc`, `cuda-nvcc`, `hip-hipcc`,
`sycl-default-icpx`, `sycl-default-acpp`, `openmp-host-gcc` and
`openacc-nvhpc-multicore`. The CUDA and HIP tests check the limit of the
grid (`-5` beyond 2^31 - 1 blocks, without reading the arrays). The
`stlpar-parunseq-gcc` configuration also checks the error reported for an
invalid file of parameters (`-6`, `test-invalid-parameters.cxx`) and the
refusal of a name containing `_generic_parallel`
(`check-refused-name.sh`). All the configurations check the libraries
built with the `generic_functions` option (`check-generic-functions.sh`):
the function of a point is compared with the `generic` interface, and a
parameter changed by `S_setParameter` is used by both kinds of functions.
The property generated with the runtime checks disabled is also built
with this option, including when the configuration file already sets
`disable_runtime_checks`.

The test of the stream is discriminating: a slow kernel writes the arguments on the
stream of the caller just before the evaluation. The values are only correct
if the property waits for this kernel. With the default stream, the test fails.

### `raja` configuration

The `raja` configuration is built alone in `build/tfel-raja` (installed in
`install/tfel-raja`, script `scratch/raja/build-tfel.sh`) and added to
`build/tfel-mmm-gpu` (`scratch/raja/reconfigure-main-tfel.sh`). The script
`scratch/raja/run-tests.sh` runs its tests with the TFEL of
`install/tfel-raja`:

- `behaviours/raja`: comparisons with the `generic` interface on the GPU
  (default policy `RAJA::cuda_exec_async<128>`), with the TDLS sub-test;
- `behaviours/raja-host`: the same comparisons on the host (OpenMP loop of
  RAJA), the libraries being compiled by nvcc;
- `behaviours/raja-options`: valid execution policies copied in the
  generated code, invalid ones refused (injected expressions, unbalanced
  templates, invalid tokens). The script fails if `mfront` is replaced by
  `/bin/true`, which accepts everything;
- `material-properties/raja` and `raja-host`: the shared test of the
  material properties, with an explicit resource of camp, the asynchronous
  evaluation and the preserved stack size on the GPU.

Results on 20/10, all green: `behaviours/raja` 711 comparisons OK
(including the preserved stack size), `behaviours/raja-host` 657
comparisons OK (the preserved stack size is only checked for the libraries
computing on the device), `raja-options` OK, material properties 170
assertions (`raja`) and 164 assertions (`raja-host`), without failure.
Logs in `scratch/raja/` (`run-tests-5.log`), untracked.

### `openacc` configurations

The `openacc-nvhpc-gpu` (`nvc++ -acc=gpu -gpu=mem:unified`) and
`openacc-nvhpc-multicore` (`nvc++ -acc=multicore`) configurations are added
to `build/tfel-mmm-gpu` (`scratch/openacc/reconfigure-main-tfel.sh`). The
script `scratch/openacc/run-tests.sh` runs their tests with the TFEL of
`install/tfel-mmm-gpu`: comparisons of the behaviours with the `generic`
interface (with the TDLS sub-test), and the shared test of the material
properties (with, on the GPU, the stack size preserved through the CUDA
driver, and on the cores of the host, the evaluation without runtime
checks). The reproducers `reproducers/nvhpc-openacc-min-reduction` and
`reproducers/nvhpc-openacc-managed-global` show the two defects of `nvc++`
25.11 worked around by the backend and by the configuration. The GPU
configuration relies on the heterogeneous memory management of the
driver: it is not run with `uvm_disable_hmm=1`.

Results on 20/10, all green: behaviours `openacc-nvhpc-gpu` 657
comparisons OK (with the TDLS sub-test, 1666 s), `openacc-nvhpc-multicore`
all the comparisons OK (1667 s), material properties 171 assertions on the
GPU (including the stack size lowered through the driver, status `-6`) and
all the assertions on the cores of the host (with the evaluation without
runtime checks). Logs in `scratch/openacc/` (`run-tests-*.log`), untracked.

### Sweeps

The test behaviours of MFront (`tfel/mfront/tests`) are much more numerous
than the behaviours of the tests of the interface. Sweeps cover all of them:

| Tool | What is checked | Result |
|---|---|---|
| `scripts/compile-sweep.sh` | compilation of all the test behaviours, models and properties, for nvcc, clang CUDA (which refuses any call to a host function from GPU code), AdaptiveCpp and nvc++ | 541 successful compilations out of 608 files (60 invalid, 7 fragments), no blocking user code, since the evening of 06/10 (see below) |
| `scripts/runtime-sweep.py` | execution of all the test behaviours and models, compared with the `generic` interface, with three loading amplitudes | all the remaining differences are rounding differences, checked |
| `scripts/property-sweep.py` | execution of the 56 test material properties | conforming |
| `scripts/roundoff-check.sh`, `scripts/roundoff-triage.py` | replays a behaviour on the host with and without contraction of the operations (FMA), like nvcc on GPU | the reported differences come from behaviours that amplify the rounding |
| `scripts/clang-cuda-header-scan.sh` | compilation of each header of TFEL by clang CUDA | consistent host and GPU annotations |

### Full validation

`scratch/validate3.sh <label>` runs in sequence:

1. the reconfiguration and the installation of TFEL;
2. the complete test suite of TFEL (`ctest`): 6970 tests as of 12/10,
   all green;
3. the matrix of the `generic-parallel` tests at `-O3` for all the
   configurations (`scripts/run-generic-parallel-tests.sh`);
4. the compilation sweeps with nvcc and clang.

The results of the matrices are in `results/generic-parallel-tests-*.txt`.
The last full validation dates from 08/10
(`validate3.sh 2026-10-08-r8`: 6965 tests of TFEL, `-O3` matrix 31/31,
nvcc and clang sweeps 541 OK; matrix in
`results/generic-parallel-tests-2026-10-08-r8-O3.txt`). After the patches
of 09 and 12/10, the complete suite of TFEL was run again on 12/10
(6970 tests, all green), not the matrix. After the rebase of 13/10 on
`ff0959ece`, only TDLS and the `generic-parallel` tests of `cuda-nvcc`
were run again. After the rebase of 20/10 of all the repositories on the
official masters (TFEL `d24a5cfe8`, MGIS `a30dbf5`, mfem-mgis
`5bbd5a6981`, mfem-mgis-examples `f256bd8`, mm-opera-hpc `f3fc124`): 6980
tests of TFEL, 57 of MGIS, 360 of mfem-mgis, 25 of the examples, 9 of
mm-opera-hpc and 143 of `cases/`, all green; the `-O3` matrix
(`results/generic-parallel-tests-2026-10-20-rebase-O3.txt`), the sweeps
and the suite of mfem-mgis on the GPU are run afterwards. After the
version 3 of the functions of the material properties (21/10,
`validate3.sh 2026-10-20-step2`): 6980 tests of TFEL, `-O3` matrix 41/41
(21 configurations of the behaviours and 20 of the material properties,
`results/generic-parallel-tests-2026-10-20-step2-O3.txt`), nvcc and clang
sweeps 542 OK each, 57 tests of MGIS and 360 of mfem-mgis, all green.
After the rebase on `1dc3ac594`, the common contract of the behaviours and
the `generic_functions` option (21/10, `validate3.sh 2026-10-21-step5`):
6986 tests of TFEL, `-O3` matrix 40/41
(`results/generic-parallel-tests-2026-10-21-step5-O3.txt`), nvcc and clang
sweeps 547 OK each, 58 tests of MGIS, and 358 of the 360 tests of
mfem-mgis in the whole suite. The behaviours failed with
`stlpar-parunseq-nvhpc-gpu` (defect of nvc++, worked around afterwards),
and two tests of mfem-mgis on the GPU failed for lack of memory, and pass
when run alone.

The tracked sweeps in `results/compile-sweep-*-2026-10-06.tsv` give
the state of the morning of 06/10: 504 OK, 36 blocking user codes, 1 call
to a host function. The sweeps at 541 OK (log of the evening of 06/10 and
of 08/10) stayed in `scratch/` (`sweep-r2-*.tsv`, `sweep-r3-*.tsv`),
untracked.

### Targeted tests during development

- `scratch/dev-test.sh <configuration>` runs a test with the development
  tree (MFront from `scratch/build-dev2`, headers from `tfel-dev`),
  without installing TFEL. `DEVTEST_KIND=material-properties` selects the
  material properties.
- `scratch/targeted-test.sh <configuration> <pattern> <behaviour>...` only builds
  the given behaviours.
- `reproducers/`: small standalone programs that check a precise
  point, for example the rotations of the blocks of the tangent operator, the
  character strings in GPU code, the cost of a synchronous call, the
  static variables in CUDA code, or defects of hypre, of AmgX
  and of cuDSS. See `reproducers/README.md`.

## MGIS

MGIS is built in `build/mgis` against the TFEL of the prototype. Its
tests (`ctest`, 58 tests on 21/10) include the original ones and:

| Test | What is checked |
|---|---|
| `GenericParallelLoadTest` | loading, cross refusals of the point and batch functions, memory layout of the copies of the structures of TFEL (size, alignment and offset of each member, initialization of the statistics of the behaviours and of the material properties), report function loaded, and refused for a `generic` library |
| `GenericParallelIntegrateTest` | batches on CPU (`stlpar`) compared with the point by point computation: integration, ranges, partial failures, initialize functions, post-processings, 2D and 3D rotations, finite strain, speed of sound, all the kinds of external variables, model, execution policies, memory locations and resources; asynchronous computation refused for a library without report function |
| `GenericParallelGenericFunctionsTest` | libraries also containing the functions of the `generic` interface (`generic_functions` option): both kinds of functions loaded from the same library, integration of a behaviour and rotations of an orthotropic behaviour point by point compared with the batch, value of a material property at a point compared with the batch |
| `GenericParallelEvaluateTest` | material properties on CPU: loading, values, bounds and policies, locations, invalid arguments, execution policy; refusal of a version of the functions which is not supported, of an invalid memory space and of a library without the symbol of the version (symbols of `GenericParallelInvalidMaterialProperties.cxx`) |
| `GenericParallelCUDATest` | same checks on GPU (`cuda`), data in managed memory, memory resource, non-default stream, model and material properties; uniform value kept on the host refused by `make_view`, given by `make_unchecked_view` |
| `MaterialDataViewTest` | views: `FieldView` and its contracts, `getPoints`; integration, initialize functions and post-processings on views in four layouts (point by point, component by component, with or without a stride larger than needed), identical bit for bit to the manager, with a workspace, with the synchronous and asynchronous policies and through `integrate_debug`, each computation starting from a copy of the state before the computation; the eight integration types, predictions included; speed of sound; model; stored and dissipated energies; empty sets; view without blocks of the tangent; `make_view`; host operations; `update`, `revert`, `updateValues` (with and without context); variables of a state view, blocks of the tangent, extraction; expected refusals, for both interfaces |
| `MaterialDataManagerLayoutTest` | manager stored component by component identical bit for bit to the point by point layout: integration on ranges, empty sets, thread pool (errors included), policies, `revert`; initialization and post-processing in all their forms; end of step of the material properties and external variables (external arrays, local arrays, uniform values, update policies); conversions of finite strain on views and on the manager, in 3D and in plane strain, Piola-Kirchhoff stress compared with its definition; views of the variables, *strided coalesced* views; `MaterialFunctionManager`; `save` and `restore`; `updateInputValues` (inputs updated, arrays of the states not copied) |
| `GenericParallelIntegrateTest`, `GenericParallelEvaluateTest` (additions) | rotations and `evaluate` on views, identical bit for bit to the overloads on arrays, for all the layouts of the source and of the destination, in place or not, matrices given in the five ways, without policy, synchronous and asynchronous |
| `FieldOperationsCUDATest` | CUDA operations (copies, transpositions, fills), complete cycle on GPU (arrays of the caller component by component in GPU memory only, integration and end of step on the GPU, compared with the host), integration, initialization, post-processing, `update`, `revert` and extraction directly on a `MaterialDataView` in GPU memory, rotations and `evaluate` on views in GPU memory, asynchronous computations with a managed memory pool, mass density kept on the host refused |
| `FieldOperationsParallelSTLTest` | the checks of `FieldOperationsCUDATest` (shared body `FieldOperationsTest.hxx`) with the `stlpar` backend on the host, except the two refusals that only concern a device (mixed copy, mass density kept on the host): operations of `ParallelSTLFieldOperations` (parallel algorithms of C++), arrays of the caller in the memory of the host, stresses of the cycle within 1e-12 of the integration point by point (difference measured: 0, same compiler); copies and fills with the policies `par`, `seq` and `unseq` |
| `GenericParallelNVHPCGPUTest`, `FieldOperationsParallelSTLNVHPCGPUTest` | `GenericParallelCUDATest` and the checks of `FieldOperationsCUDATest` for the `stlpar` backend compiled by `nvc++ -stdpar=gpu` (libraries and tests compiled by nvc++, option `enable-generic-parallel-stlpar-nvhpc-gpu-tests`): data in managed memory and in the memory of the device only, end of step by the parallel algorithms on the GPU (policies `par_unseq` and `par`); stresses of the cycle within 1e-11 of the host (difference measured: 1.2e-12; 9.6e-13 with nvcc, checked within 1e-12); in float, Norton within 5e-5 (stresses) and 3e-3 (tangent) of the host, differences of nvc++ on the host as on the GPU (3.1e-5 and 1.5e-3 at the onset of the flow) |
| `AllocationTest` | no heap allocation in the calling thread for the integration, `update` and `revert` of the managers (two layouts) and of the views, the rotations and `evaluate`, for both interfaces |
| `RotationMatrixCopyTest` | copies and moves of the rotation matrices, local and external storage, rotation after destruction of the original |
| `ProfilingTest` (addition) | synchronization function of `setProfilingSynchronization` called at the start and at the end of each active section, never for an inactive section |

The tests of the C and Fortran bindings are built by the CMake target `check`,
not by `build-tests`: without it, `ctest` runs old
executables again. The same holds for the libraries generated by MFront
for the tests: after a change of MFront, `ctest` alone loads the old
libraries (refused by MGIS when the version of the functions changed).

The build `build/mgis-full` adds HDF5 and the C, Fortran and
Python bindings (96 tests), including the tests of the layout in each binding
(`LayoutTest`: layout of the manager, invalid layout refused in C,
blocks of the tangent in Fortran in both layouts, strides of the
numpy arrays and strided arrays refused in Python). The script `scripts/check-mgis-gpu-headers.sh` compiles the
GPU headers of MGIS (copy operations and memory resources) with the
CUDA, HIP, SYCL, Kokkos, OpenMP offload, OpenACC, `nvc++ -stdpar=gpu` and
RAJA toolchains of the prototype, and runs them on the GPU
(`reproducers/mgis-gpu-headers/check.cxx`): seven successes on 20/10 for
all but OpenACC, then OpenACC alone (`nvc++ -acc=gpu`) the same day. RAJA
is skipped if TFEL is built without the `raja` configuration.

The CUDA tests require the CMake option `enable-generic-parallel-cuda-tests`,
the tests of `nvc++ -stdpar=gpu` the option
`enable-generic-parallel-stlpar-nvhpc-gpu-tests` (and `MGIS_NVHPC_CXX_COMPILER`
if `nvc++` is not in the `PATH`). The continuous integration has no `nvc++`:
it does not configure this option, which is only checked on the
development machine. The header `ParallelSTLFieldOperations.hxx` also
compiles with `clang` and `libc++` (`-fexperimental-library`), checked by
hand on 19/10.
The batch libraries of the tests are built by the functions of
`mgis/cmake/modules/generic-parallel.cmake`.

## mfem-mgis

mfem-mgis is built in `build/mfem-mgis` (or in a separate build per
work package, `build/mfem-mgis-<package>`, with the targets `BehaviourTest` and
`GenericParallelBehaviour-*`, outside of `all` and of `build-tests`). Its suite
(293 tests as of 13/10/2026 on `gpu-proto` 898f63edf9, of which 77 have a name
matching the pattern `cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos`; 244 and
58 on 12/10) is run in three ways, which must all pass (the
dated figures below are those of the runs noted in the log):

- as it is: `ctest -j 3 -E
  "cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"` (216 tests, 216/216 on
  13/10 on `gpu-proto` 6b9ebf5bcc in 432 s; 186/186 on 12/10 on
  45afaa0653), then the variants that use the GPU, `ctest -j 2 -R
  "cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"` (77 tests, 77/77 on 13/10
  in 226 s). The variants that compute on the GPU share the lock
  `RESOURCE_LOCK gpu`: without it, hypre reserves a pool on the GPU for
  each process of two simultaneous tests and the allocation fails (seen on
  13/10 with `ExecutionSettingsTest`, fixed by 898f63edf9 then
  a836c34207: the lock is given to any test whose name ends with
  `-cuda`, `-hip` or `-gpu`);
- with all the computations forced to the kernels of the batch path, through the
  variables `MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch`,
  `MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels`,
  `MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels` and
  `MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels` (182/182 on 12/10, GPU
  variants excluded);
- on the GPU, with the variable `MFEM_MGIS_DEVICE=cuda` (equivalent to
  the option `--device cuda` of each test), with or without the
  previous variables, one test at a time (`ctest -j 1 -E "cuda|gpu|GPU|debug"`,
  182 tests, 180/182 on 12/10 in 624 s without the previous variables,
  179/182 in 652 s with them): simultaneous tests
  share the GPU and fail. The script
  `scripts/mfem-mgis-gpu-suite.sh [default|kernels]` (kernels forced without argument) runs this suite
  with the ordered assembly of the stiffness matrices
  (`MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered`), which makes the results
  identical from one run to the next: 184/186 on 13/10 in 668 s with the
  kernels forced (mfem-mgis 54ff8e6c8a), only the two tests that require
  the absence of a device failing. The
  tests with 4 processes also require
  `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1`: otherwise hypre reserves 4 GiB per
  process on the GPU (16 GB). Only `DeviceForallTest-cpu` and
  `DeviceOptionTest-none`, which require the absence of a device,
  then fail. `FixedPointSolverTest-kernels-cpu`, which keeps
  its option `--device cpu`, and `UnitTestParameters`, which does not initialize
  mfem-mgis, then do not go through the GPU. Neither do the variants
  `ExecutionSettingsTest-cpu*` and `ExecutionSettingsTest-refused-*`:
  with the option `--execution`, the variable `MFEM_MGIS_DEVICE` is
  ignored (warning). `FixedPointSolverTest`
  fails from time to time on the GPU with the kernels forced: once
  in 26 runs on 11/10, then in the suite of 12/10 and once in 20
  isolated runs. It is always the acceleration `Crossed2Deltabis`: on
  one step, it falls into a cycle where the relative residual stagnates around 1e-10
  up to the limit of 1,000 iterations, whereas it converges in 71
  iterations in the other runs. The sums done on the GPU
  are not reproducible from one run to the next (even the number
  of Newton iterations varies, 8 or 9); the port of the algorithm of
  MTest is faithful (same histories, same formula, same collinearity
  threshold); on CPU, where the computations are reproducible, the test
  always passes. Sources checked on 12/10 (histories compared over
  several runs, hypre on the GPU or on the host): the assembly of
  the tangent by kernels (atomic additions in the sparse matrix;
  differences from the first step, wherever hypre is),
  then the action of the tangent without assembly; with the tangent
  assembled on the host, the first 478 lines of iterations are
  identical from one run to the next. The second source is in fact the
  same: in the sub-test without assembly, the elastic operator K0 of the
  preconditioner is assembled by kernels with atomic additions
  (the action by kernels and the reduction of MFEM are deterministic).
  The ordered assembly (parameter `TangentOperatorAssembly` set to `Ordered`,
  variable `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`) removes these differences:
  with a device, element matrices stored, then added by
  one thread per value of the matrix, in a fixed order; without a device,
  the host path already adds the elements one after the other (cost
  measured on the beam r4: `results/ordered-assembly-2026-10-12`). With it, four runs of
  `FixedPointSolverTest` on the GPU, kernels forced, give histories
  identical bit for bit (the two variants of the diagnosis), against three
  different histories in three runs with the atomic assembly
  (`scratch/night-2026-10-12/determinisme-ordonne`, untracked). The test
  `FixedPointSolverTest-kernels-ordered-cuda` therefore always passes; the
  ordered cases of `BatchIntegrationTest` check that two identical computations
  give the same results bit for bit.

The execution settings (`--execution cpu|gpu|mixed`, see
[plan-execution.md](plan-execution.md)) are checked by the 33 variants
of `ExecutionSettingsTest`: without option, `--execution cpu` (with the
placement variables, which must be ignored, and with `--device
debug`), `--execution gpu`, `--execution mixed` (default placement,
placement variables ignored, 16 placements of the phases), and 10 refusals of invalid or contradictory options (expected
message). The other variants compute the uniaxial tension of a cube,
whose stress is known, where the settings place each phase,
and check the refusal of the behaviours and of the solvers that compute elsewhere
(`tests/ExecutionSettingsTest.cxx` of mfem-mgis). On 13/10, on
b5e7d12fc4: 230/230 tests without `cuda` or `gpu` in their name, 62/62 with.

The compilation for an AMD GPU (`scripts/build-hip-amd.sh`, MFEM 4.10
compiled with HIP for gfx90a) builds mfem-mgis in
`build/mfem-mgis-hip` and the examples in `build/examples-hip`, without
the tests that require a GPU, then runs the others: 213/213 for
mfem-mgis and 15/15 for the examples on 13/10, build redone from
scratch (mfem-mgis 7db41904c0, examples a48701a). Nothing runs on an
AMD GPU there.

The same script with `HIP_ARCH=gfx1036 GPU_TESTS=ON` also runs the tests
on the AMD GPU integrated in the processor (13 and 14/10, and on 20/10 with
`HYPRE_HIP=ON`). With ROCm 10.1
(packages of TheRock), the writes of a finished kernel to device memory
are lost when this GPU stays idle about 300 µs before the next operation
(`reproducers/hip-idle-write-loss`, reported as
[ROCm/rocm-systems#13469](https://github.com/ROCm/rocm-systems/issues/13469)).
The tests therefore run with the kernels and the copies serialised
(`AMD_SERIALIZE_KERNEL=3`, `AMD_SERIALIZE_COPY=3`):
ExecutionSettingsTest 33/33. The kernels of two finite strain behaviours
in 3D compiled at -O2 or -O3 for this GPU never end, a defect of the code
generation of the compiler of ROCm 10.1
(`reproducers/hip-finite-strain-3d-hang`): the
examples rve and ssna303_3d_hypre stayed blocked at the first launch of a
kernel of their library of behaviours. The script compiles the libraries
of behaviours at -O1 for gfx1036 only. On 20/10, with this option, the
examples pass (26/26, the eleven tests on the GPU included:
ssna303_3d_hypre-gpu in 91 s, rve-gpu in 1.3 s). mfem-mgis passes 269 tests
out of 279: the tests of the displacement gradient (six) and the two
AddBehaviourIntegratorTest need a TFEL providing the displacement
gradient, which the TFEL of this environment does not;
TwoLayerCube-1-PreconditionerOverrideP stops on the host (Newton stagnating
at a relative residual of 2.7e-12); HypreMGRTestP-gpu, which uses no
behaviour, exceeds its 300 s. mox2 failed once on 13/10, then passed four
times out of four, and again on 20/10.

The option `enable-kokkos` (initialization of Kokkos by mfem-mgis) is
validated in two separate builds: `build/mfem-mgis-kokkos-off`
(without Kokkos: `ctest -R KokkosOptionTest`, 6 tests, option refused and
variable ignored with a warning) and `build/mfem-mgis-kokkos` (with
the Kokkos of the environment, built for CUDA, installed in
`install/mfem-mgis-kokkos`: `ctest -R "KokkosOptionTest|-kokkos"`, 7
tests including 4 on the GPU). On the client side, mm-opera-hpc
(`build/mm-opera-hpc-kokkos`) and `cases/` (`build/cases-kokkos`, with
`MM_OPERA_HPC_POLYCRYSTAL_WITH_KOKKOS=ON`) are built against this
installation without looking for Kokkos; their tests `beam-gpu-kokkos` and
`polycrystal-batch-kokkos` (with its comparison) use the behaviour generated
by the kokkos backend.

The heterogeneous memory management (HMM) hides here the host data
read by mistake on the GPU: the mutation that gives the kernel of a model
at the points the host pointers of its inputs is not detected with HMM.

The incomplete factorization ILU(k) of hypre 3.2.0 with level 1 or more
crashes on the GPU (reproducer `reproducers/hypre-ilu-gpu`): mfem-mgis
takes level 0 there by default and refuses the higher levels. The
ILUT factorization (`HypreILUType` 1) works there. The
ParaSails preconditioner, which only computes on the host, is refused when hypre computes
on the GPU. On the host, `mfem::HypreParaSails::SetOperator` reads the
communicator of the previous operator, already destroyed: mfem-mgis gives
the new operator before the call (mfem-mgis 52bfb818b0, reproducer
`reproducers/hypre-parasails-operator`).

The examples of mfem-mgis-examples and the cases of mm-opera-hpc are validated
in `cases/` (build `build/cases-post`, linked to an installation of
mfem-mgis and to the builds of the examples and of mm-opera-hpc): each
case keeps the check of the example (reference or analytical solution)
and runs point by point, by kernels, by batches with each host backend
and on the GPU. On the CPU, `ctest -E
"cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"` (the Kokkos tests
use the GPU); on the GPU, `ctest -R cuda`, the tests that share
an output file or the GPU being serialized by locks.

## Single precision

Each library is also built and tested in `float`, with the
twin spack environments `mmm-float-cpu` and `mmm-double-cpu`, and
`mmm-float-gpu` (MFEM and hypre in single with CUDA). Results: MGIS
49/49, and 51/51 with CUDA; mfem-mgis 223/223 on the CPU (with the
ParaSails tests, hypre using its internal LAPACK, and the tests of the Newton
stagnation criterion) and 280/280 on the GPU (ten tests remain in double
only: UMFPACK, GSLIB, Umpire, cuDSS and CUDA, absent from the MFEM in
float); examples 21/21;
mm-opera-hpc 6/6; cases of `cases/` 134 in float (143 in double). Measurements and
diagnosis: [plan-precision.md](plan-precision.md),
`results/precision-beam-*`, `results/diagnostic-float-gpu-2026-10-14`.

## Compiler defects

The defects encountered are worked around in a targeted way, with a
reproducer for each one: `results/nvhpc-bugs/README.md`. The script
`results/nvhpc-bugs/run-nvcxx-reproducers.sh` replays them, to check
whether they are fixed in a new version of a compiler.
