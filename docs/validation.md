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
- `scripts/env.sh` gives the complete environment of all the scripts. It
  exports `OMP_PROC_BIND=spread` and `OMP_PLACES=threads` (Kokkos,
  AdaptiveCpp): the OpenMP runtime of each process of a parallel test then
  binds its main thread to the first processor, and the processes of the
  test, with those of the other parallel tests running at the same time,
  share one core. Before running the parallel tests, unset these variables
  (`unset OMP_PROC_BIND OMP_PLACES`), as `scripts/mfem-mgis-gpu-suite.sh`
  and the measurement scripts of `cases/` do: otherwise the tests are tens
  of times slower or exceed their time limit, on an idle machine.
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
  (`-6`, option `--invalid-parameters`), and line of the error (option
  `--invalid-parameters-message`, invalid number of values on the third
  line);
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
invalid file of parameters (`-6`, `test-invalid-parameters.cxx`, error
on the third line of the file and given in the message), the reading of
a file of parameters which does not change `errno`
(`test-parameters-errno.cxx`, `generic` interface, subnormal value) and
the refusal of a name containing `_generic_parallel`
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
when run alone. With the fixes found by this validation (same day): 6986
tests of TFEL, `-O3` tests of `stlpar-parunseq-gcc` and
`stlpar-parunseq-nvhpc-gpu` 4/4
(`results/generic-parallel-tests-2026-10-21-fixes-O3.txt`), 58 tests of
MGIS and 360 of mfem-mgis, all green.
After the rebase on `1dd066c7e`, the execution policies and the `arcane`
backend (23/10, `validate3.sh 2026-10-23-arcane-rebase`): 6993 tests of
TFEL, `-O3` matrix 43/43
(`results/generic-parallel-tests-2026-10-23-arcane-rebase-O3.txt`), nvcc
and clang sweeps 553 OK each, 58 tests of MGIS, 360 of mfem-mgis, 25 of the
examples, 9 of mm-opera-hpc and 143 cases, all green. The rebase on
`4d50821c9` (fix of the runtime path found by this validation) gives the
same tree; the one on `d489171ee` only adds the changes of master
(documentation, docstrings of the Python bindings, not built by the
prototype). After the rebase on `9b810aa15` (merges of #1060 and #1071,
six patches already upstream): 6993 tests of TFEL, `-O3` tests of the
`arcane`, `cuda-nvcc`, `stlpar-parunseq-gcc` and `openmp-host-gcc`
configurations 8/8
(`results/generic-parallel-tests-2026-10-23-rebase-1060-1071-O3.txt`), 58
tests of MGIS, 360 of mfem-mgis, 25 of the examples, 9 of mm-opera-hpc and
143 cases, all green. GPU suite of mfem-mgis (`scripts/mfem-mgis-gpu-suite.sh`, 265
tests): 262 with the kernels forced, 261 with the default computations,
the two tests requiring the absence of a device failing as expected;
`TwoLayerCube-1-MatrixFreeExecution` (added on 17/10), which checks that
`--execution` gives no warning about the variables of placement, failed
because the suite sets them: it is now excluded from the suite.
`PETScSolversTest` exceeded its time limit once with the default
computations, without any output, and passes three times out of three
alone.
After the rebase of 09/11 on the masters (TFEL `60a42dd12`, MGIS
`b1afeb3`; mfem-mgis, mm-opera-hpc and the examples unchanged): 7831
tests of TFEL, 64 tests of MGIS in double and 60 in float, 372 of mfem-mgis, 25 of
the examples, 9 of mm-opera-hpc and 143 cases, all green; GPU suite of
mfem-mgis with the kernels forced: 262 of 264, the two tests requiring the
absence of a device failing as expected; `-O3` matrix
(`results/generic-parallel-tests-2026-11-09-rebase-O3.txt`): all the
configurations green except `openacc-nvhpc-gpu` and
`stlpar-parunseq-nvhpc-gpu`, which failed while test suites of mfem-mgis
were using the GPU and pass alone
(`results/generic-parallel-tests-2026-11-10-rebase-nvhpc-gpu-rerun-O3.txt`);
see [JOURNAL.md](../JOURNAL.md), 09/11/2026.

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
| `ProfilingTest` (addition) | synchronization function of `setProfilingSynchronization` called at the start and at the end of each active section, never for an inactive section; markers of `setProfilingMarkers` called with the name of each active section, in order, after the synchronization of its start and after the one of its end, never for an inactive section nor for a section begun before they were set, the end being reported by the function set at the start; `ProfilingRange` reported in order with the sections, its name being built only then, and absent from the profiling tree |

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
  give the same results bit for bit. On the A100 GPUs of Jean Zay
  (26/10/2026), the case with the regularization of Faltus and the
  prediction still differed by about one unit in the last place in 5 runs
  out of 6: the prediction computes `B -= K X` with the assembled matrix,
  and MFEM 4.10 makes this product on a GPU with cuSPARSE
  (`CUSPARSE_SPMV_CSR_ALG1`), which does not give the same result at each
  run. The matrices of the ordered assembly use the kernel of MFEM instead
  (one thread per row, `UseGPUSparse(false)`): 15 runs out of 15
  identical (patch 0278 of mfem-mgis).

The layouts of the data at the integration points (parameter
`QuadratureDataLayout`, see [mgis.md](mgis.md#views-and-layouts)) add
passes to this suite. Two environment variables, read when the parameter is
absent, select them:

- `MFEM_MGIS_QUADRATURE_DATA_LAYOUT=ComponentByComponent`, with the kernels
  forced (variables above), or with `MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch`
  only (the finite elements are then computed as each test selects, element
  by element by default);
- `MFEM_MGIS_QUADRATURE_DATA_PADDING=3`: three NaN values follow the values
  of each point in each array. A code that computes the position of a value
  instead of using the strides reads NaNs, and its test fails.

On the GPU, `scripts/mfem-mgis-gpu-suite.sh` runs with the same variables.
Its `default` mode also checks the automatic layout: the materials
integrated by batches on the GPU are then stored component by component.
On 10/11/2026, counting the reruns of the tests fixed during the runs (283
tests without device in their name, 122 of them on several processes):
283/283 by default on the CPU, 282/283 component by component with the
kernels forced, 282/283 component by component by batches, 283/283 with the
padding. The failure of both component by component passes is expected:
`TwoLayerCube-1-MatrixFreeExecution` checks that `--execution` gives no
warning, and the forced variables give one. On the GPU (264 tests):
262/264 component by component with the kernels forced,
262/264 with the padding in the `default` mode, the two tests
that require the absence of a device failing as expected; the 104 tests
whose name selects a device: 104/104 by default, 104/104 component by
component with the kernels forced, 104/104 with the padding. The
runs found defects, all fixed: kernels launched with more than 64 elements
per block (limit of the third dimension of a CUDA block), views of zero or
one point refused in the component by component layout, forces in the
global frame of orthotropic materials stored point by point, and tests
reading arrays without their strides. A review of the changes then found
that the kernels using one thread per point refused the elements of more
than 128 points, which the former kernels accepted: they now process them
by chunks of 128 points (`IntegrationRuleOrderTest` uses a rule of 216
points). The 52 tests touched by the last fixes (mfem-mgis 9d05e2b71e)
pass again by default, component by component with the kernels forced and
with the padding. On the prototype rebuilt with these commits, the tests of
the examples pass on the CPU (25/25), as those of mm-opera-hpc (9/9)
and of the cases (109/109 on the CPU, 37/37 on the GPU, where the
examples and mm-opera-hpc are run).

`BatchIntegrationTest` compares the layouts: 152 comparisons of the
component by component layout and 152 of a layout per array with a padding
of 2, against the point by point layout (gradients, forces and stiffness
matrices, by kernels and element by element, F-bar, regularization of
Faltus et al.). All are identical to the bit. Mutations that ignore one
stride each are all detected by this test: inner forces (475 failed
checks), action of the tangent operator without assembly (10), gradients
(988), F-bar (34), Faltus et al. (36), computations element by element
(288).

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

GPU-aware MPI (option `--gpu-aware-mpi` of `mfem_mgis::initialize`) is
validated by `ctest -R "GPUAwareMPIOptionTest|gpu-aware-mpi"`: 12 tests.
Four run without a GPU (option refused, variable ignored with a warning).
Six check the option and the variable on the GPU; they are skipped if the
MPI library can't read the memory of the GPUs. Two complete cases run on
two processes: TwoLayerCube with HyprePCG, and the parallel solvers of
AmgX, which then report their communicator `MPI_DIRECT` ("Using CUDA-Aware
MPI (GPU Direct) communicator"). The Open MPI 5.0.10 of the environment
reads the memory of the GPUs. With one GPU, the two processes share it:
the exchanges between GPUs (NVLink, GPUDirect RDMA) are not run here, nor
the allocators of Umpire created on another GPU than the GPU 0.

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

The line search of the Newton solver and of NewtonK0 (option `LineSearch`,
patches 0290 and 0291 of 14/11/2026) is validated by `LineSearchTest`
(options, decrease criterion, refused trials and failures of the search on
scalar equations), `AdaptiveLinearSolverToleranceTest` (norm of the linear
model at a reduced step, compared to the one computed on vectors),
`NewtonLineSearchTest` and `FixedPointSolverTest`. The last two use the test
behaviour `SaturatingElasticity`: on a cube loaded close to the saturation
then unloaded in one time step, the Newton iterations are the ones of a
scalar equation; without line search they oscillate between two saturated
states, with it they converge with the expected step lengths. On 14/11/2026:

- on the CI (run 43545707075), the whole suite on the CPU, in double (258
  tests on 258) and in float (254 on 254), and the tests of mm-opera-hpc (12
  on 12 in both precisions);
- on Adastra (separate build on a node of the MI300A partition,
  `cases/polycrystal/adastra-fm/ls-tests.slurm`, jobs 5978719 and 5978819):
  the 11 tests of the line search, and 84 other tests using the non linear
  solvers.

Options disabled, the iterations are unchanged: the polycrystal of 100 grains
on 96 cores gives the same outputs bit for bit as the previous executable
(`ls-identity.slurm`, jobs 5978684 and 5978982: macroscopic values, averages
per grain and lines of the iterations of the logs), with NewtonK0 and each
choice of Eisenstat and Walker, with the Newton solver, and on the raw mesh
whose integration failures move the estimates (139 moved estimates and 21
failed resolutions in 5 s simulated). The suite on a GPU was run afterwards
(next paragraph).

The prediction by linear extrapolation (strategy
`LINEAR_EXTRAPOLATION_PREDICTION`, patch 0292 of 14/11/2026) is validated by
`LinearExtrapolationPredictionTest`, with the Newton solver, the fixed point
on the elastic operator and NewtonK0. A resolution limited to zero iteration
leaves the prediction in the unknowns: on a cube in its elastic domain, moved
by a displacement in t², the test compares the predicted values (previous
increment weighted by the ratio of the time increments, exact imposed
values), the reference of the relative criterion (the one of the default
prediction) and the ratio of the residual of the prediction to this
reference, which is analytical; it checks that the extrapolation is made
again after `revert`, that a second resolution of a time step keeps the
solution of the first one, and that another strategy releases the increment.
With a displacement proportional to the time, no iteration is made from the
second time step; on the plastic cube, the solution is the one of the default
prediction with fewer iterations (3 against 15). Eight mutations of the
library each make the test fail (job 5979297). On 14/11/2026:

- on the CI (run 43554531170): the whole suite on the CPU, in double (262
  tests of mfem-mgis, 14 of mm-opera-hpc) and in float (258 and 14), all
  passed;
- on Adastra (`cases/polycrystal/adastra-fm/pred-tests.slurm`, job 5979268):
  the 15 tests of the prediction and of the line search, and 84 other tests
  using the non linear solvers;
- on one MI300A (`ls-gpu-tests.slurm`, settings of
  `scripts/mfem-mgis-gpu-suite.sh`), the tests using the non linear solvers:
  with the line search only (job 5979191), 95 out of 95 with the kernels
  forced, 95 out of 95 with the computations selected by each test, and the
  test selecting a device; with the extrapolation (job 5979353), 98 out of 98
  in both modes (the same tests and three variants of
  `LinearExtrapolationPredictionTest`), and the test selecting a device.

Without the option, the outputs of the polycrystal of 100 grains remain the
ones of the executable preceding the line search, bit for bit, on the four
cases above (`ls-identity.slurm`, job 5979298): NewtonK0 with each choice of
Eisenstat and Walker, the Newton solver, and the raw mesh (3214 lines of
iterations, 139 estimates moved after a failure of the behaviour, 21 failed
resolutions). The refusal with the non linear solver of PETSc selected
globally is not tested.

Since the patch 0293 of 14/11/2026, a resolution of this strategy starts
again from the default prediction when the behaviours can't be integrated at
the extrapolated estimate (before, it failed and the time step was rejected).
`LinearExtrapolationPredictionTest` checks it with the
`FailingIsotropicLinearElasticity` behaviour, whose integration fails beyond
a trace of the strain of 1e-3, on a cube whose imposed displacement increases
then decreases: at the time step of the reversal, the behaviour is integrated
at the default prediction (trace -4.4e-4) and not at the extrapolated
estimate (-1.72e-3). The first estimate is then the default prediction, bit
for bit, the warning is printed once, and the reference, the iterations and
the solution are the ones of the default strategy; with half of the time
step, the extrapolated estimate is integrated and kept (ratio of its residual
to the reference 11/3, analytical); at the next time step, the extrapolation
of the increment of the time step solved from the default prediction is the
solution. With an imposed stress beyond the limit of the behaviour, the
extrapolated estimate is integrated and the following ones are not: the
resolution fails without starting again. Eighteen mutations of the library
each make the test fail in at least one of its variants
(`pred-mutations.slurm`, job 5979634): the eight ones of the extrapolation
and ten of the fallback (no fallback, default prediction or imposed values
not restored, prediction operators taken for an integration, fallback after
any failure, no warning, failure not reported by each of the three
integration paths, outcome of a later integration kept). On 14/11/2026:

- on the CI (run 43560551512): the whole suite on the CPU, in double (262
  tests of mfem-mgis, 14 of mm-opera-hpc) and in float (258 and 14), all
  passed;
- on Adastra (`pred-tests.slurm`, job 5979633): the 15 tests of the
  prediction and of the line search, and 84 other tests using the non linear
  solvers, all passed;
- on one MI300A (`ls-gpu-tests.slurm`, job 5979655): the tests using the non
  linear solvers, 98 out of 98 with the kernels forced, 98 out of 98 with the
  computations selected by each test, and the test selecting a device;
- without the option, the outputs of the polycrystal of 100 grains are the
  ones of the executable preceding the line search, bit for bit, on the four
  cases (`ls-identity.slurm`, job 5979653);
- with the option, when the fallback is not triggered, the outputs are the
  ones of the executable preceding the fallback, bit for bit (`.res`,
  averages per grain, lines of the iterations): 400 grains on 96 cores, steps
  of 1 s, line search and extrapolation (job 5979654 against job 5979328;
  Slurm step 1308 s against 1307 s); the polycrystal of 5 grains of the tests
  (four computations); and the raw mesh of 100 grains, with the extrapolation
  alone and with the line search (`pred-raw.slurm`, jobs 5979781 and 5979782
  against 5979784 and 5979785). On 4 MI300A (job 5979656, 2 runs), the loop
  takes 103 and 102 s, against 102 and 101 s before the fallback.

The fallback was not triggered in any computation of the polycrystal: it is
only exercised by the unit test. On the raw mesh of 100 grains (5 steps of
1 s, 96 cores), where the behaviour fails to integrate at some estimates of
the resolutions, the line search and the extrapolation accept the 5 steps
without any failure (85 lines of iterations, loop 216 s), against 140 steps,
21 failed resolutions, 1307 lines and 1695 s without them; with the
extrapolation alone, 2 attempts fail (7 steps, 145 lines, 16 estimates moved
after an integration failure, 304 s).

After the review of the patches of 13 and 14/11/2026, tests were added to
three of them, and to the option `--single-material` of the polycrystal of
mm-opera-hpc:

- `IntegrationPointsWeightsTest` (patch 0288): on a cube made of two
  materials whose elements have different sizes, the sum of the weights of
  the integration points kept by the problem is the volume of each material
  (3/8 and 5/8), each weight is the one computed from the transformation of
  its element, a second call returns the same array, the integrals of
  thermodynamic forces which are uniform, then linear in x, are the
  analytical ones, and the weights are not computed again once the mesh is
  modified;
- `BatchIntegrationTest` (patch 0289): material axes given by partial
  quadrature functions are modified in place, in 3D and in plane strain;
  the integration by batches keeps the rotation matrices of the former axes
  until `Material::invalidateRotationMatrices` is called, then gives the
  stresses of the integration point by point (the stresses change by 136 %
  with the axes in 3D and by 61 % in plane strain, the two integrations
  then differ by 0); the kept matrices, transposed ones included, are the
  ones of the current axes;
- `FixedPointSolverTest` (patch 0291): a line search of NewtonK0 which finds
  no acceptable step makes the resolution fail (reductions 1, 2, 4, 5, 6, 8
  and 10 of the steps, then a search failing after 10 reductions); with a
  second, elastic, material integrated by batches even when the first one
  fails, a trial is refused although the norm of its residual would be
  accepted;
- `UniaxialPolycrystalTest-SingleMaterial` and the two tests
  `-SingleMaterial-Comparison` of mm-opera-hpc (patch 0046): one material
  for all the grains gives the main output file and the averages of the
  stresses per grain of one material per grain (differences of 1e-16 on the
  stresses on one and two processes, identical averages; tolerances of 1e-9
  and 2e-5, 1 % in single precision).

Eleven mutations of the library and of the polycrystal each make these
tests fail (`review-mutations.slurm`, jobs 5980689 and 5980798): four of
the kept weights (doubled, of another element, computed again at each call,
ignored by the integrals), three of the kept rotation matrices
(invalidation without effect, transposed matrices not invalidated, matrices
computed at each call), two of the line search of NewtonK0 (resolution
going on after a failed search, integration failure of a trial ignored),
two of the single material of the polycrystal (all the points in the first
grain, averages of the next grain). On 14/11/2026:

- on Adastra (`review-tests.slurm`, jobs 5980637 and 5980787, sources of the
  heads of the local repositories): the 10 tests of these executables and
  of `MaterialArraysTest`, 118 other tests using the materials, the
  rotations and the non linear solvers, and the 11 tests of the polycrystal
  by ctest, all passed;
- on one MI300A (`ls-gpu-tests.slurm`, job 5980897): the tests using the non
  linear solvers, with `IntegrationPointsWeightsTest` and
  `MaterialArraysTest`, 101 out of 101 with the kernels forced, 101 out of
  101 with the computations selected by each test, and the test selecting a
  device;
- on the CI (run 43571405244): the whole suite on the CPU, in double (264
  tests of mfem-mgis, 17 of mm-opera-hpc) and in float (260 and 17), all
  passed;
- without the options, the outputs of the polycrystal of 100 grains are the
  ones of the executable preceding the line search, bit for bit, on the four
  cases (`ls-identity.slurm`, job 5980664);
- with one material for all the grains, the line search and the
  extrapolation, and with one material per grain, the outputs are the ones
  of the executable of the measurements, bit for bit (jobs 5980665 and
  5980666; 100 grains, 96 cores, 3 steps).

With the option `enable-roctx` (patch 0294 of 14/11/2026), the profiling
sections are declared as ranges of roctx to the profilers of ROCm, through
the markers of MGIS (`setProfilingMarkers`, patch 0104). On Adastra, MGIS,
mfem-mgis and the polycrystal are built without the option, with it, and
with the roctx of roctracer instead of the one of rocprofiler-sdk
(`roctx-tests.slurm`, job 5982989, sources of the commits):

- `ProfilingMarkersTest`: no marker on the host; on a MI300A, markers with
  the option and none without it; same profiling tree in every case;
- `ProfilingMarkersTest-rocprofv3-hip`: under `rocprofv3 --marker-trace`,
  the three steps and the nine sections of the test are in the trace of the
  markers, with both libraries of roctx;
- the tests of the polycrystal of mm-opera-hpc (11 out of 11) and nine other
  tests of mfem-mgis on the host (42 out of 42; the executables of 4 other
  tests, matched by the names of the list, were not built), without the
  option;
- without the option, the outputs of the polycrystal of 100 grains (96
  cores, 3 steps, one material, line search and extrapolation) are the ones
  of the executable of the measurements, bit for bit: `.res` of 274 bytes,
  averages per grain, 120 lines of iterations (`ls-identity.slurm`, job
  5983076). The executable built with the option gives them too on the
  cores, where it sets no marker (job 5982469, with the build preceding
  the formatting of the sources of the commits; the same check of the
  final build, job 5983077, was running when the connection to Adastra was
  lost, and its result is not known here);
- with the option and without profiler, on 4 MI300A (400 grains, 20 steps of
  1 s, final configuration; `roctx-ab.slurm`, job 5982471, three executables
  in turn, 3 runs each): loop of 100.96 s on average with the ranges,
  101.74 s without the option and 101.60 s with the executable of the
  measurements, for a dispersion of 0.3 to 0.9 % between the runs of one
  executable. The whole computation opens about 1,100 ranges per process;
- under `rocprofv3 --runtime-trace --marker-trace` (`profile-trace.slurm`,
  jobs 5982472 and 5983078, the 4 processes traced): loop of 105.81 and
  106.31 s, against 101.86 and 101.74 s without profiler on the same
  node (+3.9 and +4.5 %); 2.0 GB of traces, written in 30 s at the end
  of the program. The Perfetto trace of a process holds 1,451,795 slices,
  the sum of the rows of its tables (`pftrace-tracks.py`): none was
  discarded;
- on the CI (run 43585510008): the whole suite on the CPU, in double (61
  tests of MGIS, 265 of mfem-mgis, 17 of mm-opera-hpc) and in float (61,
  261 and 17), all passed, with `ProfilingTest` and
  `ProfilingMarkersTest-cpu`; the builds with CUDA and with HIP, without the
  option.

Two processes out of four are traced by `trace-ranks.sh` on the polycrystal
of the tests, and one under `rocprof-sys-sample`, whose trace also holds
the ranges of roctx, the calls of MPI and the sampled call stacks of the
host. A buffer of 4 GB for the trace of Perfetto
(`--perfetto-buffer-size 4194304`) is refused by Perfetto, and the program
then hangs at its end: the default buffer of 1 GB is kept.

In the trace of 400 grains (`trace-zones.py`, mean of the 4 processes), a
kernel runs during 66.3 % of the loop over the time steps (74.2 % in the
second half of the steps). The tangent operator (`NS::getJacobian`)
and the linear solver (`MFEM::Mult(r,c)`) take 76.6 % of the loop and
67.9 % of the time during which the GPU waits; the integration of
the behaviour (`FPS::integrate`) keeps the GPU busy 92.1 % of its time;
the macroscopic stresses, the post-processings and the end of the time
step, computed on the host, take 9.0 % of the loop.

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
