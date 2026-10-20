---
title: The `generic-parallel` interface for material properties
author: Thomas Helfer
date: 28/04/2026
lang: en-EN
link-citations: true
colorlinks: true
geometry:
  - margin=2cm
papersize: a4
figPrefixTemplate: "$$i$$"
tblPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
---

This page describes the `generic-parallel` interface for material
properties.

The `generic-parallel` interface has been designed to be compatible with
current parallel programming models (see the `Backends` section). The
same backends are available in the [`generic-parallel` interface for
behaviours](generic-parallel-behaviour-interface.html), which describes
the notions shared by both interfaces (memory spaces, queues, stack size
of the threads of the devices, configuration files delivered with
`TFEL`).

# Prototypes of the generated functions

The `generic-parallel` interface for material properties generates two
functions. Those functions can be used:

- directly by including the generated headers and linking with the
  generated library.
- indirectly by loading the generated function using system calls like
  `dlopen`/`dlsym` on `Unix` systems or `LoadLibrary`/`GetProcAddress`
  on `Windows`.

In addition, the `generic-parallel` interface also exports:

- symbols describing the material property (names of the inputs, list of
  parameters, default values, etc.),
- the `<name>_generic_parallel_backend` symbol, which gives the backend of
  the library (`cuda`, `hip`, `stlpar`, `kokkos`, `raja`, `sycl`, `openmp`
  or `openacc`), and the `<name>_generic_parallel_api_version` symbol (an
  `unsigned short`), equal to the `MFRONT_GPMP_API_VERSION` macro of
  `MFront/GenericParallelMaterialProperty/MaterialProperty.h`, which is
  incremented at each incompatible change of the generated functions,
- the `<name>_generic_parallel_memory_space` symbol (a `const char*`),
  which tells where the arguments and the output must be located (`host`,
  `device` or `any`, as for the behaviours). The name of the backend does
  not suffice: for instance, the `stlpar` backend computes on the host or
  on a GPU depending on the compiler. This memory space is set when the
  library is compiled (see `MFront/GenericParallel/MemorySpace.hxx`),
- the `<name>_generic_parallel_report` function, which analyses the
  statistics of an evaluation (see below),
- a function allowing to modify at runtime the value of the parameters
  of the material property, if modification of parameters is allowed
  (see the `parameters_as_static_variables` DSL option).

## First prototype supporting data with strides

The first generated function matches the following prototype:

~~~~{.cxx}
void (*)(mfront_gpmp_OutputStatus* const,    // output status
         mfront_gmp_real* const,             // output values
         const mfront_gmp_size_type,         // output stride
         const mfront_gmp_real* const* const,  // values of the arguments
         const mfront_gmp_size_type* const,  // strides of the arguments
         const mfront_gmp_size_type,         // number of arguments
         const mfront_gmp_size_type,         // number of points
         const mfront_gmp_OutOfBoundsPolicy,  // out of bounds policy
         const mfront_gpmp_ExecutionOptions* const);  // execution options
~~~~

The output status (`mfront_gpmp_OutputStatus`) is described below. The
`mfront_gmp_OutOfBoundsPolicy` enumeration type is described in the page
dedicated to the [`generic` interface for material
properties](generic-material-property-interface.html). The
`mfront_gmp_real` type is by default an alias to double precision
floating point numbers.

If the number of points is null, no computation is made.

A null stride means that the associated argument is uniform. The stride
of the output can be null only if the number of points is one and if all
the strides of the arguments are null. Otherwise, the status `-7` is
reported.

All the backends treat the case where all the strides are equal to one,
and the case of a uniform output (one point, all the strides being null),
as special cases, with an optimized implementation.

## Second prototype

The second generated function matches the following prototype:

~~~~{.cxx}
void (*)(mfront_gpmp_OutputStatus* const,    // output status
         mfront_gmp_real* const,             // output values
         const mfront_gmp_real* const* const,  // values of the arguments
         const mfront_gmp_size_type,         // number of arguments
         const mfront_gmp_size_type,         // number of points
         const mfront_gmp_OutOfBoundsPolicy,  // out of bounds policy
         const mfront_gpmp_ExecutionOptions* const);  // execution options
~~~~

This function assumes that the values of the arguments and the values of
the output are stored contiguously in memory. In other words, this
prototype is equivalent to the the first prototype when all strides are
equal to one.

## Execution options

The last argument of both functions gives the options of the execution
(`mfront_gpmp_ExecutionOptions`, declared in
`MFront/GenericParallelMaterialProperty/MaterialProperty.h`), or is null
for the default options. It has the same meaning on all the backends, as
the execution options of the behaviours:

- `queue`: queue on which the computations are submitted: `cudaStream_t` or
  `hipStream_t` (the default stream if null), pointer to a `sycl::queue`
  (required by the `SYCL` backend), pointer to an instance of the Kokkos
  execution space (default instance if null), pointer to a resource of camp
  matching the execution policy of RAJA (default stream of the current
  device with `CUDA` and `HIP`, default resource otherwise, if null). It is
  not used by the C++ parallel algorithms, `OpenMP` and `OpenACC`. The
  evaluation is
  ordered after the computations previously submitted to the queue, which
  may compute the arguments. With the `SYCL` backend, an out-of-order queue
  gives no order between the kernels: an asynchronous evaluation then
  requires an in-order queue (`sycl::property::queue::in_order`), and is
  refused with an out-of-order queue (status `-5`). A synchronous evaluation
  synchronizes an out-of-order queue before the kernel.
- `statistics`: statistics of an asynchronous evaluation, or null for a
  synchronous evaluation (see below).
- `preserve_stack_size`: if not null, the stack size of the threads of the
  device is never modified and a runtime error (status `-6`) is reported
  if it is too small. Otherwise, as for the behaviours, the `CUDA` and
  `HIP` backends, the `Kokkos` backend on a `CUDA` or `HIP` execution
  space, and the `RAJA` backend with an execution policy of `CUDA` or
  `HIP`, increase the stack size if the kernel requires it. The other
  backends increase it if the `device_stack_size` option is given.

If no statistics are given, the evaluation is synchronous: the function
returns once the values are computed, and the output status reports the
result.

Otherwise, the evaluation may be asynchronous, as the integration of the
behaviours with statistics: the function returns once the computations
are submitted to the queue (`CUDA`, `HIP`, `Kokkos` and `RAJA` backends,
`SYCL` backend if no value is copied), and the output status only reports
invalid arguments and runtime errors. The statistics
(`mfront_gpmp_Statistics`) are initialized by
`mfront_gpmp_initializeStatistics` and must be accessible by the device as
the values. With the C++ parallel algorithms, `OpenMP` and `OpenACC`, which
merge the
statuses of the points by a reduction, the statistics are also updated by
the host. With the automatic data transfers of the `SYCL` backend,
statistics which are not accessible by the device are copied, as the
values. Once the computations are completed, the
`<name>_generic_parallel_report` function reads them on the host and
reports the result as a synchronous evaluation does:

~~~~{.cxx}
// statistics in managed memory
mfront_gpmp_Statistics* statistics = nullptr;
cudaMallocManaged(&statistics, sizeof(mfront_gpmp_Statistics));
mfront_gpmp_initializeStatistics(statistics);
const auto options =
    mfront_gpmp_ExecutionOptions{.queue = stream, .statistics = statistics};
// submission of the computations
Inconel600_YoungModulus(&output, E, 1, args.data(), args_strides.data(), 1,
                        npoints, policy, &options);
// ... other computations submitted to the stream ...
cudaStreamSynchronize(stream);
Inconel600_YoungModulus_generic_parallel_report(&output, statistics, policy);
~~~~

A synchronous evaluation calls the same report function on statistics
held by blocks of memory reused from one call to another (`CUDA`, `HIP`,
`SYCL`, `Kokkos` and `RAJA` backends) or by a local variable (C++ parallel
algorithms, `OpenMP` and `OpenACC`).

The statistics have the same structure as the ones of the behaviours: the
number of points where the evaluation failed, the first of them, the
number of points where a bound is violated with the warning policy, and
the first check which failed (see below). The kernels compute the status
of each point and merge it in the statistics only when an event occurs:
with integer atomic operations (`CUDA`, `HIP`, `SYCL`, `Kokkos` and `RAJA`
backends), or by a reduction (C++ parallel algorithms, `OpenMP` and
`OpenACC`), as
for the behaviours. With the atomic operations, the reporting of the
errors costs nothing when the evaluation succeeds. A reduction combines
the statuses of all the points and copies its result, a few integers, to
the host at each call. If the runtime checks are disabled
(`disable_runtime_checks` DSL option), the statistics are never updated.

## Output status

The output status (`mfront_gpmp_OutputStatus`) summarizes the evaluation
of all the points. Its `status` and `bounds_status` members follow the
conventions of the `generic` interface for material properties:

- `0`: the result has been evaluated at all points;
- `1`: an argument or the result is out of its bounds at some points,
  with the warning policy. The result is evaluated at all points but must
  be used with caution. `bounds_status` gives the rank of the variable
  (starting at one, the output being after the arguments);
- `-1`: an argument or the result is out of its physical bounds, or out
  of its bounds with the strict policy, at some points. The result is not
  evaluated at those points. `bounds_status` gives the opposite of the
  rank of the variable;
- `-4`: the result is not finite (`nan` or infinite) at some points.

The following statuses are specific to this interface:

- `-5`: invalid arguments (number of points, number of arguments,
  negative stride, no queue given to the `SYCL` backend, asynchronous
  evaluation with an out-of-order `SYCL` queue);
- `-6`: runtime error of the backend (allocation, data transfer,
  submission or execution of a kernel, stack size of the threads of the
  device), or invalid parameters. The exceptions thrown by the runtimes
  of the backends based on C++ programming models (`SYCL`, `Kokkos`, `RAJA`,
  C++ parallel algorithms, `OpenMP`, `OpenACC`) are reported with this
  status, as the
  errors of the `CUDA` and `HIP` runtimes. The status `-2` of the
  `generic` interface (exception thrown by the evaluation) is not used:
  the kernels do not throw;
- `-7`: uniform output (null stride) with several points or with non
  uniform arguments.

As for the behaviours, the output status also gives:

- `number_of_failures`: the number of points where the evaluation failed
  (status `-1` or `-4`);
- `first_failure_index`: the index of the first of them, or the maximal
  value of its type. The `generic` interface of the same material
  property, called at this point, gives the detailed error;
- `number_of_unreliable_results`: the number of points where a bound is
  violated with the warning policy.

When several checks fail, the first one in the following order is
reported: invalid result, physical bounds, then bounds, the variables
being taken in order, the output after the arguments. The error message
describes this check, with the number of points concerned.

If the runtime checks are disabled (`disable_runtime_checks` DSL option),
the number of points, the number of arguments, the signs of the strides,
the parameters, the bounds and the results are not checked. The uniform
outputs are still checked (status `-7`).

An empty set of points is accepted without any check, as an empty batch
by the behaviours: nothing is read, written or submitted, and the status
is `0`. The arrays and the queue may then be null.

The value of `errno` is restored before returning. It is not used to
detect the errors: it is thread-local and does not exist on GPUs.

# Backends

Backends are associated with parallel programming models.

## `CUDA` backend

### Example of usage

Let `UO2_ShearModulus.mfront` be an implementation of a material
property computing the shear modulus of uranium dioxide, which depends
on the temperature and on the porosity. This file can be compiled as
follows:

~~~~{.bash}
$ mfront --obuild --configuration-file=config-cuda-XXXX.json \
         --interface=generic-parallel UO2_ShearModulus.mfront
The following library has been built :
- libGenericParallelUO2-cuda.so :  UO2_ShearModulus UO2_ShearModulus2
~~~~

where the configuration file `config-cuda-XXXX.json` provide the
required information to call a `CUDA` compiler (`nvcc` or `clang++`)
with the appropriate flags. This file may also contain the options
associated with the `CUDA` backend, see below. Depending on the parallel
configurations selected during the `cmake` configuration steps, two
configuration files may be installed: `config-cuda-nvcc.json` and
`config-cuda-clang.json` (see the list of the [configuration files
delivered with
`TFEL`](generic-parallel-behaviour-interface.html#configuration-files-delivered-with-tfel)).

`UO2_ShearModulus` implements the first prototype, while
`UO2_ShearModulus2` implements the second.

The following code shows how to call the function `UO2_ShearModulus`. In
this example, the porosity is assumed to be uniform and has thus a
stride equal to zero. For the sake of simplicity, memory allocation on
the host and on the device are handled by the [Thrust
library](https://developer.nvidia.com/thrust).

~~~~{.cxx}
thrust::host_vector<double> G;
thrust::host_vector<double> T = {300, 500, 300, 800};
thrust::host_vector<double> f = {0.1};

thrust::device_vector<double> d_T(T);
thrust::device_vector<double> d_f(f);
thrust::device_vector<double> d_G(T.size());

auto output = mfront_gpmp_OutputStatus{};
const auto policy = GENERIC_MATERIALPROPERTY_NONE_POLICY;
const auto args = std::array<double *, 2u>{thrust::raw_pointer_cast(d_T.data()),
                                     thrust::raw_pointer_cast(d_f.data())};
const auto args_stride = std::array<mfront_gmp_size_type, 2u>{1, 0};
UO2_ShearModulus(&output, thrust::raw_pointer_cast(d_G.data()), 1,
                 args.data(), args_stride.data(), 2, 4, policy, nullptr);
G = d_G;
~~~~

The following remarks can be made:

- the user is responsible for consistently passing pointers to memory
  allocated on the device to the `UO2_ShearModulus` function, including
  uniform values. Shared memory can also be used.
- The [Thrust library](https://developer.nvidia.com/thrust) is not used
  by the `CUDA` backend. It is only used in the example for simplicity.
- The functions generated by the `CUDA` backend and the presented
  example can be compiled with both `nvcc` and `clang++`.

#### Example of configuration file for `nvcc`

The following configuration file exemplifies how to use the `nvcc`
compiler to build the source code generated by the `CUDA` backend

~~~~{.json}
compilation_options : {
  cuda : {
    compiler: "/usr/local/cuda-12.8/bin/nvcc",
    compilation_flags: {"-O2 -std=c++20 -arch=native",
       "-diag-suppress 20012,20208,2810,3189 -diag-error 20011,20014",
       "-Xcompiler -fPIC --expt-relaxed-constexpr"}
  }
}
~~~~

Those flags are the ones of the `cuda-nvcc` configuration. `nvcc` only
emits a warning (numbers 20011 and 20014) when device code calls a host
function, and the resulting library gives wrong results or crashes:
`-diag-error 20011,20014` turns this warning into an error.

#### Example of configuration file for `clang++`

The following configuration file exemplifies how to use the `clang++`
compiler to build the source code generated by the `CUDA` backend:

~~~~{.json}
compilation_options : {
  cuda : {
    compiler: "clang++",
    compilation_flags: {"-O3 -std=c++20 -march=native",
      "-x cuda  --cuda-compile-host-device --cuda-gpu-arch=sm_86",
      "--cuda-path=/usr/local/cuda-12.8/ -fPIC -DPIC"}
  }
},
linking_options: {
  linker_flags: "-L/usr/local/cuda-12.8/lib64 -lcudart -ldl -lrt -pthread"
}
~~~~

The `cuda-clang` configuration delivered with `TFEL` uses `-O2` and
`--cuda-gpu-arch=native`, and the `CUDA`
toolkit given by the `TFEL_CLANG_CUDA_PATH` variable of `cmake` (see the
[`TFEL/Config` page](tfel-config.html)): a version of `clang` only
supports the versions of `CUDA` released before it.

### Options

The only option available is `number_of_threads_per_block`, the number
of threads per block. By default, the value used is \(128\) (\(64\) in
the configuration files delivered with `TFEL`). This number is rounded up
at runtime to a multiple of the number of threads per warp, which is
typically \(32\) or \(64\). Typical values are thus \(64\), \(128\)
and \(256\). This option has the same meaning as for the behaviours.

### Internal details

When required, the strides are passed by value to the kernel, as kernel
parameters, which `CUDA` stores in constant memory.

If the material property exposes parameters, those are stored in a
global variable on the host (the CPU). For efficiency, the values of the
parameters are also passed by value to the kernel when evaluating the
material property. This mechanism allows to keep the flexibility provided by
parameters: this flexibility is required to perform sensitivity
analyses, uncertainty propagation studies or perform a new
identification. Note that parameters' handling can still be disabled by
setting the `parameters_as_static_variables` DSL option to `true`.

The errors are reported in the statistics of the evaluation, held by managed
memory. Each thread computes the status of its point and only accesses the
statistics, with integer atomic operations, if an event occurs, minimizing
the cost of error handling. For a synchronous evaluation, the statistics are
taken from a pool of blocks allocated at the first call on a device and
reused by the following calls (see
`mfront::generic_parallel::StatusBlockPool`), so that no memory is allocated
per call. The `HIP`, `SYCL`, `Kokkos` and `RAJA` backends handle the errors
in the same way, with the memory accessible by the host and the device of
the backend (managed memory, shared USM memory, `Kokkos::SharedSpace`,
managed memory of the resources of camp). Those blocks are never freed,
except with the `Kokkos` backend, which frees them when Kokkos is finalized.
The C++ parallel algorithms, `OpenMP` and `OpenACC` merge the statuses of
the points by a reduction.

Before launching a kernel, the stack size of the threads is increased if
the kernel requires it, as for the behaviours (see the
`preserve_stack_size` execution option).

On output, the user knows the number of points where the evaluation
failed and the first of them. For instance, the `UO2_ShearModulus`
function may report that the temperature is negative at two points, the
first one being the point 12. For material properties, errors are
associated with violation of bounds and physical bounds, and with
invalid results. Those tests are disabled by setting the
`disable_runtime_checks` DSL option to `true`. Note that this option also
disables other checks performed on the host (notably regarding the
number of arguments passed).

## `HIP` backend

The `HIP` backend is similar to the `CUDA` backend.

### Example of usage

~~~~{.cxx}
auto output = mfront_gpmp_OutputStatus{};
auto E = createManagedData(4);
const auto T = createManagedData({300, 500, 300, 800});
const auto policy = mfront_gmp_OutOfBoundsPolicy{};
const auto args = std::array<const double *, 1u>{T.get()};
const auto args_strides = std::array<mfront_gmp_size_type, 1u>{1};
Inconel600_YoungModulus(&output, E.get(), 1, args.data(),
                        args_strides.data(), 1, 4, policy, nullptr);
~~~~

where `createManagedData` is a simple wrapper around `hipMallocManaged`
returning an `std::unique_ptr` for RAII.

## Parallel STL backend (`stlpar`)

### Example of usage

Let `UO2_ShearModulus.mfront` be an implementation of a material
property computing the shear modulus of uranium dioxide, which depends
on the temperature and on the porosity. This file can be compiled as
follows:

~~~~{.bash}
$ mfront --obuild --configuration-file=config-stlpar-parunseq-gcc.json \
         --interface=generic-parallel UO2_ShearModulus.mfront
The following library has been built :
- libGenericParallelUO2-stlpar.so :  UO2_ShearModulus UO2_ShearModulus2
~~~~

where the configuration file `config-stlpar-parunseq-gcc.json` provide the required
information to compile the generated source files and also selects the
execution policy to be used (see below).

`UO2_ShearModulus` implements the first prototype, while
`UO2_ShearModulus2` implements the second.

~~~~{.cxx}
auto G = std::vector<double>(4);
const auto T = std::vector<double>{300, 500, 300, 800};
const auto f = std::vector<double>{0.1};
auto output = mfront_gpmp_OutputStatus{};
const auto policy = GENERIC_MATERIALPROPERTY_NONE_POLICY;
const auto args = std::array<const double *, 2u>{T.data(), f.data()};
const auto args_strides = std::array<mfront_gmp_size_type, 2u>{1, 0};
UO2_ShearModulus(&output, G.data(), 1, args.data(), args_strides.data(), 2, 4,
                 policy, nullptr);
~~~~

#### Example of configuration file for `g++`

The following configuration file exemplifies how to use the `g++`
compiler to build the source code generated by the `stlpar` backend:

~~~~{.json}
interfaces_options: {
  generic-parallel: {
    backend: {stlpar: {execution_policy: "parallel_unsequenced_policy"}}
  }
},
compilation_options : {
  cxx : {
    compiler: "g++",
    compilation_flags: "-O2 -std=c++20 -march=native" 
  }
},
linking_options : {
  linker_flags : "-ltbb -Wl,-z,nodelete"
}
~~~~

The library is marked as never unloaded (`-Wl,-z,nodelete`): unloading
TBB while its threads are alive crashes the process. The
`stlpar-parunseq-gcc` configuration delivered with `TFEL` uses the
optimisation flags of the libraries generated by the `generic` interface
(`-O2 -ftree-vectorize -DNDEBUG -DTFEL_NO_RUNTIME_CHECK_BOUNDS`, see the
[page of the `generic-parallel` interface for
behaviours](generic-parallel-behaviour-interface.html)).

#### Example of configuration file for `nvc++`

The following configuration file exemplifies how to use NVIDIA's `nvc++`
compiler to build the source code generated by the `stlpar` backend and
execute the computation of the material property on the device:

~~~~{.json}
interfaces_options: {
  generic-parallel: {
    backend: {stlpar: {execution_policy: "parallel_unsequenced_policy",
                       device_stack_size: 16384}}
  }
},
compilation_options : {
  cxx : {
    compiler: "nvc++",
    compilation_flags: "-O2 -stdpar=gpu -gpu=mem:managed,nordc -std=c++20 -march=native" 
  }
},
linking_options : {
  linker_flags: "-stdpar=gpu -gpu=mem:managed,nordc"
}
~~~~

Those options are the ones of the `stlpar-parunseq-nvhpc-gpu`
configuration. The memory mode shall be given explicitly
(`-gpu=mem:managed`) and the device code of each source file is compiled
without relocatable device code (`-gpu=nordc`): see the page of the
[`generic-parallel` interface for
behaviours](generic-parallel-behaviour-interface.html) for the reasons.

### Options

Two options are available:

- `execution_policy`, which can have one of the following values:
  - `sequenced_policy` (or `seq`)
  - `unsequenced_policy` (or `unseq`)
  - `parallel_policy` (or `par`)
  - `parallel_unsequenced_policy` (or `par_unseq`, default)
- `device_stack_size`: minimal stack size, in bytes, of the threads of the
  GPU when the computations are offloaded by `nvc++ -stdpar=gpu` (16384 in
  the `stlpar-parunseq-nvhpc-gpu` configuration).

The exact meaning of those policies are implementations defined and may
depend on compiler flags used. See [this
page](https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t)
for details.

### Notes

- `sequenced_policy` denotes a standard sequential loop.
- `unsequenced_policy` seems generally associated to vectorization.
- `parallel_policy` and `parallel_unsequenced_policy` generally rely on
  a multithreading approach. This note applies to the current
  implementations provided by the `Visual Studio`, `gcc`, `clang++` or
  `icpx` compilers and the standard libraries they rely on.
- The compiler flag `-stdpar=gpu` provided by the `nvhpc` compiler has
  two effects:
  - all memory allocations made by standard containers (including
    `std::vector`) are managed: eventual memory transfers between the
    host (the CPU) and the device (the GPU) are done transparently.
  - the parallel algorithms are executed on the device (the GPU).
- `libstdc++` implements the parallel STL algorithms on the top of the
  Threading Building Blocks library and requires an explicit link to
  this library. The configurations linking TBB mark the library as never
  unloaded (`-Wl,-z,nodelete`).

## `SYCL` backend

[`SYCL`](https://khronos.org/sycl) is an open cross-platform abstraction
layer that enables code for heterogeneous and offload processors. `SYCL`
allows mixing code on CPUs (leveraging multithreading), GPUs and
accelerators.

The `SYCL` backend has been tested with the following compilers:
Intel's `icpx` and AdaptiveCpp's `acpp`.

### Example of usage

~~~~{.cxx}
//! \brief type of the memory allocator
using usm_allocator_type =
    sycl::usm_allocator<double, sycl::usm::alloc::host>;
//! \brief a simple alias for a vectors
using usm_vector = std::vector<double, usm_allocator_type>;
//
auto output = mfront_gpmp_OutputStatus{};
//
auto Q = sycl::queue{sycl::gpu_selector_v};
auto usm_allocator = usm_allocator_type{Q};
//
auto E = usm_vector(4 * 3, Q);
const auto T = usm_vector({300, -300, 500, -500, 300, -300, 800, -800}, Q);
//
auto d_E = to_device(E, Q); // function transferring data to the device
auto d_T = to_device(T, Q);
//
const auto policy = mfront_gmp_OutOfBoundsPolicy{};
const auto args = std::array<const double *, 1u>{d_T.get()};
const auto args_strides = std::array<mfront_gmp_size_type, 1u>{2};
const auto options = mfront_gpmp_ExecutionOptions{.queue = &Q};
Inconel600_YoungModulus(&output, d_E.get(), 3, args.data(),
                        args_strides.data(), 1, 4, policy, &options);
~~~~

### Queue

`SYCL` has no default queue: the queue must be given by the execution
options. It gives the device and the context of the unified shared memory
of the calling code. An asynchronous evaluation requires an in-order queue
(`sycl::property::queue::in_order`). With an out-of-order queue, a
synchronous evaluation waits for the computations already submitted to
the queue before the kernel is launched.

### Backend options

This backend has two options:

- `data_transfer_policy` (see below);
- `device_stack_size`: minimal stack size, in bytes, of the threads of an
  NVIDIA GPU used through the `CUDA` runtime of AdaptiveCpp (8192 in the
  `sycl-default-acpp` configuration).

#### Data transfer policy

Two policies regarding data transfer between host and device are
available:

- `automatic` (or `auto`, default): the arrays which are not accessible
  by the device are copied in device memory: the arguments, the output and
  the statistics of an asynchronous evaluation. An array is accessible if
  it is allocated by `SYCL` (unified shared memory of any kind: device,
  shared or host) or if the device is the host. The output is copied too,
  so that the points where the evaluation fails keep their values. The
  output and the statistics are copied back once the kernel is completed:
  the evaluation is then synchronous. This policy is shared with the
  behaviours (`MFront/GenericParallel/SYCLDeviceCopies.hxx`).
- `none`: this policy assumes that the data are accessible by the device.

The errors of the kernels of a synchronous evaluation are delivered to
the asynchronous handler of the queue (`sycl::event::wait_and_throw`), as
for the behaviours: they are reported with the status `-6` if this handler
throws them, as the handlers of most calling codes do.

> **Providing data on the device**
> 
> Data copy in `SYCL` are asynchronous operations. Thus, if the
> data are provided on the device, the user is responsible for ensuring
> that the data are effectively copied before calling `MFront`'s generated
> functions.

## `Kokkos` backend

The `Kokkos` backend evaluates the material property with
`Kokkos::parallel_for` in the execution space given by the
`execution_space` option. The arguments and the output must be accessible
from this execution space, and `Kokkos` must be initialized by the calling
code. The queue of the execution options is a pointer to an instance of
the execution space (the default instance if null). The source files are
compiled with the compiler and the flags used by `Kokkos` (see the
`config-kokkos.json` and `config-kokkos-host.json` configuration files).
As for the behaviours, `Kokkos` must be a shared library, built without
relocatable device code, so that the calling code and the generated
libraries share the same instance of `Kokkos`.

The statistics of the synchronous evaluations are held by memory
allocated in `Kokkos::SharedSpace`, as for the behaviours. This memory
space is provided by all the backends of `Kokkos` except the `OpenACC`
and `OpenMPTarget` backends: a static assertion of the generated sources
reports this requirement. With a `CUDA` or a `HIP` execution space, the
stack size of the threads is increased if the kernel requires it, as for
the behaviours.

### Options

- `execution_space`: `Kokkos::DefaultExecutionSpace` by default,
  `Kokkos::DefaultHostExecutionSpace` in the `kokkos-host` configuration.
  The execution space is copied in the generated code: it must be the name
  of a type, made of identifiers separated by `::`.

## `RAJA` backend

The `RAJA` backend evaluates the material property with `RAJA::forall` and
the execution policy given by the `execution_policy` option. The arguments
and the output must be accessible from this execution policy. The queue of
the execution options is a pointer to a resource of camp matching the
execution policy (`camp::resources::Cuda`, `camp::resources::Hip` or
`camp::resources::Host`). If null, the default stream of the current device
is used with `CUDA` and `HIP`, as for the behaviours, and the default
resource otherwise. The evaluation is asynchronous with an asynchronous
execution policy (such as the default one on GPUs). The source files
are compiled by `nvcc` with extended lambdas if RAJA is built with `CUDA`,
by the `C++` compiler of `TFEL` otherwise, and linked with the shared
library of RAJA (see the `config-raja.json` and `config-raja-host.json`
configuration files and the section of the behaviours).

The statistics of the synchronous evaluations are held by managed memory
allocated by the resource of camp, one block per device, kept until the
end of the process. With an execution policy of `CUDA` or `HIP`, the stack
size of the threads is increased if the kernel requires it, as for the
behaviours. The memory space of the library is `device` for the execution
policies of `CUDA` and `HIP`, `host` otherwise. The execution policies of
OpenMP target and `SYCL` are refused at compile time, as for the
behaviours.

### Options

- `execution_policy`: `mfront::generic_parallel::raja::DefaultExecutionPolicy`
  by default, `mfront::generic_parallel::raja::DefaultHostExecutionPolicy`
  in the `raja-host` configuration file (see
  `MFront/GenericParallel/RAJAExecutionPolicy.hxx`). The execution policy is
  copied in the generated code: it must be the name of a type, made of
  qualified identifiers, integers, the boolean literals and template
  arguments.

## `OpenMP` backend

The `OpenMP` backend evaluates the material property on the default device
(`target teams distribute parallel for` construct) or with the threads of
the host (`parallel for` construct) depending on the `offload` option
(`true` by default). When the computations are offloaded, the arguments and
the output must be accessible from the device, for instance in managed
memory with `nvc++ -mp=gpu -gpu=mem:managed` (see the
`config-openmp-offload-nvhpc-gpu.json` and `config-openmp-host-gcc.json`
configuration files). The `device_stack_size` option has the same meaning
as for the `stlpar` backend (16384 in the `openmp-offload-nvhpc-gpu`
configuration). See the page of the [`generic-parallel` interface for
behaviours](generic-parallel-behaviour-interface.html) for the constraints
of `nvc++ -mp=gpu`.

## `OpenACC` backend

The `OpenACC` backend evaluates the material property with a `parallel
loop` construct, on the device selected by the compiler and the OpenACC
runtime: a GPU with `nvc++ -acc=gpu`, the cores of the host with
`nvc++ -acc=multicore` (see the `config-openacc-nvhpc-gpu.json` and
`config-openacc-nvhpc-multicore.json` configuration files). On a GPU, the
arguments and the output must be accessible from it: the
`openacc-nvhpc-gpu` configuration uses the unified memory mode of `nvc++`
(`-gpu=mem:unified`), see the page of the behaviours. The memory space of
the library is deduced from the compilation options, as for the behaviours.
The evaluation is synchronous: the queue of the execution options is not
used. The statuses of the points are merged by a reduction; the indices are
reduced as `unsigned long`, the minimum reductions of OpenACC on
`unsigned long long` being wrong on NVIDIA GPUs with `nvc++` 25.11. The
`device_stack_size` option sets the stack size of the threads of the GPU
through the CUDA driver (16384 in the `openacc-nvhpc-gpu` configuration).
See the page of the [`generic-parallel` interface for
behaviours](generic-parallel-behaviour-interface.html).

# Tests

Unlike other material property interfaces, the tests for the
`generic-parallel` interface are performed in the installation directory.
This allows to test various configurations that may use various `C++`
compilers, a feature unsupported by the `cmake` build system.

The list of configurations tested must be specified by defining a list
named `generic-parallel-configurations` when invoking `cmake` as follows:

~~~~
cmake -Dgeneric-parallel-configurations="cuda-nvcc;cuda-clang"
~~~~

The same list selects the configuration files installed with `TFEL` and
the tests of the `generic-parallel` interface for behaviours.

The following tests configurations are currently supported:

- `cuda-nvcc`: compilation of `CUDA` files using the `nvcc` compiler.
- `cuda-clang`: compilation of `CUDA` files using the `clang++`
  compiler.
- `hip-hipcc`: compilation of `HIP` files using the `hipcc` compiler.
- `hip-clang`: compilation of `HIP` files using the `clang++` compiler,
  only on the NVIDIA platform (`HIP_PLATFORM=nvidia`). On the AMD
  platform, this configuration installs no configuration file and
  declares no test, but `HIP` and `clang++` are still required.
- `stlpar-parunseq-gcc`, `stlpar-parunseq-clang`, `stlpar-parunseq-icpx`:
  parallel algorithms of the standard library with the `par_unseq`
  execution policy on the threads of the host (`g++`, `clang++` and `icpx`
  with TBB).
- `stlpar-parunseq-nvhpc-gpu`: parallel algorithms offloaded to an NVIDIA
  GPU by `nvc++ -stdpar=gpu`.
- `sycl-default-icpx`: compilation of `SYCL` files using the `icpx` compiler.
- `sycl-default-acpp`: compilation of `SYCL` files using the `acpp`
  compiler on the CPU.
- `sycl-gpu-acpp`: compilation of `SYCL` files using the `acpp`
  compiler on the GPU using the `generic` target.
- `kokkos`: compilation with the compiler and flags of `Kokkos`, tests
  on the default execution space and on the default host execution space.
- `raja`: compilation by `nvcc` (RAJA built with `CUDA`) or the `C++`
  compiler of `TFEL`, tests with the default execution policy and with the
  default host execution policy.
- `openmp-host-gcc`: compilation with `g++ -fopenmp`, computations on the
  threads of the host.
- `openmp-offload-nvhpc-gpu`: compilation with `nvc++ -mp=gpu`,
  computations offloaded to an NVIDIA GPU.
- `openacc-nvhpc-gpu`: compilation with `nvc++ -acc=gpu`, computations
  performed by an NVIDIA GPU.
- `openacc-nvhpc-multicore`: compilation with `nvc++ -acc=multicore`,
  computations performed by the cores of the host.

### Notes

For configurations based on HIP, the variable `HIP_PLATFORM` must be
passed to `cmake` or defined in the environment. The value of this variable
must be either `nvidia` or `amd`. If `HIP_PLATFORM` is `nvidia`, the
`CUDA` platform must be available. The value of `HIP_PLATFORM` is cached,
so that a reconfiguration does not depend on the environment.
