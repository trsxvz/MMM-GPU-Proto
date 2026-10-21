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

# Generated functions

For a material property `S` (the name of the function of the `generic`
interface, `<material>_<law>` or `<law>`), the `generic-parallel`
interface generates two functions, declared by the generated header:

~~~~{.cxx}
// evaluation of the material property on a batch of points
int S_generic_parallel(mfront_gpmp_OutputStatus* const,
                       const mfront_gpmp_BatchData* const);
// analysis of the statistics of an evaluation (see below)
int S_generic_parallel_report(mfront_gpmp_OutputStatus* const,
                              const mfront_gpmp_Statistics* const);
~~~~

Both functions return the status of the output status. The data
structures are declared in
`MFront/GenericParallelMaterialProperty/MaterialProperty.h`, with the
types of the pointers to those functions
(`mfront_gpmp_BatchMaterialPropertyPtr` and `mfront_gpmp_ReportPtr`).
The name of a material property can't contain `_generic_parallel`.

Those functions can be used:

- directly by including the generated headers and linking with the
  generated library.
- indirectly by loading the generated function using system calls like
  `dlopen`/`dlsym` on `Unix` systems or `LoadLibrary`/`GetProcAddress`
  on `Windows`.

In addition, the `generic-parallel` interface also exports:

- the symbols describing the material property exported by the `generic`
  interface (names of the inputs, list of parameters, default values,
  etc.), and the `S_mfront_interface` symbol, equal to `GenericParallel`.
  The library does not contain the function `S` of the `generic`
  interface;
- the `S_generic_parallel_api_version` symbol (an `unsigned short`),
  equal to the `MFRONT_GPMP_API_VERSION` macro (`3`), shifted by 1000 for
  a material property generated with another numerical type than `double`
  (`real_type` DSL option). This symbol tells
  that the library contains the functions of batches of points. It only
  changes with an incompatible change of those functions: a calling code
  must refuse any version it does not know;
- the `S_generic_parallel_backend` symbol, which gives the backend of the
  library (`cuda`, `hip`, `stlpar`, `kokkos`, `raja`, `sycl`, `openmp` or
  `openacc`), for information;
- the `S_generic_parallel_memory_space` symbol (a `const char*`), which
  tells where the values and the statuses of the points must be located
  (`host`, `device` or `any`, as for the behaviours). The statistics must
  also be readable by the host (see below). The
  name of the backend does not suffice: for instance, the `stlpar` backend
  computes on the host or on a GPU depending on the compiler. This memory
  space is set when the library is compiled (see
  `MFront/GenericParallel/MemorySpace.hxx`);
- a function allowing to modify at runtime the value of the parameters
  of the material property, if modification of parameters is allowed
  (see the `parameters_as_static_variables` DSL option).

## Data of a batch

The data of a batch are given by the `mfront_gpmp_BatchData` structure:

~~~~{.cxx}
typedef struct mfront_gpmp_BatchData {
  mfront_gmp_size_type number_of_points;
  mfront_gpmp_Field output;                 // output values
  const mfront_gpmp_ConstField* arguments;  // array read on the host
  mfront_gmp_size_type number_of_arguments;
  mfront_gmp_OutOfBoundsPolicy out_of_bounds_policy;
  int* statuses;                            // optional, code of each point
  mfront_gpmp_ExecutionOptions execution;
} mfront_gpmp_BatchData;
~~~~

The output and each argument are given by a pointer to their first value
(`values`) and by the stride between two points (`point_stride`): the
value at the `i`-th point is `values[i * point_stride]`. A null stride
means that the argument is uniform. The stride of the output can only be
null for a batch of one point. The array of the arguments is read on the
host, and may be null if the material property has no argument. The
`mfront_gmp_OutOfBoundsPolicy` enumeration type is described in the page
dedicated to the [`generic` interface for material
properties](generic-material-property-interface.html). The
`mfront_gmp_real` type is by default an alias to double precision
floating point numbers.

The evaluation at each point is the evaluation of the `generic`
interface: same checks, same status. When several checks fail at a point,
the first check making the evaluation fail is retained, otherwise the
last bound violated with the warning policy, as with the `generic`
interface. At a point where the evaluation fails (status `-1`), the
output is `NaN`. At a point where the result is not finite (status
`-4`), the output is the computed value. If the `statuses` array is
given, the status of each point is written in it: `0`, `1`, `-1` or
`-4`, with the meaning of the status of the `generic` interface. It is
located in the memory space of the library, as the values, and is not
written if an error of the batch is detected before the submission of
the computations.

Contrary to the `generic` interface, `errno` is not used to detect
errors: it is thread-local and does not exist on GPUs. A result which is
not finite (`nan` or infinite) is reported by the status `-4`. The value
of `errno` is restored before returning.

An empty batch is accepted without any check, as an empty batch by the
behaviours: nothing is read, written or submitted, and the status is
`0`. The arrays and the queue may then be null.

All the backends treat a batch of one point, and the case where all the
strides are equal to one, with a kernel without strides.

## Execution options

The `execution` member of the data of the batch gives the options of the
execution (`mfront_gpmp_ExecutionOptions`). It has the same meaning on
all the backends, as the execution options of the behaviours:

- `queue`: queue on which the computations are submitted: `cudaStream_t` or
  `hipStream_t` (the default stream if null), pointer to a `sycl::queue`
  (required by the `SYCL` backend), pointer to an instance of the Kokkos
  execution space (default instance if null), pointer to a resource of camp
  matching the execution policy of RAJA (default stream of the current
  device with `CUDA` and `HIP`, default resource otherwise, if null). It is
  not used by the C++ parallel algorithms, `OpenMP` and `OpenACC`. The
  evaluation is ordered after the computations previously submitted to the
  queue, which may compute the arguments. With the `SYCL` backend, an
  out-of-order queue gives no order between the kernels: an evaluation
  with statistics then requires an in-order queue
  (`sycl::property::queue::in_order`), and is refused with an out-of-order
  queue (status `-5`). An evaluation without statistics synchronizes an
  out-of-order queue before the kernel.
- `statistics`: statistics of the evaluation, or null (see below).
- `preserve_stack_size`: if not null, the stack size of the threads of the
  device is never modified and a runtime error (status `-7`) is reported
  if it is too small. Otherwise, as for the behaviours, the `CUDA` and
  `HIP` backends, the `Kokkos` backend on a `CUDA` or `HIP` execution
  space, and the `RAJA` backend with an execution policy of `CUDA` or
  `HIP`, increase the stack size if the kernel requires it. The other
  backends increase it if the `device_stack_size` option is given.

If no statistics are given, the evaluation is synchronous: the function
returns once the values are computed, and the output status reports the
result.

Otherwise, on all the backends, the output status only reports the errors
of the batch (statuses `-5`, `-6` and `-7`), and the
`S_generic_parallel_report` function gives the result once the
computations are completed. The backends with a queue (`CUDA`, `HIP`,
`SYCL`, `Kokkos` and `RAJA`) may return once the computations are
submitted to the queue: the calling code synchronizes its queue before
calling the report function. The C++ parallel algorithms, `OpenMP` and
`OpenACC` return once the computations are completed. The statistics
(`mfront_gpmp_Statistics`) are initialized by
`mfront_gpmp_initializeStatistics`. They must be accessible by the host
and by the device where the computations are performed, for instance in
managed or shared memory if the memory space of the library is `device`:
the kernels update them, and the report function reads them on the host
(the C++ parallel algorithms, `OpenMP` and `OpenACC` also merge the
statuses of the points in them on the host). With the automatic data
transfers of the `SYCL`
backend, statistics which are not accessible by the device are copied,
as the values, and the evaluation is then synchronous:

~~~~{.cxx}
// statistics in managed memory
mfront_gpmp_Statistics* statistics = nullptr;
cudaMallocManaged(&statistics, sizeof(mfront_gpmp_Statistics));
mfront_gpmp_initializeStatistics(statistics);
auto d = mfront_gpmp_BatchData{/* ... */};
d.execution = mfront_gpmp_ExecutionOptions{.queue = stream,
                                           .statistics = statistics};
// submission of the computations: only the errors of the batch are
// reported here
auto output = mfront_gpmp_OutputStatus{};
Inconel600_YoungModulus_generic_parallel(&output, &d);
// ... other computations submitted to the stream ...
cudaStreamSynchronize(stream);
const auto status =
    Inconel600_YoungModulus_generic_parallel_report(&output, statistics);
~~~~

The report function is only meaningful if the evaluation returned `0`.
A synchronous evaluation calls the same report function on statistics
held by blocks of memory reused from one call to another (`CUDA`, `HIP`,
`SYCL`, `Kokkos` and `RAJA` backends) or by a local variable (C++ parallel
algorithms, `OpenMP` and `OpenACC`): the output status is the same.

The statistics have the same structure as the ones of the behaviours: the
number of points where the evaluation failed, the number of points where
a bound is violated with the warning policy, and two keys giving the
first point of each kind and the check retained at this point. The keys
are opaque: the minimum of the keys gives the first point, whatever the
order of the computations, and only the report function of the same
library decodes them. A structure of statistics is used for one material
property and one batch. The kernels compute the status of each point and
merge it in the statistics only when an event occurs: with integer atomic
operations (`CUDA`, `HIP`, `SYCL`, `Kokkos` and `RAJA` backends), or by a
reduction (C++ parallel algorithms, `OpenMP` and `OpenACC`), as for the
behaviours. With the atomic operations, the reporting of the errors costs
nothing when the evaluation succeeds. A reduction combines the statuses of
all the points and copies its result, a few integers, to the host at each
call.

## Output status

The output status (`mfront_gpmp_OutputStatus`) summarizes the evaluation
of the batch. Its `status` member is the error of the batch, if any,
otherwise the status of the first point where the evaluation failed,
otherwise the status of the first point where the result is unreliable,
otherwise `0`. The statuses of the points follow the conventions of the
`generic` interface for material properties:

- `0`: the result has been evaluated at all points;
- `1`: an argument or the result is out of its bounds at some points,
  with the warning policy. The result is evaluated at all points but must
  be used with caution;
- `-1`: an argument or the result is out of its physical bounds, or out
  of its bounds with the strict policy, at some points. The output is
  `NaN` at those points;
- `-4`: the result is not finite (`nan` or infinite) at some points.

The errors of the batch are specific to this interface. Nothing is then
computed:

- `-5`: invalid arguments: no output status (returned value only), no
  data of the batch, no output values, no
  array of arguments or no values for an argument, invalid number of
  arguments, uniform output with several points, more than
  \(2^{48}-1\) points, more blocks of threads than the grid of the device
  can hold with the `CUDA` and `HIP` backends (\(2^{31}-1\) blocks), no
  queue given to the `SYCL` backend, statistics with an out-of-order
  `SYCL` queue, or no statistics given to the report function;
- `-6`: invalid parameters;
- `-7`: runtime error of the backend: allocation, data transfer,
  submission or execution of a kernel, stack size of the threads of the
  device, or exceptions thrown on the host by the runtimes of the backends
  based on C++ programming models (`SYCL`, `Kokkos`, `RAJA`, C++ parallel
  algorithms, `OpenMP`, `OpenACC`). No exception goes through the
  generated functions: an unexpected exception also gives `-7`.

The output status also gives:

- `number_of_failures`: the number of points where the evaluation failed
  (status `-1` or `-4`);
- `first_failure_index`: the index of the first of them, or the maximal
  value of its type. The `generic` interface of the same material
  property, called at this point, gives the detailed error;
- `number_of_unreliable_results`: the number of points where a bound is
  violated with the warning policy;
- `first_unreliable_result_index`: the index of the first of them, or the
  maximal value of its type;
- `bounds_status`: for the statuses `-1` and `1`, the rank of the
  variable of the check retained at the reported point (starting at one,
  the output being after the arguments), negative if the evaluation failed
  (physical bound, or bound with the strict policy), positive otherwise
  (bound with the warning policy), as with the `generic` interface. It is
  null otherwise;
- `msg`: a message describing the reported point (check, index, numbers of
  points). It never gives the values, which the `generic` interface gives
  at this point. Its text is not part of the interface.

The following differences with the `generic` interface are deliberate:
there is no status `-2` (exception) or `-3` (`errno`) at a point (the
kernels do not throw and `errno` does not exist on GPUs), the bounds
status of a point where the result is not finite is null (the `generic`
interface keeps the bounds status of a bound violated with the warning
policy), the message does not give the values, and no `MTest` file is
written when the evaluation fails.

A material property calling another material property (`@MaterialLaw`)
uses, in device code, the implementation of the called material property
generated in the `mfront::device` namespace (see the [`generic-parallel`
interface for behaviours](generic-parallel-behaviour-interface.html)).
This implementation returns `NaN` if a physical bound of the called
material property is not satisfied, whereas the `generic` interface does
not check those bounds: on a device, the evaluation then fails with the
status `-4`, while the `generic` interface computes the value.

If the runtime checks are disabled (`disable_runtime_checks` DSL option),
as with the `generic` interface, the number of arguments, the parameters,
the bounds and the results are not checked. The checks of the structure
of the batch (pointers, uniform output, number of points, queue of the
`SYCL` backend) remain: they cost little on the host and avoid invalid
memory accesses. The statuses of the points, if requested, are then null,
and the report function returns `0`.

## Library also containing the functions of the `generic` interface

With the `generic_functions` option of the interface, the generated
library also contains the function of the `generic` interface, `S`,
which evaluates the material property at a point:

~~~~{.cxx}
interfaces_options: {
  generic-parallel: {
    generic_functions: true
  }
}
~~~~

The library keeps its name (`GenericParallel<library>-<backend>`) and
exports `Generic` as interface (`S_mfront_interface` symbol): it can be
used as a library of the `generic` interface, by `MTest` or `MGIS` for
instance, and its functions of batches of points are detected by the
`S_generic_parallel_api_version` symbol. The library has one set of
parameters, changed by `S_setParameter`, used by the function of a point
and by the functions of batches of points.

The function of a point is the one of the `generic` interface: same
checks, same statuses (including `-2` for an exception and `-3` for an
error of the `C` library), same error messages, the bounds checks being
disabled by the `NO_GENERIC_BOUNDS_CHECK` macro. It is written in a
separate source file, compiled by the `C++` compiler of the configuration
(the compiler of the backend for the backends based on `C++` programming
models, the `C++` compiler for the `CUDA` and `HIP` backends): its results
are the ones of the `generic` interface up to the rounding differences
between compilers. With the compilers which translate host and device
code in a single pass (`nvc++ -stdpar=gpu`, `nvc++ -mp=gpu`,
`nvc++ -acc=gpu`, generic compiler of AdaptiveCpp), the device code is
selected in the whole translation unit (`TFEL_DEVICE_CODE` macro): a
material property called by the material property (`@MaterialLaw`) is
then the one of device code, which returns `NaN` out of its physical
bounds, as for the functions of batches of points.

The configuration files of `icpx` delivered with `TFEL` request a strict
floating point model (`-fp-model=precise`, the default model of `icpx`
being fast), which also restores the support of `errno` by the
mathematical functions, so that the statuses `-3` are reported as with
the `generic` interface. With `nvc++`, `errno` is only set by the
mathematical functions with the `-Mnobuiltin` option, which the
configuration files delivered with `TFEL` do not give (it may slow down
the computations on the host): the functions of a point then do not
report the statuses `-3`.

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
- libGenericParallelUO2-cuda.so :  UO2_ShearModulus_generic_parallel
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

The following code shows how to call the function
`UO2_ShearModulus_generic_parallel`. In this example, the porosity is
assumed to be uniform and has thus a stride equal to zero. For the sake of simplicity, memory allocation on
the host and on the device are handled by the [Thrust
library](https://developer.nvidia.com/thrust).

~~~~{.cxx}
thrust::host_vector<double> G;
thrust::host_vector<double> T = {300, 500, 300, 800};
thrust::host_vector<double> f = {0.1};

thrust::device_vector<double> d_T(T);
thrust::device_vector<double> d_f(f);
thrust::device_vector<double> d_G(T.size());

const auto args = std::array<mfront_gpmp_ConstField, 2u>{
    mfront_gpmp_ConstField{.values = thrust::raw_pointer_cast(d_T.data()),
                           .point_stride = 1},
    mfront_gpmp_ConstField{.values = thrust::raw_pointer_cast(d_f.data()),
                           .point_stride = 0}};
const auto d = mfront_gpmp_BatchData{
    .number_of_points = 4,
    .output = mfront_gpmp_Field{.values = thrust::raw_pointer_cast(d_G.data()),
                                .point_stride = 1},
    .arguments = args.data(),
    .number_of_arguments = args.size(),
    .out_of_bounds_policy = GENERIC_MATERIALPROPERTY_NONE_POLICY};
auto output = mfront_gpmp_OutputStatus{};
const auto status = UO2_ShearModulus_generic_parallel(&output, &d);
G = d_G;
~~~~

The following remarks can be made:

- the user is responsible for consistently passing pointers to memory
  allocated on the device to the `UO2_ShearModulus_generic_parallel`
  function, including
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
auto E = createManagedData(4);
const auto T = createManagedData({300, 500, 300, 800});
const auto args = std::array<mfront_gpmp_ConstField, 1u>{
    mfront_gpmp_ConstField{.values = T.get(), .point_stride = 1}};
const auto d = mfront_gpmp_BatchData{
    .number_of_points = 4,
    .output = mfront_gpmp_Field{.values = E.get(), .point_stride = 1},
    .arguments = args.data(),
    .number_of_arguments = args.size(),
    .out_of_bounds_policy = GENERIC_MATERIALPROPERTY_NONE_POLICY};
auto output = mfront_gpmp_OutputStatus{};
Inconel600_YoungModulus_generic_parallel(&output, &d);
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
- libGenericParallelUO2-stlpar.so :  UO2_ShearModulus_generic_parallel
~~~~

where the configuration file `config-stlpar-parunseq-gcc.json` provide the required
information to compile the generated source files and also selects the
execution policy to be used (see below).

~~~~{.cxx}
auto G = std::vector<double>(4);
const auto T = std::vector<double>{300, 500, 300, 800};
const auto f = std::vector<double>{0.1};
const auto args = std::array<mfront_gpmp_ConstField, 2u>{
    mfront_gpmp_ConstField{.values = T.data(), .point_stride = 1},
    mfront_gpmp_ConstField{.values = f.data(), .point_stride = 0}};
const auto d = mfront_gpmp_BatchData{
    .number_of_points = 4,
    .output = mfront_gpmp_Field{.values = G.data(), .point_stride = 1},
    .arguments = args.data(),
    .number_of_arguments = args.size(),
    .out_of_bounds_policy = GENERIC_MATERIALPROPERTY_NONE_POLICY};
auto output = mfront_gpmp_OutputStatus{};
UO2_ShearModulus_generic_parallel(&output, &d);
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
// strided argument and output
const auto args = std::array<mfront_gpmp_ConstField, 1u>{
    mfront_gpmp_ConstField{.values = d_T.get(), .point_stride = 2}};
const auto d = mfront_gpmp_BatchData{
    .number_of_points = 4,
    .output = mfront_gpmp_Field{.values = d_E.get(), .point_stride = 3},
    .arguments = args.data(),
    .number_of_arguments = args.size(),
    .out_of_bounds_policy = GENERIC_MATERIALPROPERTY_NONE_POLICY,
    .statuses = nullptr,
    .execution = mfront_gpmp_ExecutionOptions{.queue = &Q}};
Inconel600_YoungModulus_generic_parallel(&output, &d);
~~~~

### Queue

`SYCL` has no default queue: the queue must be given by the execution
options. It gives the device and the context of the unified shared memory
of the calling code. An evaluation with statistics requires an in-order
queue (`sycl::property::queue::in_order`). With an out-of-order queue, an
evaluation without statistics waits for the computations already
submitted to the queue before the kernel is launched.

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
  by the device are copied in device memory: the arguments, the output,
  the statuses of the points and the statistics. An array is accessible if
  it is allocated by `SYCL` (unified shared memory of any kind: device,
  shared or host) or if the device is the host. The output is copied in
  both directions, so that the values between two values of a strided
  output are kept. The output, the statuses and the statistics are copied
  back once the kernel is completed: the evaluation is then synchronous
  if an array is copied. This policy is shared with the
  behaviours (`MFront/GenericParallel/SYCLDeviceCopies.hxx`).
- `none`: this policy assumes that the data are accessible by the device.

The errors of the kernels of a synchronous evaluation are delivered to
the asynchronous handler of the queue (`sycl::event::wait_and_throw`), as
for the behaviours: they are reported with the status `-7` if this handler
throws them, as the handlers of most calling codes do.

> **Providing data on the device**
> 
> Data copy in `SYCL` are asynchronous operations. The generated
> functions submit their commands after the commands previously submitted
> to the given queue: an in-order queue orders them, and the generated
> functions wait for the end of the commands previously submitted to an
> out-of-order queue. If the data are provided by commands submitted to
> another queue, the user is responsible for ensuring that the data are
> effectively copied before calling `MFront`'s generated functions.

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

Each configuration runs:

- the tests of the semantics of the generated functions shared by all the
  backends whose arrays are accessible by the host
  (`test-suite/material-properties/generic-parallel/tests.hxx`): statuses
  of the points, failures and unreliable results in the same batch,
  invalid results, errors of the batch, uniform and strided values,
  report function;
- the tests specific to the backend (queues, asynchronous evaluations,
  stack size of the threads of the device, memory of the device);
- the comparison of the evaluation on batches of points with the
  evaluation of the `generic` interface at each point (`compare.cxx`):
  values, statuses of the points and summary, for the three out of bounds
  policies, contiguous, strided and uniform values, and synchronous and
  asynchronous evaluations. Both libraries are compiled by the tested
  configuration. The values are compared exactly when the computations
  are performed on the host by a compiler which neither contracts nor
  reorders the floating point operations by default (`g++`, `clang++`,
  `nvc++` with OpenACC on the cores of the host, `acpp` on the host), and
  up to a relative tolerance of \(10^{-13}\) otherwise (`icpx` and
  computations on a device). The comparisons include a material property with bounds and
  physical bounds on two arguments and on its output
  (`BoundsOrderTest.mfront`), which checks the retained check;
- the same comparisons for libraries built with the `generic_functions`
  option (`check-generic-functions.sh`): the library must export
  `Generic` as interface, its function of a point must give the statuses,
  the bounds statuses, the error numbers, the messages and the values of
  the `generic` interface, and a parameter changed by its `setParameter`
  function must be used by all its functions.

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
