---
title:  Overview of the `TFEL/Config` library
author: Thomas Helfer
date: 12/06/2026
lang: en-EN
numbersections: true
link-citations: true
colorlinks: true
figPrefixTemplate: "$$i$$"
tblPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
lstPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
header-includes:
bibliography: bibliography.bib
csl: iso690-numeric-en.csl
---

The `TFEL/Config` library is mostly used internally by `TFEL` and
provides some information related to the project.

# Device code

The `TFEL/Config/TFELConfig.hxx` header defines the following macros:

- `TFEL_HOST_DEVICE` declares a function callable on the host and in
  device code (`__host__ __device__` with `CUDA` and `HIP`, nothing
  otherwise).
- `TFEL_DEVICE_CODE` is defined when device code is compiled: device
  compilation pass of `CUDA` and `HIP` sources, device compilation of
  `SYCL` sources. Some compilers translate host and device code in a
  single pass: `nvc++ -stdpar=gpu`, `nvc++ -mp=gpu` and the generic
  compiler of AdaptiveCpp. With those compilers, `TFEL_DEVICE_CODE` is
  defined in the whole translation unit.
- `TFEL_NO_REPORT_CONTRACT_VIOLATION` is defined when contract violations
  can't be reported, in particular in device code.

`TFEL_DEVICE_CODE` allows to select the code which can't be executed on
the device, for instance when throwing an exception, while keeping the
same behaviour on the host.

# `getDefaultSubstitutions`: TFEL's default substitutions

The `getDefaultSubstitutions` function returns some default
substitutions used by `MFront` when reading configuration files and
`tfel-check` when analysing `.check` files.

### General substitutions

- `TFEL_CXX_STANDARD`: the `C++` standard used by the project.

### Supported compilers

- `TFEL_GXX_COMPILER`: path the `g++` compiler. The value associated
  with this key is determined as follows, by order of priority:

  1. by passing the variable `TFEL_GXX_COMPILER` during the invocation of
  `cmake`.
  2. by using `g++` as the main compiler to compile `TFEL`
  3. by searching `g++` using `cmake`'s `find_program` function.
- `TFEL_CLANGXX_COMPILER`: path to `LLVM` `clang++` compiler. The value
  associated with this key is determined as follows, by order of
  priority:

  1. by passing the variable `TFEL_CLANGXX_COMPILER` during the invocation of
  `cmake`,
  2. by using `clang++` as the main compiler to compile `TFEL`,
  3. by searching `clang++` using `cmake`'s `find_program` function.
- `TFEL_ICPX_COMPILER`: path to INTEL's `icpx` compiler. The value
  associated with this key is determined as follows, by order of
  priority:

  1. by passing the variable `TFEL_ICPX_COMPILER` during the invocation of
  `cmake`,
  2. by using `icpx` as the main compiler to compile `TFEL`,
  3. by searching `icpx` using `cmake`'s `find_program` function.
- `TFEL_ACPP_COMPILER`: path to AdaptativeCpp `acpp` compiler. The value
  associated with this key is determined as follows, by order of
  priority:

  1. by passing the variable `TFEL_ACPP_COMPILER` during the invocation of
  `cmake`,
  2. by searching `acpp` using `cmake`'s `find_program` function.
- `TFEL_NVHPC_COMPILER`: path to NVIDIA's `nvhpc` compiler. The value
  associated with this key is determined as follows, by order of
  priority:

  1. by passing the variable `TFEL_NVHPC_COMPILER` during the invocation of
  `cmake`,
  2. by using `nvhpc` as the main compiler to compile `TFEL`,
  3. by searching `nvhpc` using `cmake`'s `find_program` function.
- `TFEL_NVCC_COMPILER`: path to NVIDIA's `nvcc` as found by
  `FindCUDAToolkit`.
- `TFEL_HIPCC_COMPILER`: path to AMD's `hipcc` compiler. The value
  associated with this key is determined as follows, by order of
  priority:

  1. by passing the variable `TFEL_HIPCC_COMPILER` during the invocation
    of `cmake`,
  2. by searching `hipcc` using `cmake`'s `find_program` function.

### Substitutions related to the `HIPCC` compiler

- `TFEL_HIPCC_FLAGS`: additional flags passed to the `hipcc` compiler,
  which depend on the `HIP` platform. On the NVIDIA platform, those flags
  are passed to `nvcc` (`-x cu -arch=native --compiler-options=-fPIC
  --expt-relaxed-constexpr`, diagnostics, and the include directory of
  `HIP` if it is not a directory searched by default). On the AMD
  platform, this substitution is `-fPIC`.
- `TFEL_HIP_INCLUDE_DIR`: paths to `HIP` headers

### Substitutions related to the `CUDA` platform

- `TFEL_CUDA_PATH`: path to NVIDIA's `CUDA` toolkit as found by
  `CMake`'s `FindCUDAToolkit` module.
- `TFEL_CUDA_LIBRARY_DIR`: path to the libraries of the NVIDIA's `CUDA`
  toolkit as found by `CMake`'s `FindCUDAToolkit` module.
- `TFEL_CLANG_CUDA_PATH`: path to the `CUDA` toolkit used by `clang`
  (`cuda-clang` and `hip-clang` configurations), given by the
  `TFEL_CLANG_CUDA_PATH` variable of `cmake`. By default, the toolkit of
  `nvcc` is used. A version of `clang` only supports the versions of
  `CUDA` released before it.
- `TFEL_CLANG_CUDA_LIBRARY_DIR`: path to the libraries of this toolkit.

### Substitutions related to `Kokkos`

Those substitutions are defined if the `kokkos` configuration is listed
in the `generic-parallel-configurations` variable of `cmake`. They are
extracted from the `Kokkos::kokkoscore` target of the `Kokkos`
installation found by `cmake`.

- `TFEL_KOKKOS_CXX_COMPILER`: compiler used to build `Kokkos` (for
  instance `nvcc_wrapper`).
- `TFEL_KOKKOS_CXX_FLAGS`: compilation flags of `Kokkos` (include
  directories, `OpenMP` flags, and `-ccbin` for `nvcc_wrapper`).
- `TFEL_KOKKOS_LINKER_FLAGS`: linker flags of `Kokkos` (`-lkokkoscore`
  and its directory, and `-cudart shared` if `Kokkos` uses `CUDA`).
