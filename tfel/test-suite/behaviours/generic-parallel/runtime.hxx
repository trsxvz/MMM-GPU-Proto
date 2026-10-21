/*!
 * \file   test-suite/behaviours/generic-parallel/runtime.hxx
 * \brief  This file selects the runtime of the comparisons of the
 * `generic-parallel` interface to the `generic` interface (memory holding
 * the data, synchronization, queue), shared by the comparisons of the
 * behaviours (`compare.cxx`) and of the material properties
 * (`test-suite/material-properties/generic-parallel/compare.cxx`).
 *
 * The data are allocated in managed memory if this file is compiled by a
 * CUDA compiler or by a HIP compiler. The HIP runtime is used if the HIP
 * language is enabled or if the `MFRONT_GPB_COMPARE_USE_HIP` macro is
 * defined (required on NVIDIA GPUs). If the `MFRONT_GPB_COMPARE_USE_KOKKOS`
 * macro is defined, Kokkos is initialized and the data are allocated in its
 * shared memory space (if any). If the `MFRONT_GPB_COMPARE_USE_SYCL` macro
 * is defined, the data are allocated in shared USM memory (or in host memory
 * if the `MFRONT_GPB_COMPARE_USE_HOST_MEMORY` macro is defined) and the
 * computations are performed on the device selected by the default SYCL
 * selector. Otherwise, the data are allocated in host memory.
 *
 * \date   20/10/2026
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_TESTS_RUNTIME_HXX
#define LIB_MFRONT_GENERICPARALLEL_TESTS_RUNTIME_HXX

#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <algorithm>
#if defined(MFRONT_GPB_COMPARE_USE_KOKKOS)
#include <Kokkos_Core.hpp>
#elif defined(MFRONT_GPB_COMPARE_USE_SYCL)
#include <sycl/sycl.hpp>
#elif defined(__HIP__) || defined(MFRONT_GPB_COMPARE_USE_HIP)
#include <hip/hip_runtime.h>
#define GPU_API(f) hip##f
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#define GPU_API(f) cuda##f
#endif

[[noreturn]] static void fail(const std::string& msg) {
  std::fprintf(stderr, "%s\n", msg.c_str());
  std::exit(EXIT_FAILURE);
}

#if defined(MFRONT_GPB_COMPARE_USE_KOKKOS)

//! \brief initialization of Kokkos
struct Runtime {
  Kokkos::ScopeGuard guard;
};

//! \brief memory space of the data
using MemorySpace = std::conditional_t<Kokkos::has_shared_space,
                                       Kokkos::SharedSpace,
                                       Kokkos::HostSpace>;

template <typename T>
static T* allocate(const std::size_t n) {
  return static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(
      std::max<std::size_t>(n, 1) * sizeof(T)));
}

static void deallocate(void* const p) { Kokkos::kokkos_free<MemorySpace>(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

static void synchronize() { Kokkos::fence(); }

//! \return the queue passed to the behaviour (default execution space)
static void* queue() { return nullptr; }

#elif defined(MFRONT_GPB_COMPARE_USE_SYCL)

//! \brief nothing to initialize
struct Runtime {};

/*!
 * \return the SYCL queue used by all the computations. The asynchronous
 * computations require an in-order queue.
 */
static ::sycl::queue& getQueue() {
  static auto q = ::sycl::queue{::sycl::property::queue::in_order{}};
  return q;
}

#if defined(MFRONT_GPB_COMPARE_USE_HOST_MEMORY)

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      static_cast<T*>(std::malloc(std::max<std::size_t>(n, 1) * sizeof(T)));
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

static void deallocate(void* const p) { std::free(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

#else /* defined(MFRONT_GPB_COMPARE_USE_HOST_MEMORY) */

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      ::sycl::malloc_shared<T>(std::max<std::size_t>(n, 1), getQueue());
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

static void deallocate(void* const p) { ::sycl::free(p, getQueue()); }

template <typename T>
static void prefetch(T* const p, const std::size_t n) {
  if (n != 0) {
    getQueue().prefetch(p, n * sizeof(T));
  }
}

#endif /* defined(MFRONT_GPB_COMPARE_USE_HOST_MEMORY) */

static void synchronize() { getQueue().wait_and_throw(); }

//! \return the queue passed to the behaviour
static void* queue() { return &getQueue(); }

/*!
 * \return an in-order queue on another device than the one of `getQueue`,
 * or a null pointer if there is no other device
 */
static void* otherQueue() {
  static auto q = []() -> std::optional<::sycl::queue> {
    const auto d = getQueue().get_device();
    for (const auto& o : ::sycl::device::get_devices()) {
      if (o != d) {
        return ::sycl::queue{o, ::sycl::property::queue::in_order{}};
      }
    }
    return std::nullopt;
  }();
  return q.has_value() ? &(*q) : nullptr;
}

#elif defined(GPU_API)

//! \brief nothing to initialize
struct Runtime {};

static void check(const GPU_API(Error_t) e) {
  if (e != GPU_API(Success)) {
    fail(std::string{"GPU error: "} + GPU_API(GetErrorString)(e));
  }
}

template <typename T>
static T* allocate(const std::size_t n) {
  T* p = nullptr;
  check(GPU_API(MallocManaged)(&p, std::max<std::size_t>(n, 1) * sizeof(T)));
  return p;
}

static void deallocate(void* const p) { check(GPU_API(Free)(p)); }

//! \brief move the data to the device before the computation
template <typename T>
static void prefetch(T* const p, const std::size_t n) {
  if (n == 0) {
    return;
  }
  int device = 0;
  check(GPU_API(GetDevice)(&device));
#if defined(__HIPCC__) || defined(__HIP__)
  check(hipMemPrefetchAsync(p, n * sizeof(T), device, 0));
#else
  auto location = cudaMemLocation{};
  location.type = cudaMemLocationTypeDevice;
  location.id = device;
  check(cudaMemPrefetchAsync(p, n * sizeof(T), location, 0, 0));
#endif
}

static void synchronize() { check(GPU_API(DeviceSynchronize)()); }

//! \return the queue passed to the behaviour (default stream)
static void* queue() { return nullptr; }

#else

//! \brief nothing to initialize
struct Runtime {};

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      static_cast<T*>(std::malloc(std::max<std::size_t>(n, 1) * sizeof(T)));
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

static void deallocate(void* const p) { std::free(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

static void synchronize() {}

//! \return the queue passed to the behaviour (not used)
static void* queue() { return nullptr; }

#endif

#if !defined(MFRONT_GPB_COMPARE_USE_SYCL)
//! \return a queue on another device (only available with SYCL)
static void* otherQueue() { return nullptr; }
#endif /* !defined(MFRONT_GPB_COMPARE_USE_SYCL) */

static void* load(const std::string& lib, const std::string& f) {
  auto* const l = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (l == nullptr) {
    fail(dlerror());
  }
  auto* const s = dlsym(l, f.c_str());
  if (s == nullptr) {
    fail(dlerror());
  }
  return s;
}

#endif /* LIB_MFRONT_GENERICPARALLEL_TESTS_RUNTIME_HXX */
