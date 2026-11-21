/*!
 * \file   test-suite/material-properties/generic-parallel/runtime.hxx
 * \brief  This file selects the runtime of the tests of the
 * `generic-parallel` interface for material properties (memory holding the
 * data, synchronization, execution policy), shared by the comparison to the
 * `generic` interface (`compare.cxx`) and by the test of the runtime checks
 * (`test-no-runtime-checks.cxx`).
 *
 * The data are allocated in managed memory if this file is compiled by a
 * CUDA compiler or by a HIP compiler. The HIP runtime is used if the HIP
 * language is enabled or if the `MFRONT_GPMP_COMPARE_USE_HIP` macro is
 * defined (required on NVIDIA GPUs). If the `MFRONT_GPMP_COMPARE_USE_SYCL`
 * macro is defined, the data are allocated in shared USM memory (or in host
 * memory if the `MFRONT_GPMP_COMPARE_USE_HOST_MEMORY` macro is defined) and
 * the computations are performed on the device selected by the default SYCL
 * selector. If the `MFRONT_GPMP_COMPARE_USE_KOKKOS` macro is defined,
 * Kokkos is initialized and the data are allocated in its shared memory
 * space (if any). If the `MFRONT_GPMP_COMPARE_USE_ARCANE` macro is
 * defined, the accelerator runtime of Arccore is initialized, the
 * computations are launched on a queue whose execution policy is given by
 * the `MFRONT_GPMP_COMPARE_ARCANE_QUEUE` environment variable (see
 * `getArcaneExecutionPolicy`) and the data are allocated in the unified
 * memory of Arccore. The `MFRONT_GPMP_COMPARE_USE_RAJA` macro makes the
 * execution policies of the `raja` backend available (see
 * `makeExecutionPolicy`). Otherwise, the data are allocated in host memory.
 *
 * \date   20/10/2026
 */

#ifndef LIB_MFRONT_GENERICPARALLELMATERIALPROPERTY_TESTS_RUNTIME_HXX
#define LIB_MFRONT_GENERICPARALLELMATERIALPROPERTY_TESTS_RUNTIME_HXX

#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <algorithm>
#include <type_traits>
#if defined(MFRONT_GPMP_COMPARE_USE_KOKKOS)
#include <Kokkos_Core.hpp>
#include "MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx"
#elif defined(MFRONT_GPMP_COMPARE_USE_SYCL)
#include <sycl/sycl.hpp>
#include "MFront/GenericParallel/SyclBackendExecutionPolicy.hxx"
#elif defined(MFRONT_GPMP_COMPARE_USE_ARCANE)
#include <string_view>
#include <arccore/common/accelerator/Runner.h>
#include <arccore/accelerator/AcceleratorInitializer.h>
#include <arccore/common/MemoryUtils.h>
#include <arccore/common/IMemoryAllocator.h>
#include <arccore/common/AllocatedMemoryInfo.h>
#include <arccore/common/MemoryAllocationArgs.h>
#include "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx"
#elif defined(__HIP__) || defined(MFRONT_GPMP_COMPARE_USE_HIP)
#include <hip/hip_runtime.h>
#include "MFront/GenericParallel/HIPBackendExecutionPolicy.hxx"
#define GPU_API(f) hip##f
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#include "MFront/GenericParallel/CUDABackendExecutionPolicy.hxx"
#define GPU_API(f) cuda##f
#endif
#if defined(MFRONT_GPMP_COMPARE_USE_RAJA)
#include "MFront/GenericParallel/RAJABackendExecutionPolicy.hxx"
#endif
#include "MFront/GenericParallel/BackendExecutionPolicy.h"
#include "MFront/GenericParallel/ParallelSTLBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/OpenMPBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/OpenACCBackendExecutionPolicy.hxx"

[[noreturn]] static void fail(const std::string& msg) {
  std::fprintf(stderr, "%s\n", msg.c_str());
  std::exit(EXIT_FAILURE);
}

#if defined(MFRONT_GPMP_COMPARE_USE_KOKKOS)

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

template <typename T>
static void deallocate(T* const p, const std::size_t) {
  Kokkos::kokkos_free<MemorySpace>(p);
}

static void synchronize() { Kokkos::fence(); }

/*!
 * \return the execution policy of the synchronous computations: the default
 * policy (default instance of the execution space)
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  return nullptr;
}

#elif defined(MFRONT_GPMP_COMPARE_USE_SYCL)

//! \brief nothing to initialize
struct Runtime {};

/*!
 * \return the SYCL queue used by all the computations. The asynchronous
 * computations require an in-order queue. The device is the CPU if
 * `TFEL_SYCL_DEVICE` is `0`, a GPU if it is `1`, and the default device
 * otherwise.
 */
static ::sycl::queue& getQueue() {
#if defined(TFEL_SYCL_DEVICE) && (TFEL_SYCL_DEVICE == 0)
  static auto q = ::sycl::queue{::sycl::cpu_selector_v,
                                ::sycl::property::queue::in_order{}};
#elif defined(TFEL_SYCL_DEVICE) && (TFEL_SYCL_DEVICE == 1)
  static auto q = ::sycl::queue{::sycl::gpu_selector_v,
                                ::sycl::property::queue::in_order{}};
#else
  static auto q = ::sycl::queue{::sycl::property::queue::in_order{}};
#endif
  return q;
}

#if defined(MFRONT_GPMP_COMPARE_USE_HOST_MEMORY)

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      static_cast<T*>(std::malloc(std::max<std::size_t>(n, 1) * sizeof(T)));
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

template <typename T>
static void deallocate(T* const p, const std::size_t) {
  std::free(p);
}

#else /* defined(MFRONT_GPMP_COMPARE_USE_HOST_MEMORY) */

template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p =
      ::sycl::malloc_shared<T>(std::max<std::size_t>(n, 1), getQueue());
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

template <typename T>
static void deallocate(T* const p, const std::size_t) {
  ::sycl::free(p, getQueue());
}

#endif /* defined(MFRONT_GPMP_COMPARE_USE_HOST_MEMORY) */

static void synchronize() { getQueue().wait_and_throw(); }

/*!
 * \return the execution policy of the synchronous computations, on the
 * queue of `getQueue` (the `sycl` backend has no default policy)
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  static const auto p =
      mfront::generic_parallel::makeSyclBackendExecutionPolicy(getQueue());
  return &p;
}

#elif defined(MFRONT_GPMP_COMPARE_USE_ARCANE)

/*!
 * \return the execution policy of the queue of the tests, given by the
 * `MFRONT_GPMP_COMPARE_ARCANE_QUEUE` environment variable: `cuda`, `hip` or
 * `sycl` (queue of the GPU, the default, if Arccore is built with CUDA, HIP
 * or SYCL), `thread` (the default otherwise) or `sequential`
 */
static Arcane::Accelerator::eExecutionPolicy getArcaneExecutionPolicy() {
  using Arcane::Accelerator::eExecutionPolicy;
  const auto* const e = std::getenv("MFRONT_GPMP_COMPARE_ARCANE_QUEUE");
#if defined(ARCCORE_HAS_CUDA)
  const auto q = std::string_view{e == nullptr ? "cuda" : e};
#elif defined(ARCCORE_HAS_HIP)
  const auto q = std::string_view{e == nullptr ? "hip" : e};
#elif defined(ARCCORE_HAS_SYCL)
  const auto q = std::string_view{e == nullptr ? "sycl" : e};
#else
  const auto q = std::string_view{e == nullptr ? "thread" : e};
#endif
  if (q == "cuda") {
    return eExecutionPolicy::CUDA;
  }
  if (q == "hip") {
    return eExecutionPolicy::HIP;
  }
  if (q == "sycl") {
    return eExecutionPolicy::SYCL;
  }
  if (q == "thread") {
    return eExecutionPolicy::Thread;
  }
  if (q == "sequential") {
    return eExecutionPolicy::Sequential;
  }
  fail("invalid queue '" + std::string{q} + "'");
}

/*!
 * \brief initialization of the accelerator runtime of Arccore (4 threads
 * for the `Thread` queues), runner and queue of the tests, and execution
 * policy on this queue. The members are destroyed before the runtime is
 * finalized.
 */
struct Runtime {
  Runtime()
      : runner(getArcaneExecutionPolicy()),
        queue(Arcane::Accelerator::makeQueue(this->runner)),
        policy(this->queue) {
    Runtime::current = this;
  }
  Runtime(Runtime&&) = delete;
  Runtime(const Runtime&) = delete;
  Runtime& operator=(Runtime&&) = delete;
  Runtime& operator=(const Runtime&) = delete;
  ~Runtime() { Runtime::current = nullptr; }
  Arcane::Accelerator::AcceleratorInitializer initializer{true, 4};
  Arcane::Accelerator::Runner runner;
  Arcane::Accelerator::RunQueue queue;
  mfront::generic_parallel::ArcaneBackendExecutionPolicy policy;
  //! \brief runtime of the program
  static inline Runtime* current = nullptr;
};

//! \return the queue of the tests
static Arcane::Accelerator::RunQueue& getQueue() {
  if (Runtime::current == nullptr) {
    fail("the runtime of Arccore is not initialized");
  }
  return Runtime::current->queue;
}

//! \return the allocator of the unified memory of Arccore
static Arcane::IMemoryAllocator* getUnifiedMemoryAllocator() {
  return Arcane::MemoryUtils::getAllocator(
      Arcane::eMemoryResource::UnifiedMemory);
}

//! \return the size in bytes of a block of `n` values (at least one)
template <typename T>
static Arcane::Int64 getAllocationSize(const std::size_t n) {
  return static_cast<Arcane::Int64>(std::max<std::size_t>(n, 1) * sizeof(T));
}

//! \brief data in the unified memory of Arccore, accessible by all the queues
template <typename T>
static T* allocate(const std::size_t n) {
  auto* const p = static_cast<T*>(
      getUnifiedMemoryAllocator()
          ->allocate(Arcane::MemoryAllocationArgs{}, getAllocationSize<T>(n))
          .baseAddress());
  if (p == nullptr) {
    fail("allocation failed");
  }
  return p;
}

//! \brief the allocators of Arccore take the size of the block to free it
template <typename T>
static void deallocate(T* const p, const std::size_t n) {
  getUnifiedMemoryAllocator()->deallocate(
      Arcane::MemoryAllocationArgs{},
      Arcane::AllocatedMemoryInfo(p, getAllocationSize<T>(n)));
}

static void synchronize() { getQueue().barrier(); }

/*!
 * \return the execution policy of the synchronous computations, on the
 * queue of the tests (the `arcane` backend has no default policy)
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  return &(Runtime::current->policy);
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

template <typename T>
static void deallocate(T* const p, const std::size_t) {
  check(GPU_API(Free)(p));
}

static void synchronize() { check(GPU_API(DeviceSynchronize)()); }

/*!
 * \return the execution policy of the synchronous computations: the default
 * policy (default stream)
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  return nullptr;
}

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

template <typename T>
static void deallocate(T* const p, const std::size_t) {
  std::free(p);
}

static void synchronize() {}

/*!
 * \return the execution policy of the synchronous computations: the default
 * policy
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  return nullptr;
}

#endif

/*!
 * \return an execution policy whose identifier is the one expected by the
 * tested library (`<name>_generic_parallel_execution_policy` symbol), with
 * the given common options, on the queue of the tests with SYCL and Arcane
 * and on the default queue otherwise. The program fails if the identifier is
 * not the one of a policy available in this program.
 * \param[in] identifier: identifier of the policy
 * \param[in] o: common options
 */
static std::unique_ptr<mfront::generic_parallel::AbstractBackendExecutionPolicy>
makeExecutionPolicy(
    const std::string& identifier,
    const mfront::generic_parallel::BackendExecutionCommonPolicy& o) {
  namespace mgp = mfront::generic_parallel;
  if (identifier == mgp::ParallelSTLBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::ParallelSTLBackendExecutionPolicy>(o);
  }
  if (identifier == mgp::OpenMPBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::OpenMPBackendExecutionPolicy>(o);
  }
  if (identifier == mgp::OpenACCBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::OpenACCBackendExecutionPolicy>(o);
  }
#if defined(MFRONT_GPMP_COMPARE_USE_KOKKOS)
  using Default = Kokkos::DefaultExecutionSpace;
  using Host = Kokkos::DefaultHostExecutionSpace;
  if (identifier == mgp::KokkosBackendExecutionPolicy<Default>::identifier()) {
    return std::make_unique<mgp::KokkosBackendExecutionPolicy<Default>>(
        mgp::KokkosBackendExecutionOptions<Default>{.common_options = o});
  }
  if (identifier == mgp::KokkosBackendExecutionPolicy<Host>::identifier()) {
    return std::make_unique<mgp::KokkosBackendExecutionPolicy<Host>>(
        mgp::KokkosBackendExecutionOptions<Host>{.common_options = o});
  }
#elif defined(MFRONT_GPMP_COMPARE_USE_SYCL)
  if (identifier == mgp::SyclBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::SyclBackendExecutionPolicy>(getQueue(), o);
  }
#elif defined(MFRONT_GPMP_COMPARE_USE_ARCANE)
  if (identifier == mgp::ArcaneBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::ArcaneBackendExecutionPolicy>(getQueue(), o);
  }
#elif defined(__HIP__) || defined(MFRONT_GPMP_COMPARE_USE_HIP)
  if (identifier == mgp::HIPBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::HIPBackendExecutionPolicy>(
        mgp::HIPBackendExecutionOptions{.common_options = o});
  }
#elif defined(__CUDACC__)
  if (identifier == mgp::CUDABackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::CUDABackendExecutionPolicy>(
        mgp::CUDABackendExecutionOptions{.common_options = o});
  }
#endif
#if defined(MFRONT_GPMP_COMPARE_USE_RAJA)
  using DeviceResource = mgp::raja::Resource<mgp::raja::DefaultExecutionPolicy>;
  using HostResource =
      mgp::raja::Resource<mgp::raja::DefaultHostExecutionPolicy>;
  if (identifier ==
      mgp::RAJABackendExecutionPolicy<DeviceResource>::identifier()) {
    return std::make_unique<mgp::RAJABackendExecutionPolicy<DeviceResource>>(
        mgp::RAJABackendExecutionOptions<DeviceResource>{.common_options = o});
  }
  if (identifier ==
      mgp::RAJABackendExecutionPolicy<HostResource>::identifier()) {
    return std::make_unique<mgp::RAJABackendExecutionPolicy<HostResource>>(
        mgp::RAJABackendExecutionOptions<HostResource>{.common_options = o});
  }
#endif
  fail("no execution policy '" + identifier + "' in this program");
}

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

#endif /* LIB_MFRONT_GENERICPARALLELMATERIALPROPERTY_TESTS_RUNTIME_HXX */
