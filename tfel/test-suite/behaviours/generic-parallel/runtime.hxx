/*!
 * \file   test-suite/behaviours/generic-parallel/runtime.hxx
 * \brief  This file selects the runtime of the comparisons of the
 * `generic-parallel` interface to the `generic` interface (memory holding
 * the data, synchronization, execution policy), shared by the comparisons of
 * the behaviours (`compare.cxx`) and of the material properties
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
 * selector. Otherwise, the data are allocated in host memory. The
 * `MFRONT_GPB_COMPARE_USE_RAJA` macro makes the execution policies of the
 * `raja` backend available (see `makeExecutionPolicy`). If the
 * `MFRONT_GPB_COMPARE_USE_ARCANE` macro is defined, the accelerator runtime
 * of Arccore is initialized, the computations are launched on a queue whose
 * execution policy is given by the `MFRONT_GPB_COMPARE_ARCANE_QUEUE`
 * environment variable (see `getArcaneExecutionPolicy`) and the data are
 * allocated in managed memory if Arccore is built with CUDA.
 *
 * \date   20/10/2026
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_TESTS_RUNTIME_HXX
#define LIB_MFRONT_GENERICPARALLEL_TESTS_RUNTIME_HXX

#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string>
#include <algorithm>
#if defined(MFRONT_GPB_COMPARE_USE_KOKKOS)
#include <Kokkos_Core.hpp>
#include "MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx"
#elif defined(MFRONT_GPB_COMPARE_USE_SYCL)
#include <sycl/sycl.hpp>
#include "MFront/GenericParallel/SyclBackendExecutionPolicy.hxx"
#elif defined(MFRONT_GPB_COMPARE_USE_ARCANE)
#include <string_view>
#include <arccore/common/accelerator/Runner.h>
#include <arccore/accelerator/AcceleratorInitializer.h>
#include <arccore/accelerator/AcceleratorUtils.h>
#include "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx"
#if defined(ARCCORE_HAS_CUDA)
#include <cuda_runtime.h>
#define GPU_API(f) cuda##f
#endif
#elif defined(__HIP__) || defined(MFRONT_GPB_COMPARE_USE_HIP)
#include <hip/hip_runtime.h>
#include "MFront/GenericParallel/HIPBackendExecutionPolicy.hxx"
#define GPU_API(f) hip##f
#elif defined(__CUDACC__)
#include <cuda_runtime.h>
#include "MFront/GenericParallel/CUDABackendExecutionPolicy.hxx"
#define GPU_API(f) cuda##f
#endif
#if defined(MFRONT_GPB_COMPARE_USE_RAJA)
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

/*!
 * \return the execution policy of the synchronous computations: the default
 * policy (default instance of the execution space)
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  return nullptr;
}

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

/*!
 * \return the execution policy of the synchronous computations, on the
 * queue of `getQueue` (the `sycl` backend has no default policy)
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  static const auto p =
      mfront::generic_parallel::makeSyclBackendExecutionPolicy(getQueue());
  return &p;
}

/*!
 * \return an execution policy on an in-order queue on another device than
 * the one of `getQueue`, or a null pointer if there is no other device
 */
static const mfront_gp_AbstractBackendExecutionPolicy*
otherDeviceExecutionPolicy() {
  static auto p = []()
      -> std::optional<mfront::generic_parallel::SyclBackendExecutionPolicy> {
    const auto d = getQueue().get_device();
    for (const auto& o : ::sycl::device::get_devices()) {
      if (o != d) {
        return mfront::generic_parallel::makeSyclBackendExecutionPolicy(
            ::sycl::queue{o, ::sycl::property::queue::in_order{}});
      }
    }
    return std::nullopt;
  }();
  return p.has_value() ? &(*p) : nullptr;
}

#elif defined(MFRONT_GPB_COMPARE_USE_ARCANE)

/*!
 * \return the execution policy of the queue of the tests, given by the
 * `MFRONT_GPB_COMPARE_ARCANE_QUEUE` environment variable: `cuda`, `thread`
 * or `sequential` (`cuda` by default if Arccore is built with CUDA, `thread`
 * otherwise)
 */
static Arcane::Accelerator::eExecutionPolicy getArcaneExecutionPolicy() {
  using Arcane::Accelerator::eExecutionPolicy;
  const auto* const e = std::getenv("MFRONT_GPB_COMPARE_ARCANE_QUEUE");
#if defined(ARCCORE_HAS_CUDA)
  const auto q = std::string_view{e == nullptr ? "cuda" : e};
#else
  const auto q = std::string_view{e == nullptr ? "thread" : e};
#endif
  if (q == "cuda") {
    return eExecutionPolicy::CUDA;
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

#if defined(GPU_API)

static void check(const GPU_API(Error_t) e) {
  if (e != GPU_API(Success)) {
    fail(std::string{"GPU error: "} + GPU_API(GetErrorString)(e));
  }
}

//! \brief data in managed memory, accessible by all the queues
template <typename T>
static T* allocate(const std::size_t n) {
  T* p = nullptr;
  check(GPU_API(MallocManaged)(&p, std::max<std::size_t>(n, 1) * sizeof(T)));
  return p;
}

static void deallocate(void* const p) { check(GPU_API(Free)(p)); }

#else /* defined(GPU_API) */

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

#endif /* defined(GPU_API) */

template <typename T>
static void prefetch(T* const, const std::size_t) {}

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

static void deallocate(void* const p) { std::free(p); }

template <typename T>
static void prefetch(T* const, const std::size_t) {}

static void synchronize() {}

/*!
 * \return the execution policy of the synchronous computations: the default
 * policy
 */
static const mfront_gp_AbstractBackendExecutionPolicy* executionPolicy() {
  return nullptr;
}

#endif

#if defined(MFRONT_GPB_COMPARE_USE_ARCANE)
/*!
 * \return true if no computation is pending on the queue of the tests
 * (checked on the native stream of a CUDA queue, the commands of the
 * other queues are completed when launched)
 */
static bool isQueueIdle() {
#if defined(ARCCORE_HAS_CUDA)
  if (getQueue().executionPolicy() ==
      Arcane::Accelerator::eExecutionPolicy::CUDA) {
    const auto stream =
        Arcane::Accelerator::AcceleratorUtils::toCudaNativeStream(getQueue());
    return cudaStreamQuery(stream) == cudaSuccess;
  }
#endif
  return true;
}

/*!
 * \return true if the computations are performed on a GPU by the CUDA
 * runtime (CUDA queue)
 */
[[maybe_unused]] static bool computesOnGPU() {
  return getQueue().executionPolicy() ==
         Arcane::Accelerator::eExecutionPolicy::CUDA;
}
#elif defined(GPU_API)
/*!
 * \return true if no computation is pending on the default stream, used by
 * the default execution policy of the tests
 */
static bool isQueueIdle() {
  return GPU_API(StreamQuery)(nullptr) == GPU_API(Success);
}

/*!
 * \return true: the libraries computing on the host are distinguished by
 * their memory space
 */
[[maybe_unused]] static bool computesOnGPU() { return true; }
#else
/*!
 * \return true: the runtime can't tell if computations are pending (only
 * checked with CUDA and HIP)
 */
static bool isQueueIdle() { return true; }
#endif

#if !defined(MFRONT_GPB_COMPARE_USE_SYCL)
//! \return a policy on another device (only available with SYCL)
static const mfront_gp_AbstractBackendExecutionPolicy*
otherDeviceExecutionPolicy() {
  return nullptr;
}
#endif /* !defined(MFRONT_GPB_COMPARE_USE_SYCL) */

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
#if defined(MFRONT_GPB_COMPARE_USE_KOKKOS)
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
#elif defined(MFRONT_GPB_COMPARE_USE_SYCL)
  if (identifier == mgp::SyclBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::SyclBackendExecutionPolicy>(getQueue(), o);
  }
#elif defined(MFRONT_GPB_COMPARE_USE_ARCANE)
  if (identifier == mgp::ArcaneBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::ArcaneBackendExecutionPolicy>(getQueue(), o);
  }
#elif defined(__HIP__) || defined(MFRONT_GPB_COMPARE_USE_HIP)
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
#if defined(MFRONT_GPB_COMPARE_USE_RAJA)
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

/*!
 * \return an execution policy of another backend than the one whose
 * identifier is given, refused by the tested library
 * \param[in] identifier: identifier of the policy of the tested library
 */
static std::unique_ptr<mfront::generic_parallel::AbstractBackendExecutionPolicy>
makeOtherBackendExecutionPolicy(const std::string& identifier) {
  namespace mgp = mfront::generic_parallel;
  if (identifier == mgp::ParallelSTLBackendExecutionPolicy::identifier()) {
    return std::make_unique<mgp::OpenMPBackendExecutionPolicy>();
  }
  return std::make_unique<mgp::ParallelSTLBackendExecutionPolicy>();
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

#endif /* LIB_MFRONT_GENERICPARALLEL_TESTS_RUNTIME_HXX */
