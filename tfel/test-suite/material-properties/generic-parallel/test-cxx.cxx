/*!
 * \file   test-cxx.cxx
 * \brief
 * \author Thomas Helfer
 * \date 22/04/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of theses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <type_traits>
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
#include <Kokkos_Core.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#endif /* KOKKOS_ENABLE_CUDA */
#include "MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx"
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */
#ifdef MFRONT_GPMP_TEST_USE_RAJA
#include "MFront/GenericParallel/RAJAExecutionPolicy.hxx"
#include "MFront/GenericParallel/RAJABackendExecutionPolicy.hxx"
#endif /* MFRONT_GPMP_TEST_USE_RAJA */
#ifdef MFRONT_GPMP_TEST_USE_ARCANE
#include <cstdlib>
#include <arccore/common/accelerator/Runner.h>
#include <arccore/accelerator/AcceleratorInitializer.h>
#include "MFront/GenericParallel/ArcaneBackendExecutionPolicy.hxx"
#if defined(ARCCORE_HAS_CUDA)
#include <cuda_runtime.h>
#endif /* defined(ARCCORE_HAS_CUDA) */
#endif /* MFRONT_GPMP_TEST_USE_ARCANE */
#include "MFront/GenericParallel/ParallelSTLBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/OpenMPBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/OpenACCBackendExecutionPolicy.hxx"
#ifdef __NVCOMPILER_OPENACC_GPU
#include <dlfcn.h>
#include <openacc.h>
#include "MFront/GenericParallel/OpenACCStackSize.hxx"
#endif /* __NVCOMPILER_OPENACC_GPU */

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

#include "MFront/GenericMaterialProperty/MaterialProperty.h"

// This test is shared by the backends based on the C++ programming model
// (`stlpar`, `openmp`, `openacc`, `kokkos`, `raja`, `arcane`), selected by
// the `MFRONT_GPMP_TEST_BACKEND` macro (`stlpar` by default). The arrays are
// allocated by `std::vector` (on the host, or in managed memory with
// `nvc++ -gpu=mem:managed`), in the shared memory space of Kokkos if the
// `MFRONT_GPMP_TEST_USE_KOKKOS` macro is defined, in the managed memory
// of the resource of camp of the execution policy of RAJA if the
// `MFRONT_GPMP_TEST_USE_RAJA` macro is defined, or in managed memory if the
// `MFRONT_GPMP_TEST_USE_ARCANE` macro is defined and Arccore is built with
// CUDA. With Arcane, the execution policy of the queue of the tests is given
// by the `MFRONT_GPMP_TEST_ARCANE_QUEUE` environment variable (see
// `ArcaneTestRuntime`).
#ifndef MFRONT_GPMP_TEST_BACKEND
#define MFRONT_GPMP_TEST_BACKEND stlpar
#endif /* MFRONT_GPMP_TEST_BACKEND */
#define MFRONT_GPMP_TEST_STRINGIFY2(x) #x
#define MFRONT_GPMP_TEST_STRINGIFY(x) MFRONT_GPMP_TEST_STRINGIFY2(x)
// clang-format off: the spaces would appear in the name of the header
#define MFRONT_GPMP_TEST_HEADER(name) \
  MFRONT_GPMP_TEST_STRINGIFY(name-generic-parallel-MFRONT_GPMP_TEST_BACKEND.hxx)
// clang-format on

#include MFRONT_GPMP_TEST_HEADER(PoissonRatioTest)
#include MFRONT_GPMP_TEST_HEADER(BoundsOrderTest)
#include MFRONT_GPMP_TEST_HEADER(Inconel600_YoungModulus)
#include MFRONT_GPMP_TEST_HEADER(Inconel600_YoungModulus_qt)
#include MFRONT_GPMP_TEST_HEADER(UO2_YoungModulus_Martin1989)
#include MFRONT_GPMP_TEST_HEADER(LinearDataInterpolationTest)
#include MFRONT_GPMP_TEST_HEADER(MaterialLawCallTest)

#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
// execution space of the material properties (`execution_space` option of
// the backend)
#ifndef MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE
#define MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE Kokkos::DefaultExecutionSpace
#endif /* MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE */
//! \brief allocator in the shared memory space of Kokkos
template <typename T>
struct TestAllocator {
  using value_type = T;
  TestAllocator() = default;
  template <typename T2>
  TestAllocator(const TestAllocator<T2>&) noexcept {}
  T* allocate(const std::size_t n) {
    return static_cast<T*>(
        Kokkos::kokkos_malloc<Kokkos::SharedSpace>(n * sizeof(T)));
  }
  void deallocate(T* const p, const std::size_t) noexcept {
    Kokkos::kokkos_free<Kokkos::SharedSpace>(p);
  }
  template <typename T2>
  bool operator==(const TestAllocator<T2>&) const noexcept {
    return true;
  }
};
#elif defined(MFRONT_GPMP_TEST_USE_RAJA)
// execution policy of the material properties (`execution_policy` option of
// the backend)
#ifndef MFRONT_GPMP_TEST_RAJA_EXECUTION_POLICY
#define MFRONT_GPMP_TEST_RAJA_EXECUTION_POLICY \
  mfront::generic_parallel::raja::DefaultExecutionPolicy
#endif /* MFRONT_GPMP_TEST_RAJA_EXECUTION_POLICY */
//! \return the resource of camp used by the material properties with the
//! default policy (default stream of the current device with CUDA)
static auto getTestResource() {
  return mfront::generic_parallel::raja::getDefaultResource<
      MFRONT_GPMP_TEST_RAJA_EXECUTION_POLICY>();
}
/*!
 * \return a resource given explicitly by the execution policy: a resource
 * on a stream created here with CUDA (not the default stream, see
 * `test15`), the default resource otherwise
 */
static auto makeExplicitTestResource() {
  using Resource = decltype(getTestResource());
#if defined(RAJA_CUDA_ACTIVE)
  if constexpr (std::is_same_v<Resource, camp::resources::Cuda>) {
    auto stream = cudaStream_t{};
    if (cudaStreamCreate(&stream) != cudaSuccess) {
      throw std::runtime_error("cudaStreamCreate failed");
    }
    return camp::resources::Cuda::CudaFromStream(stream);
  }
#endif
  return getTestResource();
}
//! \brief allocator in the managed memory of the resource of camp
template <typename T>
struct TestAllocator {
  using value_type = T;
  TestAllocator() = default;
  template <typename T2>
  TestAllocator(const TestAllocator<T2>&) noexcept {}
  T* allocate(const std::size_t n) {
    auto r = getTestResource();
    return r.template allocate<T>(n, camp::resources::MemoryAccess::Managed);
  }
  void deallocate(T* const p, const std::size_t) noexcept {
    auto r = getTestResource();
    r.deallocate(p, camp::resources::MemoryAccess::Managed);
  }
  template <typename T2>
  bool operator==(const TestAllocator<T2>&) const noexcept {
    return true;
  }
};
#elif defined(MFRONT_GPMP_TEST_USE_ARCANE)
/*!
 * \brief accelerator runtime of Arccore (4 threads for the `Thread` queues),
 * runner and queue of the tests, and execution policy on this queue. The
 * execution policy of the queue is given by the
 * `MFRONT_GPMP_TEST_ARCANE_QUEUE` environment variable: `cuda`, `thread` or
 * `sequential` (`cuda` by default if Arccore is built with CUDA, `thread`
 * otherwise). The members are destroyed before the runtime is finalized.
 */
struct ArcaneTestRuntime {
  //! \return the execution policy of the queue of the tests
  static Arcane::Accelerator::eExecutionPolicy getExecutionPolicy() {
    using Arcane::Accelerator::eExecutionPolicy;
    const auto* const e = std::getenv("MFRONT_GPMP_TEST_ARCANE_QUEUE");
#if defined(ARCCORE_HAS_CUDA)
    const auto q = std::string_view{e == nullptr ? "cuda" : e};
#else
    const auto q = std::string_view{e == nullptr ? "thread" : e};
#endif
    if (q == "cuda") {
      return eExecutionPolicy::CUDA;
    }
    if (q == "sequential") {
      return eExecutionPolicy::Sequential;
    }
    if (q != "thread") {
      throw std::runtime_error("invalid queue '" + std::string{q} + "'");
    }
    return eExecutionPolicy::Thread;
  }
  ArcaneTestRuntime()
      : runner(getExecutionPolicy()),
        queue(Arcane::Accelerator::makeQueue(this->runner)),
        policy(this->queue) {
    ArcaneTestRuntime::current = this;
  }
  ArcaneTestRuntime(ArcaneTestRuntime&&) = delete;
  ArcaneTestRuntime(const ArcaneTestRuntime&) = delete;
  ArcaneTestRuntime& operator=(ArcaneTestRuntime&&) = delete;
  ArcaneTestRuntime& operator=(const ArcaneTestRuntime&) = delete;
  ~ArcaneTestRuntime() { ArcaneTestRuntime::current = nullptr; }
  Arcane::Accelerator::AcceleratorInitializer initializer{true, 4};
  Arcane::Accelerator::Runner runner;
  Arcane::Accelerator::RunQueue queue;
  mfront::generic_parallel::ArcaneBackendExecutionPolicy policy;
  //! \brief runtime of the tests, created by `main`
  static inline ArcaneTestRuntime* current = nullptr;
};
#if defined(ARCCORE_HAS_CUDA)
//! \brief allocator in managed memory, accessible by all the queues
template <typename T>
struct TestAllocator {
  using value_type = T;
  TestAllocator() = default;
  template <typename T2>
  TestAllocator(const TestAllocator<T2>&) noexcept {}
  T* allocate(const std::size_t n) {
    T* p = nullptr;
    if (cudaMallocManaged(&p, n * sizeof(T)) != cudaSuccess) {
      throw std::bad_alloc();
    }
    return p;
  }
  void deallocate(T* const p, const std::size_t) noexcept { cudaFree(p); }
  template <typename T2>
  bool operator==(const TestAllocator<T2>&) const noexcept {
    return true;
  }
};
#else  /* defined(ARCCORE_HAS_CUDA) */
template <typename T>
using TestAllocator = std::allocator<T>;
#endif /* defined(ARCCORE_HAS_CUDA) */
#else  /* MFRONT_GPMP_TEST_USE_KOKKOS */
template <typename T>
using TestAllocator = std::allocator<T>;
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */

//! \brief array of values given to the material properties
template <typename T>
using test_vector = std::vector<T, TestAllocator<T>>;

#ifdef __NVCOMPILER_OPENACC_GPU
/*!
 * \brief get or set the stack size of the threads of the GPU of OpenACC
 * through the CUDA driver, on the context of OpenACC (the OpenACC runtime of
 * `nvc++` does not load the CUDA runtime)
 * \return true on success
 * \param[in,out] s: stack size, in bytes
 * \param[in] set: if true, the stack size is set, otherwise it is read
 */
static bool accessOpenACCStackSize(std::size_t& s, const bool set) {
  // `CU_LIMIT_STACK_SIZE` in `CUlimit` and `CUDA_SUCCESS` in `CUresult`
  constexpr int stack_size_limit = 0;
  constexpr int success = 0;
  using PushCurrent = int (*)(void*);
  using PopCurrent = int (*)(void**);
  using GetLimit = int (*)(std::size_t*, int);
  using SetLimit = int (*)(int, std::size_t);
  auto* const context = acc_get_current_cuda_context();
  auto* const h = ::dlopen("libcuda.so.1", RTLD_LAZY | RTLD_NOLOAD);
  if ((context == nullptr) || (h == nullptr)) {
    if (h != nullptr) {
      ::dlclose(h);
    }
    return false;
  }
  const auto push =
      reinterpret_cast<PushCurrent>(::dlsym(h, "cuCtxPushCurrent_v2"));
  const auto pop =
      reinterpret_cast<PopCurrent>(::dlsym(h, "cuCtxPopCurrent_v2"));
  const auto get_limit =
      reinterpret_cast<GetLimit>(::dlsym(h, "cuCtxGetLimit"));
  const auto set_limit =
      reinterpret_cast<SetLimit>(::dlsym(h, "cuCtxSetLimit"));
  auto r = (push != nullptr) && (pop != nullptr) && (get_limit != nullptr) &&
           (set_limit != nullptr) && (push(context) == success);
  if (r) {
    r = set ? (set_limit(stack_size_limit, s) == success)
            : (get_limit(&s, stack_size_limit) == success);
    void* previous = nullptr;
    r = (pop(&previous) == success) && r;
  }
  ::dlclose(h);
  return r;
}
#endif /* __NVCOMPILER_OPENACC_GPU */

//! \brief name of the backend
#define MFRONT_GPMP_TEST_BACKEND_NAME \
  MFRONT_GPMP_TEST_STRINGIFY(MFRONT_GPMP_TEST_BACKEND)

#include "batch.hxx"

/*!
 * \return the execution of a synchronous evaluation: default policy, or
 * policy on the queue of the tests with Arcane (the `arcane` backend has no
 * default policy)
 */
static TestExecution testExecution() {
#if defined(MFRONT_GPMP_TEST_USE_ARCANE)
  return TestExecution{.policy = &(ArcaneTestRuntime::current->policy)};
#else
  return TestExecution{};
#endif
}  // end of testExecution

/*!
 * \return an execution policy of the tested backend: instance of the
 * execution space with Kokkos, resource of the tests with RAJA, queue of the
 * tests with Arcane, common options only otherwise
 * \param[in] o: common options
 */
static std::unique_ptr<mfront::generic_parallel::AbstractBackendExecutionPolicy>
makeTestPolicy(
    const mfront::generic_parallel::BackendExecutionCommonPolicy& o) {
  using namespace mfront::generic_parallel;
#if defined(MFRONT_GPMP_TEST_USE_KOKKOS)
  using ExecutionSpace = MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE;
  return std::make_unique<KokkosBackendExecutionPolicy<ExecutionSpace>>(
      KokkosBackendExecutionOptions<ExecutionSpace>{.common_options = o});
#elif defined(MFRONT_GPMP_TEST_USE_RAJA)
  using Resource = decltype(getTestResource());
  return std::make_unique<RAJABackendExecutionPolicy<Resource>>(
      RAJABackendExecutionOptions<Resource>{.common_options = o,
                                            .resource = getTestResource()});
#elif defined(MFRONT_GPMP_TEST_USE_ARCANE)
  return std::make_unique<ArcaneBackendExecutionPolicy>(
      ArcaneTestRuntime::current->queue, o);
#else
  const auto backend = std::string_view{MFRONT_GPMP_TEST_BACKEND_NAME};
  if (backend == "openmp") {
    return std::make_unique<OpenMPBackendExecutionPolicy>(o);
  }
  if (backend == "openacc") {
    return std::make_unique<OpenACCBackendExecutionPolicy>(o);
  }
  return std::make_unique<ParallelSTLBackendExecutionPolicy>(o);
#endif
}  // end of makeTestPolicy

#include "tests.hxx"

//! \brief tests specific to the backends based on the C++ programming model
struct BackendTest : public tfel::tests::TestCase {
  BackendTest()
      : tfel::tests::TestCase("MFront/MaterialProperty/GenericParallel",
                              "BackendTest") {}  // end of BackendTest

  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
#if defined(MFRONT_GPMP_TEST_USE_ARCANE)
    this->test4();
#endif
    return this->result;
  }  // end of execute

 private:
  [[nodiscard]] static double fct(const double TK) noexcept {
    const double TC = TK - 273.15;
    return (-3.1636e-3 * TC * TC - 3.8654 * TC + 2.1421e+4) * 1e7;
  }  // end of fct

  void test1() {
    // execution policy given explicitly: instance of the execution space
    // with Kokkos, resource of camp with RAJA (on a stream created by the
    // test with CUDA), queue created by the test with Arcane, common options
    // only otherwise
    constexpr auto eps = double{1e-2};
#if defined(MFRONT_GPMP_TEST_USE_KOKKOS)
    using ExecutionSpace = MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE;
    const auto execution_policy =
        mfront::generic_parallel::makeKokkosBackendExecutionPolicy(
            mfront::generic_parallel::KokkosBackendExecutionOptions<
                ExecutionSpace>{.execution_space = ExecutionSpace {}});
#elif defined(MFRONT_GPMP_TEST_USE_RAJA)
    auto resource = makeExplicitTestResource();
    const auto execution_policy =
        mfront::generic_parallel::makeRAJABackendExecutionPolicy(
            mfront::generic_parallel::RAJABackendExecutionOptions<
                decltype(resource)>{.resource = resource});
#elif defined(MFRONT_GPMP_TEST_USE_ARCANE)
    const auto execution_policy =
        mfront::generic_parallel::makeArcaneBackendExecutionPolicy(
            Arcane::Accelerator::makeQueue(ArcaneTestRuntime::current->runner));
#else
    const auto execution_policy_holder = makeTestPolicy({});
    const auto& execution_policy = *execution_policy_holder;
#endif
    const auto execution = TestExecution{.policy = &execution_policy};
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    const auto T = test_vector<double>{300, 500, 300, 800};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto r = evaluate(Inconel600_YoungModulus_generic_parallel, output,
                            E.data(), 1, args, contiguous<1>(), 4,
                            mfront_gmp_OutOfBoundsPolicy{}, execution);
    TFEL_TESTS_CHECK_EQUAL(r, 0);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
#if defined(MFRONT_GPMP_TEST_USE_RAJA) && defined(RAJA_CUDA_ACTIVE)
    // destruction of the stream created by makeExplicitTestResource. The
    // generic lambda makes the discarded branch dependent: it is not
    // checked with a resource of the host, which has no stream.
    [this](const auto& res) {
      using Resource = std::remove_cvref_t<decltype(res)>;
      if constexpr (std::is_same_v<Resource, camp::resources::Cuda>) {
        TFEL_TESTS_ASSERT(cudaStreamDestroy(res.get_stream()) == cudaSuccess);
      }
    }(resource);
#endif
  }

  void test2() {
    // asynchronous evaluation: the statistics are given by the calling code,
    // the output status only reports the errors of the batch, and the result
    // is given by the report function once the computations are completed
    // (the backends computing on the host compute synchronously). The second
    // temperature is out of its physical bounds.
    constexpr auto eps = double{1e-2};
    auto statistics = test_vector<mfront_gpmp_Statistics>(1);
    mfront_gpmp_initializeStatistics(statistics.data());
#if defined(MFRONT_GPMP_TEST_USE_KOKKOS)
    auto space = MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE{};
    const auto execution_policy =
        mfront::generic_parallel::makeKokkosBackendExecutionPolicy(
            mfront::generic_parallel::KokkosBackendExecutionOptions<
                MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE>{
                .common_options = {.asynchronous = true},
                .execution_space = space});
#elif defined(MFRONT_GPMP_TEST_USE_RAJA)
    auto resource = getTestResource();
    const auto execution_policy =
        mfront::generic_parallel::makeRAJABackendExecutionPolicy(
            mfront::generic_parallel::RAJABackendExecutionOptions<
                decltype(resource)>{.common_options = {.asynchronous = true},
                                    .resource = resource});
#elif defined(MFRONT_GPMP_TEST_USE_ARCANE)
    auto queue = ArcaneTestRuntime::current->queue;
    const auto execution_policy =
        mfront::generic_parallel::makeArcaneBackendExecutionPolicy(
            queue, {.asynchronous = true});
#else
    const auto execution_policy_holder = makeTestPolicy({.asynchronous = true});
    const auto& execution_policy = *execution_policy_holder;
#endif
    const auto execution = TestExecution{.statistics = statistics.data(),
                                         .policy = &execution_policy};
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(4);
    const auto T = test_vector<double>{300, -1, 500, 800};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto r = evaluate(Inconel600_YoungModulus_generic_parallel, output,
                            E.data(), 1, args, contiguous<1>(), 4,
                            mfront_gmp_OutOfBoundsPolicy{}, execution);
    TFEL_TESTS_CHECK_EQUAL(r, 0);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
#if defined(MFRONT_GPMP_TEST_USE_KOKKOS)
    space.fence();
#elif defined(MFRONT_GPMP_TEST_USE_RAJA)
    resource.wait();
#elif defined(MFRONT_GPMP_TEST_USE_ARCANE)
    queue.barrier();
#endif
    const auto r2 = Inconel600_YoungModulus_generic_parallel_report(
        &output, statistics.data());
    TFEL_TESTS_CHECK_EQUAL(r2, -1);
    TFEL_TESTS_CHECK_EQUAL(output.status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.bounds_status, -1);
    TFEL_TESTS_CHECK_EQUAL(output.number_of_failures, mfront_gmp_size_type{1});
    TFEL_TESTS_CHECK_EQUAL(output.first_failure_index, mfront_gmp_size_type{1});
    TFEL_TESTS_ASSERT(std::isnan(E[1]));
    for (const auto i : {0, 2, 3}) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
  }

  void test3() {
    // stack size of the threads of the device preserved by the execution
    // policy. With a CUDA execution space of Kokkos, an execution policy of
    // RAJA for CUDA or OpenACC on an NVIDIA GPU, a runtime error is reported
    // if it is too small for the kernel. On the host, the option has no
    // effect.
    constexpr auto eps = double{1e-2};
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(2);
    const auto T = test_vector<double>{300, 500};
    const auto args = std::array<const double *, 1u>{T.data()};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto preserve_policy = makeTestPolicy({.preserve_stack_size = true});
    const auto preserve = TestExecution{.policy = preserve_policy.get()};
    const auto f = Inconel600_YoungModulus_generic_parallel;
#if defined(MFRONT_GPMP_TEST_USE_KOKKOS) && defined(KOKKOS_ENABLE_CUDA)
    if constexpr (std::is_same_v<MFRONT_GPMP_TEST_KOKKOS_EXECUTION_SPACE,
                                 Kokkos::Cuda>) {
      auto initial_stack_size = std::size_t{};
      TFEL_TESTS_ASSERT(cudaDeviceGetLimit(&initial_stack_size,
                                           cudaLimitStackSize) == cudaSuccess);
      // the kernels require at least 1024 bytes (see
      // mfront::generic_parallel::gpu::ensureStackSize)
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize, 256) ==
                        cudaSuccess);
      TFEL_TESTS_CHECK_EQUAL(evaluate(f, output, E.data(), 1, args,
                                      contiguous<1>(), 2, policy, preserve),
                             -7);
      TFEL_TESTS_ASSERT(std::string(output.msg).starts_with(
          "the stack size of the threads of the device must be at least"));
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize,
                                           initial_stack_size) == cudaSuccess);
    }
#endif /* defined(MFRONT_GPMP_TEST_USE_KOKKOS) && \
          defined(KOKKOS_ENABLE_CUDA) */
#if defined(MFRONT_GPMP_TEST_USE_RAJA) && defined(RAJA_CUDA_ACTIVE)
    if constexpr (mfront::generic_parallel::raja::isDeviceExecutionPolicy<
                      MFRONT_GPMP_TEST_RAJA_EXECUTION_POLICY>) {
      auto initial_stack_size = std::size_t{};
      TFEL_TESTS_ASSERT(cudaDeviceGetLimit(&initial_stack_size,
                                           cudaLimitStackSize) == cudaSuccess);
      // the kernels require at least 1024 bytes (see
      // mfront::generic_parallel::gpu::ensureStackSize)
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize, 256) ==
                        cudaSuccess);
      TFEL_TESTS_CHECK_EQUAL(evaluate(f, output, E.data(), 1, args,
                                      contiguous<1>(), 2, policy, preserve),
                             -7);
      TFEL_TESTS_ASSERT(std::string(output.msg).starts_with(
          "the stack size of the threads of the device must be at least"));
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize,
                                           initial_stack_size) == cudaSuccess);
    }
#endif /* defined(MFRONT_GPMP_TEST_USE_RAJA) && defined(RAJA_CUDA_ACTIVE) */
#if defined(MFRONT_GPMP_TEST_USE_ARCANE) && defined(ARCCORE_HAS_CUDA)
    if (ArcaneTestRuntime::current->queue.executionPolicy() ==
        Arcane::Accelerator::eExecutionPolicy::CUDA) {
      auto initial_stack_size = std::size_t{};
      TFEL_TESTS_ASSERT(cudaDeviceGetLimit(&initial_stack_size,
                                           cudaLimitStackSize) == cudaSuccess);
      // the kernels require at least 1024 bytes (see
      // mfront::generic_parallel::gpu::ensureStackSize)
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize, 256) ==
                        cudaSuccess);
      TFEL_TESTS_CHECK_EQUAL(evaluate(f, output, E.data(), 1, args,
                                      contiguous<1>(), 2, policy, preserve),
                             -7);
      TFEL_TESTS_ASSERT(std::string(output.msg)
                            .starts_with("the stack size of the threads of the "
                                         "device must be at least"));
      TFEL_TESTS_ASSERT(cudaDeviceSetLimit(cudaLimitStackSize,
                                           initial_stack_size) == cudaSuccess);
    }
#endif /* defined(MFRONT_GPMP_TEST_USE_ARCANE) && defined(ARCCORE_HAS_CUDA) */
#ifdef __NVCOMPILER_OPENACC_GPU
    // stack size of the threads of the GPU of OpenACC, set through the CUDA
    // driver on the context of OpenACC (see `OpenACCStackSize.hxx`)
    {
      auto initial_stack_size = std::size_t{};
      TFEL_TESTS_ASSERT(accessOpenACCStackSize(initial_stack_size, false));
      // a size larger than the current one is refused if the stack size
      // must be preserved, the current one is accepted
      TFEL_TESTS_ASSERT(!mfront::generic_parallel::openacc::setMinimalStackSize(
          initial_stack_size + 1, true));
      TFEL_TESTS_ASSERT(mfront::generic_parallel::openacc::setMinimalStackSize(
          initial_stack_size, true));
      // the kernels require 16384 bytes (`device_stack_size` option of the
      // `openacc-nvhpc-gpu` configuration)
      auto small_stack_size = std::size_t{256};
      TFEL_TESTS_ASSERT(accessOpenACCStackSize(small_stack_size, true));
      TFEL_TESTS_CHECK_EQUAL(evaluate(f, output, E.data(), 1, args,
                                      contiguous<1>(), 2, policy, preserve),
                             -7);
      TFEL_TESTS_ASSERT(std::string(output.msg).starts_with(
          "unable to set the stack size of the threads"));
      TFEL_TESTS_ASSERT(accessOpenACCStackSize(initial_stack_size, true));
    }
#endif /* __NVCOMPILER_OPENACC_GPU */
    // the stack size has been increased by the previous evaluations
    const auto r =
        evaluate(f, output, E.data(), 1, args, contiguous<1>(), 2, policy,
                 preserve);
    TFEL_TESTS_CHECK_EQUAL(r, 0);
    TFEL_TESTS_CHECK_EQUAL(output.status, 0);
    for (std::size_t i = 0; i != E.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(E[i] - fct(T[i])) < eps);
    }
  }

#if defined(MFRONT_GPMP_TEST_USE_ARCANE)
  void test4() {
    // the `arcane` backend has no default policy: a null policy and a policy
    // holding a null queue are invalid arguments
    auto output = mfront_gpmp_OutputStatus{};
    auto E = test_vector<double>(2);
    const auto T = test_vector<double>{300, 500};
    const auto args = std::array<const double*, 1u>{T.data()};
    const auto policy = mfront_gmp_OutOfBoundsPolicy{};
    const auto f = Inconel600_YoungModulus_generic_parallel;
    TFEL_TESTS_CHECK_EQUAL(
        evaluate(f, output, E.data(), 1, args, contiguous<1>(), 2, policy,
                 TestExecution{}),
        -5);
    TFEL_TESTS_ASSERT(std::string_view(output.msg)
                          .ends_with("the arcane backend requires a queue"));
    const auto null_queue =
        mfront::generic_parallel::makeArcaneBackendExecutionPolicy(
            Arcane::Accelerator::RunQueue{});
    TFEL_TESTS_CHECK_EQUAL(
        evaluate(f, output, E.data(), 1, args, contiguous<1>(), 2, policy,
                 TestExecution{.policy = &null_queue}),
        -5);
    TFEL_TESTS_ASSERT(
        std::string_view(output.msg)
            .ends_with("the queue of the execution policy is null"));
  }
#endif /* defined(MFRONT_GPMP_TEST_USE_ARCANE) */
};

TFEL_TESTS_GENERATE_PROXY(BackendTest, "BackendTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main(const int argc, const char *const *argv) {
  if (argc != 2) {
    std::cerr << "no test name given on the command line (argc: " << argc
              << ")";
    return EXIT_FAILURE;
  }
#ifdef MFRONT_GPMP_TEST_USE_KOKKOS
  auto kokkos = Kokkos::ScopeGuard{};
#endif /* MFRONT_GPMP_TEST_USE_KOKKOS */
#ifdef MFRONT_GPMP_TEST_USE_ARCANE
  const auto arcane = ArcaneTestRuntime{};
#endif /* MFRONT_GPMP_TEST_USE_ARCANE */
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput(std::string{argv[1]} + ".xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
