/*!
 * \file   src/caller.cxx
 * \brief  Creates the policies, as MGIS or a user code, and checks the
 * answers of the generated functions.
 */

#include <cstdio>
#include <cstdint>
#include "Loader.h"
#include "Caller.h"
#include "MFront/GenericParallel/BatchData.h"
#include "MFront/GenericParallel/CUDABackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/ParallelSTLBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/SyclBackendExecutionPolicy.hxx"
#include "MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx"

namespace mgp = mfront::generic_parallel;

namespace {

  //! \brief result of a call of a generated function
  struct Call {
    int status = 0;
    bool all_points_computed = false;
    int computed_points = -1;
    mfront_gp_ReproducerReport report = {-1, -1, -1, -1, -1};
  };

  Call call(const mfront_gp_GeneratedFunction f,
            const mgp::AbstractBackendExecutionPolicy* const p,
            const bool with_statistics) {
    auto r = Call{};
    int statuses[4] = {0, 0, 0, 0};
    auto s = mfront_gp_Statistics{-1};
    const auto d = mfront_gp_BatchData{
        4, statuses, with_statistics ? &s : nullptr, p, &(r.report)};
    r.status = f(&d);
    r.all_points_computed = true;
    for (const auto st : statuses) {
      r.all_points_computed = r.all_points_computed && (st == 1);
    }
    r.computed_points = s.number_of_computed_points;
    return r;
  }

  struct Checks {
    int failures = 0;
    void check(const bool b, const char* const what) {
      std::printf("  %-66s %s\n", what, b ? "ok" : "FAILED");
      if (!b) {
        ++(this->failures);
      }
    }
  };

  template <typename F>
  F getFunction(const reproducer_library l, const char* const n) {
    const auto f = reinterpret_cast<F>(reproducer_symbol(l, n));
    if (f == nullptr) {
      std::fprintf(stderr, "no function '%s'\n", n);
    }
    return f;
  }

  int checkGeneratedLibrary(const char* const path, const bool c_types_only) {
    auto checks = Checks{};
    const auto l = reproducer_load(path);
    if (l == nullptr) {
      reproducer_print_error(path);
      return 1;
    }
    const auto cuda =
        getFunction<mfront_gp_GeneratedFunction>(l, "reproducer_cuda");
    const auto stlpar =
        getFunction<mfront_gp_GeneratedFunction>(l, "reproducer_stlpar");
    const auto sycl =
        getFunction<mfront_gp_GeneratedFunction>(l, "reproducer_sycl");
    const auto kokkos_cuda =
        getFunction<mfront_gp_GeneratedFunction>(l, "reproducer_kokkos_cuda");
    if ((cuda == nullptr) || (stlpar == nullptr) || (sycl == nullptr) ||
        (kokkos_cuda == nullptr)) {
      return 1;
    }
    const auto stream = reinterpret_cast<cudaStream_t>(std::intptr_t{0x1234});
    const auto cuda_policy = mgp::makeCUDABackendExecutionPolicy(
        {.common_options = {.preserve_stack_size = true}, .stream = stream});
    const auto cuda_asynchronous_policy = mgp::makeCUDABackendExecutionPolicy(
        {.common_options = {.asynchronous = true}, .stream = stream});
    const auto stlpar_asynchronous_policy =
        mgp::makeParallelSTLBackendExecutionPolicy({.asynchronous = true});
    std::printf("cuda (class, C types)\n");
    {
      const auto r = call(cuda, nullptr, false);
      checks.check((r.status == 0) && r.all_points_computed &&
                       (r.report.default_policy == 1) && (r.report.queue == 0),
                   "null policy: default policy, default stream");
    }
    {
      const auto r = call(cuda, &cuda_policy, false);
      checks.check((r.status == 0) && r.all_points_computed &&
                       (r.report.default_policy == 0) &&
                       (r.report.preserve_stack_size == 1) &&
                       (r.report.asynchronous == 0) &&
                       (r.report.queue == 0x1234),
                   "policy: identifier, stream and options read");
    }
    {
      const auto r = call(cuda, &cuda_asynchronous_policy, true);
      checks.check((r.status == 0) && (r.report.asynchronous == 1) &&
                       (r.computed_points == 4),
                   "asynchronous call with statistics");
    }
    checks.check(call(cuda, &cuda_asynchronous_policy, false).status == -5,
                 "asynchronous call without statistics refused (-5)");
    checks.check(call(cuda, &stlpar_asynchronous_policy, true).status == -5,
                 "policy of the stlpar backend refused (-5)");
    std::printf("stlpar (class, common options only)\n");
    {
      const auto r = call(stlpar, nullptr, false);
      checks.check((r.status == 0) && r.all_points_computed &&
                       (r.report.default_policy == 1),
                   "null policy: default policy");
    }
    {
      const auto r = call(stlpar, &stlpar_asynchronous_policy, true);
      checks.check((r.status == 0) && r.all_points_computed &&
                       (r.report.asynchronous == 1) && (r.computed_points == 4),
                   "asynchronous call: computed, statistics filled");
    }
    checks.check(call(stlpar, &cuda_policy, false).status == -5,
                 "policy of the cuda backend refused (-5)");
    if (!c_types_only) {
      std::printf("sycl (class holding a C++ object)\n");
      {
        const auto queue = sycl::queue{3};
        const auto sycl_policy = mgp::makeSyclBackendExecutionPolicy(
            queue, {.preserve_stack_size = true});
        const auto r = call(sycl, &sycl_policy, false);
        checks.check((r.status == 0) && r.all_points_computed &&
                         (r.report.preserve_stack_size == 1) &&
                         (r.report.queue == 3),
                     "policy: identifier, queue and options read");
        // the queue of the caller, the one of the policy and the copy
        // made by the generated function
        checks.check(r.report.queue_owners == 3,
                     "queue copied by the generated function (3 owners)");
        checks.check(queue.use_count() == 2,
                     "copy destroyed by the generated function (2 owners)");
      }
      checks.check(call(sycl, nullptr, false).status == -5,
                   "null policy refused, the queue is required (-5)");
      checks.check(call(sycl, &cuda_policy, false).status == -5,
                   "policy of the cuda backend refused (-5)");
      std::printf(
          "kokkos (class template on an execution space without "
          "visibility)\n");
      const auto kokkos_cuda_policy = mgp::makeKokkosBackendExecutionPolicy(
          mgp::KokkosBackendExecutionOptions<Kokkos::Cuda>{
              .common_options = {.preserve_stack_size = true},
              .execution_space = Kokkos::Cuda{7}});
      const auto kokkos_serial_policy =
          mgp::makeKokkosBackendExecutionPolicy<Kokkos::Serial>();
      {
        const auto r = call(kokkos_cuda, nullptr, false);
        checks.check((r.status == 0) && r.all_points_computed &&
                         (r.report.default_policy == 1) &&
                         (r.report.queue == 0),
                     "null policy: default instance");
      }
      {
        const auto r = call(kokkos_cuda, &kokkos_cuda_policy, false);
        checks.check((r.status == 0) && r.all_points_computed &&
                         (r.report.preserve_stack_size == 1) &&
                         (r.report.queue == 7),
                     "policy: identifier, instance and options read");
      }
      checks.check(call(kokkos_cuda, &kokkos_serial_policy, false).status == -5,
                   "policy of another execution space refused (-5)");
      checks.check(call(kokkos_cuda, &cuda_policy, false).status == -5,
                   "policy of the cuda backend refused (-5)");
      checks.check(call(cuda, &kokkos_cuda_policy, false).status == -5,
                   "kokkos policy refused by the cuda backend (-5)");
    }
    return checks.failures;
  }

  //! \brief for information: answers of `dynamic_cast`, not checked
  void printDynamicCastResults(const char* const path) {
    using F = int (*)(const mfront_gp_AbstractBackendExecutionPolicy*);
    const auto l = reproducer_load(path);
    if (l == nullptr) {
      reproducer_print_error(path);
      return;
    }
    const auto cuda = getFunction<F>(l, "reproducer_dynamic_cast_cuda");
    const auto kokkos_cuda =
        getFunction<F>(l, "reproducer_dynamic_cast_kokkos_cuda");
    if ((cuda == nullptr) || (kokkos_cuda == nullptr)) {
      return;
    }
    const auto cuda_policy = mgp::makeCUDABackendExecutionPolicy();
    const auto kokkos_cuda_policy =
        mgp::makeKokkosBackendExecutionPolicy<Kokkos::Cuda>();
    const auto kokkos_serial_policy =
        mgp::makeKokkosBackendExecutionPolicy<Kokkos::Serial>();
    const auto print = [](const char* const what, const bool b) {
      std::printf("  %-66s %s\n", what, b ? "yes" : "no");
    };
    std::printf("for information, dynamic_cast (not used by the design)\n");
    print("cuda policy recognised", cuda(&cuda_policy));
    print("kokkos policy recognised (class template)",
          kokkos_cuda(&kokkos_cuda_policy));
    print("kokkos policy of another execution space refused",
          !kokkos_cuda(&kokkos_serial_policy));
  }

}  // end of namespace

extern "C" REPRODUCER_CALLER_VISIBILITY int reproducer_run(
    const char* const generated,
    const char* const dynamic_cast_library,
    const int c_types_only) {
#ifdef MFRONT_GENERICPARALLEL_IDENTIFIER_FIRST
  std::printf("getIdentifier declared before the destructor\n");
#else
  std::printf("getIdentifier declared after the destructor\n");
#endif
  std::printf("generated library: %s\n", generated);
  const auto failures = checkGeneratedLibrary(generated, c_types_only != 0);
  if (dynamic_cast_library != nullptr) {
    printDynamicCastResults(dynamic_cast_library);
  }
  std::printf("%s\n", failures == 0 ? "all the checks succeeded"
                                    : "SOME CHECKS FAILED");
  return failures;
}
