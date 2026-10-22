/*!
 * \file   MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx
 * \brief  Policy of the `kokkos` backend: a class template on the execution
 * space, whose identifier is made of the name of the backend and of
 * `ExecutionSpace::name()`. Defined in this header only.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_KOKKOSBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_KOKKOSBACKENDEXECUTIONPOLICY_HXX

#include <string>
#include "StandIns.hxx"
#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  //! \brief options of the `kokkos` backend
  template <typename ExecutionSpace>
  struct KokkosBackendExecutionOptions {
    BackendExecutionCommonPolicy common_options;
    //! \brief instance of the execution space, the default one by default
    ExecutionSpace execution_space = ExecutionSpace{};
  };  // end of struct KokkosBackendExecutionOptions

  //! \brief policy of the `kokkos` backend
  template <typename ExecutionSpace>
  struct KokkosBackendExecutionPolicy final : BackendExecutionPolicyBase {
    static const char* identifier() noexcept {
      static const auto id = std::string{"kokkos:"} + ExecutionSpace::name();
      return id.c_str();
    }
    explicit KokkosBackendExecutionPolicy(
        const KokkosBackendExecutionOptions<ExecutionSpace>& opts)
        : BackendExecutionPolicyBase(opts.common_options),
          execution_space(opts.execution_space) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
    //! \brief instance of the execution space
    ExecutionSpace execution_space;
  };  // end of struct KokkosBackendExecutionPolicy

  template <typename ExecutionSpace>
  KokkosBackendExecutionPolicy<ExecutionSpace> makeKokkosBackendExecutionPolicy(
      const KokkosBackendExecutionOptions<ExecutionSpace>& opts = {}) {
    return KokkosBackendExecutionPolicy<ExecutionSpace>{opts};
  }

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_KOKKOSBACKENDEXECUTIONPOLICY_HXX */
