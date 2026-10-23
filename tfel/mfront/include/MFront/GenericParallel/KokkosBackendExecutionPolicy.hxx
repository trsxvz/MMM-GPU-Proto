/*!
 * \file mfront/include/MFront/GenericParallel/KokkosBackendExecutionPolicy.hxx
 * \brief  This file declares the execution policy of the `kokkos` backend of
 * the `generic-parallel` interfaces: a class template on the execution
 * space, whose identifier is made of the name of the backend and of the
 * name of the execution space.
 * \date   22/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_KOKKOSBACKENDEXECUTIONPOLICY_HXX
#define LIB_MFRONT_GENERICPARALLEL_KOKKOSBACKENDEXECUTIONPOLICY_HXX

#include <string>
#include <Kokkos_Core.hpp>
#include "MFront/GenericParallel/BackendExecutionPolicyBase.hxx"

namespace mfront::generic_parallel {

  //! \brief options of the execution policy of the `kokkos` backend
  template <typename ExecutionSpace>
  struct KokkosBackendExecutionOptions {
    //! \brief common options
    BackendExecutionCommonPolicy common_options;
    /*!
     * \brief instance of the execution space on which the kernels are
     * launched, the default instance by default (Kokkos must then be
     * initialized when the options are built)
     */
    ExecutionSpace execution_space = ExecutionSpace{};
  };  // end of struct KokkosBackendExecutionOptions

  /*!
   * \brief execution policy of the `kokkos` backend
   * \tparam ExecutionSpace: execution space, which must be the one selected
   * by the `execution_space` option of the backend
   */
  template <typename ExecutionSpace>
  struct KokkosBackendExecutionPolicy final : BackendExecutionPolicyBase {
    /*!
     * \return the identifier of the policy, for instance "kokkos:Cuda". The
     * string is never destroyed, so that it lives as long as the program,
     * as the identifiers of the other policies (string literals). Contrary
     * to them, this identifier is not a constant expression (the name of
     * the execution space is not one with all the versions of Kokkos): the
     * symbol exporting it is initialized when the library is loaded.
     */
    static const char* identifier() noexcept {
      static const auto* const id =
          new std::string{std::string{"kokkos:"} + ExecutionSpace::name()};
      return id->c_str();
    }
    /*!
     * \brief constructor
     * \param[in] opts: options
     */
    explicit KokkosBackendExecutionPolicy(
        const KokkosBackendExecutionOptions<ExecutionSpace>& opts = {})
        : BackendExecutionPolicyBase(opts.common_options),
          execution_space(opts.execution_space) {}
    const char* getIdentifier() const noexcept override { return identifier(); }
    //! \brief instance of the execution space
    ExecutionSpace execution_space;
  };  // end of struct KokkosBackendExecutionPolicy

  /*!
   * \return an execution policy of the `kokkos` backend
   * \param[in] opts: options
   */
  template <typename ExecutionSpace>
  KokkosBackendExecutionPolicy<ExecutionSpace> makeKokkosBackendExecutionPolicy(
      const KokkosBackendExecutionOptions<ExecutionSpace>& opts = {}) {
    return KokkosBackendExecutionPolicy<ExecutionSpace>{opts};
  }  // end of makeKokkosBackendExecutionPolicy

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_KOKKOSBACKENDEXECUTIONPOLICY_HXX */
