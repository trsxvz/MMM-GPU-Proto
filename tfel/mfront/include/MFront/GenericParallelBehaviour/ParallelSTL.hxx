/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/ParallelSTL.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with the parallel algorithms of the C++ standard library.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_PARALLELSTL_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_PARALLELSTL_HXX

#include <cstddef>
#include <iterator>
#include <numeric>
#include <execution>
#include "MFront/GenericParallel/CountingIterator.hxx"
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::stlpar {

  //! \brief a simple alias
  using CountingIterator =
      mfront::generic_parallel::CountingIterator<size_type>;

  static_assert(std::random_access_iterator<CountingIterator>);

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \tparam ExecutionPolicy: execution policy
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the parallel algorithms are offloaded to an
   * NVIDIA GPU by `nvc++ -stdpar=gpu` (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <typename PointFunction, typename ExecutionPolicy>
  int apply(OutputStatus& s,
            const size_type n,
            const real rdt0,
            const PointFunction& f,
            ExecutionPolicy&& policy,
            const std::size_t device_stack_size = 0) {
    if (!mfront::generic_parallel::setMinimalCUDAStackSize(device_stack_size)) {
      return report(s, -3, "unable to set the stack size of the threads");
    }
    const auto r = std::transform_reduce(
        policy, CountingIterator{0}, CountingIterator{n},
        BatchStatistics{.rdt = rdt0}, CombineStatistics{}, f);
    return report(s, r);
  }  // end of apply

  namespace internals {

    /*!
     * \return the primitive of the backend, as expected by
     * `mfront::gpb::integrateBatch` and the similar functions
     * \param[in,out] s: output status
     * \param[in] policy: execution policy
     * \param[in] device_stack_size: minimal stack size of the threads of the
     * device (see `apply`)
     */
    template <typename ExecutionPolicy>
    auto makeApply(OutputStatus& s,
                   const ExecutionPolicy policy,
                   const std::size_t device_stack_size) {
      return [&s, policy, device_stack_size](const size_type n, void* const,
                                             const real rdt0, const auto& f) {
        return stlpar::apply(s, n, rdt0, f, policy, device_stack_size);
      };
    }  // end of makeApply

  }  // end of namespace internals

  /*!
   * \brief integrate the behaviour at the points of a batch (see
   * `mfront::gpb::integrateBatch`).
   *
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename ParametersType,
            typename ExecutionPolicy>
  int integrate(OutputStatus& s,
                const BatchData& d,
                const ParametersType& parameters,
                const tfel::material::OutOfBoundsPolicy p,
                ExecutionPolicy&& policy,
                const std::size_t device_stack_size = 0) {
    return integrateBatch<Behaviour, PointTraits>(
        s, d, parameters, p,
        internals::makeApply(s, policy, device_stack_size));
  }  // end of integrate

  /*!
   * \brief execute an initialize function at the points of a batch (see
   * `mfront::gpb::executeInitializeFunctionOnBatch`).
   *
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] inputs: values of the inputs of the initialize function
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename InitializeFunction,
            size_type number_of_inputs,
            typename ExecutionPolicy,
            typename ParametersType>
  int executeInitializeFunction(OutputStatus& s,
                                const BatchData& d,
                                const ConstField& inputs,
                                const ParametersType& parameters,
                                const tfel::material::OutOfBoundsPolicy p,
                                ExecutionPolicy&& policy,
                                const std::size_t device_stack_size = 0) {
    return executeInitializeFunctionOnBatch<
        Behaviour, PointTraits, InitializeFunction, number_of_inputs>(
        s, d, inputs, parameters, p,
        internals::makeApply(s, policy, device_stack_size));
  }  // end of executeInitializeFunction

  /*!
   * \brief execute a post-processing at the points of a batch (see
   * `mfront::gpb::executePostProcessingOnBatch`).
   *
   * \param[out] s: output status
   * \param[out] outputs: values of the outputs of the post-processing
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename PostProcessing,
            bool use_initial_state,
            size_type number_of_outputs,
            typename ExecutionPolicy,
            typename ParametersType>
  int executePostProcessing(OutputStatus& s,
                            const Field& outputs,
                            const BatchData& d,
                            const ParametersType& parameters,
                            const tfel::material::OutOfBoundsPolicy p,
                            ExecutionPolicy&& policy,
                            const std::size_t device_stack_size = 0) {
    return executePostProcessingOnBatch<Behaviour, PointTraits, PostProcessing,
                                        use_initial_state, number_of_outputs>(
        s, outputs, d, parameters, p,
        internals::makeApply(s, policy, device_stack_size));
  }  // end of executePostProcessing

  /*!
   * \brief rotate the values of the points of a batch of an orthotropic
   * behaviour (see `mfront::gpb::rotateBatch`).
   *
   * \param[out] s: output status
   * \param[in] d: rotation data
   * \param[in] policy: execution policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename PointRotation, typename ExecutionPolicy>
  int rotate(OutputStatus& s,
             const RotationData& d,
             ExecutionPolicy&& policy,
             const std::size_t device_stack_size = 0) {
    return rotateBatch<PointRotation>(
        s, d, internals::makeApply(s, policy, device_stack_size));
  }  // end of rotate

}  // end of namespace mfront::gpb::stlpar

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_PARALLELSTL_HXX */
