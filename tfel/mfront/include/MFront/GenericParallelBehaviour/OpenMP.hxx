/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/OpenMP.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with OpenMP, on the host or on a device (target offloading).
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENMP_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENMP_HXX

#include <cstddef>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::openmp {

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * If the computations are offloaded, the fields of the batch must be
   * accessible from the device, for instance allocated in unified or
   * managed memory (`nvc++ -gpu=mem:managed`).
   *
   * \tparam offload: if true, the computations are offloaded to the default
   * device (`target` construct). Otherwise, they are performed by the
   * threads of the host (`parallel` construct) and the code does not depend
   * on the support of offloading by the OpenMP runtime.
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are offloaded to an NVIDIA
   * GPU through the CUDA runtime (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <bool offload, typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            const real rdt0,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    // the statistics are reduced member by member
    auto number_of_failures = size_type{};
    auto number_of_unreliable_results = size_type{};
    auto rdt = rdt0;
    if constexpr (offload) {
      if (!mfront::generic_parallel::setMinimalCUDAStackSize(
              device_stack_size)) {
        return report(s, -3, "unable to set the stack size of the threads");
      }
#pragma omp target teams distribute parallel for map(to : f)        \
    reduction(+ : number_of_failures, number_of_unreliable_results) \
    reduction(min : rdt)
      for (size_type i = 0; i < n; ++i) {
        const auto r = f(i);
        number_of_failures += r.number_of_failures;
        number_of_unreliable_results += r.number_of_unreliable_results;
        rdt = r.rdt < rdt ? r.rdt : rdt;
      }
    } else {
      static_cast<void>(device_stack_size);
#pragma omp parallel for reduction(+ : number_of_failures,           \
                                       number_of_unreliable_results) \
    reduction(min : rdt)
      for (size_type i = 0; i < n; ++i) {
        const auto r = f(i);
        number_of_failures += r.number_of_failures;
        number_of_unreliable_results += r.number_of_unreliable_results;
        rdt = r.rdt < rdt ? r.rdt : rdt;
      }
    }
    return report(s, BatchStatistics{.number_of_failures = number_of_failures,
                                     .number_of_unreliable_results =
                                         number_of_unreliable_results,
                                     .rdt = rdt});
  }  // end of apply

  namespace internals {

    /*!
     * \return the primitive of the backend, as expected by
     * `mfront::gpb::integrateBatch` and the similar functions
     * \tparam offload: if true, the computations are offloaded
     * \param[in,out] s: output status
     * \param[in] device_stack_size: minimal stack size of the threads of the
     * device (see `apply`)
     */
    template <bool offload>
    auto makeApply(OutputStatus& s, const std::size_t device_stack_size) {
      return [&s, device_stack_size](const size_type n, void* const,
                                     const real rdt0, const auto& f) {
        return openmp::apply<offload>(s, n, rdt0, f, device_stack_size);
      };
    }  // end of makeApply

  }  // end of namespace internals

  /*!
   * \brief integrate the behaviour at the points of a batch (see
   * `mfront::gpb::integrateBatch`).
   *
   * \tparam offload: if true, the computations are offloaded
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            bool offload,
            typename ParametersType>
  int integrate(OutputStatus& s,
                const BatchData& d,
                const ParametersType& parameters,
                const tfel::material::OutOfBoundsPolicy p,
                const std::size_t device_stack_size = 0) {
    return integrateBatch<Behaviour, PointTraits>(
        s, d, parameters, p,
        internals::makeApply<offload>(s, device_stack_size));
  }  // end of integrate

  /*!
   * \brief execute an initialize function at the points of a batch (see
   * `mfront::gpb::executeInitializeFunctionOnBatch`).
   *
   * \tparam offload: if true, the computations are offloaded
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] inputs: values of the inputs of the initialize function
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename InitializeFunction,
            size_type number_of_inputs,
            bool offload,
            typename ParametersType>
  int executeInitializeFunction(OutputStatus& s,
                                const BatchData& d,
                                const ConstField& inputs,
                                const ParametersType& parameters,
                                const tfel::material::OutOfBoundsPolicy p,
                                const std::size_t device_stack_size = 0) {
    return executeInitializeFunctionOnBatch<Behaviour, PointTraits,
                                            InitializeFunction,
                                            number_of_inputs>(
        s, d, inputs, parameters, p,
        internals::makeApply<offload>(s, device_stack_size));
  }  // end of executeInitializeFunction

  /*!
   * \brief execute a post-processing at the points of a batch (see
   * `mfront::gpb::executePostProcessingOnBatch`).
   *
   * \tparam offload: if true, the computations are offloaded
   * \param[out] s: output status
   * \param[out] outputs: values of the outputs of the post-processing
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename PostProcessing,
            bool use_initial_state,
            size_type number_of_outputs,
            bool offload,
            typename ParametersType>
  int executePostProcessing(OutputStatus& s,
                            const Field& outputs,
                            const BatchData& d,
                            const ParametersType& parameters,
                            const tfel::material::OutOfBoundsPolicy p,
                            const std::size_t device_stack_size = 0) {
    return executePostProcessingOnBatch<Behaviour, PointTraits, PostProcessing,
                                        use_initial_state, number_of_outputs>(
        s, outputs, d, parameters, p,
        internals::makeApply<offload>(s, device_stack_size));
  }  // end of executePostProcessing

  /*!
   * \brief rotate the values of the points of a batch of an orthotropic
   * behaviour (see `mfront::gpb::rotateBatch`).
   *
   * \tparam offload: if true, the computations are offloaded
   * \param[out] s: output status
   * \param[in] d: rotation data
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename PointRotation, bool offload>
  int rotate(OutputStatus& s,
             const RotationData& d,
             const std::size_t device_stack_size = 0) {
    return rotateBatch<PointRotation>(
        s, d, internals::makeApply<offload>(s, device_stack_size));
  }  // end of rotate

}  // end of namespace mfront::gpb::openmp

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_OPENMP_HXX */
