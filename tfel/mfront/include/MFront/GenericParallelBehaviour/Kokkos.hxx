/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/Kokkos.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with Kokkos.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX

#include <type_traits>
#include <Kokkos_Core.hpp>
#include "MFront/GenericParallelBehaviour/Integrate.hxx"
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include "MFront/GenericParallelBehaviour/GPU.hxx"
#endif

namespace mfront::gpb::kokkos {

  namespace internals {

    /*!
     * \brief functor computing the statistics of a batch with
     * `Kokkos::parallel_reduce`.
     * \tparam PointFunction: function object applied at one point
     */
    template <typename PointFunction>
    struct Reduction {
      //! \brief type of the reduced value
      using value_type = BatchStatistics;
      //! \brief apply the function at the `i`-th point
      KOKKOS_FUNCTION void operator()(const size_type i, value_type& s) const {
        s = combine(s, this->f(i));
      }
      //! \brief combine two partial results
      KOKKOS_FUNCTION void join(value_type& s, const value_type& s2) const {
        s = combine(s, s2);
      }
      //! \brief initialize a partial result
      KOKKOS_FUNCTION void init(value_type& s) const {
        s = value_type{.rdt = this->initial_rdt};
      }
      //! \brief function object applied at one point
      PointFunction f;
      //! \brief initial value of the time step scaling factor
      real initial_rdt;
    };

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)

    /*!
     * \brief kernel which is never launched: its attributes give the stack
     * size required to apply a function at one point. Kokkos does not
     * increase the stack size of the threads, which is too small by default
     * for some behaviours.
     */
    template <typename PointFunction>
    __global__ void stackSizeProbe(const PointFunction f) {
      f(0);
    }  // end of stackSizeProbe

#endif /* defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) */

    /*!
     * \brief prepare the execution of the given function in the given
     * execution space.
     * \return true on success
     * \param[out] s: output status
     */
    template <typename ExecutionSpace, typename PointFunction>
    bool prepare(OutputStatus& s) {
#if defined(KOKKOS_ENABLE_CUDA)
      constexpr auto device = std::is_same_v<ExecutionSpace, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
      constexpr auto device = std::is_same_v<ExecutionSpace, Kokkos::HIP>;
#else
      constexpr auto device = false;
#endif
      if constexpr (device) {
        const auto e =
            gpu::internals::ensureStackSize(stackSizeProbe<PointFunction>);
        if (e != decltype(e){}) {
          report(s, -3, "unable to set the stack size of the threads");
          return false;
        }
      }
      static_cast<void>(s);
      return true;
    }  // end of prepare

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * Kokkos must have been initialized by the calling code. This function
   * waits for the completion of the computations.
   *
   * \tparam ExecutionSpace: execution space
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] queue: instance of the execution space (default instance if
   * null)
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \return the status
   */
  template <typename ExecutionSpace, typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            void* const queue,
            const real rdt0,
            const PointFunction& f) {
    if (!Kokkos::is_initialized()) {
      return report(s, -3, "Kokkos is not initialized");
    }
    if (n == 0) {
      return s.status;
    }
    if (!internals::prepare<ExecutionSpace, PointFunction>(s)) {
      return s.status;
    }
    const auto space = queue == nullptr
                           ? ExecutionSpace{}
                           : *static_cast<const ExecutionSpace*>(queue);
    const auto range =
        Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_type>>(space,
                                                                          0, n);
    // the reduction in a scalar is synchronous
    auto r = BatchStatistics{};
    Kokkos::parallel_reduce("mfront::gpb::kokkos::apply", range,
                            internals::Reduction<PointFunction>{f, rdt0}, r);
    return report(s, r);
  }  // end of apply

  namespace internals {

    /*!
     * \return the primitive of the backend, as expected by
     * `mfront::gpb::integrateBatch` and the similar functions
     * \tparam ExecutionSpace: execution space
     * \param[in,out] s: output status
     */
    template <typename ExecutionSpace>
    auto makeApply(OutputStatus& s) {
      return [&s](const size_type n, void* const queue, const real rdt0,
                  const auto& f) {
        return kokkos::apply<ExecutionSpace>(s, n, queue, rdt0, f);
      };
    }  // end of makeApply

  }  // end of namespace internals

  /*!
   * \brief integrate the behaviour at the points of a batch (see
   * `mfront::gpb::integrateBatch`).
   *
   * \tparam ExecutionSpace: execution space
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename ExecutionSpace,
            typename ParametersType>
  int integrate(OutputStatus& s,
                const BatchData& d,
                const ParametersType& parameters,
                const tfel::material::OutOfBoundsPolicy p) {
    return integrateBatch<Behaviour, PointTraits>(
        s, d, parameters, p, internals::makeApply<ExecutionSpace>(s));
  }  // end of integrate

  /*!
   * \brief execute an initialize function at the points of a batch (see
   * `mfront::gpb::executeInitializeFunctionOnBatch`).
   *
   * \tparam ExecutionSpace: execution space
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] inputs: values of the inputs of the initialize function
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename InitializeFunction,
            size_type number_of_inputs,
            typename ExecutionSpace,
            typename ParametersType>
  int executeInitializeFunction(OutputStatus& s,
                                const BatchData& d,
                                const ConstField& inputs,
                                const ParametersType& parameters,
                                const tfel::material::OutOfBoundsPolicy p) {
    return executeInitializeFunctionOnBatch<
        Behaviour, PointTraits, InitializeFunction, number_of_inputs>(
        s, d, inputs, parameters, p, internals::makeApply<ExecutionSpace>(s));
  }  // end of executeInitializeFunction

  /*!
   * \brief execute a post-processing at the points of a batch (see
   * `mfront::gpb::executePostProcessingOnBatch`).
   *
   * \tparam ExecutionSpace: execution space
   * \param[out] s: output status
   * \param[out] outputs: values of the outputs of the post-processing
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename PostProcessing,
            bool use_initial_state,
            size_type number_of_outputs,
            typename ExecutionSpace,
            typename ParametersType>
  int executePostProcessing(OutputStatus& s,
                            const Field& outputs,
                            const BatchData& d,
                            const ParametersType& parameters,
                            const tfel::material::OutOfBoundsPolicy p) {
    return executePostProcessingOnBatch<Behaviour, PointTraits, PostProcessing,
                                        use_initial_state, number_of_outputs>(
        s, outputs, d, parameters, p, internals::makeApply<ExecutionSpace>(s));
  }  // end of executePostProcessing

  /*!
   * \brief rotate the values of the points of a batch of an orthotropic
   * behaviour (see `mfront::gpb::rotateBatch`).
   *
   * \tparam PointRotation: function object rotating the values of a point
   * (see `RotationExecutor`)
   * \tparam ExecutionSpace: execution space
   * \param[out] s: output status
   * \param[in] d: rotation data
   * \return the status
   */
  template <typename PointRotation, typename ExecutionSpace>
  int rotate(OutputStatus& s, const RotationData& d) {
    return rotateBatch<PointRotation>(s, d,
                                      internals::makeApply<ExecutionSpace>(s));
  }  // end of rotate

}  // end of namespace mfront::gpb::kokkos

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_KOKKOS_HXX */
