/*!
 * \file   mfront/include/MFront/GenericParallel/Batch.hxx
 * \brief  This file declares the statistics, the execution and the report
 * of the computations on a batch of points, shared by the
 * `generic-parallel` interfaces for material properties, behaviours and
 * models, and by the primitives of their backends.
 * \date   09/11/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BATCH_HXX
#define LIB_MFRONT_GENERICPARALLEL_BATCH_HXX

#include <bit>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <exception>
#include <string_view>
#include "TFEL/Config/TFELConfig.hxx"
#include "MFront/GenericParallel/StatisticsKeys.hxx"
#include "MFront/GenericParallel/BackendExecutionPolicy.h"

/*!
 * The primitives of the backends (`apply`, see for instance
 * `MFront/GenericParallel/ParallelSTL.hxx`) apply a function at the points
 * of a batch. This function returns the statistics of a point
 * (`BatchStatistics`). The primitives combine them, by atomic operations or
 * by reductions, and report the result through a handler of the statistics
 * given by the interface (`mfront::gpmp::StatisticsHandler` for the
 * material properties), which provides:
 *
 * - the type `OutputStatus` of the output status and the type `Statistics`
 *   of the statistics updated by the devices and given by the calling code
 *   (C structures of the interface);
 * - `initial()`: the statistics of a batch without event, neutral element of
 *   the reductions;
 * - `hasEvents()`: false if the function applied at the points never reports
 *   an event, in which case the statistics are neither combined nor
 *   reported;
 * - `initialize(Statistics&)`: initialization of the statistics updated by
 *   the devices;
 * - `report(OutputStatus&, const Statistics&)`: report of the statistics
 *   updated by the devices, once the computations are completed;
 * - `report(OutputStatus&, Statistics*, const BatchStatistics&)`: report of
 *   the statistics combined by a reduction, or accumulation of those
 *   statistics in the statistics given by the calling code, if any.
 */

namespace mfront::generic_parallel {

  //! \brief type of the indices and of the numbers of points
  using size_type = std::size_t;

  /*!
   * \brief statistics of the computations at one point or at a group of
   * points of a batch, combined by parallel reductions. The keys are built
   * by `makeKey`.
   *
   * The time step scaling factor is a statistics of the computations which
   * may reduce the time step (see `hasTimeStepScalingFactor`): the material
   * properties leave it to its neutral value.
   *
   * \tparam real: numeric type of the time step scaling factor
   */
  template <typename real>
  struct BatchStatistics {
    //! \brief combine two statistics
    TFEL_HOST_DEVICE friend constexpr BatchStatistics combine(
        const BatchStatistics& a, const BatchStatistics& b) noexcept {
      return {.number_of_failures = a.number_of_failures + b.number_of_failures,
              .first_failure = a.first_failure < b.first_failure
                                   ? a.first_failure
                                   : b.first_failure,
              .number_of_unreliable_results = a.number_of_unreliable_results +
                                              b.number_of_unreliable_results,
              .first_unreliable_result =
                  a.first_unreliable_result < b.first_unreliable_result
                      ? a.first_unreliable_result
                      : b.first_unreliable_result,
              .rdt = a.rdt < b.rdt ? a.rdt : b.rdt};
    }
    //! \brief number of points where the computations failed
    size_type number_of_failures = 0;
    //! \brief key of the first failure (see `makeKey`)
    atomic_size_type first_failure = no_event;
    //! \brief number of points where the result is unreliable
    size_type number_of_unreliable_results = 0;
    //! \brief key of the first unreliable result (see `makeKey`)
    atomic_size_type first_unreliable_result = no_event;
    //! \brief time step scaling factor
    real rdt = 1;
  };

  //! \brief function object combining two statistics
  struct CombineStatistics {
    //! \return the combination of two statistics (see `combine`)
    template <typename real>
    TFEL_HOST_DEVICE constexpr BatchStatistics<real> operator()(
        const BatchStatistics<real>& a,
        const BatchStatistics<real>& b) const noexcept {
      return combine(a, b);
    }
  };

  /*!
   * \return an unsigned integer whose ordering is the same as the ordering
   * of the given floating point values. This allows to compute the minimum
   * of floating point values with an integer atomic operation, which is how
   * the devices update the time step scaling factor.
   */
  TFEL_HOST_DEVICE constexpr atomic_size_type toOrderedKey(
      const double v) noexcept {
    const auto b = std::bit_cast<atomic_size_type>(v);
    constexpr auto sign = atomic_size_type{1} << 63;
    return (b & sign) ? ~b : (b | sign);
  }  // end of toOrderedKey

  /*!
   * \brief true if the statistics updated by the devices hold a time step
   * scaling factor, false otherwise, which is the case of the statistics of
   * the material properties
   */
  template <typename Statistics>
  inline constexpr bool hasTimeStepScalingFactor =
      requires(Statistics& s) { s.rdt; };

  /*!
   * \brief update the statistics shared by the threads with the statistics
   * of a point or of a group of points, by atomic operations. The shared
   * statistics are only accessed if an event occurred: a failure, an
   * unreliable result, or a time step scaling factor lower than the initial
   * one.
   * \tparam Atomics: atomic operations of the backend, given by the static
   * functions `add` and `min` on values of type `atomic_size_type`
   * \param[in,out] s: shared statistics
   * \param[in] r: statistics of a point or of a group of points
   * \param[in] rdt0: ordered key of the initial time step scaling factor
   * (see `toOrderedKey`)
   * \return true if the shared statistics were updated
   */
  template <typename Atomics, typename Statistics, typename real>
  TFEL_HOST_DEVICE bool updateStatistics(Statistics& s,
                                         const BatchStatistics<real>& r,
                                         const atomic_size_type rdt0) {
    auto updated = false;
    if (r.number_of_failures != 0) {
      Atomics::add(&(s.number_of_failures),
                   static_cast<atomic_size_type>(r.number_of_failures));
      Atomics::min(&(s.first_failure), r.first_failure);
      updated = true;
    }
    if (r.number_of_unreliable_results != 0) {
      Atomics::add(
          &(s.number_of_unreliable_results),
          static_cast<atomic_size_type>(r.number_of_unreliable_results));
      Atomics::min(&(s.first_unreliable_result), r.first_unreliable_result);
      updated = true;
    }
    if constexpr (hasTimeStepScalingFactor<Statistics>) {
      const auto rdt = toOrderedKey(r.rdt);
      if (rdt < rdt0) {
        Atomics::min(&(s.rdt), rdt);
        updated = true;
      }
    }
    return updated;
  }  // end of updateStatistics

  /*!
   * \brief block of statistics of the synchronous computations on a device,
   * reused from one call to another (see `StatusBlockPool`).
   *
   * The block lives in pinned host memory, accessible by the host and by the
   * device, and the statistics in device memory: the atomic operations of
   * the kernels are thus performed in device memory, which every system
   * supports. The atomic operations on host memory are not reliable: the
   * PCI Express bus only provides the atomic addition, exchange and
   * compare-and-swap operations, and the other ones are not performed on
   * the AMD GPUs connected by PCI Express without page migration (XNACK),
   * where managed memory is host memory. Managed memory may also not be
   * accessed by the host while kernels run (devices without concurrent
   * access to managed memory, for instance under Windows).
   *
   * A thread updating the statistics also sets `events`, by a plain store
   * (see `signalEvents`). Once the computations are completed, the host
   * reads the statistics only if `events` is set: a call without event
   * performs neither copy nor additional synchronization. If `events` is
   * not set, the statistics hold their neutral value (see
   * `getNeutralStatistics`) and no transfer of the block is in progress:
   * the copies of the block, performed after an event only, are
   * synchronous. A block whose `events` is set, after an event or at its
   * allocation, is reset to the neutral statistics before its use.
   *
   * \tparam Statistics: statistics updated by the devices
   */
  template <typename Statistics>
  struct StatisticsBlock {
    //! \brief statistics updated by the device, in device memory
    Statistics* statistics;
    //! \brief non zero if the statistics may differ from their neutral value
    int events;
  };

  //! \brief signal that the statistics of a block were updated
  TFEL_HOST_DEVICE inline void signalEvents(int* const events) noexcept {
    if (events != nullptr) {
      *static_cast<volatile int*>(events) = 1;
    }
  }  // end of signalEvents

  /*!
   * \return true if the statistics of a block may differ from their neutral
   * value. The flag is read again at each call: it is written by the device.
   */
  template <typename Statistics>
  bool holdsEvents(const StatisticsBlock<Statistics>& b) noexcept {
    return *static_cast<const volatile int*>(&(b.events)) != 0;
  }  // end of holdsEvents

  //! \brief state that the statistics of a block hold their neutral value
  template <typename Statistics>
  void clearEvents(StatisticsBlock<Statistics>& b) noexcept {
    *static_cast<volatile int*>(&(b.events)) = 0;
  }  // end of clearEvents

  /*!
   * \return the neutral value of the statistics of a block: the initial
   * statistics of a call, without event, except the time step scaling
   * factor, which is set to the largest key. The devices only write a time
   * step scaling factor lower than the initial one of the call (see
   * `updateStatistics`), which is restored by `mergeBlockStatistics`.
   * \param[in] s: initial statistics of a call
   */
  template <typename Statistics>
  Statistics getNeutralStatistics(Statistics s) noexcept {
    if constexpr (hasTimeStepScalingFactor<Statistics>) {
      s.rdt = no_event;
    }
    return s;
  }  // end of getNeutralStatistics

  /*!
   * \return the statistics of a call, given its initial statistics and the
   * statistics of the block, read once the computations are completed
   * \param[in] initial: initial statistics of the call
   * \param[in] s: statistics of the block
   */
  template <typename Statistics>
  Statistics mergeBlockStatistics(const Statistics& initial,
                                  Statistics s) noexcept {
    if constexpr (hasTimeStepScalingFactor<Statistics>) {
      if (initial.rdt < s.rdt) {
        s.rdt = initial.rdt;
      }
    }
    return s;
  }  // end of mergeBlockStatistics

  /*!
   * \brief execution of the computations on a batch, as given by the
   * calling code: the policy, checked by the backend against its own policy
   * (see `mfront::generic_parallel::checkBackendExecutionPolicy`), and the
   * statistics
   * \tparam Statistics: statistics updated by the devices
   */
  template <typename Statistics>
  struct BatchExecution {
    //! \brief execution policy, or null for the default policy
    const mfront_gp_AbstractBackendExecutionPolicy* policy;
    //! \brief statistics, or null
    Statistics* statistics;
  };  // end of struct BatchExecution

  /*!
   * \brief set the status and the error message of an output status
   * \param[out] s: output status
   * \param[in] status: status
   * \param[in] msg: error message
   * \return the status
   */
  template <typename OutputStatus>
  int report(OutputStatus& s,
             const int status,
             const std::string_view msg) noexcept {
    constexpr auto msize = sizeof(s.msg) - 1;
    const auto n = msg.size() < msize ? msg.size() : msize;
    s.status = status;
    std::memcpy(s.msg, msg.data(), n);
    s.msg[n] = '\0';
    return status;
  }  // end of report

  /*!
   * \brief call the body of a generated function of batches, so that no
   * exception goes through its C interface: an exception is reported as an
   * error of the backend (status `-7`). The output status is initialized
   * before the call, and `errno` is restored on return: it may be modified
   * by the runtimes of the backends and by the computations on the host.
   * \return the status
   * \param[out] s: output status
   * \param[in] initialize: function initializing the output status
   * \param[in] f: body of the function, returning the status
   */
  template <typename OutputStatus, typename Initialize, typename Function>
  int handleExceptions(OutputStatus* const s,
                       const Initialize& initialize,
                       Function&& f) noexcept {
    if (s == nullptr) {
      return -5;
    }
    struct ErrnoGuard {
      const int value = errno;
      ~ErrnoGuard() { errno = this->value; }
    };
    [[maybe_unused]] const auto errno_guard = ErrnoGuard{};
    initialize(*s);
    try {
      return f();
    } catch (std::exception& e) {
      return report(*s, -7, e.what());
    } catch (...) {
      return report(*s, -7, "unknown exception");
    }
  }  // end of handleExceptions

}  // end of namespace mfront::generic_parallel

#endif /* LIB_MFRONT_GENERICPARALLEL_BATCH_HXX */
