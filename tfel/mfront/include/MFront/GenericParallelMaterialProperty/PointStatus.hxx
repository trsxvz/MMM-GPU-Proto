/*!
 * \file   MFront/GenericParallelMaterialProperty/PointStatus.hxx
 * \brief  This file declares the status of the evaluation of a material
 * property at a point, and its merge in the statistics of an evaluation.
 * \date   08/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELMATERIALPROPERTY_POINTSTATUS_HXX
#define LIB_MFRONT_GENERICPARALLELMATERIALPROPERTY_POINTSTATUS_HXX

#include "TFEL/Config/TFELConfig.hxx"
#include "MFront/GenericParallelMaterialProperty/MaterialProperty.h"

namespace mfront::gpmp {

  //! \brief value of an index or of an identifier of a check meaning none
  inline constexpr unsigned long long no_event = ~0ull;

  /*!
   * \brief status of the evaluation of a material property at a point,
   * computed by the kernels. With atomic operations, it is merged in the
   * statistics of the evaluation only if an event occurred, so that the
   * cost of the reporting is null when the evaluation succeeds. The
   * reductions merge the statuses of all the points.
   */
  struct PointStatus {
    //! \brief true if the evaluation failed
    bool failure = false;
    //! \brief true if a bound is violated with the warning policy
    bool unreliable = false;
    /*!
     * \brief identifier of the first check which failed (see the
     * `first_violation` member of `mfront_gpmp_Statistics`)
     */
    unsigned long long violation = no_event;
    /*!
     * \brief report a check which failed
     * \param[in] c: identifier of the check
     * \param[in] f: true if the evaluation fails, false if the result is
     * computed but must be used with caution
     */
    TFEL_HOST_DEVICE constexpr void report(const unsigned long long c,
                                           const bool f) noexcept {
      this->violation = (c < this->violation) ? c : this->violation;
      if (f) {
        this->failure = true;
      } else {
        this->unreliable = true;
      }
    }  // end of report
    //! \return true if an event occurred
    [[nodiscard]] TFEL_HOST_DEVICE constexpr bool hasEvent() const noexcept {
      return this->violation != no_event;
    }  // end of hasEvent
  };   // end of struct PointStatus

  /*!
   * \return the statistics associated with the status of a point, as merged
   * by the backends based on reductions
   * \param[in] i: index of the point
   * \param[in] s: status of the point
   */
  [[nodiscard]] TFEL_HOST_DEVICE constexpr Statistics makeStatistics(
      const unsigned long long i, const PointStatus& s) noexcept {
    return Statistics{
        .number_of_failures = s.failure ? 1ull : 0ull,
        .first_failure_index = s.failure ? i : no_event,
        .number_of_unreliable_results =
            ((!s.failure) && (s.unreliable)) ? 1ull : 0ull,
        .first_violation = s.violation};
  }  // end of makeStatistics

  /*!
   * \return the merge of two statistics
   * \param[in] a: first statistics
   * \param[in] b: second statistics
   */
  [[nodiscard]] TFEL_HOST_DEVICE constexpr Statistics merge(
      const Statistics& a, const Statistics& b) noexcept {
    return Statistics{
        .number_of_failures = a.number_of_failures + b.number_of_failures,
        .first_failure_index = (a.first_failure_index < b.first_failure_index)
                                   ? a.first_failure_index
                                   : b.first_failure_index,
        .number_of_unreliable_results =
            a.number_of_unreliable_results + b.number_of_unreliable_results,
        .first_violation = (a.first_violation < b.first_violation)
                               ? a.first_violation
                               : b.first_violation};
  }  // end of merge

}  // end of namespace mfront::gpmp

#endif /* LIB_MFRONT_GENERICPARALLELMATERIALPROPERTY_POINTSTATUS_HXX */
