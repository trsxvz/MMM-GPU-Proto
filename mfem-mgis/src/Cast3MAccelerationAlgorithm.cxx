/*!
 * \file   src/Cast3MAccelerationAlgorithm.cxx
 * \brief  This file implements the acceleration of the fixed point
 * iterations of Cast3M
 * \date   10/10/2026
 */

#include <cmath>
#include <limits>
#include <utility>
#include <string>
#include <vector>
#include "MFEMMGIS/Parameters.hxx"
#include "AccelerationAlgorithms.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief acceleration of Cast3M, as `mtest::CastemAccelerationAlgorithm`:
   * every `period` iterations from the `trigger`-th one, the new estimate is
   * the combination of the last three estimates whose combination of the
   * associated residuals has the smallest norm (or of the last two estimates
   * if the last residuals are almost colinear)
   */
  struct Cast3MAccelerationAlgorithm final : AccelerationAlgorithm {
    /*!
     * \param[in] t: first accelerated iteration
     * \param[in] p: number of iterations between two accelerations
     */
    Cast3MAccelerationAlgorithm(const size_type t, const size_type p)
        : trigger(t), period(p) {}
    void initialize(const size_type n, ScalarProduct p) override {
      this->dot = std::move(p);
      for (auto* const v : {&(this->u0), &(this->u1), &(this->u2), &(this->r0),
                            &(this->r1), &(this->r2), &(this->n0), &(this->n1),
                            &(this->tmp0), &(this->tmp1)}) {
        v->UseDevice(true);
        v->SetSize(static_cast<int>(n));
        *v = real{0};
      }
    }
    void execute(mfem::Vector& x,
                 const mfem::Vector&,
                 const mfem::Vector& r,
                 const real r_eps,
                 const size_type i) override {
      const auto eps = 100 * r_eps * std::numeric_limits<real>::epsilon();
      std::swap(this->u0, this->u1);
      std::swap(this->u1, this->u2);
      std::swap(this->r0, this->r1);
      std::swap(this->r1, this->r2);
      this->u2 = x;
      this->r2 = r;
      if ((i < this->trigger) || ((i - this->trigger) % this->period != 0)) {
        return;
      }
      subtract(this->r1, this->r0, this->tmp0);
      subtract(this->r2, this->r0, this->tmp1);
      const auto nr0 = std::sqrt(this->dot(this->tmp0, this->tmp0));
      if (nr0 <= eps) {
        return;
      }
      this->n0 = this->tmp0;
      this->n0 /= nr0;
      const auto ntmp1 = this->dot(this->tmp1, this->n0);
      this->tmp1.Add(-ntmp1, this->n0);
      const auto nr1 = std::sqrt(this->dot(this->tmp1, this->tmp1));
      if (nr1 > std::abs(ntmp1) / 10) {
        this->n1 = this->tmp1;
        this->n1 /= nr1;
        const auto p0 = -this->dot(this->r0, this->n0);
        const auto p1 = -this->dot(this->r0, this->n1);
        // the combination (1 - c1 - c2) r0 + c1 r1 + c2 r2 is the projection
        // of the null vector on the plane of the residuals: the same
        // combination is applied to the estimates
        const auto c2 = p1 / nr1;
        const auto c1 = (p0 - ntmp1 * c2) / nr0;
        add(1 - c2 - c1, this->u0, c1, this->u1, x);
        x.Add(c2, this->u2);
      } else {
        // the last residuals are almost colinear
        const auto c0 = -this->dot(this->r0, this->n0) / nr0;
        add(1 - c0, this->u0, c0, this->u1, x);
      }
    }

   private:
    //! \brief first accelerated iteration
    const size_type trigger;
    //! \brief number of iterations between two accelerations
    const size_type period;
    //! \brief scalar product
    ScalarProduct dot;
    //! \brief last three estimates
    mfem::Vector u0, u1, u2;
    //! \brief residuals associated with the last three estimates
    mfem::Vector r0, r1, r2;
    //! \brief workspace
    mfem::Vector n0, n1, tmp0, tmp1;
  };  // end of struct Cast3MAccelerationAlgorithm

  std::unique_ptr<AccelerationAlgorithm> makeCast3MAccelerationAlgorithm(
      Context& ctx, const Parameters& opts) noexcept {
    if (!checkParameters(ctx, opts,
                         std::vector<std::string>{"Trigger", "Period"})) {
      return {};
    }
    const auto otrigger = get_if<int>(ctx, opts, "Trigger", 4);
    const auto operiod = get_if<int>(ctx, opts, "Period", 2);
    if (isInvalid(otrigger) || isInvalid(operiod)) {
      return {};
    }
    if (*otrigger < 3) {
      return ctx.registerErrorMessage(
          "Cast3M: the trigger must be greater than or equal to 3");
    }
    if (*operiod <= 0) {
      return ctx.registerErrorMessage(
          "Cast3M: the period must be strictly positive");
    }
    return std::make_unique<Cast3MAccelerationAlgorithm>(
        static_cast<size_type>(*otrigger), static_cast<size_type>(*operiod));
  }  // end of makeCast3MAccelerationAlgorithm

}  // end of namespace mfem_mgis::internals
