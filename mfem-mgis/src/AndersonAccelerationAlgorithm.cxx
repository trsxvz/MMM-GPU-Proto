/*!
 * \file   src/AndersonAccelerationAlgorithm.cxx
 * \brief  This file implements Anderson's acceleration of the fixed point
 * iterations
 * \date   10/10/2026
 */

#include <vector>
#include <utility>
#include <string>
#include <optional>
#include <algorithm>
#include "mfem/linalg/densemat.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "AccelerationAlgorithms.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief Anderson's method (type II): the new estimate is the combination
   * of the images `g_k` of the last iterates by the fixed point map,
   * `g_i - Σ γ_k Δg_k`, whose coefficients minimize the norm of the
   * combination of the residuals `f_k`, `|| f_i - Σ γ_k Δf_k ||`, `Δ` being
   * the differences between successive iterations. The residuals are the
   * fixed point residuals `g_k - x_k` (as `UAnderson` in MTest) or the
   * residuals of the problem at `x_k` (as `FAnderson`). The least squares
   * problem is solved by its normal equations, whose matrix (scalar
   * products of the `Δf_k`) is updated incrementally: two scalar products
   * per stored iteration and per iteration. The histories are updated at
   * each iteration, the estimate is accelerated every `period` iterations
   * from the second one.
   */
  struct AndersonAccelerationAlgorithm final : AccelerationAlgorithm {
    /*!
     * \param[in] d: maximum number of differences used
     * \param[in] p: number of iterations between two accelerations
     * \param[in] f: if true, the residuals of the problem are used,
     * otherwise the fixed point residuals
     * \param[in] r: relative regularization of the normal equations
     */
    AndersonAccelerationAlgorithm(const size_type d,
                                  const size_type p,
                                  const bool f,
                                  const real r)
        : depth(d),
          period(p),
          forces(f),
          regularization(r),
          gram(static_cast<int>(d)) {
      this->df.resize(d);
      this->dg.resize(d);
    }
    void initialize(const size_type n, ScalarProduct p) override {
      this->dot = std::move(p);
      this->count = 0;
      this->next = 0;
      this->has_previous_iterate = false;
      for (auto* const v : {&(this->f_previous), &(this->g_previous)}) {
        v->UseDevice(true);
        v->SetSize(static_cast<int>(n));
      }
      for (size_type k = 0; k != this->depth; ++k) {
        for (auto* const v : {&(this->df[k]), &(this->dg[k])}) {
          v->UseDevice(true);
          v->SetSize(static_cast<int>(n));
        }
      }
    }
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector& r,
                 const real,
                 const real,
                 const size_type i) override {
      // `x` is the image of the previous iterate
      const auto& f = this->forces ? r : dx;
      if (this->has_previous_iterate) {
        // the new differences overwrite the oldest ones when the history is
        // full
        const auto k = this->next;
        subtract(f, this->f_previous, this->df[k]);
        subtract(x, this->g_previous, this->dg[k]);
        this->next = (this->next + 1) % this->depth;
        this->count = std::min(this->count + 1, this->depth);
        for (size_type l = 0; l != this->count; ++l) {
          const auto s = this->slot(l);
          const auto v = this->dot(this->df[k], this->df[s]);
          this->gram(static_cast<int>(k), static_cast<int>(s)) = v;
          this->gram(static_cast<int>(s), static_cast<int>(k)) = v;
        }
      }
      this->f_previous = f;
      this->g_previous = x;
      this->has_previous_iterate = true;
      if ((this->count == 0) || (i < 2) || ((i - 2) % this->period != 0)) {
        return;
      }
      const auto m = static_cast<int>(this->count);
      auto A = mfem::DenseMatrix(m, m);
      auto b = mfem::Vector(m);
      auto dmax = real{0};
      for (int l = 0; l != m; ++l) {
        const auto si = static_cast<int>(this->slot(l));
        for (int j = 0; j != m; ++j) {
          A(l, j) = this->gram(si, static_cast<int>(this->slot(j)));
        }
        b(l) = this->dot(this->df[si], f);
        dmax = std::max(dmax, A(l, l));
      }
      if (dmax <= 0) {
        return;
      }
      for (int l = 0; l != m; ++l) {
        A(l, l) += this->regularization * dmax;
      }
      auto gamma = mfem::Vector(m);
      auto inverse = mfem::DenseMatrixInverse(A);
      inverse.Mult(b, gamma);
      for (int l = 0; l != m; ++l) {
        x.Add(-gamma(l), this->dg[this->slot(l)]);
      }
    }

   private:
    //! \return the slot of the `l`-th stored difference, the oldest first
    [[nodiscard]] size_type slot(const size_type l) const noexcept {
      return (this->next + this->depth - this->count + l) % this->depth;
    }
    //! \brief maximum number of differences used
    const size_type depth;
    //! \brief number of iterations between two accelerations
    const size_type period;
    //! \brief if true, the residuals of the problem are used
    const bool forces;
    //! \brief relative regularization of the normal equations
    const real regularization;
    //! \brief scalar products of the differences of the fixed point residuals
    mfem::DenseMatrix gram;
    //! \brief differences of the fixed point residuals
    std::vector<mfem::Vector> df;
    //! \brief differences of the images of the iterates
    std::vector<mfem::Vector> dg;
    //! \brief previous fixed point residual
    mfem::Vector f_previous;
    //! \brief previous image
    mfem::Vector g_previous;
    //! \brief scalar product
    ScalarProduct dot;
    //! \brief number of stored differences
    size_type count = 0;
    //! \brief slot of the next difference
    size_type next = 0;
    //! \brief if true, `f_previous` and `g_previous` are set
    bool has_previous_iterate = false;
  };  // end of struct AndersonAccelerationAlgorithm

  /*!
   * \return Anderson's method
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   * \param[in] name: name of the algorithm, for the error messages
   * \param[in] default_period: default number of iterations between two
   * accelerations
   * \param[in] residuals: if not empty, residuals used (`Unknowns` or
   * `Forces`), which can't be changed by the options
   */
  static std::unique_ptr<AccelerationAlgorithm> makeAnderson(
      Context& ctx,
      const Parameters& opts,
      const std::string& name,
      const int default_period,
      const std::string& residuals) noexcept {
    auto options =
        std::vector<std::string>{"Depth", "Period", "Regularization"};
    if (residuals.empty()) {
      options.push_back("Residuals");
    }
    if (!checkParameters(ctx, opts, options)) {
      return {};
    }
    const auto odepth = get_if<int>(ctx, opts, "Depth", 4);
    const auto operiod = get_if<int>(ctx, opts, "Period", default_period);
    const auto oregularization =
        get_if<double>(ctx, opts, "Regularization", 1e-12);
    const auto oresiduals = residuals.empty()
                                ? get_if<std::string>(ctx, opts, "Residuals",
                                                      std::string{"Unknowns"})
                                : std::optional<std::string>{residuals};
    if (isInvalid(odepth) || isInvalid(operiod) || isInvalid(oregularization) ||
        isInvalid(oresiduals)) {
      return {};
    }
    if (*odepth <= 0) {
      return ctx.registerErrorMessage(name +
                                      ": the depth must be strictly positive");
    }
    if (*operiod <= 0) {
      return ctx.registerErrorMessage(name +
                                      ": the period must be strictly positive");
    }
    if (*oregularization <= 0) {
      return ctx.registerErrorMessage(
          name + ": the regularization must be strictly positive");
    }
    if ((*oresiduals != "Unknowns") && (*oresiduals != "Forces")) {
      return ctx.registerErrorMessage(name + ": invalid residuals '" +
                                      *oresiduals +
                                      "' (expected 'Unknowns' or 'Forces')");
    }
    return std::make_unique<AndersonAccelerationAlgorithm>(
        static_cast<size_type>(*odepth), static_cast<size_type>(*operiod),
        *oresiduals == "Forces", *oregularization);
  }  // end of makeAnderson

  std::unique_ptr<AccelerationAlgorithm> makeAndersonAccelerationAlgorithm(
      Context& ctx, const Parameters& opts) noexcept {
    return makeAnderson(ctx, opts, "Anderson", 1, "");
  }  // end of makeAndersonAccelerationAlgorithm

  std::unique_ptr<AccelerationAlgorithm> makeUAndersonAccelerationAlgorithm(
      Context& ctx, const Parameters& opts) noexcept {
    return makeAnderson(ctx, opts, "UAnderson", 2, "Unknowns");
  }  // end of makeUAndersonAccelerationAlgorithm

  std::unique_ptr<AccelerationAlgorithm> makeFAndersonAccelerationAlgorithm(
      Context& ctx, const Parameters& opts) noexcept {
    return makeAnderson(ctx, opts, "FAnderson", 2, "Forces");
  }  // end of makeFAndersonAccelerationAlgorithm

}  // end of namespace mfem_mgis::internals
