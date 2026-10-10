/*!
 * \file   src/AndersonAccelerationAlgorithm.cxx
 * \brief  This file implements Anderson's acceleration of the fixed point
 * iterations
 * \date   10/10/2026
 */

#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include "mfem/linalg/densemat.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "AccelerationAlgorithms.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief Anderson's method (type II): the new estimate is the combination
   * of the images `g_k` of the last iterates by the fixed point map,
   * `g_i - Σ γ_k Δg_k`, whose coefficients minimize the norm of the
   * combination of the fixed point residuals `f_k = g_k - x_k`,
   * `|| f_i - Σ γ_k Δf_k ||`, `Δ` being the differences between successive
   * iterations. The least squares problem is solved by its normal
   * equations, whose matrix (scalar products of the `Δf_k`) is updated
   * incrementally: two scalar products per stored iteration and per
   * iteration.
   */
  struct AndersonAccelerationAlgorithm final : AccelerationAlgorithm {
    /*!
     * \param[in] d: maximum number of differences used
     * \param[in] r: relative regularization of the normal equations
     */
    AndersonAccelerationAlgorithm(const size_type d, const real r)
        : depth(d), regularization(r), gram(static_cast<int>(d)) {
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
                 const mfem::Vector&,
                 const real,
                 const size_type) override {
      // `x` is the image of the previous iterate and `dx` the fixed point
      // residual
      if (this->has_previous_iterate) {
        // the new differences overwrite the oldest ones when the history is
        // full
        const auto k = this->next;
        subtract(dx, this->f_previous, this->df[k]);
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
      this->f_previous = dx;
      this->g_previous = x;
      this->has_previous_iterate = true;
      if (this->count == 0) {
        return;
      }
      const auto m = static_cast<int>(this->count);
      auto A = mfem::DenseMatrix(m, m);
      auto b = mfem::Vector(m);
      auto dmax = real{0};
      for (int i = 0; i != m; ++i) {
        const auto si = static_cast<int>(this->slot(i));
        for (int j = 0; j != m; ++j) {
          A(i, j) = this->gram(si, static_cast<int>(this->slot(j)));
        }
        b(i) = this->dot(this->df[si], dx);
        dmax = std::max(dmax, A(i, i));
      }
      if (dmax <= 0) {
        return;
      }
      for (int i = 0; i != m; ++i) {
        A(i, i) += this->regularization * dmax;
      }
      auto gamma = mfem::Vector(m);
      auto inverse = mfem::DenseMatrixInverse(A);
      inverse.Mult(b, gamma);
      for (int i = 0; i != m; ++i) {
        x.Add(-gamma(i), this->dg[this->slot(i)]);
      }
    }

   private:
    //! \return the slot of the `l`-th stored difference, the oldest first
    [[nodiscard]] size_type slot(const size_type l) const noexcept {
      return (this->next + this->depth - this->count + l) % this->depth;
    }
    //! \brief maximum number of differences used
    const size_type depth;
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

  std::unique_ptr<AccelerationAlgorithm> makeAndersonAccelerationAlgorithm(
      Context& ctx, const Parameters& opts) noexcept {
    if (!checkParameters(ctx, opts,
                         std::vector<std::string>{"Depth", "Regularization"})) {
      return {};
    }
    const auto odepth = get_if<int>(ctx, opts, "Depth", 4);
    const auto oregularization =
        get_if<double>(ctx, opts, "Regularization", 1e-12);
    if (isInvalid(odepth) || isInvalid(oregularization)) {
      return {};
    }
    if (*odepth <= 0) {
      return ctx.registerErrorMessage(
          "Anderson: the depth must be strictly positive");
    }
    if (*oregularization < 0) {
      return ctx.registerErrorMessage(
          "Anderson: the regularization must be positive");
    }
    return std::make_unique<AndersonAccelerationAlgorithm>(
        static_cast<size_type>(*odepth), *oregularization);
  }  // end of makeAndersonAccelerationAlgorithm

}  // end of namespace mfem_mgis::internals
