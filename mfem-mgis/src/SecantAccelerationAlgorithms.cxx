/*!
 * \file   src/SecantAccelerationAlgorithms.cxx
 * \brief  This file implements the acceleration algorithms of the fixed
 * point iterations of MTest based on secant and delta squared methods
 * (`Secant`, `CrossedSecant`, `AlternateSecant`, `IronsTuck`,
 * `Steffensen`, `CrossedDelta2`, `Crossed2Delta`, `Crossed2Deltabis`,
 * `AlternateDelta2` and `Alternate2Delta`)
 * \date   10/10/2026
 */

#include <cmath>
#include <array>
#include <limits>
#include <string>
#include <vector>
#include <utility>
#include "MGIS/Contract.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "AccelerationAlgorithms.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief base class of the algorithms of this file, ported from MTest
   * (`mtest::*AccelerationAlgorithm`). With the notations of MTest, at
   * iteration `n`, `G(X_n)` is the new estimate given to `execute` and
   * `Delta X_n = G(X_n) - X_n` the fixed point residual (the correction).
   * \tparam N: number of vectors kept by the algorithm
   */
  template <std::size_t N>
  struct HistoryAccelerationAlgorithmBase : AccelerationAlgorithm {
    //! \param[in] t: first accelerated iteration
    explicit HistoryAccelerationAlgorithmBase(const size_type t) : trigger(t) {}
    void initialize(const size_type n, ScalarProduct p) override {
      this->dot = std::move(p);
      for (auto& v : this->vectors) {
        v.UseDevice(true);
        v.SetSize(static_cast<int>(n));
        v = real{0};
      }
    }

   protected:
    /*!
     * \return the threshold below which a scalar quantity of the scale of
     * `eps` is considered null, as in MTest
     * \param[in] eps: scale
     */
    [[nodiscard]] static real getThreshold(const real eps) noexcept {
      return 100 * eps * std::numeric_limits<real>::epsilon();
    }
    //! \brief first accelerated iteration
    const size_type trigger;
    //! \brief scalar product
    ScalarProduct dot;
    //! \brief vectors kept by the algorithm
    std::array<mfem::Vector, N> vectors;
  };

  /*!
   * \brief secant method on the residuals (forces):
   * `G(X_n) - a (G(X_n) - G(X_{n-1}))`, `a` minimizing the norm of the
   * corresponding combination of the residuals
   */
  struct SecantAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<5> {
    using HistoryAccelerationAlgorithmBase<5>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector&,
                 const mfem::Vector& r,
                 const real r_eps,
                 const real,
                 const size_type i) override {
      auto& [u0, u1, r0, r1, dr] = this->vectors;
      std::swap(u0, u1);
      std::swap(r0, r1);
      r1 = r;
      u1 = x;
      subtract(r1, r0, dr);
      if (i >= this->trigger) {
        const auto nr2_dr = this->dot(dr, dr);
        if (nr2_dr > getThreshold(r_eps)) {
          const auto a = this->dot(r1, dr) / nr2_dr;
          x.Add(-a, u1);
          x.Add(a, u0);
        }
      }
    }
  };

  /*!
   * \brief crossed secant method: `G(X_n) - a Delta X_n` with
   * `a = (dG . dDelta) / |dDelta|²`
   */
  struct CrossedSecantAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<6> {
    using HistoryAccelerationAlgorithmBase<6>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, r0, r1, du, dr] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(r0, r1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du);
      subtract(r1, r0, dr);
      if (i >= this->trigger) {
        const auto nr2_dr = this->dot(dr, dr);
        if (nr2_dr > eps * eps) {
          x.Add(-this->dot(du, dr) / nr2_dr, r1);
        }
      }
    }
  };

  /*!
   * \brief alternate secant method: `G(X_n) - a dG` with
   * `a = (dDelta . Delta X_n) / |dDelta|²`
   */
  struct AlternateSecantAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<6> {
    using HistoryAccelerationAlgorithmBase<6>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, r0, r1, du, dr] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(r0, r1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du);
      subtract(r1, r0, dr);
      if (i >= this->trigger) {
        const auto nr2_dr = this->dot(dr, dr);
        if (nr2_dr > eps * eps) {
          x.Add(-this->dot(dr, r1) / nr2_dr, du);
        }
      }
    }
  };

  /*!
   * \brief method of Irons and Tuck, applied every two iterations:
   * `G(X_n) - a Delta X_n` with `a = (Delta X_n . dDelta) / |dDelta|²`
   */
  struct IronsTuckAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<3> {
    using HistoryAccelerationAlgorithmBase<3>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [r0, r1, dr] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(r0, r1);
      r1 = dx;
      if ((i >= this->trigger) && ((i - this->trigger) % 2 == 0)) {
        subtract(r1, r0, dr);
        const auto nr2_dr = this->dot(dr, dr);
        if (nr2_dr > eps * eps) {
          x.Add(-this->dot(r1, dr) / nr2_dr, r1);
        }
      }
    }
  };

  /*!
   * \brief method of Steffensen, component by component, applied every two
   * iterations. The components are combined on the host.
   */
  struct SteffensenAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<3> {
    using HistoryAccelerationAlgorithmBase<3>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector&,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, u2] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(u1, u2);
      u2 = x;
      if ((i < this->trigger) || ((i - this->trigger) % 2 != 0)) {
        return;
      }
      const auto* const pu0 = u0.HostRead();
      const auto* const pu1 = u1.HostRead();
      const auto* const pu2 = u2.HostRead();
      auto* const px = x.HostReadWrite();
      for (int id = 0; id != x.Size(); ++id) {
        const auto du2 = pu2[id] - pu1[id];
        const auto du1 = pu1[id] - pu0[id];
        if ((std::abs(du2) > eps) && (std::abs(du1) > eps)) {
          const auto i1 = 1 / du2;
          const auto i2 = 1 / du1;
          if (std::abs(i1 - i2) > eps) {
            px[id] = pu1[id] + 1 / (i1 - i2);
          }
        }
      }
    }
  };

  /*!
   * \brief crossed delta squared method: crossed secant at the second
   * iteration, then `G(X_n) - a dDelta_n` with `a = (dG . d²Delta) /
   * |d²Delta|²`
   */
  struct CrossedDelta2AccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<8> {
    using HistoryAccelerationAlgorithmBase<8>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, r0, r1, dr0, dr1, du, d2r] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(r0, r1);
      std::swap(dr0, dr1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du);
      subtract(r1, r0, dr1);
      subtract(dr1, dr0, d2r);
      if (i < this->trigger) {
        return;
      }
      if (i == 2) {
        const auto nr2_dr = this->dot(dr1, dr1);
        if (nr2_dr > eps * eps) {
          x.Add(-this->dot(du, dr1) / nr2_dr, r1);
        }
      } else {
        const auto nr2_d2r = this->dot(d2r, d2r);
        if (nr2_d2r > eps * eps) {
          x.Add(-this->dot(du, d2r) / nr2_d2r, dr1);
        }
      }
    }
  };

  /*!
   * \brief crossed method of order 2: crossed secant at the second
   * iteration, then `G(X_n) - l1 Delta X_n - l2 Delta X_{n-1}`, `l1` and
   * `l2` minimizing `|dG - l1 dDelta_n - l2 dDelta_{n-1}|`, unless the last
   * differences are almost colinear
   */
  struct Crossed2DeltaAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<7> {
    using HistoryAccelerationAlgorithmBase<7>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, r0, r1, dr0, dr1, du] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(r0, r1);
      std::swap(dr0, dr1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du);
      subtract(r1, r0, dr1);
      if (i < this->trigger) {
        return;
      }
      const auto nr2_dr1 = this->dot(dr1, dr1);
      if (i != 2) {
        const auto nr2_dr0 = this->dot(dr0, dr0);
        const auto dr1_dr0 = this->dot(dr1, dr0);
        const auto det = nr2_dr1 * nr2_dr0 - dr1_dr0 * dr1_dr0;
        if ((dr1_dr0 * dr1_dr0) / (nr2_dr1 * nr2_dr0) < 0.99) {
          const auto b1 = this->dot(dr1, du);
          const auto b2 = this->dot(dr0, du);
          const auto l1 = (nr2_dr0 * b1 - dr1_dr0 * b2) / det;
          const auto l2 = (nr2_dr1 * b2 - dr1_dr0 * b1) / det;
          x.Add(-l1, r1);
          x.Add(-l2, r0);
          return;
        }
      }
      if (nr2_dr1 > eps * eps) {
        x.Add(-this->dot(du, dr1) / nr2_dr1, r1);
      }
    }
  };

  /*!
   * \brief variant of the crossed method of order 2 using the differences
   * of the accelerated estimates (`Crossed2Deltabis` of MTest)
   */
  struct Crossed2DeltabisAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<12> {
    using HistoryAccelerationAlgorithmBase<
        12>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, x0, x1, x2, dx0, dx1, dx2, r0, r1, du, dr] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(x0, x1);
      std::swap(x1, x2);
      std::swap(dx0, dx1);
      std::swap(dx1, dx2);
      std::swap(r0, r1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du);
      subtract(r1, r0, dr);
      if (i >= this->trigger) {
        const auto ma = this->dot(dr, dr);
        auto crossed_secant = [&] {
          if (ma > eps * eps) {
            x.Add(-this->dot(du, dr) / ma, r1);
          }
        };
        if (i == 2) {
          crossed_secant();
        } else {
          // `dudx` is stored in `dx2`, set at the end of this method
          auto& dudx = dx2;
          subtract(du, dx0, dudx);
          const auto mc = this->dot(dudx, dudx);
          const auto mb = this->dot(dr, dudx);
          const auto det = ma * mc - mb * mb;
          if ((mb * mb) / (ma * mc) < 0.99) {
            const auto b1 = this->dot(dr, du);
            const auto b2 = this->dot(dudx, du);
            const auto l1 = (mc * b1 - mb * b2) / det;
            const auto l2 = (ma * b2 - mb * b1) / det;
            // u1 - x0 is computed in `dudx`, no longer used
            subtract(u1, x0, dudx);
            x.Add(-l1, r1);
            x.Add(-l2, dudx);
          } else {
            crossed_secant();
          }
        }
      }
      x2 = x;
      subtract(x2, x1, dx2);
    }
  };

  /*!
   * \brief alternate delta squared method: alternate secant at the second
   * iteration, then `G(X_n) - a d²G` with `a = (d²Delta . Delta X_n) /
   * |d²Delta|²`
   */
  struct AlternateDelta2AccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<10> {
    using HistoryAccelerationAlgorithmBase<
        10>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, du0, du1, dr0, dr1, r0, r1, d2u, d2r] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(du0, du1);
      std::swap(dr0, dr1);
      std::swap(r0, r1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du1);
      subtract(du1, du0, d2u);
      subtract(r1, r0, dr1);
      subtract(dr1, dr0, d2r);
      if (i < this->trigger) {
        return;
      }
      if (i == 2) {
        const auto nr2_dr = this->dot(dr1, dr1);
        if (nr2_dr > eps * eps) {
          x.Add(-this->dot(dr1, r1) / nr2_dr, du1);
        }
      } else {
        const auto nr2_d2r = this->dot(d2r, d2r);
        if (nr2_d2r > eps * eps) {
          x.Add(-this->dot(d2r, r1) / nr2_d2r, d2u);
        }
      }
    }
  };

  /*!
   * \brief alternate method of order 2: alternate secant at the second
   * iteration, then `G(X_n) - l1 dG_n - l2 dG_{n-1}`, `l1` and `l2`
   * minimizing `|Delta X_n - l1 dDelta_n - l2 dDelta_{n-1}|`, unless the
   * last differences are almost colinear
   */
  struct Alternate2DeltaAccelerationAlgorithm final
      : HistoryAccelerationAlgorithmBase<8> {
    using HistoryAccelerationAlgorithmBase<8>::HistoryAccelerationAlgorithmBase;
    void execute(mfem::Vector& x,
                 const mfem::Vector& dx,
                 const mfem::Vector&,
                 const real,
                 const real u_eps,
                 const size_type i) override {
      auto& [u0, u1, r0, r1, du0, du1, dr0, dr1] = this->vectors;
      const auto eps = getThreshold(u_eps);
      std::swap(u0, u1);
      std::swap(r0, r1);
      std::swap(du0, du1);
      std::swap(dr0, dr1);
      r1 = dx;
      u1 = x;
      subtract(u1, u0, du1);
      subtract(r1, r0, dr1);
      if (i < this->trigger) {
        return;
      }
      const auto nr2_dr1 = this->dot(dr1, dr1);
      if (i != 2) {
        const auto nr2_dr0 = this->dot(dr0, dr0);
        const auto dr1_dr0 = this->dot(dr1, dr0);
        const auto det = nr2_dr1 * nr2_dr0 - dr1_dr0 * dr1_dr0;
        if ((dr1_dr0 * dr1_dr0) / (nr2_dr1 * nr2_dr0) < 0.99) {
          const auto b1 = this->dot(dr1, r1);
          const auto b2 = this->dot(dr0, r1);
          const auto l1 = (nr2_dr0 * b1 - dr1_dr0 * b2) / det;
          const auto l2 = (nr2_dr1 * b2 - dr1_dr0 * b1) / det;
          x.Add(-l1, du1);
          x.Add(-l2, du0);
          return;
        }
      }
      if (nr2_dr1 > eps * eps) {
        x.Add(-this->dot(dr1, r1) / nr2_dr1, du1);
      }
    }
  };

  /*!
   * \return an algorithm of this file, whose only option is `Trigger`, the
   * first accelerated iteration
   * \tparam Algorithm: algorithm
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   * \param[in] name: name of the algorithm
   * \param[in] default_trigger: default value of the trigger
   * \param[in] minimum_trigger: minimum value of the trigger
   */
  template <typename Algorithm>
  static std::unique_ptr<AccelerationAlgorithm> makeTriggeredAlgorithm(
      Context& ctx,
      const Parameters& opts,
      const char* const name,
      const int default_trigger,
      const int minimum_trigger) noexcept {
    if (!checkParameters(ctx, opts, std::vector<std::string>{"Trigger"})) {
      return {};
    }
    const auto otrigger = get_if<int>(ctx, opts, "Trigger", default_trigger);
    if (isInvalid(otrigger)) {
      return {};
    }
    if (*otrigger < minimum_trigger) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the trigger must be greater than or equal to " +
          std::to_string(minimum_trigger));
    }
    return std::make_unique<Algorithm>(static_cast<size_type>(*otrigger));
  }  // end of makeTriggeredAlgorithm

  void declareSecantAccelerationAlgorithms(
      AccelerationAlgorithmFactory& f) noexcept {
    auto ctx = Context{};
    auto check = [&ctx](const bool b) {
      if (!b) {
        auto eh = mgis::ContractViolationHandler{};
        eh.registerErrorMessage(ctx.getErrorMessage().c_str());
      }
    };
    // name, default trigger and minimum trigger as in MTest
    check(f.add(ctx, "Secant", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<SecantAccelerationAlgorithm>(c, o, "Secant",
                                                                 3, 3);
    }));
    check(f.add(ctx, "CrossedSecant", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<CrossedSecantAccelerationAlgorithm>(
          c, o, "CrossedSecant", 2, 2);
    }));
    check(f.add(ctx, "AlternateSecant", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<AlternateSecantAccelerationAlgorithm>(
          c, o, "AlternateSecant", 2, 2);
    }));
    check(f.add(ctx, "IronsTuck", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<IronsTuckAccelerationAlgorithm>(
          c, o, "IronsTuck", 2, 2);
    }));
    check(f.add(ctx, "Steffensen", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<SteffensenAccelerationAlgorithm>(
          c, o, "Steffensen", 3, 3);
    }));
    check(f.add(ctx, "CrossedDelta2", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<CrossedDelta2AccelerationAlgorithm>(
          c, o, "CrossedDelta2", 3, 2);
    }));
    check(f.add(ctx, "Crossed2Delta", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<Crossed2DeltaAccelerationAlgorithm>(
          c, o, "Crossed2Delta", 2, 2);
    }));
    check(f.add(ctx, "Crossed2Deltabis", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<Crossed2DeltabisAccelerationAlgorithm>(
          c, o, "Crossed2Deltabis", 2, 2);
    }));
    check(f.add(ctx, "AlternateDelta2", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<AlternateDelta2AccelerationAlgorithm>(
          c, o, "AlternateDelta2", 3, 2);
    }));
    check(f.add(ctx, "Alternate2Delta", [](Context& c, const Parameters& o) {
      return makeTriggeredAlgorithm<Alternate2DeltaAccelerationAlgorithm>(
          c, o, "Alternate2Delta", 2, 2);
    }));
  }  // end of declareSecantAccelerationAlgorithms

}  // end of namespace mfem_mgis::internals
