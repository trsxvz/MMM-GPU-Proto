/*!
 * \file   include/TFEL/Material/LogarithmicStrainHandler.ixx
 * \brief
 * \author Thomas Helfer
 * \date   25/06/2017
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_MATERIAL_LOGARITHMICSTRAINHANDLER_IXX
#define LIB_TFEL_MATERIAL_LOGARITHMICSTRAINHANDLER_IXX

#include <cmath>
#include "TFEL/Math/st2tost2.hxx"
#include "TFEL/Raise.hxx"
#include "TFEL/ContractViolation.hxx"
#include "TFEL/Math/General/Abs.hxx"
#include "TFEL/Math/General/IEEE754.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"

namespace tfel::material {

  namespace internals {

    /*!
     * \return the divided difference of the function
     * \(f\left(x\right)=\log\left(1+x\right)/2\) at \(a\) and \(b\),
     * computed without cancellation, even when \(a\) and \(b\) are close:
     *
     * \[
     * f\left[a,b\right]=\frac{\log\left(1+u\right)}{2\,u\,\left(1+b\right)}
     * \quad\text{with}\quad u=\frac{a-b}{1+b}
     * \]
     *
     * These are the divided differences of the Hencky strain with respect to
     * the eigenvalues of \(C-I\).
     * \param[in] a: first point
     * \param[in] b: second point
     */
    template <typename real>
    TFEL_HOST_DEVICE real halfLog1pDividedDifference(const real a,
                                                     const real b) noexcept {
      const auto y = 1 + b;
      const auto u = (a - b) / y;
      if (tfel::math::ieee754::fpclassify(u) == FP_ZERO) {
        return 1 / (2 * y);
      }
      return std::log1p(u) / (2 * u * y);
    }  // end of halfLog1pDividedDifference

    /*!
     * \return the second order divided difference of the function
     * \(f\left(x\right)=\log\left(1+x\right)/2\) at \(a\), \(b\) and
     * \(c\), which may be equal, computed without cancellation. With \(b\)
     * the median value, \(u=\left(a-b\right)/\left(1+b\right)\) and
     * \(v=\left(c-b\right)/\left(1+b\right)\):
     *
     * \[
     * f\left[a,b,c\right]=\frac{1}{\left(1+b\right)^{2}}
     * \sum_{k=2}^{\infty}\frac{\left(-1\right)^{k+1}}{2\,k}\,h_{k-2}\left(u,v\right)
     * \quad\text{with}\quad h_{n}\left(u,v\right)=\sum_{i=0}^{n}u^{i}\,v^{n-i}
     * \]
     *
     * This series is used when \(u\) and \(v\) are small. Otherwise, the
     * difference of the first order divided differences is divided by
     * \(a-c\), the largest gap.
     * \param[in] a: first point
     * \param[in] b: second point
     * \param[in] c: third point
     */
    template <typename real>
    TFEL_HOST_DEVICE real halfLog1pSecondDividedDifference(real a,
                                                           real b,
                                                           real c) noexcept {
      const auto swap = [](real& x, real& y) {
        const auto tmp = x;
        x = y;
        y = tmp;
      };
      // a <= b <= c
      if (a > b) {
        swap(a, b);
      }
      if (b > c) {
        swap(b, c);
      }
      if (a > b) {
        swap(a, b);
      }
      const auto y = 1 + b;
      const auto u = (a - b) / y;
      const auto v = (c - b) / y;
      if ((tfel::math::abs(u) < real{0.05}) &&
          (tfel::math::abs(v) < real{0.05})) {
        // the terms beyond k = 15 are smaller than 1e-18
        auto r = real{0};
        auto h = real{1};
        auto vn = real{1};
        for (int k = 2; k != 16; ++k) {
          r += ((k % 2 == 0) ? -h : h) / (2 * k);
          vn *= v;
          h = u * h + vn;
        }
        return r / (y * y);
      }
      return (halfLog1pDividedDifference(a, b) -
              halfLog1pDividedDifference(b, c)) /
             (a - c);
    }  // end of halfLog1pSecondDividedDifference

  }  // end of namespace internals

  TFEL_HOST_DEVICE LogarithmicStrainHandlerBase::LogarithmicStrainHandlerBase(
      const Setting c)
      : s(c) {}  // end of LogarithmicStrainHandlerBase

  TFEL_HOST_DEVICE inline void
  LogarithmicStrainHandlerBase::checkLagrangianSetting() const {
    if (this->s != LAGRANGIAN) {
#ifdef TFEL_DEVICE_CODE
      tfel::reportContractViolation(
          "LogarithmicStrainHandlerBase::checkLagrangianSetting: "
          "setting is not set to LAGRANGIAN");
#else  /* TFEL_DEVICE_CODE */
      tfel::raise(
          "LogarithmicStrainHandlerBase::checkLagrangianSetting: "
          "setting is not set to LAGRANGIAN");
#endif /* TFEL_DEVICE_CODE */
    }
  }  // end of checkLagrangianSetting

  TFEL_HOST_DEVICE inline void
  LogarithmicStrainHandlerBase::checkEulerianSetting() const {
    if (this->s != EULERIAN) {
#ifdef TFEL_DEVICE_CODE
      tfel::reportContractViolation(
          "LogarithmicStrainHandlerBase::checkEulerianSetting: "
          "setting is not set to EULERIAN");
#else  /* TFEL_DEVICE_CODE */
      tfel::raise(
          "LogarithmicStrainHandlerBase::checkEulerianSetting: "
          "setting is not set to EULERIAN");
#endif /* TFEL_DEVICE_CODE */
    }
  }  // end of checkEulerianSetting

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<1u, StressType>::LogarithmicStrainHandler(
      const Setting c, const DeformationGradient& F1) noexcept
      : LogarithmicStrainHandlerBase(c), F(F1) {
    this->e = {std::log(this->F[0]), std::log(this->F[1]),
               std::log(this->F[2])};
  }  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<1u, StressType>::LogarithmicStrainHandler(
      const Setting c,
      const DeformationGradient& F1,
      const StrainStensor& E1) noexcept
      : LogarithmicStrainHandlerBase(c), F(F1) {
    // log(F) = log(1 + 2 E) / 2, without cancellation for small strains
    this->e = {std::log1p(2 * E1[0]) / 2, std::log1p(2 * E1[1]) / 2,
               std::log1p(2 * E1[2]) / 2};
  }  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<1u, StressType>::updateAxialDeformationGradient(
      const real Fzz) {
    this->F[1] = Fzz;
  }  // end of updateAxialDeformationGradient

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::StrainStensor
      LogarithmicStrainHandler<1u, StressType>::getHenckyLogarithmicStrain()
          const {
    return {this->e[0], this->e[1], this->e[2]};
  }  // end of getHenckyLogarithmicStrain

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<1u, StressType>::getHenckyLogarithmicStrain(
      real* const elog) const {
    elog[0] = this->e[0];
    elog[1] = this->e[1];
    elog[2] = this->e[2];
  }  // end of getHenckyLogarithmicStrain

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::StressStensor
      LogarithmicStrainHandler<1u, StressType>::
          convertFromSecondPiolaKirchhoffStress(const StressStensor& S) const {
    return {S[0] * (this->F[0] * this->F[0]), S[1] * (this->F[1] * this->F[1]),
            S[2] * (this->F[2] * this->F[2])};
  }  // end of convertFromSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void LogarithmicStrainHandler<1u, StressType>::
      convertFromSecondPiolaKirchhoffStress(stress* const ST) const {
    ST[0] *= this->F[0] * this->F[0];
    ST[1] *= this->F[1] * this->F[1];
    ST[2] *= this->F[2] * this->F[2];
  }  // end of convertFromSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE typename LogarithmicStrainHandler<1u,
                                                     StressType>::StressStensor
  LogarithmicStrainHandler<1u, StressType>::convertToSecondPiolaKirchhoffStress(
      const StressStensor& T) const {
    return {T[0] / (this->F[0] * this->F[0]), T[1] / (this->F[1] * this->F[1]),
            T[2] / (this->F[2] * this->F[2])};
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<1u, StressType>::convertToSecondPiolaKirchhoffStress(
      stress* const TS) const {
    TS[0] /= this->F[0] * this->F[0];
    TS[1] /= this->F[1] * this->F[1];
    TS[2] /= this->F[2] * this->F[2];
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::StressStensor
      LogarithmicStrainHandler<1u, StressType>::convertToCauchyStress(
          const StressStensor& T) const {
    return T / tfel::math::det(this->F);
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<1u, StressType>::convertToCauchyStress(
      stress* const T) const {
    const auto iJ = 1 / tfel::math::det(this->F);
    T[0] *= iJ;
    T[1] *= iJ;
    T[2] *= iJ;
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::StressStensor
      LogarithmicStrainHandler<1u, StressType>::convertFromCauchyStress(
          const StressStensor& sig) const {
    return sig * tfel::math::det(this->F);
  }  // end of convertFromCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<1u, StressType>::convertFromCauchyStress(
      stress* const T) const {
    const auto J = tfel::math::det(this->F);
    T[0] *= J;
    T[1] *= J;
    T[2] *= J;
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::TangentOperator
      LogarithmicStrainHandler<1u, StressType>::convertToMaterialTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    const real iC[3] = {1 / ((this->F[0]) * (this->F[0])),
                        1 / ((this->F[1]) * (this->F[1])),
                        1 / ((this->F[2]) * (this->F[2]))};
    auto Kr = TangentOperator();
    Kr(0, 0) = (Ks(0, 0) - 2 * T[0]) * iC[0] * iC[0];
    Kr(0, 1) = Ks(0, 1) * iC[0] * iC[1];
    Kr(0, 2) = Ks(0, 2) * iC[0] * iC[2];
    Kr(1, 0) = Ks(1, 0) * iC[0] * iC[1];
    Kr(1, 1) = (Ks(1, 1) - 2 * T[1]) * iC[1] * iC[1];
    Kr(1, 2) = Ks(1, 2) * iC[1] * iC[2];
    Kr(2, 0) = Ks(2, 0) * iC[0] * iC[2];
    Kr(2, 1) = Ks(2, 1) * iC[1] * iC[2];
    Kr(2, 2) = (Ks(2, 2) - 2 * T[2]) * iC[2] * iC[2];
    return Kr;
  }

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::TangentOperator
      LogarithmicStrainHandler<1u, StressType>::convertToSpatialTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    auto Kr = TangentOperator();
    Kr(0, 0) = (Ks(0, 0) - 2 * T[0]);
    Kr(0, 1) = Ks(0, 1);
    Kr(0, 2) = Ks(0, 2);
    Kr(1, 0) = Ks(1, 0);
    Kr(1, 1) = (Ks(1, 1) - 2 * T[1]);
    Kr(1, 2) = Ks(1, 2);
    Kr(2, 0) = Ks(2, 0);
    Kr(2, 1) = Ks(2, 1);
    Kr(2, 2) = (Ks(2, 2) - 2 * T[2]);
    return Kr;
  }  // end of convertToSpatialTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<1u, StressType>::TangentOperator
      LogarithmicStrainHandler<1u, StressType>::
          convertToCauchyStressTruesdellRateTangentModuli(
              const TangentOperator& Ks, const StressStensor& T) const {
    const auto iJ = 1 / tfel::math::det(this->F);
    auto Kr = TangentOperator();
    Kr(0, 0) = (Ks(0, 0) - 2 * T[0]) * iJ;
    Kr(0, 1) = Ks(0, 1) * iJ;
    Kr(0, 2) = Ks(0, 2) * iJ;
    Kr(1, 0) = Ks(1, 0) * iJ;
    Kr(1, 1) = (Ks(1, 1) - 2 * T[1]) * iJ;
    Kr(1, 2) = Ks(1, 2) * iJ;
    Kr(2, 0) = Ks(2, 0) * iJ;
    Kr(2, 1) = Ks(2, 1) * iJ;
    Kr(2, 2) = (Ks(2, 2) - 2 * T[2]) * iJ;
    return Kr;
  }  // end of convertToCauchyStressTruesdellRateTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE void LogarithmicStrainHandler<1u, StressType>::
      convertToCauchyStressTruesdellRateTangentModuli(
          stress* const K, const stress* const T) const {
    const auto iJ = 1 / tfel::math::det(this->F);
    tfel::math::ST2toST2View<1u, stress> Kr(K);
    auto transpose = [&Kr]() {
      std::swap(Kr(0, 1), Kr(1, 0));
      std::swap(Kr(0, 2), Kr(2, 0));
      std::swap(Kr(1, 2), Kr(2, 1));
    };
    transpose();
    Kr(0, 0) = (Kr(0, 0) - 2 * T[0]) * iJ;
    Kr(0, 1) *= iJ;
    Kr(0, 2) *= iJ;
    Kr(1, 0) *= iJ;
    Kr(1, 1) = (Kr(1, 1) - 2 * T[1]) * iJ;
    Kr(1, 2) *= iJ;
    Kr(2, 0) *= iJ;
    Kr(2, 1) *= iJ;
    Kr(2, 2) = (Kr(2, 2) - 2 * T[2]) * iJ;
    transpose();
  }  // end of convertToCauchyStressTruesdellRateTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE const typename LogarithmicStrainHandler<1u, StressType>::
      DeformationGradient&
      LogarithmicStrainHandler<1u, StressType>::getDeformationGradient() const {
    return this->F;
  }  // end of getDeformationGradient

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<2u, StressType>::LogarithmicStrainHandler(
      const Setting c, const DeformationGradient& F1, const bool b)
      : LogarithmicStrainHandler(Builder(c, F1, b), c, F1) {
  }  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<2u, StressType>::LogarithmicStrainHandler(
      const Setting c,
      const DeformationGradient& F1,
      const StrainStensor& E1,
      const bool b)
      : LogarithmicStrainHandler(Builder(c, F1, E1, b), c, F1) {
  }  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE LogarithmicStrainHandler<2u, StressType>::Builder::Builder(
      const Setting c, const DeformationGradient& F1, const bool b) {
    using namespace tfel::math;
    using stensor = stensor<2u, real>;
    const auto fl = [](const real x) { return std::log1p(x - 1) / 2; };
    const auto C = computeRightCauchyGreenTensor(F1);
    std::tie(this->vp, this->m) =
        C.template computeEigenVectors<stensor::FSESJACOBIEIGENSOLVER>();
    this->svp = map([](const real x) { return x - 1; }, this->vp);
    if (b) {
      this->e = map(fl, this->vp);
    } else {
      this->e[0] = fl(this->vp(0));
      this->e[1] = fl(this->vp(1));
      this->e[2] = real(0);
    }
    this->computeDerivative(c, F1);
  }

  template <typename StressType>
  TFEL_HOST_DEVICE LogarithmicStrainHandler<2u, StressType>::Builder::Builder(
      const Setting c,
      const DeformationGradient& F1,
      const StrainStensor& E1,
      const bool b) {
    using namespace tfel::math;
    using stensor = stensor<2u, real>;
    // eigenvalues of C - I = 2 E and half their logarithm plus one, without
    // the cancellation of C - I for small strains
    const auto fl = [](const real x) { return std::log1p(x) / 2; };
    std::tie(this->svp, this->m) =
        stensor{2 * E1}
            .template computeEigenVectors<stensor::FSESJACOBIEIGENSOLVER>();
    this->vp = map([](const real x) { return 1 + x; }, this->svp);
    if (b) {
      this->e = map(fl, this->svp);
    } else {
      this->e[0] = fl(this->svp(0));
      this->e[1] = fl(this->svp(1));
      this->e[2] = real(0);
    }
    this->computeDerivative(c, F1);
  }

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::Builder::computeDerivative(
      const Setting c, const DeformationGradient& F1) {
    using namespace tfel::math;
    const auto dfl = [](const real x) { return 1 / (2 * x); };
    const auto d = map(dfl, this->vp);
    // p is one half of the tensor defined by Miehe. The divided difference
    // of the logarithm is computed without cancellation from the
    // eigenvalues of C - I, whose differences are of the order of the
    // rounding errors when two eigenvalues are close: no special case is
    // needed for equal eigenvalues
    const auto theta =
        internals::halfLog1pDividedDifference(this->svp[0], this->svp[1]);
    const auto N = getNTensors(m);
    const auto M = (c == LAGRANGIAN) ? N : getEulerianMTensors(m, F1);
    this->p = (d[0] * (N(0) ^ M(0)) + d[1] * (N(1) ^ M(1))) / 4;
    this->p(2, 2) = (d[2] / 4) * ((N(2))(2)) * ((M(2))(2));
    this->p += (theta / 2) * (N(3) ^ M(3));
  }

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<2u, StressType>::LogarithmicStrainHandler(
      Builder&& src, const Setting c, const DeformationGradient& F1)
      : LogarithmicStrainHandlerBase(c),
        F(F1),
        p(src.p),
        m(src.m),
        vp(src.vp),
        svp(src.svp),
        e(src.e) {}  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::StrainStensor
      LogarithmicStrainHandler<2u, StressType>::getHenckyLogarithmicStrain()
          const {
    return StrainStensor::computeIsotropicFunction(this->e, this->m);
  }  // end of getHenckyLogarithmicStrain

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::getHenckyLogarithmicStrain(
      real* const elog) const {
    constexpr auto cste = tfel::math::Cste<real>::sqrt2;
    const auto el = StrainStensor::computeIsotropicFunction(this->e, this->m);
    elog[0] = el[0];
    elog[1] = el[1];
    elog[2] = el[2];
    elog[3] = el[3] * cste;
  }  // end of getHenckyLogarithmicStrain

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::StressStensor
      LogarithmicStrainHandler<2u, StressType>::
          convertFromSecondPiolaKirchhoffStress(const StressStensor& S) const {
    this->checkLagrangianSetting();
    return (S | tfel::math::invert(p)) / 2;
  }  // end of convertFromSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void LogarithmicStrainHandler<2u, StressType>::
      convertFromSecondPiolaKirchhoffStress(stress* const ST) const {
    this->checkLagrangianSetting();
    tfel::math::stensor<2u, stress> S;
    S.importTab(ST);
    const auto T = this->convertFromSecondPiolaKirchhoffStress(S);
    T.exportTab(ST);
  }  // end of convertFromSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE typename LogarithmicStrainHandler<2u,
                                                     StressType>::StressStensor
  LogarithmicStrainHandler<2u, StressType>::convertToSecondPiolaKirchhoffStress(
      const StressStensor& T) const {
    this->checkLagrangianSetting();
    return 2 * (T | p);
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::convertToSecondPiolaKirchhoffStress(
      stress* const TS) const {
    this->checkLagrangianSetting();
    tfel::math::stensor<2u, stress> T;
    T.importTab(TS);
    const auto S = this->convertToSecondPiolaKirchhoffStress(T);
    S.exportTab(TS);
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::StressStensor
      LogarithmicStrainHandler<2u, StressType>::convertFromCauchyStress(
          const StressStensor& sig) const {
    if (this->s == LAGRANGIAN) {
      const auto S =
          tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(sig,
                                                                      this->F);
      return this->convertFromSecondPiolaKirchhoffStress(S);
    } else {
      return (sig | tfel::math::invert(p)) * tfel::math::det(this->F) / 2;
    }
  }  // end of ::convertFromCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::convertFromCauchyStress(
      stress* const sT) const {
    StressStensor sig;
    sig.importTab(sT);
    this->convertFromCauchyStress(sig).exportTab(sT);
  }  // end of convertFromCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::StressStensor
      LogarithmicStrainHandler<2u, StressType>::convertToCauchyStress(
          const StressStensor& T) const {
    if (this->s == LAGRANGIAN) {
      const auto S = this->convertToSecondPiolaKirchhoffStress(T);
      return tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(
          S, this->F);
    }
    return (2 * (T | p)) / tfel::math::det(this->F);
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::convertToCauchyStress(
      stress* const Ts) const {
    StressStensor T;
    T.importTab(Ts);
    this->convertToCauchyStress(T).exportTab(Ts);
  }  // end of convertToCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::updateAxialDeformationGradient(
      const real Fzz) {
    this->F[2] = Fzz;
  }  // end of updateAxialDeformationGradient

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::convertTangentModuli(
      TangentOperator& Kr,
      const TangentOperator& Ks,
      const StressStensor& T,
      const tfel::math::tvector<4u, tfel::math::stensor<2u, real>>& N,
      const tfel::math::tvector<4u, tfel::math::stensor<2u, real>>& M) const {
    using namespace tfel::math;
    // compute the derivative of the Hencky strain with respect to C
    const auto f = map([](const real x) { return -2 / (x * x); }, this->vp);
    // second order divided differences, without cancellation
    const tvector<2u, real> xsi = {
        internals::halfLog1pSecondDividedDifference(this->svp[0], this->svp[1],
                                                    this->svp[1]),
        internals::halfLog1pSecondDividedDifference(this->svp[1], this->svp[0],
                                                    this->svp[0])};
    const tvector<4u, real> dzeta = {(T | N(0)) / 2, (T | N(1)) / 2, T[2],
                                     (T | N(3)) / 2};
#if (!defined __INTEL_COMPILER) && (!defined __PGI)
    Kr = 4 * transpose(p) * Ks * p +
         (f[0] * dzeta(0) * (M(0) ^ M(0)) + f[1] * dzeta(1) * (M(1) ^ M(1)) +
          f[2] * dzeta(2) * (M(2) ^ M(2))) /
             4;
#else
    const st2tost2<2u, real> tp = transpose(p);
    Kr = 4 * tp * Ks * p +
         (f[0] * dzeta(0) * (M(0) ^ M(0)) + f[1] * dzeta(1) * (M(1) ^ M(1)) +
          f[2] * dzeta(2) * (M(2) ^ M(2))) /
             4;
#endif
    Kr +=
        2 * xsi[0] *
        (dzeta(3) * ((M(3) ^ M(1)) + (M(1) ^ M(3))) + dzeta(1) * (M(3) ^ M(3)));
    Kr +=
        2 * xsi[1] *
        (dzeta(3) * ((M(3) ^ M(0)) + (M(0) ^ M(3))) + dzeta(0) * (M(3) ^ M(3)));
  }  // end of convertTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE tfel::math::tvector<
      4u,
      tfel::math::
          stensor<2u, typename LogarithmicStrainHandler<2u, StressType>::real>>
  LogarithmicStrainHandler<2u, StressType>::getNTensors(
      const tfel::math::tmatrix<3u, 3u, real>& m) {
    using tvector = tfel::math::tvector<3u, real>;
    using stensor = tfel::math::stensor<2u, real>;
    const tvector v[2] = {m.template column_view<0u>(),
                          m.template column_view<1u>()};
    auto N = tfel::math::tvector<4u, stensor>();
    N(0) = stensor::buildFromVectorsSymmetricDiadicProduct(v[0], v[0]);
    N(1) = stensor::buildFromVectorsSymmetricDiadicProduct(v[1], v[1]);
    N(2) = {0, 0, 2, 0};
    N(3) = stensor::buildFromVectorsSymmetricDiadicProduct(v[0], v[1]);
    return N;
  }  // end of getNTensors

  template <typename StressType>
  TFEL_HOST_DEVICE tfel::math::tvector<
      4u,
      tfel::math::
          stensor<2u, typename LogarithmicStrainHandler<2u, StressType>::real>>
  LogarithmicStrainHandler<2u, StressType>::getEulerianMTensors(
      const tfel::math::tmatrix<3u, 3u, real>& m,
      const DeformationGradient& F) {
    using tvector = tfel::math::tvector<3u, real>;
    using stensor = tfel::math::stensor<2u, real>;
    constexpr auto zero = real{};
    const tvector v[2] = {m.template column_view<0u>(),
                          m.template column_view<1u>()};
    tfel::math::tvector<3u, real> n[2];
    n[0][0] = F[0] * v[0][0] + F[3] * v[0][1];
    n[0][1] = F[4] * v[0][0] + F[1] * v[0][1];
    n[0][2] = zero;
    n[1][0] = F[0] * v[1][0] + F[3] * v[1][1];
    n[1][1] = F[4] * v[1][0] + F[1] * v[1][1];
    n[1][2] = zero;
    auto M = tfel::math::tvector<4u, stensor>();
    M(0) = stensor::buildFromVectorsSymmetricDiadicProduct(n[0], n[0]);
    M(1) = stensor::buildFromVectorsSymmetricDiadicProduct(n[1], n[1]);
    M(2) = {zero, zero, 2 * F[2] * F[2], zero};
    M(3) = stensor::buildFromVectorsSymmetricDiadicProduct(n[0], n[1]);
    return M;
  }  // end of getEulerianMTensors

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::TangentOperator
      LogarithmicStrainHandler<2u, StressType>::convertToMaterialTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    this->checkLagrangianSetting();
    // real work starts here
    const auto N = getNTensors(this->m);
    TangentOperator Kr;
    this->convertTangentModuli(Kr, Ks, T, N, N);
    return Kr;
  }

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::TangentOperator
      LogarithmicStrainHandler<2u, StressType>::convertToSpatialTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    if (this->s == LAGRANGIAN) {
      using FSTOBase = FiniteStrainBehaviourTangentOperatorBase;
      const auto F0 = tfel::math::tensor<2u, real>::Id();
      const auto sig = this->convertToCauchyStress(T);
      const auto Cse = this->convertToMaterialTangentModuli(Ks, T);
      return convert<FSTOBase::SPATIAL_MODULI, FSTOBase::DS_DEGL>(Cse, F0,
                                                                  this->F, sig);
    } else {
      const auto N = getNTensors(this->m);
      const auto M = getEulerianMTensors(this->m, this->F);
      TangentOperator Kr;
      this->convertTangentModuli(Kr, Ks, T, N, M);
      return Kr;
    }
  }  // end of convertToSpatialTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::TangentOperator
      LogarithmicStrainHandler<2u, StressType>::
          convertToCauchyStressTruesdellRateTangentModuli(
              const TangentOperator& Ks, const StressStensor& T) const {
    if (this->s == LAGRANGIAN) {
      using FSTOBase = FiniteStrainBehaviourTangentOperatorBase;
      const auto F0 = tfel::math::tensor<2u, real>::Id();
      const auto sig = this->convertToCauchyStress(T);
      const auto Cse = this->convertToMaterialTangentModuli(Ks, T);
      auto Kr = convert<FSTOBase::SPATIAL_MODULI, FSTOBase::DS_DEGL>(
          Cse, F0, this->F, sig);
      Kr /= tfel::math::det(this->F);
      return Kr;
    }
    const auto N = getNTensors(this->m);
    const auto M = getEulerianMTensors(this->m, this->F);
    TangentOperator Kr;
    this->convertTangentModuli(Kr, Ks, T, N, M);
    Kr /= tfel::math::det(this->F);
    return Kr;
  }  // end of convertToCauchyStressTruesdellRateTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE void LogarithmicStrainHandler<2u, StressType>::
      convertToCauchyStressTruesdellRateTangentModuli(
          stress* const K, const stress* const T) const {
    tfel::math::ST2toST2View<2u, stress> k(K);
    auto transpose = [&k]() {
      std::swap(k(0, 1), k(1, 0));
      std::swap(k(0, 2), k(2, 0));
      std::swap(k(0, 3), k(3, 0));
      std::swap(k(1, 2), k(2, 1));
      std::swap(k(1, 3), k(3, 1));
      std::swap(k(2, 3), k(3, 2));
    };
    auto to_tfel = [&k, transpose] {
      constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
      k(0, 3) *= cste;
      k(1, 3) *= cste;
      k(2, 3) *= cste;
      k(3, 0) *= cste;
      k(3, 1) *= cste;
      k(3, 2) *= cste;
      k(3, 3) *= 2;
      transpose();
    };
    auto to_abaqus = [&k, &transpose] {
      constexpr const auto icste = tfel::math::Cste<real>::isqrt2;
      k(0, 3) *= icste;
      k(1, 3) *= icste;
      k(2, 3) *= icste;
      k(3, 0) *= icste;
      k(3, 1) *= icste;
      k(3, 2) *= icste;
      k(3, 3) /= 2;
      transpose();
    };
    to_tfel();
    tfel::math::stensor<2u, stress> t;
    t.importTab(T);
    k = this->convertToCauchyStressTruesdellRateTangentModuli(k, t);
    to_abaqus();
  }  // end of convertToCauchyStressTruesdellRateTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<2u, StressType>::TangentOperator
      LogarithmicStrainHandler<2u, StressType>::convertToAbaqusTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    using FSTOBase = FiniteStrainBehaviourTangentOperatorBase;
    const auto F0 = tfel::math::tensor<2u, real>::Id();
    const auto sig = this->convertToCauchyStress(T);
    if (this->s == LAGRANGIAN) {
      const auto Cse = this->convertToMaterialTangentModuli(Ks, T);
      auto Kr =
          convert<FSTOBase::ABAQUS, FSTOBase::DS_DEGL>(Cse, F0, this->F, sig);
      Kr /= tfel::math::det(this->F);
      return Kr;
    }
    const auto N = getNTensors(this->m);
    const auto M = getEulerianMTensors(this->m, this->F);
    TangentOperator Kr;
    this->convertTangentModuli(Kr, Ks, T, N, M);
    return convert<FSTOBase::ABAQUS, FSTOBase::SPATIAL_MODULI>(Kr, F0, this->F,
                                                               sig);
  }  // end of convertToAbaqusTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<2u, StressType>::convertToAbaqusTangentModuli(
      stress* const K, const stress* const T) const {
    tfel::math::ST2toST2View<2u, stress> k(K);
    auto transpose = [&k]() {
      std::swap(k(0, 1), k(1, 0));
      std::swap(k(0, 2), k(2, 0));
      std::swap(k(0, 3), k(3, 0));
      std::swap(k(1, 2), k(2, 1));
      std::swap(k(1, 3), k(3, 1));
      std::swap(k(2, 3), k(3, 2));
    };
    auto to_tfel = [&k, transpose] {
      constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
      k(0, 3) *= cste;
      k(1, 3) *= cste;
      k(2, 3) *= cste;
      k(3, 0) *= cste;
      k(3, 1) *= cste;
      k(3, 2) *= cste;
      k(3, 3) *= 2;
      transpose();
    };
    auto to_abaqus = [&k, &transpose] {
      constexpr const auto icste = tfel::math::Cste<real>::isqrt2;
      k(0, 3) *= icste;
      k(1, 3) *= icste;
      k(2, 3) *= icste;
      k(3, 0) *= icste;
      k(3, 1) *= icste;
      k(3, 2) *= icste;
      k(3, 3) /= 2;
      transpose();
    };
    to_tfel();
    tfel::math::stensor<2u, stress> t;
    t.importTab(T);
    k = this->convertToAbaqusTangentModuli(k, t);
    to_abaqus();
  }  // end of convertToAbaqusTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE const typename LogarithmicStrainHandler<2u, StressType>::
      DeformationGradient&
      LogarithmicStrainHandler<2u, StressType>::getDeformationGradient() const {
    return this->F;
  }  // end of getDeformationGradient

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<3u, StressType>::LogarithmicStrainHandler(
      const Setting c, const DeformationGradient& F1)
      : LogarithmicStrainHandler(Builder(c, F1), c, F1) {
  }  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<3u, StressType>::LogarithmicStrainHandler(
      const Setting c, const DeformationGradient& F1, const StrainStensor& E1)
      : LogarithmicStrainHandler(Builder(c, F1, E1), c, F1) {
  }  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE LogarithmicStrainHandler<3u, StressType>::Builder::Builder(
      const Setting c, const DeformationGradient& F1) {
    using namespace tfel::math;
    using stensor = stensor<3u, real>;
    const auto fl = [](const real x) { return std::log1p(x - 1) / 2; };
    const auto C = computeRightCauchyGreenTensor(F1);
    std::tie(this->vp, this->m) =
        C.template computeEigenVectors<stensor::FSESJACOBIEIGENSOLVER>();
    this->svp = map([](const real x) { return x - 1; }, this->vp);
    this->e = map(fl, this->vp);
    this->computeDerivative(c, F1);
  }

  template <typename StressType>
  TFEL_HOST_DEVICE LogarithmicStrainHandler<3u, StressType>::Builder::Builder(
      const Setting c, const DeformationGradient& F1, const StrainStensor& E1) {
    using namespace tfel::math;
    using stensor = stensor<3u, real>;
    // eigenvalues of C - I = 2 E and half their logarithm plus one, without
    // the cancellation of C - I for small strains
    std::tie(this->svp, this->m) =
        stensor{2 * E1}
            .template computeEigenVectors<stensor::FSESJACOBIEIGENSOLVER>();
    this->vp = map([](const real x) { return 1 + x; }, this->svp);
    this->e = map([](const real x) { return std::log1p(x) / 2; }, this->svp);
    this->computeDerivative(c, F1);
  }

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::Builder::computeDerivative(
      const Setting c, const DeformationGradient& F1) {
    using namespace tfel::math;
    const auto dfl = [](const real x) { return 1 / (2 * x); };
    const auto d = map(dfl, this->vp);
    // p is one half of the tensor defined by Miehe. The divided differences
    // of the logarithm are computed without cancellation from the
    // eigenvalues of C - I, whose differences are of the order of the
    // rounding errors when two eigenvalues are close: no special case is
    // needed for equal eigenvalues
    const auto N = getNTensors(m);
    const auto M = (c == LAGRANGIAN) ? N : getEulerianMTensors(m, F1);
    this->p = (d[0] * (N(0, 0) ^ M(0, 0)) + d[1] * (N(1, 1) ^ M(1, 1)) +
               d[2] * (N(2, 2) ^ M(2, 2))) /
              4;
    for (size_type i = 0; i != 3; ++i) {
      for (size_type j = 0; j != 3; ++j) {
        if (i == j) {
          continue;
        }
        const auto theta =
            internals::halfLog1pDividedDifference(this->svp[i], this->svp[j]);
        this->p += (theta / 4) * (N(i, j) ^ M(i, j));
      }
    }
  }

  template <typename StressType>
  TFEL_HOST_DEVICE
  LogarithmicStrainHandler<3u, StressType>::LogarithmicStrainHandler(
      Builder&& src, const Setting c, const DeformationGradient& F1)
      : LogarithmicStrainHandlerBase(c),
        F(F1),
        p(src.p),
        m(src.m),
        vp(src.vp),
        svp(src.svp),
        e(src.e) {}  // end of LogarithmicStrainHandler

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::StrainStensor
      LogarithmicStrainHandler<3u, StressType>::getHenckyLogarithmicStrain()
          const {
    return StrainStensor::computeIsotropicFunction(this->e, this->m);
  }  // end of getHenckyLogarithmicStrain

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::getHenckyLogarithmicStrain(
      real* const elog) const {
    constexpr auto cste = tfel::math::Cste<real>::sqrt2;
    const auto el = StrainStensor::computeIsotropicFunction(this->e, this->m);
    elog[0] = el[0];
    elog[1] = el[1];
    elog[2] = el[2];
    elog[3] = el[3] * cste;
    elog[4] = el[4] * cste;
    elog[5] = el[5] * cste;
  }  // end of getHenckyLogarithmicStrain

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::StressStensor
      LogarithmicStrainHandler<3u, StressType>::
          convertFromSecondPiolaKirchhoffStress(const StressStensor& S) const {
    this->checkLagrangianSetting();
    return (S | tfel::math::invert(p)) / 2;
  }  // end of convertFromSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void LogarithmicStrainHandler<3u, StressType>::
      convertFromSecondPiolaKirchhoffStress(stress* const ST) const {
    this->checkLagrangianSetting();
    tfel::math::stensor<3u, stress> S;
    S.importTab(ST);
    const auto T = this->convertFromSecondPiolaKirchhoffStress(S);
    T.exportTab(ST);
  }  // end of convertFromSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE typename LogarithmicStrainHandler<3u,
                                                     StressType>::StressStensor
  LogarithmicStrainHandler<3u, StressType>::convertToSecondPiolaKirchhoffStress(
      const StressStensor& T) const {
    this->checkLagrangianSetting();
    return 2 * (T | p);
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::convertToSecondPiolaKirchhoffStress(
      stress* const TS) const {
    this->checkLagrangianSetting();
    tfel::math::stensor<3u, stress> T;
    T.importTab(TS);
    const auto S = this->convertToSecondPiolaKirchhoffStress(T);
    S.exportTab(TS);
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::StressStensor
      LogarithmicStrainHandler<3u, StressType>::convertToCauchyStress(
          const StressStensor& T) const {
    if (this->s == LAGRANGIAN) {
      const auto S = this->convertToSecondPiolaKirchhoffStress(T);
      return tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(
          S, this->F);
    } else {
      return (2 * (T | p)) / tfel::math::det(this->F);
    }
  }  // end of convertToSecondPiolaKirchhoffStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::convertToCauchyStress(
      stress* const Ts) const {
    StressStensor T;
    T.importTab(Ts);
    this->convertToCauchyStress(T).exportTab(Ts);
  }  // end of convertToCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::StressStensor
      LogarithmicStrainHandler<3u, StressType>::convertFromCauchyStress(
          const StressStensor& sig) const {
    if (this->s == LAGRANGIAN) {
      const auto S =
          tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(sig,
                                                                      this->F);
      return this->convertFromSecondPiolaKirchhoffStress(S);
    } else {
      return (sig | tfel::math::invert(p)) * tfel::math::det(this->F) / 2;
    }
  }  // end of convertFromCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::convertFromCauchyStress(
      stress* const sT) const {
    StressStensor sig;
    sig.importTab(sT);
    this->convertFromCauchyStress(sig).exportTab(sT);
  }  // end of convertFromCauchyStress

  template <typename StressType>
  TFEL_HOST_DEVICE tfel::math::tmatrix<
      3u,
      3u,
      tfel::math::
          stensor<3u, typename LogarithmicStrainHandler<3u, StressType>::real>>
  LogarithmicStrainHandler<3u, StressType>::getNTensors(
      const tfel::math::tmatrix<3u, 3u, real>& m) {
    using namespace tfel::math;
    using tvector = tfel::math::tvector<3u, real>;
    using stensor = tfel::math::stensor<3u, real>;
    const tvector v[3] = {m.template column_view<0u>(),
                          m.template column_view<1u>(),
                          m.template column_view<2u>()};
    auto N = tmatrix<3u, 3u, stensor>{};
    for (size_type i = 0; i != 3; ++i) {
      for (size_type j = 0; j != 3; ++j) {
        N(i, j) = stensor::buildFromVectorsSymmetricDiadicProduct(v[i], v[j]);
      }
    }
    return N;
  }  // end of getNTensors

  template <typename StressType>
  TFEL_HOST_DEVICE tfel::math::tmatrix<
      3u,
      3u,
      tfel::math::
          stensor<3u, typename LogarithmicStrainHandler<3u, StressType>::real>>
  LogarithmicStrainHandler<3u, StressType>::getEulerianMTensors(
      const tfel::math::tmatrix<3u, 3u, real>& m,
      const DeformationGradient& F) {
    using namespace tfel::math;
    using tvector = tfel::math::tvector<3u, real>;
    using stensor = tfel::math::stensor<3u, real>;
    auto apply = [&F](const tvector& v2) -> tvector {
      return {F[0] * v2[0] + F[3] * v2[1] + F[5] * v2[2],
              F[4] * v2[0] + F[1] * v2[1] + F[7] * v2[2],
              F[6] * v2[0] + F[8] * v2[1] + F[2] * v2[2]};
    };
    const tvector v[3] = {m.template column_view<0u>(),
                          m.template column_view<1u>(),
                          m.template column_view<2u>()};
    const tvector n[3] = {apply(v[0]), apply(v[1]), apply(v[2])};
    auto M = tmatrix<3u, 3u, stensor>{};
    for (size_type i = 0; i != 3; ++i) {
      for (size_type j = 0; j != 3; ++j) {
        M(i, j) = stensor::buildFromVectorsSymmetricDiadicProduct(n[i], n[j]);
      }
    }
    return M;
  }  // end of getEulerianMTensors

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::convertTangentModuli(
      TangentOperator& Kr,
      const TangentOperator& Ks,
      const StressStensor& T,
      const tfel::math::tmatrix<3u, 3u, tfel::math::stensor<3u, real>>& N,
      const tfel::math::tmatrix<3u, 3u, tfel::math::stensor<3u, real>>& M)
      const {
    using tfel::math::tmatrix;
    const auto lk = [](const size_type i, const size_type j) -> size_type {
      if (i == 0) {
        return (j == 1) ? 2 : 1;
      }
      if (i == 1) {
        return (j == 0) ? 2 : 0;
      }
      return (j == 0) ? 1 : 0;
    };
    const auto f = map([](const real x) { return -2 / (x * x); }, this->vp);
    // second order divided differences, without cancellation
    const auto xsi = [this] {
      auto r = tmatrix<3u, 3u, real>{};
      for (size_type i = 0; i != 3; ++i) {
        for (size_type j = 0; j != 3; ++j) {
          r(i, j) = (i == j) ? real{}
                             : internals::halfLog1pSecondDividedDifference(
                                   this->svp[i], this->svp[j], this->svp[j]);
        }
      }
      return r;
    }();
    const auto eta = internals::halfLog1pSecondDividedDifference(
        this->svp[0], this->svp[1], this->svp[2]);
    const auto dzeta = [&N, &T] {
      auto r = tmatrix<3u, 3u, real>{};
      for (size_type i = 0; i != 3; ++i) {
        for (size_type j = 0; j != 3; ++j) {
          r(i, j) = (T | N(i, j)) / 2;
        }
      }
      return r;
    }();
#if (!defined __INTEL_COMPILER) && (!defined __PGI)
    Kr =
        4 * transpose(p) * Ks * p + (f[0] * dzeta(0, 0) * (M(0, 0) ^ M(0, 0)) +
                                     f[1] * dzeta(1, 1) * (M(1, 1) ^ M(1, 1)) +
                                     f[2] * dzeta(2, 2) * (M(2, 2) ^ M(2, 2))) /
                                        4;
#else
    const tfel::math::st2tost2<3u, real> tp = transpose(p);
    Kr = 4 * tp * Ks * p + (f[0] * dzeta(0, 0) * (M(0, 0) ^ M(0, 0)) +
                            f[1] * dzeta(1, 1) * (M(1, 1) ^ M(1, 1)) +
                            f[2] * dzeta(2, 2) * (M(2, 2) ^ M(2, 2))) /
                               4;
#endif
    for (size_type i = 0; i != 3; ++i) {
      for (size_type j = 0; j != 3; ++j) {
        if (i == j) {
          continue;
        }
        const auto k = lk(i, j);
        Kr += 2 * eta * dzeta(i, j) * (M(i, k) ^ M(j, k)) +
              2 * xsi(i, j) *
                  (dzeta(i, j) * ((M(i, j) ^ M(j, j)) + (M(j, j) ^ M(i, j))) +
                   dzeta(j, j) * (M(i, j) ^ M(i, j)));
      }
    }
  }

  template <typename StressType>
  typename LogarithmicStrainHandler<3u, StressType>::TangentOperator
      TFEL_HOST_DEVICE
      LogarithmicStrainHandler<3u, StressType>::convertToMaterialTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    this->checkLagrangianSetting();
    // real work starts here
    const auto N = getNTensors(this->m);
    TangentOperator Kr;
    this->convertTangentModuli(Kr, Ks, T, N, N);
    return Kr;
  }

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::TangentOperator
      LogarithmicStrainHandler<3u, StressType>::convertToSpatialTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    if (this->s == LAGRANGIAN) {
      using FSTOBase = FiniteStrainBehaviourTangentOperatorBase;
      const auto F0 = tfel::math::tensor<3u, real>::Id();
      const auto sig = this->convertToCauchyStress(T);
      const auto Cse = this->convertToMaterialTangentModuli(Ks, T);
      return convert<FSTOBase::SPATIAL_MODULI, FSTOBase::DS_DEGL>(Cse, F0,
                                                                  this->F, sig);
    } else {
      const auto N = getNTensors(this->m);
      const auto M = getEulerianMTensors(this->m, this->F);
      TangentOperator Kr;
      this->convertTangentModuli(Kr, Ks, T, N, M);
      return Kr;
    }
  }  // end of convertToSpatialTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::TangentOperator
      LogarithmicStrainHandler<3u, StressType>::
          convertToCauchyStressTruesdellRateTangentModuli(
              const TangentOperator& Ks, const StressStensor& T) const {
    if (this->s == LAGRANGIAN) {
      using FSTOBase = FiniteStrainBehaviourTangentOperatorBase;
      const auto F0 = tfel::math::tensor<3u, real>::Id();
      const auto sig = this->convertToCauchyStress(T);
      const auto Cse = this->convertToMaterialTangentModuli(Ks, T);
      auto Kr = convert<FSTOBase::SPATIAL_MODULI, FSTOBase::DS_DEGL>(
          Cse, F0, this->F, sig);
      Kr /= tfel::math::det(this->F);
      return Kr;
    } else {
      const auto N = getNTensors(this->m);
      const auto M = getEulerianMTensors(this->m, this->F);
      TangentOperator Kr;
      this->convertTangentModuli(Kr, Ks, T, N, M);
      Kr /= tfel::math::det(this->F);
      return Kr;
    }
  }  // end of convertToCauchyStressTruesdellRateTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE void LogarithmicStrainHandler<3u, StressType>::
      convertToCauchyStressTruesdellRateTangentModuli(
          stress* const K, const stress* const T) const {
    tfel::math::ST2toST2View<3u, stress> k(K);
    auto transpose = [&k]() {
      std::swap(k(0, 1), k(1, 0));
      std::swap(k(0, 2), k(2, 0));
      std::swap(k(0, 3), k(3, 0));
      std::swap(k(0, 4), k(4, 0));
      std::swap(k(0, 5), k(5, 0));
      std::swap(k(1, 2), k(2, 1));
      std::swap(k(1, 3), k(3, 1));
      std::swap(k(1, 4), k(4, 1));
      std::swap(k(1, 5), k(5, 1));
      std::swap(k(2, 3), k(3, 2));
      std::swap(k(2, 4), k(4, 2));
      std::swap(k(2, 5), k(5, 2));
      std::swap(k(3, 4), k(4, 3));
      std::swap(k(3, 5), k(5, 3));
      std::swap(k(4, 5), k(5, 4));
    };
    auto to_tfel = [&k, &transpose] {
      constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
      k(0, 3) *= cste;
      k(1, 3) *= cste;
      k(2, 3) *= cste;
      k(0, 4) *= cste;
      k(1, 4) *= cste;
      k(2, 4) *= cste;
      k(0, 5) *= cste;
      k(1, 5) *= cste;
      k(2, 5) *= cste;
      k(3, 0) *= cste;
      k(3, 1) *= cste;
      k(3, 2) *= cste;
      k(4, 0) *= cste;
      k(4, 1) *= cste;
      k(4, 2) *= cste;
      k(5, 0) *= cste;
      k(5, 1) *= cste;
      k(5, 2) *= cste;
      k(3, 3) *= 2;
      k(3, 4) *= 2;
      k(3, 5) *= 2;
      k(4, 3) *= 2;
      k(4, 4) *= 2;
      k(4, 5) *= 2;
      k(5, 3) *= 2;
      k(5, 4) *= 2;
      k(5, 5) *= 2;
      transpose();
    };
    auto to_abaqus = [&k, &transpose] {
      constexpr const auto icste = tfel::math::Cste<real>::isqrt2;
      k(0, 3) *= icste;
      k(1, 3) *= icste;
      k(2, 3) *= icste;
      k(0, 4) *= icste;
      k(1, 4) *= icste;
      k(2, 4) *= icste;
      k(0, 5) *= icste;
      k(1, 5) *= icste;
      k(2, 5) *= icste;
      k(3, 0) *= icste;
      k(3, 1) *= icste;
      k(3, 2) *= icste;
      k(4, 0) *= icste;
      k(4, 1) *= icste;
      k(4, 2) *= icste;
      k(5, 0) *= icste;
      k(5, 1) *= icste;
      k(5, 2) *= icste;
      k(3, 3) /= 2;
      k(3, 4) /= 2;
      k(3, 5) /= 2;
      k(4, 3) /= 2;
      k(4, 4) /= 2;
      k(4, 5) /= 2;
      k(5, 3) /= 2;
      k(5, 4) /= 2;
      k(5, 5) /= 2;
      transpose();
    };
    to_tfel();
    tfel::math::stensor<3u, stress> t;
    t.importTab(T);
    k = this->convertToCauchyStressTruesdellRateTangentModuli(k, t);
    to_abaqus();
  }  // end of convertToCauchyStressTruesdellRateTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE
      typename LogarithmicStrainHandler<3u, StressType>::TangentOperator
      LogarithmicStrainHandler<3u, StressType>::convertToAbaqusTangentModuli(
          const TangentOperator& Ks, const StressStensor& T) const {
    using FSTOBase = FiniteStrainBehaviourTangentOperatorBase;
    const auto F0 = tfel::math::tensor<3u, real>::Id();
    const auto sig = this->convertToCauchyStress(T);
    if (this->s == LAGRANGIAN) {
      const auto Cse = this->convertToMaterialTangentModuli(Ks, T);
      auto Kr =
          convert<FSTOBase::ABAQUS, FSTOBase::DS_DEGL>(Cse, F0, this->F, sig);
      Kr /= tfel::math::det(this->F);
      return Kr;
    } else {
      const auto N = getNTensors(this->m);
      const auto M = getEulerianMTensors(this->m, this->F);
      TangentOperator Kr;
      this->convertTangentModuli(Kr, Ks, T, N, M);
      return convert<FSTOBase::ABAQUS, FSTOBase::SPATIAL_MODULI>(Kr, F0,
                                                                 this->F, sig);
    }
  }  // end of convertToAbaqusTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE void
  LogarithmicStrainHandler<3u, StressType>::convertToAbaqusTangentModuli(
      stress* const K, const stress* const T) const {
    tfel::math::ST2toST2View<3u, stress> k(K);
    auto transpose = [&k]() {
      std::swap(k(0, 1), k(1, 0));
      std::swap(k(0, 2), k(2, 0));
      std::swap(k(0, 3), k(3, 0));
      std::swap(k(0, 4), k(4, 0));
      std::swap(k(0, 5), k(5, 0));
      std::swap(k(1, 2), k(2, 1));
      std::swap(k(1, 3), k(3, 1));
      std::swap(k(1, 4), k(4, 1));
      std::swap(k(1, 5), k(5, 1));
      std::swap(k(2, 3), k(3, 2));
      std::swap(k(2, 4), k(4, 2));
      std::swap(k(2, 5), k(5, 2));
      std::swap(k(3, 4), k(4, 3));
      std::swap(k(3, 5), k(5, 3));
      std::swap(k(4, 5), k(5, 4));
    };
    auto to_tfel = [&k, &transpose] {
      constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
      k(0, 3) *= cste;
      k(1, 3) *= cste;
      k(2, 3) *= cste;
      k(0, 4) *= cste;
      k(1, 4) *= cste;
      k(2, 4) *= cste;
      k(0, 5) *= cste;
      k(1, 5) *= cste;
      k(2, 5) *= cste;
      k(3, 0) *= cste;
      k(3, 1) *= cste;
      k(3, 2) *= cste;
      k(4, 0) *= cste;
      k(4, 1) *= cste;
      k(4, 2) *= cste;
      k(5, 0) *= cste;
      k(5, 1) *= cste;
      k(5, 2) *= cste;
      k(3, 3) *= 2;
      k(3, 4) *= 2;
      k(3, 5) *= 2;
      k(4, 3) *= 2;
      k(4, 4) *= 2;
      k(4, 5) *= 2;
      k(5, 3) *= 2;
      k(5, 4) *= 2;
      k(5, 5) *= 2;
      transpose();
    };
    auto to_abaqus = [&k, &transpose] {
      constexpr const auto icste = tfel::math::Cste<real>::isqrt2;
      k(0, 3) *= icste;
      k(1, 3) *= icste;
      k(2, 3) *= icste;
      k(0, 4) *= icste;
      k(1, 4) *= icste;
      k(2, 4) *= icste;
      k(0, 5) *= icste;
      k(1, 5) *= icste;
      k(2, 5) *= icste;
      k(3, 0) *= icste;
      k(3, 1) *= icste;
      k(3, 2) *= icste;
      k(4, 0) *= icste;
      k(4, 1) *= icste;
      k(4, 2) *= icste;
      k(5, 0) *= icste;
      k(5, 1) *= icste;
      k(5, 2) *= icste;
      k(3, 3) /= 2;
      k(3, 4) /= 2;
      k(3, 5) /= 2;
      k(4, 3) /= 2;
      k(4, 4) /= 2;
      k(4, 5) /= 2;
      k(5, 3) /= 2;
      k(5, 4) /= 2;
      k(5, 5) /= 2;
      transpose();
    };
    to_tfel();
    tfel::math::stensor<3u, stress> t;
    t.importTab(T);
    k = this->convertToAbaqusTangentModuli(k, t);
    to_abaqus();
  }  // end of convertToAbaqusTangentModuli

  template <typename StressType>
  TFEL_HOST_DEVICE const typename LogarithmicStrainHandler<3u, StressType>::
      DeformationGradient&
      LogarithmicStrainHandler<3u, StressType>::getDeformationGradient() const {
    return this->F;
  }  // end of getDeformationGradient

}  // end of namespace tfel::material

#endif /* LIB_TFEL_MATERIAL_LOGARITHMICSTRAINHANDLER_IXX */
