/*!
 * \file   DisplacementGradientTest.cxx
 * \brief  This file tests the strain measures computed from the displacement
 * gradient \(H = F - I\): the
 * `computeGreenLagrangeTensorFromDisplacementGradient` function and the
 * `LogarithmicStrainHandler` class built from the deformation gradient and the
 * Green-Lagrange strain.
 *
 * In double, the results are the ones obtained from the deformation
 * gradient, within the cancellation of the computation from the deformation
 * gradient (about 1e-16 divided by the norm of the displacement gradient).
 * In single precision, the strains computed from the displacement gradient
 * keep its relative precision, whereas the ones computed from the
 * deformation gradient only know it within about 6e-8 in absolute value.
 * The stresses and the tangent operators are also checked in single
 * precision when two or three principal values of the displacement
 * gradient are equal, against a reference computed in long double, and
 * the material tangent operator, for two equal stretches and an
 * anisotropic tangent operator, against a finite difference.
 * \date   14/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifdef NDEBUG
#undef NDEBUG
#endif /* NDEBUG */

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/st2tost2.hxx"
#include "TFEL/Material/LogarithmicStrainHandler.hxx"

struct DisplacementGradientTest final : public tfel::tests::TestCase {
  DisplacementGradientTest()
      : tfel::tests::TestCase("TFEL/Material", "DisplacementGradientTest") {
  }  // end of DisplacementGradientTest
  tfel::tests::TestResult execute() override {
    for (const auto a : {1e-5, 1e-3, 0.1}) {
      this->checkDouble<1u>(a);
      this->checkDouble<2u>(a);
      this->checkDouble<3u>(a);
    }
    this->checkFloat<1u>();
    this->checkFloat<2u>();
    this->checkFloat<3u>();
    this->checkFloatEqualPrincipalValues<2u>();
    this->checkFloatEqualPrincipalValues<3u>();
    this->checkTangentOperatorWithEqualStretches();
    return this->result;
  }

 private:
  //! \return a displacement gradient whose components are of order `a`
  template <unsigned short N, typename real>
  static tfel::math::tensor<N, real> getDisplacementGradient(const double a) {
    constexpr double v[9] = {1, -0.3, 0.45, 0.2, -0.15, 0.1, 0.35, -0.25, 0.05};
    auto H = tfel::math::tensor<N, real>{};
    for (unsigned short i = 0; i != H.size(); ++i) {
      H[i] = static_cast<real>(a * v[i]);
    }
    return H;
  }
  //! \return the relative gap between two objects
  template <typename T1, typename T2>
  static double gap(const T1& a, const T2& b) {
    auto m = 0.;
    auto n = 0.;
    for (unsigned short i = 0; i != a.size(); ++i) {
      const auto ai = static_cast<double>(*(a.begin() + i));
      const auto bi = static_cast<double>(*(b.begin() + i));
      m = std::max(m, std::abs(ai - bi));
      n = std::max(n, std::abs(bi));
    }
    return n == 0 ? m : m / n;
  }
  /*!
   * \brief in double, the strains, the stresses and the tangent operators
   * obtained from the displacement gradient are the ones obtained from the
   * deformation gradient
   */
  template <unsigned short N>
  void checkDouble(const double a) {
    using namespace tfel::math;
    using LSH = tfel::material::LogarithmicStrainHandler<N, double>;
    const auto H = getDisplacementGradient<N, double>(a);
    const auto F = eval(tensor<N, double>::Id() + H);
    // the cancellation of the computation from F
    const auto eps = 1e-14 / a;
    const auto E = computeGreenLagrangeTensorFromDisplacementGradient(H);
    TFEL_TESTS_ASSERT(gap(E, computeGreenLagrangeTensor(F)) < eps);
    auto T = stensor<N, double>{};
    for (unsigned short i = 0; i != T.size(); ++i) {
      T[i] = 1e8 * (0.3 + 0.1 * i);
    }
    auto K = st2tost2<N, double>{};
    for (unsigned short i = 0; i != K.size(); ++i) {
      *(K.begin() + i) = 1e11 * std::cos(static_cast<double>(i));
    }
    for (const auto s : {LSH::LAGRANGIAN, LSH::EULERIAN}) {
      const auto lF = LSH(s, F);
      const auto lH = LSH(s, F, E);
      TFEL_TESTS_ASSERT(gap(lH.getHenckyLogarithmicStrain(),
                            lF.getHenckyLogarithmicStrain()) < eps);
      TFEL_TESTS_ASSERT(gap(lH.convertToCauchyStress(T),
                            lF.convertToCauchyStress(T)) < 1e-12);
      if (s == LSH::LAGRANGIAN) {
        TFEL_TESTS_ASSERT(gap(lH.convertToMaterialTangentModuli(K, T),
                              lF.convertToMaterialTangentModuli(K, T)) < eps);
      } else {
        TFEL_TESTS_ASSERT(gap(lH.convertToSpatialTangentModuli(K, T),
                              lF.convertToSpatialTangentModuli(K, T)) < eps);
      }
    }
  }
  /*!
   * \brief in single precision, the strains computed from the displacement
   * gradient keep its relative precision for small strains
   */
  template <unsigned short N>
  void checkFloat() {
    using namespace tfel::math;
    using LSHd = tfel::material::LogarithmicStrainHandler<N, double>;
    using LSHf = tfel::material::LogarithmicStrainHandler<N, float>;
    constexpr auto a = 1e-5;
    const auto Hd = getDisplacementGradient<N, double>(a);
    const auto Hf = getDisplacementGradient<N, float>(a);
    const auto Fd = eval(tensor<N, double>::Id() + Hd);
    const auto Ff = eval(tensor<N, float>::Id() + Hf);
    // references, computed in double from the displacement gradient
    const auto Ed = computeGreenLagrangeTensorFromDisplacementGradient(Hd);
    const auto ld = LSHd(LSHd::LAGRANGIAN, Fd, Ed).getHenckyLogarithmicStrain();
    const auto Ef = computeGreenLagrangeTensorFromDisplacementGradient(Hf);
    TFEL_TESTS_ASSERT(gap(Ef, Ed) < 1e-6);
    TFEL_TESTS_ASSERT(
        gap(LSHf(LSHf::LAGRANGIAN, Ff, Ef).getHenckyLogarithmicStrain(), ld) <
        1e-6);
    // from the deformation gradient, the strains are far less precise
    TFEL_TESTS_ASSERT(gap(computeGreenLagrangeTensor(Ff), Ed) > 1e-4);
  }
  /*!
   * \return the symmetric displacement gradient whose principal values are
   * `d`, rotated by `theta` around the third axis and, in 3D, by
   * `theta / 2` around the first axis
   */
  template <unsigned short N>
  static tfel::math::tensor<N, long double> getRotatedDisplacementGradient(
      const std::array<long double, 3>& d, const long double theta) {
    const auto c = std::cos(theta), s = std::sin(theta);
    const auto c2 = (N == 3) ? std::cos(theta / 2) : 1.L;
    const auto s2 = (N == 3) ? std::sin(theta / 2) : 0.L;
    // R = Rz(theta) Rx(theta / 2)
    const long double R[3][3] = {{c, -s * c2, s * s2},  //
                                 {s, c * c2, -c * s2},  //
                                 {0, s2, c2}};
    long double H[3][3];
    for (int i = 0; i != 3; ++i) {
      for (int j = 0; j != 3; ++j) {
        H[i][j] = 0;
        for (int k = 0; k != 3; ++k) {
          H[i][j] += R[i][k] * d[k] * R[j][k];
        }
      }
    }
    auto t = tfel::math::tensor<N, long double>{};
    t[0] = H[0][0];
    t[1] = H[1][1];
    t[2] = H[2][2];
    t[3] = H[0][1];
    t[4] = H[1][0];
    if constexpr (N == 3) {
      t[5] = H[0][2];
      t[6] = H[2][0];
      t[7] = H[1][2];
      t[8] = H[2][1];
    }
    return t;
  }
  /*!
   * \brief stresses and tangent operators in single precision when two or
   * three principal values of the displacement gradient are equal (their
   * float values differ by rounding errors), compared with a reference
   * computed in long double from the same displacement gradient: the
   * divided differences of the logarithm are computed without cancellation
   */
  template <unsigned short N>
  void checkFloatEqualPrincipalValues() {
    for (const long double a : {1e-5L, 1e-3L, 0.1L}) {
      // two equal principal values (in the plane in 2D) and three
      const auto cases = {
          std::array<long double, 3>{a, (N == 2) ? a : -0.3L * a, -0.3L * a},
          std::array<long double, 3>{a, a, a}};
      for (const auto& d : cases) {
        for (const auto theta : {0.L, 0.3L, 1.1L, 2.5L}) {
          this->checkFloatStressesAndTangentOperators<N>(
              getRotatedDisplacementGradient<N>(d, theta));
        }
      }
    }
  }
  //! \brief see `checkFloatEqualPrincipalValues`
  template <unsigned short N>
  void checkFloatStressesAndTangentOperators(
      const tfel::math::tensor<N, long double>& Hl) {
    using namespace tfel::math;
    using LSHl = tfel::material::LogarithmicStrainHandler<N, long double>;
    using LSHf = tfel::material::LogarithmicStrainHandler<N, float>;
    auto Hf = tensor<N, float>{};
    for (unsigned short i = 0; i != Hf.size(); ++i) {
      Hf[i] = static_cast<float>(Hl[i]);
    }
    // the reference is computed from the rounded displacement gradient
    auto Hr = tensor<N, long double>{};
    for (unsigned short i = 0; i != Hf.size(); ++i) {
      Hr[i] = Hf[i];
    }
    const auto Fl = eval(tensor<N, long double>::Id() + Hr);
    const auto El = computeGreenLagrangeTensorFromDisplacementGradient(Hr);
    const auto Ff = eval(tensor<N, float>::Id() + Hf);
    const auto Ef = computeGreenLagrangeTensorFromDisplacementGradient(Hf);
    // dual stress of an isotropic elastic behaviour, plus a part which is
    // not coaxial with the strain (otherwise, the divided differences of
    // the logarithm do not contribute to the stresses), and tangent
    // operator of an isotropic elastic behaviour
    constexpr auto lambda = 1.2e11L;
    constexpr auto mu = 8e10L;
    constexpr long double v[6] = {0.3, -0.2, 0.1, 0.4, 0.25, -0.15};
    const auto e = LSHl(LSHl::LAGRANGIAN, Fl, El).getHenckyLogarithmicStrain();
    auto Tl =
        eval(lambda * trace(e) * stensor<N, long double>::Id() + 2 * mu * e);
    const auto ne = std::sqrt(e | e);
    for (unsigned short i = 0; i != Tl.size(); ++i) {
      Tl[i] += 2 * mu * ne * v[i];
    }
    const auto Kl = eval(lambda * st2tost2<N, long double>::IxI() +
                         2 * mu * st2tost2<N, long double>::Id());
    auto Tf = stensor<N, float>{};
    for (unsigned short i = 0; i != Tf.size(); ++i) {
      Tf[i] = static_cast<float>(Tl[i]);
    }
    auto Kf = st2tost2<N, float>{};
    for (unsigned short i = 0; i != Kf.size(); ++i) {
      *(Kf.begin() + i) = static_cast<float>(*(Kl.begin() + i));
    }
    // a few times the precision of float
    constexpr auto eps = 2e-5;
    const auto lL = LSHl(LSHl::LAGRANGIAN, Fl, El);
    const auto lE = LSHl(LSHl::EULERIAN, Fl, El);
    // handlers built from the deformation gradient and the Green-Lagrange
    // strain, and from the deformation gradient only: the conversions only
    // depend on the eigenvalues of C, known with the precision of float in
    // both cases
    for (const auto from_strain : {true, false}) {
      const auto fL = from_strain ? LSHf(LSHf::LAGRANGIAN, Ff, Ef)
                                  : LSHf(LSHf::LAGRANGIAN, Ff);
      TFEL_TESTS_ASSERT(gap(fL.convertToSecondPiolaKirchhoffStress(Tf),
                            lL.convertToSecondPiolaKirchhoffStress(Tl)) < eps);
      TFEL_TESTS_ASSERT(gap(fL.convertToMaterialTangentModuli(Kf, Tf),
                            lL.convertToMaterialTangentModuli(Kl, Tl)) < eps);
      const auto fE =
          from_strain ? LSHf(LSHf::EULERIAN, Ff, Ef) : LSHf(LSHf::EULERIAN, Ff);
      TFEL_TESTS_ASSERT(gap(fE.convertToCauchyStress(Tf),
                            lE.convertToCauchyStress(Tl)) < eps);
      TFEL_TESTS_ASSERT(gap(fE.convertToSpatialTangentModuli(Kf, Tf),
                            lE.convertToSpatialTangentModuli(Kl, Tl)) < eps);
    }
  }
  /*!
   * \brief material tangent operator of the handler built from the
   * deformation gradient, for two equal stretches, an anisotropic tangent
   * operator and a stress not coaxial with the strain, compared with a
   * centered finite difference of the second Piola-Kirchhoff stress with
   * respect to the Green-Lagrange strain, in long double
   */
  void checkTangentOperatorWithEqualStretches() {
    using namespace tfel::math;
    using LSH = tfel::material::LogarithmicStrainHandler<3u, long double>;
    auto K = st2tost2<3u, long double>{};
    for (unsigned short i = 0; i != 6; ++i) {
      for (unsigned short j = 0; j != 6; ++j) {
        K(i, j) = (i == j ? 2e11L : 0.L) +
                  1e10L * std::cos(static_cast<long double>(i + 2 * j));
      }
    }
    const auto T0 =
        stensor<3u, long double>{3e8L, -1e8L, 2e8L, 1.5e8L, -0.5e8L, 0.7e8L};
    // dual stress of the strain, not coaxial with it
    const auto dual = [&K, &T0](const LSH& l) {
      return eval(T0 + K * l.getHenckyLogarithmicStrain());
    };
    // handler associated with the right stretch tensor U = (I + 2 E)^(1/2)
    const auto handler = [](const stensor<3u, long double>& E) {
      const auto C = eval(stensor<3u, long double>::Id() + 2 * E);
      const auto U = C.computeIsotropicFunction(
          [](const long double x) { return std::sqrt(x); });
      constexpr auto icste = Cste<long double>::isqrt2;
      auto F = tensor<3u, long double>{};
      F[0] = U[0];
      F[1] = U[1];
      F[2] = U[2];
      F[3] = F[4] = U[3] * icste;
      F[5] = F[6] = U[4] * icste;
      F[7] = F[8] = U[5] * icste;
      return LSH(LSH::LAGRANGIAN, F);
    };
    // two equal stretches
    auto F = tensor<3u, long double>::Id();
    F[0] = 1.05L;
    F[1] = 0.98L;
    F[2] = 0.98L;
    const auto l = LSH(LSH::LAGRANGIAN, F);
    const auto Kr = l.convertToMaterialTangentModuli(K, dual(l));
    const auto E0 = computeGreenLagrangeTensor(F);
    constexpr auto h = 1e-7L;
    auto Kn = st2tost2<3u, long double>{};
    for (unsigned short j = 0; j != 6; ++j) {
      auto Ep = E0;
      auto Em = E0;
      Ep[j] += h;
      Em[j] -= h;
      const auto lp = handler(Ep);
      const auto lm = handler(Em);
      const auto dS = eval((lp.convertToSecondPiolaKirchhoffStress(dual(lp)) -
                            lm.convertToSecondPiolaKirchhoffStress(dual(lm))) /
                           (2 * h));
      for (unsigned short i = 0; i != 6; ++i) {
        Kn(i, j) = dS[i];
      }
    }
    // 9e-8 measured, the precision of the finite difference; 1.4e-4 with
    // the former special case of two equal eigenvalues
    TFEL_TESTS_ASSERT(gap(Kr, Kn) < 1e-6);
  }
};

TFEL_TESTS_GENERATE_PROXY(DisplacementGradientTest, "DisplacementGradientTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("DisplacementGradientTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
