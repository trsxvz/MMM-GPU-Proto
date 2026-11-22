/*!
 * \file  tests/Material/FiniteStrainBehaviourTangentOperator17.cxx
 * \brief conversions of the finite strain tangent operators: every
 * conversion available is compared to the tangent operators of a
 * neo-Hookean behaviour, in simple and double precision
 * \date 22/11/2026
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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/st2tost2.hxx"
#include "TFEL/Math/t2tost2.hxx"
#include "TFEL/Math/t2tot2.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"

struct FiniteStrainBehaviourTangentOperator17 final
    : public tfel::tests::TestCase {
  FiniteStrainBehaviourTangentOperator17()
      : tfel::tests::TestCase("TFEL/Material",
                              "FiniteStrainBehaviourTangentOperator17") {
  }  // end of FiniteStrainBehaviourTangentOperator17
  tfel::tests::TestResult execute() override {
    this->check<1u, double>(1.e-12);
    this->check<2u, double>(1.e-12);
    this->check<3u, double>(1.e-12);
    this->check<1u, float>(2.e-4f);
    this->check<2u, float>(2.e-4f);
    this->check<3u, float>(2.e-4f);
    return this->result;
  }  // end of execute

 private:
  using TangentOperator =
      tfel::material::FiniteStrainBehaviourTangentOperatorBase;
  //! \brief first Lamé coefficient
  static constexpr double l0 = 1.09465e+11;
  //! \brief shear modulus
  static constexpr double m0 = 5.6391e+10;
  //! \brief type of a tangent operator
  template <TangentOperator::Flag f, unsigned short N, typename real>
  using TangentOperatorType = typename tfel::material::
      FiniteStrainBehaviourTangentOperatorType<f, N, real>::type;
  //! \return the Cauchy stress of a neo-Hookean behaviour
  template <unsigned short N, typename real>
  static tfel::math::stensor<N, real> getCauchyStress(
      const tfel::math::tensor<N, real>& F) {
    const auto J = tfel::math::det(F);
    const auto B = tfel::math::computeLeftCauchyGreenTensor(F);
    const auto Id = tfel::math::stensor<N, real>::Id();
    return (real(l0) * std::log(J) * Id + real(m0) * (B - Id)) / J;
  }
  /*!
   * \return a tangent operator of a neo-Hookean behaviour, deduced from its
   * spatial moduli by a fixed chain of elementary conversions
   * \param[in] F0: deformation gradient at the beginning of the time step
   * \param[in] F1: deformation gradient at the end of the time step
   */
  template <TangentOperator::Flag f, unsigned short N, typename real>
  static TangentOperatorType<f, N, real> getTangentOperator(
      const tfel::math::tensor<N, real>& F0,
      const tfel::math::tensor<N, real>& F1) {
    using namespace tfel::material;
    using st2tost2 = tfel::math::st2tost2<N, real>;
    // operator from which the operator requested is deduced
    constexpr auto from = [] {
      switch (f) {
        case TangentOperator::DTAU_DDF:
        case TangentOperator::DSIG_DF:
          return TangentOperator::DTAU_DF;
        case TangentOperator::DSIG_DDF:
        case TangentOperator::DPK1_DF:
          return TangentOperator::DSIG_DF;
        case TangentOperator::DS_DC:
        case TangentOperator::DS_DF:
          return TangentOperator::DS_DEGL;
        case TangentOperator::C_TAU_JAUMANN:
          return TangentOperator::ABAQUS;
        default:
          break;
      }
      return TangentOperator::SPATIAL_MODULI;
    }();
    if constexpr (f == TangentOperator::SPATIAL_MODULI) {
      const auto J = tfel::math::det(F1);
      const auto m = real(m0) - real(l0) * std::log(J);
      return real(l0) * st2tost2::IxI() + 2 * m * st2tost2::Id();
    } else {
      return convert<f, from>(getTangentOperator<from, N, real>(F0, F1), F0, F1,
                              getCauchyStress(F1));
    }
  }
  /*!
   * \brief compare a tangent operator to its expected value
   * \param[in] K: tangent operator
   * \param[in] Kr: expected value
   * \param[in] from: tangent operator converted
   * \param[in] to: result of the conversion
   * \param[in] eps: criterion value, relative to the shear modulus
   */
  template <typename TangentOperatorType2, typename real>
  void compare(const TangentOperatorType2& K,
               const TangentOperatorType2& Kr,
               const TangentOperator::Flag from,
               const TangentOperator::Flag to,
               const real eps) {
    using namespace tfel::material;
    auto p = Kr.begin();
    for (const auto& v : K) {
      const auto e = std::abs(v - *p);
      if (!(e < eps * real(m0))) {
        std::cout << convertFiniteStrainBehaviourTangentOperatorFlagToString(
                         from)
                  << " -> "
                  << convertFiniteStrainBehaviourTangentOperatorFlagToString(to)
                  << " ("
                  << tfel::math::getSpaceDimension<TangentOperatorType2>()
                  << "D): " << v << " " << *p << " " << e << std::endl;
      }
      TFEL_TESTS_ASSERT(e < eps * real(m0));
      ++p;
    }
  }  // end of compare
  /*!
   * \brief check the conversion of a tangent operator
   * \param[in] F0: deformation gradient at the beginning of the time step
   * \param[in] F1: deformation gradient at the end of the time step
   * \param[in] eps: criterion value, relative to the shear modulus
   */
  template <TangentOperator::Flag from,
            TangentOperator::Flag to,
            unsigned short N,
            typename real>
  void checkConversion(const tfel::math::tensor<N, real>& F0,
                       const tfel::math::tensor<N, real>& F1,
                       const real eps) {
    using namespace tfel::material;
    const auto K = convert<to, from>(getTangentOperator<from, N, real>(F0, F1),
                                     F0, F1, getCauchyStress(F1));
    this->compare(K, getTangentOperator<to, N, real>(F0, F1), from, to, eps);
  }  // end of checkConversion
  /*!
   * \brief check that the conversion of the derivative of the dual stress
   * of the Hencky strain with respect to the Hencky strain (DT_DELOG) to
   * DPK1_DF is the composition of its conversion to DS_DEGL and of the
   * conversion of DS_DEGL to DPK1_DF
   * \param[in] F0: deformation gradient at the beginning of the time step
   * \param[in] F1: deformation gradient at the end of the time step
   * \param[in] eps: criterion value, relative to the shear modulus
   */
  template <unsigned short N, typename real>
  void checkLogarithmicStrainConversion(const tfel::math::tensor<N, real>& F0,
                                        const tfel::math::tensor<N, real>& F1,
                                        const real eps) {
    using namespace tfel::material;
    using st2tost2 = tfel::math::st2tost2<N, real>;
    const auto s = getCauchyStress(F1);
    const st2tost2 Ks =
        real(l0) * st2tost2::IxI() + 2 * real(m0) * st2tost2::Id();
    const auto K = convert<TangentOperator::DPK1_DF, TangentOperator::DT_DELOG>(
        Ks, F0, F1, s);
    const auto Kr = convert<TangentOperator::DPK1_DF, TangentOperator::DS_DEGL>(
        convert<TangentOperator::DS_DEGL, TangentOperator::DT_DELOG>(Ks, F0, F1,
                                                                     s),
        F0, F1, s);
    this->compare(K, Kr, TangentOperator::DT_DELOG, TangentOperator::DPK1_DF,
                  eps);
  }  // end of checkLogarithmicStrainConversion
  /*!
   * \brief check the conversions for various deformation gradients
   * \param[in] eps: criterion value, relative to the shear modulus
   */
  template <unsigned short N, typename real>
  void check(const real eps) {
    using tensor = tfel::math::tensor<N, real>;
    // small strain, moderate strain, rotation of 0.7 rad around the third
    // axis followed by a stretch, large strain
    const real v0[9u] = {1.01, 0.995,  1.02,  0.012, -0.007,
                         0.02, -0.011, 0.004, -0.009};
    const real v1[9u] = {1.03, 0.98,   1.09,   0.03, -0.012,
                         0.04, -0.028, -0.015, 0.005};
    const real v2[9u] = {0.841326, 0.726599, 0.95, -0.708640, 0.612024,
                         0.,       0.,       0.,   0.};
    const real v3[9u] = {1.35, 0.75, 1.2, 0.22, -0.15, 0.3, -0.18, 0.12, 0.25};
    const auto F0s = {tensor::Id(), tensor::Id(), tensor{v0},
                      tensor{v0},   tensor{v1},   tensor{v2}};
    const auto F1s = {tensor::Id(), tensor{v1}, tensor{v1},
                      tensor{v2},   tensor{v3}, tensor{v3}};
    auto pF0 = F0s.begin();
    for (const auto& F1 : F1s) {
      const auto& F0 = *pF0;
      ++pF0;
      this->checkConversion<TangentOperator::DSIG_DF, TangentOperator::DSIG_DDF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DF, TangentOperator::DTAU_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DF, TangentOperator::DTAU_DDF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DF, TangentOperator::DS_DEGL,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DF,
                            TangentOperator::C_TRUESDELL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DF, TangentOperator::ABAQUS,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DF, TangentOperator::DPK1_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DSIG_DDF, TangentOperator::DSIG_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DF, TangentOperator::DSIG_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DF, TangentOperator::DTAU_DDF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DF, TangentOperator::DS_DEGL,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DF,
                            TangentOperator::SPATIAL_MODULI, N, real>(F0, F1,
                                                                      eps);
      this->checkConversion<TangentOperator::DTAU_DF,
                            TangentOperator::C_TRUESDELL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DF,
                            TangentOperator::C_TAU_JAUMANN, N, real>(F0, F1,
                                                                     eps);
      this->checkConversion<TangentOperator::DTAU_DF, TangentOperator::ABAQUS,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DF, TangentOperator::DPK1_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DDF, TangentOperator::DSIG_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DDF, TangentOperator::DTAU_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DDF, TangentOperator::DS_DEGL,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DDF,
                            TangentOperator::C_TRUESDELL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DDF, TangentOperator::ABAQUS,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DTAU_DDF, TangentOperator::DPK1_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DF, TangentOperator::DTAU_DF, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC, TangentOperator::DSIG_DF, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC, TangentOperator::DTAU_DDF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC, TangentOperator::DS_DF, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC, TangentOperator::DS_DEGL, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC,
                            TangentOperator::C_TRUESDELL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC, TangentOperator::ABAQUS, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DC, TangentOperator::DPK1_DF, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL, TangentOperator::DSIG_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL, TangentOperator::DTAU_DDF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL, TangentOperator::DS_DF, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL, TangentOperator::DS_DC, N,
                            real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL,
                            TangentOperator::SPATIAL_MODULI, N, real>(F0, F1,
                                                                      eps);
      this->checkConversion<TangentOperator::DS_DEGL,
                            TangentOperator::C_TRUESDELL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL, TangentOperator::ABAQUS,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DS_DEGL, TangentOperator::DPK1_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::DSIG_DF, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::DTAU_DF, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::DTAU_DDF, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::DS_DEGL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::C_TRUESDELL, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::C_TAU_JAUMANN, N, real>(F0, F1,
                                                                     eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::ABAQUS, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::SPATIAL_MODULI,
                            TangentOperator::DPK1_DF, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::C_TRUESDELL,
                            TangentOperator::DSIG_DF, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::C_TRUESDELL,
                            TangentOperator::SPATIAL_MODULI, N, real>(F0, F1,
                                                                      eps);
      this->checkConversion<TangentOperator::C_TAU_JAUMANN,
                            TangentOperator::DTAU_DF, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::C_TAU_JAUMANN,
                            TangentOperator::SPATIAL_MODULI, N, real>(F0, F1,
                                                                      eps);
      this->checkConversion<TangentOperator::C_TAU_JAUMANN,
                            TangentOperator::ABAQUS, N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::ABAQUS, TangentOperator::DSIG_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::ABAQUS, TangentOperator::DTAU_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::ABAQUS,
                            TangentOperator::SPATIAL_MODULI, N, real>(F0, F1,
                                                                      eps);
      this->checkConversion<TangentOperator::ABAQUS,
                            TangentOperator::C_TAU_JAUMANN, N, real>(F0, F1,
                                                                     eps);
      this->checkConversion<TangentOperator::DPK1_DF, TangentOperator::DSIG_DF,
                            N, real>(F0, F1, eps);
      this->checkConversion<TangentOperator::DPK1_DF, TangentOperator::DTAU_DF,
                            N, real>(F0, F1, eps);
      this->checkLogarithmicStrainConversion<N, real>(F0, F1, eps);
    }
  }  // end of check
};

TFEL_TESTS_GENERATE_PROXY(FiniteStrainBehaviourTangentOperator17,
                          "FiniteStrainBehaviourTangentOperator17");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("FiniteStrainBehaviourTangentOperator17.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
