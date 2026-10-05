/*!
 * \file   tests/Math/ChangeBasisTest.cxx
 * \brief  This file tests the rotation of vectors and square matrices, which
 * must be consistent with the rotation of symmetric tensors.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TFEL/Math/tvector.hxx"
#include "TFEL/Math/tmatrix.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

struct ChangeBasisTest final : public tfel::tests::TestCase {
  ChangeBasisTest() : tfel::tests::TestCase("TFEL/Math", "ChangeBasisTest") {}
  tfel::tests::TestResult execute() override {
    this->check<2>();
    this->check<3>();
    return this->result;
  }  // end of execute

 private:
  //! \return a rotation matrix (rotation around the z axis in 2D)
  template <unsigned short N>
  static tfel::math::rotation_matrix<double> getRotationMatrix() {
    using namespace tfel::math;
    const auto a = 0.3;
    const auto rz = rotation_matrix<double>{std::cos(a),
                                            -std::sin(a),
                                            0,  //
                                            std::sin(a),
                                            std::cos(a),
                                            0,  //
                                            0,
                                            0,
                                            1};
    if constexpr (N == 2) {
      return rz;
    } else {
      const auto b = 0.7;
      const auto rx = rotation_matrix<double>{1, 0,           0,             //
                                              0, std::cos(b), -std::sin(b),  //
                                              0, std::sin(b), std::cos(b)};
      return rz * rx;
    }
  }
  //! \return the symmetric tensor v x v
  template <unsigned short N>
  static tfel::math::stensor<N, double> square(
      const tfel::math::tvector<N, double>& v) {
    constexpr auto cste = tfel::math::Cste<double>::sqrt2;
    if constexpr (N == 2) {
      return {v(0) * v(0), v(1) * v(1), 0, cste * v(0) * v(1)};
    } else {
      return {v(0) * v(0),        v(1) * v(1),        v(2) * v(2),
              cste * v(0) * v(1), cste * v(0) * v(2), cste * v(1) * v(2)};
    }
  }
  //! \brief check the rotations in the given space dimension
  template <unsigned short N>
  void check() {
    using namespace tfel::math;
    constexpr auto eps = 1e-14;
    const auto r = getRotationMatrix<N>();
    auto a = tvector<N, double>{};
    auto b = tvector<N, double>{};
    for (unsigned short i = 0; i != N; ++i) {
      a(i) = 1 + i;
      b(i) = 2 - 3 * i;
    }
    // vectors: (r.a) x (r.a) must be the rotated symmetric tensor a x a
    const auto ra = change_basis(a, r);
    const auto s = change_basis(square(a), r);
    const auto s2 = square(ra);
    for (unsigned short i = 0; i != StensorDimeToSize<N>::value; ++i) {
      TFEL_TESTS_ASSERT(std::abs(s(i) - s2(i)) < eps);
    }
    // the norm is preserved
    TFEL_TESTS_ASSERT(std::abs(norm(ra) - norm(a)) < eps);
    // matrices: the rotation of a x b is (r.a) x (r.b)
    auto m = tmatrix<N, N, double>{};
    for (unsigned short i = 0; i != N; ++i) {
      for (unsigned short j = 0; j != N; ++j) {
        m(i, j) = a(i) * b(j);
      }
    }
    const auto rm = change_basis(m, r);
    const auto rb = change_basis(b, r);
    for (unsigned short i = 0; i != N; ++i) {
      for (unsigned short j = 0; j != N; ++j) {
        TFEL_TESTS_ASSERT(std::abs(rm(i, j) - ra(i) * rb(j)) < eps);
      }
    }
  }  // end of check
};

TFEL_TESTS_GENERATE_PROXY(ChangeBasisTest, "ChangeBasisTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ChangeBasisTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
