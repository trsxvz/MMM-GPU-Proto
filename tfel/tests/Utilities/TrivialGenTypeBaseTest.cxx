/*!
 * \file  tests/Utilities/TrivialGenTypeBaseTest.cxx
 * \brief This file tests the `TrivialGenTypeBase` class.
 * \date  05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>
#include "TFEL/Utilities/TrivialGenTypeBase.hxx"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

//! \brief a trivially copyable type
struct Point {
  double x;
  double y;
};

using TrivialGenType = tfel::utilities::TrivialGenTypeBase<
    tfel::meta::GenerateTypeList<int, Point, std::array<double, 4>, double*>::
        type>;

static_assert(std::is_trivially_copyable_v<TrivialGenType>);

//! \return true if both values are equal
static bool areEqual(const double a, const double b) {
  return std::abs(a - b) < 1e-14;
}

struct TrivialGenTypeBaseTest final : public tfel::tests::TestCase {
  TrivialGenTypeBaseTest()
      : tfel::tests::TestCase("TFEL/Utilities", "TrivialGenTypeBaseTest") {
  }  // end of TrivialGenTypeBaseTest
  tfel::tests::TestResult execute() override {
    auto g = TrivialGenType{};
    TFEL_TESTS_ASSERT(g.empty());
    TFEL_TESTS_CHECK_EQUAL(g.getTypeIndex(), 0);
    g = 21;
    TFEL_TESTS_ASSERT(g.is<int>());
    TFEL_TESTS_ASSERT(!g.is<Point>());
    TFEL_TESTS_CHECK_EQUAL(g.get<int>(), 21);
    TFEL_TESTS_CHECK_EQUAL(g.getTypeIndex(), 1);
    g = Point{1, 2};
    TFEL_TESTS_ASSERT(g.is<Point>());
    TFEL_TESTS_CHECK_EQUAL(g.getTypeIndex(), 2);
    TFEL_TESTS_ASSERT(areEqual(g.get<Point>().y, 2));
    // conversion operator
    const Point& p = g;
    TFEL_TESTS_ASSERT(areEqual(p.x, 1));
    // copy
    const auto g2 = g;
    TFEL_TESTS_ASSERT(g2.is<Point>());
    TFEL_TESTS_ASSERT(areEqual(g2.get<Point>().x, 1));
    g.set(std::array<double, 4>{1, 2, 3, 4});
    TFEL_TESTS_ASSERT((g.is<std::array<double, 4>>()));
    TFEL_TESTS_ASSERT(areEqual(g.get<std::array<double, 4>>()[3], 4));
    // the copy is not modified
    TFEL_TESTS_ASSERT(g2.is<Point>());
    auto v = 3.;
    g = &v;
    TFEL_TESTS_ASSERT(g.is<double*>());
    *(g.get<double*>()) = 4;
    TFEL_TESTS_ASSERT(areEqual(v, 4));
    // access to a value of another type
    TFEL_TESTS_CHECK_THROW(static_cast<void>(g.get<int>()),
                           tfel::utilities::GenTypeCastError);
    return this->result;
  }  // end of execute()
};

TFEL_TESTS_GENERATE_PROXY(TrivialGenTypeBaseTest, "TrivialGenTypeBase");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("TrivialGenTypeBaseTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
