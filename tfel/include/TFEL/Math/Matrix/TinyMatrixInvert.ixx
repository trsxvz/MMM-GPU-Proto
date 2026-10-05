/*!
 * \file  include/TFEL/Math/Matrix/TinyMatrixInvert.ixx
 * \brief
 * \author Thomas Helfer
 * \brief 08 mars 2013
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_MATH_TINYMATRIXINVERTIXX
#define LIB_TFEL_MATH_TINYMATRIXINVERTIXX

namespace tfel::math {

  template <unsigned short N, typename T>
  TFEL_HOST_DEVICE bool TinyMatrixInvert<N, T>::exe(
      tfel::math::tmatrix<N, N, T>& m, const T eps) {
    constexpr auto zero = T(0);
    constexpr auto one = T(1);
    auto t = tmatrix<N, N, T>{m};
    auto p = TinyPermutation<N>{};
    auto r = TinyMatrixSolve<N, T, true>::decomp(t, p, eps);
    for (unsigned short i = 0; (r) && (i != N); ++i) {
      auto e = tvector<N, T>{zero};
      e(i) = one;
      r = TinyMatrixSolve<N, T, true>::back_substitute(t, p, e);
      for (unsigned short j = 0; j != N; ++j) {
        m(j, i) = e(j);
      }
    }
    if (!r) {
      m = tmatrix<N, N, T>{std::numeric_limits<T>::quiet_NaN()};
    }
    return r;
  }  // end of TinyMatrixInvert::exe

}  // end of namespace tfel::math

#endif /* LIB_TFEL_MATH_TINYMATRIXINVERTIXX */
