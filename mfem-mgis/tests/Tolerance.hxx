/*!
 * \file   Tolerance.hxx
 * \brief  This file declares functions adapting the tolerances of the tests,
 * written for the double numerical type, to the numerical type of mfem-mgis
 * \date   13/10/2026
 */

#ifndef LIB_MFEM_MGIS_TESTS_TOLERANCE_HXX
#define LIB_MFEM_MGIS_TESTS_TOLERANCE_HXX

#include <limits>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::unit_tests {

  /*!
   * \return a tolerance written for the double numerical type, bounded below
   * by ten times the machine epsilon of the numerical type of mfem-mgis: it
   * is unchanged in double (for a tolerance greater than 2.3e-15), and about
   * 1.2e-6 in float
   * \param[in] t: tolerance for the double numerical type
   */
  constexpr mfem_mgis::real tolerance(const double t) noexcept {
    constexpr auto e =
        10 *
        static_cast<double>(std::numeric_limits<mfem_mgis::real>::epsilon());
    return static_cast<mfem_mgis::real>(t > e ? t : e);
  }  // end of tolerance

  /*!
   * \return the tolerance of a test for the numerical type of mfem-mgis: the
   * first argument for double, the second one otherwise. The tolerance for
   * float is given explicitly when the error depends on the computation
   * (roundings accumulated by the assembly and the solvers, convergence
   * criteria of the solvers and of the behaviours bounded for this type).
   * \param[in] for_double: tolerance for the double numerical type
   * \param[in] for_float: tolerance for the float numerical type
   */
  constexpr mfem_mgis::real tolerance(const double for_double,
                                      const double for_float) noexcept {
    if constexpr (sizeof(mfem_mgis::real) == sizeof(double)) {
      return static_cast<mfem_mgis::real>(for_double);
    } else {
      return static_cast<mfem_mgis::real>(for_float);
    }
  }  // end of tolerance

}  // end of namespace mfem_mgis::unit_tests

#endif /* LIB_MFEM_MGIS_TESTS_TOLERANCE_HXX */
