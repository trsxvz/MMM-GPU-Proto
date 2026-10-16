/*!
 * \file   SequentialLinearSolver.hxx
 * \brief  This file declares the linear solver used by the sequential tests
 * in place of a direct solver
 * \date   13/10/2026
 */

#ifndef LIB_MFEM_MGIS_TESTS_SEQUENTIALLINEARSOLVER_HXX
#define LIB_MFEM_MGIS_TESTS_SEQUENTIALLINEARSOLVER_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis::unit_tests {

  /*!
   * \brief linear solver of the sequential tests: UMFPACK if MFEM provides
   * it, otherwise (MFEM in single precision, which refuses SuiteSparse) the
   * GMRES solver of MFEM with the tightest relative tolerance, preconditioned
   * by the Jacobi smoother (without it, GMRES may stagnate in float, the
   * diagonal terms of the eliminated degrees of freedom, equal to one,
   * being much smaller than the other ones)
   */
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
  inline constexpr const char* sequential_direct_solver = "UMFPackSolver";
#else  /* MFEM_MGIS_HAVE_SUITESPARSE */
  inline constexpr const char* sequential_direct_solver = "GMRESSolver";
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */

  //! \brief if true, `sequential_direct_solver` is a direct solver
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
  inline constexpr bool has_sequential_direct_solver = true;
#else  /* MFEM_MGIS_HAVE_SUITESPARSE */
  inline constexpr bool has_sequential_direct_solver = false;
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */

  //! \return the parameters of the linear solver of the sequential tests
  [[nodiscard]] inline Parameters getSequentialDirectSolverParameters() {
#ifdef MFEM_MGIS_HAVE_SUITESPARSE
    return {};
#else  /* MFEM_MGIS_HAVE_SUITESPARSE */
    return {{"VerbosityLevel", 0},
            {"RelativeTolerance", 1e-12},
            {"AbsoluteTolerance", 0.},
            {"MaximumNumberOfIterations", 10000},
            {"Preconditioner", Parameters{{"Name", "OperatorJacobiSmoother"}}}};
#endif /* MFEM_MGIS_HAVE_SUITESPARSE */
  }  // end of getSequentialDirectSolverParameters

}  // end of namespace mfem_mgis::unit_tests

#endif /* LIB_MFEM_MGIS_TESTS_SEQUENTIALLINEARSOLVER_HXX */
