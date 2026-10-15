/*!
 * \file   include/MFEMMGIS/LinearSolverConvergenceReport.hxx
 * \brief  This file declares the `LinearSolverConvergenceReport` class
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_LINEARSOLVERCONVERGENCEREPORT_HXX
#define LIB_MFEM_MGIS_LINEARSOLVERCONVERGENCEREPORT_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief interface of the linear solvers which report the outcome of their
   * last resolution. The functions `hasConverged` and
   * `getNumberOfIterationsAtConvergence` (see `SolverUtilities.hxx`) use it
   * for the solvers which are neither iterative solvers of MFEM nor Krylov
   * solvers of hypre: a linear solver wrapping another library implements
   * this interface, so that these functions need no change when a library
   * is added.
   */
  struct MFEM_MGIS_EXPORT LinearSolverConvergenceReport {
    /*!
     * \return true if the last resolution converged (iterative solvers) or
     * succeeded (direct solvers)
     */
    [[nodiscard]] virtual bool hasConverged() const noexcept = 0;
    /*!
     * \return the number of iterations of the last resolution, if
     * meaningful (iterative solvers)
     */
    [[nodiscard]] virtual std::optional<int> getNumberOfIterations()
        const noexcept = 0;
    //! \brief destructor
    virtual ~LinearSolverConvergenceReport();
  };  // end of struct LinearSolverConvergenceReport

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LINEARSOLVERCONVERGENCEREPORT_HXX */
