/*!
 * \file   MFEMMGIS/NonLinearSolvers/AccelerationAlgorithm.hxx
 * \brief  This file declares the `AccelerationAlgorithm` class, the
 * interface of the acceleration algorithms of the fixed point iterations
 * (see `FixedPointSolver`), and their factory
 * \date   10/10/2026
 */

#ifndef LIB_MFEMMGIS_NONLINEARSOLVERS_ACCELERATIONALGORITHM_HXX
#define LIB_MFEMMGIS_NONLINEARSOLVERS_ACCELERATIONALGORITHM_HXX

#include <map>
#include <memory>
#include <string>
#include <functional>
#include <string_view>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief interface of the acceleration algorithms of the fixed point
   * iterations `x_{i+1} = x_i - S r(x_i)`, `r` being the residual and `S` an
   * approximation of the inverse of the tangent operator (for instance the
   * inverse of the elastic operator).
   *
   * The algorithms follow the ones of MTest (TFEL): they are given, at each
   * iteration, the new estimate `x_{i+1}`, the correction `x_{i+1} - x_i`
   * and the residual `r(x_i)`, and may replace the new estimate by an
   * accelerated one. They only use linear combinations and scalar products
   * of vectors (`mfem::Vector`, on the device if the vectors are), the
   * scalar products being reduced over all the processes by the function
   * given to `initialize`.
   */
  struct MFEM_MGIS_EXPORT AccelerationAlgorithm {
    //! \brief scalar product reduced over all the processes
    using ScalarProduct =
        std::function<real(const mfem::Vector&, const mfem::Vector&)>;
    /*!
     * \brief initializes the algorithm at the beginning of a resolution:
     * the histories of the previous resolutions are discarded
     * \param[in] n: number of unknowns of the current process
     * \param[in] dot: scalar product
     */
    virtual void initialize(const size_type, ScalarProduct) = 0;
    /*!
     * \brief accelerates the fixed point iterations
     * \param[in, out] x: new estimate `x_{i+1}`, replaced by the accelerated
     * estimate if the algorithm applies at this iteration
     * \param[in] dx: correction `x_{i+1} - x_i`
     * \param[in] r: residual at the previous estimate `x_i`
     * \param[in] r_eps: scale of the residual (the convergence criterion on
     * the residual), below which differences of residuals are considered
     * null
     * \param[in] u_eps: scale of the unknowns (an estimate of the variation
     * of the unknowns when the residual reaches its convergence criterion),
     * below which differences of unknowns are considered null
     * \param[in] i: iteration number, starting at 1
     */
    virtual void execute(mfem::Vector&,
                         const mfem::Vector&,
                         const mfem::Vector&,
                         const real,
                         const real,
                         const size_type) = 0;
    //! \brief destructor
    virtual ~AccelerationAlgorithm();
  };  // end of struct AccelerationAlgorithm

  /*!
   * \brief factory of the acceleration algorithms. The algorithms are
   * selected by their names and set by their options:
   *
   * - `Anderson`: Anderson's method on the images of the fixed point map,
   *   options `Depth` (number of previous iterations used, 4 by default),
   *   `Period` (number of iterations between two accelerations, 1 by
   *   default), `Residuals` (`Unknowns`, by default, the fixed point
   *   residuals, or `Forces`, the residuals of the problem) and
   *   `Regularization` (relative regularization of the least squares
   *   problem, 1e-12 by default);
   * - `UAnderson` and `FAnderson`: Anderson's method on the fixed point
   *   residuals and on the residuals of the problem, accelerating every two
   *   iterations by default, as MTest (same options as `Anderson`, except
   *   `Residuals`);
   * - `Cast3M`: port of the acceleration of Cast3M (as MTest), options
   *   `Trigger` (first accelerated iteration, 4 by default) and `Period`
   *   (number of iterations between two accelerations, 2 by default);
   * - `Secant`, `CrossedSecant`, `AlternateSecant`, `IronsTuck`,
   *   `Steffensen`, `CrossedDelta2`, `Crossed2Delta`, `Crossed2Deltabis`,
   *   `AlternateDelta2` and `Alternate2Delta`: ports of the algorithms of
   *   MTest of the same names, option `Trigger` (first accelerated
   *   iteration, with the default values of MTest).
   */
  struct MFEM_MGIS_EXPORT AccelerationAlgorithmFactory {
    //! \brief generator of an acceleration algorithm from its options
    using Generator = std::function<std::unique_ptr<AccelerationAlgorithm>(
        Context&, const Parameters&)>;
    //! \return the unique instance of the factory
    static AccelerationAlgorithmFactory& get() noexcept;
    /*!
     * \brief declares a new acceleration algorithm
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the algorithm
     * \param[in] g: generator
     */
    [[nodiscard]] bool add(Context&, std::string_view, Generator) noexcept;
    /*!
     * \return the acceleration algorithm of the given name, set by the given
     * options, or a null pointer on error
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the algorithm
     * \param[in] opts: options
     */
    [[nodiscard]] std::unique_ptr<AccelerationAlgorithm> generate(
        Context&, std::string_view, const Parameters&) const noexcept;

   private:
    //! \brief default constructor, declaring the default algorithms
    AccelerationAlgorithmFactory() noexcept;
    //! \brief destructor
    ~AccelerationAlgorithmFactory() noexcept;
    //! \brief generators
    std::map<std::string, Generator, std::less<>> generators;
  };  // end of struct AccelerationAlgorithmFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEARSOLVERS_ACCELERATIONALGORITHM_HXX */
