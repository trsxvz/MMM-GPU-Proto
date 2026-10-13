/*!
 * \file   src/AccelerationAlgorithms.hxx
 * \brief  This file declares the generators of the acceleration algorithms
 * declared by default in the `AccelerationAlgorithmFactory` (internal
 * header)
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_ACCELERATIONALGORITHMS_HXX
#define LIB_MFEM_MGIS_SRC_ACCELERATIONALGORITHMS_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearSolvers/AccelerationAlgorithm.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return Anderson's method (see `AccelerationAlgorithmFactory`)
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   */
  [[nodiscard]] std::unique_ptr<AccelerationAlgorithm>
  makeAndersonAccelerationAlgorithm(Context& ctx,
                                    const Parameters& opts) noexcept;
  /*!
   * \return Anderson's method on the fixed point residuals, accelerating
   * every two iterations by default, as `UAnderson` in MTest
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   */
  [[nodiscard]] std::unique_ptr<AccelerationAlgorithm>
  makeUAndersonAccelerationAlgorithm(Context& ctx,
                                     const Parameters& opts) noexcept;
  /*!
   * \return Anderson's method on the residuals of the problem, accelerating
   * every two iterations by default, as `FAnderson` in MTest
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   */
  [[nodiscard]] std::unique_ptr<AccelerationAlgorithm>
  makeFAndersonAccelerationAlgorithm(Context& ctx,
                                     const Parameters& opts) noexcept;
  /*!
   * \brief declares the algorithms based on secant and delta squared methods
   * (`Secant`, `CrossedSecant`, ...)
   * \param[in, out] f: factory
   */
  void declareSecantAccelerationAlgorithms(
      AccelerationAlgorithmFactory& f) noexcept;
  /*!
   * \return the acceleration of Cast3M (see `AccelerationAlgorithmFactory`)
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   */
  [[nodiscard]] std::unique_ptr<AccelerationAlgorithm>
  makeCast3MAccelerationAlgorithm(Context& ctx,
                                  const Parameters& opts) noexcept;

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_ACCELERATIONALGORITHMS_HXX */
