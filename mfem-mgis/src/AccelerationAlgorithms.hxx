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
  makeAndersonAccelerationAlgorithm(Context&, const Parameters&) noexcept;
  /*!
   * \return the acceleration of Cast3M (see `AccelerationAlgorithmFactory`)
   * \param[in, out] ctx: execution context
   * \param[in] opts: options
   */
  [[nodiscard]] std::unique_ptr<AccelerationAlgorithm>
  makeCast3MAccelerationAlgorithm(Context&, const Parameters&) noexcept;

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_ACCELERATIONALGORITHMS_HXX */
