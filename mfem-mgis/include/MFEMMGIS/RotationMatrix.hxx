/*!
 * \file   include/MFEMMGIS/RotationMatrix.hxx
 * \brief  This file declares the `RotationMatrix2D`, `RotationMatrix3D` and
 * `MaterialAxis3D` aliases
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_ROTATIONMATRIX_HXX
#define LIB_MFEM_MGIS_ROTATIONMATRIX_HXX

#include <array>
#include <tuple>
#include <variant>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  /*!
   * \brief definition of the rotation matrix in 2D
   *
   * The rotation matrix is defined by:
   * - a constant matrix, computed once for all
   * - the definition of the first direction of orthotropy
   *
   * A direction given by a partial quadrature function is shared with the
   * material (`Material::setRotationMatrix`), which keeps the rotation
   * matrices computed from its values: if these values are modified
   * afterwards, `Material::invalidateRotationMatrices` must be called.
   */
  using RotationMatrix2D =
      std::variant<std::monostate,
                   std::array<real, 9u>,
                   std::shared_ptr<PartialQuadratureFunction>>;
  /*!
   * \brief a direction of orthotropy in 3D, uniform or given at each
   * integration point by a partial quadrature function (see
   * `RotationMatrix3D`)
   */
  using MaterialAxis3D =
      std::variant<std::array<real, 3u>,
                   std::shared_ptr<PartialQuadratureFunction>>;
  /*!
   * \brief definition of the rotation matrix in 3D
   *
   * The rotation matrix is defined by:
   * - a constant matrix, computed once for all
   * - the two directions of orthotropy which can be given either as a
   * partial quadrature function or as a fixed direction.
   *
   * A direction given by a partial quadrature function is shared with the
   * material (`Material::setRotationMatrix`), which keeps the rotation
   * matrices computed from its values: if these values are modified
   * afterwards, `Material::invalidateRotationMatrices` must be called.
   */
  using RotationMatrix3D = std::variant<std::monostate,
                                        std::array<real, 9u>,
                                        std::array<MaterialAxis3D, 2u>>;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ROTATIONMATRIX_HXX */
