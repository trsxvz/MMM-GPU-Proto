/*!
 * \file   mfront/include/MFront/GenericBehaviour/RotateTangentOperatorBlock.hxx
 * \brief  This file declares the `rotateTangentOperatorBlock` function.
 * \date   08/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_ROTATETANGENTOPERATORBLOCK_HXX
#define LIB_MFRONT_GENERICBEHAVIOUR_ROTATETANGENTOPERATORBLOCK_HXX

#include <type_traits>
#include "TFEL/Config/TFELConfig.hxx"
#include "TFEL/Math/tmatrix.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/tensor.hxx"

namespace mfront::gb {

  namespace internals {

    /*!
     * \return the number of components of a scalar, a vector, a symmetric
     * tensor or a tensor
     * \tparam Type: type of the object
     */
    template <typename Type>
    TFEL_HOST_DEVICE constexpr unsigned short getNumberOfComponents() noexcept {
      if constexpr (std::is_arithmetic_v<Type>) {
        return 1;
      } else {
        return static_cast<unsigned short>(Type{}.size());
      }
    }  // end of getNumberOfComponents

    /*!
     * \return the matrix of the rotation of the objects of the given type by
     * `tfel::math::change_basis`: its k-th column is the rotated k-th basis
     * object. This matrix is orthogonal in the storage of TFEL.
     * \tparam Type: type of the object (scalar, vector, symmetric tensor or
     * tensor)
     * \param[in] r: rotation matrix
     */
    template <typename Type, typename real>
    TFEL_HOST_DEVICE constexpr auto getRotationOperator(
        const tfel::math::rotation_matrix<real>& r) noexcept {
      constexpr auto s = getNumberOfComponents<Type>();
      auto q = tfel::math::tmatrix<s, s, real>{real{0}};
      if constexpr (std::is_arithmetic_v<Type>) {
        q(0, 0) = 1;
      } else {
        for (unsigned short k = 0; k != s; ++k) {
          auto e = Type{real{0}};
          e[k] = 1;
          const auto re = tfel::math::change_basis(e, r);
          for (unsigned short i = 0; i != s; ++i) {
            q(i, k) = re[i];
          }
        }
      }
      return q;
    }  // end of getRotationOperator

  }  // end of namespace internals

  /*!
   * \brief rotate a tangent operator block, i.e. the derivative of an object
   * of type `Type1` with respect to an object of type `Type2`, stored by
   * rows. Each object is a scalar, a vector, a symmetric tensor or a tensor.
   *
   * The rotated block is \f$Q_{1}\,K\,Q_{2}^{T}\f$, where \f$Q_{i}\f$ is the
   * matrix of the rotation of the objects of type `Type_i` by
   * `tfel::math::change_basis`. This function rotates the blocks for which
   * TFEL has no type which can be rotated (derivative of a vector with
   * respect to a symmetric tensor or a tensor, and conversely, derivative
   * of a tensor with respect to a symmetric tensor).
   *
   * \tparam Type1: type of the derivated object
   * \tparam Type2: type of the object with respect to which the derivative
   * is taken
   * \param[out] dest: rotated block, which may be the source
   * \param[in] src: block
   * \param[in] r: rotation matrix
   */
  template <typename Type1, typename Type2, typename real>
  TFEL_HOST_DEVICE void rotateTangentOperatorBlock(
      real* const dest,
      const real* const src,
      const tfel::math::rotation_matrix<real>& r) noexcept {
    constexpr auto s1 = internals::getNumberOfComponents<Type1>();
    constexpr auto s2 = internals::getNumberOfComponents<Type2>();
    // the numeric type can't be deduced from `rotation_matrix`, which is an
    // alias
    const auto q1 = internals::getRotationOperator<Type1, real>(r);
    const auto q2 = internals::getRotationOperator<Type2, real>(r);
    // the block is read before being written, so that the destination may
    // be the source
    auto k = tfel::math::tmatrix<s1, s2, real>{};
    for (unsigned short a = 0; a != s1; ++a) {
      for (unsigned short b = 0; b != s2; ++b) {
        k(a, b) = src[a * s2 + b];
      }
    }
    // q1 * k * transpose(q2)
    for (unsigned short a = 0; a != s1; ++a) {
      for (unsigned short b = 0; b != s2; ++b) {
        auto v = real{0};
        for (unsigned short c = 0; c != s1; ++c) {
          for (unsigned short d = 0; d != s2; ++d) {
            v += q1(a, c) * k(c, d) * q2(b, d);
          }
        }
        dest[a * s2 + b] = v;
      }
    }
  }  // end of rotateTangentOperatorBlock

}  // end of namespace mfront::gb

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_ROTATETANGENTOPERATORBLOCK_HXX */
