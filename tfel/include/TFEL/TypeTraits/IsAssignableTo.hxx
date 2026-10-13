/*!
 * \file   include/TFEL/TypeTraits/IsAssignableTo.hxx
 * \brief  This file implements the IsAssignableTo class traits and
 * its partial specialisation for each standard numeric types.
 *
 * \author Thomas Helfer
 * \date   28 jun 2006
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_TYPETRAITS_ISASSIGNABLETO_HXX
#define LIB_TFEL_TYPETRAITS_ISASSIGNABLETO_HXX

#include <type_traits>
#include "TFEL/Config/TFELConfig.hxx"
#include "TFEL/TypeTraits/Promote.hxx"
#include "TFEL/Math/Forward/Complex.hxx"

namespace tfel::typetraits::internals {

  /*!
   * \return if a value of type `X` is assignable to an object of type `Y`:
   * either `X` is promoted to `Y`, or both are floating-point types, the
   * value being then converted, possibly to a narrower type (`double` to
   * `float`). The latter allows the code written with `double` literals to
   * be compiled for `float` objects.
   * \tparam X: type of the value
   * \tparam Y: type of the object
   */
  template <typename X, typename Y>
  constexpr bool isStandardNumericTypeAssignableTo() {
    return std::is_same_v<typename Promote<X, Y>::type, Y> ||
           (std::is_floating_point_v<X> && std::is_floating_point_v<Y>);
  }  // end of isStandardNumericTypeAssignableTo

}  // end of namespace tfel::typetraits::internals

/*!
 * \def TFEL_MATH_IS_ASSIGNABLE_
 * \brief An helper macro to specialise IsAssignableTo for each
 * standard numeric types.
 * \author Thomas Helfer
 * \date   28 jun 2006
 */
#define TFEL_MATH_IS_ASSIGNABLE_(X, Y)                        \
  /*!                                                         \
   * \brief Partial specialisation for X and Y                \
   * \see   IsAssignableTo                                    \
   */                                                         \
  template <>                                                 \
  struct IsAssignableTo<X, Y> {                               \
    /*! \brief result */                                      \
    static constexpr bool value =                             \
        internals::isStandardNumericTypeAssignableTo<X, Y>(); \
    /*! \brief result */                                      \
    static constexpr bool cond =                              \
        internals::isStandardNumericTypeAssignableTo<X, Y>(); \
  }

/*!
 * \def TFEL_MATH_IS_ASSIGNABLE
 * \brief An helper macro to specialise IsAssignableTo for each
 * standard numeric types.
 * \author Thomas Helfer
 * \date   28 jun 2006
 */
#define TFEL_MATH_IS_ASSIGNABLE(X)                                     \
  TFEL_MATH_IS_ASSIGNABLE_(X, unsigned short);                         \
  TFEL_MATH_IS_ASSIGNABLE_(X, unsigned int);                           \
  TFEL_MATH_IS_ASSIGNABLE_(X, long unsigned int);                      \
  TFEL_MATH_IS_ASSIGNABLE_(X, short);                                  \
  TFEL_MATH_IS_ASSIGNABLE_(X, int);                                    \
  TFEL_MATH_IS_ASSIGNABLE_(X, long int);                               \
  TFEL_MATH_IS_ASSIGNABLE_(X, float);                                  \
  TFEL_MATH_IS_ASSIGNABLE_(X, double);                                 \
  TFEL_MATH_IS_ASSIGNABLE_(X, long double);                            \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<unsigned short>);    \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<unsigned int>);      \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<long unsigned int>); \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<short>);             \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<int>);               \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<long int>);          \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<float>);             \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<double>);            \
  TFEL_MATH_IS_ASSIGNABLE_(X, tfel::math::Complex<long double>)

namespace tfel::typetraits {

  /*!
   * \class IsAssignableTo
   * \brief Traits class which says if an object a of type A is
   * assignable to an object b of type B, ie if b = a is valid. For the
   * standard numeric types, a value is assignable if its type is promoted
   * to the type of the object, or if both types are floating-point types
   * (`double` to `float`, with a conversion).
   *
   * \tparam  A, first type
   * \tparam  B, second type
   * \return bool cond, true if b = a is valid, false otherwise.
   *
   * \see Promote
   * \author Thomas Helfer
   * \date   28 jun 2006
   */
  template <typename A, typename B>
  struct IsAssignableTo {
    //! \brief result
    static constexpr bool value = false;
    //! \brief result
    static constexpr bool cond = false;
  };

  TFEL_MATH_IS_ASSIGNABLE(unsigned short);
  TFEL_MATH_IS_ASSIGNABLE(unsigned int);
  TFEL_MATH_IS_ASSIGNABLE(long unsigned int);
  TFEL_MATH_IS_ASSIGNABLE(short);
  TFEL_MATH_IS_ASSIGNABLE(int);
  TFEL_MATH_IS_ASSIGNABLE(long int);
  TFEL_MATH_IS_ASSIGNABLE(float);
  TFEL_MATH_IS_ASSIGNABLE(double);
  TFEL_MATH_IS_ASSIGNABLE(long double);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<unsigned short>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<unsigned int>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<long unsigned int>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<short>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<int>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<long int>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<float>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<double>);
  TFEL_MATH_IS_ASSIGNABLE(tfel::math::Complex<long double>);

  /*!
   * \brief an helper function around `IsAssignableTo`
   * \tparam  A, first type
   * \tparam  B, second type
   */
  template <typename A, typename B>
  TFEL_HOST_DEVICE constexpr bool isAssignableTo() {
    return IsAssignableTo<std::remove_cv_t<A>, B>::cond;
  }  // end of isAssignableTo

}  // end of namespace tfel::typetraits

#endif /* LIB_TFEL_TYPETRAITS_ISASSIGNABLETO_HXX */
