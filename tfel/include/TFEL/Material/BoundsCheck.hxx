/*!
 * \file   include/TFEL/Material/BoundsCheck.hxx
 * \brief  This file declares the  class
 * \author Thomas Helfer
 * \date 06/03/2007
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_BOUNDSCHECK_HXX
#define LIB_TFEL_BOUNDSCHECK_HXX

#include <string>
#include <string_view>

#include "TFEL/Config/TFELConfig.hxx"
#include "TFEL/Math/qt.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/TypeTraits/BaseType.hxx"
#include "TFEL/Material/OutOfBoundsPolicy.hxx"

// fixing a bug on current glibc++ cygwin versions (19/08/2015)
#if defined __CYGWIN__ && (!defined _GLIBCXX_USE_C99)
#include <sstream>
namespace std {
  template <typename T>
  std::string to_string(const T& v) {
    std::ostringstream s;
    s << v;
    return s.str();
  }
}  // namespace std
#endif /* defined __CYGWIN__ &&  (!defined _GLIBCXX_USE_C99) */

namespace tfel::material {

  /*!
   * \brief a helper class
   *
   * All the checks return a boolean which is false if a bound is violated
   * and if this violation is not tolerated by the out of bounds policy.
   *
   * On the host, an exception is thrown in this case, so the returned value
   * is always true. In device code, where exceptions are not available, the
   * returned value must be used to report the failure. Warnings are not
   * displayed in device code. The compilers translating host and device
   * code in a single pass follow the rules of device code in the whole
   * translation unit (see `TFEL_DEVICE_CODE`).
   *
   * The names of the variables are passed as C strings: constructing a
   * `std::string_view` from a C string is not supported in device code by
   * some compilers (nvcc silently removes the code).
   */
  struct BoundsCheckBase {
    /*!
     * \brief convert to string
     * \param[in] value: value to be converted
     */
    template <typename T>
    static std::string convert_to_string(const T& value) {
      if constexpr (tfel::math::ImmutableQuantityConcept<T>) {
        return std::to_string(tfel::math::base_type_cast(value));
      } else {
        return std::to_string(value);
      }
    }

    /*!
     * \param[in] n : variable name
     * \param[in] v : value
     * \param[in] b : lower bound value
     */
    [[noreturn]] TFEL_VISIBILITY_EXPORT static void
    throwOutOfLowerBoundsException(const std::string_view,
                                   const std::string_view,
                                   const std::string_view);
    /*!
     * \param[in] n : variable name
     * \param[in] v : value
     * \param[in] b : upper bound value
     */
    [[noreturn]] TFEL_VISIBILITY_EXPORT static void
    throwOutOfUpperBoundsException(const std::string_view,
                                   const std::string_view,
                                   const std::string_view);
    /*!
     * \param[in] n  : variable name
     * \param[in] v  : value
     * \param[in] lb : lower bound value
     * \param[in] ub : upper bound value
     */
    [[noreturn]] TFEL_VISIBILITY_EXPORT static void throwOutOfBoundsException(
        const std::string_view,
        const std::string_view,
        const std::string_view,
        const std::string_view);
    /*!
     * \param[in] n : variable name
     * \param[in] v : value
     * \param[in] b : lower bound value
     */
    TFEL_VISIBILITY_EXPORT static void displayOutOfLowerBoundsWarning(
        const std::string_view, const std::string_view, const std::string_view);
    /*!
     * \param[in] n : variable name
     * \param[in] v : value
     * \param[in] b : upper bound value
     */
    TFEL_VISIBILITY_EXPORT static void displayOutOfUpperBoundsWarning(
        const std::string_view, const std::string_view, const std::string_view);
    /*!
     * \param[in] n  : variable name
     * \param[in] v  : value
     * \param[in] lb : lower bound value
     * \param[in] ub : upper bound value
     */
    TFEL_VISIBILITY_EXPORT static void displayOutOfBoundsWarning(
        const std::string_view,
        const std::string_view,
        const std::string_view,
        const std::string_view);
    template <typename T>
    TFEL_HOST_DEVICE static bool lowerBoundCheck(
        const char* const name,
        const T value,
        const T lBound,
        const OutOfBoundsPolicy p = Strict) {
      if (value < lBound) {
        return treatOutOfLowerBound(name, -1, value, lBound, p);
      }
      return true;
    }
    template <tfel::math::QuantityConcept QuantityType>
    TFEL_HOST_DEVICE static bool lowerBoundCheck(
        const char* const name,
        const QuantityType& value,
        const tfel::math::base_type<QuantityType> lBound,
        const OutOfBoundsPolicy p = Strict) {
      return lowerBoundCheck(name, tfel::math::base_type_cast(value), lBound,
                             p);
    }
    template <typename T>
    TFEL_HOST_DEVICE static bool upperBoundCheck(
        const char* const name,
        const T value,
        const T uBound,
        const OutOfBoundsPolicy p = Strict) {
      if (value > uBound) {
        return treatOutOfUpperBound(name, -1, value, uBound, p);
      }
      return true;
    }
    template <tfel::math::QuantityConcept QuantityType>
    TFEL_HOST_DEVICE static bool upperBoundCheck(
        const char* const name,
        const QuantityType& value,
        const tfel::math::base_type<QuantityType> uBound,
        const OutOfBoundsPolicy p = Strict) {
      return upperBoundCheck(name, tfel::math::base_type_cast(value), uBound,
                             p);
    }
    template <typename T>
    TFEL_HOST_DEVICE static bool lowerAndUpperBoundsChecks(
        const char* const name,
        const T value,
        const T lBound,
        const T uBound,
        const OutOfBoundsPolicy p = Strict) {
      if ((value < lBound) || (value > uBound)) {
        return treatOutOfBounds(name, -1, value, lBound, uBound, p);
      }
      return true;
    }
    template <tfel::math::QuantityConcept QuantityType>
    TFEL_HOST_DEVICE static bool lowerAndUpperBoundsChecks(
        const char* const name,
        const QuantityType& value,
        const tfel::math::base_type<QuantityType> lBound,
        const tfel::math::base_type<QuantityType> uBound,
        const OutOfBoundsPolicy p = Strict) {
      return lowerAndUpperBoundsChecks(name, tfel::math::base_type_cast(value),
                                       lBound, uBound, p);
    }

   protected:
    /*!
     * \return the name of the variable, or the name of the component of the
     * variable if the given index is positive or null
     * \param[in] n: variable name
     * \param[in] i: component index
     */
    static std::string getVariableName(const char* const n, const int i) {
      if (i < 0) {
        return std::string{n};
      }
      return std::string{n} + "(" + std::to_string(i) + ")";
    }
    /*!
     * \brief treat a value below its lower bound
     * \param[in] n: variable name
     * \param[in] i: component index, or -1 for a scalar variable
     * \param[in] v: value
     * \param[in] b: lower bound
     * \param[in] p: out of bounds policy
     */
    template <typename T>
    TFEL_HOST_DEVICE static bool treatOutOfLowerBound(
        [[maybe_unused]] const char* const n,
        [[maybe_unused]] const int i,
        [[maybe_unused]] const T v,
        [[maybe_unused]] const T b,
        const OutOfBoundsPolicy p) {
      if (p == None) {
        return true;
      }
#ifdef TFEL_DEVICE_CODE
      return p != Strict;
#else  /* TFEL_DEVICE_CODE */
      const auto name = getVariableName(n, i);
      if (p == Strict) {
        throwOutOfLowerBoundsException(name, convert_to_string(v),
                                       convert_to_string(b));
      }
      displayOutOfLowerBoundsWarning(name, convert_to_string(v),
                                     convert_to_string(b));
      return true;
#endif /* TFEL_DEVICE_CODE */
    }
    /*!
     * \brief treat a value above its upper bound
     * \param[in] n: variable name
     * \param[in] i: component index, or -1 for a scalar variable
     * \param[in] v: value
     * \param[in] b: upper bound
     * \param[in] p: out of bounds policy
     */
    template <typename T>
    TFEL_HOST_DEVICE static bool treatOutOfUpperBound(
        [[maybe_unused]] const char* const n,
        [[maybe_unused]] const int i,
        [[maybe_unused]] const T v,
        [[maybe_unused]] const T b,
        const OutOfBoundsPolicy p) {
      if (p == None) {
        return true;
      }
#ifdef TFEL_DEVICE_CODE
      return p != Strict;
#else  /* TFEL_DEVICE_CODE */
      const auto name = getVariableName(n, i);
      if (p == Strict) {
        throwOutOfUpperBoundsException(name, convert_to_string(v),
                                       convert_to_string(b));
      }
      displayOutOfUpperBoundsWarning(name, convert_to_string(v),
                                     convert_to_string(b));
      return true;
#endif /* TFEL_DEVICE_CODE */
    }
    /*!
     * \brief treat a value out of its bounds
     * \param[in] n: variable name
     * \param[in] i: component index, or -1 for a scalar variable
     * \param[in] v: value
     * \param[in] lb: lower bound
     * \param[in] ub: upper bound
     * \param[in] p: out of bounds policy
     */
    template <typename T>
    TFEL_HOST_DEVICE static bool treatOutOfBounds(
        [[maybe_unused]] const char* const n,
        [[maybe_unused]] const int i,
        [[maybe_unused]] const T v,
        [[maybe_unused]] const T lb,
        [[maybe_unused]] const T ub,
        const OutOfBoundsPolicy p) {
      if (p == None) {
        return true;
      }
#ifdef TFEL_DEVICE_CODE
      return p != Strict;
#else  /* TFEL_DEVICE_CODE */
      const auto name = getVariableName(n, i);
      if (p == Strict) {
        throwOutOfBoundsException(name, convert_to_string(v),
                                  convert_to_string(lb), convert_to_string(ub));
      }
      displayOutOfBoundsWarning(name, convert_to_string(v),
                                convert_to_string(lb), convert_to_string(ub));
      return true;
#endif /* TFEL_DEVICE_CODE */
    }
  };  // end of struct BoundsCheckBase

  /*!
   * \brief checks of the bounds of scalars and of the components of
   * symmetric tensors
   * \tparam N: space dimension
   */
  template <unsigned short N>
  struct BoundsCheck : public BoundsCheckBase {
    static_assert((N == 1u) || (N == 2u) || (N == 3u));
    using BoundsCheckBase::lowerAndUpperBoundsChecks;
    using BoundsCheckBase::lowerBoundCheck;
    using BoundsCheckBase::upperBoundCheck;
    template <typename T>
    TFEL_HOST_DEVICE static bool lowerBoundCheck(
        const char* const name,
        const tfel::math::stensor<N, T>& s,
        const tfel::math::base_type<T> lBound,
        const OutOfBoundsPolicy p = Strict) {
      auto r = true;
      for (unsigned short i = 0; i != s.size(); ++i) {
        const auto v = tfel::math::base_type_cast(s(i));
        if (v < lBound) {
          r = treatOutOfLowerBound(name, static_cast<int>(i), v, lBound, p) &&
              r;
        }
      }
      return r;
    }
    template <typename T>
    TFEL_HOST_DEVICE static bool upperBoundCheck(
        const char* const name,
        const tfel::math::stensor<N, T>& s,
        const tfel::math::base_type<T> uBound,
        const OutOfBoundsPolicy p = Strict) {
      auto r = true;
      for (unsigned short i = 0; i != s.size(); ++i) {
        const auto v = tfel::math::base_type_cast(s(i));
        if (v > uBound) {
          r = treatOutOfUpperBound(name, static_cast<int>(i), v, uBound, p) &&
              r;
        }
      }
      return r;
    }
    template <typename T>
    TFEL_HOST_DEVICE static bool lowerAndUpperBoundsChecks(
        const char* const name,
        const tfel::math::stensor<N, T>& s,
        const tfel::math::base_type<T> lBound,
        const tfel::math::base_type<T> uBound,
        const OutOfBoundsPolicy p = Strict) {
      auto r = true;
      for (unsigned short i = 0; i != s.size(); ++i) {
        const auto v = tfel::math::base_type_cast(s(i));
        if ((v < lBound) || (v > uBound)) {
          r = treatOutOfBounds(name, static_cast<int>(i), v, lBound, uBound,
                               p) &&
              r;
        }
      }
      return r;
    }
  };  // end of struct BoundsCheck

}  // end of namespace tfel::material

#endif /* LIB_TFEL_BOUNDSCHECK_HXX */
