/*!
 * \file  include/TFEL/Math/Tensor/TensorViewFromStensor.hxx
 * \brief The TensorViewFromStensor describes how a stensor can be interpreted
 * as a tensor \author Thomas Helfer \brief 04 juil. 2013 \copyright Copyright
 * (C) 2006-2018 CEA/DEN, EDF R&D. All rights reserved. This project is publicly
 * released under either the GNU GPL Licence with linking exception or the
 * CECILL-A licence. A copy of thoses licences are delivered with the sources of
 * TFEL. CEA or EDF may also distribute this project under specific licensing
 * conditions.
 */

#ifndef LIB_TFEL_MATH_TENSORVIEWFROMSTENSOR_HXX
#define LIB_TFEL_MATH_TENSORVIEWFROMSTENSOR_HXX

#include "TFEL/Math/General/EmptyRunTimeProperties.hxx"
#include "TFEL/Math/General/ConstExprMathFunctions.hxx"
#include "TFEL/Math/Tensor/TensorConcept.hxx"
#include "TFEL/Math/Stensor/StensorConcept.hxx"
#include "TFEL/Math/Forward/tensor.hxx"

namespace tfel::math {

  //! \brief tensor view from a symmetric tensor
  template <StensorConcept>
  struct TensorViewFromStensorExpr {
  };  // end of struct TensorViewFromStensorExpr

  template <unsigned short N, typename ValueType, StensorConcept T>
  struct Expr<tensor<N, ValueType>, TensorViewFromStensorExpr<T>>
      : public TensorConceptBase<
            Expr<tensor<N, ValueType>, TensorViewFromStensorExpr<T>>>,
        public ExprBase {
    typedef EmptyRunTimeProperties RunTimeProperties;

    [[nodiscard]] constexpr auto getRunTimeProperties() const noexcept {
      return RunTimeProperties();
    }

    // NOLINTBEGIN(google-explicit-constructor)
    TFEL_HOST_DEVICE constexpr Expr(T s_) noexcept : s(s_) {}  // end of Expr
    // NOLINTEND(google-explicit-constructor)

    TFEL_HOST_DEVICE constexpr Expr(const Expr&) = default;
    TFEL_HOST_DEVICE constexpr Expr(Expr&&) = default;
    TFEL_HOST_DEVICE constexpr Expr& operator=(const Expr&) = default;
    TFEL_HOST_DEVICE constexpr Expr& operator=(Expr&&) = default;

    [[nodiscard]] TFEL_HOST_DEVICE constexpr ValueType operator[](
        const unsigned short i) const {
      return this->operator()(i);
    }

    [[nodiscard]] TFEL_HOST_DEVICE constexpr ValueType operator()(
        const unsigned short i) const {
      // component of the symmetric tensor associated with each component of
      // the tensor, the off-diagonal components being divided by the square
      // root of two. A switch statement is avoided: nvc++ 25.11 (-mp=gpu,
      // -O2 and above) fails on it when this operator is inlined in a product
      // of tensors (internal compiler error).
      constexpr unsigned short indices[9] = {0, 1, 2, 3, 3, 4, 4, 5, 5};
      if (i >= TensorDimeToSize<N>::value) {
        return internals::reportInvalidTensorIndex<TensorInvalidIndexException,
                                                   ValueType>();
      }
      const auto v = this->s(indices[i]);
      return (i < 3) ? v : v * Cste<ValueType>::isqrt2;
    }  // end of operator() const

   protected:
    //! The stensor object
    ArgumentStorage<T> s;

   private:
    using stype = std::decay_t<T>;
    using traits = MathObjectTraits<stype>;

    static_assert((N == 1u) || (N == 2u) || (N == 3u));
    static_assert(traits::dime == N);
    static_assert(std::is_same_v<ValueType, typename traits::NumType>);
  };  // end of struct Expr

  template <typename T>
  struct TensorViewFromStensor {
    using stype = std::decay_t<T>;
    using traits = MathObjectTraits<stype>;
    static const unsigned short N = traits::dime;
    typedef Expr<tensor<N, typename traits::NumType>,
                 TensorViewFromStensorExpr<T>>
        type;
  };  // end of struct TensorViewFromStensor

}  // end of namespace tfel::math

#endif /* LIB_TFEL_MATH_TENSORVIEWFROMSTENSOR_HXX */
