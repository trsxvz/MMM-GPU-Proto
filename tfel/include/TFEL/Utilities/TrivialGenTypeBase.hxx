/*!
 * \file   include/TFEL/Utilities/TrivialGenTypeBase.hxx
 * \brief  This file implements the `TrivialGenTypeBase` class.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_UTILITIES_TRIVIALGENTYPEBASE_HXX
#define LIB_TFEL_UTILITIES_TRIVIALGENTYPEBASE_HXX

#include <utility>
#include <type_traits>
#include "TFEL/Raise.hxx"
#include "TFEL/Config/TFELConfig.hxx"
#include "TFEL/Metaprogramming/TypeList.hxx"
#include "TFEL/Utilities/GenTypeBase.hxx"
#include "TFEL/Utilities/GenTypeCastError.hxx"

namespace tfel::utilities::internals {

  /*!
   * \brief a union storing a value of one of the given types.
   * \tparam Types: stored types
   */
  template <typename... Types>
  union TrivialGenTypeStorage;

  //! \brief explicit specialisation ending the recursion
  template <>
  union TrivialGenTypeStorage<> {};

  /*!
   * \brief partial specialisation storing a value of the first type or a
   * value of one of the other types.
   */
  template <typename Type, typename... Types>
  union TrivialGenTypeStorage<Type, Types...> {
    //! \brief default constructor: no value is stored, no member is active
    TFEL_HOST_DEVICE constexpr TrivialGenTypeStorage() noexcept {}
    /*!
     * \brief constructor of a value of the first type
     * \param[in] args: arguments of the constructor of the value
     */
    template <typename... Arguments>
    TFEL_HOST_DEVICE constexpr TrivialGenTypeStorage(
        std::integral_constant<unsigned short, 0>, Arguments&&... args)
        : value(std::forward<Arguments>(args)...) {}
    /*!
     * \brief constructor of a value of the `I`-th type
     * \param[in] args: arguments of the constructor of the value
     */
    template <unsigned short I, typename... Arguments>
    TFEL_HOST_DEVICE constexpr TrivialGenTypeStorage(
        std::integral_constant<unsigned short, I>, Arguments&&... args)
      requires(I != 0)
        : others(std::integral_constant<unsigned short, I - 1>{},
                 std::forward<Arguments>(args)...) {}
    //! \brief value of the first type
    Type value;
    //! \brief storage of the values of the other types
    TrivialGenTypeStorage<Types...> others;
  };

  /*!
   * \return the value of the `I`-th type stored in the given union
   * \param[in] s: storage
   */
  template <unsigned short I, typename Storage>
  TFEL_HOST_DEVICE constexpr auto& getTrivialGenTypeValue(Storage& s) noexcept {
    if constexpr (I == 0) {
      return s.value;
    } else {
      return getTrivialGenTypeValue<I - 1>(s.others);
    }
  }  // end of getTrivialGenTypeValue

  //! \brief build the storage associated with a list of types
  template <typename List, typename... Types>
  struct TrivialGenTypeStorageFromTypeList {
    //! \brief a simple alias
    using type = typename TrivialGenTypeStorageFromTypeList<
        typename List::Next,
        Types...,
        typename List::Current>::type;
  };

  //! \brief partial specialisation ending the recursion
  template <typename... Types>
  struct TrivialGenTypeStorageFromTypeList<tfel::meta::TLE, Types...> {
    static_assert((std::is_trivially_copyable_v<Types> && ...),
                  "all the types must be trivially copyable");
    static_assert((std::is_trivially_destructible_v<Types> && ...),
                  "all the types must be trivially destructible");
    //! \brief a simple alias
    using type = TrivialGenTypeStorage<Types...>;
  };

}  // end of namespace tfel::utilities::internals

namespace tfel::utilities {

  /*!
   * \brief a class holding a value of one of the types of the given list, or
   * no value.
   *
   * This class has the same interface as `GenTypeBase`, but is restricted to
   * trivially copyable types. Contrary to `GenTypeBase`, it does not rely on
   * `std::variant`, whose implementations are not supported in device code
   * by some compilers. This class is trivially copyable and can be used in
   * device code.
   *
   * \tparam List: list of the types
   */
  template <typename List>
  struct TrivialGenTypeBase
      : public internals::GenTypeSpecializedMethods<TrivialGenTypeBase<List>,
                                                    List> {
    static_assert(tfel::meta::TLElementsAreUnique<List>::cond);
    //! \brief a variable checking if the given type is valid
    template <typename T>
    static constexpr bool isValidType =
        tfel::meta::TLCountNbrOfT<std::decay_t<T>, List>::value == 1;
    //! \brief default constructor, no value is stored
    TrivialGenTypeBase() = default;
    //! \brief move constructor
    TrivialGenTypeBase(TrivialGenTypeBase&&) = default;
    //! \brief copy constructor
    TrivialGenTypeBase(const TrivialGenTypeBase&) = default;
    //! \brief constructor from a value
    template <typename T1>
    TFEL_HOST_DEVICE explicit(false) TrivialGenTypeBase(T1&& value)
      requires(isValidType<T1>)
    {
      this->set(std::forward<T1>(value));
    }
    //! \brief move assignement
    TrivialGenTypeBase& operator=(TrivialGenTypeBase&&) = default;
    //! \brief assignement operator
    TrivialGenTypeBase& operator=(const TrivialGenTypeBase&) = default;
    //! \brief assignement operator from a value
    template <typename T1>
    TFEL_HOST_DEVICE TrivialGenTypeBase& operator=(T1&& value)
      requires(isValidType<T1>)
    {
      this->set(std::forward<T1>(value));
      return *this;
    }
    //! \return true if no value is stored
    [[nodiscard]] TFEL_HOST_DEVICE bool empty() const noexcept {
      return this->type_index == 0;
    }
    /*!
     * \brief copy another object
     * \param[in] src: object to be copied
     */
    [[nodiscard]] bool copy(const TrivialGenTypeBase& src) noexcept {
      *this = src;
      return true;
    }
    /*!
     * \brief set the value
     * \param[in] src: value
     */
    template <typename T1>
    TFEL_HOST_DEVICE void set(T1&& src)
      requires(isValidType<T1>)
    {
      constexpr auto i = position<T1>;
      this->storage = Storage(std::integral_constant<unsigned short, i>{},
                              std::forward<T1>(src));
      this->type_index = i + 1;
    }
    //! \return true if a value of the given type is stored
    template <typename T1>
    [[nodiscard]] TFEL_HOST_DEVICE bool is() const noexcept
      requires(isValidType<T1>)
    {
      return this->type_index == position<T1> + 1;
    }
    /*!
     * \return the stored value
     * \note the type of the value is only checked on the host (an exception
     * is thrown if the value is not of the given type), except with the
     * compilers translating host and device code in a single pass, which
     * follow the rules of device code (see `TFEL_DEVICE_CODE`).
     */
    template <typename T1>
    [[nodiscard]] TFEL_HOST_DEVICE const T1& get() const
      requires(isValidType<T1>)
    {
#ifndef TFEL_DEVICE_CODE
      tfel::raise_if<GenTypeCastError>(!this->template is<T1>());
#endif /* TFEL_DEVICE_CODE */
      return internals::getTrivialGenTypeValue<position<T1>>(this->storage);
    }
    /*!
     * \return the stored value
     * \note the type of the value is only checked on the host (an exception
     * is thrown if the value is not of the given type), except with the
     * compilers translating host and device code in a single pass, which
     * follow the rules of device code (see `TFEL_DEVICE_CODE`).
     */
    template <typename T1>
    [[nodiscard]] TFEL_HOST_DEVICE T1& get()
      requires(isValidType<T1>)
    {
#ifndef TFEL_DEVICE_CODE
      tfel::raise_if<GenTypeCastError>(!this->template is<T1>());
#endif /* TFEL_DEVICE_CODE */
      return internals::getTrivialGenTypeValue<position<T1>>(this->storage);
    }
    /*!
     * \return the index of the type of the stored value, starting from one,
     * or zero if no value is stored
     */
    [[nodiscard]] TFEL_HOST_DEVICE unsigned short getTypeIndex()
        const noexcept {
      return this->type_index;
    }

   protected:
    //! \brief remove the stored value
    TFEL_HOST_DEVICE void clear() noexcept { this->type_index = 0; }
    //! \brief store a default constructed value of the given type
    template <typename T1>
    TFEL_HOST_DEVICE void set_uninitialised()
      requires(isValidType<T1>)
    {
      this->set(T1());
    }

   private:
    //! \brief position of the given type in the list
    template <typename T1>
    static constexpr auto position = static_cast<unsigned short>(
        tfel::meta::TLFindEltPos<std::decay_t<T1>, List>::value);
    //! \brief a simple alias
    using Storage =
        typename internals::TrivialGenTypeStorageFromTypeList<List>::type;
    //! \brief stored value
    Storage storage;
    //! \brief index of the type of the stored value, zero if empty
    unsigned short type_index = 0;
  };

}  // end of namespace tfel::utilities

#endif /* LIB_TFEL_UTILITIES_TRIVIALGENTYPEBASE_HXX */
