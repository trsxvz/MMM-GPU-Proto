/*!
 * \file   MFEMMGIS/MPI.hxx
 * \brief  This file declares some MPI utilities
 * \author Thomas Helfer
 * \date   06/02/2026
 */

#ifndef LIB_MFEM_MGIS_MPI_HXX
#define LIB_MFEM_MGIS_MPI_HXX

#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <concepts>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  //! \brief MPI data type associated with the type T
  template <typename T>
  inline const auto mpi_type = []() {
    if constexpr (std::is_enum_v<T>) {
      return mpi_type<std::underlying_type_t<T>>;
    } else if constexpr (std::same_as<T, bool>) {
      return MPI_CXX_BOOL;
    } else if constexpr (std::same_as<T, int>) {
      return MPI_INT;
    } else if constexpr (std::same_as<T, char>) {
      return MPI_CHAR;
    } else if constexpr (std::same_as<T, signed char>) {
      return MPI_SIGNED_CHAR;
    } else if constexpr (std::same_as<T, unsigned char>) {
      return MPI_UNSIGNED_CHAR;
    } else if constexpr (std::same_as<T, short>) {
      return MPI_SHORT;
    } else if constexpr (std::same_as<T, unsigned short>) {
      return MPI_UNSIGNED_SHORT;
    } else if constexpr (std::same_as<T, unsigned int>) {
      return MPI_UNSIGNED;
    } else if constexpr (std::same_as<T, long>) {
      return MPI_LONG;
    } else if constexpr (std::same_as<T, unsigned long>) {
      return MPI_UNSIGNED_LONG;
    } else if constexpr (std::same_as<T, long long>) {
      return MPI_LONG_LONG_INT;
    } else if constexpr (std::same_as<T, unsigned long long>) {
      return MPI_UNSIGNED_LONG_LONG;
    } else if constexpr (std::same_as<T, float>) {
      return MPI_FLOAT;
    } else if constexpr (std::same_as<T, double>) {
      return MPI_DOUBLE;
    } else if constexpr (std::same_as<T, long double>) {
      return MPI_LONG_DOUBLE;
    } else if constexpr (std::same_as<T, int8_t>) {
      return MPI_INT8_T;
    } else if constexpr (std::same_as<T, int16_t>) {
      return MPI_INT16_T;
    } else if constexpr (std::same_as<T, int32_t>) {
      return MPI_INT32_T;
    } else if constexpr (std::same_as<T, int64_t>) {
      return MPI_INT64_T;
    } else if constexpr (std::same_as<T, uint8_t>) {
      return MPI_UINT8_T;
    } else if constexpr (std::same_as<T, uint16_t>) {
      return MPI_UINT16_T;
    } else if constexpr (std::same_as<T, uint32_t>) {
      return MPI_UINT32_T;
    } else if constexpr (std::same_as<T, uint64_t>) {
      return MPI_UINT64_T;
    } else {
      []<bool flag = false>() {
        static_assert(flag, "no MPI type associated with the given type");
      }
      ();
    }
  }();

#endif /* MFEM_USE_MPI */

  /*!
   * \brief a simple reduction for boolean values
   *
   * \param[in] m: mesh discretization
   * \param[in] b: boolean value in the current process
   * \return true if the value is true on all processes
   *
   * \note if the computations are sequential, no MPI call is made
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isTrueOnAllProcesses(
      const MeshDiscretization& m, const bool b) noexcept;

  /*!
   * \brief check the validity of a value on all processes
   * \return if the given value is valid on all processes
   *
   * \param[in] m: mesh discretization
   * \param[in] v: value in the current process
   *
   * \note if the computations are sequential, no MPI call is made
   */
  template <typename T>
  [[nodiscard]] bool isValidOnAllProcesses(const MeshDiscretization& m,
                                           const T& v) noexcept;
  /*!
   * \brief exception thrown on all the processes when a computation failed
   * on one of them (see `raiseOnFailureOnAnyProcess`): the processes can
   * handle it and continue together, unlike an exception thrown by a single
   * process
   */
  struct MFEM_MGIS_EXPORT CollectiveFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
    //! \brief destructor
    ~CollectiveFailure() noexcept override;
  };
  /*!
   * \brief throws a `CollectiveFailure` exception on all the processes if a
   * computation failed on one of them, before the collective operations
   * following it
   *
   * \param[in] m: mesh discretization
   * \param[in, out] ctx: execution context of the computation
   * \param[in] success: result of the computation in the current process
   * \param[in] what: description of the computation
   *
   * \note if the computations are sequential, no MPI call is made
   */
  MFEM_MGIS_EXPORT void raiseOnFailureOnAnyProcess(const MeshDiscretization&,
                                                   Context&,
                                                   const bool,
                                                   std::string_view);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/MPI.ixx"

#endif /* LIB_MFEM_MGIS_MPI_HXX */
