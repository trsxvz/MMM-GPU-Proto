/*!
 * \file   include/TFEL/Utilities/ConstantObject.hxx
 * \brief  This file declares the `getConstantObject` function, which gives
 * access to an object built at compile time in host and device code.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_UTILITIES_CONSTANTOBJECT_HXX
#define LIB_TFEL_UTILITIES_CONSTANTOBJECT_HXX

#include "TFEL/Config/TFELConfig.hxx"

namespace tfel::utilities::internals {

  //! \brief constant object used in host code
  template <typename T>
  inline constexpr T constant_object{};

#if (defined __CUDACC__) || (defined __HIP__)
  /*!
   * \brief constant object used in device code by the `CUDA` and `HIP`
   * compilers. Function scope static variables can't be used: their device
   * copy is not initialized by `nvcc` if a constructor is called.
   */
  template <typename T>
  __device__ constexpr T device_constant_object{};
#endif /* (defined __CUDACC__) || (defined __HIP__) */

}  // end of namespace tfel::utilities::internals

namespace tfel::utilities {

  /*!
   * \return a reference to an object of type `T` built at compile time by
   * its default constructor, which can be used in host and device code.
   *
   * This function returns a copy of this object when compiling for GPUs
   * with `nvc++ -stdpar=gpu`, `nvc++ -mp=gpu` or `nvc++ -acc=gpu`, which do
   * not allow device code to access global variables (unless the unified
   * memory mode is used). The returned value shall thus be bound to a
   * constant reference, for instance:
   *
   * \code{.cpp}
   * const auto& o = getConstantObject<T>();
   * \endcode
   *
   * \tparam T: type of the object, which must be default constructible at
   * compile time
   */
#if (defined __NVCOMPILER_STDPAR_GPU) || (defined __NVCOMPILER_OPENMP_GPU) || \
    (defined __NVCOMPILER_OPENACC_GPU)
  template <typename T>
  TFEL_HOST_DEVICE constexpr T getConstantObject() noexcept {
    return T{};
  }  // end of getConstantObject
#else /* (defined __NVCOMPILER_STDPAR_GPU) || ... */
  template <typename T>
  TFEL_HOST_DEVICE constexpr const T& getConstantObject() noexcept {
#if (defined __CUDA_ARCH__) || (defined __HIP_DEVICE_COMPILE__)
    return internals::device_constant_object<T>;
#else
    return internals::constant_object<T>;
#endif
  }  // end of getConstantObject
#endif /* (defined __NVCOMPILER_STDPAR_GPU) || ... */

}  // end of namespace tfel::utilities

#endif /* LIB_TFEL_UTILITIES_CONSTANTOBJECT_HXX */
