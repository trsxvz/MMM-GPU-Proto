/*!
 * \file   src/AmgXKernels.hxx
 * \brief  This file declares the scaling and the conversions of the vectors
 * exchanged with the solvers of AmgX (internal header)
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_AMGXKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_AMGXKERNELS_HXX

#include "mfem/general/array.hpp"
#include "mfem/linalg/vector.hpp"

namespace mfem_mgis::internals {

  /*!
   * \return the maximum of the absolute values of the components of a
   * vector (zero for an empty vector, the components which are not a
   * number being ignored), computed on the device if the vector uses it
   * \param[in] v: vector
   * \param[in, out] w: workspace of the reduction
   */
  double getMaximumAbsoluteValue(const mfem::Vector&, mfem::Array<double>&);
  /*!
   * \brief multiplies a vector by a factor, on the device if the vector
   * uses it
   * \param[out] d: result, resized
   * \param[in] s: vector
   * \param[in] f: factor
   */
  void scale(mfem::Vector&, const mfem::Vector&, const double);
  /*!
   * \brief converts a vector in double precision, multiplied by a factor,
   * to single precision, on the device if the vector uses it
   * \param[out] d: values in single precision, resized
   * \param[in] s: vector
   * \param[in] f: factor
   */
  void convertToSinglePrecision(mfem::Array<float>&,
                                const mfem::Vector&,
                                const double = 1);
  /*!
   * \brief converts values in single precision, multiplied by a factor, to
   * a vector in double precision, on the device if the vector uses it
   * \param[out] d: vector, of the size of the values
   * \param[in] s: values in single precision
   * \param[in] f: factor
   */
  void convertToDoublePrecision(mfem::Vector&,
                                const mfem::Array<float>&,
                                const double = 1);

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_AMGXKERNELS_HXX */
