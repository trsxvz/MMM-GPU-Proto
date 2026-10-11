/*!
 * \file   src/PostProcessingsComputation.hxx
 * \brief  This file declares the choice of the computation of the
 * post-processings (internal header)
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_POSTPROCESSINGSCOMPUTATION_HXX
#define LIB_MFEM_MGIS_SRC_POSTPROCESSINGSCOMPUTATION_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return true if the computations on functions defined at the integration
   * points used by the post-processings (integrals, L2 projections, values
   * at the nodes) are made by kernels, as given by the environment variable
   * `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION`: `ElementByElement` (default)
   * or `Kernels`. An invalid value is an error.
   */
  [[nodiscard]] bool usePostProcessingsKernels(attributes::MayThrow);
  /*!
   * \brief counts a computation made by kernels (see
   * `getNumberOfPostProcessingsComputationsByKernels`)
   */
  void countPostProcessingComputationByKernels() noexcept;

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_POSTPROCESSINGSCOMPUTATION_HXX */
