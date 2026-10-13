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
   * at the nodes) are made by kernels, as given by the execution settings if
   * the option `--execution` is given (`Kernels` if they place the
   * post-processings on the GPU; the environment variable
   * `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION` is then ignored, with a
   * warning), by this environment variable otherwise (`ElementByElement` or
   * `Kernels`, `ElementByElement` if it is not defined). An invalid value of
   * the variable is an error (exception).
   */
  [[nodiscard]] bool usePostProcessingsKernels(attributes::MayThrow);
  /*!
   * \brief counts a computation made by kernels (see
   * `getNumberOfPostProcessingsComputationsByKernels`)
   */
  void countPostProcessingComputationByKernels() noexcept;

}  // end of namespace mfem_mgis::internals

#endif /* LIB_MFEM_MGIS_SRC_POSTPROCESSINGSCOMPUTATION_HXX */
