/*!
 * \file   include/MFEMMGIS/PostProcessingsComputation.hxx
 * \brief  This file declares the counter of the computations made by
 * kernels for the post-processings
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_POSTPROCESSINGSCOMPUTATION_HXX
#define LIB_MFEM_MGIS_POSTPROCESSINGSCOMPUTATION_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \return the number of computations on functions defined at the
   * integration points made by kernels since the start of the program:
   * integrals and measures (`computeIntegral`, `computeMeasure`), L2
   * projections and implicit gradient regularizations, values at the nodes
   * (`updateGridFunction`). The computations are made by kernels, if the
   * kernels support them, when the option `--execution` places the
   * post-processings on the GPU (see `ExecutionSettings`; the environment
   * variable `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION` is then ignored, with
   * a warning), or, without this option, when this environment variable is
   * `Kernels`; otherwise, they are made element by element. This counter
   * makes the choice observable, for the tests.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type
  getNumberOfPostProcessingsComputationsByKernels() noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_POSTPROCESSINGSCOMPUTATION_HXX */
