/*!
 * \file   src/Batch/GradientsKernels.hxx
 * \brief  This file declares the kernels computing the gradients of the
 * behaviours (internal header)
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_GRADIENTSKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_GRADIENTSKERNELS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"

namespace mfem_mgis::batch::internals {

  //! \brief arguments of the kernels computing the gradients
  struct GradientsKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    //! \brief derivatives of the unknowns (see `UnknownsDerivatives`)
    const real* derivatives = nullptr;
    //! \brief gradients of the material
    real* gradients = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief number of components of the unknowns
    int vdim = 0;
    //! \brief dimension of the space
    int dim = 0;
    //! \brief number of integration points of the material
    int number_of_points = 0;
    //! \brief number of components of the gradients
    int gradients_size = 0;
    //! \brief layout of the gradients: true if stored point by point
    bool point_by_point = true;
    //! \brief macroscopic gradients
    real macroscopic_gradients[9] = {};
  };

  /*!
   * \brief launches the kernel computing the gradients of all the points of
   * a material
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device (if a device is
   * enabled), and the pointers of the arguments are located on the device
   */
  void launchGradientsKernel(const GradientsKind,
                             const GradientsKernelArguments&,
                             const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_GRADIENTSKERNELS_HXX */
