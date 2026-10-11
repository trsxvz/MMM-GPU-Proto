/*!
 * \file   src/Batch/InnerForcesKernels.hxx
 * \brief  This file declares the kernels computing the inner forces of the
 * standard behaviour integrators (mechanics, heat transfer, micromorphic
 * damage) (internal header)
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_INNERFORCESKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_INNERFORCESKERNELS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"

namespace mfem_mgis::batch::internals {

  //! \brief arguments of the kernels computing the inner forces
  struct InnerForcesKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    /*!
     * \brief gradients of the shape functions on the reference element
     * (`nq x dim x nd`)
     */
    const real* shape_functions_gradients = nullptr;
    /*!
     * \brief values of the shape functions on the reference element (`nq x
     * nd`), used if the gradients are the values of the unknowns
     */
    const real* shape_functions = nullptr;
    //! \brief jacobians of the element transformations (`nq x dim x dim x ne`)
    const real* jacobians = nullptr;
    //! \brief determinants of the jacobians (`nq x ne`)
    const real* determinants = nullptr;
    //! \brief weights of the integration points on the reference element
    const real* weights = nullptr;
    //! \brief thermodynamic forces
    const real* thermodynamic_forces = nullptr;
    //! \brief stride between two points of the thermodynamic forces
    long long point_stride = 0;
    //! \brief stride between two components of the thermodynamic forces
    long long component_stride = 0;
    /*!
     * \brief values element by element (`nd x vdim x ne`, `vdim` being the
     * number of components of the unknowns)
     */
    real* values = nullptr;
    /*!
     * \brief if true, the values are stored by element of the material
     * (`nd x vdim x number_of_elements`) instead of by element of the mesh
     */
    bool by_material_element = false;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief dimension of the space
    int dim = 0;
  };

  /*!
   * \brief launches the kernel adding the inner forces of the elements of
   * a material
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device (if a device is
   * enabled), and the pointers of the arguments are located on the device
   */
  void launchInnerForcesKernel(const GradientsKind,
                               const InnerForcesKernelArguments&,
                               const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_INNERFORCESKERNELS_HXX */
