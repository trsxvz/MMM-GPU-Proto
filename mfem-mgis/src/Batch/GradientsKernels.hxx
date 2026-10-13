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
    /*!
     * \brief derivatives of the unknowns with respect to the reference
     * coordinates (see `UnknownsDerivatives`)
     */
    const real* derivatives = nullptr;
    /*!
     * \brief values of the unknowns at the points of all the elements
     * (layout of `UnknownsDerivatives::point_values`), required if the
     * gradients depend on the values of the unknowns (see
     * `gradientsRequireUnknownsValues`)
     */
    const real* unknowns_values = nullptr;
    /*!
     * \brief inverses of the jacobians of the element transformations (`nq x
     * dim x dim x ne`, see `ElementsJacobians::invJ`)
     */
    const real* inverse_jacobians = nullptr;
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
   * \brief arguments of the kernel computing the derivatives of the unknowns
   * with respect to the reference coordinates (see `UnknownsDerivatives`)
   */
  struct UnknownsDerivativesKernelArguments {
    /*!
     * \brief values of the unknowns element by element, in the layout of the
     * element restriction of MFEM (`node + nd * (component + vdim *
     * element)`)
     */
    const real* element_values = nullptr;
    /*!
     * \brief gradients of the shape functions on the reference element
     * (`nq x dim x nd`, layout of `mfem::DofToQuad` in the `FULL` mode)
     */
    const real* shape_functions_gradients = nullptr;
    /*!
     * \brief values of the shape functions on the reference element
     * (`nq x nd`, layout of `mfem::DofToQuad` in the `FULL` mode), required
     * if the values at the points are computed
     */
    const real* shape_functions = nullptr;
    //! \brief derivatives (layout of `UnknownsDerivatives::values`)
    real* derivatives = nullptr;
    /*!
     * \brief values at the points (layout of
     * `UnknownsDerivatives::point_values`), not computed if null
     */
    real* point_values = nullptr;
    //! \brief number of elements
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief number of components of the unknowns
    int vdim = 0;
    //! \brief dimension of the space
    int dim = 0;
  };

  /*!
   * \brief arguments of the kernel gathering the values of the unknowns at
   * the integration points of a material (see `getUnknownsValues`)
   */
  struct UnknownsValuesKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    //! \brief values at the points of all the elements
    const real* point_values = nullptr;
    //! \brief values at the points of the material
    real* values = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief number of components of the unknowns
    int vdim = 0;
  };

  /*!
   * \brief launches the kernel gathering the values of the unknowns at the
   * integration points of a material: one thread per point
   * \param[in] a: arguments
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchUnknownsValuesKernel(const UnknownsValuesKernelArguments& a,
                                  const bool use_device);

  /*!
   * \brief launches the kernel computing the derivatives of the unknowns
   * with respect to the reference coordinates at all the points of all the
   * elements, and their values if required: one thread per point of an
   * element
   * \param[in] a: arguments
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchUnknownsDerivativesKernel(
      const UnknownsDerivativesKernelArguments& a, const bool use_device);

  /*!
   * \brief launches the kernel computing the gradients of all the points of
   * a material
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchGradientsKernel(const GradientsKind k,
                             const GradientsKernelArguments& a,
                             const bool use_device);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_GRADIENTSKERNELS_HXX */
