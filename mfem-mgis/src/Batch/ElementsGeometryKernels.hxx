/*!
 * \file   src/Batch/ElementsGeometryKernels.hxx
 * \brief  This file declares the kernel computing the jacobians of the
 * element transformations of a mesh without nodes (internal header)
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_ELEMENTSGEOMETRYKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_ELEMENTSGEOMETRYKERNELS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief arguments of the kernel computing the jacobians of the element
   * transformations of a mesh without nodes from the coordinates of the
   * vertices
   */
  struct JacobiansKernelArguments {
    /*!
     * \brief vertices of the elements (`node + nd * element`), in the order
     * of the degrees of freedom of the finite element of the transformation
     */
    const int* element_vertices = nullptr;
    /*!
     * \brief coordinates of the vertices (`vertex + nv * component`, layout
     * of `mfem::Mesh::GetVertices`)
     */
    const real* vertices = nullptr;
    /*!
     * \brief gradients of the shape functions of the finite element of the
     * transformation on the reference element (`nq x dim x nd`, layout of
     * `mfem::DofToQuad` in the `FULL` mode)
     */
    const real* shape_functions_gradients = nullptr;
    /*!
     * \brief jacobians (`nq x dim x dim x ne`, layout of
     * `mfem::GeometricFactors::J`)
     */
    real* jacobians = nullptr;
    //! \brief determinants of the jacobians (`nq x ne`)
    real* determinants = nullptr;
    //! \brief number of elements
    int number_of_elements = 0;
    //! \brief number of vertices
    int number_of_vertices = 0;
    //! \brief number of vertices per element
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief dimension of the space (1, 2 or 3)
    int dim = 0;
  };

  /*!
   * \brief launches the kernel computing the jacobians of the element
   * transformations of a mesh without nodes and their determinants: one
   * thread per point of an element
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device (if a device is
   * enabled), and the pointers of the arguments are located on the device
   */
  void launchJacobiansKernel(const JacobiansKernelArguments&, const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_ELEMENTSGEOMETRYKERNELS_HXX */
