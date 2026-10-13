/*!
 * \file   src/Batch/ElementsGeometryKernels.hxx
 * \brief  This file declares the kernel computing the determinants and the
 * inverses of the jacobians of the
 * element transformations (internal header)
 * \date   10/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_ELEMENTSGEOMETRYKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_ELEMENTSGEOMETRYKERNELS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief arguments of the kernel computing the determinants and the
   * inverses of the jacobians of the element
   * transformations from the coordinates of their nodes (the nodes of the
   * mesh, or its vertices)
   */
  struct JacobiansKernelArguments {
    /*!
     * \brief nodes of the elements (`node + nodes_per_element * element`),
     * in the order of the degrees of freedom of the finite element of the
     * transformation
     */
    const int* element_nodes = nullptr;
    /*!
     * \brief coordinates of the nodes: the component `c` of the node `n` is
     * `coordinates[n * node_stride + c * component_stride]`
     */
    const real* coordinates = nullptr;
    //! \brief stride between two nodes in the coordinates
    int node_stride = 0;
    //! \brief stride between two components in the coordinates
    int component_stride = 0;
    /*!
     * \brief gradients of the shape functions of the finite element of the
     * transformation on the reference element (`nq x dim x
     * nodes_per_element`, layout of `mfem::DofToQuad` in the `FULL` mode)
     */
    const real* shape_functions_gradients = nullptr;
    /*!
     * \brief inverses of the jacobians (`nq x dim x dim x ne`, layout of
     * `mfem::GeometricFactors::J`)
     */
    real* inverse_jacobians = nullptr;
    //! \brief determinants of the jacobians (`nq x ne`)
    real* determinants = nullptr;
    //! \brief number of elements
    int number_of_elements = 0;
    //! \brief number of nodes per element
    int nodes_per_element = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief dimension of the space (1, 2 or 3)
    int dim = 0;
  };

  /*!
   * \brief launches the kernel computing the jacobians of the element
   * transformations, their determinants and their inverses (see
   * `invertJacobian`, the jacobians themselves are not stored): one thread
   * per point of an element
   * \param[in] a: arguments
   * \param[in] use_device: if true, the kernel runs on the device (if a device
   * is enabled), and the pointers of the arguments are located on the device
   */
  void launchJacobiansKernel(const JacobiansKernelArguments& a,
                             const bool use_device);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_ELEMENTSGEOMETRYKERNELS_HXX */
