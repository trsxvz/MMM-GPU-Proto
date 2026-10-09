/*!
 * \file   src/Batch/StiffnessMatricesKernels.hxx
 * \brief  This file declares the kernels computing the stiffness matrices
 * of the elements of the standard mechanical behaviour integrators
 * (internal header)
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_STIFFNESSMATRICESKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_STIFFNESSMATRICESKERNELS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"

namespace mfem_mgis::batch::internals {

  //! \brief arguments of the kernels computing the stiffness matrices
  struct StiffnessMatricesKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    /*!
     * \brief gradients of the shape functions on the reference element
     * (`nq x dim x nd`)
     */
    const real* shape_functions_gradients = nullptr;
    //! \brief jacobians of the element transformations (`nq x dim x dim x ne`)
    const real* jacobians = nullptr;
    //! \brief determinants of the jacobians (`nq x ne`)
    const real* determinants = nullptr;
    //! \brief weights of the integration points on the reference element
    const real* weights = nullptr;
    //! \brief tangent operator blocks
    const real* tangent_operator = nullptr;
    //! \brief stride between two points of the tangent operator blocks
    long long point_stride = 0;
    //! \brief stride between two components of the tangent operator blocks
    long long component_stride = 0;
    //! \brief matrices element by element
    real* matrices = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief number of components of the gradients
    int gradients_size = 0;
  };

  //! \brief destination of the kernel: assembled sparse matrix
  struct SparseMatrixOutput {
    //! \brief number of rows
    int number_of_rows = 0;
    //! \brief offsets of the rows
    const int* I = nullptr;
    //! \brief columns
    const int* J = nullptr;
    //! \brief values
    real* values = nullptr;
    //! \brief degrees of freedom of the elements
    const int* element_dofs = nullptr;
    /*!
     * \brief number of entries which are not in the pattern of the matrix,
     * or whose row is not in the matrix
     */
    int* missing = nullptr;
  };

  /*!
   * \brief launches the kernel adding the stiffness matrices of the
   * elements of a material to an assembled sparse matrix
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments (`matrices` is not used)
   * \param[in] o: assembled matrix
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchStiffnessMatricesKernel(const GradientsKind,
                                     const StiffnessMatricesKernelArguments&,
                                     const SparseMatrixOutput&,
                                     const bool);
  /*!
   * \brief sets values to zero
   * \param[out] v: values
   * \param[in] n: number of values
   * \param[in] d: if true, the kernel runs on the device
   */
  void zeroValues(real* const, const int, const bool);

  /*!
   * \brief launches the kernel adding the stiffness matrices of the
   * elements of a material
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device (if a device is
   * enabled), and the pointers of the arguments are located on the device
   */
  void launchStiffnessMatricesKernel(const GradientsKind,
                                     const StiffnessMatricesKernelArguments&,
                                     const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_STIFFNESSMATRICESKERNELS_HXX */
