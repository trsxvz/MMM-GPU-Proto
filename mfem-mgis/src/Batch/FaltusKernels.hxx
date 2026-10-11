/*!
 * \file   src/Batch/FaltusKernels.hxx
 * \brief  This file declares the kernels of the regularization of Faltus et
 * al. (internal header, see `MFEMMGIS/Batch/Faltus.hxx`)
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_FALTUSKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_FALTUSKERNELS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "StiffnessMatricesKernels.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief arguments of the kernel computing the deformation gradients at
   * the center of the elements of a material
   */
  struct FaltusCentersKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    /*!
     * \brief derivatives of the unknowns at the center of all the elements
     * of the mesh (see `UnknownsDerivatives`), null for null unknowns
     */
    const real* derivatives = nullptr;
    /*!
     * \brief inverses of the jacobians at the center of all the elements of
     * the mesh (see `ElementsJacobians::invJ`)
     */
    const real* inverse_jacobians = nullptr;
    /*!
     * \brief deformation gradients `F0` at the center of the elements of
     * the material (`c + gs * k` for the component `c` of the `k`-th
     * element of the material)
     */
    real* centers = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief macroscopic gradients
    real macroscopic_gradients[9] = {};
  };

  /*!
   * \brief geometric data of the elements of a material used by the kernels
   * of the regularization: the operator `D = B - B0` at an integration
   * point, `B` and `B0` being the operators of the gradients at the point
   * and at the center of the element
   */
  struct FaltusGeometryKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    /*!
     * \brief gradients of the shape functions at the integration points of
     * the reference element (`nq x dim x nd`)
     */
    const real* shape_functions_gradients = nullptr;
    /*!
     * \brief inverses of the jacobians at the integration points of all the
     * elements (see `ElementsJacobians::invJ`)
     */
    const real* inverse_jacobians = nullptr;
    //! \brief determinants of the jacobians at the integration points
    const real* determinants = nullptr;
    //! \brief weights of the integration points
    const real* weights = nullptr;
    /*!
     * \brief gradients of the shape functions at the center of the
     * reference element (`dim x nd`)
     */
    const real* center_shape_functions_gradients = nullptr;
    /*!
     * \brief inverses of the jacobians at the center of all the elements
     * (see `ElementsJacobians::invJ`)
     */
    const real* center_inverse_jacobians = nullptr;
    //! \brief penalization coefficient `α`
    real alpha = 0;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
  };

  /*!
   * \brief arguments of the kernel adding the contribution of the
   * regularization to the inner forces of the elements of a material
   */
  struct FaltusInnerForcesKernelArguments {
    //! \brief geometric data
    FaltusGeometryKernelArguments geometry;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    //! \brief deformation gradients at the integration points
    const real* gradients = nullptr;
    //! \brief layout of the gradients: true if stored point by point
    bool point_by_point = true;
    //! \brief number of integration points of the material
    int number_of_points = 0;
    //! \brief deformation gradients `F0` at the center of the elements
    const real* centers = nullptr;
    /*!
     * \brief inner forces, element by element, in the layout of the element
     * restriction of MFEM (`node + nd * (component + vdim * e)`, `e` being
     * the element of the mesh)
     */
    real* values = nullptr;
  };

  /*!
   * \brief launches the kernel computing the deformation gradients at the
   * center of the elements of a material: one thread per element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFaltusCentersKernel(const GradientsKind,
                                 const FaltusCentersKernelArguments&,
                                 const bool);
  /*!
   * \brief launches the kernel adding `Σ_q w_q Dᵀ α (F_q - F0)` to the inner
   * forces of the elements of a material: one thread per element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFaltusInnerForcesKernel(const GradientsKind,
                                     const FaltusInnerForcesKernelArguments&,
                                     const bool);
  /*!
   * \brief launches the kernel adding `α Σ_q w_q Dᵀ D` to the stiffness
   * matrices of the elements of a material, stored element by element
   * (`i + n * (j + n * e)` for the row `i` and the column `j` of the
   * element `e` of the mesh): one thread per pair of nodes of an element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in, out] m: matrices element by element
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFaltusMatricesKernel(const GradientsKind,
                                  const FaltusGeometryKernelArguments&,
                                  real* const,
                                  const bool);
  /*!
   * \brief launches the kernel adding `α Σ_q w_q Dᵀ D` to an assembled
   * sparse matrix: one thread per node of an element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] o: assembled matrix
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFaltusMatricesKernel(const GradientsKind,
                                  const FaltusGeometryKernelArguments&,
                                  const SparseMatrixOutput&,
                                  const bool);
  /*!
   * \brief launches the kernel adding the diagonals of `α Σ_q w_q Dᵀ D` to
   * values element by element: one thread per node of an element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in, out] v: values element by element
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFaltusDiagonalsKernel(const GradientsKind,
                                   const FaltusGeometryKernelArguments&,
                                   real* const,
                                   const bool);
  /*!
   * \brief launches the kernel adding the action of `α Σ_q w_q Dᵀ D` on
   * values given element by element: one thread per element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] x: values element by element
   * \param[in, out] y: action, element by element
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFaltusActionKernel(const GradientsKind,
                                const FaltusGeometryKernelArguments&,
                                const real* const,
                                real* const,
                                const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_FALTUSKERNELS_HXX */
