/*!
 * \file   src/Batch/FBarKernels.hxx
 * \brief  This file declares the kernels of the FBar formulation (internal
 * header, see `MFEMMGIS/Batch/FBar.hxx`)
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_FBARKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_FBARKERNELS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "StiffnessMatricesKernels.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief arguments of the kernel computing the data of the FBar
   * formulation at the center of the elements of a material
   */
  struct FBarCentersKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    /*!
     * \brief derivatives of the unknowns at the center of all the elements
     * of the mesh (see `UnknownsDerivatives`), null for null unknowns
     */
    const real* derivatives = nullptr;
    //! \brief jacobians at the center of all the elements of the mesh
    const real* jacobians = nullptr;
    //! \brief determinants of the jacobians at the center of the elements
    const real* determinants = nullptr;
    /*!
     * \brief gradients of the shape functions at the center of the
     * reference element (`dim x nd`)
     */
    const real* shape_functions_gradients = nullptr;
    //! \brief determinant `J0` at the center of the elements of the material
    real* centers = nullptr;
    /*!
     * \brief vectors `b` of the elements of the material (`node + nd *
     * (component + vdim * k)` for the `k`-th element of the material)
     */
    real* b = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief macroscopic gradients
    real macroscopic_gradients[9] = {};
  };

  /*!
   * \brief arguments of the kernel computing the coefficients `α` and the
   * gradients `F̄` at the integration points of a material
   */
  struct FBarGradientsKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    //! \brief data at the center of the elements of the material
    const real* centers = nullptr;
    /*!
     * \brief gradients of the material (deformation gradients computed by
     * the kernel of the gradients, replaced by `F̄`), null for null unknowns
     */
    real* gradients = nullptr;
    //! \brief deformation gradients `F`, point by point
    real* deformation_gradients = nullptr;
    //! \brief coefficients `α`
    real* coefficients = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of integration points per element
    int nq = 0;
    //! \brief number of integration points of the material
    int number_of_points = 0;
    //! \brief layout of the gradients: true if stored point by point
    bool point_by_point = true;
    //! \brief macroscopic gradients (for null unknowns)
    real macroscopic_gradients[9] = {};
  };

  /*!
   * \brief arguments of the kernel computing the tangent operator `K'` and
   * the pseudo-stresses `α / (ie J0) K F` at the integration points of a
   * material
   */
  struct FBarTangentOperatorKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief offsets of the integration points of the elements
    const int* offsets = nullptr;
    //! \brief data at the center of the elements of the material
    const real* centers = nullptr;
    //! \brief deformation gradients `F`, point by point
    const real* deformation_gradients = nullptr;
    //! \brief coefficients `α`
    const real* coefficients = nullptr;
    //! \brief tangent operator blocks of the behaviour
    const real* tangent_operator = nullptr;
    //! \brief stride between two points of the tangent operator blocks
    long long point_stride = 0;
    //! \brief stride between two components of the tangent operator blocks
    long long component_stride = 0;
    //! \brief tangent operator `K'`, point by point
    real* modified_tangent_operator = nullptr;
    //! \brief pseudo-stresses, point by point
    real* pseudo_stresses = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of integration points per element
    int nq = 0;
  };

  /*!
   * \brief arguments of the kernels adding the matrices of rank one `a ⊗ b`
   * of the elements of a material
   */
  struct FBarRankOneKernelArguments {
    //! \brief elements of the material
    const int* elements = nullptr;
    //! \brief vectors `a`, by element of the material
    const real* a = nullptr;
    //! \brief vectors `b`, by element of the material
    const real* b = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom of an element
    int n = 0;
  };

  /*!
   * \brief launches the kernel computing the data of the FBar formulation
   * at the center of the elements of a material: one thread per element
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarCentersKernel(const GradientsKind,
                               const FBarCentersKernelArguments&,
                               const bool);
  /*!
   * \brief launches the kernel computing the coefficients `α` and the
   * gradients `F̄`: one thread per integration point
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarGradientsKernel(const GradientsKind,
                                 const FBarGradientsKernelArguments&,
                                 const bool);
  /*!
   * \brief launches the kernel computing the tangent operator `K'` and the
   * pseudo-stresses: one thread per integration point
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarTangentOperatorKernel(
      const GradientsKind,
      const FBarTangentOperatorKernelArguments&,
      const bool);
  /*!
   * \brief adds the matrices `a ⊗ b` to matrices stored element by element,
   * column by column (see `addStiffnessMatrices`): one thread per entry
   * \param[in] a: arguments
   * \param[in, out] m: matrices of all the elements of the mesh
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarRankOneMatricesKernel(const FBarRankOneKernelArguments&,
                                       real* const,
                                       const bool);
  /*!
   * \brief adds the matrices `a ⊗ b` to an assembled sparse matrix: one
   * thread per row of an element
   * \param[in] a: arguments
   * \param[in] o: assembled matrix
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarRankOneMatricesKernel(const FBarRankOneKernelArguments&,
                                       const SparseMatrixOutput&,
                                       const bool);
  /*!
   * \brief adds the diagonals of the matrices `a ⊗ b`: one thread per
   * degree of freedom of an element
   * \param[in] a: arguments
   * \param[in, out] diagonals: values element by element
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarRankOneDiagonalsKernel(const FBarRankOneKernelArguments&,
                                        real* const,
                                        const bool);
  /*!
   * \brief adds the action of the matrices `a ⊗ b`, `a (b · x)`: one thread
   * per element
   * \param[in] a: arguments
   * \param[in] x: values element by element
   * \param[in, out] y: action element by element
   * \param[in] d: if true, the kernel runs on the device
   */
  void launchFBarRankOneActionKernel(const FBarRankOneKernelArguments&,
                                     const real* const,
                                     real* const,
                                     const bool);

}  // end of namespace mfem_mgis::batch::internals

#endif /* LIB_MFEM_MGIS_SRC_BATCH_FBARKERNELS_HXX */
