/*!
 * \file   include/MFEMMGIS/Batch/FBar.hxx
 * \brief  This file declares the computations of the FBar formulation by
 * kernels of MFEM, for all the integration points of a material at once
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_FBAR_HXX
#define LIB_MFEM_MGIS_BATCH_FBAR_HXX

#include "mfem/linalg/vector.hpp"
#include "MGIS/Context.hxx"
#include "MGIS/FieldView.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"

namespace mfem {
  // forward declaration
  class IntegrationRule;
}  // end of namespace mfem

namespace mfem_mgis {
  // forward declarations
  struct Material;
  struct PartialQuadratureSpace;
}  // end of namespace mfem_mgis

namespace mfem_mgis::batch {

  /*!
   * \brief data of the FBar formulation of a material, computed by kernels.
   *
   * The gradient of the behaviour is `F̄ = α F`, `F` being the deformation
   * gradient at the integration point, `α = (J0 / J)^(1/ie)`, `J = det F`,
   * `J0 = det F0`, `F0` being the deformation gradient at the center of the
   * element, `ie` being 3 in 3D and 2 in plane strain (where `F̄33 = 1`).
   * The inner forces are the ones of the standard finite strain
   * formulation. The stiffness matrix of an element is
   *
   * `Σ_q w_q det J_q B_qᵀ (K_q dF̄/dF B_q + K_q dF̄/dF0 B0)`,
   *
   * `K_q` being the tangent operator blocks of the behaviour, `dF̄/dF = α I +
   * F ⊗ dα/dF`, `dF̄/dF0 = F ⊗ dα/dF0`, `dα/dF = -α / (ie J) dJ/dF` and
   * `dα/dF0 = α / (ie J0) dJ0/dF0`, `B0` being the operator `B` at the
   * center of the element. It is computed as the stiffness matrix of the
   * standard formulation with the tangent operator `K'_q = α K_q + (K_q F)
   * ⊗ dα/dF`, plus the matrix of rank one `a ⊗ b`, `a = Σ_q w_q det J_q α
   * / (ie J0) B_qᵀ K_q F`, `b = B0ᵀ dJ0/dF0` (derivative of `J0` with
   * respect to the unknowns of the element). In plane strain, the third
   * component of `F` and of the derivatives of the determinants are
   * ignored, as by the FBar behaviour integrators (`F̄33 = 1`).
   */
  struct MFEM_MGIS_EXPORT FBarData {
    /*!
     * \brief deformation gradients `F` at the integration points of the
     * material, in the global frame, point by point
     */
    mfem::Vector deformation_gradients;
    //! \brief coefficients `α` at the integration points of the material
    mfem::Vector coefficients;
    /*!
     * \brief determinant `J0` and its derivative `dJ0/dF0` at the center of
     * each element of the material (`1 + gradients size` values per
     * element)
     */
    mfem::Vector centers;
    /*!
     * \brief vectors `b` of the elements, in the layout of the element
     * restriction of MFEM, for all the elements of the mesh
     */
    mfem::Vector b;
    /*!
     * \brief vectors `a` of the elements, in the layout of the element
     * restriction of MFEM, for all the elements of the mesh
     */
    mfem::Vector a;
    /*!
     * \brief tangent operator `K'` at the integration points of the
     * material, point by point (see `computeFBarTangentOperator`)
     */
    mfem::Vector tangent_operator;
    //! \brief workspace: `α / (ie J0) K F` at the integration points
    mfem::Vector pseudo_stresses;
    /*!
     * \brief if true, the tangent operator `K'` and the vectors `a` match
     * the current tangent operator blocks of the behaviour
     */
    bool has_tangent_operator = false;
  };

  /*!
   * \return the integration rule giving the center of the elements of a
   * material (rule of order 0 of their geometry), or a null pointer if the
   * material has no element
   * \param[in] qs: quadrature space of the material
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::IntegrationRule*
  getElementsCenterIntegrationRule(const PartialQuadratureSpace&);
  /*!
   * \return true if the FBar formulation is supported by the kernels for
   * the given kind of gradients (deformation gradient in 3D or in 2D)
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isFBarSupportedByKernels(
      const GradientsKind) noexcept;
  /*!
   * \brief computes the data of the FBar formulation at the center of the
   * elements of a material and at its integration points, and replaces its
   * gradients at the end of the time step, computed by `computeGradients`,
   * by the gradients `F̄` of the FBar formulation
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] data: data of the FBar formulation
   * \param[in, out] m: material
   * \param[in] d0: derivatives of the unknowns at the center of the
   * elements
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool applyFBar(Context&,
                                                FBarData&,
                                                Material&,
                                                const UnknownsDerivatives&,
                                                const ElementsGeometry&,
                                                const GradientsKind) noexcept;
  /*!
   * \brief initializes the data of the FBar formulation of a material for
   * null unknowns (deformation gradients equal to the macroscopic
   * gradients, `α = 1`), before the first computation of the gradients
   * (prediction of the first time step)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] data: data of the FBar formulation
   * \param[in] m: material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool initializeFBar(
      Context&,
      FBarData&,
      const Material&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;
  /*!
   * \brief computes the tangent operator `K'` at the integration points of
   * a material and the vectors `a` of its elements, from the tangent
   * operator blocks of its behaviour (see `FBarData`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks in the global frame
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeFBarTangentOperator(
      Context&,
      FBarData&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;
  /*!
   * \brief adds the stiffness matrices of the elements of a material in the
   * FBar formulation to matrices stored element by element (see
   * `addStiffnessMatrices`), from the data computed by
   * `computeFBarTangentOperator`
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: matrices element by element
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarStiffnessMatrices(
      Context&,
      mfem::Vector&,
      const FBarData&,
      const PartialQuadratureSpace&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;
  /*!
   * \brief adds the stiffness matrices of the elements of a material in the
   * FBar formulation to an assembled sparse matrix (see
   * `addStiffnessMatrices`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: assembled matrix
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarStiffnessMatrices(
      Context&,
      const AssembledStiffnessMatrix&,
      const FBarData&,
      const PartialQuadratureSpace&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;
  /*!
   * \brief adds the diagonals of the stiffness matrices of the elements of
   * a material in the FBar formulation (see `addStiffnessMatricesDiagonals`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] d: values element by element
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarStiffnessMatricesDiagonals(
      Context&,
      mfem::Vector&,
      const FBarData&,
      const PartialQuadratureSpace&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;
  /*!
   * \brief adds the action of the stiffness matrices of the elements of a
   * material in the FBar formulation, without assembling them (see
   * `addTangentOperatorAction`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] y: action, element by element
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] x: values element by element
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarTangentOperatorAction(
      Context&,
      mfem::Vector&,
      const FBarData&,
      const PartialQuadratureSpace&,
      const mfem::Vector&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_FBAR_HXX */
