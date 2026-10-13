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
   * ⊗ dα/dF`, applied by the kernels of the stiffness matrices without being
   * stored (only `v = K_q F` and `c = dα/dF` are stored, see
   * `TangentOperatorModification`), plus the matrix of rank one `a ⊗ b`,
   * `a = Σ_q w_q det J_q α / (ie J0) B_qᵀ K_q F`, `b = B0ᵀ dJ0/dF0`
   * (derivative of `J0` with respect to the unknowns of the element). The
   * functions computing the stiffness matrices must be given the tangent
   * operator blocks from which `v` and `c` are computed (checked by
   * `tangent_operator_blocks`). In plane strain, the third component of `F`
   * and of the derivatives of the determinants are ignored in `v` and `c`,
   * as by the FBar behaviour integrators (`F̄33 = 1`); `J` and `J0` are
   * computed with the macroscopic `F33`, one by default.
   *
   * The data are computed with the gradients by kernels (`applyFBar`): the
   * tangent operator computed by kernels is the one of the unknowns of the
   * last computation of the gradients (or of the beginning of the time step
   * after a revert, see `revertFBarData`), whatever the unknowns given to
   * the assembly of the jacobian, and it can't be computed if the gradients
   * have been computed element by element (see `valid`, set by the
   * behaviour integrator; `computeFBarTangentOperator` does not check it,
   * the caller must).
   */
  struct MFEM_MGIS_EXPORT FBarData {
    /*!
     * \brief deformation gradients `F` at the integration points of the
     * material, in the global frame, point by point
     */
    mfem::Vector deformation_gradients;
    //! \brief coefficients `α` at the integration points of the material
    mfem::Vector coefficients;
    //! \brief determinant `J0` at the center of each element of the material
    mfem::Vector centers;
    /*!
     * \brief vectors `b` of the elements of the material (`node + nd *
     * (component + vdim * k)` for the `k`-th element of the material)
     */
    mfem::Vector b;
    //! \brief vectors `a` of the elements of the material (as `b`)
    mfem::Vector a;
    /*!
     * \brief vectors `v = K F` of the tangent operator `K' = α K + v ⊗ c`
     * at the integration points of the material, point by point (see
     * `computeFBarTangentOperator`)
     */
    mfem::Vector tangent_v;
    //! \brief vectors `c = dα/dF` of the tangent operator `K'`
    mfem::Vector tangent_c;
    /*!
     * \brief tangent operator blocks of the behaviour from which `v` and `c`
     * are computed: the functions computing the stiffness matrices check
     * that they are given the same blocks
     */
    const real* tangent_operator_blocks = nullptr;
    //! \brief workspace: `α / (ie J0) K F` at the integration points
    mfem::Vector pseudo_stresses;
    /*!
     * \brief if true, the tangent operator `K'` and the vectors `a` match
     * the current tangent operator blocks of the behaviour
     */
    bool has_tangent_operator = false;
    /*!
     * \brief false if the gradients have been computed element by element
     * since the last computation of the data by kernels: the data do not
     * match the current unknowns, and the tangent operator can't be computed
     * by kernels. It is set by the behaviour integrator:
     * `computeFBarTangentOperator` does not check it, the caller must.
     */
    bool valid = true;
    //! \brief copies of the data at the beginning of the time step
    struct {
      //! \brief deformation gradients `F`
      mfem::Vector deformation_gradients;
      //! \brief coefficients `α`
      mfem::Vector coefficients;
      //! \brief data at the center of the elements
      mfem::Vector centers;
      //! \brief vectors `b`
      mfem::Vector b;
      //! \brief validity of the data (see `FBarData::valid`)
      bool valid = true;
      //! \brief if true, the copies are set (see `updateFBarData`)
      bool defined = false;
    } beginning_of_time_step;
  };

  /*!
   * \return the integration rule giving the center of the elements of a
   * material (rule of order 0 of their geometry), or a null pointer if the
   * material has no element
   * \param[in] qs: quadrature space of the material
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::IntegrationRule*
  getElementsCenterIntegrationRule(const PartialQuadratureSpace& qs);
  /*!
   * \return true if the FBar formulation is supported by the kernels for
   * the given kind of gradients (deformation gradient in 3D or in 2D)
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isFBarSupportedByKernels(
      const GradientsKind k) noexcept;
  /*!
   * \brief computes the data of the FBar formulation at the center of the
   * elements of a material and at its integration points, and replaces its
   * gradients at the end of the time step, computed by `computeGradients`,
   * by the gradients `F̄` of the FBar formulation
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] data: data of the FBar formulation
   * \param[in, out] m: material
   * \param[in] d0: derivatives of the unknowns at the center of the
   * elements
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool applyFBar(Context& ctx,
                                                FBarData& data,
                                                Material& m,
                                                const UnknownsDerivatives& d0,
                                                const ElementsGeometry& g0,
                                                const GradientsKind k) noexcept;
  /*!
   * \brief initializes the data of the FBar formulation of a material for
   * null unknowns (deformation gradients equal to the macroscopic
   * gradients, `α = 1`), before the first computation of the gradients
   * (prediction of the first time step)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] data: data of the FBar formulation
   * \param[in] m: material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool initializeFBar(
      Context& ctx,
      FBarData& data,
      const Material& m,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;
  /*!
   * \brief keeps the data of the FBar formulation, computed for the
   * converged unknowns of the time step, as the data at the beginning of the
   * next time step (see `revertFBarData`). The tangent operator `K'` must
   * then be computed again (`has_tangent_operator` is reset).
   * \param[in, out] data: data of the FBar formulation
   */
  MFEM_MGIS_EXPORT void updateFBarData(FBarData& data);
  /*!
   * \brief restores the data of the FBar formulation at the beginning of the
   * time step, kept by `updateFBarData`, or, before the first update, removes
   * the data (the data at the beginning of the first time step are computed
   * for null unknowns, see `initializeFBar`, which supposes that the initial
   * unknowns are null). The tangent operator `K'` must then be computed
   * again (`has_tangent_operator` is reset).
   * \param[in, out] data: data of the FBar formulation
   */
  MFEM_MGIS_EXPORT void revertFBarData(FBarData& data);
  /*!
   * \brief computes the vectors `v` and `c` of the tangent operator `K'` at
   * the integration points of a material and the vectors `a` of its
   * elements, from the tangent operator blocks of its behaviour (see
   * `FBarData`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks in the global frame
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeFBarTangentOperator(
      Context& ctx,
      FBarData& data,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;
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
   * \param[in] K: tangent operator blocks of the behaviour in the global
   * frame (see `computeFBarTangentOperator`)
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarStiffnessMatrices(
      Context& ctx,
      mfem::Vector& m,
      const FBarData& data,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the stiffness matrices of the elements of a material in the
   * FBar formulation to an assembled sparse matrix (see
   * `addStiffnessMatrices`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: assembled matrix
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks of the behaviour in the global
   * frame (see `computeFBarTangentOperator`)
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarStiffnessMatrices(
      Context& ctx,
      const AssembledStiffnessMatrix& m,
      const FBarData& data,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the diagonals of the stiffness matrices of the elements of
   * a material in the FBar formulation (see `addStiffnessMatricesDiagonals`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] d: values element by element
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks of the behaviour in the global
   * frame (see `computeFBarTangentOperator`)
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarStiffnessMatricesDiagonals(
      Context& ctx,
      mfem::Vector& d,
      const FBarData& data,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the action of the stiffness matrices of the elements of a
   * material in the FBar formulation, without assembling them (see
   * `addTangentOperatorAction`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] y: action, element by element
   * \param[in] data: data of the FBar formulation
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks of the behaviour in the global
   * frame (see `computeFBarTangentOperator`)
   * \param[in] x: values element by element
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFBarTangentOperatorAction(
      Context& ctx,
      mfem::Vector& y,
      const FBarData& data,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const mfem::Vector& x,
      const ElementsGeometry& g,
      const GradientsKind k) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_FBAR_HXX */
