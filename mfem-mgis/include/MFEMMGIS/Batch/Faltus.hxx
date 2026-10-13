/*!
 * \file   include/MFEMMGIS/Batch/Faltus.hxx
 * \brief  This file declares the computations of the regularization of
 * Faltus et al. by kernels of MFEM, for all the elements of a material at
 * once
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_FALTUS_HXX
#define LIB_MFEM_MGIS_BATCH_FALTUS_HXX

#include "mfem/linalg/vector.hpp"
#include "MGIS/Context.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"

namespace mfem_mgis {
  // forward declarations
  struct Material;
  struct PartialQuadratureSpace;
}  // end of namespace mfem_mgis

namespace mfem_mgis::batch {

  /*!
   * \brief data of the regularization of Faltus et al. of a material,
   * computed by kernels.
   *
   * The regularization adds to the inner forces of an element
   *
   * `Σ_q w_q det J_q D_qᵀ α (F_q - F0)`,
   *
   * and to its stiffness matrix
   *
   * `α Σ_q w_q det J_q D_qᵀ D_q`,
   *
   * `F_q` being the deformation gradient at the integration point (gradient
   * of the behaviour), `F0` the deformation gradient at the center of the
   * element, `D_q = B_q - B0`, `B_q` and `B0` being the operators of the
   * gradients at the integration point and at the center of the element,
   * and `α` the penalization coefficient. The stiffness matrix does not
   * depend on the unknowns: it is computed on the fly by the kernels, from
   * the geometric data of the integration points and of the center of the
   * elements.
   *
   * `F0` is computed with the gradients by kernels (see
   * `computeFaltusCenters`): the inner forces computed by kernels use the
   * gradients and `F0` of the last computation of the gradients (or of the
   * beginning of the time step after a revert, see `revertFaltusData`), and
   * can't be computed if the gradients have been computed element by element
   * (see `valid`). Before the first computation of the gradients
   * (prediction of the first time step), `F0` is computed for null unknowns,
   * as the data of the FBar formulation: the initial unknowns are supposed
   * to be null.
   */
  struct MFEM_MGIS_EXPORT FaltusData {
    /*!
     * \brief deformation gradients `F0` at the center of the elements of the
     * material (`c + gs * k` for the component `c` of the `k`-th element of
     * the material, `gs` being the number of components of the gradients)
     */
    mfem::Vector centers;
    /*!
     * \brief false if the gradients have been computed element by element
     * since the last computation of `F0` by kernels, or if this computation
     * failed after the checks of the kind of the gradients and of the size
     * of `centers` (`valid` is left unchanged by these two failures): `F0`
     * does not match the gradients, and the inner forces can't be computed
     * by kernels
     */
    bool valid = true;
    //! \brief copies of the data at the beginning of the time step
    struct {
      //! \brief deformation gradients `F0`
      mfem::Vector centers;
      //! \brief validity of the data (see `FaltusData::valid`)
      bool valid = true;
      //! \brief if true, the copies are set (see `updateFaltusData`)
      bool defined = false;
    } beginning_of_time_step;
  };

  /*!
   * \return true if the regularization of Faltus et al. is supported by the
   * kernels for the given kind of gradients (deformation gradient in 3D or
   * in 2D)
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isFaltusSupportedByKernels(
      const GradientsKind k) noexcept;
  /*!
   * \brief computes the deformation gradients `F0` at the center of the
   * elements of a material
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] data: data of the regularization
   * \param[in] m: material (macroscopic gradients)
   * \param[in] d0: derivatives of the unknowns at the center of the
   * elements, null for null unknowns (before the first computation of the
   * gradients)
   * \param[in] g0: geometric data of the center of the elements (see
   * `getElementsCenterIntegrationRule`)
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeFaltusCenters(
      Context& ctx,
      FaltusData& data,
      const Material& m,
      const UnknownsDerivatives* const d0,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;
  /*!
   * \brief keeps the data of the regularization, computed for the converged
   * unknowns of the time step, as the data at the beginning of the next time
   * step (see `revertFaltusData`)
   * \param[in, out] data: data of the regularization
   */
  MFEM_MGIS_EXPORT void updateFaltusData(FaltusData& data);
  /*!
   * \brief restores the data of the regularization at the beginning of the
   * time step, kept by `updateFaltusData`, or, before the first update,
   * removes the data (computed again for null unknowns)
   * \param[in, out] data: data of the regularization
   */
  MFEM_MGIS_EXPORT void revertFaltusData(FaltusData& data);
  /*!
   * \brief adds the contribution of the regularization to the inner forces
   * of the elements of a material, in the layout of the element restriction
   * of MFEM (`nd x vdim x ne`), on the device if a device is enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] f: values element by element
   * \param[in] data: data of the regularization
   * \param[in] m: material (gradients at the end of the time step)
   * \param[in] alpha: penalization coefficient
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFaltusInnerForces(
      Context& ctx,
      mfem::Vector& f,
      const FaltusData& data,
      const Material& m,
      const real alpha,
      const ElementsGeometry& g,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the contribution of the regularization to the stiffness
   * matrices of the elements of a material, stored element by element (see
   * `addStiffnessMatrices`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: matrices element by element
   * \param[in] qs: quadrature space of the material
   * \param[in] alpha: penalization coefficient
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFaltusStiffnessMatrices(
      Context& ctx,
      mfem::Vector& m,
      const PartialQuadratureSpace& qs,
      const real alpha,
      const ElementsGeometry& g,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the contribution of the regularization to an assembled
   * sparse matrix (see `addStiffnessMatrices`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: assembled matrix
   * \param[in] qs: quadrature space of the material
   * \param[in] alpha: penalization coefficient
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFaltusStiffnessMatrices(
      Context& ctx,
      const AssembledStiffnessMatrix& m,
      const PartialQuadratureSpace& qs,
      const real alpha,
      const ElementsGeometry& g,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the diagonals of the contribution of the regularization to
   * the stiffness matrices of the elements of a material (see
   * `addStiffnessMatricesDiagonals`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] d: values element by element
   * \param[in] qs: quadrature space of the material
   * \param[in] alpha: penalization coefficient
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFaltusStiffnessMatricesDiagonals(
      Context& ctx,
      mfem::Vector& d,
      const PartialQuadratureSpace& qs,
      const real alpha,
      const ElementsGeometry& g,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;
  /*!
   * \brief adds the action of the contribution of the regularization to the
   * stiffness matrices of the elements of a material, without assembling
   * them (see `addTangentOperatorAction`)
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] y: action, element by element
   * \param[in] qs: quadrature space of the material
   * \param[in] x: values element by element
   * \param[in] alpha: penalization coefficient
   * \param[in] g: geometric data of the integration rule of the material
   * \param[in] g0: geometric data of the center of the elements
   * \param[in] k: kind of the gradients
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addFaltusTangentOperatorAction(
      Context& ctx,
      mfem::Vector& y,
      const PartialQuadratureSpace& qs,
      const mfem::Vector& x,
      const real alpha,
      const ElementsGeometry& g,
      const ElementsGeometry& g0,
      const GradientsKind k) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_FALTUS_HXX */
