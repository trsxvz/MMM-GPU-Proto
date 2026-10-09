/*!
 * \file   include/MFEMMGIS/Batch/StiffnessMatrices.hxx
 * \brief  This file declares the computation of the stiffness matrices of
 * the elements of the standard mechanical behaviour integrators by kernels
 * of MFEM, for all the elements of a material at once
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_STIFFNESSMATRICES_HXX
#define LIB_MFEM_MGIS_BATCH_STIFFNESSMATRICES_HXX

#include <memory>
#include "mfem/general/array.hpp"
#include "mfem/linalg/vector.hpp"
#include "MGIS/FieldView.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"

namespace mfem {
  // forward declarations
  class SparseMatrix;
  class FiniteElementSpace;
}  // end of namespace mfem

namespace mfem_mgis {
  // forward declaration
  struct PartialQuadratureSpace;
}  // end of namespace mfem_mgis

namespace mfem_mgis::batch {

  /*!
   * \brief adds the stiffness matrices of the elements of a material, `Σ Bᵀ
   * K B w det J` summed over the integration points, to matrices stored
   * element by element, on the device if a device is enabled. The matrix of
   * an element has `n = nd * vdim` rows and columns, the degrees of freedom
   * being ordered as in the element restriction of MFEM (native ordering of
   * the nodes, `node + nd * component`). It is stored column by column, as
   * a `mfem::DenseMatrix`, at the offset `n * n * e` for the element `e` of
   * the mesh.
   * \return true on success
   * \param[in, out] m: matrices element by element
   * \param[in] qs: quadrature space of the material. All its elements must
   * use the integration rule of the geometric data.
   * \param[in] K: tangent operator blocks in the global frame at all the
   * points of the material, row by row, located on the device if a device
   * is enabled
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients, which defines the operator `B`
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addStiffnessMatrices(
      Context&,
      mfem::Vector&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;

  /*!
   * \brief sparse matrix in which the stiffness matrices of the elements are
   * assembled by kernels (see `addStiffnessMatrices`)
   */
  struct MFEM_MGIS_EXPORT AssembledStiffnessMatrix {
    /*!
     * \brief matrix, finalized, on the local vector of the finite element
     * space, whose pattern contains all the entries of the stiffness
     * matrices of the elements (see `makeAssembledStiffnessMatrix`)
     */
    mfem::SparseMatrix* matrix = nullptr;
    /*!
     * \brief degrees of freedom of all the elements of the mesh (see
     * `computeElementsDegreesOfFreedom`)
     */
    const mfem::Array<int>* element_dofs = nullptr;
  };

  /*!
   * \brief computes the degrees of freedom of all the elements of the mesh,
   * in the ordering of `mfem::FiniteElementSpace::GetElementVDofs` (`nd *
   * vdim` values per element), which is the ordering of the stiffness
   * matrices of the elements (see `addStiffnessMatrices`)
   * \return true on success. The degrees of freedom with an orientation
   * (negative indices) are refused.
   * \param[in, out] ctx: execution context
   * \param[out] dofs: degrees of freedom, readable on the device
   * \param[in] fes: finite element space
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeElementsDegreesOfFreedom(
      Context&, mfem::Array<int>&, const mfem::FiniteElementSpace&) noexcept;
  /*!
   * \brief creates the sparse matrix in which the stiffness matrices of all
   * the elements are assembled (see `AssembledStiffnessMatrix`), on the
   * local vector of the finite element space: its pattern contains the
   * couplings between the degrees of freedom of each element, and its
   * values are null. In each row, the diagonal entry comes first, then the
   * other columns, sorted: the ordering of hypre, which
   * `mfem::HypreParMatrix` imposes on the rows of a local matrix that it
   * uses (the columns of a sorted row would not stay sorted).
   * \return the matrix, or a null pointer on error
   * \param[in, out] ctx: execution context
   * \param[in] fes: finite element space
   * \param[in] dofs: degrees of freedom of all the elements (see
   * `computeElementsDegreesOfFreedom`)
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::unique_ptr<mfem::SparseMatrix>
  makeAssembledStiffnessMatrix(Context&,
                               const mfem::FiniteElementSpace&,
                               const mfem::Array<int>&) noexcept;
  /*!
   * \brief sets the values of a sparse matrix to zero, on the device if a
   * device is enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] m: matrix (finalized)
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool zeroSparseMatrixValues(
      Context&, mfem::SparseMatrix&) noexcept;
  /*!
   * \brief adds the stiffness matrices of the elements of a material to an
   * assembled sparse matrix, on the device if a device is enabled, without
   * storing them: the position of each entry is searched in its row (by
   * bisection for the ordering of `makeAssembledStiffnessMatrix` or sorted
   * rows, linearly otherwise), and the additions are atomic (several
   * materials may share degrees of freedom)
   * \return true on success, false if an entry is not in the pattern of the
   * matrix
   * \param[in, out] m: assembled matrix
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks in the global frame
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients, which defines the operator `B`
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addStiffnessMatrices(
      Context&,
      const AssembledStiffnessMatrix&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_STIFFNESSMATRICES_HXX */
