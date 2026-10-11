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
#include <cstddef>
#include <cstdint>
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
   * \brief scaling and update of rank one of the tangent operator blocks at
   * each integration point of a material, `K' = α K + v ⊗ c` (tangent
   * operator of the FBar formulation, see `FBarData`), applied by the
   * kernels when they multiply the tangent operator blocks by a vector:
   * `K' x = α K x + v (c · x)`, the modified blocks not being stored. Only
   * supported for the deformation gradients. The arrays are given for all
   * the points of the material, in the order of its quadrature space (one
   * coefficient, and as many values of `v` and of `c` as components of the
   * gradients, point by point), located on the device if a device is
   * enabled; the modification is not applied if `coefficients` is null.
   */
  struct MFEM_MGIS_EXPORT TangentOperatorModification {
    //! \brief coefficients `α` at the points of the material
    const real* coefficients = nullptr;
    //! \brief vectors `v` at the points of the material, point by point
    const real* v = nullptr;
    //! \brief vectors `c` at the points of the material, point by point
    const real* c = nullptr;
  };

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
   * \param[in] mod: modification of the tangent operator blocks, if any
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addStiffnessMatrices(
      Context&,
      mfem::Vector&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind,
      const TangentOperatorModification& = {}) noexcept;

  /*!
   * \brief adds the diagonals of the stiffness matrices of the elements of a
   * material to values stored element by element, in the layout of the
   * element restriction of MFEM (`node + nd * (component + dim * element)`,
   * as `addInnerForces`), on the device if a device is enabled
   * \return true on success
   * \param[in, out] d: values element by element
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks in the global frame
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients
   * \param[in] mod: modification of the tangent operator blocks, if any
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addStiffnessMatricesDiagonals(
      Context&,
      mfem::Vector&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind,
      const TangentOperatorModification& = {}) noexcept;
  /*!
   * \brief adds the action of the stiffness matrices of the elements of a
   * material, computed from the tangent operator blocks without assembling
   * them, `Σ Bᵀ K B x w det J`, on values given element by element, by a
   * kernel which computes the increments of the gradients and of the
   * thermodynamic forces at the points of each element without storing
   * them, on the device if a device is enabled. The values are given in the
   * layout of the element restriction of MFEM (`node + nd * (component +
   * dim * element)`, as `addInnerForces`), for all the elements of the mesh.
   * \return true on success
   * \param[in, out] y: action, element by element
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks in the global frame
   * \param[in] x: values element by element
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients
   * \param[in] mod: modification of the tangent operator blocks, if any
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addTangentOperatorAction(
      Context&,
      mfem::Vector&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const mfem::Vector&,
      const ElementsGeometry&,
      const GradientsKind,
      const TangentOperatorModification& = {}) noexcept;

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
    /*!
     * \brief positions of the entries of the stiffness matrices of all the
     * elements in the rows of `matrix` (see
     * `computeStiffnessMatricesPositions`), or a null pointer or an empty
     * array: the entries are then searched in their rows. Otherwise, the
     * array must have `n * n` values per element of the mesh, `n` being the
     * number of degrees of freedom of an element. The position of an entry
     * is read first, and used only if it designates the column of the entry
     * in its row, the entry being searched otherwise (and counted in
     * `searched_entries`): positions computed for another matrix, with the
     * same number of values, give the same results, more slowly.
     */
    const mfem::Array<std::uint16_t>* positions = nullptr;
    /*!
     * \brief if not null, the number of entries whose position, read in
     * `positions`, does not designate their column, and which are searched
     * in their rows, is added to this counter
     */
    std::size_t* searched_entries = nullptr;
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
   * \brief computes, on the device if a device is enabled, the positions of
   * the entries of the stiffness matrices of all the elements in the rows of
   * an assembled sparse matrix, relative to the first entry of each row:
   * `n * n` values per element, `n` being its number of degrees of freedom,
   * the value of the entry `(i, j)` of the element `e` being at the index
   * `(n * e + i) * n + j`. The assembly of the stiffness matrices then reads
   * the position of each entry instead of searching it in its row (see
   * `AssembledStiffnessMatrix::positions`). An entry which is not in the
   * matrix, or whose relative position does not fit in 16 bits, is given
   * the largest value of `std::uint16_t` (it is then searched). The
   * positions take `2 n²` bytes per element, of the order of the size of
   * the matrix for linear elements. If the mesh has no element, or if the
   * number of positions does not fit in an `int`, the positions are not
   * computed and `positions` is left empty.
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] positions: positions
   * \param[in] fes: finite element space
   * \param[in] dofs: degrees of freedom of all the elements (see
   * `computeElementsDegreesOfFreedom`)
   * \param[in] m: assembled matrix, finalized
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool computeStiffnessMatricesPositions(
      Context&,
      mfem::Array<std::uint16_t>&,
      const mfem::FiniteElementSpace&,
      const mfem::Array<int>&,
      const mfem::SparseMatrix&) noexcept;
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
   * storing them: the position of each entry is read in
   * `AssembledStiffnessMatrix::positions`, or searched in its row (by
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
   * \param[in] mod: modification of the tangent operator blocks, if any
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addStiffnessMatrices(
      Context&,
      const AssembledStiffnessMatrix&,
      const PartialQuadratureSpace&,
      const mgis::FieldView<const real>&,
      const ElementsGeometry&,
      const GradientsKind,
      const TangentOperatorModification& = {}) noexcept;

}  // end of namespace mfem_mgis::batch

#endif /* LIB_MFEM_MGIS_BATCH_STIFFNESSMATRICES_HXX */
