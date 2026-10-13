/*!
 * \file   src/Batch/StiffnessMatricesKernels.hxx
 * \brief  This file declares the kernels computing the stiffness matrices
 * of the elements of the standard behaviour integrators (mechanics, heat
 * transfer, micromorphic damage) (internal header)
 * \date   09/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_BATCH_STIFFNESSMATRICESKERNELS_HXX
#define LIB_MFEM_MGIS_SRC_BATCH_STIFFNESSMATRICESKERNELS_HXX

#include <limits>
#include <cstddef>
#include <cstdint>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "KernelsUtilities.hxx"

namespace mfem_mgis::batch {

  // forward declaration
  struct AssembledStiffnessMatrix;

}  // end of namespace mfem_mgis::batch

namespace mfem_mgis::batch::internals {

  /*!
   * \brief position given to the entries which are not in the matrix, or
   * whose position relative to the first entry of their row does not fit in
   * 16 bits (see `computeStiffnessMatricesPositions`)
   */
  inline constexpr int unknown_position =
      std::numeric_limits<std::uint16_t>::max();

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
    /*!
     * \brief values of the shape functions on the reference element
     * (`nq x nd`), used by the derivatives of the thermodynamic forces with
     * respect to the unknowns (stationary heat transfer behaviours) and by
     * the gradients which are the values of the unknowns (transient term of
     * the heat equation)
     */
    const real* shape_functions = nullptr;
    /*!
     * \brief inverses of the jacobians of the element transformations (`nq x
     * dim x dim x ne`, see `ElementsJacobians::invJ`)
     */
    const real* inverse_jacobians = nullptr;
    //! \brief determinants of the jacobians (`nq x ne`)
    const real* determinants = nullptr;
    //! \brief weights of the integration points on the reference element
    const real* weights = nullptr;
    /*!
     * \brief tangent operator blocks: derivatives of the thermodynamic forces
     * with respect to the gradients, followed, for the heat transfer
     * behaviours, by their derivatives with respect to the temperature
     */
    const real* tangent_operator = nullptr;
    //! \brief stride between two points of the tangent operator blocks
    long long point_stride = 0;
    //! \brief stride between two components of the tangent operator blocks
    long long component_stride = 0;
    /*!
     * \brief coefficients `α` of the modification of the tangent operator
     * blocks at the points, `K' = α K + v ⊗ c` (see
     * `TangentOperatorModification`), null if the blocks are not modified
     */
    const real* modification_coefficients = nullptr;
    //! \brief vectors `v` of the modification, point by point
    const real* modification_v = nullptr;
    //! \brief vectors `c` of the modification, point by point
    const real* modification_c = nullptr;
    //! \brief matrices element by element
    real* matrices = nullptr;
    //! \brief number of elements of the material
    int number_of_elements = 0;
    //! \brief number of degrees of freedom per element and per component
    int nd = 0;
    //! \brief number of integration points per element
    int nq = 0;
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
     * \brief positions of the entries of the stiffness matrices of the
     * elements relative to the first entry of their rows (see
     * `computeStiffnessMatricesPositions`), or a null pointer
     */
    const std::uint16_t* positions = nullptr;
    /*!
     * \brief number of entries which are not in the pattern of the matrix,
     * or whose row is not in the matrix
     */
    int* missing = nullptr;
    /*!
     * \brief number of entries whose position, given by `positions`, does
     * not designate their column (they are searched in their rows), or a
     * null pointer. The entries whose position is unknown are not counted
     * (see `unknown_position`).
     */
    int* searched = nullptr;
    /*!
     * \return the position of an entry of the stiffness matrix of an element
     * in the values of the matrix, or `end` if the entry is not in the
     * matrix: the position read in `positions` if it designates the column
     * of the entry, the position searched in the row otherwise (see
     * `findColumn`, the entry being then counted in `searched`)
     * \param[in] e: element (in the mesh)
     * \param[in] n: number of degrees of freedom of an element
     * \param[in] i: row of the entry in the stiffness matrix of the element
     * \param[in] j: column of the entry in the stiffness matrix of the
     * element
     * \param[in] begin: position of the first entry of the row of the entry
     * \param[in] end: position after the last entry of the row of the entry
     * \param[in] c: column of the entry in the matrix
     */
    MFEM_HOST_DEVICE int findPosition(const std::size_t e,
                                      const std::size_t n,
                                      const std::size_t i,
                                      const std::size_t j,
                                      const int begin,
                                      const int end,
                                      const int c) const {
      if (this->positions != nullptr) {
        const auto o = static_cast<int>(this->positions[(n * e + i) * n + j]);
        if ((o < end - begin) && (this->J[begin + o] == c)) {
          return begin + o;
        }
        if ((this->searched != nullptr) && (o != unknown_position)) {
          ::AtomicAdd(this->searched[0], 1);
        }
      }
      return findColumn(this->J, begin, end, c);
    }
  };

  /*!
   * \return true if an assembled matrix can receive the stiffness matrices
   * of the elements, false otherwise (an error message is then registered):
   * finalized matrix, degrees of freedom of all the elements, positions of
   * the entries of their stiffness matrices (if any) of the expected sizes
   * \param[in, out] ctx: execution context
   * \param[in] m: assembled matrix
   * \param[in] ne: number of elements of the mesh
   * \param[in] n: number of degrees of freedom of an element
   * \param[in] f: calling function (error messages)
   */
  [[nodiscard]] bool checkAssembledStiffnessMatrix(
      Context&,
      const AssembledStiffnessMatrix&,
      const std::int64_t,
      const std::int64_t,
      const char* const);
  /*!
   * \brief fills the destination of the kernels adding the stiffness
   * matrices of the elements to an assembled matrix, checked by
   * `checkAssembledStiffnessMatrix`
   * \param[out] o: destination of the kernels
   * \param[out] counters: counters of the entries which are not in the
   * pattern of the matrix (`missing`) and of the entries searched in their
   * rows (`searched`), set to zero
   * \param[in] m: assembled matrix
   * \param[in] use_device: if true, the destination is read on the device
   */
  void getSparseMatrixOutput(SparseMatrixOutput&,
                             mfem::Array<int>&,
                             const AssembledStiffnessMatrix&,
                             const bool);
  /*!
   * \brief reads the counters filled by a kernel adding the stiffness
   * matrices of the elements to an assembled matrix: the number of entries
   * searched in their rows is added to
   * `AssembledStiffnessMatrix::searched_entries`, if given
   * \return true if all the entries are in the pattern of the matrix, false
   * otherwise (an error message is then registered)
   * \param[in, out] ctx: execution context
   * \param[in] counters: counters (see `getSparseMatrixOutput`)
   * \param[in] m: assembled matrix
   * \param[in] f: calling function (error message)
   */
  [[nodiscard]] bool checkSparseMatrixOutputCounters(
      Context&,
      const mfem::Array<int>&,
      const AssembledStiffnessMatrix&,
      const char* const);
  /*!
   * \brief launches the kernel computing the positions of the entries of the
   * stiffness matrices of all the elements in the rows of an assembled
   * matrix (see `computeStiffnessMatricesPositions`): one thread per row of
   * the stiffness matrix of an element
   * \param[in] o: assembled matrix (`values`, `positions`, `missing` and
   * `searched` are not used)
   * \param[in] ne: number of elements of the mesh
   * \param[in] n: number of degrees of freedom of an element
   * \param[out] positions: positions (`ne * n * n` values)
   * \param[in] use_device: if true, the kernel runs on the device
   */
  void launchStiffnessMatricesPositionsKernel(const SparseMatrixOutput&,
                                              const int,
                                              const int,
                                              std::uint16_t* const,
                                              const bool);

  /*!
   * \brief launches the kernel adding stiffness matrices stored element by
   * element to an assembled sparse matrix (see
   * `addStiffnessMatricesByEntries`): one thread per value of the matrix,
   * which adds the contributions to this value in a fixed order, without
   * atomic additions, after checking their columns (the contributions whose
   * column is not the one of the value are counted in `missing` and not
   * added)
   * \param[in] o: assembled matrix (`positions` and `searched` are not used)
   * \param[in] nnz: number of values of the matrix
   * \param[in] offsets: offsets of the values (see `EntriesContributions`)
   * \param[in] contributions: indices of the entries of the stiffness
   * matrices of the elements contributing to each value (see
   * `EntriesContributions`)
   * \param[in] k: stiffness matrices of all the elements, stored element by
   * element, column by column
   * \param[in] n: number of degrees of freedom of an element
   * \param[in] use_device: if true, the kernel runs on the device
   */
  void launchStiffnessMatricesByEntriesKernel(const SparseMatrixOutput&,
                                              const int,
                                              const int* const,
                                              const int* const,
                                              const real* const,
                                              const int,
                                              const bool);

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
   * \brief launches the kernel adding the diagonals of the stiffness
   * matrices of the elements of a material to values stored element by
   * element, in the layout of the element restriction of MFEM (`node + nd *
   * (component + vdim * element)`, `vdim` being the number of components
   * of the unknowns)
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments (`matrices` is not used)
   * \param[in, out] d: values element by element
   * \param[in] use_device: if true, the kernel runs on the device
   */
  void launchStiffnessMatricesDiagonalsKernel(
      const GradientsKind,
      const StiffnessMatricesKernelArguments&,
      real* const,
      const bool);
  /*!
   * \brief launches the kernel adding the action of the stiffness matrices of
   * the elements of a material, `Σ Bᵀ K B x w det J`, on values given
   * element by element, in the layout of the element restriction of MFEM
   * (`node + nd * (component + vdim * element)`, `vdim` being the number
   * of components of the unknowns)
   * \param[in] k: kind of the gradients
   * \param[in] a: arguments (`matrices` is not used)
   * \param[in] x: values element by element
   * \param[in, out] y: action, element by element
   * \param[in] use_device: if true, the kernel runs on the device
   */
  void launchTangentOperatorActionKernel(
      const GradientsKind,
      const StiffnessMatricesKernelArguments&,
      const real* const,
      real* const,
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
