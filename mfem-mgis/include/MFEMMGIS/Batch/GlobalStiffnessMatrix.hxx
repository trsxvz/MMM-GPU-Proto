/*!
 * \file   include/MFEMMGIS/Batch/GlobalStiffnessMatrix.hxx
 * \brief  This file declares the assembly of the global matrix of a
 * parallel finite element space from the local matrices of the processes,
 * whose pattern and exchanges between processes are computed once
 * \date   11/10/2026
 */

#ifndef LIB_MFEM_MGIS_BATCH_GLOBALSTIFFNESSMATRIX_HXX
#define LIB_MFEM_MGIS_BATCH_GLOBALSTIFFNESSMATRIX_HXX

#include "mfem/config/config.hpp"
#include "MFEMMGIS/Config.hxx"

#ifdef MFEM_USE_MPI

#include <memory>
#include <vector>
#include "mfem/linalg/hypre.hpp"

namespace mfem {
  // forward declarations
  class SparseMatrix;
  class ParFiniteElementSpace;
}  // end of namespace mfem

namespace mfem_mgis::batch {

  /*!
   * \brief global matrix of a parallel finite element space, `Pᵀ A P`, `A`
   * being the block diagonal matrix of the local matrices of the processes
   * (on the local vectors of the space) and `P` the prolongation of the true
   * degrees of freedom (`mfem::ParFiniteElementSpace::Dof_TrueDof_Matrix`),
   * as built by `mfem::ParNonlinearForm::GetGradient`.
   *
   * When each row of `P` has a single entry, equal to one (conforming
   * spaces), each entry `(i, j)` of a local matrix is added to the entry
   * `(T(i), T(j))` of the global matrix, `T` being the global numbering of
   * the true degrees of freedom. The first assembly computes the pattern of
   * the global matrix (product of hypre, on a matrix of ones, so that no
   * entry is lost), the position of each entry of the local matrix in the
   * part of the global matrix owned by the process, or its place in the
   * message sent to the process owning its row, and the positions of the
   * received entries. These positions give the list of the contributions to
   * each entry of the global matrix: entries of the local matrix, in
   * increasing order, then received values, in the order of the processes.
   * The following assemblies only copy, exchange and add the values, as
   * long as, on every process, the local matrix has the same number of
   * rows, the same number of entries and the same arrays of the pattern
   * (rows and columns), and the finite element space is not updated:
   * otherwise, the setup is made again by all the processes. The contents
   * of the arrays of the pattern are not compared: they must not be
   * modified between two assemblies.
   *
   * Each value of the global matrix is the sum of its contributions, added
   * in the order of their list by a kernel run where hypre computes: on the
   * host, or on the GPU, where the values of the local matrix are read and
   * the ones of the global matrix written without copy to the host. The
   * values are the same, bit for bit, on the host and on the GPU. They are
   * the ones of the product `Pᵀ A P` of hypre up to round-off errors (sums
   * in another order). The values exchanged between the processes are given
   * to MPI on the GPU with GPU-aware MPI (`mfem::Device::GetGPUAwareMPI`),
   * on the host otherwise. The setup is made on the host, from copies of
   * the prolongation and of the pattern of the global matrix: on a GPU, the
   * matrices are not moved to the host, and the lists of the contributions
   * are copied to the GPU without being kept on the host.
   *
   * The global matrix is the same object at each assembly. Each assembly
   * frees the copies of its former values kept in it on a GPU: its
   * transposed blocks, built by MFEM for the transposed products
   * (`mfem::HypreParMatrix::ResetTranspose`), and the copies of hypre
   * (sorted columns and triangular solves of the Gauss-Seidel relaxations).
   */
  struct MFEM_MGIS_EXPORT GlobalStiffnessMatrix {
    /*!
     * \return true if the global matrices of a finite element space can be
     * assembled by this class: each row of the prolongation has a single
     * entry, equal to one, on all the processes (collective call)
     * \param[in] fes: finite element space
     */
    [[nodiscard]] static bool isSupported(
        const mfem::ParFiniteElementSpace& fes) noexcept;
    /*!
     * \brief constructor
     * \param[in] f: finite element space (see `isSupported`), which must
     * outlive this object
     */
    explicit GlobalStiffnessMatrix(const mfem::ParFiniteElementSpace& f);
    //! \brief move constructor
    GlobalStiffnessMatrix(GlobalStiffnessMatrix&&) noexcept;
    //! \brief move assignment
    GlobalStiffnessMatrix& operator=(GlobalStiffnessMatrix&&) noexcept;
    //! \brief destructor
    ~GlobalStiffnessMatrix();
    /*!
     * \brief assembles the global matrix from the local matrix of the
     * process (collective call)
     * \return the global matrix, owned by this object, whose address does
     * not change while the setup is not made again, or a null pointer on
     * error (on all the processes, the error message being registered on
     * the processes where it occurred, and the message `failure on another
     * process` on the other ones). An exception thrown between two
     * collective calls of the setup (for example a failed allocation) is not
     * shared with the other processes.
     * \param[in, out] ctx: execution context
     * \param[in] A: local matrix, finalized, on the local vector of the
     * finite element space
     */
    [[nodiscard]] mfem::HypreParMatrix* assemble(
        Context& ctx, const mfem::SparseMatrix& A) noexcept;

   private:
    /*!
     * \brief computes the pattern of the global matrix and the exchanges
     * between the processes
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] A: local matrix
     */
    [[nodiscard]] bool setup(Context&, const mfem::SparseMatrix&);
    //! \brief finite element space
    const mfem::ParFiniteElementSpace* fes;
    //! \brief global matrix
    std::unique_ptr<mfem::HypreParMatrix> matrix;
    //! \brief rows of the pattern of the local matrix of the setup
    const int* rows = nullptr;
    //! \brief columns of the pattern of the local matrix of the setup
    const int* columns = nullptr;
    //! \brief sequence of the finite element space at the setup
    long sequence = -1;
    //! \brief number of rows of the local matrix
    int number_of_rows = -1;
    //! \brief number of entries of the local matrix
    int number_of_entries = -1;
    //! \brief number of entries of the diagonal block of the global matrix
    int diagonal_block_entries = 0;
    /*!
     * \brief number of entries of the off diagonal block of the global
     * matrix
     */
    int off_diagonal_block_entries = 0;
    /*!
     * \brief offsets of the contributions to each entry of the global
     * matrix: entries of the diagonal block, then of the off diagonal block
     */
    mfem::Array<int> offsets;
    /*!
     * \brief contributions to the entries of the global matrix: for each
     * entry, the entries of the local matrix added to it, in increasing
     * order, then, from the number of entries of the local matrix, the
     * received values added to it, in increasing order
     */
    mfem::Array<int> contributions;
    //! \brief entries of the local matrix sent, ordered by process
    mfem::Array<int> sent_entries;
    //! \brief values exchanged with a process
    struct Message {
      //! \brief rank of the process
      int process;
      //! \brief offset of the values in the values sent or received
      int offset;
      //! \brief number of values
      int size;
    };
    //! \brief messages sent to the processes owning rows of the local matrix
    std::vector<Message> sent_messages;
    //! \brief messages received from the other processes
    std::vector<Message> received_messages;
    //! \brief requests of the exchanges (workspace)
    std::vector<MPI_Request> requests;
    //! \brief values sent (workspace)
    mfem::Vector send_buffer;
    //! \brief values received (workspace)
    mfem::Vector receive_buffer;
  };

}  // end of namespace mfem_mgis::batch

#endif /* MFEM_USE_MPI */

#endif /* LIB_MFEM_MGIS_BATCH_GLOBALSTIFFNESSMATRIX_HXX */
