/*!
 * \file   src/Batch/GlobalStiffnessMatrix.cxx
 * \brief  This file implements the assembly of the global matrix of a
 * parallel finite element space from the local matrices of the processes
 * \date   11/10/2026
 */

#include "MFEMMGIS/Batch/GlobalStiffnessMatrix.hxx"

#ifdef MFEM_USE_MPI

#include <tuple>
#include <limits>
#include <string>
#include <utility>
#include <algorithm>
#include <exception>
#include "mfem/general/forall.hpp"
#include "mfem/general/communication.hpp"
#include "mfem/linalg/handle.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/pfespace.hpp"
#include "GlobalStiffnessMatrixKernels.hxx"

namespace mfem_mgis::batch {

  //! \brief tag of the messages of the assemblies
  constexpr int message_tag = 23571;

  //! \return true if a value is exactly one
  static bool isOne(const real v) { return (v >= real{1}) && (v <= real{1}); }

  //! \brief copy on the host of a block of a matrix of hypre
  struct BlockCopy {
    //! \brief offsets of the rows
    std::vector<HYPRE_Int> rows;
    //! \brief columns of the entries
    std::vector<HYPRE_Int> columns;
    //! \brief values of the entries, if they were asked
    std::vector<HYPRE_Complex> values;
  };

  /*!
   * \return a copy on the host of the pattern of a block of a matrix of
   * hypre, and of its values if they are asked, made from the memory where
   * the block is (`hypre_CSRMatrixMemoryLocation`, which MFEM keeps up to
   * date when it moves the arrays of a matrix): the matrix is neither moved
   * nor copied as a whole, which `mfem::HypreParMatrix::HostRead` does on a
   * GPU, keeping the copy as long as the matrix lives
   * \param[in] block: block
   * \param[in] values: if true, the values are copied
   */
  static BlockCopy copyBlock(const hypre_CSRMatrix* const block,
                             const bool values) {
    const auto location = hypre_CSRMatrixMemoryLocation(block);
    const auto nr = hypre_CSRMatrixNumRows(block);
    const auto nnz = hypre_CSRMatrixNumNonzeros(block);
    auto c = BlockCopy{};
    c.rows.assign(static_cast<std::size_t>(nr) + 1, 0);
    if (hypre_CSRMatrixI(block) != nullptr) {
      hypre_TMemcpy(c.rows.data(), hypre_CSRMatrixI(block), HYPRE_Int, nr + 1,
                    HYPRE_MEMORY_HOST, location);
    }
    if (nnz != 0) {
      c.columns.resize(static_cast<std::size_t>(nnz));
      hypre_TMemcpy(c.columns.data(), hypre_CSRMatrixJ(block), HYPRE_Int, nnz,
                    HYPRE_MEMORY_HOST, location);
      if (values) {
        c.values.resize(static_cast<std::size_t>(nnz));
        hypre_TMemcpy(c.values.data(), hypre_CSRMatrixData(block),
                      HYPRE_Complex, nnz, HYPRE_MEMORY_HOST, location);
      }
    }
    return c;
  }  // end of copyBlock

  /*!
   * \return the global numbers of the true degrees of freedom of the local
   * degrees of freedom of a finite element space, read in a copy of the
   * prolongation, or an empty vector if a row of the prolongation has not a
   * single entry equal to one (on this process)
   * \param[in] fes: finite element space
   */
  static std::vector<HYPRE_BigInt> getGlobalTrueDegreesOfFreedom(
      const mfem::ParFiniteElementSpace& fes) {
    const hypre_ParCSRMatrix* const p = *(fes.Dof_TrueDof_Matrix());
    const auto diag = copyBlock(hypre_ParCSRMatrixDiag(p), true);
    const auto offd = copyBlock(hypre_ParCSRMatrixOffd(p), true);
    const auto* const cmap = hypre_ParCSRMatrixColMapOffd(p);
    const auto first = hypre_ParCSRMatrixFirstColDiag(p);
    const auto n = diag.rows.size() - 1;
    auto T = std::vector<HYPRE_BigInt>(n);
    for (std::size_t r = 0; r != n; ++r) {
      const auto d = static_cast<std::size_t>(diag.rows[r]);
      const auto o = static_cast<std::size_t>(offd.rows[r]);
      const auto nd = diag.rows[r + 1] - diag.rows[r];
      const auto no = offd.rows[r + 1] - offd.rows[r];
      if ((nd == 1) && (no == 0) && isOne(diag.values[d])) {
        T[r] = first + diag.columns[d];
      } else if ((nd == 0) && (no == 1) && isOne(offd.values[o])) {
        T[r] = cmap[offd.columns[o]];
      } else {
        return {};
      }
    }
    return T;
  }  // end of getGlobalTrueDegreesOfFreedom

  /*!
   * \brief positions of the entries of the part of a global matrix owned by
   * the process, searched row by row: the positions of the entries of the
   * last row asked are kept by column, so that the entries of a row of a
   * local matrix, which follow each other and belong to the same row of
   * the global matrix, are found without searching this row again
   */
  struct GlobalPositionsFinder {
    /*!
     * \brief constructor
     * \param[in] m: global matrix, which must outlive this object. The
     * patterns of its blocks are copied to the host (see `copyBlock`).
     */
    explicit GlobalPositionsFinder(const hypre_ParCSRMatrix* const m)
        : A(m),
          diag(copyBlock(hypre_ParCSRMatrixDiag(m), false)),
          offd(copyBlock(hypre_ParCSRMatrixOffd(m), false)),
          first_row(hypre_ParCSRMatrixFirstRowIndex(m)),
          last_row(hypre_ParCSRMatrixLastRowIndex(m)),
          first_column(hypre_ParCSRMatrixFirstColDiag(m)),
          last_column(hypre_ParCSRMatrixLastColDiag(m)),
          diagonal_block(static_cast<std::size_t>(
              hypre_CSRMatrixNumCols(hypre_ParCSRMatrixDiag(m)))),
          off_diagonal_block(static_cast<std::size_t>(
              hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(m)))) {
    }  // end of GlobalPositionsFinder
    /*!
     * \return the position of the entry `(I, J)` (global numbers) in the
     * values of the part of the global matrix owned by the process
     * (positions in the diagonal block, then in the off diagonal block from
     * the number of entries of the diagonal block), or -1 if the entry is
     * not in the pattern of the matrix or if its row is not owned by the
     * process
     * \param[in] I: row
     * \param[in] J: column
     */
    [[nodiscard]] int find(const HYPRE_BigInt I, const HYPRE_BigInt J) {
      if ((I < this->first_row) || (I > this->last_row)) {
        return -1;
      }
      const auto r = static_cast<HYPRE_Int>(I - this->first_row);
      if (r != this->row) {
        this->row = r;
        this->setRow(this->diagonal_block, this->diag, 0);
        this->setRow(this->off_diagonal_block, this->offd,
                     static_cast<HYPRE_Int>(this->diag.columns.size()));
      }
      if ((J >= this->first_column) && (J <= this->last_column)) {
        return this->getPosition(
            this->diagonal_block,
            static_cast<std::size_t>(J - this->first_column));
      }
      const auto* const cmap = hypre_ParCSRMatrixColMapOffd(this->A);
      const auto nc = this->off_diagonal_block.size();
      const auto* const c = std::lower_bound(cmap, cmap + nc, J);
      if ((c == cmap + nc) || (*c != J)) {
        return -1;
      }
      return this->getPosition(this->off_diagonal_block,
                               static_cast<std::size_t>(c - cmap));
    }  // end of find

   private:
    //! \brief position of the entry of a column in a row of a block
    struct Position {
      //! \brief row of the block in which the entry was found
      HYPRE_Int row = -1;
      //! \brief position of the entry in the values of the global matrix
      int position = -1;
    };
    /*!
     * \brief keeps the positions of the entries of the current row of a
     * block
     * \param[out] positions: positions by column of the block
     * \param[in] block: pattern of the block
     * \param[in] offset: position of the first value of the block
     */
    void setRow(std::vector<Position>& positions,
                const BlockCopy& block,
                const HYPRE_Int offset) const {
      const auto r = static_cast<std::size_t>(this->row);
      for (auto p = block.rows[r]; p != block.rows[r + 1]; ++p) {
        positions[static_cast<std::size_t>(
            block.columns[static_cast<std::size_t>(p)])] = {
            this->row, static_cast<int>(offset + p)};
      }
    }  // end of setRow
    /*!
     * \return the position of the entry of a column of a block in the
     * current row, or -1
     * \param[in] positions: positions by column of the block
     * \param[in] c: column of the block
     */
    [[nodiscard]] int getPosition(const std::vector<Position>& positions,
                                  const std::size_t c) const {
      return positions[c].row == this->row ? positions[c].position : -1;
    }  // end of getPosition
    //! \brief global matrix
    const hypre_ParCSRMatrix* A;
    //! \brief pattern of the diagonal block
    BlockCopy diag;
    //! \brief pattern of the off diagonal block
    BlockCopy offd;
    //! \brief first row owned by the process
    HYPRE_BigInt first_row;
    //! \brief last row owned by the process
    HYPRE_BigInt last_row;
    //! \brief first column of the diagonal block
    HYPRE_BigInt first_column;
    //! \brief last column of the diagonal block
    HYPRE_BigInt last_column;
    //! \brief current row, whose positions are kept
    HYPRE_Int row = -1;
    //! \brief positions of the entries of the diagonal block, by column
    std::vector<Position> diagonal_block;
    //! \brief positions of the entries of the off diagonal block, by column
    std::vector<Position> off_diagonal_block;
  };  // end of struct GlobalPositionsFinder

  /*!
   * \brief frees the copies of the values of a block of a matrix that hypre
   * keeps in it when it computes on a GPU: the block with sorted columns
   * and the one of the triangular solves (Gauss-Seidel relaxations), which
   * hypre only builds if they do not exist. Without this call, hypre would
   * use the values of the previous assembly.
   *
   * These fields are internal to hypre. They were listed by reading hypre
   * 3.2.0 (`seq_mv/csr_matrix.h`, `parcsr_mv/par_csr_matrix.h`), to be read
   * again with another version: the other data kept in a matrix are its
   * transposed blocks (`diagT`, `offdT`: see `ResetTranspose` in
   * `assemble`), the description of its pattern for the sparse library of
   * the GPU (`mat_data`, without values), and two fields that the matrices
   * assembled here do not have (`bdiaginv`, inverses of the diagonal blocks
   * of the multigrid reduction, and `proc_ordering`, order of the ordered
   * Gauss-Seidel relaxation). `GlobalStiffnessMatrixTest` applies products
   * and smoothers of hypre after a second assembly to detect such a copy.
   * \param[in, out] m: block of a matrix of hypre
   */
  static void freeValuesKeptByHypre([[maybe_unused]] hypre_CSRMatrix* const m) {
#if defined(HYPRE_USING_CUSPARSE) || defined(HYPRE_USING_ROCSPARSE) || \
    defined(HYPRE_USING_ONEMKLSPARSE)
    hypre_TFree(hypre_CSRMatrixSortedData(m), HYPRE_MEMORY_DEVICE);
    hypre_TFree(hypre_CSRMatrixSortedJ(m), HYPRE_MEMORY_DEVICE);
    hypre_CsrsvDataDestroy(hypre_CSRMatrixCsrsvData(m));
    hypre_CSRMatrixCsrsvData(m) = nullptr;
#endif
  }  // end of freeValuesKeptByHypre

  /*!
   * \return the pattern of the global matrix of a local matrix: the product
   * `Pᵀ A P` of hypre on a matrix of ones, the sums of positive values
   * being never null. The matrix of ones and the operands of the product
   * are freed before the tables of the setup are built.
   * \param[in] fes: finite element space
   * \param[in] A: local matrix
   */
  static std::unique_ptr<mfem::HypreParMatrix> makeGlobalPattern(
      const mfem::ParFiniteElementSpace& fes, const mfem::SparseMatrix& A) {
    auto ones = mfem::SparseMatrix(A);
    ones = real{1};
    auto dA = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto Ph = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto G = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    dA.MakeSquareBlockDiag(fes.GetComm(), fes.GlobalVSize(),
                           fes.GetDofOffsets(), &ones);
    Ph.ConvertFrom(fes.Dof_TrueDof_Matrix());
    G.MakePtAP(dA, Ph);
    G.SetOperatorOwner(false);
    return std::unique_ptr<mfem::HypreParMatrix>(G.As<mfem::HypreParMatrix>());
  }  // end of makeGlobalPattern

  /*!
   * \brief copies indices in an array read by the kernels. When the kernels
   * run on a GPU, the array is allocated for the device and the indices are
   * copied to it directly: an array written on the host by MFEM then read
   * on the device keeps both copies, twice its size on a GPU which shares
   * the memory of the host.
   * \param[out] a: array
   * \param[in] v: indices
   * \param[in] use_device: if true, the kernels run on the device
   */
  static void setIndices(mfem::Array<int>& a,
                         const std::vector<int>& v,
                         const bool use_device) {
    const auto n = static_cast<int>(v.size());
    if (!use_device) {
      a.SetSize(n);
      std::copy(v.begin(), v.end(), a.HostWrite());
      return;
    }
    a.SetSize(n, mfem::Device::GetDeviceMemoryType());
    if (n != 0) {
      a.GetMemory().CopyFromHost(v.data(), n);
    }
  }  // end of setIndices

  /*!
   * \brief builds the lists of the contributions to the entries of the
   * global matrix: for each entry, the entries of the local matrix, then
   * the received values (numbered from the number of entries of the local
   * matrix), in increasing order
   * \param[out] offsets: offsets of the contributions to each entry
   * \param[out] contributions: contributions
   * \param[in] positions: positions in the global matrix of the entries of
   * the local matrix, or -1
   * \param[in] received_positions: positions in the global matrix of the
   * received values, or -1
   * \param[in] ng: number of entries of the global matrix
   * \param[in] use_device: if true, the kernels run on the device
   */
  static void setContributions(mfem::Array<int>& offsets,
                               mfem::Array<int>& contributions,
                               const std::vector<int>& positions,
                               const std::vector<int>& received_positions,
                               const std::size_t ng,
                               const bool use_device) {
    auto o = std::vector<int>(ng + 1, 0);
    for (const auto* const ps : {&positions, &received_positions}) {
      for (const auto p : *ps) {
        if (p >= 0) {
          ++o[static_cast<std::size_t>(p) + 1];
        }
      }
    }
    for (std::size_t p = 0; p != ng; ++p) {
      o[p + 1] += o[p];
    }
    auto c = std::vector<int>(static_cast<std::size_t>(o[ng]));
    auto next = std::vector<int>(o.begin(), o.end() - 1);
    auto contribution = 0;
    for (const auto* const ps : {&positions, &received_positions}) {
      for (const auto p : *ps) {
        if (p >= 0) {
          const auto t = next[static_cast<std::size_t>(p)]++;
          c[static_cast<std::size_t>(t)] = contribution;
        }
        ++contribution;
      }
    }
    setIndices(offsets, o, use_device);
    setIndices(contributions, c, use_device);
  }  // end of setContributions

  bool GlobalStiffnessMatrix::isSupported(
      const mfem::ParFiniteElementSpace& fes) noexcept {
    auto supported = int{0};
    try {
      // the construction of the prolongation is collective: it is built on
      // all the processes, including the ones without degree of freedom
      static_cast<void>(fes.Dof_TrueDof_Matrix());
      supported = ((fes.GetVSize() == 0) ||
                   (!getGlobalTrueDegreesOfFreedom(fes).empty()))
                      ? 1
                      : 0;
    } catch (...) {
      supported = 0;
    }
    auto all = int{0};
    MPI_Allreduce(&supported, &all, 1, MPI_INT, MPI_MIN, fes.GetComm());
    return all == 1;
  }  // end of isSupported

  GlobalStiffnessMatrix::GlobalStiffnessMatrix(
      const mfem::ParFiniteElementSpace& f)
      : fes(&f) {}  // end of GlobalStiffnessMatrix

  GlobalStiffnessMatrix::GlobalStiffnessMatrix(
      GlobalStiffnessMatrix&&) noexcept = default;

  GlobalStiffnessMatrix& GlobalStiffnessMatrix::operator=(
      GlobalStiffnessMatrix&&) noexcept = default;

  GlobalStiffnessMatrix::~GlobalStiffnessMatrix() = default;

  bool GlobalStiffnessMatrix::setup(Context& ctx, const mfem::SparseMatrix& A) {
    const auto comm = this->fes->GetComm();
    auto nprocs = int{};
    auto rank = int{};
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
    // checks made before the collective calls, whose result is shared by
    // all the processes: no process returns while the others wait
    const auto T = getGlobalTrueDegreesOfFreedom(*(this->fes));
    const auto valid = (A.Finalized()) &&
                       (A.Height() == this->fes->GetVSize()) &&
                       (A.Width() == this->fes->GetVSize()) &&
                       (T.size() == static_cast<std::size_t>(A.Height()));
    auto all_valid = int{0};
    const auto local_valid = valid ? 1 : 0;
    MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN, comm);
    if (all_valid == 0) {
      return ctx.registerErrorMessage(
          valid ? "GlobalStiffnessMatrix: failure on another process"
                : "GlobalStiffnessMatrix: unsupported finite element space or "
                  "local matrix of a wrong size");
    }
    // from here, all the processes take part in all the collective calls;
    // an entry not found is reported at the end
    auto found = true;
    // the tables read by the kernels are located where hypre computes
    const auto use_device = isHypreUsingGPU();
    this->matrix = makeGlobalPattern(*(this->fes), A);
    const hypre_ParCSRMatrix* const g = *(this->matrix);
    const auto nd = hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixDiag(g));
    const auto no = hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixOffd(g));
    auto too_many_entries =
        static_cast<long long>(nd) + no >= std::numeric_limits<int>::max();
    this->diagonal_block_entries = static_cast<int>(nd);
    this->off_diagonal_block_entries = static_cast<int>(no);
    // first row owned by each process
    auto first_rows =
        std::vector<HYPRE_BigInt>(static_cast<std::size_t>(nprocs) + 1);
    const auto first_row = hypre_ParCSRMatrixFirstRowIndex(g);
    MPI_Allgather(&first_row, 1, HYPRE_MPI_BIG_INT, first_rows.data(), 1,
                  HYPRE_MPI_BIG_INT, comm);
    first_rows.back() = hypre_ParCSRMatrixGlobalNumRows(g);
    const auto owner = [&first_rows](const HYPRE_BigInt I) {
      const auto it =
          std::upper_bound(first_rows.begin(), first_rows.end() - 1, I);
      return static_cast<int>(it - first_rows.begin()) - 1;
    };
    // destinations of the entries of the local matrix
    auto finder = GlobalPositionsFinder(g);
    const auto* const Ai = A.HostReadI();
    const auto* const Aj = A.HostReadJ();
    const auto n = A.NumNonZeroElems();
    // positions of the entries of the local matrix whose rows are owned by
    // the process in the values of the global matrix: positions in the
    // diagonal block, then, from its number of entries, in the off diagonal
    // block; -1 for the entries sent to other processes
    auto positions = std::vector<int>(static_cast<std::size_t>(n), -1);
    auto sent_to =
        std::vector<std::vector<int>>(static_cast<std::size_t>(nprocs));
    auto pairs = std::vector<std::vector<HYPRE_BigInt>>(
        static_cast<std::size_t>(nprocs));
    for (int r = 0; r != A.Height(); ++r) {
      const auto I = T[static_cast<std::size_t>(r)];
      const auto o = owner(I);
      for (auto k = Ai[r]; k != Ai[r + 1]; ++k) {
        const auto J = T[static_cast<std::size_t>(Aj[k])];
        if (o == rank) {
          const auto p = finder.find(I, J);
          found = found && (p >= 0);
          positions[static_cast<std::size_t>(k)] = p;
        } else {
          sent_to[static_cast<std::size_t>(o)].push_back(k);
          pairs[static_cast<std::size_t>(o)].push_back(I);
          pairs[static_cast<std::size_t>(o)].push_back(J);
        }
      }
    }
    // exchanges: the entries sent, ordered by process
    const auto np = static_cast<std::size_t>(nprocs);
    auto send_counts = std::vector<int>(np);
    auto send_offsets = std::vector<int>(np);
    auto sent_entries_values = std::vector<int>{};
    auto sent_pairs = std::vector<HYPRE_BigInt>{};
    for (std::size_t q = 0; q != np; ++q) {
      send_offsets[q] = static_cast<int>(sent_entries_values.size());
      send_counts[q] = static_cast<int>(sent_to[q].size());
      sent_entries_values.insert(sent_entries_values.end(), sent_to[q].begin(),
                                 sent_to[q].end());
      sent_pairs.insert(sent_pairs.end(), pairs[q].begin(), pairs[q].end());
    }
    const auto sent = static_cast<int>(sent_entries_values.size());
    setIndices(this->sent_entries, sent_entries_values, use_device);
    auto receive_counts = std::vector<int>(np);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, receive_counts.data(), 1,
                 MPI_INT, comm);
    auto receive_offsets = std::vector<int>(np);
    auto received = 0;
    for (std::size_t q = 0; q != np; ++q) {
      receive_offsets[q] = received;
      received += receive_counts[q];
    }
    // messages of the assemblies: the processes with which values are
    // exchanged
    this->sent_messages.clear();
    this->received_messages.clear();
    for (std::size_t q = 0; q != np; ++q) {
      if (send_counts[q] != 0) {
        this->sent_messages.push_back(
            {static_cast<int>(q), send_offsets[q], send_counts[q]});
      }
      if (receive_counts[q] != 0) {
        this->received_messages.push_back(
            {static_cast<int>(q), receive_offsets[q], receive_counts[q]});
      }
    }
    // global numbers of the entries sent: one pair per entry, so that the
    // counts and the offsets of the messages are the ones of the entries
    auto received_pairs =
        std::vector<HYPRE_BigInt>(2 * static_cast<std::size_t>(received));
    auto pair_type = MPI_Datatype{};
    MPI_Type_contiguous(2, HYPRE_MPI_BIG_INT, &pair_type);
    MPI_Type_commit(&pair_type);
    MPI_Alltoallv(sent_pairs.data(), send_counts.data(), send_offsets.data(),
                  pair_type, received_pairs.data(), receive_counts.data(),
                  receive_offsets.data(), pair_type, comm);
    MPI_Type_free(&pair_type);
    // positions of the received entries in the values of the global matrix
    auto received_positions =
        std::vector<int>(static_cast<std::size_t>(received));
    for (std::size_t i = 0; i != static_cast<std::size_t>(received); ++i) {
      const auto p =
          finder.find(received_pairs[2 * i], received_pairs[2 * i + 1]);
      found = found && (p >= 0);
      received_positions[i] = p;
    }
    // a contribution is the index of an entry of the local matrix, or the
    // one of a received value from the number of entries of the local matrix
    too_many_entries =
        too_many_entries || (static_cast<long long>(n) + received >
                             std::numeric_limits<int>::max());
    if (!too_many_entries) {
      setContributions(
          this->offsets, this->contributions, positions, received_positions,
          static_cast<std::size_t>(nd) + static_cast<std::size_t>(no),
          use_device);
    }
    this->send_buffer.UseDevice(true);
    this->send_buffer.SetSize(sent);
    this->receive_buffer.UseDevice(true);
    this->receive_buffer.SetSize(received);
    this->requests.reserve(this->sent_messages.size() +
                           this->received_messages.size());
    this->number_of_rows = A.Height();
    this->number_of_entries = n;
    if (too_many_entries) {
      return ctx.registerErrorMessage(
          "GlobalStiffnessMatrix: too many entries in the global matrix");
    }
    if (!found) {
      return ctx.registerErrorMessage(
          "GlobalStiffnessMatrix: entry not found in the global matrix");
    }
    return true;
  }  // end of setup

  mfem::HypreParMatrix* GlobalStiffnessMatrix::assemble(
      Context& ctx, const mfem::SparseMatrix& A) noexcept {
    try {
      const auto comm = this->fes->GetComm();
      // a new setup is made by all the processes if the local matrix or the
      // finite element space changed on one of them (the arrays of the
      // pattern of the local matrix are compared, not their contents)
      const auto* const Ai = A.HostReadI();
      const auto* const Aj = A.HostReadJ();
      const auto changed = (this->matrix == nullptr) ||
                           (A.Height() != this->number_of_rows) ||
                           (A.NumNonZeroElems() != this->number_of_entries) ||
                           (Ai != this->rows) || (Aj != this->columns) ||
                           (this->fes->GetSequence() != this->sequence);
      auto local_changed = changed ? 1 : 0;
      auto any_changed = int{0};
      MPI_Allreduce(&local_changed, &any_changed, 1, MPI_INT, MPI_MAX, comm);
      if (any_changed == 1) {
        // collective: all the processes take part in the setup, even if it
        // fails on one of them, and all of them forget the previous one
        const auto success = this->setup(ctx, A) ? 1 : 0;
        auto all = int{0};
        MPI_Allreduce(&success, &all, 1, MPI_INT, MPI_MIN, comm);
        if (all == 0) {
          this->matrix.reset();
          if (success == 1) {
            std::ignore = ctx.registerErrorMessage(
                "GlobalStiffnessMatrix::assemble: failure on another "
                "process");
          }
          return nullptr;
        }
        this->rows = Ai;
        this->columns = Aj;
        this->sequence = this->fes->GetSequence();
      }
      // the values are computed where hypre computes
      const auto use_device = isHypreUsingGPU();
      this->matrix->HypreReadWrite();
      hypre_ParCSRMatrix* const g = *(this->matrix);
      auto* const diag = hypre_ParCSRMatrixDiag(g);
      auto* const offd = hypre_ParCSRMatrixOffd(g);
      const auto* const a = A.ReadData(use_device);
      const auto sent = this->send_buffer.Size();
      const auto received = this->receive_buffer.Size();
      // exchange of the values of the rows owned by other processes, given
      // to MPI on the device with GPU-aware MPI, the kernel filling the
      // values sent being finished (as MFEM does in `GroupCommunicator`)
      const auto mpi_on_device = use_device && mfem::Device::GetGPUAwareMPI();
      const auto* s = static_cast<const real*>(nullptr);
      if (sent != 0) {
        internals::launchGlobalStiffnessMatrixSentValuesKernel(
            this->send_buffer.Write(use_device),
            this->sent_entries.Read(use_device), a, sent, use_device);
        if (mpi_on_device) {
          MFEM_STREAM_SYNC;
        }
        s = mpi_on_device ? this->send_buffer.Read()
                          : this->send_buffer.HostRead();
      }
      auto* r = static_cast<real*>(nullptr);
      if (received != 0) {
        r = mpi_on_device ? this->receive_buffer.Write()
                          : this->receive_buffer.HostWrite();
      }
      const auto mpi_real = mfem::MPITypeMap<real>::mpi_type;
      this->requests.clear();
      for (const auto& m : this->received_messages) {
        MPI_Irecv(r + m.offset, m.size, mpi_real, m.process, message_tag, comm,
                  &(this->requests.emplace_back()));
      }
      for (const auto& m : this->sent_messages) {
        MPI_Isend(s + m.offset, m.size, mpi_real, m.process, message_tag, comm,
                  &(this->requests.emplace_back()));
      }
      MPI_Waitall(static_cast<int>(this->requests.size()),
                  this->requests.data(), MPI_STATUSES_IGNORE);
      auto k = internals::GlobalStiffnessMatrixKernelArguments{};
      k.offsets = this->offsets.Read(use_device);
      k.contributions = this->contributions.Read(use_device);
      k.local_values = a;
      if (received != 0) {
        k.received_values = this->receive_buffer.Read(use_device);
      }
      k.diagonal_block = hypre_CSRMatrixData(diag);
      k.off_diagonal_block = hypre_CSRMatrixData(offd);
      k.local_entries = this->number_of_entries;
      k.diagonal_block_entries = this->diagonal_block_entries;
      k.off_diagonal_block_entries = this->off_diagonal_block_entries;
      internals::launchGlobalStiffnessMatrixKernel(k, use_device);
      // copies of the former values kept in the matrix: its transposed
      // blocks, built by MFEM for the transposed products on a GPU, and
      // the copies of hypre
      this->matrix->ResetTranspose();
      if (use_device) {
        freeValuesKeptByHypre(diag);
        freeValuesKeptByHypre(offd);
      }
      return this->matrix.get();
    } catch (std::exception& e) {
      std::ignore = ctx.registerErrorMessage(
          std::string{"GlobalStiffnessMatrix::assemble: "} + e.what());
    } catch (...) {
      std::ignore = ctx.registerErrorMessage(
          "GlobalStiffnessMatrix::assemble: unknown exception");
    }
    return nullptr;
  }  // end of assemble

}  // end of namespace mfem_mgis::batch

#endif /* MFEM_USE_MPI */
