/*!
 * \file   src/Batch/GlobalStiffnessMatrix.cxx
 * \brief  This file implements the assembly of the global matrix of a
 * parallel finite element space from the local matrices of the processes
 * \date   11/10/2026
 */

#include "MFEMMGIS/Batch/GlobalStiffnessMatrix.hxx"

#ifdef MFEM_USE_MPI

#include <tuple>
#include <atomic>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <exception>
#include "mfem/general/forall.hpp"
#include "mfem/general/communication.hpp"
#include "mfem/linalg/handle.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/pfespace.hpp"
#include "HostThreads.hxx"
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
    std::unique_ptr<HYPRE_Int[]> columns;
    //! \brief values of the entries, if they were asked
    std::unique_ptr<HYPRE_Complex[]> values;
    //! \brief number of entries
    HYPRE_Int size = 0;
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
      c.size = nnz;
      c.columns = std::make_unique_for_overwrite<HYPRE_Int[]>(
          static_cast<std::size_t>(nnz));
      // the array is first written by the threads: the copy from a GPU
      // takes several times longer when it also maps its pages
      internals::forEachThreadRange(
          nnz, [&c](const HYPRE_Int begin, const HYPRE_Int end) {
            std::fill(c.columns.get() + begin, c.columns.get() + end, 0);
          });
      hypre_TMemcpy(c.columns.get(), hypre_CSRMatrixJ(block), HYPRE_Int, nnz,
                    HYPRE_MEMORY_HOST, location);
      if (values) {
        c.values = std::make_unique_for_overwrite<HYPRE_Complex[]>(
            static_cast<std::size_t>(nnz));
        hypre_TMemcpy(c.values.get(), hypre_CSRMatrixData(block), HYPRE_Complex,
                      nnz, HYPRE_MEMORY_HOST, location);
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
   * \brief copy on the host of the pattern of the part of a global matrix
   * owned by the process, in which the positions of its entries are
   * searched (see `GlobalPositionsFinder`)
   */
  struct GlobalPattern {
    /*!
     * \brief constructor
     * \param[in] m: global matrix, which must outlive this object. The
     * patterns of its blocks are copied to the host (see `copyBlock`).
     */
    explicit GlobalPattern(const hypre_ParCSRMatrix* const m)
        : diag(copyBlock(hypre_ParCSRMatrixDiag(m), false)),
          offd(copyBlock(hypre_ParCSRMatrixOffd(m), false)),
          column_map(hypre_ParCSRMatrixColMapOffd(m)),
          first_row(hypre_ParCSRMatrixFirstRowIndex(m)),
          last_row(hypre_ParCSRMatrixLastRowIndex(m)),
          first_column(hypre_ParCSRMatrixFirstColDiag(m)),
          last_column(hypre_ParCSRMatrixLastColDiag(m)),
          diagonal_block_columns(static_cast<std::size_t>(
              hypre_CSRMatrixNumCols(hypre_ParCSRMatrixDiag(m)))),
          off_diagonal_block_columns(static_cast<std::size_t>(
              hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(m)))) {
    }  // end of GlobalPattern
    //! \brief pattern of the diagonal block
    const BlockCopy diag;
    //! \brief pattern of the off diagonal block
    const BlockCopy offd;
    //! \brief global numbers of the columns of the off diagonal block
    const HYPRE_BigInt* const column_map;
    //! \brief first row owned by the process
    const HYPRE_BigInt first_row;
    //! \brief last row owned by the process
    const HYPRE_BigInt last_row;
    //! \brief first column of the diagonal block
    const HYPRE_BigInt first_column;
    //! \brief last column of the diagonal block
    const HYPRE_BigInt last_column;
    //! \brief number of columns of the diagonal block
    const std::size_t diagonal_block_columns;
    //! \brief number of columns of the off diagonal block
    const std::size_t off_diagonal_block_columns;
  };  // end of struct GlobalPattern

  /*!
   * \brief positions of the entries of the part of a global matrix owned by
   * the process, searched row by row: the entries of the selected row are
   * kept by column, so that the entries of a row of a local matrix, which
   * follow each other and belong to the same row of the global matrix, are
   * found without searching this row again. The pattern is shared: each
   * thread searching entries has its own finder.
   */
  struct GlobalPositionsFinder {
    /*!
     * \brief constructor
     * \param[in] p: pattern of the global matrix, which must outlive this
     * object
     */
    explicit GlobalPositionsFinder(const GlobalPattern& p)
        : pattern(p),
          diagonal_block(p.diagonal_block_columns, -1),
          off_diagonal_block(p.off_diagonal_block_columns, -1) {
    }  // end of GlobalPositionsFinder
    /*!
     * \brief selects the row of the entries searched by `find`
     * \param[in] I: row (global number)
     */
    void setRow(const HYPRE_BigInt I) {
      if (I == this->row) {
        return;
      }
      const auto& p = this->pattern;
      this->row = I;
      if ((I < p.first_row) || (I > p.last_row)) {
        this->diagonal = {};
        this->off_diagonal = {};
        return;
      }
      const auto r = static_cast<std::size_t>(I - p.first_row);
      this->diagonal = keepEntries(this->diagonal_block, p.diag, r);
      this->off_diagonal = keepEntries(this->off_diagonal_block, p.offd, r);
    }  // end of setRow
    /*!
     * \return the position of the entry of a column of the selected row in
     * the values of the part of the global matrix owned by the process
     * (positions in the diagonal block, then in the off diagonal block from
     * the number of entries of the diagonal block), or -1 if the entry is
     * not in the pattern of the matrix or if the row is not owned by the
     * process
     * \param[in] J: column (global number)
     */
    [[nodiscard]] int find(const HYPRE_BigInt J) const {
      const auto& p = this->pattern;
      if ((J >= p.first_column) && (J <= p.last_column)) {
        const auto c = static_cast<std::size_t>(J - p.first_column);
        return getPosition(this->diagonal_block[c], this->diagonal, 0);
      }
      const auto* const cmap = p.column_map;
      const auto nc = p.off_diagonal_block_columns;
      const auto* const c = std::lower_bound(cmap, cmap + nc, J);
      if ((c == cmap + nc) || (*c != J)) {
        return -1;
      }
      return getPosition(
          this->off_diagonal_block[static_cast<std::size_t>(c - cmap)],
          this->off_diagonal, p.diag.size);
    }  // end of find

   private:
    //! \brief entries of a row of a block
    struct Row {
      //! \brief first entry
      HYPRE_Int begin = 0;
      //! \brief entry following the last one
      HYPRE_Int end = 0;
    };
    /*!
     * \brief keeps the entries of a row of a block
     * \return the entries of the row
     * \param[out] entries: last entry met in each column of the block
     * \param[in] block: pattern of the block
     * \param[in] r: row
     */
    static Row keepEntries(std::vector<HYPRE_Int>& entries,
                           const BlockCopy& block,
                           const std::size_t r) {
      const auto row = Row{block.rows[r], block.rows[r + 1]};
      for (auto e = row.begin; e != row.end; ++e) {
        entries[static_cast<std::size_t>(
            block.columns[static_cast<std::size_t>(e)])] = e;
      }
      return row;
    }  // end of keepEntries
    /*!
     * \return the position of an entry kept for a column if it belongs to
     * the selected row, or -1: the entry kept for a column is the last one
     * met in this column, in any row
     * \param[in] e: entry kept for the column
     * \param[in] row: entries of the selected row
     * \param[in] offset: position of the first value of the block
     */
    [[nodiscard]] static int getPosition(const HYPRE_Int e,
                                         const Row& row,
                                         const HYPRE_Int offset) {
      return ((e >= row.begin) && (e < row.end)) ? static_cast<int>(offset + e)
                                                 : -1;
    }  // end of getPosition
    //! \brief pattern of the global matrix
    const GlobalPattern& pattern;
    //! \brief selected row (global number), none at first
    HYPRE_BigInt row = -1;
    //! \brief entries of the selected row in the diagonal block
    Row diagonal;
    //! \brief entries of the selected row in the off diagonal block
    Row off_diagonal;
    /*!
     * \brief last entry met in each column of the diagonal block (index in
     * the block), or -1
     */
    std::vector<HYPRE_Int> diagonal_block;
    //! \brief last entry met in each column of the off diagonal block, or -1
    std::vector<HYPRE_Int> off_diagonal_block;
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
    // matrix of ones of the pattern of `A`, built where hypre computes: on a
    // GPU, its arrays are copied from the ones of `A` on the GPU, which the
    // kernels of the assembly use, and its values are set by a kernel,
    // without going through the host
    const auto use_device = isHypreUsingGPU();
    auto ones =
        mfem::SparseMatrix(A, true,
                           use_device ? mfem::Device::GetDeviceMemoryType()
                                      : mfem::MemoryType::PRESERVE);
    internals::launchGlobalStiffnessMatrixOnesKernel(
        ones.WriteData(use_device), ones.NumNonZeroElems(), use_device);
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
   * \param[in] n: number of indices
   * \param[in] use_device: if true, the kernels run on the device
   */
  static void setIndices(mfem::Array<int>& a,
                         const int* const v,
                         const int n,
                         const bool use_device) {
    if (!use_device) {
      a.SetSize(n);
      std::copy(v, v + n, a.HostWrite());
      return;
    }
    a.SetSize(n, mfem::Device::GetDeviceMemoryType());
    if (n != 0) {
      a.GetMemory().CopyFromHost(v, n);
    }
  }  // end of setIndices

  /*!
   * \brief rows of a local matrix added to each row of the part of the
   * global matrix owned by the process, in increasing order. The entries of
   * two rows of the global matrix have no common contribution from the
   * local matrix: the threads treat different rows of the global matrix
   * (see `forEachThreadRange`).
   */
  struct LocalRows {
    /*!
     * \brief constructor
     * \param[in] T: global numbers of the true degrees of freedom of the
     * rows of the local matrix
     * \param[in] p: pattern of the global matrix
     */
    LocalRows(const std::vector<HYPRE_BigInt>& T, const GlobalPattern& p)
        : offsets(p.diag.rows.size(), 0) {
      const auto owned = [&p](const HYPRE_BigInt I) {
        return (I >= p.first_row) && (I <= p.last_row);
      };
      for (const auto I : T) {
        if (owned(I)) {
          ++(this->offsets[static_cast<std::size_t>(I - p.first_row) + 1]);
        }
      }
      std::partial_sum(this->offsets.begin(), this->offsets.end(),
                       this->offsets.begin());
      this->rows.resize(static_cast<std::size_t>(this->offsets.back()));
      auto next =
          std::vector<int>(this->offsets.begin(), this->offsets.end() - 1);
      for (std::size_t r = 0; r != T.size(); ++r) {
        if (owned(T[r])) {
          const auto i = static_cast<std::size_t>(T[r] - p.first_row);
          this->rows[static_cast<std::size_t>(next[i]++)] = static_cast<int>(r);
        }
      }
    }  // end of LocalRows
    //! \return the number of rows of the global matrix
    [[nodiscard]] int size() const {
      return static_cast<int>(this->offsets.size()) - 1;
    }
    //! \brief offsets of the local rows of each row of the global matrix
    std::vector<int> offsets;
    //! \brief local rows
    std::vector<int> rows;
  };  // end of struct LocalRows

  /*!
   * \brief builds the lists of the contributions to the entries of the
   * global matrix: for each entry, the entries of the local matrix, then
   * the received values (numbered from the number of entries of the local
   * matrix), in increasing order
   * \param[out] offsets: offsets of the contributions to each entry
   * \param[out] contributions: contributions
   * \param[in] rows: rows of the local matrix of each row of the global
   * matrix
   * \param[in] Ai: offsets of the rows of the local matrix
   * \param[in] positions: positions in the global matrix of the entries of
   * the local matrix, or -1
   * \param[in] received_positions: positions in the global matrix of the
   * received values, or -1
   * \param[in] n: number of entries of the local matrix
   * \param[in] ng: number of entries of the global matrix
   * \param[in] use_device: if true, the kernels run on the device
   */
  static void setContributions(mfem::Array<int>& offsets,
                               mfem::Array<int>& contributions,
                               const LocalRows& rows,
                               const int* const Ai,
                               const int* const positions,
                               const std::vector<int>& received_positions,
                               const int n,
                               const std::size_t ng,
                               const bool use_device) {
    // calls `f(k, p)` for each entry `k` of the local matrix added to an
    // entry `p` of the global matrix, in increasing order for each `p`
    const auto forEachLocalContribution = [&rows, Ai,
                                           positions](const auto& f) {
      internals::forEachThreadRange(
          rows.size(),
          [&rows, Ai, positions, &f](const int begin, const int end) {
            for (auto i = rows.offsets[static_cast<std::size_t>(begin)];
                 i != rows.offsets[static_cast<std::size_t>(end)]; ++i) {
              const auto r = rows.rows[static_cast<std::size_t>(i)];
              for (auto k = Ai[r]; k != Ai[r + 1]; ++k) {
                if (positions[k] >= 0) {
                  f(k, static_cast<std::size_t>(positions[k]));
                }
              }
            }
          });
    };
    // number of contributions to each entry, kept in the next offset
    auto o = std::make_unique_for_overwrite<int[]>(ng + 1);
    internals::forEachThreadRange(
        ng + 1, [&o](const std::size_t begin, const std::size_t end) {
          std::fill(o.get() + begin, o.get() + end, 0);
        });
    forEachLocalContribution(
        [&o](const int, const std::size_t p) { ++o[p + 1]; });
    for (const auto p : received_positions) {
      if (p >= 0) {
        ++o[static_cast<std::size_t>(p) + 1];
      }
    }
    // first contribution to each entry, kept in the next offset too
    auto size = 0;
    for (std::size_t p = 0; p != ng; ++p) {
      const auto count = o[p + 1];
      o[p + 1] = size;
      size += count;
    }
    // the contributions are visited by increasing index, which sorts the
    // ones of each entry, whose next offset ends as the offset of the next
    // entry
    auto c =
        std::make_unique_for_overwrite<int[]>(static_cast<std::size_t>(size));
    forEachLocalContribution([&o, &c](const int k, const std::size_t p) {
      c[static_cast<std::size_t>(o[p + 1]++)] = k;
    });
    auto contribution = n;
    for (const auto p : received_positions) {
      if (p >= 0) {
        c[static_cast<std::size_t>(o[static_cast<std::size_t>(p) + 1]++)] =
            contribution;
      }
      ++contribution;
    }
    setIndices(offsets, o.get(), static_cast<int>(ng + 1), use_device);
    setIndices(contributions, c.get(), size, use_device);
  }  // end of setContributions

  bool GlobalStiffnessMatrix::isSupported(
      const mfem::ParFiniteElementSpace& fes) noexcept {
    // a failure, for instance an allocation of hypre on the GPU, is not an
    // unsupported space: its message is written on the standard error
    // stream, as the error stream of mfem-mgis may be disabled on this
    // process
    const auto report = [&fes](const char* const what) {
      auto rank = int{};
      MPI_Comm_rank(fes.GetComm(), &rank);
      std::cerr << "process " << rank
                << ": GlobalStiffnessMatrix::isSupported: " << what
                << std::endl;
    };
    auto supported = int{0};
    try {
      // the construction of the prolongation is collective: it is built on
      // all the processes, including the ones without degree of freedom
      static_cast<void>(fes.Dof_TrueDof_Matrix());
      supported = ((fes.GetVSize() == 0) ||
                   (!getGlobalTrueDegreesOfFreedom(fes).empty()))
                      ? 1
                      : 0;
    } catch (const std::exception& e) {
      report(e.what());
      supported = 0;
    } catch (...) {
      report("unknown exception");
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
    auto found = std::atomic<bool>{true};
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
    const auto pattern = GlobalPattern(g);
    const auto local_rows = LocalRows(T, pattern);
    const auto* const Ai = A.HostReadI();
    const auto* const Aj = A.HostReadJ();
    const auto n = A.NumNonZeroElems();
    // positions of the entries of the local matrix whose rows are owned by
    // the process in the values of the global matrix: positions in the
    // diagonal block, then, from its number of entries, in the off diagonal
    // block; -1 for the entries sent to other processes
    auto positions =
        std::make_unique_for_overwrite<int[]>(static_cast<std::size_t>(n));
    internals::forEachThreadRange(
        local_rows.size(), [&](const int begin, const int end) {
          auto finder = GlobalPositionsFinder(pattern);
          auto all = true;
          for (auto i = local_rows.offsets[static_cast<std::size_t>(begin)];
               i != local_rows.offsets[static_cast<std::size_t>(end)]; ++i) {
            const auto r = local_rows.rows[static_cast<std::size_t>(i)];
            finder.setRow(T[static_cast<std::size_t>(r)]);
            for (auto k = Ai[r]; k != Ai[r + 1]; ++k) {
              const auto p = finder.find(T[static_cast<std::size_t>(Aj[k])]);
              all = all && (p >= 0);
              positions[static_cast<std::size_t>(k)] = p;
            }
          }
          if (!all) {
            found = false;
          }
        });
    // rows owned by other processes: their entries are sent, ordered by
    // process. The place of the entries of a row among the ones sent to its
    // owner follows from the rows before it: the rows are then written by
    // the threads.
    struct SentRow {
      //! \brief row of the local matrix
      int row;
      //! \brief process owning the row
      int owner;
      //! \brief first entry of the row in the entries sent to its owner
      int offset;
    };
    const auto np = static_cast<std::size_t>(nprocs);
    auto send_counts = std::vector<int>(np, 0);
    auto sent_rows = std::vector<SentRow>{};
    for (int r = 0; r != A.Height(); ++r) {
      const auto o = owner(T[static_cast<std::size_t>(r)]);
      if (o != rank) {
        auto& count = send_counts[static_cast<std::size_t>(o)];
        sent_rows.push_back({r, o, count});
        count += Ai[r + 1] - Ai[r];
      }
    }
    auto send_offsets = std::vector<int>(np);
    std::exclusive_scan(send_counts.begin(), send_counts.end(),
                        send_offsets.begin(), 0);
    const auto sent = send_offsets.back() + send_counts.back();
    auto sent_entries_values = std::vector<int>(static_cast<std::size_t>(sent));
    auto sent_pairs =
        std::vector<HYPRE_BigInt>(2 * static_cast<std::size_t>(sent));
    internals::forEachThreadRange(
        sent_rows.size(), [&](const std::size_t begin, const std::size_t end) {
          for (auto i = begin; i != end; ++i) {
            const auto r = sent_rows[i].row;
            auto e = static_cast<std::size_t>(
                send_offsets[static_cast<std::size_t>(sent_rows[i].owner)] +
                sent_rows[i].offset);
            for (auto k = Ai[r]; k != Ai[r + 1]; ++k, ++e) {
              positions[static_cast<std::size_t>(k)] = -1;
              sent_entries_values[e] = k;
              sent_pairs[2 * e] = T[static_cast<std::size_t>(r)];
              sent_pairs[2 * e + 1] = T[static_cast<std::size_t>(Aj[k])];
            }
          }
        });
    setIndices(this->sent_entries, sent_entries_values.data(), sent,
               use_device);
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
    internals::forEachThreadRange(
        static_cast<std::size_t>(received),
        [&](const std::size_t begin, const std::size_t end) {
          auto finder = GlobalPositionsFinder(pattern);
          auto all = true;
          for (auto i = begin; i != end; ++i) {
            finder.setRow(received_pairs[2 * i]);
            const auto p = finder.find(received_pairs[2 * i + 1]);
            all = all && (p >= 0);
            received_positions[i] = p;
          }
          if (!all) {
            found = false;
          }
        });
    // a contribution is the index of an entry of the local matrix, or the
    // one of a received value from the number of entries of the local matrix
    too_many_entries =
        too_many_entries || (static_cast<long long>(n) + received >
                             std::numeric_limits<int>::max());
    if (!too_many_entries) {
      setContributions(
          this->offsets, this->contributions, local_rows, Ai, positions.get(),
          received_positions, n,
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
