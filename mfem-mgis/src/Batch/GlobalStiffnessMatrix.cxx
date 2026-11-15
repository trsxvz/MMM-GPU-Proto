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

  /*!
   * \return the global numbers of the true degrees of freedom of the local
   * degrees of freedom of a finite element space, read in the prolongation,
   * or an empty vector if a row of the prolongation has not a single entry
   * equal to one (on this process). The prolongation is read on the host,
   * then given back to hypre where it computes.
   * \param[in] fes: finite element space
   */
  static std::vector<HYPRE_BigInt> getGlobalTrueDegreesOfFreedom(
      const mfem::ParFiniteElementSpace& fes) {
    auto& P = *(fes.Dof_TrueDof_Matrix());
    P.HostRead();
    const hypre_ParCSRMatrix* const p = P;
    const auto* const diag = hypre_ParCSRMatrixDiag(p);
    const auto* const offd = hypre_ParCSRMatrixOffd(p);
    const auto* const di = hypre_CSRMatrixI(diag);
    const auto* const dj = hypre_CSRMatrixJ(diag);
    const auto* const dv = hypre_CSRMatrixData(diag);
    const auto* const oi = hypre_CSRMatrixI(offd);
    const auto* const oj = hypre_CSRMatrixJ(offd);
    const auto* const ov = hypre_CSRMatrixData(offd);
    const auto* const cmap = hypre_ParCSRMatrixColMapOffd(p);
    const auto first = hypre_ParCSRMatrixFirstColDiag(p);
    const auto n = static_cast<int>(hypre_CSRMatrixNumRows(diag));
    auto T = std::vector<HYPRE_BigInt>(static_cast<std::size_t>(n));
    for (int r = 0; r != n; ++r) {
      const auto nd = di[r + 1] - di[r];
      const auto no = oi[r + 1] - oi[r];
      if ((nd == 1) && (no == 0) && isOne(dv[di[r]])) {
        T[static_cast<std::size_t>(r)] = first + dj[di[r]];
      } else if ((nd == 0) && (no == 1) && isOne(ov[oi[r]])) {
        T[static_cast<std::size_t>(r)] = cmap[oj[oi[r]]];
      } else {
        T.clear();
        break;
      }
    }
    P.HypreRead();
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
     * \param[in] m: global matrix, readable on the host, which must outlive
     * this object
     */
    explicit GlobalPositionsFinder(const hypre_ParCSRMatrix* const m)
        : A(m),
          diag(hypre_ParCSRMatrixDiag(m)),
          offd(hypre_ParCSRMatrixOffd(m)),
          first_row(hypre_ParCSRMatrixFirstRowIndex(m)),
          last_row(hypre_ParCSRMatrixLastRowIndex(m)),
          first_column(hypre_ParCSRMatrixFirstColDiag(m)),
          last_column(hypre_ParCSRMatrixLastColDiag(m)),
          diagonal_block(
              static_cast<std::size_t>(hypre_CSRMatrixNumCols(diag))),
          off_diagonal_block(
              static_cast<std::size_t>(hypre_CSRMatrixNumCols(offd))) {
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
        const auto nd = hypre_CSRMatrixNumNonzeros(this->diag);
        this->setRow(this->diagonal_block, this->diag, 0);
        this->setRow(this->off_diagonal_block, this->offd, nd);
      }
      if ((J >= this->first_column) && (J <= this->last_column)) {
        return this->getPosition(
            this->diagonal_block,
            static_cast<std::size_t>(J - this->first_column));
      }
      const auto* const cmap = hypre_ParCSRMatrixColMapOffd(this->A);
      const auto nc = hypre_CSRMatrixNumCols(this->offd);
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
     * \param[in] block: block
     * \param[in] offset: position of the first value of the block
     */
    void setRow(std::vector<Position>& positions,
                const hypre_CSRMatrix* const block,
                const HYPRE_Int offset) const {
      const auto* const i = hypre_CSRMatrixI(block);
      const auto* const j = hypre_CSRMatrixJ(block);
      for (auto p = i[this->row]; p != i[this->row + 1]; ++p) {
        positions[static_cast<std::size_t>(j[p])] = {
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
    //! \brief diagonal block
    const hypre_CSRMatrix* diag;
    //! \brief off diagonal block
    const hypre_CSRMatrix* offd;
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
    // pattern of the global matrix: product of hypre on a matrix of ones,
    // the sums of positive values being never null
    auto ones = mfem::SparseMatrix(A);
    ones = real{1};
    auto dA = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto Ph = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto G = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    dA.MakeSquareBlockDiag(comm, this->fes->GlobalVSize(),
                           this->fes->GetDofOffsets(), &ones);
    Ph.ConvertFrom(this->fes->Dof_TrueDof_Matrix());
    G.MakePtAP(dA, Ph);
    G.SetOperatorOwner(false);
    this->matrix.reset(G.As<mfem::HypreParMatrix>());
    this->matrix->HostRead();
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
    auto sent = 0;
    for (std::size_t q = 0; q != np; ++q) {
      send_offsets[q] = sent;
      send_counts[q] = static_cast<int>(sent_to[q].size());
      sent += send_counts[q];
    }
    this->sent_entries.SetSize(sent);
    auto* const se = this->sent_entries.HostWrite();
    auto sent_pairs = std::vector<HYPRE_BigInt>{};
    sent_pairs.reserve(2 * static_cast<std::size_t>(sent));
    for (std::size_t q = 0; q != np; ++q) {
      std::copy(sent_to[q].begin(), sent_to[q].end(), se + send_offsets[q]);
      sent_pairs.insert(sent_pairs.end(), pairs[q].begin(), pairs[q].end());
    }
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
    // pairs of global numbers (two values per entry)
    auto pair_counts = std::vector<int>(np);
    auto pair_offsets = std::vector<int>(np);
    auto received_pair_counts = std::vector<int>(np);
    auto received_pair_offsets = std::vector<int>(np);
    for (std::size_t q = 0; q != np; ++q) {
      pair_counts[q] = 2 * send_counts[q];
      pair_offsets[q] = 2 * send_offsets[q];
      received_pair_counts[q] = 2 * receive_counts[q];
      received_pair_offsets[q] = 2 * receive_offsets[q];
    }
    auto received_pairs =
        std::vector<HYPRE_BigInt>(2 * static_cast<std::size_t>(received));
    MPI_Alltoallv(sent_pairs.data(), pair_counts.data(), pair_offsets.data(),
                  HYPRE_MPI_BIG_INT, received_pairs.data(),
                  received_pair_counts.data(), received_pair_offsets.data(),
                  HYPRE_MPI_BIG_INT, comm);
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
      // contributions to each entry of the global matrix: the entries of
      // the local matrix, then the received values, in increasing order
      const auto ng =
          static_cast<std::size_t>(nd) + static_cast<std::size_t>(no);
      this->offsets.SetSize(static_cast<int>(ng) + 1);
      auto* const o = this->offsets.HostWrite();
      std::fill(o, o + ng + 1, 0);
      for (const auto* const ps : {&positions, &received_positions}) {
        for (const auto p : *ps) {
          if (p >= 0) {
            ++o[p + 1];
          }
        }
      }
      for (std::size_t p = 0; p != ng; ++p) {
        o[p + 1] += o[p];
      }
      this->contributions.SetSize(o[ng]);
      auto* const c = this->contributions.HostWrite();
      auto next = std::vector<int>(o, o + ng);
      auto contribution = 0;
      for (const auto* const ps : {&positions, &received_positions}) {
        for (const auto p : *ps) {
          if (p >= 0) {
            c[next[static_cast<std::size_t>(p)]++] = contribution;
          }
          ++contribution;
        }
      }
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
