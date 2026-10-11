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
#include "mfem/general/communication.hpp"
#include "mfem/linalg/handle.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/pfespace.hpp"

namespace mfem_mgis::batch {

  //! \return true if a value is exactly one
  static bool isOne(const real v) { return (v >= real{1}) && (v <= real{1}); }

  /*!
   * \return the global numbers of the true degrees of freedom of the local
   * degrees of freedom of a finite element space, read in the prolongation,
   * or an empty vector if a row of the prolongation has not a single entry
   * equal to one (on this process)
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
        return {};
      }
    }
    return T;
  }  // end of getGlobalTrueDegreesOfFreedom

  /*!
   * \return the position of the entry `(I, J)` (global numbers) in the
   * values of the part of a global matrix owned by the process (positions
   * in the diagonal block, then in the off diagonal block from the number
   * of entries of the diagonal block), or -1 if the entry is not in the
   * pattern of the matrix or if its row is not owned by the process
   * \param[in] A: global matrix
   * \param[in] I: row
   * \param[in] J: column
   */
  static int findGlobalPosition(const hypre_ParCSRMatrix* const A,
                                const HYPRE_BigInt I,
                                const HYPRE_BigInt J) {
    const auto first_row = hypre_ParCSRMatrixFirstRowIndex(A);
    const auto last_row = hypre_ParCSRMatrixLastRowIndex(A);
    if ((I < first_row) || (I > last_row)) {
      return -1;
    }
    const auto r = static_cast<HYPRE_Int>(I - first_row);
    const auto* const diag = hypre_ParCSRMatrixDiag(A);
    const auto first_column = hypre_ParCSRMatrixFirstColDiag(A);
    const auto last_column = hypre_ParCSRMatrixLastColDiag(A);
    if ((J >= first_column) && (J <= last_column)) {
      const auto c = static_cast<HYPRE_Int>(J - first_column);
      const auto* const i = hypre_CSRMatrixI(diag);
      const auto* const j = hypre_CSRMatrixJ(diag);
      for (auto p = i[r]; p != i[r + 1]; ++p) {
        if (j[p] == c) {
          return static_cast<int>(p);
        }
      }
      return -1;
    }
    const auto* const offd = hypre_ParCSRMatrixOffd(A);
    const auto* const cmap = hypre_ParCSRMatrixColMapOffd(A);
    const auto nc = hypre_CSRMatrixNumCols(offd);
    const auto* const c = std::lower_bound(cmap, cmap + nc, J);
    if ((c == cmap + nc) || (*c != J)) {
      return -1;
    }
    const auto cc = static_cast<HYPRE_Int>(c - cmap);
    const auto* const i = hypre_CSRMatrixI(offd);
    const auto* const j = hypre_CSRMatrixJ(offd);
    const auto nd = hypre_CSRMatrixNumNonzeros(diag);
    for (auto p = i[r]; p != i[r + 1]; ++p) {
      if (j[p] == cc) {
        return static_cast<int>(nd + p);
      }
    }
    return -1;
  }  // end of findGlobalPosition

  bool GlobalStiffnessMatrix::isSupported(
      const mfem::ParFiniteElementSpace& fes) noexcept {
    auto supported = int{0};
    try {
      // the construction of the prolongation is collective: it is built on
      // all the processes, including the ones without degree of freedom
      static_cast<void>(fes.Dof_TrueDof_Matrix());
      supported = ((!isHypreUsingGPU()) &&
                   ((fes.GetVSize() == 0) ||
                    (!getGlobalTrueDegreesOfFreedom(fes).empty())))
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
    const auto too_many_entries =
        static_cast<long long>(nd) + no > std::numeric_limits<int>::max();
    this->diagonal_block_entries = static_cast<int>(nd);
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
    const auto* const Ai = A.HostReadI();
    const auto* const Aj = A.HostReadJ();
    const auto n = A.NumNonZeroElems();
    this->positions.assign(static_cast<std::size_t>(n), -1);
    auto sent = std::vector<std::vector<int>>(static_cast<std::size_t>(nprocs));
    auto pairs = std::vector<std::vector<HYPRE_BigInt>>(
        static_cast<std::size_t>(nprocs));
    for (int r = 0; r != A.Height(); ++r) {
      const auto I = T[static_cast<std::size_t>(r)];
      const auto o = owner(I);
      for (auto k = Ai[r]; k != Ai[r + 1]; ++k) {
        const auto J = T[static_cast<std::size_t>(Aj[k])];
        if (o == rank) {
          const auto p = findGlobalPosition(g, I, J);
          found = found && (p >= 0);
          this->positions[static_cast<std::size_t>(k)] = p;
        } else {
          sent[static_cast<std::size_t>(o)].push_back(k);
          pairs[static_cast<std::size_t>(o)].push_back(I);
          pairs[static_cast<std::size_t>(o)].push_back(J);
        }
      }
    }
    // exchanges
    const auto np = static_cast<std::size_t>(nprocs);
    this->send_counts.assign(np, 0);
    this->send_offsets.assign(np, 0);
    this->sent_entries.clear();
    auto sent_pairs = std::vector<HYPRE_BigInt>{};
    for (std::size_t q = 0; q != np; ++q) {
      this->send_offsets[q] = static_cast<int>(this->sent_entries.size());
      this->send_counts[q] = static_cast<int>(sent[q].size());
      this->sent_entries.insert(this->sent_entries.end(), sent[q].begin(),
                                sent[q].end());
      sent_pairs.insert(sent_pairs.end(), pairs[q].begin(), pairs[q].end());
    }
    this->receive_counts.assign(np, 0);
    MPI_Alltoall(this->send_counts.data(), 1, MPI_INT,
                 this->receive_counts.data(), 1, MPI_INT, comm);
    this->receive_offsets.assign(np, 0);
    auto received = 0;
    for (std::size_t q = 0; q != np; ++q) {
      this->receive_offsets[q] = received;
      received += this->receive_counts[q];
    }
    // pairs of global numbers (two values per entry)
    auto pair_counts = std::vector<int>(np);
    auto pair_offsets = std::vector<int>(np);
    auto received_pair_counts = std::vector<int>(np);
    auto received_pair_offsets = std::vector<int>(np);
    for (std::size_t q = 0; q != np; ++q) {
      pair_counts[q] = 2 * this->send_counts[q];
      pair_offsets[q] = 2 * this->send_offsets[q];
      received_pair_counts[q] = 2 * this->receive_counts[q];
      received_pair_offsets[q] = 2 * this->receive_offsets[q];
    }
    auto received_pairs =
        std::vector<HYPRE_BigInt>(2 * static_cast<std::size_t>(received));
    MPI_Alltoallv(sent_pairs.data(), pair_counts.data(), pair_offsets.data(),
                  HYPRE_MPI_BIG_INT, received_pairs.data(),
                  received_pair_counts.data(), received_pair_offsets.data(),
                  HYPRE_MPI_BIG_INT, comm);
    this->received_positions.resize(static_cast<std::size_t>(received));
    for (std::size_t i = 0; i != static_cast<std::size_t>(received); ++i) {
      const auto p = findGlobalPosition(g, received_pairs[2 * i],
                                        received_pairs[2 * i + 1]);
      found = found && (p >= 0);
      this->received_positions[i] = p;
    }
    this->send_buffer.resize(this->sent_entries.size());
    this->receive_buffer.resize(static_cast<std::size_t>(received));
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
      this->matrix->HostReadWrite();
      hypre_ParCSRMatrix* const g = *(this->matrix);
      auto* const dv = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(g));
      auto* const ov = hypre_CSRMatrixData(hypre_ParCSRMatrixOffd(g));
      const auto nd = this->diagonal_block_entries;
      std::fill(dv, dv + nd, real{0});
      std::fill(ov, ov + hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixOffd(g)),
                real{0});
      const auto add = [dv, ov, nd](const int p, const real v) {
        if (p < nd) {
          dv[p] += v;
        } else {
          ov[p - nd] += v;
        }
      };
      const auto* const a = A.HostReadData();
      for (std::size_t k = 0; k != this->positions.size(); ++k) {
        if (this->positions[k] >= 0) {
          add(this->positions[k], a[k]);
        }
      }
      for (std::size_t i = 0; i != this->sent_entries.size(); ++i) {
        this->send_buffer[i] =
            a[static_cast<std::size_t>(this->sent_entries[i])];
      }
      const auto mpi_real = mfem::MPITypeMap<real>::mpi_type;
      MPI_Alltoallv(this->send_buffer.data(), this->send_counts.data(),
                    this->send_offsets.data(), mpi_real,
                    this->receive_buffer.data(), this->receive_counts.data(),
                    this->receive_offsets.data(), mpi_real, comm);
      for (std::size_t i = 0; i != this->received_positions.size(); ++i) {
        add(this->received_positions[i], this->receive_buffer[i]);
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
