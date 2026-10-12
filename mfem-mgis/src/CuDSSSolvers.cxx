/*!
 * \file   src/CuDSSSolvers.cxx
 * \brief  This file implements the solvers of the cuDSS library of NVIDIA
 * \date   12/10/2026
 */

#include "CuDSSSolvers.hxx"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS)

#include <utility>
#include <cuda_runtime_api.h>
#include "mfem/general/error.hpp"

namespace mfem_mgis::internals {

  bool canGPUReadPageableHostMemory() {
    auto device = int{};
    auto pageable_memory_access = int{};
    if (cudaGetDevice(&device) != cudaSuccess) {
      return false;
    }
    if (cudaDeviceGetAttribute(&pageable_memory_access,
                               cudaDevAttrPageableMemoryAccess,
                               device) != cudaSuccess) {
      return false;
    }
    return pageable_memory_access != 0;
  }  // end of canGPUReadPageableHostMemory

  CuDSSDirectSolver::CuDSSDirectSolver(const MPI_Comm c,
                                       const CuDSSSolverSettings& s,
                                       const bool r)
      : communicator(c), settings(s), reuse_analysis(r) {}

  CuDSSDirectSolver::Layout CuDSSDirectSolver::getLayout(
      const mfem::HypreParMatrix& A) {
    const auto* const m = static_cast<hypre_ParCSRMatrix*>(A);
    return {.number_of_rows = A.Height(),
            .number_of_columns = A.Width(),
            .global_number_of_rows = hypre_ParCSRMatrixGlobalNumRows(m),
            .first_row = hypre_ParCSRMatrixFirstRowIndex(m)};
  }  // end of getLayout

  void CuDSSDirectSolver::makeSolver() {
    // the previous solver and its device memory are released first
    this->solver.reset();
    this->solver = std::make_unique<mfem::CuDSSSolver>(this->communicator);
    // the view type depends on the matrix type
    this->solver->SetMatrixSymType(this->settings.matrix_type);
    this->solver->SetMatrixViewType(this->settings.view);
    this->solver->SetReorderingReuse(this->reuse_analysis);
  }  // end of makeSolver

  void CuDSSDirectSolver::SetOperator(const mfem::Operator& op) {
    const auto* const A = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    MFEM_VERIFY(A != nullptr,
                "CuDSSSolver: the operator must be a HypreParMatrix");
    const auto l = getLayout(*A);
    // collective decision: a new solver of MFEM is created on all processes
    // or on none
    const auto same_layout = this->layout == l;
    const auto same =
        this->reuse_analysis
            ? this->pattern.isSameOnAllProcesses(*A, this->communicator,
                                                 same_layout)
            : isTrueOnAllProcesses(same_layout, this->communicator);
    // the layout and the pattern are forgotten until the factorization
    // succeeds: if it fails by an exception, the next matrix is given to a
    // new solver
    auto previous = std::exchange(this->pattern, HypreMatrixPatternRecord{});
    this->layout.reset();
    if (!same) {
      this->makeSolver();
    }
    this->height = op.Height();
    this->width = op.Width();
    this->matrix = nullptr;
    this->solver->SetOperator(op);
    this->matrix = A;
    this->layout = l;
    if (!this->reuse_analysis) {
      // MFEM analyses each matrix
      ++(this->number_of_analyses);
      return;
    }
    if (same) {
      this->pattern = std::move(previous);
      return;
    }
    ++(this->number_of_analyses);
    this->pattern.record(*A);
  }  // end of SetOperator

  void CuDSSDirectSolver::checkSolutions(
      const mfem::Array<const mfem::Vector*>& X,
      const mfem::Array<mfem::Vector*>& Y) const {
    // relative bound of the residuals
    constexpr auto tolerance = real{1e-3};
    auto succeeded_on_all_processes = true;
    for (int i = 0; i != Y.Size(); ++i) {
      const auto& x = *(X[i]);
      const auto& y = *(Y[i]);
      // computed where the vectors are
      this->residual.SetSize(x.Size());
      this->residual.UseDevice(true);
      this->matrix->Mult(y, this->residual);
      this->residual -= x;
      const auto r2 = mfem::InnerProduct(this->communicator, this->residual,
                                         this->residual);
      const auto b2 = mfem::InnerProduct(this->communicator, x, x);
      // false for values which are not finite
      succeeded_on_all_processes =
          succeeded_on_all_processes && (r2 <= tolerance * tolerance * b2);
    }
    this->succeeded = succeeded_on_all_processes;
  }  // end of checkSolutions

  void CuDSSDirectSolver::Mult(const mfem::Vector& x, mfem::Vector& y) const {
    MFEM_VERIFY(this->solver != nullptr, "no operator given");
    this->solver->Mult(x, y);
    auto X = mfem::Array<const mfem::Vector*>(1);
    auto Y = mfem::Array<mfem::Vector*>(1);
    X[0] = &x;
    Y[0] = &y;
    this->checkSolutions(X, Y);
  }  // end of Mult

  void CuDSSDirectSolver::ArrayMult(const mfem::Array<const mfem::Vector*>& X,
                                    mfem::Array<mfem::Vector*>& Y) const {
    MFEM_VERIFY(this->solver != nullptr, "no operator given");
    this->solver->ArrayMult(X, Y);
    this->checkSolutions(X, Y);
  }  // end of ArrayMult

  bool CuDSSDirectSolver::hasSucceeded() const noexcept {
    return this->succeeded;
  }  // end of hasSucceeded

  bool CuDSSDirectSolver::reusesAnalysis() const noexcept {
    return this->reuse_analysis;
  }  // end of reusesAnalysis

  size_type CuDSSDirectSolver::getNumberOfAnalyses() const noexcept {
    return this->number_of_analyses;
  }  // end of getNumberOfAnalyses

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS) */
