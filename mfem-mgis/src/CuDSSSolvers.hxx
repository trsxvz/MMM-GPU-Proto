/*!
 * \file   src/CuDSSSolvers.hxx
 * \brief  This file declares the solvers of the cuDSS library of NVIDIA
 * (direct solver on the GPU) built by the factory of linear solvers
 * (internal header)
 * \date   12/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_CUDSSSOLVERS_HXX
#define LIB_MFEM_MGIS_SRC_CUDSSSOLVERS_HXX

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS)

#include <memory>
#include <optional>
#include <mpi.h>
#include "mfem/linalg/cudss.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "HypreMatrixPattern.hxx"

namespace mfem_mgis::internals {

  //! \brief settings of a cuDSS solver
  struct CuDSSSolverSettings {
    //! \brief type of the matrices
    mfem::CuDSSSolver::MatType matrix_type =
        mfem::CuDSSSolver::MatType::NONSYMMETRIC;
    //! \brief part of the matrices read by cuDSS
    mfem::CuDSSSolver::MatViewType view = mfem::CuDSSSolver::MatViewType::FULL;
  };

  /*!
   * \return true if the current GPU reads the pageable memory of the host
   * (heterogeneous memory management, or address translation services).
   * When hypre computes on the host, `mfem::CuDSSSolver` (MFEM 4.10) gives
   * cuDSS host arrays through a copy between device arrays, which is only
   * valid on such a GPU.
   */
  [[nodiscard]] bool canGPUReadPageableHostMemory();

  /*!
   * \brief cuDSS solver of the factory of linear solvers.
   *
   * `mfem::CuDSSSolver` (MFEM 4.10) stops the program when given a matrix
   * whose local size differs from the one of the previous matrix, keeps the
   * descriptors of the right-hand side and of the solution of its first
   * matrix, and the reuse of its analysis (ordering and symbolic
   * factorization) is chosen before its first matrix, after which only the
   * number of coefficients of the next matrices is checked (MFEM stops the
   * program if it changes). A new `mfem::CuDSSSolver` is thus created on
   * all processes when the layout of the matrix (local numbers of rows and
   * columns, global number of rows, first row) changed on one process and,
   * if the analysis is reused, when its pattern changed on one process.
   * With the reuse, the analysis of the previous matrix is reused when the
   * new one has the same pattern on all processes, as the matrices of the
   * successive iterations of the non linear solvers.
   *
   * cuDSS does not report the failures of the factorization (null pivots,
   * matrix declared positive definite which is not), which give wrong
   * solutions, finite or not: each solution is checked by its relative
   * residual, `|b - A x| <= 1e-3 |b|`, computed where the vectors are. The
   * result is given by `hasSucceeded`, read by `mfem_mgis::hasConverged`.
   * The bound only rejects wrong solutions: the relative residual of a
   * stable factorization is of the order of the precision times the
   * condition number of the matrix. As for the iterative solvers of MFEM,
   * the operator must exist as long as the solver is used.
   */
  struct CuDSSDirectSolver final : public mfem::Solver,
                                   public LinearSolverConvergenceReport {
    /*!
     * \brief constructor
     * \param[in] c: communicator
     * \param[in] s: settings
     * \param[in] r: reuse the analysis for the matrices of the same pattern
     */
    CuDSSDirectSolver(MPI_Comm c, const CuDSSSolverSettings& s, const bool r);
    //! \brief set the operator, a parallel matrix of hypre
    void SetOperator(const mfem::Operator&) override;
    //! \brief solve
    void Mult(const mfem::Vector&, mfem::Vector&) const override;
    //! \brief solve for several right-hand sides
    void ArrayMult(const mfem::Array<const mfem::Vector*>&,
                   mfem::Array<mfem::Vector*>&) const override;
    //! \return true if the residuals of the last solutions are small
    [[nodiscard]] bool hasSucceeded() const noexcept;
    //! \return `hasSucceeded()`
    [[nodiscard]] bool hasConverged() const noexcept override;
    //! \return an empty value: a direct solver makes no iteration
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override;
    //! \return true if the analysis is reused
    [[nodiscard]] bool reusesAnalysis() const noexcept;
    //! \return the number of analyses made
    [[nodiscard]] size_type getNumberOfAnalyses() const noexcept;

   private:
    //! \brief layout of a matrix
    struct Layout {
      //! \brief comparison operator
      bool operator==(const Layout&) const = default;
      //! \brief local number of rows
      int number_of_rows = 0;
      //! \brief local number of columns
      int number_of_columns = 0;
      //! \brief global number of rows
      HYPRE_BigInt global_number_of_rows = 0;
      //! \brief global index of the first local row
      HYPRE_BigInt first_row = 0;
    };
    //! \return the layout of a matrix
    static Layout getLayout(const mfem::HypreParMatrix&);
    //! \brief create a new solver of MFEM, the previous one being released
    void makeSolver();
    /*!
     * \brief check the solutions by their relative residuals
     * \param[in] X: right-hand sides
     * \param[in] Y: solutions
     */
    void checkSolutions(const mfem::Array<const mfem::Vector*>&,
                        const mfem::Array<mfem::Vector*>&) const;
    //! \brief communicator
    const MPI_Comm communicator;
    //! \brief settings
    const CuDSSSolverSettings settings;
    //! \brief reuse of the analysis
    const bool reuse_analysis;
    //! \brief solver of MFEM
    std::unique_ptr<mfem::CuDSSSolver> solver;
    //! \brief operator
    const mfem::HypreParMatrix* matrix = nullptr;
    //! \brief residual of a solution
    mutable mfem::Vector residual;
    //! \brief layout of the matrix of the solver, if any
    std::optional<Layout> layout;
    //! \brief pattern of the last analysed matrix, if the analysis is reused
    HypreMatrixPatternRecord pattern;
    //! \brief number of analyses made
    size_type number_of_analyses = 0;
    //! \brief if true, the residuals of the last solutions are small
    mutable bool succeeded = true;
  };  // end of struct CuDSSDirectSolver

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS) */

#endif /* LIB_MFEM_MGIS_SRC_CUDSSSOLVERS_HXX */
