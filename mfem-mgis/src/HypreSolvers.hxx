/*!
 * \file   src/HypreSolvers.hxx
 * \brief  This file declares the solvers and the preconditioners of hypre
 * which are not wrapped by MFEM: its BiCGSTAB, COGMRES and LGMRES Krylov
 * solvers, its hybrid solver and its factorized sparse approximate inverse
 * (internal header, see `HypreLinearSolvers.cxx`)
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_HYPRESOLVERS_HXX
#define LIB_MFEM_MGIS_SRC_HYPRESOLVERS_HXX

#include "mfem/config/config.hpp"

#ifdef MFEM_USE_MPI

#include <vector>
#include <optional>
#include <functional>
#include "mfem/linalg/hypre.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"

namespace mfem_mgis::internals {

  //! \brief functions of the interface of a Krylov solver of hypre
  struct HypreKrylovSolverFunctions {
    //! \brief creation of the solver
    HYPRE_Int (*create)(MPI_Comm, HYPRE_Solver*);
    //! \brief destruction of the solver
    HYPRE_Int (*destroy)(HYPRE_Solver);
    //! \brief setup of the solver
    HYPRE_PtrToParSolverFcn setup;
    //! \brief resolution
    HYPRE_PtrToParSolverFcn solve;
    //! \brief relative tolerance
    HYPRE_Int (*set_tolerance)(HYPRE_Solver, HYPRE_Real);
    //! \brief absolute tolerance
    HYPRE_Int (*set_absolute_tolerance)(HYPRE_Solver, HYPRE_Real);
    //! \brief maximum number of iterations
    HYPRE_Int (*set_maximum_number_of_iterations)(HYPRE_Solver, HYPRE_Int);
    //! \brief print level
    HYPRE_Int (*set_print_level)(HYPRE_Solver, HYPRE_Int);
    //! \brief preconditioner
    HYPRE_Int (*set_preconditioner)(HYPRE_Solver,
                                    HYPRE_PtrToParSolverFcn,
                                    HYPRE_PtrToParSolverFcn,
                                    HYPRE_Solver);
    //! \brief number of iterations of the last resolution
    HYPRE_Int (*get_number_of_iterations)(HYPRE_Solver, HYPRE_Int*);
    //! \brief final relative residual norm of the last resolution
    HYPRE_Int (*get_final_relative_residual_norm)(HYPRE_Solver, HYPRE_Real*);
    //! \brief convergence of the last resolution
    HYPRE_Int (*get_converged)(HYPRE_Solver, HYPRE_Int*);
    //! \brief dimension of the Krylov space, if any
    HYPRE_Int (*set_krylov_dimension)(HYPRE_Solver, HYPRE_Int) = nullptr;
  };

  //! \return the functions of the BiCGSTAB solver of hypre
  [[nodiscard]] HypreKrylovSolverFunctions getHypreBiCGSTABFunctions();
  //! \return the functions of the COGMRES solver of hypre
  [[nodiscard]] HypreKrylovSolverFunctions getHypreCOGMRESFunctions();
  //! \return the functions of the LGMRES solver of hypre
  [[nodiscard]] HypreKrylovSolverFunctions getHypreLGMRESFunctions();

  /*!
   * \brief Krylov solver of hypre which is not wrapped by MFEM, described by
   * the functions of its interface (see `HypreKrylovSolverFunctions`), and
   * used as the Krylov solvers of hypre wrapped by MFEM (`mfem::HypreGMRES`
   * for instance): its preconditioner is a solver of hypre, which is given
   * the operators. A resolution which does not converge is not an error
   * (the errors of hypre are ignored): its convergence is given by the flag
   * of hypre and by the code returned by hypre (see `hasConverged`).
   */
  struct HypreKrylovSolver final : public mfem::HypreSolver,
                                   public LinearSolverConvergenceReport {
    /*!
     * \param[in] comm: communicator
     * \param[in] f: functions of the solver
     */
    HypreKrylovSolver(MPI_Comm comm, const HypreKrylovSolverFunctions& f);
    HypreKrylovSolver(HypreKrylovSolver&&) = delete;
    HypreKrylovSolver(const HypreKrylovSolver&) = delete;
    HypreKrylovSolver& operator=(HypreKrylovSolver&&) = delete;
    HypreKrylovSolver& operator=(const HypreKrylovSolver&) = delete;
    //! \return the solver of hypre
    operator HYPRE_Solver() const override;
    //! \return the setup function of hypre
    HYPRE_PtrToParSolverFcn SetupFcn() const override;
    //! \return the solve function of hypre
    HYPRE_PtrToParSolverFcn SolveFcn() const override;
    /*!
     * \brief set the operator, which must be a matrix of hypre, and give it
     * to the preconditioner, as `mfem::HypreGMRES::SetOperator`
     * \param[in] op: operator
     */
    void SetOperator(const mfem::Operator& op) override;
    /*!
     * \brief solve the system, as `mfem::HypreSolver::Mult`, the code
     * returned by hypre being kept (see `hasConverged`)
     * \param[in] b: right hand side
     * \param[in, out] x: solution
     */
    void Mult(const mfem::HypreParVector& b,
              mfem::HypreParVector& x) const override;
    using mfem::HypreSolver::Mult;
    /*!
     * \brief set the preconditioner
     * \param[in] p: preconditioner (not owned)
     */
    void SetPreconditioner(mfem::HypreSolver& p);
    /*!
     * \brief set the relative tolerance
     * \param[in] tol: tolerance
     */
    void SetTol(const real tol);
    /*!
     * \brief set the absolute tolerance
     * \param[in] tol: tolerance
     */
    void SetAbsTol(const real tol);
    /*!
     * \brief set the maximum number of iterations
     * \param[in] n: maximum number of iterations
     */
    void SetMaxIter(const int n);
    /*!
     * \brief set the print level of hypre
     * \param[in] l: print level
     */
    void SetPrintLevel(const int l);
    /*!
     * \brief set the dimension of the Krylov space
     * \return false if the solver has no Krylov space dimension
     * \param[in] n: dimension
     */
    [[nodiscard]] bool SetKDim(const int n);
    [[nodiscard]] bool hasConverged() const noexcept override;
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override;
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override;
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override;
    //! \brief destructor
    ~HypreKrylovSolver() override;

   private:
    //! \brief functions of the solver
    const HypreKrylovSolverFunctions functions;
    //! \brief solver of hypre
    HYPRE_Solver solver = nullptr;
    //! \brief preconditioner, if any (not owned)
    mfem::HypreSolver* preconditioner = nullptr;
    //! \brief code returned by hypre at the last resolution
    mutable HYPRE_Int error_code = 0;
    //! \brief relative tolerance
    real relative_tolerance = 1e-6;
    //! \brief maximum number of iterations
    int maximum_number_of_iterations = 100;
  };  // end of struct HypreKrylovSolver

  //! \brief a setting of a solver of hypre, applied to its handle
  using HypreSetting = std::function<void(HYPRE_Solver)>;

  /*!
   * \brief hybrid solver of hypre: a Krylov solver (conjugate gradient by
   * default) preconditioned by the diagonal scaling, switching to the
   * BoomerAMG preconditioner when the convergence is too slow. The solver of
   * hypre is created again, with the settings, for each new operator.
   */
  struct HypreHybridSolver final : public mfem::HypreSolver,
                                   public LinearSolverConvergenceReport {
    /*!
     * \param[in] s: settings
     * \param[in] tol: relative tolerance
     * \param[in] n: maximum number of iterations of each phase
     */
    HypreHybridSolver(std::vector<HypreSetting> s, const real tol, const int n);
    HypreHybridSolver(HypreHybridSolver&&) = delete;
    HypreHybridSolver(const HypreHybridSolver&) = delete;
    HypreHybridSolver& operator=(HypreHybridSolver&&) = delete;
    HypreHybridSolver& operator=(const HypreHybridSolver&) = delete;
    //! \return the solver of hypre
    operator HYPRE_Solver() const override;
    //! \return the setup function of hypre
    HYPRE_PtrToParSolverFcn SetupFcn() const override;
    //! \return the solve function of hypre
    HYPRE_PtrToParSolverFcn SolveFcn() const override;
    /*!
     * \brief set the operator, which must be a matrix of hypre: the hybrid
     * solver of hypre is built again
     * \param[in] op: operator
     */
    void SetOperator(const mfem::Operator& op) override;
    [[nodiscard]] bool hasConverged() const noexcept override;
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override;
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override;
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override;
    //! \brief destructor
    ~HypreHybridSolver() override;

   private:
    //! \brief create the solver of hypre and apply the settings
    void create();
    //! \brief solver of hypre
    HYPRE_Solver solver = nullptr;
    //! \brief settings
    const std::vector<HypreSetting> settings;
    //! \brief relative tolerance
    real relative_tolerance;
    //! \brief maximum number of iterations of each phase
    const int maximum_number_of_iterations;
  };  // end of struct HypreHybridSolver

  /*!
   * \brief factorized sparse approximate inverse of hypre (FSAI), applied
   * as a preconditioner (one iteration from a null initial guess). The
   * solver of hypre is created again, with the settings, for each new
   * operator.
   */
  struct HypreFSAI final : public mfem::HypreSolver {
    //! \param[in] s: settings
    explicit HypreFSAI(std::vector<HypreSetting> s);
    HypreFSAI(HypreFSAI&&) = delete;
    HypreFSAI(const HypreFSAI&) = delete;
    HypreFSAI& operator=(HypreFSAI&&) = delete;
    HypreFSAI& operator=(const HypreFSAI&) = delete;
    //! \return the solver of hypre
    operator HYPRE_Solver() const override;
    //! \return the setup function of hypre
    HYPRE_PtrToParSolverFcn SetupFcn() const override;
    //! \return the solve function of hypre
    HYPRE_PtrToParSolverFcn SolveFcn() const override;
    /*!
     * \brief set the operator, which must be a matrix of hypre: the
     * factorized sparse approximate inverse of hypre is built again
     * \param[in] op: operator
     */
    void SetOperator(const mfem::Operator& op) override;
    //! \brief destructor
    ~HypreFSAI() override;

   private:
    //! \brief create the solver of hypre and apply the settings
    void create();
    //! \brief solver of hypre
    HYPRE_Solver solver = nullptr;
    //! \brief settings
    const std::vector<HypreSetting> settings;
  };  // end of struct HypreFSAI

  /*!
   * \brief multigrid reduction of hypre (MGR), applied as a preconditioner
   * (one cycle from a null initial guess). The unknowns are labelled by
   * their component (point markers, whatever the ordering of the unknowns);
   * at each level of reduction, the components which are kept (coarse
   * points) are given, the other ones (fine points) being relaxed and
   * eliminated; the coarsest level is solved by the BoomerAMG of hypre
   * (default of MGR). The solver of hypre is created again, with the
   * settings, for each new operator.
   */
  struct HypreMGR final : public mfem::HypreSolver {
    /*!
     * \param[in] s: settings
     * \param[in] c: components kept at each level of reduction
     * \param[in] vdim: number of components of the unknowns
     * \param[in] by_nodes: if true, the unknowns are ordered by nodes
     * \param[in] r: F-relaxation of each level, if not empty (it then
     * replaces the F-relaxation given by the settings)
     */
    HypreMGR(std::vector<HypreSetting> s,
             std::vector<std::vector<HYPRE_Int>> c,
             const int vdim,
             const bool by_nodes,
             std::vector<HYPRE_Int> r = {});
    HypreMGR(HypreMGR&&) = delete;
    HypreMGR(const HypreMGR&) = delete;
    HypreMGR& operator=(HypreMGR&&) = delete;
    HypreMGR& operator=(const HypreMGR&) = delete;
    //! \return the solver of hypre
    operator HYPRE_Solver() const override;
    //! \return the setup function of hypre
    HYPRE_PtrToParSolverFcn SetupFcn() const override;
    //! \return the solve function of hypre
    HYPRE_PtrToParSolverFcn SolveFcn() const override;
    /*!
     * \brief set the operator, which must be a matrix of hypre whose local
     * size is a multiple of the number of components: the multigrid
     * reduction of hypre is built again
     * \param[in] op: operator
     */
    void SetOperator(const mfem::Operator& op) override;
    //! \brief destructor
    ~HypreMGR() override;

   private:
    //! \brief create the solver of hypre and apply the settings
    void create();
    //! \brief solver of hypre
    HYPRE_Solver solver = nullptr;
    //! \brief settings
    const std::vector<HypreSetting> settings;
    //! \brief components kept at each level of reduction
    std::vector<std::vector<HYPRE_Int>> coarse_components;
    //! \brief number of components kept at each level of reduction
    std::vector<HYPRE_Int> number_of_coarse_components;
    //! \brief pointers to the components kept at each level
    std::vector<HYPRE_Int*> coarse_components_pointers;
    //! \brief F-relaxation of each level, if not empty
    std::vector<HYPRE_Int> level_relaxation_types;
    /*!
     * \brief component of each local unknown, read by hypre at the setup
     * (hypre keeps the pointer)
     */
    std::vector<HYPRE_Int> point_markers;
    //! \brief number of components of the unknowns
    const int number_of_components;
    //! \brief if true, the unknowns are ordered by nodes
    const bool ordered_by_nodes;
  };  // end of struct HypreMGR

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */

#endif /* LIB_MFEM_MGIS_SRC_HYPRESOLVERS_HXX */
