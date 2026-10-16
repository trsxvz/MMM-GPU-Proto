/*!
 * \file   src/STRUMPACKSolvers.cxx
 * \brief  This file implements the module of the solvers and of the
 * preconditioner of STRUMPACK (see `LinearSolverModules.hxx`): its sparse
 * direct solver, in double or mixed precision, which compresses its frontal
 * matrices (block low-rank, HSS, HODLR, ZFP) on demand, and the
 * preconditioner given by an approximate factorization
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_STRUMPACK)

#include <memory>
#include <functional>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <utility>
#include <optional>
#include "mfem/linalg/strumpack.hpp"
#include "mfem/fem/pfespace.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief command line options of STRUMPACK, kept alive as long as the
   * solver which reads them
   */
  struct STRUMPACKCommandLine {
    /*!
     * \brief add an option
     * \param[in] n: name of the option, without the leading dashes
     * \param[in] v: value (empty for an option without value)
     */
    void add(const std::string& n, const std::string& v) {
      this->arguments.push_back("--" + n);
      if (!v.empty()) {
        this->arguments.push_back(v);
      }
    }
    //! \return the number of arguments
    [[nodiscard]] int argc() const {
      return static_cast<int>(this->arguments.size()) + 1;
    }
    //! \return the arguments, the first one being the name of the program
    [[nodiscard]] char** argv() {
      this->pointers.clear();
      this->pointers.push_back(this->program.data());
      for (auto& a : this->arguments) {
        this->pointers.push_back(a.data());
      }
      this->pointers.push_back(nullptr);
      return this->pointers.data();
    }

   private:
    //! \brief name of the program
    std::string program = "mfem-mgis";
    //! \brief arguments
    std::vector<std::string> arguments;
    //! \brief pointers to the arguments
    std::vector<char*> pointers;
  };  // end of struct STRUMPACKCommandLine

  /*!
   * \brief add an option given by a parameter to the command line options of
   * STRUMPACK: a boolean gives the option alone (a flag, such as
   * `sp_enable_replace_tiny_pivots`) when true and nothing when false, an
   * integer, a floating-point value (written with all its digits) or a
   * string gives the option and its value, an empty string giving the
   * option alone
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in, out] c: command line options
   * \param[in] n: name of the option
   * \param[in] p: parameter
   */
  [[nodiscard]] static bool addSTRUMPACKOption(Context& ctx,
                                               STRUMPACKCommandLine& c,
                                               const std::string& n,
                                               const Parameter& p) {
    if (n.empty() || n.starts_with('-')) {
      return ctx.registerErrorMessage(
          "invalid option '" + n +
          "' of STRUMPACK (the options are given without the leading "
          "dashes)");
    }
    if (mfem_mgis::is<bool>(p)) {
      const auto ov = mfem_mgis::get<bool>(ctx, p);
      if (isInvalid(ov)) {
        return false;
      }
      if (*ov) {
        c.add(n, "");
      }
      return true;
    }
    if (mfem_mgis::is<int>(p)) {
      const auto ov = mfem_mgis::get<int>(ctx, p);
      if (isInvalid(ov)) {
        return false;
      }
      c.add(n, std::to_string(*ov));
      return true;
    }
    if (mfem_mgis::is<double>(p)) {
      const auto ov = mfem_mgis::get<double>(ctx, p);
      if (isInvalid(ov)) {
        return false;
      }
      auto os = std::ostringstream{};
      os << std::setprecision(17) << *ov;
      c.add(n, os.str());
      return true;
    }
    if (mfem_mgis::is<std::string>(p)) {
      const auto ov = mfem_mgis::get<std::string>(ctx, p);
      if (isInvalid(ov)) {
        return false;
      }
      c.add(n, *ov);
      return true;
    }
    return ctx.registerErrorMessage(
        "invalid value of the option '" + n +
        "' of STRUMPACK (a boolean, an integer, a floating-point value or a "
        "string is expected)");
  }  // end of addSTRUMPACKOption

  /*!
   * \brief solver of MFEM giving access to the number of iterations of the
   * solver of STRUMPACK, which MFEM does not report
   */
  template <typename SolverType>
  struct STRUMPACKSolverWithIterations final : public SolverType {
    using SolverType::SolverType;
    /*!
     * \return the number of iterations of the Krylov solver (or of the
     * iterative refinement) of STRUMPACK during the last resolution, null
     * without Krylov solver
     */
    [[nodiscard]] int getNumberOfIterations() const {
      return this->solver_->Krylov_iterations();
    }
    //! \return the maximum number of iterations of STRUMPACK
    [[nodiscard]] int getMaximumNumberOfIterations() const {
      return this->solver_->options().maxit();
    }
  };  // end of struct STRUMPACKSolverWithIterations

  /*!
   * \brief solver of STRUMPACK: the matrices of hypre are copied to the
   * format of STRUMPACK, and the options of STRUMPACK are given by its
   * command line options, then by the settings of mfem-mgis (tolerances,
   * verbosity, reuse of the reordering). The solver of MFEM is built again,
   * with the command line options and the settings, when the global size of
   * the operator changes: the reuse of the reordering requires the same
   * sparsity pattern.
   *
   * STRUMPACK reports a success when its Krylov solver (used with the
   * compression of the frontal matrices) or its iterative refinement stops
   * at the maximum number of iterations: such a resolution is reported as
   * not converged.
   */
  template <typename SolverType>
  struct STRUMPACKDirectSolver final : public mfem::Solver,
                                       public LinearSolverConvergenceReport {
    //! \brief a setting of the solver
    using Setting = std::function<void(SolverType&)>;
    /*!
     * \param[in] comm: communicator
     * \param[in] c: command line options
     * \param[in] s: settings
     */
    STRUMPACKDirectSolver(MPI_Comm comm,
                          STRUMPACKCommandLine c,
                          std::vector<Setting> s)
        : communicator(comm),
          command_line(std::move(c)),
          settings(std::move(s)) {}
    void SetOperator(const mfem::Operator& op) override {
      const auto* const a = dynamic_cast<const mfem::HypreParMatrix*>(&op);
      if (a == nullptr) {
        mgis::raise(
            "STRUMPACKDirectSolver::SetOperator: the operator is not a matrix "
            "of hypre");
      }
      // the global size is the same on all processes, which take the same
      // decision (the creation of the solver of MFEM is collective)
      const auto n = a->GetGlobalNumRows();
      if ((this->solver == nullptr) || (n != this->global_size)) {
        this->create();
        this->global_size = n;
      }
      auto m = std::make_unique<mfem::STRUMPACKRowLocMatrix>(op);
      this->solver->SetOperator(*m);
      this->matrix = std::move(m);
      this->height = op.Height();
      this->width = op.Width();
    }
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
      if (this->solver == nullptr) {
        mgis::raise("STRUMPACKDirectSolver::Mult: the operator is not set");
      }
      this->solver->Mult(x, y);
      this->resolved = true;
    }
    [[nodiscard]] bool hasConverged() const noexcept override {
      return this->resolved && (this->solver->getNumberOfIterations() <
                                this->solver->getMaximumNumberOfIterations());
    }
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override {
      if (!this->resolved) {
        return {};
      }
      return this->solver->getNumberOfIterations();
    }

   private:
    //! \brief create the solver of MFEM
    void create() {
      this->solver.reset();
      this->resolved = false;
      this->solver =
          std::make_unique<STRUMPACKSolverWithIterations<SolverType>>(
              this->communicator, this->command_line.argc(),
              this->command_line.argv());
      this->solver->SetFromCommandLine();
      for (const auto& setting : this->settings) {
        setting(*(this->solver));
      }
    }
    //! \brief communicator
    const MPI_Comm communicator;
    //! \brief command line options, read by the solver
    STRUMPACKCommandLine command_line;
    //! \brief settings
    const std::vector<Setting> settings;
    /*!
     * \brief matrix in the format of STRUMPACK, which STRUMPACK copies: it
     * is kept since the solver of MFEM keeps a pointer to it
     */
    std::unique_ptr<mfem::STRUMPACKRowLocMatrix> matrix;
    //! \brief solver of MFEM
    std::unique_ptr<STRUMPACKSolverWithIterations<SolverType>> solver;
    //! \brief global size of the operator
    HYPRE_BigInt global_size = 0;
    //! \brief if true, the solver of MFEM performed a resolution
    mutable bool resolved = false;
  };  // end of struct STRUMPACKDirectSolver

  /*!
   * \return the solver of STRUMPACK
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] params: options (see `declareSTRUMPACKLinearSolvers`)
   * \param[in] preconditioner: if true, the solver is an approximate
   * factorization applied once (block low-rank compression and no Krylov
   * solver of STRUMPACK by default)
   */
  template <typename SolverType>
  [[nodiscard]] static std::unique_ptr<STRUMPACKDirectSolver<SolverType>>
  makeSTRUMPACKSolver(Context& ctx,
                      FiniteElementSpace<true>& fespace,
                      const Parameters& params,
                      const bool preconditioner) {
    using Problem = AbstractNonLinearEvolutionProblem;
    using Setting = typename STRUMPACKDirectSolver<SolverType>::Setting;
    auto options = std::vector<std::string>{"Options", "ReuseAnalysis",
                                            Problem::SolverVerbosityLevel};
    if (!preconditioner) {
      options.insert(options.end(), {Problem::SolverRelativeTolerance,
                                     Problem::SolverAbsoluteTolerance,
                                     Problem::SolverMaximumNumberOfIterations});
    }
    if (!checkParameters(ctx, params, options)) {
      return {};
    }
    auto c = STRUMPACKCommandLine{};
    if (preconditioner) {
      c.add("sp_compression", "blr");
      c.add("sp_Krylov_solver", "direct");
    }
    // the options of STRUMPACK given by the user come last, and override
    // the previous ones
    const auto oopts = get_if<Parameters>(ctx, params, "Options", Parameters{});
    if (isInvalid(oopts)) {
      return {};
    }
    for (const auto& [n, p] : *oopts) {
      if (!addSTRUMPACKOption(ctx, c, n, p)) {
        return {};
      }
    }
    // the settings of mfem-mgis, applied after the command line options
    auto settings = std::vector<Setting>{};
    const auto olvl =
        get_if<int>(ctx, params, Problem::SolverVerbosityLevel, 0);
    const auto oreuse = get_if<bool>(ctx, params, "ReuseAnalysis", false);
    if (isInvalid(olvl) || isInvalid(oreuse)) {
      return {};
    }
    settings.push_back([v = *olvl > 0, r = *oreuse](SolverType& s) {
      s.SetPrintFactorStatistics(v);
      s.SetPrintSolveStatistics(v);
      s.SetReorderingReuse(r);
    });
    if (contains(params, Problem::SolverRelativeTolerance)) {
      const auto ov =
          get<double>(ctx, params, Problem::SolverRelativeTolerance);
      if (isInvalid(ov)) {
        return {};
      }
      settings.push_back([v = static_cast<double>(boundRelativeTolerance(*ov))](
                             SolverType& s) { s.SetRelTol(v); });
    }
    if (contains(params, Problem::SolverAbsoluteTolerance)) {
      const auto ov =
          get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
      if (isInvalid(ov)) {
        return {};
      }
      settings.push_back([v = *ov](SolverType& s) { s.SetAbsTol(v); });
    }
    if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
      const auto ov =
          get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
      if (isInvalid(ov)) {
        return {};
      }
      settings.push_back([v = *ov](SolverType& s) { s.SetMaxIter(v); });
    }
    return std::make_unique<STRUMPACKDirectSolver<SolverType>>(
        fespace.GetComm(), std::move(c), std::move(settings));
  }  // end of makeSTRUMPACKSolver

  void declareSTRUMPACKLinearSolvers(LinearSolverFactory<true>& f) {
    using Problem = AbstractNonLinearEvolutionProblem;
    const auto description = [](std::string d) {
      return LinearSolverDescription{
          .library = "STRUMPACK",
          .kind = LinearSolverKind::DIRECT,
          .placement = LinearSolverPlacement::HOST,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = {"Options", "ReuseAnalysis", Problem::SolverVerbosityLevel,
                      Problem::SolverRelativeTolerance,
                      Problem::SolverAbsoluteTolerance,
                      Problem::SolverMaximumNumberOfIterations},
          .description = std::move(d)};
    };
    const auto add = [&f, &description]<typename SolverType>(
                         const std::string& n, std::string d) {
      declare(
          f, n, description(std::move(d)),
          [](Context& ctx, FiniteElementSpace<true>& fespace,
             const Parameters& params) -> LinearSolverHandler {
            auto s =
                makeSTRUMPACKSolver<SolverType>(ctx, fespace, params, false);
            if (s == nullptr) {
              return {};
            }
            return LinearSolverHandler{
                std::move(s), std::unique_ptr<LinearSolverPreconditioner>{}};
          });
    };
    add.template operator()<mfem::STRUMPACKSolver>(
        "STRUMPACKSolver",
        "sparse multifrontal solver of STRUMPACK, with optional compression "
        "of the frontal matrices (options of STRUMPACK)");
#if STRUMPACK_VERSION_MAJOR >= 7
    add.template operator()<mfem::STRUMPACKMixedPrecisionSolver>(
        "STRUMPACKMixedPrecisionSolver",
        "sparse multifrontal solver of STRUMPACK factorizing in single "
        "precision, with iterative refinement in double precision");
#endif /* STRUMPACK_VERSION_MAJOR >= 7 */
  }  // end of declareSTRUMPACKLinearSolvers

  void declareSTRUMPACKPreconditioners(PreconditionerFactory<true>& f) {
    using Problem = AbstractNonLinearEvolutionProblem;
    declare(f, "STRUMPACKPreconditioner",
            LinearSolverDescription{
                .library = "STRUMPACK",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Options", "ReuseAnalysis",
                            Problem::SolverVerbosityLevel},
                .description =
                    "approximate factorization of STRUMPACK (block low-rank "
                    "compression by default), applied once"},
            [](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& opts)
                -> std::unique_ptr<LinearSolverPreconditioner> {
              return makeSTRUMPACKSolver<mfem::STRUMPACKSolver>(ctx, fespace,
                                                                opts, true);
            });
  }  // end of declareSTRUMPACKPreconditioners

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_STRUMPACK) */
