/*!
 * \file   src/PETScLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of PETSc (see `LinearSolverModules.hxx`): a Krylov
 * solver (`KSP`) and a preconditioner (`PC`) of PETSc, configured by the
 * options database of PETSc
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)

#include <atomic>
#include <memory>
#include <tuple>
#include <string>
#include <vector>
#include <iomanip>
#include <optional>
#include <sstream>
#include <utility>
#include <string_view>
#include "petscsys.h"
#include "petscksp.h"
#include "mfem/linalg/petsc.hpp"
#include "mfem/fem/pfespace.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "PETScSupport.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  /*!
   * \brief initialize PETSc, unless it is already initialized (by
   * mfem-mgis, or by the program, which then finalizes it)
   */
  static void initializePETSc() {
    auto initialized = PetscBool{};
    if ((PetscInitialized(&initialized) == 0) && (initialized == PETSC_TRUE)) {
      return;
    }
    initializePETScOnDemand();
  }  // end of initializePETSc

  /*!
   * \return a prefix of the options database of PETSc which is unique in
   * the process, so that the options of two solvers do not interfere
   * \param[in] kind: kind of object (`ksp`, `pc`)
   */
  [[nodiscard]] static std::string makePETScOptionsPrefix(
      const std::string_view kind) {
    static auto count = std::atomic<std::size_t>{};
    return "mfem_mgis_" + std::string{kind} + std::to_string(count++) + "_";
  }  // end of makePETScOptionsPrefix

  /*!
   * \return the value of an option of PETSc given by a parameter, or an
   * empty value on error: a boolean (`true` or `false`), an integer, a
   * floating-point value (written with all its digits) or a string (an empty
   * string being an option without value)
   * \param[in, out] ctx: execution context
   * \param[in] n: name of the option
   * \param[in] p: parameter
   */
  [[nodiscard]] static std::optional<std::string> getPETScOptionValue(
      Context& ctx, const std::string& n, const Parameter& p) {
    if (mfem_mgis::is<bool>(p)) {
      const auto ov = mfem_mgis::get<bool>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return std::string{*ov ? "true" : "false"};
    }
    if (mfem_mgis::is<int>(p)) {
      const auto ov = mfem_mgis::get<int>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return std::to_string(*ov);
    }
    if (mfem_mgis::is<double>(p)) {
      const auto ov = mfem_mgis::get<double>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      auto os = std::ostringstream{};
      os << std::setprecision(17) << *ov;
      return os.str();
    }
    if (mfem_mgis::is<std::string>(p)) {
      const auto ov = mfem_mgis::get<std::string>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return std::string{*ov};
    }
    std::ignore = ctx.registerErrorMessage(
        "invalid value of the option '" + n +
        "' of PETSc (a boolean, an integer, a floating-point value or a "
        "string is expected)");
    return {};
  }  // end of getPETScOptionValue

  /*!
   * \brief options of a solver of PETSc, written in the options database of
   * PETSc with a prefix which is unique to the solver, and removed from it
   * when the solver is destroyed. This class is the first base class of the
   * solvers of PETSc, so that the prefix is known when the object of PETSc
   * is created and that the options are removed after its destruction.
   *
   * PETSc ignores the options it does not know (for instance a misspelt
   * option): the options which PETSc never read are reported when a solver
   * which has been used is destroyed.
   */
  struct PETScOptions {
    /*!
     * \param[in] c: communicator
     * \param[in] kind: kind of object (`ksp`, `pc`)
     */
    PETScOptions(MPI_Comm c, const std::string_view kind)
        : prefix(makePETScOptionsPrefix(kind)), communicator(c) {}
    PETScOptions(PETScOptions&&) = delete;
    PETScOptions(const PETScOptions&) = delete;
    PETScOptions& operator=(PETScOptions&&) = delete;
    PETScOptions& operator=(const PETScOptions&) = delete;
    /*!
     * \brief write options in the options database of PETSc
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] options: options, given without the leading dash and the
     * prefix (for instance `{"ksp_type", "cg"}`)
     */
    [[nodiscard]] bool setOptions(Context& ctx, const Parameters& options) {
      for (const auto& [n, p] : options) {
        if (n.empty() || n.starts_with('-')) {
          return ctx.registerErrorMessage(
              "invalid option '" + n +
              "' of PETSc (the options are given without the leading dash)");
        }
        const auto ovalue = getPETScOptionValue(ctx, n, p);
        if (isInvalid(ovalue)) {
          return false;
        }
        auto option = "-" + this->prefix + n;
        const auto* const value = ovalue->empty() ? nullptr : ovalue->c_str();
        if (PetscOptionsSetValue(nullptr, option.c_str(), value) != 0) {
          return ctx.registerErrorMessage("the option '" + n +
                                          "' can't be given to PETSc");
        }
        this->options.push_back(std::move(option));
      }
      return true;
    }
    /*!
     * \brief destructor, reporting the options which PETSc never read, if
     * the solver has been used, and removing the options from the options
     * database (unless PETSc has been finalized)
     */
    ~PETScOptions() {
      auto finalized = PetscBool{};
      if ((PetscFinalized(&finalized) != 0) || (finalized == PETSC_TRUE)) {
        return;
      }
      auto rank = int{};
      MPI_Comm_rank(this->communicator, &rank);
      for (const auto& o : this->options) {
        auto read = PetscBool{};
        if ((this->used) && (rank == 0) &&
            (PetscOptionsUsed(nullptr, o.c_str() + 1, &read) == 0) &&
            (read == PETSC_FALSE)) {
          getErrorStream() << "warning: the option '"
                           << o.substr(this->prefix.size() + 1)
                           << "' was never read by PETSc\n";
        }
        std::ignore = PetscOptionsClearValue(nullptr, o.c_str());
      }
    }
    //! \brief prefix of the options
    const std::string prefix;

   protected:
    //! \brief mark the solver as used (called by the `Mult` methods)
    void markAsUsed() const noexcept { this->used = true; }

   private:
    //! \brief communicator of the solver
    const MPI_Comm communicator;
    //! \brief options written in the options database, with their prefix
    std::vector<std::string> options;
    //! \brief if true, the solver has been used
    mutable bool used = false;
  };  // end of struct PETScOptions

  /*!
   * \brief Krylov solver of PETSc (`KSP`) configured by the options database
   * of PETSc. The operators are converted to the AIJ format of PETSc at each
   * call to `SetOperator`, so that the preconditioners of PETSc can be built
   * (they are thus built again at each call). A preconditioner of mfem-mgis
   * is seen by PETSc as a shell preconditioner, a preconditioner of PETSc
   * (`PETScPreconditioner`) becomes the preconditioner of the Krylov solver.
   *
   * The `iterative_mode` member of MFEM is given to PETSc at the
   * construction only (the solvers of mfem-mgis start from a null initial
   * guess). A resolution stopped by the maximum number of iterations
   * without norm (option `ksp_norm_type none`) is reported as converged, as
   * PETSc does.
   */
  struct PETScKrylovSolver final : public PETScOptions,
                                   public mfem::PetscLinearSolver,
                                   public LinearSolverConvergenceReport {
    /*!
     * \param[in] comm: communicator
     */
    explicit PETScKrylovSolver(MPI_Comm comm)
        : PETScOptions(comm, "ksp"),
          mfem::PetscLinearSolver(comm, this->prefix, false) {}
    /*!
     * \brief set a preconditioner of mfem-mgis: a preconditioner of PETSc
     * (`PETScPreconditioner`) becomes the preconditioner of the Krylov
     * solver, the others are seen by PETSc as shell preconditioners. This
     * method hides the one of MFEM, which is not virtual.
     * \param[in] p: preconditioner
     */
    void SetPreconditioner(mfem::Solver& p) {
      if (dynamic_cast<mfem::PetscPreconditioner*>(&p) == nullptr) {
        this->shell_preconditioner = &p;
      }
      mfem::PetscLinearSolver::SetPreconditioner(p);
    }
    void SetOperator(const mfem::Operator& op) override {
      // unlike the Krylov solvers of MFEM, the solver of PETSc does not give
      // the operator to a shell preconditioner (the operator of a
      // preconditioner of PETSc is set by PETSc)
      if (this->shell_preconditioner != nullptr) {
        this->shell_preconditioner->SetOperator(op);
      }
      mfem::PetscLinearSolver::SetOperator(op);
    }
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      this->markAsUsed();
      mfem::PetscLinearSolver::Mult(b, x);
    }
    [[nodiscard]] bool hasConverged() const noexcept override {
      // the accessors of MFEM are not constant
      return const_cast<PETScKrylovSolver*>(this)->GetConverged() > 0;
    }
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override {
      return const_cast<PETScKrylovSolver*>(this)->GetNumIterations();
    }
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override {
      auto rtol = PetscReal{};
      if (KSPGetTolerances(KSP(*this), &rtol, nullptr, nullptr, nullptr) != 0) {
        return {};
      }
      return static_cast<real>(rtol);
    }
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override {
      this->SetRelTol(tol);
      return true;
    }

   private:
    //! \brief shell preconditioner, if any (not owned)
    mfem::Solver* shell_preconditioner = nullptr;
  };  // end of struct PETScKrylovSolver

  //! \brief preconditioner of PETSc (`PC`) configured by the options database
  struct PETScPreconditioner final : public PETScOptions,
                                     public mfem::PetscPreconditioner {
    /*!
     * \param[in] comm: communicator
     */
    explicit PETScPreconditioner(MPI_Comm comm)
        : PETScOptions(comm, "pc"),
          mfem::PetscPreconditioner(comm, this->prefix) {}
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override {
      this->markAsUsed();
      mfem::PetscPreconditioner::Mult(b, x);
    }
  };  // end of struct PETScPreconditioner

  //! \return the options of the Krylov solver of PETSc
  [[nodiscard]] static std::vector<std::string> getPETScKrylovSolverOptions() {
    auto options = getIterativeSolverParametersList();
    options.insert(options.end(), {"Options", "Preconditioner"});
    return options;
  }  // end of getPETScKrylovSolverOptions

  /*!
   * \return the generator of the Krylov solver of PETSc
   * \param[in] d: description of the solver
   */
  [[nodiscard]] static LinearSolverFactory<true>::Generator
  buildPETScKrylovSolverGenerator(const LinearSolverDescription d) {
    return [d](Context& ctx, FiniteElementSpace<true>& fespace,
               const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      const auto oopts =
          get_if<Parameters>(ctx, params, "Options", Parameters{});
      if (isInvalid(oopts)) {
        return {};
      }
      initializePETSc();
      auto s = std::make_unique<PETScKrylovSolver>(fespace.GetComm());
      if (!s->setOptions(ctx, *oopts)) {
        return {};
      }
      auto& ksp = *s;
      // the options of mfem-mgis, overridden by the ones of PETSc
      if (contains(params, Problem::SolverRelativeTolerance)) {
        const auto ov =
            get<double>(ctx, params, Problem::SolverRelativeTolerance);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetRelTol(boundRelativeTolerance(*ov));
      }
      if (contains(params, Problem::SolverAbsoluteTolerance)) {
        const auto ov =
            get<double>(ctx, params, Problem::SolverAbsoluteTolerance);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetAbsTol(*ov);
      }
      if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
        const auto ov =
            get<int>(ctx, params, Problem::SolverMaximumNumberOfIterations);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetMaxIter(*ov);
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto ov = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(ov)) {
          return {};
        }
        ksp.SetPrintLevel(*ov);
      }
      // a preconditioner of mfem-mgis is seen by PETSc as a shell
      // preconditioner
      auto oprec = makePreconditioner<true>(ctx, fespace, "PETScLinearSolver",
                                            d, params);
      if (isInvalid(oprec)) {
        return {};
      }
      if (*oprec != nullptr) {
        // PETSc would replace the preconditioner by the one given by these
        // options
        for (const auto& [n, p] : *oopts) {
          if (n.starts_with("pc_")) {
            std::ignore = ctx.registerErrorMessage(
                "the option '" + n +
                "' of the preconditioner of PETSc can't be given to a solver "
                "with a preconditioner of mfem-mgis (the options of the "
                "'PETScPreconditioner' preconditioner are given by its "
                "'Options' option)");
            return {};
          }
        }
        s->SetPreconditioner(**oprec);
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildPETScKrylovSolverGenerator

  /*!
   * \return the preconditioner of PETSc
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] opts: options (`Options`, the options of PETSc)
   */
  [[nodiscard]] static std::unique_ptr<LinearSolverPreconditioner>
  makePETScPreconditioner(Context& ctx,
                          FiniteElementSpace<true>& fespace,
                          const Parameters& opts) {
    if (!checkParameters(ctx, opts, {"Options"})) {
      return {};
    }
    const auto oopts = get_if<Parameters>(ctx, opts, "Options", Parameters{});
    if (isInvalid(oopts)) {
      return {};
    }
    initializePETSc();
    auto p = std::make_unique<PETScPreconditioner>(fespace.GetComm());
    if (!p->setOptions(ctx, *oopts)) {
      return {};
    }
    return p;
  }  // end of makePETScPreconditioner

  void declarePETScLinearSolvers(LinearSolverFactory<true>& f) {
    const auto d = LinearSolverDescription{
        .library = "PETSc",
        .kind = LinearSolverKind::ITERATIVE,
        .placement = LinearSolverPlacement::HOST,
        .preconditioners = AcceptedPreconditioners::ANY,
        .options = getPETScKrylovSolverOptions(),
        .description =
            "Krylov solver of PETSc (KSP), configured by the "
            "options of PETSc given by the 'Options' parameter"};
    declare(f, "PETScLinearSolver", d, buildPETScKrylovSolverGenerator(d));
  }  // end of declarePETScLinearSolvers

  void declarePETScPreconditioners(PreconditionerFactory<true>& f) {
    declare(f, "PETScPreconditioner",
            LinearSolverDescription{
                .library = "PETSc",
                .kind = LinearSolverKind::PRECONDITIONER,
                .placement = LinearSolverPlacement::HOST,
                .preconditioners = AcceptedPreconditioners::NONE,
                .options = {"Options"},
                .description = "preconditioner of PETSc (PC), configured by "
                               "the options of PETSc given by the 'Options' "
                               "parameter"},
            makePETScPreconditioner);
  }  // end of declarePETScPreconditioners

}  // end of namespace mfem_mgis::internals

#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC) */
