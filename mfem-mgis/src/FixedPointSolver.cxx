/*!
 * \file   src/FixedPointSolver.cxx
 * \brief  This file implements the `FixedPointSolver` class
 * \date   10/10/2026
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <iomanip>
#include <utility>
#include <optional>
#include <algorithm>
#include "mfem/linalg/sparsemat.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearSolvers/FixedPointSolver.hxx"

namespace mfem_mgis {

  /*!
   * \return a copy of an assembled operator, or of the assembled version of
   * an operator applied without assembly (see `AssembledOperatorProvider`)
   * \param[in] op: operator
   */
  [[nodiscard]] static std::unique_ptr<mfem::Operator> copyAssembledOperator(
      const mfem::Operator& op) {
    const auto* a = &op;
    if (const auto* const p = dynamic_cast<const AssembledOperatorProvider*>(a);
        p != nullptr) {
      a = &(p->getAssembledOperator());
    }
    if (const auto* const m = dynamic_cast<const mfem::SparseMatrix*>(a);
        m != nullptr) {
      return std::make_unique<mfem::SparseMatrix>(*m);
    }
#ifdef MFEM_USE_MPI
    if (const auto* const m = dynamic_cast<const mfem::HypreParMatrix*>(a);
        m != nullptr) {
      return std::make_unique<mfem::HypreParMatrix>(*m);
    }
#endif /* MFEM_USE_MPI */
    raise(
        "FixedPointSolver: the elastic operator is neither a sparse matrix "
        "nor a parallel matrix of hypre");
  }  // end of copyAssembledOperator

  namespace {

    /*!
     * \brief preconditioner applying a linear solver whose setup is kept:
     * the operators given to `SetOperator` are ignored
     */
    struct KeptLinearSolver final : mfem::Solver {
      //! \param[in] s: linear solver
      explicit KeptLinearSolver(mfem::Solver& s)
          : mfem::Solver(s.Height(), s.Width()), solver(s) {}
      void SetOperator(const mfem::Operator&) override {}
      void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
        this->solver.Mult(x, y);
      }

     private:
      //! \brief linear solver
      mfem::Solver& solver;
    };

  }  // end of namespace

  /*!
   * \brief builds the functions of the fixed point solver which depend on
   * the type of the problem
   * \param[out] integrate: integration of the behaviours
   * \param[out] generate: generator of linear solvers
   * \param[in] p: non linear evolution problem
   */
  template <bool parallel>
  static void initializeFixedPointSolverFunctions(
      std::function<bool(const mfem::Vector&, const IntegrationType)>&
          integrate,
      std::function<LinearSolverHandler(
          Context&, std::string_view, const Parameters&)>& generate,
      NonLinearEvolutionProblemImplementation<parallel>& p) {
    integrate = [&p](const mfem::Vector& u, const IntegrationType it) {
      return p.integrate(u, it, {});
    };
    generate = [&p](Context& ctx, std::string_view n,
                    const Parameters& opts) -> LinearSolverHandler {
      auto& fespace = p.getFiniteElementDiscretization()
                          .template getFiniteElementSpace<parallel>();
      return LinearSolverFactory<parallel>::getFactory().generate(
          ctx, n, fespace, opts);
    };
  }  // end of initializeFixedPointSolverFunctions

#ifdef MFEM_USE_MPI

  FixedPointSolver::FixedPointSolver(
      NonLinearEvolutionProblemImplementation<true>& p)
      : NonLinearSolverBase(p) {
    // the behaviours are integrated by this class, with the required
    // integration type
    this->nue_actions.clear();
    initializeFixedPointSolverFunctions(this->integrate_behaviours,
                                        this->generate_linear_solver, p);
  }  // end of FixedPointSolver

#endif /* MFEM_USE_MPI */

  FixedPointSolver::FixedPointSolver(
      NonLinearEvolutionProblemImplementation<false>& p)
      : NonLinearSolverBase(p) {
    this->nue_actions.clear();
    initializeFixedPointSolverFunctions(this->integrate_behaviours,
                                        this->generate_linear_solver, p);
  }  // end of FixedPointSolver

  bool FixedPointSolver::setSolverParameters(
      Context& ctx, const Parameters& params) noexcept {
    const auto options =
        std::vector<std::string>{"Method",
                                 "ElasticOperatorComputation",
                                 "ElasticOperatorUpdate",
                                 "ElasticOperatorLinearSolver",
                                 "Acceleration",
                                 "PredictionLinearSolver",
                                 "TangentLinearSolver",
                                 "AdaptiveLinearSolverTolerance"};
    auto others = Parameters{};
    for (const auto& [k, v] : params) {
      if (std::find(options.begin(), options.end(), k) == options.end()) {
        if (!others.insert(ctx, k, v)) {
          return false;
        }
      }
    }
    // reads a value among the allowed ones
    auto read = [&ctx, &params](const char* const n,
                                const std::vector<std::string>& values)
        -> std::optional<std::string> {
      const auto ov = get<std::string>(ctx, params, n);
      if (isInvalid(ov)) {
        return {};
      }
      if (std::find(values.begin(), values.end(), *ov) == values.end()) {
        auto msg = "FixedPointSolver: invalid value '" + *ov +
                   "' for option '" + std::string{n} + "' (expected";
        for (const auto& value : values) {
          msg += " '" + value + "'";
        }
        return ctx.registerErrorMessage(msg + ")");
      }
      return *ov;
    };
    // the options are read and checked before being applied: a refused call
    // leaves the solver unchanged
    auto new_method = this->method;
    auto new_by_integration = this->elastic_operator_by_integration;
    auto new_once = this->elastic_operator_assembled_once;
    auto new_linear_solver = std::optional<LinearSolverHandler>{};
    auto new_acceleration = std::unique_ptr<AccelerationAlgorithm>{};
    if (contains(params, "Method")) {
      const auto om =
          read("Method", {"ElasticOperator", "ModifiedNewton", "NewtonK0"});
      if (isInvalid(om)) {
        return false;
      }
      new_method = (*om == "ElasticOperator")  ? Method::ELASTIC_OPERATOR
                   : (*om == "ModifiedNewton") ? Method::MODIFIED_NEWTON
                                               : Method::NEWTON_K0;
    }
    if (contains(params, "ElasticOperatorComputation")) {
      const auto oc =
          read("ElasticOperatorComputation", {"Prediction", "Integration"});
      if (isInvalid(oc)) {
        return false;
      }
      new_by_integration = *oc == "Integration";
    }
    if (contains(params, "ElasticOperatorUpdate")) {
      const auto ou = read("ElasticOperatorUpdate", {"EachResolution", "Once"});
      if (isInvalid(ou)) {
        return false;
      }
      new_once = *ou == "Once";
    }
    // reads the name and the options of an object built by a factory; the
    // value `None` removes the object and takes no option
    auto read_factory_argument = [&ctx, &params](const char* const n)
        -> std::optional<std::pair<std::string, Parameters>> {
      auto oa = getFactoryArgument(ctx, params, n);
      if (isInvalid(oa)) {
        return {};
      }
      if ((oa->first == "None") && (!oa->second.empty())) {
        return ctx.registerErrorMessage(
            "FixedPointSolver: no option expected "
            "with the value 'None' of option '" +
            std::string{n} + "'");
      }
      return oa;
    };
    auto remove_linear_solver = false;
    if (contains(params, "ElasticOperatorLinearSolver")) {
      const auto oa = read_factory_argument("ElasticOperatorLinearSolver");
      if (isInvalid(oa)) {
        return false;
      }
      remove_linear_solver = oa->first == "None";
      if (!remove_linear_solver) {
        try {
          auto ls = this->generate_linear_solver(ctx, oa->first, oa->second);
          if (isInvalid(ls)) {
            return false;
          }
          ls.linear_solver->iterative_mode = false;
          new_linear_solver = std::move(ls);
        } catch (...) {
          return registerExceptionInErrorBacktrace(ctx);
        }
      }
    }
    auto remove_acceleration = false;
    if (contains(params, "Acceleration")) {
      const auto oa = read_factory_argument("Acceleration");
      if (isInvalid(oa)) {
        return false;
      }
      remove_acceleration = oa->first == "None";
      if (!remove_acceleration) {
        new_acceleration = AccelerationAlgorithmFactory::get().generate(
            ctx, oa->first, oa->second);
        if (new_acceleration == nullptr) {
          return false;
        }
      }
    }
    // the iterative solvers of the predictions and of the tangent operator
    // are given by the factory of the linear solvers
    auto new_prediction = std::optional<LinearSolverHandler>{};
    if (contains(params, "PredictionLinearSolver")) {
      const auto oa = getFactoryArgument(ctx, params, "PredictionLinearSolver");
      if (isInvalid(oa)) {
        return false;
      }
      const auto& [n, opts] = *oa;
      if (n == "Problem") {
        if (!opts.empty()) {
          return ctx.registerErrorMessage(
              "FixedPointSolver: no option expected with the value 'Problem' "
              "of option 'PredictionLinearSolver'");
        }
        new_prediction = LinearSolverHandler{};
      } else {
        auto s = this->generateIterativeSolver(
            ctx, "PredictionLinearSolver",
            n == "ElasticOperator" ? "CGSolver" : n, opts, 1e-8);
        if (isInvalid(s)) {
          return false;
        }
        new_prediction = std::move(s);
      }
    }
    auto new_tangent = std::optional<LinearSolverHandler>{};
    if (contains(params, "TangentLinearSolver")) {
      const auto oa = getFactoryArgument(ctx, params, "TangentLinearSolver");
      if (isInvalid(oa)) {
        return false;
      }
      auto s = this->generateIterativeSolver(ctx, "TangentLinearSolver",
                                             oa->first, oa->second, 1e-4);
      if (isInvalid(s)) {
        return false;
      }
      new_tangent = std::move(s);
    } else if ((new_method == Method::NEWTON_K0) &&
               (this->tangent_linear_solver.linear_solver == nullptr)) {
      auto s = this->generateIterativeSolver(ctx, "TangentLinearSolver",
                                             "CGSolver", {}, 1e-4);
      if (isInvalid(s)) {
        return false;
      }
      new_tangent = std::move(s);
    }
    auto new_adaptive_tolerance = this->adaptive_linear_solver_tolerance;
    if (contains(params, "AdaptiveLinearSolverTolerance")) {
      const auto ot = AdaptiveLinearSolverTolerance::read(ctx, params);
      if (isInvalid(ot)) {
        return false;
      }
      new_adaptive_tolerance = *ot;
    }
    if ((new_method != Method::NEWTON_K0) &&
        (new_adaptive_tolerance.isEnabled())) {
      return ctx.registerErrorMessage(
          "FixedPointSolver: the adaptive tolerance of the linear solver "
          "(option 'AdaptiveLinearSolverTolerance') is only used by the "
          "'NewtonK0' method (set it to 'None' to change the method)");
    }
    const auto has_linear_solver =
        (!remove_linear_solver) &&
        (new_linear_solver.has_value() ||
         (this->elastic_operator_linear_solver.linear_solver != nullptr));
    if ((new_method == Method::NEWTON_K0) && (!has_linear_solver)) {
      return ctx.registerErrorMessage(
          "FixedPointSolver: the 'NewtonK0' method requires a linear solver "
          "of the elastic operator (option 'ElasticOperatorLinearSolver')");
    }
    const auto has_acceleration =
        (!remove_acceleration) &&
        ((new_acceleration != nullptr) || (this->acceleration != nullptr));
    if ((new_method == Method::NEWTON_K0) && has_acceleration) {
      return ctx.registerErrorMessage(
          "FixedPointSolver: the 'NewtonK0' method can't be accelerated");
    }
    if ((new_method == Method::MODIFIED_NEWTON) && has_linear_solver) {
      return ctx.registerErrorMessage(
          "FixedPointSolver: the linear solver of the elastic operator "
          "(option 'ElasticOperatorLinearSolver') is not used by the "
          "'ModifiedNewton' method");
    }
    if (!NonLinearSolverBase::setSolverParameters(ctx, others)) {
      return false;
    }
    // the preconditioners applying the linear solver of the elastic
    // operator are removed before this linear solver is replaced
    const auto new_elastic_operator_linear_solver =
        new_linear_solver.has_value() || remove_linear_solver;
    if (new_elastic_operator_linear_solver || new_prediction.has_value()) {
      this->prediction_preconditioner.reset();
    }
    if (new_elastic_operator_linear_solver || new_tangent.has_value()) {
      this->tangent_preconditioner.reset();
    }
    // the elastic operator is computed again if its computation or its
    // linear solver changed
    if ((new_method != this->method) ||
        (new_by_integration != this->elastic_operator_by_integration) ||
        new_elastic_operator_linear_solver) {
      this->elastic_operator.reset();
    }
    this->method = new_method;
    this->adaptive_linear_solver_tolerance = new_adaptive_tolerance;
    this->elastic_operator_by_integration = new_by_integration;
    this->elastic_operator_assembled_once = new_once;
    if (new_linear_solver.has_value()) {
      this->elastic_operator_linear_solver = std::move(*new_linear_solver);
    }
    if (remove_linear_solver) {
      this->elastic_operator_linear_solver = LinearSolverHandler{};
    }
    if (new_acceleration != nullptr) {
      this->acceleration = std::move(new_acceleration);
    }
    if (remove_acceleration) {
      this->acceleration.reset();
    }
    if (new_prediction.has_value()) {
      this->prediction_linear_solver = std::move(*new_prediction);
    }
    if (new_tangent.has_value()) {
      this->tangent_linear_solver = std::move(*new_tangent);
    }
    return true;
  }  // end of setSolverParameters

  LinearSolverHandler FixedPointSolver::generateIterativeSolver(
      Context& ctx,
      std::string_view o,
      std::string_view n,
      const Parameters& opts,
      const real tol) const noexcept {
    if (contains(opts, "Preconditioner")) {
      return ctx.registerErrorMessage(
          "FixedPointSolver: the iterative solver of option '" +
          std::string{o} +
          "' is preconditioned by the linear solver of the elastic operator "
          "(option 'Preconditioner' not allowed)");
    }
    try {
      // default values of the tolerances and of the number of iterations
      auto options = Parameters{{"RelativeTolerance", tol},
                                {"AbsoluteTolerance", real{0}},
                                {"MaximumNumberOfIterations", 1000}};
      for (const auto& [k, v] : opts) {
        options.replaceOrInsert(k, v);
      }
      auto s = this->generate_linear_solver(ctx, n, options);
      if (isInvalid(s)) {
        return {};
      }
      if (dynamic_cast<mfem::IterativeSolver*>(s.linear_solver.get()) ==
          nullptr) {
        return ctx.registerErrorMessage(
            "FixedPointSolver: '" + std::string{n} + "' (option '" +
            std::string{o} +
            "') is not an iterative solver of MFEM (expected for instance "
            "'CGSolver', 'GMRESSolver' or 'BiCGSTABSolver')");
      }
      s.linear_solver->iterative_mode = false;
      return s;
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
  }  // end of generateIterativeSolver

  mfem::Solver& FixedPointSolver::getPreconditionedSolver(
      const LinearSolverHandler& s, std::unique_ptr<mfem::Solver>& p) const {
    if (p == nullptr) {
      auto& is = dynamic_cast<mfem::IterativeSolver&>(*(s.linear_solver));
      p = std::make_unique<KeptLinearSolver>(
          *(this->elastic_operator_linear_solver.linear_solver));
      is.SetPreconditioner(*p);
    }
    return *(s.linear_solver);
  }  // end of getPreconditionedSolver

  LinearSolver* FixedPointSolver::getElasticPredictionLinearSolver() noexcept {
    if ((this->prediction_linear_solver.linear_solver == nullptr) ||
        (this->method == Method::MODIFIED_NEWTON) ||
        (this->elastic_operator_linear_solver.linear_solver == nullptr) ||
        (this->elastic_operator == nullptr)) {
      return nullptr;
    }
    try {
      return &(this->getPreconditionedSolver(this->prediction_linear_solver,
                                             this->prediction_preconditioner));
    } catch (...) {
      this->prediction_preconditioner.reset();
      return nullptr;
    }
  }  // end of getElasticPredictionLinearSolver

  bool FixedPointSolver::integrate(const mfem::Vector& u,
                                   const IntegrationType it) const {
    auto profiler =
        this->ctx_ptr != nullptr
            ? this->ctx_ptr->startNewProfiling(
                  "FPS::integrate", this->ctx_ptr->isProfilingEnabled())
            : mgis::ProfilingSection{};
    if (!this->integrate_behaviours(u, it)) {
      return false;
    }
    if (this->nue_actions.empty()) {
      return true;
    }
    return this->processNewUnknownsEstimate(u);
  }  // end of integrate

  bool FixedPointSolver::setElasticOperator(const mfem::Vector& u) const {
    if ((this->elastic_operator == nullptr) ||
        (!this->elastic_operator_assembled_once)) {
      const auto it = this->elastic_operator_by_integration
                          ? IntegrationType::INTEGRATION_ELASTIC_OPERATOR
                          : IntegrationType::PREDICTION_ELASTIC_OPERATOR;
      // the actions associated with a new estimate are performed by the
      // integration which follows (see `Mult`)
      auto profiler_integration =
          this->ctx_ptr != nullptr
              ? this->ctx_ptr->startNewProfiling(
                    "FPS::integrate", this->ctx_ptr->isProfilingEnabled())
              : mgis::ProfilingSection{};
      if (!this->integrate_behaviours(u, it)) {
        return false;
      }
      // the operator given by the problem is modified by the next
      // assemblies: it is copied. The linear solver is set up on the new
      // operator before the previous one is destroyed; on failure, the
      // elastic operator is computed again at the next resolution.
      auto K0 = copyAssembledOperator(this->getJacobian(u));
      if (this->elastic_operator_linear_solver.linear_solver != nullptr) {
        auto profiler = this->ctx_ptr != nullptr
                            ? this->ctx_ptr->startNewProfiling(
                                  "NS::setLinearSolverOperator",
                                  this->ctx_ptr->isProfilingEnabled())
                            : mgis::ProfilingSection{};
        try {
          this->elastic_operator_linear_solver.linear_solver->SetOperator(*K0);
        } catch (...) {
          this->elastic_operator.reset();
          throw;
        }
      }
      this->elastic_operator = std::move(K0);
    }
    if (this->elastic_operator_linear_solver.linear_solver == nullptr) {
      // the prediction of the problem gives its own operator to its linear
      // solver
      auto profiler = this->ctx_ptr != nullptr
                          ? this->ctx_ptr->startNewProfiling(
                                "NS::setLinearSolverOperator",
                                this->ctx_ptr->isProfilingEnabled())
                          : mgis::ProfilingSection{};
      this->prec->SetOperator(*(this->elastic_operator));
    }
    return true;
  }  // end of setElasticOperator

  mfem::Solver& FixedPointSolver::getCorrectionSolver() const {
    if (this->method == Method::NEWTON_K0) {
      return this->getPreconditionedSolver(this->tangent_linear_solver,
                                           this->tangent_preconditioner);
    }
    if ((this->method == Method::ELASTIC_OPERATOR) &&
        (this->elastic_operator_linear_solver.linear_solver != nullptr)) {
      return *(this->elastic_operator_linear_solver.linear_solver);
    }
    return *(this->prec);
  }  // end of getCorrectionSolver

  void FixedPointSolver::Mult(const mfem::Vector&, mfem::Vector& x) const {
    auto profiler = this->ctx_ptr != nullptr
                        ? this->ctx_ptr->startNewProfiling(
                              "FPS::Mult", this->ctx_ptr->isProfilingEnabled())
                        : mgis::ProfilingSection{};
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");
    this->iterations_information.clear();
    this->converged = false;
    this->final_iter = 0;
    this->final_norm = std::numeric_limits<real>::max();
    auto& log = this->getLogStream();
    const auto shall_print = this->shallPrint();
    // operator of the corrections and first integration. The failures of
    // the assembly and of the setup of the linear solvers are handled as by
    // the Newton solver (see `callAndHandleFailures`)
    // integration type at each new estimate: the tangent operator is only
    // required by the `NewtonK0` method
    const auto it_type =
        (this->method == Method::NEWTON_K0)
            ? IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR
            : IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR;
    auto* solver_ptr = static_cast<mfem::Solver*>(nullptr);
    // relative tolerance of the iterative solver of the tangent operator
    // adapted to the decrease of the residual (`NewtonK0`), restored at the
    // end of the resolution
    auto adaptive_tolerance =
        std::optional<AdaptiveLinearSolverToleranceUpdater>{};
    const auto initialized = this->callAndHandleFailures([this, &x, it_type,
                                                          &solver_ptr,
                                                          &adaptive_tolerance] {
      if (this->method != Method::MODIFIED_NEWTON) {
        if (!(this->setElasticOperator(x) && this->integrate(x, it_type))) {
          return false;
        }
        solver_ptr = &(this->getCorrectionSolver());
        if ((this->method == Method::NEWTON_K0) &&
            (this->adaptive_linear_solver_tolerance.isEnabled())) {
          adaptive_tolerance.emplace(this->adaptive_linear_solver_tolerance,
                                     *solver_ptr);
        }
        return true;
      }
      if (!this->integrate(
              x, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR)) {
        return false;
      }
      auto& J = this->getJacobian(x);
      auto profiler_setup = this->ctx_ptr != nullptr
                                ? this->ctx_ptr->startNewProfiling(
                                      "NS::setLinearSolverOperator",
                                      this->ctx_ptr->isProfilingEnabled())
                                : mgis::ProfilingSection{};
      this->prec->SetOperator(J);
      solver_ptr = &(this->getCorrectionSolver());
      return true;
    });
    if (!initialized) {
      return;
    }
    auto& solver = *solver_ptr;
    const auto n = x.Size();
    auto r = mfem::Vector(n);
    auto c = mfem::Vector(n);
    auto xn = mfem::Vector(n);
    auto dx = mfem::Vector(n);
    for (auto* const v : {&r, &c, &xn, &dx}) {
      v->UseDevice(true);
    }
    auto update_residual = [this, &r](const mfem::Vector& u) {
      auto profiler_residual =
          this->ctx_ptr != nullptr
              ? this->ctx_ptr->startNewProfiling(
                    "NS::computeResidual", this->ctx_ptr->isProfilingEnabled())
              : mgis::ProfilingSection{};
      this->computeResidual(r, u);
      return this->Norm(r);
    };
    auto norm = update_residual(x);
    if (!this->reference_residual_norm.has_value()) {
      this->reference_residual_norm = norm;
    }
    this->initial_norm = *(this->reference_residual_norm);
    const auto goal = std::max(
        this->rel_tol * (*(this->reference_residual_norm)), this->abs_tol);
    if (this->acceleration != nullptr) {
      this->acceleration->initialize(
          static_cast<size_type>(n),
          [this](const mfem::Vector& a, const mfem::Vector& b) {
            return this->Dot(a, b);
          });
    }
    // scale of the variations of the unknowns at convergence, estimated
    // from the first correction (see `AccelerationAlgorithm::execute`)
    auto unknowns_scale = real{};
    // relative tolerance of the last linear solve, if adapted, and norm of
    // its residual, if computed (see `AdaptiveLinearSolverTolerance`)
    auto linear_solver_tolerance = std::optional<real>{};
    auto linear_residual_norm = std::optional<real>{};
    // information on an iteration (see `NewtonSolver`)
    auto add_iteration_information = [this, &solver, &linear_solver_tolerance,
                                      &linear_residual_norm](const real nr,
                                                             const bool ls_cv) {
      auto iteration = Parameters{};
      iteration.replaceOrInsert("Norm", nr);
      const auto ls_iterations = getNumberOfIterationsAtConvergence(solver);
      if (isValid(ls_iterations)) {
        auto ls = Parameters{};
        ls.replaceOrInsert("NumberOfIterations", ls_iterations.value());
        ls.replaceOrInsert("HasConverged", ls_cv);
        if (linear_solver_tolerance.has_value()) {
          ls.replaceOrInsert("RelativeTolerance", *linear_solver_tolerance);
        }
        if (linear_residual_norm.has_value()) {
          ls.replaceOrInsert("ResidualNorm", *linear_residual_norm);
        }
        iteration.replaceOrInsert("LinearSolver", ls);
      }
      this->iterations_information.push_back(iteration);
    };
    auto it = size_type{};
    while (true) {
      if (!std::isfinite(norm)) {
        if (shall_print) {
          log << "Fixed point iteration " << std::setw(2) << it
              << ": the norm of the residual is not finite\n";
        }
        break;
      }
      if (shall_print) {
        log << "Fixed point iteration " << std::setw(2) << it
            << " : ||r|| = " << norm
            << ", ||r||/||r_0|| = " << norm / (*(this->reference_residual_norm))
            << '\n';
      }
      this->Monitor(static_cast<int>(it), norm, r, x);
      if (norm <= goal) {
        this->converged = true;
        break;
      }
      if (it >= static_cast<size_type>(this->max_iter)) {
        break;
      }
      // correction
      if (adaptive_tolerance.has_value()) {
        linear_solver_tolerance = adaptive_tolerance->update(it, norm, goal);
        if (shall_print) {
          log << "Fixed point iteration " << std::setw(2) << it
              << " : relative tolerance of the linear solver "
              << *linear_solver_tolerance << '\n';
        }
      }
      const auto solved = this->callAndHandleFailures(
          [this, &solver, &r, &c, &x, &adaptive_tolerance] {
            auto* J = static_cast<mfem::Operator*>(nullptr);
            if (this->method == Method::NEWTON_K0) {
              // the tangent operator is given to the iterative solver, whose
              // preconditioner is kept
              J = &(this->getJacobian(x));
              auto profiler_setup =
                  this->ctx_ptr != nullptr
                      ? this->ctx_ptr->startNewProfiling(
                            "NS::setLinearSolverOperator",
                            this->ctx_ptr->isProfilingEnabled())
                      : mgis::ProfilingSection{};
              solver.SetOperator(*J);
            }
            auto profiler_solve = this->ctx_ptr != nullptr
                                      ? this->ctx_ptr->startNewProfiling(
                                            "MFEM::Mult(r,c)",
                                            this->ctx_ptr->isProfilingEnabled())
                                      : mgis::ProfilingSection{};
            c = real{0};
            solver.Mult(r, c);
            clearHypreErrors();
            if ((J != nullptr) && (adaptive_tolerance.has_value()) &&
                (adaptive_tolerance->requiresLinearResidualNorm())) {
              // norm of the residual of the linear solve, J c - r
              auto lr = mfem::Vector(r.Size());
              lr.UseDevice(true);
              J->Mult(c, lr);
              lr -= r;
              adaptive_tolerance->setLinearResidualNorm(this->Norm(lr));
            }
            return true;
          });
      if (adaptive_tolerance.has_value()) {
        linear_residual_norm = adaptive_tolerance->getLinearResidualNorm();
      }
      if (!solved) {
        break;
      }
      if ((!this->discardLinearSolverFailure) && (!hasConverged(solver))) {
        add_iteration_information(norm, false);
        break;
      }
      if (it == 0) {
        unknowns_scale = this->rel_tol * this->Norm(c);
      }
      subtract(x, c, xn);
      if (this->acceleration != nullptr) {
        auto profiler_acceleration =
            this->ctx_ptr != nullptr
                ? this->ctx_ptr->startNewProfiling(
                      "FPS::accelerate", this->ctx_ptr->isProfilingEnabled())
                : mgis::ProfilingSection{};
        subtract(xn, x, dx);
        this->acceleration->execute(xn, dx, r, goal, unknowns_scale, it + 1);
      }
      // integration at the new estimate: on failure, the estimate is moved
      // towards the previous one
      auto integrated = this->integrate(xn, it_type);
      for (int i = 0; (!integrated) && (i != 10); ++i) {
        if (shall_print) {
          log << "Fixed point iteration " << std::setw(2) << it
              << ": integration failure, the estimate is moved towards the "
                 "previous one\n";
        }
        xn += x;
        xn *= real{1} / 2;
        integrated = this->integrate(xn, it_type);
      }
      if (!integrated) {
        break;
      }
      x = xn;
      norm = update_residual(x);
      add_iteration_information(norm, true);
      ++it;
    }
    this->final_iter = static_cast<int>(it);
    this->final_norm = norm;
  }  // end of Mult

  FixedPointSolver::~FixedPointSolver() noexcept = default;

}  // end of namespace mfem_mgis
