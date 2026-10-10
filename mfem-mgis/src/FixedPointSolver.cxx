/*!
 * \file   src/FixedPointSolver.cxx
 * \brief  This file implements the `FixedPointSolver` class
 * \date   10/10/2026
 */

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

  /*!
   * \return the name and the options of an object built by a factory, given
   * by its name or by a dictionary with a unique entry associating its name
   * with its options (see `extractFactoryArgument`)
   * \param[in, out] ctx: execution context
   * \param[in] params: parameters
   * \param[in] n: name of the parameter
   */
  [[nodiscard]] static std::optional<std::pair<std::string, Parameters>>
  getFactoryArgument(Context& ctx,
                     const Parameters& params,
                     const char* const n) noexcept {
    const auto op = params.get(ctx, n);
    if (isInvalid(op)) {
      return {};
    }
    if (is<std::string>(*op)) {
      const auto oname = get<std::string>(ctx, *op);
      if (isInvalid(oname)) {
        return {};
      }
      return std::pair<std::string, Parameters>{*oname, Parameters{}};
    }
    if (is<Parameters>(*op)) {
      const auto oparams = get<Parameters>(ctx, *op);
      if (isInvalid(oparams)) {
        return {};
      }
      return extractFactoryArgument(ctx, *oparams);
    }
    return ctx.registerErrorMessage(
        "FixedPointSolver: the option '" + std::string{n} +
        "' must be a name or a dictionary associating a name with options");
  }  // end of getFactoryArgument

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
    const auto options = std::vector<std::string>{
        "Method", "ElasticOperatorComputation", "ElasticOperatorUpdate",
        "ElasticOperatorLinearSolver", "Acceleration"};
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
    if (contains(params, "Method")) {
      const auto om = read("Method", {"ElasticOperator", "ModifiedNewton"});
      if (isInvalid(om)) {
        return false;
      }
      this->method = (*om == "ElasticOperator") ? Method::ELASTIC_OPERATOR
                                                : Method::MODIFIED_NEWTON;
    }
    if (contains(params, "ElasticOperatorComputation")) {
      const auto oc =
          read("ElasticOperatorComputation", {"Prediction", "Integration"});
      if (isInvalid(oc)) {
        return false;
      }
      this->elastic_operator_by_integration = *oc == "Integration";
    }
    if (contains(params, "ElasticOperatorUpdate")) {
      const auto ou = read("ElasticOperatorUpdate", {"EachResolution", "Once"});
      if (isInvalid(ou)) {
        return false;
      }
      this->elastic_operator_assembled_once = *ou == "Once";
    }
    if (contains(params, "ElasticOperatorLinearSolver")) {
      const auto oa =
          getFactoryArgument(ctx, params, "ElasticOperatorLinearSolver");
      if (isInvalid(oa)) {
        return false;
      }
      try {
        auto s = this->generate_linear_solver(ctx, oa->first, oa->second);
        if (isInvalid(s)) {
          return false;
        }
        s.linear_solver->iterative_mode = false;
        this->elastic_operator_linear_solver = std::move(s);
        this->elastic_operator.reset();
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    if (contains(params, "Acceleration")) {
      const auto oa = getFactoryArgument(ctx, params, "Acceleration");
      if (isInvalid(oa)) {
        return false;
      }
      auto a = AccelerationAlgorithmFactory::get().generate(ctx, oa->first,
                                                            oa->second);
      if (a == nullptr) {
        return false;
      }
      this->acceleration = std::move(a);
    }
    return NonLinearSolverBase::setSolverParameters(ctx, others);
  }  // end of setSolverParameters

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
    return this->processNewUnknownsEstimate(u);
  }  // end of integrate

  bool FixedPointSolver::setElasticOperator(const mfem::Vector& u) const {
    if ((this->elastic_operator == nullptr) ||
        (!this->elastic_operator_assembled_once)) {
      const auto it = this->elastic_operator_by_integration
                          ? IntegrationType::INTEGRATION_ELASTIC_OPERATOR
                          : IntegrationType::PREDICTION_ELASTIC_OPERATOR;
      if (!this->integrate(u, it)) {
        return false;
      }
      // the operator given by the problem is modified by the next
      // assemblies: it is copied
      this->elastic_operator = copyAssembledOperator(this->getJacobian(u));
      if (this->elastic_operator_linear_solver.linear_solver != nullptr) {
        auto profiler = this->ctx_ptr != nullptr
                            ? this->ctx_ptr->startNewProfiling(
                                  "NS::setLinearSolverOperator",
                                  this->ctx_ptr->isProfilingEnabled())
                            : mgis::ProfilingSection{};
        this->elastic_operator_linear_solver.linear_solver->SetOperator(
            *(this->elastic_operator));
      }
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
    // operator of the corrections and first integration
    if (this->method == Method::ELASTIC_OPERATOR) {
      if ((!this->setElasticOperator(x)) ||
          (!this->integrate(
              x, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR))) {
        return;
      }
    } else {
      if (!this->integrate(
              x, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR)) {
        return;
      }
      auto profiler_setup = this->ctx_ptr != nullptr
                                ? this->ctx_ptr->startNewProfiling(
                                      "NS::setLinearSolverOperator",
                                      this->ctx_ptr->isProfilingEnabled())
                                : mgis::ProfilingSection{};
      this->prec->SetOperator(this->getJacobian(x));
    }
    auto& solver = this->getCorrectionSolver();
    const auto n = x.Size();
    auto r = mfem::Vector(n);
    auto c = mfem::Vector(n);
    auto xn = mfem::Vector(n);
    auto dx = mfem::Vector(n);
    for (auto* const v : {&r, &c, &xn, &dx}) {
      v->UseDevice(true);
    }
    auto update_residual = [this, &r](const mfem::Vector& u) {
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
    auto it = size_type{};
    while (true) {
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
      {
        auto profiler_solve =
            this->ctx_ptr != nullptr
                ? this->ctx_ptr->startNewProfiling(
                      "MFEM::Mult(r,c)", this->ctx_ptr->isProfilingEnabled())
                : mgis::ProfilingSection{};
        c = real{0};
        solver.Mult(r, c);
        clearHypreErrors();
      }
      if ((!this->discardLinearSolverFailure) && (!hasConverged(solver))) {
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
      auto integrated =
          this->integrate(xn, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR);
      for (int i = 0; (!integrated) && (i != 10); ++i) {
        if (shall_print) {
          log << "Fixed point iteration " << std::setw(2) << it
              << ": integration failure, the estimate is moved towards the "
                 "previous one\n";
        }
        xn += x;
        xn *= real{1} / 2;
        integrated = this->integrate(
            xn, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR);
      }
      if (!integrated) {
        break;
      }
      x = xn;
      norm = update_residual(x);
      auto iteration = Parameters{};
      iteration.replaceOrInsert("Norm", norm);
      this->iterations_information.push_back(iteration);
      ++it;
    }
    this->final_iter = static_cast<int>(it);
    this->final_norm = norm;
  }  // end of Mult

  FixedPointSolver::~FixedPointSolver() noexcept = default;

}  // end of namespace mfem_mgis
