/*!
 * \file   src/GinkgoLinearSolvers.cxx
 * \brief  This file implements the module of the linear solvers and of the
 * preconditioners of Ginkgo (see `LinearSolverModules.hxx`). Each entry is
 * described by a configuration of Ginkgo (see `gko::config::parse`) and
 * built on the classes of MFEM wrapping Ginkgo (`mfem::Ginkgo`), which
 * convert the matrices and the vectors of MFEM.
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#ifdef MFEM_MGIS_HAVE_GINKGO

#include <array>
#include <memory>
#include <tuple>
#include <string>
#include <vector>
#include <cstdint>
#include <cstring>
#include <optional>
#include <utility>
#include <algorithm>
#include <exception>
#include <string_view>
#include <ginkgo/ginkgo.hpp>
#include "mfem/general/device.hpp"
#include "mfem/linalg/ginkgo.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverConvergenceReport.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "LinearSolverModules.hxx"

// the distributed solvers of Ginkgo require Ginkgo built with MPI
#if defined(MFEM_USE_MPI) && GINKGO_BUILD_MPI
#define MFEM_MGIS_GINKGO_DISTRIBUTED
#endif

namespace mfem_mgis::internals {

  //! \brief executor of Ginkgo, as wrapped by MFEM
  using GinkgoExecutor = mfem::Ginkgo::GinkgoExecutor;
  //! \brief configuration of Ginkgo
  using GinkgoConfiguration = gko::config::pnode;
  //! \brief map of a configuration of Ginkgo
  using GinkgoConfigurationMap = gko::config::pnode::map_type;

  //! \brief name of the preconditioner given to a solver of Ginkgo
  static constexpr const char* ginkgo_preconditioner_name =
      "mfem_mgis_preconditioner";
  //! \brief name of the stopping criteria given to a solver of Ginkgo
  static constexpr const char* ginkgo_criteria_name = "mfem_mgis_criteria";

  /*!
   * \return the executor of Ginkgo matching the device of MFEM: the CUDA or
   * HIP executor when MFEM computes on a GPU, the OpenMP executor when MFEM
   * computes with OpenMP (and Ginkgo supports it), the reference executor
   * otherwise, the host executor of the GPU executors being the reference
   * executor. This is the choice of the constructor of `GinkgoExecutor`
   * taking the device of MFEM when MFEM is built with MPI; this constructor
   * is not used here since it requires a device of MFEM, whose destructor
   * finalizes hypre. The executor is shared by the solvers and the
   * preconditioners alive, so that Ginkgo copies no data between executors
   * (the factories of linear solvers are not thread-safe).
   */
  [[nodiscard]] static std::shared_ptr<GinkgoExecutor> getGinkgoExecutor() {
    static auto shared = std::weak_ptr<GinkgoExecutor>{};
    if (auto e = shared.lock(); e != nullptr) {
      return e;
    }
    using mfem::Backend;
    using mfem::Device;
    auto e = std::shared_ptr<GinkgoExecutor>{};
    if (Device::Allows(Backend::CUDA_MASK)) {
      e = std::make_shared<GinkgoExecutor>(GinkgoExecutor::CUDA,
                                           GinkgoExecutor::REFERENCE);
    } else if (Device::Allows(Backend::HIP_MASK)) {
      e = std::make_shared<GinkgoExecutor>(GinkgoExecutor::HIP,
                                           GinkgoExecutor::REFERENCE);
    } else if ((Device::Allows(Backend::OMP_MASK)) &&
               (std::strcmp(gko::version_info::get().omp_version.tag,
                            "not compiled") != 0)) {
      e = std::make_shared<GinkgoExecutor>(GinkgoExecutor::OMP);
    } else {
      e = std::make_shared<GinkgoExecutor>(GinkgoExecutor::REFERENCE);
    }
    shared = e;
    return e;
  }  // end of getGinkgoExecutor

  /*!
   * \return the types used by Ginkgo: the floating-point type of MFEM, and
   * the index types of the matrices converted by MFEM (the ones of hypre in
   * parallel)
   */
  template <bool parallel>
  [[nodiscard]] static gko::config::type_descriptor getGinkgoTypes() {
    if constexpr (parallel) {
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
      return gko::config::make_type_descriptor<
          mfem::real_t, mfem::Ginkgo::gko_hypre_int,
          mfem::Ginkgo::gko_hypre_bigint>();
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
    }
    return gko::config::make_type_descriptor<mfem::real_t, int>();
  }  // end of getGinkgoTypes

  /*!
   * \return the configuration of Ginkgo given by a parameter, or an empty
   * value on error: a boolean, an integer, a floating-point value, a
   * string, a list or a map of parameters
   * \param[in, out] ctx: execution context
   * \param[in] n: key of the parameter
   * \param[in] p: parameter
   */
  [[nodiscard]] static std::optional<GinkgoConfiguration>
  convertToGinkgoConfiguration(Context& ctx,
                               const std::string& n,
                               const Parameter& p) {
    // the integers are tested before the floating-point values, since
    // `is<double>` is also true for an integer
    if (mfem_mgis::is<bool>(p)) {
      const auto ov = mfem_mgis::get<bool>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return GinkgoConfiguration{static_cast<bool>(*ov)};
    }
    if (mfem_mgis::is<int>(p)) {
      const auto ov = mfem_mgis::get<int>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return GinkgoConfiguration{static_cast<std::int64_t>(*ov)};
    }
    if (mfem_mgis::is<double>(p)) {
      const auto ov = mfem_mgis::get<double>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return GinkgoConfiguration{static_cast<double>(*ov)};
    }
    if (mfem_mgis::is<std::string>(p)) {
      const auto ov = mfem_mgis::get<std::string>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      return GinkgoConfiguration{std::string{*ov}};
    }
    if (mfem_mgis::is<std::vector<Parameter>>(p)) {
      const auto ov = mfem_mgis::get<std::vector<Parameter>>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      auto values = GinkgoConfiguration::array_type{};
      for (const auto& v : *ov) {
        auto oc = convertToGinkgoConfiguration(ctx, n, v);
        if (isInvalid(oc)) {
          return {};
        }
        values.push_back(std::move(*oc));
      }
      return GinkgoConfiguration{values};
    }
    if (mfem_mgis::is<Parameters>(p)) {
      const auto ov = mfem_mgis::get<Parameters>(ctx, p);
      if (isInvalid(ov)) {
        return {};
      }
      auto values = GinkgoConfigurationMap{};
      for (const auto& [k, v] : *ov) {
        auto oc = convertToGinkgoConfiguration(ctx, k, v);
        if (isInvalid(oc)) {
          return {};
        }
        values.emplace(k, std::move(*oc));
      }
      return GinkgoConfiguration{values};
    }
    std::ignore = ctx.registerErrorMessage(
        "invalid value of the key '" + n +
        "' of the configuration of Ginkgo (a boolean, an integer, a "
        "floating-point value, a string, a list or a map is expected)");
    return {};
  }  // end of convertToGinkgoConfiguration

  /*!
   * \brief add the keys given by the `Configuration` parameter to a
   * configuration of Ginkgo
   * \return true on success, false if a key is already defined or on error
   * \param[in, out] ctx: execution context
   * \param[in, out] c: configuration
   * \param[in] params: parameters of the solver or of the preconditioner
   */
  [[nodiscard]] static bool addConfiguration(Context& ctx,
                                             GinkgoConfigurationMap& c,
                                             const Parameters& params) {
    if (!contains(params, "Configuration")) {
      return true;
    }
    const auto oc = get<Parameters>(ctx, params, "Configuration");
    if (isInvalid(oc)) {
      return false;
    }
    for (const auto& [k, v] : *oc) {
      if (c.contains(k)) {
        return ctx.registerErrorMessage(
            "the key '" + k +
            "' of the configuration of Ginkgo is already defined");
      }
      auto ov = convertToGinkgoConfiguration(ctx, k, v);
      if (isInvalid(ov)) {
        return false;
      }
      c.emplace(k, std::move(*ov));
    }
    return true;
  }  // end of addConfiguration

  /*!
   * \return true if a configuration of Ginkgo, or one of its nested
   * configurations, uses the L1 smoother of the Schwarz preconditioner,
   * which requires a matrix whose diagonal blocks are sorted (see
   * `mfem::Ginkgo::SchwarzPreconditioner`)
   * \param[in] c: configuration
   */
  [[nodiscard]] static bool usesL1Smoother(const GinkgoConfiguration& c) {
    using tag_t = GinkgoConfiguration::tag_t;
    if (c.get_tag() == tag_t::map) {
      for (const auto& [k, v] : c.get_map()) {
        if ((k == "l1_smoother") && (v.get_tag() == tag_t::boolean) &&
            (v.get_boolean())) {
          return true;
        }
        if (usesL1Smoother(v)) {
          return true;
        }
      }
    } else if (c.get_tag() == tag_t::array) {
      for (const auto& v : c.get_array()) {
        if (usesL1Smoother(v)) {
          return true;
        }
      }
    }
    return false;
  }  // end of usesL1Smoother

  /*!
   * \return true if a type of Ginkgo is the one of an iterative solver,
   * which reports its convergence
   * \param[in] type: type
   */
  [[nodiscard]] static bool isGinkgoIterativeSolver(const std::string& type) {
    static const auto types = std::array<const char*, 14>{
        "solver::Cg",        "solver::Fcg",      "solver::PipeCg",
        "solver::Bicg",      "solver::Bicgstab", "solver::Cgs",
        "solver::Gmres",     "solver::CbGmres",  "solver::Gcr",
        "solver::Idr",       "solver::Minres",   "solver::Ir",
        "solver::Chebyshev", "solver::Multigrid"};
    return std::find(types.begin(), types.end(), type) != types.end();
  }  // end of isGinkgoIterativeSolver

  /*!
   * \return the factory of Ginkgo described by a configuration, or a null
   * pointer on error
   * \param[in, out] ctx: execution context
   * \param[in] c: configuration
   * \param[in] r: registry of the objects referenced by the configuration
   * \param[in] e: executor
   */
  template <bool parallel>
  [[nodiscard]] static std::shared_ptr<gko::LinOpFactory> parseConfiguration(
      Context& ctx,
      const GinkgoConfigurationMap& c,
      const gko::config::registry& r,
      const std::shared_ptr<GinkgoExecutor>& e) {
    try {
      auto f = gko::config::parse(GinkgoConfiguration{c}, r,
                                  getGinkgoTypes<parallel>())
                   .on(e->GetExecutor());
      if (f == nullptr) {
        std::ignore = ctx.registerErrorMessage(
            "invalid configuration of Ginkgo (no factory built)");
      }
      return f;
    } catch (std::exception& ex) {
      std::ignore = ctx.registerErrorMessage(
          "invalid configuration of Ginkgo (" + std::string{ex.what()} + ")");
    }
    return {};
  }  // end of parseConfiguration

  /*!
   * \brief check a configuration of Ginkgo key by key, the nested
   * configurations being checked first. Ginkgo terminates the program on an
   * error in a configuration whose other keys have not been read yet (the
   * destructor of its checker of the keys throws during the propagation of
   * the error): each key is thus checked alone, with the keys selecting the
   * type of the object, so that errors are reported by exceptions.
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in] c: configuration
   * \param[in] e: executor
   */
  template <bool parallel>
  [[nodiscard]] static bool checkConfiguration(
      Context& ctx,
      const GinkgoConfiguration& c,
      const std::shared_ptr<GinkgoExecutor>& e) {
    using tag_t = GinkgoConfiguration::tag_t;
    if (c.get_tag() == tag_t::array) {
      for (const auto& v : c.get_array()) {
        if (!checkConfiguration<parallel>(ctx, v, e)) {
          return false;
        }
      }
      return true;
    }
    if (c.get_tag() != tag_t::map) {
      return true;
    }
    const auto& m = c.get_map();
    for (const auto& [k, v] : m) {
      // the stopping criteria are not described by configurations of linear
      // operators: they are checked with their parent
      if (k == "criteria") {
        continue;
      }
      if (!checkConfiguration<parallel>(ctx, v, e)) {
        return false;
      }
    }
    // the maps without type (the simplified form of the stopping criteria)
    // are checked with their parent
    if (!m.contains("type")) {
      return true;
    }
    static const auto type_keys = std::array<const char*, 6>{
        "type",      "value_type",   "index_type", "global_index_type",
        "isai_type", "reverse_apply"};
    auto base = GinkgoConfigurationMap{};
    for (const auto& k : type_keys) {
      if (const auto p = m.find(k); p != m.end()) {
        base.emplace(p->first, p->second);
      }
    }
    for (const auto& [k, v] : m) {
      if (base.contains(k)) {
        continue;
      }
      auto one = base;
      one.emplace(k, v);
      if (parseConfiguration<parallel>(ctx, one, gko::config::registry{}, e) ==
          nullptr) {
        return false;
      }
    }
    return true;
  }  // end of checkConfiguration

  /*!
   * \brief preconditioner of Ginkgo built by a factory of Ginkgo: the base
   * class of the preconditioners of MFEM wrapping Ginkgo converts the
   * matrix given to `SetOperator`
   */
  struct GinkgoFactoryPreconditioner final
      : public mfem::Ginkgo::GinkgoPreconditioner {
    /*!
     * \param[in] e: executor
     * \param[in] f: factory
     */
    GinkgoFactoryPreconditioner(GinkgoExecutor& e,
                                std::shared_ptr<gko::LinOpFactory> f)
        : mfem::Ginkgo::GinkgoPreconditioner(e) {
      this->precond_gen = std::move(f);
    }
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
    /*!
     * \param[in] e: executor
     * \param[in] comm: communicator
     * \param[in] f: factory of a distributed preconditioner
     */
    GinkgoFactoryPreconditioner(GinkgoExecutor& e,
                                MPI_Comm comm,
                                std::shared_ptr<gko::LinOpFactory> f)
        : mfem::Ginkgo::GinkgoPreconditioner(e, comm) {
      this->precond_gen = std::move(f);
    }
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
  };  // end of struct GinkgoFactoryPreconditioner

  /*!
   * \brief preconditioner (or direct solver) of Ginkgo. Given to a solver of
   * Ginkgo, its factory is used, the preconditioner being generated by Ginkgo
   * for each operator. Given to another solver, the preconditioner is
   * generated by `SetOperator`, which requires an assembled matrix.
   */
  struct GinkgoConfiguredPreconditioner final : public mfem::Solver {
    /*!
     * \param[in] e: executor
     * \param[in] p: preconditioner of MFEM wrapping Ginkgo
     * \param[in] l1: if true, the preconditioner is the Schwarz
     * preconditioner of Ginkgo with the L1 smoother
     */
    GinkgoConfiguredPreconditioner(
        std::shared_ptr<GinkgoExecutor> e,
        std::unique_ptr<mfem::Ginkgo::GinkgoPreconditioner> p,
        const bool l1)
        : shared_executor(std::move(e)),
          preconditioner(std::move(p)),
          l1_smoother(l1) {}
    void SetOperator(const mfem::Operator& op) override {
      this->height = op.Height();
      this->width = op.Width();
      this->preconditioner->SetOperator(op);
    }
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
      this->preconditioner->Mult(x, y);
    }
    //! \return the factory of the preconditioner
    [[nodiscard]] std::shared_ptr<gko::LinOpFactory> getFactory() const {
      return this->preconditioner->GetFactory();
    }
    /*!
     * \return true if the preconditioner is the Schwarz preconditioner of
     * Ginkgo with the L1 smoother, which requires a matrix whose diagonal
     * block is sorted
     */
    [[nodiscard]] bool usesL1Smoother() const noexcept {
      return this->l1_smoother;
    }

   private:
    //! \brief executor, kept alive (see `getGinkgoExecutor`)
    const std::shared_ptr<GinkgoExecutor> shared_executor;
    //! \brief preconditioner of MFEM wrapping Ginkgo
    const std::unique_ptr<mfem::Ginkgo::GinkgoPreconditioner> preconditioner;
    //! \brief if true, the Schwarz preconditioner with the L1 smoother is used
    const bool l1_smoother;
  };  // end of struct GinkgoConfiguredPreconditioner

  /*!
   * \brief stopping criteria of a solver of Ginkgo given by the options of
   * mfem-mgis. The default values are the ones of the wrappers of MFEM (no
   * tolerance, ten iterations).
   */
  struct GinkgoStoppingCriteria {
    //! \brief tolerance relative to the norm of the initial residual
    double relative_tolerance = 0;
    //! \brief absolute tolerance
    double absolute_tolerance = 0;
    //! \brief maximum number of iterations
    int maximum_number_of_iterations = 10;
    /*!
     * \brief if true, the residual norm computed by the solver (the
     * preconditioned residual norm of MINRES) is used instead of the norm of
     * the residual, which MINRES does not compute
     */
    bool implicit_residual_norm = false;
  };

  /*!
   * \return the factory of the stopping criteria of Ginkgo: maximum number
   * of iterations, relative and absolute tolerances
   * \param[in] c: stopping criteria
   * \param[in] e: executor
   */
  [[nodiscard]] static std::shared_ptr<gko::stop::CriterionFactory>
  makeGinkgoStoppingCriteria(const GinkgoStoppingCriteria& c,
                             const std::shared_ptr<const gko::Executor>& e) {
    using real_t = mfem::real_t;
    auto iterations = gko::stop::Iteration::build()
                          .with_max_iters(static_cast<gko::size_type>(
                              c.maximum_number_of_iterations))
                          .on(e);
    const auto make = [&e, &c]<typename Criterion>() {
      auto relative =
          Criterion::build()
              .with_baseline(gko::stop::mode::initial_resnorm)
              .with_reduction_factor(static_cast<real_t>(c.relative_tolerance))
              .on(e);
      auto absolute =
          Criterion::build()
              .with_baseline(gko::stop::mode::absolute)
              .with_reduction_factor(static_cast<real_t>(c.absolute_tolerance))
              .on(e);
      return std::make_pair(
          std::shared_ptr<gko::stop::CriterionFactory>(std::move(relative)),
          std::shared_ptr<gko::stop::CriterionFactory>(std::move(absolute)));
    };
    const auto [relative, absolute] =
        c.implicit_residual_norm
            ? make.template
              operator()<gko::stop::ImplicitResidualNorm<real_t>>()
            : make.template operator()<gko::stop::ResidualNorm<real_t>>();
    return gko::stop::Combined::build()
        .with_criteria(std::move(iterations), relative, absolute)
        .on(e);
  }  // end of makeGinkgoStoppingCriteria

  /*!
   * \brief linear solver of Ginkgo described by a configuration of Ginkgo,
   * built on the base class of the iterative solvers of MFEM wrapping Ginkgo
   * (conversion of the matrices and of the vectors, convergence logger).
   *
   * A preconditioner of Ginkgo (`GinkgoConfiguredPreconditioner`) is given
   * by its factory, Ginkgo generating the preconditioner for each operator.
   * Another preconditioner is wrapped by MFEM
   * (`mfem::Ginkgo::MFEMPreconditioner`), the wrapper reading the size of the
   * preconditioner: the factory of the solver is then built again for each
   * operator.
   *
   * The stopping criteria given by the options of mfem-mgis are given to
   * Ginkgo through its registry, so that the relative tolerance can be
   * changed (see `setRelativeTolerance`).
   */
  template <bool parallel>
  struct GinkgoLinearSolver final : public mfem::Ginkgo::GinkgoIterativeSolver,
                                    public LinearSolverConvergenceReport {
    /*!
     * \param[in] e: executor
     * \param[in] c: configuration, without stopping criteria if they are
     * given by `sc`
     * \param[in] sc: stopping criteria given by the options of mfem-mgis
     * \param[in] k: key of the preconditioner in the configuration
     * (`preconditioner`, or `solver` for the iterative refinement)
     */
    GinkgoLinearSolver(std::shared_ptr<GinkgoExecutor> e,
                       GinkgoConfigurationMap c,
                       std::optional<GinkgoStoppingCriteria> sc,
                       std::string k)
      requires(!parallel)
        : mfem::Ginkgo::GinkgoIterativeSolver(*e, false),
          shared_executor(std::move(e)),
          configuration(std::move(c)),
          stopping_criteria(std::move(sc)),
          preconditioner_key(std::move(k)) {}
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
    /*!
     * \param[in] e: executor
     * \param[in] comm: communicator
     * \param[in] c: configuration, without stopping criteria if they are
     * given by `sc`
     * \param[in] sc: stopping criteria given by the options of mfem-mgis
     * \param[in] k: key of the preconditioner in the configuration
     * (`preconditioner`, or `solver` for the iterative refinement)
     */
    GinkgoLinearSolver(std::shared_ptr<GinkgoExecutor> e,
                       MPI_Comm comm,
                       GinkgoConfigurationMap c,
                       std::optional<GinkgoStoppingCriteria> sc,
                       std::string k)
      requires(parallel)
        : mfem::Ginkgo::GinkgoIterativeSolver(*e, comm, false),
          shared_executor(std::move(e)),
          communicator(comm),
          configuration(std::move(c)),
          stopping_criteria(std::move(sc)),
          preconditioner_key(std::move(k)) {}
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
    /*!
     * \brief set the preconditioner, if any, and build the factory of the
     * solver. With a preconditioner which is not one of Ginkgo, the
     * configuration is checked with a placeholder of the preconditioner.
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] p: preconditioner (not owned), or a null pointer
     */
    [[nodiscard]] bool setup(Context& ctx, mfem::Solver* const p) {
      if (const auto* const gp =
              dynamic_cast<const GinkgoConfiguredPreconditioner*>(p);
          gp != nullptr) {
        this->preconditioner_factory = gp->getFactory();
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
        this->needs_sorted_diagonal_mat = gp->usesL1Smoother();
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
      } else if (p != nullptr) {
        this->shell_preconditioner = p;
        this->needs_wrapped_vecs = true;
        this->sub_op_needs_wrapped_vecs = true;
      }
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
      // a Schwarz preconditioner described by the configuration of the
      // solver
      if (usesL1Smoother(GinkgoConfiguration{this->configuration})) {
        this->needs_sorted_diagonal_mat = true;
      }
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
      return this->buildFactory(ctx);
    }
    /*!
     * \brief only accept assembled operators (sparse matrices of MFEM and
     * matrices of hypre), for the solvers which transpose the operator or
     * use vectors of several columns
     */
    void requireAssembledOperators() noexcept {
      this->assembled_operators_only = true;
    }
    void SetOperator(const mfem::Operator& op) override {
      if ((this->assembled_operators_only) &&
          (dynamic_cast<const mfem::SparseMatrix*>(&op) == nullptr)
#ifdef MFEM_USE_MPI
          && (dynamic_cast<const mfem::HypreParMatrix*>(&op) == nullptr)
#endif /* MFEM_USE_MPI */
      ) {
        raise(
            "GinkgoLinearSolver::SetOperator: this solver of Ginkgo requires "
            "an assembled operator");
      }
      if (this->shell_preconditioner != nullptr) {
        // unlike the Krylov solvers of MFEM, the solvers of Ginkgo do not
        // give the operator to their preconditioner
        this->shell_preconditioner->SetOperator(op);
        // the wrapper reads the size of the preconditioner: it is built
        // again for each operator, by all the processes in parallel
        this->wrapPreconditioner();
        auto ctx = Context{};
        if (!this->buildFactory(ctx)) {
          raise("GinkgoLinearSolver::SetOperator: " + ctx.getErrorMessage());
        }
      }
      mfem::Ginkgo::GinkgoIterativeSolver::SetOperator(op);
    }
    [[nodiscard]] bool hasConverged() const noexcept override {
      return this->GetConverged() == 1;
    }
    [[nodiscard]] std::optional<int> getNumberOfIterations()
        const noexcept override {
      return this->GetNumIterations();
    }
    [[nodiscard]] std::optional<real> getRelativeTolerance()
        const noexcept override {
      if (!this->stopping_criteria.has_value()) {
        return {};
      }
      return static_cast<real>(this->stopping_criteria->relative_tolerance);
    }
    [[nodiscard]] bool setRelativeTolerance(const real tol) noexcept override {
      if (!this->stopping_criteria.has_value()) {
        return false;
      }
      try {
        const auto previous = this->stopping_criteria->relative_tolerance;
        this->stopping_criteria->relative_tolerance = tol;
        auto ctx = Context{};
        if (!this->buildFactory(ctx)) {
          this->stopping_criteria->relative_tolerance = previous;
          return false;
        }
        // the solver generated for the current operator, if any
        if (auto* const s =
                dynamic_cast<gko::solver::IterativeBase*>(this->solver.get());
            s != nullptr) {
          s->set_stop_criterion_factory(makeGinkgoStoppingCriteria(
              *(this->stopping_criteria), this->executor));
        }
        return true;
      } catch (...) {
      }
      return false;
    }

   private:
    //! \brief wrap the preconditioner which is not one of Ginkgo
    void wrapPreconditioner() {
      using MFEMPreconditioner = mfem::Ginkgo::MFEMPreconditioner;
      auto& e = *(this->shared_executor);
      if constexpr (parallel) {
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
        this->wrapper = std::make_unique<MFEMPreconditioner>(
            e, *(this->shell_preconditioner), this->communicator);
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
      } else {
        this->wrapper = std::make_unique<MFEMPreconditioner>(
            e, *(this->shell_preconditioner));
      }
    }
    /*!
     * \brief build the factory of the solver from its configuration, its
     * stopping criteria and its preconditioner, given through the registry of
     * Ginkgo. A preconditioner which is not one of Ginkgo is replaced by a
     * placeholder until it is wrapped.
     * \return true on success
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] bool buildFactory(Context& ctx) {
      auto r = gko::config::registry{};
      auto c = this->configuration;
      if (this->stopping_criteria.has_value()) {
        r.emplace(ginkgo_criteria_name,
                  makeGinkgoStoppingCriteria(*(this->stopping_criteria),
                                             this->executor));
        if (!c.emplace("criteria", GinkgoConfiguration{ginkgo_criteria_name})
                 .second) {
          return ctx.registerErrorMessage(
              "the stopping criteria are given both by the options of "
              "mfem-mgis and by the configuration of Ginkgo");
        }
      }
      const auto generated_key = "generated_" + this->preconditioner_key;
      if (((this->preconditioner_factory != nullptr) ||
           (this->shell_preconditioner != nullptr)) &&
          ((c.contains(this->preconditioner_key)) ||
           (c.contains(generated_key)))) {
        return ctx.registerErrorMessage(
            "the preconditioner is given both by the option 'Preconditioner' "
            "and by the configuration of Ginkgo");
      }
      if (this->preconditioner_factory != nullptr) {
        r.emplace(ginkgo_preconditioner_name, this->preconditioner_factory);
        c.emplace(this->preconditioner_key,
                  GinkgoConfiguration{ginkgo_preconditioner_name});
      } else if (this->shell_preconditioner != nullptr) {
        if (this->wrapper != nullptr) {
          r.emplace(ginkgo_preconditioner_name,
                    this->wrapper->GetGeneratedPreconditioner());
        } else {
          r.emplace(ginkgo_preconditioner_name,
                    std::shared_ptr<gko::LinOp>(
                        gko::matrix::Identity<mfem::real_t>::create(
                            this->executor, 1)));
        }
        c.emplace(generated_key,
                  GinkgoConfiguration{ginkgo_preconditioner_name});
      }
      auto f = parseConfiguration<parallel>(ctx, c, r, this->shared_executor);
      if (f == nullptr) {
        return false;
      }
      this->solver_gen = std::move(f);
      return true;
    }
    //! \brief executor, kept alive (see `getGinkgoExecutor`)
    const std::shared_ptr<GinkgoExecutor> shared_executor;
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
    //! \brief communicator
    MPI_Comm communicator = MPI_COMM_NULL;
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
    //! \brief configuration of the solver, without its preconditioner
    const GinkgoConfigurationMap configuration;
    //! \brief stopping criteria given by the options of mfem-mgis, if any
    std::optional<GinkgoStoppingCriteria> stopping_criteria;
    //! \brief key of the preconditioner in the configuration
    const std::string preconditioner_key;
    //! \brief factory of the preconditioner of Ginkgo, if any
    std::shared_ptr<gko::LinOpFactory> preconditioner_factory;
    //! \brief preconditioner which is not one of Ginkgo, if any (not owned)
    mfem::Solver* shell_preconditioner = nullptr;
    //! \brief wrapper of this preconditioner
    std::unique_ptr<mfem::Ginkgo::MFEMPreconditioner> wrapper;
    //! \brief if true, only assembled operators are accepted
    bool assembled_operators_only = false;
  };  // end of struct GinkgoLinearSolver

  /*!
   * \brief description of a linear solver of Ginkgo registered in the
   * factory of linear solvers
   */
  struct GinkgoSolverEntry {
    //! \brief name of the entry
    std::string name;
    //! \brief type of the solver in the configurations of Ginkgo
    std::string type;
    //! \brief description of the entry
    std::string description;
    /*!
     * \brief if true, the stopping criteria use the residual norm computed
     * by the solver (see `GinkgoStoppingCriteria`)
     */
    bool implicit_residual_norm = false;
    //! \brief if true, the `KDim` option sets the dimension of the Krylov space
    bool krylov_dimension = false;
    //! \brief if true, the solver is available in parallel
    bool distributed = true;
    /*!
     * \brief if true, the solver only accepts the preconditioners of Ginkgo
     * (the solvers which apply the transposed preconditioner, or whose
     * workspace is not supported by the vectors wrapped by MFEM) and
     * assembled operators
     */
    bool ginkgo_preconditioners_only = false;
    /*!
     * \brief if true, the solver applies the transposed preconditioner, which
     * restricts the preconditioners of Ginkgo it accepts (see
     * `getTransposableGinkgoPreconditioners`)
     */
    bool transposed_preconditioner = false;
  };

  /*!
   * \brief read the stopping criteria given by the options of mfem-mgis
   * (relative and absolute tolerances, maximum number of iterations), unless
   * the configuration of Ginkgo defines them
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] sc: stopping criteria, empty if the configuration of Ginkgo
   * defines them
   * \param[in] c: configuration
   * \param[in] params: parameters of the solver
   * \param[in] implicit: if true, the residual norm computed by the solver
   * is used
   */
  [[nodiscard]] static bool getStoppingCriteria(
      Context& ctx,
      std::optional<GinkgoStoppingCriteria>& sc,
      const GinkgoConfigurationMap& c,
      const Parameters& params,
      const bool implicit) {
    using Problem = AbstractNonLinearEvolutionProblem;
    const auto given =
        contains(params, Problem::SolverRelativeTolerance) ||
        contains(params, Problem::SolverAbsoluteTolerance) ||
        contains(params, Problem::SolverMaximumNumberOfIterations);
    if (c.contains("criteria")) {
      if (given) {
        return ctx.registerErrorMessage(
            "the stopping criteria are given both by the options of "
            "mfem-mgis and by the configuration of Ginkgo");
      }
      sc.reset();
      return true;
    }
    const auto ortol =
        get_if<double>(ctx, params, Problem::SolverRelativeTolerance, 0.);
    const auto oatol =
        get_if<double>(ctx, params, Problem::SolverAbsoluteTolerance, 0.);
    const auto omaxit =
        get_if<int>(ctx, params, Problem::SolverMaximumNumberOfIterations, 10);
    if (isInvalid(ortol) || isInvalid(oatol) || isInvalid(omaxit)) {
      return false;
    }
    if (*omaxit < 0) {
      return ctx.registerErrorMessage(
          "invalid maximum number of iterations (must not be negative)");
    }
    sc = GinkgoStoppingCriteria{
        .relative_tolerance = boundRelativeTolerance(*ortol),
        .absolute_tolerance = *oatol,
        .maximum_number_of_iterations = *omaxit,
        .implicit_residual_norm = implicit};
    return true;
  }  // end of getStoppingCriteria

  /*!
   * \return the preconditioners of Ginkgo which can be transposed, as
   * required by the solvers applying the transposed preconditioner
   */
  [[nodiscard]] static std::vector<std::string>
  getTransposableGinkgoPreconditioners() {
    return {"GinkgoJacobi", "GinkgoILU", "GinkgoILUISAI", "GinkgoIC",
            "GinkgoICISAI"};
  }  // end of getTransposableGinkgoPreconditioners

  /*!
   * \brief check that the preconditioner given to a solver applying the
   * transposed preconditioner can be transposed
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in] n: name of the solver
   * \param[in] params: parameters of the solver
   */
  [[nodiscard]] static bool checkTransposablePreconditioner(
      Context& ctx, const std::string& n, const Parameters& params) {
    if (!contains(params, "Preconditioner")) {
      return true;
    }
    const auto op = get<Parameters>(ctx, params, "Preconditioner");
    if (isInvalid(op)) {
      return false;
    }
    const auto oname = get<std::string>(ctx, *op, "Name");
    if (isInvalid(oname)) {
      return false;
    }
    const auto names = getTransposableGinkgoPreconditioners();
    if ((*oname == "None") ||
        (std::find(names.begin(), names.end(), *oname) != names.end())) {
      return true;
    }
    auto msg = "the linear solver '" + n +
               "' requires a preconditioner which can be transposed:";
    for (const auto& name : names) {
      msg += " " + name;
    }
    return ctx.registerErrorMessage(msg);
  }  // end of checkTransposablePreconditioner

  /*!
   * \return the generator of a linear solver of Ginkgo
   * \param[in] entry: entry
   * \param[in] d: description of the entry
   */
  template <bool parallel>
  [[nodiscard]] static typename LinearSolverFactory<parallel>::Generator
  buildGinkgoSolverGenerator(const GinkgoSolverEntry entry,
                             const LinearSolverDescription d) {
    return [entry, d](Context& ctx, FiniteElementSpace<parallel>& fespace,
                      const Parameters& params) -> LinearSolverHandler {
      using Problem = AbstractNonLinearEvolutionProblem;
      if (!checkParameters(ctx, params, d.options)) {
        return {};
      }
      auto c = GinkgoConfigurationMap{};
      if (!entry.type.empty()) {
        c.emplace("type", GinkgoConfiguration{entry.type});
      } else if (!contains(params, "Configuration")) {
        std::ignore = ctx.registerErrorMessage(
            "no configuration of Ginkgo given to '" + entry.name + "'");
        return {};
      }
      if ((entry.krylov_dimension) && (contains(params, "KDim"))) {
        const auto odim = get<int>(ctx, params, "KDim");
        if (isInvalid(odim)) {
          return {};
        }
        if (*odim < 1) {
          std::ignore = ctx.registerErrorMessage(
              "invalid Krylov dimension (must be strictly positive)");
          return {};
        }
        c.emplace("krylov_dim",
                  GinkgoConfiguration{static_cast<std::int64_t>(*odim)});
      }
      if (!addConfiguration(ctx, c, params)) {
        return {};
      }
      if (!c.contains("type")) {
        std::ignore =
            ctx.registerErrorMessage("the configuration of Ginkgo given to '" +
                                     entry.name + "' has no type");
        return {};
      }
      // the convergence of the other types (direct solvers, preconditioners)
      // is not reported by Ginkgo
      const auto& type = c.at("type");
      if ((type.get_tag() != GinkgoConfiguration::tag_t::string) ||
          (!isGinkgoIterativeSolver(type.get_string()))) {
        std::ignore = ctx.registerErrorMessage(
            "the configuration of Ginkgo given to '" + entry.name +
            "' does not describe an iterative solver of Ginkgo");
        return {};
      }
      auto e = getGinkgoExecutor();
      if (!checkConfiguration<parallel>(ctx, GinkgoConfiguration{c}, e)) {
        return {};
      }
      auto sc = std::optional<GinkgoStoppingCriteria>{};
      if (!getStoppingCriteria(ctx, sc, c, params,
                               entry.implicit_residual_norm)) {
        return {};
      }
      // the iterative refinement takes its preconditioner as inner solver
      const auto key =
          ((type.get_tag() == GinkgoConfiguration::tag_t::string) &&
           (type.get_string() == "solver::Ir"))
              ? "solver"
              : "preconditioner";
      auto s = std::unique_ptr<GinkgoLinearSolver<parallel>>{};
      if constexpr (parallel) {
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
        s = std::make_unique<GinkgoLinearSolver<true>>(
            e, fespace.GetComm(), std::move(c), std::move(sc), key);
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
      } else {
        s = std::make_unique<GinkgoLinearSolver<false>>(e, std::move(c),
                                                        std::move(sc), key);
      }
      if (entry.ginkgo_preconditioners_only) {
        s->requireAssembledOperators();
      }
      if (contains(params, Problem::SolverVerbosityLevel)) {
        const auto ov = get<int>(ctx, params, Problem::SolverVerbosityLevel);
        if (isInvalid(ov)) {
          return {};
        }
        s->SetPrintLevel(*ov);
      }
      auto oprec =
          makePreconditioner<parallel>(ctx, fespace, entry.name, d, params);
      if (isInvalid(oprec)) {
        return {};
      }
      if ((entry.transposed_preconditioner) &&
          (!checkTransposablePreconditioner(ctx, entry.name, params))) {
        return {};
      }
      if (!s->setup(ctx, oprec->get())) {
        return {};
      }
      return LinearSolverHandler{std::move(s), std::move(*oprec)};
    };
  }  // end of buildGinkgoSolverGenerator

  /*!
   * \brief description of a preconditioner (or of a direct solver) of
   * Ginkgo registered in the factory of preconditioners
   */
  struct GinkgoPreconditionerEntry {
    //! \brief name of the entry
    std::string name;
    //! \brief configuration (empty for the generic entry)
    GinkgoConfigurationMap configuration;
    //! \brief description of the entry
    std::string description;
    /*!
     * \brief if true, the preconditioner is local: it is the local solver of
     * the Schwarz preconditioner of Ginkgo in parallel
     */
    bool local = true;
    /*!
     * \brief keys added to the configuration when the `Configuration` option
     * does not define them
     */
    GinkgoConfigurationMap defaults = {};
  };

  /*!
   * \return the preconditioner of Ginkgo described by an entry and its
   * options, or a null pointer on error
   * \param[in, out] ctx: execution context
   * \param[in] fespace: finite element space
   * \param[in] entry: entry
   * \param[in] opts: options (`Configuration`, and `L1Smoother` for the
   * local preconditioners in parallel)
   */
  template <bool parallel>
  [[nodiscard]] static std::unique_ptr<GinkgoConfiguredPreconditioner>
  makeGinkgoPreconditioner(
      Context& ctx,
      [[maybe_unused]] FiniteElementSpace<parallel>& fespace,
      const GinkgoPreconditionerEntry& entry,
      const Parameters& opts) {
    const auto schwarz = parallel && entry.local;
    auto options = std::vector<std::string>{"Configuration"};
    if (schwarz) {
      options.push_back("L1Smoother");
    }
    if (!checkParameters(ctx, opts, options)) {
      return {};
    }
    auto c = entry.configuration;
    if (!addConfiguration(ctx, c, opts)) {
      return {};
    }
    for (const auto& [k, v] : entry.defaults) {
      c.emplace(k, v);
    }
    if (!c.contains("type")) {
      std::ignore =
          ctx.registerErrorMessage("the configuration of Ginkgo given to '" +
                                   entry.name + "' has no type");
      return {};
    }
    // the L1 smoother requires the sorting of the diagonal blocks of the
    // matrix made by the Schwarz preconditioner of MFEM (`L1Smoother`
    // option of the local preconditioners) or by the solvers of Ginkgo
    if (usesL1Smoother(GinkgoConfiguration{c})) {
      std::ignore = ctx.registerErrorMessage(
          "the L1 smoother of the Schwarz preconditioner of Ginkgo can't be "
          "given by the configuration of '" +
          entry.name +
          "' (see the option 'L1Smoother' of the local preconditioners)");
      return {};
    }
    auto e = getGinkgoExecutor();
    // in parallel, the local preconditioners are described with the types of
    // the local matrices
    const auto checked =
        schwarz ? checkConfiguration<false>(ctx, GinkgoConfiguration{c}, e)
                : checkConfiguration<parallel>(ctx, GinkgoConfiguration{c}, e);
    if (!checked) {
      return {};
    }
    if constexpr (parallel) {
#ifdef MFEM_MGIS_GINKGO_DISTRIBUTED
      if (schwarz) {
        const auto ol1 = get_if<bool>(ctx, opts, "L1Smoother", false);
        if (isInvalid(ol1)) {
          return {};
        }
        // the local solver of the Schwarz preconditioner of MFEM, which
        // sorts the diagonal block of the matrix required by the L1 smoother
        auto f = parseConfiguration<false>(ctx, c, gko::config::registry{}, e);
        if (f == nullptr) {
          return {};
        }
        auto local = GinkgoFactoryPreconditioner(*e, std::move(f));
        auto p = std::make_unique<mfem::Ginkgo::SchwarzPreconditioner>(
            *e, fespace.GetComm(), local, *ol1);
        return std::make_unique<GinkgoConfiguredPreconditioner>(e, std::move(p),
                                                                *ol1);
      }
      auto f = parseConfiguration<true>(ctx, c, gko::config::registry{}, e);
      if (f == nullptr) {
        return {};
      }
      auto p = std::make_unique<GinkgoFactoryPreconditioner>(
          *e, fespace.GetComm(), std::move(f));
      return std::make_unique<GinkgoConfiguredPreconditioner>(e, std::move(p),
                                                              false);
#else  /* MFEM_MGIS_GINKGO_DISTRIBUTED */
      return {};
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
    } else {
      auto f = parseConfiguration<false>(ctx, c, gko::config::registry{}, e);
      if (f == nullptr) {
        return {};
      }
      auto p = std::make_unique<GinkgoFactoryPreconditioner>(*e, std::move(f));
      return std::make_unique<GinkgoConfiguredPreconditioner>(e, std::move(p),
                                                              false);
    }
  }  // end of makeGinkgoPreconditioner

  /*!
   * \return a configuration of Ginkgo of the given type
   * \param[in] type: type of the object of Ginkgo
   */
  [[nodiscard]] static GinkgoConfigurationMap makeConfiguration(
      const std::string& type) {
    return GinkgoConfigurationMap{{"type", GinkgoConfiguration{type}}};
  }  // end of makeConfiguration

  /*!
   * \return the configuration of the incomplete sparse approximate inverse
   * \param[in] isai_type: type of the inverse (`lower`, `upper`)
   */
  [[nodiscard]] static GinkgoConfiguration makeISAIConfiguration(
      const std::string& isai_type) {
    return GinkgoConfiguration{GinkgoConfigurationMap{
        {"type", GinkgoConfiguration{"preconditioner::Isai"}},
        {"isai_type", GinkgoConfiguration{isai_type}}}};
  }  // end of makeISAIConfiguration

  /*!
   * \return the default keys of the configuration of the sparse direct
   * solver of Ginkgo: the LU factorization
   */
  [[nodiscard]] static GinkgoConfigurationMap getDirectSolverDefaults() {
    return GinkgoConfigurationMap{
        {"factorization",
         GinkgoConfiguration{makeConfiguration("factorization::Lu")}}};
  }  // end of getDirectSolverDefaults

  //! \return the preconditioners of Ginkgo
  [[nodiscard]] static std::vector<GinkgoPreconditionerEntry>
  getGinkgoPreconditionerEntries() {
    auto ilu_isai = makeConfiguration("preconditioner::Ilu");
    ilu_isai.emplace("l_solver", makeISAIConfiguration("lower"));
    ilu_isai.emplace("u_solver", makeISAIConfiguration("upper"));
    auto ic_isai = makeConfiguration("preconditioner::Ic");
    ic_isai.emplace("l_solver", makeISAIConfiguration("lower"));
    auto multigrid = GinkgoConfigurationMap{};
    multigrid.emplace(
        "mg_level",
        GinkgoConfiguration{GinkgoConfiguration::array_type{
            GinkgoConfiguration{makeConfiguration("multigrid::Pgm")}}});
    // one cycle: the simplified form of the stopping criteria (a map) is not
    // supported by the configuration of the multigrid of Ginkgo 1.11
    multigrid.emplace(
        "criteria",
        GinkgoConfiguration{GinkgoConfiguration::array_type{
            GinkgoConfiguration{GinkgoConfigurationMap{
                {"type", GinkgoConfiguration{"Iteration"}},
                {"max_iters", GinkgoConfiguration{std::int64_t{1}}}}}}});
    return {
        {.name = "GinkgoJacobi",
         .configuration = makeConfiguration("preconditioner::Jacobi"),
         .description = "block Jacobi preconditioner of Ginkgo"},
        {.name = "GinkgoILU",
         .configuration = makeConfiguration("preconditioner::Ilu"),
         .description = "incomplete LU factorization of Ginkgo"},
        {.name = "GinkgoILUISAI",
         .configuration = std::move(ilu_isai),
         .description = "incomplete LU factorization of Ginkgo, whose "
                        "triangular factors are applied by incomplete sparse "
                        "approximate inverses"},
        {.name = "GinkgoIC",
         .configuration = makeConfiguration("preconditioner::Ic"),
         .description = "incomplete Cholesky factorization of Ginkgo"},
        {.name = "GinkgoICISAI",
         .configuration = std::move(ic_isai),
         .description = "incomplete Cholesky factorization of Ginkgo, whose "
                        "triangular factors are applied by incomplete sparse "
                        "approximate inverses"},
        {.name = "GinkgoSOR",
         .configuration = makeConfiguration("preconditioner::Sor"),
         .description = "successive over-relaxation of Ginkgo"},
        {.name = "GinkgoGaussSeidel",
         .configuration = makeConfiguration("preconditioner::GaussSeidel"),
         .description = "Gauss-Seidel preconditioner of Ginkgo"},
        {.name = "GinkgoDirect",
         .configuration = makeConfiguration("solver::Direct"),
         .description = "sparse LU factorization of Ginkgo",
         .defaults = getDirectSolverDefaults()},
        {.name = "GinkgoMultigrid",
         .configuration = makeConfiguration("solver::Multigrid"),
         .description = "one cycle of the algebraic multigrid of Ginkgo "
                        "(parallel graph match aggregation)",
         .local = false,
         .defaults = std::move(multigrid)},
        {.name = "GinkgoPreconditioner",
         .configuration = {},
         .description = "preconditioner of Ginkgo described by the "
                        "configuration given by the 'Configuration' "
                        "parameter",
         .local = false}};
  }  // end of getGinkgoPreconditionerEntries

  //! \return the linear solvers of Ginkgo
  [[nodiscard]] static std::vector<GinkgoSolverEntry> getGinkgoSolverEntries() {
    return {{.name = "GinkgoCG",
             .type = "solver::Cg",
             .description = "conjugate gradient of Ginkgo"},
            {.name = "GinkgoFCG",
             .type = "solver::Fcg",
             .description = "flexible conjugate gradient of Ginkgo"},
            // its workspace has vectors of several columns, which the
            // vectors wrapped by MFEM (in parallel, or with a preconditioner
            // of another library) do not support
            {.name = "GinkgoPipeCG",
             .type = "solver::PipeCg",
             .description = "pipelined conjugate gradient of Ginkgo",
             .distributed = false,
             .ginkgo_preconditioners_only = true},
            {.name = "GinkgoBiCGSTAB",
             .type = "solver::Bicgstab",
             .description = "BiCGSTAB of Ginkgo"},
            {.name = "GinkgoCGS",
             .type = "solver::Cgs",
             .description = "conjugate gradient squared of Ginkgo"},
            {.name = "GinkgoGMRES",
             .type = "solver::Gmres",
             .description = "GMRES of Ginkgo",
             .krylov_dimension = true},
            {.name = "GinkgoGCR",
             .type = "solver::Gcr",
             .description = "generalized conjugate residual of Ginkgo",
             .krylov_dimension = true},
            // the distributed version fails with the vectors wrapped by MFEM
            {.name = "GinkgoIDR",
             .type = "solver::Idr",
             .description = "induced dimension reduction of Ginkgo",
             .distributed = false},
            // MINRES does not compute the residual, whose norm would be
            // computed by an additional product at each iteration
            {.name = "GinkgoMINRES",
             .type = "solver::Minres",
             .description = "MINRES of Ginkgo",
             .implicit_residual_norm = true},
            {.name = "GinkgoIR",
             .type = "solver::Ir",
             .description = "iterative refinement of Ginkgo (stationary "
                            "iterations with the preconditioner as inner "
                            "solver)"},
            {.name = "GinkgoBiCG",
             .type = "solver::Bicg",
             .description = "BiCG of Ginkgo",
             .distributed = false,
             .ginkgo_preconditioners_only = true,
             .transposed_preconditioner = true},
            {.name = "GinkgoCBGMRES",
             .type = "solver::CbGmres",
             .description = "GMRES of Ginkgo storing its Krylov basis in a "
                            "reduced precision",
             .krylov_dimension = true,
             .distributed = false},
            {.name = "GinkgoSolver",
             .type = "",
             .description = "solver of Ginkgo described by the configuration "
                            "given by the 'Configuration' parameter"}};
  }  // end of getGinkgoSolverEntries

  template <bool parallel>
  void declareGinkgoLinearSolvers(LinearSolverFactory<parallel>& f) {
#ifndef MFEM_MGIS_GINKGO_DISTRIBUTED
    if constexpr (parallel) {
      return;
    }
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
    for (const auto& entry : getGinkgoSolverEntries()) {
      if ((parallel) && (!entry.distributed)) {
        continue;
      }
      auto options = getIterativeSolverParametersList();
      options.insert(options.end(), {"Preconditioner", "Configuration"});
      if (entry.krylov_dimension) {
        options.push_back("KDim");
      }
      const auto d = LinearSolverDescription{
          .library = "Ginkgo",
          .kind = LinearSolverKind::ITERATIVE,
          .placement = LinearSolverPlacement::MFEM,
          .preconditioners = entry.ginkgo_preconditioners_only
                                 ? AcceptedPreconditioners::SAME_LIBRARY
                                 : AcceptedPreconditioners::ANY,
          .options = std::move(options),
          .description = entry.description};
      declare(f, entry.name, d, buildGinkgoSolverGenerator<parallel>(entry, d));
    }
    if constexpr (!parallel) {
      // the direct solver of Ginkgo, which is not distributed
      const auto d = LinearSolverDescription{
          .library = "Ginkgo",
          .kind = LinearSolverKind::DIRECT,
          .placement = LinearSolverPlacement::MFEM,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = {"Configuration"},
          .description = "sparse LU factorization of Ginkgo"};
      declare(
          f, "GinkgoDirect", d,
          [](Context& ctx, FiniteElementSpace<false>& fespace,
             const Parameters& params) -> LinearSolverHandler {
            const auto entry = GinkgoPreconditionerEntry{
                .name = "GinkgoDirect",
                .configuration = makeConfiguration("solver::Direct"),
                .description = {},
                .defaults = getDirectSolverDefaults()};
            auto s =
                makeGinkgoPreconditioner<false>(ctx, fespace, entry, params);
            if (s == nullptr) {
              return {};
            }
            return LinearSolverHandler{
                std::move(s), std::unique_ptr<LinearSolverPreconditioner>{}};
          });
    }
  }  // end of declareGinkgoLinearSolvers

  template <bool parallel>
  void declareGinkgoPreconditioners(PreconditionerFactory<parallel>& f) {
#ifndef MFEM_MGIS_GINKGO_DISTRIBUTED
    if constexpr (parallel) {
      return;
    }
#endif /* MFEM_MGIS_GINKGO_DISTRIBUTED */
    for (auto& entry : getGinkgoPreconditionerEntries()) {
      auto options = std::vector<std::string>{"Configuration"};
      if ((parallel) && (entry.local)) {
        options.push_back("L1Smoother");
      }
      const auto d = LinearSolverDescription{
          .library = "Ginkgo",
          .kind = LinearSolverKind::PRECONDITIONER,
          .placement = LinearSolverPlacement::MFEM,
          .preconditioners = AcceptedPreconditioners::NONE,
          .options = std::move(options),
          .description = entry.description};
      declare(f, entry.name, d,
              [entry](Context& ctx, FiniteElementSpace<parallel>& fespace,
                      const Parameters& opts)
                  -> std::unique_ptr<LinearSolverPreconditioner> {
                return makeGinkgoPreconditioner<parallel>(ctx, fespace, entry,
                                                          opts);
              });
    }
  }  // end of declareGinkgoPreconditioners

#ifdef MFEM_USE_MPI
  template void declareGinkgoLinearSolvers<true>(LinearSolverFactory<true>&);
  template void declareGinkgoPreconditioners<true>(
      PreconditionerFactory<true>&);
#endif /* MFEM_USE_MPI */
  template void declareGinkgoLinearSolvers<false>(LinearSolverFactory<false>&);
  template void declareGinkgoPreconditioners<false>(
      PreconditionerFactory<false>&);

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_MGIS_HAVE_GINKGO */
