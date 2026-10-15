/*!
 * \file   src/LinearSolverFactory.cxx
 * \brief  This file implements the `LinearSolverFactory` class. The linear
 * solvers are registered by the modules of each library (see
 * `LinearSolverModules.hxx`).
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <string>
#include <utility>
#include <exception>
#include "mfem/config/config.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "ExecutionSettings.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis {

  namespace {

    /*!
     * \brief implementation of `LinearSolverFactory<parallel>::add`
     * \param[in, out] ctx: execution context
     * \param[in, out] entries: registered linear solvers
     * \param[in] n: name of the linear solver
     * \param[in] e: description and generator of the linear solver
     */
    template <typename Entries, typename Entry>
    [[nodiscard]] bool addLinearSolver(Context& ctx,
                                       Entries& entries,
                                       const std::string_view n,
                                       Entry e) noexcept {
      if (entries.find(n) != entries.end()) {
        return ctx.registerErrorMessage(
            "LinearSolverFactory::add: a linear solver called '" +
            std::string{n} + "' has already been declared");
      }
      entries.insert({std::string{n}, std::move(e)});
      return true;
    }  // end of addLinearSolver

    /*!
     * \brief implementation of `LinearSolverFactory<parallel>::generate`
     * \param[in, out] ctx: execution context
     * \param[in] entries: registered linear solvers
     * \param[in] n: name of the linear solver
     * \param[in] fespace: finite element space
     * \param[in] params: parameters of the linear solver
     */
    template <bool parallel, typename Entries>
    LinearSolverHandler generateLinearSolver(
        Context& ctx,
        const Entries& entries,
        const std::string_view n,
        FiniteElementSpace<parallel>& fespace,
        const Parameters& params) {
      const auto pe = entries.find(n);
      if (pe == entries.end()) {
#ifdef MFEM_USE_MPI
        if constexpr (!parallel) {
          if (LinearSolverFactory<true>::getFactory().getDescription(n) !=
              nullptr) {
            return ctx.registerErrorMessage(
                "LinearSolverFactory<false>::generate: the linear solver '" +
                std::string{n} + "' is only available in parallel");
          }
        } else {
          if (LinearSolverFactory<false>::getFactory().getDescription(n) !=
              nullptr) {
            return ctx.registerErrorMessage(
                "LinearSolverFactory<true>::generate: the linear solver '" +
                std::string{n} + "' is only available in sequential");
          }
        }
#endif /* MFEM_USE_MPI */
        const auto hint = internals::getUnavailableLibraryHint(n);
        if (!hint.empty()) {
          return ctx.registerErrorMessage(
              std::string{"LinearSolverFactory<"} +
              (parallel ? "true" : "false") +
              ">::generate: no linear solver called '" + std::string{n} +
              "' declared: " + hint);
        }
        auto msg = std::string{"LinearSolverFactory<"} +
                   (parallel ? "true" : "false") +
                   ">::generate: no linear solver called '" + std::string{n} +
                   "' declared. Here is the list of available linear "
                   "solvers:";
        for (const auto& e : entries) {
          msg += " " + e.first;
        }
        return ctx.registerErrorMessage(msg);
      }
      const auto& e = pe->second;
      if ((e.description.placement.has_value()) &&
          (!internals::checkLinearSolverExecutionLocation(
              ctx, n, *(e.description.placement)))) {
        return {};
      }
      // a preconditioner given to a linear solver which accepts none is
      // refused, `None` (no preconditioner) being left to the solver
      if ((e.description.preconditioners == AcceptedPreconditioners::NONE) &&
          (contains(params, "Preconditioner"))) {
        const auto opr = get<Parameters>(ctx, params, "Preconditioner");
        if (isInvalid(opr)) {
          return {};
        }
        const auto oname = get<std::string>(ctx, *opr, "Name");
        if (isInvalid(oname)) {
          return {};
        }
        if (*oname != "None") {
          return ctx.registerErrorMessage("the linear solver '" +
                                          std::string{n} +
                                          "' accepts no preconditioner");
        }
      }
      try {
        return e.generator(ctx, fespace, params);
      } catch (std::exception& ex) {
        return ctx.registerErrorMessage(
            std::string{"LinearSolverFactory<"} +
            (parallel ? "true" : "false") +
            ">::generate: error while generating linear solver '" +
            std::string{n} + "'\n" + ex.what());
      }
    }  // end of generateLinearSolver

    //! \return the names of the registered linear solvers
    template <typename Entries>
    std::vector<std::string> getLinearSolverNames(const Entries& entries) {
      auto names = std::vector<std::string>{};
      names.reserve(entries.size());
      for (const auto& e : entries) {
        names.push_back(e.first);
      }
      return names;
    }  // end of getLinearSolverNames

    /*!
     * \return the description of the linear solvers registered without
     * description: they compute where the device of MFEM computes, as all
     * the linear solvers did before the descriptions
     */
    LinearSolverDescription getDefaultLinearSolverDescription() {
      return LinearSolverDescription{
          .library = "",
          .kind = LinearSolverKind::ITERATIVE,
          .placement = LinearSolverPlacement::MFEM,
          .preconditioners = AcceptedPreconditioners::ANY,
          .options = {},
          .description = ""};
    }  // end of getDefaultLinearSolverDescription

  }  // end of namespace

#ifdef MFEM_USE_MPI

  LinearSolverFactory<true>& LinearSolverFactory<true>::getFactory() {
    static LinearSolverFactory<true> factory;
    return factory;
  }  // end of getFactory

  bool LinearSolverFactory<true>::add(Context& ctx,
                                      std::string_view n,
                                      LinearSolverDescription d,
                                      Generator g) noexcept {
    return addLinearSolver(ctx, this->entries, n,
                           Entry{std::move(d), std::move(g)});
  }  // end of add

  bool LinearSolverFactory<true>::add(Context& ctx,
                                      std::string_view n,
                                      Generator g) noexcept {
    return this->add(ctx, n, getDefaultLinearSolverDescription(), std::move(g));
  }  // end of add

  LinearSolverHandler LinearSolverFactory<true>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<true>& fespace,
      const Parameters& params) const {
    return generateLinearSolver<true>(ctx, this->entries, n, fespace, params);
  }  // end of generate

  const LinearSolverDescription* LinearSolverFactory<true>::getDescription(
      std::string_view n) const noexcept {
    const auto pe = this->entries.find(n);
    return pe == this->entries.end() ? nullptr : &(pe->second.description);
  }  // end of getDescription

  std::vector<std::string> LinearSolverFactory<true>::getNames() const {
    return getLinearSolverNames(this->entries);
  }  // end of getNames

  LinearSolverFactory<true>::LinearSolverFactory() {
    internals::declareLinearSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<true>::~LinearSolverFactory() = default;

#endif /* MFEM_USE_MPI */

  LinearSolverFactory<false>& LinearSolverFactory<false>::getFactory() {
    static LinearSolverFactory<false> factory;
    return factory;
  }  // end of getFactory

  bool LinearSolverFactory<false>::add(Context& ctx,
                                       std::string_view n,
                                       LinearSolverDescription d,
                                       Generator g) noexcept {
    return addLinearSolver(ctx, this->entries, n,
                           Entry{std::move(d), std::move(g)});
  }  // end of add

  bool LinearSolverFactory<false>::add(Context& ctx,
                                       std::string_view n,
                                       Generator g) noexcept {
    return this->add(ctx, n, getDefaultLinearSolverDescription(), std::move(g));
  }  // end of add

  LinearSolverHandler LinearSolverFactory<false>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<false>& fespace,
      const Parameters& params) const {
    return generateLinearSolver<false>(ctx, this->entries, n, fespace, params);
  }  // end of generate

  const LinearSolverDescription* LinearSolverFactory<false>::getDescription(
      std::string_view n) const noexcept {
    const auto pe = this->entries.find(n);
    return pe == this->entries.end() ? nullptr : &(pe->second.description);
  }  // end of getDescription

  std::vector<std::string> LinearSolverFactory<false>::getNames() const {
    return getLinearSolverNames(this->entries);
  }  // end of getNames

  LinearSolverFactory<false>::LinearSolverFactory() {
    internals::declareLinearSolvers(*this);
  }  // end of LinearSolverFactory

  LinearSolverFactory<false>::~LinearSolverFactory() = default;

}  // end of namespace mfem_mgis
