/*!
 * \file   src/PreconditionerFactory.cxx
 * \brief  This file implements the `PreconditionerFactory` class
 * \date   15/10/2026
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
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "ExecutionSettings.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis {

  template <bool parallel>
  PreconditionerFactory<parallel>&
  PreconditionerFactory<parallel>::getFactory() {
    static PreconditionerFactory factory;
    return factory;
  }  // end of getFactory

  template <bool parallel>
  bool PreconditionerFactory<parallel>::add(Context& ctx,
                                            std::string_view n,
                                            LinearSolverDescription d,
                                            Generator g) noexcept {
    if (this->entries.find(n) != this->entries.end()) {
      return ctx.registerErrorMessage(
          "PreconditionerFactory::add: a preconditioner called '" +
          std::string{n} + "' has already been declared");
    }
    this->entries.insert({std::string{n}, Entry{std::move(d), std::move(g)}});
    return true;
  }  // end of add

  template <bool parallel>
  std::unique_ptr<LinearSolverPreconditioner>
  PreconditionerFactory<parallel>::generate(
      Context& ctx,
      std::string_view n,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& opts) const {
    const auto pe = this->entries.find(n);
    if (pe == this->entries.end()) {
#ifdef MFEM_USE_MPI
      if constexpr (!parallel) {
        if (PreconditionerFactory<true>::getFactory().getDescription(n) !=
            nullptr) {
          return ctx.registerErrorMessage(
              "setLinearSolverPreconditioner: the preconditioner '" +
              std::string{n} + "' is only available in parallel");
        }
      } else {
        if (PreconditionerFactory<false>::getFactory().getDescription(n) !=
            nullptr) {
          return ctx.registerErrorMessage(
              "setLinearSolverPreconditioner: the preconditioner '" +
              std::string{n} + "' is only available in sequential");
        }
      }
#endif /* MFEM_USE_MPI */
      const auto hint = internals::getUnavailableLibraryHint(n);
      if (!hint.empty()) {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: unsupported preconditioner '" +
            std::string{n} + "': " + hint);
      }
      auto msg = "setLinearSolverPreconditioner: unsupported preconditioner '" +
                 std::string{n} +
                 "'. Here is the list of available "
                 "preconditioners:";
      for (const auto& e : this->entries) {
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
    try {
      return e.generator(ctx, fespace, opts);
    } catch (std::exception& ex) {
      return ctx.registerErrorMessage(
          "PreconditionerFactory::generate: error while generating "
          "preconditioner '" +
          std::string{n} + "'\n" + ex.what());
    }
  }  // end of generate

  template <bool parallel>
  const LinearSolverDescription*
  PreconditionerFactory<parallel>::getDescription(
      std::string_view n) const noexcept {
    const auto pe = this->entries.find(n);
    return pe == this->entries.end() ? nullptr : &(pe->second.description);
  }  // end of getDescription

  template <bool parallel>
  std::vector<std::string> PreconditionerFactory<parallel>::getNames() const {
    auto names = std::vector<std::string>{};
    names.reserve(this->entries.size());
    for (const auto& e : this->entries) {
      names.push_back(e.first);
    }
    return names;
  }  // end of getNames

  template <bool parallel>
  PreconditionerFactory<parallel>::PreconditionerFactory() {
    internals::declarePreconditioners(*this);
  }  // end of PreconditionerFactory

  template <bool parallel>
  PreconditionerFactory<parallel>::~PreconditionerFactory() = default;

#ifdef MFEM_USE_MPI
  template struct PreconditionerFactory<true>;
#endif /* MFEM_USE_MPI */
  template struct PreconditionerFactory<false>;

}  // end of namespace mfem_mgis
