/*!
 * \file   src/AccelerationAlgorithm.cxx
 * \brief  This file implements the `AccelerationAlgorithm` and the
 * `AccelerationAlgorithmFactory` classes
 * \date   10/10/2026
 */

#include <string>
#include <tuple>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/ErrorBacktrace.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "AccelerationAlgorithms.hxx"
#include "MFEMMGIS/NonLinearSolvers/AccelerationAlgorithm.hxx"

namespace mfem_mgis {

  AccelerationAlgorithm::~AccelerationAlgorithm() = default;

  AccelerationAlgorithmFactory& AccelerationAlgorithmFactory::get() noexcept {
    static AccelerationAlgorithmFactory f;
    return f;
  }  // end of get

  AccelerationAlgorithmFactory::AccelerationAlgorithmFactory() noexcept {
    auto ctx = Context{};
    auto check = [&ctx](const bool b) {
      if (!b) {
        auto eh = mgis::ContractViolationHandler{};
        eh.registerErrorMessage(ctx.getErrorMessage().c_str());
      }
    };
    check(this->add(ctx, "Anderson",
                    internals::makeAndersonAccelerationAlgorithm));
    check(this->add(ctx, "Cast3M", internals::makeCast3MAccelerationAlgorithm));
  }  // end of AccelerationAlgorithmFactory

  bool AccelerationAlgorithmFactory::add(Context& ctx,
                                         std::string_view n,
                                         Generator g) noexcept {
    if (!g) {
      return ctx.registerErrorMessage("invalid generator");
    }
    if (this->generators.contains(n)) {
      return ctx.registerErrorMessage("an acceleration algorithm named '" +
                                      std::string{n} +
                                      "' has already been declared");
    }
    this->generators.insert({std::string{n}, std::move(g)});
    return true;
  }  // end of add

  std::unique_ptr<AccelerationAlgorithm> AccelerationAlgorithmFactory::generate(
      Context& ctx, std::string_view n, const Parameters& opts) const noexcept {
    const auto p = this->generators.find(n);
    if (p == this->generators.end()) {
      return ctx.registerErrorMessage("no acceleration algorithm named '" +
                                      std::string{n} + "' declared");
    }
    try {
      return p->second(ctx, opts);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of generate

  AccelerationAlgorithmFactory::~AccelerationAlgorithmFactory() noexcept =
      default;

}  // end of namespace mfem_mgis
