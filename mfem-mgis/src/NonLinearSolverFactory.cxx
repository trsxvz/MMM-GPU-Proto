/*!
 * \file   NonLinearSolverFactory.cxx
 * \brief  This file implements the `NonLinearSolverFactory` class
 * \author Thomas Helfer
 * \date   20/09/2026
 */

#include <string>
#include <string_view>
#include "MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx"
#include "MFEMMGIS/NonLinearSolvers/FixedPointSolver.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.hxx"
#include "KINSOLSolver.hxx"
#include "PETScSNESSolver.hxx"

namespace mfem_mgis {

  AbstractNonLinearSolverGenerator::
      ~AbstractNonLinearSolverGenerator() noexcept = default;

  NonLinearSolverFactory& NonLinearSolverFactory::get() noexcept {
    static NonLinearSolverFactory f;
    return f;
  }

  /*!
   * \return why an optional nonlinear solver is not declared (its library
   * is not available), or an empty string
   * \param[in] n: name of the nonlinear solver
   */
  static std::string getUnavailableNonLinearSolverHint(
      const std::string_view n) {
#ifndef MFEM_MGIS_HAVE_SUNDIALS
    if (n == "KINSOL") {
      return " (mfem-mgis does not use SUNDIALS: option enable-sundials, MFEM "
             "compiled with SUNDIALS)";
    }
#endif /* MFEM_MGIS_HAVE_SUNDIALS */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC))
    if (n == "SNES") {
      return " (mfem-mgis does not use PETSc: option enable-petsc, MFEM "
             "compiled with PETSc and MPI)";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC)) */
    static_cast<void>(n);
    return {};
  }  // end of getUnavailableNonLinearSolverHint

  NonLinearSolverFactory::NonLinearSolverFactory() noexcept {
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    auto g = make_unique<StandardNonLinearSolverGenerator<NewtonSolver>>(ctx) |
             or_die;
    this->add(ctx, "Newton", std::move(g)) | or_die;
    auto g2 =
        make_unique<StandardNonLinearSolverGenerator<FixedPointSolver>>(ctx) |
        or_die;
    this->add(ctx, "FixedPoint", std::move(g2)) | or_die;
#ifdef MFEM_MGIS_HAVE_SUNDIALS
    auto g3 =
        make_unique<StandardNonLinearSolverGenerator<internals::KINSOLSolver>>(
            ctx) |
        or_die;
    this->add(ctx, "KINSOL", std::move(g3)) | or_die;
#endif /* MFEM_MGIS_HAVE_SUNDIALS */
#if defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC)
    auto g4 = make_unique<
                  StandardNonLinearSolverGenerator<internals::PETScSNESSolver>>(
                  ctx) |
              or_die;
    this->add(ctx, "SNES", std::move(g4)) | or_die;
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_PETSC) */
  }  // end of NonLinearSolverFactory

  bool NonLinearSolverFactory::add(
      Context& ctx,
      std::string_view n,
      std::unique_ptr<AbstractNonLinearSolverGenerator> g) noexcept {
    if (g.get() == nullptr) {
      return ctx.registerErrorMessage("invalid generator");
    }
    if (this->generators.contains(n)) {
      return ctx.registerErrorMessage("a nonlinear solver named '" +
                                      std::string{n} +
                                      "' has already been declared");
    }
    this->generators.insert({std::string{n}, std::move(g)});
    return true;
  }  // end of add

#ifdef MFEM_USE_MPI

  std::unique_ptr<AbstractNonLinearSolver> NonLinearSolverFactory::generate(
      Context& ctx,
      std::string_view n,
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& parameters) const noexcept {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      return ctx.registerErrorMessage("no nonlinear solver named '" +
                                      std::string{n} + "' declared" +
                                      getUnavailableNonLinearSolverHint(n));
    }
    auto ptr = pg->second->operator()(ctx, p, parameters);
    if (isInvalid(ptr)) {
      return ctx.registerErrorMessage("generation of nonlinear solver '" +
                                      std::string{n} + "' failed");
    }
    return ptr;
  }  // end of generate

#endif /* MFEM_USE_MPI */

  std::unique_ptr<AbstractNonLinearSolver> NonLinearSolverFactory::generate(
      Context& ctx,
      std::string_view n,
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& parameters) const noexcept {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      return ctx.registerErrorMessage("no nonlinear solver named '" +
                                      std::string{n} + "' declared" +
                                      getUnavailableNonLinearSolverHint(n));
    }
    auto ptr = pg->second->operator()(ctx, p, parameters);
    if (isInvalid(ptr)) {
      return ctx.registerErrorMessage("generation of nonlinear solver '" +
                                      std::string{n} + "' failed");
    }
    return ptr;
  }  // end of generate

  NonLinearSolverFactory::~NonLinearSolverFactory() = default;

}  // end of namespace mfem_mgis