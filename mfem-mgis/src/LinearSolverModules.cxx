/*!
 * \file   src/LinearSolverModules.cxx
 * \brief  This file implements the functions shared by the modules of the
 * linear solvers and of the preconditioners, and calls these modules
 * \date   15/10/2026
 */

#include <tuple>
#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <string_view>
#include "mfem/config/config.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Parameters.hxx"
#include "LinearSolverModules.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return true if the parameters of a preconditioner select no
   * preconditioner (name `None`, without options), false otherwise, or an
   * empty value on error. `generatePreconditioner` returns a null
   * preconditioner both in this case and on error: its callers must test
   * this case before.
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters of the preconditioner
   */
  [[nodiscard]] static std::optional<bool> isNoPreconditioner(
      Context& ctx, const Parameters& p) {
    const auto oname = get<std::string>(ctx, p, "Name");
    if (isInvalid(oname)) {
      return {};
    }
    if (*oname != "None") {
      return false;
    }
    if (contains(p, "Options")) {
      return ctx.registerErrorMessage(
          "setLinearSolverPreconditioner: no options expected for "
          "preconditioner 'None'");
    }
    if (!checkParameters(ctx, p, std::vector<std::string>{"Name"})) {
      return {};
    }
    return true;
  }  // end of isNoPreconditioner

  std::string getUnavailableLibraryHint(const std::string_view n) {
#ifndef MFEM_MGIS_HAVE_AMGX
    if ((n == "AmgX") || (n == "AmgXPCG")) {
      return "mfem-mgis is not compiled with AmgX (option enable-amgx)";
    }
#endif /* MFEM_MGIS_HAVE_AMGX */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_TRILINOS))
    if ((n == "BelosSolver") || (n == "Amesos2Solver") || (n == "Ifpack2") ||
        (n == "MueLu") || (n == "FROSch") || (n == "Amesos2Preconditioner")) {
      return "mfem-mgis is not compiled with Trilinos (option "
             "enable-trilinos, in parallel)";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_MGIS_HAVE_TRILINOS)) */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS))
    if ((n == "MUMPSSolver") || (n == "MUMPSPreconditioner")) {
      return "MFEM is not compiled with MUMPS";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)) */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS))
    if (n == "CuDSSSolver") {
      return "MFEM is not compiled with cuDSS";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_USE_CUDSS)) */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_USE_SUPERLU))
    if (n == "SuperLUSolver") {
      return "MFEM is not compiled with SuperLU_DIST";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_USE_SUPERLU)) */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_USE_STRUMPACK))
    if (n.starts_with("STRUMPACK")) {
      return "MFEM is not compiled with STRUMPACK";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_USE_STRUMPACK)) */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC))
    if ((n == "PETScLinearSolver") || (n == "PETScPreconditioner")) {
      return "MFEM is not compiled with PETSc";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_USE_PETSC)) */
#ifndef MFEM_USE_GINKGO
    if (n.starts_with("Ginkgo")) {
      return "MFEM is not compiled with Ginkgo";
    }
#endif /* MFEM_USE_GINKGO */
#ifndef MFEM_USE_MKL_PARDISO
    if (n == "MKLPardisoSolver") {
      return "MFEM is not compiled with the PARDISO solver of oneMKL";
    }
#endif /* MFEM_USE_MKL_PARDISO */
#if !(defined(MFEM_USE_MPI) && defined(MFEM_USE_MKL_CPARDISO))
    if (n == "MKLClusterSparseSolver") {
      return "MFEM is not compiled with the Cluster Sparse Solver of oneMKL";
    }
#endif /* !(defined(MFEM_USE_MPI) && defined(MFEM_USE_MKL_CPARDISO)) */
#ifndef MFEM_MGIS_HAVE_MAGMA_SPARSE
    if (n.starts_with("MAGMA")) {
      return "mfem-mgis is not compiled with the sparse part of MAGMA (MFEM "
             "compiled with MAGMA in double precision, and the library "
             "magma_sparse)";
    }
#endif /* MFEM_MGIS_HAVE_MAGMA_SPARSE */
#ifndef MFEM_USE_SUITESPARSE
    if ((n == "UMFPackSolver") || (n == "KLUSolver") || (n == "UMFPack") ||
        (n == "KLU")) {
      return "MFEM is not compiled with SuiteSparse";
    }
#endif /* MFEM_USE_SUITESPARSE */
    static_cast<void>(n);
    return {};
  }  // end of getUnavailableLibraryHint

  template <bool parallel>
  std::unique_ptr<LinearSolverPreconditioner> generatePreconditioner(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const Parameters& p) {
    if (!checkParameters(ctx, p, std::vector<std::string>{"Name", "Options"})) {
      return {};
    }
    const auto oname = get<std::string>(ctx, p, "Name");
    if (isInvalid(oname)) {
      return {};
    }
    if (*oname == "None") {
      if (contains(p, "Options")) {
        return ctx.registerErrorMessage(
            "setLinearSolverPreconditioner: "
            "no options expected for preconditioner 'None'");
      }
      return {};
    }
    const auto oopts = get_if<Parameters>(ctx, p, "Options", Parameters{});
    if (isInvalid(oopts)) {
      return {};
    }
    return PreconditionerFactory<parallel>::getFactory().generate(
        ctx, *oname, fespace, *oopts);
  }  // end of generatePreconditioner

  template <bool parallel>
  std::optional<std::unique_ptr<LinearSolverPreconditioner>> makePreconditioner(
      Context& ctx,
      FiniteElementSpace<parallel>& fespace,
      const std::string_view n,
      const LinearSolverDescription& d,
      const Parameters& params) {
    const char* const Preconditioner = "Preconditioner";
    if (!contains(params, Preconditioner)) {
      return std::unique_ptr<LinearSolverPreconditioner>{};
    }
    const auto opr = get<Parameters>(ctx, params, Preconditioner);
    if (isInvalid(opr)) {
      return {};
    }
    const auto onone = isNoPreconditioner(ctx, *opr);
    if (isInvalid(onone)) {
      return {};
    }
    if (*onone) {
      return std::unique_ptr<LinearSolverPreconditioner>{};
    }
    if (d.preconditioners == AcceptedPreconditioners::NONE) {
      std::ignore =
          ctx.registerErrorMessage("the linear solver '" + std::string{n} +
                                   "' accepts no preconditioner");
      return {};
    }
    if (d.preconditioners == AcceptedPreconditioners::SAME_LIBRARY) {
      const auto oname = get<std::string>(ctx, *opr, "Name");
      if (isInvalid(oname)) {
        return {};
      }
      const auto* const pd =
          PreconditionerFactory<parallel>::getFactory().getDescription(*oname);
      if ((pd != nullptr) && (pd->library != d.library)) {
        std::ignore = ctx.registerErrorMessage(
            "the linear solver '" + std::string{n} +
            "' only accepts the preconditioners of " + d.library + " ('" +
            *oname + "' is a preconditioner of " + pd->library + ")");
        return {};
      }
    }
    auto prec = generatePreconditioner<parallel>(ctx, fespace, *opr);
    if (prec == nullptr) {
      return {};
    }
    return prec;
  }  // end of makePreconditioner

  template <bool parallel>
  void declareLinearSolvers(LinearSolverFactory<parallel>& f) {
    declareMFEMLinearSolvers(f);
#ifdef MFEM_USE_GINKGO
    declareGinkgoLinearSolvers(f);
#endif /* MFEM_USE_GINKGO */
#ifdef MFEM_MGIS_HAVE_AMGX
    declareAmgXLinearSolvers(f);
#endif /* MFEM_MGIS_HAVE_AMGX */
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      declareHypreLinearSolvers(f);
#ifdef MFEM_USE_MUMPS
      declareMUMPSLinearSolvers(f);
#endif /* MFEM_USE_MUMPS */
#ifdef MFEM_USE_CUDSS
      declareCuDSSLinearSolvers(f);
#endif /* MFEM_USE_CUDSS */
#ifdef MFEM_USE_SUPERLU
      declareSuperLULinearSolvers(f);
#endif /* MFEM_USE_SUPERLU */
#ifdef MFEM_USE_STRUMPACK
      declareSTRUMPACKLinearSolvers(f);
#endif /* MFEM_USE_STRUMPACK */
#ifdef MFEM_USE_PETSC
      declarePETScLinearSolvers(f);
#endif /* MFEM_USE_PETSC */
#ifdef MFEM_USE_MKL_CPARDISO
      declareMKLClusterSparseLinearSolvers(f);
#endif /* MFEM_USE_MKL_CPARDISO */
#ifdef MFEM_MGIS_HAVE_TRILINOS
      declareTrilinosLinearSolvers(f);
#endif /* MFEM_MGIS_HAVE_TRILINOS */
    }
#endif /* MFEM_USE_MPI */
#ifdef MFEM_USE_MKL_PARDISO
    if constexpr (!parallel) {
      declareMKLPardisoLinearSolvers(f);
    }
#endif /* MFEM_USE_MKL_PARDISO */
#ifdef MFEM_MGIS_HAVE_MAGMA_SPARSE
    if constexpr (!parallel) {
      declareMAGMALinearSolvers(f);
    }
#endif /* MFEM_MGIS_HAVE_MAGMA_SPARSE */
  }  // end of declareLinearSolvers

  template <bool parallel>
  void declarePreconditioners(PreconditionerFactory<parallel>& f) {
    declareMFEMPreconditioners(f);
    declareNodalBlockJacobiPreconditioner(f);
#ifdef MFEM_USE_GINKGO
    declareGinkgoPreconditioners(f);
#endif /* MFEM_USE_GINKGO */
#ifdef MFEM_MGIS_HAVE_AMGX
    declareAmgXPreconditioners(f);
#endif /* MFEM_MGIS_HAVE_AMGX */
#ifdef MFEM_MGIS_HAVE_MAGMA_SPARSE
    declareMAGMAPreconditioners(f);
#endif /* MFEM_MGIS_HAVE_MAGMA_SPARSE */
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      declareHyprePreconditioners(f);
#ifdef MFEM_USE_MUMPS
      declareMUMPSPreconditioners(f);
#endif /* MFEM_USE_MUMPS */
#ifdef MFEM_USE_STRUMPACK
      declareSTRUMPACKPreconditioners(f);
#endif /* MFEM_USE_STRUMPACK */
#ifdef MFEM_USE_PETSC
      declarePETScPreconditioners(f);
#endif /* MFEM_USE_PETSC */
#ifdef MFEM_MGIS_HAVE_TRILINOS
      declareTrilinosPreconditioners(f);
#endif /* MFEM_MGIS_HAVE_TRILINOS */
    }
#endif /* MFEM_USE_MPI */
  }  // end of declarePreconditioners

#ifdef MFEM_USE_MPI
  template std::unique_ptr<LinearSolverPreconditioner>
  generatePreconditioner<true>(Context&,
                               FiniteElementSpace<true>&,
                               const Parameters&);
  template std::optional<std::unique_ptr<LinearSolverPreconditioner>>
  makePreconditioner<true>(Context&,
                           FiniteElementSpace<true>&,
                           std::string_view,
                           const LinearSolverDescription&,
                           const Parameters&);
  template void declareLinearSolvers<true>(LinearSolverFactory<true>&);
  template void declarePreconditioners<true>(PreconditionerFactory<true>&);
#endif /* MFEM_USE_MPI */
  template std::unique_ptr<LinearSolverPreconditioner>
  generatePreconditioner<false>(Context&,
                                FiniteElementSpace<false>&,
                                const Parameters&);
  template std::optional<std::unique_ptr<LinearSolverPreconditioner>>
  makePreconditioner<false>(Context&,
                            FiniteElementSpace<false>&,
                            std::string_view,
                            const LinearSolverDescription&,
                            const Parameters&);
  template void declareLinearSolvers<false>(LinearSolverFactory<false>&);
  template void declarePreconditioners<false>(PreconditionerFactory<false>&);

}  // end of namespace mfem_mgis::internals
