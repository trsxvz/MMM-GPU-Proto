/*!
 * \file   src/TrilinosBackend.cxx
 * \brief  This file implements the interface of the library of mfem-mgis
 * built against Trilinos (see `TrilinosBackend.hxx`), compiled by the
 * compiler of Kokkos: the solvers are built for the node of Tpetra on the
 * GPU or on the host (see `TrilinosTpetraBackend.hxx`)
 * \date   15/10/2026
 */

#include <cctype>
#include <algorithm>
#include "TrilinosTpetraBackend.hxx"

namespace mfem_mgis::trilinos {

  double Operator::getRelativeTolerance() const { return -1; }

  bool Operator::setRelativeTolerance(const double) { return false; }

  Operator::~Operator() = default;

  bool hasGPUNode() noexcept { return TpetraBackend<DefaultNode>::gpu; }

  bool hasHostNode() noexcept { return !TpetraBackend<HostNode>::gpu; }

  std::unique_ptr<Operator> makeIfpack2Preconditioner(const std::string& type,
                                                      const Options& options,
                                                      const bool gpu) {
    if (gpu) {
      return TpetraBackend<DefaultNode>::makeIfpack2Preconditioner(type,
                                                                   options);
    }
    return TpetraBackend<HostNode>::makeIfpack2Preconditioner(type, options);
  }  // end of makeIfpack2Preconditioner

  std::unique_ptr<Operator> makeMueLuPreconditioner(const Options& options,
                                                    NodesCoordinates c,
                                                    const bool gpu) {
    if (gpu) {
      return TpetraBackend<DefaultNode>::makeMueLuPreconditioner(options,
                                                                 std::move(c));
    }
    return TpetraBackend<HostNode>::makeMueLuPreconditioner(options,
                                                            std::move(c));
  }  // end of makeMueLuPreconditioner

  std::unique_ptr<Operator> makeFROSchPreconditioner(const Options& options,
                                                     NodesCoordinates c,
                                                     const bool) {
    // FROSch 17.2.1 runs some loops in the default execution space of
    // Kokkos, whatever the node (FROSch_ExtractSubmatrices_def.hpp,
    // FROSch_CoarseSpace_def.hpp): it computes with the default node
    return TpetraBackend<DefaultNode>::makeFROSchPreconditioner(options,
                                                                std::move(c));
  }  // end of makeFROSchPreconditioner

  std::unique_ptr<Operator> makeAmesos2Solver(const std::string& solver,
                                              const Options& options,
                                              const bool gpu) {
    // Tacho (`Tacho` or `Amesos2_Tacho`) computes in the default execution
    // space of Kokkos, whatever the node (Amesos2_Tacho_decl.hpp)
    auto name = solver;
    std::ranges::transform(name, name.begin(),
                           [](const unsigned char c) noexcept {
                             return static_cast<char>(std::tolower(c));
                           });
    if ((gpu) || (name.find("tacho") != std::string::npos)) {
      return TpetraBackend<DefaultNode>::makeAmesos2Solver(solver, options);
    }
    return TpetraBackend<HostNode>::makeAmesos2Solver(solver, options);
  }  // end of makeAmesos2Solver

  bool isAmesos2SolverAvailable(const std::string& solver) {
    return Amesos2::query(solver);
  }  // end of isAmesos2SolverAvailable

  std::unique_ptr<Operator> makeBelosSolver(const std::string& solver,
                                            const Options& options,
                                            const BelosSettings& settings,
                                            const Operator* const p,
                                            ExternalPreconditioner e,
                                            const bool gpu) {
    // the solver computes on the node of its preconditioner of Trilinos
    if ((p != nullptr) ? p->computesOnTheGPU() : gpu) {
      return TpetraBackend<DefaultNode>::makeBelosSolver(
          solver, options, settings, p, std::move(e));
    }
    return TpetraBackend<HostNode>::makeBelosSolver(solver, options, settings,
                                                    p, std::move(e));
  }  // end of makeBelosSolver

}  // end of namespace mfem_mgis::trilinos
