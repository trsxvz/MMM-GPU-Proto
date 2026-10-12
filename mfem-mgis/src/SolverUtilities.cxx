/*!
 * \file   src/SolverUtilities.cxx
 * \brief  This file implements the functions declared in
 * `MFEMMGIS/Utilities/SolverUtilities.hxx`
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#include <optional>
#include "mfem/linalg/solvers.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "AmgXSolvers.hxx"
#include "MUMPSSolvers.hxx"

namespace mfem_mgis {

  std::vector<std::string> getIterativeSolverParametersList() {
    using Problem = AbstractNonLinearEvolutionProblem;
    return {Problem::SolverVerbosityLevel, Problem::SolverRelativeTolerance,
            Problem::SolverAbsoluteTolerance,
            Problem::SolverMaximumNumberOfIterations};
  }  // end of getIterativeSolverParametersList

  template <typename SolverType>
  static void setSolverParametersImplementation(attributes::Throwing,
                                                SolverType& s,
                                                const Parameters& params) {
    using Problem = AbstractNonLinearEvolutionProblem;
    checkParameters(throwing, params, getIterativeSolverParametersList());
    // the values are read before being applied: an invalid value leaves the
    // solver unchanged
    auto read = [&params]<typename ValueType>(const char* const n) {
      return contains(params, n)
                 ? std::optional<ValueType>{get<ValueType>(throwing, params, n)}
                 : std::optional<ValueType>{};
    };
    const auto verbosity =
        read.template operator()<int>(Problem::SolverVerbosityLevel);
    const auto rtol =
        read.template operator()<double>(Problem::SolverRelativeTolerance);
    const auto atol =
        read.template operator()<double>(Problem::SolverAbsoluteTolerance);
    const auto maximum_number_of_iterations =
        read.template operator()<int>(Problem::SolverMaximumNumberOfIterations);
    if (verbosity.has_value()) {
      s.SetPrintLevel(*verbosity);
    }
    if (rtol.has_value()) {
      s.SetRelTol(*rtol);
    }
    if (atol.has_value()) {
      s.SetAbsTol(*atol);
    }
    if (maximum_number_of_iterations.has_value()) {
      s.SetMaxIter(*maximum_number_of_iterations);
    }
  }  // end of setSolverParametersImplementation

  bool setSolverParameters(Context& ctx,
                           IterativeSolver& s,
                           const Parameters& params) noexcept {
    try {
      setSolverParametersImplementation(throwing, s, params);
    } catch (...) {
      return mgis::registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of setSolverParameters

#ifdef MFEM_USE_PETSC
  bool setSolverParameters(Context& ctx,
                           mfem::PetscNonlinearSolver& s,
                           const Parameters& params) noexcept {
    try {
      setSolverParametersImplementation(throwing, s, params);
    } catch (...) {
      return mgis::registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }    // end of setSolverParameters
#endif /* MFEM_USE_PETSC */

  void setSolverParameters(attributes::Throwing,
                           IterativeSolver& s,
                           const Parameters& params) {
    setSolverParametersImplementation(throwing, s, params);
  }  // end of setSolverParameters

#ifdef MFEM_USE_PETSC
  void setSolverParameters(attributes::Throwing,
                           mfem::PetscNonlinearSolver& s,
                           const Parameters& params) {
    setSolverParametersImplementation(throwing, s, params);
  }    // end of setSolverParameters
#endif /* MFEM_USE_PETSC */

  template <typename SolverType>
  inline constexpr auto hasGetMaxIter = requires(const SolverType& l) {
    { l.GetMaxIter() } -> std::same_as<int>;
  };

  void clearHypreErrors() noexcept {
#ifdef MFEM_USE_MPI
    HYPRE_ClearAllErrors();
#endif /* MFEM_USE_MPI */
  }  // end of clearHypreErrors

  bool hasConverged(const LinearSolver& ls) noexcept {
    if (const auto* isolver = dynamic_cast<const IterativeSolver*>(&ls);
        isolver != nullptr) {
      return isolver->GetConverged();
    }
#ifdef MFEM_MGIS_HAVE_AMGX
    if (const auto* aptr = dynamic_cast<const internals::AmgXSolver*>(&ls);
        aptr != nullptr) {
      return aptr->hasConverged();
    }
#endif /* MFEM_MGIS_HAVE_AMGX */
#ifdef MFEM_USE_MPI
    auto check_hypre_solver_convergence =
        []<typename SolverType>(const SolverType* ptr) {
          if constexpr (hasGetMaxIter<SolverType>) {
            int niter;
            ptr->GetNumIterations(niter);
            return niter < ptr->GetMaxIter();
          } else {
            return true;
          }
        };
    if (const auto* hptr = dynamic_cast<const mfem::HyprePCG*>(&ls);
        hptr != nullptr) {
      return check_hypre_solver_convergence(hptr);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreGMRES*>(&ls);
        hptr != nullptr) {
      return check_hypre_solver_convergence(hptr);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreFGMRES*>(&ls);
        hptr != nullptr) {
      return check_hypre_solver_convergence(hptr);
    }
#endif /* MFEM_USE_MPI */
    return true;
  }  // end of hasConverged

  std::optional<int> getNumberOfIterationsAtConvergence(
      const LinearSolver& ls) noexcept {
    if (const auto* isolver = dynamic_cast<const IterativeSolver*>(&ls);
        isolver != nullptr) {
      return isolver->GetNumIterations();
    }
#ifdef MFEM_MGIS_HAVE_AMGX
    if (const auto* aptr = dynamic_cast<const internals::AmgXSolver*>(&ls);
        aptr != nullptr) {
      return aptr->getNumberOfIterations();
    }
#endif /* MFEM_MGIS_HAVE_AMGX */
#ifdef MFEM_USE_MPI
    auto get_hypre_solver_iterations =
        []<typename SolverType>(const SolverType* ptr) {
          int niter;
          ptr->GetNumIterations(niter);
          return niter;
        };
    if (const auto* hptr = dynamic_cast<const mfem::HyprePCG*>(&ls);
        hptr != nullptr) {
      return get_hypre_solver_iterations(hptr);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreGMRES*>(&ls);
        hptr != nullptr) {
      return get_hypre_solver_iterations(hptr);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreFGMRES*>(&ls);
        hptr != nullptr) {
      return get_hypre_solver_iterations(hptr);
    }
#endif /* MFEM_USE_MPI */
    // reached if the solver is not an iterative solver
    return {};
  }  // end of getNumberOfIterationsAtConvergence

  std::optional<size_type> getNumberOfMUMPSAnalyses(
      [[maybe_unused]] const LinearSolver& ls) noexcept {
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS)
    if (const auto* const s =
            dynamic_cast<const internals::MUMPSSolverReusingAnalysis*>(&ls);
        s != nullptr) {
      return s->getNumberOfAnalyses();
    }
#endif /* defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) */
    return {};
  }  // end of getNumberOfMUMPSAnalyses

  namespace {

    /*!
     * \brief access to the relative tolerance of the iterative solvers of
     * MFEM, which have no getter: the pointer to this protected member is
     * formed in a derived class, then applied to any iterative solver
     */
    struct IterativeSolverRelativeTolerance : IterativeSolver {
      //! \return the relative tolerance of an iterative solver
      static real get(const IterativeSolver& s) noexcept {
        return s.*(&IterativeSolverRelativeTolerance::rel_tol);
      }
    };

  }  // end of namespace

  std::optional<real> getRelativeTolerance(const LinearSolver& ls) noexcept {
    if (const auto* isolver = dynamic_cast<const IterativeSolver*>(&ls);
        isolver != nullptr) {
      return IterativeSolverRelativeTolerance::get(*isolver);
    }
#ifdef MFEM_USE_MPI
    auto tol = HYPRE_Real{};
    if (const auto* hptr = dynamic_cast<const mfem::HyprePCG*>(&ls);
        hptr != nullptr) {
      HYPRE_PCGGetTol(*hptr, &tol);
      return static_cast<real>(tol);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreGMRES*>(&ls);
        hptr != nullptr) {
      HYPRE_GMRESGetTol(*hptr, &tol);
      return static_cast<real>(tol);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreFGMRES*>(&ls);
        hptr != nullptr) {
      HYPRE_FlexGMRESGetTol(*hptr, &tol);
      return static_cast<real>(tol);
    }
#endif /* MFEM_USE_MPI */
    return {};
  }  // end of getRelativeTolerance

  bool setRelativeTolerance(LinearSolver& ls, const real tol) noexcept {
    if (auto* isolver = dynamic_cast<IterativeSolver*>(&ls);
        isolver != nullptr) {
      isolver->SetRelTol(tol);
      return true;
    }
#ifdef MFEM_USE_MPI
    if (auto* hptr = dynamic_cast<mfem::HyprePCG*>(&ls); hptr != nullptr) {
      hptr->SetTol(tol);
      return true;
    }
    if (auto* hptr = dynamic_cast<mfem::HypreGMRES*>(&ls); hptr != nullptr) {
      hptr->SetTol(tol);
      return true;
    }
    if (auto* hptr = dynamic_cast<mfem::HypreFGMRES*>(&ls); hptr != nullptr) {
      hptr->SetTol(tol);
      return true;
    }
#endif /* MFEM_USE_MPI */
    return false;
  }  // end of setRelativeTolerance
}  // end of namespace mfem_mgis
