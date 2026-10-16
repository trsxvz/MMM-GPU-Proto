#pragma once

#include <string>
#include <string_view>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cassert>
#include <algorithm>

#include <sys/resource.h>
#include <mpi.h>

#include "mfem.hpp"

#include "MFEMMGIS/Parameter.hxx"

#include "Setup.hxx"

inline bool contains(const std::vector<std::string>& list,
                     const std::string& value) {
  return std::find(list.begin(), list.end(), value) != list.end();
}

// UMFPackSolver is not listed since it is only available for sequential
// computations
inline const std::vector<std::string> direct_solvers = {"MUMPSSolver",
                                                        "CuDSSSolver"};

inline const std::vector<std::string> iterative_solvers = {
    "CGSolver",  "GMRESSolver", "BiCGSTABSolver", "MINRESSolver",
    "SLISolver", "HyprePCG",    "HypreGMRES",     "HypreFGMRES"};

// MueLu requires mfem-mgis built with Trilinos (the rigid body modes are
// given by mfem-mgis for the mechanics), with a Krylov solver of MFEM
inline const std::vector<std::string> preconditioners = {
    "HypreBoomerAMG", "HypreDiagScale", "HypreEuclid", "HypreILU",
    "HypreParaSails", "MueLu"};

inline long get_memory_checkpoint() {
  rusage obj;
  [[maybe_unused]] const auto r = getrusage(RUSAGE_SELF, &obj);
  assert((r == 0) && "error: getrusage has failed");
  long res = 0;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
}

inline void print_memory_footprint(std::string msg) {
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6;
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
}

template <typename Implementation>
inline void print_mesh_information(Implementation& impl) {
  using mfem_mgis::Profiler::Utils::Message;
  using mfem_mgis::Profiler::Utils::sum;
  Message("INFO: print_mesh_information");
  auto mesh = impl.getFiniteElementSpace().GetMesh();
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);
  auto& fespace = impl.getFiniteElementSpace();
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);
  Message("INFO: number of elements -> ", numbers_of_elements);
  Message("INFO: Number of finite element unknowns: ", unknowns);
}

template <typename Problem>
inline void add_post_processings(mfem_mgis::attributes::MayAbort,
                                 mfem_mgis::Context& ctx,
                                 Problem& p,
                                 std::string msg,
                                 std::string field_name) {
  auto or_die = ctx.getFatalFailureHandler();
  p.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", msg}, {"OutputFieldName", field_name}}) |
      or_die;
}

/*!
 * \return true if a symmetry of the matrices given to the direct solvers is
 * valid: empty (default symmetry), Unsymmetric, SymmetricIndefinite or
 * SymmetricPositiveDefinite
 * \param[in] symmetry: symmetry
 */
inline bool isValidMatrixSymmetry(const std::string_view symmetry) {
  return symmetry.empty() || (symmetry == "Unsymmetric") ||
         (symmetry == "SymmetricIndefinite") ||
         (symmetry == "SymmetricPositiveDefinite");
}  // end of isValidMatrixSymmetry

/*!
 * \return the parameters of the direct solvers of mfem-mgis (MUMPSSolver,
 * CuDSSSolver) for the symmetry of the matrices given by the option -msMc,
 * or for the default symmetry if this option is not given. The symmetric
 * modes give the lower triangle of the matrices to the solver without
 * checking their symmetry: a non symmetric matrix would slow down or break
 * the convergence of the Newton method.
 * \param[in] symmetry: value of the option -msMc, checked by
 * `isValidMatrixSymmetry`
 * \param[in] default_symmetry: default symmetry
 * \param[in] reuse_analysis: reuse the analysis of the solver for the
 * matrices of the same pattern (option -mra)
 */
inline mfem_mgis::Parameters getDirectSolverParameters(
    const std::string_view symmetry,
    const std::string_view default_symmetry,
    const bool reuse_analysis) {
  const auto s = symmetry.empty() ? default_symmetry : symmetry;
  auto parameters = mfem_mgis::Parameters{};
  if (s == "SymmetricIndefinite") {
    parameters = mfem_mgis::Parameters{{"Symmetric", true}};
  } else if (s == "SymmetricPositiveDefinite") {
    parameters =
        mfem_mgis::Parameters{{"Symmetric", true}, {"PositiveDefinite", true}};
  } else if ((!s.empty()) && (s != "Unsymmetric")) {
    mfem_mgis::abort(
        "Invalid symmetry of the matrices given to the direct solver: " +
        std::string{s} +
        " (expected Unsymmetric, SymmetricIndefinite or "
        "SymmetricPositiveDefinite)");
  }
  if (reuse_analysis) {
    // only given if requested, for the versions of mfem-mgis without this
    // parameter
    parameters.insert(mfem_mgis::attributes::Throwing{}, "ReuseAnalysis", true);
  }
  return parameters;
}  // end of getDirectSolverParameters

template <typename Problem>
inline static void setLinearSolver(mfem_mgis::attributes::MayAbort,
                                   mfem_mgis::Context& ctx,
                                   Problem& p,
                                   const std::string& physics_type,
                                   const TestParameters& param,
                                   const int verbosity = 0,
                                   const mfem_mgis::real Tol = 1e-9) {
  CatchTimeSection(ctx, "set_linear_solver");
  auto or_die = ctx.getFatalFailureHandler();
  std::string solver;
  std::string precond;

  if (physics_type == "heat_transfer") {
    solver = param.solver_thermo;
    precond = param.precond_thermo;
  } else if (physics_type == "mechanics") {
    solver = param.solver_meca;
    precond = param.precond_meca;
  } else {
    mfem_mgis::abort("unknown physics (" + physics_type + ")");
  }
  // options of the mechanical problem without effect on its linear solver
  if (physics_type == "mechanics") {
    if ((*param.amg_strategy_meca != '\0') &&
        ((!contains(iterative_solvers, solver)) ||
         (precond != "HypreBoomerAMG"))) {
      mfem_mgis::Profiler::Utils::Message(
          "warning: the option -amgMc is ignored (no BoomerAMG preconditioner "
          "for the mechanical problem)");
    }
    if ((*param.matrix_symmetry_meca != '\0') &&
        (!contains(direct_solvers, solver)) &&
        (std::string_view{param.nonlinear_solver_meca} != "NewtonK0")) {
      mfem_mgis::Profiler::Utils::Message(
          "warning: the option -msMc is ignored (the linear solver of the "
          "mechanical problem is not a direct solver)");
    }
  }

  if (contains(iterative_solvers, solver)) {
    auto solverParameters = mfem_mgis::Parameters{
        {"VerbosityLevel", verbosity}, {"MaximumNumberOfIterations", 10000}};
    if (solver.starts_with("Hypre")) {
      // relative tolerance of the Hypre solvers
      solverParameters.insert(ctx, "Tolerance", Tol) | or_die;
    } else {
      // the tolerances of the MFEM solvers are zero by default
      solverParameters.insert(ctx, "RelativeTolerance", Tol) | or_die;
      solverParameters.insert(ctx, "AbsoluteTolerance", 0.) | or_die;
    }

    if (!precond.empty()) {
      if (!contains(preconditioners, precond)) {
        mfem_mgis::abort("Invalid preconditioner: " + precond);
      }

      auto options = mfem_mgis::Parameters{};
      if (precond == "MueLu") {
        // the parameters of MueLu, without verbosity level
        const auto smoother = std::string_view{param.muelu_smoother};
        if (smoother == "jacobi") {
          const auto jacobi = mfem_mgis::Parameters{
              {"relaxation: type", "Jacobi"},
              {"relaxation: sweeps", 2},
              {"relaxation: damping factor", param.muelu_damping}};
          options.insert(
              ctx, "Parameters",
              mfem_mgis::Parameters{{"smoother: type", "RELAXATION"},
                                    {"smoother: params", jacobi}}) |
              or_die;
        } else if (smoother == "chebyshev") {
          // the default estimate of the largest eigenvalue (10 iterations)
          // may be too low, which breaks the positivity of the
          // preconditioner
          const auto chebyshev = mfem_mgis::Parameters{
              {"chebyshev: degree", 2},
              {"chebyshev: eigenvalue max iterations", 30},
              {"chebyshev: boost factor", 1.2}};
          options.insert(
              ctx, "Parameters",
              mfem_mgis::Parameters{{"smoother: type", "CHEBYSHEV"},
                                    {"smoother: params", chebyshev}}) |
              or_die;
        } else if (smoother != "default") {
          mfem_mgis::abort("invalid smoother of MueLu: " +
                           std::string{smoother});
        }
      } else {
        options.insert(ctx, "VerbosityLevel", verbosity) | or_die;
      }
      if ((physics_type == "mechanics") && (precond == "HypreBoomerAMG") &&
          (*param.amg_strategy_meca != '\0')) {
        options.insert(ctx, "Strategy", param.amg_strategy_meca) | or_die;
      }

      auto preconditioner =
          mfem_mgis::Parameters{{"Name", precond}, {"Options", options}};

      solverParameters.insert(ctx, "Preconditioner", preconditioner) | or_die;
    }

    p.setLinearSolver(ctx, solver, solverParameters) | or_die;
  }

  else if (contains(direct_solvers, solver)) {
    const auto solverParameters =
        (physics_type == "mechanics")
            ? getDirectSolverParameters(param.matrix_symmetry_meca,
                                        "Unsymmetric", param.reuse_analysis)
            : getDirectSolverParameters("", "", param.reuse_analysis);
    p.setLinearSolver(ctx, solver, solverParameters) | or_die;
  }

  else {
    mfem_mgis::abort("Unknown solver type: " + solver);
  }
}
