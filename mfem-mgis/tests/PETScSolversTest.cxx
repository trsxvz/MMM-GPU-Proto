/*!
 * \file   PETScSolversTest.cxx
 * \brief  Test of the Krylov solver (`PETScLinearSolver`) and of the
 * preconditioner (`PETScPreconditioner`) of PETSc built by the factories,
 * on the Laplace equation on a cube, in parallel (the generic test
 * `LinearSolverFactoriesTest` solves it with the default options):
 *
 * - the options of PETSc given by the `Options` parameter reach the solver
 *   (types of the Krylov solver and of the preconditioner, tolerances,
 *   options of each type of value and option without value), and the
 *   options of two solvers do not interfere;
 * - the options of mfem-mgis set the tolerances, the options of PETSc
 *   overriding them, and the relative tolerance can be read and changed
 *   (see `getRelativeTolerance` and `setRelativeTolerance`);
 * - a resolution stopped by the maximum number of iterations is reported
 *   as not converged;
 * - the options of a solver are removed from the options database of PETSc
 *   when the solver is destroyed;
 * - the algebraic multigrid of PETSc (GAMG) reduces the number of
 *   iterations compared to a Jacobi preconditioner;
 * - a preconditioner of mfem-mgis is used as a shell preconditioner, which
 *   receives the operators of different sizes given to the solver;
 * - the preconditioner of PETSc preconditions the Krylov solvers of PETSc
 *   and of MFEM;
 * - invalid options are refused: an invalid value, a leading dash, and the
 *   options of the preconditioner of PETSc given to a solver with a
 *   preconditioner of mfem-mgis.
 * \date   15/10/2026
 */

#include <cmath>
#include <string>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "mfem.hpp"
#include "petscksp.h"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "LaplaceSystem.hxx"

using mfem_mgis::tests::LaplaceSystem;

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief outcome of a resolution
struct Resolution {
  //! \brief if true, the resolution converged
  bool converged = false;
  //! \brief number of iterations
  int iterations = 0;
  //! \brief norm of the residual relative to the right hand side
  double residual = 1;
};

/*!
 * \return true if a resolution converged in at least one iteration with a
 * relative residual lower than 1e-8
 * \param[in] r: resolution
 */
static bool isSolved(const Resolution& r) {
  return r.converged && (r.iterations > 0) && (r.residual < 1e-8);
}  // end of isSolved

/*!
 * \return true if a multigrid preconditioner is active: the number of
 * iterations is less than two thirds of the one of the Jacobi
 * preconditioner (about 10 iterations against 25 on this mesh)
 * \param[in] r: resolution with the multigrid preconditioner
 * \param[in] rj: resolution with the Jacobi preconditioner
 */
static bool isFaster(const Resolution& r, const Resolution& rj) {
  return (r.iterations > 0) && (3 * r.iterations < 2 * rj.iterations);
}  // end of isFaster

/*!
 * \return the linear solver of the factory, throws on failure
 * \param[in] s: system
 * \param[in] name: name of the linear solver
 * \param[in] parameters: parameters of the linear solver
 */
static mfem_mgis::LinearSolverHandler generate(
    LaplaceSystem<true>& s,
    const std::string& name,
    const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<true>::getFactory().generate(
      ctx, name, s.fespace, parameters);
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  return ls;
}  // end of generate

/*!
 * \return the error message of the generation of a linear solver which
 * shall fail, or an empty string if the generation succeeds
 * \param[in] s: system
 * \param[in] parameters: parameters of `PETScLinearSolver`
 */
static std::string getGenerationError(LaplaceSystem<true>& s,
                                      const mfem_mgis::Parameters& parameters) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  const auto ls = LinearSolverFactory<true>::getFactory().generate(
      ctx, "PETScLinearSolver", s.fespace, parameters);
  return isInvalid(ls) ? ctx.getErrorMessage() : "";
}  // end of getGenerationError

//! \return the resolution of the system by a linear solver
static Resolution solve(LaplaceSystem<true>& s, mfem::Solver& ls) {
  using namespace mfem_mgis;
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.SetOperator(s.A);
  ls.Mult(s.B, x);
  auto r = Resolution{};
  r.converged = hasConverged(ls);
  r.iterations = getNumberOfIterationsAtConvergence(ls).value_or(-1);
  r.residual = s.getRelativeResidual(x);
  if (rank == 0) {
    std::cout << r.iterations << " iterations, residual " << r.residual
              << (r.converged ? "" : " (not converged)") << '\n';
  }
  return r;
}  // end of solve

//! \return the Krylov solver of PETSc of a linear solver of the factory
static KSP getKSP(const mfem_mgis::LinearSolverHandler& ls) {
  const auto* const p =
      dynamic_cast<const mfem::PetscLinearSolver*>(ls.linear_solver.get());
  if (p == nullptr) {
    throw std::runtime_error("not a linear solver of PETSc");
  }
  return KSP(*p);
}  // end of getKSP

//! \return the type of the Krylov solver of PETSc
static std::string getKSPType(const mfem_mgis::LinearSolverHandler& ls) {
  auto type = KSPType{};
  KSPGetType(getKSP(ls), &type);
  return type == nullptr ? "" : type;
}  // end of getKSPType

//! \return the type of a preconditioner of PETSc
static std::string getPCType(PC pc) {
  auto type = PCType{};
  PCGetType(pc, &type);
  return type == nullptr ? "" : type;
}  // end of getPCType

//! \return the type of the preconditioner of the Krylov solver of PETSc
static std::string getPCType(const mfem_mgis::LinearSolverHandler& ls) {
  auto pc = PC{};
  KSPGetPC(getKSP(ls), &pc);
  return getPCType(pc);
}  // end of getPCType

//! \brief tolerances of a Krylov solver of PETSc
struct Tolerances {
  //! \brief relative tolerance
  double relative = 0;
  //! \brief absolute tolerance
  double absolute = 0;
  //! \brief maximum number of iterations
  int maximum_number_of_iterations = 0;
};

//! \return the tolerances of the Krylov solver of PETSc
static Tolerances getTolerances(const mfem_mgis::LinearSolverHandler& ls) {
  auto rtol = PetscReal{};
  auto atol = PetscReal{};
  auto dtol = PetscReal{};
  auto maxits = PetscInt{};
  KSPGetTolerances(getKSP(ls), &rtol, &atol, &dtol, &maxits);
  return {.relative = static_cast<double>(rtol),
          .absolute = static_cast<double>(atol),
          .maximum_number_of_iterations = static_cast<int>(maxits)};
}  // end of getTolerances

/*!
 * \return true if an option is in the options database of PETSc
 * \param[in] prefix: prefix of the option
 * \param[in] n: name of the option, without the leading dash
 */
static bool hasOption(const std::string& prefix, const std::string& n) {
  auto found = PetscBool{};
  const auto option = "-" + n;
  PetscOptionsHasName(nullptr, prefix.c_str(), option.c_str(), &found);
  return found == PETSC_TRUE;
}  // end of hasOption

//! \return true if two values are equal up to a relative precision of 1e-12
static bool areClose(const double a, const double b) {
  return std::abs(a - b) <= 1e-12 * std::abs(b);
}  // end of areClose

/*!
 * \brief checks the options and the tolerances
 * \param[in] s: system
 */
static void checkOptions(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  // options of each type of value, and an option without value
  auto gmres = generate(
      s, "PETScLinearSolver",
      Parameters{{"Options", Parameters{{"ksp_type", "gmres"},
                                        {"ksp_rtol", 1e-11},
                                        {"ksp_max_it", 700},
                                        {"ksp_gmres_restart", 40},
                                        {"ksp_error_if_not_converged", true},
                                        {"ksp_initial_guess_nonzero", ""}}}});
  check(isSolved(solve(s, *(gmres.linear_solver))), "GMRES: resolution");
  check(getKSPType(gmres) == "gmres", "GMRES: type of solver");
  auto restart = PetscInt{};
  KSPGMRESGetRestart(getKSP(gmres), &restart);
  check(restart == 40, "GMRES: restart");
  auto error = PetscBool{};
  KSPGetErrorIfNotConverged(getKSP(gmres), &error);
  check(error == PETSC_TRUE, "GMRES: boolean option");
  auto nonzero = PetscBool{};
  KSPGetInitialGuessNonzero(getKSP(gmres), &nonzero);
  check(nonzero == PETSC_TRUE, "GMRES: option without value");
  const auto t = getTolerances(gmres);
  check(areClose(t.relative, 1e-11) && (t.maximum_number_of_iterations == 700),
        "GMRES: tolerances");
  // options of mfem-mgis, overridden by the ones of PETSc
  auto tolerances =
      generate(s, "PETScLinearSolver",
               Parameters{{"RelativeTolerance", 1e-9},
                          {"AbsoluteTolerance", 1e-30},
                          {"MaximumNumberOfIterations", 300},
                          {"Options", Parameters{{"ksp_max_it", 200}}}});
  static_cast<void>(solve(s, *(tolerances.linear_solver)));
  const auto t2 = getTolerances(tolerances);
  check(areClose(t2.relative, 1e-9) && areClose(t2.absolute, 1e-30) &&
            (t2.maximum_number_of_iterations == 200),
        "options of mfem-mgis");
  // access to the relative tolerance
  const auto rtol = getRelativeTolerance(*(tolerances.linear_solver));
  check(rtol.has_value() && areClose(*rtol, 1e-9),
        "relative tolerance not read");
  check(setRelativeTolerance(*(tolerances.linear_solver), 1e-7) &&
            areClose(getTolerances(tolerances).relative, 1e-7),
        "relative tolerance not set");
  // the options of a solver are removed when it is destroyed
  auto prefix = static_cast<const char*>(nullptr);
  KSPGetOptionsPrefix(getKSP(tolerances), &prefix);
  const auto removed_prefix = std::string{prefix == nullptr ? "" : prefix};
  check(hasOption(removed_prefix, "ksp_max_it"), "options of a solver");
  tolerances = LinearSolverHandler{};
  check(!hasOption(removed_prefix, "ksp_max_it"),
        "options of a destroyed solver");
  // a resolution stopped by the maximum number of iterations
  auto stopped =
      generate(s, "PETScLinearSolver",
               Parameters{{"RelativeTolerance", 1e-10},
                          {"Options", Parameters{{"ksp_type", "cg"},
                                                 {"pc_type", "jacobi"},
                                                 {"ksp_max_it", 2}}}});
  const auto rs = solve(s, *(stopped.linear_solver));
  check((!rs.converged) && (rs.iterations == 2),
        "resolution stopped by the maximum number of iterations");
}  // end of checkOptions

/*!
 * \brief checks the preconditioners
 * \param[in] s: system
 */
static void checkPreconditioners(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  // types of the solver and of the preconditioner, two solvers alive
  auto jacobi =
      generate(s, "PETScLinearSolver",
               Parameters{{"Options", Parameters{{"ksp_type", "cg"},
                                                 {"pc_type", "jacobi"},
                                                 {"ksp_rtol", 1e-10}}}});
  auto gamg =
      generate(s, "PETScLinearSolver",
               Parameters{{"Options", Parameters{{"ksp_type", "cg"},
                                                 {"pc_type", "gamg"},
                                                 {"ksp_rtol", 1e-10}}}});
  const auto rj = solve(s, *(jacobi.linear_solver));
  const auto rg = solve(s, *(gamg.linear_solver));
  check(getKSPType(jacobi) == "cg", "CG, Jacobi: type of solver");
  check(getPCType(jacobi) == "jacobi", "CG, Jacobi: type of preconditioner");
  check(getPCType(gamg) == "gamg", "CG, GAMG: type of preconditioner");
  check(isSolved(rj) && isSolved(rg), "CG: resolution");
  check(isFaster(rg, rj),
        "CG, GAMG: the multigrid does not reduce the number of iterations");
  // a preconditioner of mfem-mgis as a shell preconditioner, with operators
  // of different sizes
  auto amg = generate(
      s, "PETScLinearSolver",
      Parameters{{"Options", Parameters{{"ksp_type", "cg"}}},
                 {"RelativeTolerance", 1e-10},
                 {"Preconditioner",
                  Parameters{{"Name", "HypreBoomerAMG"},
                             {"Options", Parameters{{"VerbosityLevel", 0}}}}}});
  const auto ramg = solve(s, *(amg.linear_solver));
  check(getPCType(amg) == "shell", "CG, BoomerAMG: shell preconditioner");
  check(isSolved(ramg) && isFaster(ramg, rj), "CG, BoomerAMG: resolution");
  auto mesh =
      mfem::Mesh::MakeCartesian3D(12, 12, 12, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  auto s2 = LaplaceSystem<true>(mesh, fec);
  check(isSolved(solve(s2, *(amg.linear_solver))),
        "CG, BoomerAMG: resolution on a smaller mesh");
  // the preconditioner of PETSc with the Krylov solver of PETSc
  const auto petsc_gamg = Parameters{
      {"Name", "PETScPreconditioner"},
      {"Options", Parameters{{"Options", Parameters{{"pc_type", "gamg"}}}}}};
  auto ksp = generate(s, "PETScLinearSolver",
                      Parameters{{"Options", Parameters{{"ksp_type", "cg"}}},
                                 {"RelativeTolerance", 1e-10},
                                 {"Preconditioner", petsc_gamg}});
  const auto rksp = solve(s, *(ksp.linear_solver));
  check(getPCType(ksp) == "gamg",
        "CG, PETScPreconditioner: type of preconditioner");
  check(isSolved(rksp) && isFaster(rksp, rj),
        "CG, PETScPreconditioner: resolution");
  // the preconditioner of PETSc with a Krylov solver of MFEM
  auto cg = generate(s, "CGSolver",
                     Parameters{{"RelativeTolerance", 1e-10},
                                {"MaximumNumberOfIterations", 500},
                                {"Preconditioner", petsc_gamg}});
  const auto rcg = solve(s, *(cg.linear_solver));
  const auto* const pc =
      dynamic_cast<const mfem::PetscPreconditioner*>(cg.preconditioner.get());
  check((pc != nullptr) && (getPCType(PC(*pc)) == "gamg"),
        "CGSolver, PETScPreconditioner: type of preconditioner");
  check(isSolved(rcg) && isFaster(rcg, rj),
        "CGSolver, PETScPreconditioner: resolution");
}  // end of checkPreconditioners

/*!
 * \brief checks that invalid options are refused
 * \param[in] s: system
 */
static void checkInvalidOptions(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  const auto refused = [&s](const Parameters& parameters,
                            const std::string& msg) {
    return getGenerationError(s, parameters).find(msg) != std::string::npos;
  };
  check(refused(Parameters{{"Options", Parameters{{"ksp_type", Parameters{}}}}},
                "invalid value of the option"),
        "invalid value of an option not refused");
  check(refused(Parameters{{"Options", Parameters{{"-ksp_type", "cg"}}}},
                "without the leading dash"),
        "option with a leading dash not refused");
  check(refused(Parameters{{"Options", Parameters{{"pc_type", "gamg"}}},
                           {"Preconditioner",
                            Parameters{{"Name", "HypreBoomerAMG"}}}},
                "can't be given to a solver with a preconditioner"),
        "option of the preconditioner of PETSc not refused");
}  // end of checkInvalidOptions

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto mesh =
      mfem::Mesh::MakeCartesian3D(16, 16, 16, mfem::Element::HEXAHEDRON);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    auto s = LaplaceSystem<true>(mesh, fec);
    checkOptions(s);
    checkPreconditioners(s);
    checkInvalidOptions(s);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
