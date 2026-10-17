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
 * - on the linear elasticity of a beam, the rigid body modes given to GAMG
 *   (unknowns ordered by vector dimension) reduce the number of iterations,
 *   for the Krylov solvers of PETSc and of MFEM, and are not given with the
 *   unknowns ordered by nodes; when the preconditioner of PETSc
 *   preconditions the Krylov solver of PETSc, its own option applies;
 * - a preconditioner of mfem-mgis is used as a shell preconditioner, which
 *   receives the operators of different sizes given to the solver;
 * - the preconditioner of PETSc preconditions the Krylov solvers of PETSc
 *   and of MFEM;
 * - invalid options are refused: an invalid value, a leading dash, and the
 *   options of the preconditioner of PETSc given to a solver with a
 *   preconditioner of mfem-mgis;
 * - the external packages of PETSc, if PETSc is built with them: the sparse
 *   direct solvers MUMPS, PaStiX, SuperLU_DIST and STRUMPACK (LU
 *   factorization of PETSc), and the Krylov solvers of HPDDM (GCRODR,
 *   recycling a subspace between two resolutions);
 * - the matrix given to the Krylov solver of PETSc is the matrix of hypre
 *   (products by a vector compared), in the AIJ format of the GPU when MFEM
 *   computes on a GPU for which PETSc is built (the solvers of PETSc then
 *   compute on the GPU), of the host otherwise;
 * - when MFEM computes on the host, the solvers of PETSc do not initialize
 *   the GPU.
 * \date   15/10/2026
 */

#include <cmath>
#include <array>
#include <algorithm>
#include <string>
#include <string_view>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "mfem.hpp"
#include "petscksp.h"
#include "petscdevice.h"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "LaplaceSystem.hxx"
#include "ElasticitySystem.hxx"

using mfem_mgis::tests::ElasticitySystem;
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
template <typename System>
static mfem_mgis::LinearSolverHandler generate(
    System& s,
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
template <typename System>
static Resolution solve(System& s, mfem::Solver& ls) {
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

/*!
 * \brief check that the matrix given to the Krylov solver of PETSc is the
 * matrix of hypre: their products by a vector are compared
 * \param[in] s: system
 * \param[in] name: name of the system
 */
template <typename System>
static void checkMatrixCopy(System& s, const std::string& name) {
  using namespace mfem_mgis;
  auto petsc = [](const PetscErrorCode e, const char* const what) {
    if (e != 0) {
      throw std::runtime_error(std::string{what} + " failed");
    }
  };
  auto ls = generate(s, "PETScLinearSolver", Parameters{});
  ls.linear_solver->SetOperator(s.A);
  auto a = Mat{};
  petsc(KSPGetOperators(getKSP(ls), &a, nullptr), "KSPGetOperators");
  // the solvers of PETSc compute on the GPU when MFEM computes on a GPU for
  // which PETSc is built
#if defined(MFEM_USE_CUDA) && defined(PETSC_HAVE_CUDA)
  const auto device = mfem::Device::Allows(mfem::Backend::CUDA_MASK);
#elif defined(MFEM_USE_HIP) && defined(PETSC_HAVE_HIP)
  const auto device = mfem::Device::Allows(mfem::Backend::HIP_MASK);
#else
  const auto device = false;
#endif
  auto type = MatType{};
  petsc(MatGetType(a, &type), "MatGetType");
  const auto t = std::string_view{type};
  check((t.ends_with("cusparse") || t.ends_with("hipsparse")) == device,
        name + ": the matrix given to PETSc is of type " + std::string{t});
  const auto n = s.B.Size();
  auto x = mfem::Vector(n);
  x.UseDevice(true);
  x.Randomize(1 + rank);
  auto y = mfem::Vector(n);
  y.UseDevice(true);
  s.A.Mult(x, y);
  auto px = Vec{};
  auto py = Vec{};
  petsc(MatCreateVecs(a, &px, &py), "MatCreateVecs");
  auto* v = static_cast<PetscScalar*>(nullptr);
  petsc(VecGetArrayWrite(px, &v), "VecGetArrayWrite");
  const auto* const hx = x.HostRead();
  std::copy(hx, hx + n, v);
  petsc(VecRestoreArrayWrite(px, &v), "VecRestoreArrayWrite");
  petsc(MatMult(a, px, py), "MatMult");
  // maximum difference and maximum value of the products
  auto e = std::array<double, 2>{0, 0};
  const auto* w = static_cast<const PetscScalar*>(nullptr);
  petsc(VecGetArrayRead(py, &w), "VecGetArrayRead");
  const auto* const hy = y.HostRead();
  for (int i = 0; i != n; ++i) {
    e[0] = std::max(e[0], std::abs(static_cast<double>(w[i] - hy[i])));
    e[1] = std::max(e[1], std::abs(static_cast<double>(hy[i])));
  }
  petsc(VecRestoreArrayRead(py, &w), "VecRestoreArrayRead");
  petsc(VecDestroy(&px), "VecDestroy");
  petsc(VecDestroy(&py), "VecDestroy");
  MPI_Allreduce(MPI_IN_PLACE, e.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  // the sums of the products differ by their order (hypre sums the
  // diagonal entry first, PETSc sorts the columns)
  check((e[1] > 0) && (e[0] <= 1e-12 * e[1]),
        name + ": the matrix given to PETSc differs from the matrix of hypre");
}  // end of checkMatrixCopy

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
 * \return true if PETSc is built with an external package
 * \param[in] n: name of the package
 */
static bool hasExternalPackage(const char* const n) {
  auto has = PetscBool{PETSC_FALSE};
  if (PetscHasExternalPackage(n, &has) != PETSC_SUCCESS) {
    return false;
  }
  return has == PETSC_TRUE;
}  // end of hasExternalPackage

/*!
 * \brief checks the external packages of PETSc which PETSc is built with:
 * sparse direct solvers and Krylov solvers of HPDDM
 * \param[in] s: system
 */
static void checkExternalPackages(LaplaceSystem<true>& s) {
  using namespace mfem_mgis;
  for (const auto* const n : {"mumps", "pastix", "superlu_dist", "strumpack"}) {
    if (!hasExternalPackage(n)) {
      if (rank == 0) {
        std::cout << "PETSc, " << n << ": skipped (not available)\n";
      }
      continue;
    }
    auto direct = generate(
        s, "PETScLinearSolver",
        Parameters{{"Options", Parameters{{"ksp_type", "preonly"},
                                          {"pc_type", "lu"},
                                          {"pc_factor_mat_solver_type", n}}}});
    const auto r = solve(s, *(direct.linear_solver));
    if (rank == 0) {
      std::cout << "PETSc, LU factorization of " << n << ": residual "
                << r.residual << '\n';
    }
    check(r.converged && (r.residual < 1e-10),
          std::string{"PETSc, LU factorization of "} + n + ": resolution");
  }
  if (hasExternalPackage("hpddm")) {
    auto gcrodr =
        generate(s, "PETScLinearSolver",
                 Parameters{{"Options", Parameters{{"ksp_type", "hpddm"},
                                                   {"ksp_hpddm_type", "gcrodr"},
                                                   {"ksp_hpddm_recycle", 10},
                                                   {"pc_type", "jacobi"},
                                                   {"ksp_rtol", 1e-10}}}});
    const auto r1 = solve(s, *(gcrodr.linear_solver));
    const auto r2 = solve(s, *(gcrodr.linear_solver));
    if (rank == 0) {
      std::cout << "PETSc, GCRODR of HPDDM: " << r1.iterations << " then "
                << r2.iterations << " iterations\n";
    }
    check(isSolved(r1) && isSolved(r2) && (r2.iterations < r1.iterations),
          "PETSc, GCRODR of HPDDM: no recycling of the subspace");
  } else if (rank == 0) {
    std::cout << "PETSc, hpddm: skipped (not available)\n";
  }
}  // end of checkExternalPackages

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
  // the shell preconditioner applied once (preonly), to the vectors of PETSc
  // wrapping the right hand side and the solution given to the solver: the
  // result is the one of the preconditioner applied by itself
  {
    auto once = generate(
        s, "PETScLinearSolver",
        Parameters{
            {"Options", Parameters{{"ksp_type", "preonly"}}},
            {"Preconditioner",
             Parameters{{"Name", "HypreBoomerAMG"},
                        {"Options", Parameters{{"VerbosityLevel", 0}}}}}});
    const auto n = s.B.Size();
    auto x1 = mfem::Vector(n);
    x1.UseDevice(true);
    x1 = 0.;
    once.linear_solver->SetOperator(s.A);
    once.linear_solver->Mult(s.B, x1);
    auto x2 = mfem::Vector(n);
    x2.UseDevice(true);
    once.preconditioner->Mult(s.B, x2);
    // `mfem::ParNormlp` reads the host array without synchronizing it in
    // MFEM 4.10 (mfem/mfem#5586)
    const auto norm = std::sqrt(mfem::InnerProduct(MPI_COMM_WORLD, x2, x2));
    x2 -= x1;
    const auto difference =
        std::sqrt(mfem::InnerProduct(MPI_COMM_WORLD, x2, x2));
    check((norm > 0) && (difference <= 1e-10 * norm),
          "preonly, BoomerAMG: the shell preconditioner gives another result "
          "than the preconditioner");
  }
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
 * \brief checks the rigid body modes given to the algebraic multigrid of
 * PETSc (GAMG) on the linear elasticity: with the unknowns ordered by vector
 * dimension, they reduce the number of iterations compared to the same
 * resolution without them (option `RigidBodyModes`), whereas they are not
 * given with the unknowns ordered by nodes (same iterations with and
 * without the option)
 * \param[in] e1: system, unknowns ordered by vector dimension
 * \param[in] e2: system, unknowns ordered by nodes
 */
static void checkRigidBodyModes(ElasticitySystem<true>& e1,
                                ElasticitySystem<true>& e2) {
  using namespace mfem_mgis;
  // the conjugate gradients stop on the norm of the preconditioned
  // residual, the residual being then about 1e-7
  const auto solved = [](const Resolution& r) {
    return r.converged && (r.residual < 1e-5);
  };
  const auto parameters = [](const bool rbm) {
    return Parameters{
        {"Options", Parameters{{"ksp_type", "cg"}, {"pc_type", "gamg"}}},
        {"RelativeTolerance", 1e-8},
        {"MaximumNumberOfIterations", 1000},
        {"RigidBodyModes", rbm}};
  };
  auto with = generate(e1, "PETScLinearSolver", parameters(true));
  const auto r1 = solve(e1, *(with.linear_solver));
  auto without = generate(e1, "PETScLinearSolver", parameters(false));
  const auto r2 = solve(e1, *(without.linear_solver));
  check(solved(r1) && solved(r2) && isFaster(r1, r2),
        "CG, GAMG: the rigid body modes do not reduce the number of "
        "iterations");
  auto with2 = generate(e2, "PETScLinearSolver", parameters(true));
  const auto r3 = solve(e2, *(with2.linear_solver));
  auto without2 = generate(e2, "PETScLinearSolver", parameters(false));
  const auto r4 = solve(e2, *(without2.linear_solver));
  check(solved(r3) && (r3.iterations == r4.iterations),
        "CG, GAMG: rigid body modes given with the unknowns ordered by "
        "nodes");
  // the preconditioner of PETSc with the Krylov solver of MFEM
  auto cg = generate(
      e1, "CGSolver",
      Parameters{
          {"RelativeTolerance", 1e-8},
          {"MaximumNumberOfIterations", 1000},
          {"Preconditioner",
           Parameters{
               {"Name", "PETScPreconditioner"},
               {"Options",
                Parameters{{"Options", Parameters{{"pc_type", "gamg"}}}}}}}});
  const auto r5 = solve(e1, *(cg.linear_solver));
  check(solved(r5) && isFaster(r5, r2),
        "CGSolver, PETScPreconditioner (GAMG): the rigid body modes do not "
        "reduce the number of iterations");
  // the preconditioner of PETSc with the Krylov solver of PETSc, which
  // share their object of PETSc: the option of the preconditioner applies
  const auto nested = [](const bool solver_rbm, const bool pc_rbm) {
    return Parameters{
        {"Options", Parameters{{"ksp_type", "cg"}}},
        {"RelativeTolerance", 1e-8},
        {"MaximumNumberOfIterations", 1000},
        {"RigidBodyModes", solver_rbm},
        {"Preconditioner",
         Parameters{{"Name", "PETScPreconditioner"},
                    {"Options",
                     Parameters{{"Options", Parameters{{"pc_type", "gamg"}}},
                                {"RigidBodyModes", pc_rbm}}}}}};
  };
  auto n1 = generate(e1, "PETScLinearSolver", nested(false, true));
  const auto r6 = solve(e1, *(n1.linear_solver));
  check(solved(r6) && (r6.iterations == r1.iterations),
        "CG, PETScPreconditioner (GAMG): the rigid body modes of the "
        "preconditioner are not used");
  auto n2 = generate(e1, "PETScLinearSolver", nested(true, false));
  const auto r7 = solve(e1, *(n2.linear_solver));
  check(solved(r7) && (r7.iterations == r2.iterations),
        "CG, PETScPreconditioner (GAMG): rigid body modes given despite the "
        "option of the preconditioner");
}  // end of checkRigidBodyModes

/*!
 * \brief checks the fields of the field split preconditioner of PETSc, one
 * per component of the unknowns (option `Fields` set to `Components`), for
 * both orderings of the unknowns: an additive field split (block Jacobi by
 * component) with GAMG on each component solves the linear elasticity, with
 * close numbers of iterations for both orderings, and with three fields
 * \param[in] e1: system, unknowns ordered by vector dimension
 * \param[in] e2: system, unknowns ordered by nodes
 */
static void checkFieldSplit(ElasticitySystem<true>& e1,
                            ElasticitySystem<true>& e2) {
  using namespace mfem_mgis;
  // the conjugate gradients stop on the norm of the preconditioned
  // residual
  const auto solved = [](const Resolution& r) {
    return r.converged && (r.residual < 1e-5);
  };
  const auto options = Parameters{{"pc_type", "fieldsplit"},
                                  {"pc_fieldsplit_type", "additive"},
                                  {"fieldsplit_0_pc_type", "gamg"},
                                  {"fieldsplit_1_pc_type", "gamg"},
                                  {"fieldsplit_2_pc_type", "gamg"}};
  auto ksp_options = options;
  ksp_options.insert(may_throw, Parameters{{"ksp_type", "cg"}});
  const auto parameters = Parameters{{"Options", ksp_options},
                                     {"RelativeTolerance", 1e-8},
                                     {"MaximumNumberOfIterations", 1000},
                                     {"Fields", "Components"}};
  const auto number_of_fields = [](const LinearSolverHandler& ls) {
    auto pc = PC{};
    auto n = PetscInt{};
    auto* sub = static_cast<KSP*>(nullptr);
    if ((KSPGetPC(getKSP(ls), &pc) != 0) ||
        (PCFieldSplitGetSubKSP(pc, &n, &sub) != 0)) {
      return PetscInt{-1};
    }
    PetscFree(sub);
    return n;
  };
  auto s1 = generate(e1, "PETScLinearSolver", parameters);
  const auto r1 = solve(e1, *(s1.linear_solver));
  check(solved(r1) && (number_of_fields(s1) == 3),
        "CG, field split by component: resolution (byVDIM)");
  auto s2 = generate(e2, "PETScLinearSolver", parameters);
  const auto r2 = solve(e2, *(s2.linear_solver));
  check(solved(r2) && (number_of_fields(s2) == 3),
        "CG, field split by component: resolution (byNODES)");
  check(
      std::abs(r1.iterations - r2.iterations) <= std::max(3, r1.iterations / 5),
      "CG, field split by component: the ordering of the unknowns changes "
      "the number of iterations");
  // the preconditioner of PETSc with the Krylov solver of MFEM
  auto cg = generate(
      e2, "CGSolver",
      Parameters{
          {"RelativeTolerance", 1e-8},
          {"MaximumNumberOfIterations", 1000},
          {"Preconditioner",
           Parameters{{"Name", "PETScPreconditioner"},
                      {"Options", Parameters{{"Options", options},
                                             {"Fields", "Components"}}}}}});
  const auto r3 = solve(e2, *(cg.linear_solver));
  check(solved(r3), "CGSolver, PETScPreconditioner (field split): resolution");
}  // end of checkFieldSplit

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
  check(refused(Parameters{{"Fields", "Rows"}},
                "invalid value 'Rows' of the option 'Fields'"),
        "invalid value of the option 'Fields' not refused");
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
    checkExternalPackages(s);
    // linear elasticity of a beam
    auto beam = mfem::Mesh::MakeCartesian3D(16, 4, 4, mfem::Element::HEXAHEDRON,
                                            4., 1., 1.);
    auto e1 = ElasticitySystem<true>(beam, fec, mfem::Ordering::byVDIM);
    auto e2 = ElasticitySystem<true>(beam, fec, mfem::Ordering::byNODES);
    checkRigidBodyModes(e1, e2);
    checkFieldSplit(e1, e2);
    checkMatrixCopy(s, "Laplace equation");
    checkMatrixCopy(e1, "linear elasticity");
    // the solvers of PETSc compute on the host: when MFEM does too, they
    // must not initialize the GPU, which opens a context per process
    if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK)) {
      check((!PetscDeviceInitialized(PETSC_DEVICE_CUDA)) &&
                (!PetscDeviceInitialized(PETSC_DEVICE_HIP)),
            "the solvers of PETSc initialized the GPU");
    }
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
