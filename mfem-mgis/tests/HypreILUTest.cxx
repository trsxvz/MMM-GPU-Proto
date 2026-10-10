/*!
 * \file   tests/HypreILUTest.cxx
 * \brief  Test of the incomplete factorisation of hypre (preconditioner
 * `HypreILU` of the factory of linear solvers): its level of fill is kept
 * after a new operator (`mfem::HypreILU::SetOperator` restores the default
 * options of MFEM, a level of fill of one), and, when hypre computes on a
 * GPU, where the factorisations with a level of fill greater than zero
 * fail, the level of fill is zero by default and a greater level is refused.
 * \date   10/10/2026
 */

#include <cstdlib>
#include <iostream>
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pgridfunc.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \return the numbers of iterations of two resolutions by the flexible
 * GMRES of hypre preconditioned by its incomplete factorisation, the
 * operator being given again before the second one, or an empty vector if
 * the solver is refused
 * \param[in] fespace: finite element space
 * \param[in] A: matrix
 * \param[in] B: right hand side
 * \param[in] ilu: options of the incomplete factorisation
 */
static std::vector<int> solve(mfem_mgis::FiniteElementSpace<true>& fespace,
                              mfem::HypreParMatrix& A,
                              const mfem::Vector& B,
                              const mfem_mgis::Parameters& ilu) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto& factory = LinearSolverFactory<true>::getFactory();
  auto h = factory.generate(
      ctx, "HypreFGMRES", fespace,
      Parameters{{"Tolerance", 1e-12},
                 {"MaximumNumberOfIterations", 1000},
                 {"VerbosityLevel", 0},
                 {"Preconditioner",
                  Parameters{{"Name", "HypreILU"}, {"Options", ilu}}}});
  if (isInvalid(h)) {
    if (fespace.GetMyRank() == 0) {
      std::cout << "refused: " << ctx.getErrorMessage() << '\n';
    }
    return {};
  }
  auto n = std::vector<int>{};
  for (int i = 0; i != 2; ++i) {
    auto X = mfem::Vector(B.Size());
    X = 0.0;
    h.linear_solver->SetOperator(A);
    h.linear_solver->Mult(B, X);
    const auto oit = getNumberOfIterationsAtConvergence(*(h.linear_solver));
    n.push_back(oit.has_value() ? *oit : -1);
  }
  return n;
}  // end of solve

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  // Laplacian on a cube, homogeneous Dirichlet conditions
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = Mesh<true>(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fespace = FiniteElementSpace<true>(&pmesh, &fec);
  auto boundaries = mfem::Array<int>(pmesh.bdr_attributes.Max());
  boundaries = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  auto one = mfem::ConstantCoefficient(1.);
  auto a = BilinearForm<true>(&fespace);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.Assemble();
  auto b = LinearForm<true>(&fespace);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  b.Assemble();
  auto x = GridFunction<true>(&fespace);
  x = 0.0;
  auto A = mfem::HypreParMatrix{};
  auto B = mfem::Vector{};
  auto X = mfem::Vector{};
  a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
  const auto verbose = fespace.GetMyRank() == 0;
  auto print = [verbose](const char* const name, const std::vector<int>& n) {
    if (verbose) {
      std::cout << name << ':';
      for (const auto i : n) {
        std::cout << ' ' << i;
      }
      std::cout << '\n';
    }
  };
  const auto n0 = solve(fespace, A, B, {{"HypreILULevelOfFill", 0}});
  print("level of fill 0", n0);
  const auto nd = solve(fespace, A, B, {});
  print("default level of fill", nd);
  // the level of fill is kept after a new operator
  check((n0.size() == 2) && (n0[0] > 0) && (n0[1] == n0[0]),
        "level of fill 0: numbers of iterations");
  check((nd.size() == 2) && (nd[0] > 0) && (nd[1] == nd[0]),
        "default level of fill: numbers of iterations");
  const auto n2 = solve(fespace, A, B, {{"HypreILULevelOfFill", 2}});
  print("level of fill 2", n2);
  if (isHypreUsingGPU()) {
    // ILU(0) by default, greater levels of fill refused
    check(nd == n0, "the default level of fill is not zero on a GPU");
    check(n2.empty(), "a level of fill of 2 is not refused on a GPU");
  } else {
    check((n2.size() == 2) && (n2[0] > 0) && (n2[1] == n2[0]),
          "level of fill 2: numbers of iterations");
    check((n2.size() == 2) && (n2[0] < n0[0]),
          "the level of fill 2 does not reduce the number of iterations");
  }
  if (verbose) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
