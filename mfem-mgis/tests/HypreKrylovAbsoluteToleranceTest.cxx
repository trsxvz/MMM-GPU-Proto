/*!
 * \file   HypreKrylovAbsoluteToleranceTest.cxx
 * \brief  Test of the absolute tolerance of the Krylov solvers of hypre
 * built by the linear solver factory
 * \date   09/10/2026
 */

#include <cstdlib>
#include <iostream>
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pgridfunc.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

/*!
 * \return the number of iterations of a solver of the factory, or -1 on
 * failure
 * \param[in] fespace: finite element space
 * \param[in] A: matrix
 * \param[in] B: right hand side
 * \param[in] name: name of the solver
 * \param[in] params: parameters of the solver
 */
static int solve(mfem_mgis::FiniteElementSpace<true>& fespace,
                 mfem::HypreParMatrix& A,
                 const mfem::Vector& B,
                 const char* const name,
                 const mfem_mgis::Parameters& params) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto& factory = LinearSolverFactory<true>::getFactory();
  auto h = factory.generate(ctx, name, fespace, params);
  if (isInvalid(h)) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return -1;
  }
  auto X = mfem::Vector(B.Size());
  X = 0.0;
  h.linear_solver->SetOperator(A);
  h.linear_solver->Mult(B, X);
  const auto oit = getNumberOfIterationsAtConvergence(*(h.linear_solver));
  if (!oit.has_value()) {
    return -1;
  }
  return *oit;
}  // end of solve

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto smesh =
      mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON, 1, 1, 1);
  auto pmesh = mfem_mgis::Mesh<true>(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fespace = mfem_mgis::FiniteElementSpace<true>(&pmesh, &fec, 1);
  // Poisson problem, zero on the whole boundary
  auto one = mfem::ConstantCoefficient(1);
  auto a = mfem_mgis::BilinearForm<true>(&fespace);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  a.Assemble();
  auto b = mfem_mgis::LinearForm<true>(&fespace);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  b.Assemble();
  auto x = mfem_mgis::GridFunction<true>(&fespace);
  x = 0.0;
  auto boundaries = mfem::Array<int>(pmesh.bdr_attributes.Max());
  boundaries = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  auto A = mfem::HypreParMatrix{};
  auto B = mfem::Vector{};
  auto X = mfem::Vector{};
  a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
  //
  using mfem_mgis::Parameters;
  auto success = true;
  for (const auto* const name : {"HyprePCG", "HypreGMRES", "HypreFGMRES"}) {
    const auto relative = Parameters{{"Tolerance", 1e-12},
                                     {"MaximumNumberOfIterations", 1000},
                                     {"VerbosityLevel", 0}};
    // an absolute tolerance much larger than the initial residual: the
    // solver must stop at once
    auto absolute = relative;
    absolute.insert(mfem_mgis::may_throw,
                    {{"AbsoluteTolerance", 1e30}});
    const auto n = solve(fespace, A, B, name, relative);
    const auto na = solve(fespace, A, B, name, absolute);
    if (fespace.GetMyRank() == 0) {
      std::cout << name << ": " << n << " iterations, " << na
                << " iterations with a large absolute tolerance\n";
    }
    if ((n <= 1) || (na < 0) || (na >= n)) {
      success = false;
    }
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
