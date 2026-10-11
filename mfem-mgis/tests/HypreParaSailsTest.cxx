/*!
 * \file   tests/HypreParaSailsTest.cxx
 * \brief  Test of the ParaSails preconditioner of hypre (preconditioner
 * `HypreParaSails` of the factory of linear solvers): two resolutions with
 * two operators, the first one being destroyed before the second one is
 * given, as the non linear solvers do (`mfem::HypreParaSails::SetOperator`
 * reads the communicator of the previous operator). When hypre computes on
 * a GPU, where ParaSails is not supported, the preconditioner is refused.
 * \date   11/10/2026
 */

#include <cstdlib>
#include <memory>
#include <string>
#include <iostream>
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "mfem/fem/plinearform.hpp"
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
 * \return the matrix of a diffusion problem on a cube, with homogeneous
 * Dirichlet conditions, owned by the caller
 * \param[in] fespace: finite element space
 * \param[in] ess_tdofs: degrees of freedom of the Dirichlet conditions
 * \param[in] k: conductivity
 */
static std::unique_ptr<mfem::HypreParMatrix> assemble(
    mfem_mgis::FiniteElementSpace<true>& fespace,
    const mfem::Array<int>& ess_tdofs,
    const mfem_mgis::real k) {
  auto c = mfem::ConstantCoefficient(k);
  auto a = mfem_mgis::BilinearForm<true>(&fespace);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(c));
  a.Assemble();
  a.Finalize();
  auto A = std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  delete A->EliminateRowsCols(ess_tdofs);
  return A;
}  // end of assemble

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = Mesh<true>(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fespace = FiniteElementSpace<true>(&pmesh, &fec);
  auto boundaries = mfem::Array<int>(pmesh.bdr_attributes.Max());
  boundaries = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  auto one = mfem::ConstantCoefficient(1.);
  auto b = LinearForm<true>(&fespace);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  b.Assemble();
  auto B = std::unique_ptr<mfem::HypreParVector>(b.ParallelAssemble());
  B->SetSubVector(ess_tdofs, 0.);
  const auto verbose = fespace.GetMyRank() == 0;
  auto ctx = Context{};
  auto& factory = LinearSolverFactory<true>::getFactory();
  auto h = factory.generate(
      ctx, "HypreGMRES", fespace,
      Parameters{{"Tolerance", 1e-12},
                 {"MaximumNumberOfIterations", 1000},
                 {"VerbosityLevel", 0},
                 {"Preconditioner",
                  Parameters{{"Name", "HypreParaSails"},
                             {"Options", Parameters{{"VerbosityLevel", 0}}}}}});
  if (isHypreUsingGPU()) {
    check(isInvalid(h), "ParaSails is not refused on a GPU");
    if (verbose && isInvalid(h)) {
      std::cout << "refused: " << ctx.getErrorMessage() << '\n';
    }
  } else {
    check(!isInvalid(h), "ParaSails is refused");
    if (!isInvalid(h)) {
      auto A = std::unique_ptr<mfem::HypreParMatrix>{};
      for (int i = 0; i != 2; ++i) {
        // the new operator is assembled before the destruction of the
        // previous one, so that nothing is allocated in the memory freed
        // before the new operator is given
        A = assemble(fespace, ess_tdofs, 1 + i);
        auto X = mfem::Vector(B->Size());
        X = 0.0;
        h.linear_solver->SetOperator(*A);
        h.linear_solver->Mult(*B, X);
        const auto oit = getNumberOfIterationsAtConvergence(*(h.linear_solver));
        check(hasConverged(*(h.linear_solver)) && oit.has_value() && (*oit > 0),
              "resolution " + std::to_string(i) + ": no convergence");
        if (verbose) {
          std::cout << "resolution " << i << ": "
                    << (oit.has_value() ? *oit : -1) << " iterations\n";
        }
      }
    }
  }
  if (verbose) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
