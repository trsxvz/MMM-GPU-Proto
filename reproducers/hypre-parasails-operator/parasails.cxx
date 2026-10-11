// ParaSails preconditioner of hypre (mfem::HypreParaSails) preconditioning
// the GMRES of hypre, for two successive operators, the first one being
// destroyed before the second one is given, as a non linear solver does
// when the tangent operator is assembled again.
//
// mfem::HypreParaSails::SetOperator reads the communicator of the previous
// operator to create a new preconditioner of hypre: the first operator
// being destroyed, the communicator read is undefined.
//
// usage: parasails <device> <mode>
// - device: cpu or cuda;
// - mode: mfem for mfem::HypreParaSails, kept for mfem::HypreParaSails with
//   the first operator kept alive, fixed for a preconditioner which points
//   to the new operator before calling mfem::HypreParaSails::SetOperator.
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string_view>
#include "mfem/general/communication.hpp"
#include "mfem/general/device.hpp"
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/lininteg.hpp"
#include "mfem/linalg/hypre.hpp"

// ParaSails preconditioner reading the communicator of the new operator
struct FixedParaSails : mfem::HypreParaSails {
  using mfem::HypreParaSails::HypreParaSails;
  void SetOperator(const mfem::Operator& op) override {
    if (A != nullptr) {
      A = dynamic_cast<const mfem::HypreParMatrix*>(&op);
    }
    mfem::HypreParaSails::SetOperator(op);
  }
};

int main(int argc, char** argv) {
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  mfem::Device device(argc > 1 ? argv[1] : "cpu");
  const auto mode = std::string_view{argc > 2 ? argv[2] : "mfem"};
  if ((mode != "mfem") && (mode != "kept") && (mode != "fixed")) {
    std::cerr << "invalid mode '" << mode << "'\n";
    return EXIT_FAILURE;
  }
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fes = mfem::ParFiniteElementSpace(&pmesh, &fec);
  auto ess_bdr = mfem::Array<int>(pmesh.bdr_attributes.Max());
  ess_bdr = 1;
  auto ess = mfem::Array<int>{};
  fes.GetEssentialTrueDofs(ess_bdr, ess);
  auto one = mfem::ConstantCoefficient(1.);
  auto b = mfem::ParLinearForm(&fes);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  b.Assemble();
  // right-hand side, zero on the boundary
  auto B = std::unique_ptr<mfem::HypreParVector>(b.ParallelAssemble());
  B->SetSubVector(ess, 0.);
  auto X = mfem::HypreParVector(*B);
  // operator of a diffusion problem whose conductivity is k, owned by the
  // caller
  auto assemble = [&](const double k) {
    auto c = mfem::ConstantCoefficient(k);
    auto a = mfem::ParBilinearForm(&fes);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(c));
    a.Assemble();
    a.Finalize();
    auto A = std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
    delete A->EliminateRowsCols(ess);
    return A;
  };
  auto mfem_parasails = mfem::HypreParaSails(MPI_COMM_WORLD);
  auto fixed_parasails = FixedParaSails(MPI_COMM_WORLD);
  mfem::HypreSolver& prec = (mode == "fixed")
                                ? static_cast<mfem::HypreSolver&>(fixed_parasails)
                                : mfem_parasails;
  auto solver = mfem::HypreGMRES(MPI_COMM_WORLD);
  solver.SetTol(1e-12);
  solver.SetMaxIter(500);
  solver.SetPreconditioner(prec);
  auto A = std::unique_ptr<mfem::HypreParMatrix>{};
  auto kept = std::unique_ptr<mfem::HypreParMatrix>{};
  for (int i = 0; i != 2; ++i) {
    // the new operator is assembled before the destruction of the previous
    // one (kept alive in the mode kept), so that nothing is allocated in
    // the memory freed before the new operator is given
    auto next = assemble(1. + i);
    if (mode == "kept") {
      kept = std::move(A);
    }
    A = std::move(next);
    X = 0.;
    solver.SetOperator(*A);
    solver.Mult(*B, X);
    int n = 0;
    solver.GetNumIterations(n);
    if (mfem::Mpi::Root()) {
      std::cout << mode << " resolution " << i << ": " << n << " iterations\n";
    }
  }
  return EXIT_SUCCESS;
}
