// mfem::HypreParaSails preconditioning the conjugate gradient of a
// mfem::NewtonSolver applied to a mfem::ParNonlinearForm: only classes of
// MFEM are used.
//
// At each Newton iteration, ParNonlinearForm::GetGradient destroys the
// previous gradient before assembling the new one, then NewtonSolver gives
// the new gradient to the linear solver, which gives it to the
// preconditioner. HypreParaSails::SetOperator then reads the communicator of
// the destroyed gradient.
//
// The problem is linear (diffusion), but the loose tolerance of the
// conjugate gradient makes the Newton solver iterate several times.
//
// usage: newton <mode>
// - mode: mfem for mfem::HypreParaSails, fixed for a preconditioner which
//   reads the communicator of the new operator, as mfem::HypreEuclid does.
#include <cstdlib>
#include <iostream>
#include <string_view>
#include "mfem/general/communication.hpp"
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/fem/pnonlinearform.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/lininteg.hpp"
#include "mfem/linalg/hypre.hpp"
#include "mfem/linalg/solvers.hpp"

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
  const auto mode = std::string_view{argc > 1 ? argv[1] : "mfem"};
  if ((mode != "mfem") && (mode != "fixed")) {
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
  auto B = mfem::Vector(fes.GetTrueVSize());
  b.ParallelAssemble(B);
  B.SetSubVector(ess, 0.);
  auto X = mfem::Vector(fes.GetTrueVSize());
  X = 0.;
  auto form = mfem::ParNonlinearForm(&fes);
  form.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  form.SetEssentialTrueDofs(ess);
  auto mfem_parasails = mfem::HypreParaSails(MPI_COMM_WORLD);
  auto fixed_parasails = FixedParaSails(MPI_COMM_WORLD);
  mfem::Solver& prec = (mode == "fixed")
                           ? static_cast<mfem::Solver&>(fixed_parasails)
                           : mfem_parasails;
  auto cg = mfem::CGSolver(MPI_COMM_WORLD);
  cg.SetRelTol(1e-2);
  cg.SetMaxIter(500);
  cg.SetPreconditioner(prec);
  auto newton = mfem::NewtonSolver(MPI_COMM_WORLD);
  newton.SetSolver(cg);
  newton.SetOperator(form);
  newton.SetRelTol(1e-10);
  newton.SetMaxIter(20);
  newton.Mult(B, X);
  if (mfem::Mpi::Root()) {
    std::cout << mode << " newton: " << newton.GetNumIterations()
              << " iterations, "
              << (newton.GetConverged() ? "converged" : "not converged")
              << '\n';
  }
  return newton.GetConverged() ? EXIT_SUCCESS : EXIT_FAILURE;
}
