// Incomplete factorisation of hypre (ILU) preconditioning its flexible
// GMRES, on the CPU or on a GPU (hypre built with CUDA): three resolutions,
// the operator being given again before each one.
//
// usage: ilu <device> <mode> <level of fill>
// - device: cpu or cuda;
// - mode: 0 for mfem::HypreILU, which creates a new preconditioner of hypre
//   with the default options of MFEM (level of fill of one) at each new
//   operator, 1 for a preconditioner kept from one operator to the other
//   (setup again by hypre), whose level of fill is kept;
// - level of fill (mode 1 only).
#include <cstdlib>
#include <iostream>
#include "mfem/general/communication.hpp"
#include "mfem/general/device.hpp"
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pgridfunc.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/lininteg.hpp"
#include "mfem/linalg/hypre.hpp"
struct KeptILU : mfem::HypreSolver {
  HYPRE_Solver ilu = nullptr;
  explicit KeptILU(const int fill) { HYPRE_ILUCreate(&ilu); HYPRE_ILUSetMaxIter(ilu, 1); HYPRE_ILUSetTol(ilu, 0.); HYPRE_ILUSetLevelOfFill(ilu, fill); }
  ~KeptILU() override { HYPRE_ILUDestroy(ilu); }
  operator HYPRE_Solver() const override { return ilu; }
  void SetOperator(const mfem::Operator& op) override {
    A = const_cast<mfem::HypreParMatrix*>(dynamic_cast<const mfem::HypreParMatrix*>(&op));
    height = A->Height(); width = A->Width();
    setup_called = 0;
    delete X; delete B; B = X = nullptr;
    auxB.Delete(); auxB.Reset(); auxX.Delete(); auxX.Reset();
  }
  HYPRE_PtrToParSolverFcn SetupFcn() const override { return (HYPRE_PtrToParSolverFcn)HYPRE_ILUSetup; }
  HYPRE_PtrToParSolverFcn SolveFcn() const override { return (HYPRE_PtrToParSolverFcn)HYPRE_ILUSolve; }
};
int main(int argc, char** argv) {
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  mfem::Device device(argc > 1 ? argv[1] : "cpu");
  const int kept = argc > 2 ? std::atoi(argv[2]) : 0;
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fes = mfem::ParFiniteElementSpace(&pmesh, &fec);
  auto ess_bdr = mfem::Array<int>(pmesh.bdr_attributes.Max());
  ess_bdr = 1;
  auto ess = mfem::Array<int>{};
  fes.GetEssentialTrueDofs(ess_bdr, ess);
  auto one = mfem::ConstantCoefficient(1.);
  auto a = mfem::ParBilinearForm(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.Assemble();
  auto b = mfem::ParLinearForm(&fes);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
  b.Assemble();
  auto x = mfem::ParGridFunction(&fes);
  x = 0.;
  auto A = mfem::HypreParMatrix{};
  auto B = mfem::Vector{};
  auto X = mfem::Vector{};
  a.FormLinearSystem(ess, x, b, A, X, B);
  auto mfem_ilu = mfem::HypreILU();
  auto kept_ilu = KeptILU(argc > 3 ? std::atoi(argv[3]) : 0);
  mfem::HypreSolver& prec = kept ? static_cast<mfem::HypreSolver&>(kept_ilu) : mfem_ilu;
  auto solver = mfem::HypreFGMRES(MPI_COMM_WORLD);
  solver.SetTol(1e-12);
  solver.SetMaxIter(500);
  solver.SetPreconditioner(prec);
  for (int i = 0; i != 3; ++i) {
    X = 0.;
    solver.SetOperator(A);
    solver.Mult(B, X);
    int n = 0;
    solver.GetNumIterations(n);
    if (mfem::Mpi::Root()) {
      std::cout << (kept ? "kept" : "mfem") << " resolution " << i << ": " << n << " iterations\n";
    }
  }
  return EXIT_SUCCESS;
}
