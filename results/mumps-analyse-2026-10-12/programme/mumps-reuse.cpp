// Cost of the symbolic analysis of MUMPS for repeated factorizations of a
// matrix with the same pattern (elasticity on the mesh of ex8 of
// mfem-mgis-examples, linear elements, byNODES ordering like ex8).
//
// Usage: mumps-reuse mesh.msh [number of factorizations] [ordering] [BLR
// tolerance]. ordering: automatic, amd, amf, pord, metis, parmetis, scotch
// or ptscotch. BLR tolerance: 0 for none.
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <memory>
#include "mfem.hpp"

static auto ordering = mfem::MUMPSSolver::ReorderingStrategy::AUTOMATIC;
static double blr_tolerance = 0;

static double factorize(mfem::HypreParMatrix& A, const bool reuse,
                        const int n, const mfem::Vector& b, mfem::Vector& x) {
  mfem::StopWatch sw;
  auto s = mfem::MUMPSSolver(A.GetComm());
  s.SetPrintLevel(0);
  s.SetMatrixSymType(mfem::MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
  s.SetReorderingReuse(reuse);
  s.SetReorderingStrategy(ordering);
  if (blr_tolerance > 0) {
    s.SetBLRTol(blr_tolerance);
  }
  auto total = 0.;
  auto first = 0.;
  for (int i = 0; i != n; ++i) {
    // values change, pattern kept
    A *= (1 + 0.01 * i) / (i == 0 ? 1 : (1 + 0.01 * (i - 1)));
    sw.Clear();
    sw.Start();
    s.SetOperator(A);
    sw.Stop();
    if (i == 0) {
      first = sw.RealTime();
    } else {
      total += sw.RealTime();
    }
    s.Mult(b, x);
  }
  if (mfem::Mpi::Root()) {
    std::cout << (reuse ? "reuse   " : "no reuse") << ": first " << first
              << " s, next ones " << total / (n - 1) << " s each\n";
  }
  return total;
}

int main(int argc, char** argv) {
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " mesh [n]\n";
    return EXIT_FAILURE;
  }
  const int n = argc > 2 ? std::atoi(argv[2]) : 10;
  if (argc > 3) {
    using R = mfem::MUMPSSolver::ReorderingStrategy;
    const auto* const o = argv[3];
    ordering = std::strcmp(o, "amd") == 0        ? R::AMD
               : std::strcmp(o, "amf") == 0      ? R::AMF
               : std::strcmp(o, "pord") == 0     ? R::PORD
               : std::strcmp(o, "metis") == 0    ? R::METIS
               : std::strcmp(o, "parmetis") == 0 ? R::PARMETIS
               : std::strcmp(o, "scotch") == 0   ? R::SCOTCH
               : std::strcmp(o, "ptscotch") == 0 ? R::PTSCOTCH
                                                 : R::AUTOMATIC;
  }
  if (argc > 4) {
    blr_tolerance = std::atof(argv[4]);
  }
  auto smesh = mfem::Mesh(argv[1], 1, 1);
  auto mesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  smesh.Clear();
  const auto dim = mesh.Dimension();
  auto fec = mfem::H1_FECollection(1, dim);
  auto fes = mfem::ParFiniteElementSpace(&mesh, &fec, dim,
                                         mfem::Ordering::byNODES);
  const auto size = fes.GlobalTrueVSize();
  if (mfem::Mpi::Root()) {
    std::cout << "unknowns: " << size << '\n';
  }
  // clamp all the boundaries (the mesh of ex8 has several parts)
  auto ess_bdr = mfem::Array<int>(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  auto lambda = mfem::ConstantCoefficient(1.);
  auto mu = mfem::ConstantCoefficient(1.);
  auto a = mfem::ParBilinearForm(&fes);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
  a.Assemble();
  auto ones = mfem::VectorConstantCoefficient(mfem::Vector({1., 1., 1.}));
  auto f = mfem::ParLinearForm(&fes);
  f.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(ones));
  f.Assemble();
  auto u = mfem::ParGridFunction(&fes);
  u = 0.;
  auto A = mfem::HypreParMatrix{};
  auto X = mfem::Vector{};
  auto B = mfem::Vector{};
  a.FormLinearSystem(ess_tdofs, u, f, A, X, B);
  factorize(A, false, n, B, X);
  factorize(A, true, n, B, X);
  // error of the last solution
  auto R = mfem::Vector(B.Size());
  A.Mult(X, R);
  R -= B;
  const auto r = std::sqrt(mfem::InnerProduct(A.GetComm(), R, R) /
                           mfem::InnerProduct(A.GetComm(), B, B));
  if (mfem::Mpi::Root()) {
    std::cout << "relative residual of the last solution: " << r << '\n';
  }
  return EXIT_SUCCESS;
}
