// Factorizations of the elasticity matrix of ex8 of mfem-mgis-examples
// (linear elements, byNODES ordering, all boundaries clamped) by MUMPS on
// the host and by cuDSS on the GPU, with and without the reuse of the
// analysis, for matrices with the same pattern and different values.
//
// Usage: cudss-bench mesh.msh [number of factorizations] [device] [hypre]
// hypre: device (default) or host (hypre computes on the host, as with the
// option --hypre-on-host of mfem-mgis)
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include "mfem.hpp"

//! \return the relative residual of the solution of A x = b
static double residual(const mfem::HypreParMatrix& A, const mfem::Vector& b,
                       const mfem::Vector& x) {
  auto r = mfem::Vector(b.Size());
  r.UseDevice(true);
  A.Mult(x, r);
  r -= b;
  return std::sqrt(mfem::InnerProduct(A.GetComm(), r, r) /
                   mfem::InnerProduct(A.GetComm(), b, b));
}

/*!
 * \brief time the factorizations of `n` matrices `(1 + 0.01 i) A`
 * \param[in] name: name of the solver
 * \param[in] s: solver
 * \param[in] A: matrix, scaled back at the end
 */
//! \brief progress message of all the processes
static void progress(const std::string& what) {
  std::cerr << "rank " << mfem::Mpi::WorldRank() << ": " << what << std::endl;
}

static void run(const std::string& name, mfem::Solver& s,
                mfem::HypreParMatrix& A, const int n, const mfem::Vector& b,
                mfem::Vector& x) {
  progress(name);
  mfem::StopWatch sw;
  auto first = 0.;
  auto next = 0.;
  auto solve = 0.;
  auto scale = 1.;
  for (int i = 0; i != n; ++i) {
    const auto f = 1 + 0.01 * i;
    A *= f / scale;
    scale = f;
    MPI_Barrier(A.GetComm());
    sw.Clear();
    sw.Start();
    s.SetOperator(A);
    progress(name + ": operator " + std::to_string(i));
    MPI_Barrier(A.GetComm());
    sw.Stop();
    (i == 0 ? first : next) += sw.RealTime();
    sw.Clear();
    sw.Start();
    x = 0.;
    s.Mult(b, x);
    MPI_Barrier(A.GetComm());
    sw.Stop();
    solve += sw.RealTime();
  }
  const auto r = residual(A, b, x);
  A *= 1 / scale;
  if (mfem::Mpi::Root()) {
    std::cout << name << ": first " << first << " s, next ones "
              << next / (n - 1) << " s each, solve " << solve / n
              << " s, relative residual " << r << '\n';
  }
}

int main(int argc, char** argv) {
  mfem::Mpi::Init(argc, argv);
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " mesh [n] [device]\n";
    return EXIT_FAILURE;
  }
  const int n = argc > 2 ? std::atoi(argv[2]) : 6;
  const auto hypre_on_host = (argc > 4) && (std::string{argv[4]} == "host");
  if (hypre_on_host) {
    mfem::Hypre::configure_runtime_policy_from_mfem = false;
  }
  auto device = mfem::Device(argc > 3 ? argv[3] : "cuda");
  mfem::Hypre::Init();
  if (hypre_on_host) {
    HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
    HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
  }
  if (mfem::Mpi::Root()) {
    device.Print();
    std::cout << "hypre on the " << (hypre_on_host ? "host" : "device")
              << '\n';
  }
  progress("mesh");
  auto smesh = mfem::Mesh(argv[1], 1, 1);
  progress("serial mesh");
  auto mesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  progress("parallel mesh");
  smesh.Clear();
  const auto dim = mesh.Dimension();
  auto fec = mfem::H1_FECollection(1, dim);
  auto fes = mfem::ParFiniteElementSpace(&mesh, &fec, dim,
                                         mfem::Ordering::byNODES);
  if (mfem::Mpi::Root()) {
    std::cout << "unknowns: " << fes.GlobalTrueVSize() << '\n';
  }
  auto ess_bdr = mfem::Array<int>(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  auto lambda = mfem::ConstantCoefficient(1.);
  auto mu = mfem::ConstantCoefficient(1.);
  auto a = mfem::ParBilinearForm(&fes);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
  a.Assemble();
  progress("assembly");
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
  progress("linear system");
  {
    auto s = mfem::MUMPSSolver(A.GetComm());
    s.SetPrintLevel(0);
    s.SetMatrixSymType(mfem::MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
    run("MUMPS, no reuse", s, A, n, B, X);
  }
  {
    auto s = mfem::MUMPSSolver(A.GetComm());
    s.SetPrintLevel(0);
    s.SetMatrixSymType(mfem::MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
    s.SetReorderingReuse(true);
    run("MUMPS, reuse", s, A, n, B, X);
  }
  for (const auto reuse : {false, true}) {
    for (const auto spd : {false, true}) {
      auto s = mfem::CuDSSSolver(A.GetComm());
      if (spd) {
        s.SetMatrixSymType(
            mfem::CuDSSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
        s.SetMatrixViewType(mfem::CuDSSSolver::MatViewType::LOWER);
      }
      s.SetReorderingReuse(reuse);
      run(std::string{"cuDSS, "} + (spd ? "SPD lower" : "general") +
              (reuse ? ", reuse" : ", no reuse"),
          s, A, n, B, X);
    }
  }
  return EXIT_SUCCESS;
}
