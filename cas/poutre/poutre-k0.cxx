/*!
 * \file   cas/poutre/poutre-k0.cxx
 * \brief  Resolution of the elastic operator of the beam, the operator K0
 * of the NewtonK0 method, by a linear solver of the factory of mfem-mgis:
 * same discretization as the beam (finite element space of mfem-mgis,
 * unknowns ordered by nodes), same elastic properties as its behaviour,
 * same essential boundary conditions (clamped end, transverse displacement
 * at the other end). The setup of the solver (`SetOperator`) and the
 * resolutions are timed separately, the times being the largest ones of
 * the processes. The setup of the preconditioners of hypre is done at the
 * first resolution: the first resolution and the next ones are reported
 * separately.
 */

#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "solver-options.hxx"

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh_file = "beam-tet.mesh";
  const char* solver = "cg-amg";
  const char* amg_options = "";
  const char* amgx_configuration_file = "";
  const char* amgx_options = "";
  auto amgx_single_precision = false;
  auto refinements = 3;
  auto resolutions = 3;
  auto tolerance = real{1e-8};
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file.");
  args.AddOption(&refinements, "-r", "--refinements",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&solver, "-s", "--solver",
                 "Linear solver: cg-amg (conjugate gradient of MFEM and one "
                 "cycle of BoomerAMG), pcg-amg (conjugate gradient of hypre "
                 "and one cycle of BoomerAMG), cg-amgx (conjugate gradient "
                 "of MFEM and one cycle of AmgX), amgx-pcg (conjugate "
                 "gradient of AmgX) or cg-jacobi.");
  args.AddOption(&amg_options, "-amgo", "--amg-options",
                 "Options of BoomerAMG, as Name=value separated by commas.");
  args.AddOption(&amgx_configuration_file, "-amgxc",
                 "--amgx-configuration-file",
                 "Configuration file of the cycle of AmgX (cg-amgx).");
  args.AddOption(&amgx_single_precision, "-amgxs",
                 "--amgx-single-precision", "-no-amgxs",
                 "--no-amgx-single-precision",
                 "Computations of AmgX in single precision (cg-amgx, "
                 "amgx-pcg; tolerance of amgx-pcg of at least 1e-6).");
  args.AddOption(&amgx_options, "-amgxo", "--amgx-options",
                 "Parameters of AmgX completing the default configuration "
                 "of mfem-mgis (one V-cycle of the classical multigrid), as "
                 "name=value separated by commas with the names of AmgX, "
                 "for instance strength_threshold=0.8,interpolator=MULTIPASS "
                 "(cg-amgx, amgx-pcg).");
  args.AddOption(&resolutions, "-n", "--resolutions",
                 "Number of resolutions after the setup.");
  args.AddOption(&tolerance, "-lt", "--linear-tolerance",
                 "Relative tolerance of the conjugate gradient.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto ctx = Context{};
  // discretization of the beam (see poutre.cxx)
  auto fed = FiniteElementDiscretization{
      ctx, Parameters{{"MeshFileName", mesh_file},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", 1},
                      {"UnknownsSize", 3},
                      {"NumberOfUniformRefinements", refinements},
                      {"Parallel", true}}};
  auto& pmesh = fed.getMesh<true>();
  auto& fespace = fed.getFiniteElementSpace<true>();
  const auto dim = pmesh.Dimension();
  const auto root = fespace.GetMyRank() == 0;
  // essential degrees of freedom of the beam: clamped end (boundary 1),
  // transverse displacement imposed at the other end (boundary 2)
  auto ess_tdofs = mfem::Array<int>{};
  {
    auto clamped = mfem::Array<int>(pmesh.bdr_attributes.Max());
    clamped = 0;
    clamped[0] = 1;
    fespace.GetEssentialTrueDofs(clamped, ess_tdofs);
    auto end = mfem::Array<int>(pmesh.bdr_attributes.Max());
    end = 0;
    end[1] = 1;
    auto end_tdofs = mfem::Array<int>{};
    fespace.GetEssentialTrueDofs(end, end_tdofs, 1);
    ess_tdofs.Append(end_tdofs);
  }
  // linear elasticity, with the elastic properties of the behaviour of the
  // beam (IsotropicLinearHardeningPlasticity), the transverse displacement
  // imposed at the end giving the right hand side, as for the elastic
  // prediction
  const auto young = real{70e9};
  const auto poisson = real{0.34};
  auto lambda = mfem::ConstantCoefficient(young * poisson /
                                          ((1 + poisson) * (1 - 2 * poisson)));
  auto mu = mfem::ConstantCoefficient(young / (2 * (1 + poisson)));
  auto a = mfem::ParBilinearForm(&fespace);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
  a.Assemble();
  auto b = mfem::ParLinearForm(&fespace);
  b.Assemble();
  auto x = mfem::ParGridFunction(&fespace);
  x = 0.0;
  {
    auto u = mfem::Vector(dim);
    u = 0.0;
    u(1) = 0.05;
    auto uc = mfem::VectorConstantCoefficient(u);
    auto end = mfem::Array<int>(pmesh.bdr_attributes.Max());
    end = 0;
    end[1] = 1;
    x.ProjectBdrCoefficient(uc, end);
  }
  auto A = mfem::HypreParMatrix{};
  auto B = mfem::Vector{};
  auto X = mfem::Vector{};
  a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
  // linear solver
  const auto amg = parseBoomerAMGOptions(amg_options);
  const auto amgx =
      parseAmgXOptions(amgx_configuration_file, std::string{amgx_options},
                       amgx_single_precision);
  const auto s = std::string{solver};
  const auto cg = [tolerance](const Parameters& p) {
    return Parameters{{"RelativeTolerance", tolerance},
                      {"AbsoluteTolerance", 0.},
                      {"MaximumNumberOfIterations", 10000},
                      {"Preconditioner", p}};
  };
  auto& factory = LinearSolverFactory<true>::getFactory();
  auto ls = [&]() -> LinearSolverHandler {
    if (s == "cg-amg") {
      return factory.generate(
          ctx, "CGSolver", fespace,
          cg(Parameters{{"Name", "HypreBoomerAMG"}, {"Options", amg}}));
    } else if (s == "pcg-amg") {
      return factory.generate(
          ctx, "HyprePCG", fespace,
          Parameters{{"Tolerance", tolerance},
                     {"MaximumNumberOfIterations", 10000},
                     {"Preconditioner", Parameters{{"Name", "HypreBoomerAMG"},
                                                   {"Options", amg}}}});
    } else if (s == "cg-amgx") {
      return factory.generate(
          ctx, "CGSolver", fespace,
          cg(Parameters{{"Name", "AmgX"}, {"Options", amgx}}));
    } else if (s == "amgx-pcg") {
      return factory.generate(ctx, "AmgXPCG", fespace,
                              amgxPCGOptions(amgx, tolerance, 10000));
    } else if (s == "cg-jacobi") {
      return factory.generate(
          ctx, "CGSolver", fespace,
          cg(Parameters{{"Name", "OperatorJacobiSmoother"}}));
    }
    return ctx.registerErrorMessage("invalid linear solver '" + s + "'");
  }();
  if (isInvalid(ls)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  // the kernels of a GPU are asynchronous: the GPU is synchronized before
  // reading the clock; the time is the largest one of the processes
  using clock = std::chrono::steady_clock;
  const auto gpu = mfem::Device::Allows(mfem::Backend::CUDA_MASK);
  const auto synchronize = [gpu] {
    if (gpu) {
      MFEM_DEVICE_SYNC;
    }
  };
  const auto elapsed = [&synchronize, &fespace](const clock::time_point& t0) {
    synchronize();
    auto t = std::chrono::duration<double>(clock::now() - t0).count();
    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MAX,
                  fespace.GetComm());
    return t;
  };
  synchronize();
  MPI_Barrier(MPI_COMM_WORLD);
  auto t0 = clock::now();
  ls.linear_solver->SetOperator(A);
  const auto setup = elapsed(t0);
  if (root) {
    std::cout << "dofs " << A.GetGlobalNumRows() << " | solver " << s
              << " | setup " << std::setprecision(4) << setup << " s\n";
  }
  auto norm = [&fespace](const mfem::Vector& v) {
    return std::sqrt(mfem::InnerProduct(fespace.GetComm(), v, v));
  };
  const auto nB = norm(B);
  auto success = true;
  for (int i = 0; i != resolutions + 1; ++i) {
    X.UseDevice(true);
    X = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = clock::now();
    ls.linear_solver->Mult(B, X);
    const auto t = elapsed(t0);
    const auto converged = hasConverged(*(ls.linear_solver));
    const auto iterations =
        getNumberOfIterationsAtConvergence(*(ls.linear_solver));
    auto r = mfem::Vector(B.Size());
    r.UseDevice(true);
    A.Mult(X, r);
    r -= B;
    const auto residual = norm(r) / nB;
    success = success && converged;
    if (root) {
      std::cout << (i == 0 ? "first resolution" : "resolution") << " | time "
                << std::setprecision(4) << t << " s | iterations "
                << (iterations.has_value() ? std::to_string(*iterations)
                                           : std::string{"-"})
                << " | residual " << residual
                << (converged ? "" : " | not converged") << '\n';
    }
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
