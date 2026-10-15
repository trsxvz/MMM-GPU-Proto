/*!
 * \file   NodalBlockJacobiTest.cxx
 * \brief  Test of the Jacobi preconditioner by nodal blocks
 * (`NodalBlockJacobi`) on a cantilever beam in linear elasticity, with the
 * unknowns ordered by nodes and by vector dimension, in sequential and in
 * parallel (`--parallel`):
 *
 * - the conjugate gradient of MFEM needs the same number of iterations with
 *   each available backend of the batched linear algebra of MFEM (native
 *   kernels, cuBLAS or hipBLAS, MAGMA), with the inverses or the LU factors
 *   of the blocks, and with both orderings of the unknowns;
 * - with blocks of size one, the preconditioner is the Jacobi smoother of
 *   MFEM; the nodal blocks need fewer iterations (in double precision; in
 *   single precision, the numbers of iterations are only compared up to a
 *   quarter); the GMRES solver of MFEM
 *   converges with the default backend (MAGMA, if available, which must be
 *   synchronized with the kernels of MFEM);
 * - in sequential, the preconditioner of the block diagonal part of the
 *   matrix is its inverse: the conjugate gradient converges in one
 *   iteration (check of the extraction of the blocks);
 * - invalid options and unavailable backends are refused.
 * \date   15/10/2026
 */

#include <array>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/PreconditionerFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

//! \brief number of failed checks
static int failures = 0;
//! \brief rank of the process
static int rank = 0;
//! \brief parallel computations
static bool parallel = false;
//! \brief relative tolerance of the conjugate gradient
static constexpr double tolerance =
    std::is_same_v<mfem::real_t, double> ? 1e-10 : 1e-5;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief cantilever beam in linear elasticity, clamped on the face x = 0
template <bool is_parallel>
struct ElasticitySystem {
  using Matrix =
      std::conditional_t<is_parallel, mfem::HypreParMatrix, mfem::SparseMatrix>;
  /*!
   * \param[in] smesh: sequential mesh (partitioned in parallel)
   * \param[in] fec: finite element collection
   * \param[in] o: ordering of the unknowns
   */
  ElasticitySystem(mfem::Mesh& smesh,
                   mfem::FiniteElementCollection& fec,
                   const mfem::Ordering::Type o)
      : mesh(makeMesh(smesh)), fespace(mesh.get(), &fec, 3, o), a(&fespace) {
    auto lambda = mfem::ConstantCoefficient(1.5);
    auto mu = mfem::ConstantCoefficient(1);
    this->a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
    this->a.Assemble();
    auto f = mfem::Vector(3);
    f = 0.;
    f(2) = -1;
    auto fc = mfem::VectorConstantCoefficient(f);
    auto b = mfem_mgis::LinearForm<is_parallel>(&(this->fespace));
    b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(fc));
    b.Assemble();
    auto x = mfem_mgis::GridFunction<is_parallel>(&(this->fespace));
    x = 0.;
    auto boundaries = mfem::Array<int>(this->mesh->bdr_attributes.Max());
    boundaries = 0;
    boundaries[4] = 1;
    auto ess_tdofs = mfem::Array<int>{};
    this->fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
    auto X = mfem::Vector{};
    auto rhs = mfem::Vector{};
    this->a.FormLinearSystem(ess_tdofs, x, b, this->A, X, rhs);
    this->B = rhs;
  }
  std::unique_ptr<mfem_mgis::Mesh<is_parallel>> mesh;
  mfem_mgis::FiniteElementSpace<is_parallel> fespace;
  mfem_mgis::BilinearForm<is_parallel> a;
  Matrix A;
  mfem::Vector B;

 private:
  static std::unique_ptr<mfem_mgis::Mesh<is_parallel>> makeMesh(
      mfem::Mesh& smesh) {
    if constexpr (is_parallel) {
      return std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, smesh);
    } else {
      return std::make_unique<mfem::Mesh>(smesh);
    }
  }
};

/*!
 * \return the number of iterations of the conjugate gradient of MFEM with
 * a preconditioner, if it converged
 * \param[in] s: system
 * \param[in] a: matrix
 * \param[in] p: parameters of the preconditioner
 * \param[in] name: name of the resolution
 * \param[in] solver: name of the Krylov solver
 */
template <bool is_parallel>
static std::optional<int> solve(
    ElasticitySystem<is_parallel>& s,
    const typename ElasticitySystem<is_parallel>::Matrix& a,
    const mfem_mgis::Parameters& p,
    const std::string& name,
    const char* const solver = "CGSolver") {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto ls = LinearSolverFactory<is_parallel>::getFactory().generate(
      ctx, solver, s.fespace,
      Parameters{{"RelativeTolerance", tolerance},
                 {"AbsoluteTolerance", 0.},
                 {"MaximumNumberOfIterations", 5000},
                 {"Preconditioner", p}});
  if (isInvalid(ls)) {
    throw std::runtime_error(ctx.getErrorMessage());
  }
  auto x = mfem::Vector(s.B.Size());
  x.UseDevice(true);
  x = 0.;
  ls.linear_solver->SetOperator(a);
  ls.linear_solver->Mult(s.B, x);
  const auto converged = hasConverged(*(ls.linear_solver));
  const auto iterations =
      getNumberOfIterationsAtConvergence(*(ls.linear_solver)).value_or(-1);
  if (rank == 0) {
    std::cout << name << ": " << iterations << " iterations"
              << (converged ? "" : " (not converged)") << '\n';
  }
  if (!converged) {
    return {};
  }
  return iterations;
}  // end of solve

//! \return the parameters of the preconditioner
static mfem_mgis::Parameters makeNodalBlockJacobi(
    const mfem_mgis::Parameters& options) {
  return mfem_mgis::Parameters{{"Name", "NodalBlockJacobi"},
                               {"Options", options}};
}  // end of makeNodalBlockJacobi

//! \brief computations in double precision
static constexpr bool double_precision = std::is_same_v<mfem::real_t, double>;

/*!
 * \return true if two numbers of iterations are close: equal up to one in
 * double precision, up to a quarter in single precision, where the
 * rounding errors change the numbers of iterations of the variants of the
 * same preconditioner (inverse or LU factors, ordering of the unknowns)
 */
static bool areClose(const std::optional<int>& a, const std::optional<int>& b) {
  if ((!a.has_value()) || (!b.has_value())) {
    return false;
  }
  const auto d = std::abs(*a - *b);
  return double_precision ? d <= 1 : 4 * d <= std::max(*a, *b);
}  // end of areClose

/*!
 * \return the number of iterations with the nodal blocks
 * \param[in] s: system
 * \param[in] ordering: name of the ordering
 */
template <bool is_parallel>
static std::optional<int> run(ElasticitySystem<is_parallel>& s,
                              const std::string& ordering) {
  using namespace mfem_mgis;
  using Batched = mfem::BatchedLinAlg;
  const auto reference = solve(s, s.A, makeNodalBlockJacobi(Parameters{}),
                               "NodalBlockJacobi, " + ordering);
  check(reference.has_value(), ordering + ": resolution");
  const auto backends = std::array<std::pair<const char*, Batched::Backend>, 3>{
      {{"Native", Batched::NATIVE},
       {"GPUBLAS", Batched::GPU_BLAS},
       {"MAGMA", Batched::MAGMA}}};
  for (const auto& [name, backend] : backends) {
    if (!Batched::IsAvailable(backend)) {
      if (rank == 0) {
        std::cout << "backend " << name << ": not available\n";
      }
      continue;
    }
    for (const auto* const factorization : {"Inverse", "LU"}) {
      const auto it =
          solve(s, s.A,
                makeNodalBlockJacobi(Parameters{
                    {"Backend", name}, {"Factorization", factorization}}),
                std::string{"NodalBlockJacobi ("} + name + ", " +
                    factorization + "), " + ordering);
      check(areClose(it, reference),
            ordering + ": the backend " + name + " with the factorization " +
                factorization + " changes the number of iterations");
    }
  }
  // GMRES of MFEM, which gives to the preconditioner the vectors it has
  // just computed: with the backends computing on their own queue (MAGMA),
  // the device must be synchronized
  for (const auto bs : {3, 1}) {
    const auto it =
        solve(s, s.A, makeNodalBlockJacobi(Parameters{{"BlockSize", bs}}),
              "GMRESSolver + NodalBlockJacobi (blocks of size " +
                  std::to_string(bs) + "), " + ordering,
              "GMRESSolver");
    const auto jacobi =
        solve(s, s.A, Parameters{{"Name", "DSmoother"}},
              "GMRESSolver + DSmoother, " + ordering, "GMRESSolver");
    // in single precision, the restarted GMRES stagnates close to its
    // tolerance: only its convergence is checked
    check(it.has_value() && jacobi.has_value() &&
              ((bs == 3) || (!double_precision) || areClose(it, jacobi)),
          ordering + ": GMRES with the nodal blocks of size " +
              std::to_string(bs));
  }
  // blocks of size one: Jacobi preconditioner
  const auto scalar =
      solve(s, s.A, makeNodalBlockJacobi(Parameters{{"BlockSize", 1}}),
            "NodalBlockJacobi (blocks of size one), " + ordering);
  const auto jacobi = solve(s, s.A, Parameters{{"Name", "DSmoother"}},
                            "DSmoother, " + ordering);
  check(areClose(scalar, jacobi),
        ordering + ": blocks of size one differ from the Jacobi smoother");
  if constexpr (double_precision) {
    check(
        reference.has_value() && scalar.has_value() && (*reference < *scalar),
        ordering + ": the nodal blocks do not reduce the number of iterations");
  }
  // block diagonal part of the matrix: the preconditioner is its inverse
  if constexpr (!is_parallel) {
    const auto n = s.A.Height();
    const auto nnodes = n / 3;
    const auto by_nodes = s.fespace.GetOrdering() == mfem::Ordering::byNODES;
    const auto node = [nnodes, by_nodes](const int i) {
      return by_nodes ? i % nnodes : i / 3;
    };
    auto d = mfem::SparseMatrix(s.A);
    const auto* const I = d.HostReadI();
    const auto* const J = d.HostReadJ();
    auto* const values = d.HostReadWriteData();
    for (int i = 0; i != n; ++i) {
      for (int k = I[i]; k != I[i + 1]; ++k) {
        if (node(i) != node(J[k])) {
          values[k] = 0;
        }
      }
    }
    const auto it =
        solve(s, d, makeNodalBlockJacobi(Parameters{}),
              "NodalBlockJacobi, block diagonal matrix, " + ordering);
    check(it.has_value() && (*it <= 1),
          ordering +
              ": the preconditioner is not the inverse of the block "
              "diagonal part");
  }
  return reference;
}  // end of run

//! \brief checks the refusals of invalid options
template <bool is_parallel>
static void checkRefusals(ElasticitySystem<is_parallel>& s) {
  using namespace mfem_mgis;
  const auto refused = [&s](const Parameters& options, const std::string& msg) {
    auto ctx = Context{};
    const auto p = PreconditionerFactory<is_parallel>::getFactory().generate(
        ctx, "NodalBlockJacobi", s.fespace, options);
    const auto error = std::string{ctx.getErrorMessage()};
    check((p == nullptr) && (error.find(msg) != std::string::npos),
          "invalid options not refused ('" + error + "')");
  };
  refused(Parameters{{"Backend", "cuSOLVER"}}, "invalid backend");
  refused(Parameters{{"Factorization", "Cholesky"}}, "invalid factorization");
  refused(Parameters{{"BlockSize", 0}}, "invalid block size");
  refused(Parameters{{"Damping", -1.}}, "invalid damping factor");
  if (!mfem::BatchedLinAlg::IsAvailable(mfem::BatchedLinAlg::MAGMA)) {
    refused(Parameters{{"Backend", "MAGMA"}}, "is not available");
  }
}  // end of checkRefusals

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "parallel computations");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
#ifdef MFEM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* MFEM_USE_MPI */
  auto mesh = mfem::Mesh::MakeCartesian3D(12, 3, 3, mfem::Element::HEXAHEDRON,
                                          4., 1., 1.);
  auto fec = mfem::H1_FECollection(1, 3);
  try {
    const auto test = [&mesh, &fec]<bool is_parallel>() {
      auto s1 =
          ElasticitySystem<is_parallel>(mesh, fec, mfem::Ordering::byNODES);
      auto s2 =
          ElasticitySystem<is_parallel>(mesh, fec, mfem::Ordering::byVDIM);
      const auto i1 = run(s1, "byNODES");
      const auto i2 = run(s2, "byVDIM");
      check(areClose(i1, i2), "the orderings change the number of iterations");
      checkRefusals(s1);
    };
    if (parallel) {
#ifdef MFEM_USE_MPI
      test.template operator()<true>();
#else
      throw std::runtime_error("MFEM is not built with MPI");
#endif
    } else {
      test.template operator()<false>();
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
