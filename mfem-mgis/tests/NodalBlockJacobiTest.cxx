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
 * - its application to the right hand side matches the inverses of the
 *   nodal blocks computed on the host, from the matrix or from its diagonal
 *   block in parallel, with each backend, factorization and a damping;
 * - an operator giving its nodal blocks (as the tangent operator applied
 *   without assembly) gives the same preconditioner as its matrix, and
 *   blocks of size one are refused for it on all the processes;
 * - in parallel, the same checks with a process owning no element;
 * - invalid options and unavailable backends are refused.
 * \date   15/10/2026
 */

#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
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
#include "MFEMMGIS/MatrixFreeTangentOperator.hxx"
#include "ElasticitySystem.hxx"

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

using mfem_mgis::tests::ElasticitySystem;

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

/*!
 * \brief checks the application of the preconditioner to the right hand
 * side against a reference computed on the host: the inverses of the nodal
 * blocks of the matrix (of its diagonal block in parallel), computed by
 * `mfem::DenseMatrix::Invert`, with each available backend, with the
 * inverses and the LU factors of the blocks, and with a damping factor
 * \param[in] s: system
 * \param[in] ordering: name of the ordering of the unknowns
 */
template <bool is_parallel>
static void checkApplication(ElasticitySystem<is_parallel>& s,
                             const std::string& ordering) {
  using namespace mfem_mgis;
  using Batched = mfem::BatchedLinAlg;
  // reference, on the host
  auto local = mfem::SparseMatrix{};
  const mfem::SparseMatrix* a = nullptr;
  if constexpr (is_parallel) {
    s.A.GetDiag(local);
    a = &local;
  } else {
    a = &(s.A);
  }
  const auto n = a->Height();
  const auto nb = n / 3;
  const auto by_nodes = s.fespace.GetOrdering() == mfem::Ordering::byNODES;
  const auto unknown = [nb, by_nodes](const int node, const int c) {
    return by_nodes ? c * nb + node : 3 * node + c;
  };
  auto reference = mfem::Vector(n);
  {
    const auto* const I = a->HostReadI();
    const auto* const J = a->HostReadJ();
    const auto* const A = a->HostReadData();
    const auto* const b = s.B.HostRead();
    auto* const z = reference.HostWrite();
    for (int node = 0; node != nb; ++node) {
      auto block = mfem::DenseMatrix(3);
      block = 0.;
      for (int c = 0; c != 3; ++c) {
        const auto i = unknown(node, c);
        for (int k = I[i]; k != I[i + 1]; ++k) {
          for (int d = 0; d != 3; ++d) {
            if (J[k] == unknown(node, d)) {
              block(c, d) += A[k];
            }
          }
        }
      }
      block.Invert();
      for (int c = 0; c != 3; ++c) {
        auto v = mfem::real_t{0};
        for (int d = 0; d != 3; ++d) {
          v += block(c, d) * b[unknown(node, d)];
        }
        z[unknown(node, c)] = v;
      }
    }
  }
  const auto rmax = reference.Normlinf();
  const auto bound = double_precision ? 1e-12 : 1e-4;
  const auto backends = std::array<std::pair<const char*, Batched::Backend>, 3>{
      {{"Native", Batched::NATIVE},
       {"GPUBLAS", Batched::GPU_BLAS},
       {"MAGMA", Batched::MAGMA}}};
  for (const auto& [name, backend] : backends) {
    if (!Batched::IsAvailable(backend)) {
      continue;
    }
    for (const auto* const factorization : {"Inverse", "LU"}) {
      for (const auto damping : {1., 0.5}) {
        auto ctx = Context{};
        auto p = PreconditionerFactory<is_parallel>::getFactory().generate(
            ctx, "NodalBlockJacobi", s.fespace,
            Parameters{{"Backend", name},
                       {"Factorization", factorization},
                       {"Damping", damping}});
        if (p == nullptr) {
          throw std::runtime_error(ctx.getErrorMessage());
        }
        p->SetOperator(s.A);
        auto z = mfem::Vector(n);
        z.UseDevice(true);
        p->Mult(s.B, z);
        z.Add(-damping, reference);
        // a process owning no unknown has a null reference
        const auto e = rmax > 0 ? z.Normlinf() / rmax : z.Normlinf();
        check(e < bound,
              ordering + ": wrong application of the preconditioner (backend " +
                  name + ", " + factorization + ", damping " +
                  std::to_string(damping) + ", relative error " +
                  std::to_string(e) + ")");
      }
    }
  }
}  // end of checkApplication

/*!
 * \brief matrix giving its nodal blocks, extracted from the matrix (from its
 * diagonal block in parallel), as the tangent operator applied without
 * assembly does (see `mfem_mgis::NodalBlocksProvider`)
 */
template <bool is_parallel>
struct NodalBlocksOfMatrix final : mfem::Operator,
                                   mfem_mgis::NodalBlocksProvider {
  //! \param[in] s: system
  explicit NodalBlocksOfMatrix(ElasticitySystem<is_parallel>& s)
      : mfem::Operator(s.A.Height(), s.A.Width()), system(s) {}
  void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    this->system.A.Mult(x, y);
  }
  void AssembleNodalBlocks(mfem::Vector& b) const override {
    auto local = mfem::SparseMatrix{};
    const mfem::SparseMatrix* a = nullptr;
    if constexpr (is_parallel) {
      this->system.A.GetDiag(local);
      a = &local;
    } else {
      a = &(this->system.A);
    }
    const auto n = a->Height();
    const auto nb = n / 3;
    const auto by_nodes =
        this->system.fespace.GetOrdering() == mfem::Ordering::byNODES;
    const auto node = [nb, by_nodes](const int i) {
      return by_nodes ? i % nb : i / 3;
    };
    const auto component = [nb, by_nodes](const int i) {
      return by_nodes ? i / nb : i % 3;
    };
    const auto* const I = a->HostReadI();
    const auto* const J = a->HostReadJ();
    const auto* const A = a->HostReadData();
    b.SetSize(3 * n);
    auto* const pb = b.HostWrite();
    std::fill(pb, pb + 3 * n, mfem::real_t{0});
    // the entry (c1, c2) of the block of a node is the entry of the
    // component c1 of this node in the column c2
    for (int i = 0; i != n; ++i) {
      for (int k = I[i]; k != I[i + 1]; ++k) {
        if (node(J[k]) == node(i)) {
          pb[component(J[k]) * n + i] += A[k];
        }
      }
    }
  }

 private:
  //! \brief system
  ElasticitySystem<is_parallel>& system;
};  // end of struct NodalBlocksOfMatrix

/*!
 * \brief checks the preconditioner of an operator giving its nodal blocks:
 * its application to the right hand side is the one of the preconditioner of
 * the matrix, and blocks of size one are refused on all the processes (an
 * error on some processes only would block the other ones)
 * \param[in] s: system
 * \param[in] ordering: name of the ordering of the unknowns
 */
template <bool is_parallel>
static void checkNodalBlocksProvider(ElasticitySystem<is_parallel>& s,
                                     const std::string& ordering) {
  using namespace mfem_mgis;
  const auto op = NodalBlocksOfMatrix<is_parallel>(s);
  const auto apply = [&s](const mfem::Operator& a, const Parameters& options) {
    auto ctx = Context{};
    auto p = PreconditionerFactory<is_parallel>::getFactory().generate(
        ctx, "NodalBlockJacobi", s.fespace, options);
    if (p == nullptr) {
      throw std::runtime_error(ctx.getErrorMessage());
    }
    p->SetOperator(a);
    auto z = mfem::Vector(a.Height());
    z.UseDevice(true);
    p->Mult(s.B, z);
    return z;
  };
  auto z = apply(s.A, Parameters{});
  const auto zmax = z.Normlinf();
  z -= apply(op, Parameters{});
  check(z.Normlinf() <= (double_precision ? 1e-14 : 1e-6) * zmax,
        ordering +
            ": the blocks given by the operator differ from the ones of the "
            "matrix");
  auto refused = false;
  try {
    std::ignore = apply(op, Parameters{{"BlockSize", 1}});
  } catch (std::exception& e) {
    refused =
        std::string{e.what()}.find("number of components") != std::string::npos;
  }
  check(refused, ordering +
                     ": blocks of size one not refused for an operator "
                     "giving its nodal blocks");
}  // end of checkNodalBlocksProvider

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
  refused(Parameters{{"BlockSize", 2}}, "invalid block size");
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
      checkApplication(s1, "byNODES");
      checkApplication(s2, "byVDIM");
      checkNodalBlocksProvider(s1, "byNODES");
      checkNodalBlocksProvider(s2, "byVDIM");
      checkRefusals(s1);
      if constexpr (is_parallel) {
        // all the elements on the first process, the other ones owning no
        // unknown: the nodal blocks, which are in the diagonal blocks of the
        // processes, do not change
        const auto partitioning =
            std::vector<int>(static_cast<std::size_t>(mesh.GetNE()), 0);
        for (const auto o : {mfem::Ordering::byNODES, mfem::Ordering::byVDIM}) {
          auto s =
              ElasticitySystem<is_parallel>(mesh, fec, o, partitioning.data());
          const auto name =
              std::string{o == mfem::Ordering::byNODES ? "byNODES" : "byVDIM"} +
              ", empty process";
          check((rank == 0) || (s.A.Height() == 0),
                name + ": the process owns unknowns");
          check(areClose(run(s, name), i1),
                name + ": the partitioning changes the number of iterations");
          checkApplication(s, name);
          checkNodalBlocksProvider(s, name);
        }
      }
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
