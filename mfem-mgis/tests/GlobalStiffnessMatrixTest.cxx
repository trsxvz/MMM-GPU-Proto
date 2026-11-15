/*!
 * \file   GlobalStiffnessMatrixTest.cxx
 * \brief  Test of the assembly of the global matrix of a parallel finite
 * element space from the local matrices (`batch::GlobalStiffnessMatrix`):
 * same matrix as the product `Pᵀ A P` of hypre, at the first assembly and
 * after a change of the values (the essential degrees of freedom being
 * eliminated from the global matrix in between, as by the non linear
 * problems), for scalar and vector spaces of orders 1 and 2; a local matrix
 * with another pattern or other arrays on one process, or with its columns
 * sorted in other arrays, makes all the processes compute the pattern
 * again; a non conforming space is refused. A checksum of the entries of
 * each global matrix is printed: it does not depend on the order of the
 * entries, and it is the same when hypre computes on the host and on the
 * GPU (device given by the MFEM_MGIS_DEVICE environment variable), the
 * values being the same bit for bit (see GlobalStiffnessMatrixTest.cmake).
 * Last, products and smoothers of hypre are applied after a second assembly
 * in the same global matrix: they must use the new values (on a GPU, MFEM
 * and hypre keep in the matrix data computed from its values). To be run on
 * several processes.
 * \date   11/10/2026
 */

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/handle.hpp"
#include "mfem/linalg/hypre.hpp"
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pfespace.hpp"
#include "Tolerance.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"
#include "MFEMMGIS/Batch/GlobalStiffnessMatrix.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \brief sets the values of a local matrix from a seed
 * \param[in, out] A: matrix
 * \param[in] seed: seed
 */
static void setValues(mfem::SparseMatrix& A, const int seed) {
  auto* const v = A.HostWriteData();
  for (int k = 0; k != A.NumNonZeroElems(); ++k) {
    v[k] = std::sin(0.37 * k + seed) + 0.1 * seed;
  }
}  // end of setValues

/*!
 * \return the maximum relative difference between the products of two
 * global matrices by vectors
 * \param[in] K: first matrix
 * \param[in] R: second matrix
 */
static mfem_mgis::real compareMatrices(const mfem::HypreParMatrix& K,
                                       const mfem::HypreParMatrix& R) {
  using mfem_mgis::real;
  auto d = real{};
  for (int s = 0; s != 3; ++s) {
    auto x = mfem::Vector(R.Width());
    auto* const xv = x.HostWrite();
    const auto first = R.ColPart()[0];
    for (int i = 0; i != x.Size(); ++i) {
      xv[i] = std::cos(0.11 * static_cast<real>(first + i) + s);
    }
    auto y = mfem::Vector(R.Height());
    auto yr = mfem::Vector(R.Height());
    K.Mult(x, y);
    R.Mult(x, yr);
    auto ld = real{};
    auto lm = real{};
    const auto* const yv = y.HostRead();
    const auto* const rv = yr.HostRead();
    for (int i = 0; i != y.Size(); ++i) {
      ld = std::max(ld, std::abs(yv[i] - rv[i]));
      lm = std::max(lm, std::abs(rv[i]));
    }
    auto gd = real{};
    auto gm = real{};
    const auto mpi_real = mfem::MPITypeMap<real>::mpi_type;
    MPI_Allreduce(&ld, &gd, 1, mpi_real, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lm, &gm, 1, mpi_real, MPI_MAX, MPI_COMM_WORLD);
    d = std::max(d, gd / gm);
  }
  return d;
}  // end of compareMatrices

/*!
 * \return a checksum of the entries of a global matrix (global row, global
 * column and bits of the value of each entry), which does not depend on the
 * order of the entries nor on the process owning them
 * \param[in] K: matrix
 */
static std::uint64_t computeChecksum(const mfem::HypreParMatrix& K) {
  // mixing function of splitmix64
  const auto mix = [](std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  };
  K.HostRead();
  const hypre_ParCSRMatrix* const k = K;
  const auto first_row = hypre_ParCSRMatrixFirstRowIndex(k);
  const auto first_column = hypre_ParCSRMatrixFirstColDiag(k);
  const auto* const cmap = hypre_ParCSRMatrixColMapOffd(k);
  auto checksum = std::uint64_t{};
  const auto add = [&checksum, &mix, first_row](const hypre_CSRMatrix* const b,
                                                const auto& column) {
    const auto* const i = hypre_CSRMatrixI(b);
    const auto* const j = hypre_CSRMatrixJ(b);
    const auto* const v = hypre_CSRMatrixData(b);
    for (HYPRE_Int r = 0; r != hypre_CSRMatrixNumRows(b); ++r) {
      for (auto p = i[r]; p != i[r + 1]; ++p) {
        auto bits = std::uint64_t{};
        std::memcpy(&bits, v + p, sizeof(v[p]));
        checksum ^= mix(mix(mix(static_cast<std::uint64_t>(first_row + r)) ^
                            static_cast<std::uint64_t>(column(j[p]))) ^
                        bits);
      }
    }
  };
  add(hypre_ParCSRMatrixDiag(k),
      [first_column](const HYPRE_Int c) { return first_column + c; });
  add(hypre_ParCSRMatrixOffd(k), [cmap](const HYPRE_Int c) { return cmap[c]; });
  K.HypreRead();
  auto all = std::uint64_t{};
  MPI_Allreduce(&checksum, &all, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD);
  return all;
}  // end of computeChecksum

/*!
 * \brief compares the global matrices assembled by
 * `batch::GlobalStiffnessMatrix` and by hypre
 * \param[in] mesh: mesh
 * \param[in] order: order of the finite elements
 * \param[in] vdim: number of components of the unknowns
 * \param[in] name: name of the case
 */
static void testCase(mfem::ParMesh& mesh,
                     const int order,
                     const int vdim,
                     const std::string& name) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto fec = mfem::H1_FECollection(order, mesh.Dimension());
  auto fes = mfem::ParFiniteElementSpace(&mesh, &fec, vdim);
  check(batch::GlobalStiffnessMatrix::isSupported(fes),
        name + ": conforming space not supported");
  auto dofs = mfem::Array<int>{};
  check(batch::computeElementsDegreesOfFreedom(ctx, dofs, fes),
        name + ": degrees of freedom of the elements");
  const auto A = batch::makeAssembledStiffnessMatrix(ctx, fes, dofs);
  check(A != nullptr, name + ": local matrix");
  if (A == nullptr) {
    return;
  }
  // essential degrees of freedom of the first boundary, eliminated from the
  // global matrix after each assembly, as by the non linear problems
  auto marker = mfem::Array<int>(mesh.bdr_attributes.Max());
  marker = 0;
  marker[0] = 1;
  auto essential_dofs = mfem::Array<int>{};
  fes.GetEssentialTrueDofs(marker, essential_dofs);
  auto g = batch::GlobalStiffnessMatrix(fes);
  // assembles the global matrix of a local matrix, and compares it to the
  // product of hypre
  const auto assemble =
      [&fes, &g, &name, &essential_dofs](
          const mfem::SparseMatrix& M,
          const std::string& what) -> const mfem::HypreParMatrix* {
    auto c = Context{};
    auto dA = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto Ph = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto R = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
    auto copy = mfem::SparseMatrix(M);
    dA.MakeSquareBlockDiag(fes.GetComm(), fes.GlobalVSize(),
                           fes.GetDofOffsets(), &copy);
    Ph.ConvertFrom(fes.Dof_TrueDof_Matrix());
    R.MakePtAP(dA, Ph);
    auto* const K = g.assemble(c, M);
    check(K != nullptr,
          name + ", " + what + ": assembly (" + c.getErrorMessage() + ")");
    if (K == nullptr) {
      return nullptr;
    }
    const auto d = compareMatrices(*K, *(R.As<mfem::HypreParMatrix>()));
    const auto checksum = computeChecksum(*K);
    if (mfem::Mpi::Root()) {
      std::cout << name << ", " << what << ": maximum relative difference " << d
                << '\n'
                << name << ", " << what << ": checksum " << std::hex << checksum
                << std::dec << '\n';
    }
    // in float, the relative difference reaches 1e-7 (sums in another
    // order)
    check(d < mfem_mgis::unit_tests::tolerance(1e-14, 1e-6),
          name + ", " + what + ": global matrix");
    K->EliminateBC(essential_dofs, mfem::Operator::DIAG_ONE);
    return K;
  };
  const auto* address = static_cast<const mfem::HypreParMatrix*>(nullptr);
  for (const auto seed : {1, 2}) {
    setValues(*A, seed);
    const auto* const K = assemble(*A, "values " + std::to_string(seed));
    check((address == nullptr) || (address == K),
          name + ": address of the global matrix changed");
    address = K;
  }
  // a local matrix of another pattern on one process makes all the
  // processes compute the pattern again
  auto other = std::make_unique<mfem::SparseMatrix>(*A);
  if (mfem::Mpi::Root()) {
    other = std::make_unique<mfem::SparseMatrix>(A->Height());
    for (int i = 0; i != A->Height(); ++i) {
      other->Set(i, i, 1 + i);
    }
    other->Finalize();
  } else {
    setValues(*other, 3);
  }
  static_cast<void>(assemble(*other, "another pattern on one process"));
  setValues(*A, 4);
  static_cast<void>(assemble(*A, "first pattern again"));
  // the same entries in other arrays, the columns being sorted, on all the
  // processes: the arrays of the pattern are compared, not their contents
  auto sorted = mfem::SparseMatrix(*A);
  sorted.SortColumnIndices();
  setValues(sorted, 5);
  static_cast<void>(assemble(sorted, "columns sorted in other arrays"));
  // other arrays on one process only (the other ones giving the arrays of
  // the previous assembly): new setup on all the processes
  setValues(*A, 6);
  static_cast<void>(assemble(*A, "first pattern on all the processes"));
  auto copy = mfem::SparseMatrix(*A);
  setValues(copy, 7);
  setValues(*A, 7);
  static_cast<void>(
      assemble(mfem::Mpi::Root() ? copy : *A, "other arrays on one process"));
}  // end of testCase

/*!
 * \brief sets the values of a local matrix of symmetric pattern: a
 * symmetric matrix with a dominant diagonal
 * \param[in, out] A: matrix
 * \param[in] d: factor of the diagonal
 */
static void setDiagonallyDominantValues(mfem::SparseMatrix& A,
                                        const mfem_mgis::real d) {
  const auto* const I = A.HostReadI();
  const auto* const J = A.HostReadJ();
  auto* const v = A.HostWriteData();
  for (int r = 0; r != A.Height(); ++r) {
    for (auto k = I[r]; k != I[r + 1]; ++k) {
      v[k] = (J[k] == r) ? d * (100 + r % 7) : -(1 + (r + J[k]) % 3) / 10.;
    }
  }
}  // end of setDiagonallyDominantValues

/*!
 * \brief checks that the products and the smoothers of hypre applied to the
 * global matrix after a second assembly use its new values: same results
 * as with the product of hypre, a new matrix. On a GPU, MFEM and hypre keep
 * in a matrix data computed from its values: the transposed blocks of the
 * transposed product, the copies used by the triangular solves of the
 * Gauss-Seidel smoother. Each operation is applied after the first
 * assembly, so that these data exist when the values change.
 * \param[in] mesh: mesh
 */
static void testOperationsAfterAssembly(mfem::ParMesh& mesh) {
  using namespace mfem_mgis;
  using Smoother = mfem::HypreSmoother;
  const auto name = std::string{"operations after a second assembly"};
  auto ctx = Context{};
  auto fec = mfem::H1_FECollection(1, mesh.Dimension());
  auto fes = mfem::ParFiniteElementSpace(&mesh, &fec, 3);
  auto dofs = mfem::Array<int>{};
  check(batch::computeElementsDegreesOfFreedom(ctx, dofs, fes),
        name + ": degrees of freedom of the elements");
  const auto A = batch::makeAssembledStiffnessMatrix(ctx, fes, dofs);
  check(A != nullptr, name + ": local matrix");
  if (A == nullptr) {
    return;
  }
  auto b = mfem::Vector(fes.GetTrueVSize());
  b.UseDevice(true);
  auto* const bv = b.HostWrite();
  for (int i = 0; i != b.Size(); ++i) {
    bv[i] = std::cos(0.11 * static_cast<real>(fes.GetMyTDofOffset() + i));
  }
  // an operation applied to `b`. The Chebyshev smoother is not among them:
  // on a GPU, its estimate of the largest eigenvalue starts from random
  // values which change at each call
  struct Operation {
    //! \brief name
    const char* name;
    //! \brief computes the result of the operation for a matrix
    std::function<void(const mfem::HypreParMatrix&, mfem::Vector&)> apply;
  };
  // two sweeps of a smoother of hypre from zero
  const auto smooth = [&b](const Smoother::Type type) {
    return [&b, type](const mfem::HypreParMatrix& K, mfem::Vector& x) {
      x = real{0};
      auto smoother = Smoother(K, type, 2);
      smoother.Mult(b, x);
    };
  };
  const auto operations = std::vector<Operation>{
      {"product",
       [&b](const mfem::HypreParMatrix& K, mfem::Vector& x) { K.Mult(b, x); }},
      {"transposed product", [&b](const mfem::HypreParMatrix& K,
                                  mfem::Vector& x) { K.MultTranspose(b, x); }},
      {"Gauss-Seidel smoother", smooth(Smoother::GS)},
      {"l1 Gauss-Seidel smoother", smooth(Smoother::l1GS)},
      {"l1-Jacobi smoother", smooth(Smoother::l1Jacobi)}};
  // results of all the operations for a matrix
  const auto apply = [&b, &operations](const mfem::HypreParMatrix& K) {
    auto results = std::vector<mfem::Vector>(operations.size());
    for (std::size_t i = 0; i != operations.size(); ++i) {
      results[i].UseDevice(true);
      results[i].SetSize(b.Size());
      operations[i].apply(K, results[i]);
    }
    return results;
  };
  // maximum relative difference between two vectors
  const auto compare = [](const mfem::Vector& x, const mfem::Vector& y) {
    auto ld = real{};
    auto lm = real{};
    const auto* const xv = x.HostRead();
    const auto* const yv = y.HostRead();
    for (int i = 0; i != x.Size(); ++i) {
      ld = std::max(ld, std::abs(xv[i] - yv[i]));
      lm = std::max(lm, std::abs(yv[i]));
    }
    auto gd = real{};
    auto gm = real{};
    const auto mpi_real = mfem::MPITypeMap<real>::mpi_type;
    MPI_Allreduce(&ld, &gd, 1, mpi_real, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lm, &gm, 1, mpi_real, MPI_MAX, MPI_COMM_WORLD);
    return gd / gm;
  };
  auto g = batch::GlobalStiffnessMatrix(fes);
  setDiagonallyDominantValues(*A, 1);
  auto* K = g.assemble(ctx, *A);
  check(K != nullptr, name + ": first assembly");
  if (K == nullptr) {
    return;
  }
  const auto first = apply(*K);
  setDiagonallyDominantValues(*A, 2);
  K = g.assemble(ctx, *A);
  check(K != nullptr, name + ": second assembly");
  if (K == nullptr) {
    return;
  }
  const auto second = apply(*K);
  auto dA = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
  auto Ph = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
  auto R = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
  auto copy = mfem::SparseMatrix(*A);
  dA.MakeSquareBlockDiag(fes.GetComm(), fes.GlobalVSize(), fes.GetDofOffsets(),
                         &copy);
  Ph.ConvertFrom(fes.Dof_TrueDof_Matrix());
  R.MakePtAP(dA, Ph);
  const auto reference = apply(*(R.As<mfem::HypreParMatrix>()));
  for (std::size_t i = 0; i != operations.size(); ++i) {
    const auto operation = name + ", " + operations[i].name;
    const auto d = compare(second[i], reference[i]);
    const auto d1 = compare(first[i], reference[i]);
    if (mfem::Mpi::Root()) {
      std::cout << operation << ": maximum relative difference " << d
                << " with the product of hypre (" << d1
                << " for the first values)\n";
    }
    check(d < mfem_mgis::unit_tests::tolerance(1e-12, 1e-5),
          operation + ": values of the first assembly used");
    check(d1 > 0.1, operation + ": the two assemblies give close results");
  }
}  // end of testOperationsAfterAssembly

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh_file = "cube-two-materials.mesh";
  auto args = mfem::OptionsParser(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "mesh");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto serial = mfem::Mesh(mesh_file, 1, 1);
  serial.UniformRefinement();
  auto mesh = mfem::ParMesh(MPI_COMM_WORLD, serial);
  for (const auto order : {1, 2}) {
    for (const auto vdim : {1, 3}) {
      testCase(mesh, order, vdim,
               "order " + std::to_string(order) + ", " + std::to_string(vdim) +
                   " components");
    }
  }
  testOperationsAfterAssembly(mesh);
  // non conforming space: the prolongation has several entries per row
  {
    auto nc = mfem::Mesh(mesh_file, 1, 1);
    nc.EnsureNCMesh();
    auto refined = mfem::Array<int>{};
    refined.Append(0);
    nc.GeneralRefinement(refined);
    auto pnc = mfem::ParMesh(MPI_COMM_WORLD, nc);
    auto fec = mfem::H1_FECollection(1, pnc.Dimension());
    auto fes = mfem::ParFiniteElementSpace(&pnc, &fec, 3);
    check(!batch::GlobalStiffnessMatrix::isSupported(fes),
          "non conforming space not refused");
  }
  auto all = 0;
  MPI_Allreduce(&failures, &all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (mfem::Mpi::Root()) {
    std::cout << (all == 0 ? "success" : "failure") << '\n';
  }
  return all == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
