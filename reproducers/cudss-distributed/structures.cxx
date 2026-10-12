// Distributed factorization by cuDSS (mfem::CuDSSSolver of MFEM 4.10, cuDSS
// 0.7.1) of a matrix whose diagonal coefficients are equal to 4 and whose
// off-diagonal ones are equal to -1, the numbers of rows of the processes
// being given on the command line.
//
// usage: mpirun -n <p> structures <structure> <rows of process 0> ...
// - structure: tridiagonal, or next (the row i is coupled with the row
//   f + 2 + ((i + 1) mod 2), f being the first row of the next process, as
//   the coupling matrices of the test DirectSolverReuseAnalysisTest of
//   mfem-mgis; the next process must have four rows at least).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "mfem.hpp"

int main(int argc, char** argv) {
  mfem::Mpi::Init(argc, argv);
  const auto rank = mfem::Mpi::WorldRank();
  const auto size = mfem::Mpi::WorldSize();
  if (argc != size + 2) {
    std::cerr << "usage: " << argv[0] << " structure <rows of each process>\n";
    return EXIT_FAILURE;
  }
  const auto tridiagonal = std::string{argv[1]} == "tridiagonal";
  auto device = mfem::Device("cuda");
  mfem::Hypre::Init();
  auto first = HYPRE_BigInt{};
  auto global_size = HYPRE_BigInt{};
  auto next = HYPRE_BigInt{};
  for (int r = 0; r != size; ++r) {
    const auto nr = std::atoi(argv[r + 2]);
    if (r == (rank + 1) % size) {
      next = global_size;
    }
    if (r < rank) {
      first += nr;
    }
    global_size += nr;
  }
  const auto n = std::atoi(argv[rank + 2]);
  auto rows = std::vector<int>{0};
  auto columns = std::vector<HYPRE_BigInt>{};
  auto values = std::vector<mfem::real_t>{};
  for (int l = 0; l != n; ++l) {
    const auto i = first + l;
    if (tridiagonal && (i > 0)) {
      columns.push_back(i - 1);
      values.push_back(-1);
    }
    columns.push_back(i);
    values.push_back(4);
    if (tridiagonal && (i + 1 < global_size)) {
      columns.push_back(i + 1);
      values.push_back(-1);
    }
    if (!tridiagonal) {
      columns.push_back(next + 2 + ((i + 1) % 2));
      values.push_back(-1);
    }
    rows.push_back(static_cast<int>(columns.size()));
  }
  HYPRE_BigInt partition[2] = {first, first + n};
  auto A = mfem::HypreParMatrix(MPI_COMM_WORLD, n, global_size, global_size,
                                rows.data(), columns.data(), values.data(),
                                partition, partition);
  auto s = mfem::CuDSSSolver(MPI_COMM_WORLD);
  s.SetOperator(A);
  auto b = mfem::Vector(n);
  b = 1.;
  auto x = mfem::Vector(n);
  x = 0.;
  s.Mult(b, x);
  auto r = mfem::Vector(n);
  A.Mult(x, r);
  r -= b;
  const auto e = std::sqrt(mfem::InnerProduct(MPI_COMM_WORLD, r, r) /
                           mfem::InnerProduct(MPI_COMM_WORLD, b, b));
  if (mfem::Mpi::Root()) {
    std::cout << "relative residual " << e << '\n';
  }
  return EXIT_SUCCESS;
}
