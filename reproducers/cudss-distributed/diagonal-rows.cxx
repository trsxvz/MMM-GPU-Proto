// Distributed factorization by cuDSS (mfem::CuDSSSolver of MFEM 4.10, cuDSS
// 0.7.1) of a matrix whose rows are tridiagonal, except the second half of
// the rows of the last process, which only have a diagonal coefficient (as
// the rows of imposed degrees of freedom when the elimination removes their
// other coefficients).
//
// usage: mpirun -n <p> diagonal-rows <coupled> <symmetric> <n>
// - coupled: 1 to couple the last row of each process with the first row of
//   the next one, 0 for a block diagonal matrix;
// - symmetric: 1 for a symmetric pattern (the rows next to the diagonal
//   rows are not coupled with them: the diagonal rows are isolated vertices
//   of the graph of the matrix), 0 otherwise;
// - n: number of rows per process.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "mfem.hpp"

int main(int argc, char** argv) {
  mfem::Mpi::Init(argc, argv);
  if (argc < 4) {
    std::cerr << "usage: " << argv[0] << " coupled symmetric n\n";
    return EXIT_FAILURE;
  }
  auto device = mfem::Device("cuda");
  mfem::Hypre::Init();
  const auto coupled = std::atoi(argv[1]) != 0;
  const auto symmetric = std::atoi(argv[2]) != 0;
  const auto n = std::atoi(argv[3]);
  const auto rank = mfem::Mpi::WorldRank();
  const auto last = rank + 1 == mfem::Mpi::WorldSize();
  const auto first = HYPRE_BigInt{n} * rank;
  const auto global_size = HYPRE_BigInt{n} * mfem::Mpi::WorldSize();
  auto is_diagonal_row = [&](const int l) { return last && (l >= n / 2); };
  auto rows = std::vector<int>{0};
  auto columns = std::vector<HYPRE_BigInt>{};
  auto values = std::vector<mfem::real_t>{};
  auto add = [&](const HYPRE_BigInt j, const mfem::real_t v) {
    columns.push_back(j);
    values.push_back(v);
  };
  for (int l = 0; l != n; ++l) {
    const auto i = first + l;
    if (!is_diagonal_row(l) && ((l > 0) || (coupled && (rank > 0)))) {
      add(i - 1, -1);
    }
    add(i, 4);
    const auto next_is_diagonal_row = is_diagonal_row(l + 1);
    if (!is_diagonal_row(l) && ((l + 1 < n) || (coupled && !last)) &&
        !(symmetric && next_is_diagonal_row)) {
      add(i + 1, -1);
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
