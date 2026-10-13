/*!
 * \file   tests/DirectSolverReuseAnalysisTest.cxx
 * \brief  Test of the `ReuseAnalysis` parameter of the direct solvers of the
 * factory of linear solvers, `MUMPSSolver` (by default) or `CuDSSSolver`
 * (option `--solver`, refused without a device on a CUDA GPU, and when
 * hypre computes on the host on a GPU which does not read the pageable
 * memory of the host): a sequence
 * of operators whose pattern is
 * kept or changes, each one being destroyed after the next one is given,
 * as the non linear solvers do. Each resolution is checked by its
 * residual, with and without the reuse of the analysis, for symmetric and
 * non symmetric matrices, and, with the reuse, the number of analyses made
 * by the solver is checked after each operator: a new analysis is expected if
 * and only if the pattern changed on one of the processes, including the
 * case where only the global indices of the off-diagonal columns change on
 * all processes.
 * \date   12/10/2026
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <optional>
#include <iostream>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "MFEMMGIS/Config.hxx"
#include "Tolerance.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

#ifdef MFEM_USE_CUDA
#include <cuda_runtime_api.h>
#endif /* MFEM_USE_CUDA */

//! \brief number of failed checks
static int failures = 0;

/*!
 * \brief tolerance on the relative residual of the resolutions. In float,
 * it reaches 9.7e-8
 */
static constexpr auto direct_solver_residual_tolerance =
    mfem_mgis::unit_tests::tolerance(1e-10, 1e-6);

//! \return true if the cuDSS solver is expected to be refused
static bool isCuDSSSolverRefused() {
  if (!mfem::Device::Allows(mfem::Backend::CUDA_MASK)) {
    return true;
  }
  if (mfem_mgis::isHypreUsingGPU()) {
    return false;
  }
#ifdef MFEM_USE_CUDA
  // when hypre computes on the host, the GPU must read the pageable memory
  // of the host
  auto device = int{};
  auto pageable_memory_access = int{};
  if ((cudaGetDevice(&device) != cudaSuccess) ||
      (cudaDeviceGetAttribute(&pageable_memory_access,
                              cudaDevAttrPageableMemoryAccess,
                              device) != cudaSuccess)) {
    return true;
  }
  return pageable_memory_access == 0;
#else  /* MFEM_USE_CUDA */
  return true;
#endif /* MFEM_USE_CUDA */
}  // end of isCuDSSSolverRefused

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \brief an operator of the sequence and the matrix holding the local block
 * of the block diagonal operators
 */
struct TestOperator {
  //! \brief local block of a block diagonal operator, used by `A`
  std::unique_ptr<mfem::SparseMatrix> local_block;
  //! \brief operator, destroyed before `local_block`
  std::unique_ptr<mfem::HypreParMatrix> A;
};

/*!
 * \return the matrix of a diffusion problem on a cube, with homogeneous
 * Dirichlet conditions
 * \param[in] fespace: finite element space
 * \param[in] ess_tdofs: degrees of freedom of the Dirichlet conditions
 * \param[in] k: conductivity
 */
static TestOperator assemble(mfem_mgis::FiniteElementSpace<true>& fespace,
                             const mfem::Array<int>& ess_tdofs,
                             const mfem_mgis::real k) {
  auto c = mfem::ConstantCoefficient(k);
  auto a = mfem_mgis::BilinearForm<true>(&fespace);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(c));
  a.Assemble();
  a.Finalize();
  auto A = std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  delete A->EliminateRowsCols(ess_tdofs);
  return {nullptr, std::move(A)};
}  // end of assemble

/*!
 * \return a block diagonal operator whose local block is the local
 * diagonal block of the given matrix, multiplied by `k`, or, on the last
 * process if `diagonal_on_last_process` is true, a diagonal matrix
 * \param[in] A: matrix
 * \param[in] k: factor
 * \param[in] diagonal_on_last_process: diagonal block on the last process
 */
static TestOperator makeBlockDiagonal(const mfem::HypreParMatrix& A,
                                      const mfem_mgis::real k,
                                      const bool diagonal_on_last_process) {
  auto r = TestOperator{};
  const auto last = mfem::Mpi::WorldRank() + 1 == mfem::Mpi::WorldSize();
  if (diagonal_on_last_process && last) {
    auto d = mfem::Vector(A.Height());
    d = 2 * k;
    r.local_block = std::make_unique<mfem::SparseMatrix>(d);
  } else {
    auto d = mfem::SparseMatrix{};
    A.GetDiag(d);
    r.local_block = std::make_unique<mfem::SparseMatrix>(d);
    *(r.local_block) *= k;
  }
  r.A = std::make_unique<mfem::HypreParMatrix>(
      A.GetComm(), A.GetGlobalNumRows(), A.GetRowStarts(), r.local_block.get());
  return r;
}  // end of makeBlockDiagonal

//! \brief number of rows of a process, given its rank
using RowsPerProcess = int (*)(const int);

//! \return four rows for all processes
static int fourRows(const int) { return 4; }

//! \return the global index of the first row of a process
static HYPRE_BigInt getFirstRow(const RowsPerProcess rows, const int rank) {
  auto first = HYPRE_BigInt{};
  for (int r = 0; r != rank; ++r) {
    first += rows(r);
  }
  return first;
}  // end of getFirstRow

/*!
 * \return a matrix given by the global column index of the off-diagonal
 * coefficient of each row (the diagonal coefficients are equal to 4 and the
 * off-diagonal ones to -1)
 * \param[in] coupled: global column index of the off-diagonal coefficient
 * of a row, given its global index
 * \param[in] rows: number of rows of each process
 */
template <typename CoupledColumn>
static TestOperator makeCouplingMatrix(const CoupledColumn& coupled,
                                       const RowsPerProcess rows = fourRows) {
  const auto rank = mfem::Mpi::WorldRank();
  const auto n = rows(rank);
  const auto first = getFirstRow(rows, rank);
  const auto size = getFirstRow(rows, mfem::Mpi::WorldSize());
  auto offsets = std::vector<int>{0};
  auto columns = std::vector<HYPRE_BigInt>{};
  auto values = std::vector<mfem_mgis::real>{};
  for (int l = 0; l != n; ++l) {
    columns.push_back(first + l);
    values.push_back(4);
    columns.push_back(coupled(first + l));
    values.push_back(-1);
    offsets.push_back(static_cast<int>(columns.size()));
  }
  const auto partition = std::array<HYPRE_BigInt, 2>{first, first + n};
  return {nullptr,
          std::make_unique<mfem::HypreParMatrix>(
              MPI_COMM_WORLD, n, size, size, offsets.data(), columns.data(),
              values.data(), partition.data(), partition.data())};
}  // end of makeCouplingMatrix

/*!
 * \return the symmetric matrix `4 I - C - Cᵀ`, where `C` couples each
 * unknown `k` of the first half of the unknowns with the unknown `h + k` of
 * the second half, or, if `swapped` is true, with the unknown `h + (k xor
 * 1)`, `h` being the half of the number of unknowns. On two processes, the
 * two matrices only differ by the column indices of their off-diagonal
 * blocks; on one process, by the column indices of their diagonal blocks.
 * \param[in] swapped: couplings of the swapped pairs of unknowns
 */
static TestOperator makeMatching(const bool swapped) {
  const auto h = HYPRE_BigInt{2} * mfem::Mpi::WorldSize();
  return makeCouplingMatrix([h, swapped](const HYPRE_BigInt i) {
    const auto k = (i < h) ? i : i - h;
    const auto l = swapped ? (k ^ 1) : k;
    return (i < h) ? h + l : l;
  });
}  // end of makeMatching

/*!
 * \return a non symmetric matrix coupling the row `l` of each process with
 * the column `offset + ((l + 1) mod 2)` of the next process (of the same
 * process, outside the diagonal, on one process). The matrices of two
 * offsets only differ by the global indices of their off-diagonal columns
 * on two processes or more; on one process, by the column indices of their
 * diagonal blocks.
 * \param[in] offset: offset of the coupled columns, 0 or 2
 * \param[in] rows: number of rows of each process, two at least
 */
static TestOperator makeNextProcessCoupling(
    const HYPRE_BigInt offset, const RowsPerProcess rows = fourRows) {
  const auto size = mfem::Mpi::WorldSize();
  const auto next = getFirstRow(rows, (mfem::Mpi::WorldRank() + 1) % size);
  return makeCouplingMatrix(
      [next, offset](const HYPRE_BigInt i) {
        return next + offset + ((i + 1) % 2);
      },
      rows);
}  // end of makeNextProcessCoupling

/*!
 * \return the number of rows of a process: four, except on two processes
 * or more, where the last process has six
 */
static int sixRowsOnTheLastProcess(const int rank) {
  const auto size = mfem::Mpi::WorldSize();
  return ((size > 1) && (rank + 1 == size)) ? 6 : 4;
}  // end of sixRowsOnTheLastProcess

/*!
 * \return the number of rows of a process: four, except on three processes
 * or more, where a row of the second process is moved to the last one, the
 * layout of the first process being unchanged (cuDSS 0.7.1 fails to
 * analyse the coupling matrix when two rows are moved)
 */
static int rowsMovedToTheLastProcess(const int rank) {
  const auto size = mfem::Mpi::WorldSize();
  if (size < 3) {
    return 4;
  }
  return (rank == 1) ? 3 : ((rank + 1 == size) ? 5 : 4);
}  // end of rowsMovedToTheLastProcess

//! \brief an operator of the sequence
struct Step {
  //! \brief description
  const char* description;
  //! \brief if true, the operator is symmetric
  bool symmetric;
  /*!
   * \brief if the pattern of the operator differs from the one of the
   * previous operator of the sequence on one process at least, empty if
   * not checked
   */
  std::optional<bool> pattern_changed;
  //! \brief if true, some rows of the operator only have a diagonal entry
  bool diagonal_rows = false;
};

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  // name of the solver
  const char* solver_name = "MUMPSSolver";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&solver_name, "-s", "--solver",
                 "direct solver: MUMPSSolver or CuDSSSolver");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    finalize();
    return EXIT_FAILURE;
  }
  const auto solver_type = std::string{solver_name};
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = Mesh<true>(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fespace = FiniteElementSpace<true>(&pmesh, &fec);
  auto boundaries = mfem::Array<int>(pmesh.bdr_attributes.Max());
  boundaries = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  const auto verbose = fespace.GetMyRank() == 0;
  if ((solver_type == "CuDSSSolver") && (isCuDSSSolverRefused())) {
    auto ctx = Context{};
    const auto h = LinearSolverFactory<true>::getFactory().generate(
        ctx, solver_name, fespace, Parameters{{"ReuseAnalysis", true}});
    check(isInvalid(h), "cuDSS accepted");
    if (verbose) {
      std::cout << "refused: " << ctx.getErrorMessage() << '\n'
                << (failures == 0 ? "success" : "failure") << '\n';
    }
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  const auto reference = assemble(fespace, ess_tdofs, 1);
  // sequence of operators. On one process, the block diagonal operator is
  // the diffusion operator, stored differently: the analysis may be reused
  // or not. On three processes or more, the operator whose rows are moved
  // to the last process does not change the layout of the first process,
  // and is the previous one on one and two processes. The last operator
  // changes the number of rows of the last process on two processes or
  // more (the global number of rows changes on all processes).
  const auto parallel = mfem::Mpi::WorldSize() > 1;
  const auto block_diagonal_pattern_changed =
      parallel ? std::optional<bool>{true} : std::optional<bool>{};
  const auto rows_moved = mfem::Mpi::WorldSize() > 2;
  const auto steps = std::array<Step, 12>{
      Step{"diffusion", true, true},
      Step{"diffusion, values changed", true, false},
      Step{"diffusion, values changed again", true, false},
      Step{"block diagonal", true, block_diagonal_pattern_changed},
      Step{"diagonal on the last process", true, true, true},
      Step{"diagonal on the last process, values changed", true, false, true},
      Step{"matching", true, true},
      Step{"swapped matching", true, true},
      Step{"coupling with the next process", false, true},
      Step{"coupling with the next process, other columns", false, true},
      Step{"coupling with the next process, rows moved to the last process",
           false, rows_moved},
      Step{"coupling with the next process, rows added to the last process",
           false, true}};
  auto make_operator = [&](const std::size_t i) {
    if (i == 0) {
      return assemble(fespace, ess_tdofs, 1);
    } else if (i == 1) {
      return assemble(fespace, ess_tdofs, 3);
    } else if (i == 2) {
      return assemble(fespace, ess_tdofs, 7);
    } else if (i == 3) {
      return makeBlockDiagonal(*(reference.A), 1, false);
    } else if (i == 4) {
      return makeBlockDiagonal(*(reference.A), 1, true);
    } else if (i == 5) {
      return makeBlockDiagonal(*(reference.A), 5, true);
    } else if (i == 6) {
      return makeMatching(false);
    } else if (i == 7) {
      return makeMatching(true);
    } else if (i == 8) {
      return makeNextProcessCoupling(0);
    } else if (i == 9) {
      return makeNextProcessCoupling(2);
    } else if (i == 10) {
      return makeNextProcessCoupling(2, rowsMovedToTheLastProcess);
    }
    return makeNextProcessCoupling(0, sixRowsOnTheLastProcess);
  };
  for (const auto reuse : {false, true}) {
    for (const auto symmetric : {false, true}) {
      auto ctx = Context{};
      auto& factory = LinearSolverFactory<true>::getFactory();
      auto h = factory.generate(ctx, solver_name, fespace,
                                Parameters{{"Symmetric", symmetric},
                                           {"PositiveDefinite", symmetric},
                                           {"ReuseAnalysis", reuse}});
      check(!isInvalid(h),
            solver_type + " is refused: " + ctx.getErrorMessage());
      if (isInvalid(h)) {
        continue;
      }
      const auto& solver = *(h.linear_solver);
      const auto label = std::string{reuse ? "reuse" : "no reuse"} +
                         (symmetric ? ", symmetric" : ", non symmetric");
      check(getNumberOfAnalyses(solver).has_value() == reuse,
            label + ": number of analyses " +
                (reuse ? "not available" : "available"));
      auto op = TestOperator{};
      auto analyses = size_type{};
      for (std::size_t i = 0; i != steps.size(); ++i) {
        if (symmetric && !steps[i].symmetric) {
          continue;
        }
        // the distributed analysis of cuDSS 0.7.1 fails when some rows only
        // have a diagonal entry (observed with two processes on one GPU);
        // the matrices of mfem-mgis keep the entries of the rows of the
        // imposed degrees of freedom
        if ((solver_type == "CuDSSSolver") && parallel &&
            steps[i].diagonal_rows) {
          continue;
        }
        const auto step = label + ", " + steps[i].description;
        // the previous operator is destroyed after the new one is given
        auto previous = make_operator(i);
        std::swap(op, previous);
        h.linear_solver->SetOperator(*(op.A));
        // right-hand side depending on the global indices of the rows
        const auto first_row = op.A->GetRowStarts()[0];
        auto B = mfem::Vector(op.A->Height());
        for (int r = 0; r != B.Size(); ++r) {
          B[r] = static_cast<mfem_mgis::real>(1 + first_row + r);
        }
        auto X = mfem::Vector(B.Size());
        X = 0.0;
        h.linear_solver->Mult(B, X);
        auto R = mfem::Vector(B.Size());
        op.A->Mult(X, R);
        R -= B;
        const auto r = std::sqrt(mfem::InnerProduct(MPI_COMM_WORLD, R, R));
        const auto nb = std::sqrt(mfem::InnerProduct(MPI_COMM_WORLD, B, B));
        check(r < direct_solver_residual_tolerance * nb,
              step + ": relative residual " + std::to_string(r / nb));
        if (verbose) {
          std::cout << step << ": relative residual " << r / nb;
        }
        if (reuse) {
          const auto n = getNumberOfAnalyses(solver).value_or(0);
          if (const auto& changed = steps[i].pattern_changed;
              changed.has_value()) {
            const auto expected = analyses + (*changed ? 1 : 0);
            check(n == expected, step + ": " + std::to_string(n) +
                                     " analyses instead of " +
                                     std::to_string(expected));
          }
          if (verbose) {
            std::cout << ", " << n << " analyses";
          }
          analyses = n;
        }
        if (verbose) {
          std::cout << '\n';
        }
      }
    }
  }
  if (solver_type == "CuDSSSolver") {
    // cuDSS does not report the failures of the factorization: the wrong
    // solution given for a matrix declared positive definite which is not
    // is detected by its residual and reported by hasConverged
    auto ctx = Context{};
    auto h = LinearSolverFactory<true>::getFactory().generate(
        ctx, solver_name, fespace,
        Parameters{{"Symmetric", true}, {"PositiveDefinite", true}});
    check(!isInvalid(h), solver_type + " is refused: " + ctx.getErrorMessage());
    if (!isInvalid(h)) {
      const auto op = makeBlockDiagonal(*(reference.A), -1, false);
      h.linear_solver->SetOperator(*(op.A));
      auto B = mfem::Vector(op.A->Height());
      B = 1.0;
      auto X = mfem::Vector(B.Size());
      X = 0.0;
      h.linear_solver->Mult(B, X);
      check(!hasConverged(*(h.linear_solver)),
            "failure of the factorization not reported");
      if (verbose) {
        std::cout << "negative definite matrix declared positive definite: "
                  << (hasConverged(*(h.linear_solver)) ? "not " : "")
                  << "reported\n";
      }
    }
  }
  // a value which is not a boolean is refused
  auto ctx = Context{};
  check(isInvalid(LinearSolverFactory<true>::getFactory().generate(
            ctx, solver_name, fespace, Parameters{{"ReuseAnalysis", 1}})),
        "invalid value of ReuseAnalysis accepted");
  if (verbose) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
