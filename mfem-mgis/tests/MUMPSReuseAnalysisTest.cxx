/*!
 * \file   tests/MUMPSReuseAnalysisTest.cxx
 * \brief  Test of the `ReuseAnalysis` parameter of the MUMPS solver of the
 * factory of linear solvers: a sequence of operators whose pattern is
 * kept or changes, each one being destroyed after the next one is given,
 * as the non linear solvers do. Each resolution is checked by its
 * residual, with and without the reuse of the analysis, for symmetric and
 * non symmetric matrices, and, with the reuse, the number of analyses made
 * by MUMPS is checked after each operator: a new analysis is expected if
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
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

//! \brief number of failed checks
static int failures = 0;

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

/*!
 * \return a matrix with four rows per process, given by the global column
 * index of the off-diagonal coefficient of each row (the diagonal
 * coefficients are equal to 4 and the off-diagonal ones to -1)
 * \param[in] coupled: global column index of the off-diagonal coefficient
 * of a row, given its global index
 */
template <typename CoupledColumn>
static TestOperator makeCouplingMatrix(const CoupledColumn& coupled) {
  constexpr auto n = int{4};
  const auto first = HYPRE_BigInt{n} * mfem::Mpi::WorldRank();
  const auto size = HYPRE_BigInt{n} * mfem::Mpi::WorldSize();
  auto rows = std::vector<int>{0};
  auto columns = std::vector<HYPRE_BigInt>{};
  auto values = std::vector<mfem_mgis::real>{};
  for (int l = 0; l != n; ++l) {
    columns.push_back(first + l);
    values.push_back(4);
    columns.push_back(coupled(first + l));
    values.push_back(-1);
    rows.push_back(static_cast<int>(columns.size()));
  }
  const auto partition = std::array<HYPRE_BigInt, 2>{first, first + n};
  return {nullptr,
          std::make_unique<mfem::HypreParMatrix>(
              MPI_COMM_WORLD, n, size, size, rows.data(), columns.data(),
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
 */
static TestOperator makeNextProcessCoupling(const HYPRE_BigInt offset) {
  const auto size = mfem::Mpi::WorldSize();
  const auto next = HYPRE_BigInt{4} * ((mfem::Mpi::WorldRank() + 1) % size);
  return makeCouplingMatrix([next, offset](const HYPRE_BigInt i) {
    return next + offset + ((i + 1) % 2);
  });
}  // end of makeNextProcessCoupling

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
};

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = Mesh<true>(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  auto fespace = FiniteElementSpace<true>(&pmesh, &fec);
  auto boundaries = mfem::Array<int>(pmesh.bdr_attributes.Max());
  boundaries = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  const auto verbose = fespace.GetMyRank() == 0;
  const auto reference = assemble(fespace, ess_tdofs, 1);
  // sequence of operators. On one process, the block diagonal operator is
  // the diffusion operator, stored differently: the analysis may be reused
  // or not.
  const auto parallel = mfem::Mpi::WorldSize() > 1;
  const auto block_diagonal_pattern_changed =
      parallel ? std::optional<bool>{true} : std::optional<bool>{};
  const auto steps = std::array<Step, 10>{
      Step{"diffusion", true, true},
      Step{"diffusion, values changed", true, false},
      Step{"diffusion, values changed again", true, false},
      Step{"block diagonal", true, block_diagonal_pattern_changed},
      Step{"diagonal on the last process", true, true},
      Step{"diagonal on the last process, values changed", true, false},
      Step{"matching", true, true},
      Step{"swapped matching", true, true},
      Step{"coupling with the next process", false, true},
      Step{"coupling with the next process, other columns", false, true}};
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
    }
    return makeNextProcessCoupling(2);
  };
  for (const auto reuse : {false, true}) {
    for (const auto symmetric : {false, true}) {
      auto ctx = Context{};
      auto& factory = LinearSolverFactory<true>::getFactory();
      auto h = factory.generate(ctx, "MUMPSSolver", fespace,
                                Parameters{{"Symmetric", symmetric},
                                           {"PositiveDefinite", symmetric},
                                           {"ReuseAnalysis", reuse}});
      check(!isInvalid(h), "MUMPS is refused: " + ctx.getErrorMessage());
      if (isInvalid(h)) {
        continue;
      }
      const auto& solver = *(h.linear_solver);
      const auto label = std::string{reuse ? "reuse" : "no reuse"} +
                         (symmetric ? ", symmetric" : ", non symmetric");
      check(getNumberOfMUMPSAnalyses(solver).has_value() == reuse,
            label + ": number of analyses " +
                (reuse ? "not available" : "available"));
      auto op = TestOperator{};
      auto analyses = size_type{};
      for (std::size_t i = 0; i != steps.size(); ++i) {
        if (symmetric && !steps[i].symmetric) {
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
        check(r < 1e-10 * nb,
              step + ": relative residual " + std::to_string(r / nb));
        if (verbose) {
          std::cout << step << ": relative residual " << r / nb;
        }
        if (reuse) {
          const auto n = getNumberOfMUMPSAnalyses(solver).value_or(0);
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
  // a value which is not a boolean is refused
  auto ctx = Context{};
  check(isInvalid(LinearSolverFactory<true>::getFactory().generate(
            ctx, "MUMPSSolver", fespace, Parameters{{"ReuseAnalysis", 1}})),
        "invalid value of ReuseAnalysis accepted");
  if (verbose) {
    std::cout << (failures == 0 ? "success" : "failure") << '\n';
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
