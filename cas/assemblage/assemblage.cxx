/*!
 * \file   assemblage.cxx
 * \brief  Time of the assembly by kernels of the stiffness matrices of all
 * the elements in the local sparse matrix of the tangent operator, with the
 * positions of the entries in the rows of the matrix computed once
 * (`batch::computeStiffnessMatricesPositions`) and without them (entries
 * searched in their rows), in the same program. Mechanics in small or finite
 * strain (vector unknowns) or heat transfer (scalar unknowns), arbitrary
 * tangent operator blocks, one process or several (`--parallel`), on the
 * host or on a device. The values of the matrix may be written to a file,
 * to compare two versions of mfem-mgis (one process). In parallel, the
 * construction of the global matrix from the local ones is also timed, as
 * made by mfem-mgis (`assembleTangentOperatorByKernels`): block diagonal
 * matrix, product `Pᵀ A P` of hypre, elimination of the degrees of freedom
 * of a clamped end (boundary 1), and, if mfem-mgis provides it, the
 * assembly by `batch::GlobalStiffnessMatrix` (pattern and exchanges kept);
 * the times are the maxima over the processes.
 */

#include <cmath>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include "mfem/general/backends.hpp"
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/linalg/handle.hpp"
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/pfespace.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Batch/ElementsGeometry.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"
#ifdef MFEM_MGIS_HAS_GLOBAL_STIFFNESS_MATRIX
#include "MFEMMGIS/Batch/GlobalStiffnessMatrix.hxx"
#endif /* MFEM_MGIS_HAS_GLOBAL_STIFFNESS_MATRIX */

//! \return the median of a list of times
static double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const auto n = v.size();
  return (n % 2 == 1) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
}  // end of median

/*!
 * \return the maximum of a value over the processes
 * \param[in] v: value
 * \param[in] parallel: if false, the value is returned
 */
static double maximumOverProcesses(const double v, const bool parallel) {
  if (!parallel) {
    return v;
  }
  auto r = v;
  MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return r;
}  // end of maximumOverProcesses

//! \brief waits for the end of the computations launched on a GPU
static void synchronize() {
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                           mfem::Backend::HIP_MASK)) {
    MFEM_DEVICE_SYNC;
  }
}  // end of synchronize

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "beam-tet.mesh";
  const char* library = "src/libBehaviour.so";
  const char* physics = "mechanics";
  const char* values_file = "";
  auto refinements = 3;
  auto order = 1;
  auto repetitions = 10;
  auto tolerance = real{1e-13};
  auto parallel = false;
  mfem::OptionsParser args(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "Mesh file.");
  args.AddOption(&library, "-l", "--library", "Behaviour library.");
  args.AddOption(&physics, "-t", "--physics",
                 "mechanics (IsotropicLinearHardeningPlasticity), "
                 "finite-strain (SaintVenantKirchhoffElasticity) or thermal "
                 "(StationaryNonLinearHeatTransfer).");
  args.AddOption(&values_file, "-ov", "--output-values",
                 "File in which the values of the matrix assembled without "
                 "the positions are written (raw binary values, without "
                 "header, in the order of the entries of the matrix).");
  args.AddOption(&refinements, "-r", "--refinements",
                 "Number of uniform refinements of the mesh.");
  args.AddOption(&order, "-o", "--order", "Order of the finite elements.");
  args.AddOption(&repetitions, "-n", "--repetitions",
                 "Number of assemblies with and without the positions.");
  args.AddOption(&tolerance, "-tol", "--tolerance",
                 "Relative tolerance of the comparison of the two matrices.");
  args.AddOption(&parallel, "-p", "--parallel", "-no-p", "--no-parallel",
                 "Parallel computation, timing also the construction of the "
                 "global matrix.");
  args.Parse();
  if ((!args.Good()) || (repetitions < 1)) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  if (parallel && (!std::string_view{values_file}.empty())) {
    std::cerr << "the values of the matrix are written in sequential only\n";
    return EXIT_FAILURE;
  }
  const auto thermal = std::string_view{physics} == "thermal";
  const auto finite_strain = std::string_view{physics} == "finite-strain";
  if ((!thermal) && (!finite_strain) &&
      (std::string_view{physics} != "mechanics")) {
    std::cerr << "invalid physics '" << physics << "'\n";
    return EXIT_FAILURE;
  }
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem =
      construct<NonLinearEvolutionProblem>(
          ctx, Parameters{{"MeshFileName", mesh},
                          {"FiniteElementFamily", "H1"},
                          {"FiniteElementOrder", order},
                          {"UnknownsSize", thermal ? 1 : 3},
                          {"NumberOfUniformRefinements", refinements},
                          {"Hypothesis", "Tridimensional"},
                          {"Parallel", parallel}}) |
      or_die;
  auto& fed = problem.getFiniteElementDiscretization();
  const auto& fes =
      parallel ? static_cast<const mfem::FiniteElementSpace&>(
                     fed.getFiniteElementSpace<true>())
               : fed.getFiniteElementSpace<false>();
  const auto attributes = parallel ? fed.getMesh<true>().attributes
                                   : fed.getMesh<false>().attributes;
  // behaviour integrators on all the element attributes, arbitrary tangent
  // operator blocks
  auto integrators = std::vector<AbstractBehaviourIntegrator*>{};
  auto geometry = batch::ElementsGeometrySet{};
  for (const auto m : attributes) {
    if (thermal) {
      problem.addBehaviourIntegrator(ctx, "StationaryNonLinearHeatTransfer",
                                     m, library,
                                     "StationaryNonLinearHeatTransfer") |
          or_die;
    } else {
      problem.addBehaviourIntegrator(
          ctx, "Mechanics", m, library,
          finite_strain ? "SaintVenantKirchhoffElasticity"
                        : "IsotropicLinearHardeningPlasticity") |
          or_die;
    }
    auto& mat = problem.getMaterial(ctx, m, 0) | or_die;
    if (mat.n == 0) {
      // no element of this material on this process
      continue;
    }
    mat.synchronizeOnHost();
    for (std::size_t i = 0; i != mat.K.size(); ++i) {
      mat.K[i] = std::cos(real(i));
    }
    const auto* const r =
        mat.getPartialQuadratureSpace().getUniformIntegrationRule();
    if (r == nullptr) {
      std::cerr << "no uniform integration rule\n";
      return EXIT_FAILURE;
    }
    if (!geometry.contains(r)) {
      batch::computeElementsGeometry(ctx, geometry[r], fes, *r) | or_die;
    }
    integrators.push_back(&(problem.getBehaviourIntegrator(ctx, m, 0) | or_die));
  }
  auto dofs = mfem::Array<int>{};
  batch::computeElementsDegreesOfFreedom(ctx, dofs, fes) | or_die;
  auto matrix = batch::makeAssembledStiffnessMatrix(ctx, fes, dofs);
  if (matrix == nullptr) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  using clock = std::chrono::steady_clock;
  const auto elapsed = [](const clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(clock::now() - t0)
        .count();
  };
  // positions computed twice (allocations and first accesses at the first
  // computation), the second computation being timed
  auto positions = mfem::Array<std::uint16_t>{};
  auto positions_time = 0.;
  for (int i = 0; i != 2; ++i) {
    const auto t0 = clock::now();
    batch::computeStiffnessMatricesPositions(ctx, positions, fes, dofs,
                                             *matrix) |
        or_die;
    synchronize();
    positions_time = elapsed(t0);
  }
  const auto no_positions = mfem::Array<std::uint16_t>{};
  // entries searched in their rows with the positions (none expected)
  auto searched_entries = std::size_t{};
  // assembly of the stiffness matrices of all the elements, values set to
  // zero first (as `MultiMaterialNonLinearIntegrator`); the reading of the
  // counters of the entries by `addStiffnessMatricesByKernels` waits for the
  // end of the kernels
  const auto assemble = [&](const mfem::Array<std::uint16_t>& p) {
    const auto t = clock::now();
    batch::zeroSparseMatrixValues(ctx, *matrix) | or_die;
    const auto a = batch::AssembledStiffnessMatrix{
        .matrix = matrix.get(),
        .element_dofs = &dofs,
        .positions = &p,
        .searched_entries = &searched_entries};
    for (auto* const bi : integrators) {
      bi->addStiffnessMatricesByKernels(ctx, a, geometry) | or_die;
    }
    return elapsed(t);
  };
  const auto nnz = matrix->NumNonZeroElems();
  // first assemblies (allocations on the device), and values compared
  assemble(no_positions);
  const auto reference = std::vector<real>(
      matrix->HostReadData(), matrix->HostReadData() + nnz);
  if (!std::string_view{values_file}.empty()) {
    auto f = std::ofstream(values_file, std::ios::binary);
    f.write(reinterpret_cast<const char*>(reference.data()),
            static_cast<std::streamsize>(reference.size() * sizeof(real)));
    if (!f) {
      std::cerr << "can't write the values in '" << values_file << "'\n";
      return EXIT_FAILURE;
    }
  }
  assemble(positions);
  auto difference = real{};
  auto maximum = real{};
  {
    const auto* const v = matrix->HostReadData();
    for (int i = 0; i != nnz; ++i) {
      difference = std::max(difference, std::abs(v[i] - reference[i]));
      maximum = std::max(maximum, std::abs(reference[i]));
    }
  }
  // assemblies with and without the positions, alternated
  auto with = std::vector<double>{};
  auto without = std::vector<double>{};
  for (int i = 0; i != repetitions; ++i) {
    if (i % 2 == 0) {
      without.push_back(assemble(no_positions));
      with.push_back(assemble(positions));
    } else {
      with.push_back(assemble(positions));
      without.push_back(assemble(no_positions));
    }
  }
  // construction of the global matrix, as by
  // `NonLinearEvolutionProblemImplementation::assembleTangentOperatorByKernels`
  auto global = std::vector<double>{};
  auto elimination = std::vector<double>{};
  // assembly by `batch::GlobalStiffnessMatrix`: first assembly (pattern and
  // exchanges), following assemblies, relative difference with the product
  // of hypre
  auto kept_setup = -1.;
  auto kept = std::vector<double>{};
  auto kept_difference = -1.;
  if (parallel) {
    const auto& pfes = fed.getFiniteElementSpace<true>();
    const auto& pmesh = fed.getMesh<true>();
    auto marker = mfem::Array<int>(pmesh.bdr_attributes.Max());
    marker = 0;
    marker[0] = 1;
    auto essential_dofs = mfem::Array<int>{};
    pfes.GetEssentialTrueDofs(marker, essential_dofs);
    auto* const P = pfes.Dof_TrueDof_Matrix();
    for (int i = 0; i != repetitions; ++i) {
      assemble(positions);
      MPI_Barrier(MPI_COMM_WORLD);
      const auto t0 = clock::now();
      auto A = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
      auto dA = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
      auto Ph = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
      dA.MakeSquareBlockDiag(pfes.GetComm(), pfes.GlobalVSize(),
                             pfes.GetDofOffsets(), matrix.get());
      Ph.ConvertFrom(P);
      A.MakePtAP(dA, Ph);
      synchronize();
      MPI_Barrier(MPI_COMM_WORLD);
      global.push_back(elapsed(t0));
      const auto t1 = clock::now();
      A.As<mfem::HypreParMatrix>()->EliminateBC(essential_dofs,
                                                mfem::Operator::DIAG_ONE);
      synchronize();
      MPI_Barrier(MPI_COMM_WORLD);
      elimination.push_back(elapsed(t1));
    }
#ifdef MFEM_MGIS_HAS_GLOBAL_STIFFNESS_MATRIX
    if (batch::GlobalStiffnessMatrix::isSupported(pfes)) {
      auto g = batch::GlobalStiffnessMatrix(pfes);
      assemble(positions);
      MPI_Barrier(MPI_COMM_WORLD);
      const auto t0 = clock::now();
      g.assemble(ctx, *matrix) | or_die;
      MPI_Barrier(MPI_COMM_WORLD);
      kept_setup = elapsed(t0);
      for (int i = 0; i != repetitions; ++i) {
        assemble(positions);
        MPI_Barrier(MPI_COMM_WORLD);
        const auto t1 = clock::now();
        g.assemble(ctx, *matrix) | or_die;
        MPI_Barrier(MPI_COMM_WORLD);
        kept.push_back(elapsed(t1));
      }
      // comparison with the product of hypre, by products with a vector
      auto dA = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
      auto Ph = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
      auto R = mfem::OperatorHandle(mfem::Operator::Hypre_ParCSR);
      dA.MakeSquareBlockDiag(pfes.GetComm(), pfes.GlobalVSize(),
                             pfes.GetDofOffsets(), matrix.get());
      Ph.ConvertFrom(P);
      R.MakePtAP(dA, Ph);
      auto* const K = g.assemble(ctx, *matrix) | or_die;
      auto x = mfem::Vector(pfes.GetTrueVSize());
      auto* const xv = x.HostWrite();
      for (int i = 0; i != x.Size(); ++i) {
        xv[i] = std::cos(0.3 * i);
      }
      auto y = mfem::Vector(pfes.GetTrueVSize());
      auto yr = mfem::Vector(pfes.GetTrueVSize());
      K->Mult(x, y);
      R.As<mfem::HypreParMatrix>()->Mult(x, yr);
      auto d = 0.;
      auto m = 0.;
      for (int i = 0; i != y.Size(); ++i) {
        d = std::max(d, std::abs(y(i) - yr(i)));
        m = std::max(m, std::abs(yr(i)));
      }
      kept_difference = maximumOverProcesses(d, true) /
                        maximumOverProcesses(m, true);
    }
#endif /* MFEM_MGIS_HAS_GLOBAL_STIFFNESS_MATRIX */
  }
  const auto ne = fes.GetNE();
  const auto n = (ne == 0) ? 0 : dofs.Size() / ne;
  const auto without_time = maximumOverProcesses(median(without), parallel);
  const auto with_time = maximumOverProcesses(median(with), parallel);
  // checks made on all the processes
  const auto failed = maximumOverProcesses(
      ((positions.Size() == 0) || (searched_entries != 0) ||
       (!(maximum > 0)) || (difference > tolerance * maximum))
          ? 1
          : 0,
      parallel);
  auto rank = 0;
  auto size = 1;
  if (parallel) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }
  if (rank != 0) {
    return failed > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
  }
  if (parallel) {
    std::cout << std::setprecision(4) << size
              << " processes, maxima over the processes, first process:\n"
              << "assembly without positions: median " << without_time
              << " ms\n"
              << "assembly with positions: median " << with_time << " ms\n"
              << "block diagonal matrix and product Pt A P: median "
              << median(global) << " ms\n"
              << "elimination of the essential degrees of freedom: median "
              << median(elimination) << " ms\n";
    if (!kept.empty()) {
      std::cout << "global matrix with the pattern and the exchanges kept: "
                << "first assembly " << kept_setup << " ms, median "
                << median(kept) << " ms, relative difference with the "
                << "product of hypre " << kept_difference << '\n';
    }
  }
  std::cout << std::setprecision(4) << "elements " << ne
            << ", degrees of freedom per element " << n << ", entries " << nnz
            << ", positions " << positions.Size() * sizeof(std::uint16_t) / 1e6
            << " MB (" << positions_time << " ms)\n"
            << "assembly without positions: median "
            << median(without) << " ms, minimum "
            << *std::min_element(without.begin(), without.end()) << " ms\n"
            << "assembly with positions: median " << median(with)
            << " ms, minimum " << *std::min_element(with.begin(), with.end())
            << " ms\n"
            << "entries searched with the positions " << searched_entries
            << '\n'
            << "ratio of the medians " << median(without) / median(with)
            << ", maximum difference " << difference << " (maximum value "
            << maximum << ")\n";
  if (failed > 0) {
    std::cerr << "assemblage: no positions, entries searched or different "
                 "matrices\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}  // end of main
