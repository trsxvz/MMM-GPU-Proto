/*!
 * \file   GradientsKernelsTest.cxx
 * \brief  Test of the computations of the gradients of the behaviours, of
 * the inner forces and of the stiffness matrices by kernels
 * (`AbstractBehaviourIntegrator::computeGradientsByKernels`,
 * `addInnerForcesByKernels` and `addStiffnessMatricesByKernels`, the
 * stiffness matrices being stored element by element or assembled in a
 * sparse matrix): same values as the computations element by element and
 * the assembly by MFEM, entries outside of the pattern of the sparse
 * matrix detected, positions of the entries in the rows of the sparse
 * matrix (`batch::computeStiffnessMatricesPositions`) checked and used,
 * including positions of another matrix, wrong positions and positions of
 * a wrong size, assembly row by row in a fixed order
 * (`batch::addStiffnessMatricesByRows`, same values at each assembly), for
 * small and finite strain behaviours, isotropic and
 * orthotropic, in 3D (hexahedra and tetrahedra, orders 1 and 2, two
 * materials with macroscopic gradients) and under the plane strain
 * hypothesis (the plane stress hypothesis uses the same operator, but no
 * behaviour of the tests supports it), also with the FBar formulation and
 * the regularization of Faltus et al. (hexahedra of orders 1 and 2 and
 * quadrangles on distorted meshes, tetrahedra), on the device selected by
 * the `--device` option.
 * \date   09/10/2026
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cstdint>
#include <limits>
#include <optional>
#include <algorithm>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/restriction.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/Batch/Gradients.hxx"
#include "MFEMMGIS/Batch/InnerForces.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

//! \brief description of a case
struct Case {
  //! \brief name of the case
  std::string name;
  //! \brief mesh
  const char* mesh;
  //! \brief modelling hypothesis
  const char* hypothesis;
  //! \brief materials and their behaviours
  std::vector<std::pair<mfem_mgis::size_type, std::string>> materials;
  //! \brief order of the finite elements
  int order = 1;
  //! \brief number of uniform refinements
  int refinement = 1;
  //! \brief if true, non zero macroscopic gradients are given
  bool macroscopic_gradients = false;
  /*!
   * \brief regularization of the finite strain behaviours: `FBar`,
   * `Faltus2026` or none if empty
   */
  std::string regularization = {};
  /*!
   * \brief if true, the vertices of the mesh of the unit cube (square) are
   * moved by a non linear transformation (see `distortMesh`)
   */
  bool distorted_mesh = false;
};

/*!
 * \brief moves the vertices of a mesh of the unit cube or of the unit
 * square by a non linear transformation keeping its faces (sides) in their
 * planes (lines): the elements are not affine and have different shapes
 * \param[in, out] m: mesh
 */
static void distortMesh(mfem::Mesh& m) {
  m.Transform([](const mfem::Vector& x, mfem::Vector& y) {
    using mfem_mgis::real;
    const auto dim = x.Size();
    const auto c = [&x, dim](const int i) { return i < dim ? x[i] : real{0}; };
    y = x;
    y[0] += real(0.1) * c(0) * (1 - c(0)) * (c(1) + real(0.5) * c(2));
    y[1] += real(0.1) * c(1) * (1 - c(1)) * (real(0.3) + c(0) * c(2));
    if (dim == 3) {
      y[2] += real(0.1) * c(2) * (1 - c(2)) * (c(0) + real(0.7) * c(1));
    }
  });
}  // end of distortMesh

/*!
 * \brief checks the positions of the entries of the stiffness matrices of
 * all the elements in the rows of a matrix (see
 * `batch::computeStiffnessMatricesPositions`): the position of an entry is
 * the one of its column in its row, relatively to the first entry of the
 * row, or the largest value of `std::uint16_t` if the entry is not in the
 * matrix
 * \return the number of entries which are not in the matrix
 * \param[in] name: name of the case
 * \param[in] positions: positions
 * \param[in] dofs: degrees of freedom of all the elements
 * \param[in] ne: number of elements
 * \param[in] m: matrix
 */
static int checkPositions(const std::string& name,
                          const mfem::Array<std::uint16_t>& positions,
                          const mfem::Array<int>& dofs,
                          const int ne,
                          const mfem::SparseMatrix& m) {
  constexpr auto unknown = std::numeric_limits<std::uint16_t>::max();
  const auto n = static_cast<std::size_t>(dofs.Size() / ne);
  if (static_cast<std::size_t>(positions.Size()) !=
      n * n * static_cast<std::size_t>(ne)) {
    check(false, name + ": number of positions");
    return 0;
  }
  const auto* const p = positions.HostRead();
  const auto* const d = dofs.HostRead();
  const auto* const I = m.HostReadI();
  const auto* const J = m.HostReadJ();
  auto missing = 0;
  auto wrong = 0;
  for (std::size_t e = 0; e != static_cast<std::size_t>(ne); ++e) {
    for (std::size_t i = 0; i != n; ++i) {
      const auto row = d[n * e + i];
      for (std::size_t j = 0; j != n; ++j) {
        const auto* const b = J + I[row];
        const auto* const f = std::find(b, J + I[row + 1], d[n * e + j]);
        const auto v = p[(n * e + i) * n + j];
        if (f == J + I[row + 1]) {
          ++missing;
          wrong += (v == unknown) ? 0 : 1;
        } else {
          wrong += (v == f - b) ? 0 : 1;
        }
      }
    }
  }
  check(wrong == 0, name + ": " + std::to_string(wrong) + " wrong positions");
  return missing;
}  // end of checkPositions

/*!
 * \brief compare the gradients computed by kernels and element by element
 * \param[in, out] ctx: execution context
 * \param[in] library: behaviour library
 * \param[in] c: case
 */
static void testCase(mfem_mgis::Context& ctx,
                     const char* const library,
                     const Case& c) {
  using namespace mfem_mgis;
  auto or_die = ctx.getFatalFailureHandler();
  const auto tridimensional = std::string{c.hypothesis} == "Tridimensional";
  const auto dim = tridimensional ? 3 : 2;
  auto problem =
      construct<NonLinearEvolutionProblem>(
          ctx, Parameters{{"MeshFileName", c.mesh},
                          {"FiniteElementFamily", "H1"},
                          {"FiniteElementOrder", c.order},
                          {"UnknownsSize", dim},
                          {"NumberOfUniformRefinements", c.refinement},
                          {"Hypothesis", c.hypothesis},
                          {"Parallel", false}}) |
      or_die;
  if (c.distorted_mesh) {
    distortMesh(problem.getFiniteElementDiscretization().getMesh<false>());
  }
  for (const auto& [m, b] : c.materials) {
    if (c.regularization.empty()) {
      problem.addBehaviourIntegrator(ctx, "Mechanics", m, library, b) | or_die;
    } else if (c.regularization == "Faltus2026") {
      // penalization coefficient of the order of the arbitrary tangent
      // operator blocks and thermodynamic forces of the comparisons, the
      // differences of the computations being proportional to it: the
      // computation element by element adds the contributions of `B_q` and
      // of `B0` to the stiffness matrices separately, with cancellation
      // errors (for instance on linear tetrahedra, where `B_q = B0`), and
      // the rounding errors of the deformation gradients computed element
      // by element and by kernels are multiplied by the coefficient in the
      // inner forces
      problem.addBehaviourIntegrator(
          ctx, "Mechanics", m, library, b,
          {{"Regularization",
            Parameters{{"Faltus2026",
                        Parameters{{"PenalizationCoefficient", 10.}}}}}}) |
          or_die;
    } else {
      problem.addBehaviourIntegrator(
          ctx, "Mechanics", m, library, b,
          {{"Regularization", Parameters{{c.regularization, Parameters{}}}}}) |
          or_die;
    }
    auto& mat = problem.getMaterial(ctx, m, 0) | or_die;
    if (mat.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
      // a rotation which is not a symmetry of the material, nor its own
      // inverse
      const auto cs = std::cos(0.3);
      const auto sn = std::sin(0.3);
      const auto r = std::array<real, 9u>{cs, sn, 0, -sn, cs, 0, 0, 0, 1};
      if (tridimensional) {
        mat.setRotationMatrix(RotationMatrix3D{r});
      } else {
        mat.setRotationMatrix(RotationMatrix2D{r});
      }
    }
  }
  auto& fes =
      problem.getFiniteElementDiscretization().getFiniteElementSpace<false>();
  // displacement with a non uniform gradient
  auto gf = mfem::GridFunction(&fes);
  auto f = mfem::VectorFunctionCoefficient(
      dim, [dim](const mfem::Vector& x, mfem::Vector& v) {
        v[0] = 1e-3 * (x[0] * x[0] + 0.3 * x[1]);
        v[1] = 1e-3 * (0.5 * x[0] * x[1] - 0.2 * x[1] * x[1]);
        if (dim == 3) {
          v[0] += 2e-4 * x[2];
          v[1] += 1e-4 * x[0] * x[2];
          v[2] = -1e-3 * (0.2 * x[1] * x[1] + 0.4 * x[2] * x[0]);
        }
      });
  gf.ProjectCoefficient(f);
  const auto& u = static_cast<const mfem::Vector&>(gf);
  auto derivatives = batch::UnknownsDerivativesSet{};
  auto geometry = batch::ElementsGeometrySet{};
  for (const auto& [m, b] : c.materials) {
    auto& bi = problem.getBehaviourIntegrator(ctx, m, 0) | or_die;
    auto& mat = problem.getMaterial(ctx, m, 0) | or_die;
    const auto name = c.name + ", material " + std::to_string(m) + ", " + b;
    check(bi.supportsKernels(), name + ": kernels not supported");
    if (c.macroscopic_gradients) {
      auto g = std::vector<real>(mat.s1.gradients_stride);
      for (std::size_t i = 0; i != g.size(); ++i) {
        g[i] = 1e-3 * (1 + real(i));
      }
      if (mat.b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
        for (std::size_t i = 0; i != 3; ++i) {
          g[i] += 1;
        }
      }
      mat.setMacroscopicGradients(g);
    }
    // gradients computed element by element, on the host
    mat.synchronizeOnHost();
    auto vdofs = mfem::Array<int>{};
    auto ue = mfem::Vector{};
    for (int i = 0; i != fes.GetNE(); ++i) {
      if (fes.GetAttribute(i) != static_cast<int>(m)) {
        continue;
      }
      fes.GetElementVDofs(i, vdofs);
      u.GetSubVector(vdofs, ue);
      check(bi.computeGradients(*(fes.GetFE(i)),
                                *(fes.GetElementTransformation(i)), ue),
            name + ": computation of the gradients element by element");
    }
    const auto reference =
        std::vector<real>(mat.s1.gradients.begin(), mat.s1.gradients.end());
    // gradients computed by kernels
    const auto* const r =
        mat.getPartialQuadratureSpace().getUniformIntegrationRule();
    check(r != nullptr, name + ": no uniform integration rule");
    if (r == nullptr) {
      continue;
    }
    // the integration rule of the material and the additional rules of the
    // regularizations (center of the elements)
    auto rules = bi.getAdditionalIntegrationRulesForKernels();
    rules.push_back(r);
    for (const auto* const ri : rules) {
      if (!derivatives.contains(ri)) {
        check(batch::computeUnknownsDerivatives(ctx, derivatives[ri], fes, u,
                                                *ri),
              name + ": computation of the derivatives of the unknowns (" +
                  ctx.getErrorMessage() + ")");
      }
      if (!geometry.contains(ri)) {
        check(batch::computeElementsGeometry(ctx, geometry[ri], fes, *ri),
              name + ": computation of the geometric data (" +
                  ctx.getErrorMessage() + ")");
      }
    }
    check(bi.computeGradientsByKernels(ctx, derivatives, geometry),
          name + ": computation of the gradients by kernels");
    mat.synchronizeOnHostForReading();
    auto d = real{};
    auto vmax = real{};
    for (std::size_t i = 0; i != reference.size(); ++i) {
      d = std::max(d, std::abs(mat.s1.gradients[i] - reference[i]));
      vmax = std::max(vmax, std::abs(reference[i]));
    }
    std::cout << name << ": " << reference.size()
              << " values, maximum difference " << d << " (maximum value "
              << vmax << ")\n";
    check((!reference.empty()) && (vmax != 0) && (d <= 1e-13 * vmax),
          name + ": gradients");
    // inner forces computed from arbitrary thermodynamic forces
    mat.synchronizeOnHost();
    for (std::size_t i = 0; i != mat.s1.thermodynamic_forces.size(); ++i) {
      mat.s1.thermodynamic_forces[i] = std::sin(real(i));
    }
    const auto& g = geometry[r];
    auto forces = mfem::Vector{};
    forces.UseDevice(true);
    forces.SetSize(g.nd * dim * g.ne);
    forces = real{0};
    check(bi.addInnerForcesByKernels(ctx, forces, geometry),
          name + ": computation of the inner forces by kernels (" +
              ctx.getErrorMessage() + ")");
    const auto* const fv = forces.HostRead();
    // the kernels read the thermodynamic forces on the device: phase of
    // read-only accesses on the host (see `MaterialArrays`)
    mat.synchronizeOnHostForReading();
    auto df = real{};
    auto fmax = real{};
    auto Fe = mfem::Vector{};
    for (int i = 0; i != fes.GetNE(); ++i) {
      if (fes.GetAttribute(i) != static_cast<int>(m)) {
        continue;
      }
      // inner forces of the thermodynamic forces, and the contribution of
      // the regularization of Faltus et al. (see `updateResidual`)
      fes.GetElementVDofs(i, vdofs);
      u.GetSubVector(vdofs, ue);
      bi.updateResidual(Fe, *(fes.GetFE(i)), *(fes.GetElementTransformation(i)),
                        ue);
      check(Fe.Size() == g.nd * dim, name + ": size of the inner forces");
      for (int j = 0; j != Fe.Size(); ++j) {
        const auto v = fv[j + static_cast<std::size_t>(g.nd) * dim * i];
        df = std::max(df, std::abs(v - Fe[j]));
        fmax = std::max(fmax, std::abs(Fe[j]));
      }
    }
    std::cout << name << ": inner forces, maximum difference " << df
              << " (maximum value " << fmax << ")\n";
    check((fmax != 0) && (df <= 1e-13 * fmax), name + ": inner forces");
    // stiffness matrices computed from arbitrary tangent operator blocks
    mat.synchronizeOnHost();
    for (std::size_t i = 0; i != mat.K.size(); ++i) {
      mat.K[i] = std::cos(real(i));
    }
    const auto ndofs = g.nd * dim;
    auto ke = mfem::Vector{};
    ke.UseDevice(true);
    ke.SetSize(ndofs * ndofs * g.ne);
    ke = real{0};
    check(bi.addStiffnessMatricesByKernels(ctx, ke, geometry),
          name + ": computation of the stiffness matrices by kernels (" +
              ctx.getErrorMessage() + ")");
    const auto* const kv = ke.HostRead();
    mat.synchronizeOnHostForReading();
    auto dk = real{};
    auto kmax = real{};
    auto Ke = mfem::DenseMatrix{};
    for (int i = 0; i != fes.GetNE(); ++i) {
      if (fes.GetAttribute(i) != static_cast<int>(m)) {
        continue;
      }
      fes.GetElementVDofs(i, vdofs);
      u.GetSubVector(vdofs, ue);
      bi.updateJacobian(Ke, *(fes.GetFE(i)), *(fes.GetElementTransformation(i)),
                        ue);
      check((Ke.Height() == ndofs) && (Ke.Width() == ndofs),
            name + ": size of the stiffness matrix");
      const auto* const k = Ke.GetData();
      for (int j = 0; j != ndofs * ndofs; ++j) {
        const auto v = kv[j + static_cast<std::size_t>(ndofs) * ndofs * i];
        dk = std::max(dk, std::abs(v - k[j]));
        kmax = std::max(kmax, std::abs(k[j]));
      }
    }
    std::cout << name << ": stiffness matrices, maximum difference " << dk
              << " (maximum value " << kmax << ")\n";
    check((kmax != 0) && (dk <= 1e-13 * kmax), name + ": stiffness matrices");
    // stiffness matrices assembled by kernels in a sparse matrix, compared to
    // the assembly by MFEM of the matrices computed element by element (the
    // pattern of the matrix being the one built by MFEM)
    auto reference_matrix = mfem::SparseMatrix(fes.GetVSize());
    for (int i = 0; i != fes.GetNE(); ++i) {
      if (fes.GetAttribute(i) != static_cast<int>(m)) {
        continue;
      }
      fes.GetElementVDofs(i, vdofs);
      u.GetSubVector(vdofs, ue);
      bi.updateJacobian(Ke, *(fes.GetFE(i)), *(fes.GetElementTransformation(i)),
                        ue);
      reference_matrix.AddSubMatrix(vdofs, vdofs, Ke, 0);
    }
    reference_matrix.Finalize(0);
    // the assembled matrix keeps the ordering of MFEM (search of the
    // entries by bisection, then linearly), and is sorted to be compared to
    // the reference
    auto assembled_matrix = mfem::SparseMatrix(reference_matrix);
    const auto mfem_ordering = mfem::SparseMatrix(reference_matrix);
    reference_matrix.SortColumnIndices();
    auto element_dofs = mfem::Array<int>{};
    check(batch::computeElementsDegreesOfFreedom(ctx, element_dofs, fes),
          name + ": degrees of freedom of the elements (" +
              ctx.getErrorMessage() + ")");
    check(batch::zeroSparseMatrixValues(ctx, assembled_matrix),
          name + ": initialisation of the assembled matrix (" +
              ctx.getErrorMessage() + ")");
    check(bi.addStiffnessMatricesByKernels(
              ctx,
              batch::AssembledStiffnessMatrix{.matrix = &assembled_matrix,
                                              .element_dofs = &element_dofs},
              geometry),
          name + ": assembly of the stiffness matrices by kernels (" +
              ctx.getErrorMessage() + ")");
    assembled_matrix.SortColumnIndices();
    const auto* const av = assembled_matrix.HostReadData();
    const auto* const rv = reference_matrix.HostReadData();
    auto da = real{};
    auto amax = real{};
    for (int i = 0; i != reference_matrix.NumNonZeroElems(); ++i) {
      da = std::max(da, std::abs(av[i] - rv[i]));
      amax = std::max(amax, std::abs(rv[i]));
    }
    std::cout << name << ": assembled matrix, "
              << reference_matrix.NumNonZeroElems()
              << " values, maximum difference " << da << " (maximum value "
              << amax << ")\n";
    check((amax != 0) && (da <= 1e-13 * amax), name + ": assembled matrix");
    // stiffness matrices stored element by element, added row by row in a
    // fixed order (`batch::addStiffnessMatricesByRows`), in the pattern of
    // all the elements: same values as the atomic additions, and the same
    // values at each assembly; rows of another matrix, stiffness matrices of
    // a wrong size and degrees of freedom outside of the matrix refused
    {
      const auto all =
          batch::makeAssembledStiffnessMatrix(ctx, fes, element_dofs);
      auto rows = batch::ElementsRows{};
      check((all != nullptr) &&
                batch::computeElementsRows(ctx, rows, element_dofs, ndofs,
                                           fes.GetVSize()),
            name + ": rows of the stiffness matrices of the elements (" +
                ctx.getErrorMessage() + ")");
      auto sorted = rows.offsets.Size() == fes.GetVSize() + 1;
      for (int r = 0; sorted && (r != fes.GetVSize()); ++r) {
        sorted = std::is_sorted(rows.entries.HostRead() + rows.offsets[r],
                                rows.entries.HostRead() + rows.offsets[r + 1]);
      }
      check(sorted, name +
                        ": rows of the stiffness matrices of the elements "
                        "sorted");
      if (all != nullptr) {
        const auto atomic = [&](mfem::SparseMatrix& matrix) {
          auto c2 = Context{};
          return batch::zeroSparseMatrixValues(c2, matrix) &&
                 bi.addStiffnessMatricesByKernels(
                     c2,
                     batch::AssembledStiffnessMatrix{
                         .matrix = &matrix, .element_dofs = &element_dofs},
                     geometry);
        };
        const auto by_rows = [&](mfem::SparseMatrix& matrix) {
          auto c2 = Context{};
          return batch::zeroSparseMatrixValues(c2, matrix) &&
                 batch::addStiffnessMatricesByRows(
                     c2,
                     batch::AssembledStiffnessMatrix{
                         .matrix = &matrix, .element_dofs = &element_dofs},
                     rows, ke);
        };
        auto a1 = mfem::SparseMatrix(*all);
        auto a2 = mfem::SparseMatrix(*all);
        auto r1 = mfem::SparseMatrix(*all);
        auto r2 = mfem::SparseMatrix(*all);
        check(atomic(a1) && atomic(a2),
              name + ": atomic assembly in the matrix of all the elements");
        check(by_rows(r1) && by_rows(r2), name + ": assembly row by row");
        const auto nnz = all->NumNonZeroElems();
        const auto* const va1 = a1.HostReadData();
        const auto* const va2 = a2.HostReadData();
        const auto* const vr1 = r1.HostReadData();
        const auto* const vr2 = r2.HostReadData();
        auto dr = real{};
        for (int i = 0; i != nnz; ++i) {
          dr = std::max(dr, std::abs(vr1[i] - va1[i]));
        }
        const auto same_atomic = std::equal(va1, va1 + nnz, va2);
        const auto same_rows = std::equal(vr1, vr1 + nnz, vr2);
        std::cout << name << ": assembly row by row, maximum difference to "
                  << "the atomic additions " << dr << ", same values at each "
                  << "assembly: " << (same_rows ? "yes" : "no")
                  << " (atomic additions: " << (same_atomic ? "yes" : "no")
                  << ")\n";
        check(dr <= 1e-13 * amax, name + ": assembly row by row");
        check(same_rows,
              name + ": assembly row by row, same values at each assembly");
        auto c3 = Context{};
        auto other_rows = batch::ElementsRows{};
        check(batch::computeElementsRows(c3, other_rows, element_dofs, ndofs,
                                         fes.GetVSize() + 1) &&
                  !batch::addStiffnessMatricesByRows(
                      c3,
                      batch::AssembledStiffnessMatrix{
                          .matrix = &r1, .element_dofs = &element_dofs},
                      other_rows, ke),
              name + ": rows of another matrix not detected");
        auto small = mfem::Vector(ke.Size() - 1);
        small = real{0};
        check(!batch::addStiffnessMatricesByRows(
                  c3,
                  batch::AssembledStiffnessMatrix{
                      .matrix = &r1, .element_dofs = &element_dofs},
                  rows, small),
              name + ": stiffness matrices of a wrong size not detected");
        check(!batch::computeElementsRows(c3, other_rows, element_dofs, ndofs,
                                          fes.GetVSize() - 1),
              name + ": degree of freedom outside of the matrix not detected");
        // with other materials, the entries of their elements are not in the
        // pattern of the matrix of the elements of the material
        auto other_materials = false;
        for (int i = 0; i != fes.GetNE(); ++i) {
          other_materials =
              other_materials || (fes.GetAttribute(i) != static_cast<int>(m));
        }
        if (other_materials) {
          auto partial = mfem::SparseMatrix(mfem_ordering);
          check(!by_rows(partial),
                name + ": entries outside of the pattern not detected");
        }
      }
    }
    // positions of the entries in the rows of the matrix, read instead of
    // searched: positions of the matrix (no entry searched), positions of
    // the matrix created by `makeAssembledStiffnessMatrix` (other ordering
    // of the rows) and wrong positions (a position is used only if it
    // designates the column of the entry, the other entries being searched
    // and counted), positions of a wrong size (refused). Returns the number
    // of entries searched, or nothing if the assembly failed.
    const auto assemble_with_positions =
        [&](const mfem::Array<std::uint16_t>& positions,
            const std::string& what) -> std::optional<std::size_t> {
      auto matrix = mfem::SparseMatrix(mfem_ordering);
      auto c2 = Context{};
      check(batch::zeroSparseMatrixValues(c2, matrix),
            name + ": initialisation of the assembled matrix");
      auto searched = std::size_t{};
      if (!bi.addStiffnessMatricesByKernels(
              c2,
              batch::AssembledStiffnessMatrix{.matrix = &matrix,
                                              .element_dofs = &element_dofs,
                                              .positions = &positions,
                                              .searched_entries = &searched},
              geometry)) {
        return {};
      }
      matrix.SortColumnIndices();
      const auto* const pv = matrix.HostReadData();
      auto dp = real{};
      for (int i = 0; i != reference_matrix.NumNonZeroElems(); ++i) {
        dp = std::max(dp, std::abs(pv[i] - rv[i]));
      }
      std::cout << name << ": assembled matrix, " << what
                << ", maximum difference " << dp << ", " << searched
                << " entries searched\n";
      check(dp <= 1e-13 * amax, name + ": assembled matrix, " + what);
      return searched;
    };
    auto positions = mfem::Array<std::uint16_t>{};
    check(batch::computeStiffnessMatricesPositions(ctx, positions, fes,
                                                   element_dofs, mfem_ordering),
          name + ": positions of the entries (" + ctx.getErrorMessage() + ")");
    const auto missing = checkPositions(name, positions, element_dofs,
                                        fes.GetNE(), mfem_ordering);
    // the matrix contains the entries of the elements of the material only
    auto other_elements = false;
    for (int i = 0; i != fes.GetNE(); ++i) {
      other_elements =
          other_elements || (fes.GetAttribute(i) != static_cast<int>(m));
    }
    check((missing != 0) == other_elements,
          name + ": positions of the entries outside of the matrix");
    check(assemble_with_positions(positions, "positions of the matrix") ==
              std::size_t{0},
          name + ": assembly with the positions of the matrix");
    {
      const auto other =
          batch::makeAssembledStiffnessMatrix(ctx, fes, element_dofs);
      auto other_positions = mfem::Array<std::uint16_t>{};
      check((other != nullptr) &&
                batch::computeStiffnessMatricesPositions(
                    ctx, other_positions, fes, element_dofs, *other),
            name + ": positions of the entries in another matrix (" +
                ctx.getErrorMessage() + ")");
      const auto searched = assemble_with_positions(
          other_positions, "positions of another matrix");
      check(searched.has_value() && (*searched != 0),
            name + ": assembly with the positions of another matrix");
    }
    {
      auto wrong_positions = mfem::Array<std::uint16_t>(positions.Size());
      auto* const wp = wrong_positions.HostWrite();
      for (int i = 0; i != wrong_positions.Size(); ++i) {
        wp[i] = static_cast<std::uint16_t>((7 * i) % 13);
      }
      const auto searched =
          assemble_with_positions(wrong_positions, "wrong positions");
      check(searched.has_value() && (*searched != 0),
            name + ": assembly with wrong positions");
      auto empty_positions = mfem::Array<std::uint16_t>{};
      check(assemble_with_positions(empty_positions, "no positions") ==
                std::size_t{0},
            name + ": assembly without positions");
      wrong_positions.SetSize(positions.Size() - 1);
      check(!assemble_with_positions(wrong_positions, "positions").has_value(),
            name + ": positions of a wrong size not detected");
    }
    // action and diagonal of the tangent operator computed by kernels
    // without assembly, compared to the product by the assembled matrix and
    // to its diagonal
    {
      const auto* const restriction =
          fes.GetElementRestriction(mfem::ElementDofOrdering::NATIVE);
      auto x = mfem::Vector(fes.GetVSize());
      x.UseDevice(true);
      {
        auto* const xv = x.HostWrite();
        for (int i = 0; i != x.Size(); ++i) {
          xv[i] = std::sin(0.7 * real(i)) + 0.1 * std::cos(real(i));
        }
      }
      auto xe = mfem::Vector(restriction->Height());
      xe.UseDevice(true);
      restriction->Mult(x, xe);
      auto ye = mfem::Vector(restriction->Height());
      ye.UseDevice(true);
      ye = real{0};
      check(bi.addTangentOperatorActionByKernels(ctx, ye, xe, geometry),
            name + ": action of the tangent operator (" +
                ctx.getErrorMessage() + ")");
      auto y = mfem::Vector(fes.GetVSize());
      y.UseDevice(true);
      y = real{0};
      restriction->AddMultTranspose(ye, y);
      auto de = mfem::Vector(restriction->Height());
      de.UseDevice(true);
      de = real{0};
      check(bi.addTangentOperatorDiagonalByKernels(ctx, de, geometry),
            name + ": diagonal of the tangent operator (" +
                ctx.getErrorMessage() + ")");
      auto diag = mfem::Vector(fes.GetVSize());
      diag.UseDevice(true);
      diag = real{0};
      restriction->AddMultTranspose(de, diag);
      auto yref = mfem::Vector(fes.GetVSize());
      reference_matrix.Mult(x, yref);
      auto dref = mfem::Vector{};
      reference_matrix.GetDiag(dref);
      const auto* const yv = y.HostRead();
      const auto* const yr = yref.HostRead();
      const auto* const dv = diag.HostRead();
      const auto* const dr = dref.HostRead();
      auto dy = real{};
      auto ymax = real{};
      auto dd = real{};
      auto dmax = real{};
      for (int i = 0; i != y.Size(); ++i) {
        dy = std::max(dy, std::abs(yv[i] - yr[i]));
        ymax = std::max(ymax, std::abs(yr[i]));
        dd = std::max(dd, std::abs(dv[i] - dr[i]));
        dmax = std::max(dmax, std::abs(dr[i]));
      }
      std::cout << name << ": action of the tangent operator, maximum "
                << "difference " << dy << " (maximum value " << ymax
                << "), diagonal, maximum difference " << dd
                << " (maximum value " << dmax << ")\n";
      check((ymax != 0) && (dy <= 1e-13 * ymax),
            name + ": action of the tangent operator");
      check((dmax != 0) && (dd <= 1e-13 * dmax),
            name + ": diagonal of the tangent operator");
    }
    // entries outside of the pattern of the matrix are detected
    auto diagonal_matrix = mfem::SparseMatrix(fes.GetVSize());
    for (int i = 0; i != fes.GetVSize(); ++i) {
      diagonal_matrix.Set(i, i, real{0});
    }
    diagonal_matrix.Finalize(0);
    auto c2 = Context{};
    check(!bi.addStiffnessMatricesByKernels(
              c2,
              batch::AssembledStiffnessMatrix{.matrix = &diagonal_matrix,
                                              .element_dofs = &element_dofs},
              geometry),
          name + ": entries outside of the pattern of the matrix not detected");
    auto c3 = Context{};
    check(!bi.addStiffnessMatricesByKernels(
              c3,
              batch::AssembledStiffnessMatrix{.matrix = &diagonal_matrix,
                                              .element_dofs = &element_dofs,
                                              .positions = &positions},
              geometry),
          name +
              ": entries outside of the pattern of the matrix not detected "
              "with the positions of another matrix");
  }
  // pattern of the matrix in which the stiffness matrices are assembled: the
  // one built by MFEM from the matrices of all the elements, the diagonal
  // entry being the first of each row
  auto element_dofs = mfem::Array<int>{};
  check(batch::computeElementsDegreesOfFreedom(ctx, element_dofs, fes),
        c.name + ": degrees of freedom of the elements (" +
            ctx.getErrorMessage() + ")");
  const auto m = batch::makeAssembledStiffnessMatrix(ctx, fes, element_dofs);
  check(m != nullptr, c.name + ": creation of the assembled matrix (" +
                          ctx.getErrorMessage() + ")");
  if (m == nullptr) {
    return;
  }
  auto reference = mfem::SparseMatrix(fes.GetVSize());
  auto vdofs = mfem::Array<int>{};
  auto Ke = mfem::DenseMatrix{};
  for (int i = 0; i != fes.GetNE(); ++i) {
    fes.GetElementVDofs(i, vdofs);
    Ke.SetSize(vdofs.Size());
    Ke = real{0};
    reference.AddSubMatrix(vdofs, vdofs, Ke, 0);
  }
  reference.Finalize(0);
  reference.SortColumnIndices();
  const auto nnz = reference.NumNonZeroElems();
  auto same_pattern = (m->Height() == reference.Height()) &&
                      (m->Width() == reference.Width()) &&
                      (m->NumNonZeroElems() == nnz);
  if (same_pattern) {
    const auto* const I = m->HostReadI();
    const auto* const J = m->HostReadJ();
    const auto* const v = m->HostReadData();
    for (int r = 0; r != m->Height(); ++r) {
      same_pattern = same_pattern && (I[r] != I[r + 1]) && (J[I[r]] == r);
    }
    same_pattern = same_pattern &&
                   std::all_of(v, v + nnz, [](const real x) { return x == 0; });
    auto sorted = mfem::SparseMatrix(*m);
    sorted.SortColumnIndices();
    const auto* const sI = sorted.HostReadI();
    const auto* const sJ = sorted.HostReadJ();
    const auto* const rI = reference.HostReadI();
    const auto* const rJ = reference.HostReadJ();
    same_pattern = same_pattern && std::equal(sI, sI + m->Height() + 1, rI) &&
                   std::equal(sJ, sJ + nnz, rJ);
  }
  std::cout << c.name << ": pattern of the assembled matrix, " << nnz
            << " values\n";
  check(same_pattern, c.name + ": pattern of the assembled matrix");
  auto positions = mfem::Array<std::uint16_t>{};
  check(batch::computeStiffnessMatricesPositions(ctx, positions, fes,
                                                 element_dofs, *m),
        c.name + ": positions of the entries (" + ctx.getErrorMessage() + ")");
  check(checkPositions(c.name, positions, element_dofs, fes.GetNE(), *m) == 0,
        c.name + ": positions of the entries of the assembled matrix");
}  // end of testCase

/*!
 * \brief checks how the multi-material integrator keeps the positions of the
 * entries of the stiffness matrices (see
 * `MultiMaterialNonLinearIntegrator::assembleStiffnessMatrixByKernels`):
 * computed once for the matrices created by `makeAssembledStiffnessMatrix`,
 * once for these matrices with sorted columns (as by a direct solver), even
 * when the two kinds of matrices are alternated, computed again after an
 * assembly in another layout of the rows, at most three times, then
 * abandoned; never computed with a device which is not the host. The
 * ordered assembly gives the same values at each assembly, close to the ones
 * of the atomic additions.
 * \param[in, out] ctx: execution context
 * \param[in] library: behaviour library
 * \param[in] mesh: mesh of a cube (hexahedra)
 */
static void testIntegratorPositions(mfem_mgis::Context& ctx,
                                    const char* const library,
                                    const char* const mesh) {
  using namespace mfem_mgis;
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", mesh},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 3},
                                     {"NumberOfUniformRefinements", 1},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", false}}) |
                 or_die;
  auto mi = MultiMaterialNonLinearIntegrator(
      problem.getFiniteElementDiscretizationPointer(),
      Hypothesis::TRIDIMENSIONAL);
  mi.addBehaviourIntegrator(ctx, "Mechanics", 1, library,
                            "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& mat = mi.getMaterial(ctx, 1, 0) | or_die;
  mat.synchronizeOnHost();
  for (std::size_t i = 0; i != mat.K.size(); ++i) {
    mat.K[i] = std::cos(real(i));
  }
  const auto& fes =
      problem.getFiniteElementDiscretization().getFiniteElementSpace<false>();
  const auto m1 = mi.makeAssembledStiffnessMatrix(ctx, fes);
  const auto m2 = mi.makeAssembledStiffnessMatrix(ctx, fes);
  check((m1 != nullptr) && (m2 != nullptr),
        "integrator positions: creation of the matrices");
  if ((m1 == nullptr) || (m2 == nullptr)) {
    return;
  }
  // matrix in the layout of MFEM (neither sorted nor diagonal first)
  auto m3 = mfem::SparseMatrix(fes.GetVSize());
  {
    auto vdofs = mfem::Array<int>{};
    auto Ke = mfem::DenseMatrix{};
    for (int i = 0; i != fes.GetNE(); ++i) {
      fes.GetElementVDofs(i, vdofs);
      Ke.SetSize(vdofs.Size());
      Ke = real{0};
      m3.AddSubMatrix(vdofs, vdofs, Ke, 0);
    }
    m3.Finalize(0);
  }
  const auto host = !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  const auto assemble = [&ctx, &mi, &fes](mfem::SparseMatrix& m) {
    check(mi.assembleStiffnessMatrixByKernels(ctx, fes, m),
          "integrator positions: assembly (" + ctx.getErrorMessage() + ")");
  };
  // expected numbers of computations and positions kept, for the matrices
  // whose columns are not sorted, then sorted
  const auto expect = [&mi, host](const std::string& what, const int n1,
                                  const bool p1, const int n2, const bool p2) {
    const auto c1 = mi.getNumberOfStiffnessMatricesPositionsComputations(false);
    const auto c2 = mi.getNumberOfStiffnessMatricesPositionsComputations(true);
    std::cout << "integrator positions, " << what << ": " << c1 << " and " << c2
              << " computations\n";
    check((c1 == (host ? n1 : 0)) &&
              (mi.hasStiffnessMatricesPositions(false) == (host && p1)) &&
              (c2 == (host ? n2 : 0)) &&
              (mi.hasStiffnessMatricesPositions(true) == (host && p2)),
          "integrator positions: " + what);
  };
  assemble(*m1);
  const auto reference = std::vector<real>(
      m1->HostReadData(), m1->HostReadData() + m1->NumNonZeroElems());
  expect("first assembly", 1, true, 0, false);
  assemble(*m2);
  assemble(*m1);
  expect("two matrices of the same layout", 1, true, 0, false);
  m1->SortColumnIndices();
  for (int i = 0; i != 3; ++i) {
    assemble(*m1);
    assemble(*m2);
  }
  expect("sorted and unsorted matrices alternated", 1, true, 1, true);
  // another layout of the unsorted rows: positions computed again, at most
  // three times, then abandoned
  assemble(m3);
  expect("matrix in another layout", 1, false, 1, true);
  assemble(*m2);
  assemble(m3);
  assemble(*m2);
  expect("matrices in two layouts alternated", 3, true, 1, true);
  assemble(m3);
  assemble(*m2);
  expect("positions abandoned", 3, false, 1, true);
  // the assembled values do not depend on the positions
  const auto* const v = m2->HostReadData();
  auto same =
      static_cast<std::size_t>(m2->NumNonZeroElems()) == reference.size();
  for (std::size_t i = 0; same && (i != reference.size()); ++i) {
    same = std::abs(v[i] - reference[i]) <=
           1e-14 * std::max(std::abs(reference[i]), real{1});
  }
  check(same, "integrator positions: assembled values");
  // ordered assembly (see `setOrderedStiffnessMatricesAssembly`), in the
  // three layouts of the rows: the same values at each assembly, close to
  // the ones of the atomic additions
  check(mi.getNumberOfOrderedStiffnessMatricesAssemblies() == 0,
        "integrator positions: no ordered assembly by default");
  mi.setOrderedStiffnessMatricesAssembly(true);
  for (auto* const m : {m2.get(), &m3, m1.get()}) {
    assemble(*m);
    const auto* const first = m->HostReadData();
    const auto values = std::vector<real>(first, first + m->NumNonZeroElems());
    assemble(*m);
    check(std::equal(values.begin(), values.end(), m->HostReadData()),
          "integrator positions: ordered assembly, same values at each "
          "assembly");
  }
  const auto* const ov = m2->HostReadData();
  auto close = true;
  for (std::size_t i = 0; close && (i != reference.size()); ++i) {
    close = std::abs(ov[i] - reference[i]) <=
            1e-14 * std::max(std::abs(reference[i]), real{1});
  }
  check(close, "integrator positions: ordered assembly, assembled values");
  check(mi.getNumberOfOrderedStiffnessMatricesAssemblies() == 6,
        "integrator positions: number of ordered assemblies");
  mi.setOrderedStiffnessMatricesAssembly(false);
}  // end of testIntegratorPositions

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* cube = "cube.mesh";
  const char* tetrahedra = "beam-tet.mesh";
  const char* two_materials = "cube_2mat_per.mesh";
  const char* square = "cube2d.msh";
  const char* library = "src/libBehaviourTest.so";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&cube, "-m", "--mesh", "mesh of a cube (hexahedra)");
  args.AddOption(&tetrahedra, "-mt", "--tetrahedral-mesh",
                 "tetrahedral mesh (material 1)");
  args.AddOption(&two_materials, "-m2", "--two-materials-mesh",
                 "mesh of a cube made of two materials (1 and 2)");
  args.AddOption(&square, "-ms", "--square-mesh",
                 "two-dimensional mesh of a square (material 1)");
  args.AddOption(&library, "-l", "--library", "behaviour library");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto ctx = Context{};
  const auto plasticity = std::string{"IsotropicLinearHardeningPlasticity"};
  const auto elasticity = std::string{"IsotropicLinearElasticity"};
  const auto svk = std::string{"SaintVenantKirchhoffElasticity"};
  const auto orthotropic = std::string{"OrthotropicElasticity"};
  auto cases = std::vector<Case>{};
  for (const auto order : {1, 2}) {
    const auto o = ", order " + std::to_string(order);
    cases.push_back({.name = "hexahedra" + o,
                     .mesh = cube,
                     .hypothesis = "Tridimensional",
                     .materials = {{1, plasticity}},
                     .order = order});
    cases.push_back({.name = "hexahedra, finite strain" + o,
                     .mesh = cube,
                     .hypothesis = "Tridimensional",
                     .materials = {{1, svk}},
                     .order = order});
    cases.push_back({.name = "tetrahedra" + o,
                     .mesh = tetrahedra,
                     .hypothesis = "Tridimensional",
                     .materials = {{1, orthotropic}},
                     .order = order,
                     .refinement = 0});
    cases.push_back({.name = "plane strain" + o,
                     .mesh = square,
                     .hypothesis = "PlaneStrain",
                     .materials = {{1, plasticity}},
                     .order = order,
                     .refinement = 0});
  }
  cases.push_back({.name = "two materials, macroscopic gradients",
                   .mesh = two_materials,
                   .hypothesis = "Tridimensional",
                   .materials = {{1, plasticity}, {2, elasticity}},
                   .refinement = 0,
                   .macroscopic_gradients = true});
  cases.push_back({.name = "finite strain, macroscopic gradients",
                   .mesh = cube,
                   .hypothesis = "Tridimensional",
                   .materials = {{1, svk}},
                   .macroscopic_gradients = true});
  cases.push_back({.name = "plane strain, finite strain",
                   .mesh = square,
                   .hypothesis = "PlaneStrain",
                   .materials = {{1, svk}},
                   .refinement = 0});
  // regularizations of the finite strain behaviours, on distorted meshes
  // (elements which are not affine and have different shapes)
  for (const auto* const regularization : {"FBar", "Faltus2026"}) {
    const auto rn = std::string{regularization} + ", ";
    for (const auto order : {1, 2}) {
      cases.push_back({.name = rn + "hexahedra, distorted mesh, order " +
                               std::to_string(order),
                       .mesh = cube,
                       .hypothesis = "Tridimensional",
                       .materials = {{1, svk}},
                       .order = order,
                       .regularization = regularization,
                       .distorted_mesh = true});
    }
    // quadratic tetrahedra: on linear tetrahedra, the deformation gradient
    // is uniform in each element, and the regularizations change nothing
    cases.push_back({.name = rn + "tetrahedra, order 2",
                     .mesh = tetrahedra,
                     .hypothesis = "Tridimensional",
                     .materials = {{1, svk}},
                     .order = 2,
                     .refinement = 0,
                     .regularization = regularization});
    cases.push_back({.name = rn + "plane strain, distorted mesh",
                     .mesh = square,
                     .hypothesis = "PlaneStrain",
                     .materials = {{1, svk}},
                     .refinement = 0,
                     .regularization = regularization,
                     .distorted_mesh = true});
  }
  for (const auto& c : cases) {
    testCase(ctx, library, c);
  }
  testIntegratorPositions(ctx, library, cube);
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
