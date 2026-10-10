/*!
 * \file   src/Batch/StiffnessMatrices.cxx
 * \brief  This file implements the computation of the stiffness matrices of
 * the elements of the standard mechanical behaviour integrators by kernels
 * (host part: checks and arguments of the kernels)
 * \date   09/10/2026
 */

#include <tuple>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/linalg/sparsemat.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/StiffnessMatrices.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "StiffnessMatricesKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES

  /*!
   * \brief checks the arguments of `addStiffnessMatrices` and fills the
   * arguments of the kernels which do not depend on the destination
   * \return true on success, false if there is nothing to do or on error
   * (see `error`)
   * \param[in, out] ctx: execution context
   * \param[out] a: arguments of the kernels
   * \param[out] error: true on error
   * \param[in] qs: quadrature space of the material
   * \param[in] K: tangent operator blocks
   * \param[in] g: geometric data
   * \param[in] k: kind of the gradients
   * \param[in] use_device: if true, the pointers are located on the device
   * \param[in] name: name of the calling function, for the error messages
   * \param[in] diagonals: if true, the kernels compute the diagonals of the
   * stiffness matrices (one thread per node of an element), otherwise the
   * stiffness matrices (one thread per pair of nodes)
   */
  static bool getStiffnessMatricesKernelArguments(
      Context& ctx,
      internals::StiffnessMatricesKernelArguments& a,
      bool& error,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k,
      const bool use_device,
      const std::string& name,
      const bool diagonals) {
    error = false;
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return false;
    }
    error = true;
    if ((g.rule == nullptr) || (g.maps == nullptr) ||
        (g.jacobians == nullptr) ||
        (qs.getUniformIntegrationRule() != g.rule)) {
      return ctx.registerErrorMessage(
          name +
          ": the elements of the material do not use the "
          "integration rule of the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    if (g.dim != dim) {
      return ctx.registerErrorMessage(
          name + ": unexpected dimension of the space (" +
          std::to_string(g.dim) + " instead of " + std::to_string(dim) + ")");
    }
    // one thread per node (diagonals) or per pair of nodes (matrices) of
    // an element
    const auto threads_per_element =
        diagonals ? std::int64_t{g.nd} : std::int64_t{g.nd} * g.nd;
    if (std::int64_t{elements.Size()} * threads_per_element >
        std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage(name +
                                      ": too many elements for the kernels");
    }
    const auto gsize = getGradientsSize(k);
    if ((K.getNumberOfPoints() != qs.getNumberOfIntegrationPoints()) ||
        (K.getNumberOfComponents() != getTangentOperatorBlocksSize(k))) {
      return ctx.registerErrorMessage(
          name +
          ": unexpected number of points or of components "
          "of the tangent operator blocks");
    }
    a.number_of_elements = elements.Size();
    a.nd = g.nd;
    a.nq = g.nq;
    a.gradients_size = static_cast<int>(gsize);
    a.tangent_operator = K.data();
    a.point_stride = static_cast<long long>(K.getPointStride());
    a.component_stride = static_cast<long long>(K.getComponentStride());
    a.elements = elements.Read(use_device);
    a.offsets = offsets.Read(use_device);
    a.shape_functions_gradients = g.maps->G.Read(use_device);
    a.shape_functions = g.maps->B.Read(use_device);
    a.jacobians = g.jacobians->J.Read(use_device);
    a.determinants = g.jacobians->detJ.Read(use_device);
    a.weights = g.weights.Read(use_device);
    error = false;
    return true;
  }  // end of getStiffnessMatricesKernelArguments

#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

  bool addStiffnessMatrices(Context& ctx,
                            mfem::Vector& m,
                            const PartialQuadratureSpace& qs,
                            const mgis::FieldView<const real>& K,
                            const ElementsGeometry& g,
                            const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "addStiffnessMatrices: the kernels computing the stiffness matrices "
        "are not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::StiffnessMatricesKernelArguments{};
      auto error = false;
      if (!getStiffnessMatricesKernelArguments(ctx, a, error, qs, K, g, k,
                                               use_device,
                                               "addStiffnessMatrices", false)) {
        return !error;
      }
      const auto ndofs = std::int64_t{g.nd} * getUnknownsSize(k);
      if (static_cast<std::int64_t>(m.Size()) != ndofs * ndofs * g.ne) {
        return ctx.registerErrorMessage(
            "addStiffnessMatrices: unexpected size of the matrices element by "
            "element");
      }
      a.matrices = m.ReadWrite(use_device);
      internals::launchStiffnessMatricesKernel(k, a, use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"addStiffnessMatrices: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addStiffnessMatrices

  bool addStiffnessMatricesDiagonals(Context& ctx,
                                     mfem::Vector& d,
                                     const PartialQuadratureSpace& qs,
                                     const mgis::FieldView<const real>& K,
                                     const ElementsGeometry& g,
                                     const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(d);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "addStiffnessMatricesDiagonals: the kernels computing the stiffness "
        "matrices are not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::StiffnessMatricesKernelArguments{};
      auto error = false;
      if (!getStiffnessMatricesKernelArguments(
              ctx, a, error, qs, K, g, k, use_device,
              "addStiffnessMatricesDiagonals", true)) {
        return !error;
      }
      const auto size =
          std::int64_t{g.nd} * getUnknownsSize(k) * std::int64_t{g.ne};
      if (static_cast<std::int64_t>(d.Size()) != size) {
        return ctx.registerErrorMessage(
            "addStiffnessMatricesDiagonals: unexpected size of the values "
            "element by element");
      }
      internals::launchStiffnessMatricesDiagonalsKernel(
          k, a, d.ReadWrite(use_device), use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addStiffnessMatricesDiagonals: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addStiffnessMatricesDiagonals: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addStiffnessMatricesDiagonals

  bool addTangentOperatorAction(Context& ctx,
                                mfem::Vector& y,
                                const PartialQuadratureSpace& qs,
                                const mgis::FieldView<const real>& K,
                                const mfem::Vector& x,
                                const ElementsGeometry& g,
                                const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(y);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(x);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "addTangentOperatorAction: the kernels computing the stiffness "
        "matrices are not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::StiffnessMatricesKernelArguments{};
      auto error = false;
      // one thread per element: the check of the number of threads made for
      // the diagonals (one thread per node) is sufficient
      if (!getStiffnessMatricesKernelArguments(
              ctx, a, error, qs, K, g, k, use_device,
              "addTangentOperatorAction", true)) {
        return !error;
      }
      const auto size =
          std::int64_t{g.nd} * getUnknownsSize(k) * std::int64_t{g.ne};
      if ((static_cast<std::int64_t>(x.Size()) != size) ||
          (static_cast<std::int64_t>(y.Size()) != size)) {
        return ctx.registerErrorMessage(
            "addTangentOperatorAction: unexpected size of the values element "
            "by element");
      }
      internals::launchTangentOperatorActionKernel(
          k, a, x.Read(use_device), y.ReadWrite(use_device), use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addTangentOperatorAction: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addTangentOperatorAction: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addTangentOperatorAction

  bool computeElementsDegreesOfFreedom(
      Context& ctx,
      mfem::Array<int>& dofs,
      const mfem::FiniteElementSpace& fes) noexcept {
    try {
      const auto ne = fes.GetNE();
      if (ne == 0) {
        dofs.SetSize(0);
        return true;
      }
      const auto* const fe = fes.GetTypicalFE();
      const auto ndofs =
          static_cast<std::int64_t>(fe->GetDof()) * fes.GetVDim();
      if (ndofs * ne > std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "computeElementsDegreesOfFreedom: too many values for MFEM");
      }
      dofs.SetSize(static_cast<int>(ndofs * ne));
      auto* const d = dofs.HostWrite();
      auto vdofs = mfem::Array<int>{};
      for (int e = 0; e != ne; ++e) {
        fes.GetElementVDofs(e, vdofs);
        if (vdofs.Size() != ndofs) {
          return ctx.registerErrorMessage(
              "computeElementsDegreesOfFreedom: elements with different "
              "numbers of degrees of freedom are not supported");
        }
        for (int i = 0; i != vdofs.Size(); ++i) {
          if (vdofs[i] < 0) {
            return ctx.registerErrorMessage(
                "computeElementsDegreesOfFreedom: degrees of freedom with an "
                "orientation are not supported");
          }
          d[ndofs * e + i] = vdofs[i];
        }
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeElementsDegreesOfFreedom: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeElementsDegreesOfFreedom: unknown exception");
    }
    return true;
  }  // end of computeElementsDegreesOfFreedom

  std::unique_ptr<mfem::SparseMatrix> makeAssembledStiffnessMatrix(
      Context& ctx,
      const mfem::FiniteElementSpace& fes,
      const mfem::Array<int>& dofs) noexcept {
    try {
      const auto n = fes.GetVSize();
      const auto ne = fes.GetNE();
      const auto nde = (ne == 0) ? 0 : dofs.Size() / ne;
      if (nde * ne != dofs.Size()) {
        std::ignore = ctx.registerErrorMessage(
            "makeAssembledStiffnessMatrix: invalid number of degrees of "
            "freedom of the elements");
        return {};
      }
      const auto* const d = dofs.HostRead();
      // elements of each degree of freedom
      auto offsets = std::vector<int>(static_cast<std::size_t>(n) + 1, 0);
      for (int i = 0; i != dofs.Size(); ++i) {
        if ((d[i] < 0) || (d[i] >= n)) {
          std::ignore = ctx.registerErrorMessage(
              "makeAssembledStiffnessMatrix: invalid degree of freedom");
          return {};
        }
        ++offsets[static_cast<std::size_t>(d[i]) + 1];
      }
      std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
      auto elements = std::vector<int>(static_cast<std::size_t>(offsets[n]));
      auto positions = std::vector<int>(offsets.begin(), offsets.end() - 1);
      for (int e = 0; e != ne; ++e) {
        for (int j = 0; j != nde; ++j) {
          const auto r = static_cast<std::size_t>(d[nde * e + j]);
          elements[static_cast<std::size_t>(positions[r]++)] = e;
        }
      }
      // columns of each row: degrees of freedom of the elements of the row,
      // without repetitions, the diagonal first and the others sorted
      auto rows = std::vector<int>(static_cast<std::size_t>(n) + 1, 0);
      auto columns = std::vector<int>{};
      auto last_row = std::vector<int>(static_cast<std::size_t>(n), -1);
      for (int r = 0; r != n; ++r) {
        const auto begin = columns.size();
        for (auto p = offsets[r]; p != offsets[r + 1]; ++p) {
          const auto e = elements[static_cast<std::size_t>(p)];
          for (int j = 0; j != nde; ++j) {
            const auto c = d[nde * e + j];
            if (last_row[static_cast<std::size_t>(c)] != r) {
              last_row[static_cast<std::size_t>(c)] = r;
              columns.push_back(c);
            }
          }
        }
        const auto first = columns.begin() + static_cast<std::ptrdiff_t>(begin);
        std::sort(first, columns.end());
        const auto diagonal = std::lower_bound(first, columns.end(), r);
        if ((diagonal != columns.end()) && (*diagonal == r)) {
          std::rotate(first, diagonal, diagonal + 1);
        }
        if (columns.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
          std::ignore = ctx.registerErrorMessage(
              "makeAssembledStiffnessMatrix: too many values for MFEM");
          return {};
        }
        rows[static_cast<std::size_t>(r) + 1] =
            static_cast<int>(columns.size());
      }
      // arrays allocated by the memory manager of MFEM (as
      // `mfem::Transpose`), owned by the matrix, whose values are allocated
      // and set to zero by its constructor
      const auto nnz = static_cast<int>(columns.size());
      int* const I = mfem::Memory<int>(n + 1);
      int* const J = mfem::Memory<int>(nnz);
      std::copy(rows.begin(), rows.end(), I);
      std::copy(columns.begin(), columns.end(), J);
      return std::make_unique<mfem::SparseMatrix>(I, J, nullptr, n, n, true,
                                                  true, false);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of makeAssembledStiffnessMatrix

  bool zeroSparseMatrixValues(Context& ctx, mfem::SparseMatrix& m) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    return ctx.registerErrorMessage(
        "zeroSparseMatrixValues: the kernels are not compiled (MFEM compiled "
        "with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      if (!m.Finalized()) {
        return ctx.registerErrorMessage(
            "zeroSparseMatrixValues: the matrix is not finalized");
      }
      const auto use_device = mfem::Device::IsEnabled();
      internals::zeroValues(m.WriteData(use_device), m.NumNonZeroElems(),
                            use_device);
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"zeroSparseMatrixValues: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "zeroSparseMatrixValues: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of zeroSparseMatrixValues

  bool addStiffnessMatrices(Context& ctx,
                            const AssembledStiffnessMatrix& m,
                            const PartialQuadratureSpace& qs,
                            const mgis::FieldView<const real>& K,
                            const ElementsGeometry& g,
                            const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return ctx.registerErrorMessage(
        "addStiffnessMatrices: the kernels computing the stiffness matrices "
        "are not compiled (MFEM compiled with HIP)");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      if ((m.matrix == nullptr) || (m.element_dofs == nullptr) ||
          (!m.matrix->Finalized())) {
        return ctx.registerErrorMessage(
            "addStiffnessMatrices: invalid assembled matrix");
      }
      const auto ndofs = std::int64_t{g.nd} * getUnknownsSize(k);
      if (static_cast<std::int64_t>(m.element_dofs->Size()) != ndofs * g.ne) {
        return ctx.registerErrorMessage(
            "addStiffnessMatrices: unexpected number of degrees of freedom "
            "of the elements");
      }
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::StiffnessMatricesKernelArguments{};
      auto error = false;
      if (!getStiffnessMatricesKernelArguments(ctx, a, error, qs, K, g, k,
                                               use_device,
                                               "addStiffnessMatrices", false)) {
        return !error;
      }
      auto missing = mfem::Array<int>(1);
      missing.HostWrite()[0] = 0;
      auto o = internals::SparseMatrixOutput{};
      o.number_of_rows = m.matrix->Height();
      o.I = m.matrix->ReadI(use_device);
      o.J = m.matrix->ReadJ(use_device);
      o.values = m.matrix->ReadWriteData(use_device);
      o.element_dofs = m.element_dofs->Read(use_device);
      o.missing = missing.ReadWrite(use_device);
      internals::launchStiffnessMatricesKernel(k, a, o, use_device);
      if (missing.HostRead()[0] != 0) {
        return ctx.registerErrorMessage(
            "addStiffnessMatrices: " + std::to_string(missing.HostRead()[0]) +
            " entries of the stiffness matrices are not in the pattern of the "
            "assembled matrix");
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"addStiffnessMatrices: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addStiffnessMatrices

}  // end of namespace mfem_mgis::batch
