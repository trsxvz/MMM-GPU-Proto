/*!
 * \file   src/Batch/StiffnessMatrices.cxx
 * \brief  This file implements the computation of the stiffness matrices of
 * the elements of the standard mechanical behaviour integrators by kernels
 * (host part: checks and arguments of the kernels)
 * \date   09/10/2026
 */

#include <limits>
#include <string>
#include <cstdint>
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
   */
  static bool getStiffnessMatricesKernelArguments(
      Context& ctx,
      internals::StiffnessMatricesKernelArguments& a,
      bool& error,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k,
      const bool use_device) {
    error = false;
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    if (elements.Size() == 0) {
      return false;
    }
    error = true;
    if ((g.rule == nullptr) || (g.maps == nullptr) || (g.factors == nullptr) ||
        (qs.getUniformIntegrationRule() != g.rule)) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: the elements of the material do not use the "
          "integration rule of the geometric data");
    }
    const auto dim = getSpaceDimension(k);
    if (g.dim != dim) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unexpected dimension of the space (" +
          std::to_string(g.dim) + " instead of " + std::to_string(dim) + ")");
    }
    const auto gsize = getGradientsSize(k);
    if ((K.getNumberOfPoints() != qs.getNumberOfIntegrationPoints()) ||
        (K.getNumberOfComponents() != gsize * gsize)) {
      return ctx.registerErrorMessage(
          "addStiffnessMatrices: unexpected number of points or of components "
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
    a.jacobians = g.factors->J.Read(use_device);
    a.determinants = g.factors->detJ.Read(use_device);
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
                                               use_device)) {
        return !error;
      }
      const auto ndofs = std::int64_t{g.nd} * getSpaceDimension(k);
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
      const auto ndofs = std::int64_t{g.nd} * getSpaceDimension(k);
      if (static_cast<std::int64_t>(m.element_dofs->Size()) != ndofs * g.ne) {
        return ctx.registerErrorMessage(
            "addStiffnessMatrices: unexpected number of degrees of freedom "
            "of the elements");
      }
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::StiffnessMatricesKernelArguments{};
      auto error = false;
      if (!getStiffnessMatricesKernelArguments(ctx, a, error, qs, K, g, k,
                                               use_device)) {
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
