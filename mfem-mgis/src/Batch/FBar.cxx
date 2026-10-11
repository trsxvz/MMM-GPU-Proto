/*!
 * \file   src/Batch/FBar.cxx
 * \brief  This file implements the computations of the FBar formulation by
 * kernels of MFEM (see `MFEMMGIS/Batch/FBar.hxx`)
 * \date   11/10/2026
 */

#include <limits>
#include <string>
#include <cstdint>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "mfem/fem/intrules.hpp"
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Batch/InnerForces.hxx"
#include "MFEMMGIS/Batch/FBar.hxx"
#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES
#include "FBarKernels.hxx"
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

namespace mfem_mgis::batch {

  const mfem::IntegrationRule* getElementsCenterIntegrationRule(
      const PartialQuadratureSpace& qs) {
    const auto& elements = qs.getElements();
    if (elements.Size() == 0) {
      return nullptr;
    }
    const auto& fes =
        getFiniteElementSpace(qs.getFiniteElementDiscretization());
    const auto e = elements.HostRead()[0];
    return &(mfem::IntRules.Get(fes.GetFE(e)->GetGeomType(), 0));
  }  // end of getElementsCenterIntegrationRule

  /*!
   * \brief copies the values of an array of the memory of MFEM, on the
   * device if a device is enabled
   * \param[out] d: destination
   * \param[in] s: source
   */
  static void copyFBarArray(mfem::Vector& d, const mfem::Vector& s) {
    d.UseDevice(true);
    d.SetSize(s.Size());
    if (s.Size() != 0) {
      d = s;
    }
  }  // end of copyFBarArray

  void updateFBarData(FBarData& data) {
    auto& b = data.beginning_of_time_step;
    copyFBarArray(b.deformation_gradients, data.deformation_gradients);
    copyFBarArray(b.coefficients, data.coefficients);
    copyFBarArray(b.centers, data.centers);
    copyFBarArray(b.b, data.b);
    b.valid = data.valid;
    b.defined = true;
    data.has_tangent_operator = false;
  }  // end of updateFBarData

  void revertFBarData(FBarData& data) {
    const auto& b = data.beginning_of_time_step;
    if (b.defined) {
      copyFBarArray(data.deformation_gradients, b.deformation_gradients);
      copyFBarArray(data.coefficients, b.coefficients);
      copyFBarArray(data.centers, b.centers);
      copyFBarArray(data.b, b.b);
      data.valid = b.valid;
    } else {
      data.deformation_gradients.SetSize(0);
      data.coefficients.SetSize(0);
      data.centers.SetSize(0);
      data.b.SetSize(0);
      data.valid = true;
    }
    data.has_tangent_operator = false;
  }  // end of revertFBarData

  bool isFBarSupportedByKernels(const GradientsKind k) noexcept {
    return (k == GradientsKind::DEFORMATION_GRADIENT_3D) ||
           (k == GradientsKind::DEFORMATION_GRADIENT_2D);
  }  // end of isFBarSupportedByKernels

#ifdef MFEM_MGIS_HAVE_DEVICE_SOURCES

  /*!
   * \brief resizes an array of the memory of MFEM, used on the device if a
   * device is enabled
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] v: array
   * \param[in] s: size
   */
  static bool resizeFBarArray(Context& ctx,
                              mfem::Vector& v,
                              const std::int64_t s) {
    if (s > std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage("FBar: array too large for MFEM");
    }
    v.UseDevice(true);
    v.SetSize(static_cast<int>(s));
    return true;
  }  // end of resizeFBarArray

  /*!
   * \brief computes the data at the center of the elements and at the
   * integration points (see `applyFBar` and `initializeFBar`)
   * \param[in] d0: derivatives at the center, null for null unknowns
   * \param[in] gradients: if true, the gradients of the material are
   * replaced by `F̄`
   */
  static bool computeFBarData(Context& ctx,
                              FBarData& data,
                              const Material& m,
                              real* const gradients,
                              const bool point_by_point,
                              const UnknownsDerivatives* const d0,
                              const ElementsGeometry& g0,
                              const GradientsKind k,
                              const char* const name) {
    if (!isFBarSupportedByKernels(k)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the FBar formulation requires a deformation gradient in 3D or in "
          "2D");
    }
    const auto& qs = m.getPartialQuadratureSpace();
    const auto& elements = qs.getElements();
    const auto& offsets = qs.getElementOffsets();
    const auto gsize = static_cast<std::int64_t>(getGradientsSize(k));
    const auto np = std::int64_t{qs.getNumberOfIntegrationPoints()};
    const auto vdim = getUnknownsSize(k);
    data.has_tangent_operator = false;
    data.valid = true;
    if ((!resizeFBarArray(ctx, data.deformation_gradients, np * gsize)) ||
        (!resizeFBarArray(ctx, data.coefficients, np)) ||
        (!resizeFBarArray(ctx, data.centers, elements.Size())) ||
        (!resizeFBarArray(ctx, data.b,
                          std::int64_t{elements.Size()} * g0.nd * vdim))) {
      return false;
    }
    if (elements.Size() == 0) {
      return true;
    }
    const auto* const r0 = getElementsCenterIntegrationRule(qs);
    if ((g0.rule != r0) || (g0.nq != 1) || (g0.maps == nullptr) ||
        (g0.jacobians == nullptr) || (g0.dim != getSpaceDimension(k))) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": invalid geometric data of the center of the elements");
    }
    if ((d0 != nullptr) && ((d0->rule != r0) || (d0->nq != 1) ||
                            (d0->vdim != vdim) || (d0->ne != g0.ne))) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": invalid derivatives of the unknowns at the center of the "
          "elements");
    }
    const auto* const r = qs.getUniformIntegrationRule();
    if (r == nullptr) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the elements of the material do not use the same integration "
          "rule");
    }
    const auto nq = r->GetNPoints();
    if (std::int64_t{elements.Size()} * nq != np) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": unexpected number of points");
    }
    const auto use_device = mfem::Device::IsEnabled();
    const auto macroscopic_gradients = m.getMacroscopicGradients();
    // data at the center of the elements
    auto ca = internals::FBarCentersKernelArguments{};
    ca.elements = elements.Read(use_device);
    ca.derivatives = (d0 != nullptr) ? d0->values.Read(use_device) : nullptr;
    ca.jacobians = g0.jacobians->J.Read(use_device);
    ca.determinants = g0.jacobians->detJ.Read(use_device);
    ca.shape_functions_gradients = g0.maps->G.Read(use_device);
    ca.centers = data.centers.Write(use_device);
    ca.b = data.b.Write(use_device);
    ca.number_of_elements = elements.Size();
    ca.nd = g0.nd;
    for (std::int64_t c = 0; c != gsize; ++c) {
      ca.macroscopic_gradients[c] = macroscopic_gradients[c];
    }
    internals::launchFBarCentersKernel(k, ca, use_device);
    // data at the integration points
    auto ga = internals::FBarGradientsKernelArguments{};
    ga.elements = elements.Read(use_device);
    ga.offsets = offsets.Read(use_device);
    ga.centers = data.centers.Read(use_device);
    ga.gradients = gradients;
    ga.deformation_gradients = data.deformation_gradients.Write(use_device);
    ga.coefficients = data.coefficients.Write(use_device);
    ga.number_of_elements = elements.Size();
    ga.nq = nq;
    ga.number_of_points = static_cast<int>(np);
    ga.point_by_point = point_by_point;
    for (std::int64_t c = 0; c != gsize; ++c) {
      ga.macroscopic_gradients[c] = macroscopic_gradients[c];
    }
    internals::launchFBarGradientsKernel(k, ga, use_device);
    return true;
  }  // end of computeFBarData

  /*!
   * \return the arguments of the kernels adding the matrices of rank one,
   * or an error if the tangent operator of the FBar formulation is not
   * computed or was computed from other tangent operator blocks than `K`
   */
  static bool getFBarRankOneKernelArguments(
      Context& ctx,
      internals::FBarRankOneKernelArguments& a,
      const FBarData& data,
      const PartialQuadratureSpace& qs,
      const mgis::FieldView<const real>& K,
      const ElementsGeometry& g,
      const GradientsKind k,
      const bool use_device,
      const char* const name,
      const std::int64_t threads_per_element) {
    const auto n = std::int64_t{g.nd} * getUnknownsSize(k);
    const auto& elements = qs.getElements();
    const auto size = n * elements.Size();
    const auto np = std::int64_t{qs.getNumberOfIntegrationPoints()};
    const auto gsize = static_cast<std::int64_t>(getGradientsSize(k));
    if ((!data.has_tangent_operator) ||
        (static_cast<std::int64_t>(data.a.Size()) != size) ||
        (static_cast<std::int64_t>(data.b.Size()) != size) ||
        (static_cast<std::int64_t>(data.coefficients.Size()) != np) ||
        (static_cast<std::int64_t>(data.tangent_v.Size()) != np * gsize) ||
        (static_cast<std::int64_t>(data.tangent_c.Size()) != np * gsize)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the tangent operator of the FBar formulation is not computed");
    }
    if (K.data() != data.tangent_operator_blocks) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the tangent operator blocks are not the ones from which the "
          "tangent operator of the FBar formulation is computed");
    }
    if (std::int64_t{elements.Size()} * threads_per_element >
        std::numeric_limits<int>::max()) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": too many elements for the kernels");
    }
    a.elements = elements.Read(use_device);
    a.a = data.a.Read(use_device);
    a.b = data.b.Read(use_device);
    a.number_of_elements = elements.Size();
    a.n = static_cast<int>(n);
    return true;
  }  // end of getFBarRankOneKernelArguments

  /*!
   * \return the modification of the tangent operator blocks of the behaviour
   * giving the tangent operator `K'` of the FBar formulation
   * \param[in] data: data of the FBar formulation
   * \param[in] use_device: if true, the arrays are read on the device
   */
  static TangentOperatorModification getFBarTangentOperatorModification(
      const FBarData& data, const bool use_device) {
    return {.coefficients = data.coefficients.Read(use_device),
            .v = data.tangent_v.Read(use_device),
            .c = data.tangent_c.Read(use_device)};
  }  // end of getFBarTangentOperatorModification

#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */

  //! \brief error reported if the kernels are not compiled
  [[maybe_unused]] static bool reportFBarKernelsNotCompiled(
      Context& ctx, const char* const name) {
    return ctx.registerErrorMessage(
        std::string{name} +
        ": the kernels of the FBar formulation are not compiled (MFEM "
        "compiled with HIP)");
  }  // end of reportFBarKernelsNotCompiled

  bool applyFBar(Context& ctx,
                 FBarData& data,
                 Material& m,
                 const UnknownsDerivatives& d0,
                 const ElementsGeometry& g0,
                 const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(data);
    static_cast<void>(m);
    static_cast<void>(d0);
    static_cast<void>(g0);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx, "applyFBar");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      if (m.s1.gradients_stride != getGradientsSize(k)) {
        return ctx.registerErrorMessage(
            "applyFBar: unexpected number of components of the gradients");
      }
      const auto use_device = mfem::Device::IsEnabled();
      auto* const gradients = m.getArrays().readWrite(
          MaterialArrays::StateArray::GRADIENTS, ets, use_device);
      return computeFBarData(
          ctx, data, m, gradients,
          m.getArrays().array_layout == mgis::FieldLayout::POINT_BY_POINT, &d0,
          g0, k, "applyFBar");
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"applyFBar: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage("applyFBar: unknown exception");
    }
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of applyFBar

  bool initializeFBar(Context& ctx,
                      FBarData& data,
                      const Material& m,
                      const ElementsGeometry& g0,
                      const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(data);
    static_cast<void>(m);
    static_cast<void>(g0);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx, "initializeFBar");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      return computeFBarData(ctx, data, m, nullptr, true, nullptr, g0, k,
                             "initializeFBar");
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"initializeFBar: "} +
                                      e.what());
    } catch (...) {
      return ctx.registerErrorMessage("initializeFBar: unknown exception");
    }
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of initializeFBar

  bool computeFBarTangentOperator(Context& ctx,
                                  FBarData& data,
                                  const PartialQuadratureSpace& qs,
                                  const mgis::FieldView<const real>& K,
                                  const ElementsGeometry& g,
                                  const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(data);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx, "computeFBarTangentOperator");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      data.has_tangent_operator = false;
      data.tangent_operator_blocks = nullptr;
      const auto& elements = qs.getElements();
      const auto& offsets = qs.getElementOffsets();
      const auto np = std::int64_t{qs.getNumberOfIntegrationPoints()};
      const auto gsize = static_cast<std::int64_t>(getGradientsSize(k));
      if ((!isFBarSupportedByKernels(k)) ||
          (static_cast<std::int64_t>(data.coefficients.Size()) != np) ||
          (data.centers.Size() != elements.Size())) {
        return ctx.registerErrorMessage(
            "computeFBarTangentOperator: the data of the FBar formulation are "
            "not computed");
      }
      if ((K.getNumberOfPoints() != qs.getNumberOfIntegrationPoints()) ||
          (K.getNumberOfComponents() != getTangentOperatorBlocksSize(k))) {
        return ctx.registerErrorMessage(
            "computeFBarTangentOperator: unexpected number of points or of "
            "components of the tangent operator blocks");
      }
      const auto n = std::int64_t{g.nd} * getUnknownsSize(k);
      if ((!resizeFBarArray(ctx, data.tangent_v, np * gsize)) ||
          (!resizeFBarArray(ctx, data.tangent_c, np * gsize)) ||
          (!resizeFBarArray(ctx, data.pseudo_stresses, np * gsize)) ||
          (!resizeFBarArray(ctx, data.a, n * elements.Size()))) {
        return false;
      }
      data.a = real{0};
      if (elements.Size() != 0) {
        const auto* const r = qs.getUniformIntegrationRule();
        if ((r == nullptr) || (g.rule != r)) {
          return ctx.registerErrorMessage(
              "computeFBarTangentOperator: the elements of the material do "
              "not use the integration rule of the geometric data");
        }
        const auto use_device = mfem::Device::IsEnabled();
        auto a = internals::FBarTangentOperatorKernelArguments{};
        a.elements = elements.Read(use_device);
        a.offsets = offsets.Read(use_device);
        a.centers = data.centers.Read(use_device);
        a.deformation_gradients = data.deformation_gradients.Read(use_device);
        a.coefficients = data.coefficients.Read(use_device);
        a.tangent_operator = K.data();
        a.point_stride = static_cast<long long>(K.getPointStride());
        a.component_stride = static_cast<long long>(K.getComponentStride());
        a.v = data.tangent_v.Write(use_device);
        a.c = data.tangent_c.Write(use_device);
        a.pseudo_stresses = data.pseudo_stresses.Write(use_device);
        a.number_of_elements = elements.Size();
        a.nq = r->GetNPoints();
        internals::launchFBarTangentOperatorKernel(k, a, use_device);
        // vectors `a`: inner forces of the pseudo-stresses
        const auto s = std::span<const real>(
            data.pseudo_stresses.Read(use_device),
            static_cast<std::size_t>(data.pseudo_stresses.Size()));
        if (!addMaterialInnerForces(ctx, data.a, qs,
                                    mgis::makePointByPointFieldView(
                                        s, qs.getNumberOfIntegrationPoints(),
                                        static_cast<size_type>(gsize),
                                        static_cast<size_type>(gsize)),
                                    g, k)) {
          return false;
        }
      }
      data.tangent_operator_blocks = K.data();
      data.has_tangent_operator = true;
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"computeFBarTangentOperator: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "computeFBarTangentOperator: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of computeFBarTangentOperator

  bool addFBarStiffnessMatrices(Context& ctx,
                                mfem::Vector& m,
                                const FBarData& data,
                                const PartialQuadratureSpace& qs,
                                const mgis::FieldView<const real>& K,
                                const ElementsGeometry& g,
                                const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    static_cast<void>(data);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx, "addFBarStiffnessMatrices");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::FBarRankOneKernelArguments{};
      // one thread per entry of the matrices
      const auto n = std::int64_t{g.nd} * getUnknownsSize(k);
      if (!getFBarRankOneKernelArguments(ctx, a, data, qs, K, g, k, use_device,
                                         "addFBarStiffnessMatrices", n * n)) {
        return false;
      }
      if (!addStiffnessMatrices(
              ctx, m, qs, K, g, k,
              getFBarTangentOperatorModification(data, use_device))) {
        return false;
      }
      if (a.number_of_elements != 0) {
        internals::launchFBarRankOneMatricesKernel(a, m.ReadWrite(use_device),
                                                   use_device);
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFBarStiffnessMatrices: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFBarStiffnessMatrices: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addFBarStiffnessMatrices

  bool addFBarStiffnessMatrices(Context& ctx,
                                const AssembledStiffnessMatrix& m,
                                const FBarData& data,
                                const PartialQuadratureSpace& qs,
                                const mgis::FieldView<const real>& K,
                                const ElementsGeometry& g,
                                const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(m);
    static_cast<void>(data);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx, "addFBarStiffnessMatrices");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::FBarRankOneKernelArguments{};
      // one thread per row of the matrices
      if (!getFBarRankOneKernelArguments(
              ctx, a, data, qs, K, g, k, use_device, "addFBarStiffnessMatrices",
              std::int64_t{g.nd} * getUnknownsSize(k))) {
        return false;
      }
      if (!addStiffnessMatrices(
              ctx, m, qs, K, g, k,
              getFBarTangentOperatorModification(data, use_device))) {
        return false;
      }
      if (a.number_of_elements != 0) {
        auto missing = mfem::Array<int>(1);
        missing.HostWrite()[0] = 0;
        auto o = internals::SparseMatrixOutput{};
        o.number_of_rows = m.matrix->Height();
        o.I = m.matrix->ReadI(use_device);
        o.J = m.matrix->ReadJ(use_device);
        o.values = m.matrix->ReadWriteData(use_device);
        o.element_dofs = m.element_dofs->Read(use_device);
        o.missing = missing.ReadWrite(use_device);
        internals::launchFBarRankOneMatricesKernel(a, o, use_device);
        if (missing.HostRead()[0] != 0) {
          return ctx.registerErrorMessage(
              "addFBarStiffnessMatrices: " +
              std::to_string(missing.HostRead()[0]) +
              " entries of the stiffness matrices are not in the pattern of "
              "the assembled matrix");
        }
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFBarStiffnessMatrices: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFBarStiffnessMatrices: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addFBarStiffnessMatrices

  bool addFBarStiffnessMatricesDiagonals(Context& ctx,
                                         mfem::Vector& d,
                                         const FBarData& data,
                                         const PartialQuadratureSpace& qs,
                                         const mgis::FieldView<const real>& K,
                                         const ElementsGeometry& g,
                                         const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(d);
    static_cast<void>(data);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(g);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx,
                                        "addFBarStiffnessMatricesDiagonals");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::FBarRankOneKernelArguments{};
      // one thread per degree of freedom
      if (!getFBarRankOneKernelArguments(
              ctx, a, data, qs, K, g, k, use_device,
              "addFBarStiffnessMatricesDiagonals",
              std::int64_t{g.nd} * getUnknownsSize(k))) {
        return false;
      }
      if (!addStiffnessMatricesDiagonals(
              ctx, d, qs, K, g, k,
              getFBarTangentOperatorModification(data, use_device))) {
        return false;
      }
      if (a.number_of_elements != 0) {
        internals::launchFBarRankOneDiagonalsKernel(a, d.ReadWrite(use_device),
                                                    use_device);
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFBarStiffnessMatricesDiagonals: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFBarStiffnessMatricesDiagonals: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addFBarStiffnessMatricesDiagonals

  bool addFBarTangentOperatorAction(Context& ctx,
                                    mfem::Vector& y,
                                    const FBarData& data,
                                    const PartialQuadratureSpace& qs,
                                    const mgis::FieldView<const real>& K,
                                    const mfem::Vector& x,
                                    const ElementsGeometry& g,
                                    const GradientsKind k) noexcept {
#ifndef MFEM_MGIS_HAVE_DEVICE_SOURCES
    static_cast<void>(y);
    static_cast<void>(data);
    static_cast<void>(qs);
    static_cast<void>(K);
    static_cast<void>(x);
    static_cast<void>(g);
    static_cast<void>(k);
    return reportFBarKernelsNotCompiled(ctx, "addFBarTangentOperatorAction");
#else  /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::FBarRankOneKernelArguments{};
      // one thread per element
      if (!getFBarRankOneKernelArguments(ctx, a, data, qs, K, g, k, use_device,
                                         "addFBarTangentOperatorAction", 1)) {
        return false;
      }
      if (!addTangentOperatorAction(
              ctx, y, qs, K, x, g, k,
              getFBarTangentOperatorModification(data, use_device))) {
        return false;
      }
      if (a.number_of_elements != 0) {
        internals::launchFBarRankOneActionKernel(
            a, x.Read(use_device), y.ReadWrite(use_device), use_device);
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(
          std::string{"addFBarTangentOperatorAction: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage(
          "addFBarTangentOperatorAction: unknown exception");
    }
    return true;
#endif /* MFEM_MGIS_HAVE_DEVICE_SOURCES */
  }  // end of addFBarTangentOperatorAction

}  // end of namespace mfem_mgis::batch
