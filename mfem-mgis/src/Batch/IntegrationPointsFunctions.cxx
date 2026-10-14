/*!
 * \file   src/Batch/IntegrationPointsFunctions.cxx
 * \brief  This file implements the computations by kernels on functions
 * defined at the integration points of a material (host part: checks and
 * arguments of the kernels)
 * \date   11/10/2026
 */

#include <limits>
#include <string>
#include <cstdint>
#include <exception>
#include "mfem/general/device.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/Batch/IntegrationPointsFunctions.hxx"
#include "IntegrationPointsFunctionsKernels.hxx"

namespace mfem_mgis::batch {

  bool gatherValues(Context& ctx,
                    mfem::Vector& v,
                    const ImmutablePartialQuadratureFunctionView& f) noexcept {
    try {
      const auto& qs = f.getPartialQuadratureSpace();
      const auto& elements = qs.getElements();
      const auto& offsets = qs.getElementOffsets();
      const auto n = static_cast<std::int64_t>(elements.Size());
      const auto* const ir = qs.getUniformIntegrationRule();
      if ((n != 0) && (ir == nullptr)) {
        return ctx.registerErrorMessage(
            "gatherValues: the elements of the partial quadrature space do "
            "not use the same integration rule");
      }
      const auto nq = (n == 0) ? std::int64_t{0} : ir->GetNPoints();
      const auto nc = static_cast<std::int64_t>(f.getNumberOfComponents());
      if (nc * nq * n > std::numeric_limits<int>::max()) {
        return ctx.registerErrorMessage(
            "gatherValues: too many values for a vector of MFEM");
      }
      v.SetSize(static_cast<int>(nc * nq * n));
      v.UseDevice(true);
      auto* const pv = v.HostWrite();
      const auto* const po = offsets.HostRead();
      for (std::int64_t k = 0; k != n; ++k) {
        for (std::int64_t q = 0; q != nq; ++q) {
          const auto values =
              f.getIntegrationPointValues(static_cast<size_type>(po[k] + q));
          for (std::int64_t c = 0; c != nc; ++c) {
            pv[c + nc * (q + nq * k)] = values[c];
          }
        }
      }
    } catch (std::exception& e) {
      return ctx.registerErrorMessage(std::string{"gatherValues: "} + e.what());
    } catch (...) {
      return ctx.registerErrorMessage("gatherValues: unknown exception");
    }
    return true;
  }  // end of gatherValues

  /*!
   * \brief checks the geometric data of a partial quadrature space and the
   * size of the values of a function, and fills the arguments of the
   * kernels
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] a: arguments of the kernels
   * \param[in] qs: partial quadrature space
   * \param[in] v: values of the function, or null
   * \param[in] nc: number of components of the function
   * \param[in] g: geometric data
   * \param[in] use_device: if true, the data are read on the device
   * \param[in] name: name of the calling function, for the messages
   */
  static bool makeKernelArguments(
      Context& ctx,
      internals::IntegrationPointsFunctionsKernelArguments& a,
      const PartialQuadratureSpace& qs,
      const mfem::Vector* const v,
      const size_type nc,
      const ElementsGeometry& g,
      const bool use_device,
      const char* const name) {
    const auto& elements = qs.getElements();
    if ((g.rule == nullptr) || (g.maps == nullptr) ||
        (g.jacobians == nullptr) ||
        (qs.getUniformIntegrationRule() != g.rule)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": the elements of the partial quadrature space do not use the "
          "integration rule of the geometric data");
    }
    if (nc == 0) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": invalid number of components");
    }
    const auto size =
        std::int64_t{g.nq} * static_cast<std::int64_t>(nc) * elements.Size();
    if ((v != nullptr) && (static_cast<std::int64_t>(v->Size()) != size)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": unexpected size of the values of the function");
    }
    a.elements = elements.Read(use_device);
    a.points.shape_functions = g.maps->B.Read(use_device);
    a.points.determinants = g.jacobians->detJ.Read(use_device);
    a.points.weights = g.weights.Read(use_device);
    a.points.nq = g.nq;
    a.values = (v != nullptr) ? v->Read(use_device) : nullptr;
    a.nc = static_cast<int>(nc);
    a.number_of_elements = elements.Size();
    a.nd = g.nd;
    return true;
  }  // end of makeKernelArguments

  /*!
   * \return an error message reporting the exception being handled (to be
   * called in a `catch` block only)
   * \param[in] name: name of the calling function
   */
  static std::string getExceptionMessage(const char* const name) {
    try {
      throw;
    } catch (std::exception& e) {
      return std::string{name} + ": " + e.what();
    } catch (...) {
    }
    return std::string{name} + ": unknown exception";
  }  // end of getExceptionMessage

  std::optional<real> computeIntegral(Context& ctx,
                                      const PartialQuadratureSpace& qs,
                                      const mfem::Vector* const v,
                                      const ElementsGeometry& g) noexcept {
    if (qs.getElements().Size() == 0) {
      return real{0};
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::IntegrationPointsFunctionsKernelArguments{};
      if (!makeKernelArguments(ctx, a, qs, v, 1, g, use_device,
                               "computeIntegral")) {
        return {};
      }
      // geometric data of the elements of the space only
      a.compact_geometry = (g.ne == qs.getElements().Size());
      auto w = mfem::Array<double>{};
      return static_cast<real>(
          internals::launchIntegralKernel(a, w, use_device));
    } catch (...) {
      return ctx.registerErrorMessage(getExceptionMessage("computeIntegral"));
    }
  }  // end of computeIntegral

  bool addL2ProjectionRightHandSide(Context& ctx,
                                    mfem::Vector& r,
                                    const PartialQuadratureSpace& qs,
                                    const mfem::Vector& v,
                                    const size_type nc,
                                    const size_type c,
                                    const ElementsGeometry& g) noexcept {
    const auto* const name = "addL2ProjectionRightHandSide";
    if (qs.getElements().Size() == 0) {
      return true;
    }
    if (c >= nc) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": invalid component");
    }
    if (static_cast<std::int64_t>(r.Size()) != std::int64_t{g.nd} * g.ne) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": unexpected size of the values element by element");
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::IntegrationPointsFunctionsKernelArguments{};
      if (!makeKernelArguments(ctx, a, qs, &v, nc, g, use_device, name)) {
        return false;
      }
      internals::launchL2ProjectionRightHandSideKernel(
          a, static_cast<int>(c), r.ReadWrite(use_device), use_device);
    } catch (...) {
      return ctx.registerErrorMessage(getExceptionMessage(name));
    }
    return true;
  }  // end of addL2ProjectionRightHandSide

  bool computeLocalProjections(Context& ctx,
                               mfem::Vector& r,
                               const PartialQuadratureSpace& qs,
                               const mfem::Vector& v,
                               const size_type nc,
                               const ElementsGeometry& g,
                               const ElementsGeometry& m) noexcept {
    const auto* const name = "computeLocalProjections";
    if (qs.getElements().Size() == 0) {
      return true;
    }
    if ((m.rule == nullptr) || (m.maps == nullptr) ||
        (m.jacobians == nullptr) || (m.nd != g.nd) || (m.ne != g.ne)) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": inconsistent geometric data of the mass matrices");
    }
    if (g.nd > maximumNumberOfNodesOfLocalProjections) {
      return ctx.registerErrorMessage(
          std::string{name} + ": too many degrees of freedom per element (" +
          std::to_string(g.nd) + ", at most " +
          std::to_string(maximumNumberOfNodesOfLocalProjections) + ")");
    }
    if (static_cast<std::int64_t>(r.Size()) !=
        std::int64_t{g.nd} * static_cast<std::int64_t>(nc) * g.ne) {
      return ctx.registerErrorMessage(
          std::string{name} +
          ": unexpected size of the values element by element");
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      auto a = internals::IntegrationPointsFunctionsKernelArguments{};
      if (!makeKernelArguments(ctx, a, qs, &v, nc, g, use_device, name)) {
        return false;
      }
      auto am = internals::IntegrationRuleKernelArguments{};
      am.shape_functions = m.maps->B.Read(use_device);
      am.determinants = m.jacobians->detJ.Read(use_device);
      am.weights = m.weights.Read(use_device);
      am.nq = m.nq;
      internals::launchLocalProjectionsKernel(a, am, r.ReadWrite(use_device),
                                              use_device);
    } catch (...) {
      return ctx.registerErrorMessage(getExceptionMessage(name));
    }
    return true;
  }  // end of computeLocalProjections

  bool computeMassMatrices(Context& ctx,
                           mfem::Vector& m,
                           const ElementsGeometry& g,
                           const mfem::ElementDofOrdering o,
                           const bool add) noexcept {
    const auto* const name = "computeMassMatrices";
    if ((g.rule == nullptr) || (g.maps == nullptr) ||
        (g.jacobians == nullptr)) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": invalid geometric data");
    }
    if (static_cast<std::int64_t>(m.Size()) !=
        std::int64_t{g.nd} * g.nd * g.ne) {
      return ctx.registerErrorMessage(std::string{name} +
                                      ": unexpected size of the mass matrices");
    }
    if (g.ne == 0) {
      return true;
    }
    try {
      const auto use_device = mfem::Device::IsEnabled();
      // native index of each degree of freedom in the lexicographic ordering,
      // copied from the finite element: the arrays of the finite elements,
      // used by MFEM on the host, are not read on the device (the debug
      // device of MFEM would protect their memory on the host)
      auto permutation = mfem::Array<int>{};
      const int* p = nullptr;
      if (o == mfem::ElementDofOrdering::LEXICOGRAPHIC) {
        const auto* const tfe =
            dynamic_cast<const mfem::TensorBasisElement*>(g.maps->FE);
        if (tfe == nullptr) {
          return ctx.registerErrorMessage(
              std::string{name} +
              ": the lexicographic ordering requires a tensor basis");
        }
        const auto& dof_map = tfe->GetDofMap();
        if (dof_map.Size() != 0) {
          if (dof_map.Size() != g.nd) {
            return ctx.registerErrorMessage(std::string{name} +
                                            ": invalid lexicographic ordering");
          }
          permutation.SetSize(dof_map.Size());
          for (int i = 0; i != dof_map.Size(); ++i) {
            permutation[i] = dof_map[i];
          }
          p = permutation.Read(use_device);
        }
      }
      auto a = internals::IntegrationRuleKernelArguments{};
      a.shape_functions = g.maps->B.Read(use_device);
      a.determinants = g.jacobians->detJ.Read(use_device);
      a.weights = g.weights.Read(use_device);
      a.nq = g.nq;
      auto* const r = add ? m.ReadWrite(use_device) : m.Write(use_device);
      internals::launchMassMatricesKernel(a, p, g.nd, g.ne, r, add, use_device);
    } catch (...) {
      return ctx.registerErrorMessage(getExceptionMessage(name));
    }
    return true;
  }  // end of computeMassMatrices

}  // end of namespace mfem_mgis::batch
