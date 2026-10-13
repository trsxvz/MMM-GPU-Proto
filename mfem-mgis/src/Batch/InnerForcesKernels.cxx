/*!
 * \file   src/Batch/InnerForcesKernels.cxx
 * \brief  This file implements the kernels computing the inner forces of the
 * standard behaviour integrators (mechanics, heat transfer, micromorphic
 * damage). It contains kernels of MFEM: it is compiled as CUDA or HIP code when
 * MFEM is compiled with CUDA or HIP.
 * \date   09/10/2026
 */

#include <cstddef>
#include "mfem/general/forall.hpp"
#include "KernelsUtilities.hxx"
#include "InnerForcesKernels.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief adds the contribution of a point to the inner forces of a node,
   * `w Bᵀ σ`, `B` being the operator of the gradients (see
   * `computeGradientsAtPoint` in `GradientsKernels.cxx`)
   * \tparam K: kind of the gradients
   * \tparam Accessor: type of the accessor to the components of the
   * thermodynamic forces
   * \param[in, out] f: inner forces of the node
   * \param[in] dN: derivatives of the shape function of the node
   * \param[in] Nv: value of the shape function of the node (transient term
   * of the heat equation, micromorphic damage)
   * \param[in] s: thermodynamic force (accessor)
   * \param[in] w: weight of the point multiplied by the determinant of the
   * jacobian
   */
  template <GradientsKind K, typename Accessor>
  MFEM_HOST_DEVICE inline void addInnerForces(real* const f,
                                              const real* const dN,
                                              const real Nv,
                                              const Accessor& s,
                                              const real w) {
    constexpr auto icste = getInverseOfSquareRootOfTwo();
    if constexpr (K == GradientsKind::SMALL_STRAIN_3D) {
      f[0] += w * (dN[0] * s(0) + (dN[1] * s(3) + dN[2] * s(4)) * icste);
      f[1] += w * (dN[1] * s(1) + (dN[0] * s(3) + dN[2] * s(5)) * icste);
      f[2] += w * (dN[2] * s(2) + (dN[0] * s(4) + dN[1] * s(5)) * icste);
    } else if constexpr (K == GradientsKind::SMALL_STRAIN_2D) {
      f[0] += w * (dN[0] * s(0) + dN[1] * s(3) * icste);
      f[1] += w * (dN[1] * s(1) + dN[0] * s(3) * icste);
    } else if constexpr (K == GradientsKind::DEFORMATION_GRADIENT_3D) {
      f[0] += w * (dN[0] * s(0) + dN[1] * s(3) + dN[2] * s(5));
      f[1] += w * (dN[0] * s(4) + dN[1] * s(1) + dN[2] * s(7));
      f[2] += w * (dN[0] * s(6) + dN[1] * s(8) + dN[2] * s(2));
    } else if constexpr (K == GradientsKind::DEFORMATION_GRADIENT_2D) {
      f[0] += w * (dN[0] * s(0) + dN[1] * s(3));
      f[1] += w * (dN[0] * s(4) + dN[1] * s(1));
    } else if constexpr (isTemperatureGradient<K>()) {
      // heat flux
      for (int c = 0; c != getDimension<K>(); ++c) {
        f[0] += w * dN[c] * s(c);
      }
    } else if constexpr (isMicromorphicDamage<K>()) {
      // dual forces of the micromorphic damage and of its gradient
      auto v = s(0) * Nv;
      for (int c = 0; c != getDimension<K>(); ++c) {
        v += s(1 + c) * dN[c];
      }
      f[0] += w * v;
    } else {
      // opposite of the rate of the enthalpy
      static_assert(isTemperature<K>());
      f[0] += w * Nv * s(0);
    }
  }  // end of addInnerForces

  /*!
   * \brief launches the kernel adding the inner forces of the elements of a
   * material: one thread per element, which loads the inverse of the
   * jacobian once per point (see `ElementsJacobians`)
   * \tparam K: kind of the gradients
   */
  template <GradientsKind K>
  static void launchInnerForcesKernel(const InnerForcesKernelArguments& a,
                                      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a] MFEM_HOST_DEVICE(int k) {
          // dimension and number of components of the unknowns
          constexpr auto N = getDimension<K>();
          constexpr auto V = getUnknownsSize<K>();
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto nq = static_cast<std::size_t>(a.nq);
          const auto nd = static_cast<std::size_t>(a.nd);
          const auto ke =
              a.by_material_element ? static_cast<std::size_t>(k) : e;
          auto* const fe = a.values + nd * V * ke;
          for (std::size_t q = 0; q != nq; ++q) {
            const auto detJ = a.determinants[q + nq * e];
            real iJ[N * N];
            loadInverseJacobian<N>(iJ, a.inverse_jacobians + q + nq * N * N * e,
                                   nq);
            const auto o = static_cast<long long>(a.offsets[k]) +
                           static_cast<long long>(q);
            const auto* const s0 = a.thermodynamic_forces + o * a.point_stride;
            const auto cs = a.component_stride;
            auto s = [s0, cs](const int c) { return s0[c * cs]; };
            const auto w = a.weights[q] * detJ;
            for (std::size_t node = 0; node != nd; ++node) {
              real dN[N];
              computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, q,
                                            nq, node, iJ);
              const auto Nv =
                  getShapeFunctionValue<K>(a.shape_functions, q, nq, node);
              real f[V] = {};
              addInnerForces<K>(f, dN, Nv, s, w);
              for (std::size_t c = 0; c != V; ++c) {
                fe[node + nd * c] += f[c];
              }
            }
          }
        });
  }  // end of launchInnerForcesKernel

  void launchInnerForcesKernel(const GradientsKind k,
                               const InnerForcesKernelArguments& a,
                               const bool use_device) {
    switch (k) {
      case GradientsKind::SMALL_STRAIN_3D:
        launchInnerForcesKernel<GradientsKind::SMALL_STRAIN_3D>(a, use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchInnerForcesKernel<GradientsKind::SMALL_STRAIN_2D>(a, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchInnerForcesKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
            a, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchInnerForcesKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
            a, use_device);
        break;
      case GradientsKind::TEMPERATURE_GRADIENT_3D:
        launchInnerForcesKernel<GradientsKind::TEMPERATURE_GRADIENT_3D>(
            a, use_device);
        break;
      case GradientsKind::TEMPERATURE_GRADIENT_2D:
        launchInnerForcesKernel<GradientsKind::TEMPERATURE_GRADIENT_2D>(
            a, use_device);
        break;
      case GradientsKind::TEMPERATURE_3D:
        launchInnerForcesKernel<GradientsKind::TEMPERATURE_3D>(a, use_device);
        break;
      case GradientsKind::TEMPERATURE_2D:
        launchInnerForcesKernel<GradientsKind::TEMPERATURE_2D>(a, use_device);
        break;
      case GradientsKind::MICROMORPHIC_DAMAGE_3D:
        launchInnerForcesKernel<GradientsKind::MICROMORPHIC_DAMAGE_3D>(
            a, use_device);
        break;
      case GradientsKind::MICROMORPHIC_DAMAGE_2D:
        launchInnerForcesKernel<GradientsKind::MICROMORPHIC_DAMAGE_2D>(
            a, use_device);
        break;
    }
  }  // end of launchInnerForcesKernel

}  // end of namespace mfem_mgis::batch::internals
