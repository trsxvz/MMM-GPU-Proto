/*!
 * \file   src/Batch/InnerForcesKernels.cxx
 * \brief  This file implements the kernels computing the inner forces of
 * the standard mechanical behaviour integrators. It contains kernels of
 * MFEM: it is compiled as CUDA code when MFEM is compiled with CUDA.
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
   * \param[in, out] f: inner forces of the node
   * \param[in] dN: derivatives of the shape function of the node
   * \param[in] s: thermodynamic force (accessor)
   * \param[in] w: weight of the point
   */
  template <GradientsKind K, typename Accessor>
  MFEM_HOST_DEVICE inline void addInnerForces(real* const f,
                                              const real* const dN,
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
    } else {
      static_assert(K == GradientsKind::DEFORMATION_GRADIENT_2D);
      f[0] += w * (dN[0] * s(0) + dN[1] * s(3));
      f[1] += w * (dN[0] * s(4) + dN[1] * s(1));
    }
  }  // end of addInnerForces

  /*!
   * \brief launches the kernel adding the inner forces of the elements of a
   * material: one thread per element, which inverts the jacobian once per
   * point (a thread per node of an element inverted it for each node)
   * \tparam K: kind of the gradients
   * \tparam N: dimension
   */
  template <GradientsKind K, int N>
  static void launchInnerForcesKernel(const InnerForcesKernelArguments& a,
                                      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a] MFEM_HOST_DEVICE(int k) {
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto nq = static_cast<std::size_t>(a.nq);
          const auto nd = static_cast<std::size_t>(a.nd);
          auto* const fe = a.values + nd * N * e;
          for (std::size_t q = 0; q != nq; ++q) {
            const auto detJ = a.determinants[q + nq * e];
            real iJ[N * N];
            invertJacobian<N>(iJ, a.jacobians + q + nq * N * N * e, nq, detJ);
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
              real f[N] = {};
              addInnerForces<K>(f, dN, s, w);
              for (std::size_t c = 0; c != N; ++c) {
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
        launchInnerForcesKernel<GradientsKind::SMALL_STRAIN_3D, 3>(a,
                                                                   use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchInnerForcesKernel<GradientsKind::SMALL_STRAIN_2D, 2>(a,
                                                                   use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchInnerForcesKernel<GradientsKind::DEFORMATION_GRADIENT_3D, 3>(
            a, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchInnerForcesKernel<GradientsKind::DEFORMATION_GRADIENT_2D, 2>(
            a, use_device);
        break;
    }
  }  // end of launchInnerForcesKernel

}  // end of namespace mfem_mgis::batch::internals
