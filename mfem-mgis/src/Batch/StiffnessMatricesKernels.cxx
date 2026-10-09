/*!
 * \file   src/Batch/StiffnessMatricesKernels.cxx
 * \brief  This file implements the kernels computing the stiffness matrices
 * of the elements of the standard mechanical behaviour integrators. It
 * contains kernels of MFEM: it is compiled as CUDA code when MFEM is
 * compiled with CUDA.
 * \date   09/10/2026
 */

#include <cstddef>
#include "mfem/general/forall.hpp"
#include "KernelsUtilities.hxx"
#include "StiffnessMatricesKernels.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief launches the kernel adding the stiffness matrices of the
   * elements of a material: one thread per pair of nodes of an element,
   * computing the `N x N` block `Σ B_aᵀ K B_b w det J`
   * \tparam K: kind of the gradients
   * \tparam N: dimension
   */
  template <GradientsKind K, int N>
  static void launchStiffnessMatricesKernel(
      const StiffnessMatricesKernelArguments& a, const bool use_device) {
    const auto n = a.number_of_elements * a.nd * a.nd;
    mfem::forall_switch(use_device, n, [a] MFEM_HOST_DEVICE(int i) {
      const auto nd = static_cast<std::size_t>(a.nd);
      const auto nq = static_cast<std::size_t>(a.nq);
      const auto gs = static_cast<std::size_t>(a.gradients_size);
      const auto k = static_cast<std::size_t>(i) / (nd * nd);
      const auto r = static_cast<std::size_t>(i) - k * nd * nd;
      const auto na = r % nd;
      const auto nb = r / nd;
      const auto e = static_cast<std::size_t>(a.elements[k]);
      real m[N][N] = {};
      for (std::size_t q = 0; q != nq; ++q) {
        const auto detJ = a.determinants[q + nq * e];
        real iJ[N * N];
        invertJacobian<N>(iJ, a.jacobians + q + nq * N * N * e, nq, detJ);
        real dNa[N];
        real dNb[N];
        computePhysicalDerivatives<N>(dNa, a.shape_functions_gradients, q, nq,
                                      na, iJ);
        computePhysicalDerivatives<N>(dNb, a.shape_functions_gradients, q, nq,
                                      nb, iJ);
        real Ba[9][3];
        real Bb[9][3];
        computeB<K>(Ba, dNa);
        computeB<K>(Bb, dNb);
        const auto o =
            static_cast<long long>(a.offsets[k]) + static_cast<long long>(q);
        const auto* const t = a.tangent_operator + o * a.point_stride;
        const auto cs = a.component_stride;
        const auto w = a.weights[q] * detJ;
        // K B_b, then B_aᵀ (K B_b)
        for (std::size_t d = 0; d != N; ++d) {
          real KB[9];
          for (std::size_t r1 = 0; r1 != gs; ++r1) {
            KB[r1] = 0;
            for (std::size_t r2 = 0; r2 != gs; ++r2) {
              KB[r1] +=
                  t[static_cast<long long>(r1 * gs + r2) * cs] * Bb[r2][d];
            }
          }
          for (std::size_t c = 0; c != N; ++c) {
            auto v = real{0};
            for (std::size_t r1 = 0; r1 != gs; ++r1) {
              v += Ba[r1][c] * KB[r1];
            }
            m[c][d] += w * v;
          }
        }
      }
      // column-major element matrix, degrees of freedom `node + nd * c`
      const auto ndofs = nd * N;
      auto* const me = a.matrices + ndofs * ndofs * e;
      for (std::size_t c = 0; c != N; ++c) {
        for (std::size_t d = 0; d != N; ++d) {
          me[(na + nd * c) + ndofs * (nb + nd * d)] += m[c][d];
        }
      }
    });
  }  // end of launchStiffnessMatricesKernel

  void launchStiffnessMatricesKernel(const GradientsKind k,
                                     const StiffnessMatricesKernelArguments& a,
                                     const bool use_device) {
    switch (k) {
      case GradientsKind::SMALL_STRAIN_3D:
        launchStiffnessMatricesKernel<GradientsKind::SMALL_STRAIN_3D, 3>(
            a, use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchStiffnessMatricesKernel<GradientsKind::SMALL_STRAIN_2D, 2>(
            a, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchStiffnessMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_3D,
                                      3>(a, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchStiffnessMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_2D,
                                      2>(a, use_device);
        break;
    }
  }  // end of launchStiffnessMatricesKernel

}  // end of namespace mfem_mgis::batch::internals
