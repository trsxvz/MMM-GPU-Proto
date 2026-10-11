/*!
 * \file   src/Batch/ElementsGeometryKernels.cxx
 * \brief  This file implements the kernel computing the determinants and the
 * inverses of the jacobians of the element transformations. It contains kernels
 * of MFEM: it is compiled as CUDA code when MFEM is compiled with CUDA.
 * \date   10/10/2026
 */

#include <cstddef>
#include "mfem/general/forall.hpp"
#include "KernelsUtilities.hxx"
#include "ElementsGeometryKernels.hxx"

namespace mfem_mgis::batch::internals {

  void launchJacobiansKernel(const JacobiansKernelArguments& a,
                             const bool use_device) {
    const auto n = a.number_of_elements * a.nq;
    mfem::forall_switch(use_device, n, [a] MFEM_HOST_DEVICE(int i) {
      const auto nq = static_cast<std::size_t>(a.nq);
      const auto nd = static_cast<std::size_t>(a.nodes_per_element);
      const auto dim = static_cast<std::size_t>(a.dim);
      const auto ns = static_cast<std::size_t>(a.node_stride);
      const auto cs = static_cast<std::size_t>(a.component_stride);
      const auto e = static_cast<std::size_t>(i) / nq;
      const auto q = static_cast<std::size_t>(i) - e * nq;
      const auto* const nodes = a.element_nodes + nd * e;
      // `J[r + 3 * c]`: derivative of the `r`-th coordinate with respect to
      // the `c`-th reference coordinate
      real J[9] = {};
      for (std::size_t c = 0; c != dim; ++c) {
        const auto* const G = a.shape_functions_gradients + q + nq * c;
        for (std::size_t node = 0; node != nd; ++node) {
          const auto g = G[nq * dim * node];
          const auto* const x =
              a.coordinates + static_cast<std::size_t>(nodes[node]) * ns;
          for (std::size_t r = 0; r != dim; ++r) {
            J[r + 3 * c] += g * x[cs * r];
          }
        }
      }
      auto det = J[0];
      if (dim == 2) {
        det = J[0] * J[4] - J[1] * J[3];
      } else if (dim == 3) {
        det = J[0] * (J[4] * J[8] - J[5] * J[7]) -
              J[3] * (J[1] * J[8] - J[2] * J[7]) +
              J[6] * (J[1] * J[5] - J[2] * J[4]);
      }
      a.determinants[q + nq * e] = det;
      // inverse, `iJ[r + dim * c]` being its component `(r, c)`
      real Jd[9] = {};
      for (std::size_t c = 0; c != dim; ++c) {
        for (std::size_t r = 0; r != dim; ++r) {
          Jd[r + dim * c] = J[r + 3 * c];
        }
      }
      real iJ[9] = {};
      if (dim == 2) {
        invertJacobian<2>(iJ, Jd, 1, det);
      } else if (dim == 3) {
        invertJacobian<3>(iJ, Jd, 1, det);
      } else {
        iJ[0] = 1 / det;
      }
      auto* const iJe = a.inverse_jacobians + q + nq * dim * dim * e;
      for (std::size_t k = 0; k != dim * dim; ++k) {
        iJe[nq * k] = iJ[k];
      }
    });
  }  // end of launchJacobiansKernel

}  // end of namespace mfem_mgis::batch::internals
