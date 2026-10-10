/*!
 * \file   src/Batch/ElementsGeometryKernels.cxx
 * \brief  This file implements the kernel computing the jacobians of the
 * element transformations of a mesh without nodes. It contains kernels of
 * MFEM: it is compiled as CUDA code when MFEM is compiled with CUDA.
 * \date   10/10/2026
 */

#include <cstddef>
#include "mfem/general/forall.hpp"
#include "ElementsGeometryKernels.hxx"

namespace mfem_mgis::batch::internals {

  void launchJacobiansKernel(const JacobiansKernelArguments& a,
                             const bool use_device) {
    const auto n = a.number_of_elements * a.nq;
    mfem::forall_switch(use_device, n, [a] MFEM_HOST_DEVICE(int i) {
      const auto nq = static_cast<std::size_t>(a.nq);
      const auto nd = static_cast<std::size_t>(a.vertices_per_element);
      const auto nv = static_cast<std::size_t>(a.number_of_vertices);
      const auto dim = static_cast<std::size_t>(a.dim);
      const auto e = static_cast<std::size_t>(i) / nq;
      const auto q = static_cast<std::size_t>(i) - e * nq;
      const auto* const v = a.element_vertices + nd * e;
      // `J[r + 3 * c]`: derivative of the `r`-th coordinate with respect to
      // the `c`-th reference coordinate
      real J[9] = {};
      for (std::size_t c = 0; c != dim; ++c) {
        const auto* const G = a.shape_functions_gradients + q + nq * c;
        for (std::size_t node = 0; node != nd; ++node) {
          const auto g = G[nq * dim * node];
          const auto* const x = a.vertices + static_cast<std::size_t>(v[node]);
          for (std::size_t r = 0; r != dim; ++r) {
            J[r + 3 * c] += g * x[nv * r];
          }
        }
      }
      auto* const Je = a.jacobians + q + nq * dim * dim * e;
      for (std::size_t c = 0; c != dim; ++c) {
        for (std::size_t r = 0; r != dim; ++r) {
          Je[nq * (r + dim * c)] = J[r + 3 * c];
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
    });
  }  // end of launchJacobiansKernel

}  // end of namespace mfem_mgis::batch::internals
