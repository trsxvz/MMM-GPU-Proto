/*!
 * \file   src/NodalBlockJacobiKernels.cxx
 * \brief  This file implements the kernels building the nodal blocks of the
 * Jacobi preconditioner by nodal blocks. It contains kernels of MFEM: it is
 * compiled as CUDA or HIP code when MFEM is compiled with CUDA or HIP (see
 * `src/CMakeLists.txt`).
 * \date   15/10/2026
 */

#include "mfem/general/backends.hpp"
#include "mfem/general/forall.hpp"
#include "NodalBlockJacobiKernels.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return the unknown of a component of a node
   * \param[in] l: layout of the unknowns
   * \param[in] node: node
   * \param[in] c: component
   */
  MFEM_HOST_DEVICE inline int getUnknown(const NodalBlocksLayout& l,
                                         const int node,
                                         const int c) {
    return l.by_nodes ? c * l.nb + node : node * l.bs + c;
  }  // end of getUnknown

  /*!
   * \return the node of an unknown
   * \param[in] l: layout of the unknowns
   * \param[in] i: unknown
   */
  MFEM_HOST_DEVICE inline int getNode(const NodalBlocksLayout& l, const int i) {
    return l.by_nodes ? i % l.nb : i / l.bs;
  }  // end of getNode

  /*!
   * \return the component of an unknown
   * \param[in] l: layout of the unknowns
   * \param[in] i: unknown
   */
  MFEM_HOST_DEVICE inline int getComponent(const NodalBlocksLayout& l,
                                           const int i) {
    return l.by_nodes ? i / l.nb : i % l.bs;
  }  // end of getComponent

  void launchNodalBlocksExtractionKernel(const NodalBlocksLayout& l,
                                         const int* const I,
                                         const int* const J,
                                         const real* const A,
                                         real* const blocks,
                                         const bool use_device) {
    const auto n = l.nb * l.bs;
    mfem::forall_switch(use_device, n,
                        [l, I, J, A, blocks] MFEM_HOST_DEVICE(int i) {
                          const auto bs = l.bs;
                          const auto node = getNode(l, i);
                          const auto c = getComponent(l, i);
                          // the row c of the block of the node, written by this
                          // thread only
                          auto* const b = blocks + bs * bs * node;
                          for (int d = 0; d != bs; ++d) {
                            b[c + bs * d] = 0;
                          }
                          for (int k = I[i]; k != I[i + 1]; ++k) {
                            const auto j = J[k];
                            if (getNode(l, j) == node) {
                              b[c + bs * getComponent(l, j)] += A[k];
                            }
                          }
                        });
  }  // end of launchNodalBlocksExtractionKernel

  void launchNodalBlocksGatherKernel(const NodalBlocksLayout& l,
                                     const real* const columns,
                                     real* const blocks,
                                     const bool use_device) {
    const auto bs = l.bs;
    const auto n = l.nb * bs * bs;
    mfem::forall_switch(
        use_device, n, [l, columns, blocks] MFEM_HOST_DEVICE(int t) {
          const auto bs = l.bs;
          const auto c = t % bs;
          const auto d = (t / bs) % bs;
          const auto node = t / (bs * bs);
          blocks[t] = columns[d * l.nb * bs + getUnknown(l, node, c)];
        });
  }  // end of launchNodalBlocksGatherKernel

}  // end of namespace mfem_mgis::internals
