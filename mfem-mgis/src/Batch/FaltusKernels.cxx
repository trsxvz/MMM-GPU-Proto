/*!
 * \file   src/Batch/FaltusKernels.cxx
 * \brief  This file implements the kernels of the regularization of Faltus
 * et al. It contains kernels of MFEM: it is compiled as CUDA code when MFEM
 * is compiled with CUDA.
 * \date   11/10/2026
 */

#include <cstddef>
#include "mfem/general/forall.hpp"
#include "KernelsUtilities.hxx"
#include "FaltusKernels.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \brief inverse of the jacobian at the center of an element
   * \tparam N: dimension
   * \param[out] iJ0: inverse of the jacobian
   * \param[in] a: geometric data
   * \param[in] e: element of the mesh
   */
  template <int N>
  MFEM_HOST_DEVICE inline void invertCenterJacobian(
      real* const iJ0,
      const FaltusGeometryKernelArguments& a,
      const std::size_t e) {
    invertJacobian<N>(iJ0,
                      a.center_jacobians + static_cast<std::size_t>(N * N) * e,
                      1, a.center_determinants[e]);
  }  // end of invertCenterJacobian

  /*!
   * \brief inverse of the jacobian at an integration point of an element
   * \tparam N: dimension
   * \return the weight of the point, `w det J`
   * \param[out] iJ: inverse of the jacobian
   * \param[in] a: geometric data
   * \param[in] e: element of the mesh
   * \param[in] q: integration point
   */
  template <int N>
  MFEM_HOST_DEVICE inline real invertPointJacobian(
      real* const iJ,
      const FaltusGeometryKernelArguments& a,
      const std::size_t e,
      const std::size_t q) {
    const auto nq = static_cast<std::size_t>(a.nq);
    const auto detJ = a.determinants[q + nq * e];
    invertJacobian<N>(iJ, a.jacobians + q + nq * N * N * e, nq, detJ);
    return a.weights[q] * detJ;
  }  // end of invertPointJacobian

  /*!
   * \brief computes the block of a node of the operator `D = B - B0` at an
   * integration point: `D[c][v]` is the derivative of the `c`-th component
   * of `F - F0` with respect to the `v`-th component of the unknowns at the
   * node
   * \tparam K: kind of the gradients
   * \param[out] D: block of the operator
   * \param[in] a: geometric data
   * \param[in] q: integration point
   * \param[in] node: node
   * \param[in] iJ: inverse of the jacobian at the point
   * \param[in] iJ0: inverse of the jacobian at the center
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline void computeFaltusOperator(
      real (&D)[9][3],
      const FaltusGeometryKernelArguments& a,
      const std::size_t q,
      const std::size_t node,
      const real* const iJ,
      const real* const iJ0) {
    constexpr auto N = getDimension<K>();
    real dN[N];
    computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, q,
                                  static_cast<std::size_t>(a.nq), node, iJ);
    real dN0[N];
    computePhysicalDerivatives<N>(dN0, a.center_shape_functions_gradients, 0, 1,
                                  node, iJ0);
    real B0[9][3];
    computeB<K>(D, dN, real{0});
    computeB<K>(B0, dN0, real{0});
    for (int c = 0; c != getGradientsSize<K>(); ++c) {
      for (int v = 0; v != getUnknownsSize<K>(); ++v) {
        D[c][v] -= B0[c][v];
      }
    }
  }  // end of computeFaltusOperator

  /*!
   * \brief computes the difference `d = ∇N(q) - ∇N(0)` of the derivatives
   * of the shape function of a node at an integration point and at the
   * center of the element. For a deformation gradient, the block of the node
   * of the operator `D = B - B0` has the components `D[(v, j)][v] = d[j]`,
   * `(v, j)` being the component `dv/dxj` of the gradient: the block `(a,
   * b)` of `Dᵀ D` is `d_a · d_b` times the identity
   * \tparam N: dimension
   * \param[out] d: difference of the derivatives
   * \param[in] a: geometric data
   * \param[in] q: integration point
   * \param[in] node: node
   * \param[in] iJ: inverse of the jacobian at the point
   * \param[in] iJ0: inverse of the jacobian at the center
   */
  template <int N>
  MFEM_HOST_DEVICE inline void computeFaltusDerivativesDifference(
      real* const d,
      const FaltusGeometryKernelArguments& a,
      const std::size_t q,
      const std::size_t node,
      const real* const iJ,
      const real* const iJ0) {
    real dN0[N];
    computePhysicalDerivatives<N>(d, a.shape_functions_gradients, q,
                                  static_cast<std::size_t>(a.nq), node, iJ);
    computePhysicalDerivatives<N>(dN0, a.center_shape_functions_gradients, 0, 1,
                                  node, iJ0);
    for (int c = 0; c != N; ++c) {
      d[c] -= dN0[c];
    }
  }  // end of computeFaltusDerivativesDifference

  /*!
   * \return the coefficient `α Σ_q w_q d_a · d_b` of the block of the nodes
   * `a` and `b` of the matrix `α Σ_q w_q Dᵀ D` of an element (see
   * `computeFaltusDerivativesDifference`)
   * \tparam K: kind of the gradients (deformation gradient)
   * \param[in] a: geometric data
   * \param[in] e: element of the mesh
   * \param[in] na: first node
   * \param[in] nb: second node
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline real computeFaltusNodesCoefficient(
      const FaltusGeometryKernelArguments& a,
      const std::size_t e,
      const std::size_t na,
      const std::size_t nb) {
    static_assert((K == GradientsKind::DEFORMATION_GRADIENT_3D) ||
                  (K == GradientsKind::DEFORMATION_GRADIENT_2D));
    constexpr auto N = getDimension<K>();
    real iJ0[N * N];
    invertCenterJacobian<N>(iJ0, a, e);
    auto v = real{0};
    for (std::size_t q = 0; q != static_cast<std::size_t>(a.nq); ++q) {
      real iJ[N * N];
      const auto w = invertPointJacobian<N>(iJ, a, e, q);
      real da[N];
      real db[N];
      computeFaltusDerivativesDifference<N>(da, a, q, na, iJ, iJ0);
      computeFaltusDerivativesDifference<N>(db, a, q, nb, iJ, iJ0);
      auto s = real{0};
      for (int c = 0; c != N; ++c) {
        s += da[c] * db[c];
      }
      v += w * s;
    }
    return a.alpha * v;
  }  // end of computeFaltusNodesCoefficient

  template <GradientsKind K>
  static void launchFaltusCentersKernel(const FaltusCentersKernelArguments& a,
                                        const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a] MFEM_HOST_DEVICE(int k) {
          constexpr auto N = getDimension<K>();
          constexpr auto V = getUnknownsSize<K>();
          constexpr auto gs = getGradientsSize<K>();
          const auto e = static_cast<std::size_t>(a.elements[k]);
          real iJ[N * N];
          invertJacobian<N>(iJ,
                            a.jacobians + static_cast<std::size_t>(N * N) * e,
                            1, a.determinants[e]);
          computeDeformationGradientAtCenter<K>(
              a.centers +
                  static_cast<std::size_t>(gs) * static_cast<std::size_t>(k),
              a.macroscopic_gradients,
              (a.derivatives != nullptr)
                  ? a.derivatives + static_cast<std::size_t>(V * N) * e
                  : nullptr,
              iJ);
        });
  }  // end of launchFaltusCentersKernel

  template <GradientsKind K>
  static void launchFaltusInnerForcesKernel(
      const FaltusInnerForcesKernelArguments& a, const bool use_device) {
    mfem::forall_switch(
        use_device, a.geometry.number_of_elements, [a] MFEM_HOST_DEVICE(int k) {
          constexpr auto N = getDimension<K>();
          constexpr auto V = getUnknownsSize<K>();
          constexpr auto gs = getGradientsSize<K>();
          const auto& g = a.geometry;
          const auto e = static_cast<std::size_t>(g.elements[k]);
          const auto nd = static_cast<std::size_t>(g.nd);
          constexpr auto sgs = static_cast<std::size_t>(gs);
          const auto np = static_cast<std::size_t>(a.number_of_points);
          const auto* const F0 = a.centers + sgs * static_cast<std::size_t>(k);
          auto* const fe = a.values + nd * V * e;
          real iJ0[N * N];
          invertCenterJacobian<N>(iJ0, g, e);
          for (std::size_t q = 0; q != static_cast<std::size_t>(g.nq); ++q) {
            real iJ[N * N];
            const auto w = invertPointJacobian<N>(iJ, g, e, q);
            const auto o = static_cast<std::size_t>(a.offsets[k]) + q;
            // α (F - F0)
            real dF[9];
            for (int c = 0; c != gs; ++c) {
              const auto sc = static_cast<std::size_t>(c);
              const auto F = a.point_by_point ? a.gradients[o * sgs + sc]
                                              : a.gradients[sc * np + o];
              dF[c] = g.alpha * (F - F0[c]);
            }
            for (std::size_t node = 0; node != nd; ++node) {
              real D[9][3];
              computeFaltusOperator<K>(D, g, q, node, iJ, iJ0);
              for (int v = 0; v != V; ++v) {
                auto f = real{0};
                for (int c = 0; c != gs; ++c) {
                  f += D[c][v] * dF[c];
                }
                fe[node + nd * static_cast<std::size_t>(v)] += w * f;
              }
            }
          }
        });
  }  // end of launchFaltusInnerForcesKernel

  template <GradientsKind K>
  static void launchFaltusMatricesKernel(const FaltusGeometryKernelArguments& a,
                                         real* const m,
                                         const bool use_device) {
    // one thread per pair of nodes: the block of the pair is the coefficient
    // of the nodes times the identity
    const auto nd = a.nd;
    mfem::forall_switch(
        use_device, a.number_of_elements * nd * nd,
        [a, m, nd] MFEM_HOST_DEVICE(int t) {
          constexpr auto V = getUnknownsSize<K>();
          const auto k = t / (nd * nd);
          const auto r = t - k * nd * nd;
          const auto nb = static_cast<std::size_t>(r / nd);
          const auto na = static_cast<std::size_t>(r) - nb * nd;
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto snd = static_cast<std::size_t>(nd);
          const auto n = snd * V;
          const auto s = computeFaltusNodesCoefficient<K>(a, e, na, nb);
          // row `na + nd * v`, column `nb + nd * v`, stored column by column
          for (std::size_t v = 0; v != V; ++v) {
            m[n * n * e + (na + snd * v) + n * (nb + snd * v)] += s;
          }
        });
  }  // end of launchFaltusMatricesKernel

  template <GradientsKind K>
  static void launchFaltusMatricesKernel(const FaltusGeometryKernelArguments& a,
                                         const SparseMatrixOutput& o,
                                         const bool use_device) {
    // one thread per node of an element (rows of the node)
    const auto nd = a.nd;
    mfem::forall_switch(
        use_device, a.number_of_elements * nd,
        [a, o, nd] MFEM_HOST_DEVICE(int t) {
          constexpr auto V = getUnknownsSize<K>();
          const auto k = t / nd;
          const auto na = static_cast<std::size_t>(t - k * nd);
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto snd = static_cast<std::size_t>(nd);
          const auto* const dofs = o.element_dofs + snd * V * e;
          for (std::size_t nb = 0; nb != snd; ++nb) {
            const auto s = computeFaltusNodesCoefficient<K>(a, e, na, nb);
            for (std::size_t v = 0; v != V; ++v) {
              const auto row = dofs[na + snd * v];
              if ((row < 0) || (row >= o.number_of_rows)) {
                ::AtomicAdd(o.missing[0], 1);
                continue;
              }
              const auto begin = o.I[row];
              const auto end = o.I[row + 1];
              const auto p =
                  o.findPosition(e, snd * V, na + snd * v, nb + snd * v, begin,
                                 end, dofs[nb + snd * v]);
              if (p == end) {
                ::AtomicAdd(o.missing[0], 1);
              } else {
                ::AtomicAdd(o.values[p], s);
              }
            }
          }
        });
  }  // end of launchFaltusMatricesKernel

  template <GradientsKind K>
  static void launchFaltusDiagonalsKernel(
      const FaltusGeometryKernelArguments& a,
      real* const d,
      const bool use_device) {
    // one thread per node of an element
    const auto nd = a.nd;
    mfem::forall_switch(
        use_device, a.number_of_elements * nd,
        [a, d, nd] MFEM_HOST_DEVICE(int t) {
          constexpr auto V = getUnknownsSize<K>();
          const auto k = t / nd;
          const auto na = static_cast<std::size_t>(t - k * nd);
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto snd = static_cast<std::size_t>(nd);
          const auto s = computeFaltusNodesCoefficient<K>(a, e, na, na);
          for (std::size_t v = 0; v != V; ++v) {
            d[snd * V * e + na + snd * v] += s;
          }
        });
  }  // end of launchFaltusDiagonalsKernel

  template <GradientsKind K>
  static void launchFaltusActionKernel(const FaltusGeometryKernelArguments& a,
                                       const real* const x,
                                       real* const y,
                                       const bool use_device) {
    // one thread per element: `y_a = α Σ_q w_q d_a (Σ_b d_b · x_b)` for each
    // component of the unknowns
    mfem::forall_switch(
        use_device, a.number_of_elements, [a, x, y] MFEM_HOST_DEVICE(int k) {
          constexpr auto N = getDimension<K>();
          constexpr auto V = getUnknownsSize<K>();
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto nd = static_cast<std::size_t>(a.nd);
          const auto* const xe = x + nd * V * e;
          auto* const ye = y + nd * V * e;
          real iJ0[N * N];
          invertCenterJacobian<N>(iJ0, a, e);
          for (std::size_t q = 0; q != static_cast<std::size_t>(a.nq); ++q) {
            real iJ[N * N];
            const auto w = invertPointJacobian<N>(iJ, a, e, q);
            // g[v][j] = Σ_b d_b[j] x_b[v]
            real g[V][N] = {};
            for (std::size_t node = 0; node != nd; ++node) {
              real d[N];
              computeFaltusDerivativesDifference<N>(d, a, q, node, iJ, iJ0);
              for (std::size_t v = 0; v != V; ++v) {
                const auto xv = xe[node + nd * v];
                for (int j = 0; j != N; ++j) {
                  g[v][j] += d[j] * xv;
                }
              }
            }
            for (std::size_t node = 0; node != nd; ++node) {
              real d[N];
              computeFaltusDerivativesDifference<N>(d, a, q, node, iJ, iJ0);
              for (std::size_t v = 0; v != V; ++v) {
                auto s = real{0};
                for (int j = 0; j != N; ++j) {
                  s += d[j] * g[v][j];
                }
                ye[node + nd * v] += w * a.alpha * s;
              }
            }
          }
        });
  }  // end of launchFaltusActionKernel

  void launchFaltusCentersKernel(const GradientsKind k,
                                 const FaltusCentersKernelArguments& a,
                                 const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFaltusCentersKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, use_device);
    } else {
      launchFaltusCentersKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, use_device);
    }
  }  // end of launchFaltusCentersKernel

  void launchFaltusInnerForcesKernel(const GradientsKind k,
                                     const FaltusInnerForcesKernelArguments& a,
                                     const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFaltusInnerForcesKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, use_device);
    } else {
      launchFaltusInnerForcesKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, use_device);
    }
  }  // end of launchFaltusInnerForcesKernel

  void launchFaltusMatricesKernel(const GradientsKind k,
                                  const FaltusGeometryKernelArguments& a,
                                  real* const m,
                                  const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFaltusMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, m, use_device);
    } else {
      launchFaltusMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, m, use_device);
    }
  }  // end of launchFaltusMatricesKernel

  void launchFaltusMatricesKernel(const GradientsKind k,
                                  const FaltusGeometryKernelArguments& a,
                                  const SparseMatrixOutput& o,
                                  const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFaltusMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, o, use_device);
    } else {
      launchFaltusMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, o, use_device);
    }
  }  // end of launchFaltusMatricesKernel

  void launchFaltusDiagonalsKernel(const GradientsKind k,
                                   const FaltusGeometryKernelArguments& a,
                                   real* const d,
                                   const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFaltusDiagonalsKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, d, use_device);
    } else {
      launchFaltusDiagonalsKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, d, use_device);
    }
  }  // end of launchFaltusDiagonalsKernel

  void launchFaltusActionKernel(const GradientsKind k,
                                const FaltusGeometryKernelArguments& a,
                                const real* const x,
                                real* const y,
                                const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFaltusActionKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, x, y, use_device);
    } else {
      launchFaltusActionKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, x, y, use_device);
    }
  }  // end of launchFaltusActionKernel

}  // end of namespace mfem_mgis::batch::internals
