/*!
 * \file   src/Batch/StiffnessMatricesKernels.cxx
 * \brief  This file implements the kernels computing the stiffness matrices
 * of the elements of the standard behaviour integrators (mechanics, heat
 * transfer, micromorphic damage). It contains kernels of MFEM: it is
 * compiled as CUDA code when MFEM is compiled with CUDA.
 * \date   09/10/2026
 */

#include <vector>
#include <cstddef>
#include "mfem/general/backends.hpp"
#include "mfem/general/forall.hpp"
#include "KernelsUtilities.hxx"
#include "StiffnessMatricesKernels.hxx"

namespace mfem_mgis::batch::internals {

  //! \brief stores the blocks in the matrices of the elements
  struct ElementMatricesStorage {
    //! \brief matrices of the elements
    real* matrices;
    /*!
     * \brief adds a block of the matrix of an element
     * \tparam N: number of components of the unknowns
     * \param[in] e: element
     * \param[in] nd: number of nodes per element
     * \param[in] na: first node
     * \param[in] nb: second node
     * \param[in] m: block
     */
    template <int N>
    MFEM_HOST_DEVICE void add(const std::size_t e,
                              const std::size_t nd,
                              const std::size_t na,
                              const std::size_t nb,
                              const real (&m)[N][N]) const {
      // column-major element matrix, degrees of freedom `node + nd * c`
      const auto ndofs = nd * N;
      auto* const me = this->matrices + ndofs * ndofs * e;
      for (std::size_t c = 0; c != N; ++c) {
        for (std::size_t d = 0; d != N; ++d) {
          me[(na + nd * c) + ndofs * (nb + nd * d)] += m[c][d];
        }
      }
    }
  };

  //! \brief adds the blocks to an assembled sparse matrix
  struct SparseMatrixAssembly {
    //! \brief assembled matrix
    SparseMatrixOutput o;
    /*!
     * \brief if true, the values are added atomically (kernel). Otherwise,
     * the blocks are added by a single thread (see
     * `computeStiffnessMatricesOnHost`).
     */
    bool atomic = true;
    //! \brief see `ElementMatricesStorage::add`
    template <int N>
    MFEM_HOST_DEVICE void add(const std::size_t e,
                              const std::size_t nd,
                              const std::size_t na,
                              const std::size_t nb,
                              const real (&m)[N][N]) const {
      const auto ndofs = nd * N;
      const auto* const dofs = this->o.element_dofs + ndofs * e;
      for (std::size_t c = 0; c != N; ++c) {
        const auto row = dofs[na + nd * c];
        if ((row < 0) || (row >= this->o.number_of_rows)) {
          ::AtomicAdd(this->o.missing[0], N);
          continue;
        }
        const auto end = this->o.I[row + 1];
        for (std::size_t d = 0; d != N; ++d) {
          const auto p =
              findColumn(this->o.J, this->o.I[row], end, dofs[nb + nd * d]);
          if (p == end) {
            ::AtomicAdd(this->o.missing[0], 1);
          } else if (this->atomic) {
            ::AtomicAdd(this->o.values[p], m[c][d]);
          } else {
            this->o.values[p] += m[c][d];
          }
        }
      }
    }
  };

  //! \brief adds the diagonals of the blocks to values element by element
  struct ElementDiagonalsStorage {
    //! \brief values element by element (layout of the element restriction)
    real* values;
    //! \brief see `ElementMatricesStorage::add`
    template <int N>
    MFEM_HOST_DEVICE void add(const std::size_t e,
                              const std::size_t nd,
                              const std::size_t na,
                              const std::size_t,
                              const real (&m)[N][N]) const {
      for (std::size_t c = 0; c != N; ++c) {
        this->values[na + nd * (c + N * e)] += m[c][c];
      }
    }
  };

  /*!
   * \brief adds, at a point, the contribution of the derivatives of the
   * thermodynamic forces with respect to the unknowns (heat flux with respect
   * to the temperature) to a block of the stiffness matrix of an element,
   * `w B_aᵀ (dq/dT) N_b`
   * \tparam K: kind of the gradients
   * \param[in, out] m: block
   * \param[in] Ba: operator `B` of the first node
   * \param[in] t: tangent operator blocks of the point
   * \param[in] cs: stride between two components of the blocks
   * \param[in] Nb: value of the shape function of the second node
   * \param[in] w: weight of the point
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline void addUnknownsBlock(
      real (&m)[getUnknownsSize<K>()][getUnknownsSize<K>()],
      const real (&Ba)[9][3],
      const real* const t,
      const long long cs,
      const real Nb,
      const real w) {
    if constexpr (hasUnknownsBlock<K>()) {
      constexpr auto gs = getGradientsSize<K>();
      auto v = real{0};
      for (int r = 0; r != gs; ++r) {
        v += Ba[r][0] * t[static_cast<long long>(gs * gs + r) * cs];
      }
      m[0][0] += w * v * Nb;
    }
  }  // end of addUnknownsBlock

  /*!
   * \brief version of `launchStiffnessMatricesKernel` used without device:
   * one iteration per element, the physical derivatives of the shape functions,
   * the operators `B` and the products `K B` of the nodes being computed once
   * per point instead of once per pair of nodes. The operations on each
   * block are made in the same order as by the kernel: the results are
   * identical.
   * \tparam K: kind of the gradients
   * \tparam Output: destination of the blocks
   * \tparam diagonal: if true, only the pairs of a node with itself are
   * computed
   */
  template <GradientsKind K, typename Output, bool diagonal>
  static void computeStiffnessMatricesOnHost(
      const StiffnessMatricesKernelArguments& a, const Output& o) {
    constexpr auto N = getDimension<K>();
    constexpr auto V = getUnknownsSize<K>();
    //! \brief operator `B` of a node at a point, and its product by `K`
    struct NodeOperators {
      real B[9][3];
      real KB[9][3];
    };
    const auto nd = static_cast<std::size_t>(a.nd);
    const auto nq = static_cast<std::size_t>(a.nq);
    const auto gs = static_cast<std::size_t>(a.gradients_size);
    const auto cs = a.component_stride;
    auto operators = std::vector<NodeOperators>(nq * nd);
    auto weights = std::vector<real>(nq);
    auto tangent_operators = std::vector<const real*>(nq);
    for (int k = 0; k != a.number_of_elements; ++k) {
      const auto e = static_cast<std::size_t>(a.elements[k]);
      for (std::size_t q = 0; q != nq; ++q) {
        const auto detJ = a.determinants[q + nq * e];
        real iJ[N * N];
        invertJacobian<N>(iJ, a.jacobians + q + nq * N * N * e, nq, detJ);
        const auto o2 =
            static_cast<long long>(a.offsets[k]) + static_cast<long long>(q);
        const auto* const t = a.tangent_operator + o2 * a.point_stride;
        weights[q] = a.weights[q] * detJ;
        tangent_operators[q] = t;
        for (std::size_t node = 0; node != nd; ++node) {
          auto& op = operators[q * nd + node];
          real dN[N];
          computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, q, nq,
                                        node, iJ);
          computeB<K>(op.B, dN,
                      getShapeFunctionValue<K>(a.shape_functions, q, nq, node));
          for (std::size_t d = 0; d != V; ++d) {
            for (std::size_t r1 = 0; r1 != gs; ++r1) {
              op.KB[r1][d] = 0;
              for (std::size_t r2 = 0; r2 != gs; ++r2) {
                op.KB[r1][d] +=
                    getTangentOperatorComponent<K>(t, r1, r2, cs) * op.B[r2][d];
              }
            }
          }
        }
      }
      const auto npairs = diagonal ? nd : nd * nd;
      for (std::size_t r = 0; r != npairs; ++r) {
        const auto na = diagonal ? r : r % nd;
        const auto nb = diagonal ? r : r / nd;
        real m[V][V] = {};
        for (std::size_t q = 0; q != nq; ++q) {
          const auto& Ba = operators[q * nd + na].B;
          const auto& KB = operators[q * nd + nb].KB;
          const auto w = weights[q];
          for (std::size_t d = 0; d != V; ++d) {
            for (std::size_t c = 0; c != V; ++c) {
              auto v = real{0};
              for (std::size_t r1 = 0; r1 != gs; ++r1) {
                v += Ba[r1][c] * KB[r1][d];
              }
              m[c][d] += w * v;
            }
          }
          if constexpr (hasUnknownsBlock<K>()) {
            addUnknownsBlock<K>(m, Ba, tangent_operators[q], cs,
                                a.shape_functions[q + nq * nb], w);
          }
        }
        o.template add<V>(e, nd, na, nb, m);
      }
    }
  }  // end of computeStiffnessMatricesOnHost

  /*!
   * \brief launches the kernel adding the stiffness matrices of the
   * elements of a material: one thread per pair of nodes of an element,
   * computing the `V x V` block `Σ B_aᵀ K B_b w det J` (plus the
   * contribution of the derivatives with respect to the unknowns, see
   * `addUnknownsBlock`), `V` being the number of components of the unknowns.
   * Without device, one iteration per element (see
   * `computeStiffnessMatricesOnHost`); the devices of MFEM computing on the
   * host (`omp`, `debug`) use the kernel.
   * \tparam K: kind of the gradients
   * \tparam Output: destination of the blocks
   * \tparam diagonal: if true, only the pairs of a node with itself are
   * computed
   */
  template <GradientsKind K, typename Output, bool diagonal = false>
  static void launchStiffnessMatricesKernel(
      const StiffnessMatricesKernelArguments& a,
      const Output o,
      const bool use_device) {
    if (!use_device) {
      computeStiffnessMatricesOnHost<K, Output, diagonal>(a, o);
      return;
    }
    const auto npairs = diagonal ? a.nd : a.nd * a.nd;
    const auto n = a.number_of_elements * npairs;
    mfem::forall_switch(use_device, n, [a, o, npairs] MFEM_HOST_DEVICE(int i) {
      // dimension and number of components of the unknowns
      constexpr auto N = getDimension<K>();
      constexpr auto V = getUnknownsSize<K>();
      const auto nd = static_cast<std::size_t>(a.nd);
      const auto nq = static_cast<std::size_t>(a.nq);
      const auto gs = static_cast<std::size_t>(a.gradients_size);
      const auto np = static_cast<std::size_t>(npairs);
      const auto k = static_cast<std::size_t>(i) / np;
      const auto r = static_cast<std::size_t>(i) - k * np;
      const auto na = diagonal ? r : r % nd;
      const auto nb = diagonal ? r : r / nd;
      const auto e = static_cast<std::size_t>(a.elements[k]);
      real m[V][V] = {};
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
        computeB<K>(Ba, dNa,
                    getShapeFunctionValue<K>(a.shape_functions, q, nq, na));
        computeB<K>(Bb, dNb,
                    getShapeFunctionValue<K>(a.shape_functions, q, nq, nb));
        const auto o2 =
            static_cast<long long>(a.offsets[k]) + static_cast<long long>(q);
        const auto* const t = a.tangent_operator + o2 * a.point_stride;
        const auto cs = a.component_stride;
        const auto w = a.weights[q] * detJ;
        // K B_b, then B_aᵀ (K B_b)
        for (std::size_t d = 0; d != V; ++d) {
          real KB[9];
          for (std::size_t r1 = 0; r1 != gs; ++r1) {
            KB[r1] = 0;
            for (std::size_t r2 = 0; r2 != gs; ++r2) {
              KB[r1] +=
                  getTangentOperatorComponent<K>(t, r1, r2, cs) * Bb[r2][d];
            }
          }
          for (std::size_t c = 0; c != V; ++c) {
            auto v = real{0};
            for (std::size_t r1 = 0; r1 != gs; ++r1) {
              v += Ba[r1][c] * KB[r1];
            }
            m[c][d] += w * v;
          }
        }
        if constexpr (hasUnknownsBlock<K>()) {
          addUnknownsBlock<K>(m, Ba, t, cs, a.shape_functions[q + nq * nb], w);
        }
      }
      o.template add<V>(e, nd, na, nb, m);
    });
  }  // end of launchStiffnessMatricesKernel

  /*!
   * \brief launches the kernel for the given kind of gradients
   * \tparam Output: destination of the blocks
   * \tparam diagonal: if true, only the pairs of a node with itself are
   * computed
   */
  template <typename Output, bool diagonal = false>
  static void dispatch(const GradientsKind k,
                       const StiffnessMatricesKernelArguments& a,
                       const Output& o,
                       const bool use_device) {
    switch (k) {
      case GradientsKind::SMALL_STRAIN_3D:
        launchStiffnessMatricesKernel<GradientsKind::SMALL_STRAIN_3D, Output,
                                      diagonal>(a, o, use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchStiffnessMatricesKernel<GradientsKind::SMALL_STRAIN_2D, Output,
                                      diagonal>(a, o, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchStiffnessMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_3D,
                                      Output, diagonal>(a, o, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchStiffnessMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_2D,
                                      Output, diagonal>(a, o, use_device);
        break;
      case GradientsKind::TEMPERATURE_GRADIENT_3D:
        launchStiffnessMatricesKernel<GradientsKind::TEMPERATURE_GRADIENT_3D,
                                      Output, diagonal>(a, o, use_device);
        break;
      case GradientsKind::TEMPERATURE_GRADIENT_2D:
        launchStiffnessMatricesKernel<GradientsKind::TEMPERATURE_GRADIENT_2D,
                                      Output, diagonal>(a, o, use_device);
        break;
      case GradientsKind::TEMPERATURE_3D:
        launchStiffnessMatricesKernel<GradientsKind::TEMPERATURE_3D, Output,
                                      diagonal>(a, o, use_device);
        break;
      case GradientsKind::TEMPERATURE_2D:
        launchStiffnessMatricesKernel<GradientsKind::TEMPERATURE_2D, Output,
                                      diagonal>(a, o, use_device);
        break;
      case GradientsKind::MICROMORPHIC_DAMAGE_3D:
        launchStiffnessMatricesKernel<GradientsKind::MICROMORPHIC_DAMAGE_3D,
                                      Output, diagonal>(a, o, use_device);
        break;
      case GradientsKind::MICROMORPHIC_DAMAGE_2D:
        launchStiffnessMatricesKernel<GradientsKind::MICROMORPHIC_DAMAGE_2D,
                                      Output, diagonal>(a, o, use_device);
        break;
    }
  }  // end of dispatch

  void launchStiffnessMatricesKernel(const GradientsKind k,
                                     const StiffnessMatricesKernelArguments& a,
                                     const bool use_device) {
    dispatch(k, a, ElementMatricesStorage{a.matrices}, use_device);
  }  // end of launchStiffnessMatricesKernel

  void launchStiffnessMatricesKernel(const GradientsKind k,
                                     const StiffnessMatricesKernelArguments& a,
                                     const SparseMatrixOutput& o,
                                     const bool use_device) {
    // the blocks are added by a single thread on the host
    dispatch(k, a, SparseMatrixAssembly{.o = o, .atomic = use_device},
             use_device);
  }  // end of launchStiffnessMatricesKernel

  void launchStiffnessMatricesDiagonalsKernel(
      const GradientsKind k,
      const StiffnessMatricesKernelArguments& a,
      real* const d,
      const bool use_device) {
    dispatch<ElementDiagonalsStorage, true>(k, a, ElementDiagonalsStorage{d},
                                            use_device);
  }  // end of launchStiffnessMatricesDiagonalsKernel

  /*!
   * \brief launches the kernel adding the action of the stiffness matrices of
   * the elements of a material on values given element by element (see
   * `launchTangentOperatorActionKernel`): one thread per element, which
   * computes, at each point, the increment of the gradients `B x`, the
   * increment of the thermodynamic forces `K B x` (plus the contribution of
   * the increment of the unknowns `N x` for the heat transfer behaviours)
   * and adds `Bᵀ K B x w det J` to the values of the element
   * \tparam K: kind of the gradients
   */
  template <GradientsKind K>
  static void launchTangentOperatorActionKernel(
      const StiffnessMatricesKernelArguments& a,
      const real* const x,
      real* const y,
      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a, x, y] MFEM_HOST_DEVICE(int k) {
          // dimension and number of components of the unknowns
          constexpr auto N = getDimension<K>();
          constexpr auto V = getUnknownsSize<K>();
          constexpr auto gs = static_cast<std::size_t>(getGradientsSize<K>());
          const auto nd = static_cast<std::size_t>(a.nd);
          const auto nq = static_cast<std::size_t>(a.nq);
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto* const xe = x + nd * V * e;
          auto* const ye = y + nd * V * e;
          for (std::size_t q = 0; q != nq; ++q) {
            const auto detJ = a.determinants[q + nq * e];
            real iJ[N * N];
            invertJacobian<N>(iJ, a.jacobians + q + nq * N * N * e, nq, detJ);
            // increment of the gradients, B x
            real dg[gs] = {};
            for (std::size_t node = 0; node != nd; ++node) {
              real dN[N];
              computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, q,
                                            nq, node, iJ);
              real B[9][3];
              computeB<K>(
                  B, dN,
                  getShapeFunctionValue<K>(a.shape_functions, q, nq, node));
              for (std::size_t r = 0; r != gs; ++r) {
                for (std::size_t c = 0; c != V; ++c) {
                  dg[r] += B[r][c] * xe[node + nd * c];
                }
              }
            }
            // increment of the thermodynamic forces, K B x
            const auto o = static_cast<long long>(a.offsets[k]) +
                           static_cast<long long>(q);
            const auto* const t = a.tangent_operator + o * a.point_stride;
            const auto cs = a.component_stride;
            real ds[gs];
            for (std::size_t r1 = 0; r1 != gs; ++r1) {
              ds[r1] = 0;
              for (std::size_t r2 = 0; r2 != gs; ++r2) {
                ds[r1] +=
                    getTangentOperatorComponent<K>(t, r1, r2, cs) * dg[r2];
              }
            }
            if constexpr (hasUnknownsBlock<K>()) {
              // increment of the unknowns at the point, N x, times the
              // derivatives of the thermodynamic forces with respect to the
              // unknowns
              auto du = real{0};
              for (std::size_t node = 0; node != nd; ++node) {
                du += a.shape_functions[q + nq * node] * xe[node];
              }
              for (std::size_t r1 = 0; r1 != gs; ++r1) {
                ds[r1] += t[static_cast<long long>(gs * gs + r1) * cs] * du;
              }
            }
            // Bᵀ K B x w det J
            const auto w = a.weights[q] * detJ;
            for (std::size_t node = 0; node != nd; ++node) {
              real dN[N];
              computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, q,
                                            nq, node, iJ);
              real B[9][3];
              computeB<K>(
                  B, dN,
                  getShapeFunctionValue<K>(a.shape_functions, q, nq, node));
              for (std::size_t c = 0; c != V; ++c) {
                auto v = real{0};
                for (std::size_t r = 0; r != gs; ++r) {
                  v += B[r][c] * ds[r];
                }
                ye[node + nd * c] += w * v;
              }
            }
          }
        });
  }  // end of launchTangentOperatorActionKernel

  void launchTangentOperatorActionKernel(
      const GradientsKind k,
      const StiffnessMatricesKernelArguments& a,
      const real* const x,
      real* const y,
      const bool use_device) {
    switch (k) {
      case GradientsKind::SMALL_STRAIN_3D:
        launchTangentOperatorActionKernel<GradientsKind::SMALL_STRAIN_3D>(
            a, x, y, use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchTangentOperatorActionKernel<GradientsKind::SMALL_STRAIN_2D>(
            a, x, y, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchTangentOperatorActionKernel<
            GradientsKind::DEFORMATION_GRADIENT_3D>(a, x, y, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchTangentOperatorActionKernel<
            GradientsKind::DEFORMATION_GRADIENT_2D>(a, x, y, use_device);
        break;
      case GradientsKind::TEMPERATURE_GRADIENT_3D:
        launchTangentOperatorActionKernel<
            GradientsKind::TEMPERATURE_GRADIENT_3D>(a, x, y, use_device);
        break;
      case GradientsKind::TEMPERATURE_GRADIENT_2D:
        launchTangentOperatorActionKernel<
            GradientsKind::TEMPERATURE_GRADIENT_2D>(a, x, y, use_device);
        break;
      case GradientsKind::TEMPERATURE_3D:
        launchTangentOperatorActionKernel<GradientsKind::TEMPERATURE_3D>(
            a, x, y, use_device);
        break;
      case GradientsKind::TEMPERATURE_2D:
        launchTangentOperatorActionKernel<GradientsKind::TEMPERATURE_2D>(
            a, x, y, use_device);
        break;
      case GradientsKind::MICROMORPHIC_DAMAGE_3D:
        launchTangentOperatorActionKernel<
            GradientsKind::MICROMORPHIC_DAMAGE_3D>(a, x, y, use_device);
        break;
      case GradientsKind::MICROMORPHIC_DAMAGE_2D:
        launchTangentOperatorActionKernel<
            GradientsKind::MICROMORPHIC_DAMAGE_2D>(a, x, y, use_device);
        break;
    }
  }  // end of launchTangentOperatorActionKernel

  void zeroValues(real* const v, const int n, const bool use_device) {
    mfem::forall_switch(use_device, n,
                        [v] MFEM_HOST_DEVICE(int i) { v[i] = real{0}; });
  }  // end of zeroValues

}  // end of namespace mfem_mgis::batch::internals
