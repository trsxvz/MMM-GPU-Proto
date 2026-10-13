/*!
 * \file   src/Batch/StiffnessMatricesKernels.cxx
 * \brief  This file implements the kernels computing the stiffness matrices
 * of the elements of the standard behaviour integrators (mechanics, heat
 * transfer, micromorphic damage). It contains kernels of MFEM: it is
 * compiled as CUDA code when MFEM is compiled with CUDA.
 * \date   09/10/2026
 */

#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
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
        const auto begin = this->o.I[row];
        const auto end = this->o.I[row + 1];
        for (std::size_t d = 0; d != N; ++d) {
          const auto p =
              this->o.findPosition(e, ndofs, na + nd * c, nb + nd * d, begin,
                                   end, dofs[nb + nd * d]);
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
   * \param[in] xa: derivatives of the shape function of the first node in
   * the physical frame, followed by its value (see `getBOperatorEntries`)
   * \param[in] t: tangent operator blocks of the point
   * \param[in] cs: stride between two components of the blocks
   * \param[in] Nb: value of the shape function of the second node
   * \param[in] w: weight of the point
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline void addUnknownsBlock(
      real (&m)[getUnknownsSize<K>()][getUnknownsSize<K>()],
      const real* const xa,
      const real* const t,
      const long long cs,
      const real Nb,
      const real w) {
    if constexpr (hasUnknownsBlock<K>()) {
      constexpr auto gs = getGradientsSize<K>();
      constexpr auto entries = getBOperatorEntries<K>();
      auto v = real{0};
      for (const auto& e : entries) {
        if (e.column == 0) {
          v += (e.factor * xa[e.source]) *
               t[static_cast<long long>(gs * gs + e.row) * cs];
        }
      }
      m[0][0] += w * v * Nb;
    }
  }  // end of addUnknownsBlock

  /*!
   * \brief product of the tangent operator blocks at a point by the column
   * `d` of the block of a node of the operator `B`, `y = K B[:, d]`, or `y =
   * α K B[:, d] + v (c · B[:, d])` if the blocks are modified (see
   * `multiplyByTangentOperator`), without forming `B`: the null entries of
   * `B` are skipped (see `getBOperatorEntries`), the other terms being
   * added in the same order
   * \tparam K: kind of the gradients
   * \param[out] y: product (one value per component of the gradients)
   * \param[in] a: arguments of the kernels
   * \param[in] o: offset of the point in the material
   * \param[in] x: derivatives of the shape function of the node in the
   * physical frame, followed by its value
   * \param[in] d: column of `B` (component of the unknowns)
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline void multiplyByTangentOperatorAndB(
      real* const y,
      const StiffnessMatricesKernelArguments& a,
      const long long o,
      const real* const x,
      const std::size_t d) {
    constexpr auto gs = static_cast<std::size_t>(getGradientsSize<K>());
    constexpr auto entries = getBOperatorEntries<K>();
    const auto* const t = a.tangent_operator + o * a.point_stride;
    const auto cs = a.component_stride;
    for (std::size_t r1 = 0; r1 != gs; ++r1) {
      y[r1] = 0;
      for (const auto& e : entries) {
        if (static_cast<std::size_t>(e.column) == d) {
          const auto r2 = static_cast<std::size_t>(e.row);
          y[r1] += getTangentOperatorComponent<K>(t, r1, r2, cs) *
                   (e.factor * x[e.source]);
        }
      }
    }
    constexpr auto deformation_gradient =
        (K == GradientsKind::DEFORMATION_GRADIENT_3D) ||
        (K == GradientsKind::DEFORMATION_GRADIENT_2D);
    if (deformation_gradient && (a.modification_coefficients != nullptr)) {
      const auto so = static_cast<std::size_t>(o);
      const auto* const v = a.modification_v + gs * so;
      const auto* const c = a.modification_c + gs * so;
      const auto alpha = a.modification_coefficients[so];
      auto cx = real{0};
      for (const auto& e : entries) {
        if (static_cast<std::size_t>(e.column) == d) {
          cx += c[e.row] * (e.factor * x[e.source]);
        }
      }
      for (std::size_t r1 = 0; r1 != gs; ++r1) {
        y[r1] = alpha * y[r1] + v[r1] * cx;
      }
    }
  }  // end of multiplyByTangentOperatorAndB

  /*!
   * \return the product of the transpose of the column `c` of the block of a
   * node of the operator `B` by a vector, `B[:, c] · y`, without forming
   * `B` (see `getBOperatorEntries`)
   * \tparam K: kind of the gradients
   * \param[in] x: derivatives of the shape function of the node in the
   * physical frame, followed by its value
   * \param[in] c: column of `B` (component of the unknowns)
   * \param[in] y: vector (one value per component of the gradients)
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline real multiplyByTransposedB(const real* const x,
                                                     const std::size_t c,
                                                     const real* const y) {
    constexpr auto entries = getBOperatorEntries<K>();
    auto v = real{0};
    for (const auto& e : entries) {
      if (static_cast<std::size_t>(e.column) == c) {
        v += (e.factor * x[e.source]) * y[e.row];
      }
    }
    return v;
  }  // end of multiplyByTransposedB

  /*!
   * \brief product of the tangent operator blocks at a point by a vector,
   * `y = K x`, or `y = α K x + v (c · x)` if the blocks are modified (see
   * `TangentOperatorModification`)
   * \tparam K: kind of the gradients
   * \param[out] y: product
   * \param[in] a: arguments of the kernels
   * \param[in] o: offset of the point in the material
   * \param[in] x: vector (one value per component of the gradients)
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline void multiplyByTangentOperator(
      real* const y,
      const StiffnessMatricesKernelArguments& a,
      const long long o,
      const real* const x) {
    constexpr auto gs = static_cast<std::size_t>(getGradientsSize<K>());
    const auto* const t = a.tangent_operator + o * a.point_stride;
    const auto cs = a.component_stride;
    for (std::size_t r1 = 0; r1 != gs; ++r1) {
      y[r1] = 0;
      for (std::size_t r2 = 0; r2 != gs; ++r2) {
        y[r1] += getTangentOperatorComponent<K>(t, r1, r2, cs) * x[r2];
      }
    }
    // the tangent operator blocks are only modified for the deformation
    // gradients (FBar formulation)
    constexpr auto deformation_gradient =
        (K == GradientsKind::DEFORMATION_GRADIENT_3D) ||
        (K == GradientsKind::DEFORMATION_GRADIENT_2D);
    if (deformation_gradient && (a.modification_coefficients != nullptr)) {
      const auto so = static_cast<std::size_t>(o);
      const auto* const v = a.modification_v + gs * so;
      const auto* const c = a.modification_c + gs * so;
      const auto alpha = a.modification_coefficients[so];
      auto cx = real{0};
      for (std::size_t r = 0; r != gs; ++r) {
        cx += c[r] * x[r];
      }
      for (std::size_t r1 = 0; r1 != gs; ++r1) {
        y[r1] = alpha * y[r1] + v[r1] * cx;
      }
    }
  }  // end of multiplyByTangentOperator

  /*!
   * \brief version of `launchStiffnessMatricesKernel` used without device:
   * one iteration per element, the physical derivatives of the shape functions
   * and the products `K B` of the nodes being computed once per point instead
   * of once per pair of nodes. The operations on each block are made in the
   * same order as by the kernel: the results are identical.
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
    /*!
     * \brief derivatives of the shape function of a node at a point in the
     * physical frame, followed by its value (see `getBOperatorEntries`), and
     * the products of the tangent operator blocks by the columns of the
     * operator `B` of the node (`KB[d]` for the column `d`)
     */
    struct NodeOperators {
      real x[N + 1];
      real KB[V][getGradientsSize<K>()];
    };
    const auto nd = static_cast<std::size_t>(a.nd);
    const auto nq = static_cast<std::size_t>(a.nq);
    const auto cs = a.component_stride;
    auto operators = std::vector<NodeOperators>(nq * nd);
    auto weights = std::vector<real>(nq);
    auto tangent_operators = std::vector<const real*>(nq);
    for (int k = 0; k != a.number_of_elements; ++k) {
      const auto e = static_cast<std::size_t>(a.elements[k]);
      for (std::size_t q = 0; q != nq; ++q) {
        const auto detJ = a.determinants[q + nq * e];
        real iJ[N * N];
        loadInverseJacobian<N>(iJ, a.inverse_jacobians + q + nq * N * N * e,
                               nq);
        const auto o2 =
            static_cast<long long>(a.offsets[k]) + static_cast<long long>(q);
        const auto* const t = a.tangent_operator + o2 * a.point_stride;
        weights[q] = a.weights[q] * detJ;
        tangent_operators[q] = t;
        for (std::size_t node = 0; node != nd; ++node) {
          auto& op = operators[q * nd + node];
          computePhysicalDerivatives<N>(op.x, a.shape_functions_gradients, q,
                                        nq, node, iJ);
          op.x[N] = getShapeFunctionValue<K>(a.shape_functions, q, nq, node);
          for (std::size_t d = 0; d != V; ++d) {
            multiplyByTangentOperatorAndB<K>(op.KB[d], a, o2, op.x, d);
          }
        }
      }
      const auto npairs = diagonal ? nd : nd * nd;
      for (std::size_t r = 0; r != npairs; ++r) {
        const auto na = diagonal ? r : r % nd;
        const auto nb = diagonal ? r : r / nd;
        real m[V][V] = {};
        for (std::size_t q = 0; q != nq; ++q) {
          const auto& xa = operators[q * nd + na].x;
          const auto& KB = operators[q * nd + nb].KB;
          const auto w = weights[q];
          for (std::size_t d = 0; d != V; ++d) {
            if constexpr (diagonal) {
              // only the diagonal of the block is used
              m[d][d] += w * multiplyByTransposedB<K>(xa, d, KB[d]);
            } else {
              for (std::size_t c = 0; c != V; ++c) {
                m[c][d] += w * multiplyByTransposedB<K>(xa, c, KB[d]);
              }
            }
          }
          if constexpr (hasUnknownsBlock<K>()) {
            addUnknownsBlock<K>(m, xa, tangent_operators[q], cs,
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
        loadInverseJacobian<N>(iJ, a.inverse_jacobians + q + nq * N * N * e,
                               nq);
        // derivatives of the shape functions in the physical frame, followed
        // by their values (see `getBOperatorEntries`)
        real xa[N + 1];
        real xb[N + 1];
        computePhysicalDerivatives<N>(xa, a.shape_functions_gradients, q, nq,
                                      na, iJ);
        computePhysicalDerivatives<N>(xb, a.shape_functions_gradients, q, nq,
                                      nb, iJ);
        xa[N] = getShapeFunctionValue<K>(a.shape_functions, q, nq, na);
        xb[N] = getShapeFunctionValue<K>(a.shape_functions, q, nq, nb);
        const auto o2 =
            static_cast<long long>(a.offsets[k]) + static_cast<long long>(q);
        // tangent operator blocks of the point, read by `addUnknownsBlock`
        [[maybe_unused]] const auto* const t =
            a.tangent_operator + o2 * a.point_stride;
        [[maybe_unused]] const auto cs = a.component_stride;
        const auto w = a.weights[q] * detJ;
        // K B_b, then B_aᵀ (K B_b), column by column, without forming B
        for (std::size_t d = 0; d != V; ++d) {
          real KB[getGradientsSize<K>()];
          multiplyByTangentOperatorAndB<K>(KB, a, o2, xb, d);
          if constexpr (diagonal) {
            // only the diagonal of the block is used
            m[d][d] += w * multiplyByTransposedB<K>(xa, d, KB);
          } else {
            for (std::size_t c = 0; c != V; ++c) {
              m[c][d] += w * multiplyByTransposedB<K>(xa, c, KB);
            }
          }
        }
        if constexpr (hasUnknownsBlock<K>()) {
          addUnknownsBlock<K>(m, xa, t, cs, a.shape_functions[q + nq * nb], w);
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

  void launchStiffnessMatricesByEntriesKernel(const SparseMatrixOutput& o,
                                              const int nnz,
                                              const int* const offsets,
                                              const int* const contributions,
                                              const real* const k,
                                              const int n,
                                              const bool use_device) {
    mfem::forall_switch(
        use_device, nnz,
        [o, offsets, contributions, k, n] MFEM_HOST_DEVICE(int p) {
          const auto sn = static_cast<std::size_t>(n);
          const auto column = o.J[p];
          auto v = o.values[p];
          for (auto t = offsets[p]; t != offsets[p + 1]; ++t) {
            // entry (i, j) of the stiffness matrix of the element e, stored
            // column by column (see `ElementMatricesStorage`)
            const auto q = static_cast<std::size_t>(contributions[t]);
            const auto e = q / (sn * sn);
            const auto j = (q - e * sn * sn) / sn;
            if (o.element_dofs[sn * e + j] != column) {
              ::AtomicAdd(o.missing[0], 1);
              continue;
            }
            v += k[q];
          }
          o.values[p] = v;
        });
  }  // end of launchStiffnessMatricesByEntriesKernel

  void launchStiffnessMatricesPositionsKernel(const SparseMatrixOutput& o,
                                              const int ne,
                                              const int n,
                                              std::uint16_t* const positions,
                                              const bool use_device) {
    mfem::forall_switch(
        use_device, ne * n, [o, n, positions] MFEM_HOST_DEVICE(int t) {
          const auto sn = static_cast<std::size_t>(n);
          const auto st = static_cast<std::size_t>(t);
          const auto e = st / sn;
          const auto i = st - e * sn;
          const auto* const dofs = o.element_dofs + sn * e;
          auto* const p = positions + sn * st;
          const auto row = dofs[i];
          const auto valid = (row >= 0) && (row < o.number_of_rows);
          const auto begin = valid ? o.I[row] : 0;
          const auto end = valid ? o.I[row + 1] : 0;
          for (std::size_t j = 0; j != sn; ++j) {
            const auto c = valid ? findColumn(o.J, begin, end, dofs[j]) : end;
            p[j] = static_cast<std::uint16_t>(
                ((c != end) && (c - begin < unknown_position))
                    ? c - begin
                    : unknown_position);
          }
        });
  }  // end of launchStiffnessMatricesPositionsKernel

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
            loadInverseJacobian<N>(iJ, a.inverse_jacobians + q + nq * N * N * e,
                                   nq);
            // increment of the gradients, B x, without forming B (see
            // `getBOperatorEntries`)
            constexpr auto entries = getBOperatorEntries<K>();
            real dg[gs] = {};
            for (std::size_t node = 0; node != nd; ++node) {
              real xn[N + 1];
              computePhysicalDerivatives<N>(xn, a.shape_functions_gradients, q,
                                            nq, node, iJ);
              xn[N] = getShapeFunctionValue<K>(a.shape_functions, q, nq, node);
              for (const auto& en : entries) {
                const auto c = static_cast<std::size_t>(en.column);
                dg[en.row] += (en.factor * xn[en.source]) * xe[node + nd * c];
              }
            }
            // increment of the thermodynamic forces, K B x
            const auto o = static_cast<long long>(a.offsets[k]) +
                           static_cast<long long>(q);
            // tangent operator blocks of the point, for the derivatives with
            // respect to the unknowns
            [[maybe_unused]] const auto* const t =
                a.tangent_operator + o * a.point_stride;
            [[maybe_unused]] const auto cs = a.component_stride;
            real ds[gs];
            multiplyByTangentOperator<K>(ds, a, o, dg);
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
              real xn[N + 1];
              computePhysicalDerivatives<N>(xn, a.shape_functions_gradients, q,
                                            nq, node, iJ);
              xn[N] = getShapeFunctionValue<K>(a.shape_functions, q, nq, node);
              for (std::size_t c = 0; c != V; ++c) {
                ye[node + nd * c] += w * multiplyByTransposedB<K>(xn, c, ds);
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
