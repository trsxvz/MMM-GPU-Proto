/*!
 * \file   src/Batch/StiffnessMatricesKernels.cxx
 * \brief  This file implements the kernels computing the stiffness matrices
 * of the elements of the standard mechanical behaviour integrators. It
 * contains kernels of MFEM: it is compiled as CUDA code when MFEM is
 * compiled with CUDA.
 * \date   09/10/2026
 */

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
     * \tparam N: dimension
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

  /*!
   * \return the position of a column in a row of a sparse matrix, or `end`
   * if the column is not in the row. The first column is tested first, then
   * the other columns are searched by bisection, which is exact for sorted
   * rows and for rows whose first column is the diagonal one, the others
   * being sorted (see `makeAssembledStiffnessMatrix`). Otherwise, the row is
   * searched linearly.
   * \param[in] J: columns of the matrix
   * \param[in] begin: position of the first column of the row
   * \param[in] end: position after the last column of the row
   * \param[in] c: column
   */
  MFEM_HOST_DEVICE inline int findColumn(const int* const J,
                                         const int begin,
                                         const int end,
                                         const int c) {
    if (begin == end) {
      return end;
    }
    if (J[begin] == c) {
      return begin;
    }
    auto lo = begin + 1;
    auto hi = end;
    while (lo != hi) {
      const auto mid = lo + (hi - lo) / 2;
      if (J[mid] < c) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    if ((lo != end) && (J[lo] == c)) {
      return lo;
    }
    for (auto p = begin + 1; p != end; ++p) {
      if (J[p] == c) {
        return p;
      }
    }
    return end;
  }  // end of findColumn

  //! \brief adds the blocks to an assembled sparse matrix
  struct SparseMatrixAssembly {
    //! \brief assembled matrix
    SparseMatrixOutput o;
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
          if (p != end) {
            ::AtomicAdd(this->o.values[p], m[c][d]);
          } else {
            ::AtomicAdd(this->o.missing[0], 1);
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
   * \brief launches the kernel adding the stiffness matrices of the
   * elements of a material: one thread per pair of nodes of an element,
   * computing the `N x N` block `Σ B_aᵀ K B_b w det J`
   * \tparam K: kind of the gradients
   * \tparam N: dimension
   * \tparam Output: destination of the blocks
   * \tparam diagonal: if true, only the pairs of a node with itself are
   * computed
   */
  template <GradientsKind K, int N, typename Output, bool diagonal = false>
  static void launchStiffnessMatricesKernel(
      const StiffnessMatricesKernelArguments& a,
      const Output o,
      const bool use_device) {
    const auto npairs = diagonal ? a.nd : a.nd * a.nd;
    const auto n = a.number_of_elements * npairs;
    mfem::forall_switch(use_device, n, [a, o, npairs] MFEM_HOST_DEVICE(int i) {
      const auto nd = static_cast<std::size_t>(a.nd);
      const auto nq = static_cast<std::size_t>(a.nq);
      const auto gs = static_cast<std::size_t>(a.gradients_size);
      const auto np = static_cast<std::size_t>(npairs);
      const auto k = static_cast<std::size_t>(i) / np;
      const auto r = static_cast<std::size_t>(i) - k * np;
      const auto na = diagonal ? r : r % nd;
      const auto nb = diagonal ? r : r / nd;
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
        const auto o2 =
            static_cast<long long>(a.offsets[k]) + static_cast<long long>(q);
        const auto* const t = a.tangent_operator + o2 * a.point_stride;
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
      o.template add<N>(e, nd, na, nb, m);
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
        launchStiffnessMatricesKernel<GradientsKind::SMALL_STRAIN_3D, 3, Output,
                                      diagonal>(a, o, use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchStiffnessMatricesKernel<GradientsKind::SMALL_STRAIN_2D, 2, Output,
                                      diagonal>(a, o, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchStiffnessMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_3D, 3,
                                      Output, diagonal>(a, o, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchStiffnessMatricesKernel<GradientsKind::DEFORMATION_GRADIENT_2D, 2,
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
    dispatch(k, a, SparseMatrixAssembly{o}, use_device);
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
   * increment of the thermodynamic forces `K B x` and adds
   * `Bᵀ K B x w det J` to the values of the element
   * \tparam K: kind of the gradients
   * \tparam N: dimension
   */
  template <GradientsKind K, int N>
  static void launchTangentOperatorActionKernel(
      const StiffnessMatricesKernelArguments& a,
      const real* const x,
      real* const y,
      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a, x, y] MFEM_HOST_DEVICE(int k) {
          constexpr auto gs = static_cast<std::size_t>(getGradientsSize<K>());
          const auto nd = static_cast<std::size_t>(a.nd);
          const auto nq = static_cast<std::size_t>(a.nq);
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto* const xe = x + nd * N * e;
          auto* const ye = y + nd * N * e;
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
              computeB<K>(B, dN);
              for (std::size_t r = 0; r != gs; ++r) {
                for (std::size_t c = 0; c != N; ++c) {
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
                ds[r1] += t[static_cast<long long>(r1 * gs + r2) * cs] * dg[r2];
              }
            }
            // Bᵀ K B x w det J
            const auto w = a.weights[q] * detJ;
            for (std::size_t node = 0; node != nd; ++node) {
              real dN[N];
              computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, q,
                                            nq, node, iJ);
              real B[9][3];
              computeB<K>(B, dN);
              for (std::size_t c = 0; c != N; ++c) {
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
        launchTangentOperatorActionKernel<GradientsKind::SMALL_STRAIN_3D, 3>(
            a, x, y, use_device);
        break;
      case GradientsKind::SMALL_STRAIN_2D:
        launchTangentOperatorActionKernel<GradientsKind::SMALL_STRAIN_2D, 2>(
            a, x, y, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_3D:
        launchTangentOperatorActionKernel<
            GradientsKind::DEFORMATION_GRADIENT_3D, 3>(a, x, y, use_device);
        break;
      case GradientsKind::DEFORMATION_GRADIENT_2D:
        launchTangentOperatorActionKernel<
            GradientsKind::DEFORMATION_GRADIENT_2D, 2>(a, x, y, use_device);
        break;
    }
  }  // end of launchTangentOperatorActionKernel

  void zeroValues(real* const v, const int n, const bool use_device) {
    mfem::forall_switch(use_device, n,
                        [v] MFEM_HOST_DEVICE(int i) { v[i] = real{0}; });
  }  // end of zeroValues

}  // end of namespace mfem_mgis::batch::internals
