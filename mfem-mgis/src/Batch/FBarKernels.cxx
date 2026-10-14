/*!
 * \file   src/Batch/FBarKernels.cxx
 * \brief  This file implements the kernels of the FBar formulation. It contains
 * kernels of MFEM: it is compiled as CUDA or HIP code when MFEM is compiled
 * with CUDA or HIP, without contraction of the floating point operations into
 * fused multiply-add operations in single precision (see
 * `src/CMakeLists.txt`).
 * \date   11/10/2026
 */

#include <cmath>
#include <cstddef>
#include "mfem/general/forall.hpp"
#include "KernelsUtilities.hxx"
#include "FBarKernels.hxx"

namespace mfem_mgis::batch::internals {

  /*!
   * \return the determinant of a deformation gradient, in the order of the
   * gradients of the behaviours (`F11, F22, F33, F12, F21, F13, F31, F23,
   * F32` in 3D, `F11, F22, F33, F12, F21` in 2D)
   * \tparam K: kind of the gradients
   * \param[in] F: deformation gradient
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline real computeDeterminant(const real* const F) {
    if constexpr (K == GradientsKind::DEFORMATION_GRADIENT_3D) {
      return F[0] * (F[1] * F[2] - F[7] * F[8]) -
             F[3] * (F[4] * F[2] - F[7] * F[6]) +
             F[5] * (F[4] * F[8] - F[1] * F[6]);
    } else {
      static_assert(K == GradientsKind::DEFORMATION_GRADIENT_2D);
      return F[2] * (F[0] * F[1] - F[3] * F[4]);
    }
  }  // end of computeDeterminant

  /*!
   * \brief computes the derivative of the determinant of a deformation
   * gradient with respect to its components (cofactors), in the order of
   * the gradients
   * \tparam K: kind of the gradients
   * \param[out] dJ: derivative
   * \param[in] F: deformation gradient
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline void computeDeterminantDerivative(
      real* const dJ, const real* const F) {
    if constexpr (K == GradientsKind::DEFORMATION_GRADIENT_3D) {
      // F11 F22 F33 F12 F21 F13 F31 F23 F32
      // 0   1   2   3   4   5   6   7   8
      dJ[0] = F[1] * F[2] - F[7] * F[8];
      dJ[1] = F[0] * F[2] - F[5] * F[6];
      dJ[2] = F[0] * F[1] - F[3] * F[4];
      dJ[3] = F[7] * F[6] - F[4] * F[2];
      dJ[4] = F[5] * F[8] - F[3] * F[2];
      dJ[5] = F[4] * F[8] - F[1] * F[6];
      dJ[6] = F[3] * F[7] - F[5] * F[1];
      dJ[7] = F[3] * F[6] - F[0] * F[8];
      dJ[8] = F[5] * F[4] - F[0] * F[7];
    } else {
      static_assert(K == GradientsKind::DEFORMATION_GRADIENT_2D);
      dJ[0] = F[1] * F[2];
      dJ[1] = F[0] * F[2];
      dJ[2] = F[0] * F[1] - F[3] * F[4];
      dJ[3] = -F[4] * F[2];
      dJ[4] = -F[3] * F[2];
    }
  }  // end of computeDeterminantDerivative

  //! \return the inverse of the exponent of the FBar formulation
  template <GradientsKind K>
  MFEM_HOST_DEVICE constexpr int getFBarExponentInverse() {
    return K == GradientsKind::DEFORMATION_GRADIENT_3D ? 3 : 2;
  }  // end of getFBarExponentInverse

  /*!
   * \return `(J0 / J)^(1/ie)`
   * \tparam K: kind of the gradients
   * \param[in] J0: determinant of the deformation gradient at the center of
   * the element
   * \param[in] J: determinant of the deformation gradient at the point
   */
  template <GradientsKind K>
  MFEM_HOST_DEVICE inline real computeFBarCoefficient(const real J0,
                                                      const real J) {
    if constexpr (getFBarExponentInverse<K>() == 3) {
      return cbrt(J0 / J);
    } else {
      return sqrt(J0 / J);
    }
  }  // end of computeFBarCoefficient

  template <GradientsKind K>
  static void launchFBarCentersKernel(const FBarCentersKernelArguments& a,
                                      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a] MFEM_HOST_DEVICE(int k) {
          constexpr auto N = getDimension<K>();
          constexpr auto V = getUnknownsSize<K>();
          constexpr auto gs = getGradientsSize<K>();
          const auto e = static_cast<std::size_t>(a.elements[k]);
          // inverse of the jacobian at the center (one point per element)
          real iJ[N * N];
          loadInverseJacobian<N>(
              iJ, a.inverse_jacobians + static_cast<std::size_t>(N * N) * e, 1);
          // deformation gradient at the center
          real F0[9];
          computeDeformationGradientAtCenter<K>(
              F0, a.macroscopic_gradients,
              (a.derivatives != nullptr)
                  ? a.derivatives + static_cast<std::size_t>(V * N) * e
                  : nullptr,
              iJ);
          real dJ0[9];
          computeDeterminantDerivative<K>(dJ0, F0);
          a.centers[k] = computeDeterminant<K>(F0);
          // derivative of J0 with respect to the unknowns of the element
          const auto nd = static_cast<std::size_t>(a.nd);
          auto* const b = a.b + nd * V * static_cast<std::size_t>(k);
          for (std::size_t node = 0; node != nd; ++node) {
            real dN[N];
            computePhysicalDerivatives<N>(dN, a.shape_functions_gradients, 0, 1,
                                          node, iJ);
            real B[9][3];
            computeB<K>(B, dN, real{0});
            for (int comp = 0; comp != V; ++comp) {
              real v = 0;
              for (int i = 0; i != gs; ++i) {
                v += dJ0[i] * B[i][comp];
              }
              b[node + nd * comp] = v;
            }
          }
        });
  }  // end of launchFBarCentersKernel

  template <GradientsKind K>
  static void launchFBarGradientsKernel(const FBarGradientsKernelArguments& a,
                                        const bool use_device) {
    const auto n = a.number_of_elements * a.nq;
    mfem::forall_switch(use_device, n, [a] MFEM_HOST_DEVICE(int i) {
      constexpr auto gs = getGradientsSize<K>();
      const auto k = i / a.nq;
      const auto q = i - k * a.nq;
      const auto o = static_cast<std::size_t>(a.offsets[k] + q);
      const auto np = static_cast<std::size_t>(a.number_of_points);
      const auto sgs = static_cast<std::size_t>(gs);
      // deformation gradient at the point
      real F[9];
      for (std::size_t c = 0; c != sgs; ++c) {
        if (a.gradients == nullptr) {
          F[c] = a.macroscopic_gradients[c];
        } else {
          F[c] = a.point_by_point ? a.gradients[o * sgs + c]
                                  : a.gradients[c * np + o];
        }
      }
      const auto J0 = a.centers[k];
      const auto alpha =
          computeFBarCoefficient<K>(J0, computeDeterminant<K>(F));
      a.coefficients[o] = alpha;
      for (std::size_t c = 0; c != sgs; ++c) {
        a.deformation_gradients[o * sgs + c] = F[c];
      }
      if (a.gradients != nullptr) {
        for (std::size_t c = 0; c != sgs; ++c) {
          // in plane strain, the third component of F̄ is one
          const auto v =
              ((K == GradientsKind::DEFORMATION_GRADIENT_2D) && (c == 2))
                  ? real{1}
                  : alpha * F[c];
          if (a.point_by_point) {
            a.gradients[o * sgs + c] = v;
          } else {
            a.gradients[c * np + o] = v;
          }
        }
      }
    });
  }  // end of launchFBarGradientsKernel

  template <GradientsKind K>
  static void launchFBarTangentOperatorKernel(
      const FBarTangentOperatorKernelArguments& a, const bool use_device) {
    const auto n = a.number_of_elements * a.nq;
    mfem::forall_switch(use_device, n, [a] MFEM_HOST_DEVICE(int i) {
      constexpr auto gs = getGradientsSize<K>();
      constexpr auto ie = getFBarExponentInverse<K>();
      const auto k = i / a.nq;
      const auto q = i - k * a.nq;
      const auto o = static_cast<std::size_t>(a.offsets[k] + q);
      const auto sgs = static_cast<std::size_t>(gs);
      const auto* const F = a.deformation_gradients + o * sgs;
      const auto alpha = a.coefficients[o];
      const auto J = computeDeterminant<K>(F);
      const auto J0 = a.centers[k];
      // dα/dF; in plane strain, the third component is ignored
      real c[9];
      computeDeterminantDerivative<K>(c, F);
      for (int j = 0; j != gs; ++j) {
        c[j] *= -alpha / (ie * J);
      }
      if constexpr (K == GradientsKind::DEFORMATION_GRADIENT_2D) {
        c[2] = 0;
      }
      // v = K F
      const auto* const Kb =
          a.tangent_operator + static_cast<long long>(o) * a.point_stride;
      auto Kv = [Kb, &a](const int r, const int col) {
        return Kb[static_cast<long long>(r * gs + col) * a.component_stride];
      };
      real v[9];
      for (int r = 0; r != gs; ++r) {
        v[r] = 0;
        for (int col = 0; col != gs; ++col) {
          if ((K == GradientsKind::DEFORMATION_GRADIENT_2D) && (col == 2)) {
            continue;
          }
          v[r] += Kv(r, col) * F[col];
        }
      }
      // vectors of the modification of the tangent operator, `K' = α K + v
      // ⊗ c`, applied by the kernels of the stiffness matrices
      auto* const vo = a.v + o * sgs;
      auto* const co = a.c + o * sgs;
      for (int r = 0; r != gs; ++r) {
        vo[r] = v[r];
        co[r] = c[r];
      }
      auto* const s = a.pseudo_stresses + o * sgs;
      for (int r = 0; r != gs; ++r) {
        s[r] = alpha / (ie * J0) * v[r];
      }
    });
  }  // end of launchFBarTangentOperatorKernel

  void launchFBarCentersKernel(const GradientsKind k,
                               const FBarCentersKernelArguments& a,
                               const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFBarCentersKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, use_device);
    } else {
      launchFBarCentersKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, use_device);
    }
  }  // end of launchFBarCentersKernel

  void launchFBarGradientsKernel(const GradientsKind k,
                                 const FBarGradientsKernelArguments& a,
                                 const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFBarGradientsKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, use_device);
    } else {
      launchFBarGradientsKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, use_device);
    }
  }  // end of launchFBarGradientsKernel

  void launchFBarTangentOperatorKernel(
      const GradientsKind k,
      const FBarTangentOperatorKernelArguments& a,
      const bool use_device) {
    if (k == GradientsKind::DEFORMATION_GRADIENT_3D) {
      launchFBarTangentOperatorKernel<GradientsKind::DEFORMATION_GRADIENT_3D>(
          a, use_device);
    } else {
      launchFBarTangentOperatorKernel<GradientsKind::DEFORMATION_GRADIENT_2D>(
          a, use_device);
    }
  }  // end of launchFBarTangentOperatorKernel

  void launchFBarRankOneMatricesKernel(const FBarRankOneKernelArguments& a,
                                       real* const m,
                                       const bool use_device) {
    const auto n = a.number_of_elements * a.n * a.n;
    mfem::forall_switch(use_device, n, [a, m] MFEM_HOST_DEVICE(int t) {
      const auto nn = a.n * a.n;
      const auto k = t / nn;
      const auto r = t - k * nn;
      // row i (first index of `a`), column j, stored column by column
      const auto j = static_cast<std::size_t>(r / a.n);
      const auto i = static_cast<std::size_t>(r) - j * a.n;
      const auto e = static_cast<std::size_t>(a.elements[k]);
      const auto sk = static_cast<std::size_t>(k);
      const auto sn = static_cast<std::size_t>(a.n);
      m[sn * sn * e + i + sn * j] += a.a[sn * sk + i] * a.b[sn * sk + j];
    });
  }  // end of launchFBarRankOneMatricesKernel

  void launchFBarRankOneMatricesKernel(const FBarRankOneKernelArguments& a,
                                       const SparseMatrixOutput& o,
                                       const bool use_device) {
    const auto n = a.number_of_elements * a.n;
    mfem::forall_switch(use_device, n, [a, o] MFEM_HOST_DEVICE(int t) {
      const auto k = t / a.n;
      const auto i = t - k * a.n;
      const auto e = static_cast<std::size_t>(a.elements[k]);
      const auto sk = static_cast<std::size_t>(k);
      const auto sn = static_cast<std::size_t>(a.n);
      const auto* const dofs = o.element_dofs + sn * e;
      const auto row = dofs[i];
      if ((row < 0) || (row >= o.number_of_rows)) {
        ::AtomicAdd(o.missing[0], a.n);
        return;
      }
      const auto ai = a.a[sn * sk + static_cast<std::size_t>(i)];
      const auto begin = o.I[row];
      const auto end = o.I[row + 1];
      for (int j = 0; j != a.n; ++j) {
        const auto p =
            o.findPosition(e, sn, static_cast<std::size_t>(i),
                           static_cast<std::size_t>(j), begin, end, dofs[j]);
        if (p == end) {
          ::AtomicAdd(o.missing[0], 1);
        } else {
          ::AtomicAdd(o.values[p],
                      ai * a.b[sn * sk + static_cast<std::size_t>(j)]);
        }
      }
    });
  }  // end of launchFBarRankOneMatricesKernel

  void launchFBarRankOneDiagonalsKernel(const FBarRankOneKernelArguments& a,
                                        real* const d,
                                        const bool use_device) {
    const auto n = a.number_of_elements * a.n;
    mfem::forall_switch(use_device, n, [a, d] MFEM_HOST_DEVICE(int t) {
      const auto k = t / a.n;
      const auto i = static_cast<std::size_t>(t - k * a.n);
      const auto sn = static_cast<std::size_t>(a.n);
      const auto e = static_cast<std::size_t>(a.elements[k]);
      const auto p = sn * static_cast<std::size_t>(k) + i;
      d[sn * e + i] += a.a[p] * a.b[p];
    });
  }  // end of launchFBarRankOneDiagonalsKernel

  void launchFBarRankOneActionKernel(const FBarRankOneKernelArguments& a,
                                     const real* const x,
                                     real* const y,
                                     const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a, x, y] MFEM_HOST_DEVICE(int k) {
          const auto sn = static_cast<std::size_t>(a.n);
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto sk = static_cast<std::size_t>(k);
          real bx = 0;
          for (std::size_t j = 0; j != sn; ++j) {
            bx += a.b[sn * sk + j] * x[sn * e + j];
          }
          for (std::size_t i = 0; i != sn; ++i) {
            y[sn * e + i] += a.a[sn * sk + i] * bx;
          }
        });
  }  // end of launchFBarRankOneActionKernel

}  // end of namespace mfem_mgis::batch::internals
