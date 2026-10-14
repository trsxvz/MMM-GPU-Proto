/*!
 * \file   src/Batch/IntegrationPointsFunctionsKernels.cxx
 * \brief  This file implements the kernels computing integrals, right-hand
 * sides of L2 projections and local projections of functions defined at the
 * integration points of a material, and mass matrices of the elements. It
 * contains kernels of MFEM: it is compiled as CUDA or HIP code when MFEM is
 * compiled with CUDA or HIP.
 * \date   11/10/2026
 */

#include <cmath>
#include <cstddef>
#include "mfem/general/forall.hpp"
#include "mfem/general/reducers.hpp"
#include "MFEMMGIS/Batch/IntegrationPointsFunctions.hxx"
#include "IntegrationPointsFunctionsKernels.hxx"

namespace mfem_mgis::batch::internals {

  double launchIntegralKernel(
      const IntegrationPointsFunctionsKernelArguments& a,
      mfem::Array<double>& w,
      const bool use_device) {
    auto r = double{};
    mfem::reduce(
        a.number_of_elements, r,
        [a] MFEM_HOST_DEVICE(const int k, double& s) {
          const auto ke = static_cast<std::size_t>(k);
          const auto e =
              a.compact_geometry ? ke : static_cast<std::size_t>(a.elements[k]);
          const auto nq = static_cast<std::size_t>(a.points.nq);
          for (std::size_t q = 0; q != nq; ++q) {
            const auto detJ = a.points.determinants[q + nq * e];
            if (a.values == nullptr) {
              s += a.points.weights[q] * detJ;
            } else {
              s += a.points.weights[q] * (detJ * a.values[q + nq * ke]);
            }
          }
        },
        mfem::SumReducer<double>{}, use_device, w);
    return r;
  }  // end of launchIntegralKernel

  void launchL2ProjectionRightHandSideKernel(
      const IntegrationPointsFunctionsKernelArguments& a,
      const int c,
      real* const r,
      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a, c, r] MFEM_HOST_DEVICE(int k) {
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto nq = static_cast<std::size_t>(a.points.nq);
          const auto nd = static_cast<std::size_t>(a.nd);
          const auto nc = static_cast<std::size_t>(a.nc);
          const auto ke = static_cast<std::size_t>(k);
          auto* const re = r + nd * e;
          for (std::size_t q = 0; q != nq; ++q) {
            const auto f =
                a.values[static_cast<std::size_t>(c) + nc * (q + nq * ke)];
            // same operations as the integrators of the L2 projection
            const auto w =
                a.points.weights[q] * (a.points.determinants[q + nq * e] * f);
            for (std::size_t i = 0; i != nd; ++i) {
              re[i] += w * a.points.shape_functions[q + nq * i];
            }
          }
        });
  }  // end of launchL2ProjectionRightHandSideKernel

  void launchLocalProjectionsKernel(
      const IntegrationPointsFunctionsKernelArguments& a,
      const IntegrationRuleKernelArguments& m,
      real* const r,
      const bool use_device) {
    mfem::forall_switch(
        use_device, a.number_of_elements, [a, m, r] MFEM_HOST_DEVICE(int k) {
          constexpr auto N =
              static_cast<std::size_t>(maximumNumberOfNodesOfLocalProjections);
          const auto e = static_cast<std::size_t>(a.elements[k]);
          const auto nq = static_cast<std::size_t>(a.points.nq);
          const auto nqm = static_cast<std::size_t>(m.nq);
          const auto nd = static_cast<std::size_t>(a.nd);
          const auto nc = static_cast<std::size_t>(a.nc);
          const auto ke = static_cast<std::size_t>(k);
          // lower part of the mass matrix of the element, then of its
          // Cholesky factor (row major, N x N)
          real L[N * N];
          for (std::size_t i = 0; i != nd; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
              L[i * N + j] = real{0};
            }
          }
          // the absolute values of the determinants are used for the mass
          // matrix and for the right-hand side, so that the mass matrix is
          // positive definite for elements of negative orientation, the
          // projection being unchanged
          for (std::size_t q = 0; q != nqm; ++q) {
            const auto w = m.weights[q] * fabs(m.determinants[q + nqm * e]);
            for (std::size_t i = 0; i != nd; ++i) {
              const auto wNi = w * m.shape_functions[q + nqm * i];
              for (std::size_t j = 0; j <= i; ++j) {
                L[i * N + j] += wNi * m.shape_functions[q + nqm * j];
              }
            }
          }
          // right-hand sides, stored in the values of the element
          auto* const re = r + nd * nc * e;
          for (std::size_t c = 0; c != nc; ++c) {
            for (std::size_t i = 0; i != nd; ++i) {
              re[i + nd * c] = real{0};
            }
          }
          for (std::size_t q = 0; q != nq; ++q) {
            const auto w =
                a.points.weights[q] * fabs(a.points.determinants[q + nq * e]);
            const auto* const f = a.values + nc * (q + nq * ke);
            for (std::size_t i = 0; i != nd; ++i) {
              const auto wNi = w * a.points.shape_functions[q + nq * i];
              for (std::size_t c = 0; c != nc; ++c) {
                re[i + nd * c] += wNi * f[c];
              }
            }
          }
          // Cholesky factorization of the mass matrix
          for (std::size_t j = 0; j != nd; ++j) {
            auto d = L[j * N + j];
            for (std::size_t p = 0; p != j; ++p) {
              d -= L[j * N + p] * L[j * N + p];
            }
            const auto ljj = sqrt(d);
            L[j * N + j] = ljj;
            for (std::size_t i = j + 1; i != nd; ++i) {
              auto s = L[i * N + j];
              for (std::size_t p = 0; p != j; ++p) {
                s -= L[i * N + p] * L[j * N + p];
              }
              L[i * N + j] = s / ljj;
            }
          }
          // resolutions, component by component
          for (std::size_t c = 0; c != nc; ++c) {
            auto* const x = re + nd * c;
            for (std::size_t i = 0; i != nd; ++i) {
              auto s = x[i];
              for (std::size_t p = 0; p != i; ++p) {
                s -= L[i * N + p] * x[p];
              }
              x[i] = s / L[i * N + i];
            }
            for (std::size_t i = nd; i-- != 0;) {
              auto s = x[i];
              for (std::size_t p = i + 1; p != nd; ++p) {
                s -= L[p * N + i] * x[p];
              }
              x[i] = s / L[i * N + i];
            }
          }
        });
  }  // end of launchLocalProjectionsKernel

  void launchMassMatricesKernel(const IntegrationRuleKernelArguments& m,
                                const int* const p,
                                const int nd,
                                const int ne,
                                real* const r,
                                const bool add,
                                const bool use_device) {
    mfem::forall_switch(
        use_device, ne, [m, p, nd, r, add] MFEM_HOST_DEVICE(int k) {
          const auto e = static_cast<std::size_t>(k);
          const auto nq = static_cast<std::size_t>(m.nq);
          const auto n = static_cast<std::size_t>(nd);
          auto* const re = r + n * n * e;
          if (!add) {
            for (std::size_t i = 0; i != n * n; ++i) {
              re[i] = real{0};
            }
          }
          // native index of a degree of freedom of the mass matrix
          auto native = [p](const std::size_t i) {
            return (p != nullptr) ? static_cast<std::size_t>(p[i]) : i;
          };
          // same operations as the mass integrator of MFEM
          for (std::size_t q = 0; q != nq; ++q) {
            const auto w = m.weights[q] * m.determinants[q + nq * e];
            for (std::size_t j = 0; j != n; ++j) {
              const auto wNj = w * m.shape_functions[q + nq * native(j)];
              for (std::size_t i = 0; i != n; ++i) {
                re[i + n * j] += wNj * m.shape_functions[q + nq * native(i)];
              }
            }
          }
        });
  }  // end of launchMassMatricesKernel

}  // end of namespace mfem_mgis::batch::internals
