/*!
 * \file   src/Batch/GlobalStiffnessMatrixKernels.cxx
 * \brief  This file implements the kernels assembling the values of the
 * global matrix of a parallel finite element space from the local matrix of
 * the process and from the values received from the other processes. It
 * contains kernels of MFEM: it is compiled as CUDA or HIP code when MFEM is
 * compiled with CUDA or HIP.
 * \date   15/11/2026
 */

#include "mfem/general/forall.hpp"
#include "GlobalStiffnessMatrixKernels.hxx"

namespace mfem_mgis::batch::internals {

  void launchGlobalStiffnessMatrixOnesKernel(real* const v,
                                             const int n,
                                             const bool use_device) {
    mfem::forall_switch(use_device, n,
                        [v] MFEM_HOST_DEVICE(int i) { v[i] = real{1}; });
  }  // end of launchGlobalStiffnessMatrixOnesKernel

  void launchGlobalStiffnessMatrixSentValuesKernel(real* const s,
                                                   const int* const entries,
                                                   const real* const a,
                                                   const int n,
                                                   const bool use_device) {
    mfem::forall_switch(use_device, n, [s, entries, a] MFEM_HOST_DEVICE(int i) {
      s[i] = a[entries[i]];
    });
  }  // end of launchGlobalStiffnessMatrixSentValuesKernel

  void launchGlobalStiffnessMatrixKernel(
      const GlobalStiffnessMatrixKernelArguments& a, const bool use_device) {
    mfem::forall_switch(
        use_device, a.diagonal_block_entries + a.off_diagonal_block_entries,
        [a] MFEM_HOST_DEVICE(int p) {
          auto v = real{0};
          for (auto t = a.offsets[p]; t != a.offsets[p + 1]; ++t) {
            const auto c = a.contributions[t];
            v += (c < a.local_entries) ? a.local_values[c]
                                       : a.received_values[c - a.local_entries];
          }
          if (p < a.diagonal_block_entries) {
            a.diagonal_block[p] = v;
          } else {
            a.off_diagonal_block[p - a.diagonal_block_entries] = v;
          }
        });
  }  // end of launchGlobalStiffnessMatrixKernel

}  // end of namespace mfem_mgis::batch::internals
