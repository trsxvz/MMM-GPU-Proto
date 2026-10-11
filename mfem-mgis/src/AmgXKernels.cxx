/*!
 * \file   src/AmgXKernels.cxx
 * \brief  This file implements the scaling and the conversions of the
 * vectors exchanged with the solvers of AmgX
 * \date   10/10/2026
 * \note   compiled as CUDA code when MFEM is compiled with CUDA
 */

#include "mfem/general/forall.hpp"
#include "mfem/general/reducers.hpp"
#include "AmgXKernels.hxx"

namespace mfem_mgis::internals {

  double getMaximumAbsoluteValue(const mfem::Vector& v,
                                 mfem::Array<double>& w) {
    const auto use_device = v.UseDevice() && mfem::Device::IsEnabled();
    const auto* const pv = v.Read(use_device);
    auto r = 0.;
    mfem::reduce(
        v.Size(), r,
        [=] MFEM_HOST_DEVICE(const int i, double& m) {
          m = fmax(m, fabs(pv[i]));
        },
        mfem::MaxReducer<double>{}, use_device, w);
    return r;
  }  // end of getMaximumAbsoluteValue

  void scale(mfem::Vector& d, const mfem::Vector& s, const double f) {
    const auto n = s.Size();
    const auto use_device = s.UseDevice() && mfem::Device::IsEnabled();
    d.SetSize(n);
    d.UseDevice(s.UseDevice());
    const auto* const ps = s.Read(use_device);
    auto* const pd = d.Write(use_device);
    mfem::forall_switch(use_device, n, [=] MFEM_HOST_DEVICE(const int i) {
      pd[i] = f * ps[i];
    });
  }  // end of scale

  void convertToSinglePrecision(mfem::Array<float>& d,
                                const mfem::Vector& s,
                                const double f) {
    const auto n = s.Size();
    const auto use_device = s.UseDevice() && mfem::Device::IsEnabled();
    d.SetSize(n);
    const auto* const ps = s.Read(use_device);
    auto* const pd = d.Write(use_device);
    mfem::forall_switch(use_device, n, [=] MFEM_HOST_DEVICE(const int i) {
      pd[i] = static_cast<float>(f * ps[i]);
    });
  }  // end of convertToSinglePrecision

  void convertToDoublePrecision(mfem::Vector& d,
                                const mfem::Array<float>& s,
                                const double f) {
    const auto n = s.Size();
    const auto use_device = d.UseDevice() && mfem::Device::IsEnabled();
    d.SetSize(n);
    const auto* const ps = s.Read(use_device);
    auto* const pd = d.Write(use_device);
    mfem::forall_switch(use_device, n, [=] MFEM_HOST_DEVICE(const int i) {
      pd[i] = f * static_cast<double>(ps[i]);
    });
  }  // end of convertToDoublePrecision

}  // end of namespace mfem_mgis::internals
