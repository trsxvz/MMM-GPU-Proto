/*!
 * \file   src/Batch/MFEMFieldOperations.cxx
 * \brief  This file implements the `MFEMFieldOperations` class
 * \date   09/10/2026
 * \note   compiled as CUDA code when MFEM is compiled with CUDA
 */

#include <limits>
#include <cstddef>
#include "mfem/general/forall.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Batch/MFEMFieldOperations.hxx"

namespace mfem_mgis::batch {

  /*!
   * \return the number of values of a view, checked against the range of
   * the indices of the kernels of MFEM
   * \param[in] n: number of points
   * \param[in] nc: number of components
   */
  static int getNumberOfValues(const std::size_t n, const std::size_t nc) {
    if ((nc != 0) &&
        (n > static_cast<std::size_t>(std::numeric_limits<int>::max()) / nc)) {
      mgis::raise("MFEMFieldOperations: too many values");
    }
    return static_cast<int>(n * nc);
  }  // end of getNumberOfValues

  /*!
   * \brief copy the values of a view in another one with a kernel of MFEM
   * \param[in] use_device: if true, the kernel runs on the device
   * \param[in] d: destination
   * \param[in] s: source
   * \note nvcc does not allow the extended lambdas of the kernels in the
   * private member functions of `MFEMFieldOperations`
   */
  static void copyFieldView(const bool use_device,
                            const mgis::FieldView<real>& d,
                            const mgis::FieldView<const real>& s) {
    const auto n = d.getNumberOfPoints();
    const auto nc = d.getNumberOfComponents();
    auto* const dv = d.data();
    const auto dps = d.getPointStride();
    const auto dcs = d.getComponentStride();
    const auto* const sv = s.data();
    const auto sps = s.getPointStride();
    const auto scs = s.getComponentStride();
    const auto nv = getNumberOfValues(n, nc);
    // consecutive threads process consecutive values of the destination
    if (dcs <= dps) {
      mfem::forall_switch(use_device, nv, [=] MFEM_HOST_DEVICE(int k) {
        const auto i = static_cast<std::size_t>(k) / nc;
        const auto j = static_cast<std::size_t>(k) % nc;
        dv[i * dps + j * dcs] = sv[i * sps + j * scs];
      });
    } else {
      mfem::forall_switch(use_device, nv, [=] MFEM_HOST_DEVICE(int k) {
        const auto i = static_cast<std::size_t>(k) % n;
        const auto j = static_cast<std::size_t>(k) / n;
        dv[i * dps + j * dcs] = sv[i * sps + j * scs];
      });
    }
  }  // end of copyFieldView

  /*!
   * \brief fill a view with a kernel of MFEM (see `copyFieldView`)
   * \param[in] use_device: if true, the kernel runs on the device
   * \param[in] d: view
   * \param[in] v: value
   */
  static void fillFieldView(const bool use_device,
                            const mgis::FieldView<real>& d,
                            const real v) {
    const auto n = d.getNumberOfPoints();
    const auto nc = d.getNumberOfComponents();
    auto* const dv = d.data();
    const auto ps = d.getPointStride();
    const auto cs = d.getComponentStride();
    const auto nv = getNumberOfValues(n, nc);
    if (cs <= ps) {
      mfem::forall_switch(use_device, nv, [=] MFEM_HOST_DEVICE(int k) {
        const auto i = static_cast<std::size_t>(k) / nc;
        const auto j = static_cast<std::size_t>(k) % nc;
        dv[i * ps + j * cs] = v;
      });
    } else {
      mfem::forall_switch(use_device, nv, [=] MFEM_HOST_DEVICE(int k) {
        const auto i = static_cast<std::size_t>(k) % n;
        const auto j = static_cast<std::size_t>(k) / n;
        dv[i * ps + j * cs] = v;
      });
    }
  }  // end of fillFieldView

  MFEMFieldOperations::MFEMFieldOperations(const bool d) noexcept
      : use_device(d) {}  // end of MFEMFieldOperations

  void MFEMFieldOperations::do_copy(const mgis::FieldView<real>& d,
                                    const mgis::FieldView<const real>& s) {
    if (!this->use_device) {
      // dense copies on the host
      auto ops = mgis::HostFieldOperations{};
      ops.copy(d, s);
      return;
    }
    copyFieldView(this->use_device, d, s);
  }  // end of do_copy

  void MFEMFieldOperations::do_fill(const mgis::FieldView<real>& d,
                                    const real v) {
    if (!this->use_device) {
      auto ops = mgis::HostFieldOperations{};
      ops.fill(d, v);
      return;
    }
    fillFieldView(this->use_device, d, v);
  }  // end of do_fill

  MFEMFieldOperations::~MFEMFieldOperations() = default;

}  // end of namespace mfem_mgis::batch
