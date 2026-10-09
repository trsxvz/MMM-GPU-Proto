/*!
 * \file   MFEMFieldOperationsTest.cxx
 * \brief  Test of the copies and fills of views of MGIS by kernels of MFEM
 * (`MFEMFieldOperations`), on the host and on the device selected by the
 * `--device` option, in both layouts
 * \date   09/10/2026
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Batch/MFEMFieldOperations.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \return a view of an array, on the device or on the host
 * \param[in] v: array
 * \param[in] d: if true, view of the copy on the device
 * \param[in] n: number of points
 * \param[in] nc: number of components
 * \param[in] ps: stride between two points
 * \param[in] cs: stride between two components
 */
static mgis::FieldView<mfem_mgis::real> getView(mfem::Vector& v,
                                                const bool d,
                                                const mfem_mgis::size_type n,
                                                const mfem_mgis::size_type nc,
                                                const mfem_mgis::size_type ps,
                                                const mfem_mgis::size_type cs) {
  auto* const p = d ? v.ReadWrite() : v.HostReadWrite();
  return {std::span<mfem_mgis::real>(p, static_cast<std::size_t>(v.Size())), n,
          nc, ps, cs};
}  // end of getView

static void test(const bool d) {
  using namespace mfem_mgis;
  const auto where = std::string{d ? "device" : "host"};
  constexpr size_type n = 1000;
  constexpr size_type nc = 6;
  auto ops = batch::MFEMFieldOperations(d);
  auto a = mfem::Vector(n * nc);
  auto b = mfem::Vector(n * nc);
  a.UseDevice(true);
  b.UseDevice(true);
  auto* const ha = a.HostWrite();
  for (size_type i = 0; i != n * nc; ++i) {
    ha[i] = static_cast<real>(i);
  }
  // copy from the point by point layout to the component by component one
  ops.copy(getView(b, d, n, nc, 1, n), getView(a, d, n, nc, nc, 1));
  const auto* hb = b.HostRead();
  auto ok = true;
  for (size_type i = 0; i != n; ++i) {
    for (size_type j = 0; j != nc; ++j) {
      ok = ok && (hb[j * n + i] == static_cast<real>(i * nc + j));
    }
  }
  check(ok, where + ": copy between layouts");
  // copy of a uniform value (stride between points equal to zero) to the
  // component by component layout
  ops.copy(getView(b, d, n, nc, 1, n), getView(a, d, n, nc, 0, 1));
  hb = b.HostRead();
  ok = true;
  for (size_type i = 0; i != n; ++i) {
    for (size_type j = 0; j != nc; ++j) {
      ok = ok && (hb[j * n + i] == static_cast<real>(j));
    }
  }
  check(ok, where + ": copy of a uniform value");
  // fill of the second component in the point by point layout
  {
    auto* const p = d ? b.ReadWrite() : b.HostReadWrite();
    ops.fill(
        mgis::FieldView<real>{
            std::span<real>(p + 1, static_cast<std::size_t>(n * nc - 1)), n, 1,
            nc, 1},
        real{-2});
  }
  hb = b.HostRead();
  ok = true;
  for (size_type k = 0; k != n * nc; ++k) {
    ok = ok && ((k % nc == 1) ? (hb[k] == -2) : (hb[k] != -2));
  }
  check(ok, where + ": fill of one component");
  // copy from the component by component layout to the point by point one
  auto c = mfem::Vector(n * nc);
  c.UseDevice(true);
  auto* const hc = c.HostWrite();
  for (size_type i = 0; i != n * nc; ++i) {
    hc[i] = static_cast<real>(-1);
  }
  auto* const hb2 = b.HostWrite();
  for (size_type i = 0; i != n; ++i) {
    for (size_type j = 0; j != nc; ++j) {
      hb2[j * n + i] = static_cast<real>(100 * i + j);
    }
  }
  ops.copy(getView(c, d, n, nc, nc, 1), getView(b, d, n, nc, 1, n));
  const auto* hc2 = c.HostRead();
  ok = true;
  for (size_type i = 0; i != n; ++i) {
    for (size_type j = 0; j != nc; ++j) {
      ok = ok && (hc2[i * nc + j] == static_cast<real>(100 * i + j));
    }
  }
  check(ok, where + ": copy from the component by component layout");
  // fill of the third component in the component by component layout (the
  // values of a component are contiguous: stride between points equal to
  // one)
  {
    auto* const p = d ? b.ReadWrite() : b.HostReadWrite();
    ops.fill(mgis::FieldView<real>{std::span<real>(p + 2 * n,
                                                   static_cast<std::size_t>(n)),
                                   n, 1, 1, n},
             real{-3});
  }
  hb = b.HostRead();
  ok = true;
  for (size_type i = 0; i != n; ++i) {
    for (size_type j = 0; j != nc; ++j) {
      ok = ok && ((j == 2) ? (hb[j * n + i] == -3)
                           : (hb[j * n + i] == static_cast<real>(100 * i + j)));
    }
  }
  check(ok, where +
                ": fill of one component in the component by component "
                "layout");
  // fill of all the components of a view in the component by component
  // layout (the points are processed first by the kernels)
  ops.fill(getView(b, d, n, nc, 1, n), real{5});
  hb = b.HostRead();
  ok = true;
  for (size_type k = 0; k != n * nc; ++k) {
    ok = ok && (hb[k] == 5);
  }
  check(ok, where + ": fill in the component by component layout");
}  // end of test

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  test(false);
  test(true);
  std::cout << (failures == 0 ? "success" : "failure") << " (device "
            << (mfem::Device::IsEnabled() ? "enabled" : "disabled") << ")\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
