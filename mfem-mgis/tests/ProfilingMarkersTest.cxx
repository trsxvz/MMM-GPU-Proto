/*!
 * \file   ProfilingMarkersTest.cxx
 * \brief  Test of the profiling sections declared as ranges of roctx, the
 * annotation library of ROCm: the markers are set by `mfem_mgis::initialize`
 * if mfem-mgis is built with the option `enable-roctx` and if the device
 * selected by the `--device` option uses a GPU, and only then. With or
 * without markers, the profiling tree is the same.
 * \date   14/11/2026
 */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include "mfem/general/device.hpp"
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/vector.hpp"
#include "MGIS/Context.hxx"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Config.hxx"

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  // `yes` if the profiling sections are expected to be declared to roctx,
  // `no` otherwise
  const char* expected = "no";
  auto args = mfem::OptionsParser(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&expected, "-e", "--expected",
                 "yes if the profiling sections are expected to be declared "
                 "to roctx, no otherwise");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  const auto markers = mgis::getProfilingMarkers();
  if ((markers.begin == nullptr) != (markers.end == nullptr)) {
    std::cerr << "only one of the two markers is set\n";
    return EXIT_FAILURE;
  }
  if ((markers.begin != nullptr) != (std::strcmp(expected, "yes") == 0)) {
    std::cerr << "the profiling sections are "
              << (markers.begin != nullptr ? "" : "not ")
              << "declared to roctx\n";
    return EXIT_FAILURE;
  }
  // sections and ranges around computations of the device, with or without
  // markers and without profiler
  constexpr auto n = 1000;
  auto norm = mfem::real_t{};
  auto ctx = mgis::Context{};
  ctx.enableProfiling(true);
  for (int i = 0; i != 3; ++i) {
    const auto range = mgis::ProfilingRange{
        ctx, [i] { return "ProfilingMarkersTest: step " + std::to_string(i); }};
    CatchTimeSection(ctx, "ProfilingMarkersTest::Step");
    auto v = mfem::Vector(n);
    v.UseDevice(true);
    {
      CatchTimeSection(ctx, "ProfilingMarkersTest::Fill");
      v = mfem::real_t{2};
    }
    {
      CatchTimeSection(ctx, "ProfilingMarkersTest::Norm");
      norm = v.Norml2();
    }
  }
  if (std::abs(norm * norm - mfem::real_t{4 * n}) > mfem::real_t{n} / 1000) {
    std::cerr << "invalid norm " << norm << '\n';
    return EXIT_FAILURE;
  }
  // the ranges are not added to the profiling tree
  const auto& root = ctx.getProfilingResultTree();
  if ((root.children.size() != 1) ||
      (root.children[0]->name != "ProfilingMarkersTest::Step") ||
      (root.children[0]->calls != 3) ||
      (root.children[0]->children.size() != 2) ||
      (root.children[0]->children[0]->name != "ProfilingMarkersTest::Fill") ||
      (root.children[0]->children[0]->calls != 3) ||
      (root.children[0]->children[1]->name != "ProfilingMarkersTest::Norm") ||
      (root.children[0]->children[1]->calls != 3)) {
    std::cerr << "unexpected profiling tree\n";
    return EXIT_FAILURE;
  }
  std::cout << "profiling sections " << (markers.begin != nullptr ? "" : "not ")
            << "declared to roctx\n";
  return EXIT_SUCCESS;
}  // end of main
