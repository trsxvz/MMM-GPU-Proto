/*!
 * \file   tests/ParaviewExportResultsTest.cxx
 * \brief  This test checks the files written by the export of the unknowns
 * of a non linear evolution problem for ParaView (`ParaviewExportResults`):
 * the coordinates of the nodes and the fields are written in double
 * precision by default, in single precision with the parameter
 * `SinglePrecision`; the fields added to the export (`addField`) are written
 * in the same files as the unknowns, on the same mesh.
 * \date   18/11/2026
 */

#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string_view>
#include "mfem/general/optparser.hpp"
#include "mfem/fem/fe_coll.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/fem/gridfunc.hpp"
#include "mfem/mesh/mesh.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/ParaviewExportResults.hxx"

static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    ++failures;
    std::cerr << "check failed: " << msg << '\n';
  }
}  // end of check

/*!
 * \return the content of the file written by this process for a cycle of
 * an export, empty if it can't be read
 * \param[in] d: directory of the export
 * \param[in] c: cycle
 */
static std::string readExportedFile(const std::string& d, const int c) {
  auto cycle = std::to_string(c);
  auto rank = std::to_string(mfem_mgis::getMPIRankForDefaultCommunicator());
  const auto f = d + "/Cycle" + std::string(6 - cycle.size(), '0') + cycle +
                 "/proc" + std::string(6 - rank.size(), '0') + rank + ".vtu";
  auto in = std::ifstream(f);
  auto s = std::ostringstream{};
  s << in.rdbuf();
  check(!s.str().empty(), "file '" + f + "' not read");
  return s.str();
}  // end of readExportedFile

//! \return the number of occurrences of a pattern in a string
static std::size_t count(const std::string& s, const std::string_view p) {
  auto n = std::size_t{};
  for (auto pos = s.find(p); pos != std::string::npos;
       pos = s.find(p, pos + p.size())) {
    ++n;
  }
  return n;
}  // end of count

/*!
 * \brief adds to a problem an export of the unknowns with another field, on
 * the space of the unknowns, set to 1 at the initial time and to 2 at the
 * end of the time steps
 * \return the time step stages given to the function of the field
 * \param[in, out] ctx: execution context
 * \param[in, out] problem: problem
 * \param[in] d: directory of the export
 */
template <bool parallel>
static std::shared_ptr<std::vector<mfem_mgis::TimeStepStage>> addFieldExport(
    mfem_mgis::Context& ctx,
    mfem_mgis::NonLinearEvolutionProblem& problem,
    const std::string& d) {
  using namespace mfem_mgis;
  auto& i = problem.getImplementation<parallel>();
  auto e = std::make_unique<ParaviewExportResults<parallel>>(
      ctx, i, Parameters{{"OutputFileName", d}});
  auto f = std::make_shared<GridFunction<parallel>>(&i.getFiniteElementSpace());
  auto stages = std::make_shared<std::vector<TimeStepStage>>();
  e->addField("f", [f, stages](Context&, const TimeStepStage s) {
    stages->push_back(s);
    *f = (s == bts) ? real{1} : real{2};
    return f.get();
  });
  i.addPostProcessing(ctx, std::move(e)) | ctx.getFatalFailureHandler();
  return stages;
}  // end of addFieldExport

/*!
 * \brief checks that a field defined on another mesh than the one of the
 * export is refused (sequential computations only)
 * \param[in, out] ctx: execution context
 * \param[in, out] problem: problem
 * \param[in] d: directory of the export
 */
static void checkFieldOnAnotherMesh(
    mfem_mgis::Context& ctx,
    mfem_mgis::NonLinearEvolutionProblem& problem,
    const std::string& d) {
  using namespace mfem_mgis;
  auto& i = problem.getImplementation<false>();
  auto e = ParaviewExportResults<false>(ctx, i, {{"OutputFileName", d}});
  auto mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON);
  auto collection = mfem::H1_FECollection(1, 3);
  auto space = mfem::FiniteElementSpace(&mesh, &collection, 3);
  auto g = mfem::GridFunction(&space);
  g = real{0};
  e.addField("g", [&g](Context&, const TimeStepStage) { return &g; });
  auto c = Context{};
  check(!e.execute(c, i, 0, 1), "a field on another mesh is not refused");
}  // end of checkFieldOnAnotherMesh

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "cube.mesh";
  int parallel = 0;
  auto args = mfem::OptionsParser(argc, argv);
  args.AddOption(&mesh, "-m", "--mesh", "mesh file");
  args.AddOption(&parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", mesh},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 3},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", bool(parallel)}}) |
                 or_die;
  const auto prefix = std::string{"ParaviewExportResultsTest-"} +
                      std::to_string(getMPISizeForDefaultCommunicator());
  // the coordinates of the nodes and the displacement are written in double
  // precision by default, in single precision with the parameter
  // `SinglePrecision`
  const auto d = prefix + "-double";
  const auto s = prefix + "-single";
  problem.addPostProcessing(ctx, "ParaviewExportResults",
                            {{"OutputFileName", d}}) |
      or_die;
  problem.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", s}, {"SinglePrecision", true}}) |
      or_die;
  // the displacement and another field, in the same files
  const auto field_export = prefix + "-field";
#ifdef MFEM_USE_MPI
  const auto stages = parallel
                          ? addFieldExport<true>(ctx, problem, field_export)
                          : addFieldExport<false>(ctx, problem, field_export);
#else  /* MFEM_USE_MPI */
  const auto stages = addFieldExport<false>(ctx, problem, field_export);
#endif /* MFEM_USE_MPI */
  problem.executeInitialPostProcessings(ctx, 0) | or_die;
  problem.executePostProcessings(ctx, 0, 1) | or_die;
  check((stages->size() == 2) && (stages->front() == bts) &&
            (stages->back() == ets),
        "the function of the field is not called once per export, with the "
        "stage of the exported unknowns");
  for (const auto c : {0, 1}) {
    const auto fd = readExportedFile(d, c);
    check((count(fd, "type=\"Float64\"") == 2) &&
              (count(fd, "type=\"Float32\"") == 0),
          "cycle " + std::to_string(c) +
              ": coordinates and displacement not written in double "
              "precision");
    const auto fs = readExportedFile(s, c);
    check((count(fs, "type=\"Float32\"") == 2) &&
              (count(fs, "type=\"Float64\"") == 0),
          "cycle " + std::to_string(c) +
              ": coordinates and displacement not written in single "
              "precision");
    const auto ff = readExportedFile(field_export, c);
    check((count(ff, "<Points>") == 1) && (count(ff, "Name=\"u\"") == 1) &&
              (count(ff, "Name=\"f\"") == 1),
          "cycle " + std::to_string(c) +
              ": the displacement and the other field are not written in "
              "the same file, on one mesh");
  }
  if (!parallel) {
    checkFieldOnAnotherMesh(ctx, problem, prefix + "-other-mesh");
  }
  std::cout << (failures == 0 ? "success" : "failure") << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
