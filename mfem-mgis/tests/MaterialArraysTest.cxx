/*!
 * \file   MaterialArraysTest.cxx
 * \brief  Test of the arrays of the materials in the memory of MFEM
 * (`MaterialArrays`): sizes, use without copy by MGIS, identity of the
 * deformation gradients, validity of the copies on the host and on the
 * device selected by the `--device` option
 * \date   09/10/2026
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include "mfem/general/forall.hpp"
#include "mfem/general/optparser.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

/*!
 * \brief checks that MGIS uses an array of the memory of MFEM without copy
 * \param[in] v: array of the memory of MFEM
 * \param[in] s: view used by MGIS
 * \param[in] n: name of the array
 */
static void checkAlias(const mfem::Vector& v,
                       const std::span<mfem_mgis::real> s,
                       const std::string& n) {
  if (v.Size() == 0) {
    return;
  }
  check(static_cast<std::size_t>(v.Size()) == s.size(), n + ": size");
  check(v.HostRead() == s.data(), n + ": MGIS does not use the array");
}  // end of checkAlias

static void testMaterial(
    mfem_mgis::Context& ctx,
    std::shared_ptr<const mfem_mgis::PartialQuadratureSpace> qspace,
    const char* const library,
    const char* const behaviour) {
  using namespace mfem_mgis;
  auto b = mfem_mgis::load(ctx, library, behaviour, Hypothesis::TRIDIMENSIONAL);
  if (isInvalid(b)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  auto m = Material(qspace, std::move(b));
  auto& a = m.getArrays();
  const auto n = qspace->getNumberOfIntegrationPoints();
  const auto name = std::string{behaviour};
  // sizes and use without copy
  for (std::size_t i = 0; i != 2; ++i) {
    auto& s = i == 0 ? m.s0 : m.s1;
    check(
        static_cast<size_type>(a.gradients[i].Size()) == n * s.gradients_stride,
        name + ": size of the gradients");
    checkAlias(a.gradients[i], s.gradients, name + ": gradients");
    checkAlias(a.thermodynamic_forces[i], s.thermodynamic_forces,
               name + ": thermodynamic forces");
    checkAlias(a.internal_state_variables[i], s.internal_state_variables,
               name + ": internal state variables");
    checkAlias(a.stored_energies[i], s.stored_energies,
               name + ": stored energies");
    checkAlias(a.dissipated_energies[i], s.dissipated_energies,
               name + ": dissipated energies");
  }
  check(static_cast<size_type>(a.tangent_operator_blocks.Size()) ==
            n * m.K_stride,
        name + ": size of the tangent operator blocks");
  checkAlias(a.tangent_operator_blocks, m.K,
             name + ": tangent operator blocks");
  // identity of the deformation gradients
  if ((m.b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
      (m.b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
    for (const auto* s : {&m.s0, &m.s1}) {
      for (size_type o = 0; o != n; ++o) {
        const auto* const F = s->gradients.data() + o * s->gradients_stride;
        for (size_type c = 0; c != s->gradients_stride; ++c) {
          check(F[c] == (c < 3 ? 1 : 0),
                name + ": deformation gradient not initialized to identity");
        }
      }
    }
  }
  // validity of the copies: values written by a kernel on the device are
  // seen by MGIS on the host after `synchronizeOnHost`
  auto& v = a.thermodynamic_forces[1];
  const auto size = v.Size();
  auto* const d = v.ReadWrite();
  mfem::forall(size, [=] MFEM_HOST_DEVICE(int i) { d[i] = 3 * i + 1; });
  a.synchronizeOnHost();
  const auto& sig = m.s1.thermodynamic_forces;
  for (int i = 0; i != size; ++i) {
    if (sig[i] != 3 * i + 1) {
      check(false, name + ": values of the device not seen by MGIS");
      break;
    }
  }
  // the copy on the device is invalid: MGIS may modify the values on the
  // host, a kernel then reads them
  check(!v.GetMemory().DeviceIsValid() || !mfem::Device::IsEnabled(),
        name + ": the copy on the device is still valid");
  sig[0] = -1;
  auto* const r = v.Read();
  auto first = mfem::Vector(1);
  first.UseDevice(true);
  auto* const f = first.Write();
  mfem::forall(1, [=] MFEM_HOST_DEVICE(int) { f[0] = r[0]; });
  check(first.HostRead()[0] == -1,
        name + ": values of MGIS not seen on the device");
}  // end of testMaterial

/*!
 * \brief writes the values `3 * i + k` in an array on the device
 * \param[in] v: array
 * \param[in] k: value at the first index, different for each call, so that
 * a stale copy on the host can't hold the expected values
 */
static void writeOnTheDevice(mfem::Vector& v, const int k) {
  auto* const d = v.ReadWrite();
  mfem::forall(v.Size(), [=] MFEM_HOST_DEVICE(int i) { d[i] = 3 * i + k; });
}  // end of writeOnTheDevice

/*!
 * \brief checks the values written by `writeOnTheDevice`
 * \param[in] v: values
 * \param[in] stride: stride between two values
 * \param[in] n: number of values
 * \param[in] k: value at the first index
 * \param[in] msg: message
 */
static void checkWrittenValues(const mfem_mgis::real* const v,
                               const mfem_mgis::size_type stride,
                               const mfem_mgis::size_type n,
                               const int k,
                               const std::string& msg) {
  for (mfem_mgis::size_type i = 0; i != n; ++i) {
    if (v[i * stride] != static_cast<mfem_mgis::real>(3 * i * stride + k)) {
      check(false, msg);
      return;
    }
  }
}  // end of checkWrittenValues

/*!
 * \brief the accessors of the materials and the update of a problem see
 * the values written on the device
 */
static void testHostAccesses(mfem_mgis::Context& ctx,
                             const char* const mesh,
                             const char* const library) {
  using namespace mfem_mgis;
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", mesh},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 3},
                                     {"NumberOfUniformRefinements", 1},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", false}}) |
                 or_die;
  problem.addBehaviourIntegrator(ctx, "Mechanics", 1, library,
                                 "IsotropicLinearHardeningPlasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, 1, 0) | or_die;
  auto& isvs = m.getArrays().internal_state_variables[1];
  const auto n = m.n;
  const auto stride = m.s1.internal_state_variables_stride;
  // first internal state variable (first component of the elastic strain)
  // at the end of the time step, through the accessor
  writeOnTheDevice(isvs, 1);
  const auto of = getInternalStateVariable(ctx, m, "ElasticStrain", ets);
  if (isInvalid(of)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  checkWrittenValues(of->getValues().data(), stride, n, 1,
                     "accessor: values of the device not seen");
  // read-only accessor: the copy on the device stays valid
  writeOnTheDevice(isvs, 2);
  const auto& cm = static_cast<const Material&>(m);
  const auto ocf = getInternalStateVariable(ctx, cm, "ElasticStrain", ets);
  if (isInvalid(ocf)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  checkWrittenValues(ocf->getValues().data(), stride, n, 2,
                     "read-only accessor: values of the device not seen");
  check(!mfem::Device::IsEnabled() || isvs.GetMemory().DeviceIsValid(),
        "read-only accessor: the copy on the device is invalidated");
  // update: the values at the end of the time step written on the device
  // are copied at the beginning of the time step
  writeOnTheDevice(isvs, 3);
  if (!problem.update(ctx)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  checkWrittenValues(m.s0.internal_state_variables.data(), 1, n * stride, 3,
                     "update: values of the device not seen");
}  // end of testHostAccesses

/*!
 * \brief the arrays of the elements and of the offsets of a partial
 * quadrature space agree with its hash tables, and are readable on the
 * device
 * \param[in] s: partial quadrature space
 */
static void testQuadratureSpaceArrays(
    const mfem_mgis::PartialQuadratureSpace& s) {
  using namespace mfem_mgis;
  const auto& elements = s.getElements();
  const auto& offsets = s.getElementOffsets();
  check(static_cast<size_type>(elements.Size()) == s.getNumberOfElements(),
        "quadrature space: number of elements");
  check(offsets.Size() == elements.Size() + 1,
        "quadrature space: number of offsets");
  check(static_cast<size_type>(offsets.HostRead()[elements.Size()]) ==
            s.getNumberOfIntegrationPoints(),
        "quadrature space: last offset");
  const auto* const he = elements.HostRead();
  const auto* const ho = offsets.HostRead();
  for (int i = 0; i != elements.Size(); ++i) {
    const auto e = static_cast<size_type>(he[i]);
    check((i == 0) || (he[i - 1] < he[i]),
          "quadrature space: elements not in increasing order");
    check(s.getOffset(e) == static_cast<size_type>(ho[i]),
          "quadrature space: offset of element " + std::to_string(e));
    check(s.getNumberOfQuadraturePoints(e) ==
              static_cast<size_type>(ho[i + 1] - ho[i]),
          "quadrature space: number of points of element " + std::to_string(e));
  }
  // number of points of each element, computed on the device
  auto npoints = mfem::Array<int>(elements.Size());
  auto* const np = npoints.Write();
  const auto* const o = offsets.Read();
  mfem::forall(elements.Size(),
               [=] MFEM_HOST_DEVICE(int i) { np[i] = o[i + 1] - o[i]; });
  const auto* const hnp = npoints.HostRead();
  auto total = size_type{};
  for (int i = 0; i != npoints.Size(); ++i) {
    total += static_cast<size_type>(hnp[i]);
  }
  check(total == s.getNumberOfIntegrationPoints(),
        "quadrature space: offsets read on the device");
}  // end of testQuadratureSpaceArrays

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "cube.mesh";
  const char* library = "src/libBehaviourTest.so";
  const char* multi_material_mesh = "cube_2mat_per.mesh";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "mesh file");
  args.AddOption(&library, "-l", "--library", "behaviour library");
  args.AddOption(&multi_material_mesh, "-mm", "--multi-material-mesh",
                 "mesh with two materials");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto ctx = Context{};
  auto fed = std::make_shared<FiniteElementDiscretization>(
      ctx, Parameters{{"MeshFileName", mesh},
                      {"FiniteElementFamily", "H1"},
                      {"FiniteElementOrder", 1},
                      {"UnknownsSize", 3},
                      {"NumberOfUniformRefinements", 2},
                      {"Parallel", false}});
  auto qspace = std::make_shared<const PartialQuadratureSpace>(
      *fed, 1,
      [](const mfem::FiniteElement& e,
         const mfem::ElementTransformation& tr) noexcept
          -> const mfem::IntegrationRule& {
        return mfem::IntRules.Get(e.GetGeomType(), 2 * tr.OrderGrad(&e));
      });
  for (const auto* const b :
       {"IsotropicLinearHardeningPlasticity", "SaintVenantKirchhoffElasticity",
        "IsotropicLinearElasticity"}) {
    testMaterial(ctx, qspace, library, b);
  }
  testHostAccesses(ctx, mesh, library);
  testQuadratureSpaceArrays(*qspace);
  {
    // a material covering a part of the mesh
    auto fed2 = FiniteElementDiscretization(
        ctx, Parameters{{"MeshFileName", multi_material_mesh},
                        {"FiniteElementFamily", "H1"},
                        {"FiniteElementOrder", 1},
                        {"UnknownsSize", 3},
                        {"NumberOfUniformRefinements", 0},
                        {"Parallel", false}});
    for (const size_type m : {1, 2}) {
      const auto s = PartialQuadratureSpace(
          fed2, m,
          [](const mfem::FiniteElement& e,
             const mfem::ElementTransformation& tr) noexcept
              -> const mfem::IntegrationRule& {
            return mfem::IntRules.Get(e.GetGeomType(), 2 * tr.OrderGrad(&e));
          });
      check((s.getNumberOfElements() != 0) &&
                (s.getNumberOfElements() <
                 static_cast<size_type>(fed2.getMesh<false>().GetNE())),
            "the material does not cover a part of the mesh");
      testQuadratureSpaceArrays(s);
    }
  }
  std::cout << (failures == 0 ? "success" : "failure") << " (device "
            << (mfem::Device::IsEnabled() ? "enabled" : "disabled") << ")\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
