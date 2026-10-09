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

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  const char* mesh = "cube.mesh";
  const char* library = "src/libBehaviourTest.so";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "mesh file");
  args.AddOption(&library, "-l", "--library", "behaviour library");
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
  std::cout << (failures == 0 ? "success" : "failure") << " (device "
            << (mfem::Device::IsEnabled() ? "enabled" : "disabled") << ")\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
