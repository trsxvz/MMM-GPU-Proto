/*!
 * \file   MaterialArraysTest.cxx
 * \brief  Test of the arrays of the materials in the memory of MFEM
 * (`MaterialArrays`): sizes, use without copy by MGIS, identity of the
 * deformation gradients, layouts, validity of the copies on the host and on
 * the device selected by the `--device` option at each access of
 * mfem-mgis on the host, arrays of the elements of a partial quadrature
 * space.
 *
 * With the `debug` device of MFEM, an access on the host to an array whose
 * copy on the host is not valid stops the program (protected memory): each
 * entry point is called after writing all the arrays on the device, so
 * that a missing synchronization makes the test fail.
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
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionUtilities.hxx"
#include "MFEMMGIS/PointWiseModel.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/TimeStep.hxx"

//! \brief number of failed checks
static int failures = 0;

static void check(const bool b, const std::string& msg) {
  if (!b) {
    std::cerr << "check failed: " << msg << '\n';
    ++failures;
  }
}  // end of check

using StateArray = mfem_mgis::MaterialArrays::StateArray;

//! \return the name of the device, used in the names of the output files
static std::string getDeviceName() {
  if (mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE)) {
    return "debug";
  }
  if (mfem::Device::Allows(mfem::Backend::CUDA_MASK)) {
    return "cuda";
  }
  return "cpu";
}  // end of getDeviceName

//! \brief arrays of a state
static constexpr StateArray state_arrays[] = {
    StateArray::GRADIENTS, StateArray::THERMODYNAMIC_FORCES,
    StateArray::INTERNAL_STATE_VARIABLES, StateArray::STORED_ENERGIES,
    StateArray::DISSIPATED_ENERGIES};

//! \brief states
static constexpr mfem_mgis::TimeStepStage states[] = {mfem_mgis::bts,
                                                      mfem_mgis::ets};

/*!
 * \return the value written at the first index of an array by
 * `writeAllOnTheDevice`
 * \param[in] k: value of the call
 * \param[in] a: array
 * \param[in] s: state
 */
static int getFirstValue(const int k,
                         const StateArray a,
                         const mfem_mgis::TimeStepStage s) {
  return k + 1000 * static_cast<int>(a) + (s == mfem_mgis::ets ? 500 : 0);
}  // end of getFirstValue

/*!
 * \brief writes the values `3 * i + k` in an array on the device
 * \param[in] d: pointer to the values on the device
 * \param[in] n: number of values
 * \param[in] k: value at the first index, different for each call, so that
 * a stale copy on the host can't hold the expected values
 */
static void writeOnTheDevice(mfem_mgis::real* const d,
                             const int n,
                             const int k) {
  mfem::forall(n, [=] MFEM_HOST_DEVICE(int i) { d[i] = 3 * i + k; });
}  // end of writeOnTheDevice

/*!
 * \brief writes all the arrays of a material on the device (values
 * `3 * i + getFirstValue(k, a, s)`, and `3 * i + k + 6000` for the tangent
 * operator blocks)
 * \param[in] a: arrays
 * \param[in] k: value of the call
 */
static void writeAllOnTheDevice(mfem_mgis::MaterialArrays& a, const int k) {
  for (const auto s : states) {
    for (const auto sa : state_arrays) {
      const auto n = a.getArray(sa, s).Size();
      if (n != 0) {
        writeOnTheDevice(a.readWrite(sa, s), n, getFirstValue(k, sa, s));
      }
    }
  }
  const auto n = a.getTangentOperatorBlocks().Size();
  if (n != 0) {
    writeOnTheDevice(a.readWriteTangentOperatorBlocks(), n, k + 6000);
  }
}  // end of writeAllOnTheDevice

/*!
 * \return true if the values of an array are `3 * i * stride + k`
 * \param[in] v: values
 * \param[in] stride: stride between two values
 * \param[in] n: number of values
 * \param[in] k: value at the first index
 */
static bool checkWrittenValues(const mfem_mgis::real* const v,
                               const mfem_mgis::size_type stride,
                               const mfem_mgis::size_type n,
                               const int k) {
  for (mfem_mgis::size_type i = 0; i != n; ++i) {
    if (v[i * stride] != static_cast<mfem_mgis::real>(3 * i * stride + k)) {
      return false;
    }
  }
  return true;
}  // end of checkWrittenValues

/*!
 * \brief checks that MGIS uses an array of the memory of MFEM without copy
 * \param[in] v: array of the memory of MFEM
 * \param[in] s: view used by MGIS
 * \param[in] n: name of the array
 */
static void checkAlias(const mfem::Vector& v,
                       const std::span<mfem_mgis::real> s,
                       const std::string& n) {
  check(static_cast<std::size_t>(v.Size()) == s.size(), n + ": size");
  if (v.Size() != 0) {
    check(v.HostRead() == s.data(), n + ": MGIS does not use the array");
  }
}  // end of checkAlias

/*!
 * \return true if two arrays have the same values on the host
 * \param[in] a: first array
 * \param[in] b: second array
 */
static bool haveSameValues(const mfem::Vector& a, const mfem::Vector& b) {
  if (a.Size() != b.Size()) {
    return false;
  }
  const auto* const ha = a.HostRead();
  const auto* const hb = b.HostRead();
  for (int i = 0; i != a.Size(); ++i) {
    if (ha[i] != hb[i]) {
      return false;
    }
  }
  return true;
}  // end of haveSameValues

static void testMaterial(
    mfem_mgis::Context& ctx,
    std::shared_ptr<const mfem_mgis::PartialQuadratureSpace> qspace,
    const char* const library,
    const char* const behaviour,
    const mgis::FieldLayout layout) {
  using namespace mfem_mgis;
  auto b = mfem_mgis::load(ctx, library, behaviour, Hypothesis::TRIDIMENSIONAL);
  if (isInvalid(b)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  auto m = Material(qspace, std::move(b), layout);
  auto& a = m.getArrays();
  const auto n = qspace->getNumberOfIntegrationPoints();
  const auto pbp = layout == mgis::FieldLayout::POINT_BY_POINT;
  const auto name = std::string{behaviour} +
                    (pbp ? " (point by point)" : " (component by component)");
  check((m.layout == layout) && (m.s0.layout == layout) &&
            (m.s1.layout == layout),
        name + ": layout not given to MGIS");
  // sizes and use without copy
  for (const auto st : states) {
    auto& s = st == bts ? m.s0 : m.s1;
    check(
        static_cast<size_type>(a.getArray(StateArray::GRADIENTS, st).Size()) ==
            n * s.gradients_stride,
        name + ": size of the gradients");
    checkAlias(a.getArray(StateArray::GRADIENTS, st), s.gradients,
               name + ": gradients");
    checkAlias(a.getArray(StateArray::THERMODYNAMIC_FORCES, st),
               s.thermodynamic_forces, name + ": thermodynamic forces");
    checkAlias(a.getArray(StateArray::INTERNAL_STATE_VARIABLES, st),
               s.internal_state_variables, name + ": internal state variables");
    checkAlias(a.getArray(StateArray::STORED_ENERGIES, st), s.stored_energies,
               name + ": stored energies");
    checkAlias(a.getArray(StateArray::DISSIPATED_ENERGIES, st),
               s.dissipated_energies, name + ": dissipated energies");
    check(m.b.computesStoredEnergy ==
              (a.getArray(StateArray::STORED_ENERGIES, st).Size() != 0),
          name + ": stored energies allocated if and only if computed");
    check(m.b.computesDissipatedEnergy ==
              (a.getArray(StateArray::DISSIPATED_ENERGIES, st).Size() != 0),
          name + ": dissipated energies allocated if and only if computed");
  }
  check(static_cast<size_type>(a.getTangentOperatorBlocks().Size()) ==
            n * m.K_stride,
        name + ": size of the tangent operator blocks");
  checkAlias(a.getTangentOperatorBlocks(), m.K,
             name + ": tangent operator blocks");
  // identity of the deformation gradients
  if ((m.b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
      (m.b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
    for (auto* s : {&m.s0, &m.s1}) {
      const auto F = mgis::behaviour::getGradientsView(*s);
      auto identity = true;
      for (size_type o = 0; o != n; ++o) {
        for (size_type c = 0; c != s->gradients_stride; ++c) {
          identity = identity && (F(o, c) == (c < 3 ? 1 : 0));
        }
      }
      check(identity,
            name + ": deformation gradient not initialized to identity");
    }
  }
  // the accessors based on a stride between points require the point by
  // point layout (the six accessors, mutable and read-only)
  {
    const auto& cm = static_cast<const Material&>(m);
    auto refused = [](const auto& o) { return isInvalid(o); };
    const auto g = m.b.gradients.front().name;
    const auto t = m.b.thermodynamic_forces.front().name;
    auto c = Context{};
    check(refused(getGradient(c, m, g, ets)) != pbp,
          name + ": layout and mutable accessor of the gradients");
    check(refused(getGradient(c, cm, g, ets)) != pbp,
          name + ": layout and read-only accessor of the gradients");
    check(refused(getThermodynamicForce(c, m, t, ets)) != pbp,
          name + ": layout and mutable accessor of the thermodynamic forces");
    check(refused(getThermodynamicForce(c, cm, t, ets)) != pbp,
          name + ": layout and read-only accessor of the thermodynamic forces");
    if (!m.b.isvs.empty()) {
      const auto i = m.b.isvs.front().name;
      check(refused(getInternalStateVariable(c, m, i, ets)) != pbp,
            name +
                ": layout and mutable accessor of the internal state "
                "variables");
      check(refused(getInternalStateVariable(c, cm, i, ets)) != pbp,
            name +
                ": layout and read-only accessor of the internal state "
                "variables");
    }
    check(!refused(getStoredEnergy(c, cm, ets)) == m.b.computesStoredEnergy,
          name + ": the energies are accessible in both layouts");
  }
  // validity of the copies: values written by a kernel on the device are
  // seen by MGIS on the host after `synchronizeOnHost`
  const auto& v = a.getArray(StateArray::THERMODYNAMIC_FORCES, ets);
  const auto size = v.Size();
  writeOnTheDevice(a.readWrite(StateArray::THERMODYNAMIC_FORCES, ets), size, 1);
  a.synchronizeOnHost();
  check(checkWrittenValues(m.s1.thermodynamic_forces.data(), 1, size, 1),
        name + ": values of the device not seen by MGIS");
  // the copy on the device is invalid: MGIS may modify the values on the
  // host, a kernel then reads them
  check(!v.GetMemory().DeviceIsValid() || !mfem::Device::IsEnabled(),
        name + ": the copy on the device is still valid");
  m.s1.thermodynamic_forces[0] = -1;
  auto* const r = v.Read();
  auto first = mfem::Vector(1);
  first.UseDevice(true);
  auto* const f = first.Write();
  mfem::forall(1, [=] MFEM_HOST_DEVICE(int) { f[0] = r[0]; });
  check(first.HostRead()[0] == -1,
        name + ": values of MGIS not seen on the device");
  // read-only accessor of the stored energies (one component: both
  // layouts)
  if (m.b.computesStoredEnergy) {
    writeAllOnTheDevice(a, 9);
    auto c = Context{};
    const auto oe = getStoredEnergy(c, static_cast<const Material&>(m), ets);
    check(!isInvalid(oe) &&
              checkWrittenValues(
                  oe->getValues().data(), 1, n,
                  getFirstValue(9, StateArray::STORED_ENERGIES, ets)),
          name + ": stored energies written on the device not seen");
  }
  // update and revert, in the layout of the material: the arrays written on
  // the device are copied from a state to the other one
  const auto copied = {StateArray::GRADIENTS, StateArray::THERMODYNAMIC_FORCES,
                       StateArray::INTERNAL_STATE_VARIABLES,
                       StateArray::STORED_ENERGIES};
  writeAllOnTheDevice(a, 7);
  m.updateMaterialStates();
  for (const auto sa : copied) {
    const auto& v0 = a.getArray(sa, bts);
    check(haveSameValues(v0, a.getArray(sa, ets)) &&
              checkWrittenValues(v0.HostRead(), 1, v0.Size(),
                                 getFirstValue(7, sa, ets)),
          name + ": update of array " + std::to_string(static_cast<int>(sa)));
  }
  // the tangent operator blocks are set to zero by MGIS on the host: their
  // copy on the device, written before, must not be seen
  {
    const auto& K = a.getTangentOperatorBlocks();
    auto zero = true;
    for (int i = 0; i != K.Size(); ++i) {
      zero = zero && (K.HostRead()[i] == 0);
    }
    check(zero, name + ": tangent operator blocks after the update");
  }
  writeAllOnTheDevice(a, 8);
  m.revertMaterialStates();
  for (const auto sa : copied) {
    const auto& v1 = a.getArray(sa, ets);
    check(haveSameValues(a.getArray(sa, bts), v1) &&
              checkWrittenValues(v1.HostRead(), 1, v1.Size(),
                                 getFirstValue(8, sa, bts)),
          name + ": revert of array " + std::to_string(static_cast<int>(sa)));
  }
}  // end of testMaterial

/*!
 * \brief a material without integration point (as on a process which owns
 * no element of the material)
 */
static void testEmptyMaterial(mfem_mgis::Context& ctx,
                              const char* const library) {
  using namespace mfem_mgis;
  auto b = mfem_mgis::load(ctx, library, "SaintVenantKirchhoffElasticityV2",
                           Hypothesis::TRIDIMENSIONAL);
  if (isInvalid(b)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  for (const auto l : {mgis::FieldLayout::POINT_BY_POINT,
                       mgis::FieldLayout::COMPONENT_BY_COMPONENT}) {
    try {
      auto a = MaterialArrays(*b, 0, l);
      auto m = mgis::behaviour::MaterialDataManager(*b, 0, a.getInitializer());
      a.synchronizeOnHost();
      a.synchronizeOnHostForReading();
      auto ops = mgis::HostFieldOperations{};
      mgis::behaviour::update(m, ops);
      mgis::behaviour::revert(m, ops);
      auto empty = a.getTangentOperatorBlocks().Size() == 0;
      for (const auto st : states) {
        for (const auto sa : state_arrays) {
          empty = empty && (a.getArray(sa, st).Size() == 0);
        }
      }
      check(empty, "material without point: empty arrays");
    } catch (std::exception& e) {
      check(false, std::string{"material without point: "} + e.what());
    }
  }
}  // end of testEmptyMaterial

/*!
 * \brief each entry point of mfem-mgis accessing the arrays on the host is
 * called after writing all the arrays on the device: the accessors, the
 * integration of the behaviours, the post-processings, revert and update
 * see the values written on the device (and, with the `debug` device of
 * MFEM, do not access an invalid copy on the host)
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
  for (const auto* const pp :
       {"MeanThermodynamicForces", "ComputeResultantForceOnBoundary"}) {
    auto params =
        Parameters{{"OutputFileName", std::string{pp} + "-MaterialArraysTest-" +
                                          getDeviceName() + ".txt"}};
    if (std::string_view{pp} == "ComputeResultantForceOnBoundary") {
      params.insert(throwing, {{"Boundary", 1}});
    }
    problem.addPostProcessing(ctx, pp, params) | or_die;
  }
  auto& m = problem.getMaterial(ctx, 1, 0) | or_die;
  auto& a = m.getArrays();
  const auto n = m.n;
  const auto stride = m.s1.internal_state_variables_stride;
  const auto first = [](const int k, const TimeStepStage s) {
    return getFirstValue(k, StateArray::INTERNAL_STATE_VARIABLES, s);
  };
  // first internal state variable (first component of the elastic strain)
  // at the end of the time step, through the accessor
  writeAllOnTheDevice(a, 1);
  const auto of = getInternalStateVariable(ctx, m, "ElasticStrain", ets);
  if (isInvalid(of)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  check(checkWrittenValues(of->getValues().data(), stride, n, first(1, ets)),
        "accessor: values of the device not seen");
  // read-only accessor: the copy on the device stays valid
  writeAllOnTheDevice(a, 2);
  const auto& cm = static_cast<const Material&>(m);
  const auto ocf = getInternalStateVariable(ctx, cm, "ElasticStrain", ets);
  if (isInvalid(ocf)) {
    std::cerr << ctx.getErrorMessage() << '\n';
    ++failures;
    return;
  }
  check(checkWrittenValues(ocf->getValues().data(), stride, n, first(2, ets)),
        "read-only accessor: values of the device not seen");
  check(!mfem::Device::IsEnabled() ||
            a.getArray(StateArray::INTERNAL_STATE_VARIABLES, ets)
                .GetMemory()
                .DeviceIsValid(),
        "read-only accessor: the copy on the device is invalidated");
  // integration of the behaviours (gradients computed from the unknowns,
  // all the arrays of the end of the time step written)
  writeAllOnTheDevice(a, 3);
  auto u = mfem::Vector(problem.getUnknowns(ets).Size());
  u = real{0};
  check(problem.integrate(
            u, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR, 1),
        "integration");
  // the integration writes the arrays of the end of the time step on the
  // host: with the debug device, a missing synchronization stops the program
  // (the values computed from the values written on the device are not
  // checked)
  // post-processings, which read the arrays on the host
  writeAllOnTheDevice(a, 4);
  check(problem.executePostProcessings(ctx, 0, 1), "post-processings");
  // mean values of the thermodynamic forces, called directly
  writeAllOnTheDevice(a, 12);
  check(!isInvalid(computeMeanThermodynamicForcesValues<false>(
            ctx, problem.getImplementation<false>())),
        "mean values of the thermodynamic forces");
  // the other accessors, mutable and read-only
  const auto gs = m.s1.gradients_stride;
  const auto ts = m.s1.thermodynamic_forces_stride;
  const auto first_value = [](const int k, const StateArray sa) {
    return getFirstValue(k, sa, ets);
  };
  writeAllOnTheDevice(a, 13);
  const auto og = getGradient(ctx, m, "Strain", ets);
  check(!isInvalid(og) &&
            checkWrittenValues(og->getValues().data(), gs, n,
                               first_value(13, StateArray::GRADIENTS)),
        "mutable accessor of the gradients");
  writeAllOnTheDevice(a, 14);
  const auto ocg = getGradient(ctx, cm, "Strain", ets);
  check(!isInvalid(ocg) &&
            checkWrittenValues(ocg->getValues().data(), gs, n,
                               first_value(14, StateArray::GRADIENTS)),
        "read-only accessor of the gradients");
  writeAllOnTheDevice(a, 15);
  const auto ot = getThermodynamicForce(ctx, m, "Stress", ets);
  check(
      !isInvalid(ot) &&
          checkWrittenValues(ot->getValues().data(), ts, n,
                             first_value(15, StateArray::THERMODYNAMIC_FORCES)),
      "mutable accessor of the thermodynamic forces");
  writeAllOnTheDevice(a, 16);
  const auto oct = getThermodynamicForce(ctx, cm, "Stress", ets);
  check(
      !isInvalid(oct) &&
          checkWrittenValues(oct->getValues().data(), ts, n,
                             first_value(16, StateArray::THERMODYNAMIC_FORCES)),
      "read-only accessor of the thermodynamic forces");
  // revert: the values of the beginning of the time step written on the
  // device are copied at the end of the time step
  writeAllOnTheDevice(a, 5);
  check(problem.revert(ctx), "revert");
  check(checkWrittenValues(m.s1.internal_state_variables.data(), 1, n * stride,
                           first(5, bts)),
        "revert: values of the device not seen");
  // update: the values at the end of the time step written on the device
  // are copied at the beginning of the time step
  writeAllOnTheDevice(a, 6);
  check(problem.update(ctx), "update");
  check(checkWrittenValues(m.s0.internal_state_variables.data(), 1, n * stride,
                           first(6, ets)),
        "update: values of the device not seen");
}  // end of testHostAccesses

/*!
 * \brief the thermodynamic forces written on the device are rotated to the
 * global frame (identity rotation) on the host, and the rotation is refused
 * in the component by component layout
 */
static void testRotation(
    mfem_mgis::Context& ctx,
    std::shared_ptr<const mfem_mgis::PartialQuadratureSpace> qspace,
    const char* const library) {
  using namespace mfem_mgis;
  for (const auto l : {mgis::FieldLayout::POINT_BY_POINT,
                       mgis::FieldLayout::COMPONENT_BY_COMPONENT}) {
    auto b = mfem_mgis::load(ctx, library, "OrthotropicElasticity",
                             Hypothesis::TRIDIMENSIONAL);
    if (isInvalid(b)) {
      std::cerr << ctx.getErrorMessage() << '\n';
      ++failures;
      return;
    }
    auto m = Material(qspace, std::move(b), l);
    m.setRotationMatrix(
        RotationMatrix3D{std::array<real, 9u>{1, 0, 0, 0, 1, 0, 0, 0, 1}});
    auto& a = m.getArrays();
    writeAllOnTheDevice(a, 11);
    auto f = PartialQuadratureFunction(qspace, 6);
    auto c = Context{};
    const auto ok = rotateThermodynamicsForcesToGlobalFrame(c, f, m, ets);
    if (l == mgis::FieldLayout::POINT_BY_POINT) {
      check(ok && checkWrittenValues(
                      f.getValues().data(), 6, m.n,
                      getFirstValue(11, StateArray::THERMODYNAMIC_FORCES, ets)),
            "rotation of the thermodynamic forces written on the device");
    } else {
      check(!ok, "rotation in the component by component layout");
    }
  }
}  // end of testRotation

/*!
 * \brief the integration of a point-wise model sees the values written on
 * the device
 */
static void testPointWiseModel(mfem_mgis::Context& ctx,
                               const char* const mesh,
                               const char* const library) {
  using namespace mfem_mgis;
  auto or_die = ctx.getFatalFailureHandler();
  auto fed = make_shared<FiniteElementDiscretization>(
                 ctx, Parameters{{"MeshFileName", mesh},
                                 {"FiniteElementFamily", "H1"},
                                 {"FiniteElementOrder", 1},
                                 {"UnknownsSize", 1},
                                 {"NumberOfUniformRefinements", 0},
                                 {"Parallel", false}}) |
             or_die;
  auto qspace = make_shared<const PartialQuadratureSpace>(
                    ctx, *fed, 5,
                    [](const mfem::FiniteElement& e,
                       const mfem::ElementTransformation&) noexcept
                        -> const mfem::IntegrationRule& {
                      return mfem::IntRules.Get(e.GetGeomType(), 0);
                    }) |
                or_die;
  auto model = make_shared<PointWiseModel>(
                   ctx, qspace,
                   Parameters{{"Library", library},
                              {"Model", "UO2Shrinkage_RAPHAEL2008"},
                              {"Hypothesis", "PlaneStrain"}}) |
               or_die;
  auto& m = model->getMaterial();
  for (auto* s : {&m.s0, &m.s1}) {
    mgis::behaviour::setExternalStateVariable(ctx, *s, "Temperature", 893.15) |
        or_die;
  }
  mgis::behaviour::setExternalStateVariable(ctx, m.s0, "BurnUp_AtPercent", 0) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m.s1, "BurnUp_AtPercent", 5) |
      or_die;
  // the shrinkage is computed incrementally from its value at the beginning
  // of the time step: this value is set to one on the host, then to zero on
  // the device, so that the integration gives the expected value only if it
  // sees the values of the device
  auto& a = m.getArrays();
  a.synchronizeOnHost();
  std::fill(m.s0.internal_state_variables.begin(),
            m.s0.internal_state_variables.end(), real{1});
  writeAllOnTheDevice(a, 17);
  {
    const auto n0 =
        a.getArray(StateArray::INTERNAL_STATE_VARIABLES, bts).Size();
    auto* const d0 = a.readWrite(StateArray::INTERNAL_STATE_VARIABLES, bts);
    mfem::forall(n0, [=] MFEM_HOST_DEVICE(int i) { d0[i] = 0; });
  }
  auto ps = construct<PhysicalSystem>(ctx, *fed) | or_die;
  ps.setModel(ctx, model) | or_die;
  const auto r = ps.computeNextState(ctx, {.begin = 0, .end = 1, .dt = 1});
  check(r.first.shallContinue(), "integration of the point-wise model");
  const auto s = getInternalStateVariable(ctx, m, "Shrinkage", ets) | or_die;
  const auto Ta = real{750};
  const auto Tm = real{893.15};
  const auto A = std::max(real{5.e-3}, real{-1.26e-2 + 1.8e-5 * Tm});
  const auto s_ref = A * (std::exp(-3.e-2 * 5 * (1 + std::exp(-Ta / Tm))) - 1);
  auto ok = true;
  for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
    ok = ok && (std::abs(s.getIntegrationPointValue(i) - s_ref) < 1e-8);
  }
  check(ok, "point-wise model: shrinkage");
}  // end of testPointWiseModel

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
  const char* mesh_2d = "bar.msh";
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&mesh, "-m", "--mesh", "mesh file");
  args.AddOption(&library, "-l", "--library", "behaviour library");
  args.AddOption(&multi_material_mesh, "-mm", "--multi-material-mesh",
                 "mesh with two materials");
  args.AddOption(&mesh_2d, "-m2", "--mesh-2d",
                 "two-dimensional mesh (point-wise model on material 5)");
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
  // the second Saint-Venant Kirchhoff behaviour computes the stored energy
  for (const auto* const b :
       {"IsotropicLinearHardeningPlasticity", "SaintVenantKirchhoffElasticity",
        "SaintVenantKirchhoffElasticityV2", "IsotropicLinearElasticity"}) {
    for (const auto l : {mgis::FieldLayout::POINT_BY_POINT,
                         mgis::FieldLayout::COMPONENT_BY_COMPONENT}) {
      testMaterial(ctx, qspace, library, b, l);
    }
  }
  testEmptyMaterial(ctx, library);
  testRotation(ctx, qspace, library);
  testPointWiseModel(ctx, mesh_2d, library);
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
