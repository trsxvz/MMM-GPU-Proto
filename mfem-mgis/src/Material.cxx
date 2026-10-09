/*!
 * \file   src/Material.cxx
 * \brief  This file implements the `Material` class
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#include <limits>
#include <string>
#include <variant>
#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/RotationMatrix.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  //! \brief arrays of a state of a material
  using StateArray = MaterialArrays::StateArray;

  [[noreturn]] static std::array<real, 9u> raiseInvalidGetRotationMatrixCall(
      const RotationMatrix2D &, const RotationMatrix3D &, const size_type) {
    raise(
        "Material::getRotationMatrix: "
        "invalid call for isotropic behaviours");
  }  // end of raiseInvalidGetRotationMatrixCall

  [[noreturn]] static std::array<real, 9u> raiseUnsetRotationMatrix(
      const RotationMatrix2D &, const RotationMatrix3D &, const size_type) {
    raise(
        "Material::getRotationMatrix: "
        "unset rotation matrix");
  }  // end of raiseInvalidGetRotationMatrixCall

  [[nodiscard]] static const Behaviour &Material_dereferenceBehaviour(
      std::unique_ptr<const Behaviour> &b_ptr) {
    if (isInvalid(b_ptr)) {
      raise("invalid behaviour");
    }
    return *b_ptr;
  }

  [[nodiscard]] static const PartialQuadratureSpace &
  Material_dereferenceQuadratureSpace(
      std::shared_ptr<const PartialQuadratureSpace> s) {
    if (isInvalid(s)) {
      raise("invalid partial quadrature space");
    }
    return *s;
  }

  Material::Material(std::shared_ptr<const PartialQuadratureSpace> s,
                     std::unique_ptr<const Behaviour> b_ptr,
                     const mgis::FieldLayout l)
      : MaterialArrays(Material_dereferenceBehaviour(b_ptr),
                       Material_dereferenceQuadratureSpace(s)
                           .getNumberOfIntegrationPoints(),
                       l),
        MaterialDataManager(Material_dereferenceBehaviour(b_ptr),
                            Material_dereferenceQuadratureSpace(s)
                                .getNumberOfIntegrationPoints(),
                            this->getInitializer()),
        quadrature_space(s),
        macroscopic_gradients(this->s1.gradients_stride, real(0)),
        get_rotation_fct_ptr(this->b.symmetry ==
                                     mgis::behaviour::Behaviour::ORTHOTROPIC
                                 ? &raiseUnsetRotationMatrix
                                 : &raiseInvalidGetRotationMatrixCall),
        behaviour_ptr(std::move(b_ptr)) {
    this->allocateArrayOfTangentOperatorBlocks();
  }  // end of Material::Material

  void Material::updateMaterialStates() {
    // MGIS uses the copies of the arrays on the host
    this->synchronizeOnHost();
    auto ops = mgis::HostFieldOperations{};
    mgis::behaviour::update(*this, ops);
  }  // end of updateMaterialStates

  void Material::revertMaterialStates() {
    this->synchronizeOnHost();
    auto ops = mgis::HostFieldOperations{};
    mgis::behaviour::revert(*this, ops);
  }  // end of revertMaterialStates

  MaterialArrays &Material::getArrays() {
    return *this;
  }  // end of getArrays

  const MaterialArrays &Material::getArrays() const {
    return *this;
  }  // end of getArrays

  const PartialQuadratureSpace &Material::getPartialQuadratureSpace() const {
    return *(this->quadrature_space);
  }  // end of getPartialQuadratureSpace

  std::shared_ptr<const PartialQuadratureSpace>
  Material::getPartialQuadratureSpacePointer() const {
    return this->quadrature_space;
  }  // end of getPartialQuadratureSpacePointer

  void Material::setMacroscopicGradients(std::span<const real> g) {
    if (g.size() != this->s1.gradients_stride) {
      raise(
          "Material::setMacroscopicGradients: "
          "invalid number of components of the gradients");
    }
    std::copy(g.begin(), g.end(), this->macroscopic_gradients.begin());
  }  // end of setMacroscopicGradients

  static void checkBehaviourSymmetry(const Behaviour &b) {
    if (b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour is not orthotropic)");
    }
  }  // end of checkBehaviourSymmetry

  void Material::setRotationMatrix(const RotationMatrix2D &r) {
    checkBehaviourSymmetry(this->b);
    if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) != 2u) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour' hypothesis is not 2D)");
    }
    if (std::holds_alternative<std::array<real, 9u>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &rm, const RotationMatrix3D &,
              const size_type) { return std::get<std::array<real, 9u>>(rm); };
    } else if (std::holds_alternative<
                   std::shared_ptr<PartialQuadratureFunction>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &rm, const RotationMatrix3D &,
              const size_type o) {
            const auto &f =
                *(std::get<std::shared_ptr<PartialQuadratureFunction>>(rm));
            return mgis::behaviour::buildRotationMatrix(
                f.getIntegrationPointValues<2>(o));
          };
    } else {
      raise("Material::setRotationMatrix: unimplemented case yet");
    }
    this->r2D = r;
  }  // end of setRotationMatrix

  void Material::setRotationMatrix(const RotationMatrix3D &r) {
    checkBehaviourSymmetry(this->b);
    if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) != 3u) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour hypothesis is not 3D)");
    }
    if (std::holds_alternative<std::array<real, 9u>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &,  //
              const RotationMatrix3D &rm,
              const size_type) { return std::get<std::array<real, 9u>>(rm); };
    } else if (std::holds_alternative<std::array<MaterialAxis3D, 2u>>(r)) {
      const auto &[v1, v2] = std::get<std::array<MaterialAxis3D, 2u>>(r);
      if ((std::holds_alternative<std::array<real, 3u>>(v1)) &&
          (std::holds_alternative<std::array<real, 3u>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          return mgis::behaviour::buildRotationMatrix(
              std::get<std::array<real, 3u>>(a1),
              std::get<std::array<real, 3u>>(a2));
        };
      } else if ((std::holds_alternative<std::array<real, 3u>>(v1)) &&
                 (std::holds_alternative<std::array<real, 3u>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          return mgis::behaviour::buildRotationMatrix(
              std::get<std::array<real, 3u>>(a1),
              std::get<std::array<real, 3u>>(a2));
        };
      } else if ((std::holds_alternative<std::array<real, 3u>>(v1)) &&
                 (std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type o) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          const auto &f2 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a2));
          return mgis::behaviour::buildRotationMatrix(
              std::get<std::array<real, 3u>>(a1),
              f2.getIntegrationPointValues<3>(o));
        };
      } else if ((std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v1)) &&
                 (std::holds_alternative<std::array<real, 3u>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type o) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          const auto &f1 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a1));
          return mgis::behaviour::buildRotationMatrix(
              f1.getIntegrationPointValues<3>(o),
              std::get<std::array<real, 3u>>(a2));
        };
      } else if ((std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v1)) &&
                 (std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type o) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          const auto &f1 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a1));
          const auto &f2 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a2));
          return mgis::behaviour::buildRotationMatrix(
              f1.getIntegrationPointValues<3>(o),
              f2.getIntegrationPointValues<3>(o));
        };
      } else {
        raise("Material::setRotationMatrix: unimplemented case yet");
      }
    } else {
      raise("Material::setRotationMatrix: unimplemented case yet");
    }
    this->r3D = r;
  }  // end of setRotationMatrix

  /*!
   * \brief resizes an array of the memory of MFEM
   * \param[out] v: array
   * \param[in] s: size
   * \param[in] m: calling method
   */
  static void resizeArray(mfem::Vector &v,
                          const std::size_t s,
                          const char *const m) {
    if (s > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      raise(std::string{m} + ": array too large for MFEM");
    }
    v.SetSize(static_cast<int>(s));
  }  // end of resizeArray

  mgis::FieldView<const real> Material::computeRotationMatrices(
      mfem::Vector &r, const bool t) const {
    if (this->b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      raise(
          "Material::computeRotationMatrices: "
          "invalid call (behaviour is not orthotropic)");
    }
    // an unset rotation matrix is not uniform: the first call to
    // `getRotationMatrixAtIntegrationPoint` reports it
    const auto uniform = [this] {
      if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) == 2u) {
        return std::holds_alternative<std::array<real, 9u>>(this->r2D);
      }
      if (std::holds_alternative<std::array<real, 9u>>(this->r3D)) {
        return true;
      }
      const auto *const axes =
          std::get_if<std::array<MaterialAxis3D, 2u>>(&(this->r3D));
      return (axes != nullptr) &&
             (std::holds_alternative<std::array<real, 3u>>((*axes)[0])) &&
             (std::holds_alternative<std::array<real, 3u>>((*axes)[1]));
    }();
    const auto np = static_cast<size_type>(this->n);
    const auto nm = uniform ? size_type{1} : np;
    const auto size = 9 * static_cast<std::size_t>(nm);
    resizeArray(r, size, "Material::computeRotationMatrices");
    if (size == 0) {
      return mgis::makePointByPointFieldView(std::span<const real>{}, 0, 9, 9);
    }
    auto *const values = r.HostWrite();
    for (size_type i = 0; i != nm; ++i) {
      const auto m = this->getRotationMatrixAtIntegrationPoint(i);
      auto *const v = values + 9 * static_cast<std::size_t>(i);
      for (std::size_t j = 0; j != 3; ++j) {
        for (std::size_t k = 0; k != 3; ++k) {
          v[3 * j + k] = t ? m[3 * k + j] : m[3 * j + k];
        }
      }
    }
    const auto matrices = std::span<const real>(values, size);
    if (uniform) {
      return mgis::makeUniformFieldView(matrices, np);
    }
    return mgis::makePointByPointFieldView(matrices, np, 9, 9);
  }  // end of computeRotationMatrices

  void Material::computeValuesInGlobalFrame(mfem::Vector &v,
                                            const StateArray a,
                                            const StateSelection s) const {
    const auto gradients = a == StateArray::GRADIENTS;
    const auto *const method =
        gradients ? "Material::computeGradientsInGlobalFrame"
                  : "Material::computeThermodynamicForcesInGlobalFrame";
    if (this->b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      raise(std::string{method} +
            ": invalid call (behaviour is not orthotropic)");
    }
    const auto &sm = getStateManager(*this, s);
    const auto np = static_cast<size_type>(this->n);
    const auto nc =
        gradients ? sm.gradients_stride : sm.thermodynamic_forces_stride;
    const auto size = static_cast<std::size_t>(np) * nc;
    resizeArray(v, size, method);
    if (size == 0) {
      return;
    }
    this->synchronizeOnHostForReading(a, s);
    auto r = mfem::Vector{};
    // the gradients are rotated by the transposed matrices, from the
    // material frame to the global frame
    const auto matrices = this->computeRotationMatrices(r, gradients);
    const auto d = mgis::makePointByPointFieldView(
        std::span<real>(v.HostWrite(), size), np, nc, nc);
    if (gradients) {
      mgis::behaviour::rotateGradients(
          d, this->b, mgis::behaviour::getGradientsView(sm), matrices);
    } else {
      mgis::behaviour::rotateThermodynamicForces(
          d, this->b, mgis::behaviour::getThermodynamicForcesView(sm),
          matrices);
    }
  }  // end of computeValuesInGlobalFrame

  void Material::computeGradientsInGlobalFrame(mfem::Vector &v,
                                               const StateSelection s) const {
    this->computeValuesInGlobalFrame(v, StateArray::GRADIENTS, s);
  }  // end of computeGradientsInGlobalFrame

  void Material::computeThermodynamicForcesInGlobalFrame(
      mfem::Vector &v, const StateSelection s) const {
    this->computeValuesInGlobalFrame(v, StateArray::THERMODYNAMIC_FORCES, s);
  }  // end of computeThermodynamicForcesInGlobalFrame

  Material::~Material() = default;

  static std::optional<PartialQuadratureFunction>
  buildPartialQuadratureFunction(
      Context &ctx,
      std::shared_ptr<const PartialQuadratureSpace> qs,
      std::span<mgis::real> values,
      const std::vector<mgis::behaviour::Variable> &variables,
      const std::string_view n,
      const Hypothesis h) noexcept {
    using namespace mgis::behaviour;
    const auto ostride = getArraySize(ctx, variables, h);
    const auto ov = getVariable(ctx, variables, n);
    const auto oo = getVariableOffset(ctx, variables, n, h);
    if (isInvalid(ostride) || isInvalid(oo) || isInvalid(ov)) {
      return {};
    }
    const auto os = getVariableSize(ctx, *ov, h);
    if (isInvalid(os)) {
      return {};
    }
    return PartialQuadratureFunction::borrow(
        ctx, qs, values,
        {.data_begin = static_cast<size_type>(*oo),
         .data_size = static_cast<size_type>(*os),
         .data_stride = static_cast<size_type>(*ostride)});
  }  // end of buildPartialQuadratureFunction

  //   static PartialQuadratureFunction buildPartialQuadratureFunction(
  //       std::shared_ptr<const PartialQuadratureSpace> qs,
  //       std::span<mgis::real> values,
  //       const std::vector<mgis::behaviour::Variable> &variables,
  //       const std::string_view n,
  //       const Hypothesis h) {
  //     const auto o = getVariableOffset(variables, n, h);
  //     const auto s =
  //         getVariableSize(mgis::behaviour::getVariable(variables, n), h);
  //     return PartialQuadratureFunction(qs, StorageMode::EXTERNAL_STORAGE,
  //     values,
  //                                      o, s);
  //   }  // end of buildPartialQuadratureFunction

  static std::optional<ImmutablePartialQuadratureFunctionView>
  buildImmutablePartialQuadratureFunctionView(
      Context &ctx,
      std::shared_ptr<const PartialQuadratureSpace> qs,
      std::span<const mgis::real> values,
      const std::vector<mgis::behaviour::Variable> &variables,
      const std::string_view n,
      const Hypothesis h) {
    using namespace mgis::behaviour;
    const auto ostride = getArraySize(ctx, variables, h);
    const auto ov = getVariable(ctx, variables, n);
    const auto oo = getVariableOffset(ctx, variables, n, h);
    if (isInvalid(ostride) || isInvalid(oo) || isInvalid(ov)) {
      return {};
    }
    const auto os = getVariableSize(ctx, *ov, h);
    if (isInvalid(os)) {
      return {};
    }
    return ImmutablePartialQuadratureFunctionView(
        qs, values,
        {.data_begin = static_cast<size_type>(*oo),
         .data_size = static_cast<size_type>(*os),
         .data_stride = static_cast<size_type>(*ostride)});
  }  // end of buildImmutablePartialQuadratureFunctionView

  bool checkPointByPointLayout(AbstractErrorHandler &eh,
                               const Material &m) noexcept {
    if (m.layout != mgis::FieldLayout::POINT_BY_POINT) {
      return eh.registerErrorMessage(
          "the accessors and the evaluators of the material require the "
          "point by point layout of its arrays");
    }
    return true;
  }  // end of checkPointByPointLayout

  /*!
   * \brief synchronize an array of a material on the host, reporting the
   * exceptions in the context
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] a: array
   * \param[in] st: state
   */
  static bool synchronizeOnHost(Context &ctx,
                                Material &m,
                                const StateArray a,
                                const Material::StateSelection st) noexcept {
    try {
      m.synchronizeOnHost(a, st);
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of synchronizeOnHost

  /*!
   * \brief synchronize an array of a material on the host for reading,
   * reporting the exceptions in the context
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] a: array
   * \param[in] st: state
   */
  static bool synchronizeOnHostForReading(
      Context &ctx,
      const Material &m,
      const StateArray a,
      const Material::StateSelection st) noexcept {
    try {
      m.synchronizeOnHostForReading(a, st);
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of synchronizeOnHostForReading

  std::optional<PartialQuadratureFunction> getGradient(
      Context &ctx,
      Material &m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!checkPointByPointLayout(ctx, m)) {
      return {};
    }
    if (!synchronizeOnHost(ctx, m, StateArray::GRADIENTS, s)) {
      return {};
    }
    return buildPartialQuadratureFunction(
        ctx, m.getPartialQuadratureSpacePointer(),
        getStateManager(m, s).gradients, m.b.gradients, n, m.b.hypothesis);
  }  // end of getGradient

  std::optional<ImmutablePartialQuadratureFunctionView> getGradient(
      Context &ctx,
      const Material &m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!checkPointByPointLayout(ctx, m)) {
      return {};
    }
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::GRADIENTS, s)) {
      return {};
    }
    return buildImmutablePartialQuadratureFunctionView(
        ctx, m.getPartialQuadratureSpacePointer(),
        getStateManager(m, s).gradients, m.b.gradients, n, m.b.hypothesis);
  }  // end of getGradient

  std::optional<PartialQuadratureFunction> getThermodynamicForce(
      Context &ctx,
      Material &m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!checkPointByPointLayout(ctx, m)) {
      return {};
    }
    if (!synchronizeOnHost(ctx, m, StateArray::THERMODYNAMIC_FORCES, s)) {
      return {};
    }
    return buildPartialQuadratureFunction(
        ctx, m.getPartialQuadratureSpacePointer(),
        getStateManager(m, s).thermodynamic_forces, m.b.thermodynamic_forces, n,
        m.b.hypothesis);
  }  // end of getThermodynamicForce

  std::optional<ImmutablePartialQuadratureFunctionView> getThermodynamicForce(
      Context &ctx,
      const Material &m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!checkPointByPointLayout(ctx, m)) {
      return {};
    }
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::THERMODYNAMIC_FORCES, s)) {
      return {};
    }
    return buildImmutablePartialQuadratureFunctionView(
        ctx, m.getPartialQuadratureSpacePointer(),
        getStateManager(m, s).thermodynamic_forces, m.b.thermodynamic_forces, n,
        m.b.hypothesis);
  }  // end of getThermodynamicForce

  std::optional<PartialQuadratureFunction> getInternalStateVariable(
      Context &ctx,
      Material &m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!checkPointByPointLayout(ctx, m)) {
      return {};
    }
    if (!synchronizeOnHost(ctx, m, StateArray::INTERNAL_STATE_VARIABLES, s)) {
      return {};
    }
    return buildPartialQuadratureFunction(
        ctx, m.getPartialQuadratureSpacePointer(),
        getStateManager(m, s).internal_state_variables, m.b.isvs, n,
        m.b.hypothesis);
  }  // end of getInternalStateVariable

  std::optional<ImmutablePartialQuadratureFunctionView>
  getInternalStateVariable(Context &ctx,
                           const Material &m,
                           const std::string_view n,
                           const Material::StateSelection s) noexcept {
    if (!checkPointByPointLayout(ctx, m)) {
      return {};
    }
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::INTERNAL_STATE_VARIABLES, s)) {
      return {};
    }
    return buildImmutablePartialQuadratureFunctionView(
        ctx, m.getPartialQuadratureSpacePointer(),
        getStateManager(m, s).internal_state_variables, m.b.isvs, n,
        m.b.hypothesis);
  }  // end of getInternalStateVariable

  std::optional<PartialQuadratureFunction> getStoredEnergy(
      Context &ctx, Material &m, const Material::StateSelection s) noexcept {
    // one component: same array in both layouts
    if (!synchronizeOnHost(ctx, m, StateArray::STORED_ENERGIES, s)) {
      return {};
    }
    const auto &sm = getStateManager(m, s);
    if (!sm.b.computesStoredEnergy) {
      return ctx.registerErrorMessage(
          "behaviour does not compute the stored energy");
    }
    return PartialQuadratureFunction::borrow(
        ctx, m.getPartialQuadratureSpacePointer(), sm.stored_energies,
        {.data_begin = 0, .data_size = 1, .data_stride = 1});
  }  // end of getStoredEnergy

  std::optional<ImmutablePartialQuadratureFunctionView> getStoredEnergy(
      Context &ctx,
      const Material &m,
      const Material::StateSelection s) noexcept {
    // one component: same array in both layouts
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::STORED_ENERGIES, s)) {
      return {};
    }
    const auto &sm = getStateManager(m, s);
    if (!sm.b.computesStoredEnergy) {
      return ctx.registerErrorMessage(
          "behaviour does not compute the stored energy");
    }
    return ImmutablePartialQuadratureFunctionView(
        m.getPartialQuadratureSpacePointer(), sm.stored_energies,
        {.data_begin = 0, .data_size = 1, .data_stride = 1});
  }  // end of getStoredEnergy

  std::optional<PartialQuadratureFunction> getDissipatedEnergy(
      Context &ctx, Material &m, const Material::StateSelection s) noexcept {
    // one component: same array in both layouts
    if (!synchronizeOnHost(ctx, m, StateArray::DISSIPATED_ENERGIES, s)) {
      return {};
    }
    const auto &sm = getStateManager(m, s);
    if (!sm.b.computesDissipatedEnergy) {
      return ctx.registerErrorMessage(
          "behaviour does not compute the dissipated energy");
    }
    return PartialQuadratureFunction::borrow(
        ctx, m.getPartialQuadratureSpacePointer(), sm.dissipated_energies,
        {.data_begin = 0, .data_size = 1, .data_stride = 1});
  }  // end of getDissipatedEnergy

  std::optional<ImmutablePartialQuadratureFunctionView> getDissipatedEnergy(
      Context &ctx,
      const Material &m,
      const Material::StateSelection s) noexcept {
    // one component: same array in both layouts
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::DISSIPATED_ENERGIES, s)) {
      return {};
    }
    const auto &sm = getStateManager(m, s);
    if (!sm.b.computesDissipatedEnergy) {
      return ctx.registerErrorMessage(
          "behaviour does not compute the dissipated energy");
    }
    return ImmutablePartialQuadratureFunctionView(
        m.getPartialQuadratureSpacePointer(), sm.dissipated_energies,
        {.data_begin = 0, .data_size = 1, .data_stride = 1});
  }  // end of getDissipatedEnergy

  std::optional<real> computeStoredEnergy(
      Context &ctx,
      const AbstractBehaviourIntegrator &bi,
      const Material::StateSelection s) noexcept {
    const auto om = bi.getMaterial(ctx);
    if (isInvalid(om)) {
      return {};
    }
    if (!om->b.computesStoredEnergy) {
      return ctx.registerErrorMessage(
          "behaviour does not compute the stored energy");
    }
    const auto oe = getStoredEnergy(ctx, *om, s);
    if (isInvalid(oe)) {
      return {};
    }
    return computeIntegral<real>(bi, *oe);
  }  // end of computeStoredEnergy

  std::optional<real> computeDissipatedEnergy(
      Context &ctx,
      const AbstractBehaviourIntegrator &bi,
      const Material::StateSelection s) noexcept {
    const auto om = bi.getMaterial(ctx);
    if (isInvalid(om)) {
      return {};
    }
    if (!om->b.computesDissipatedEnergy) {
      return ctx.registerErrorMessage(
          "behaviour does not compute the dissipated energy");
    }
    const auto oe = getDissipatedEnergy(ctx, *om, s);
    if (isInvalid(oe)) {
      return {};
    }
    return computeIntegral<real>(bi, *oe);
  }  // end of computeDissipatedEnergy

}  // end of namespace mfem_mgis
