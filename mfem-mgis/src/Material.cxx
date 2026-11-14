/*!
 * \file   src/Material.cxx
 * \brief  This file implements the `Material` class
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#include <span>
#include <limits>
#include <tuple>
#include <string>
#include <cstddef>
#include <utility>
#include <variant>
#include <algorithm>
#include "mfem/general/device.hpp"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourDescription.hxx"
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

  Material::Material(const PartialQuadratureSpace& s,
                     std::unique_ptr<const Behaviour> b_ptr,
                     const mgis::behaviour::MaterialDataManagerLayout& l)
      : MaterialArrays(Material_dereferenceBehaviour(b_ptr),
                       s.getNumberOfIntegrationPoints(),
                       l),
        MaterialDataManager(Material_dereferenceBehaviour(b_ptr),
                            s.getNumberOfIntegrationPoints(),
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

  void Material::updateMaterialStates(
      const mgis::behaviour::UpdateArguments &args) {
    // MGIS uses the copies of the arrays on the host
    this->synchronizeOnHost();
    auto ops = mgis::HostFieldOperations{};
    mgis::behaviour::update(*this, args, ops);
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
    return this->quadrature_space;
  }  // end of getPartialQuadratureSpace

  void Material::setMacroscopicGradients(std::span<const real> g) {
    if (g.size() != this->s1.gradients_stride) {
      raise(
          "Material::setMacroscopicGradients: "
          "invalid number of components of the gradients");
    }
    std::copy(g.begin(), g.end(), this->macroscopic_gradients.begin());
  }  // end of setMacroscopicGradients

  std::span<const real> Material::getMacroscopicGradients() const noexcept {
    return this->macroscopic_gradients;
  }  // end of getMacroscopicGradients

  static void checkBehaviourSymmetry(const Behaviour &b) {
    if (b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour is not orthotropic)");
    }
  }  // end of checkBehaviourSymmetry

  /*!
   * \return the values of a function defining a material axis at an
   * integration point, copied in an array: the rotation matrices are built
   * from contiguous values, whatever the layout of the function
   * \param[in] f: function
   * \param[in] o: offset of the integration point
   */
  template <std::size_t N>
  static std::array<real, N> getMaterialAxisValues(
      const PartialQuadratureFunction& f, const size_type o) {
    const auto v = f(o);
    auto a = std::array<real, N>{};
    for (std::size_t i = 0; i != N; ++i) {
      a[i] = v[i];
    }
    return a;
  }  // end of getMaterialAxisValues

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
                getMaterialAxisValues<2>(f, o));
          };
    } else {
      raise("Material::setRotationMatrix: unimplemented case yet");
    }
    this->r2D = r;
    this->invalidateRotationMatrices();
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
              getMaterialAxisValues<3>(f2, o));
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
              getMaterialAxisValues<3>(f1, o),
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
              getMaterialAxisValues<3>(f1, o), getMaterialAxisValues<3>(f2, o));
        };
      } else {
        raise("Material::setRotationMatrix: unimplemented case yet");
      }
    } else {
      raise("Material::setRotationMatrix: unimplemented case yet");
    }
    this->r3D = r;
    this->invalidateRotationMatrices();
  }  // end of setRotationMatrix

  void Material::invalidateRotationMatrices() noexcept {
    this->kept_rotation_matrices_computed = {false, false};
  }  // end of invalidateRotationMatrices

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

  std::pair<const mfem::Vector *, mgis::FieldView<const real>>
  Material::getRotationMatrices(const bool t) const {
    const auto i = t ? std::size_t{1} : std::size_t{0};
    auto &r = this->kept_rotation_matrices[i];
    if (!this->kept_rotation_matrices_computed[i]) {
      std::ignore = this->computeRotationMatrices(r, t);
      this->kept_rotation_matrices_computed[i] = true;
    }
    const auto np = static_cast<size_type>(this->n);
    const auto size = static_cast<std::size_t>(r.Size());
    if (size == 0) {
      return {&r, mgis::makePointByPointFieldView(std::span<const real>{}, 0, 9,
                                                  9)};
    }
    const auto matrices = std::span<const real>(r.HostRead(), size);
    // one matrix: uniform rotation (or a single point)
    if (size == 9) {
      return {&r, mgis::makeUniformFieldView(matrices, np)};
    }
    return {&r, mgis::makePointByPointFieldView(matrices, np, 9, 9)};
  }  // end of getRotationMatrices

  void Material::computeValuesInGlobalFrame(mfem::Vector& v,
                                            const StateArray a,
                                            const StateSelection s,
                                            const mgis::FieldLayout ol) const {
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
    // the gradients are rotated by the transposed matrices, from the
    // material frame to the global frame
    const auto [r, hmatrices] = this->getRotationMatrices(gradients);
    const auto hsource = gradients
                             ? mgis::behaviour::getGradientsView(sm)
                             : mgis::behaviour::getThermodynamicForcesView(sm);
    // the values are rotated on the host, or on the copies of the arrays on
    // the device if the computations of the behaviour can't access the
    // memory of the host
    const auto host = mgis::behaviour::isMemoryLocationAccessible(
        this->b, mgis::MemoryLocation::HOST);
    if (host) {
      this->synchronizeOnHostForReading(a, s);
    }
    const auto l =
        host ? mgis::MemoryLocation::HOST : getDeviceMemoryLocation();
    const auto &array = this->getArray(a, s);
    const auto source = host ? mgis::FieldView<const real>(hsource)
                             : relocate(mgis::FieldView<const real>(hsource),
                                        this->read(a, s), array.Size());
    const auto matrices =
        host ? hmatrices : relocate(hmatrices, r->Read(), r->Size());
    auto *const values = host ? v.HostWrite() : v.Write();
    const auto d =
        mgis::makeFieldView(std::span<real>(values, size), ol, np, nc);
    if (gradients) {
      mgis::behaviour::rotateGradients(d, this->b, source, matrices, l);
    } else {
      mgis::behaviour::rotateThermodynamicForces(d, this->b, source, matrices,
                                                 l);
    }
  }  // end of computeValuesInGlobalFrame

  void Material::computeGradientsInGlobalFrame(
      mfem::Vector& v,
      const StateSelection s,
      const mgis::FieldLayout l) const {
    this->computeValuesInGlobalFrame(v, StateArray::GRADIENTS, s, l);
  }  // end of computeGradientsInGlobalFrame

  void Material::computeThermodynamicForcesInGlobalFrame(
      mfem::Vector& v,
      const StateSelection s,
      const mgis::FieldLayout l) const {
    this->computeValuesInGlobalFrame(v, StateArray::THERMODYNAMIC_FORCES, s, l);
  }  // end of computeThermodynamicForcesInGlobalFrame

  /*!
   * \return a view of a copy, readable on the device, of the values of a
   * material property, of an external state variable or of the mass
   * density, given by a view of their values on the host
   * \param[out] c: copy, stored densely in the layout of the values on the
   * host (point by point unless the components are not contiguous), one
   * point for uniform values
   * \param[in] v: view of the values on the host
   */
  static mgis::FieldView<const real> copyToDevice(
      mfem::Vector& c, const mgis::FieldView<const real>& v) {
    const auto np = v.isUniform() ? size_type{1} : v.getNumberOfPoints();
    const auto nc = v.getNumberOfComponents();
    const auto size = static_cast<std::size_t>(np) * nc;
    if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      raise("Material::makeDeviceView: array too large for MFEM");
    }
    c.UseDevice(true);
    c.SetSize(static_cast<int>(size));
    const auto l = (v.isUniform() || (nc <= 1) || v.hasContiguousComponents())
                       ? mgis::FieldLayout::POINT_BY_POINT
                       : mgis::FieldLayout::COMPONENT_BY_COMPONENT;
    const auto dst =
        mgis::makeFieldView(std::span<real>(c.HostWrite(), size), l, np, nc);
    for (size_type i = 0; i != np; ++i) {
      for (size_type j = 0; j != nc; ++j) {
        dst(i, j) = v(i, j);
      }
    }
    const auto d = std::span<const real>(c.Read(), size);
    if (v.isUniform()) {
      return mgis::makeUniformFieldView(d, v.getNumberOfPoints());
    }
    return mgis::makeFieldView(d, l, np, nc);
  }  // end of copyToDevice

  mgis::behaviour::MaterialDataView Material::makeDeviceView(
      std::vector<mfem::Vector>& inputs,
      const bool tangent,
      const bool gradients) {
    auto& arrays = this->getArrays();
    // view of the arrays on the host, whose views are replaced by views of
    // the copies on the device
    auto v = mgis::behaviour::make_unchecked_view(*this);
    v.memory_location = getDeviceMemoryLocation();
    // arrays read by the behaviour (`Read`), which are not modified by it,
    // or read and written (`ReadWrite`)
    const auto relocateArray =
        [&arrays](mgis::FieldView<real>& view, const StateArray a,
                  const TimeStepStage s, const bool output) {
          const auto& array = arrays.getArray(a, s);
          auto* const p = output ? arrays.readWrite(a, s)
                                 : const_cast<real*>(arrays.read(a, s));
          view = relocate(view, p, array.Size());
        };
    for (const auto& [view, a] :
         {std::pair{&v.s0.gradients, StateArray::GRADIENTS},
          std::pair{&v.s0.thermodynamic_forces,
                    StateArray::THERMODYNAMIC_FORCES},
          std::pair{&v.s0.internal_state_variables,
                    StateArray::INTERNAL_STATE_VARIABLES},
          std::pair{&v.s0.stored_energies, StateArray::STORED_ENERGIES},
          std::pair{&v.s0.dissipated_energies,
                    StateArray::DISSIPATED_ENERGIES}}) {
      relocateArray(*view, a, bts, false);
    }
    relocateArray(v.s1.gradients, StateArray::GRADIENTS, ets, gradients);
    for (const auto& [view, a] :
         {std::pair{&v.s1.thermodynamic_forces,
                    StateArray::THERMODYNAMIC_FORCES},
          std::pair{&v.s1.internal_state_variables,
                    StateArray::INTERNAL_STATE_VARIABLES},
          std::pair{&v.s1.stored_energies, StateArray::STORED_ENERGIES},
          std::pair{&v.s1.dissipated_energies,
                    StateArray::DISSIPATED_ENERGIES}}) {
      relocateArray(*view, a, ets, true);
    }
    if (tangent) {
      v.K = relocate(v.K, arrays.writeTangentOperatorBlocks(),
                     arrays.getTangentOperatorBlocks().Size());
    } else {
      v.K = {};
    }
    v.speed_of_sound = {};
    // material properties, external state variables and mass density
    // copied to the device
    auto views = std::vector<mgis::FieldView<const real>*>{};
    for (auto* const st : {&v.s0, &v.s1}) {
      for (auto& mp : st->material_properties) {
        views.push_back(&mp);
      }
      for (auto& esv : st->external_state_variables) {
        views.push_back(&esv);
      }
      views.push_back(&(st->mass_density));
    }
    inputs.resize(views.size());
    for (std::size_t i = 0; i != views.size(); ++i) {
      if (!views[i]->empty()) {
        *(views[i]) = copyToDevice(inputs[i], *(views[i]));
      }
    }
    return v;
  }  // end of makeDeviceView

  Material::~Material() = default;

  /*!
   * \return the specifications of a partial quadrature function viewing the
   * values of a variable of a material, whatever the layout of the array
   * holding it
   * \param[in] values: values of the array holding the variable
   * \param[in] v: view of the variable, given by MGIS
   */
  template <typename ValueType>
  [[nodiscard]] static ViewSpecifications getViewSpecifications(
      std::span<ValueType> values,
      const mgis::FieldView<ValueType>& v) noexcept {
    // the strides of an empty view may be null
    return {.data_begin = static_cast<size_type>(v.data() - values.data()),
            .data_size = static_cast<size_type>(v.getNumberOfComponents()),
            .data_stride = std::max(static_cast<size_type>(v.getPointStride()),
                                    size_type{1}),
            .data_component_stride = std::max(
                static_cast<size_type>(v.getComponentStride()), size_type{1})};
  }  // end of getViewSpecifications

  /*!
   * \return a partial quadrature function viewing the values of a variable
   * of a material
   * \param[in, out] ctx: execution context
   * \param[in] qs: quadrature space
   * \param[in] values: values of the array holding the variable
   * \param[in] get_view: function returning the view of the variable, given
   * by MGIS, which may throw if the variable does not exist
   */
  [[nodiscard]] static std::optional<PartialQuadratureFunction>
  buildPartialQuadratureFunction(
      Context& ctx,
      const PartialQuadratureSpace& qs,
      std::span<mgis::real> values,
      const std::function<mgis::FieldView<mgis::real>()>& get_view) noexcept {
    try {
      const auto v = get_view();
      return PartialQuadratureFunction::borrow(
          ctx, qs, values, getViewSpecifications(values, v));
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of buildPartialQuadratureFunction

  //! \brief immutable counterpart of `buildPartialQuadratureFunction`
  [[nodiscard]] static std::optional<ImmutablePartialQuadratureFunctionView>
  buildImmutablePartialQuadratureFunctionView(
      Context& ctx,
      const PartialQuadratureSpace& qs,
      std::span<const mgis::real> values,
      const std::function<mgis::FieldView<const mgis::real>()>&
          get_view) noexcept {
    try {
      const auto v = get_view();
      return ImmutablePartialQuadratureFunctionView(
          qs, values, getViewSpecifications(values, v));
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of buildImmutablePartialQuadratureFunctionView

  bool isIntegratedByBatchesOnDevice(const Behaviour& b) noexcept {
    return mfem::Device::Allows(mfem::Backend::DEVICE_MASK) &&
           mgis::behaviour::isGenericParallelBehaviour(b) &&
           mgis::behaviour::isMemoryLocationAccessible(
               b, getDeviceMemoryLocation());
  }  // end of isIntegratedByBatchesOnDevice

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
      Context& ctx,
      Material& m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!synchronizeOnHost(ctx, m, StateArray::GRADIENTS, s)) {
      return {};
    }
    auto& sm = getStateManager(m, s);
    return buildPartialQuadratureFunction(
        ctx, m.getPartialQuadratureSpace(), sm.gradients,
        [&sm, n] { return mgis::behaviour::getGradientView(sm, n); });
  }  // end of getGradient

  std::optional<ImmutablePartialQuadratureFunctionView> getGradient(
      Context &ctx,
      const Material &m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!synchronizeOnHostForReading(ctx, m, StateArray::GRADIENTS, s)) {
      return {};
    }
    const auto& sm = getStateManager(m, s);
    return buildImmutablePartialQuadratureFunctionView(
        ctx, m.getPartialQuadratureSpace(), sm.gradients,
        [&sm, n] { return mgis::behaviour::getGradientView(sm, n); });
  }  // end of getGradient

  std::optional<PartialQuadratureFunction> getThermodynamicForce(
      Context& ctx,
      Material& m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!synchronizeOnHost(ctx, m, StateArray::THERMODYNAMIC_FORCES, s)) {
      return {};
    }
    auto& sm = getStateManager(m, s);
    return buildPartialQuadratureFunction(
        ctx, m.getPartialQuadratureSpace(), sm.thermodynamic_forces,
        [&sm, n] { return mgis::behaviour::getThermodynamicForceView(sm, n); });
  }  // end of getThermodynamicForce

  std::optional<ImmutablePartialQuadratureFunctionView> getThermodynamicForce(
      Context& ctx,
      const Material& m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::THERMODYNAMIC_FORCES, s)) {
      return {};
    }
    const auto& sm = getStateManager(m, s);
    return buildImmutablePartialQuadratureFunctionView(
        ctx, m.getPartialQuadratureSpace(), sm.thermodynamic_forces,
        [&sm, n] { return mgis::behaviour::getThermodynamicForceView(sm, n); });
  }  // end of getThermodynamicForce

  std::optional<PartialQuadratureFunction> getInternalStateVariable(
      Context& ctx,
      Material& m,
      const std::string_view n,
      const Material::StateSelection s) noexcept {
    if (!synchronizeOnHost(ctx, m, StateArray::INTERNAL_STATE_VARIABLES, s)) {
      return {};
    }
    auto& sm = getStateManager(m, s);
    return buildPartialQuadratureFunction(
        ctx, m.getPartialQuadratureSpace(), sm.internal_state_variables,
        [&sm, n] {
          return mgis::behaviour::getInternalStateVariableView(sm, n);
        });
  }  // end of getInternalStateVariable

  std::optional<ImmutablePartialQuadratureFunctionView>
  getInternalStateVariable(Context& ctx,
                           const Material& m,
                           const std::string_view n,
                           const Material::StateSelection s) noexcept {
    if (!synchronizeOnHostForReading(
            ctx, m, StateArray::INTERNAL_STATE_VARIABLES, s)) {
      return {};
    }
    const auto& sm = getStateManager(m, s);
    return buildImmutablePartialQuadratureFunctionView(
        ctx, m.getPartialQuadratureSpace(), sm.internal_state_variables,
        [&sm, n] {
          return mgis::behaviour::getInternalStateVariableView(sm, n);
        });
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
        ctx, m.getPartialQuadratureSpace(), sm.stored_energies,
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
        m.getPartialQuadratureSpace(), sm.stored_energies,
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
        ctx, m.getPartialQuadratureSpace(), sm.dissipated_energies,
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
        m.getPartialQuadratureSpace(), sm.dissipated_energies,
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
