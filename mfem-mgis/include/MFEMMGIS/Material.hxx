/*!
 * \file   include/MFEMMGIS/Material.hxx
 * \brief  This file declares the `Material` class
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_MATERIAL_HXX
#define LIB_MFEM_MGIS_MATERIAL_HXX

#include <span>
#include <array>
#include <vector>
#include <memory>
#include <optional>
#include <string_view>

#ifdef MGIS_FUNCTION_SUPPORT
#include "MGIS/Function/EvaluatorConcept.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/MaterialArrays.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/RotationMatrix.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  // forward declarations
  struct PartialQuadratureSpace;
  struct AbstractBehaviourIntegrator;

  /*!
   * \brief a simple structure describing a material and the associated data.
   *
   * The arrays of the data at the integration points are in the memory of
   * MFEM (see `MaterialArrays`, and its rules of validity of the copies on
   * the host and on the device). They are held by a base class declared
   * before the material data manager of MGIS, so that they are created
   * before it.
   */
  struct MFEM_MGIS_EXPORT Material
      : private MaterialArrays,
        public mgis::behaviour::MaterialDataManager {
    //! \brief a simple alias used to select the state of the material
    using StateSelection = TimeStepStage;
    //! \brief beginning of the time step
    static constexpr TimeStepStage BEGINNING_OF_TIME_STEP =
        TimeStepStage::BEGINNING_OF_TIME_STEP;
    //! \brief end of the time step
    static constexpr TimeStepStage END_OF_TIME_STEP =
        TimeStepStage::END_OF_TIME_STEP;
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     * \param[in] l: layout of the arrays of the data at the integration
     * points. The component by component layout is not supported by the
     * accessors based on a stride between points (`getGradient`,
     * `getInternalStateVariable`...), nor by the behaviour integrators of
     * mfem-mgis, which use the point by point layout.
     */
    Material(std::shared_ptr<const PartialQuadratureSpace> s,
             std::unique_ptr<const Behaviour> b_ptr,
             const mgis::FieldLayout l = mgis::FieldLayout::POINT_BY_POINT);
    /*!
     * \brief set the macroscopic gradients
     * \param[in] g: macroscopic gradients
     */
    void setMacroscopicGradients(std::span<const real> g);
    /*!
     * \brief set the rotation matrix
     * \param[in] r: rotation matrix
     * \note this call is only meaningful in 2D hypotheses for orthotropic
     * behaviours.
     */
    void setRotationMatrix(const RotationMatrix2D &r);
    /*!
     * \brief set the rotation matrix
     * \param[in] r: rotation matrix
     * \note this call is only meaningful in 3D for orthotropic behaviours
     */
    void setRotationMatrix(const RotationMatrix3D &r);
    /*!
     * \brief makes the copies of the arrays on the host valid before an
     * access on the host, and invalidates their copies on the device (see
     * `MaterialArrays`). Required before a direct access to the states (`s0`,
     * `s1`) or to the tangent operator blocks (`K`) after a computation on
     * the device. The accessors (`getInternalStateVariable`, ...), the
     * integration of the behaviours, `update`, `revert` and the
     * post-processings call it.
     */
    using MaterialArrays::synchronizeOnHost;
    /*!
     * \brief makes the copies of the arrays on the host valid before a
     * read-only access on the host, keeping their copies on the device valid
     */
    using MaterialArrays::synchronizeOnHostForReading;
    /*!
     * \brief copies the state at the end of the time step on the state at
     * the beginning of the time step (see `mgis::behaviour::update`). The
     * arrays are copied by kernels of MFEM (`batch::MFEMFieldOperations`),
     * in the memory where MGIS uses them (the host).
     */
    void updateMaterialStates();
    /*!
     * \brief copies the state at the beginning of the time step on the state
     * at the end of the time step (see `mgis::behaviour::revert` and
     * `updateMaterialStates`)
     */
    void revertMaterialStates();
    //! \return the arrays of the data at the integration points
    MaterialArrays &getArrays();
    //! \return the arrays of the data at the integration points
    const MaterialArrays &getArrays() const;
    //! \return the quadrature space
    const PartialQuadratureSpace &getPartialQuadratureSpace() const;
    //! \return the quadrature space
    std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const;
    /*!
     * \brief compute the rotation matrix at an integration point
     * \return the rotation matrix for the given integration point
     * \param[in] i: offset of the integration point
     * \note this method is only valid for orthotropic behaviours
     */
    std::array<real, 9u> getRotationMatrixAtIntegrationPoint(
        const size_type i) const;
    //! \brief destructor
    ~Material();

   protected:
    /*!
     * \brief underlying quadrature space
     */
    const std::shared_ptr<const PartialQuadratureSpace> quadrature_space;
    /*!
     * \brief macroscopic gradients
     */
    std::vector<real> macroscopic_gradients;

   protected:
    //! \brief the rotation matrix in 2D
    RotationMatrix2D r2D;
    //! \brief the rotation matrix in 3D
    RotationMatrix3D r3D;
    //! \brief pointer to a function returning the rotation matrix
    std::array<real, 9u> (*get_rotation_fct_ptr)(const RotationMatrix2D &,
                                                 const RotationMatrix3D &,
                                                 const size_type);

   private:
    //! \brief copy constructor (disabled)
    Material(const Material &) = delete;
    //! \brief move constructor (disabled)
    Material(Material &&) = delete;
    //! \brief standard assignment (disabled)
    Material &operator=(const Material &) = delete;
    //! \brief move assignment (disabled)
    Material &operator=(Material &&) = delete;

    /*!
     * \brief underlying behaviour. Only stored for memory management.
     * \note The behaviour can be accessed through the `b` member which is
     * inherited from the `mgis::behaviour::MaterialDataManager` class.
     */
    const std::unique_ptr<const Behaviour> behaviour_ptr;

  };  // end of struct Material

  /*!
   * \brief get a gradient of the material
   * \return a partial quadrature function for the given gradient
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the gradient
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<PartialQuadratureFunction>
  getGradient(
      Context &ctx,
      Material &m,
      const std::string_view n,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get a gradient of the material
   * \return a partial quadrature function for the given gradient
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the gradient
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<ImmutablePartialQuadratureFunctionView>
  getGradient(
      Context &ctx,
      const Material &m,
      const std::string_view n,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get a thermodynamic force of the material
   * \return a partial quadrature function for the given thermodynamic force
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the thermodynamic force
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<PartialQuadratureFunction>
  getThermodynamicForce(
      Context &ctx,
      Material &m,
      const std::string_view n,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get a thermodynamic force of the material
   * \return a partial quadrature function for the given thermodynamic force
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the thermodynamic force
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<ImmutablePartialQuadratureFunctionView>
  getThermodynamicForce(
      Context &ctx,
      const Material &m,
      const std::string_view n,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get an internal state variable of the material
   * \return a partial quadrature function for the given state variable
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the state variable
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<PartialQuadratureFunction>
  getInternalStateVariable(
      Context &ctx,
      Material &m,
      const std::string_view n,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get an internal state variable of the material
   * \return a partial quadrature function for the given state variable
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the state variable
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<ImmutablePartialQuadratureFunctionView>
  getInternalStateVariable(
      Context &ctx,
      const Material &m,
      const std::string_view n,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get the stored energy of the material
   * \return a partial quadrature function holding the stored energy
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<PartialQuadratureFunction>
  getStoredEnergy(
      Context &ctx,
      Material &m,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get the stored energy of the material
   * \return a partial quadrature function holding the stored energy
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::optional<ImmutablePartialQuadratureFunctionView>
  getStoredEnergy(
      Context &ctx,
      const Material &m,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get the dissipated energy of the material
   * \return a partial quadrature function holding the dissipated energy
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<PartialQuadratureFunction>
  getDissipatedEnergy(
      Context &ctx,
      Material &m,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief get the dissipated energy of the material
   * \return a partial quadrature function holding the dissipated energy
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]] std::optional<ImmutablePartialQuadratureFunctionView>
  getDissipatedEnergy(
      Context &ctx,
      const Material &m,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief compute the energy stored by the whole material
   * \return the energy stored by the whole material, empty if the behaviour
   * does not compute it or on failure
   * \param[in, out] ctx: execution context
   * \param[in] bi: behaviour integrator
   * \param[in] s: selection of the state
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<real> computeStoredEnergy(
      Context &ctx,
      const AbstractBehaviourIntegrator &bi,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;
  /*!
   * \brief compute the energy dissipated by the whole material
   * \return the energy dissipated by the whole material, empty if the
   * behaviour does not compute it or on failure
   * \param[in, out] ctx: execution context
   * \param[in] bi: behaviour integrator
   * \param[in] s: selection of the state
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<real> computeDissipatedEnergy(
      Context &ctx,
      const AbstractBehaviourIntegrator &bi,
      const Material::StateSelection s = Material::END_OF_TIME_STEP) noexcept;

  /*!
   * \brief get the state of the material at a time step stage
   * \return the state of the material at the given time step stage
   * \param[in] m: material
   * \param[in] s: time step stage
   */
  [[nodiscard]] mgis::behaviour::MaterialStateManager &getStateManager(
      Material &m, const Material::StateSelection s) noexcept;

  /*!
   * \brief get the state of the material at a time step stage
   * \return the state of the material at the given time step stage
   * \param[in] m: material
   * \param[in] s: time step stage
   */
  [[nodiscard]] const mgis::behaviour::MaterialStateManager &getStateManager(
      const Material &m, const Material::StateSelection s) noexcept;

#ifdef MGIS_FUNCTION_SUPPORT

  //! \brief an evaluator returning the rotation matrix
  struct RotationMatrixEvaluator {
    /*!
     * \brief constructor
     * \param[in] m: material
     */
    inline RotationMatrixEvaluator(const Material &m) : material(m) {}
    /*!
     * \brief perform consistency checks
     * \return true on success
     * \param[in, out] ctx: error handler
     */
    inline bool check(AbstractErrorHandler &ctx) const {
      if (this->material.b.symmetry !=
          mgis::behaviour::Behaviour::ORTHOTROPIC) {
        return ctx.registerErrorMessage(
            "considered material is not orthotropic");
      }
      return true;
    }
    //! \return the quadrature space
    inline const PartialQuadratureSpace &getSpace() const {
      return this->material.getPartialQuadratureSpace();
    }
    /*!
     * \brief access operator
     * \return the rotation matrix at the given integration point
     * \param[in] i: offset of the integration point
     */
    inline std::array<real, 9u> operator()(const size_type i) const {
      return this->material.getRotationMatrixAtIntegrationPoint(i);
    }

   private:
    //! \brief underlying material
    const Material &material;
  };

  /*!
   * \brief get the quadrature space of an evaluator
   * \return the quadrature space
   * \param[in] e: evaluator
   */
  inline const PartialQuadratureSpace &getSpace(
      const RotationMatrixEvaluator &e) {
    return e.getSpace();
  }  // end of getSpace

  /*!
   * \brief perform consistency checks
   * \return true on success
   * \param[in, out] eh: error handler
   * \param[in] e: evaluator
   */
  bool check(AbstractErrorHandler &eh, const RotationMatrixEvaluator &e);

  /*!
   * \brief return the number of components
   * \param[in] e: evaluator
   * \return the number of components of a rotation matrix
   */
  constexpr mgis::size_type getNumberOfComponents(
      const RotationMatrixEvaluator &e) noexcept;

  static_assert(mgis::function::EvaluatorConcept<RotationMatrixEvaluator>);
  static_assert(!mgis::function::FunctionConcept<RotationMatrixEvaluator>);

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Material.ixx"

#endif /* LIB_MFEM_MGIS_MATERIAL_HXX */
