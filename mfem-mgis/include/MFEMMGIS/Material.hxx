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

#include "MGIS/FieldView.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/MaterialDataView.hxx"
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
    //! \return the macroscopic gradients
    [[nodiscard]] std::span<const real> getMacroscopicGradients()
        const noexcept;
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
     * \brief makes the copies of the arrays on the host valid at the
     * beginning of a phase of accesses on the host, and invalidates their
     * copies on the device (see `MaterialArrays` for the complete rules).
     * Required before a direct access to the states (`s0`, `s1`) or to the
     * tangent operator blocks (`K`) after a computation on the device. It is
     * called by the mutable accessors (`getInternalStateVariable`, ...), by
     * the integration of the behaviours on the host (`integrate` and `solve`
     * of the non linear evolution problems, `computeNextState` of the
     * point-wise models), by `updateMaterialStates` and by
     * `revertMaterialStates`.
     *
     * \note a behaviour integrated on the device (see
     * `BehaviourIntegratorBase::integratesBatchesOnDevice`) is integrated
     * without synchronization on the host: between the first integration of
     * a time step and the end of the resolution, the copies on the host of
     * its arrays are outdated, and a direct access to its states must be
     * preceded by a synchronization. Its update and its revert are made on
     * the device and end with a synchronization on the host, as for a
     * behaviour integrated on the host.
     *
     * \note the arrays of a material used by another one as external
     * storage (material properties or external state variables given with
     * the `EXTERNAL_STORAGE` mode) are not synchronized by the other
     * material: the code which integrates the other material must
     * synchronize this one first.
     */
    using MaterialArrays::synchronizeOnHost;
    /*!
     * \brief makes the copies of the arrays on the host valid at the
     * beginning of a phase of read-only accesses on the host, keeping their
     * copies on the device valid. It is called by the read-only accessors,
     * before the post-processings, by the functions computing mean values
     * or resultant forces from the thermodynamic forces, and after the
     * integration of a point-wise model on the device (see
     * `PointWiseModel::computeNextState`).
     */
    using MaterialArrays::synchronizeOnHostForReading;
    /*!
     * \brief copies the state at the end of the time step on the state at
     * the beginning of the time step (see `mgis::behaviour::update`), after
     * the synchronization of the arrays on the host, where MGIS uses them.
     * The arrays are copied by `mgis::HostFieldOperations`.
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
    /*!
     * \brief computes the rotation matrices from the global frame to the
     * material frame at the integration points (orthotropic behaviours), in
     * the form expected by the rotation functions of MGIS taking views (see
     * `mgis::behaviour::rotateGradients`): one matrix for an uniform
     * rotation, one matrix per integration point otherwise, 9 values per
     * matrix.
     * \return a view of the rotation matrices of all the integration points
     * (with a null stride between points for an uniform rotation), on the
     * copy on the host of the given array
     * \param[out] r: array of the rotation matrices, resized
     * \param[in] t: if true, the transposed matrices are computed, i.e. the
     * rotation matrices from the material frame to the global frame
     * \note the material axes given by partial quadrature functions are read
     * at each call
     */
    [[nodiscard]] mgis::FieldView<const real> computeRotationMatrices(
        mfem::Vector &, const bool = false) const;
    /*!
     * \brief computes the gradients of a state in the global frame at all
     * the integration points (orthotropic behaviours), by one rotation of
     * all the points (see `mgis::behaviour::rotateGradients`, which uses the
     * rotation function of the batches for a behaviour generated by the
     * `generic-parallel` interface), on the host, or on the copies of the
     * arrays on the device if the computations of the behaviour can't
     * access the memory of the host: the values must then be read through
     * the accessors of `mfem::Vector` which keep track of the valid copies
     * (`HostRead`, `Read`), not through `GetData` or `operator[]`
     * \param[out] v: values, stored point by point whatever the layout of
     * the arrays of the material (resized)
     * \param[in] s: state
     */
    void computeGradientsInGlobalFrame(
        mfem::Vector &, const StateSelection = END_OF_TIME_STEP) const;
    /*!
     * \brief computes the thermodynamic forces of a state in the global
     * frame at all the integration points (orthotropic behaviours), by one
     * rotation of all the points (see `computeGradientsInGlobalFrame`)
     * \param[out] v: values, stored point by point (resized)
     * \param[in] s: state
     */
    void computeThermodynamicForcesInGlobalFrame(
        mfem::Vector &, const StateSelection = END_OF_TIME_STEP) const;
    /*!
     * \return a view of MGIS on the copies of the arrays of the material on
     * the device, for an integration of its behaviour on the device by
     * `mgis::behaviour::integrate` with the batch execution policy of the
     * view (see `isIntegratedByBatchesOnDevice`). The arrays of the state at
     * the beginning of the time step are read, the ones of the state at the
     * end of the time step are read and written, except its gradients,
     * which are only read unless the third argument is true. The material
     * properties, the external state variables and the mass density, given
     * on the host, are copied to the device.
     * \param[in, out] inputs: copies on the device of the material properties,
     * of the external state variables and of the mass density, used by the
     * view
     * \param[in] tangent: if true, the view writes the tangent operator
     * blocks, otherwise it has no tangent operator blocks
     * \param[in] gradients: if true, the gradients at the end of the time
     * step are read and written (for instance rotated in place)
     */
    [[nodiscard]] mgis::behaviour::MaterialDataView makeDeviceView(
        std::vector<mfem::Vector>&, const bool, const bool);
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
    /*!
     * \brief computes the values of an array of a state in the global frame
     * (see `computeGradientsInGlobalFrame`)
     * \param[out] v: values
     * \param[in] a: array (gradients or thermodynamic forces)
     * \param[in] s: state
     */
    void computeValuesInGlobalFrame(mfem::Vector &,
                                    const StateArray,
                                    const StateSelection) const;
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
   * \return true if a behaviour is integrated by batches on the device
   * enabled in MFEM (see `Material::makeDeviceView`): a device is enabled,
   * the behaviour is generated by the `generic-parallel` interface and its
   * computations can access the memory of the device (see
   * `getDeviceMemoryLocation`)
   * \param[in] b: behaviour
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isIntegratedByBatchesOnDevice(
      const Behaviour&) noexcept;

  /*!
   * \return true if the arrays of the material are in the point by point
   * layout, required by the accessors and the evaluators based on a stride
   * between points; otherwise, an error message is registered
   * \param[in, out] eh: error handler
   * \param[in] m: material
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool checkPointByPointLayout(
      AbstractErrorHandler &, const Material &) noexcept;

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
