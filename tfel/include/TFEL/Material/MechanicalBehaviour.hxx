/*!
 * \file   include/TFEL/Material/MechanicalBehaviour.hxx
 * \brief  This file declares the MechanicalBehaviour class
 * \author Thomas Helfer
 * \date 22/09/2006
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_TFEL_MATERIAL_MECHANICALBEHAVIOUR_HXX
#define LIB_TFEL_MATERIAL_MECHANICALBEHAVIOUR_HXX

#include <utility>
#include <concepts>
#include "TFEL/Config/TFELTypes.hxx"
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperatorBase.hxx"

namespace tfel::material {

  /*!
   * \class MechanicalBehaviourBase
   * \brief This class declares enums for mechanical behaviours.
   * \author Thomas Helfer
   * \date 30/07/2013
   */
  struct MechanicalBehaviourBase {
    /*!
     * An indication of the type of the behaviour treated
     * If the behaviour is a small strain standard behaviour, then:
     * - the only driving variable is the total strain 'eto' (symmetric tensor)
     * - the only thermodynamic force is the stress    'sig' (symmetric tensor)
     * If the behaviour is a cohesive zone model, then:
     * - the only driving variable is the opening displacement (tvector)
     * - the only thermodynamic force is the traction (tvector)
     */
    enum BehaviourType {
      GENERALBEHAVIOUR = 0,
      STANDARDSTRAINBASEDBEHAVIOUR = 1,
      STANDARDFINITESTRAINBEHAVIOUR = 2,
      COHESIVEZONEMODEL = 3
    };  // end of enum MFrontBehaviourType
    /*!
     * An indication of the kinematic of the behaviour treated
     */
    enum Kinematic {
      UNDEFINEDKINEMATIC = 0,
      SMALLSTRAINKINEMATIC = 1,
      COHESIVEZONEKINEMATIC = 2,
      FINITESTRAINKINEMATIC_F_CAUCHY = 3,
      FINITESTRAINKINEMATIC_ETO_PK1 = 4
    };  // end of enum MFrontBehaviourType
    /*!
     * \brief return values of the integrate method
     */
    enum ExitStatus {
      FAILURE = 0,        //<! Integration failed
      SUCCESS = 1,        //<! Integration is a success
      UNRELIABLE_RESULTS  //<! Integration succeed, but one or more
                          //   internal criteria show that the
                          //   results may be inaccurate
    };
    /*!
     * \brief a small structure representing the result of the behaviour
     * integration
     * See Issue #296 for the rationale behind this class.
     * <https://github.com/thelfer/tfel/issues/296>
     */
    struct IntegrationResult {
      //! \brief default constructor
      TFEL_HOST_DEVICE constexpr IntegrationResult() noexcept
          : status(SUCCESS) {}
      //! \brief move constructor
      constexpr IntegrationResult(IntegrationResult&&) noexcept = default;
      //! \brief default constructor
      constexpr IntegrationResult(const IntegrationResult&) noexcept = default;
      //! \brief move assignement
      constexpr IntegrationResult& operator=(IntegrationResult&&) noexcept =
          default;
      //! \brief standard assignement
      constexpr IntegrationResult& operator=(
          const IntegrationResult&) noexcept = default;
      // NOLINTBEGIN(google-explicit-constructor)
      /*!
       * \brief constructor from an ExitStatus
       */
      TFEL_HOST_DEVICE constexpr explicit(false)
          IntegrationResult(const ExitStatus s) noexcept
          : status(s) {}  // end of IntegrationResult
      /*!
       * \brief constructor from a boolean value
       */
      TFEL_HOST_DEVICE constexpr explicit(false)
          IntegrationResult(const bool b) noexcept
          : status(b ? SUCCESS : FAILURE) {}  // end of IntegrationResult

      //! \brief convertion operator to the ExistStatus enumeration
      TFEL_HOST_DEVICE constexpr explicit(false) operator ExitStatus()
          const noexcept {
        return this->status;
      }
      // NOLINTEND(google-explicit-constructor)
      //! \brief destructor
      inline ~IntegrationResult() noexcept = default;

     private:
      ExitStatus status;
    };  // end of enum IntegrationResult
    //! \brief list of possible tangent operator type
    enum SMType {
      ELASTIC,
      SECANTOPERATOR,
      TANGENTOPERATOR,
      CONSISTENTTANGENTOPERATOR,
      NOSTIFFNESSREQUESTED
    };  // end of enum SMType
  };    // end of struct MechanicalBehaviourBase

  /*!
   * a trait class describing which tangent operators can be
   * computed by a mechanical behaviour.
   */
  template <MechanicalBehaviourBase::BehaviourType btype>
  struct TangentOperatorTraits {
    /*!
     * Most mechanical behaviour can only compute one tangent
     * operator. Finite strain beaviours are a noticeable exception.
     */
    enum SMFlag { STANDARDTANGENTOPERATOR };  // end of enum Flag
  };                                          // end of struct

  /*!
   * a trait class describing which tangent operators can be
   * computed by a mechanical behaviour.
   */
  template <>
  struct TangentOperatorTraits<
      MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR>
      : public FiniteStrainBehaviourTangentOperatorBase {
    /*!
     * Most mechanical behaviour can only compute one tangent
     * operator. Finite strain beaviours are a noticeable exception.
     */
    typedef FiniteStrainBehaviourTangentOperatorBase::Flag SMFlag;
  };  // end of struct

  /*!
   * \class MechanicalBehaviour
   * \brief This class is the base class of mechanical behaviours. It
   * defines the types and constants shared by all the behaviours.
   *
   * A mechanical behaviour shall provide the following methods (see the
   * `MechanicalBehaviourConcept` concept):
   *
   * - `computePredictionOperator(smflag, smt)`, which only computes a
   *   prediction operator;
   * - `integrate(smflag, smt)`, which determines the values of the
   *   internal state variables at the end of the time step and computes the
   *   tangent operator if requested;
   * - `getMinimalTimeStepScalingFactor()`, which returns the minimal
   *   scaling factor to be used if the integration failed;
   * - `computeAPrioriTimeStepScalingFactor(rdt)` and
   *   `computeAPosterioriTimeStepScalingFactor(rdt)`, which return a pair
   *   containing a boolean (saying if the integration shall be performed,
   *   or considered a success) and a time step scaling factor proposed by
   *   the behaviour, used to increase the time step if the integration
   *   succeeded, or to decrease it if the integration failed or if the
   *   results are not reliable. The argument `rdt` is the time step scaling
   *   factor proposed by the calling code (a priori) or the current one (a
   *   posteriori). The interfaces use the `hasAPrioriTimeStepScalingFactor`
   *   and `hasAPosterioriTimeStepScalingFactor` members of the
   *   `MechanicalBehaviourTraits` class to know if the behaviour gives such
   *   a factor: if not, the behaviour may return `1`.
   *
   * Those methods are not virtual: behaviours are never manipulated through
   * a reference to this class and virtual functions are not supported in
   * device code by some parallel programming models (SYCL).
   *
   * \param H, modelling hypothesis
   * \param NumType, numerical type.
   * \param bool use_qt, says if one shall use quantities.
   * \author Thomas Helfer
   * \date 28/07/2006
   */
  template <MechanicalBehaviourBase::BehaviourType btype,
            ModellingHypothesis::Hypothesis H,
            typename NumType,
            bool use_qt>
  struct MechanicalBehaviour : public TangentOperatorTraits<btype>,
                               public MechanicalBehaviourBase {
    //! \brief available tangent operator
    using SMFlag = typename TangentOperatorTraits<btype>::SMFlag;
    /*!
     * \brief dimension of the space for the the given modelling hypothesis
     */
    static constexpr unsigned short N =
        ModellingHypothesisToSpaceDimension<H>::value;
    //! \brief a simple alias
    using TFELTypes = tfel::config::Types<N, NumType, use_qt>;
  };

  /*!
   * \brief concept satisfied by the classes implementing a mechanical
   * behaviour (see the documentation of the `MechanicalBehaviour` class).
   */
  template <typename BehaviourType>
  concept MechanicalBehaviourConcept =
      requires(BehaviourType& b,
               const BehaviourType& cb,
               const typename BehaviourType::SMFlag smflag,
               const typename BehaviourType::SMType smt,
               const typename BehaviourType::TFELTypes::real rdt) {
        {
          b.computePredictionOperator(smflag, smt)
        } -> std::same_as<typename BehaviourType::IntegrationResult>;
        {
          b.integrate(smflag, smt)
        } -> std::same_as<typename BehaviourType::IntegrationResult>;
        {
          cb.getMinimalTimeStepScalingFactor()
        } -> std::convertible_to<typename BehaviourType::TFELTypes::real>;
        {
          cb.computeAPrioriTimeStepScalingFactor(rdt)
        } -> std::convertible_to<
            std::pair<bool, typename BehaviourType::TFELTypes::real>>;
        {
          cb.computeAPosterioriTimeStepScalingFactor(rdt)
        } -> std::convertible_to<
            std::pair<bool, typename BehaviourType::TFELTypes::real>>;
      };

}  // end of namespace tfel::material

#endif /* LIB_TFEL_MATERIAL_MECHANICALBEHAVIOUR_HXX */
