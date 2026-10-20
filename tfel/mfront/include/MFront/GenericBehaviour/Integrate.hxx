/*!
 * \file   mfront/include/MFront/GenericBehaviour/Integrate.hxx
 * \brief
 * \author Thomas Helfer
 * \date 03/07/2018
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_INTEGRATE_HXX
#define LIB_MFRONT_GENERICBEHAVIOUR_INTEGRATE_HXX

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include "TFEL/Raise.hxx"
#include "TFEL/Material/MechanicalBehaviour.hxx"
#include "TFEL/Math/General/IEEE754.hxx"
#include "TFEL/Math/t2tot2.hxx"
#include "TFEL/Math/t2tost2.hxx"
#include "TFEL/Math/st2tost2.hxx"
#include "TFEL/Material/OutOfBoundsPolicy.hxx"
#include "TFEL/Material/MechanicalBehaviourTraits.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include "MFront/GenericBehaviour/Types.hxx"
#include "MFront/GenericBehaviour/RotateTangentOperatorBlock.hxx"
#include "MFront/GenericBehaviour/BehaviourData.h"

namespace mfront::gb {

  /*!
   * \brief enumeration of the stress measure supported by the generic behaviour
   * interface for finite strain behaviours
   */
  enum struct StressMeasure {
    PK1,                    //!< \brief First Piola-Kirchhoff stress
    PK2,                    //!< \brief Second Piola-Kirchhoff stress
    CAUCHY,                 //!< \brief Cauchy stress
    INVALID_STRESS_MEASURE  //!< \brief invalid value
  };
  /*!
   * \brief offset added to the encoded stress measure (second value of the
   * `K` array of the `mfront_gb_BehaviourData` structure) when the gradients
   * given to a finite strain behaviour are the displacement gradient
   * \(H = F - I\) instead of the deformation gradient \(F\).
   *
   * The deformation gradient being close to the identity, its components
   * are only known within the precision of the numerical type near 1 (about
   * 6e-8 in single precision): the displacement gradient keeps the relative
   * precision of small strains. The strain measures (Green-Lagrange, Hencky)
   * are then computed from \(H\) without cancellation.
   *
   * A library generated before this option reports an invalid stress
   * measure for an offset value.
   */
  inline constexpr int displacement_gradient_stress_measure_offset = 10;

  /*!
   * \return true if the given value (second value of the `K` array of the
   * `mfront_gb_BehaviourData` structure) states that the gradients of a
   * finite strain behaviour are the displacement gradient
   * \param[in] v: encoded value
   */
  TFEL_HOST_DEVICE constexpr bool decodeDisplacementGradient(
      const real v) noexcept {
    return v > real{displacement_gradient_stress_measure_offset} - real{0.5};
  }  // end of decodeDisplacementGradient

  /*!
   * \return the stress measure encoded by the given value (second value of
   * the `K` array of the `mfront_gb_BehaviourData` structure), whether it is
   * offset by `displacement_gradient_stress_measure_offset` or not
   * \param[in] ev: encoded value
   */
  TFEL_HOST_DEVICE constexpr StressMeasure decodeStressMeasure(
      const real ev) noexcept {
    const auto v =
        decodeDisplacementGradient(ev)
            ? ev - real{displacement_gradient_stress_measure_offset}
            : ev;
    if (v < real{0.5}) {
      return StressMeasure::CAUCHY;
    } else if (v < real{1.5}) {
      return StressMeasure::PK2;
    } else if (v < real{2.5}) {
      return StressMeasure::PK1;
    }
    return StressMeasure::INVALID_STRESS_MEASURE;
  }  // end of decodeStressMeasure

  /*!
   * \return the stress measure selected by the caller of the behaviour for
   * finite strain behaviours.
   * \param[in] K: initial values of the consistent tangent operator.
   */
  inline StressMeasure getStressMeasure(const real* const K) {
    return decodeStressMeasure(K[1]);
  }  // end of getStressMeasure

  //! \brief a simple alias
  using FiniteStrainTangentOperator =
      tfel::material::FiniteStrainBehaviourTangentOperatorBase::Flag;

  /*!
   * \return the tangent operator encoded by the given values
   * \param[in] type: type of computation (first value of the `K` array of the
   * `mfront_gb_BehaviourData` structure)
   * \param[in] v: encoded tangent operator (third value of the `K` array)
   * \note FiniteStrainTangentOperator::C_TRUESDELL is returned to signal an
   * error.
   */
  TFEL_HOST_DEVICE constexpr FiniteStrainTangentOperator
  decodeFiniteStrainTangentOperator(const real type, const real v) noexcept {
    if ((type > real{-0.5}) && (type < real{0.5})) {
      // no stiffness requested, returned value is meaningless
      return FiniteStrainTangentOperator::DSIG_DF;
    }
    if (v < real{0.5}) {
      return FiniteStrainTangentOperator::DSIG_DF;
    } else if (v < real{1.5}) {
      return FiniteStrainTangentOperator::DS_DEGL;
    } else if (v < real{2.5}) {
      return FiniteStrainTangentOperator::DPK1_DF;
    } else if (v < real{3.5}) {
      return FiniteStrainTangentOperator::DTAU_DDF;
    }
    return FiniteStrainTangentOperator::C_TRUESDELL;
  }  // end of decodeFiniteStrainTangentOperator

  /*!
   * \return the tangent operator selected by the caller of the behaviour for
   * finite strain behaviours.
   * \param[in] K: initial values of the consistent tangent operator.
   * \note FiniteStrainTangentOperator::C_TRUESDELL is returned to signal an
   * error.
   */
  inline FiniteStrainTangentOperator getTangentOperator(const real* const K) {
    return decodeFiniteStrainTangentOperator(K[0], K[2]);
  }  // end of getTangentOperator

  /*!
   * \brief options of a computation, encoded by the first values of the `K`
   * array of the `mfront_gb_BehaviourData` structure:
   *
   * - the first value gives the type of computation (prediction operator,
   *   integration with or without tangent operator, type of operator). A
   *   value greater than 50 (offset of 100) also requests the speed of
   *   sound;
   * - the second and third values give the stress measure and the tangent
   *   operator of finite strain behaviours. The stress measure is offset by
   *   `displacement_gradient_stress_measure_offset` when the gradients are
   *   the displacement gradient.
   */
  struct ComputationOptions {
    /*!
     * \brief constructor from the type of computation only, used for the
     * behaviours which are not finite strain behaviours
     * \param[in] o0: first value
     */
    TFEL_HOST_DEVICE constexpr explicit ComputationOptions(
        const real o0) noexcept
        : ComputationOptions(o0, real{0}, real{0}) {}
    /*!
     * \brief constructor
     * \param[in] o0: first value
     * \param[in] o1: second value
     * \param[in] o2: third value
     */
    TFEL_HOST_DEVICE constexpr ComputationOptions(const real o0,
                                                  const real o1,
                                                  const real o2) noexcept
        : speed_of_sound(o0 > real{50}),
          type(o0 > real{50} ? o0 - real{100} : o0),
          stress_measure(decodeStressMeasure(o1)),
          displacement_gradient(decodeDisplacementGradient(o1)),
          tangent_operator(decodeFiniteStrainTangentOperator(
              o0 > real{50} ? o0 - real{100} : o0, o2)) {}
    //! \return true if only the prediction operator is requested
    TFEL_HOST_DEVICE constexpr bool isPredictionRequested() const noexcept {
      return this->type < real{-0.25};
    }
    //! \return true if a tangent operator is computed after integration
    TFEL_HOST_DEVICE constexpr bool isTangentOperatorRequested()
        const noexcept {
      return this->type > real{0.5};
    }
    //! \return the type of the prediction operator
    template <typename Behaviour>
    TFEL_HOST_DEVICE constexpr typename Behaviour::SMType
    getPredictionOperatorType() const noexcept {
      if (real{-1.5} < this->type) {
        return Behaviour::ELASTIC;
      } else if (real{-2.5} < this->type) {
        return Behaviour::SECANTOPERATOR;
      }
      return Behaviour::TANGENTOPERATOR;
    }
    //! \return the type of the operator computed after integration
    template <typename Behaviour>
    TFEL_HOST_DEVICE constexpr typename Behaviour::SMType
    getTangentOperatorType() const noexcept {
      if (this->type < real{0.5}) {
        return Behaviour::NOSTIFFNESSREQUESTED;
      } else if (this->type < real{1.5}) {
        return Behaviour::ELASTIC;
      } else if (this->type < real{2.5}) {
        return Behaviour::SECANTOPERATOR;
      } else if (this->type < real{3.5}) {
        return Behaviour::TANGENTOPERATOR;
      }
      return Behaviour::CONSISTENTTANGENTOPERATOR;
    }
    //! \brief true if the speed of sound is requested
    bool speed_of_sound;
    //! \brief type of computation (without the offset of the speed of sound)
    real type;
    //! \brief stress measure (finite strain behaviours)
    StressMeasure stress_measure;
    /*!
     * \brief true if the gradients are the displacement gradient \(H\)
     * rather than the deformation gradient \(F = I + H\) (finite strain
     * behaviours)
     */
    bool displacement_gradient;
    //! \brief tangent operator (finite strain behaviours)
    FiniteStrainTangentOperator tangent_operator;
  };

  /*!
   * \brief export the tangent operator used by some generic behaviours.
   * \tparam StressType: numeric type used
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v, const StressType K) {
    tfel::math::map<StressType>(v) = K;
  }  // end of exportTangentOperator

  /*!
   * \brief export the tangent operator used by some generic behaviours.
   * \tparam StressType: numeric type used
   * \tparam N: size of the array containing the tangent operator
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType, unsigned short N>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v,
      const tfel::math::tvector<N, StressType>& K) {
    tfel::math::map<tfel::math::tvector<N, StressType>>(v) = K;
  }  // end of exportTangentOperator

  /*!
   * \brief export the tangent operator used by cohesive zone models.
   * \tparam StressType: numeric type used
   * \tparam N: space dimension
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType, unsigned short N>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v,
      const tfel::math::tmatrix<N, N, StressType>& K) {
    tfel::math::map<tfel::math::tmatrix<N, N, StressType>>(v) = K;
  }  // end of exportTangentOperator

  /*!
   * \brief export the tangent operator used by some generic behaviours.
   * \tparam StressType: numeric type used
   * \tparam N: space dimension
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType, unsigned short N>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v,
      const tfel::math::t2tost2<N, StressType>& K) {
    tfel::math::map<tfel::math::t2tost2<N, StressType>>(v) = K;
  }  // end of exportTangentOperator

  /*!
   * \brief export the tangent operator used by some generic behaviours.
   * \tparam StressType: numeric type used
   * \tparam N: space dimension
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType, unsigned short N>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v,
      const tfel::math::t2tot2<N, StressType>& K) {
    tfel::math::map<tfel::math::t2tot2<N, StressType>>(v) = K;
  }  // end of exportTangentOperator

  /*!
   * \brief export the tangent operator used by standard small
   * strain behaviours.
   * \tparam StressType: numeric type used
   * \tparam N: space dimension
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType, unsigned short N>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v,
      const tfel::math::st2tost2<N, StressType>& K) {
    tfel::math::map<tfel::math::st2tost2<N, StressType>>(v) = K;
  }  // end of exportTangentOperator

  /*!
   * \brief export the tangent operator used by finite strain
   * strain behaviours.
   * \tparam StressType: numeric type used
   * \tparam N: space dimension
   * \param[out] v: exported values
   * \param[in] K: computed tangent operator
   */
  template <typename StressType, unsigned short N>
  TFEL_HOST_DEVICE void exportTangentOperator(
      tfel::math::base_type<StressType>* const v,
      const tfel::material::FiniteStrainBehaviourTangentOperator<N, StressType>&
          K) {
    if (K.template is<tfel::math::t2tot2<N, StressType>>()) {
      const auto& k = K.template get<tfel::math::t2tot2<N, StressType>>();
      tfel::math::map<tfel::math::t2tot2<N, StressType>>(v) = k;
    } else if (K.template is<tfel::math::t2tot2<N, StressType>*>()) {
      const auto& k = *(K.template get<tfel::math::t2tot2<N, StressType>*>());
      tfel::math::map<tfel::math::t2tot2<N, StressType>>(v) = k;
    } else if (K.template is<tfel::math::t2tost2<N, StressType>>()) {
      const auto& k = K.template get<tfel::math::t2tost2<N, StressType>>();
      tfel::math::map<tfel::math::t2tost2<N, StressType>>(v) = k;
    } else if (K.template is<tfel::math::t2tost2<N, StressType>*>()) {
      const auto& k = *(K.template get<tfel::math::t2tost2<N, StressType>*>());
      tfel::math::map<tfel::math::t2tost2<N, StressType>>(v) = k;
    } else if (K.template is<tfel::math::st2tost2<N, StressType>>()) {
      const auto& k = K.template get<tfel::math::st2tost2<N, StressType>>();
      tfel::math::map<tfel::math::st2tost2<N, StressType>>(v) = k;
    } else if (K.template is<tfel::math::st2tost2<N, StressType>*>()) {
      const auto& k = *(K.template get<tfel::math::st2tost2<N, StressType>*>());
      tfel::math::map<tfel::math::st2tost2<N, StressType>>(v) = k;
    } else {
#ifndef TFEL_DEVICE_CODE
      tfel::raise(
          "mfront::gb::exportTangentOperator: "
          "unsupported tangent operator type");
#endif /* TFEL_DEVICE_CODE */
    }
  }  // end of exportTangentOperator

  /*!
   * \brief report an error. The message is truncated if needed.
   * \param[in] d: behaviour data
   * \param[in] e: error message
   */
  TFEL_HOST_DEVICE inline void reportError(mfront_gb_BehaviourData& d,
                                           const char* const e) noexcept {
    constexpr std::size_t bsize = 511;
    if (d.error_message == nullptr) {
      return;
    }
    auto n = std::size_t{};
    for (; (n != bsize) && (e[n] != '\0'); ++n) {
      d.error_message[n] = e[n];
    }
    d.error_message[n] = '\0';
  }  // end of reportError

  /*!
   * \brief a simple function to report integration failure when
   * an exception is thrown.
   * \param[in] d: behaviour data
   */
  inline void reportFailureByException(mfront_gb_BehaviourData& d) {
    try {
      throw;
    } catch (std::exception& e) {
      reportError(d, e.what());
    } catch (...) {
      reportError(d, "unknown exception");
    }
  }  // end of reportFailureByException

  /*!
   * \brief call the given function. On the host, the exceptions thrown by the
   * function are reported as errors.
   * \return the value returned by the function, or -1 if an exception has
   * been thrown
   * \param[in,out] d: behaviour data
   * \param[in] f: function
   * \param[in] on_exception: function called if an exception is thrown
   */
  template <typename Function, typename ExceptionHandler>
  TFEL_HOST_DEVICE int callAndReportExceptions(
      mfront_gb_BehaviourData& d,
      const Function& f,
      [[maybe_unused]] const ExceptionHandler& on_exception) {
#ifdef TFEL_DEVICE_CODE
    static_cast<void>(d);
    return f();
#else  /* TFEL_DEVICE_CODE */
    try {
      return f();
    } catch (...) {
      reportFailureByException(d);
      on_exception();
    }
    return -1;
#endif /* TFEL_DEVICE_CODE */
  }  // end of callAndReportExceptions

  /*!
   * \brief call the given function. On the host, the exceptions thrown by the
   * function are reported as errors.
   * \return the value returned by the function, or -1 if an exception has
   * been thrown
   * \param[in,out] d: behaviour data
   * \param[in] f: function
   */
  template <typename Function>
  TFEL_HOST_DEVICE int callAndReportExceptions(mfront_gb_BehaviourData& d,
                                               const Function& f) {
    return callAndReportExceptions(d, f, [] {});
  }  // end of callAndReportExceptions

  /*!
   * \brief compute the prediction operator
   * \tparam Behaviour: class implementing the behaviour
   * \param[in,out] b: behaviour
   * \param[in,out] d: behaviour data
   * \param[in] o: computation options
   * \param[in] f: choice of consistent tangent operator
   */
  template <typename Behaviour>
  TFEL_HOST_DEVICE int computePredictionOperator(
      Behaviour& b,
      mfront_gb_BehaviourData& d,
      const ComputationOptions& o,
      const typename Behaviour::SMFlag f) {
    const auto smt = o.template getPredictionOperatorType<Behaviour>();
    if (!b.computePredictionOperator(f, smt)) {
      return -1;
    }
    exportTangentOperator(d.K, b.getTangentOperator());
    return 1;
  }  // end of computePredictionOperator

  //! structure in charge of calling the computeInternalEnergy method
  struct InternalEnergyComputer {
    /*!
     * \brief call the `computeInternalEnergy` method
     * \param[in] d: behaviour data
     * \param[in] b: behaviour
     */
    template <typename Behaviour>
    TFEL_HOST_DEVICE static void exe(mfront_gb_BehaviourData& d,
                                     const Behaviour& b) {
      using stress = typename Behaviour::stress;
      auto e = stress{*(d.s0.stored_energy)};
      b.computeInternalEnergy(e);
      tfel::math::map<stress>(*(d.s1.stored_energy)) = e;
    }
  };  // end of struct InternalEnergyComputer
  //! structure in charge of calling the computeDissipatedEnergy method
  struct DissipatedEnergyComputer {
    /*!
     * \brief call the `computeDissipatedEnergy` method
     * \param[in] d: behaviour data
     * \param[in] b: behaviour
     */
    template <typename Behaviour>
    TFEL_HOST_DEVICE static void exe(mfront_gb_BehaviourData& d,
                                     const Behaviour& b) {
      using stress = typename Behaviour::stress;
      auto e = stress{*(d.s0.dissipated_energy)};
      b.computeDissipatedEnergy(e);
      tfel::math::map<stress>(*(d.s1.dissipated_energy)) = e;
    }  // end of exe
  };   // end of struct DissipatedEnergyComputer
  //! place holder for tag dispatching
  struct DoNothingEnergyComputer {
    //! \brief empty function
    template <typename Behaviour>
    TFEL_HOST_DEVICE static void exe(mfront_gb_BehaviourData&,
                                     const Behaviour&) {}
  };  // end of struct DoNothingEnergyComputer

  /*!
   * \brief integrate the behaviour over a time step.
   *
   * This function can be called in device code.
   *
   * \tparam Behaviour: class implementing the behaviour
   * \tparam ConstructorArguments: types of the additional arguments of the
   * constructor of the behaviour
   * \param[in,out] d: behaviour data. The first values of the `K` array are
   * not interpreted as computation options.
   * \param[in] o: computation options
   * \param[in] f: choice of consistent tangent operator
   * \param[in] p: out of bounds policy
   * \param[in] args: additional arguments of the constructor of the
   * behaviour (the parameters of the behaviour for the `generic-parallel`
   * interface)
   */
  template <typename Behaviour, typename... ConstructorArguments>
  TFEL_HOST_DEVICE int integrate(mfront_gb_BehaviourData& d,
                                 const ComputationOptions& o,
                                 const typename Behaviour::SMFlag f,
                                 const tfel::material::OutOfBoundsPolicy p,
                                 const ConstructorArguments&... args) {
    static_assert(tfel::material::MechanicalBehaviourConcept<Behaviour>);
    using MTraits = tfel::material::MechanicalBehaviourTraits<Behaviour>;
    using IEnergyComputer =
        typename std::conditional<MTraits::hasComputeInternalEnergy,
                                  InternalEnergyComputer,
                                  DoNothingEnergyComputer>::type;
    using DEnergyComputer =
        typename std::conditional<MTraits::hasComputeDissipatedEnergy,
                                  DissipatedEnergyComputer,
                                  DoNothingEnergyComputer>::type;
    using speed = typename Behaviour::speed;
    using massdensity = typename Behaviour::massdensity;
    using behaviour_real_type = typename Behaviour::real;
    auto&& rdt = tfel::math::map<behaviour_real_type>(d.rdt);
    // time step scaling factor reported if an exception is thrown. The
    // behaviour is built in the protected region, since its constructor may
    // throw on the host (for instance while building shared data).
    auto rdt_on_exception = behaviour_real_type{rdt};
    const auto integrate = [&]() -> int {
      Behaviour b(d, args...);
      rdt_on_exception = b.getMinimalTimeStepScalingFactor();
      b.setOutOfBoundsPolicy(p);
      if (!b.initialize()) {
        reportError(d, "behaviour initialisation failed");
        return -1;
      }
      if (!b.checkBounds()) {
        // reached when `TFEL_DEVICE_CODE` is defined (device code, and the
        // whole translation unit with the compilers translating host and
        // device code in a single pass): otherwise, out of bounds are
        // reported by exceptions
        reportError(d, "out of bounds");
        rdt = b.getMinimalTimeStepScalingFactor();
        return -1;
      }
      if (o.isPredictionRequested()) {
        if (o.speed_of_sound) {
          tfel::math::map<speed>(d.speed_of_sound) =
              b.computeSpeedOfSound(massdensity(*(d.s0.mass_density)));
        }
        if constexpr (!MTraits::hasPredictionOperator) {
          reportError(d, "prediction operator is not implemented");
          return -1;
        } else {
          return computePredictionOperator(b, d, o, f);
        }
      }
      const auto smt = o.template getTangentOperatorType<Behaviour>();
      if constexpr (!MTraits::hasConsistentTangentOperator) {
        if (smt != Behaviour::NOSTIFFNESSREQUESTED) {
          reportError(d, "tangent operator is not implemented");
          return -1;
        }
      }
      const auto tsf = b.computeAPrioriTimeStepScalingFactor(rdt);
      rdt = tsf.second;
      if (!tsf.first) {
        return -1;
      }
      if (b.integrate(f, smt) == Behaviour::FAILURE) {
        rdt = b.getMinimalTimeStepScalingFactor();
        return -1;
      }
      const auto atsf = b.computeAPosterioriTimeStepScalingFactor(rdt);
      if (rdt > atsf.second) {
        rdt = atsf.second;
      }
      if (!atsf.first) {
        return -1;
      }
      b.exportStateData(d.s1);
      if constexpr (MTraits::hasConsistentTangentOperator) {
        if (o.isTangentOperatorRequested()) {
          exportTangentOperator(d.K, b.getTangentOperator());
        }
      }
      IEnergyComputer::exe(d, b);
      DEnergyComputer::exe(d, b);
      if (o.speed_of_sound) {
        tfel::math::map<speed>(d.speed_of_sound) =
            b.computeSpeedOfSound(massdensity(*(d.s1.mass_density)));
      }
      return rdt < behaviour_real_type{0.99} ? 0 : 1;
    };
    return callAndReportExceptions(
        d, integrate, [&rdt, &rdt_on_exception] { rdt = rdt_on_exception; });
  }  // end of integrate

  /*!
   * \brief integrate the behaviour over a time step
   * \tparam Behaviour: class implementing the behaviour
   * \param[in,out] d: behaviour data. The computation options are given by
   * the first value of the `K` array.
   * \param[in] f: choice of consistent tangent operator
   * \param[in] p: out of bounds policy
   */
  template <typename Behaviour>
  int integrate(mfront_gb_BehaviourData& d,
                const typename Behaviour::SMFlag f,
                const tfel::material::OutOfBoundsPolicy p) {
    return integrate<Behaviour>(d, ComputationOptions{d.K[0]}, f, p);
  }  // end of integrate

  /*!
   * \brief function object calling a method of a behaviour through a
   * pointer to member.
   *
   * \note calls through pointers to members are not supported in device
   * code by some compilers (`nvc++` 25.11 silently ignores them). The code
   * generated for devices uses function objects calling the methods by their
   * names.
   *
   * \tparam m: pointer to the method
   */
  template <auto m>
  struct MethodCaller {
    //! \brief call the method of the given behaviour
    template <typename Behaviour, typename... Arguments>
    TFEL_HOST_DEVICE void operator()(Behaviour& b, Arguments&&... args) const {
      (b.*m)(std::forward<Arguments>(args)...);
    }
  };

  /*!
   * \brief execute the given initialize function.
   *
   * This function can be called in device code.
   *
   * \tparam Behaviour: class describing the initialize function.
   * \tparam InitializeFunction: type of the function object calling the
   * method implementing the initialize function.
   * \tparam ConstructorArguments: types of the additional arguments of the
   * constructor of the behaviour
   * \param[in] f: function object calling the method implementing the
   * initialize function, as `f(b, initialize_variables)` where `b` is the
   * behaviour.
   * \param[in,out] d: behaviour data.
   * \param[in] initialize_variables: pointer to the values of the
   * initialize function variables.
   * \param[in] p: out of bounds policy.
   * \param[in] args: additional arguments of the constructor of the
   * behaviour
   */
  template <typename Behaviour,
            typename InitializeFunction,
            typename... ConstructorArguments>
  TFEL_HOST_DEVICE int executeInitializeFunction(
      const InitializeFunction& f,
      mfront_gb_BehaviourData& d,
      const real* const initialize_variables,
      const tfel::material::OutOfBoundsPolicy p,
      const ConstructorArguments&... args) {
    // the behaviour is built from the state at the beginning of the time
    // step: a copy of the behaviour data, whose gradients, material
    // properties and external state variables at the end of the time step
    // are the ones at the beginning of the time step
    auto d0 = d;
    d0.s1.gradients = d.s0.gradients;
    d0.s1.material_properties = d.s0.material_properties;
    d0.s1.external_state_variables = d.s0.external_state_variables;
    const auto execute = [&]() -> int {
      Behaviour b(d0, args...);
      b.setOutOfBoundsPolicy(p);
      if (!b.initialize()) {
        reportError(d, "behaviour initialisation failed");
        return -1;
      }
      f(b, initialize_variables);
      b.exportStateData(d.s1);
      return EXIT_SUCCESS;
    };
    return callAndReportExceptions(d, execute);
  }  // end of executeInitializeFunction

  /*!
   * \brief execute the given initialize function.
   * \tparam Behaviour: class describing the initialize function.
   * \tparam m: method implementing the initialize function.
   * \param[in,out] d: behaviour data.
   * \param[in] initialize_variables: pointer to the values of the
   * initialize function variables.
   * \param[in] p: out of bounds policy.
   */
  template <typename Behaviour, void (Behaviour::*m)(const real* const)>
  int executeInitializeFunction(mfront_gb_BehaviourData& d,
                                const real* const initialize_variables,
                                const tfel::material::OutOfBoundsPolicy p) {
    return executeInitializeFunction<Behaviour>(MethodCaller<m>{}, d,
                                                initialize_variables, p);
  }  // end of executeInitializeFunction

  /*!
   * \brief execute the given post-processing.
   *
   * The behaviour is initialized with the state at the end of the time step
   * (thermodynamic forces, internal state variables, gradients and external
   * state variables).
   *
   * This function can be called in device code.
   *
   * \tparam Behaviour: class describing the post-processing.
   * \tparam use_initial_state: boolean stating if the postprocessing uses
   * the initial state of the material.
   * \tparam PostProcessing: type of the function object calling the method
   * implementing the post-processing.
   * \tparam ConstructorArguments: types of the additional arguments of the
   * constructor of the behaviour
   * \param[in] f: function object calling the method implementing the
   * post-processing, as `f(b, post_processing_variables, initial_state)`
   * where `b` is the behaviour.
   * \param[out] post_processing_variables: pointer to the values of the
   * post-processing variables.
   * \param[in,out] d: behaviour data.
   * \param[in] p: out of bounds policy.
   * \param[in] args: additional arguments of the constructor of the
   * behaviour
   */
  template <typename Behaviour,
            bool use_initial_state,
            typename PostProcessing,
            typename... ConstructorArguments>
  TFEL_HOST_DEVICE int executePostProcessing(
      const PostProcessing& f,
      real* const post_processing_variables,
      mfront_gb_BehaviourData& d,
      const tfel::material::OutOfBoundsPolicy p,
      const ConstructorArguments&... args) {
    // the behaviour is built from the thermodynamic forces and the internal
    // state variables at the end of the time step: a copy of the behaviour
    // data, whose thermodynamic forces and internal state variables at the
    // beginning of the time step are the ones at the end of the time step
    auto d1 = d;
    d1.s0.thermodynamic_forces = d.s1.thermodynamic_forces;
    d1.s0.internal_state_variables = d.s1.internal_state_variables;
    const auto execute = [&]() -> int {
      Behaviour b(d1, args...);
      b.setOutOfBoundsPolicy(p);
      // This call sets the values of the gradients and external state
      // variables at the end of the time step
      b.updateExternalStateVariables();
      if (!b.initialize()) {
        reportError(d, "behaviour initialisation failed");
        return -1;
      }
      if constexpr (use_initial_state) {
        // create an object containing the intial state
        // Here, we can't initialize an `Behaviour::BehaviourData` as it would
        // not initialize the gradients
        Behaviour initial_state(d, args...);
        f(b, post_processing_variables, initial_state);
      } else {
        // the initial state is not used
        f(b, post_processing_variables, b);
      }
      return EXIT_SUCCESS;
    };
    return callAndReportExceptions(d, execute);
  }  // end of executePostProcessing

  /*!
   * \brief execute the given post-processing.
   * \tparam Behaviour: class describing the post-processing.
   * \tparam m: method implementing the post-processing.
   * \tparam use_initial_state: boolean stating if the postprocessing uses
   * the initial state of the material.
   * \param[out] post_processing_variables: pointer to the values of the
   * post-processing variables.
   * \param[in,out] d: behaviour data.
   * \param[in] p: out of bounds policy.
   */
  template <typename Behaviour,
            void (Behaviour::*m)(real* const,
                                 const typename Behaviour::BehaviourData&),
            const bool use_initial_state>
  int executePostProcessing(real* const post_processing_variables,
                            mfront_gb_BehaviourData& d,
                            const tfel::material::OutOfBoundsPolicy p) {
    return executePostProcessing<Behaviour, use_initial_state>(
        MethodCaller<m>{}, post_processing_variables, d, p);
  }  // end of executePostProcessing

}  // end of namespace mfront::gb

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_INTEGRATE_HXX */
