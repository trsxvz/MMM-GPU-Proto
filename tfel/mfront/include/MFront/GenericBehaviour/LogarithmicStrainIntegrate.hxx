/*!
 * \file   mfront/include/MFront/GenericBehaviour/LogarithmicStrainIntegrate.hxx
 * \brief  This file declares the functions integrating the small strain
 * behaviours written in the Hencky (logarithmic) strain framework, shared
 * by the `generic` and `generic-parallel` interfaces.
 * \author Thomas Helfer
 * \date   31/01/2022
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_LOGARITHMICSTRAININTEGRATE_HXX
#define LIB_MFRONT_GENERICBEHAVIOUR_LOGARITHMICSTRAININTEGRATE_HXX

#include <cmath>
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/st2tost2.hxx"
#include "TFEL/FSAlgorithm/copy.hxx"
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include "TFEL/Material/LogarithmicStrainHandler.hxx"
#include "MFront/GenericBehaviour/GenericBehaviourTraits.hxx"
#include "MFront/GenericBehaviour/Integrate.hxx"
#include "MFront/GenericBehaviour/StressMeasures.hxx"

namespace mfront::gb::logarithmic_strain {

  namespace internals {

    //! \brief a simple alias
    template <unsigned short N>
    using LogarithmicStrainHandler =
        tfel::material::LogarithmicStrainHandler<N, real>;

    /*!
     * \brief copy the deformation gradient from the gradients
     * \tparam Behaviour: class implementing the behaviour
     * \param[out] F: deformation gradient
     * \param[in] gradients: gradients
     */
    template <typename Behaviour, unsigned short N>
    TFEL_HOST_DEVICE void getDeformationGradient(tfel::math::tensor<N, real>& F,
                                                 const real* const gradients) {
      using namespace tfel::material;
      using Traits = GenericBehaviourTraits<Behaviour>;
      constexpr auto TensorSize =
          ModellingHypothesisToTensorSize<Traits::hypothesis>::value;
      tfel::fsalgo::copy<TensorSize>::exe(gradients, F.begin());
    }  // end of getDeformationGradient

    /*!
     * \brief in plane stress, update the axial component of the deformation
     * gradient from the axial strain, which is computed by the behaviour and
     * stored in the internal state variables.
     * \return false if the axial strain is not computed by the behaviour
     * \tparam Behaviour: class implementing the behaviour
     * \param[in,out] lgh: logarithmic strain handler
     * \param[in,out] d: behaviour data
     * \param[in] isvs: internal state variables
     */
    template <typename Behaviour, unsigned short N>
    TFEL_HOST_DEVICE bool updateAxialDeformationGradient(
        [[maybe_unused]] LogarithmicStrainHandler<N>& lgh,
        [[maybe_unused]] mfront_gb_BehaviourData& d,
        [[maybe_unused]] const real* const isvs) {
      using namespace tfel::material;
      using Traits = GenericBehaviourTraits<Behaviour>;
      constexpr auto hypothesis = Traits::hypothesis;
      if constexpr ((hypothesis == ModellingHypothesis::PLANESTRESS) ||
                    (hypothesis == ModellingHypothesis::
                                       AXISYMMETRICALGENERALISEDPLANESTRESS)) {
        if constexpr (Traits::has_axial_strain_offset) {
          lgh.updateAxialDeformationGradient(
              std::exp(isvs[Traits::axial_strain_offset]));
        } else {
          reportError(d,
                      "The behaviour does not declare the axial strain"
                      "as an internal state variable");
          return false;
        }
      }
      return true;
    }  // end of updateAxialDeformationGradient

    /*!
     * \brief convert the given stress to the dual stress of the Hencky strain
     * \return false if the stress measure is invalid
     * \param[out] T: dual stress of the Hencky strain
     * \param[out] s: Cauchy stress
     * \param[in,out] d: behaviour data
     * \param[in] v: values of the stress
     * \param[in] F: deformation gradient
     * \param[in] lgh: logarithmic strain handler
     * \param[in] sm: stress measure
     */
    template <unsigned short N>
    TFEL_HOST_DEVICE bool convertToDualStress(
        tfel::math::stensor<N, real>& T,
        tfel::math::stensor<N, real>& s,
        mfront_gb_BehaviourData& d,
        const real* const v,
        const tfel::math::tensor<N, real>& F,
        const LogarithmicStrainHandler<N>& lgh,
        const StressMeasure sm) {
      if (sm == StressMeasure::CAUCHY) {
        s = tfel::math::ConstStensorView<N, real>(v);
      } else if (sm == StressMeasure::PK1) {
        s = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(
            tfel::math::ConstTensorView<N, real>(v), F);
      } else if (sm == StressMeasure::PK2) {
        s = tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(
            tfel::math::ConstStensorView<N, real>(v), F);
      } else {
        reportError(d, "invalid choice for the stress measure");
        return false;
      }
      T = lgh.convertFromCauchyStress(s);
      return true;
    }  // end of convertToDualStress

    /*!
     * \brief convert the tangent operator computed by the behaviour (the
     * derivative of the dual stress of the Hencky strain with respect to the
     * Hencky strain) to the requested tangent operator
     * \return false if the tangent operator is invalid
     * \param[out] v: values of the tangent operator
     * \param[in,out] d: behaviour data
     * \param[in] K: tangent operator computed by the behaviour
     * \param[in] lgh: logarithmic strain handler associated with the
     * deformation gradient `F`
     * \param[in] T: dual stress of the Hencky strain
     * \param[in] F0: deformation gradient at the beginning of the time step
     * \param[in] F: current deformation gradient
     * \param[in] s: Cauchy stress
     * \param[in] smf: requested tangent operator
     */
    template <unsigned short N>
    TFEL_HOST_DEVICE bool convertTangentOperator(
        real* const v,
        mfront_gb_BehaviourData& d,
        const tfel::math::st2tost2<N, real>& K,
        const LogarithmicStrainHandler<N>& lgh,
        const tfel::math::stensor<N, real>& T,
        const tfel::math::tensor<N, real>& F0,
        const tfel::math::tensor<N, real>& F,
        const tfel::math::stensor<N, real>& s,
        const FiniteStrainTangentOperator smf) {
      using tfel::material::convert;
      if ((smf == FiniteStrainTangentOperator::DSIG_DF) ||
          (smf == FiniteStrainTangentOperator::DTAU_DDF)) {
        const auto Cs = lgh.convertToSpatialTangentModuli(K, T);
        const auto Dt =
            convert<FiniteStrainTangentOperator::DTAU_DF,
                    FiniteStrainTangentOperator::SPATIAL_MODULI>(Cs, F0, F, s);
        if (smf == FiniteStrainTangentOperator::DSIG_DF) {
          tfel::math::map<tfel::math::t2tost2<N, real>>(v) =
              convert<FiniteStrainTangentOperator::DSIG_DF,
                      FiniteStrainTangentOperator::DTAU_DF>(Dt, F0, F, s);
        } else {
          tfel::math::map<tfel::math::t2tost2<N, real>>(v) =
              convert<FiniteStrainTangentOperator::DTAU_DDF,
                      FiniteStrainTangentOperator::DTAU_DF>(Dt, F0, F, s);
        }
      } else if (smf == FiniteStrainTangentOperator::DS_DEGL) {
        tfel::math::map<tfel::math::st2tost2<N, real>>(v) =
            lgh.convertToMaterialTangentModuli(K, T);
      } else if (smf == FiniteStrainTangentOperator::DPK1_DF) {
        const auto Cse = lgh.convertToMaterialTangentModuli(K, T);
        tfel::math::map<tfel::math::t2tot2<N, real>>(v) =
            convert<FiniteStrainTangentOperator::DPK1_DF,
                    FiniteStrainTangentOperator::DS_DEGL>(Cse, F0, F, s);
      } else {
        reportError(d, "invalid choice for consistent tangent operator");
        return false;
      }
      return true;
    }  // end of convertTangentOperator

    /*!
     * \return the setting of the logarithmic strain handler associated with
     * the given tangent operator
     * \param[in] smf: tangent operator
     */
    template <unsigned short N>
    TFEL_HOST_DEVICE constexpr auto getSetting(
        const FiniteStrainTangentOperator smf) noexcept {
      return (smf == FiniteStrainTangentOperator::DSIG_DF)
                 ? LogarithmicStrainHandler<N>::EULERIAN
                 : LogarithmicStrainHandler<N>::LAGRANGIAN;
    }  // end of getSetting

  }  // end of namespace internals

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
   * initialize function (see `mfront::gb::executeInitializeFunction`).
   * \param[in,out] d: behaviour data.
   * \param[in] initialize_variables: pointer to the values of the
   * initialize function variables.
   * \param[in] o: computation options (only the stress measure is used).
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
      const ComputationOptions& o,
      const tfel::material::OutOfBoundsPolicy p,
      const ConstructorArguments&... args) {
    using namespace tfel::material;
    using Traits = GenericBehaviourTraits<Behaviour>;
    constexpr auto N =
        ModellingHypothesisToSpaceDimension<Traits::hypothesis>::value;
    using LogarithmicStrainHandler = internals::LogarithmicStrainHandler<N>;
    auto F0 = tfel::math::tensor<N, real>{};
    internals::getDeformationGradient<Behaviour>(F0, d.s0.gradients);
    auto lgh0 =
        LogarithmicStrainHandler(LogarithmicStrainHandler::LAGRANGIAN, F0);
    auto e0 = lgh0.getHenckyLogarithmicStrain();
    auto s0 = tfel::math::stensor<N, real>{};
    auto T0 = tfel::math::stensor<N, real>{};
    if ((!internals::updateAxialDeformationGradient<Behaviour>(
            lgh0, d, d.s0.internal_state_variables)) ||
        (!internals::convertToDualStress(T0, s0, d, d.s0.thermodynamic_forces,
                                         F0, lgh0, o.stress_measure))) {
      return -1;
    }
    // the initialize function uses the state at the beginning of the time
    // step: the dual stress it computes is converted with the logarithmic
    // strain handler at the beginning of the time step
    auto T1 = T0;
    auto dc = d;
    dc.s0.gradients = e0.begin();
    dc.s1.gradients = e0.begin();
    dc.s0.thermodynamic_forces = T0.begin();
    dc.s1.thermodynamic_forces = T1.begin();
    const auto r = mfront::gb::executeInitializeFunction<Behaviour>(
        f, dc, initialize_variables, p, args...);
    if (r == -1) {
      return r;
    }
    const auto s1 = lgh0.convertToCauchyStress(T1);
    if (!mfront::gb::internals::convertFromCauchyStress(
            d.s1.thermodynamic_forces, d, s1, F0, o.stress_measure)) {
      return -1;
    }
    return r;
  }  // end of executeInitializeFunction

  /*!
   * \brief execute the given initialize function
   * \tparam Behaviour: class describing the initialize function.
   * \tparam m: method implementing the initialize function.
   * \param[in,out] d: behaviour data. The stress measure is given by the
   * second value of the `K` array.
   * \param[in] initialize_variables: pointer to the values of the
   * initialize function variables.
   * \param[in] p: out of bounds policy.
   */
  template <typename Behaviour, void (Behaviour::*m)(const real* const)>
  int executeInitializeFunction(mfront_gb_BehaviourData& d,
                                const real* const initialize_variables,
                                const tfel::material::OutOfBoundsPolicy p) {
    return executeInitializeFunction<Behaviour>(
        MethodCaller<m>{}, d, initialize_variables,
        ComputationOptions{d.K[0], d.K[1], d.K[2]}, p);
  }  // end of executeInitializeFunction

  /*!
   * \brief integrate the behaviour using the Hencky strain measure.
   *
   * This function can be called in device code.
   *
   * \tparam Behaviour: class implementing the behaviour
   * \tparam ConstructorArguments: types of the additional arguments of the
   * constructor of the behaviour
   * \param[in,out] d: behaviour data
   * \param[in] o: computation options
   * \param[in] p: out of bounds policy
   * \param[in] args: additional arguments of the constructor of the
   * behaviour
   */
  template <typename Behaviour, typename... ConstructorArguments>
  TFEL_HOST_DEVICE int integrate(mfront_gb_BehaviourData& d,
                                 const ComputationOptions& o,
                                 const tfel::material::OutOfBoundsPolicy p,
                                 const ConstructorArguments&... args) {
    using namespace tfel::material;
    using Traits = GenericBehaviourTraits<Behaviour>;
    constexpr auto N =
        ModellingHypothesisToSpaceDimension<Traits::hypothesis>::value;
    using LogarithmicStrainHandler = internals::LogarithmicStrainHandler<N>;
    const auto sm = o.stress_measure;
    if (sm == StressMeasure::INVALID_STRESS_MEASURE) {
      reportError(d, "invalid choice for the stress measure");
      return -1;
    }
    const auto smf = o.tangent_operator;
    if (smf == FiniteStrainTangentOperator::C_TRUESDELL) {
      reportError(d, "invalid choice for consistent tangent operator");
      return -1;
    }
    auto K = tfel::math::st2tost2<N, real>{};
    auto F0 = tfel::math::tensor<N, real>{};
    auto F1 = tfel::math::tensor<N, real>{};
    internals::getDeformationGradient<Behaviour>(F0, d.s0.gradients);
    internals::getDeformationGradient<Behaviour>(F1, d.s1.gradients);
    const auto setting = internals::getSetting<N>(smf);
    auto lgh0 = LogarithmicStrainHandler(setting, F0);
    auto lgh1 = LogarithmicStrainHandler(setting, F1);
    auto e0 = lgh0.getHenckyLogarithmicStrain();
    auto e1 = lgh1.getHenckyLogarithmicStrain();
    auto s0 = tfel::math::stensor<N, real>{};
    auto T0 = tfel::math::stensor<N, real>{};
    auto T1 = tfel::math::stensor<N, real>{};
    if ((!internals::updateAxialDeformationGradient<Behaviour>(
            lgh0, d, d.s0.internal_state_variables)) ||
        (!internals::convertToDualStress(T0, s0, d, d.s0.thermodynamic_forces,
                                         F0, lgh0, sm))) {
      return -1;
    }
    auto dc = d;
    dc.s0.gradients = e0.begin();
    dc.s1.gradients = e1.begin();
    dc.s0.thermodynamic_forces = T0.begin();
    dc.s1.thermodynamic_forces = T1.begin();
    dc.K = K.begin();
    const auto r = mfront::gb::integrate<Behaviour>(
        dc, o, Behaviour::STANDARDTANGENTOPERATOR, p, args...);
    if (r == -1) {
      return r;
    }
    if (o.isPredictionRequested()) {
      if (!internals::convertTangentOperator(d.K, d, K, lgh0, T0, F0, F0, s0,
                                             smf)) {
        return -1;
      }
      return r;
    }
    // the axial strain at the end of the time step is computed by the
    // behaviour
    if (!internals::updateAxialDeformationGradient<Behaviour>(
            lgh1, d, d.s1.internal_state_variables)) {
      return -1;
    }
    const auto s1 = lgh1.convertToCauchyStress(T1);
    if (sm == StressMeasure::CAUCHY) {
      tfel::math::map<tfel::math::stensor<N, real>>(d.s1.thermodynamic_forces) =
          s1;
    } else if (sm == StressMeasure::PK1) {
      tfel::math::map<tfel::math::tensor<N, real>>(d.s1.thermodynamic_forces) =
          tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1, F1);
    } else {
      tfel::math::map<tfel::math::stensor<N, real>>(d.s1.thermodynamic_forces) =
          tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1, F1);
    }
    if (o.isTangentOperatorRequested()) {
      if (!internals::convertTangentOperator(d.K, d, K, lgh1, T1, F0, F1, s1,
                                             smf)) {
        return -1;
      }
    }
    return r;
  }  // end of integrate

  /*!
   * \brief integrate the behaviour using the Hencky strain measure.
   * \tparam Behaviour: class implementing the behaviour
   * \param[in,out] d: behaviour data. The computation options are given by
   * the three first values of the `K` array.
   * \param[in] p: out of bounds policy
   */
  template <typename Behaviour>
  int integrate(mfront_gb_BehaviourData& d,
                const tfel::material::OutOfBoundsPolicy p) {
    return integrate<Behaviour>(d, ComputationOptions{d.K[0], d.K[1], d.K[2]},
                                p);
  }  // end of integrate

  /*!
   * \brief execute the given post-processing.
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
   * post-processing (see `mfront::gb::executePostProcessing`).
   * \param[out] post_processing_variables: pointer to the values of the
   * post-processing variables.
   * \param[in,out] d: behaviour data.
   * \param[in] o: computation options (the stress measure and the tangent
   * operator, which defines the setting of the logarithmic strain handler,
   * are used).
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
      const ComputationOptions& o,
      const tfel::material::OutOfBoundsPolicy p,
      const ConstructorArguments&... args) {
    using namespace tfel::material;
    using Traits = GenericBehaviourTraits<Behaviour>;
    constexpr auto N =
        ModellingHypothesisToSpaceDimension<Traits::hypothesis>::value;
    using LogarithmicStrainHandler = internals::LogarithmicStrainHandler<N>;
    const auto sm = o.stress_measure;
    const auto smf = o.tangent_operator;
    if (smf == FiniteStrainTangentOperator::C_TRUESDELL) {
      reportError(d, "invalid choice for consistent tangent operator");
      return -1;
    }
    auto F0 = tfel::math::tensor<N, real>{};
    auto F1 = tfel::math::tensor<N, real>{};
    internals::getDeformationGradient<Behaviour>(F0, d.s0.gradients);
    internals::getDeformationGradient<Behaviour>(F1, d.s1.gradients);
    const auto setting = internals::getSetting<N>(smf);
    auto lgh0 = LogarithmicStrainHandler(setting, F0);
    auto lgh1 = LogarithmicStrainHandler(setting, F1);
    auto e0 = lgh0.getHenckyLogarithmicStrain();
    auto e1 = lgh1.getHenckyLogarithmicStrain();
    auto s0 = tfel::math::stensor<N, real>{};
    auto s1 = tfel::math::stensor<N, real>{};
    auto T0 = tfel::math::stensor<N, real>{};
    auto T1 = tfel::math::stensor<N, real>{};
    if ((!internals::updateAxialDeformationGradient<Behaviour>(
            lgh0, d, d.s0.internal_state_variables)) ||
        (!internals::updateAxialDeformationGradient<Behaviour>(
            lgh1, d, d.s1.internal_state_variables)) ||
        (!internals::convertToDualStress(T0, s0, d, d.s0.thermodynamic_forces,
                                         F0, lgh0, sm)) ||
        (!internals::convertToDualStress(T1, s1, d, d.s1.thermodynamic_forces,
                                         F1, lgh1, sm))) {
      return -1;
    }
    auto dc = d;
    dc.s0.gradients = e0.begin();
    dc.s1.gradients = e1.begin();
    dc.s0.thermodynamic_forces = T0.begin();
    dc.s1.thermodynamic_forces = T1.begin();
    const auto r =
        mfront::gb::executePostProcessing<Behaviour, use_initial_state>(
            f, post_processing_variables, dc, p, args...);
    return r;
  }  // end of executePostProcessing

  /*!
   * \brief execute the given post-processing
   * \tparam Behaviour: class describing the post-processing.
   * \tparam m: method implementing the post-processing.
   * \tparam use_initial_state: boolean stating if the postprocessing uses
   * the initial state of the material.
   * \param[out] post_processing_variables: pointer to the values of the
   * post-processing variables.
   * \param[in,out] d: behaviour data. The computation options are given by
   * the three first values of the `K` array.
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
        MethodCaller<m>{}, post_processing_variables, d,
        ComputationOptions{d.K[0], d.K[1], d.K[2]}, p);
  }  // end of executePostProcessing

}  // end of namespace mfront::gb::logarithmic_strain

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_LOGARITHMICSTRAININTEGRATE_HXX */
