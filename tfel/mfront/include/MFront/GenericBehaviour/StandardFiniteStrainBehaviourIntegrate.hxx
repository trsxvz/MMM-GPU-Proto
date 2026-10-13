/*!
 * \file
 * mfront/include/MFront/GenericBehaviour/StandardFiniteStrainBehaviourIntegrate.hxx
 * \brief  This file declares the functions integrating the standard finite
 * strain behaviours, shared by the `generic` and `generic-parallel`
 * interfaces.
 * \author Thomas Helfer
 * \date   31/01/2022
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_STANDARDFINITESTRAINBEHAVIOURINTEGRATE_HXX
#define LIB_MFRONT_GENERICBEHAVIOUR_STANDARDFINITESTRAINBEHAVIOURINTEGRATE_HXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/FSAlgorithm/copy.hxx"
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include "MFront/GenericBehaviour/GenericBehaviourTraits.hxx"
#include "MFront/GenericBehaviour/Integrate.hxx"
#include "MFront/GenericBehaviour/StressMeasures.hxx"

namespace mfront::gb::finite_strain {

  namespace internals {

    /*!
     * \brief get the deformation gradient. In plane stress, the axial
     * component of the deformation gradient is computed by the behaviour and
     * stored in the internal state variables.
     * \return false if the axial component of the deformation gradient is
     * not computed by the behaviour (plane stress only)
     * \tparam Behaviour: class implementing the behaviour
     * \param[out] F: deformation gradient
     * \param[in,out] d: behaviour data
     * \param[in] gradients: gradients
     * \param[in] isvs: internal state variables
     */
    template <typename Behaviour, unsigned short N>
    TFEL_HOST_DEVICE bool getDeformationGradient(
        tfel::math::tensor<N, real>& F,
        mfront_gb_BehaviourData& d,
        const real* const gradients,
        [[maybe_unused]] const real* const isvs) {
      using namespace tfel::material;
      using Traits = GenericBehaviourTraits<Behaviour>;
      constexpr auto hypothesis = Traits::hypothesis;
      constexpr auto TensorSize =
          ModellingHypothesisToTensorSize<hypothesis>::value;
      tfel::fsalgo::copy<TensorSize>::exe(gradients, F.begin());
      if constexpr ((hypothesis == ModellingHypothesis::PLANESTRESS) ||
                    (hypothesis == ModellingHypothesis::
                                       AXISYMMETRICALGENERALISEDPLANESTRESS)) {
        if constexpr (Traits::has_axial_deformation_gradient_offset) {
          const auto Fzz = isvs[Traits::axial_deformation_gradient_offset];
          if constexpr (hypothesis == ModellingHypothesis::PLANESTRESS) {
            F[2] += Fzz;
          } else {
            F[1] += Fzz;
          }
        } else {
          reportError(d,
                      "the axial deformation gradient is not computed by the "
                      "behaviour");
          return false;
        }
      }
      return true;
    }  // end of getDeformationGradient

    /*!
     * \brief convert the given stress to the Cauchy stress
     * \return false if the stress measure is invalid
     * \param[out] s: Cauchy stress
     * \param[in,out] d: behaviour data
     * \param[in] v: values of the stress
     * \param[in] F: deformation gradient
     * \param[in] sm: stress measure
     */
    template <unsigned short N>
    TFEL_HOST_DEVICE bool convertToCauchyStress(
        tfel::math::stensor<N, real>& s,
        mfront_gb_BehaviourData& d,
        const real* const v,
        const tfel::math::tensor<N, real>& F,
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
      return true;
    }  // end of convertToCauchyStress
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
    const auto sm = o.stress_measure;
    if (sm == StressMeasure::INVALID_STRESS_MEASURE) {
      reportError(d, "invalid choice for the stress measure");
      return -1;
    }
    if (sm == StressMeasure::CAUCHY) {
      return mfront::gb::executeInitializeFunction<Behaviour>(
          f, d, initialize_variables, p, args...);
    }
    // the behaviour works with the Cauchy stress: the initialize function is
    // executed with a copy of the behaviour data referring to local Cauchy
    // stresses, and the stress it computes is converted. The initialize
    // function uses the state at the beginning of the time step.
    auto F0 = tfel::math::tensor<N, real>{};
    auto s0 = tfel::math::stensor<N, real>{};
    if ((!internals::getDeformationGradient<Behaviour>(
            F0, d, d.s0.gradients, d.s0.internal_state_variables)) ||
        (!internals::convertToCauchyStress(s0, d, d.s0.thermodynamic_forces, F0,
                                           sm))) {
      return -1;
    }
    auto s1 = s0;
    auto dc = d;
    dc.s0.thermodynamic_forces = s0.begin();
    dc.s1.thermodynamic_forces = s1.begin();
    const auto r = mfront::gb::executeInitializeFunction<Behaviour>(
        f, dc, initialize_variables, p, args...);
    if (r == -1) {
      return r;
    }
    if (!mfront::gb::internals::convertFromCauchyStress(
            d.s1.thermodynamic_forces, d, s1, F0, sm)) {
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
   * \brief integrate a finite strain behaviour.
   *
   * The behaviour computes the Cauchy stress. The stress given by the
   * caller, if it is not the Cauchy stress, is converted before and after
   * the integration. The tangent operator requested by the caller is
   * computed by the behaviour.
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
    if (sm == StressMeasure::CAUCHY) {
      return mfront::gb::integrate<Behaviour>(d, o, smf, p, args...);
    }
    // the behaviour works with the Cauchy stress: the behaviour is
    // integrated with a copy of the behaviour data referring to local Cauchy
    // stresses, and the stress at the end of the time step is converted
    auto F0 = tfel::math::tensor<N, real>{};
    auto s0 = tfel::math::stensor<N, real>{};
    if ((!internals::getDeformationGradient<Behaviour>(
            F0, d, d.s0.gradients, d.s0.internal_state_variables)) ||
        (!internals::convertToCauchyStress(s0, d, d.s0.thermodynamic_forces, F0,
                                           sm))) {
      return -1;
    }
    auto s1 = s0;
    auto dc = d;
    dc.s0.thermodynamic_forces = s0.begin();
    dc.s1.thermodynamic_forces = s1.begin();
    const auto r = mfront::gb::integrate<Behaviour>(dc, o, smf, p, args...);
    if ((r == -1) || (o.isPredictionRequested())) {
      return r;
    }
    auto F1 = tfel::math::tensor<N, real>{};
    if ((!internals::getDeformationGradient<Behaviour>(
            F1, d, d.s1.gradients, d.s1.internal_state_variables)) ||
        (!mfront::gb::internals::convertFromCauchyStress(
            d.s1.thermodynamic_forces, d, s1, F1, sm))) {
      return -1;
    }
    return r;
  }  // end of integrate

  /*!
   * \brief integrate a finite strain behaviour.
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
   * \param[in] o: computation options (only the stress measure is used).
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
    const auto sm = o.stress_measure;
    if (sm == StressMeasure::INVALID_STRESS_MEASURE) {
      reportError(d, "invalid choice for the stress measure");
      return -1;
    }
    if (sm == StressMeasure::CAUCHY) {
      return mfront::gb::executePostProcessing<Behaviour, use_initial_state>(
          f, post_processing_variables, d, p, args...);
    }
    // the behaviour works with the Cauchy stress: the post-processing is
    // executed with a copy of the behaviour data referring to local Cauchy
    // stresses
    auto F0 = tfel::math::tensor<N, real>{};
    auto F1 = tfel::math::tensor<N, real>{};
    auto s0 = tfel::math::stensor<N, real>{};
    auto s1 = tfel::math::stensor<N, real>{};
    if ((!internals::getDeformationGradient<Behaviour>(
            F0, d, d.s0.gradients, d.s0.internal_state_variables)) ||
        (!internals::getDeformationGradient<Behaviour>(
            F1, d, d.s1.gradients, d.s1.internal_state_variables)) ||
        (!internals::convertToCauchyStress(s0, d, d.s0.thermodynamic_forces, F0,
                                           sm)) ||
        (!internals::convertToCauchyStress(s1, d, d.s1.thermodynamic_forces, F1,
                                           sm))) {
      return -1;
    }
    auto dc = d;
    dc.s0.thermodynamic_forces = s0.begin();
    dc.s1.thermodynamic_forces = s1.begin();
    return mfront::gb::executePostProcessing<Behaviour, use_initial_state>(
        f, post_processing_variables, dc, p, args...);
  }  // end of executePostProcessing

  /*!
   * \brief execute the given post-processing
   * \tparam Behaviour: class describing the post-processing.
   * \tparam m: method implementing the post-processing.
   * \tparam use_initial_state: boolean stating if the postprocessing uses
   * the initial state of the material.
   * \param[out] post_processing_variables: pointer to the values of the
   * post-processing variables.
   * \param[in,out] d: behaviour data. The stress measure is given by the
   * second value of the `K` array.
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

}  // end of namespace mfront::gb::finite_strain

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_STANDARDFINITESTRAINBEHAVIOURINTEGRATE_HXX \
        */
