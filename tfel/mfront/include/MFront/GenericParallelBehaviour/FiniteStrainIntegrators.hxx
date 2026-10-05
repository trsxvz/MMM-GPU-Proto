/*!
 * \file
 * mfront/include/MFront/GenericParallelBehaviour/FiniteStrainIntegrators.hxx
 * \brief  This file declares the integration strategies of finite strain
 * behaviours for the `generic-parallel` interface.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_FINITESTRAININTEGRATORS_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_FINITESTRAININTEGRATORS_HXX

#include "MFront/GenericBehaviour/StandardFiniteStrainBehaviourIntegrate.hxx"
#include "MFront/GenericBehaviour/GreenLagrangeStrainIntegrate.hxx"
#include "MFront/GenericBehaviour/LogarithmicStrainIntegrate.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb {

  namespace internals {

    /*!
     * \brief base class of the integration strategies of finite strain
     * behaviours. The gradients are the deformation gradient. The
     * thermodynamic forces are the stress measure selected by the computation
     * options (Cauchy stress, first or second Piola-Kirchhoff stress) and the
     * tangent operator is the one selected by the computation options.
     */
    struct FiniteStrainIntegratorBase {
      /*!
       * \return the sizes of the data of a point. The generated maximal sizes
       * of the thermodynamic forces and of the tangent operator are the sizes
       * of a tensor and of a fourth order tensor.
       * \tparam PointTraits: description of the data of a point
       * \param[in] o: computation options
       */
      template <typename PointTraits>
      TFEL_HOST_DEVICE static constexpr PointDataSizes getPointDataSizes(
          const ComputationOptions& o) noexcept {
        using mfront::gb::FiniteStrainTangentOperator;
        using mfront::gb::StressMeasure;
        auto sizes = PointDataSizes{PointTraits::value};
        // size of a tensor
        const auto T = sizes.gradients;
        // size of a symmetric tensor (3, 5 and 9 for tensors in 1D, 2D and
        // 3D)
        const auto S = T == 9 ? size_type{6} : (T == 5 ? size_type{4} : T);
        sizes.thermodynamic_forces =
            o.stress_measure == StressMeasure::PK1 ? T : S;
        if (o.tangent_operator == FiniteStrainTangentOperator::DS_DEGL) {
          sizes.tangent_operator = S * S;
        } else if (o.tangent_operator == FiniteStrainTangentOperator::DPK1_DF) {
          sizes.tangent_operator = T * T;
        } else {
          sizes.tangent_operator = S * T;
        }
        return sizes;
      }  // end of getPointDataSizes
    };

  }  // end of namespace internals

  /*!
   * \brief integration strategy of standard finite strain behaviours, based
   * on `mfront::gb::finite_strain::integrate`.
   */
  struct FiniteStrainIntegrator : internals::FiniteStrainIntegratorBase {
    /*!
     * \brief integrate the behaviour at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour, typename ParametersType>
    TFEL_HOST_DEVICE static int integrate(
        mfront::gb::BehaviourData& d,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::finite_strain::integrate<Behaviour>(d, o, p,
                                                             parameters);
    }  // end of integrate
    /*!
     * \brief execute an initialize function at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam InitializeFunction: function object calling the method
     * implementing the initialize function
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] values: values of the inputs of the initialize function
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour,
              typename InitializeFunction,
              typename ParametersType>
    TFEL_HOST_DEVICE static int executeInitializeFunction(
        mfront::gb::BehaviourData& d,
        const real* const values,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::finite_strain::executeInitializeFunction<Behaviour>(
          InitializeFunction{}, d, values, o, p, parameters);
    }  // end of executeInitializeFunction
    /*!
     * \brief execute a post-processing at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam PostProcessing: function object calling the method
     * implementing the post-processing
     * \tparam use_initial_state: boolean stating if the post-processing uses
     * the initial state of the material
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[out] outputs: values of the outputs of the post-processing
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour,
              typename PostProcessing,
              bool use_initial_state,
              typename ParametersType>
    TFEL_HOST_DEVICE static int executePostProcessing(
        real* const outputs,
        mfront::gb::BehaviourData& d,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::finite_strain::executePostProcessing<
          Behaviour, use_initial_state>(PostProcessing{}, outputs, d, o, p,
                                        parameters);
    }  // end of executePostProcessing
  };

  /*!
   * \brief integration strategy of strain based behaviours written in the
   * Green-Lagrange strain, based on
   * `mfront::gb::green_lagrange_strain::integrate`.
   */
  struct GreenLagrangeStrainIntegrator : internals::FiniteStrainIntegratorBase {
    /*!
     * \brief integrate the behaviour at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour, typename ParametersType>
    TFEL_HOST_DEVICE static int integrate(
        mfront::gb::BehaviourData& d,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::green_lagrange_strain::integrate<Behaviour>(
          d, o, p, parameters);
    }  // end of integrate
    /*!
     * \brief execute an initialize function at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam InitializeFunction: function object calling the method
     * implementing the initialize function
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] values: values of the inputs of the initialize function
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour,
              typename InitializeFunction,
              typename ParametersType>
    TFEL_HOST_DEVICE static int executeInitializeFunction(
        mfront::gb::BehaviourData& d,
        const real* const values,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::green_lagrange_strain::executeInitializeFunction<
          Behaviour>(InitializeFunction{}, d, values, o, p, parameters);
    }  // end of executeInitializeFunction
    /*!
     * \brief execute a post-processing at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam PostProcessing: function object calling the method
     * implementing the post-processing
     * \tparam use_initial_state: boolean stating if the post-processing uses
     * the initial state of the material
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[out] outputs: values of the outputs of the post-processing
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour,
              typename PostProcessing,
              bool use_initial_state,
              typename ParametersType>
    TFEL_HOST_DEVICE static int executePostProcessing(
        real* const outputs,
        mfront::gb::BehaviourData& d,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::green_lagrange_strain::executePostProcessing<
          Behaviour, use_initial_state>(PostProcessing{}, outputs, d, o, p,
                                        parameters);
    }  // end of executePostProcessing
  };

  /*!
   * \brief integration strategy of strain based behaviours written in the
   * Hencky strain, based on `mfront::gb::logarithmic_strain::integrate`.
   */
  struct LogarithmicStrainIntegrator : internals::FiniteStrainIntegratorBase {
    /*!
     * \brief integrate the behaviour at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour, typename ParametersType>
    TFEL_HOST_DEVICE static int integrate(
        mfront::gb::BehaviourData& d,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::logarithmic_strain::integrate<Behaviour>(d, o, p,
                                                                  parameters);
    }  // end of integrate
    /*!
     * \brief execute an initialize function at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam InitializeFunction: function object calling the method
     * implementing the initialize function
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[in,out] d: behaviour data
     * \param[in] values: values of the inputs of the initialize function
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour,
              typename InitializeFunction,
              typename ParametersType>
    TFEL_HOST_DEVICE static int executeInitializeFunction(
        mfront::gb::BehaviourData& d,
        const real* const values,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::logarithmic_strain::executeInitializeFunction<
          Behaviour>(InitializeFunction{}, d, values, o, p, parameters);
    }  // end of executeInitializeFunction
    /*!
     * \brief execute a post-processing at one point
     * \tparam Behaviour: class implementing the behaviour
     * \tparam PostProcessing: function object calling the method
     * implementing the post-processing
     * \tparam use_initial_state: boolean stating if the post-processing uses
     * the initial state of the material
     * \tparam ParametersType: type of the structure holding the parameters
     * \param[out] outputs: values of the outputs of the post-processing
     * \param[in,out] d: behaviour data
     * \param[in] o: computation options
     * \param[in] p: out of bounds policy
     * \param[in] parameters: values of the parameters
     */
    template <typename Behaviour,
              typename PostProcessing,
              bool use_initial_state,
              typename ParametersType>
    TFEL_HOST_DEVICE static int executePostProcessing(
        real* const outputs,
        mfront::gb::BehaviourData& d,
        const ComputationOptions& o,
        const tfel::material::OutOfBoundsPolicy p,
        const ParametersType& parameters) {
      return mfront::gb::logarithmic_strain::executePostProcessing<
          Behaviour, use_initial_state>(PostProcessing{}, outputs, d, o, p,
                                        parameters);
    }  // end of executePostProcessing
  };

}  // end of namespace mfront::gpb

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_FINITESTRAININTEGRATORS_HXX */
