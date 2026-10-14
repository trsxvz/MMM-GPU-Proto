/*!
 * \file   mfront/include/MFront/GenericBehaviour/DisplacementGradient.hxx
 * \brief  This file declares the functions computing the deformation
 * gradient and the Green-Lagrange strain from the displacement gradient
 * \(H = F - I\) given by the caller of a finite strain behaviour (see
 * `mfront::gb::displacement_gradient_stress_measure_offset`). These
 * functions are shared by the three integration headers (standard finite
 * strain behaviours, Hencky and Green-Lagrange strain frameworks).
 * \date   14/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_DISPLACEMENTGRADIENT_HXX
#define LIB_MFRONT_GENERICBEHAVIOUR_DISPLACEMENTGRADIENT_HXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/FSAlgorithm/copy.hxx"
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "MFront/GenericBehaviour/Types.hxx"

namespace mfront::gb {

  /*!
   * \brief axial component of the gradients of a finite strain behaviour in
   * plane stress (`PLANESTRESS` and `AXISYMMETRICALGENERALISEDPLANESTRESS`),
   * computed by the behaviour: the one given by the caller is not used by
   * the interface, neither in the strain given to the behaviour nor in the
   * conversions of the stresses and of the tangent operators.
   * \tparam hypothesis: modelling hypothesis
   */
  template <tfel::material::ModellingHypothesis::Hypothesis hypothesis>
  struct PlaneStressAxialComponent {
    //! \brief true in plane stress
    static constexpr bool exists =
        (hypothesis == tfel::material::ModellingHypothesis::PLANESTRESS) ||
        (hypothesis == tfel::material::ModellingHypothesis::
                           AXISYMMETRICALGENERALISEDPLANESTRESS);
    //! \brief index of the axial component, if any
    static constexpr unsigned short index =
        (hypothesis == tfel::material::ModellingHypothesis::PLANESTRESS) ? 2
                                                                         : 1;
  };

  /*!
   * \brief compute the deformation gradient \(F = I + H\) from the
   * displacement gradient \(H\) given to a finite strain behaviour (the
   * three first components of the tensors are the diagonal ones).
   *
   * \tparam hypothesis: modelling hypothesis
   * \param[out] F: deformation gradient
   * \param[in] H: displacement gradient
   */
  template <tfel::material::ModellingHypothesis::Hypothesis hypothesis>
  TFEL_HOST_DEVICE void computeDeformationGradient(
      real* const F, const real* const H) noexcept {
    using namespace tfel::material;
    constexpr auto TensorSize =
        ModellingHypothesisToTensorSize<hypothesis>::value;
    tfel::fsalgo::copy<TensorSize>::exe(H, F);
    F[0] += real{1};
    F[1] += real{1};
    F[2] += real{1};
  }  // end of computeDeformationGradient

  /*!
   * \brief compute the deformation gradient \(F = I + H\) and the
   * Green-Lagrange strain from the displacement gradient \(H\) given to a
   * finite strain behaviour. The Green-Lagrange strain is computed from
   * \(H\), without the cancellation of \(\left(F^{T}\,F - I\right)/2\) (see
   * `tfel::math::computeGreenLagrangeTensorFromDisplacementGradient`).
   *
   * In plane stress, the axial component of the Green-Lagrange strain is
   * computed by the behaviour: the one given here is zero (see
   * `PlaneStressAxialComponent`).
   *
   * \tparam hypothesis: modelling hypothesis
   * \param[out] e: Green-Lagrange strain
   * \param[out] F: deformation gradient
   * \param[in] H: displacement gradient
   */
  template <tfel::material::ModellingHypothesis::Hypothesis hypothesis,
            unsigned short N>
  TFEL_HOST_DEVICE void computeGreenLagrangeStrainFromDisplacementGradient(
      tfel::math::stensor<N, real>& e,
      tfel::math::tensor<N, real>& F,
      const real* const H) noexcept {
    using namespace tfel::material;
    using AxialComponent = PlaneStressAxialComponent<hypothesis>;
    constexpr auto TensorSize =
        ModellingHypothesisToTensorSize<hypothesis>::value;
    auto Hv = tfel::math::tensor<N, real>{};
    tfel::fsalgo::copy<TensorSize>::exe(H, Hv.begin());
    computeDeformationGradient<hypothesis>(F.begin(), H);
    e = tfel::math::computeGreenLagrangeTensorFromDisplacementGradient(Hv);
    if constexpr (AxialComponent::exists) {
      e[AxialComponent::index] = real{0};
    }
  }  // end of computeGreenLagrangeStrainFromDisplacementGradient

}  // end of namespace mfront::gb

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_DISPLACEMENTGRADIENT_HXX */
