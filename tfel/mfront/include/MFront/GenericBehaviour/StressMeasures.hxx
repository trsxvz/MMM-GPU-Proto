/*!
 * \file   mfront/include/MFront/GenericBehaviour/StressMeasures.hxx
 * \brief  This file declares the `convertFromCauchyStress` function, which
 * converts the Cauchy stress to the stress measures of the finite strain
 * behaviours of the generic interface. This function is shared by the
 * three integration headers (standard finite strain behaviours, Hencky and
 * Green-Lagrange strain frameworks).
 * \date   06/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_STRESSMEASURES_HXX
#define LIB_MFRONT_GENERICBEHAVIOUR_STRESSMEASURES_HXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include "MFront/GenericBehaviour/Integrate.hxx"

namespace mfront::gb::internals {

  /*!
   * \brief convert the Cauchy stress to the given stress measure
   * \return false if the stress measure is invalid
   * \param[out] v: values of the stress
   * \param[in,out] d: behaviour data
   * \param[in] s: Cauchy stress
   * \param[in] F: deformation gradient
   * \param[in] sm: stress measure
   */
  template <unsigned short N>
  TFEL_HOST_DEVICE bool convertFromCauchyStress(
      real* const v,
      mfront_gb_BehaviourData& d,
      const tfel::math::stensor<N, real>& s,
      const tfel::math::tensor<N, real>& F,
      const StressMeasure sm) {
    if (sm == StressMeasure::CAUCHY) {
      tfel::math::map<tfel::math::stensor<N, real>>(v) = s;
    } else if (sm == StressMeasure::PK1) {
      tfel::math::map<tfel::math::tensor<N, real>>(v) =
          tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s, F);
    } else if (sm == StressMeasure::PK2) {
      tfel::math::map<tfel::math::stensor<N, real>>(v) =
          tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s, F);
    } else {
      reportError(d, "invalid choice for the stress measure");
      return false;
    }
    return true;
  }  // end of convertFromCauchyStress

}  // end of namespace mfront::gb::internals

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_STRESSMEASURES_HXX */
