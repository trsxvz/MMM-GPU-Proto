/*!
 * \file   include/TFEL/Material/PolyCrystalsSlidingSystems.hxx
 * \brief
 * \author Thomas Helfer
 * \brief  16 mai 2013
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef TFEL_MATERIAL_POLYCRISTALSSLIDINGSYSTEMS_HXX
#define TFEL_MATERIAL_POLYCRISTALSSLIDINGSYSTEMS_HXX 1

#include <string>

#include "TFEL/Math/tvector.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/General/RandomRotationMatrix.hxx"

namespace tfel::material {

  /*!
   * \param[in] Np      : number of phases
   * \param[in] GS      : class describing the sliding system of one phase
   * \param[in] NumType : numeric type used
   */
  template <unsigned short Ng, typename GS, typename NumType>
  struct PolyCrystalsSlidingSystems {
    //! number of sliding systems of the mono crystal
    static constexpr unsigned short Nss = GS::Nss;
    //! a simple alias
    typedef NumType real;
    //! a simple alias
    typedef tfel::math::stensor<3u> StrainStensor;
    //! return the unique instance of the class
    static const PolyCrystalsSlidingSystems& getPolyCrystalsSlidingSystems(
        const std::string& = "");
    /*!
     * \brief copy constructor. The data are stored in fixed size arrays, so
     * that a copy can be used in device code (see the `@SharedData` keyword
     * of `MFront`).
     */
    PolyCrystalsSlidingSystems(const PolyCrystalsSlidingSystems&) = default;
    //! \brief assignement operator (deleted)
    PolyCrystalsSlidingSystems& operator=(const PolyCrystalsSlidingSystems&) =
        delete;
    //! tensor of directional senses, sorted by phases
    tfel::math::tvector<Ng, tfel::math::tvector<Nss, StrainStensor>> mus;

   private:
    /*!
     * default constructor (declared private~: it can only be used
     * by the getPolyCrystalsSlidingSystems method)
     * \param[in] f : file where the angles are declared
     */
    PolyCrystalsSlidingSystems(const std::string&);
  };

  /*!
   * \brief the same class as PolyCrystalsSlidingSystems but with
   * volume fraction associated to each grain (fourth item of each line of the
   * file given as an argument to the method getPolyCrystalsSlidingSystems)
   * \param[in] Np      : number of phases
   * \param[in] GS      : class describing the sliding system of one phase
   * \param[in] NumType : numeric type used
   */
  template <unsigned short Ng, typename GS, typename NumType>
  struct ExtendedPolyCrystalsSlidingSystems {
    static constexpr unsigned short Nss = GS::Nss;
    typedef NumType real;
    typedef tfel::math::stensor<3u, real> StrainStensor;

    static const ExtendedPolyCrystalsSlidingSystems&
    getPolyCrystalsSlidingSystems(const std::string&);
    //! \brief copy constructor (see `PolyCrystalsSlidingSystems`)
    ExtendedPolyCrystalsSlidingSystems(
        const ExtendedPolyCrystalsSlidingSystems&) = default;
    //! \brief assignement operator (deleted)
    ExtendedPolyCrystalsSlidingSystems& operator=(
        const ExtendedPolyCrystalsSlidingSystems&) = delete;

    tfel::math::tvector<Ng, tfel::math::tvector<Nss, StrainStensor>> mus;

    //! \brief volume fractions per phases
    tfel::math::tvector<Ng, real> volume_fractions;

   private:
    ExtendedPolyCrystalsSlidingSystems(const std::string&);
  };

}  // end of namespace tfel::material

#include "TFEL/Material/PolyCrystalsSlidingSystems.ixx"

#endif /* TFEL_MATERIAL_POLYCRISTALSSLIDINGSYSTEMS_HXX */
