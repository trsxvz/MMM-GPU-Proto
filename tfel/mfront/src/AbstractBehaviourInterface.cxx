/*!
 * \file   mfront/src/AbstractBehaviourInterface.cxx
 *
 * \brief
 * \author Thomas Helfer
 * \date 17/01/2007
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include "MFront/AbstractBehaviourInterface.hxx"

namespace mfront {

  std::string
  AbstractBehaviourInterface::getParametersSourceInBehaviourConstructor(
      const BehaviourDescription&, const Hypothesis) const {
    return "";
  }  // end of getParametersSourceInBehaviourConstructor

  bool AbstractBehaviourInterface::requiresDeviceCompatibleBehaviourClass()
      const {
    return false;
  }  // end of requiresDeviceCompatibleBehaviourClass

  std::map<std::string, std::shared_ptr<AbstractBehaviourInterface>>
  AbstractBehaviourInterface::getInterfacesOfAdditionalBehaviourConstructors(
      const BehaviourDescription&) const {
    return {};
  }  // end of getInterfacesOfAdditionalBehaviourConstructors

  AbstractBehaviourInterface::~AbstractBehaviourInterface() = default;

}  // end of namespace mfront
