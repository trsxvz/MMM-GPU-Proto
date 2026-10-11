/*!
 * \file   src/AbstractBoundaryCondition.cxx
 * \brief  This file implements the `AbstractBoundaryCondition` class
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#include "MFEMMGIS/AbstractBoundaryCondition.hxx"

namespace mfem_mgis {

  bool AbstractBoundaryCondition::isIndependentOfTheUnknowns() const noexcept {
    return false;
  }  // end of isIndependentOfTheUnknowns

  bool AbstractBoundaryCondition::addContributionToTheResidual(
      Context& ctx, mfem::Vector&) noexcept {
    return ctx.registerErrorMessage(
        "AbstractBoundaryCondition::addContributionToTheResidual: the "
        "boundary condition depends on the unknowns");
  }  // end of addContributionToTheResidual

  AbstractBoundaryCondition::~AbstractBoundaryCondition() = default;

}  // end of namespace mfem_mgis
