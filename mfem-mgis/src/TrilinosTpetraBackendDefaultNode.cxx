/*!
 * \file   src/TrilinosTpetraBackendDefaultNode.cxx
 * \brief  This file instantiates the solvers of Trilinos for the default
 * node of Tpetra (see `TrilinosTpetraBackend.hxx`)
 * \date   16/10/2026
 */

#include "TrilinosTpetraBackend.hxx"

namespace mfem_mgis::trilinos {

  template struct TpetraBackend<DefaultNode>;

}  // end of namespace mfem_mgis::trilinos
