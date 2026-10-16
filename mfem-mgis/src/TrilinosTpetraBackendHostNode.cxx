/*!
 * \file   src/TrilinosTpetraBackendHostNode.cxx
 * \brief  This file instantiates the solvers of Trilinos for the node of
 * Tpetra on the host, when the default node computes on a GPU (see
 * `TrilinosTpetraBackend.hxx`)
 * \date   16/10/2026
 */

#include "TrilinosTpetraBackend.hxx"

namespace mfem_mgis::trilinos {

#ifdef MFEM_MGIS_TRILINOS_HOST_NODE
  template struct TpetraBackend<HostNode>;
#endif /* MFEM_MGIS_TRILINOS_HOST_NODE */

}  // end of namespace mfem_mgis::trilinos
