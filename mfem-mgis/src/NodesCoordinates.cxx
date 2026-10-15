/*!
 * \file   src/NodesCoordinates.cxx
 * \brief  This file implements the computation of the coordinates of the
 * nodes of the true degrees of freedom of a finite element space (see
 * `NodesCoordinates.hxx`)
 * \date   15/10/2026
 */

#include "mfem/config/config.hpp"

#ifdef MFEM_USE_MPI

#include "mfem/fem/pfespace.hpp"
#include "mfem/fem/pgridfunc.hpp"
#include "mfem/fem/coefficient.hpp"
#include "NodesCoordinates.hxx"

namespace mfem_mgis::internals {

  NodesCoordinates getNodesCoordinates(FiniteElementSpace<true>& fespace) {
    auto c = NodesCoordinates{};
    auto* const mesh = fespace.GetParMesh();
    const auto d = mesh->SpaceDimension();
    c.dimension = d;
    auto scalar = mfem::ParFiniteElementSpace(mesh, fespace.FEColl(), 1);
    for (int k = 0; k != d; ++k) {
      auto gf = mfem::ParGridFunction(&scalar);
      auto position = mfem::FunctionCoefficient(
          [k](const mfem::Vector& x) { return x(k); });
      gf.ProjectCoefficient(position);
      auto t = mfem::Vector{};
      gf.GetTrueDofs(t);
      const auto* const pt = t.HostRead();
      c.values.insert(c.values.end(), pt, pt + t.Size());
    }
    return c;
  }  // end of getNodesCoordinates

  bool hasRigidBodyModes(const FiniteElementSpace<true>& fespace) {
    const auto d = fespace.GetParMesh()->SpaceDimension();
    return (fespace.GetVDim() == d) && ((d == 2) || (d == 3));
  }  // end of hasRigidBodyModes

}  // end of namespace mfem_mgis::internals

#endif /* MFEM_USE_MPI */
