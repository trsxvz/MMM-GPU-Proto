/*!
 * \file   cases/tools/count_unknowns.cxx
 * \brief  Number of unknowns of a vector H1 finite element space on a mesh
 * refined uniformly, for the context of the benchmark of the examples
 * (cases/benchmark_slides.py), whose programs do not all print it. The
 * mesh is not made periodic: for the periodic meshes, the number given
 * slightly exceeds the one of the periodic space.
 *
 * usage: count_unknowns mesh refinements order components
 */

#include <cstdlib>
#include <iostream>
#include "mfem.hpp"

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: " << argv[0]
              << " mesh refinements order components\n";
    return EXIT_FAILURE;
  }
  auto mesh = mfem::Mesh::LoadFromFile(argv[1], 1, 1);
  const auto refinements = std::atoi(argv[2]);
  const auto order = std::atoi(argv[3]);
  const auto components = std::atoi(argv[4]);
  for (int i = 0; i != refinements; ++i) {
    mesh.UniformRefinement();
  }
  auto collection = mfem::H1_FECollection(order, mesh.Dimension());
  auto space = mfem::FiniteElementSpace(&mesh, &collection, components);
  std::cout << space.GetTrueVSize() << '\n';
  return EXIT_SUCCESS;
}
