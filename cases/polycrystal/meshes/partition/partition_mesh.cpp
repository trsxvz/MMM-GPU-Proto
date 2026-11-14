// Partitioned mesh of the polycrystal, read by mfem-mgis with the mesh read
// mode "Restart" (one file per process, mfem::MakeParFilename):
// usage: mpirun -n <processes> partition-mesh <mesh> <refinements> <prefix>
//
// The mesh is built as by mfem-mgis in the mode "FromScratch"
// (MeshDiscretization.cxx): serial mesh read with mfem::Mesh(name, 0, 1),
// parallel mesh with the default partitioning of MFEM (METIS), then the
// uniform refinements of the parallel mesh. Each process writes its part
// after the refinements (mfem::ParMesh::ParPrint, 17 significant digits):
// the computation reads it with "Restart" and no refinement, and gets the
// same mesh, numbering and partition as with "FromScratch", without each
// process reading and partitioning the whole serial mesh.
#include <fstream>
#include <iostream>
#include <string>
#include "mfem.hpp"

int main(int argc, char *argv[]) {
  mfem::Mpi::Init(argc, argv);
  if (argc != 4) {
    if (mfem::Mpi::Root()) {
      std::cerr << "usage: partition-mesh <mesh> <refinements> <prefix>\n";
    }
    return EXIT_FAILURE;
  }
  const std::string name = argv[1];
  const int refinements = std::stoi(argv[2]);
  const std::string prefix = argv[3];
  auto smesh = mfem::Mesh(name.c_str(), 0, 1);
  // as in mfem-mgis: the serial mesh is refined only if a process would
  // have fewer than 8 elements
  int done = 0;
  while ((double(smesh.GetNE()) / mfem::Mpi::WorldSize() < 8) &&
         (done < refinements)) {
    smesh.UniformRefinement();
    ++done;
  }
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  smesh.Clear();
  for (int i = done; i < refinements; ++i) {
    pmesh.UniformRefinement();
  }
  const auto file = mfem::MakeParFilename(prefix, mfem::Mpi::WorldRank());
  std::ofstream out(file);
  out.precision(17);
  pmesh.ParPrint(out);
  if (!out) {
    std::cerr << "can't write " << file << '\n';
    return EXIT_FAILURE;
  }
  out.close();
  long long ne = pmesh.GetNE(), total = 0;
  MPI_Reduce(&ne, &total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mfem::Mpi::Root()) {
    std::cout << prefix << ": " << mfem::Mpi::WorldSize() << " parts, "
              << total << " elements\n";
  }
  return EXIT_SUCCESS;
}
