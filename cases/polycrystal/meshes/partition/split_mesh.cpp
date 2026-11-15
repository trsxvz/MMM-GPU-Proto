// Partitioned mesh of the polycrystal in any number of parts, written by a
// few processes: the parts of partition-mesh (partition_mesh.cpp) when the
// processes of the computation are too many to read the whole serial mesh
// each (6144 processes, 39.6 M tetrahedra: 10 to 15 GB per process).
//
// usage: mpirun -n <workers> split-mesh <mesh> <refinements> <parts> <prefix>
//        mpirun -n <parts> split-mesh --mfem <mesh> <refinements> <prefix>
//        mpirun -n <workers> split-mesh --check <parts> <prefix>
//
// First form. Each worker reads the serial mesh as mfem-mgis does
// (mfem::Mesh(name, 0, 1)); the first one computes the partition that
// mfem::ParMesh would compute for <parts> processes
// (Mesh::GeneratePartitioning(parts, 1), METIS) and sends it to the others.
// The worker w then builds the parts p = w modulo <workers>: the part of the
// serial mesh and its shared entities (mfem::MeshPartitioner), refined
// <refinements> times, written in the format of mfem::ParMesh::ParPrint (17
// significant digits) in <prefix><p on 6 digits>. No process of rank p of
// <parts> exists: the refinement of the shared entities, which MFEM does in
// ParMesh::UniformRefineGroups3D, is done here on the lists of the part
// (tetrahedra only). The local mesh is refined by MFEM
// (Mesh::UniformRefinement), with the same choice of the diagonals as in the
// parallel mesh (it only depends on the vertices of each tetrahedron in the
// serial mesh): the global refined mesh is the one of partition-mesh.
//
// Second form (--mfem): the same parts made by MFEM alone, one process per
// part, each reading the whole mesh (part of mfem::MeshPartitioner read by
// mfem::ParMesh, refined and written by ParPrint). It checks the first form:
// the files are identical.
//
// Third form (--check): reads the parts (tetrahedra) and checks that they
// describe one conforming mesh: each group of parts is described by all its
// parts with the same numbers of shared vertices, edges and triangles, at
// the same positions in the same order (modulo the box of the mesh when it
// is periodic), and the shared edges and triangles are edges and boundary
// faces of the local mesh. Prints the global numbers of elements, boundary
// elements, vertices and edges (the P2 nodes are the vertices and the
// edges), and a fingerprint of the elements (sum of the hashes of their
// attributes and of the positions of their vertices), the same for all the
// partitions of the same mesh.
#include <algorithm>
#include <cmath>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "mfem.hpp"

namespace {

// shared entities of a group of parts, as local vertices
struct Shared {
  std::vector<int> vertices, edges, triangles;
};

// access to the writer of the serial part of the parallel format
struct PartMesh : mfem::Mesh {
  explicit PartMesh(mfem::Mesh &&m) : mfem::Mesh(std::move(m)) {}
  using mfem::Mesh::Printer;
};

// largest resident memory of the process, in GB
double peakMemory() {
  std::ifstream status("/proc/self/status");
  std::string word;
  while (status >> word) {
    if (word == "VmHWM:") {
      double kb = 0;
      status >> kb;
      return kb / 1048576;
    }
  }
  return 0;
}

// One uniform refinement of the local mesh and of its shared entities, as
// ParMesh::UniformRefinement3D: the vertex in the middle of an edge is
// numbered after the vertices of the coarse mesh, by the rank of the edge
// in the lexicographic order of its two vertices
// (Mesh::UniformRefinement3D_base, e2v).
void refine(mfem::Mesh &mesh, std::vector<Shared> &groups) {
  const int nv = mesh.GetNV();
  std::vector<std::pair<int, int>> edges;
  edges.reserve(6 * static_cast<std::size_t>(mesh.GetNE()));
  for (int i = 0; i < mesh.GetNE(); ++i) {
    MFEM_VERIFY(mesh.GetElementType(i) == mfem::Element::TETRAHEDRON,
                "refinement of the parts: tetrahedra only");
    const int *v = mesh.GetElement(i)->GetVertices();
    for (int a = 0; a < 4; ++a) {
      for (int b = a + 1; b < 4; ++b) {
        edges.push_back(std::minmax(v[a], v[b]));
      }
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  auto middle = [&](const int a, const int b) {
    const std::pair<int, int> e = std::minmax(a, b);
    const auto p = std::lower_bound(edges.begin(), edges.end(), e);
    MFEM_VERIFY(p != edges.end() && *p == e,
                "shared edge absent from the local mesh");
    return nv + static_cast<int>(p - edges.begin());
  };
  for (auto &g : groups) {
    // each shared edge (a, b) of middle m: vertex m, edges (a, m) and (b, m)
    const std::size_t ne = g.edges.size() / 2;
    for (std::size_t i = 0; i < ne; ++i) {
      const int a = g.edges[2 * i], b = g.edges[2 * i + 1];
      const int m = middle(a, b);
      g.vertices.push_back(m);
      g.edges.push_back(b);
      g.edges.push_back(m);
      g.edges[2 * i + 1] = m;
    }
    // each shared triangle: three inner edges, four triangles
    const std::size_t nt = g.triangles.size() / 3;
    for (std::size_t i = 0; i < nt; ++i) {
      const int v0 = g.triangles[3 * i], v1 = g.triangles[3 * i + 1],
                v2 = g.triangles[3 * i + 2];
      const int m0 = middle(v0, v1), m1 = middle(v1, v2), m2 = middle(v2, v0);
      for (const int v : {m0, m1, m1, m2, m0, m2}) {
        g.edges.push_back(v);
      }
      for (const int v : {m1, m2, m0, m0, v1, m1, m2, m1, v2}) {
        g.triangles.push_back(v);
      }
      g.triangles[3 * i + 1] = m0;
      g.triangles[3 * i + 2] = m2;
    }
  }
  mesh.UniformRefinement();
  MFEM_VERIFY(static_cast<std::size_t>(mesh.GetNV()) == nv + edges.size(),
              "unexpected numbering of the refined mesh");
}

// part refined and written in the format of ParMesh::ParPrint; returns the
// number of elements
long long writePart(mfem::MeshPart &part, const int refinements,
                    const std::string &file) {
  std::vector<Shared> groups(part.my_groups.Size());
  const mfem::Table &gv =
      part.group_shared_entity_to_vertex[mfem::Geometry::POINT];
  const mfem::Table &ge =
      part.group_shared_entity_to_vertex[mfem::Geometry::SEGMENT];
  const mfem::Table &gt =
      part.group_shared_entity_to_vertex[mfem::Geometry::TRIANGLE];
  MFEM_VERIFY(part.group_shared_entity_to_vertex[mfem::Geometry::SQUARE]
                      .Size_of_connections() == 0,
              "shared quadrilaterals not handled");
  for (int g = 0; g < part.my_groups.Size(); ++g) {
    groups[g].vertices.assign(gv.GetRow(g), gv.GetRow(g) + gv.RowSize(g));
    groups[g].edges.assign(ge.GetRow(g), ge.GetRow(g) + ge.RowSize(g));
    groups[g].triangles.assign(gt.GetRow(g), gt.GetRow(g) + gt.RowSize(g));
  }
  PartMesh mesh(std::move(part.GetMesh()));
  for (int i = 0; i < refinements; ++i) {
    refine(mesh, groups);
  }
  std::ofstream out(file);
  out.precision(17);
  mesh.Printer(out, "mfem_serial_mesh_end");
  // as mfem::GroupTopology::Save
  out << "\ncommunication_groups\n";
  out << "number_of_groups " << groups.size() << "\n\n";
  out << "# number of entities in each group, followed by ranks in group\n";
  for (int g = 0; g < part.my_groups.Size(); ++g) {
    out << part.my_groups.RowSize(g);
    for (int i = 0; i < part.my_groups.RowSize(g); ++i) {
      out << ' ' << part.my_groups.GetRow(g)[i];
    }
    out << '\n';
  }
  // as mfem::ParMesh::ParPrint
  std::size_t nv = 0, ne = 0, nt = 0;
  for (const auto &g : groups) {
    nv += g.vertices.size();
    ne += g.edges.size() / 2;
    nt += g.triangles.size() / 3;
  }
  out << "\ntotal_shared_vertices " << nv << '\n';
  out << "total_shared_edges " << ne << '\n';
  out << "total_shared_faces " << nt << '\n';
  out << "\n# group 0 has no shared entities\n";
  for (std::size_t g = 1; g < groups.size(); ++g) {
    const Shared &s = groups[g];
    out << "\n# group " << g << "\nshared_vertices " << s.vertices.size()
        << '\n';
    for (const int v : s.vertices) {
      out << v << '\n';
    }
    out << "\nshared_edges " << s.edges.size() / 2 << '\n';
    for (std::size_t i = 0; i < s.edges.size(); i += 2) {
      out << s.edges[i] << ' ' << s.edges[i + 1] << '\n';
    }
    out << "\nshared_faces " << s.triangles.size() / 3 << '\n';
    for (std::size_t i = 0; i < s.triangles.size(); i += 3) {
      out << mfem::Geometry::TRIANGLE << ' ' << s.triangles[i] << ' '
          << s.triangles[i + 1] << ' ' << s.triangles[i + 2] << '\n';
    }
  }
  out << "\nmfem_mesh_end" << std::endl;
  out.close();
  MFEM_VERIFY(out, "can't write " << file);
  return mesh.GetNE();
}

// serial mesh read and partitioned as by mfem-mgis in the mode "FromScratch"
// for <parts> processes; the partition is computed by the first process
struct SerialMesh {
  mfem::Mesh mesh;
  mfem::Array<int> partitioning;
  int serial_refinements = 0;
  SerialMesh(const std::string &name, const int refinements, const int parts)
      : mesh(name.c_str(), 0, 1) {
    // as in mfem-mgis: the serial mesh is refined only if a process would
    // have fewer than 8 elements
    while ((double(mesh.GetNE()) / parts < 8) &&
           (serial_refinements < refinements)) {
      mesh.UniformRefinement();
      ++serial_refinements;
    }
    partitioning.SetSize(mesh.GetNE());
    if (mfem::Mpi::Root()) {
      int *p = mesh.GeneratePartitioning(parts, 1);
      std::copy(p, p + mesh.GetNE(), partitioning.GetData());
      delete[] p;
    }
    MPI_Bcast(partitioning.GetData(), mesh.GetNE(), MPI_INT, 0,
              MPI_COMM_WORLD);
  }
};

int split(const std::string &name, const int refinements, const int parts,
          const std::string &prefix, const bool by_mfem) {
  const int rank = mfem::Mpi::WorldRank(), workers = mfem::Mpi::WorldSize();
  double t0 = MPI_Wtime();
  SerialMesh serial(name, refinements, parts);
  const int left = refinements - serial.serial_refinements;
  mfem::MeshPartitioner partitioner(serial.mesh, parts,
                                    serial.partitioning.GetData());
  const double t1 = MPI_Wtime();
  long long ne = 0, total = 0;
  mfem::MeshPart part;
  if (by_mfem) {
    std::stringstream coarse;
    coarse.precision(17);
    partitioner.ExtractPart(rank, part);
    part.Print(coarse);
    mfem::ParMesh pmesh(MPI_COMM_WORLD, coarse, false, 0, false);
    for (int i = 0; i < left; ++i) {
      pmesh.UniformRefinement();
    }
    std::ofstream out(mfem::MakeParFilename(prefix, rank));
    out.precision(17);
    pmesh.ParPrint(out);
    out.close();
    MFEM_VERIFY(out, "can't write the part " << rank);
    ne = pmesh.GetNE();
  } else {
    for (int p = rank; p < parts; p += workers) {
      partitioner.ExtractPart(p, part);
      ne += writePart(part, left, mfem::MakeParFilename(prefix, p));
    }
  }
  const double t2 = MPI_Wtime();
  double memory = peakMemory(), max_memory = 0;
  MPI_Reduce(&ne, &total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&memory, &max_memory, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (mfem::Mpi::Root()) {
    std::cout << "reading and partition: " << t1 - t0 << " s, parts: "
              << t2 - t1 << " s, " << workers
              << " processes, largest resident memory " << max_memory
              << " GB\n"
              << prefix << ": " << parts << " parts, " << total
              << " elements\n";
  }
  return EXIT_SUCCESS;
}

// tokens of a part file, without the comment lines
struct Reader {
  std::string text;
  const char *p = nullptr;
  explicit Reader(const std::string &file) {
    std::ifstream in(file, std::ios::binary);
    std::stringstream s;
    s << in.rdbuf();
    text = s.str();
    p = text.c_str();
  }
  void skip() {
    while (true) {
      while (*p == ' ' || *p == '\n' || *p == '\t' || *p == '\r') {
        ++p;
      }
      if (*p != '#') {
        return;
      }
      while (*p != '\n' && *p != '\0') {
        ++p;
      }
    }
  }
  long integer() {
    skip();
    char *end = nullptr;
    const long v = std::strtol(p, &end, 10);
    if (end == p) {
      throw std::runtime_error("integer expected");
    }
    p = end;
    return v;
  }
  double real() {
    skip();
    char *end = nullptr;
    const double v = std::strtod(p, &end);
    if (end == p) {
      throw std::runtime_error("number expected");
    }
    p = end;
    return v;
  }
  void expect(const char *word) {
    skip();
    const std::size_t n = std::strlen(word);
    if (std::strncmp(p, word, n) != 0) {
      throw std::runtime_error(std::string("'") + word + "' expected");
    }
    p += n;
  }
};

// hash of a sequence of integers
struct Hash {
  std::uint64_t h = 1469598103934665603ULL;
  void add(const std::uint64_t x) {
    for (int i = 0; i < 8; ++i) {
      h = (h ^ ((x >> (8 * i)) & 0xff)) * 1099511628211ULL;
    }
  }
};

// a part: counts, position of each vertex, shared entities of each group
struct Part {
  int rank = 0;
  long long elements = 0, boundary = 0, vertices = 0, edges = 0, errors = 0;
  // sum of the hashes of the elements (attribute and positions of the
  // vertices): the same for all the partitions of a mesh
  std::uint64_t fingerprint = 0;
  std::vector<double> x;  // three coordinates per vertex
  std::vector<std::vector<long long>> members;
  std::vector<Shared> groups;
};

// box: lowest coordinates, opposites of the largest ones, then periods (the
// smallest distances between two positions of a vertex)
Part readPart(const std::string &file, const int rank, const int parts,
              double *box) {
  Part c;
  c.rank = rank;
  Reader r(file);
  if (r.text.empty()) {
    throw std::runtime_error("missing or empty");
  }
  r.expect("MFEM mesh v1.2");
  r.expect("dimension");
  if (r.integer() != 3) {
    throw std::runtime_error("meshes of dimension 3 only");
  }
  r.expect("elements");
  c.elements = r.integer();
  std::vector<int> elements(4 * c.elements), attributes(c.elements);
  std::vector<std::pair<int, int>> edges;
  std::vector<std::array<int, 3>> faces;
  for (long long i = 0; i < c.elements; ++i) {
    attributes[i] = r.integer();
    if (r.integer() != mfem::Geometry::TETRAHEDRON) {
      throw std::runtime_error("tetrahedra only");
    }
    int *v = &elements[4 * i];
    for (int k = 0; k < 4; ++k) {
      v[k] = r.integer();
    }
    for (int a = 0; a < 4; ++a) {
      for (int b = a + 1; b < 4; ++b) {
        edges.push_back(std::minmax(v[a], v[b]));
      }
      std::array<int, 3> f = {v[(a + 1) % 4], v[(a + 2) % 4], v[(a + 3) % 4]};
      std::sort(f.begin(), f.end());
      faces.push_back(f);
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  std::sort(faces.begin(), faces.end());
  r.expect("boundary");
  c.boundary = r.integer();
  for (long long i = 0; i < c.boundary; ++i) {
    r.integer();
    const long geom = r.integer();
    if (geom < 0 || geom >= mfem::Geometry::NumGeom) {
      throw std::runtime_error("invalid boundary element");
    }
    for (int k = 0; k < mfem::Geometry::NumVerts[geom]; ++k) {
      r.integer();
    }
  }
  r.expect("vertices");
  const long nv = r.integer();
  if (!edges.empty() && (edges.front().first < 0 || edges.back().first >= nv)) {
    throw std::runtime_error("vertex of an element out of range");
  }
  for (const auto &e : edges) {
    if (e.second >= nv) {
      throw std::runtime_error("vertex of an element out of range");
    }
  }
  c.x.assign(3 * nv, 0.);
  Hash element;
  auto add = [&](const double *p) {
    for (int d = 0; d < 3; ++d) {
      std::uint64_t bits;
      std::memcpy(&bits, &p[d], sizeof(bits));
      element.add(bits);
      box[d] = std::min(box[d], p[d]);
      box[3 + d] = std::min(box[3 + d], -p[d]);
    }
  };
  r.skip();
  if (std::strncmp(r.p, "nodes", 5) == 0) {
    // periodic mesh: positions of the vertices of each element (the mesh
    // read from gmsh by MFEM); a vertex takes the one of its first element
    r.expect("nodes");
    r.expect("FiniteElementSpace");
    r.expect("FiniteElementCollection: L2_T1_3D_P1");
    r.expect("VDim:");
    if (r.integer() != 3) {
      throw std::runtime_error("nodes: three components expected");
    }
    r.expect("Ordering:");
    const long ordering = r.integer();
    std::vector<double> nodes(12 * c.elements);
    for (auto &xi : nodes) {
      xi = r.real();
    }
    std::vector<char> seen(nv, 0);
    const long long n = 4 * c.elements;
    for (long long i = 0; i < n; ++i) {
      double p[3];
      for (int d = 0; d < 3; ++d) {
        p[d] = (ordering == 1) ? nodes[3 * i + d] : nodes[d * n + i];
      }
      if (i % 4 == 0) {
        element = Hash();
        element.add(attributes[i / 4]);
      }
      add(p);
      if (i % 4 == 3) {
        c.fingerprint += element.h;
      }
      const int v = elements[i];
      if (!seen[v]) {
        seen[v] = 1;
        std::copy(p, p + 3, &c.x[3 * v]);
        continue;
      }
      for (int d = 0; d < 3; ++d) {
        const double distance = std::abs(p[d] - c.x[3 * v + d]);
        if (distance > 1e-8 * (1 + std::abs(p[d]))) {
          box[6 + d] = std::min(box[6 + d], distance);
        }
      }
    }
  } else {
    if (r.integer() != 3) {
      throw std::runtime_error("vertices: three coordinates expected");
    }
    for (auto &xi : c.x) {
      xi = r.real();
    }
    for (long long i = 0; i < c.elements; ++i) {
      element = Hash();
      element.add(attributes[i]);
      for (int k = 0; k < 4; ++k) {
        add(&c.x[3 * elements[4 * i + k]]);
      }
      c.fingerprint += element.h;
    }
  }
  r.expect("mfem_serial_mesh_end");
  r.expect("communication_groups");
  r.expect("number_of_groups");
  const long ng = r.integer();
  c.members.resize(ng);
  c.groups.resize(ng);
  for (auto &m : c.members) {
    m.resize(r.integer());
    for (auto &q : m) {
      q = r.integer();
      if (q < 0 || q >= parts) {
        throw std::runtime_error("part of a group out of range");
      }
    }
    std::sort(m.begin(), m.end());
    if (!std::binary_search(m.begin(), m.end(), rank)) {
      throw std::runtime_error("group without the part itself");
    }
  }
  if (ng < 1 || c.members[0].size() != 1) {
    throw std::runtime_error("first group: the part alone expected");
  }
  r.expect("total_shared_vertices");
  const long total_vertices = r.integer();
  r.expect("total_shared_edges");
  const long total_edges = r.integer();
  r.expect("total_shared_faces");
  const long total_faces = r.integer();
  c.vertices = nv;
  c.edges = edges.size();
  long sum_vertices = 0, sum_edges = 0, sum_faces = 0;
  auto vertex = [&]() {
    const long v = r.integer();
    if (v < 0 || v >= nv) {
      throw std::runtime_error("shared vertex out of range");
    }
    return static_cast<int>(v);
  };
  for (long g = 1; g < ng; ++g) {
    Shared &s = c.groups[g];
    r.expect("shared_vertices");
    const long gv = r.integer();
    for (long i = 0; i < gv; ++i) {
      s.vertices.push_back(vertex());
    }
    r.expect("shared_edges");
    const long ge = r.integer();
    for (long i = 0; i < ge; ++i) {
      const int a = vertex(), b = vertex();
      s.edges.push_back(a);
      s.edges.push_back(b);
      // a shared edge is an edge of the local mesh
      const std::pair<int, int> e = std::minmax(a, b);
      if (!std::binary_search(edges.begin(), edges.end(), e)) {
        ++c.errors;
      }
    }
    r.expect("shared_faces");
    const long gf = r.integer();
    for (long i = 0; i < gf; ++i) {
      if (r.integer() != mfem::Geometry::TRIANGLE) {
        throw std::runtime_error("shared triangles only");
      }
      std::array<int, 3> f;
      for (int k = 0; k < 3; ++k) {
        f[k] = vertex();
        s.triangles.push_back(f[k]);
      }
      // a shared triangle is a face of one local tetrahedron only
      std::sort(f.begin(), f.end());
      const auto range = std::equal_range(faces.begin(), faces.end(), f);
      if (range.second - range.first != 1) {
        ++c.errors;
      }
    }
    sum_vertices += gv;
    sum_edges += ge;
    sum_faces += gf;
    // the entities of a group are counted by its part of lowest rank
    if (c.members[g][0] != rank) {
      c.vertices -= gv;
      c.edges -= ge;
    }
  }
  r.expect("mfem_mesh_end");
  if (sum_vertices != total_vertices || sum_edges != total_edges ||
      sum_faces != total_faces) {
    throw std::runtime_error("totals of the shared entities");
  }
  return c;
}

// Records of the groups of a part: size, parts, numbers of shared vertices,
// edges and triangles, and two hashes of their positions in the order of the
// file (the vertices of an edge or of a triangle in any order). The
// positions are taken on a grid of 2^26 points per period of the mesh,
// modulo the period (per size of the box in a direction without period): a
// vertex of a periodic mesh has several positions, which differ by periods
// and by the rounding. A position next to the limit of two cells of the
// grid may fall in either: the second hash uses the grid shifted by half a
// cell, and two descriptions agree if one hash agrees.
void records(const Part &c, const double *box, std::vector<long long> &out) {
  const long long cells = 1LL << 26;
  auto point = [&](const int v, const double shift) {
    Hash h;
    for (int d = 0; d < 3; ++d) {
      const bool periodic = (box[6 + d] < std::numeric_limits<double>::max());
      const double size = periodic ? box[6 + d] : -box[3 + d] - box[d];
      long long q = 0;
      if (size > 0) {
        q = static_cast<long long>(
            std::floor((c.x[3 * v + d] - box[d]) / size * cells + shift));
        if (periodic) {
          q = ((q % cells) + cells) % cells;
        }
      }
      h.add(static_cast<std::uint64_t>(q));
    }
    return h.h;
  };
  for (std::size_t g = 1; g < c.groups.size(); ++g) {
    const Shared &s = c.groups[g];
    out.push_back(c.members[g].size());
    out.insert(out.end(), c.members[g].begin(), c.members[g].end());
    out.push_back(s.vertices.size());
    out.push_back(s.edges.size() / 2);
    out.push_back(s.triangles.size() / 3);
    for (const double shift : {0.5, 0.}) {
      Hash h;
      auto entities = [&](const std::vector<int> &list, const int n) {
        std::uint64_t e[3];
        for (std::size_t i = 0; i < list.size(); i += n) {
          for (int k = 0; k < n; ++k) {
            e[k] = point(list[i + k], shift);
          }
          std::sort(e, e + n);
          for (int k = 0; k < n; ++k) {
            h.add(e[k]);
          }
        }
      };
      entities(s.vertices, 1);
      entities(s.edges, 2);
      entities(s.triangles, 3);
      out.push_back(static_cast<long long>(h.h));
    }
  }
}

int check(const int parts, const std::string &prefix) {
  const int rank = mfem::Mpi::WorldRank(), workers = mfem::Mpi::WorldSize();
  std::vector<Part> read;
  long long counts[5] = {0, 0, 0, 0, 0}, totals[5];
  unsigned long long fingerprint = 0, global_fingerprint = 0;
  double box[9], global_box[9];
  std::fill(box, box + 9, std::numeric_limits<double>::max());
  for (int p = rank; p < parts; p += workers) {
    const std::string file = mfem::MakeParFilename(prefix, p);
    try {
      Part c = readPart(file, p, parts, box);
      counts[0] += c.elements;
      counts[1] += c.boundary;
      counts[2] += c.vertices;
      counts[3] += c.edges;
      counts[4] += c.errors;
      fingerprint += c.fingerprint;
      if (c.errors != 0) {
        std::cerr << file << ": " << c.errors
                  << " shared edges or triangles absent from the local mesh\n";
      }
      read.push_back(std::move(c));
    } catch (const std::exception &e) {
      std::cerr << file << ": " << e.what() << '\n';
      ++counts[4];
    }
  }
  MPI_Allreduce(counts, totals, 5, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(box, global_box, 9, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&fingerprint, &global_fingerprint, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, MPI_COMM_WORLD);
  std::vector<long long> mine;
  for (const Part &c : read) {
    records(c, global_box, mine);
  }
  // records of all the workers on the first one
  int size = mine.size();
  std::vector<int> sizes(workers), offsets(workers + 1, 0);
  MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  for (int w = 0; w < workers; ++w) {
    offsets[w + 1] = offsets[w] + sizes[w];
  }
  std::vector<long long> all(mfem::Mpi::Root() ? offsets[workers] : 0);
  MPI_Gatherv(mine.data(), size, MPI_LONG_LONG, all.data(), sizes.data(),
              offsets.data(), MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  if (!mfem::Mpi::Root()) {
    return EXIT_SUCCESS;
  }
  // each group: as many descriptions as parts, which agree
  std::map<std::vector<long long>, std::vector<std::array<long long, 5>>>
      groups;
  for (std::size_t i = 0; i < all.size();) {
    const long long n = all[i];
    const std::vector<long long> key(all.begin() + i + 1,
                                     all.begin() + i + 1 + n);
    i += 1 + n;
    groups[key].push_back(
        {all[i], all[i + 1], all[i + 2], all[i + 3], all[i + 4]});
    i += 5;
  }
  long long inconsistent = 0, shared_faces = 0;
  for (const auto &g : groups) {
    const auto &first = g.second.front();
    bool ok = (g.second.size() == g.first.size());
    for (const auto &d : g.second) {
      ok = ok && d[0] == first[0] && d[1] == first[1] && d[2] == first[2] &&
           (d[3] == first[3] || d[4] == first[4]);
    }
    if (!ok) {
      if (++inconsistent <= 10) {
        std::cerr << "group";
        for (const auto q : g.first) {
          std::cerr << ' ' << q;
        }
        std::cerr << ": " << g.second.size()
                  << " descriptions, different or missing\n";
      }
    } else {
      shared_faces += first[2];
    }
  }
  std::cout.precision(17);
  std::cout << prefix << ": periods";
  for (int d = 0; d < 3; ++d) {
    if (global_box[6 + d] < std::numeric_limits<double>::max()) {
      std::cout << ' ' << global_box[6 + d];
    } else {
      std::cout << " none";
    }
  }
  std::cout << ", fingerprint of the elements " << std::hex
            << global_fingerprint << std::dec << '\n'
            << prefix << ": " << parts << " parts, " << totals[0]
            << " elements, " << totals[1] << " boundary elements, "
            << totals[2] << " vertices, " << totals[3] << " edges, "
            << groups.size() << " groups, " << shared_faces
            << " shared triangles, " << totals[4] << " errors in the parts, "
            << inconsistent << " inconsistent groups\n";
  if (totals[4] != 0 || inconsistent != 0) {
    std::cout << prefix << ": check FAILED\n";
    return EXIT_FAILURE;
  }
  std::cout << prefix << ": check passed\n";
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char *argv[]) {
  mfem::Mpi::Init(argc, argv);
  const std::string mode = (argc > 1) ? argv[1] : "";
  if (mode == "--check" && argc == 4) {
    return check(std::stoi(argv[2]), argv[3]);
  }
  if (mode == "--mfem" && argc == 5) {
    return split(argv[2], std::stoi(argv[3]), mfem::Mpi::WorldSize(), argv[4],
                 true);
  }
  if (mode.rfind("--", 0) != 0 && argc == 5) {
    return split(argv[1], std::stoi(argv[2]), std::stoi(argv[3]), argv[4],
                 false);
  }
  if (mfem::Mpi::Root()) {
    std::cerr
        << "usage: split-mesh <mesh> <refinements> <parts> <prefix>\n"
        << "       split-mesh --mfem <mesh> <refinements> <prefix>\n"
        << "       split-mesh --check <parts> <prefix>\n";
  }
  return EXIT_FAILURE;
}
