# Meshes of the polycrystal for the series in number of grains

Periodic RVEs of UO2 for the polycrystal of mm-opera-hpc (13/11/2026,
measurements on Adastra, see `JOURNAL.md`). The meshes themselves are not in
the repository (570 MB for 6,400 grains): `make-mesh.sh` builds them again.

- **Geometry.** Laguerre tessellation of Merope, with germs from a saturated
  RSA packing as in `mm-opera-hpc/polycrystal/mesh/1000crystals.py`. The
  mean volume of a grain is 1000, so the cube has side (1000 N)^(1/3).
  - The RVE of Merope is not a box: it is a union of whole grains, and its
    periodic boundary follows the grains.
- **Mesh.**
  - gmsh 4.15.2 (msh 2.2), mesh size 2, about 6 elements across a grain
    before the uniform refinements of the computation.
  - P2 and refined once (r1): about 33,000 unknowns per grain.
- **Orientations.** `randomVectorGeneration` of mm-opera-hpc (algorithm
  `angles`, numpy seed N), two vectors per grain.
- **Short edges.** `clean_mesh.py` removes them (see below).

| Grains | Seed | Nodes | Tetrahedra | Max aspect before / after | Shortest edge before / after |
|---|---|---|---|---|---|
| 100 | 1 | 20,785 | 103,119 | 194 / 39 | 0.0113 / 0.321 |
| 400 | 1 | 77,857 | 413,155 | 1,638 / 48 | 0.0020 / 0.302 |
| 1,600 | 1 | 296,866 | 1,662,414 | 1,914 / 164 | 0.0013 / 0.157 |
| 6,400 | 2 | 1,138,980 | 6,597,793 | 71,090 / 198 | 4.2e-5 / 0.057 |

The aspect is the longest edge over the smallest height of a tetrahedron.

For 6,400 grains, the seed 1 gave a geometry that gmsh could not build: 3
plane surfaces were not created and 4 volumes were left without elements.

## Plane surfaces that are not planar

The geometries of many grains often have a few plane surfaces that
OpenCASCADE refuses ("Could not add plane surface"), leaving the volumes
they bound without elements (14/11/2026): seeds 1 to 5 of 12,800 grains,
1 to 8 of 38,400 grains.

- **Cause.** When a face of the Laguerre tessellation is tiny and seen by
  one cell only, Merope merges its vertices into one point (warning "An
  entire face should be merged", without checks). The faces around this
  point are no longer planar: their vertices are up to 1e-4 from their
  plane (edges of 1 to 8), above the tolerance of OpenCASCADE (1e-7).
  - The number of these warnings grows with the number of grains: 0 to 5
    for 12,800 grains, 8 to 16 for 38,400 grains (8 seeds each).
- **Fix.** `fix_plane_surfaces.py in.geo out.geo` defines these surfaces
  as `Surface` instead of `Plane Surface`: OpenCASCADE builds a surface
  through the loop, as flat as the loop. The other lines are copied.
  - A plane surface is rewritten when one of its vertices is farther than
    1e-8 from the least squares plane of the vertices of its loop. This
    needs no run of gmsh: 11 s for the 301,522 plane surfaces of 38,400
    grains (17 rewritten, largest distance 7.75e-5), against 1712 s for
    a parsing of the geometry of 12,800 grains by gmsh.
  - `adastra-fm/mesh.slurm` applies it before gmsh.
- **Check.** Seed 1 of 12,800 grains, refused before (4 plane surfaces):
  with the 4 surfaces rewritten, gmsh gives no error and no volume without
  elements, 32 warnings on ill-shaped tetrahedra (42 for the seed 6,
  meshed as is).

## Short edges of the tessellation and flat tetrahedra

The Laguerre tessellation has very short edges, and their length decreases
with the number of grains. gmsh meshes them by edges of the same length,
between tetrahedra of the size of the mesh, which gives flat tetrahedra.

In these elements, the small errors of the iterates of the inexact Newton
method give large strains: the behaviour fails at a few integration points,
and the whole step is then divided. On 100 grains (r1, 5 steps of 1 s, 4
MI300A):

| Mesh | Steps | Resolutions | Newton iterations | Failed resolutions | Messages of failed integration | Total time |
|---|---|---|---|---|---|---|
| Raw | 16 | 35 | 338 | 8 | 103 | 105.3 s |
| Cleaned | 5 | 16 | 97 | 0 | 0 | 36.4 s |

The axial stress at 5 s is 74.44 MPa on the raw mesh and 74.32 MPa on the
cleaned one.

## How `clean_mesh.py` works

`clean_mesh.py in.msh eps out.msh` collapses the edges shorter than eps (0.3,
15 % of the mesh size), as the regularization of Neper removes the short
edges of a tessellation.

1. **Merging.** The nodes of the short edges are merged, together with
   their periodic images.
2. **Position.** A merged class takes the position of its most constrained
   node: first the periodic boundary (number of distinct translations to
   its images), then the grains (quadruple point, triple line, grain
   boundary).
3. **Exact periodicity.** The periodic nodes are put exactly at the position
   of their master plus the translation.
4. **Giving up bad collapses.** A collapse that reverses an element, crushes
   it (more than 99 % of its volume lost) or makes its aspect much worse is
   given up, with its periodic images, and the collapse is done again.
5. **Result.** The total volume is kept and the volume of the grains changes
   by 9e-4 at most (100 grains, `grain_volumes.py`). The nodes are
   renumbered.

## Partitioned meshes and binary meshes

In the mode "FromScratch" of mfem-mgis, each process reads the whole mesh,
builds the parallel mesh (METIS) and refines it: 74 s on 64 MI300A and
104 s on 1536 cores for 6,400 grains (14/11/2026). `partition/` builds the
same refined parallel mesh once (`partition-mesh`, MFEM, same calls as
mfem-mgis) and writes the part of each process (`ParPrint`, 17 digits),
read by the computation with `--mesh_mode Restart` and no refinement
(`adastra-fm/partition.slurm`, `FM_MESH_GPU` and `FM_MESH_CPU` of
`adastra-fm/series.slurm`): 7.7 s and 1.2 s for 6,400 grains. The
vertices are numbered again at the reading: the axial stress differs by
7e-9 relatively from the one of "FromScratch" (100 grains, CPU).

`msh22_binary.py` writes the same mesh in binary Gmsh 2.2 (read by MFEM
4.10, identical outputs on the CPU), but the reading is hardly faster
(2.67 s against 2.29 s for 400 grains): the cost is the construction of
the mesh, not the parsing of the file.
