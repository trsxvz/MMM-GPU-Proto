# Cost of the post-processings of the polycrystal, work package L8 (13/10/2026)

Polycrystal of mm-opera-hpc (5 grains, single crystal of UO2 in finite
strain, 20,553 unknowns), 10 steps, all the computations by kernels,
Newton, on the GPU (one process) and on 8 CPU processes, without and then
with the post-processings of the case (ParaView export of the results and
averages of the stresses at each step), computed on the host and then by
kernels (`MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels`). mfem-mgis
54ff8e6c8a (separate installation `install/mfem-mgis-54ff8e6c8a`, build
`build/mm-opera-hpc-54ff8e6c8a`). Script:
`cases/measure-post-processings.sh`.

Total time (seconds), and time of the post-processings
(`NLEP::PP::execute`):

| | Without | On the host | By kernels |
|---|---|---|---|
| GPU, total | 11.00 | 12.48 | 12.53 |
| GPU, post-processings | | 1.48 | 1.53 |
| 8 processes, total | 28.68 | 30.28 | 30.02 |
| 8 processes, post-processings (average) | | 0.19 | 0.18 |

- The post-processings cost 1.5 s on the GPU (12 % of the computation),
  and gain nothing from being computed by kernels: this case mostly does
  the ParaView export (writing of the files on the host), whereas the
  bubble, where the kernels divide the total time by 2.8
  (`results/post-processings-kernels-2026-10-11`), mostly does L2
  projections and values at the nodes.
- The copy of the materials back to the host before the post-processings
  (`NLEP::PP::synchronizeOnHost`) costs 0.02 s.
- Same macroscopic stresses in all the variants (relative gap at most
  2.4e-13 on the GPU, zero on 8 processes).
