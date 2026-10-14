# kokkos backend on the polycrystal (11/10/2026)

Script: `cases/measure-backends-polycrystal.sh`, with the variables
`polycrystal=build/mm-opera-hpc-kokkos/polycrystal`,
`cases=build/cases-kokkos` and `backends="cuda:cuda kokkos:kokkos"`.

Polycrystal of mm-opera-hpc (5 grains, single crystal behaviour of UO2,
orthotropic, finite strain, implicit), 10 steps, one process and the GPU,
everything by kernels. Same executable for all the computations: the one
of mm-opera-hpc built against the installation `install/mfem-mgis-kokkos`
(kokkos branch of mfem-mgis, aed33abb8a, option `enable-kokkos`). The line
"mfem-mgis" of `summary.txt` gives the version of the submodule, not that
of this installation.

Behaviours compared:

- `generic`: generic interface, integrated on the host (reference of the
  macroscopic stresses);
- cuda: generic-parallel interface, configuration `config-cuda-nvcc.json`;
- kokkos: generic-parallel interface, configuration `config-kokkos.json`
  (Kokkos 5.2.2 built for CUDA, compiled by `nvcc_wrapper`). Kokkos is
  initialised by mfem-mgis (variable `MFEM_MGIS_KOKKOS=1`), on the GPU of
  the device of MFEM.

Quiet machine (web browser only), two runs.

Macroscopic stresses: relative gap at most 2.8e-13 with respect to the
first computation with the `generic` behaviour (70 values). Same Newton
iterations (98) in the six computations.

| Library | Total (s) | Integration, 158 calls (s) | Solution (s) |
|---|---|---|---|
| `generic` (host) | 107.65 and 107.48 | 99.97 and 99.84 | 5.73 and 5.66 |
| cuda (nvcc) | 12.88 and 12.79 | 5.61 and 5.61 | 5.64 and 5.56 |
| kokkos (nvcc) | 12.97 and 13.22 | 5.62 and 5.69 | 5.70 and 5.88 |

- kokkos and cuda give the same integration time, within 1.4%. Both
  toolchains go through nvcc. As for cuda and hip
  (`results/backends-polycrystal-2026-10-11`), the compiler is what
  matters, not the programming model.
- The gap on the total (0.1 and 0.4 s) comes mostly from the solution,
  which does not depend on the behaviour. The initialisation of Kokkos,
  done by `mfem_mgis::initialize`, precedes the solution of the problem
  and is not counted in the total.
- The program had to initialise Kokkos. This is now the role of
  mfem-mgis: no change of mm-opera-hpc was needed.
