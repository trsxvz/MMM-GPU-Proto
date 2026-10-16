# MueLu: factories written with Kokkos and rigid body modes (16/10/2026)

With the rigid body modes, MueLu needed about twice the iterations with its
factories written with Kokkos (`use kokkos refactor` set to true, the
default of MueLu on all the nodes of Tpetra but the serial one) as with
the other ones: 22 against 12 on the linear elasticity of 421,443
unknowns on the GPU, 78 against 11 per resolution on the nonlinear beam
(`../beam-muelu-2026-10-16`). This directory finds the cause: the local
QR factorization of the tentative prolongator written with Kokkos
(`TentativePFactory_kokkos` of Trilinos 17.2.1) does not compute the Q
factor of the nullspace when the nullspace has more than one vector.

`cases/measure-muelu-factories.sh`: linear elasticity of the benchmark
`LinearSolversBenchmark` of mfem-mgis, conjugate gradient of MFEM
preconditioned by MueLu with a Jacobi smoother (two sweeps damped by 0.6),
relative tolerance 1e-8, one process: on the GPU at 421,443 unknowns
(`gpu-*`), on the host with the serial node of Tpetra at 56,355 unknowns
(`host-*`). The variants add parameters of MueLu (`summary.txt`). The logs
`stats-*` were run again with `VERBOSITY=medium` (aggregates and levels).
The logs come from three runs of the same commands (GPU reserved), whose
iterations agree to one iteration.

| variant | parameters | GPU | host |
|---|---|---|---|
| `classic` | `use kokkos refactor` false | 12 | 15 |
| `kokkos` | `use kokkos refactor` true | 22 to 23 | 25 |
| `kokkos-host-aggregation` | the same, aggregation on the host | 22 | 23 |
| `classic-kokkos-aggregation` | false, aggregation written with Kokkos | 12 | 22 |
| `kokkos-random`, `kokkos-phase1-serial` | true, other aggregation options | 22 to 23 | 23 to 25 |
| `*-unfiltered`, `kokkos-no-lumping` | filtered matrix not used, or not lumped | unchanged | unchanged |
| `classic-noqr`, `kokkos-host-aggregation-noqr` | no local QR (`tentative: calculate qr` false) | the same results for both (not converged) | the same results for both |

- On the GPU, the factories written with Kokkos with the aggregation of
  the host build the same hierarchy as the other factories: same number of
  aggregates (5,203 then 240), same sizes and numbers of nonzeros of the
  coarse operators, same damping factor of the prolongator on the finest
  level (0.500648). Only the values of the coarse operators differ (largest
  eigenvalue of the diagonally scaled operator of the first coarse level
  1.80 against 1.51), and the iterations double.
- Without the local QR, both give exactly the same resolution: the
  dropping and the aggregation are not the cause, the local QR is.
- `qr-check.py` replays the local QR of `LocalQRDecompFunctor`
  (`MueLu_TentativePFactory_kokkos_def.hpp` of Trilinos 17.2.1) on the
  rigid body modes of random aggregates: the Householder reflections are
  applied to Q^T on its trailing block only (`for (int j = k; j < m; j++)`),
  in increasing order, and Q^T is used as Q (its transposition is under
  `#if 0`): the result is not orthonormal and `Q R` differs from the
  nullspace by 19 to 51 percent. The same reflections applied to all the
  columns of Q^T, or accumulated in decreasing order, give the
  factorization to the rounding errors. The case of a nullspace of one
  vector (Laplace equation) has its own code, which is correct. An
  independent check compiled the functor of 17.2.1 unchanged against a
  minimal mock of Kokkos (not versioned): `P R` differs from the rigid
  body modes by 27 to 100 percent, and applying the reflections to all the
  columns of Q^T and transposing it gives the factorization to 1e-15.
- The same check found a second defect, also in Trilinos develop: for an
  aggregate with fewer unknowns than nullspace vectors (an aggregate of one
  node with the rigid body modes: 3 unknowns, 6 vectors), the functor
  writes and reads `q(i, j)` for `j < n` in a view of `m` by `m` entries
  (`m < n`). With the layout of the nodes of the host (`LayoutRight`:
  serial, OpenMP), the entries out of the view overlap the next rows, and
  the rows of the tentative prolongator of the node are wrong (`[1 0 0 0 1
  0]` instead of `[1 0 0 0 0 0]`); with the layout of CUDA, they land
  after the view, and the result is right. mfem-mgis does not use these
  factories by default.
- Trilinos develop replaced this QR by the batched QR of Kokkos Kernels
  (`KokkosBatched::SerialQR` and `SerialApplyQ`, file `MueLu_LocalQR.hpp`,
  pull request 13507 merged on 08/09/2026, presented as a refactoring),
  after the release 17.2.1 (31/08/2026). mfem-mgis sets `use kokkos
  refactor` to false by default, which avoids the defect with the releases
  up to 17.2.1.
- On the host at 56,355 unknowns, the aggregation written with Kokkos alone
  also increases the iterations (22 against 15, same number of aggregates),
  while on the GPU at 421,443 unknowns it gives the same resolution as the
  aggregation of the host. Not investigated (different sizes).
