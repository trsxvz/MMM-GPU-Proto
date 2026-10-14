# ex8: symmetric matrices given to MUMPS (11/10/2026)

Example ex8 of mfem-mgis-examples (thermomechanical plate), original mesh,
5 steps, MUMPS for the mechanics, everything by kernels. Option `-msMc` of
ex8 (mfem-mgis-examples 0782d73, then 03c347e): symmetry of the matrices of
the mechanics given to MUMPS (parameters `Symmetric` and `PositiveDefinite`
of the MUMPS solver of mfem-mgis). Build `build/examples-post-ex8`, linked
to `install/mfem-mgis-post` (standard host memory). Quiet machine.

The tangent of the mechanics of ex8 is symmetric: same Newton iterations
(29 solutions of the mechanics, hence 29 factorisations, and 35 of the
thermal) and same statistics of the fields in all cases.

## CPU, 2 processes, `generic` behaviours (one run)

| Symmetry | Total (s) | Factorisation (s) |
|---|---|---|
| `Unsymmetric` (default) | 22.05 | 7.49 |
| `SymmetricIndefinite` | 20.54 | 5.99 |
| `SymmetricPositiveDefinite` | 19.85 | 5.31 |

## GPU, 1 process, cuda behaviours (two runs)

| Symmetry | Total (s) | Factorisation (s) | Solutions, two physics (s) |
|---|---|---|---|
| `Unsymmetric` (default) | 16.50 and 16.55 | 11.37 and 11.41 | 3.54 and 3.57 |
| `SymmetricPositiveDefinite` | 12.34 and 12.48 | 7.30 and 7.33 | 3.47 and 3.56 |

- On the GPU, ex8 goes from 16.5 to 12.4 s. MUMPS (on the host) remains
  the bulk: 7.3 s of factorisations and 0.31 s of solutions (29) out of
  12.4 s. The 3.5 s of the solutions column are mostly those of the
  thermal: GMRES and BoomerAMG on the GPU, 3.15 s for 35 solutions (about
  1.5 s on 2 CPU processes), the second item after MUMPS.
- With the measurements of 11/10 (`results/ex8-behaviours-cuda-2026-10-11`,
  18.5 s) and the fix of the host memory
  (`results/host-memory-2026-10-11`, 16.9 s), ex8 entirely on the GPU goes
  down to 12.4 s.
