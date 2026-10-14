# cuDSS on the polycrystal (12/10/2026)

Polycrystal of mm-opera-hpc (5 grains, single crystal behaviour of UO2 in
finite strain, periodic cell, uniaxial tension, 20,553 unknowns, linear
elements, 10 steps), everything by kernels, cuda behaviour, one process on
the GPU (cores 0 to 3). mfem-mgis 429f76ded3 (`install/mfem-mgis-mumps`),
build `build/mm-opera-hpc-cudss` (only `uniaxial-polycrystal`). Option
`--reuse-analysis` and direct solvers `MUMPSSolver` and `CuDSSSolver`
added to the option `--linear-solver` of the polycrystal (non symmetric
matrices). Two passes (commands in `measure.sh`, macroscopic stresses in
`*.res`). The load of the machine was 3.3 in the first pass and 2 in the
second (another session was debugging on cores 12 to 15). The two passes
agree within 1 %. Total times, in seconds:

| Method | Solver | Time |
|---|---|---|
| Newton | HyprePCG and BoomerAMG (hypre on the GPU) | 11.46 and 11.44 |
| NewtonK0, Eisenstat and Walker tolerance | conjugate gradient and BoomerAMG built on K0 | 8.78 and 8.76 |
| Newton | cuDSS, analysis reused, hypre on the GPU | 19.12 and 19.11 |
| Newton | cuDSS, analysis reused, hypre on the host | 19.30 and 19.26 |

- cuDSS does not pay off on this case: one factorisation takes 132 ms
  (12.9 s for 98), against 57 ms for one solve by PCG and BoomerAMG (5.6 s
  for 98). The reuse of the analysis works (short trial over two steps:
  134 ms per factorisation with it, 222 ms without), but the numerical
  factorisation of this matrix is expensive: the periodic cell links the
  opposite faces, and the matrix is non symmetric.
- At a similar size, the mechanical matrix of ex8 (28,548 unknowns, non
  periodic, symmetric positive definite) is factorised in 24 ms
  (`results/cudss-2026-10-12`).
- Same macroscopic stresses with cuDSS and with PCG (relative gap
  2.6e-9). NewtonK0 with Eisenstat and Walker departs from Newton by about
  1e-6 on the axial stress (a few pascals on the transverse stresses,
  which are almost zero).
- The best variant stays NewtonK0 with the Eisenstat and Walker tolerance
  (8.8 s; 9.6 s on 10/10 with mfem-mgis 4b72095ee0,
  `results/L5-polycrystal-2026-10-10-0502`).
