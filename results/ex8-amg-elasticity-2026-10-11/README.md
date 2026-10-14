# ex8: BoomerAMG with the elasticity options (11/10/2026)

Example ex8 of mfem-mgis-examples (thermomechanical plate), original mesh,
5 steps, 2 CPU processes, everything by kernels for the mechanics. Build
`build/examples-post-ex8` (ex8 with the options `-oMc` and `-amgMc`, linked
to `install/mfem-mgis-post`).

| Solver of the mechanics | Coupling over the 5 steps (s) |
|---|---|
| MUMPS (reference) | 23.52 |
| PCG + BoomerAMG, elasticity options, byVDIM ordering (`-lsMc HyprePCG -pcMc HypreBoomerAMG -oMc byVDIM -amgMc Elasticity`) | 201.16 |

- Same statistics of the fields (temperatures and norms of the
  displacements) in both computations.
- The elasticity options (interpolation by rigid modes) do worse than the
  systems options tried before (about 240 iterations, 93 s, measurement of
  the afternoon of 11/10 noted in `JOURNAL.md`, with no raw output kept):
  441 PCG iterations on average per solution (29 solutions). The plate is
  thin and its materials are contrasted: BoomerAMG does not fit as is (the
  first level coarsens little: 28,548 then 22,482 rows, operator
  complexity 3.6; the explanation by the geometry and the materials is not
  verified). MUMPS remains the choice for the mechanics of ex8.
