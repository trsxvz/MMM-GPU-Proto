# Ordered assembly of the tangent (12 and 13/10/2026)

Parameter `TangentOperatorAssembly` of mfem-mgis (`Atomic` by default,
`Ordered`), variable `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`. With a device,
the stiffness matrices of the elements are added to the sparse matrix by
atomic additions, in an order that changes from one run to the next: the
last digits of the results change too, and the number of Newton iterations
may vary (32 or 31 tangents below with the atomic assembly). The ordered
assembly stores the matrices of the elements, then adds them in a fixed
order. Without a device, the host path already adds the elements one after
the other: both values of the parameter give the same assembly there.

Beam r4 of the prototype (111,843 unknowns, von Mises with linear
hardening, 10 steps), everything by kernels, cuda behaviour, one process on
the GPU (cores 0 to 3), hypre on the GPU. Separate installation
`install/mfem-mgis-ordered`, build `build/cases-ordered` (only the `beam`
executable). Commands in `measure.sh`.

## Two versions of the ordered assembly on the GPU

- First version (mfem-mgis 456c298dd6): one thread per row of the matrix,
  which looks for the column of each of its entries in the row. Nsight
  Systems profile: 7.1 ms per assembly for this kernel, against 0.48 ms for
  the kernel of the atomic additions. With a single thread per row
  (111,843 threads), each thread does in sequence nearly 200 searches and
  scattered reads.
- Version kept (mfem-mgis f57072220f): one thread per value of the matrix
  (about 9 million threads), which adds in a fixed order the two to three
  contributions to this value. The list of the contributions is computed
  once per layout of the rows, from the positions of the entries (existing
  kernel), then checked column by column at each assembly. Same values as
  the first version.

## Times (seconds)

Total time of the loop over the steps, two runs per variant; cumulated
times of the assembly of the tangent (`NS::getJacobian`) and of the
assembly of the prediction (`PRED::assemble`, which contains the first
assembly, hence the computation of the contributions).

| Method | Assembly | Total | Tangent | Prediction |
|---|---|---|---|---|
| Newton, PCG and BoomerAMG | atomic | 3.09 and 3.09 | 0.137 (32) and 0.146 (33) | 0.104 and 0.106 |
| Newton, PCG and BoomerAMG | ordered, by row | 3.20 and 3.20 | 0.354 (32) and 0.353 (32) | 0.178 and 0.175 |
| Newton, PCG and BoomerAMG | ordered, by value | 3.09 and 3.67 | 0.189 (32) and 0.212 (32) | 0.234 and 0.247 |
| NewtonK0, Eisenstat and Walker | atomic | 1.89 and 1.95 | 0.123 (51) and 0.140 (53) | 0.106 and 0.109 |
| NewtonK0, Eisenstat and Walker | ordered, by row | 2.40 and 2.61 | 0.493 (53) and 0.542 (53) | 0.173 and 0.183 |
| NewtonK0, Eisenstat and Walker | ordered, by value | 2.22 and 2.35 | 0.218 (53) and 0.234 (53) | 0.235 and 0.244 |

In parentheses, the number of assemblies of the tangent. The total times
vary by about 10% from one run to the next (other sessions on the
machine); the times per assembly are more stable.

- Per assembly of the tangent: 4.3 ms with atomics, 11 ms by row, 6 ms by
  value. The first assembly costs 0.13 s more by value (computation of the
  contributions, once per layout of the rows; 0.45 s when the entries were
  looked for on the host).
- Extra memory with a device: the matrices of the elements (8 n² bytes per
  element), the contributions (4 n² bytes per element, plus 4 bytes per
  value of the matrix) and their positions, during their computation (2 n²
  bytes per element).

## CPU

Same beam on 4 processes (cores 0 to 3), without device, openmp behaviour,
Newton, PCG and BoomerAMG. Alternated runs: atomic 17.96 and 16.36 s,
ordered 17.84 and 16.71 s (tangent 0.76 and 0.70 s for both). The first
version of the ordered assembly, which also stored the matrices on the
host, took 22.2 s (tangent 4.1 to 4.3 s, without the positions of the
entries): files `beam-cpu4-newton-Ordered-rows-*`.

## Reproducibility

`FixedPointSolverTest` on the GPU, kernels forced, ordered assembly: four
runs give bitwise identical iteration histories, with both versions; with
the atomic assembly, three runs give three different histories.
