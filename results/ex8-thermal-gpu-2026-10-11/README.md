# ex8: solver of the thermal on the GPU (11/10/2026)

Example ex8 of mfem-mgis-examples (thermomechanical plate), original mesh
(9,516 unknowns for the thermal, 28,548 for the mechanics), 5 steps,
everything by kernels. Mechanics: MUMPS with
`-msMc SymmetricPositiveDefinite` (`results/ex8-mumps-symmetric-2026-10-11`).
Build `build/examples-post-ex8`, linked to `install/mfem-mgis-post`.

After the switch of the mechanics to symmetric MUMPS, the thermal is the
second item of ex8 on the GPU. Variants of the solver of the thermal
(`-lsTh`, `-pcTh`): GMRES of hypre (default) with BoomerAMG (default),
DiagScale, ILU or ParaSails, and MUMPS on the host.

## GPU, 1 process, cuda behaviours (one run)

| Thermal | Total (s) | Thermal (s) | Solutions (35) | Setup |
|---|---|---|---|---|
| GMRES and BoomerAMG (default) | 12.62 | 4.26 | 3.49 | 0.003 |
| GMRES and DiagScale | 12.80 | 4.56 | 3.83 | 0.000 |
| MUMPS (`-lsTh MUMPSSolver`) | **11.50** | **3.03** | 0.06 | 2.17 (factorisations) |
| GMRES and ILU | stopped at 600 s | | | |
| GMRES and ParaSails | crash | | | |

The mechanics does not change (8.2 to 8.4 s, including 7.3 s of
factorisations). With MUMPS, BoomerAMG or DiagScale, same Newton iterations
and same statistics of the fields.

- MUMPS on the host is the fastest for the thermal of ex8 on the GPU: 4.26
  to 3.03 s for the thermal, 12.6 to 11.5 s for ex8. The 35
  factorisations (2.17 s) cost less than the 35 solutions of GMRES on the
  GPU (3.49 s).
- On the GPU, one solution of GMRES and BoomerAMG takes 100 ms for 9,516
  unknowns, against 42 ms on 2 CPU processes. DiagScale, with no setup,
  takes almost as long as BoomerAMG. The system is probably too small for
  the GPU: the time would come mostly from the fixed cost of each
  iteration (kernel launches, reductions); the number of iterations was not
  recorded.
- ILU: on the GPU, mfem-mgis imposes ILU(0) (the higher fill levels crash
  in hypre, `reproducers/hypre-ilu-gpu`). GMRES does not converge from the
  second Newton iteration of the first step on. The time step is then
  split, and after 600 s the computation is only at 12,700 s out of the
  20,000 s of the first step. On the CPU (2 processes, block ILU(1), the
  default of MFEM), ILU converges (next table). These attempts do not tell
  whether the failure comes from level 0 or from the GPU.
- ParaSails: two defects, fixed in mfem-mgis (branch `parasails`,
  607f895a0d, details below).

## CPU, 2 processes, generic behaviours (one run, cores 0 to 3)

| Thermal | Thermal (s) | Solutions (35) | Newton iterations |
|---|---|---|---|
| GMRES and BoomerAMG (default) | 6.04 | 1.49 | 114 |
| GMRES and block ILU(1) | 8.09 | 3.53 | 114 |
| GMRES and ParaSails, before the fix | crash at the 2nd solution | | |
| GMRES and ParaSails, after the fix | failure: time steps split down to the limit | | |

The totals (28 to 30 s) are not comparable to the 19.9 s of
`results/ex8-mumps-symmetric-2026-10-11`: these attempts compute element by
element (without the variables of mfem-mgis that select the batches and
the kernels), while that measurement computes everything by kernels. The
solutions of the thermal (1.49 s) and the factorisations of the mechanics
(5.32 s) are the same. Only the columns of the thermal are used to compare
the preconditioners.

## ParaSails

1. On the CPU as on the GPU, `mfem::HypreParaSails::SetOperator` (MFEM
   4.10, also on the main branch of MFEM) creates the new preconditioner
   of hypre with the communicator of the previous operator. The nonlinear
   solvers of mfem-mgis destroy this operator before giving the new one.
   The communicator read is then invalid, and the second solution crashes
   (`MPI_ERR_COMM` or segmentation fault). Reproducer:
   `reproducers/hypre-parasails-operator`. `HypreEuclid` reads the
   communicator of the new operator and does not have this defect. Fix in
   mfem-mgis: the preconditioner of the factory replaces the operator by
   the new one before calling this method. Possible upstream fix: one line
   of MFEM (`*A` replaced by `*new_A`).
2. ParaSails only computes on the host. When hypre computes on the GPU,
   its setup reads the matrix stored on the GPU and crashes in
   `MatrixSetRow`, from the first setup on. mfem-mgis now refuses it with
   a message (`gpu-gmres-parasails-refused.txt`).

After the fix, ParaSails (default options of MFEM) remains a bad
preconditioner for the thermal of ex8. GMRES takes 1,944 to 5,626
iterations when it converges, and otherwise stops at 10,000 iterations on
a relative residual of 0.5 to 0.7 (`np2-gmres-parasails-v1.txt`). The time
step is split down to the limit of the number of failures.

## Files

- `gpu-*.txt`: GPU, one process (`gpu-gmres-parasails.txt` before the fix,
  `gpu-gmres-parasails-refused.txt` after).
- `np2-*.txt`: CPU, 2 processes (`np2-gmres-parasails.txt` before the fix,
  `np2-gmres-parasails-fixed.txt` after, `np2-gmres-parasails-v1.txt`
  after with `-v 1`, stopped at 40 s).

Command (GPU):

```
./rjh_plate -m assemblage_hexa.msh -r 0 -ns 5 -no-pp --device cuda \
  -lU .../libGenericParallelU3SI2-cuda.so \
  -lA .../libGenericParallelALFENI-cuda.so \
  -msMc SymmetricPositiveDefinite [-pcTh HypreDiagScale|HypreILU|HypreParaSails] \
  [-lsTh MUMPSSolver]
```
