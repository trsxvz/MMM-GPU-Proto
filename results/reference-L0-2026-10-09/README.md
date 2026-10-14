# Reference times of the work package L0 (09/10/2026)

Current path of mfem-mgis (branch `gpu-proto`, 373550f0d0), MFEM 4.10 of
the environment `mmm-gpu`, on CPU (no device), one OpenMP thread per
process. Machine: AMD Ryzen 7 9800X3D (8 cores, 16 threads), 30 GB, quiet
machine. Script: `cases/measure-reference.sh`; summary per phase:
`cases/summary-profiling.py *.txt`.

## Cases

| Case | Unknowns | Processes | Steps | Newton (total) |
|---|---|---|---|---|
| beam r3, PCG + BoomerAMG | 15,795 | 4 | 10 | 22 |
| beam r3, PCG + BoomerAMG, tolerance 1e-4 | 15,795 | 4 | 10 | 25 |
| beam r3, MUMPS | 15,795 | 4 | 10 | 22 |
| beam r4, PCG + BoomerAMG, tolerance 1e-4 | 111,843 | 8 | 10 | 31 |
| beam r4, MUMPS | 111,843 | 8 | 10 | 26 |
| ssna303 3D (ex4), FGMRES + BoomerAMG | 27,945 | 4 | 10 | 26 |
| ssna303 3D (ex4) r1 | about 8 times more | 8 | 10 | 26 |
| polycrystal, 5 grains (mm-opera-hpc) | 20,553 | 4 | 10 | 98 |
| polycrystal r1 | 168,885 | 8 | 10 | 147 |
| thermomechanical (ex8), MUMPS | 9,516 | 4 | 2 | 17 |

The beam gives the same solution with all the solvers (final norm
18.38573629 at r3, 48.83704586 at r4).

## Times per phase (seconds, maximum over the processes)

In parentheses: number of calls. Integration: integration of the
behaviours (`NS::processNewUnknownsEstimate`); preconditioner: construction
of the linear solver (factorisation of MUMPS; BoomerAMG is built at the
first product, counted in the solution); prediction: elastic prediction at
the beginning of the step (assembly and solution).

| Case | Total | Integration | Residual | Tangent | Preconditioner | Solution | Prediction | End of step |
|---|---|---|---|---|---|---|---|---|
| beam r3 PCG+AMG | 2.55 | 0.09 | 0.08 | 0.32 | 0.00 | 1.32 | 0.73 | 0.00 |
| beam r3 PCG+AMG 1e-4 | 2.20 | 0.10 | 0.09 | 0.36 | 0.00 | 1.13 | 0.51 | 0.00 |
| beam r3 MUMPS | 3.42 | 0.09 | 0.08 | 0.32 | 1.32 | 0.06 | 1.54 | 0.00 |
| beam r4 PCG+AMG 1e-4 | 15.45 | 0.50 | 0.43 | 2.19 | 0.01 | 9.43 | 2.84 | 0.03 |
| beam r4 MUMPS | 65.62 | 0.41 | 0.36 | 1.84 | 35.81 | 0.79 | 26.34 | 0.03 |
| ssna303 3D | 10.35 | 3.18 | 0.31 | 1.91 | 0.00 | 2.44 | 2.37 | 0.05 |
| ssna303 3D r1 | 66.86 | 12.61 | 1.26 | 8.37 | 0.01 | 28.95 | 14.87 | 0.48 |
| polycrystal | 49.73 | 27.48 | 1.13 | 3.22 | 0.00 | 17.65 | - | 0.01 |
| polycrystal r1 | 386.71 | 139.16 | 6.57 | 26.13 | 0.13 | 214.33 | - | 0.11 |
| thermomechanical | 5.07 | 1.62 | 0.37 | 0.85 | 1.70 | 0.46 | - | 0.01 |

## Findings

- **What the GPU must target**: the linear solution (up to 61% of the time
  for the beam, 55% for the refined polycrystal) and the integration of
  the behaviours (55% of the polycrystal, 31% of ssna303). The tangent
  weighs 7 to 18%. The residual is small everywhere.
- **Prediction**: with MUMPS, it refactorises the elastic operator at each
  step (26 s out of 66 for the beam r4). This is the gain expected from
  the prediction with `K0` factorised once (plan, fixed point).
- **Comparison with the old stack** (MFEM 4.9, mfem-mgis master
  4c8a91d815, study of 04/10), same case, same Newton of mfem-mgis, beam
  r4 on 8 processes, same solutions and same iterations:
  - PCG + BoomerAMG: 15.4 s against 17.5 s (faster);
  - MUMPS: 65.6 s against 42.2 s (1.55 times slower). OpenMPI compiled
    with CUDA explains about 10% of it: 59.3 s with
    `mpirun --mca accelerator null`. The rest is not explained: MUMPS
    5.9.1 against 5.8.1, ScaLAPACK 2.2.3 against 2.2.2, OpenBLAS with
    OpenMP against without, MFEM 4.10 against 4.9. Both MUMPS are built
    without METIS or SCOTCH (AMD or PORD reordering). To try when MFEM is
    rebuilt: MUMPS `+metis +parmetis`.
- The time of OpenMPI compiled with CUDA on the CPU exchanges (about 10%
  here) is to be kept in mind for the CPU against GPU comparisons: the CPU
  measurements can be redone with `--mca accelerator null`.
