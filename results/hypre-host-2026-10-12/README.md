# hypre on the host with the device on the GPU (12/10/2026)

Option `--hypre-on-host` of `mfem_mgis::initialize` (or variable
`MFEM_MGIS_HYPRE_ON_HOST`, mfem-mgis 279b5583d4 and 3ad387b002, branch
hypre-host): the behaviours, the kernels of mfem-mgis and the computations
of MFEM stay on the GPU, hypre computes on the host (MFEM gives it its data
on the host). For the small linear systems, whose solutions on the GPU are
dominated by the kernel launches.

Build `build/examples-hypre` (mfem-mgis-examples, ex7 and ex8), linked to
`install/mfem-mgis-hypre`. GPU, one process launched by `mpirun -n 1` on
cores 4 to 7, everything by kernels, cuda behaviours, two runs per
variant. Same Newton iterations and same results in both cases.

## Complete ex7 (40 steps, 2,634 unknowns, GMRES and BoomerAMG)

| hypre | Total (s) | Solution (s) | Linear solver (s) | Jacobian (s) |
|---|---|---|---|---|
| on the GPU | 6.73 and 6.65 | 5.61 and 5.76 | 4.94 and 5.08 | 0.44 and 0.45 |
| on the host | 4.32 and 4.57 | 3.46 and 3.47 | 3.12 and 3.13 | 0.11 and 0.11 |

On the CPU (`results/examples-gpu-2026-10-12`): 8.0 s on one core and
6.3 s on two, by kernels. With hypre on the host, the GPU takes the lead.
The Jacobian is also faster: the global matrix is assembled by the kernels
of mfem-mgis instead of the triple product of hypre on the GPU.

## ex8 (5 steps, thermal with 9,516 unknowns by GMRES and BoomerAMG, mechanics by symmetric MUMPS)

| hypre | Coupling (s) | Thermal (s) | Thermal solutions (s) | Mechanics (s) |
|---|---|---|---|---|
| on the GPU | 12.70 and 12.46 | 4.01 and 3.94 | 3.18 and 3.13 | 8.65 and 8.49 |
| on the host | 11.76 and 11.40 | 2.89 and 2.71 | 1.97 and 1.85 | 8.84 and 8.65 |

The thermal gains 1.2 s: as much as MUMPS for the thermal
(`results/ex8-thermal-gpu-2026-10-11`, 11.5 s), without changing the
solver. The mechanics (MUMPS on the host) does not change.

## Tests

`ctest.txt` (mfem-mgis 3ad387b002, after review): `DeviceOptionTest`
(where hypre computes, option, variable, override, refusal of
`--hypre-device-pool-size`, device configured by `MFEM_DEVICE` or by the
program), `HypreILUTestP-cuda-hypre-on-host`, `BatchIntegrationTestP-cuda`
and `BatchIntegrationTestP-cuda-hypre-on-host` (63 s and 32 s), 18/18. Both
settings of `initialize` are needed: without the flag of MFEM (set before
the initialisation of hypre), the test where the program creates the
device fails; without the host policy (set after the device), four tests
fail, since hypre compiled for the GPU computes on the GPU by default. The
measurements were done with 279b5583d4, whose behaviour is the same for a
device created by `initialize`.
