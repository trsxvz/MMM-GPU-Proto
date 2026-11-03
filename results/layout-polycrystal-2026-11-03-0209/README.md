# Layout of the data at the integration points on the polycrystal (03/11/2026)

Polycrystal of mm-opera-hpc (5 grains, single crystal behaviour of UO2,
finite strain, tangent operator blocks of 81 components, 20,553 unknowns),
10 steps. mfem-mgis b283483f4b, script `cases/measure-layout.sh
polycrystal`, executable built against `install/mfem-mgis-layout`. One RTX
5080, 16 cores, quiet machine (load 1.1). Stiffness matrices assembled in
a fixed order.

| Computation | AoS | SoA |
|---|---|---|
| GPU (integration of the behaviour) | 11.33 s (5.52 s) | 11.30 s (5.50 s) |
| CPU, 8 processes (integration of the behaviour) | 29.08 s (14.10 s) | 28.15 s (13.89 s) |

Total times (`cases/summary-profiling.py`). The macroscopic stresses are
the same to the last bit in both layouts (`summary.txt`). The integration
of the single crystal behaviour dominates and is bound by its
computations: the layout makes no difference on this case.
