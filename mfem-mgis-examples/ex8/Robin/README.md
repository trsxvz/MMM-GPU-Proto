# Robin Thermal Convection Test

**Author:** Pierre Dubois

## Overview

This test case validates the implementation of thermal convection through a **Robin boundary condition** in MFEM-MGIS.

The problem consists of a steady-state heat conduction test in a three-dimensional bar with constant thermal conductivity. A prescribed temperature is imposed on one end of the domain, while a Robin (convective) boundary condition is applied on the opposite end:

$$
\mathbf{j} \cdot \mathbf{n} = h \left(T - T_\infty\right),
$$

where $h$ is the heat transfer coefficient and $T_\infty$ is the ambient temperature.

Because the analytical solution is known, this example provides a straightforward verification of the implementation by comparing the numerical and exact solutions. It also validates the residual and tangent contributions associated with the Robin boundary condition.

The test `robin_test`, run by `ctest` on 2 processes, fails if the L2 norm of the difference between the numerical and exact temperatures exceeds 1e-8.

The library of the behaviour is built for the CPU and, if MFEM uses a GPU, for the GPU; the one
placed by the option `--execution` of mfem-mgis is used, unless a library is given by the option
`-l`. With the `MFEM_MGIS_EXAMPLES_GPU_TESTS` option (see the main README), `robin_test-gpu` runs the
same test on the GPU (`--execution gpu`).
