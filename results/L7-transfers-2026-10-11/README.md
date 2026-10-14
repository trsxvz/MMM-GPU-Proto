# Transfers between the host and the GPU (work package L7.3, 11/10/2026)

nsys profile (`cases/transfers-profile.sh`) of the steady-state thermal
problem r5 (279,873 unknowns, 5 steps, 18 Newton iterations) and of the
beam r4 (111,843 unknowns, 10 steps, NewtonK0 and Eisenstat-Walker 1,
aggressive BoomerAMG), all the kernels, `cuda` behaviour, mfem-mgis
e03fd927f7. Copies grouped by direction and by size. Their origin is given
by gdb (breakpoints on `mfem::CuMemcpyDtoH` and `mfem::CuMemcpyHtoD`
filtered by size, `gdb-copies.txt`, thermal problem over one step).

| Case | Duration | Host to GPU | GPU to host |
|---|---|---|---|
| Thermal r5 | 1.46 s | 60 ms | 58 ms |
| Beam r4 | 1.97 s | 33 ms | 30 ms |

Thermal: the copies represent 8 % of the time.

| Copy | Count | Time | Origin |
|---|---|---|---|
| states, 3 components per point, GPU to host | 40 | 40 ms | end of step (`update`), state arrays brought back to the host |
| states, 3 components per point, host to GPU | 30 | 30 ms | next integration, arrays sent back to the GPU |
| temperature at the points, host to GPU | 92 | 21 ms | integration: external state variables at the beginning and at the end of the step |
| temperature at the points, GPU to host | 46 | 14 ms | values of the unknowns at the points copied into the external state variable |

Beam: stresses and strains (6 components) and internal state variables (7
components), 80 and 40 copies in each direction, 56 ms in total (3 %),
same origin (end of step, then next integration).

## Decisions

- End of step: the update is done on the GPU, then the arrays are brought
  back to the host and their copy on the GPU is invalidated. This is the
  documented contract of `MaterialArrays`: after an update, the state may
  be modified directly on the host (case "write after update" of
  `BatchIntegrationTest`). A read-only synchronisation, which would keep
  the copy on the GPU valid, would save the copy back (2 % on the thermal
  problem, 1.4 % on the beam) but breaks this contract: trial cancelled
  (the test fails on the `debug` and `cuda` devices). An option that would
  keep the states on the GPU at the price of explicitly synchronised
  writes on the host remains possible.
- Temperature: keeping it on the GPU (source on the GPU for the external
  state variable at the end of the step, copy at the beginning of the step
  kept from one iteration to the other, copy to the host only at the end
  of the step) would save at most 2.3 % of the thermal problem, for one
  more validity management in `BehaviourIntegratorBase`. Not done for now.
