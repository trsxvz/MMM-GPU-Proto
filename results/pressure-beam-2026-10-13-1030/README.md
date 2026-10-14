# Imposed pressure on the beam, work package L7.6 (13/10/2026)

Beam r4 of the prototype (111,843 unknowns), pressure of 1e8 Pa on the
end, 10 steps, NewtonK0 (Eisenstat and Walker tolerance, one cycle of
BoomerAMG on K0), all the computations by kernels, on the GPU (one process,
hypre on the GPU) and on 8 CPU processes. Script:
`cases/measure-pressure-beam.sh`.

- Before: executable of 11/10 (`scratch/night-2026-10-11/l76/beam-before`,
  version of mfem-mgis earlier than the work package L7.6): the pressure
  is a boundary integrator of the nonlinear form, evaluated on the host at
  each residual; its presence also forces the assembly of the tangent by
  MFEM on the host, from the matrices of the elements computed by kernels.
- After: mfem-mgis 54ff8e6c8a (separate installation
  `install/mfem-mgis-54ff8e6c8a`, build `build/cases-54ff8e6c8a`): the
  pressure, independent of the unknowns, is assembled once for a unit
  pressure and added to the residual on the GPU, without any integrator in
  the form; the tangent is still assembled by kernels, or applied without
  assembly.

Time of the loop over the steps (seconds):

| | Without pressure, before | With pressure, before | Without pressure, after | With pressure, after |
|---|---|---|---|---|
| GPU | 1.81 | 17.63 | 1.50 | 1.70 |
| 8 processes | 7.13 | 9.88 | 5.70 | 5.68 |

With the tangent applied without assembly and the pressure: 1.86 s on the
GPU.

- On the GPU, the pressure cost 15.8 s (tangent 14.3 s, assembled by MFEM
  on the host; residual 1.8 s); it now costs only 0.2 s: the computation
  with pressure is 10 times faster.
- On 8 processes, the extra cost of the pressure (2.75 s) disappears.
- Without pressure, the new version is also faster (1.81 → 1.50 s on the
  GPU, 7.13 → 5.70 s on 8 processes): work packages of 11/10 later than the
  executable of before (positions of the entries, matrices without forming
  B, inverses of the Jacobians).
- A first run, done during compilations of another session, gave wrong CPU
  times (23 s instead of 7 s); it was redone on a quiet machine.
