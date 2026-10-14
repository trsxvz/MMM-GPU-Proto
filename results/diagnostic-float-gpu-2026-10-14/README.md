# Divergence of float on the GPU, beam (14/10/2026)

Follow-up of results/precision-beam-gpu-2026-10-14-0604: on the GPU, in
float, with the Newton criterion 1e-4, one beam computation out of three
diverged, and the converged computations differed from double by 3e-4. Beam
r4 (111,843 unknowns), J2 plasticity, PCG of hypre and BoomerAMG (tolerance
1e-4), one process on the GPU, NVIDIA RTX 5080 GPU.

## Conclusion

Cause: by default, nvcc contracts multiplications and additions into fused
multiply-adds (FMA, `-fmad=true`). In the kernels of the stiffness
matrices, in float, these contractions made Newton diverge, or converge
towards less accurate solutions. Most likely: Newton drifted along the soft
modes of the beam (bending), which the criterion on the residual does not
see.

Float works close to its limit here. The tangents by kernels and element by
element differ by only 2e-7 of the largest term, but at step 8, the Newton
direction is only known to within 16 to 32% (gap between PCG and MUMPS). The
rounding of the FMA was enough to tip the computations over.

Fix in mfem-mgis (efb06c1267): in float, the GPU code of the kernels that
compute stiffness matrices (standard, FBar, Faltus et al.) is compiled
without contraction (`-fmad=false` for nvcc, `-ffp-contract=off` for
clang), by the CMake function
`mfem_mgis_disable_device_floating_point_contraction`. In double, the
contractions are kept: they did not change the convergence, and the
assembly of the beam was 13 to 16% slower without them. With nvcc, the code
compiled for the host is unchanged.

After the fix, 31 computations out of 31 converge (28 with atomic assembly,
n = 8 to 14 steps, and 3 identical computations with ordered assembly),
within 4.2e-5 at most of double. In float, the assembly costs as much as
before (83 to 85 ms for 16 assemblies, before as after). Performance
measurement redone: results/precision-beam-gpu-2026-10-14-0843 (float 1.37
times faster than double, 3 computations out of 3 converged, within 3.4e-5
at most of double). Reproducer test: `beam-gpu-r4` in cases/beam (ordered
assembly, n = 10, norm compared to double within 5e-5).

## Steps

1. Localisation (`scratch`, ordered assembly, n = 10). With PCG and
   BoomerAMG, only the element by element tangent converges. It is then
   computed and assembled on the host by the loop of MFEM (425 ms per
   assembly against 5.3 ms by kernels), while the behaviours, the gradients
   and the internal forces stay on the GPU. Gradients and internal forces
   by kernels are ruled out. All the attempts by kernels with ordered
   assembly at n = 10 diverge (PCG of hypre, conjugate gradient of MFEM,
   linear tolerance 1e-6). A computation with MUMPS (on the CPU)
   converges.
2. Tangents (`tangents.txt`). Assembled by kernels (K) and element by
   element (E) at the same state, at the end of steps 1, 5, 7 and 8 (step 8
   being not converged, after 50 iterations): maximum gap 2e-7 of the
   largest term of the matrix, Frobenius gap 4e-8, that is the rounding of
   float. Same number of terms (4,702,185). Same asymmetry from step 5 on
   (5e-8 of the maximum term; at step 1, 0 for K and 1.3e-8 for E). The
   assembly does not change the residual (bitwise identical,
   `directions-step8.txt`).
3. Newton directions (`directions-step8.txt`), at the state of the first
   iteration of step 8 of the trajectory by kernels (residual 5.70e5), on
   copies of the matrices migrated to the GPU (the original matrices give
   the same ratios within 1%). A single state measured.

   | Matrix | Solver | Iterations | True linear residual | New residual |
   |---|---|---|---|---|
   | K, float | MUMPS | | 5.4e-3 | 1.15e5 (0.20) |
   | K, float | PCG 1e-4 | 65 | 7.2e-3 | 1.64e5 (0.29) |
   | K, float | PCG 1e-6 | 102 | 7.8e-3 | 1.64e5 (0.29) |
   | E, float | MUMPS | | 3.9e-3 | 7.99e4 (0.14) |
   | E, float | PCG 1e-4 | 62 | 4.1e-3 | 9.45e4 (0.17) |
   | K, double (its trajectory) | MUMPS | | 8.9e-14 | 1.69e4 (0.030) |
   | K, double (its trajectory) | PCG 1e-4 | 64 | 1.1e-4 | 1.69e4 (0.030) |

   In float, the true linear residual levels off between 4e-3 and 9e-3,
   whatever the solver and the tolerance. The directions of PCG and of
   MUMPS differ by 32% with K and by 16% with E, against 6e-5 to 4e-4 in
   double (copies of the matrices).
   A Newton step reduces the residual by 3.5 to 7 times, against 34 times
   in double. The criterion of PCG on the 2-norm of the residual changes
   nothing. In the computations that diverge, the relative residual
   reaches its minimum at the second iteration (1.0e-4 to 1.5e-4, just
   above the criterion), then grows by 1.5 to 2 times per iteration up to
   3e-3.
4. Sensitivity (`gaps.txt`), n = 8 to 14 steps, criterion 1e-4, atomic
   assembly (not deterministic, one computation per n), gap of the norm of
   the displacements to double on the GPU at the same number of steps:

   | Configuration | Converged | Maximum gap |
   |---|---|---|
   | float, GPU, tangent by kernels | 4 out of 7 | 2.8e-4 |
   | float, element by element tangent (host), behaviours on the GPU | 7 out of 7 | 4.9e-5 |
   | float, CPU, kernels (two processes) | 7 out of 7 | 7.6e-5 |
   | float, GPU, kernels, criterion 3e-4 | 7 out of 7 | 2.4e-4 |
   | float, GPU, kernels, criterion 1e-3 | 7 out of 7 | 5.3e-4 |

   On the computations of the day kept from before the fix, with criterion
   1e-4, with atomic assembly (n = 8 to 14, measurement 0604, profiling), 5
   computations out of 11 diverge (4 out of 8 in a first series).
5. Without contraction. mfem-mgis entirely compiled with `-fmad=false`
   (`gaps-without-contraction.txt`): 7 computations out of 7 converge by
   kernels, within 3.7e-5 at most of double. Element by element, the
   results are identical to 10 digits with and without contraction: on this
   case, only the stiffness kernels matter. At the first state of step 8
   of this trajectory (`directions-step8-without-contraction.txt`): true
   linear residual 5e-4, MUMPS step that reduces the residual by 20 times,
   K and E equivalent, and rigid translations better respected by K (4.2e-8
   of the maximum term against 5.1e-8 with contraction). Fix targeted on
   the stiffness kernels (`gaps-fix.txt`, first version then final
   version): 31 computations out of 31, including 3 identical ones with
   ordered assembly.

## Method

Diagnostic program (draft, not tracked): copy of cases/beam with two hooks.
The first one acts after the solution of the steps listed in
`COMPARE_TANGENT_STEPS`. It assembles the tangent by kernels
(`GetGradient`), then element by element
(`ParNonlinearForm::GetGradient`, after
`setStiffnessMatricesComputedByKernels(false)` and
`synchronizeOnHostForReading`), and compares the two matrices term by term
on the host.

The second one stops the step `COMPARE_DIRECTIONS_STEP` after one Newton
iteration (`COMPARE_MAX_ITERATIONS=1`), then the program stops. At this
state, it computes the Newton directions with MUMPS and with PCG and
BoomerAMG configured as in mfem-mgis (system options set after the
operator), on copies migrated to the GPU, then on the original matrices.

For each direction, it computes the true linear residual (accumulated in
double on a copy of the float matrix) and the Newton residual obtained
after the integration of the behaviours. The solutions are done before the
copy to the host, which moves the matrix of hypre to the host. The assembly
is ordered in all these computations.

Other attempts: without prediction, the computation fails at step 1 in all
modes (a case that requires the prediction); cuDSS is not available in the
float environment.
