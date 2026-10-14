# Defects found and corrected

Defects of the existing code were found on the way, mostly by the
comparison tests, the GPU compilers (which refuse what the host
compiler accepts) and independent reviews. In TFEL and MGIS, each one
is corrected in its own commit and covered by an existing or a new test.
In mfem-mgis, some commits still mix several subjects and some
corrections have no test yet: they will be separated in the upstream
contributions. Some corrections change CPU results: that was their
purpose. The sources, in French, are
`docs/tfel.md`, `docs/mgis.md`, `docs/integration-amont.md`,
`docs/plan-precision.md` and `JOURNAL.md`.

Compiler defects and defects of other libraries that were worked around
rather than corrected are listed in [Pitfalls](pitfalls.md).

## TFEL and MFront

### `generic` interface of behaviours

- Status not checked after a finite strain integration (`if (r)` in the
  finite strain drivers).
- `DTAU_DDF` operator with the Hencky (logarithmic) strain: the end of
  step tangent used the beginning of step handler.
- Post-processings saw the initial internal state variables with the
  final stresses: the end of step state is now given to them.
- Post-processing without initial state: a null pointer was passed by
  reference (undefined behaviour).
- Initialize functions of the Green-Lagrange and Hencky strategies had
  their arguments inverted (compilation error).
- Speed of sound with the prediction operator: the decoding ignored the
  offset of the speed of sound.
- Pointers not restored after an exception thrown by an initialize
  function.
- Stress computed by an initialize function was only exported for the
  Cauchy stress of standard finite strain behaviours; it is now
  converted for all stress measures.
- Exceptions thrown by the constructor of the behaviour (for example a
  missing file read by `@SharedData`) escaped the generated function
  instead of being reported as a failure.
- Rotations of orthotropic behaviours: point offset forgotten for second
  order tensors that are not the first variable, scalars not copied,
  vectors and vector/vector blocks not supported (`std::exit(-1)` at run
  time), and now rotation of all the blocks of the tangent operator
  (arrays of variables, pairs of variables without a TFEL type).

### Models and material properties

- Bounds of `generic` models checked in all compilation modes; on single
  pass compilers, physical bounds of models were silently ignored.
- Generated classes of `generic` models now have internal linkage: the
  parameters of a model were shared by all its libraries loaded in a
  process (the singleton was a `STB_GNU_UNIQUE` symbol).
- A model with material properties and parameters treated as static
  variables did not compile with the `generic` interface.
- A static variable with units (`@StaticVariable temperature Tref = 300;`)
  in a model using quantities did not compile.
- Physical bounds of the arguments of material properties were not
  exported when the argument had no standard bounds.
- `generic` interface of material properties: wrong message for the
  physical upper bound ("below" instead of "over"), `errno` not restored
  when an upper bound is violated with the `Strict` policy, and
  documentation of the sign of the rank of the variable fixed to match
  the code.
- `generic-parallel` interface of material properties (existing code).
  The generated code did not compile with the runtime checks disabled.
  The SYCL device memory was never freed, the extent of strided fields
  was too large, and the values between the values of a strided output
  were overwritten (SYCL). Transfers were made outside the `try` block,
  the include guard missed the backend name, and bounds and default
  values were written with 6 digits instead of 14. The grid size of the
  CUDA kernels was wrong, and the kernels read the array of argument
  pointers, often on the caller's stack. In the kernels, `errno` is
  replaced by the detection of non finite results.
- Properties defined by data interpolation (`@Data`): the CUDA and HIP
  kernels used the name of the abscissae for their input fields.
- With clang CUDA, a property calling another law did not compile (only
  the host versions of the called laws were visible).

### DSLs, bricks and TFEL headers

- Runge-Kutta algorithms: the stiffness modifiers kept, in a static
  lambda, the name of the first behaviour processed by a call of MFront.
- DDIF2 brick: inverted test on the number of non constant fracture
  stresses, first fracture energy used in all three directions, and
  checks of the options (only one of `softening_slope`,
  `softening_slopes`, `fracture_energy` and `fracture_energies`).
- Isotropic DSLs with `use_stress_update_algorithm`: Poisson ratio
  converted to `real` instead of `stress`, which did not compile with
  `@UseQt true` (test behaviour `Plasticity8`).
- `rkCastem`: the variable `young` is now used as the error
  normalisation factor, as announced by the TFEL 4.1 release notes.
- `StandardElastoViscoPlasticity` brick: the name of
  `mfront_computeYieldRadiusAndDerivative` was doubled in the spline
  variant; the two argument declarations of the Chu-Needleman nucleation
  functions were never defined.
- `MetallicCFCGenericSlidingSystemsInteractionMatrix`: the documented
  cache of the first call, which ignored the coefficients of later
  calls, is removed (the matrix is computed at each call).
- 85 definitions without `__host__ __device__` whose declaration had
  them (accepted by nvcc, refused by clang CUDA).
- Test behaviour `DigitSeparatorTest`: uninitialised derivative of the
  flow rule.

### Build, configuration and tests

- Generated Makefiles placed the linker options before the object files:
  TBB was dropped with `--as-needed`.
- CMake: the substitution of a compiler was only recorded when CMake
  searched it; `HIP_PLATFORM` not cached; the `hip-hipcc` configuration
  passed `-isystem /usr/include`, which breaks the `#include_next` of
  libstdc++; HIP headers not given to hipcc on the NVIDIA platform;
  wrong `endif` argument in the `hip-clang` configuration.
- The `sycl-gpu-acpp` tests used configuration files only installed with
  `sycl-default-acpp`.
- CUDA test of material properties wrote host memory from a kernel
  (hidden by heterogeneous memory management, HMM).
- HIP test of material properties: `hipMallocManaged` received a number
  of reals instead of a size in bytes, and the test assumed the managed
  memory to be zero.

## MGIS

- Rotation matrix axes refused when their number of values was correct
  (inverted test), and the number of points compared to a remainder
  instead of a quotient: no `RotationMatrix2D` or `RotationMatrix3D`
  could be built.
- A material property could not be given by a `std::span` of one value.
- The failure of a point-wise initialize function or post-processing
  did not set the status to `-1`.
- A uniform material property value in local storage is stored as a
  scalar, as for external variables.
- The mass density can be given by a `std::span` of one value.
- HDF5 restoration of a multi-component external variable in local
  storage wrote outside its array.
- CUDA, HIP and Kokkos memory resources refused alignments stricter than
  their memory: a standard library pool failed above 256 bytes.
- `MaterialDataManagerInitializer` had no constructor in Python.
- Error message of the rotations of the tangent blocks spoke about
  thermodynamic forces.
- `MaterialAxisStorage` not exported (copying a rotation matrix outside
  MGIS did not link), and the copy of a locally stored axis pointed to
  the storage of the original.
- Fortran bindings of the tangent operator returned an array of wrong
  shape when the operator is not square.
- The Python binding read a strided numpy array as if it were contiguous.
- `update` and `revert` of a manager allocated an execution context at
  each call.
- Before this work, a `generic-parallel` library was loaded without error
  and its batch functions called with the arguments of one point
  (undefined behaviour); it is now refused with a clear message.
- Documentation: `update` described backwards, names of the rotations of
  finite strain behaviours, wrongly annotated parameters.

## mfem-mgis

Corrections of the existing code, without GPU code:

- `HypreParaSails`: the preconditioner is given the new operator before
  the call, working around a defect of MFEM (see below); ParaSails
  crashed everywhere.
- Settings of the ILU of hypre reapplied after each new operator: MFEM
  resets them to their defaults (the fill level was lost at the second
  operator, also on CPU).
- Absolute tolerance of the Krylov solvers of hypre: dead code for
  `HypreGMRES`, and `HypreFGMRES` called `SetTol` instead.
- Short names of the options declared by `declareDefaultOptions`
  (`-up`, `-no-up`, `petscrc_file`, and `-dev` for the device) ignored
  by `initialize`.
- An exception thrown by the tangent operator during the Newton
  correction (`noexcept`) ended in `std::terminate`; it is now a
  failure of the resolution in sequential, and a message followed by
  `MPI_Abort` in parallel.
- Exceptions of the L2 projections, declared `noexcept`, reported in the
  execution context.
- Sequential prediction: diagonal set to one at the imposed degrees of
  freedom, as in parallel.
- `SystemHypreBoomerAMG` reset two BoomerAMG settings at each
  `SetOperator`.
- CMake: several directories returned by `tfel-config`, one line of
  `env.sh` per directory; `FindMFEM.cmake` truncated the values of
  `config.mk` containing an equal sign.
- Documentation of the constructors that take an execution context.

Corrections of the GPU path of the existing code, found by running the
whole test suite on the GPU:

- Host accesses declared to MFEM (`HostRead`, `HostReadWrite`): the L2
  projection and the implicit gradient regularisation of a
  multi-component function gave zeros on the GPU; the ParaView export
  overwrote the unknowns computed on the GPU with stale host values;
  Dirichlet imposed values, field interpolation (GSLIB) and the
  analytical L2 error read without synchronisation.
- On the GPU, the `Elasticity` strategy of BoomerAMG stops the program
  in hypre: replaced by `System` with a warning.
- ILU of hypre limited to level 0, and ParaSails refused, when hypre
  computes on the GPU (crashes in hypre).

Single precision (official code that did not work with a float MFEM):

- `MPI_DOUBLE` used on `real` data (silent overflow in
  `ComputeResultantForceOnBoundary.cxx`, no compilation in
  `MeanThermodynamicForces.ixx`), and four such reductions in the tests
  ("stack smashing" in parallel).
- `double Eval(...) override` in `PartialQuadratureFunction.cxx`, and
  `numeric_limits<double>::max()` stored in a `real`.
- `FindMFEM.cmake` required SuiteSparse, which MFEM refuses in single
  precision.

## mfem-mgis-examples

- ex8 on the GPU read its results in stale host copies (comparison with
  the reference failed while the solutions were right).
- ex8: the thermal problem reads the mechanical gradients on the host;
  they must be synchronised when the mechanics is integrated on the GPU.
- ex4 now prints the cause of resolution failures.

## mm-opera-hpc

- The loop over the bubbles, on host data, used `mfem::forall`, which
  since MFEM 4.10 requires a CUDA compilation when MFEM has CUDA (and
  would launch the loop on the GPU): replaced by `mfem::hypre_forall_cpu`.
- The default bubble case (467,109 unknowns) failed even point by point
  (PCG with the diagonal stopped at 1000 iterations); a preconditioner
  option was added, and it converges with BoomerAMG.
- From the review: help `-h` without failure, unit of the bubble
  pressure, output files opened after the options and written by the
  main process only.
- Single precision: MPI types of `real`, and the Kröner loop bounded
  (tolerance `max(1e-10, 10 ε)`, at most 1000 iterations).

## Upstream libraries: pull requests

Corrections proposed to MFEM and AmgX during the work (open on
13/10/2026):

| Pull request | Defect |
|---|---|
| MFEM #5604 | `AmgXSolver`: its resources, shared by all its solvers, kept the communicator and configuration of the first one; destroying the first solver made the others fail. The fix keeps a shared resource with its own configuration and a duplicated communicator, and calls `InitSerial` once |
| MFEM #5605 | Order of the vertices of the boundary elements of a `SubMesh`, now taken from the parent boundary element |
| MFEM #5615 | `HypreParaSails::SetOperator` created the new preconditioner with the communicator of the previous operator, already freed by the nonlinear solvers (`MPI_ERR_COMM`, reproducer `reproducers/hypre-parasails-operator`) |
| AmgX #386 | Classical multigrid with CUDA 13 failed with "invalid device function": a weak copy of `cub::DeviceFor::Bulk` in `amg.cu` launched a null kernel (closes AmgX issue #380; reproducer `reproducers/amgx-classical-cuda13`) |

The rebuild of MFEM and AmgX in the spack environment waits for their
merge upstream, then in spack-packages.

## Single precision on the GPU (14/10/2026)

Source: `results/diagnostic-float-gpu-2026-10-14/README.md`.

- **Symptom.** On the GPU, in float, with a Newton criterion of 1e-4,
  the beam r4 (111,843 unknowns, J2 plasticity, PCG and BoomerAMG)
  diverged in 5 computations out of 11 with the atomic assembly (8 to 14
  steps), and the converged ones differed from double by up to 3e-4.
  With the ordered assembly and ten steps, every computation by kernels
  with an iterative solver diverged; there, only the tangent assembled
  element by element on the host converged.
- **Cause.** nvcc contracts by default multiplications and additions
  into fused multiply-add operations (FMA, `-fmad=true`). In the
  stiffness kernels, in float, these contractions made Newton diverge or
  converge to less accurate solutions, most probably by a drift along
  the soft modes of the beam (bending), which the residual criterion
  does not see.
- **Why so sensitive.** The kernel and element by element tangents only
  differ by 2e-7 of the largest term, but at step 8 the Newton direction
  is only known to within 16 to 32 % (gap between PCG and MUMPS): float
  works here near its limit.
- **Correction in mfem-mgis.** In float, the GPU code of the kernels that
  compute stiffness matrices (standard, FBar, Faltus et al.) is compiled
  without contraction (`-fmad=false` for nvcc, `-ffp-contract=off` for
  clang), through the CMake function
  `mfem_mgis_disable_device_floating_point_contraction`. In double the
  contractions are kept: they did not change convergence, and the
  assembly of the beam was 13 to 16 % slower without them.
- **Result.** 31 computations out of 31 converge, within 4.2e-5 of
  double, at the same assembly cost. Measurement redone
  (`results/precision-poutre-gpu-2026-10-14-0843`): float 1.37 times
  faster than double, same iterations, within 3.4e-5 of double. A
  reproducer test, `poutre-gpu-r4`, is in `cas/poutre`.
