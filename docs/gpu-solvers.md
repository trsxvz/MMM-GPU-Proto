# Solvers and techniques for the CPU and GPU study of mfem-mgis

Status: list drawn up on 09/10/2026, updated on 13/10/2026. The rows and
the items marked "done" are implemented in mfem-mgis; the others remain
leads. This list is used to
choose the performance studies of mfem-mgis and the packages of the
spack environment `mmm-gpu` ([deps/spack](../deps/spack/README.md)).

Scope: only what runs on CPU and on GPU, or on GPU only.
What runs only on CPU is excluded (see [Excluded](#excluded-cpu-only)).

Each item states where the feature exists, what was checked and how,
what mfem-mgis would need, and an interest out of 10 for this project.
"Checked" means read in the sources or the release notes of the quoted
version (see [Sources](#sources)), not yet tried on the machine. "To be
checked" means not confirmed.

## Contents

1. [Overview table](#overview-table)
2. [Reusing the preconditioners](#reusing-the-preconditioners)
3. [Nonlinear solvers](#nonlinear-solvers)
4. [Matrix free preconditioners](#matrix-free-preconditioners)
5. [AMG on GPU](#amg-on-gpu)
6. [Mixed precision](#mixed-precision)
7. [Direct solvers on GPU](#direct-solvers-on-gpu)
8. [Several GPUs](#several-gpus)
9. [Assembly and tangent operator](#assembly-and-tangent-operator)
10. [Coupled problems](#coupled-problems)
11. [Measurement tools](#measurement-tools)
12. [Excluded (CPU only)](#excluded-cpu-only)
13. [MFEM 4.10](#mfem-410)
14. [Mapping to the spack environment](#mapping-to-the-spack-environment)
15. [Sources](#sources)

## Overview table

The rows without a mark come from the initial list. The rows marked
"addition" complete this list. Ginkgo is installed in the environment
(MFEM is compiled with it), not yet tried: it appears in several
rows.

| Technique | CPU/GPU | Interest | Where |
|---|---|---|---|
| Reuse the preconditioner (AMG built once, on the elastic operator, or kept for several iterations or time steps) | both | 9 | done: method `NewtonK0` of the `FixedPoint` solver of mfem-mgis (AMG built once per solve on the elastic operator, or once for the whole computation with `ElasticOperatorUpdate` set to `Once`) |
| Addition: partial reuse of the AMG hierarchy (prolongators kept, coarse operators recomputed) | both | 7 | GAMG of PETSc; AmgX (`AMGX_solver_resetup`, not exposed by MFEM) |
| Fixed point accelerated by Anderson | both | 9 | done: `FixedPoint` solver of mfem-mgis, accelerations of MTest (parameter `Acceleration`: `Anderson`, `UAnderson`, `FAnderson`, `Secant`, `IronsTuck`, etc.); KINSOL not tried |
| Addition: assembly of the global tangent matrix on GPU (option c on GPU), sparsity pattern and hypre matrix reused from one iteration to the next | both | 9 | done in mfem-mgis: tangent assembled by kernels (see [below](#assembly-and-tangent-operator)) |
| Addition: choice of the integration rule (8 points instead of 27 on the first order hexahedron) | both | 8 | done: parameter `IntegrationRuleOrder` of the discretization of mfem-mgis |
| Inexact Newton (Eisenstat-Walker, `NewtonSolver::SetAdaptiveLinRtol` of MFEM) | both | 8 | done in the Newton of mfem-mgis and the method `NewtonK0` (option `AdaptiveLinearSolverTolerance`) |
| Matrix free Jacobi and Chebyshev smoothers on the diagonal (`OperatorJacobiSmoother`, `OperatorChebyshevSmoother`), required by a fully matrix free option a | both | 8 | done: preconditioners `OperatorJacobiSmoother` and `OperatorChebyshevSmoother` of the factory of mfem-mgis (MFEM) |
| Addition: matrix free Jacobi by 3×3 nodal blocks (block diagonal of the tangent) | both | 7 | mfem-mgis, like the `hooke` mini-application of MFEM; done on 15/10/2026 (`NodalBlockJacobi`, patches `0210` to `0213`) |
| Addition: matrix free "absolute value" l1-Jacobi | both | 5 | MFEM 4.9, `miniapps/diag-smoothers` |
| GPU settings of BoomerAMG exposed in mfem-mgis (PMIS, ext+i, l1-Jacobi, Chebyshev, number of levels, thresholds) | both | 8 | done: options of `HypreBoomerAMG` in mfem-mgis (see [below](#boomeramg-hypre)) |
| Addition: BoomerAMG for systems on GPU ("System" strategy, settings preserved at each `SetOperator`; "Elasticity" stops the program with hypre built for GPU) | both | 8 | done in mfem-mgis: "Elasticity" replaced by "System" with a warning, settings applied again after each `SetOperator` |
| Mixed precision: preconditioner in single precision, Krylov in double; large gain on the RTX 5080 (slow FP64), smaller on A100 or H100 | mostly GPU | 6 (8 on consumer GPU) | Ginkgo, hypre 3, AmgX |
| L-BFGS preconditioned by an AMG or p-MG cycle (`LBFGSSolver` of MFEM), like Ratel | both | 7 | to be written: `LBFGSSolver` of MFEM ignores the preconditioner; `SNESQN` of PETSc has one |
| Matrix free p-multigrid with AMG on the Q1 (order 2 and higher), or LOR (`LORSolver`) | both | 7 | MFEM |
| Addition: distributed Ginkgo (MPI matrices, two level Schwarz, distributed PGM multigrid) | both | 7 | Ginkgo 1.8 to 1.11; MFEM 4.10 wrappers |
| Addition: AmgX in mixed precision (mode `dDFI`) and 3×3 block matrices | NVIDIA GPU | 6 | AmgX; not in the MFEM wrapper |
| Line search for Newton (hook `ComputeScalingFactor`) | both | 6 | to be written in the Newton of mfem-mgis (today it only divides the correction by 2 when the behaviour fails) |
| Elasticity AMG with rigid body modes on GPU: GAMG of PETSc with near null space (Kokkos), new block path on GPU (Adams, arXiv 2606.24748); requires repairing the PETSc path of mfem-mgis | both | 6 | PETSc 3.26 (`aijkokkos`; the block path of the paper is not there). Since 15/10/2026, the solvers of PETSc compute on the host when MFEM computes on a GPU: with MFEM 4.10, the vectors of MFEM placed in the vectors of PETSc of the GPU were seen as null, and these vectors made `PetscFinalize` abort; GAMG on the GPU needs this interoperability first |
| KINSOL of SUNDIALS through MFEM (`KINSolver`: fixed point with Anderson, JFNK, line search) | both | 6 | SUNDIALS 7.9, MFEM; requires an operator whose `Mult` integrates the behaviour |
| Addition: prediction by extrapolation of the displacement increment | both | 5 | mfem-mgis |
| AmgX through `AmgXSolver` of MFEM (comparison point with BoomerAMG, MPI modes) | NVIDIA GPU | 5 | done by a wrapper specific to mfem-mgis (`AmgX`, `AmgXPCG`, see [below](#amgx-nvidia)) |
| cuDSS (NVIDIA direct solver, one or several GPUs, MPI, refactorization) | NVIDIA GPU | 5 | cuDSS; wrapper `CuDSSSolver` of MFEM 4.10; done: solver `CuDSSSolver` of mfem-mgis |
| Addition: pipelined Krylov (pipelined CG of Ginkgo, `KSPPIPECG` of PETSc) | both | 5 | Ginkgo 1.10, PETSc |
| Addition: Ratel as reference computation on the same case | both | 5 | Ratel, PETSc, libCEED |
| SuperLU_DIST and STRUMPACK on GPU (existing MFEM wrappers) | both | 4 | MFEM |
| MGR of hypre (GPU since 2.29) for monolithic coupled problems | both | 4 | hypre (setup on GPU since 2.24, without unified memory since 2.29) |
| FSAI (GPU since 2.30) and ILU of hypre on GPU | both | 4 | hypre |
| Addition: recent hypre (filtering of the couplings between components, flexible cycle, overlapping Schwarz, Chebyshev solver) | both | 4 | hypre 2.32 to 3.2 |
| Addition: adaptive time step driven by the convergence of Newton | both | 4 | mfem-mgis |
| Addition: nonlinear solvers of PETSc on GPU vectors (NGMRES, Anderson, quasi-Newton, NASM) | both | 4 | PETSc SNES: the entry `SNES` of the factory of nonlinear solvers is done (15/10/2026), on vectors of the host (see the row of GAMG) |
| Addition: MueLu of Trilinos (smoothed aggregation with rigid body modes, Kokkos) | both | 3 | Trilinos; no MFEM interface |
| Addition: recycling of Krylov subspaces between Newton iterations | both | 3 | PETSc (HPDDM, GCRO-DR); GPU to be checked |
| Addition: behaviours in single precision on consumer GPU | GPU | 3 | requires TFEL and MGIS |
| Addition: batched dense algebra of MFEM (`BatchedLinAlg`, 3×3 blocks and element matrices) and CUDA graphs (small kernels) | mostly GPU | 4 | MFEM 4.8 and later, CUDA |
| Addition (tools): NVTX profiles of hypre, Caliper, Nsight Systems | both | 7 | hypre `+gpu-profiling`, MFEM |
| Addition (tools): Umpire memory pools for MFEM and hypre | mostly GPU | 6 | Umpire |
| Addition (tools): GPU-aware MPI | GPU | 5 | MPI of the cluster, MFEM, hypre |

## Reusing the preconditioners

On GPU, building an AMG often costs more than several
solves. The Newton of mfem-mgis rebuilds it at each iteration: the
linear solver receives the new tangent
(`this->prec->SetOperator(J)` in `NewtonSolver.cxx`), and the `SetOperator` of a Krylov solver
always calls the one of its preconditioner (Krylov solvers of MFEM,
`linalg/solvers.cpp:184-193`; `HyprePCG::SetOperator`; the setup
phase of the Krylov solvers of hypre also rebuilds the preconditioner).

- **Full reuse**: build the AMG once and keep it.
  Variants: on the elastic operator (kept for the whole simulation), on
  the tangent at the start of the step (kept for one step), on the tangent kept
  for `n` iterations. The Krylov solver stays exact, only the preconditioner
  ages. This needs a preconditioner wrapper whose `SetOperator`
  rebuilds nothing (and, with the Krylov solvers of hypre, a preconditioner
  whose setup phase does nothing), plus the update policy
  in the Newton of mfem-mgis. Mechanisms exist elsewhere:
  `KSPSetReusePreconditioner` and `-snes_lag_preconditioner` of PETSc,
  `SetMaxSetupCalls` of KINSOL. Interest 9. Done: the method `NewtonK0`
  of the `FixedPoint` solver of mfem-mgis keeps an AMG built on
  the elastic operator.
- **Partial reuse** (addition): keep the structure of the
  hierarchy (aggregates, prolongators) and recompute only the
  coarse operators `Pᵀ A P` with the new values. Much lower setup
  cost, quality close to a full setup
  as long as the tangent does not change too much. Exists in:
  - GAMG of PETSc (`PCGAMGSetReuseInterpolation`);
  - AmgX (`AMGX_solver_resetup`, option `structure_reuse_levels`), not
    exposed by MFEM: `AmgXSolver::UpdateOperator` only replaces the
    coefficients of the matrix, without a new setup, which is
    a full reuse (checked in `linalg/amgxsolver.cpp` of
    MFEM 4.9);
  - MGR of hypre 3.2 ("component reuse", release notes);
  - BoomerAMG: no equivalent found.
  Interest 7.

## Nonlinear solvers

### Prerequisite: an operator that integrates the behaviour

In mfem-mgis, `Mult` (the residual) does not integrate the behaviour: the integration
is triggered by the nonlinear solver of mfem-mgis, through an action
called at each new estimate (`NonLinearSolverBase.cxx:39-58`).
The `Mult` of the form only assembles the forces already computed.
An external nonlinear solver (`LBFGSSolver`, `KINSolver`,
`mfem::NewtonSolver`, SNES of PETSc) only calls `Mult`: it
would compute the residual with outdated stresses. The only exception is the
PETSc mode, where `AssembleElementVector` integrates the behaviour itself; a failure
of the behaviour there stops the program (`MFEM_VERIFY` of
`MultiMaterialNonLinearIntegrator::AssembleElementVector`), which may
explain why the PETSc path crashed.

So all the external solvers need an operator whose `Mult`
integrates the behaviour then assembles, and which reports a failure of the behaviour without
stopping the program. This is the role of the dedicated operator of the
detailed plan (entry point B).

The Newton of mfem-mgis (`mfem_mgis::NewtonSolver`) is not the one of
MFEM: it derives from `mfem::IterativeSolver` through `NonLinearSolverBase`,
and rewrites `Mult`. The hooks of `mfem::NewtonSolver`
(`SetAdaptiveLinRtol`, `ComputeScalingFactor`) are of no use there.

### Techniques

All these solvers only use vectors, operator actions
and scalar products: they run on GPU as soon as the vectors and
the operators are there.

- **Fixed point accelerated by Anderson** (interest 9). Done in mfem-mgis:
  `FixedPoint` solver (option b of the plan), accelerations of MTest
  (parameter `Acceleration`). Another way, not tried: KINSOL
  (`KINSolver::EnableAndersonAcc`, for `KIN_FP` and `KIN_PICARD`, checked
  in `linalg/sundials.hpp` of MFEM 4.9). With `KIN_FP`, KINSOL expects the
  fixed point function `G(u) = u - K0⁻¹ F(u)`, not the residual: this needs
  an operator that applies it, or `KIN_PICARD`. The plain fixed point
  converged poorly on the polycrystal of the preliminary trials: Anderson is
  the first thing to try.
- **Inexact Newton** (interest 8): tolerance of the linear solver adapted
  to the nonlinear residual (Eisenstat-Walker). Done in the Newton of
  mfem-mgis and the method `NewtonK0`, following
  `mfem::NewtonSolver::SetAdaptiveLinRtol` (parameter
  `AdaptiveLinearSolverTolerance`).
- **L-BFGS** (interest 7): with an AMG or p-multigrid cycle as
  initial inverse, like Ratel. `LBFGSSolver` of MFEM is not suitable:
  it ignores the preconditioner (`SetPreconditioner` prints "L-BFGS
  won't use the given preconditioner", `linalg/solvers.hpp:917-920`, in
  4.9 as in 4.10). To be written; `SNESQN` of PETSc
  (`-snes_qn_type lbfgs`, `-snes_qn_scale_type jacobian`) is the
  existing reference of an L-BFGS with a preconditioned initial inverse.
- **Line search** (interest 6): the Newton of mfem-mgis today
  divides the correction by 2 only when the behaviour fails
  (`NewtonSolver::Mult`). A real line search (on the
  norm of the residual or the energy) is to be written there.
- **KINSOL** (interest 6): fixed point with Anderson, Picard, matrix free
  Newton-Krylov (`SetJFNK`), line search, reused Jacobian
  (`SetMaxSetupCalls`). The MFEM wrapper uses the CUDA vectors of
  SUNDIALS when MFEM is compiled with CUDA (checked in
  `linalg/sundials.hpp`). Requires the operator that integrates the behaviour. To be
  tried on the machine.
- **Prediction by extrapolation** (addition, interest 5): start of each step
  from the increment of the previous step (linear), or of the two
  previous ones (quadratic). Cheaper than the prediction by the tangent
  operator of mfem-mgis: no solve. Reduces the number
  of Newton iterations when the loading is regular.
- **Adaptive time step** (addition, interest 4): mfem-mgis has the framework
  (`AbstractTimeIncrementComputer`, substeps on failure), but the
  default computation goes directly to the end of the sequence. A computation of the
  step driven by the number of Newton iterations directly influences the
  total cost, hence the comparison of the variants.
- **Nonlinear solvers of PETSc** (addition, interest 4): NGMRES,
  Anderson, quasi-Newton (L-BFGS, `SNESQN`), NASM, with vectors on
  GPU. Useful as reference implementations for the items
  above, once the PETSc path of mfem-mgis is repaired (operator that
  integrates the behaviour).

## Matrix free preconditioners

Useful for option a (matrix free tangent). Option a of the plan
is preconditioned by an AMG built on the assembled elastic operator;
the smoothers below, which only need the diagonal or blocks,
allow a fully matrix free option a (alone, or inside a
p-multigrid).

- **Jacobi and Chebyshev on the diagonal** (interest 8):
  `OperatorJacobiSmoother` and `OperatorChebyshevSmoother` of MFEM
  (checked in 4.9). The diagonal comes from `AssembleGradDiagonalPA`
  (entry point of the partial assembly) or from a kernel of mfem-mgis.
  Chebyshev requires an estimate of the largest eigenvalue.
- **Jacobi by 3×3 nodal blocks** (addition, interest 7): invert, node by
  node, the 3×3 block of the diagonal of the tangent. More robust than the
  scalar Jacobi in elasticity, where the components of a node are
  strongly coupled. The `hooke` mini-application of MFEM 4.9 (matrix free
  nonlinear elasticity, on GPU) does it
  (`ElasticityDiagonalPreconditioner`, type `BlockDiagonal`): it is a
  direct model for mfem-mgis. The batched dense algebra of MFEM
  (`BatchedLinAlg`, since 4.8) can invert these blocks on GPU. Done on
  15/10/2026: the `NodalBlockJacobi` preconditioner, whose blocks are
  extracted from the assembled matrix or computed by kernels for the
  operator applied without assembly, inverted by the native kernels,
  cuBLAS or MAGMA.
- **Absolute value l1-Jacobi** (addition, interest 5): family of matrix free
  smoothers convergent by construction (MFEM 4.9,
  `miniapps/diag-smoothers`, examples in elasticity and in multigrid).
- **p-multigrid and LOR** (interest 7): `GeometricMultigrid` of MFEM with
  matrix free smoothers and an AMG on the first order space, or
  `LORSolver` (AMG on the low order discretization). Only useful
  for order 2 and higher. MFEM 4.10 improves the batched LOR assembly.

## AMG on GPU

### BoomerAMG (hypre)

- What runs on GPU (documentation of hypre): PMIS coarsening,
  direct (3), BAMG direct (15), extended (14), extended+i
  (6) and extended+e (18) interpolations, aggressive coarsening with multipass
  or two stage interpolation, Jacobi (7), l1-Jacobi (18),
  hybrid Gauss-Seidel (3, 4, 6), two stage Gauss-Seidel (11, 12) and
  Chebyshev (16) smoothers. The environment compiles hypre without unified memory: an
  option without GPU support fails instead of falling back to the CPU.
- MFEM chooses values compatible with the GPU when hypre is compiled
  for the GPU (PMIS, ext+i, l1-Jacobi, without aggressive coarsening;
  checked in `HypreBoomerAMG::SetDefaultOptions`, MFEM 4.9).
- **Settings exposed in mfem-mgis** (interest 8, done): options
  `StrengthThreshold`, `CoarseningType`, `AggressiveCoarseningLevels`,
  `InterpolationType`, `RelaxationType`, `RelaxationSweeps`, `CycleType`,
  `NumberOfCycles` and `MaximumNumberOfLevels` of `HypreBoomerAMG`, in addition
  to the strategy. Best setting measured on the beam:
  `AggressiveCoarseningLevels=1` (`results/L5-beam-2026-10-10-1034`).
- **Strategies of mfem-mgis on GPU** (addition, interest 8):
  - "System" (`SetSystemsOptions`, unknown-based approach) works on
    GPU: MFEM copies the table of the components to the GPU. It is the default
    of mfem-mgis without `Strategy` when the space has several components.
  - "Elasticity" calls `HypreBoomerAMG::SetElasticityOptions`, which
    stops the program with hypre built for GPU ("this method is not supported in
    hypre built with GPU support", checked in MFEM 4.9). This only
    concerns the byVDIM numbering: in byNODES, mfem-mgis already switches
    to "System" with a warning. Done: when hypre computes on the
    GPU, mfem-mgis also switches to "System" with a warning. The
    rigid body modes (interpolation vectors) stay on CPU with
    BoomerAMG.
  - Calling `SetSystemsOptions` again at each `SetOperator`, as
    `SystemHypreBoomerAMG` did, resets the number of aggressive
    coarsening levels to 0 and the threshold to 0.5 (`linalg/hypre.cpp:5231-5232`): exposed
    settings would be overwritten at each Newton iteration. Done:
    the class `ConfiguredHypreBoomerAMG` of mfem-mgis applies all the
    settings again, strategy then options, after each `SetOperator`.
  - `SetDefaultOptions` chooses the values for CPU or GPU at the
    construction of the AMG (`linalg/hypre.cpp:5029`): the AMG must be
    created after `mfem::Device`.
- **Filtering of the couplings between components** (addition, interest 4):
  `HYPRE_BoomerAMGSetFilterFunctions` (hypre 2.32) removes the couplings
  between components from the matrix used for the AMG: lower
  complexity, convergence to be measured in elasticity. Assumes interleaved
  unknowns (byVDIM). GPU to be checked.
- **Small systems: hypre on the host** (done on 12/10). Below about 10,000
  unknowns, a solve by GMRES and BoomerAMG on the GPU is
  dominated by the kernel launches: 100 ms for the 9,516 unknowns
  of the thermal problem of ex8, against 42 ms on two cores
  (`results/ex8-thermal-gpu-2026-10-11`). The option
  `--hypre-on-host` of `mfem_mgis::initialize` keeps the behaviours and the
  kernels on the GPU and makes hypre compute on the host: full ex7
  (2,634 unknowns) 6.7 -> 4.4 s, ex8 12.5 -> 11.6 s
  (`results/hypre-host-2026-10-12`). MUMPS on the host gives the same gain
  for the thermal problem of ex8 (`results/ex8-thermal-gpu-2026-10-11`).

### GAMG (PETSc)

Smoothed aggregation with near null space (rigid body modes): the reference
AMG for elasticity. With PETSc compiled with Kokkos or CUDA, the
matrices are on GPU (`aijkokkos`, `aijcusparse`); the share of the
setup done on GPU depends on the version, to be measured. PETSc 3.26
uses Kokkos 5.2.2 by default, the version of the environment. The PETSc path
of mfem-mgis crashed in earlier trials: to be repaired first
(see [the prerequisite of the nonlinear solvers](#prerequisite-an-operator-that-integrates-the-behaviour)).
MFEM 4.10 requires PETSc compiled with hypre. The block path on GPU of
the paper by Adams (arXiv 2606.24748) is not in PETSc 3.26.0, which only has
the `aijkokkos` and `aijcusparse` matrices on GPU. Interest 6.

State on 15/10/2026: `PETScLinearSolver` and `PETScPreconditioner` give
GAMG the block size and the rigid body modes when the unknowns are ordered
by vector dimension (15 conjugate gradient iterations instead of 51 on the
beam of the tests of mfem-mgis). With MFEM on a GPU, the solvers of PETSc
compute on the host (see [plan-solver-libraries.md](plan-solver-libraries.md)).

### AmgX (NVIDIA)

- `AmgXSolver` of MFEM: scalar matrix, double precision only
  (`AMGX_mode_dDDI`, comment "currently only double is supported"
  in `linalg/amgxsolver.hpp`), MPI modes, full reuse through
  `UpdateOperator` (coefficients replaced, no new
  setup). Its resources, shared by all its solvers, keep
  the communicator and the configuration of the first one: destroying the first
  solver makes the next ones fail; GPU 0 is imposed; the convergence
  status of AmgX is not accessible. **mfem-mgis does not use it.**
- **Wrapper specific to mfem-mgis** (done, option `enable-amgx`, AmgX given
  by `AMGX_DIR`): C interface of AmgX, preconditioner and solver `AmgX`
  (cycles of the multigrid), solver `AmgXPCG` (conjugate gradient of AmgX).
  - A single set of AmgX resources per process: they hold global
    objects (GPU memory pool, cuSPARSE and cuBLAS handles),
    released by the destruction of any resources, and keep
    the address of their communicator. They are created for the first
    solver on `MPI_COMM_WORLD` (AmgX makes its own copy of it for
    its exchanges) and destroyed with the last one; creation and destruction
    are not collective.
  - Communicator of a parallel solver: all the processes in the
    same order, or a single process (the matrix is then uploaded
    as a local matrix). Each process must own rows
    of the matrix, at most 2^31 - 1 rows in total: collective checks
    in `SetOperator`. After a failure of `SetOperator`, the solver refuses
    to solve until the next setup.
  - The solvers must be destroyed before `MPI_Finalize` (AmgX releases
    its communicators there).
  - Floors of AmgX: convergence is declared as soon as the norm of the
    residual goes below the maximum of 1e-20 and `e` times the norm of the
    initial residual (`e` is 1e-12 in double, 1e-6 in single), and
    immediately if the initial norm is below 1e-20. Without an
    absolute tolerance, `AmgXPCG` multiplies the right-hand side by a power of 2
    that brings the maximum of the absolute values of its components between
    1 and 2 (computed on all the processes; the norm itself
    would overflow at the extremes), then the solution by the inverse
    (exact operations): only the relative floor remains, and the
    relative tolerance must be at least `e`. The statistics
    printed by AmgX are then those of the scaled system.
    With an absolute tolerance (double precision only), the right-hand
    side is not modified: the absolute tolerance must be at least
    1e-20, and with an absolute tolerance alone the solve also stops
    when the norm of the residual has been multiplied by `e`. In single
    precision, the right-hand side is always scaled: a
    relative tolerance is required and an absolute tolerance refused. The
    preconditioner in single precision (default configuration) also
    scales the right-hand side, to stay within the range of
    single precision: a tolerance passed through `AmgXParameters`
    then applies to the scaled right-hand side.
  - Double copy of AmgX, until MFEM is rebuilt:
    libmfem.so contains the unpatched static AmgX of MFEM. A program
    must link libamgxsh.so (patched) before MFEM, as
    `MFEMMGISConfig.cmake` does; otherwise the classical multigrid fails
    ("Thrust failure"). A program that ends with a solver
    still alive (`std::exit`, static object) can abort or
    hang on a double release of the global objects of the two
    copies.
  - Current GPU of the process, convergence status of AmgX, errors
    reported through the return codes (without the parameter
    `exception_handling`, which on the contrary makes AmgX abort).
  - Parameters of AmgX by their name (`AmgXParameters`, for example
    `strength_threshold=0.8, interpolator=MULTIPASS`), which complete the
    default configuration (the one of MFEM).
  - Single precision (`Precision`: `Single`, mode `dFFI`, vectors
    converted on the GPU). The mixed precision `dDFI` (matrix in single,
    vectors in double) is no longer supported by AmgX 2.5
    ("Mixed precision modes not currently supported for CUDA 10.1 or
    later").
  - Still to do: 3×3 block matrices (ordering by components),
    partial reuse (`AMGX_solver_resetup`), setup without
    going through the host.
- AmgX 2.5.0 announces CUDA 13 (tested up to 13.0). Two complementary
  fixes followed on the main branch in 2026 (C++17 with
  CUDA 13.1 and later, tuples of CUDA 13): the environment uses a
  commit of this branch. With CUDA 13 on sm_120, the classical multigrid
  failed ("invalid device function"): patch of the spack
  recipe, reproducer `reproducers/amgx-classical-cuda13`.

### Ginkgo

Algebraic multigrid by aggregation (`Pgm`), distributed since 1.8. No
MFEM wrapper for the multigrid of Ginkgo (neither in 4.9 nor in 4.10): to be
configured through the API of Ginkgo. The MFEM wrappers only offer
the reference, OpenMP, CUDA and HIP executors (not SYCL). See also
[Several GPUs](#several-gpus) and [Mixed precision](#mixed-precision).

### MueLu (Trilinos)

Smoothed aggregation with rigid body modes on GPU through Kokkos. No MFEM
interface: heavy integration. Interest 3, to keep in mind as a
comparison point.

State on 15/10/2026: module of Trilinos written in mfem-mgis (MueLu with
the rigid body modes, FROSch, Ifpack2, Belos, Amesos2 with Tacho), in a
shared library compiled by the compiler of Kokkos, the unknowns ordered by
nodes being permuted for MueLu; Trilinos 17.2.1 built with CUDA
(environment `mmm-trilinos`); not yet built in mfem-mgis (see
[plan-solver-libraries.md](plan-solver-libraries.md)).

## Mixed precision

On the RTX 5080, computing in double precision is about 64 times
slower than in single. The matrix-vector products are limited by the
memory: in CSR, switching the values to single precision reduces the
bytes read by about a third (8 bytes per coefficient instead of 12,
the indices staying on 4 bytes). The dense computations (behaviours,
small blocks) are limited by the double precision compute. On
A100 or H100, the ratio is only 2.

Common idea: the outer Krylov solver and the residual stay in double, the
preconditioner works in single (or lower). The solution keeps
double precision.

- **Ginkgo**: block Jacobi with adaptive precision (`storage_opt` and
  `accuracy` of `JacobiPreconditioner` in the MFEM 4.9 wrapper,
  checked), CB-GMRES (compressed Krylov basis) and iterative refinement
  (`CBGMRESSolver`, `IRSolver` in the MFEM 4.9 wrapper), half
  precision and bfloat16 (1.9, 1.10), distributed matrices in mixed
  precision (1.10), multigrid in mixed precision.
- **hypre 3**: mixed precision functions (3.0: choice of the precision
  at run time, Krylov preconditioned in lower precision), on GPU
  since 3.2 (release notes). But the spack recipe of hypre forbids
  `precision=mixed` with CUDA (`conflicts("precision=mixed")` under
  `+cuda`), and MFEM requires a hypre in double (`hypre precision=double`
  in its recipe): the recipe would have to be changed, then the
  compatibility with MFEM checked.
- **AmgX**: mode `dDFI` (see above).
- **MFEM**: no mixed precision. MFEM is compiled entirely in
  single or in double (`precision`).

Interest 6, 8 on consumer GPU.

## Direct solvers on GPU

Interest for small cases, mixed configurations and as a reference.

- **cuDSS** (interest 5, done on 12/10): direct solver of NVIDIA, one or
  several GPUs, MPI. MFEM 4.10 wraps it (`CuDSSSolver`, with
  `SetReorderingReuse`: the symbolic analysis is kept when only the
  value of the coefficients changes, which is the case from one Newton
  iteration to the next). The spack recipe of MFEM 4.10 accepts cuDSS 0.6 to
  0.7.1; the environment installs cuDSS 0.7.1 for CUDA 13 and MFEM with
  `+cudss`. Solver `CuDSSSolver` of the factory of mfem-mgis, with
  `Symmetric`, `PositiveDefinite` and `ReuseAnalysis` (pattern compared on
  all the processes, MFEM solver recreated when the pattern or the size
  changes: MFEM only allows choosing the reuse before the first
  matrix and stops on a different size). On the elasticity matrix
  of ex8, a factorization takes 0.018 s against 0.203 s for
  MUMPS on the host; ex8 on the GPU 9.8 -> 4.4 s, then 2.6 s with the
  thermal problem also solved by cuDSS (9,516 unknowns: factorization 6 ms,
  solve 1 ms, against 52 ms by GMRES and BoomerAMG on the host)
  (`results/cudss-2026-10-12`). It does not pay off on large
  three-dimensional problems (beam r4, 111,843 unknowns: factorization 0.24 s
  against 68 ms per solve by PCG and BoomerAMG,
  `results/cudss-beam-2026-10-12`) nor on the periodic and
  nonsymmetric polycrystal (20,553 unknowns: 132 ms against 57 ms,
  `results/cudss-polycrystal-2026-10-12`). Each solution is checked by its
  residual (cuDSS does not report the failures of the factorization). The distributed
  analysis of cuDSS 0.7.1 fails for some matrices that MUMPS
  solves: rows reduced to their diagonal coefficient,
  tridiagonal matrix coupling the processes, some uneven distributions
  (`reproducers/cudss-distributed`, several processes on one GPU); a
  single process is not affected, and the diffusion matrices pass.
- **MUMPS on the host, reused analysis** (done on 12/10). In ex8 on
  the GPU, the 29 factorizations of MUMPS for the mechanics are the largest
  cost item. The parameter `ReuseAnalysis` of the solver `MUMPSSolver` of
  mfem-mgis keeps the analysis (ordering and symbolic factorization)
  as long as the pattern of the matrix does not change on any process:
  `SetReorderingReuse` of MFEM does not check the pattern, mfem-mgis
  compares it before each factorization. Factorizations −19 % on the GPU
  (one process) and −26 % on CPU (two processes), ex8 on the GPU
  11.1 -> 9.7 s, same results (`results/mumps-analysis-2026-10-12`).
  The other settings of MUMPS (orderings, BLR approximate
  factorization) gain nothing on this matrix. The same comparison of the pattern
  is also used by the cuDSS solver (above): with `SetReorderingReuse`,
  MFEM only checks the number of coefficients there (and stops if it
  changes), not their columns.
- **SuperLU_DIST** (interest 4): factorization on GPU (CUDA). MFEM 4.10
  adds `SuperLUSolver::SetDeviceOffload`.
- **STRUMPACK** (interest 4): `STRUMPACKSolverBase::EnableGPU` in MFEM
  4.9 (checked).

## Several GPUs

Useful for the study on a cluster (one process per GPU).

- **Distributed Ginkgo** (addition, interest 7): distributed matrices and vectors,
  distributed multigrid with PGM coarsening (1.8),
  two level Schwarz (1.10), pipelined CG (1.10), communications overlapped with computation
  (1.11). MFEM wrappers for the distributed solvers and preconditioners
  in MFEM 4.10 (example `examples/ginkgo/ex1p`, requires
  Ginkgo 1.11 compiled with MPI, which the environment has). Before 4.10,
  the wrapper only took sequential matrices.
- **Pipelined Krylov** (addition, interest 5): a single global reduction per
  iteration, overlapped with computation. Pipelined CG of Ginkgo, `KSPPIPECG`
  and related solvers in PETSc. Gain only with a large number of GPUs.
- **GPU-aware MPI** (tool, interest 5): exchanges directly between
  GPU memories. `Device::SetGPUAwareMPI` or variable `MFEM_GPU_AWARE_MPI`
  (since MFEM 4.9); `GroupCommunicator` on the device (MFEM
  4.10); variant `+gpu-aware-mpi` of hypre (not enabled). The OpenMPI of
  the environment is compiled with CUDA (internal hwloc without OpenCL), so
  that the code stays generic; on a single card, the effect is nil. On
  a cluster, the MPI of the cluster. Without GPU-aware MPI, the exchanges
  of the prolongation and those of hypre go through the host at each
  application: "no hidden copy" is then only
  reachable with a single process.

## Assembly and tangent operator

The plan provides three options for the tangent operator
([plan-mfem-mgis.md](plan-mfem-mgis.md#choices-made)). Additions:

- **Integration rule** (addition, interest 8): mfem-mgis integrates with
  `IntRules.Get(geom, 2 * OrderGrad)` (code of the generated integrators),
  i.e. 27 points on the first order hexahedron. A rule with 8 points divides
  by about 3.4 the cost of the behaviour, of the residual and of the tangent, on CPU
  as on GPU. It is the most direct cost lever of the hexahedral
  meshes. Done: parameter `IntegrationRuleOrder` of the
  discretization of mfem-mgis (option `--integration-rule-order` of ex4).
  Order 3 gives 8 points, which integrate exactly the elastic
  stiffness of affine first order hexahedra; orders 0 and 1 give
  one point, with spurious modes (documentation of mfem-mgis).
- **Assembly of the global matrix on GPU** (addition, interest 9): for
  option c on GPU, the element matrices computed on GPU must
  give the hypre matrix without going through the host. Two paths:
  - the sparse matrix filling of MFEM on the device
    (`ElementRestriction::FillJAndData`, assembly level `FULL` of the
    bilinear forms). It works on scalar degrees of
    freedom; the elasticity of MFEM uses it component by component.
    For a space with 3 components, to be checked, or to be done by blocks of
    components;
  - the IJ assembly of hypre on GPU (pointers in GPU memory, since
    hypre 2.20).
  The sparsity pattern, the `HypreParMatrix` and its communication scheme
  must be kept from one iteration to the next (only the values
  change), and the imposed degrees of freedom eliminated on the
  device (possible since MFEM 4.5). Otherwise, option c on GPU
  copies the tangent back to the host, or rebuilds its structure, at each
  iteration.
  Done in mfem-mgis: on the GPU, the kernels add the blocks of the
  element matrices directly into the local sparse matrix (pattern
  of the first assembly of MFEM; atomic additions, or fixed order with
  the parameter `TangentOperatorAssembly` set to `Ordered`), then the product
  `PᵀAP` and the elimination of the imposed degrees of freedom (`EliminateBC`) are
  done on the GPU. When hypre computes on the host, the global matrix
  keeps its pattern and its exchange scheme from one iteration to the next
  (`batch::GlobalStiffnessMatrix`,
  `results/global-matrix-end-to-end-2026-10-11`).
- **Element assembly** (addition, interest 6): keep the element
  matrices and apply the tangent element by element (level
  `ELEMENT` of MFEM for the bilinear forms). Intermediate between
  option a (nothing stored but `K` at the points) and option c (global
  matrix): more memory than a, no global assembly, diagonal
  and nodal blocks immediately available. `BatchedLinAlg` of MFEM (4.8) for the
  batched products; CUDA graphs to reduce the launch cost of the
  small kernels (small cases, several materials).
- The `hooke` mini-application of MFEM (matrix free nonlinear elasticity
  on GPU, Newton and CG, diagonal or block preconditioner)
  is the closest reference to option a in MFEM.

## Coupled problems

- **MGR of hypre** (interest 4): multigrid reduction for the
  monolithic systems (thermomechanics, micromorphic). Setup on GPU
  since 2.24, without unified memory since 2.29, many options since then (3.2: reuse of
  components, user-provided Schur complement, several cycles).
- Block preconditioners of MFEM (`BlockDiagonalPreconditioner`,
  Schur complement): generic, they run on GPU if the blocks are
  there.

## Measurement tools

Required by the performance studies, without being solvers.

- **Profiles** (interest 7): Nsight Systems and Nsight Compute (installed
  on the machine); NVTX markers in hypre (variant
  `+gpu-profiling`); Caliper (not yet in the environment). Without
  markers, a profile of hypre is hard to read.
- **Memory pools** (interest 6): Umpire for hypre (`+umpire`, in
  the environment) and for MFEM (`+umpire`). The setup of an AMG on
  GPU allocates many small arrays. MFEM compiled with Umpire takes
  the Umpire memory types by default: to keep in mind in the
  measurements.
- **`debug` device of MFEM**: it simulates the GPU on the CPU (kernels
  executed on the host, distinct "device" memory) and detects the
  accesses to an outdated copy. It cannot run a behaviour compiled
  for `cuda`: with it, use an `openmp` or `stlpar` behaviour and a
  location accessible by the host.
- GPU without HMM and Nsight Systems profile: check that no hidden copy
  takes place on the real GPU (already in the plan).

## Excluded (CPU only)

- Direct solvers MUMPS and UMFPACK (they stay in the environment: they are
  the current solvers of mfem-mgis, and the reference of the mixed
  configurations);
- `HypreEuclid`, `HypreParaSails`; elasticity options of BoomerAMG
  (`SetElasticityOptions`, see above);
- BDDC and FETI (PETSc), HPDDM: partial or no GPU support;
- (no longer excluded) libCEED through MFEM: backends `ceed-cpu` and
  `ceed-cuda` of MFEM, installed for comparison with the native partial
  assembly; also useful through Ratel.

## MFEM 4.10

The whole prototype uses MFEM 4.10 (decision of 09/10/2026; mfem-mgis
was used with 4.9 until then). MFEM 4.10 (01/09/2026) brings, for
this study (release notes and sources of 4.10):

- wrapper of cuDSS (`CuDSSSolver`);
- distributed solvers and preconditioners of Ginkgo;
- `SuperLUSolver::SetDeviceOffload`;
- `GroupCommunicator` on the device, with GPU-aware MPI;
- improved partial assembly kernels (diagonal of the elasticity);
- batched LOR assembly on highly connected meshes.

The facts checked on the sources of 4.9 and quoted in this document
(BoomerAMG on GPU, L-BFGS, `mfem::forall`, AmgX, Ginkgo, partial
assembly of the nonlinear forms) were checked again on those of 4.10.

## Mapping to the spack environment

| Technique | Package and variant of `mmm-gpu` |
|---|---|
| Krylov solvers of MFEM, matrix free smoothers, p-MG, LOR, `BatchedLinAlg` | `mfem@4.10.0 +cuda cxxstd=20` |
| Same kernels on multicore CPU, libCEED backend | `mfem +openmp +threadsafe +libceed`, `libceed +cuda +openmp` |
| Inexact Newton, fixed point and Anderson | done in mfem-mgis |
| Line search, preconditioned L-BFGS | to be written in mfem-mgis |
| BoomerAMG, MGR, FSAI, ILU, Krylov solvers of hypre | `hypre@3.2.0 +cuda +umpire` |
| Ginkgo, distributed included | `ginkgo@1.11.0 +cuda +openmp +mpi`, `mfem +ginkgo` |
| GAMG, SNES, pipelined Krylov | `petsc@3.26.0 +cuda +kokkos +openmp +hypre`, `kokkos-kernels`, `mfem +petsc` |
| KINSOL | `sundials +cuda +mpi`, `mfem +sundials` |
| AmgX | `amgx` (commit compatible with CUDA 13), `mfem +amgx` |
| cuDSS | `cudss@0.7.1-13 +mpi`, `mfem +cudss` |
| SuperLU_DIST, STRUMPACK on GPU | `superlu-dist +cuda`, `strumpack +cuda`, `mfem +superlu-dist +strumpack` |
| GPU-aware MPI | `openmpi@5 +cuda` (internal hwloc without OpenCL) |
| Memory pools | `umpire +cuda`, `mfem +umpire`, `hypre +umpire` |
| CPU direct solvers of mfem-mgis | `mfem +mumps +suite-sparse` |

Not yet in the environment: Caliper, hypre `+gpu-profiling`, hypre
in mixed precision, MueLu, Ratel (which requires PETSc 3.24).

## Sources

- Sources of MFEM 4.9.0 and 4.10.0 (archives of the spack recipe,
  checksums verified): `linalg/hypre.cpp` (`SetDefaultOptions`,
  `SetSystemsOptions`, `SetElasticityOptions`), `linalg/ginkgo.hpp`,
  `linalg/amgxsolver.hpp` and `.cpp`, `linalg/sundials.hpp`,
  `linalg/strumpack.hpp`, `linalg/solvers.hpp`, `linalg/cudss.hpp` (4.10),
  `miniapps/hooke`, `miniapps/diag-smoothers`, `CHANGELOG`, `INSTALL`.
- hypre: `CHANGELOG` of version 3.2.0, `src/parcsr_ls/HYPRE_parcsr_ls.h`
  (3.2.0), documentation of BoomerAMG (GPU options).
- Ginkgo: `CHANGELOG.md` of version 1.11.0.
- PETSc 3.24.6 and 3.26.0: `config/BuildSystem/config/packages/kokkos.py`
  and `kokkos-kernels.py`.
- AmgX: history of the main branch (CUDA 13 fixes of
  January and July 2026, after version 2.5.0 of December 2025).
- mfem-mgis: `src/LinearSolverFactory.cxx` (strategies of BoomerAMG),
  `src/NewtonSolver.cxx`, `src/NonLinearSolverBase.cxx`,
  `src/MultiMaterialNonLinearIntegrator.cxx`, generated integrators, read on
  09/10/2026; "done" states and references reviewed again on 13/10/2026 on
  `gpu-proto` b5e7d12fc4.
- Independent review of this document on 09/10/2026 (findings corrected).
- Spack recipes of the repository `spack-packages` at commit `3b7fdffc24`.
- Ratel: `https://gitlab.com/micromorph/ratel`, arXiv 2204.01722.
