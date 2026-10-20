# Plan for mfem-mgis

Status: choices made on 08/10/2026, first step (views API of MGIS) done
on 09/10/2026, detailed plan on 09/10/2026
([plan-mfem-mgis-detailed.md](plan-mfem-mgis-detailed.md)).
Work packages L0 to L10 are done (base master `a011d97d6f`, branch
rebased on 13/10/2026 onto master `8fac67336e`, on 16/10/2026 onto
master `20dcdf0a9b` (merge of PR #487), then on 20/10/2026 onto master
`5bbd5a6981` (merge of PR #500), repository `mfem-mgis/`;
status in [overview.md](overview.md#progress)). This plan relies on
TFEL and MGIS as described in [tfel.md](tfel.md) and
[mgis.md](mgis.md).

## Contents

1. [Spirit](#spirit)
2. [Starting point](#starting-point)
3. [Choices made](#choices-made)
4. [Target architecture](#target-architecture)
5. [Code organization](#code-organization)
6. [Portability](#portability)
7. [Solvers](#solvers)
8. [Mixed configurations](#mixed-configurations)
9. [Performance](#performance)
10. [Validation](#validation)
11. [Steps](#steps)
12. [Open questions](#open-questions)
13. [Risks](#risks)

## Spirit

The work on mfem-mgis has two goals, which are compatible:

- **explore**: measure, compare, try many variants;
- **propose** to the developers of mfem-mgis a clean, complete, flexible
  and modular code, in the spirit of TFEL and MGIS.

So everything can be chosen and nothing is excluded. The current path of
mfem-mgis remains the default path. The new paths are activated by
options. When several ways of doing exist, they are all implemented,
chosen by an option, and compared on the same cases. The API must be
ready for the experiments: better to provide too much than not enough,
as long as the code stays clean.

## Starting point

### mfem-mgis at the start (base `a011d97d6f`, 06/10/2026)

This section describes the code before the prototype. The current status
is in [overview.md](overview.md#progress).

- Element by element assembly on the host. One MFEM integrator per
  material (`MultiMaterialNonLinearIntegrator`, `AssembleElementVector`,
  `AssembleElementGrad`), then one class per hypothesis, kinematics and
  symmetry (`include/MFEMMGIS/*BehaviourIntegrator.hxx`: small and
  finite strains, plane stress and plane strain, 3D, nonlinear thermal,
  micromorphic, FBar, regularizations...).
- In each element: gradients at the Gauss points, integration of the
  behaviour point by point
  (`BehaviourIntegratorBase::performsLocalBehaviourIntegration`), then
  contributions to the residual and to the tangent matrix. This path
  requires a `generic` library: MGIS refuses a batch library on it, with
  a clear message.
- Data at the Gauss points: `PartialQuadratureFunction` (`std::vector`
  storage or external view). The material (`Material`) derives from the
  `MaterialDataManager` of MGIS.
- The `FiniteElementSpaceOrdering` parameter (`byNODES` by default, or
  `byVDIM`) chooses the order of the **degrees of freedom** of the finite
  element space. This is a different question from the layout of the
  **data at the Gauss points** (see [Choices made](#choices-made),
  point 2).
- Nonlinear solver: `Newton` only, specific to mfem-mgis (it does not
  derive from `mfem::NewtonSolver`). It gives the tangent to the linear
  solver at each iteration (`prec->SetOperator(J)`, where `prec` is the
  linear solver), which rebuilds the preconditioner at each iteration. No
  inexact Newton, no quasi-Newton, no fixed point; the correction is only
  divided by 2 if the behaviour fails (no real line search). The
  behaviour is integrated by an action of the solver, not by the `Mult`
  of the form: an external nonlinear solver would not integrate it.
- Linear solvers: see [Solvers](#solvers). Ginkgo is not connected.
- No use of the partial assembly of MFEM (`AssemblyLevel`), of
  `mfem::Device`, `ElementRestriction`, `QuadratureInterpolator` or
  `mfem::forall`.

### What MGIS now offers

- Batches of behaviours, models and material properties, on the arrays
  of the caller, without copy.
- Automatic switching of the MGIS functions that work on a
  `MaterialDataManager`.
- Execution policy (queue), checked location of the arrays, memory
  resources.
- Synchronous computation by default, or asynchronous: several materials
  can be submitted before a single synchronization.
- On 08/10, not yet: component by component layout of the data, update
  of the states on GPU (`update`, `revert`). Both are the first step of
  this plan, done on 09/10.

### Preliminary trials

Trials made outside this repository (directory `mfem-mgis-gpu`)
compared, on a beam in plasticity and 8 processes:

- Newton with PCG and AMG: 14.4 s on CPU;
- fixed point with the elastic operator and one AMG cycle: 8.1 s on CPU;
- on GPU, with a J2 behaviour written by hand in CUDA: 0.69 s, 0.34 s
  and 0.15 s for three variants.

On a polycrystal, the fixed point converges badly: Newton with GMRES and
AMG takes 22.6 s, the fixed point 67 s. No method wins everywhere.

## Choices made

1. **Entry point in MFEM: both.** The partial assembly of
   `NonlinearForm` (`AssemblyLevel::PARTIAL`: `AddMultPA`,
   `AssembleGradPA`, `AddMultGradPA`, `AssembleGradDiagonalPA`) and an
   operator specific to mfem-mgis. Choice by option, without duplicating
   the core: the batch kernels (gradients, integration by MGIS, residual,
   action of the tangent) are written once and serve both.
2. **Layout of the data at the Gauss points: both.** Point by point
   (components of one point contiguous) and component by component
   (values of one component contiguous, coalesced accesses on GPU). TFEL
   already accepts both. MGIS today only accepts point by point: it needs
   an API on views, designed first since the rest depends on it. The
   order of the degrees of freedom (`FiniteElementSpaceOrdering`) is an
   independent setting, already available.
3. **Tangent operator: three options, chosen by an option.**
   - a. matrix free: action `Σ Bᵀ K B` at the points, Krylov solver
     preconditioned by the AMG of the elastic operator;
   - b. fixed point with the elastic operator: residual only at each
     iteration;
   - c. assembled tangent: element matrices, then hypre global matrix.

   A fixed point nonlinear solver is also needed. The detailed plan adds
   a variant of option a, the element assembly (element matrices kept),
   and specifies the combinations of parameters of each option.
4. **Hypotheses, kinematics, integrators**: in the end, all of
   mfem-mgis. We start with 3D small strains. Each step remains testable
   on CPU and on GPU.
5. **Material properties and external variables**: all the ways, in the
   end (uniform, arrays computed by MFEM, `generic-parallel` properties
   evaluated by MGIS).
6. **Several materials**: handled, with the asynchronous computation of
   MGIS.
7. **Outputs, post-processings, restarts**: current behaviour by default
   (copy to the host).
8. **MFEM version: 4.10** for the whole prototype (decision of
   09/10/2026; mfem-mgis was used with 4.9). MFEM 4.10 brings cuDSS,
   distributed Ginkgo and `GroupCommunicator` on the device (see
   [gpu-solvers.md](gpu-solvers.md#mfem-410)). Compiled with CUDA,
   Ginkgo, PETSc, SUNDIALS, AmgX, cuDSS, SuperLU_DIST, STRUMPACK, MUMPS and
   SuiteSparse in `deps/spack/mmm-gpu`, with a GPU-aware OpenMPI.
9. **Fine measurements later.** Priority to a clean framework that works.
   The first measurements serve to compare the variants.
10. **Solvers**: choose the combinations that work on GPU (Krylov of
    MFEM or hypre, BoomerAMG on GPU, Jacobi), connect Ginkgo in the
    solver factory (CPU and GPU, chosen by the executor), keep the direct
    solvers on CPU for the mixed configurations.
11. **mfem-mgis-examples and mm-opera-hpc**: tracked in this repository
    from their base. Each example will receive parameters to choose CPU
    or GPU, and the other options (assembly level, tangent operator,
    solver, data layout).

## Target architecture

At each iteration, for each material, everything stays on the GPU:

1. **Gradients at the Gauss points.** The element values are extracted
   from the solution (`ElementRestriction` of MFEM), differentiated at
   the points (`QuadratureInterpolator`), then converted into gradients of
   the behaviour (strain, deformation gradient, temperature gradient),
   written directly into the end of step gradients of the material, in
   its layout.
2. **Integration of the batch** by MGIS, on these same arrays.
3. **Residual**: `Bᵀ σ w det J` per point, sum per element, then
   parallel assembly by MFEM on the GPU.
4. **Tangent operator**: one of the three options (a, b or c).
5. **Linear solvers**: hypre on GPU (PCG, GMRES, BoomerAMG), Ginkgo as
   an option (CUDA, HIP, OpenMP).
6. **End of step**: `update` and `revert` of MGIS on the GPU. mfem-mgis
   owns the memory, MGIS knows which arrays to copy. The transfers to the
   host are limited to the outputs.
7. **Post-processings**: on the host first, on the GPU next.

## Code organization

- **The current path does not change** (`AssemblyLevel::LEGACY`, element
  by element on the host). It remains the default path.
- **Common batch kernels.** Gradients, call to MGIS, residual, action of
  the tangent and diagonal are written once, generic over the gradient
  type, with `mfem::forall` (portable CUDA, HIP and CPU).
- **Two entry points** on these kernels:
  - partial assembly: the integrators of mfem-mgis implement the methods
    of `NonlinearFormIntegrator`. `AddMultPA(x, y)` computes the
    gradients, integrates the batch and adds `Bᵀ σ`. `AssembleGradPA`
    prepares the tangent operator (already computed by `AddMultPA` in the
    usual case). `AddMultGradPA(v, y)` adds `Σ Bᵀ K B v`.
    `AssembleGradDiagonalPA` gives the diagonal, for a Jacobi smoother;
  - operator specific to mfem-mgis: one `mfem::Operator` for the residual
    and one for the tangent, without going through `NonlinearForm`. It
    leaves more freedom (several materials submitted asynchronously,
    order of the kernels, fusion of kernels).
- The current classes (one per hypothesis, symmetry and kinematics) keep
  their role.
- **Data at the Gauss points** in the memory of MFEM (`mfem::Memory`),
  point by point or component by component. MGIS sees them without copy,
  through the views API to be designed. The host path reads the same
  memory.
- Rotations of orthotropic behaviours: batch rotation functions of MGIS,
  with matrices given as an array.
- Choices at run time: MFEM device (`mfem::Device`), library of the
  behaviour, entry point, data layout, option of the tangent operator,
  solver. MGIS checks that the memory of MFEM suits the library (location
  against memory space).

### Lead for later: the behaviour inside the element kernel

The planned path separates the gradients kernel, the MGIS batch and the
residual kernel. The data go through global memory between the three.
This is modular, but requires more bandwidth. Ratel, libCEED and dFEM on
the contrary do all the work of a point in a single kernel.

To try this variant, TFEL would have to offer a new output: a function
of one point, generated in a header, compiled into the kernel of the
calling code. With CUDA, calling from a kernel a function of another
library requires relocatable code. This output would also allow using
MFront behaviours in libCEED or Ratel. A large work package in TFEL, to
be compared with the separate path.

## Portability

| MFEM device | TFEL backend of the behaviour |
|---|---|
| CUDA | `cuda`, or `kokkos` with CUDA |
| HIP | `hip` |
| CPU (sequential, OpenMP) | `openmp`, `stlpar`, or `generic` library |

MFEM has no SYCL backend: the SYCL backends of TFEL serve other codes.
The HIP behaviours of TFEL run here on the NVIDIA platform; mfem-mgis
with an MFEM compiled with HIP for an AMD GPU compiles and links
(`scripts/build-hip-amd.sh`), without running on an AMD GPU.

## Solvers

### What existed in mfem-mgis at the start

Factory `src/LinearSolverFactory.cxx` of the base:

- MFEM: `CGSolver`, `GMRESSolver`, `BiCGSTABSolver`, `MINRESSolver`,
  `SLISolver` (sequential and parallel);
- parallel: `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `MUMPSSolver` (if
  MFEM has MUMPS);
- sequential: `UMFPackSolver` (if SuiteSparse);
- preconditioners: `HypreBoomerAMG`, `HypreEuclid`, `HypreILU`,
  `HypreParaSails`, `HypreDiagScale`, or none;
- PETSc: a path exists (`Config.cxx`), it crashed during earlier trials.

Added since by the prototype: solvers `FGMRESSolver`, `AmgX`, `AmgXPCG`
and `CuDSSSolver`, `ReuseAnalysis` parameter of `MUMPSSolver`,
preconditioners `OperatorJacobiSmoother`, `OperatorChebyshevSmoother` and
`AssembledOperator`, options of `HypreBoomerAMG` (see
[gpu-solvers.md](gpu-solvers.md)).

### Portability on GPU

Table of 08/10, completed on 13/10 with what was tried on the machine
(ILU, ParaSails, added solvers).

| Solver or preconditioner | GPU |
|---|---|
| Krylov of MFEM (CG, GMRES, BiCGSTAB, MINRES, SLI) | yes, if the operator and the preconditioner are on the GPU |
| Krylov of hypre (PCG, GMRES, FGMRES) | yes, with hypre compiled for CUDA or HIP |
| `HypreBoomerAMG` | yes with GPU hypre, for part of the options (PMIS coarsening, direct and extended interpolations, Jacobi, l1-Jacobi, hybrid Gauss-Seidel, Chebyshev smoothers); MFEM chooses values compatible with GPU |
| `HypreDiagScale` | yes |
| `HypreILU` | partly: with hypre 3.2.0, ILU(k) of level 1 or more crashes on the GPU; mfem-mgis takes level 0 there and refuses the others; ILUT works (`reproducers/hypre-ilu-gpu`) |
| `HypreEuclid` | no |
| `HypreParaSails` | no: refused by mfem-mgis when hypre computes on the GPU |
| `AmgX`, `AmgXPCG` (mfem-mgis) | NVIDIA GPU only |
| `CuDSSSolver` (mfem-mgis) | NVIDIA GPU only |
| `MUMPSSolver`, `UMFPackSolver` | no: direct solvers on CPU, matrix copied to the host |
| PETSc | only with PETSc compiled for GPU; path to be repaired |
| Ginkgo (`mfem::Ginkgo` wrappers) | CPU (OpenMP) and GPU (CUDA, HIP) with the same code, chosen by the executor (Ginkgo itself also has SYCL, not the MFEM wrappers); Krylov, Jacobi, ILU, IC, ISAI in the wrappers, not the AMG of Ginkgo; sequential matrices or matrix free operators; sequential matrices only in the wrappers of MFEM 4.9 (distributed in MFEM 4.10); requires MFEM compiled with Ginkgo |

The "elasticity" options of BoomerAMG (nodal coarsening, rigid body mode
vectors) do not work on GPU: `SetElasticityOptions` stops the program
with GPU hypre (checked in MFEM 4.9 and 4.10). The "Elasticity" strategy
of mfem-mgis is replaced by "System" (`SetSystemsOptions`, which works
on GPU), with a warning, when hypre computes on the GPU.

### Planned additions

The complete list of the solvers, preconditioners and techniques chosen
for the study (CPU and GPU, or GPU only), with their interest, what was
checked and the corresponding packages of the environment, is in
[gpu-solvers.md](gpu-solvers.md). The order in which they enter
mfem-mgis is given by the [detailed plan](plan-mfem-mgis-detailed.md).

Main reference: Ratel (`https://gitlab.com/micromorph/ratel`, arXiv
2204.01722). Matrix free p-multigrid with Chebyshev or Jacobi smoothers
and AMG on the assembled Q1, Newton-CG with steps solved to 1e-3, L-BFGS
with one multigrid cycle as initial inverse Hessian, line search, data at
the points kept for the tangent.

## Mixed configurations

The choices are made at run time, independently: MFEM device, library
of the behaviour, entry point, data layout, option of the tangent
operator, linear solver.

| Configuration | Data | Host and GPU copies per iteration |
|---|---|---|
| All CPU | host | none |
| All GPU | GPU | none, except the outputs |
| Behaviour on GPU, assembly and solver on CPU (heavy behaviours) | shared | gradients to the GPU, stresses and tangents to the host |
| Assembly on GPU, direct solver on CPU | GPU, matrix copied | matrix and right-hand side to the host |
| A non portable behaviour in a GPU computation | this material in shared memory, `generic` behaviour | data of this material only |

The copies of the mixes are explicit: `SHARED` location declared to
MGIS, or copies of MFEM. On the all GPU path, the absence of hidden copy
is checked in three ways: MFEM `debug` device, GPU without HMM, `nsys`
profile.

## Performance

- The data stay on the GPU from one step to the next. No managed memory
  on the all GPU path.
- On 08/10, the `generic-parallel` kernels had not yet been measured
  seriously. First measurements planned: J2 behaviour by batch against
  the hand-written J2 kernel of the preliminary trials, point by point
  layout against component by component, then heavier behaviours
  (implicit Norton, single crystal). Measured since: for example
  `results/L2-beam-2026-10-09`, `results/backends-gpu-2026-10-11` and
  `results/backends-polycrystal-2026-10-11`.
- Several materials: asynchronous submission and a single
  synchronization.
- MPI: one process per GPU (polycrystal stack), GPU-aware MPI if
  available.

## Validation

- Batch paths on CPU: same results as the current path on the test suite
  of mfem-mgis, for both entry points and both layouts.
- GPU against CPU on the examples (beam in plasticity, finite strains,
  thermal, polycrystal), with a rounding tolerance.
- Comparison of the variants: time per step and per phase, multicore CPU
  against GPU, the three options of the tangent operator on the same
  cases.

## Steps

The detailed split of these steps into work packages (files touched,
tests, completion criteria, proposed parameters) is in
[plan-mfem-mgis-detailed.md](plan-mfem-mgis-detailed.md).

1. **MGIS: API on views** for both layouts of the data at the points
   (behaviours, models, material properties; synchronous and
   asynchronous), and update of the states on GPU. Tests on CPU and
   CUDA, independent review. Done on 2026-10-09 (see
   [mgis.md](mgis.md#views-and-layouts)).
2. **MFEM 4.10 with CUDA** in the `mmm-gpu` environment, with all the
   solver libraries (few jobs at a time: `ptxas` takes 22 GB for one
   file). Done.
3. **mfem-mgis on CPU**: data at the points in the memory of MFEM, then
   batch kernels (gradients, MGIS, residual, action of the tangent), 3D
   small strains, results identical to the current path. Both entry
   points and both layouts. Done (work packages L1 to L3).
4. **Tangent operator**: the three options and the fixed point; GPU
   compatible solvers; Ginkgo in the factory. Done, except Ginkgo in the
   factory.
5. **GPU**: MFEM device, same results as on CPU up to rounding, without
   hidden copy (`debug` device, GPU without HMM, `nsys`). Done (work
   package L6).
6. **End of step without copy to the host**, several materials
   asynchronously. Partly: the end of step is done on the GPU, but the
   states then come back to the host, as required by the contract of
   `MaterialArrays` (`results/L7-transfers-2026-10-11`); mfem-mgis does
   not use the asynchronous computations of MGIS.
7. **All of mfem-mgis**: hypotheses, kinematics, orthotropic rotations,
   material properties, external variables. Done (work package L7,
   leftovers in the detailed plan).
8. **Examples** of mfem-mgis-examples and mm-opera-hpc with the CPU/GPU
   options and the other options; measurements and comparison of the
   variants. Done (work package L9).
9. **Documentation** of the prototype and preparation of the upstream
   proposals ([upstream-integration.md](upstream-integration.md)). Done;
   nothing is proposed yet.

At each step: green tests (targeted pass), independent review before
declaring the step finished, log up to date, commit, synchronization,
push, bundle.

## Open questions

Questions of 08/10, settled since:

- Owner of the memory at the points on GPU: `mfem::Memory` in the
  `PartialQuadratureFunction`, seen by MGIS. The end of a time step
  without copy to the host is settled on the MGIS side: copy on the GPU
  by `update` with operations provided by the caller. Swapping the arrays
  is ruled out: it is not equivalent to `update` (a prediction does not
  write the end of step state). Settled as follows: arrays of the
  materials in the memory of MFEM (`MaterialArrays`), end of step by
  `update` with the operations of mfem-mgis (`MFEMFieldOperations`).
- Upstream proposal: a parallel path activated by option, acceptable to
  the maintainers of mfem-mgis. Split proposed in
  [upstream-integration.md](upstream-integration.md).

## Risks

- Compilation of MFEM 4.10 with CUDA 13: duration and memory.
- mfem-mgis has only been used with MFEM 4.8 and 4.9: compatibility with
  4.10 to be checked first.
- Scope: mfem-mgis has many integrators. The batch path is generic and
  starts with the most used cases.
- Cost of building the AMG on GPU, memory.
