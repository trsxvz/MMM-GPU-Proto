# Plan: solvers and preconditioners of the major libraries in mfem-mgis

Status: survey of 15/10/2026, scope and order decided the same day (see
[Decisions taken](#decisions-taken-15102026)). Phase 0 done on
15/10/2026 (patches `patches/mfem-mgis/0186` to `0188`): descriptions,
factory of preconditioners, one module per library, generic test
`LinearSolverFactoriesTest`, documentation of the solvers in the user
guide of mfem-mgis. Phase 1 started: PETSc (patches `0189` to `0191`):
`PETScLinearSolver` and `PETScPreconditioner`, configured by the options
of PETSc, any preconditioner of mfem-mgis as shell preconditioner; Ginkgo
(patches `0192` to `0194`): 15 solvers and 10 preconditioners described
by configurations of Ginkgo, on the CPU and the GPU, local
preconditioners in the Schwarz preconditioner of Ginkgo in parallel. It
complements [gpu-solvers.md](gpu-solvers.md), organised
by technique, with a view by library: what each major library offers,
what MFEM already wraps, what mfem-mgis exposes today, and how to plug
the rest in a modular way.

## Goal

Any solver or preconditioner of the major libraries can be selected in
mfem-mgis by its parameters, as the solvers of MFEM and hypre are today,
and combined with the others where the libraries allow it. Some of them
only run on the CPU, others only on the GPU: the execution settings of
mfem-mgis (`cpu`, `gpu`, `mixed`) already refuse a solver that does not
compute where they place the linear solvers, with a clear message. The
same rule applies to the new ones.

## Current state

### mfem-mgis (branch `gpu-proto`)

| Kind | Entries | Options |
|---|---|---|
| Krylov solvers of MFEM | `CGSolver`, `GMRESSolver`, `FGMRESSolver`, `BiCGSTABSolver`, `MINRESSolver`, `SLISolver` | tolerances, iterations, verbosity, `KDim`, `Preconditioner` |
| Smoothers of MFEM (preconditioners) | `OperatorJacobiSmoother`, `OperatorChebyshevSmoother` | damping, order, power iterations |
| hypre | solvers `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `HypreBoomerAMG`; preconditioners `HypreBoomerAMG`, `HypreEuclid`, `HypreILU`, `HypreParaSails`, `HypreDiagScale` | BoomerAMG: strategy (`System`), coarsening, interpolation, relaxation, cycles, levels, threshold, aggressive coarsening; ILU: type, level of fill |
| Direct solvers | `MUMPSSolver` (parallel), `CuDSSSolver` (parallel, NVIDIA GPU), `UMFPackSolver` (sequential) | MUMPS: symmetric, positive definite, reuse of the analysis |
| AmgX (own wrapper) | solvers `AmgX`, `AmgXPCG`; preconditioner `AmgX` | full JSON configuration, precision |
| Wrapper | `AssembledOperator` | preconditioner updated every `UpdatePeriod` calls |
| PETSc | not in the factories: a global switch (`--petsc-configuration-file`) replaces the nonlinear solver by `PetscNonlinearSolver`, configured by a PETSc options file | |

The preconditioners are selected by a chain of `if` in
`getLinearSolverPreconditioner`, and the location of each solver (host,
device, as MFEM, as hypre) is deduced from its name
(`getLinearSolverPlacement`). Both must change before adding many
libraries (see [Architecture](#architecture)).

### MFEM 4.10 wrappers

| Library | Classes of MFEM 4.10 |
|---|---|
| MFEM | `CGSolver`, `GMRESSolver`, `FGMRESSolver`, `BiCGSTABSolver`, `MINRESSolver`, `SLISolver`, `OperatorJacobiSmoother`, `OperatorChebyshevSmoother`, `BlockILU`, `GSSmoother`, `DSmoother`, `AMGFSolver` (filtered AMG), `LORSolver`, multigrid classes, `OrthoSolver`, constrained solvers, `NewtonSolver`, `LBFGSSolver` |
| hypre | `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `HypreBoomerAMG`, `HypreSmoother` (Jacobi, l1-Jacobi, l1-GS, truncated l1-GS, lumped Jacobi, GS, OPFS, Chebyshev, Taubin, FIR), `HypreEuclid`, `HypreParaSails`, `HypreILU`, `HypreDiagScale`, `HypreAMS`, `HypreADS`, `HypreLOBPCG`, `HypreAME` |
| PETSc | `PetscLinearSolver`, `PetscPCGSolver`, `PetscPreconditioner`, `PetscBDDCSolver`, `PetscFieldSplitSolver`, `PetscNonlinearSolver`, `PetscODESolver` (all configured by the PETSc options database) |
| Ginkgo | `CGSolver`, `BICGSTABSolver`, `CGSSolver`, `FCGSolver`, `GMRESSolver`, `CBGMRESSolver`, `IRSolver`; `JacobiPreconditioner`, `IluPreconditioner`, `IluIsaiPreconditioner`, `IcPreconditioner`, `IcIsaiPreconditioner`, `SchwarzPreconditioner` (distributed), `MFEMPreconditioner` (any `mfem::Solver` seen by Ginkgo) |
| Direct solvers | `MUMPSSolver`, `SuperLUSolver` (SuperLU_DIST), `STRUMPACKSolver`, `STRUMPACKMixedPrecisionSolver`, `CuDSSSolver`, `PardisoSolver` and `CPardisoSolver` (oneMKL), `UMFPackSolver`, `KLUSolver` (SuiteSparse) |
| Others | `AmgXSolver`, `KINSolver` (SUNDIALS), batched dense algebra (`BatchedLinAlg`, native, cuBLAS/hipBLAS or MAGMA) |

### Libraries in the environments of the prototype

The MFEM of `mmm-gpu` (CUDA, double) has MPI, MUMPS, SuperLU_DIST
(CUDA), STRUMPACK (CUDA, SLATE), cuDSS, Ginkgo 1.11 (CUDA, OpenMP, MPI),
AmgX, PETSc 3.26 (CUDA, Kokkos, hypre; without MUMPS, HPDDM, STRUMPACK,
SuiteSparse, ML), SUNDIALS 7.9 (CUDA), SuiteSparse 7.14 (without CUDA),
libCEED and Umpire. Kokkos Kernels 5.2.2 is installed.

The MFEM of the three other environments (`mmm-double-cpu`,
`mmm-float-cpu`, `mmm-float-gpu`) only has MPI and MUMPS: Ginkgo,
PETSc, SuperLU_DIST, STRUMPACK and SuiteSparse are not there.

Available in spack and not installed: MAGMA 2.10 (CUDA and ROCm, with
its sparse part), Trilinos 17.2 (Tpetra, Belos, Ifpack2, MueLu,
Amesos2, ShyLU with FROSch), PaStiX 6.4, rocALUTION, oneMKL, HPDDM,
PSBLAS and AMG4PSBLAS. AMGCL is not in spack (header only library).

## Inventory by library

Interest from 1 to 10 for a nonlinear solid mechanics code on CPU and
GPU; effort in days of work, tests and review included.

### Libraries already partly plugged

**hypre** (CPU, NVIDIA and AMD GPU, MPI, float with hypre built for it).
To add: the smoothers of `HypreSmoother` as preconditioners (l1-Jacobi,
Chebyshev and l1-GS run on the GPU: the natural GPU smoothers); FSAI
(GPU since 2.30) and MGR (coupled problems), not wrapped by MFEM (API of
hypre); the Krylov solvers of hypre not wrapped by MFEM (BiCGSTAB,
COGMRES, LGMRES, the hybrid solver which switches from diagonal scaling
to AMG); an audit of the options of BoomerAMG (smoother types and
orders, F/C relaxation, maximum row sum, interpolation truncation,
nodal options). AMS and ADS are for H(curl) and H(div): out of scope.
Interest 8, effort 1 to 1.5.

**MUMPS** (CPU, MPI, OpenMP). To add: block low rank compression, choice
of the ordering, memory relaxation, number of threads, verbosity; the
direct solver used as a preconditioner (low rank factorization at low
accuracy). Interest 6, effort 0.5.

**cuDSS** (NVIDIA GPU, one or several GPUs, MPI). To add: the options of
the reordering, of the hybrid memory and hybrid execution modes, and
the multi-GPU mode. Interest 5, effort 0.5.

**AmgX** (NVIDIA GPU). Everything is reachable through the JSON
configuration. Possible addition: ready-made configurations by name
(classical AMG, aggregation AMG, PCG with each) to avoid writing JSON.
Interest 4, effort 0.5.

**MFEM** (CPU and GPU). To add: `BlockILU` and the sparse smoothers
(sequential or local blocks), `KLUSolver`, `AMGFSolver`, `LORSolver` and
the p-multigrid (high order), `OrthoSolver` (problems defined up to a
rigid motion). Interest 6, effort 1.5.

### Libraries wrapped by MFEM and installed, not plugged

**Ginkgo 1.11** (CPU reference and OpenMP, NVIDIA, AMD and Intel GPU,
MPI, all precisions, mixed precision). Solvers: CG, BiCG, BiCGSTAB, CGS,
FCG, GCR, GMRES, CB-GMRES (compressed basis), IDR, IR, MINRES,
pipelined CG, Chebyshev, multigrid; experimental direct solver (sparse
LU and Cholesky, CPU and GPU). Preconditioners: Jacobi (scalar and
block, adaptive precision), ILU and IC (exact, ParILU, ParIC, ParILUT,
ParICT), ISAI, SOR and Gauss-Seidel, distributed Schwarz, Pgm algebraic
multigrid (aggregation, distributed). MFEM wraps seven solvers and the
local preconditioners; the local ones need the distributed Schwarz
preconditioner in parallel; Pgm, the direct solver, GCR, IDR, MINRES,
pipelined CG and Chebyshev need the API of Ginkgo. The executor follows
the device of MFEM. Interest 8 (GPU, mixed precision, a second AMG on
the GPU), effort 2.

**PETSc 3.26** (CPU, NVIDIA GPU with cuSPARSE or Kokkos, AMD GPU with
Kokkos or hipSPARSE, MPI). Through the options database, one entry
gives access to all the Krylov solvers (CG, GMRES, FGMRES, LGMRES, GCR,
BiCGSTAB, TFQMR, MINRES, Chebyshev, pipelined variants) and all the
preconditioners (Jacobi, point block Jacobi, block Jacobi and additive
Schwarz with ILU or ICC, SOR, GAMG with the rigid body modes, BDDC,
field split, hypre, multigrid; with external packages when PETSc is
built with them: MUMPS, SuperLU_DIST, STRUMPACK, PaStiX, MKL PARDISO,
HPDDM for GenEO and Krylov recycling, ML). To plug: a linear solver and
a preconditioner whose options are given in the parameters of
mfem-mgis (translated into the options database, with a prefix per
instance), the matrices converted to the GPU types of PETSc when the
linear solvers run on the GPU, and the rigid body modes given to GAMG.
The global PETSc path of mfem-mgis (`PetscNonlinearSolver`) moves into
the factory of nonlinear solvers. Interest 9 (it brings tens of
methods at once), effort 2.

**SuperLU_DIST 9.2** (CPU, NVIDIA and AMD GPU offload, MPI; single
precision with mixed precision refinement since version 8): MFEM
`SuperLUSolver`. Interest 5, effort 0.5.

**STRUMPACK 8** (CPU, NVIDIA, AMD and Intel GPU, MPI; compression HSS,
BLR, HODLR, lossy): MFEM `STRUMPACKSolver` and
`STRUMPACKMixedPrecisionSolver`; usable as a preconditioner with a
compressed factorization. Interest 6, effort 0.5.

**SUNDIALS 7.9**: KINSOL (Newton with line search, Picard, fixed point
with Anderson acceleration) through `KINSolver`, as a nonlinear solver
of mfem-mgis. Its linear solvers (SPGMR, SPFGMR, SPBCGS, SPTFQMR, PCG and
the wrappers of KLU, SuperLU, cuSOLVER, MAGMA, Ginkgo) only serve inside
SUNDIALS. Interest 5, effort 1.

**SuiteSparse 7.14** (CPU): `KLUSolver` (MFEM), CHOLMOD (supernodal
Cholesky, GPU when built with CUDA), SPQR, ParU (parallel LU). No MFEM
wrapper for CHOLMOD and ParU. Sequential (one process, or local blocks).
Interest 4, effort 1.

### Libraries to install

**MAGMA 2.10** (NVIDIA and AMD GPU, one GPU). Two parts:
- the sparse part: Krylov solvers (CG, BiCGSTAB, GMRES, IDR, QMR, TFQMR,
  CGS) and preconditioners (Jacobi, block Jacobi, ILU with cuSPARSE or
  hipSPARSE, ParILU, ParILUT, ISAI) on one GPU, in single and double
  precision. Its authors develop Ginkgo, which covers the same methods
  with MPI: MAGMA sparse is useful as a comparison point and on one GPU;
- the dense part, batched: MFEM uses it for its batched dense algebra
  (`MFEM_USE_MAGMA`), which concerns the element matrices and the 3×3
  blocks of the kernels of mfem-mgis rather than the linear solvers.
Interest 4 (sparse), 5 (batched dense); effort 1.5 (wrapper of the
sparse part, CSR matrix on the GPU from the local hypre matrix).

**Trilinos 17.2** (CPU, NVIDIA and AMD GPU through Kokkos, MPI). Belos
(Krylov, block methods, GCRO-DR recycling), Ifpack2 (relaxations,
Chebyshev, ILU(k), ILUT, additive Schwarz, block relaxation), MueLu
(smoothed aggregation AMG with the rigid body modes, on GPU through
Kokkos), Amesos2 (KLU2, Basker, Tacho on GPU, SuperLU_DIST, MUMPS),
FROSch (two level Schwarz, GDSW coarse spaces, on GPU). No MFEM wrapper:
conversion of the hypre matrix to a Tpetra matrix, and adaptors between
`mfem::Solver` and `Tpetra::Operator`. Building Trilinos with Kokkos
for CUDA takes hours. Interest 7 (MueLu and FROSch are references for
elasticity), effort 4.

**oneMKL** (Intel and AMD CPU): PARDISO (shared memory) and the Cluster
Sparse Solver (MPI), wrapped by MFEM (`PardisoSolver`,
`CPardisoSolver`), and reachable through PETSc. Requires MFEM built with
them. Interest 5, effort 0.5 (plus the build of MFEM).

**PaStiX 6.4** (CPU, NVIDIA and AMD GPU through StarPU or PaRSEC, MPI,
low rank compression). No MFEM wrapper: through PETSc first
(`MATSOLVERPASTIX`), its own API later. Interest 4, effort 1.

**rocALUTION** (CPU OpenMP and AMD GPU, MPI). Krylov solvers (CG,
BiCGSTAB, GMRES, FGMRES, IDR, QMRCGSTAB, Chebyshev, mixed precision
defect correction), preconditioners (Jacobi, block Jacobi, Gauss-Seidel
and SGS, ILU(p), ILUT, IC, multicolored ILU and Gauss-Seidel, FSAI,
SPAI, approximate inverse Chebyshev, additive Schwarz) and algebraic
multigrids (smoothed, unsmoothed and pairwise aggregation,
Ruge-Stüben). The counterpart of AmgX on AMD GPU. Interest 5, effort 1.5.

**AMGCL** (header only; CPU OpenMP, NVIDIA GPU through CUDA or VexCL,
MPI). Smoothed aggregation and Ruge-Stüben AMG with 3×3 block values and
rigid body modes, relaxations (Jacobi, Gauss-Seidel, ILU(0), ILU(k),
ILUT, SPAI, Chebyshev), Krylov solvers (CG, BiCGSTAB, GMRES, LGMRES,
FGMRES, IDR(s)), mixed precision. Interest 5, effort 1.

**Kokkos Kernels 5.2** (installed): SPILUK, triangular solves, Gauss-
Seidel (point, cluster, two stage) on GPU, as preconditioners and
smoothers of the Krylov solvers of MFEM. Mostly reached through PETSc
(Kokkos matrices) and Trilinos. Interest 3, effort 1.

### Considered, not retained for now

- cuSPARSE and hipSPARSE (ILU(0), IC(0), triangular solves) and cuSOLVER
  sparse: covered by PETSc, Ginkgo, MAGMA and cuDSS.
- PSBLAS and AMG4PSBLAS (Fortran, AMG on GPU): to reconsider after
  Trilinos and Ginkgo.
- HSL (MA57, MA86, MA97, MI20): licence.
- Lis, ViennaCL, Paralution: little maintained or superseded
  (Paralution by rocALUTION).
- SLEPc (eigenvalues), libCEED and Ratel (operators, reference
  computations), H2Opus (through PETSc), ButterflyPACK (through
  STRUMPACK).

## Architecture

1. **Registries.** The factory of linear solvers exists; a factory of
   preconditioners is added on the same model (name, generator), and
   replaces the chain of `if`. The factory of nonlinear solvers exists.
2. **Description of each entry**, given at its registration: library,
   kind (Krylov solver, preconditioner, direct solver, multigrid),
   sequential or parallel, where it computes (host, device, as MFEM,
   as hypre), supported precisions, whether it needs an assembled
   matrix, which preconditioners it accepts (none, any `mfem::Solver`,
   those of its library), and its options (name, type, default value,
   description). The checks of the execution settings and of the
   options, the error messages and the tables of the documentation are
   derived from this description, not from the name.
3. **One module per library**: `src/<Library>Solvers.{hxx,cxx}`, with a
   function `declare<Library>Solvers` called when the library is
   available (`MFEM_USE_<LIBRARY>` for the libraries wrapped by MFEM, an
   option of mfem-mgis found by CMake for the others, as AmgX today).
   Adding a library is adding one module.
4. **Composition between libraries.** The Krylov solvers of MFEM accept
   any preconditioner. The Krylov solvers of hypre only accept those of
   hypre. PETSc accepts any `mfem::Solver` as a shell preconditioner.
   Ginkgo accepts any `mfem::Solver` through `MFEMPreconditioner`. The
   direct solvers, AmgX and cuDSS are self-contained. The description of
   each entry says what it accepts, and the factory refuses the other
   combinations with a message.
5. **Matrices.** mfem-mgis provides a `HypreParMatrix` (parallel) or a
   `SparseMatrix` (sequential), or a matrix free operator. Each module
   converts it to the format of its library (PETSc, Tpetra, Ginkgo,
   CSR on the GPU for MAGMA and rocALUTION), keeping the sparsity
   pattern when only the values change (Newton iterations). The matrix
   free operator is accepted by the Krylov solvers of MFEM, Ginkgo,
   PETSc (shell matrix) and Belos.
6. **Reuse.** The existing mechanisms (elastic operator kept by
   `NewtonK0`, `UpdatePeriod` of `AssembledOperator`, reuse of the
   analysis of the direct solvers) apply to all the entries.
7. **Tests.** One generic test runs each registered entry on the
   Laplace equation and on elasticity, sequential and parallel, on CPU
   and GPU, in double and float, with the options of its description,
   and checks the convergence and the residual. Each module adds its
   specific checks. A benchmark campaign (beam, polycrystal, cermet)
   then compares them.

## Decisions taken (15/10/2026)

Scope: the whole coverage of what MFEM 4.10 already wraps, with all the
solvers and preconditioners of each library (MFEM, hypre with FSAI, MGR
and its Krylov solvers not wrapped by MFEM, PETSc linear and nonlinear,
Ginkgo, MUMPS, SuperLU_DIST, STRUMPACK, cuDSS, UMFPACK, KLU, oneMKL
PARDISO and Cluster Sparse Solver, AmgX, KINSOL), plus MAGMA (sparse
and batched dense) and Trilinos (MueLu, Ifpack2, Belos, Amesos2 with
Tacho, FROSch). PaStiX through PETSc built with it. Not retained:
rocALUTION, AMGCL, Kokkos Kernels used directly, CHOLMOD and ParU, and
the libraries listed in [Considered, not retained](#considered-not-retained-for-now).

## Order of work

| Phase | Content | Effort |
|---|---|---|
| 0 | registries, descriptions, generic test; existing entries moved to modules | 1.5 |
| 1 | libraries wrapped by MFEM and present in `mmm-gpu` (no rebuild): PETSc (linear, then nonlinear), Ginkgo, hypre completions, SuperLU_DIST, STRUMPACK, KLU, MFEM completions, options of MUMPS and cuDSS, KINSOL | 8 |
| 2 | separate environment with the libraries which require a rebuild: oneMKL, PETSc with PaStiX, MUMPS and HPDDM, MAGMA | 3 |
| 3 | MAGMA (sparse and batched dense) | 2 |
| 4 | Trilinos (MueLu, Ifpack2, Belos, Amesos2, FROSch) | 4 |
| 5 | benchmark campaign and recommendations by case | 2 |

## Environments and continuous integration

- Phase 1 needs no rebuild: all these libraries are in the MFEM of
  `mmm-gpu`, which also computes on the CPU (device `cpu`).
- The MFEM of `mmm-gpu` is not rebuilt before the merge of the pending
  upstream pull requests (AmgX, MFEM): the libraries which change its
  build (oneMKL, PETSc with PaStiX, MUMPS and HPDDM, MAGMA, Trilinos) go
  into a separate environment.
- The continuous integration builds a reduced environment. A complete
  environment, with all the libraries, fits a nightly run rather than
  each push.
