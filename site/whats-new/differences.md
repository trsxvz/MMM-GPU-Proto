# Differences with the official versions

Each library of the prototype starts from a commit of its official
repository (its base). All changes are kept as a series of patches
against this base. This page lists the bases, what changes in each
library, how the series are organised and checked, and how they are
meant to be proposed upstream. The sources, in French, are `README.md`
and `docs/integration-amont.md`.

## Library by library

| Library | Official base | What changes | Patches |
|---|---|---|---|
| TFEL/MFront | master `ff0959ece` (12/10/2026) | `generic-parallel` interface for behaviours and models (6 backends, 14 configurations), extended interface for material properties, GPU code in TFEL and in the generated code, single precision (DSL option `real_type`), corrections of the existing code | 120 |
| MGIS | master `e8ef4b4` (13/10/2026) | batches of behaviours, models and material properties, execution policy, location of the arrays, memory resources, views and layouts, end of step on the GPU, single precision (option `MGIS_REAL_TYPE`), corrections | 70 |
| mfem-mgis | master `8fac67336e` (13/10/2026) | batch integration on CPU and GPU (data at integration points in the memory of MFEM), kernels for gradients, internal forces and tangent operator (assembled or matrix free; mechanics and stationary thermal problems), prediction by kernels, `FixedPoint` solver (fixed point, accelerations of MTest, NewtonK0), inexact Newton, BoomerAMG options, matrix free smoothers, reuse of the MUMPS analysis, cuDSS direct solver on the GPU, ordered assembly of the tangent operator (same results at each run on the GPU), execution settings, kernels compiled for CUDA or HIP (AMD GPU: compilation checked, not run), float or double precision given by MFEM, corrections | 176 |
| mfem-mgis-examples | master `600da70` of the reference repository of the examples | options of the examples (behaviour library, device, solvers, profiling), synchronisations on the host for the GPU, choice of the direct solvers, behaviour libraries for CPU and GPU chosen by `--execution`, tests on the GPU, examples in single precision | 31 |
| mm-opera-hpc | main `49147f2` of its reference repository | polycrystal: choice of the nonlinear solver (Newton, NewtonK0), of the adaptive tolerance and of AmgX, behaviour integrated on the GPU; bubble and cermet: options of mfem-mgis and of the behaviour, synchronisations on the host, profiling of the post-processings, direct solvers, behaviour libraries for CPU and GPU chosen by `--execution`, tests on the GPU, cases in single precision | 31 |

The execution settings of mfem-mgis are the option
`--execution cpu|gpu|mixed`, the placement of each phase in `mixed`
mode (`--finite-elements`, `--behaviours`, `--linear-solver-device`,
`--post-processings`), the placement of the behaviour of each material
(parameter `BehaviourExecutionLocation`), the counting of host and GPU
copies, and behaviour libraries for the CPU and the GPU. In single
precision, the tests run in float, and the stiffness kernels of the GPU
are compiled without contraction into fused multiply-add (FMA).

The details are in `docs/tfel.md` (TFEL), `docs/mgis.md` (MGIS) and
`docs/plan-mfem-mgis-detaille.md` (mfem-mgis).

## How the patch series are organised

The repository holds, for each library:

- `patches/tfel/`, `patches/mgis/`, `patches/mfem-mgis/`,
  `patches/mfem-mgis-examples/` and `patches/mm-opera-hpc/`: the series
  of patches from the base commit, produced by `git format-patch`;
- `tfel/`, `mfem-mgis/`, `mfem-mgis-examples/` and `mm-opera-hpc/`: the
  complete sources of the last commit, ready to compile, without the
  history of the official repository. MGIS is only kept as its patch
  series.

The scripts `scripts/sync-tfel.sh` and `scripts/sync-repository.sh`
regenerate a series from the local working repository and copy the
sources of its last commit. `patches/tfel/README.md` and
`patches/mgis/README.md` (in English) list the patches by subject.

A series is applied with `git am`. For TFEL, `git am --keep-cr` is
required: some files of TFEL (test behaviours) have CRLF line endings,
and some patches do not apply without this option.

## How the patch series are checked

`scripts/verify-patches.sh repository patches base branch` applies a
series on its base commit, in a temporary worktree of the local
repository, then compares the resulting tree with the tree of the
branch of the prototype. The check fails if a patch does not apply, or
if the trees differ.

The draft of the continuous integration (`docs/plan-ci-site.md`) has a
`patches` job that does the same on the upstream repositories: each
series is applied on its base, and the tree obtained is compared with
the sources at the root of the repository (`diff -r`). Simulated
locally for the five series, all apply and give back the tree of the
root.

## How the changes are to be proposed upstream

The changes are not yet proposed to the official repositories.
`docs/integration-amont.md` splits them into contributions that can
each be reviewed and tested alone. Some commits of the prototype mix
several subjects, so the contributions will be rebuilt from the final
state of the branches, file by file, rather than by picking commits.

### TFEL, seven contributions in dependency order

1. Corrections of the existing code, without GPU code. They may change
   results on CPU. Each one is covered by an existing or a new test.
2. GPU code in TFEL (mathematical and material headers). CPU results do
   not change.
3. Generated code compatible with the GPU in MFront. Results of the
   existing interfaces do not change.
4. The `generic-parallel` interface for behaviours and models, with a
   first set of backends: CUDA, HIP, C++ parallel algorithms and OpenMP
   on the host. The data contract is versioned from this contribution
   on, so that MGIS can check it.
5. Additional backends and compilers: Kokkos, SYCL, `nvc++` offloading
   (`stdpar` and OpenMP), clang CUDA and HIP. Most compiler workarounds
   go into this contribution.
6. Extensions of the interface for material properties.
7. Single precision (DSL option `real_type`): numerical type of the
   `generic` and `generic-parallel` interfaces, other interfaces refusing
   `float`, MTest checking the type, convergence criteria bounded at the
   generation. In double, the generated code only changes by its size
   symbol.

Each contribution is checked with the full TFEL test suite (`ctest`),
the `generic-parallel` comparisons for its configurations at `-O3`, no
change of the results of the existing interfaces (except the
corrections of the first contribution), and `clang-format` limited to
the modified lines.

### MGIS, eight contributions

They follow the fourth contribution of TFEL, whose contract they read:
corrections, batches of behaviours and models, memory (location of the
arrays and memory resources), material properties, layouts and copy
operations, profiling synchronisation, documentation, then single
precision (option `MGIS_REAL_TYPE`, after the seventh contribution of
TFEL). The CUDA test of
MGIS needs an NVIDIA GPU. It is disabled by default (option
`enable-generic-parallel-cuda-tests`), so that the continuous
integration of MGIS stays on CPU.

### mfem-mgis, thirteen contributions

The branch `gpu-proto` starts from the official master `8fac67336e`.
In dependency order:

1. corrections of the existing code, without GPU code;
2. features without GPU, useful on CPU (integration rule order, BoomerAMG
   options, inexact Newton, Jacobi and Chebyshev smoothers, reuse of the
   MUMPS analysis, profiling sections...);
3. build and configuration for the GPU (MFEM with CUDA or HIP, `--device`
   option, hypre memory pool, cuDSS direct solver...);
4. GPU path of the existing code: the whole test suite on the GPU, and
   the corrections found this way;
5. material data in the memory of MFEM;
6. batch integration of the behaviours, on the host and on the GPU;
7. mechanics by kernels (element geometry, gradients, internal forces,
   stiffness matrices, direct assembly into the sparse matrix, matrix
   free operator, other integrators);
8. other formulations by kernels (stationary and transient thermal
   problems, micromorphic damage, FBar, Faltus et al. regularisation);
9. the fixed point solver (`FixedPoint`);
10. post-processings by kernels;
11. execution settings (`--execution` and the placement options);
12. single precision: real type given by MFEM, MPI types, tolerances of
    the tests, and, in float, GPU stiffness kernels compiled without
    floating point contraction;
13. optional, later: Kokkos initialised by mfem-mgis, and the AmgX
    solvers once the AmgX correction and its spack recipe are merged.

Contributions 1 and 2 need neither the TFEL nor the MGIS of the
prototype. Contribution 3 needs the profiling synchronisation of MGIS,
contribution 5 the location of the arrays and the layouts of MGIS, and
from contribution 6 on, the `generic-parallel` interface of TFEL and
the batches of MGIS. Many commits of `gpu-proto` touch several
contributions, so these too will be rebuilt from the final state, piece
by piece.

Each contribution is checked with the full test suite of mfem-mgis on
CPU, sequential and on two processes; from contribution 4 on, the suite
on the GPU (`MFEM_MGIS_DEVICE=cuda`); from contribution 6 on, the suite
with the variables that force batches and kernels, on CPU and GPU; the
examples and the cases of mm-opera-hpc, with the same results on CPU
and GPU up to rounding; a test for each correction and each feature;
and `clang-format` limited to the modified lines.

### Examples and mm-opera-hpc

Their changes follow the contributions of mfem-mgis they use, and will
be proposed after them, by theme: options in each example and each
case, synchronisations on the host, choice of the nonlinear solver and
of the direct solvers, tests on the GPU, and in mm-opera-hpc, tests
compared with reference solutions and corrections found during review.

### Already proposed upstream

Corrections proposed to MFEM and AmgX during the work (open on
13/10/2026):

- MFEM #5604: shared resources of `AmgXSolver` (lifetime, sequential
  initialisation);
- MFEM #5605: order of the vertices of the boundary elements of a
  `SubMesh`, taken from the parent boundary element;
- MFEM #5615: `HypreParaSails::SetOperator` takes the communicator of
  the new operator (mfem-mgis works around this defect in its first
  contribution);
- AmgX #386: classical multigrid with CUDA 13 ("invalid device
  function", reproducer `reproducers/amgx-classical-cuda13`).
