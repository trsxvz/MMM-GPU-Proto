# Plan: choice of execution on CPU, on GPU or mixed

Status: plan of 13/10/2026 (task 9a). Option scheme chosen on 13/10/2026.
Implemented on 13/10/2026 (steps 1 to 5, see the last section):
mfem-mgis 75bfd92a54 to 1348a95a5e, examples dbf659f and the next one,
mm-opera-hpc 7416eb8 and the next ones, two independent reviews handled.

## Objective

Before this task, placing a computation on the GPU required several
independent settings, and each example exposed them in its own way:

- the `--device` option (or the `MFEM_MGIS_DEVICE` variable) configures
  the MFEM device;
- the variables `MFEM_MGIS_BEHAVIOUR_INTEGRATION`,
  `MFEM_MGIS_GRADIENTS_COMPUTATION`, `MFEM_MGIS_INNER_FORCES_COMPUTATION`,
  `MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION` and
  `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION` choose the computations by
  kernels;
- `--hypre-on-host` (or `MFEM_MGIS_HYPRE_ON_HOST`) keeps hypre on the
  host;
- the library of behaviours (`-l`) decides where the behaviours are
  integrated: a `generic-parallel` library of a GPU backend integrates
  them on the GPU, a `generic` library on the host; the examples only
  build the second one;
- ex1 and ex3 only read the environment variables, the MUMPS and PETSc
  variants of ex4 have no `-l` option, the Robin test has no option at
  all.

The goal: a single set of options, declared by mfem-mgis and understood
in the same way by all the examples, and all the possible placement
combinations, except for a real, documented technical impossibility.

## Options

Declared by `mfem_mgis::declareDefaultOptions`, read by
`mfem_mgis::initialize`:

- `--execution cpu|gpu|mixed` (`-x`):
  - `cpu`: no device, everything is computed on the host;
  - `gpu`: MFEM `gpu` device (CUDA or HIP depending on how it was
    built), everything is computed on the GPU;
  - `mixed`: `gpu` device, and each phase is placed by its sub-switch
    (GPU by default);
- in `mixed` only (refused otherwise, because contradictory):
  - `--finite-elements cpu|gpu`: MFEM data, gradients, internal forces,
    tangent;
  - `--behaviours cpu|gpu`: integration of the behaviours (also material
    by material, through the API);
  - `--linear-solver-device cpu|gpu`: hypre and the linear solvers;
  - `--post-processings cpu|gpu`;
- `--device` is kept for experts (MFEM names: `debug`, `omp`,
  `ceed-cuda`...). It is refused if it contradicts `--execution` (for
  example `--execution cpu --device cuda`).

Without `--execution`, nothing changes: the current behaviour (existing
options and variables) is kept, which leaves the tests of mfem-mgis
unchanged.

Order of priority, from the strongest to the weakest: parameter given to
the problem, environment variable, setting deduced from `--execution`,
default value. A parameter given explicitly therefore keeps control (a
test can always force a path).

Implemented differently (see [Implementation](#implementation-13102026)):
with `--execution`, the placement environment variables
(`MFEM_MGIS_DEVICE`, `MFEM_MGIS_HYPRE_ON_HOST`,
`MFEM_MGIS_BEHAVIOUR_INTEGRATION`, `MFEM_MGIS_*_COMPUTATION`) are
ignored, with a warning. The order is therefore: parameter given to the
problem, setting deduced from `--execution`, default value. Without
`--execution`, the variables apply as before.

## Mapping to the existing settings

| Phase | GPU | CPU |
|---|---|---|
| finite elements | `BehaviourIntegration=Batch`, gradients, internal forces and tangent `Kernels` on the device | `ElementByElement`, MFEM loops on the host (implemented: in `mixed` only; `Kernels` on the host in `cpu` mode) |
| behaviours | library of a GPU backend of `generic-parallel` | `generic` library (or `generic-parallel` on the host) |
| linear solver | hypre on the GPU | hypre on the host (`--hypre-on-host`) |
| post-processings | `PostProcessingsComputation=Kernels` on the device | `ElementByElement` on the host |

In `cpu` mode, the finite elements and the post-processings keep the
default path of mfem-mgis (element by element); the kernels remain
accessible through the existing parameters and variables, for the
measurements.

Implemented differently: in `cpu` mode, the finite elements are computed
by kernels on the host, which is faster, with a fallback to element by
element if the integrators do not support them. They remain element by
element when they are placed on the host in `mixed` mode. Only the
post-processings keep the element by element path on the host.

## Combinations

In `mixed` mode, the 16 combinations of the four phases are targeted.
What is already possible:

- finite elements on the host and behaviours on the GPU: the gradients
  computed element by element are written into the arrays of the
  materials, in the memory of MFEM, then read on the GPU by the batch
  integration;
- behaviours on the host and finite elements on the GPU: the arrays are
  synchronized on the host before the point by point integration, then
  read again on the GPU;
- solver on the GPU with a matrix assembled on the host: MFEM copies the
  matrix into the memory of hypre;
- post-processings on the host with a computation on the GPU:
  synchronization of the materials on the host, already done by the
  post-processings.

Real limits, to be refused with a clear message:

- a solver that only exists on one side: MUMPS, UMFPACK, PETSc (on the
  host); AmgX and cuDSS (on the GPU). The solver chosen by the
  application wins: a conflict with `--linear-solver-device` is an error,
  not a silent move;
- the placement of hypre is global to the process (fixed before its
  initialization): all the hypre solvers of a computation are on the same
  side;
- the MFEM device is global: in `mixed`, it is always configured, and
  "finite elements on the CPU" means MFEM loops on the host, with the
  vectors remaining in the memory of MFEM.

## Transfers

The transfers between the host and the GPU go through the memory of MFEM
(`Read`, `Write`, `HostRead`...) and through the synchronizations of the
materials. So that their overhead can be measured:

- a profiling section per synchronization of the materials
  (`Material::synchronizeOnHost`, `...OnDevice`), with the number of
  bytes;
- a section around the copies of the matrix to hypre when the assembly
  and the solver are not on the same side;
- a measurement table (beam r4, one combination per row: total time,
  time of the transfers) in `results/`.

Asynchronous copies: the memory of MFEM copies synchronously. An
asynchronous copy (CUDA or HIP stream, pinned memory) will only be added
if the measurement shows a gain, and only where the order of the
operations is guaranteed (for example the copy of the stresses back to
the host during the assembly of the tangent on the GPU).

## Behaviours

- A common CMake function, installed by mfem-mgis:
  `mfem_mgis_add_behaviours_library(NAME <name> SOURCES <.mfront...>
  [BACKENDS ...])`, which builds the `generic` library and, if MFEM uses
  a GPU, the one of the GPU backend of `generic-parallel` (`cuda` or
  `hip`), with a fixed naming. The implemented signature is different
  (see step 3 of the [Implementation](#implementation-13102026)).
- A function of mfem-mgis that chooses the library from the placement of
  the behaviours: `mfem_mgis::getBehavioursLibrary(ctx, name)`.
- The `-l` option of the examples forces a library (for the
  measurements).

## Steps

1. mfem-mgis: execution options and settings (`ExecutionSettings`),
   default values of the parameters deduced from the placement, refusal
   of contradictions; test of all the combinations on a small problem
   (same results as everything on the host) and of the refusals.
2. Profiling sections for the transfers.
3. CMake function for the behaviours and automatic choice of the
   library.
4. Examples of mfem-mgis-examples, then mm-opera-hpc, then ex1 and ex3:
   common options, behaviours built by the common function, environment
   variables removed from the scripts.
5. Measurement of the overhead of each combination (beam),
   documentation.

## Implementation (13/10/2026)

- Step 1: options and settings (`ExecutionSettings`,
  `getExecutionSettings`), read before the initialization of MPI and
  hypre; with `--execution`, the placement environment variables are
  ignored (warning), a parameter given to the problem keeps control. In
  `cpu` mode, the finite elements are computed by kernels on the host
  (faster: 12.2 s instead of 16.5 s on the beam, 8 processes), with a
  fallback to element by element if the integrators do not support them;
  the default batch integration has the same fallback. The
  `BehaviourExecutionLocation` parameter of `addBehaviourIntegrator`
  places the behaviours of a material.
- Step 2: copies of the arrays of the materials counted
  (`getMaterialArraysTransfers`) and displayed after the profiling table.
  The other transfers (material properties, external variables, MFEM
  vectors and matrices) are not counted: the total overhead is read in
  the times.
- Step 3: `mfem_mgis_add_behaviours_library` (installed CMake module) and
  `getBehavioursLibrary`. Implemented signature, without `NAME` or
  `BACKENDS`:
  `mfem_mgis_add_behaviours_library(<target> SOURCES <files.mfront>...
  [LIBRARIES <names>...] [INCLUDE_DIRECTORY <directory>] [DEPENDS
  <files>...] [EXCLUDE_FROM_ALL] [NO_GPU])`. The libraries are placed in
  `<directory>/cpu/` and `<directory>/gpu/`
  (`cmake/modules/MFEMMGISBehaviours.cmake` of mfem-mgis).
- Step 4: all the examples, mm-opera-hpc, ex1 and ex3; tests on GPU
  through `--execution`.
- Step 5: `results/execution-beam-2026-10-13-1432`. The copies of the
  arrays cost at most 0.42 s (from 0.2 % to 4.5 % of the computation
  depending on the placement); the cost of a mixed placement comes from
  the processor that computes. The asynchronous copies are therefore not
  done: the possible gain is small compared with the risk (order of the
  operations, pinned memory).
- Real limits (documented in the guide of mfem-mgis): solvers that only
  exist on one side; MFEM solvers on the GPU as soon as a GPU is used;
  global placement of hypre; cuDSS on the GPU with hypre on the host
  only through the expert options. With HIP, the kernels of mfem-mgis are
  compiled for HIP since task 9c: compilation and linking checked for an
  AMD GPU (gfx90a), not run on such a GPU.
- Leads: keep the states on the GPU after the end of step update (done
  on the GPU, but followed by a return to the host required by the
  contract of `MaterialArrays`, see `results/L7-transfers-2026-10-11`),
  with writes on the host synchronized explicitly (120 avoidable copies
  on the beam, 0.07 s); count the other transfers.
