# MMM-GPU-Proto

Prototype of the MFront (TFEL), MGIS and mfem-mgis chain, on CPU and on
GPU.

## In short

- **TFEL/MFront** generates, for each behaviour, model and material
  property, functions that process a whole batch of integration points.
  The same MFront file compiles for CUDA, HIP, Kokkos, SYCL, OpenMP and the
  parallel algorithms of C++. The results are systematically compared to
  those of the `generic` interface.
- **MGIS** loads these libraries and computes batches without copying the
  data. Its existing functions switch to batch computation by themselves.
  The caller chooses the memory, the execution queue, and synchronous (by
  default) or asynchronous computation. The data can be stored component
  by component, or belong to the caller and be described by views,
  including on GPU, including the end of the time step.
- **mfem-mgis** integrates the behaviours by batches, on CPU and on GPU.
  The data at the points live in the memory of MFEM. Gradients, inner
  forces, tangent operator (assembled or matrix free), prediction,
  thermal, micromorphic, FBar and post-processings are computed by
  kernels. Solvers suited to the GPU are added to Newton: accelerated
  fixed point, NewtonK0, inexact Newton, cuDSS direct solver. The whole
  test suite passes on the GPU, except the two tests that require the
  absence of a GPU. The plan followed is split into work packages
  ([docs/plan-mfem-mgis-detailed.md](docs/plan-mfem-mgis-detailed.md)).
- **mfem-mgis-examples** and **mm-opera-hpc**: all the examples and the
  cases (polycrystal, cermet, bubble) run by batches and on the GPU. On an
  RTX 5080, the beam takes 1.96 s with NewtonK0, against 9.5 s on 8 cores
  with the same method (13.6 s with Newton), and ex8 2.6 s against 19.2 s
  on 2 cores
  ([docs/overview.md](docs/overview.md#progress)).
- The whole chain also works in single precision (float), the real type
  being the one of the installed MFEM: on the beam, float is 1.11 times
  faster than double on CPU and 1.37 times on GPU
  ([docs/plan-precision.md](docs/plan-precision.md)).
- The changes are not yet proposed to the official repositories. Their
  split into contributions is described in
  [docs/upstream-integration.md](docs/upstream-integration.md).

## Differences with the official versions

| Library | Official base | What changes | Details |
|---|---|---|---|
| TFEL/MFront | master `ff0959ece` | `generic-parallel` interface of the behaviours and of the models (6 backends, 14 configurations), extended interface of the material properties, GPU code in TFEL and in the generated code, single precision (DSL option `real_type`), displacement gradient given to the finite strain behaviours, logarithmic strain handler exact in float, axial component computed by the finite strain behaviours in plane stress, fixes of the existing code; 125 patches | [docs/tfel.md](docs/tfel.md) |
| MGIS | master `e8ef4b4` | batches of behaviours, models and material properties, execution policy, location of the arrays, memory resources, views and layouts, end of the time step on GPU, single precision (option `MGIS_REAL_TYPE`), displacement gradient (option `gradient` of the finite strain behaviours), fixes; 78 patches | [docs/mgis.md](docs/mgis.md) |
| mfem-mgis | master `8fac67336e` | batch integration on CPU and GPU (data at the points in the memory of MFEM), kernels of the gradients, of the internal forces and of the tangent operator (assembled or matrix free; mechanics and stationary thermal), prediction by kernels, `FixedPoint` solver (fixed point, accelerations of MTest, NewtonK0), inexact Newton, BoomerAMG options, matrix free smoothers, reuse of the MUMPS analysis, cuDSS direct solver on the GPU, ordered assembly of the tangent operator (same results at each pass on the GPU), execution settings (`--execution cpu|gpu|mixed`, placement of each phase in `mixed` mode by `--finite-elements`, `--behaviours`, `--linear-solver-device` and `--post-processings`, placement of the behaviour of each material by the parameter `BehaviourExecutionLocation`, host/GPU copies counted, libraries of behaviours for the CPU and the GPU), kernels compiled for CUDA or HIP (AMD GPU: compilation checked, not run), float or double precision given by MFEM (tests in float; in float, stiffness kernels of the GPU without FMA contraction), displacement gradient given to the finite strain behaviours (parameter `DisplacementGradient`), linear solvers and preconditioners registered by one module per library with their descriptions (location, kind, accepted preconditioners, options), Krylov solver and preconditioner of PETSc configured by the options of PETSc, fixes; 191 patches | [docs/plan-mfem-mgis-detailed.md](docs/plan-mfem-mgis-detailed.md), [JOURNAL.md](JOURNAL.md) |
| mfem-mgis-examples | master `600da70` (repository `latug0/mfem-mgis-examples`) | options of the examples (library of behaviours, device, solvers, profiling), synchronisations on the host for the GPU, choice of the direct solvers, behaviours for the CPU and the GPU chosen by `--execution`, tests on GPU, examples in float; 31 patches | `patches/mfem-mgis-examples/` |
| mm-opera-hpc | main `49147f2` (repository `rprat-pro/mm-opera-hpc`) | polycrystal: choice of the non linear solver (Newton, NewtonK0), of the adaptive tolerance and of AmgX, behaviour integrated on the GPU; bubble and cermet: options of mfem-mgis and of the behaviour, synchronisations on the host, profiling of the post-processings, direct solvers, behaviours for the CPU and the GPU chosen by `--execution`, tests on GPU, cases in float, macroscopic displacement gradient (option `--displacement-gradient`); 34 patches | `patches/mm-opera-hpc/`; measurement of the polycrystal: [results/L5-polycrystal-2026-10-10-0502/README.md](results/L5-polycrystal-2026-10-10-0502/README.md) |

## Documentation

To read in this order:

1. [docs/overview.md](docs/overview.md): goal, role of each library,
   principles, vocabulary, progress.
2. [docs/tfel.md](docs/tfel.md): changes of TFEL and MFront, and their
   justifications.
3. [docs/mgis.md](docs/mgis.md): changes of MGIS.
4. [docs/validation.md](docs/validation.md): tests, sweeps, tools and
   results.
5. [docs/upstream-integration.md](docs/upstream-integration.md): split into
   contributions for the official repositories.
6. [docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md): plan for mfem-mgis,
   with the questions to decide, and its split into work packages
   ([docs/plan-mfem-mgis-detailed.md](docs/plan-mfem-mgis-detailed.md)).
7. [docs/plan-execution.md](docs/plan-execution.md): choice of the
   execution on CPU, on GPU or mixed (`--execution`), done.
8. [docs/plan-precision.md](docs/plan-precision.md): float or double
   precision on the whole chain, done.
9. [docs/plan-ci-site.md](docs/plan-ci-site.md): continuous integration and
   site of the prototype (decisions taken on 13/10; repository public
   since 15/10/2026, CI and site started at each push on main).
10. [docs/gpu-solvers.md](docs/gpu-solvers.md): CPU and GPU solvers,
    preconditioners and techniques retained for the study of
    mfem-mgis.
11. [docs/displacement-gradient.md](docs/displacement-gradient.md):
    displacement gradient given to the finite strain behaviours, for the
    precision of small strains in float.
12. [JOURNAL.md](JOURNAL.md): day by day progress.

The documentation of the interfaces, in English, is also in the tree of
each library: `tfel/docs/web/generic-parallel-*.md` for TFEL,
`mgis/docs/web/` for MGIS.

## Content of the repository

| Path | Content |
|---|---|
| `docs/` | documentation of the prototype; `docs/history/` keeps the outdated working documents (review of 06/10, initial design of MGIS) |
| `JOURNAL.md` | dated log of the work |
| `tfel/` | complete sources of the modified TFEL, ready to compile (copy of the last commit, without the history of TFEL) |
| `patches/tfel/` | the same changes as a patch series from the master of TFEL (`git am --keep-cr`) |
| `patches/mgis/` | changes of MGIS as a patch series from its master |
| `mfem-mgis/`, `mfem-mgis-examples/`, `mm-opera-hpc/` | sources of these repositories (copy of the last commit, without their history) |
| `patches/mfem-mgis/`, `patches/mfem-mgis-examples/`, `patches/mm-opera-hpc/` | their changes as patch series from their base |
| `build-tfel.sh` | configuration, build and installation of TFEL |
| `build-mgis.sh`, `build-mfem-mgis.sh` | build and tests of MGIS, then of mfem-mgis, of the examples, of mm-opera-hpc and of the cases of `cases/`, against the TFEL installed by `build-tfel.sh` |
| `cases/` | measurement cases (beam, thermal, polycrystal, bubble, cermet, ssna303, examples) and measurement scripts `measure-*.sh`, see `cases/README.md` |
| `scripts/` | environment (`env.sh`), synchronisation of the sources and of the patches (`sync-tfel.sh`, `sync-mgis.sh`, `sync-repository.sh`) and their check (`verify-patches.sh`), test matrix, suite of mfem-mgis on the GPU (`mfem-mgis-gpu-suite.sh`), build for an AMD GPU (`build-hip-amd.sh`), GPU headers of MGIS (`check-mgis-gpu-headers.sh`), host and GPU annotations of the headers of TFEL (`clang-cuda-header-scan.sh`, `annotate-function-templates.py`), compilation and runtime sweeps, roundoff check, memory guard of the long compilations |
| `scratch/` | tracked development tools: complete validation (`validate3.sh`), tests with the development tree (`dev-test.sh`, `targeted-test.sh`), development configurations (`b2/`) |
| `reproducers/` | small standalone programs that check a precise point |
| `results/` | results of the matrices and of the sweeps; compiler defects and their reproducers (`results/nvhpc-bugs/`) |
| `site/` | pages of the site of the prototype, in English (Sphinx, MyST); the complete site, with the documentations of TFEL, MGIS and mfem-mgis, is built by `scripts/build-site.sh` and published by `.github/workflows/pages.yml` |
| `deps/spack/` | spack environments of the prototype (`mmm-gpu`; `mmm-float-cpu` and its double twin `mmm-double-cpu`; `mmm-float-gpu`; `ci-cpu` and `ci-cpu-float`) and complementary recipes |
| `deps/spack-sycl/` | specification of the LLVM of the old environment |

The working repositories are local clones. `tfel/`, `mfem-mgis/`,
`mfem-mgis-examples/` and `mm-opera-hpc/` are also copied here, without
their history, by `scripts/sync-tfel.sh` and `scripts/sync-repository.sh`.
`tfel-dev/` (development branch of TFEL) and `mgis/` are not copied: MGIS
is tracked by its patch series.

## Build

The paths of the scripts are those of the development machine.
`scripts/env.sh` describes the environment: compilers, Kokkos,
AdaptiveCpp, HIP, CUDA.

TFEL, from the sources of this repository:

```bash
./build-tfel.sh install   # builds tfel/ and installs it in install/tfel-mmm-gpu
```

TFEL, from the official repository:

```bash
git clone <official repository of TFEL> tfel && cd tfel
git checkout -b generic-parallel ff0959ece
git am --keep-cr ../patches/tfel/*.patch   # some files of TFEL use CRLF
```

MGIS, from the official repository, against this TFEL:

```bash
git clone <official repository of MGIS> mgis && cd mgis
git checkout -b gpu-proto e8ef4b4
git am ../patches/mgis/*.patch
cd .. && source scripts/env.sh
cmake -S mgis -B build/mgis -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL" -DCMAKE_CXX_COMPILER=g++-14 \
  -Denable-generic-parallel-cuda-tests=ON \
  -DCMAKE_CUDA_COMPILER="$(which nvcc)" -DCMAKE_CUDA_HOST_COMPILER=g++-14
cmake --build build/mgis --target check -j16   # builds and runs the tests
```

The scripts of the root do these steps with the local repositories:
`./build-mgis.sh` builds MGIS (directory `mgis/`), tests it and installs
it in `install/mgis`; `./build-mfem-mgis.sh` builds, tests and installs
mfem-mgis, then builds and tests the examples, mm-opera-hpc and the cases
of `cases/` (see the header of each script). The dependencies (MFEM,
hypre, solvers, compilers) come from the spack environment `mmm-gpu`,
described in [deps/spack/README.md](deps/spack/README.md) and installed by
`deps/spack/install-mmm-gpu.sh`.

The tests and how to read them are described in
[docs/validation.md](docs/validation.md).

Compilation for an AMD GPU (compilation and link only: the integrated AMD
GPU of the development machine is accessible only to the members of the
group render): `scripts/build-hip-amd.sh` builds hypre (for the host),
MFEM compiled with HIP for gfx90a, mfem-mgis and the examples with this
MFEM, and runs the tests that do not need a GPU. ROCm is installed
without administrator rights (Python packages of ROCm, see the header of
the script).
