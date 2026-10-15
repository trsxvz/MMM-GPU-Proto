# Continuous integration and website

The repository is public since 15/10/2026. The plan and the decisions
that led to this setup are in [plan-ci-site.md](plan-ci-site.md).

## Continuous integration (`.github/workflows/ci.yml`)

Triggers: each push and pull request on `main`, every night at 1:37 UTC
(to keep the build cache up to date with the latest dependencies), and
by hand (`workflow_dispatch`, the input `precisions` selecting the
numerical types of the job `cpu`, double by default). A new run of a
pull request cancels the previous one. The runs of `main` wait for the
run in progress, and only the latest waiting run is kept.

The machines of GitHub have no GPU: the jobs `cuda` and `hip` compile
only. The tests on the GPU stay local
(`scripts/mfem-mgis-gpu-suite.sh`).

| Job | Content |
|---|---|
| `patches` | the series of patches of `patches/` applied on their upstream base commits (TFEL, MGIS, mfem-mgis, mfem-mgis-examples, mm-opera-hpc); the patched sources must be identical to the sources of the repository |
| `cpu` | matrix over the precision (`double`, `float`): dependencies of the spack environment `deps/spack/ci-cpu` or `deps/spack/ci-cpu-float`; TFEL built from `tfel/`; MGIS from its patches (option `MGIS_REAL_TYPE`); mfem-mgis, mfem-mgis-examples and mm-opera-hpc built and tested |
| `cuda` | the same chain compiled with CUDA 13 (MFEM and hypre with CUDA, `cuda-nvcc` configuration of TFEL): environment `deps/spack/ci-cuda`, MFEM with the variant `cuda_ptx_only` of the spack repository `mmm_gpu` |
| `hip` | the same chain compiled with HIP for an AMD GPU (`gfx90a`): ROCm from its Python packages, MFEM built with HIP by `scripts/build-hip-amd.sh` |
| `update_build_cache` | single update of the index of the build cache once the other jobs have pushed, on `main` only |

### Caches

- Dependencies installed by spack: build cache on GHCR
  (`ghcr.io/trsxvz/mmm-gpu-proto-buildcache`), pushed from `main` only,
  then pruned of the packages no longer used. Composite actions of
  mfem-mgis (`mfem-mgis/.github/actions/`).
- MFEM compiled with HIP: `actions/cache`, keyed by the architecture,
  the version of ROCm, the hashes of the dependencies and the script.
- Compilations of the prototype: ccache, one cache per job and
  precision.

The dependencies are pushed to the build cache even when their
installation fails or stops (step `Install the dependencies` limited to
300 minutes): the packages already built are kept, and the next run
resumes from them. The push runs with the spack environment of the job
active (`SPACK_ENV`), which declares the repository `mmm_gpu` of the
packages `mfem` and `openmpi` of the job `cuda`.

### Durations

| Run | Job | Duration |
|---|---|---|
| first run, empty build cache | `cpu` (double) | 60 min: dependencies 40.6, TFEL 3.3, MGIS 2.1, mfem-mgis 8.6, mfem-mgis-examples 2.6, mm-opera-hpc 1.1 |
| fourth run, build cache filled | `cpu` (double) | 13 min |
| fourth run, build cache filled | `cpu` (float) | 14 min |
| fourth run | `hip` | 28 min |
| fourth run, empty build cache for CUDA | `cuda` | 3 h 41 min for the dependencies (hypre and MFEM compiled with one job) |

### Memory of the runners

A runner has 16 GB. MFEM 4.10 compiled with CUDA for `sm_120` needs more
than 12 GB for `fem/lor/lor_batched.cpp` in `ptxas`, even at `-O0`: the
variant `cuda_ptx_only` of the package `mfem` of `mmm_gpu` generates PTX
only (no `ptxas`), which is enough for a CI that does not run on a GPU.
hypre and MFEM are compiled with a single job (`cicc` alone takes 8.8 GB
on this file).

### Reproducing a job locally

Each job only uses the commands of the workflow: the spack environment
of the job (`spack -e deps/spack/ci-cpu install`), then the CMake
commands of its steps, in the same order (TFEL, MGIS, mfem-mgis,
mfem-mgis-examples, mm-opera-hpc).

## Website (`.github/workflows/pages.yml`)

Published on GitHub Pages at <https://trsxvz.github.io/MMM-GPU-Proto/> at
each push on `main` (and by hand). `scripts/build-site.sh` builds:

- the pages of `site/` (Sphinx, MyST Markdown) at the root, among which
  the tab "What's new with this prototype" (differences with the
  official versions, additions, choices, fixes, pitfalls, journal);
- the web pages of TFEL and MGIS, built by their CMake target `website`
  (pandoc and pandoc-crossref), from the sources of the prototype;
- the Sphinx and Doxygen documentation of mfem-mgis.

The workflow installs the tools, applies the patches of MGIS (only its
series of patches is in the repository), runs the script and deploys
the result (environment `github-pages`, default branch only).
