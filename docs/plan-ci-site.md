# Plan: continuous integration and site of the prototype

Status: plan of 13/10/2026 (tasks 9f and 9g), decisions taken on
13/10 (see below); drafts of the CI on the CPU and of the site on
14/10, triggered by hand, not yet run. The two tasks
depend on the visibility of the repository (decision E1).

## Constraints

- The repository is private. With a free GitHub account (conditions to be
  checked at the time of the choice): about 2,000 minutes of Actions per month
  for private repositories, about 500 MB of storage for private
  packages (the spack cache of the CI of mfem-mgis weighs nearly 2 GB), and
  GitHub Pages is not available for a private repository. A public repository
  does not have these limits, but publishes all its content.
- The machines of GitHub have no GPU: the CI can compile for CUDA
  and for HIP, not run on a GPU.
- Costs measured on the CI of mfem-mgis (4 core machines): dependencies
  without cache 70 min, with cache 2 min. Here TFEL (the longest,
  compiled from the sources of the repository), MGIS, the examples and mm-opera-hpc
  are added, and for CUDA an MFEM compiled with CUDA (several hours the first
  time, see `deps/spack/README.md`).
- The spack environment mmm-gpu is too big for the CI (NVIDIA HPC SDK,
  LLVM, Kokkos, PETSc, Ginkgo, SUNDIALS...). The CI has its own reduced
  environments.

## Proposal for the CI (9f)

Three levels, each in its own job:

1. CPU, at each push: spack environment `ci-cpu` (OpenMPI, hypre,
   MFEM with MUMPS and SuiteSparse), and `ci-cpu-float` for single
   precision, TFEL from the sources of the repository (ccache),
   MGIS applied from `patches/mgis`, mfem-mgis, examples, mm-opera-hpc;
   CPU tests of each. Also checks that the patch series
   apply (`scripts/verify-patches.sh`).
2. CUDA, compilation only, every night: environment `ci-cuda` (CUDA from
   spack, MFEM with CUDA, without GPU), kernels of mfem-mgis and `cuda`
   backends of TFEL compiled, tests on CPU.
3. HIP, compilation only, every night: ROCm through Python packages and
   `scripts/build-hip-amd.sh` (MFEM HIP, mfem-mgis, examples), tests without
   GPU.

Cache: spack dependencies on GHCR as for mfem-mgis (composite
actions, single index, monthly purge), plus ccache for TFEL. The
lessons of the CI of mfem-mgis are reused: pinned versions of
spack-packages (otherwise an update of autoconf invalidates the whole cache),
no `continue-on-error` on the push to the cache, read-only tokens
for the pull requests coming from forks.

The GPU tests stay on the development machine
(`scripts/mfem-mgis-gpu-suite.sh`), unless decided otherwise in E3.

## Proposal for the site (9g)

A single GitHub Pages site for the repository (one site per repository), built and
pushed by the CI:

- home: what the prototype is, links to the three parts;
- TFEL, MGIS, mfem-mgis: their documentations built from the
  sources of the prototype (web pages of TFEL and of MGIS, Sphinx and Doxygen of
  mfem-mgis);
- a tab "What's new with this prototype": differences with the
  official versions, additions, choices and justifications, fixes
  found, pitfalls encountered, and the log as the last sub-tab. These
  pages are drawn from `README.md`, `docs/` and `JOURNAL.md`, in English
  like the documentations of the three libraries (the documents of the
  prototype were translated into English on 14/10/2026).

## Decisions to take

- E1. Visibility. (a) Public repository: free minutes and packages, Pages
  possible, all the content published. (b) Private repository with a paid account
  (Pages then requires a public site anyway). (c) Private repository, CI on
  one's own machine (self-hosted runner), no site. Recommended: (a)
  if publication is wanted; otherwise (c) for the CI, without site.
- E2. Levels of the CI: the three, or only the CPU level at the start.
  Recommended: CPU first, then CUDA and HIP at night.
- E3. GPU tests through a self-hosted runner on the development
  machine (sharing of the GPU with the other work, security: never
  for pull requests coming from forks). Recommended: no for
  now.
- E4. Language of the site: English (like the documentations of the
  libraries), the repository keeping its documents in French (they were
  then translated into English on 14/10/2026). Recommended.

## Decisions taken (13/10/2026)

E1 (a): public repository, after an audit of the whole history (done on
13/10, nothing to rewrite); the change of visibility remains to be done by the
owner of the repository. E2: CPU first, then CUDA and HIP at night. E3:
no runner on the development machine. E4: site in English.

## Draft of the CPU level (14/10/2026)

`.github/workflows/cpu.yml`, `deps/spack/ci-cpu/spack.yaml` and
`deps/spack/ci-cpu-float/spack.yaml`, not yet run. Manual
trigger only as long as the repository is private (machine minutes);
pushes, pull requests and nightly run to be added once the
repository is public.

- Job `patches`: each series of `patches/` is applied on the upstream
  repository at its base, and the resulting tree is compared with the sources of the root repository
  (`diff -r`). `scripts/verify-patches.sh` does the same check
  locally, on the nested repositories, which do not exist on the CI. Simulated
  locally for the five series: all apply and give back the tree of
  the root repository.
- Job `cpu`: dependencies through spack (packages pinned to the same commit as
  the environments of the prototype), TFEL from the sources of the root repository (host
  backends `stlpar` and `openmp` only, without TDLS), MGIS applied from
  `patches/mgis`, then mfem-mgis, the examples and mm-opera-hpc, with their
  tests on the CPU; ccache through the composite actions of mfem-mgis.
  OpenMPI 5.0.10: `OMPI_MCA_btl=self,tcp` against the collision of the
  job identifiers of concurrent `mpirun`, and oversubscription
  allowed (one test launches three processes).
- Job `cpu` simulated locally on 14/10, end to end, in an
  empty environment: same spack (1.2.2) and same commit of
  spack-packages, copy of `ci-cpu` with gcc 14, then TFEL, MGIS (49/49),
  mfem-mgis (208/208), examples (22/22) and mm-opera-hpc (6/6). The
  simulation found a defect, fixed: the behaviour library of the
  tests of mfem-mgis (target `BehaviourTest`) was not built, hence
  126 failures out of 208.
- Without a cache of the dependencies, each job of the matrix rebuilds its
  own dependencies, OpenMPI, hypre, MUMPS and MFEM (and SuiteSparse in
  double only: MFEM refuses it in float), then TFEL: several
  hours on a two core machine. The cache on GHCR (like
  mfem-mgis) comes next.
- Float (9b): the job `cpu` is a matrix over the precision (`double`,
  `float`). The environment `deps/spack/ci-cpu-float` is the twin of
  `mmm-float-cpu`, and MGIS receives `MGIS_REAL_TYPE` from the matrix. Launched
  by hand, the job only takes double by default (input
  `precisions`); both precisions will run at each push once
  the repository is public. Simulated locally in float: MGIS 49/49, mfem-mgis
  204/204, examples 21/21, mm-opera-hpc 6/6. The simulation found a
  defect of mm-opera-hpc in float: the fixed point on the lateral
  stresses of the first step of the cermet did not go below its tolerance
  of 1e4 Pa (stagnation around 1.1e4 Pa or oscillation between 1.5e4 and
  1.7e4 Pa, depending on the rounding of the whole computation). The correction of the
  lateral components of F, just below 1, rounds to zero below
  about 1.5e4 Pa of residual. The step was cut, and the axial stress
  of the cermet at 0.2 s differed by 2.7 % from double, an effect of the time step
  and not of the rounding (0.24 to 0.27 % from double at the same steps). The fixed
  point now takes at least twice its resolution, computed from
  F and the moduli (30,804 Pa here, negligible in double):
  gap to double of 0.10 to 0.16 % at 0.2 s, and the test checks the same
  instants in both precisions (`results/cermet-float-2026-10-14`).
- The behaviours generated by mfront in the tests are compiled with
  `tfel-config --oflags`, which contains `-march=native` even in a
  portable build: their results depend on the processor of the
  CI machine. On the cermet in float, the behaviours compiled by gcc 13
  and by gcc 14 give the same stresses within 1e-5; the stopping
  point of the fixed point depends on the rounding, which moves the axial stress
  by at most 0.08 %, far from the tolerance of 1 % of the test.

## Draft of the site (14/10/2026)

- Generator: Sphinx with MyST (Markdown), theme `sphinx_rtd_theme`,
  like the documentation of mfem-mgis. The pages of the prototype, in
  English, are in `site/`: home and tab "What's new with this
  prototype" (overview, differences with the official versions,
  additions, choices, fixed defects, pitfalls, milestones). They
  summarize the documents of the repository, which remain the reference.
- `scripts/build-site.sh` builds the complete site in
  `build/site/html`: these pages at the root, the web pages of TFEL and of
  MGIS (CMake target `website`, pandoc 3.9 and pandoc-crossref 0.3.24, which
  must be built with the same version of pandoc), and the
  Sphinx and Doxygen documentation of mfem-mgis. Nothing is compiled: TFEL
  and MGIS are only configured. The patched sources of MGIS must
  be in `mgis/` (the CI applies `patches/mgis`). Python packages:
  `site/requirements.txt`. Built locally without warning from the
  pages of the prototype (145 MB).
- `.github/workflows/pages.yml` installs the tools, applies the patches
  of MGIS, runs the script and publishes on GitHub Pages. Manual
  trigger only: Pages requires a public repository.

## Steps, once E1 to E4 are decided

1. Environments `ci-cpu` and `ci-cpu-float`, CPU job without cache, then
   with the GHCR cache.
2. Nightly CUDA and HIP jobs.
3. Build of the site (Sphinx, see "Draft of the site"), publication
   by the CI.
4. Documentation of the CI and of the site in `docs/`.
