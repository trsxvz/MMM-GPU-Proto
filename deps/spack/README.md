# Spack environment of the prototype

Content:

- `mmm-gpu/spack.yaml`: environment. NVIDIA HPC SDK 25.11 and its CUDA
  13.0.2 as externals; Kokkos 5.2.2, HIP 7.2.3, LLVM 20.1.8, AdaptiveCpp
  25.10.0, TBB 2023.1.0 (toolchains of TFEL); MFEM 4.10.0 with CUDA,
  OpenMP and libCEED, and the
  solver libraries of the study (hypre 3.2.0, Ginkgo 1.11.0 with
  MPI, PETSc 3.26.0 with Kokkos, SUNDIALS, AmgX, SuperLU_DIST, STRUMPACK,
  MUMPS, UMFPACK), cuDSS 0.7.1, OpenMPI 5, OpenBLAS. The choice of the
  solvers is explained in [docs/gpu-solvers.md](../../docs/gpu-solvers.md);
- `mmm-float-cpu/spack.yaml`, `mmm-double-cpu/spack.yaml` and
  `mmm-float-gpu/spack.yaml`: environments of the float precision (task
  9b, [docs/plan-precision.md](../../docs/plan-precision.md)). MFEM 4.10 in
  single precision, hypre in single precision and MUMPS `+float`, on the
  CPU (`mmm-float-cpu`) or with CUDA and Umpire (`mmm-float-gpu`, without
  AmgX nor cuDSS); `mmm-double-cpu` is the double twin of `mmm-float-cpu`,
  to compare the two precisions with the same environment (in `mmm-gpu`,
  hypre is built for CUDA). They share the spack instance, the packages
  and the externals of `mmm-gpu`, which is not modified. In float,
  hypre is built with its internal BLAS and LAPACK (`~lapack`): with an
  external LAPACK, hypre 3.2.0 in single precision calls the LAPACK
  routines of the double type on float arrays (ParaSails, Gaussian
  elimination 98, 99, 198 and 199 of BoomerAMG, Schwarz, FSAI, MGR, AIR,
  GSMG, LOBPCG, AME);
- `mmm-raja/spack.yaml`: RAJA 2026.07.0 for CUDA (sm_120) and OpenMP, for
  the `raja` backend of TFEL and MGIS (`deps/spack/install-mmm-raja.sh`,
  then `build-tfel.sh`). RAJA with OpenMP needs a camp built with OpenMP,
  while the camp of `mmm-gpu` (a dependency of Umpire) is built without
  it: adding RAJA to `mmm-gpu`, whose packages are unified, would build
  Umpire, hypre and MFEM again. The environment shares the spack instance,
  the externals and the configuration of `mmm-gpu`; only RAJA and its camp
  are built (a few seconds). RAJA 2025.12.2 and older do not support
  CUDA 13;
- `ci-cpu/spack.yaml` and `ci-cpu-float/spack.yaml`: environments of the
  CI on the CPU (`.github/workflows/ci.yml`,
  [docs/plan-ci-site.md](../../docs/plan-ci-site.md)), in double (MFEM 4.10
  with MUMPS and SuiteSparse, hypre, OpenMPI) and in single precision (as
  `mmm-float-cpu`, without SuiteSparse, refused by MFEM in float). Only
  the package repository of spack, at the same commit as the other
  environments, and the compilers of the CI machine. The job `hip` of the
  CI also uses `ci-cpu` (MPI, METIS, SuiteSparse, OpenBLAS and hypre for
  the MFEM compiled with HIP by `scripts/build-hip-amd.sh`);
- `ci-cuda/spack.yaml`: environment of the job `cuda` of the CI
  (compilation only), the specs of `ci-cpu` with hypre and MFEM compiled
  with CUDA for sm_120; CUDA 13.0 is installed by the workflow from the
  packages of NVIDIA and found by `spack external find cuda`. It uses the
  repository `mmm_gpu` (patch of MFEM for CUDA 13), and MFEM embeds the
  PTX code only (variant `cuda_ptx_only`): for `fem/lor/lor_batched.cpp`,
  ptxas needs more than 12 GB even with `-O0`, and cicc 8.8 GB and 34 min;
  two jobs reached 12.5 GB for cicc, hence hypre and MFEM are compiled
  with one job;
- `mmm-light/spack.yaml`: environment of the clusters (installation without
  access to the internet, [docs/install.md](../../docs/install.md)): only
  the libraries needed by mfem-mgis, hypre 3.2 with CUDA, Umpire and 64 bit
  global indices (`+mixedint`), MUMPS, MFEM 4.10 with CUDA, Umpire, MUMPS,
  SuiteSparse and gslib, and GNU patch (built first, for the compute nodes
  without `patch`). The compilers, CUDA, MPI, `cuda_arch` and the target
  come from the `spack-config.yaml` of the machine profile
  (`deps/machines/`), included by the environment;
- `repos/spack_repo/mmm_gpu`: complementary package repository, with
  priority over the one of spack:
  - `adaptivecpp`: AdaptiveCpp 25.10 (the `hipsycl` package of spack
    stops at 24.10, which does not compile with CUDA 13);
  - `amgx`: version `2.5.0.20260701`, commit of the main branch with two
    complementary fixes for CUDA 13, later than version 2.5.0 (C++17 with
    CUDA 13.1 and above, tuples of CUDA 13);
  - `openmpi`: internal hwloc without OpenCL (see [MPI](#mpi));
  - `mfem`: a partial build resumes without recompiling everything
    (`make config` rewrites `config/config.mk`, on which all the objects
    depend: its date is kept when its content does not change); patch
    `mfem-4.10-convection-ea-cuda13.patch` (CUDA 13 only): with
    CUDA 13 (sm_120), cicc does not
    finish `fem/integ/bilininteg_convection_ea.cpp` (more than 2 h 30 and
    17 GB, then killed for lack of memory); the 3D kernels of the element
    assembly of the convection are no longer specialised beyond order 3,
    the generic kernel of the same file takes the higher orders; variant
    `cuda_ptx_only` (PTX code only, ptxas not run at compilation) for the
    compilation checks of the CI;
  - `petsc`: complete copy of the recipe of spack-packages (commit
    3b7fdffc24), whose later fixes are not taken over. Only modification:
    the `cuda_lambda` variant of Kokkos is required only up to Kokkos 4
    (it no longer exists in Kokkos 5, where CUDA lambdas are always
    enabled);
- `mmm-solvers/spack.yaml`: environment of the phase 2 of the solver
  libraries ([docs/plan-solver-libraries.md](../../docs/plan-solver-libraries.md)):
  the specs of `mmm-gpu`, with MFEM built with the sparse direct solvers
  of oneMKL (PARDISO, Cluster Sparse Solver) and MAGMA, PETSc with HPDDM
  and PaStiX, MAGMA 2.10 and PaStiX 6.4. Only MFEM, PETSc, MUMPS, MAGMA,
  PaStiX, Scotch and oneMKL are built, the other packages of `mmm-gpu`
  being reused. PETSc is built without MUMPS (PETSc requires Fortran for
  it, and MUMPS is directly available in mfem-mgis). It uses the
  repository `mmm_solvers` (below) before `mmm_gpu`;
- `repos/spack_repo/mmm_solvers`: package repository of the phase 2,
  with priority over `mmm_gpu`, whose packages are left unchanged (their
  hashes, and so `mmm-gpu` and the CI, do not change):
  - `mfem`: derives from the `mfem` of `mmm_gpu`, with the variants
    `mkl-pardiso`, `mkl-cpardiso` (with the BLACS of oneMKL for the MPI
    of the environment) and `magma`, which the package of spack does not
    have;
  - `petsc`: derives from the `petsc` of `mmm_gpu`, with the variant
    `pastix` (PaStiX 6.4, the version supported by PETSc 3.26);
  - `intel-oneapi-mkl-solvers`: oneMKL installed with the installer of
    `intel-oneapi-mkl`, without providing BLAS, LAPACK and ScaLAPACK: in a
    unified environment, a package which can provide a virtual package
    present in the graph is its provider for the whole graph, and oneMKL
    replaced OpenBLAS for all the packages (a requirement on the virtual
    package does not help: spack applies it to every provider in the
    graph);
- `install-mmm-solvers.sh`: installation of `mmm-solvers`, as
  `install-mmm-gpu.sh`, on the cores 2-7 and 10-15 with the lowest
  priority;
- `mmm-trilinos/spack.yaml`: environment of the phase 4 of the solver
  libraries: the specs of `mmm-solvers`, plus Trilinos 17.2.1 (Tpetra,
  Belos, Ifpack2, MueLu, Amesos2 with Tacho and Basker, FROSch), with
  CUDA. Trilinos 17.2.1 requires Kokkos and Kokkos Kernels 5.2.1, where
  the other environments have 5.2.2: the packages which depend on Kokkos
  (PETSc, and MFEM through PETSc) are built again, `mmm-solvers` is left
  unchanged. Trilinos with Kokkos and CUDA requires hwloc with CUDA:
  LLVM, AdaptiveCpp and TBB, which would be built again with it and are
  not used by the solvers, are left out;
- `repos/spack_repo/mmm_trilinos`: package repository of the phase 4,
  with priority over `mmm_solvers`: `trilinos` derives from the builtin
  package, with the variant `tacho` (sparse direct solver Tacho in
  Amesos2: the variant `shylu` only enables ShyLU_DD, and the optional
  dependencies of Amesos2 are not enabled; the tests of Tacho require
  GoogleTest, which TriBITS requires as soon as Tacho is enabled, and are
  disabled), and the variants `fastilu` (ShyLU_NodeFastILU in Ifpack2)
  and `teko-tpetra` (Teko on Tpetra, without Epetra);
- `mmm-trilinos-teko/spack.yaml`: the environment `mmm-trilinos`, with
  Trilinos also built with the fast incomplete factorizations of
  ShyLU_Node in Ifpack2 (variant `fastilu`: `FAST_ILU`, `FAST_IC`,
  `FAST_ILDL`) and the block preconditioners of Teko on Tpetra (variant
  `teko-tpetra`: Teko, Thyra, ThyraTpetraAdapters, Stratimikos). Only
  Trilinos is built again; the other packages, MFEM included, keep the
  hashes of `mmm-trilinos`, which is left unchanged. The variant `teko`
  of the builtin package requires the packages of Epetra (`epetraext`,
  `ml`), removed from Trilinos 17, hence the variant `teko-tpetra` of
  `mmm_trilinos`;
- `install-mmm-trilinos-teko.sh`: installation of Trilinos in
  `mmm-trilinos-teko`, 4 jobs, under `scripts/memory-governor.sh`;
- `install-mmm-trilinos.sh`: installation of `mmm-trilinos`, as
  `install-mmm-solvers.sh`: the dependencies of MFEM and Trilinos first,
  then MFEM alone with two jobs (the file `fem/lor/lor_batched.cpp` needs
  about 22 GB for ptxas with sm_120, which leaves no room for other
  compilations on a machine of 32 GB). On the machine of the prototype,
  shared with other work, this file was compiled in the stage of MFEM to
  PTX code only, before the rest of MFEM (see below);
- `make-nvhpc-cuda-prefix.sh`: CUDA prefix with the layout of the
  standalone toolkit, built from the SDK (whose mathematical libraries are
  in a separate directory);
- `install-mmm-gpu.sh`: installation in two stages, suited to a machine
  with 30 GB of memory (see below).

## Recreate the environment

```bash
git clone --depth 1 --branch v1.2.2 https://github.com/spack/spack.git spack-mmm
git clone https://github.com/spack/spack-packages.git spack-packages-mmm
git -C spack-packages-mmm checkout 3b7fdffc24
deps/spack/make-nvhpc-cuda-prefix.sh <sdk>/Linux_x86_64/25.11 13.0 <prefix>
```

In `mmm-gpu/spack.yaml`, adapt:

- the path of the `builtin` repository (clone of spack-packages);
- the prefixes of `nvhpc` (root of the SDK) and of `cuda` (prefix built
  above);
- the external `gcc` compiler;
- `cuda_arch`: 80 for A100, 120 for Blackwell (RTX 5080). Kokkos accepts
  only one architecture per build.

Then:

```bash
export SPACK_DISABLE_LOCAL_CONFIG=1
spack-mmm/bin/spack -e deps/spack/mmm-gpu concretize -f
deps/spack/install-mmm-gpu.sh spack-mmm/bin/spack
```

`SPACK_DISABLE_LOCAL_CONFIG` isolates the environment from the
configuration of the user (`~/.spack`). The environment also has its own
cache (`misc_cache`): the cache of the user is shared with the other spack
installations of the machine, whose `builtin` repository may be at another
commit, which breaks the index of the patches.

`install-mmm-gpu.sh` first installs everything except MFEM (two packages
at a time, 4 jobs each), then MFEM alone with 2 jobs. Each phase runs
under `scripts/memory-governor.sh`: when the available memory is too low,
it suspends all the processes of the build except the largest one (and
its ancestors), which can thus finish and free its memory; everything
resumes afterwards. `ptxas` takes about 22 GB for the kernels of
`fem/lor/lor_batched.cpp` (sm_120), the other files up to 8 GB. Each
phase is restarted up to three times in case of failure (lack of memory),
MFEM with a single job; an interrupted build resumes where it stopped
(`--dont-restage` and `mfem` recipe of `mmm_gpu`).

When 22 GB can't be spared, `fem/lor/lor_batched.cpp` alone can be
compiled to PTX code only, which needs no ptxas (cicc takes about 6 GB
and 18 minutes): its kernels (batched assembly of low order refined
operators, not used by mfem-mgis) are then compiled by the driver at their
first use. This was done for the MFEM of `mmm-trilinos` (16/10/2026):

```bash
E=deps/spack/mmm-trilinos
spack -e $E install --only package -u configure /hb6luie  # stage and configure
cd $(spack -e $E location -s /hb6luie)/spack-src
# the command printed by `make -n fem/lor/lor_batched.o`, with
# -arch=compute_120 instead of -arch=sm_120
spack -e $E build-env /hb6luie -- bash -c "<command>"
# spack configures MFEM again: the object is dated in the future so that it
# is not compiled again
touch -d tomorrow fem/lor/lor_batched.o
spack -e $E install --dont-restage --only package -j 2 /hb6luie
```

## CUDA prefix

In the prefix, `bin` is a directory of links, not a link: nvcc looks for
its headers in `bin/../targets`. With a link, `bin/..` designates the CUDA
of the SDK, without the headers of the mathematical libraries
(`cublas_v2.h` not found for the packages that do not pass the directory
of the headers explicitly, such as Kokkos Kernels).

## MPI

OpenMPI is compiled with CUDA (GPU aware MPI), so that the code stays
generic, even if the machine has only one GPU. It uses its internal hwloc,
built without its OpenCL plugin: with CUDA, this hwloc finds the OpenCL
library of the toolkit and fails at link time (`cannot
determine absolute directory name of 'no/opencl/lib/'`). The `openmpi`
recipe of `mmm_gpu` adds `--disable-opencl`. The hwloc of spack with CUDA
is not used: shared with LLVM and TBB, it would rebuild them.

## State of the installation of 09/10/2026

The current installation of MFEM (hash `zxpdstv`) was finished by hand to
meet the deadlines; its content is the one given by the recipes, but not
yet by a direct build:

- `fem/integ/bilininteg_convection_ea.cpp` was compiled from the patched
  version (same content as `mfem-4.10-convection-ea-cuda13.patch`), in the
  build directory;
- the static library of gslib (`libgs.a`, in the prefix of gslib) was
  recompiled with `-fPIC` (same sources, same options): not relocatable,
  it could not go into `libmfem.so`. The original is kept outside the
  repository.

To do: concretize again and rebuild MFEM (about 3 h), for an installation
entirely produced by the recipes. gslib `+shared` is already imposed in
`mmm-gpu/spack.yaml` (see [Patched AmgX](#patched-amgx-10102026)). The
condition of the convection patch (`^cuda@13:`, added after the
installation) also changes the hash of MFEM: the same new concretization
takes it into account.

## Patched AmgX (10/10/2026)

With CUDA 13 (sm_120), the classical multigrid of the AmgX of commit
91a8413e fails ("invalid device function"): the `amgx` recipe of
`mmm_gpu` applies `amgx-cuda13-sequence-kernel.patch` (reproducer
`reproducers/amgx-classical-cuda13`, upstream PR NVIDIA/AMGX #386).
The patched AmgX (hash `xupvzsc`) is installed, but the environment keeps
its previous spack.lock, whose MFEM contains the unpatched AmgX: the
rebuild of MFEM waits until the fixes (AmgX #386, MFEM #5604) are merged
upstream then in spack-packages. Until then, mfem-mgis uses the patched
AmgX through `AMGX_DIR` (set by `scripts/env.sh`), placed before MFEM at
link time. `mmm-gpu/spack.yaml` also imposes gslib `+shared`, taken into
account at the same rebuild.

## Constraints between packages

- PETSc with `+openmp`: the headers of a Kokkos compiled with OpenMP
  require OpenMP in the compiler (`#error "You enabled Kokkos OpenMP
  support without enabling OpenMP in the compiler!"`).
- MFEM with `cxxstd=20`: Umpire 2026 (through camp) uses the concepts of
  C++20; MFEM compiles in C++17 by default.
- shared kokkos-kernels requires Kokkos `+pic`: Kokkos was rebuilt with
  this variant (same version, same options otherwise).
- hypre with 32 bit integers (`~int64`): the spack recipe of mfem-mgis
  requests `hypre+int64` by default (variant `int64` of mfem-mgis), which
  hypre refuses with CUDA ("Use +mixedint for 64-bit integer support for
  GPUs!"). An installation of mfem-mgis by its recipe with CUDA therefore
  requires `mfem-mgis~int64`. The 32 bit integers limit the total number
  of unknowns to about 2 billion, which does not hinder the prototype.
- hypre is built without `+gpu-aware-mpi` (default value): its MPI
  exchanges go through the host, even with CUDA aware OpenMPI. Variant to
  compare for the measurements on GPU with several processes.

## AMD GPU (HIP)

The MFEM compiled with HIP for an AMD GPU (gfx90a) does not come from this
environment: `scripts/build-hip-amd.sh` installs ROCm by its Python
packages, without administrator rights, and builds hypre and MFEM in
`deps/hip-amd/` (not tracked), with MPI, METIS, SuiteSparse and OpenBLAS
from `mmm-gpu`. See the header of the script.

## Drivers

CUDA 13.0 requires a driver 580 or newer. With an older driver, use the
CUDA 12.9 of the same SDK (`make-nvhpc-cuda-prefix.sh
... 12.9 ...` and external `cuda@12.9.1`).
