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
- `ci-cpu/spack.yaml` and `ci-cpu-float/spack.yaml`: environments of the
  CI on the CPU (`.github/workflows/cpu.yml`,
  [docs/plan-ci-site.md](../../docs/plan-ci-site.md)), in double (MFEM 4.10
  with MUMPS and SuiteSparse, hypre, OpenMPI) and in single precision (as
  `mmm-float-cpu`, without SuiteSparse, refused by MFEM in float). Only
  the package repository of spack, at the same commit as the other
  environments, and the compilers of the CI machine;
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
    the generic kernel of the same file takes the higher orders;
  - `petsc`: complete copy of the recipe of spack-packages (commit
    3b7fdffc24), whose later fixes are not taken over. Only modification:
    the `cuda_lambda` variant of Kokkos is required only up to Kokkos 4
    (it no longer exists in Kokkos 5, where CUDA lambdas are always
    enabled);
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
