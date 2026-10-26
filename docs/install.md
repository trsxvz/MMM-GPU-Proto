# Installation on a cluster without access to the internet

This document describes how the prototype is installed on a cluster
whose nodes have no access to the internet, with spack for the libraries
and the scripts of the prototype for TFEL, MGIS and mfem-mgis. The
cluster never uses the network: everything is downloaded on another
machine and copied. The same procedure works on a cluster with access to
the internet.

Jean Zay (IDRIS) is the first machine installed this way (profiles
`jean-zay-h100` and `jean-zay-a100`, see the last sections).

## Principle

Two machines take part:

- a **connected machine** (a workstation, a laptop) with access to the
  internet. Nothing is compiled on it: it downloads and packs. Its power
  does not matter; it needs about 5 GB of disk and git, bash and Python 3;
- the **target machine** (the cluster), which needs bash, Python 3 (for
  spack), git, tar, its compilers, CUDA and MPI.

Only files go from one to the other (`scp`, `rsync`, a disk):

| Step | Machine | Command | Result |
|---|---|---|---|
| 1 | connected | `scripts/offline/bundle.sh base WORK` | `mmm-base.tar` (about 2 GB) |
| 2 | target | `tar -xf mmm-base.tar`, `./deploy.sh setup` | spack ready without network, prototype and MGIS |
| 3 | target | `deploy.sh concretize MACHINE` | `spack.lock` (a few hundred kB) |
| 4 | connected | `scripts/offline/bundle.sh mirror WORK MACHINE spack.lock` | `mmm-mirror-MACHINE.tar` (sources) |
| 5 | target | `deploy.sh install MACHINE mmm-mirror-MACHINE.tar` | libraries of the spack environment |
| 6 | target | `deploy.sh build MACHINE` | TFEL, MGIS, mfem-mgis, examples, mm-opera-hpc, cases |

The concretization (step 3) is made on the target machine, which knows
its compilers, its operating system and its external packages: the
sources of step 4 are then exactly the ones the installation needs. The
copy of `spack.lock` back to the connected machine is the price of this
exactness.

What `mmm-base.tar` contains:

- spack (`scripts/offline/versions.sh`: 1.2.2) and its packages
  (spack-packages at the commit of the environments of the prototype),
  without their history;
- the bootstrap mirror of spack (`spack bootstrap mirror
  --binary-packages`): clingo, the solver of the concretizer, and its
  dependencies, as binaries and sources. Without it, spack can't
  concretize offline;
- the prototype: git bundle of the current commit of the connected
  machine (changes which are not committed are not packed);
- MGIS at the base of the patches of the prototype (git bundle): the
  prototype only holds the patches of MGIS (`patches/mgis`);
- `deploy.sh`, the script of the target machine.

## Deployment directory

On the target machine, everything lives in one directory ROOT:

```
ROOT/
  deploy.sh  MANIFEST  *.bundle  *.tar   files of the archives
  spack/  spack-packages/                instance of spack and its packages
  bootstrap/  bootstrap-store/           bootstrap of spack (mirror, clingo)
  mirror-MACHINE/                        sources of the packages
  misc-cache/  user-cache/               caches of spack
  MMM-GPU-Proto/                         the prototype, with mgis/ (branch gpu-proto)
```

A deployment directory serves one machine profile (recorded in
`ROOT/MACHINE`, `deploy.sh` refuses another one): its spack environment,
its `spack.lock` and the builds of the prototype are those of this
profile. Two profiles of the same cluster (for instance the H100 and A100
partitions of Jean Zay) are deployed in two directories, from the same
`mmm-base.tar`.

The configuration of spack is written in the site scope of this instance
(`ROOT/spack/etc/spack`), and the configuration of the user (`~/.spack`)
is ignored (`SPACK_DISABLE_LOCAL_CONFIG`): the deployment depends on
nothing else. `deploy.sh setup` sets:

- the builtin repository of spack to `ROOT/spack-packages` (by default,
  spack 1.2 clones it from GitHub);
- the bootstrap from `ROOT/bootstrap`, the other bootstrap sources being
  disabled, and checks it (`spack bootstrap now`).

## Step by step

On the connected machine, from a clone of the prototype (WORK is any
directory, reused from one call to the next):

```bash
scripts/offline/bundle.sh base ~/mmm-offline
scp ~/mmm-offline/mmm-base.tar cluster:ROOT/
```

On the target machine:

```bash
cd ROOT && tar -xf mmm-base.tar && ./deploy.sh setup
MMM-GPU-Proto/scripts/offline/deploy.sh concretize MACHINE
# ROOT/MMM-GPU-Proto/deps/spack/mmm-light/spack.lock goes to the connected machine
```

On the connected machine:

```bash
scp cluster:ROOT/MMM-GPU-Proto/deps/spack/mmm-light/spack.lock ~/mmm-offline/spack.lock-MACHINE
scripts/offline/bundle.sh mirror ~/mmm-offline MACHINE ~/mmm-offline/spack.lock-MACHINE
scp ~/mmm-offline/mmm-mirror-MACHINE.tar cluster:ROOT/
```

On the target machine, preferably in a job (see the jobs of the profile):

```bash
MMM-GPU-Proto/scripts/offline/deploy.sh install MACHINE ROOT/mmm-mirror-MACHINE.tar
MMM-GPU-Proto/scripts/offline/deploy.sh build MACHINE
```

`deploy.sh build` runs `build-tfel.sh`, `build-mgis.sh` and
`build-mfem-mgis.sh` with the profile `MACHINE`: builds, tests (target
`check`) and installations in `ROOT/MMM-GPU-Proto/build` and
`ROOT/MMM-GPU-Proto/install`. The arguments after MACHINE go to
`build-mfem-mgis.sh` (`mfem-mgis`, `examples`, `opera`, `cases`).

Afterwards, `source ROOT/MMM-GPU-Proto/scripts/env.sh` with
`PROTO_MACHINE=MACHINE` gives the environment of the prototype (compilers,
CUDA, MPI, TFEL).

### Update of the prototype

A new commit of the prototype does not require the 2 GB archive again:

```bash
scripts/offline/bundle.sh update ~/mmm-offline      # connected machine: mmm-update.tar
cd ROOT && tar -xf mmm-update.tar && ./deploy.sh setup   # target machine
```

`deploy.sh setup` moves the prototype to the packed commit (it refuses if
the clone of the target has changes which are not committed) and
recreates `mgis/` from its base and the patches of the new commit. If
the spack environment of the machine changed, the steps 3 to 5 are run
again; spack only builds the packages which changed.

### Binary cache

With `PROTO_BUILDCACHE` set to a directory, `deploy.sh install` copies the
installed packages in a binary cache of spack (a few files per package),
with the `spack.lock` used. On a filesystem with few inodes, or whose
unused files are deleted after some time, the installation is kept this
way, and `deploy.sh restore MACHINE PROTO_BUILDCACHE` reinstalls it in a
new deployment without compiling.

## Spack environment `mmm-light`

`deps/spack/mmm-light/spack.yaml` holds the libraries needed by
mfem-mgis on CPU and GPU, nothing more: hypre 3.2 (CUDA, Umpire), MUMPS,
MFEM 4.10 (CUDA, Umpire, MUMPS, SuiteSparse, gslib), with the recipes of
the repository `mmm_gpu` (MFEM patched for CUDA 13). The solver libraries
of the study of the workstation (PETSc, Ginkgo, SUNDIALS, SuperLU_DIST,
STRUMPACK, AmgX, cuDSS, Trilinos) and the toolchains of the other
backends of TFEL (Kokkos, RAJA, LLVM, AdaptiveCpp, HIP) are not part of it.

hypre is built with `+mixedint`: 64 bit global indices (`HYPRE_BigInt`)
and 32 bit local indices (`HYPRE_Int`). The total number of unknowns can
exceed 2^31 with CUDA, where hypre refuses `+int64`; the number of
nonzero entries of the local matrices of each process stays below 2^31
(about 26 million unknowns per process in 3D elasticity with linear
elements). MFEM detects the type of the indices of hypre. The coarse mesh
is partitioned by METIS with 32 bit indices: mfem-mgis refines it in
sequential only until each process has 8 elements, the following
refinements being parallel.

The file holds nothing specific to a machine: the compilers, CUDA, MPI,
the architecture of the GPUs and the target of the processors come from
the profile of the machine (`include` of its `spack-config.yaml`).

## Machine profiles

A machine is described by a directory `deps/machines/MACHINE` (see
`deps/machines/README.md`):

| File | Role |
|---|---|
| `machine.sh` | sourced by `scripts/env.sh`: modules, compilers (`PROTO_CC`, `PROTO_CXX`, `PROTO_FC`), CUDA (`PROTO_CUDA`), spack (`PROTO_SPACK`, `PROTO_SPACK_ENVIRONMENT`), options of MGIS and mfem-mgis |
| `tfel.sh` | sourced by `build-tfel.sh`: configurations `generic-parallel` and options of CMake of TFEL |
| `spack-config.yaml` | included by `mmm-light`: external packages, `cuda_arch`, target |
| `jobs/` | jobs of the batch scheduler of the machine |

`PROTO_MACHINE` selects the profile; `workstation`, the development
machine, is the default value, so that its scripts and paths are the
ones of before.

## Jean Zay (profile `jean-zay-h100`)

- **Nodes**: partition `gpu_p6`, 4 H100 80 GB (`sm_90`) and 2 Xeon
  Platinum 8468 (Sapphire Rapids) per node, RHEL 9.6, NVIDIA driver
  595.91 (CUDA 13). The login nodes are Cascade Lake, with the same
  system (see the compilers below). The login nodes reach some sites through a proxy of IDRIS, the
  compute nodes none: the procedure above does not use this proxy.
- **External packages** (modules after `module load arch/h100`):
  gcc 14.2.0, CUDA 13.0.3, Open MPI 4.1.8 with CUDA, UCX and the PMI of
  Slurm, CMake 3.31.4. The libraries are built for `x86_64_v4` (AVX-512).
- **Storage**: the WORK and STORE spaces have quotas per project, shared
  by its members; the SCRATCH space has very large quotas, and its files
  not read nor modified for 30 days are deleted. The deployment goes in
  `$SCRATCH/mmm`, and the binary cache (`PROTO_BUILDCACHE`, a few hundred
  files) in `$WORK`: after a deletion, `deploy.sh restore` reinstalls the
  libraries in a few minutes.
- **Compilers of the H100 partition**: the modules loaded after
  `module load arch/h100` are built for the Sapphire Rapids processors of
  this partition. Their gcc stops with "internal compiler error: Illegal
  instruction" on the login nodes and on the other partitions (Cascade
  Lake). Everything that compiles therefore runs on the partition
  `compil_h100` (compilation node of the H100 partition, no hours
  consumed) or `gpu_p6`. The concretization, which compiles nothing, runs
  on a login node.
- **Installation of the libraries**: job `jobs/install.slurm` on the
  partition `compil_h100`, submitted from ROOT:

  ```bash
  cd $SCRATCH/mmm && tar -xf mmm-mirror-jean-zay-h100.tar
  PROTO_BUILDCACHE=$WORK/mmm-buildcache sbatch -A <project>@h100 \
    MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/install.slurm
  ```

- **Build of the prototype**: job `jobs/build.slurm` on a node of
  `gpu_p6` (one GPU). MFront compiles the code it generates with
  `-arch=native` (nvcc) and `-march=native` (gcc): TFEL, MGIS and
  mfem-mgis are built on the nodes where they run, whose GPU is then the
  one of the computations.

  ```bash
  cd $SCRATCH/mmm && sbatch -A <project>@h100 \
    MMM-GPU-Proto/deps/machines/jean-zay-h100/jobs/build.slurm
  ```

- **TFEL**: configurations `cuda-nvcc` and `openmp-host-gcc`; the tests of
  TFEL are built only with `PROTO_TFEL_TESTING=ON`.
- **Tests on the GPU**: job `jobs/tests-gpu.slurm` (one H100): the suite of
  mfem-mgis on the GPU (`scripts/mfem-mgis-gpu-suite.sh`, modes kernels
  and default), the tests selecting a device, and the tests of the
  examples, of mm-opera-hpc and of the cases.
- **Measurements on several GPUs**: job `jobs/measure-multi-gpu.slurm`
  (one whole node, 4 H100 and 96 cores, QoS dev, one hour at most):
  `cases/measure-multi-gpu.sh`, strong scaling of the beam on 1, 2 and 4
  GPUs (with and without `--gpu-aware-mpi`) and on 24, 48 and 96 processes
  on the CPU. The jobs of the measurement and the tests run `mpirun`
  (Open MPI of IDRIS) in the allocation of the job; each process takes the
  GPU of its rank on the node.

## Jean Zay, A100 partition (profile `jean-zay-a100`)

Same organisation as `jean-zay-h100`, with the modules of `arch/a100`:

- **Nodes**: partition `gpu_p5`, 8 A100 80 GB (`sm_80`) and 2 AMD EPYC
  with 64 cores per node (no AVX-512: target `x86_64_v3`).
- **External packages**: gcc 14.2.0, CUDA 13.0.3, Open MPI 4.1.8 with CUDA,
  PSM2 (Omni-Path) and the PMIx of Slurm, CMake 3.31.4. These compilers
  also run on the login nodes and on the CPU partition: the libraries are
  installed by `jobs/install.slurm` on `cpu_p1` (QoS dev, `<project>@cpu`).
- **Build, tests, measurements**: jobs `jobs/build.slurm`,
  `jobs/tests-gpu.slurm` (one A100) and `jobs/measure-multi-gpu.slurm`
  (one whole node: 1, 2, 4 and 8 GPUs, 32, 64 and 128 processes on the
  CPU), with `<project>@a100`.
- **Deployment**: in its own directory, for instance `$SCRATCH/mmm-a100`.
