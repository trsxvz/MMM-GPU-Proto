# Machine profiles

Each directory describes a machine on which the prototype is built. The
variable `PROTO_MACHINE` selects the profile (`workstation` by default):

- `workstation`: the development machine of the prototype (RTX 5080,
  NVIDIA HPC SDK 25.11, oneAPI, spack environment `mmm-gpu`). Its files
  hold the paths and the options which were written in `scripts/env.sh`
  and `build-tfel.sh` before the profiles;
- `jean-zay-h100`: Jean Zay (IDRIS), partition `gpu_p6` (H100), spack
  environment `mmm-light`, installed without access to the internet
  (`docs/install.md`);
- `jean-zay-a100`: Jean Zay, partition `gpu_p5` (A100), the same way.

## Files of a profile

- `machine.sh` (required), sourced by `scripts/env.sh`, which defines
  `P`, the root of the prototype. It defines:
  - `PROTO_CC`, `PROTO_CXX`, `PROTO_FC`: C, C++ and Fortran compilers;
  - `PROTO_SPACK`: spack executable; `PROTO_SPACK_ENVIRONMENT`: directory
    of the spack environment of the libraries (MFEM, hypre, MPI);
    `PROTO_SPACK_VIEW`: its view;
  - `PROTO_CUDA`: CUDA toolkit;
  - `PROTO_TFEL_BUILD`, `PROTO_TFEL_INSTALL`: build and installation
    directories of TFEL;
  - `PROTO_MGIS_OPTIONS` (array): options of CMake of MGIS;
  - optionally `PROTO_AMGX` (`ON` by default; `OFF` builds mfem-mgis
    without AmgX), `PROTO_MPI_PACKAGE` (`openmpi` by default: the spack
    package of MPI) and `PROTO_JOBS` (parallel jobs of the builds);
  - `PATH`, `LD_LIBRARY_PATH`, `CMAKE_PREFIX_PATH` and the variables of
    the toolchains (modules of the machine, for instance).
- `tfel.sh` (required by `build-tfel.sh`), sourced by `build-tfel.sh`. It
  defines `BUILD` and `INSTALL` (directories), `CONFIGURATIONS` (the
  `generic-parallel` configurations of TFEL) and `TFEL_OPTIONS` (array of
  options of CMake, compilers included).
- `spack-config.yaml` (machines using `deps/spack/mmm-light`): included
  by the environment. External packages (compilers, CUDA, MPI, CMake),
  `cuda_arch`, target of the processors, directories of the builds.
- `jobs/`: jobs of the batch scheduler of the machine.

## Adding a machine

1. Copy `jean-zay-h100` to `deps/machines/<name>`.
2. In `spack-config.yaml`, give the prefixes of the compilers, CUDA, MPI
   and CMake of the machine (`module show`, `which`), its target
   (`x86_64_v3`, `x86_64_v4`, ...) and `cuda_arch` (80 for A100, 90 for
   H100, 120 for Blackwell).
3. In `machine.sh`, load the modules of the machine and give the path of
   CUDA; in `tfel.sh`, choose the configurations of TFEL.
4. Adapt the jobs (partitions, accounts, resources).
5. Follow `docs/install.md` with `MACHINE=<name>`.
