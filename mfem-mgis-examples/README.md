# mfem-mgis-examples

This repository provides user oriented material for the external `mfem-mgis` [repository](https://github.com/thelfer/mfem-mgis).
You will find here meshes and a collection of use cases for the `mfem-mgis` library.

## How to install

### Installing mfem-mgis

The examples of this branch follow the prototype branch of mfem-mgis, which provides the CMake
function `mfem_mgis_add_behaviours_library`, the option `--execution` and the function
`mfem_mgis::getBehavioursLibrary` used by all the examples. This branch of mfem-mgis must be built
from its sources, as described by its installation guide. It requires a version of MGIS providing
`MGIS/FieldOperations.hxx`. The Spack package of mfem-mgis (`spack install mfem-mgis@master`)
installs a version without these functions, with which the configuration of the examples fails.
Spack may still provide the other dependencies (MFEM, TFEL...).

The examples build their behaviour libraries with MFront: when CMake configures the build,
`mfront` must be found in the `PATH` (or given by the CMake variable `MFRONT`).

Then build and run the examples, `MFEMMGIS_DIR` giving the installation of mfem-mgis:

```
cd mfem-mgis-examples
cmake -B build -DCMAKE_BUILD_TYPE=Release -DMFEMMGIS_DIR=<prefix of mfem-mgis>/share/mfem-mgis/cmake
cmake --build build -j 4
ctest --test-dir build
```

The tests compute the complete simulations, except in the Debug and Coverage builds, where they
only compute their beginning to be faster. The `MFEM_MGIS_EXAMPLES_TEST_MODE` option, set to `full`
or `restricted` (`auto` by default), selects the test mode explicitly:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DMFEM_MGIS_EXAMPLES_TEST_MODE=restricted
```

The `MFEM_MGIS_EXAMPLES_GPU_TESTS` option (off by default) adds the tests of the examples on a GPU
(option `--execution gpu`, or `mixed` when a phase stays on the CPU), when MFEM is compiled with
CUDA or HIP and the behaviour libraries of the example are built for the GPU (TFEL providing the
configuration file `config-cuda-nvcc.json`, resp. `config-hip-hipcc.json`, of MFront; CMake prints
why they are not built). A GPU is then required to run them.

All the examples accept the options of mfem-mgis, in particular `--execution cpu|gpu|mixed`, which
places the computations on the CPU, on the GPU, or each phase (finite elements, behaviours, linear
solvers, post-processings) on one of them. Each example builds its behaviour libraries for the CPU
and, if MFEM uses a GPU, for the GPU (CMake function `mfem_mgis_add_behaviours_library` of
mfem-mgis), and uses the one placed by `--execution`, unless a library is given by its option `-l` (`-lU` and
`-lA` for ex8).

For an installation on a supercomputer without internet please follow the procedure described here for MFEM-MGIS installation: https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html#installation-guide-on-topaze-ccrt-of-mfem-mgis-examples
This procedure installs `mfem-mgis@master`: the prototype branch of mfem-mgis must then be built
from its sources in this environment.

## Test case description

| Name | Description | Directory
|--|--|--|
| TensileTest | Cyclic tension-compression of a unit cube made of a plastic material with linear isotropic hardening. | ex1 |
| Ssna303     | 2D (plane strain) tensile test on a notched beam with a finite-strain plastic behaviour, described in the [tutorial](https://thelfer.github.io/mfem-mgis/user_guide/tutorial.html) of mfem-mgis. | ex2 |
| TwoLayerCube | Periodic cube made of two elastic layers under an imposed macroscopic strain. | ex3 |
| Ssna303_3d  | 3D tensile test on a notched beam with a finite-strain plastic behaviour, solved with MUMPS, hypre or PETSc. | ex4 |
| Satoh       | Thermo-elastic plate in plane strain, clamped on its left and right boundaries and subjected to a parabolic temperature profile. | ex5 |
| Rve-elastic | Periodic Representative Volume Element (RVE) made of two materials with a Saint Venant-Kirchhoff hyperelastic behaviour under an imposed macroscopic deformation gradient: a two-layer cube, or spherical inclusions. | ex6 |
| Mox2        | Representative Volume Element (RVE) of a mixed oxide fuel: elastic inclusions in a viscoplastic matrix under an imposed macroscopic strain. More information in ex7/README.md | ex7 |
| RJH plate   | Thermomechanical simulation of one ring of a RJH fuel assembly. The fuel is U3Si2, the cladding and stiffeners are ALFENI. The heat transfer is non-linear with convective boundary conditions. The mechanics is in finite strain with thermal expansion, irradiation creep, swelling and plasticity. Both are strongly coupled at each time step. More information in ex8/README.md | ex8 |
