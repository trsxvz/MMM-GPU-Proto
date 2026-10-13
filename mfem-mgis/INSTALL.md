# Installation guide

The installation guide is available at
<https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html>.

In short, `mfem-mgis` is packaged in [Spack](https://spack.io/), which also
installs all its dependencies. These commands set up Spack:

```sh
git clone --depth=2 --branch=v1.2.2 https://github.com/spack/spack.git
source spack/share/spack/setup-env.sh
spack repo update builtin --branch develop
```

The last command selects the `develop` branch of the Spack packages, because
their releases do not provide `mfem-mgis` yet.

The `@master` version of `mfem-mgis`, its development version, is currently
recommended, as it includes many fixes:

```sh
spack install mfem-mgis@master
```

The last release (1.0.4) is installed by:

```sh
spack install mfem-mgis
```

These commands install the upstream versions of `mfem-mgis`, without the
features of this prototype (computations by batches and by kernels,
execution settings, computations on GPUs). The prototype is built from its
sources, with versions of MGIS and TFEL which are not the `master` ones:
see the section *Building the prototype* of the installation guide
([docs/source/installation_guide/installation_guide.rst](docs/source/installation_guide/installation_guide.rst)),
which also describes the builds for GPUs (MFEM compiled with CUDA or HIP,
AmgX, cuDSS, Kokkos).
