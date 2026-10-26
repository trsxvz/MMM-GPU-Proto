# MMM-GPU-Proto

A prototype of the chain MFront (TFEL), MGIS and mfem-mgis on CPU and GPU.
MFront generates, for each behaviour, model and material property, code
which processes a whole batch of integration points on several backends
(CUDA, HIP, Kokkos, SYCL, OpenMP, parallel algorithms of C++). MGIS loads
these libraries and integrates batches without copying the data.
mfem-mgis integrates the behaviours by batches on CPU and GPU, its data at
the integration points living in the memory of MFEM, with the finite
element computations made by kernels.

## Documentation of the libraries

Built from the sources of the prototype:

- <a href="tfel/index.html">TFEL/MFront</a>, including the
  `generic-parallel` interface of the behaviours and of the models;
- <a href="mgis/index.html">MGIS</a>;
- <a href="mfem-mgis/index.html">mfem-mgis</a> and its
  <a href="mfem-mgis/doxygen/index.html">API documentation</a> (Doxygen).

## What's new with this prototype

The differences with the official versions, the additions, the choices and
their justifications, the defects found, the pitfalls met and the
milestones of the work.

```{toctree}
:maxdepth: 2

whats-new/index
```

## Installation on a cluster

How the prototype is installed with spack on a cluster without access
to the internet (Jean Zay, H100 partition, for a start).

```{toctree}
:maxdepth: 2

install
```
