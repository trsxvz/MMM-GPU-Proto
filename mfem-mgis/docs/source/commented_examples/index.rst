==================
Commented Examples
==================

The examples come from two repositories.

mfem-mgis-examples
==================

These examples are in https://github.com/latug0/mfem-mgis-examples.

``ctest`` runs their tests. The CMake option ``MFEM_MGIS_EXAMPLES_TEST_MODE``
selects their size:

- ``full`` runs the complete simulations.
- ``restricted`` only runs their beginning or a smaller case. It is much
  faster.
- ``auto``, the default, selects ``restricted`` in the ``Debug`` and
  ``Coverage`` builds and ``full`` otherwise.

Some tests are the same in both modes.

The CMake option ``MFEM_MGIS_EXAMPLES_GPU_TESTS`` (``OFF`` by default) also
registers the tests of the examples on a GPU (option ``--execution gpu``, or
``--execution mixed`` when a phase stays on the CPU), when MFEM is compiled
with CUDA or HIP and the behaviour libraries of the example are built for the
GPU (``mfront`` must find the configuration file ``config-cuda-nvcc.json``,
resp. ``config-hip-hipcc.json``, of TFEL; CMake prints why they are not
built). A GPU is then required to run them.

The options and the tests described in these pages (``--execution``, behaviour
libraries for the CPU and the GPU, tests on a GPU) are the ones of the
prototype branch of the examples, which requires the prototype branch of
mfem-mgis, built from its sources.

.. toctree::
   :maxdepth: 1

   basic_tests
   rve_elastic_inclusions
   rve_mox
   rjh_plate

mm-opera-hpc
============

These examples are in https://github.com/rprat-pro/mm-opera-hpc. They were
developed in the OperaHPC project. As for mfem-mgis-examples, the options and
the tests described in these pages are the ones of the prototype branch of
mm-opera-hpc, which requires the prototype branch of mfem-mgis.

.. toctree::
   :maxdepth: 1

   bubbles
   polycrystal
   cermet
