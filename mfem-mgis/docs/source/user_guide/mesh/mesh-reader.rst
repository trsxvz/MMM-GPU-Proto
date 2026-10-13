===========
Mesh Reader
===========

``MFEM/MGIS`` uses the classic sequential reading of meshes generated
via gmsh by specifying the `MFEM` input format with the ``-format
msh22`` option.

.. note::
  For more information: https://mfem.org/mesh-formats/

Keywords
^^^^^^^^

+--------------------------------+-------------------------------------------------------------------+
| Key                            | Description                                                       |
+================================+===================================================================+
| ``MeshFileName``               | Define the mesh file name.                                        |
+--------------------------------+-------------------------------------------------------------------+
| ``MeshReadMode``               | ``FromScratch`` (default): the whole mesh is read, then           |
|                                | partitioned; ``Restart``: each process reads its part of a mesh   |
|                                | split beforehand (``<MeshFileName><rank on 6 digits>``).          |
|                                | ``Restart`` is refused in a sequential computation and with a     |
|                                | ``.med`` file.                                                    |
+--------------------------------+-------------------------------------------------------------------+
| ``NumberOfUniformRefinements`` | Number of uniform refinements of the mesh (0 by default).         |
+--------------------------------+-------------------------------------------------------------------+
| ``Materials``                  | Names associated with the material identifiers.                   |
+--------------------------------+-------------------------------------------------------------------+
| ``Boundaries``                 | Names associated with the boundary identifiers.                   |
+--------------------------------+-------------------------------------------------------------------+
| ``Points``, ``PointsSets``     | Points and sets of points declared by the user (see               |
|                                | :doc:`curves`).                                                   |
+--------------------------------+-------------------------------------------------------------------+
| ``GeneralVerbosityLevel``      | Verbosity level.                                                  |
+--------------------------------+-------------------------------------------------------------------+
 
.. code-block::

   mfem_mgis::Parameters{{"MeshFileName", "mesh-explorer.mesh."},
                         {"MeshReadMode", "Restart"}, ... }

MED files
^^^^^^^^^

When ``MFEM`` is built with the support of the ``MED`` format
(``MFEM_USE_MED``), the meshes in ``.med`` files are read directly. If a
file with the same name and the extension ``.per`` exists, it describes
the periodicity of the mesh.

Split Mesh
^^^^^^^^^^

If your mesh is too large (memory limit), you can also use MFEM's option
of splitting the mesh into smaller meshes and reading them in parallel.

To split the mesh, we use the ``mesh-explorer`` tool
(https://mfem.org/meshing-miniapps/#mesh-explorer) and you can specify
the parallel reader using the keyword: ``MeshReadMode``.


Install Mesh-explorer
---------------------

Here are two ways to install mesh-explorer. We recommend installing via
spack, but we will also add the installation procedure with cmake.

Spack version:

.. code-block::

  spack install mfem+miniapps  

The binary file is then located at: 

.. code-block::

  `spack location -i mfem`/share/mfem/miniapps/meshing/mesh-explorer

.. note::

  To make sure the mfem used by mfem-mgis provides the miniapps, install mfem-mgis with ``spack install mfem-mgis@master ^mfem+miniapps``.


CMake version:

.. code-block::

   git clone https://github.com/mfem/mfem.git
   mkdir build && cd build
   spack load hypre
   cmake ../mfem/ -DMFEM_USE_MPI=ON -DMFEM_ENABLE_MINIAPPS=ON
   make -j 8

The binary file is then located at: 

.. code-block::

  build/miniapps/meshing/mesh-explorer


Mesh generation
---------------

The mesh explorer works interactively. Here is a series of options for splitting a gmsh mesh into `128` meshes named `output-mesh.XXXXXX`

Run:

.. code-block::

  mesh-explorer --mesh mesh.msh


Create partitioning:

.. code-block::

  press p > to generate a partitioning
  press 1 > choose the partitioning method
  press 128 > choose the number of MPI processes

Save mesh in parallel format:

.. code-block::

  press T > choose to save in parallel format
  press output-mesh. > prefix name for output files
  press 6 > number of digits


Load your mesh with MFEM-MGIS
-----------------------------

.. note::

  Please, verify that the split mesh chosen is split into 128 files if you want to run it on 128 mpi processes.

Example: 

- One example is available here: tests/ParallelReadMode.cxx

.. code-block:: c++

    auto ctx = mfem_mgis::Context{};
    auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      ctx, mfem_mgis::Parameters{{"MeshFileName", "output-mesh."},
                            {"MeshReadMode", "Restart"},
                            {"FiniteElementFamily", "H1"},
                            {"FiniteElementOrder", p.order},
                            {"UnknownsSize", dim},
                            {"NumberOfUniformRefinements", p.parallel ? 2 : 0},
                            {"Parallel", p.parallel}});
